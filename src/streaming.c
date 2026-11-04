/* streaming.c: each read's scores, made on any thread (ri_stream_prepare_read) and applied in read order
 * (ri_stream_commit_read); see streaming.h. */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif
#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 700
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "streaming.h"
#include "khash.h"
#include "onsprof.h"
#include "calls.h"
#include "slope_range.h"

struct ri_stream_state_s {
	ri_stream_config_t cfg;          /* its grid is the current slope grid */
	uint64_t total_reads;            /* reads committed so far */
	ri_slopes_t *slopes;             /* NULL unless --stream-slope-chains */
	uint32_t *tax_id;                /* per signature, its taxon as a number (UINT32_MAX: "*"), and the taxa */
	uint32_t  n_tax;
};

/* One thread's buffers for ri_stream_prepare_read, kept from read to read. */
struct ri_stream_scratch_s {
	ri_calls_scratch_t *calls;       /* the read p-value's and the near-neighbor check's */
	ri_band_buf_t       band;        /* the band counters */
	ri_slopes_buf_t     slopes;      /* the slope refits */
	uint32_t           *tax_slot;    /* per taxon, its place in the read's per-taxon table, valid when tax_read is the */
	uint32_t           *tax_read;    /* read's number in this scratch */
	uint32_t            n_read;
};

/* One taxon's best signatures in a read */
typedef struct {
	const char *taxon;
	double      top1;           /* best signature score (count * containment) */
	double      best_count;     /* largest count among the taxon's scored signatures */
	uint32_t    best_rid;       /* the signature with that count */
	uint32_t    n_sigs;         /* the taxon's signatures with a score */
	/* --stream-nonoverlap-gap */
	uint32_t    p_best;         /* the largest non-overlapping hit count over the taxon's signatures (ties: more
	                               hits) */
	uint32_t    p_x;            /* hits in that band */
	uint32_t    p_second;       /* the largest non-overlapping hit count among the taxon's other signatures */
	double      q_p;            /* the smallest band p-value (ri_band_pvalue) over the taxon's signatures */
	uint32_t    q_rid;          /* that signature, its non-overlapping hit count, its non-overlapping hits and its
	                               length */
	uint32_t    q_nonoverlap;
	uint32_t    q_units;
	uint32_t    q_len;
} per_taxon_t;

/* One read's scores, which depend only on the read and the slope grid of its minibatch (ri_stream_prepare_read),
 * applied to the stream state in read order by ri_stream_commit_read. */
struct ri_read_scores_s {
	per_taxon_t    *per_tax;     /* NULL: no taxon scored */
	uint32_t        n_per, read_events, n_slopes;
	int             scored;      /* the read p-value ran (cfg.calls and --stream-nonoverlap-gap) */
	ri_calls_set_t *sets;        /* the read's signature sets for calls.c, one per scored taxon */
	int             call_e, call_set;
	double          call_p;
	int             checked;     /* near-neighbor check: -3 skipped (target already called), -2 not run, -1 no band,
	                              * 0 dropped, 1 kept */
	double          near_p;      /* the near-neighbor p-value when checked is 0 or 1 */
	ri_slopes_read_t slope;      /* the sliding-window slope range's refits */
};

static void *xalloc(void *p)
{
	if (!p) { fprintf(stderr, "[ERROR] stream layer: out of memory\n"); exit(1); }
	return p;
}

KHASH_MAP_INIT_STR(taxid, uint32_t)

void ri_stream_config_init(ri_stream_config_t *cfg)
{
	if (!cfg) return;
	memset(cfg, 0, sizeof(*cfg));
	cfg->emit_fd                  = -1;
	cfg->emit_fp                  = NULL;
	cfg->grid.min                 = 1.0f;  /* slope 1 only */
	cfg->grid.max                 = 1.0f;
	cfg->grid.step                = 0.0f;
	cfg->nonoverlap_gap           = 0;
	cfg->slope_chains             = 0;     /* sliding-window slope range off */
	cfg->slope_alpha              = 1e-7;
	cfg->slope_width              = 10.0f;
	cfg->slope_width_sd           = 0;
	cfg->slope_cover              = 1.0f;
}

ri_stream_state_t *ri_stream_init(const ri_stream_config_t *cfg, const ri_idx_t *ri)
{
	ri_stream_state_t *s;
	if (!cfg || !ri) return NULL;
	s = (ri_stream_state_t*)calloc(1, sizeof(*s));
	if (!s) return NULL;
	s->cfg = *cfg;
	s->total_reads = 0;

	/* Lazily open emit_fd if requested. */
	if (s->cfg.emit_fd >= 0) {
		s->cfg.emit_fp = fdopen(s->cfg.emit_fd, "w");
		if (s->cfg.emit_fp) setvbuf(s->cfg.emit_fp, NULL, _IOLBF, 0);
	}

	{   /* each signature's taxon as a number, so that a read finds its taxa without comparing names */
		khash_t(taxid) *h = kh_init(taxid);
		uint32_t rid;
		int absent;
		s->tax_id = (uint32_t*)xalloc(malloc(sizeof(uint32_t) * (ri->n_seq ? ri->n_seq : 1)));
		for (rid = 0; rid < ri->n_seq; ++rid) {
			const char *t = ri_sig_taxon(ri, rid);
			khint_t k;
			if (!t || strcmp(t, "*") == 0) { s->tax_id[rid] = UINT32_MAX; continue; }
			k = kh_put(taxid, h, t, &absent);
			if (absent) kh_val(h, k) = s->n_tax++;
			s->tax_id[rid] = kh_val(h, k);
		}
		kh_destroy(taxid, h);
	}
	if (s->cfg.slope_chains > 0)
		s->slopes = ri_slopes_init(&s->cfg.grid, s->cfg.slope_chains, s->cfg.slope_alpha, s->cfg.slope_width,
		                           s->cfg.slope_width_sd, s->cfg.slope_cover, s->cfg.nonoverlap_gap, s->cfg.calls, s->cfg.emit_fp);
	return s;
}

void ri_stream_destroy(ri_stream_state_t *s)
{
	if (!s) return;
	ri_slopes_destroy(s->slopes);
	free(s->tax_id);
	if (s->cfg.emit_fp) {
		fflush(s->cfg.emit_fp);
		/* Don't fclose: the FD is owned by the caller. */
	}
	free(s);
}

/* containment(rid) = distinct_query_hashes / sig_emitted_positions[rid]. */
static double sig_containment(const ri_idx_t *ri, const ri_sig_hits_t *s)
{
	uint32_t denom;
	if (!ri || !s) return 0.0;
	if (!ri->sig_emitted_positions || s->rid >= ri->n_seq) return 0.0;
	denom = ri->sig_emitted_positions[s->rid];
	if (denom == 0) return 0.0;
	return (double)s->distinct_query_hashes / (double)denom;
}

/* The read's record on fd 3, one JSON line per read with a scored taxon. */
static void emit_read(ri_stream_state_t *s, uint64_t read_id, const per_taxon_t *per_tax, uint32_t n_per,
                          uint32_t read_events)
{
	uint32_t i;
	if (!s || !s->cfg.emit_fp) return;
	fprintf(s->cfg.emit_fp, "{\"kind\":\"read\",\"read_id\":%llu", (unsigned long long)read_id);
	{   /* per taxon with a score: [taxon, best signature score (count * containment), best count, its rid, signatures scored] */
		fputs(",\"per_taxon\":[", s->cfg.emit_fp);
		for (i = 0; per_tax && i < n_per; ++i)
			fprintf(s->cfg.emit_fp, "%s[\"%s\",%.6f,%.0f,%u,%u]", i ? "," : "", per_tax[i].taxon,
			        per_tax[i].top1, per_tax[i].best_count, per_tax[i].best_rid, per_tax[i].n_sigs);
		fputc(']', s->cfg.emit_fp);
	}
	/* --stream-nonoverlap-gap: the read's length in events (largest hit position + 1); per taxon: [taxon, largest
	 * non-overlapping hit count, hits in its band, second largest, then the signature with the smallest band p-value:
	 * rid, its non-overlapping hit count, its non-overlapping hits over all bands, its length] */
	if (s->cfg.nonoverlap_gap > 0) {
		fprintf(s->cfg.emit_fp, ",\"read_events\":%u,\"nonoverlap\":[", read_events);
		for (i = 0; per_tax && i < n_per; ++i)
			fprintf(s->cfg.emit_fp, "%s[\"%s\",%u,%u,%u,%u,%u,%u,%u]", i ? "," : "", per_tax[i].taxon,
			        per_tax[i].p_best, per_tax[i].p_x, per_tax[i].p_second, per_tax[i].q_rid, per_tax[i].q_nonoverlap,
			        per_tax[i].q_units, per_tax[i].q_len);
		fputc(']', s->cfg.emit_fp);
	}
	if (s->slopes) ri_slopes_emit(s->slopes, s->cfg.emit_fp);
	fputs("}\n", s->cfg.emit_fp);
}

/* The read's signature sets as calls.c takes them: the nonoverlap fields of its read event. */
static void read_sets(ri_calls_set_t *sets, const per_taxon_t *per_tax, uint32_t n_per)
{
	uint32_t i;
	for (i = 0; i < n_per; ++i) {
		sets[i].taxon = per_tax[i].taxon;
		sets[i].nonoverlap_max   = per_tax[i].p_best;
		sets[i].rid   = per_tax[i].q_rid;
		sets[i].nonoverlap   = per_tax[i].q_nonoverlap;
		sets[i].units = per_tax[i].q_units;
		sets[i].len   = per_tax[i].q_len;
	}
}

ri_stream_scratch_t *ri_stream_scratch_init(const ri_stream_state_t *s)
{
	ri_stream_scratch_t *x = (ri_stream_scratch_t*)xalloc(calloc(1, sizeof(ri_stream_scratch_t)));
	if (s && s->cfg.calls) x->calls = (ri_calls_scratch_t*)xalloc(ri_calls_scratch_init(s->cfg.calls));
	if (s && s->slopes) ri_slopes_buf_reserve(&x->slopes, 1024);   /* the default hit cap; grows if a signature has more */
	if (s) {
		x->tax_slot = (uint32_t*)xalloc(calloc(s->n_tax + 1, sizeof(uint32_t)));
		x->tax_read = (uint32_t*)xalloc(calloc(s->n_tax + 1, sizeof(uint32_t)));
	}
	return x;
}

void ri_stream_scratch_destroy(ri_stream_state_t *s, ri_stream_scratch_t *x)
{
	if (!x) return;
	if (x->calls) ri_calls_scratch_destroy(s ? s->cfg.calls : NULL, x->calls);
	ri_band_buf_free(&x->band);
	ri_slopes_buf_free(&x->slopes);
	free(x->tax_slot); free(x->tax_read);
	free(x);
}

ri_read_scores_t *ri_stream_prepare_read(const ri_stream_state_t *s, const ri_idx_t *ri, const ri_read_t *rd,
                                         ri_stream_scratch_t *x)
{
	ri_read_scores_t *sc;
	per_taxon_t *per_tax;
	uint32_t m, n_per = 0, read_events = 0;

	if (!s || !ri || !rd) return NULL;
	sc = (ri_read_scores_t*)xalloc(calloc(1, sizeof(ri_read_scores_t)));
	sc->call_e = sc->call_set = -1;
	sc->call_p = 1.0;
	sc->checked = -2;
	sc->near_p = 1.0;
	if (rd->n_sig_hits == 0) return sc;

	#ifdef ONS_PROF
	ons_prof_read_t *pr = ons_prof_cur;
	uint64_t pt0 = pr ? ons_now() : 0, pt1;
	#endif

	/* The read's length in events for the band p-values: its largest hit position + 1. */
	if (s->cfg.nonoverlap_gap > 0) {
		uint32_t i;
		for (m = 0; m < rd->n_sig_hits; ++m)
			for (i = 0; i < rd->sig_hits[m].n_hits; ++i)
				if (rd->sig_hits[m].hits[i].q_pos + 1 > read_events)
					read_events = rd->sig_hits[m].hits[i].q_pos + 1;
	}
	sc->read_events = read_events;

	/* Per taxon: its best signature score (count * containment) and count, and with --stream-nonoverlap-gap the
	 * non-overlapping hit counts and band p-values of its signatures. */
	per_tax = (per_taxon_t*)xalloc(calloc(s->n_tax + 1, sizeof(per_taxon_t)));
	++x->n_read;
	ri_band_pv_read(&x->band, read_events);
	for (m = 0; m < rd->n_sig_hits; ++m) {
		const ri_sig_hits_t *sh = &rd->sig_hits[m];
		uint32_t tid = sh->rid < ri->n_seq ? s->tax_id[sh->rid] : UINT32_MAX;
		double containment, count, primary;
		uint32_t d_p = 0, d_x = 0, d_units = 0, d_xmax = 0;
		int    idx;
		if (tid == UINT32_MAX) continue;
		/* The signature's count: the fullest diagonal band. With --stream-nonoverlap-gap it comes from the same
		 * sorted pass as the non-overlapping hit count. */
		if (s->cfg.nonoverlap_gap > 0)
			d_p = ri_band_nonoverlap(sh, &s->cfg.grid, s->cfg.nonoverlap_gap, &x->band, &d_x, &d_units, &d_xmax);
		count = s->cfg.nonoverlap_gap > 0 ? (double)d_xmax : (double)ri_band_fullest(sh, &s->cfg.grid);
		if (count <= 0.0) continue;
		containment = sig_containment(ri, sh);
		if (containment <= 0.0) continue;
		primary = count * containment;

		idx = x->tax_read[tid] == x->n_read ? (int)x->tax_slot[tid] : -1;
		if (idx < 0) {
			x->tax_read[tid] = x->n_read;
			x->tax_slot[tid] = n_per;
			per_tax[n_per].taxon = ri_sig_taxon(ri, sh->rid);
			per_tax[n_per].top1  = primary;
			per_tax[n_per].best_count = count;
			per_tax[n_per].best_rid   = sh->rid;
			per_tax[n_per].n_sigs     = 1;
			++n_per;
		} else {
			++per_tax[idx].n_sigs;
			if (count > per_tax[idx].best_count) {
				per_tax[idx].best_count = count;
				per_tax[idx].best_rid   = sh->rid;
			}
			if (primary > per_tax[idx].top1) per_tax[idx].top1 = primary;
		}
		if (s->cfg.nonoverlap_gap > 0) {
			per_taxon_t *t = &per_tax[idx < 0 ? (int)n_per - 1 : idx];
			uint32_t hx = d_x, units = d_units, len, p = d_p;
			double   pv;
			if (idx < 0) t->q_p = 2.0;
			if (p > t->p_best || (p == t->p_best && hx > t->p_x)) {
				if (t->p_best > t->p_second) t->p_second = t->p_best;
				t->p_best = p;
				t->p_x    = hx;
			} else if (p > t->p_second) {
				t->p_second = p;
			}
			len = ri->seq[sh->rid].len;
			pv = ri_band_pvalue_buf(&x->band, p, units, len);
			if (pv < t->q_p || (pv == t->q_p && p > t->q_nonoverlap)) {
				t->q_p     = pv;
				t->q_rid   = sh->rid;
				t->q_nonoverlap   = p;
				t->q_units = units;
				t->q_len   = len;
			}
		}
	}

	#ifdef ONS_PROF
	if (pr) { pt1 = ons_now(); pr->t_bands += pt1 - pt0; pt0 = pt1; pr->n_taxa = n_per; }
	#endif

	if (n_per == 0) {
		free(per_tax);
		return sc;
	}
	{   /* the array was sized for every taxon of the index; the read keeps it until its commit */
		per_taxon_t *shrunk = (per_taxon_t*)realloc(per_tax, sizeof(per_taxon_t) * n_per);
		if (shrunk) per_tax = shrunk;
	}
	sc->per_tax = per_tax;
	sc->n_per = n_per;

	/* The read's p-value and call (calls.c), for the read and target calls, the near-neighbor check and the
	 * sliding-window slope range, which takes the chain of a confident call. The grid changes only between minibatches,
	 * so it is the one the read ran on. */
	if (s->cfg.calls && s->cfg.nonoverlap_gap > 0) {
		sc->scored = 1;
		sc->n_slopes = ri_grid_n_slopes(&s->cfg.grid);
		sc->sets = (ri_calls_set_t*)xalloc(malloc(sizeof(ri_calls_set_t) * n_per));
		read_sets(sc->sets, per_tax, n_per);
		sc->call_e = ri_read_pvalue(s->cfg.calls, x->calls, sc->sets, n_per, read_events, sc->n_slopes, &sc->call_p,
		                            &sc->call_set);
		#ifdef ONS_PROF
		if (pr) { pt1 = ons_now(); pr->t_pvalue += pt1 - pt0; pt0 = pt1; pr->call = sc->call_e >= 0; }
		#endif
		if (s->slopes)
			ri_slopes_prepare(s->slopes, &x->slopes, x->calls, ri, rd, sc->call_e,
			                  sc->call_e >= 0 ? sc->sets[sc->call_set].rid : 0, read_events, &sc->slope);
		#ifdef ONS_PROF
		if (pr) { pt1 = ons_now(); pr->t_slopes += pt1 - pt0; pt0 = pt1; }
		#endif
		/* the near-neighbor check of a read call at the top cutoff on its signature's band, skipped once the read's
		 * target is called (ri_calls_done) */
		if (sc->call_e >= 0 && ri_calls_neighbor_wanted(s->cfg.calls, sc->call_p) && ri_calls_done(s->cfg.calls, sc->call_e))
			sc->checked = -3;
		else if (sc->call_e >= 0 && ri_calls_neighbor_wanted(s->cfg.calls, sc->call_p)) {
			uint32_t rid = sc->sets[sc->call_set].rid;
			nc_band_t band;
			for (m = 0; m < rd->n_sig_hits && rd->sig_hits[m].rid != rid; ++m) {}
			sc->checked = -1;
			if (m < rd->n_sig_hits && rd->sig_hits[m].n_hits > 0 &&
			    ri_band_box(&rd->sig_hits[m], &s->cfg.grid, s->cfg.nonoverlap_gap, &x->band, &band)) {
				double q = ri_calls_neighbor_check(s->cfg.calls, x->calls, rid, rd->levels, rd->n_levels, &band,
				                                   sc->call_p);
				if (q >= 0) {    /* kept: the read still counts at the loosest cutoff */
					sc->near_p = q;
					sc->checked = ri_calls_neighbor_wanted(s->cfg.calls, q > sc->call_p ? q : sc->call_p);
				}
			}
		}
		#ifdef ONS_PROF
		if (pr) { pr->t_neighbor += ons_now() - pt0; pr->checked = sc->checked; }
		#endif
	}
	return sc;
}

void ri_stream_commit_read(ri_stream_state_t *s, const ri_read_t *rd, ri_read_scores_t *sc, uint64_t read_id)
{
	if (!s || !rd || !sc) goto done;
	if (sc->per_tax) {
		/* The sliding-window slope range takes the chain of a confident call. */
		if (sc->scored && s->slopes) ri_slopes_commit(s->slopes, sc->call_e, sc->call_p, &sc->slope);
		/* The read's record, then its row and the target-call tests. */
		emit_read(s, read_id, sc->per_tax, sc->n_per, sc->read_events);
		if (sc->scored)
			ri_calls_record(s->cfg.calls, read_id, rd->read_name, sc->call_e >= 0 ? &sc->sets[sc->call_set] : NULL,
			                sc->call_e, sc->call_p, sc->read_events, sc->n_slopes, sc->checked, sc->near_p);
	}
	++s->total_reads;

done:
	if (sc) {
		free(sc->per_tax);
		free(sc->sets);
		free(sc);
	}
}

void ri_stream_end_batch(ri_stream_state_t *s)
{
	if (s && s->slopes) ri_slopes_end_batch(s->slopes, &s->cfg.grid, s->total_reads);
}
