/* slope_range.c: the sliding-window slope range; see slope_range.h. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "slope_range.h"
#include "utils.h"

#define SLOPE_WIDTH_SDS 4.0   /* --stream-slope-width-sd: W = 4 SD, as 10 events = 4 x 2.54 on R9 */

struct ri_slopes_s {
	uint32_t     n, head, held;  /* ring size, next slot, chains held */
	ri_chain_t  *ring;
	double      *tmp;            /* n values, for the medians */
	uint64_t     confident;      /* chains added */
	uint32_t     changes;        /* grid changes */
	uint32_t     n_slopes;       /* slopes of the current grid */
	ri_slope_grid_t search;      /* the grid at the start, on which every read's chain is judged */
	uint32_t     n_search;
	double       alpha;          /* settings: --stream-slope-alpha, -width, -width-sd and --stream-nonoverlap-gap */
	float        width;
	int          width_sd;
	float        cover;          /* --slope-grid */
	uint32_t     gap;
	const ri_calls_t *calls;     /* the read p-value of a chain */
	FILE        *emit_fp;        /* the slopes lines of fd 3; NULL: none */
	int          read_entry;     /* this read's call by read p-value (-1: none) and its p */
	double       read_p;
	int          read_chain;     /* 1 when the read's chain went into the ring: its entry and the chain */
	int          chain_entry;
	ri_chain_t   chain;
	uint64_t     n_called, n_other_sig;  /* reads with a call; of them, the most-hits signature was not the call's */
	uint64_t     n_from_most;    /* chains that only the most-hits signature gave */
};

static void *xalloc(void *p)
{
	if (!p) { fprintf(stderr, "[ERROR] sliding-window slope range: out of memory\n"); exit(1); }
	return p;
}

void ri_slopes_buf_reserve(ri_slopes_buf_t *x, uint32_t n)
{
	if (n <= x->cap_n) return;
	x->cap_n = n;
	x->cap_h = 16;
	while (x->cap_h < n * 2u) x->cap_h <<= 1;
	x->a     = (uint64_t*)xalloc(realloc(x->a, sizeof(uint64_t) * n));
	x->qq    = (double*)xalloc(realloc(x->qq, sizeof(double) * n));
	x->rr    = (double*)xalloc(realloc(x->rr, sizeof(double) * n));
	x->keep  = (uint8_t*)xalloc(realloc(x->keep, n));
	x->keys  = (int64_t*)xalloc(realloc(x->keys, sizeof(int64_t) * x->cap_h));
	x->dec   = (uint32_t*)xalloc(realloc(x->dec, sizeof(uint32_t) * x->cap_h));
	x->lastq = (uint32_t*)xalloc(realloc(x->lastq, sizeof(uint32_t) * x->cap_h));
}

void ri_slopes_buf_free(ri_slopes_buf_t *x)
{
	free(x->a); free(x->keys); free(x->dec); free(x->lastq); free(x->qq); free(x->rr); free(x->keep);
}

/* The signature's hits in x->a, sorted by read position; returns its non-overlapping hits over all bands (the n of
 * the read p-value, as ri_band_nonoverlap counts them). */
static uint32_t sorted_hits(ri_slopes_buf_t *x, const ri_sig_hits_t *s, uint32_t gap)
{
	uint32_t i, n = s->n_hits, u = 1;
	uint64_t last;
	ri_slopes_buf_reserve(x, n);
	for (i = 0; i < n; ++i) x->a[i] = ((uint64_t)s->hits[i].q_pos << 32) | s->hits[i].ref_pos;
	radix_sort_64(x->a, x->a + n);   /* in place; glibc's qsort allocates above 1 KB */
	for (i = 1, last = x->a[0] >> 32; i < n; ++i)
		if ((x->a[i] >> 32) - last >= gap) { ++u; last = x->a[i] >> 32; }
	return u;
}

/* Least-squares line q = c + b r through the points with keep[i] set (all when keep is NULL); returns the SD of the
 * residuals (numpy's std: mean removed, divided by the count). */
static double line_fit(const double *r, const double *q, const uint8_t *keep, uint32_t n, double *b, double *c)
{
	double mr = 0.0, mq = 0.0, sxx = 0.0, sxy = 0.0, me = 0.0, v = 0.0;
	uint32_t i, m = 0;
	for (i = 0; i < n; ++i) if (!keep || keep[i]) { mr += r[i]; mq += q[i]; ++m; }
	mr /= m; mq /= m;
	for (i = 0; i < n; ++i) if (!keep || keep[i]) { sxx += (r[i] - mr) * (r[i] - mr); sxy += (r[i] - mr) * (q[i] - mq); }
	*b = sxy / sxx;
	*c = mq - *b * mr;
	for (i = 0; i < n; ++i) if (!keep || keep[i]) me += q[i] - (*c + *b * r[i]);
	me /= m;
	for (i = 0; i < n; ++i) if (!keep || keep[i]) { double e = q[i] - (*c + *b * r[i]) - me; v += e * e; }
	return sqrt(v / m);
}

/* The chain of one signature: the band with the most non-overlapping hits over
 * the grid, ties to the first slope and then the lowest band (*n_nonoverlap is that
 * count, *units the hits over all bands);
 * least squares through the hits in it and its two neighbours, then once more without the hits off the line by more
 * than max(3 SD, 1 event) if at least 3 hits with an r-span stay. Returns 0 with fewer than 3 hits or no r-span. */
static int chain_fit(ri_slopes_buf_t *x, const ri_sig_hits_t *s, const ri_slope_grid_t *g, uint32_t gap,
                     ri_chain_t *ch, uint32_t *n_nonoverlap, uint32_t *units)
{
	uint32_t i, j, k, n = s->n_hits, m, h = 16, n_slopes = ri_grid_n_slopes(g), best = 0, best_k = 0, kept = 0;
	int64_t  best_key = 0;
	double   lo = g->min > 0.0f ? (double)g->min : 1.0, step = (double)g->step;
	double   b, c, sd, lim, rmin, rmax, kmin = 0.0, kmax = 0.0;
	int      unit;

	if (n < 3) return 0;
	*units = sorted_hits(x, s, gap);
	while (h < n * 2u) h <<= 1;   /* hash slots for this signature, within the reserved cap_h */
	for (k = 0; k < n_slopes; ++k) {
		uint32_t top = 0;
		int64_t  top_key = 0;
		b = lo + step * (double)k;
		unit = fabs(b - 1.0) < 1e-6;
		for (j = 0; j < h; ++j) x->keys[j] = INT64_MIN;
		for (i = 0; i < n; ++i) {   /* per band: non-overlapping hits so far and the last counted position */
			int64_t  q = (int64_t)(x->a[i] >> 32), rp = (int64_t)(x->a[i] & 0xffffffffULL);
			int64_t  key = unit ? (q - rp) / 10 : (int64_t)(((double)q - b * (double)rp) / 10.0);
			uint32_t idx = (uint32_t)((uint64_t)(key * 0x9E3779B97F4A7C15ULL) & (h - 1));
			while (x->keys[idx] != INT64_MIN && x->keys[idx] != key) idx = (idx + 1) & (h - 1);
			if (x->keys[idx] == INT64_MIN) { x->keys[idx] = key; x->dec[idx] = 1; x->lastq[idx] = (uint32_t)q; }
			else if ((uint32_t)q - x->lastq[idx] >= gap) { ++x->dec[idx]; x->lastq[idx] = (uint32_t)q; }
		}
		for (j = 0; j < h; ++j)
			if (x->keys[j] != INT64_MIN && (x->dec[j] > top || (x->dec[j] == top && x->keys[j] < top_key))) {
				top = x->dec[j];
				top_key = x->keys[j];
			}
		if (top > best) { best = top; best_k = k; best_key = top_key; }
	}

	b = lo + step * (double)best_k;
	unit = fabs(b - 1.0) < 1e-6;
	for (i = 0, m = 0; i < n; ++i) {
		int64_t q = s->hits[i].q_pos, rp = s->hits[i].ref_pos;
		int64_t key = unit ? (q - rp) / 10 : (int64_t)(((double)q - b * (double)rp) / 10.0);
		if (key >= best_key - 1 && key <= best_key + 1) { x->qq[m] = (double)q; x->rr[m] = (double)rp; ++m; }
	}
	if (m < 3) return 0;
	rmin = rmax = x->rr[0];
	for (i = 1; i < m; ++i) { if (x->rr[i] < rmin) rmin = x->rr[i]; if (x->rr[i] > rmax) rmax = x->rr[i]; }
	if (!(rmax > rmin)) return 0;
	sd = line_fit(x->rr, x->qq, NULL, m, &b, &c);
	lim = 3.0 * sd > 1.0 ? 3.0 * sd : 1.0;
	for (i = 0; i < m; ++i) {
		x->keep[i] = fabs(x->qq[i] - (c + b * x->rr[i])) <= lim;
		if (!x->keep[i]) continue;
		if (kept == 0 || x->rr[i] < kmin) kmin = x->rr[i];
		if (kept == 0 || x->rr[i] > kmax) kmax = x->rr[i];
		++kept;
	}
	if (kept >= 3 && kmax > kmin) sd = line_fit(x->rr, x->qq, x->keep, m, &b, &c);
	ch->slope = b;
	ch->span  = rmax - rmin;
	ch->sd    = sd;
	*n_nonoverlap = best;
	return 1;
}

ri_slopes_t *ri_slopes_init(const ri_slope_grid_t *search, uint32_t n_chains, double alpha, float width, int width_sd,
                            float cover, uint32_t gap, const ri_calls_t *calls, FILE *emit_fp)
{
	ri_slopes_t *r = (ri_slopes_t*)xalloc(calloc(1, sizeof(ri_slopes_t)));
	r->n = n_chains;
	r->ring = (ri_chain_t*)xalloc(calloc(r->n, sizeof(ri_chain_t)));
	r->tmp = (double*)xalloc(malloc(r->n * sizeof(double)));
	r->n_slopes = r->n_search = ri_grid_n_slopes(search);
	r->search = *search;
	r->alpha = alpha;
	r->width = width;
	r->width_sd = width_sd;
	r->cover = cover;
	r->gap = gap;
	r->calls = calls;
	r->emit_fp = emit_fp;
	r->read_entry = -1;
	fprintf(stderr, "[M::%s] slope grid from the last %u chains with read p-value <= %g on the search grid %.6f:%.6f:%.6f "
	        "(%u slopes); W %s\n", __func__, r->n, alpha, search->min, search->max, search->step, r->n_slopes,
	        width_sd ? "4 x median residual SD" : "--stream-slope-width");
	if (emit_fp)
		fprintf(emit_fp, "{\"kind\":\"slopes\",\"total_reads\":0,\"grid\":\"%.6f:%.6f:%.6f\",\"slopes\":%u,"
		        "\"chains\":0,\"confident\":0}\n", search->min, search->max, search->step, r->n_slopes);
	return r;
}

void ri_slopes_destroy(ri_slopes_t *r)
{
	if (!r) return;
	fprintf(stderr, "[M::%s] slope grid changed %u times; %llu chains, %llu of them only from the most-hits signature; "
	        "of %llu reads with a call, the most-hits signature was another in %llu\n", __func__, r->changes,
	        (unsigned long long)r->confident, (unsigned long long)r->n_from_most, (unsigned long long)r->n_called,
	        (unsigned long long)r->n_other_sig);
	free(r->ring); free(r->tmp);
	free(r);
}

/* One signature refitted on the search grid, with its read p-value there (calls.c, on that signature alone; 2 when it
 * has no chain) and its entry. */
static ri_chain_t slopes_try(const ri_slopes_t *r, ri_slopes_buf_t *x, ri_calls_scratch_t *cx, const ri_idx_t *ri,
                             const ri_sig_hits_t *sig, uint32_t read_events, int *e)
{
	ri_chain_t ch;
	ri_calls_set_t one;
	uint32_t n_nonoverlap, units;
	*e = -1;
	if (!chain_fit(x, sig, &r->search, r->gap, &ch, &n_nonoverlap, &units)) { ch.p = 2.0; return ch; }
	one.taxon = ri_sig_taxon(ri, sig->rid);
	one.nonoverlap_max = one.nonoverlap = n_nonoverlap;
	one.rid = sig->rid;
	one.units = units;
	one.len = ri->seq[sig->rid].len;
	*e = ri_read_pvalue(r->calls, cx, &one, 1, read_events, r->n_search, &ch.p, NULL);
	return ch;
}

/* Per read, with its call by read p-value (call_e, the signature call_rid it came from). Two signatures are
 * refitted on the search grid: the call's (the pre-pass's choice) and the one with the most non-overlapping hits,
 * which does not depend on the grid, for reads whose band lies outside it. ri_slopes_commit puts the one with the
 * smaller read p-value there into the ring when that p is at most slope_alpha. */
void ri_slopes_prepare(const ri_slopes_t *r, ri_slopes_buf_t *x, ri_calls_scratch_t *cx, const ri_idx_t *ri,
                       const ri_read_t *rd, int call_e, uint32_t call_rid, uint32_t read_events, ri_slopes_read_t *out)
{
	uint32_t m, u, best_u, most = 0, call = UINT32_MAX, cand = UINT32_MAX;
	uint32_t rid = call_e >= 0 ? call_rid : 0;
	out->c1.p = out->c2.p = 2.0;
	out->e1 = out->e2 = -1;
	if (call_e >= 0)
		for (m = 0; m < rd->n_sig_hits; ++m)
			if (rd->sig_hits[m].rid == rid) { call = m; break; }

	/* The signature with the most non-overlapping hits: start with the one with the most hits, then try only
	 * those with more hits than the best count so far (the count is at most the hits). */
	for (m = 0; m < rd->n_sig_hits; ++m)
		if (rd->sig_hits[m].n_hits > most && strcmp(ri_sig_taxon(ri, rd->sig_hits[m].rid), "*") != 0) {
			most = rd->sig_hits[m].n_hits;
			cand = m;
		}
	if (cand == UINT32_MAX) return;
	best_u = sorted_hits(x, &rd->sig_hits[cand], r->gap);
	for (m = 0; m < rd->n_sig_hits; ++m) {
		const ri_sig_hits_t *y = &rd->sig_hits[m];
		if (m == cand || y->n_hits <= best_u || strcmp(ri_sig_taxon(ri, y->rid), "*") == 0) continue;
		u = sorted_hits(x, y, r->gap);
		if (u > best_u) { best_u = u; cand = m; }
	}

	if (call != UINT32_MAX) out->c1 = slopes_try(r, x, cx, ri, &rd->sig_hits[call], read_events, &out->e1);
	if (cand != call) out->c2 = slopes_try(r, x, cx, ri, &rd->sig_hits[cand], read_events, &out->e2);
	out->cand = 1;
	out->called = call != UINT32_MAX;
	out->other = cand != call;
}

void ri_slopes_commit(ri_slopes_t *r, int call_e, double call_p, const ri_slopes_read_t *sr)
{
	r->read_entry = call_e;
	r->read_p = call_p;
	r->read_chain = 0;
	if (!sr->cand) return;
	if (sr->called) {
		++r->n_called;
		r->n_other_sig += sr->other;
	}
	r->chain = sr->c2.p < sr->c1.p ? sr->c2 : sr->c1;
	if (!(r->chain.p <= r->alpha)) return;
	r->chain_entry = sr->c2.p < sr->c1.p ? sr->e2 : sr->e1;
	r->n_from_most += !(sr->c1.p <= r->alpha);
	r->ring[r->head] = r->chain;
	r->head = (r->head + 1) % r->n;
	if (r->held < r->n) ++r->held;
	++r->confident;
	r->read_chain = 1;
}

/* The read record's fields: the slopes of the grid the read ran on (the read p-value's multiplier), its call [entry, read
 * p-value], and its chain [entry, read p-value on the search grid, slope, r-span, residual SD] when it went into the
 * ring. */
void ri_slopes_emit(const ri_slopes_t *r, FILE *fp)
{
	fprintf(fp, ",\"slopes\":%u,\"read_pvalue\":[\"%s\",%.6e]", r->n_slopes,
	        r->read_entry < 0 ? "*" : ri_calls_entry(r->calls, r->read_entry), r->read_p);
	if (r->read_chain)
		fprintf(fp, ",\"chain\":[\"%s\",%.6e,%.6f,%.0f,%.4f]", ri_calls_entry(r->calls, r->chain_entry),
		        r->chain.p, r->chain.slope, r->chain.span, r->chain.sd);
}

static int cmp_double(const void *a, const void *b)
{
	double x = *(const double*)a, y = *(const double*)b;
	return (x > y) - (x < y);
}

/* The middle value, or the mean of the two middle values (Python's statistics.median). */
static double median_of(double *x, uint32_t n)
{
	qsort(x, n, sizeof(double), cmp_double);
	return (n & 1) ? x[n / 2] : 0.5 * (x[n / 2 - 1] + x[n / 2]);
}

void ri_slopes_end_batch(ri_slopes_t *r, ri_slope_grid_t *grid, uint64_t total_reads)
{
	double smin, smax, span, sd = 0.0, w, step, lo, hi;
	float  flo, fhi, fstep;
	char   text[96];
	uint32_t i;

	if (!r || r->held < r->n) return;
	smin = smax = r->ring[0].slope;
	for (i = 1; i < r->n; ++i) {
		if (r->ring[i].slope < smin) smin = r->ring[i].slope;
		if (r->ring[i].slope > smax) smax = r->ring[i].slope;
	}
	if (r->cover < 1.0f) {   /* the middle of the ring */
		uint32_t drop = (uint32_t)((1.0 - r->cover) / 2.0 * r->n);
		for (i = 0; i < r->n; ++i) r->tmp[i] = r->ring[i].slope;
		qsort(r->tmp, r->n, sizeof(double), cmp_double);
		smin = r->tmp[drop];
		smax = r->tmp[r->n - 1 - drop];
	}
	for (i = 0; i < r->n; ++i) r->tmp[i] = r->ring[i].span;
	span = median_of(r->tmp, r->n);
	for (i = 0; i < r->n; ++i) r->tmp[i] = r->ring[i].sd;
	sd = median_of(r->tmp, r->n);
	w = r->width_sd ? SLOPE_WIDTH_SDS * sd : (double)r->width;
	if (!(w > 0.0) || !(span > 0.0)) return;
	step = w / span;
	lo = floor(smin * 1e6) / 1e6;
	hi = lo + (ceil((smax - lo) / step) + 0.5) * step;
	snprintf(text, sizeof(text), "%.6f:%.6f:%.6f", lo, hi, step);
	if (sscanf(text, "%f:%f:%f", &flo, &fhi, &fstep) != 3 || !(fstep > 0.0f)) return;
	if (flo == grid->min && fhi == grid->max && fstep == grid->step) return;
	grid->min  = flo;
	grid->max  = fhi;
	grid->step = fstep;
	r->n_slopes = ri_grid_n_slopes(grid);
	++r->changes;
	if (r->emit_fp)
		fprintf(r->emit_fp, "{\"kind\":\"slopes\",\"total_reads\":%llu,\"grid\":\"%s\",\"slopes\":%u,\"chains\":%u,"
		        "\"confident\":%llu,\"slope_min\":%.6f,\"slope_max\":%.6f,\"span_median\":%.1f,\"sd_median\":%.4f,"
		        "\"width\":%.4f}\n", (unsigned long long)total_reads, text, r->n_slopes, r->n,
		        (unsigned long long)r->confident, smin, smax, span, sd, w);
	fprintf(stderr, "[M::%s] after %llu reads: slope grid %s (%u slopes) from the last %u chains: slopes %.4f-%.4f, "
	        "median span %.1f, W %.3f\n", __func__, (unsigned long long)total_reads, text, r->n_slopes, r->n,
	        smin, smax, span, w);
}
