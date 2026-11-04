/* calls.c: read calls (read p-value) and target calls; see calls.h. */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif

#include "calls.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "khash.h"
#include "neighbor_check.h"
#include "utils.h"

KHASH_MAP_INIT_STR(calls_entry, int)

#define CALLS_T_MAX 1e-3   /* the loosest per-read cutoff: read calls get the near-neighbor check at p <= 1e-3 */

struct ri_calls_s {
	const ri_idx_t *ri;
	int       n_entry;            /* K, the panel entries */
	char    **entry;              /* their names, in panel order */
	uint32_t *n_sigs;             /* the index's signatures per entry */
	khash_t(calls_entry) *of;     /* index taxon -> entry */
	int       n_alpha;
	double    alpha[RI_CALLS_MAX_ALPHA];
	double    alpha_max;          /* the loosest cutoff: the near-neighbor check runs on its calls */
	int       averaged;           /* 1: the target call averaged over cutoffs in [t_min, CALLS_T_MAX] (default) */
	double    t_min;
	double  **lst;                /* averaged test, [rule][entry]: p of the calls so far, sorted, until the target call */
	uint32_t *lst_n, *lst_m;
	double   *stat;               /* averaged test: ln M at the target call */
	double    alpha_run;
	double    threshold;          /* ln(K / alpha_run) */
	nc_db_t  *nc;                 /* near-neighbor check regions (--neighbor-check); NULL = off */
	uint32_t  nc_n;               /* near-neighbor checks run and their time, summed over the threads' scratch */
	double    nc_s;
	int       n_rule;             /* 1: read calls; 2: also the calls the near-neighbor check keeps */
	uint32_t *calls;              /* [rule][alpha][entry]: read calls so far */
	uint32_t *first;              /* reads so far at the target call; 0: not called */
	float    *first_s;            /* seconds since detection began at the target call */
	FILE     *read_fp, *target_fp;
	double    t0;
};

/* One thread's scratch for the read p-value (each entry's set) and the near-neighbor check (nc_keep's buffer). */
struct ri_calls_scratch_s {
	int32_t  *pick;
	uint32_t *pick_count;
	void     *nc_buf;
	uint32_t  nc_n;
	double    nc_s;
};

static double log_pmf(uint64_t k, uint64_t n, double q)
{
	return lgamma((double)n + 1.0) - lgamma((double)k + 1.0) - lgamma((double)(n - k) + 1.0)
	     + (double)k * log(q) + (double)(n - k) * log1p(-q);
}

/* log sum_{j >= m} pmf(j) for m at or above the mode, where the terms fall from m upward. */
static double log_upper(uint64_t m, uint64_t n, double q)
{
	double t = 1.0, s = 1.0, r = q / (1.0 - q);
	uint64_t k;
	for (k = m; k < n; ++k) {
		t *= (double)(n - k) / (double)(k + 1) * r;
		s += t;
		if (t < s * 1e-17) break;
	}
	return log_pmf(m, n, q) + log(s);
}

/* log sum_{j <= m} pmf(j) for m at or below the mode, where the terms fall from m downward. */
static double log_lower(uint64_t m, uint64_t n, double q)
{
	double t = 1.0, s = 1.0, r = (1.0 - q) / q;
	uint64_t k;
	for (k = m; k > 0; --k) {
		t *= (double)k / (double)(n - k + 1) * r;
		s += t;
		if (t < s * 1e-17) break;
	}
	return log_pmf(m, n, q) + log(s);
}

/* Each tail is summed directly when it is the smaller side of the mean, else as 1 minus the other tail. */
double ri_log_binom_sf(uint64_t m, uint64_t n, double q)
{
	if (m == 0) return 0.0;
	if (m > n || q <= 0.0) return -INFINITY;
	if (q >= 1.0) return 0.0;
	if ((double)m > (double)n * q) return log_upper(m, n, q);
	return log1p(-exp(log_lower(m - 1, n, q)));
}

double ri_log_binom_cdf(uint64_t m, uint64_t n, double q)
{
	if (m >= n || q <= 0.0) return 0.0;
	if (q >= 1.0) return -INFINITY;
	if ((double)m < (double)n * q) return log_lower(m, n, q);
	return log1p(-exp(log_upper(m + 1, n, q)));
}

/* Target call at one cutoff: log LR of c >= 1 calls in n reads, rate log-uniform on [p0, 1] against p0, that is
 * betaln(c, n - c + 1) + log I_{1-p0}(n - c + 1, c) - log ln(1/p0) - c log p0 - (n - c) log(1 - p0), with
 * B(c, n - c + 1) = 1 / (c C(n, c)) and I_{1-p0}(n - c + 1, c) = P[Bin(n, p0) <= c - 1]. */
double ri_calls_log_lr(uint64_t c, uint64_t n, double p0)
{
	return ri_log_binom_cdf(c - 1, n, p0) - log((double)c) - log_pmf(c, n, p0) - log(log(1.0 / p0));
}

/* Target call averaged over cutoffs: ln of M = (1/L) int LR_t(c(t), n) dt/t over t from t_min to
 * t_max, L = ln(t_max / t_min), c(t) the calls with p <= t (p sorted, each at least t_min). c(t) is constant between the
 * p values, so the integral is taken piecewise there, in spans of at most 1 / (c + 1) in ln t with 4-point
 * Gauss-Legendre; below the smallest p, c(t) = 0 adds nothing. */
double ri_calls_log_mixture(const double *p, uint32_t n_p, uint64_t n, double t_min, double t_max)
{
	static const double gx[4] = { -0.86113631159405257522, -0.33998104358485626481,
	                              0.33998104358485626481, 0.86113631159405257522 };
	static const double gw[4] = { 0.34785484513745385737, 0.65214515486254614263,
	                              0.65214515486254614263, 0.34785484513745385737 };
	double top = -INFINITY, sum = 0.0;
	uint32_t k = 0, c, j, g, m;
	while (k < n_p) {
		double a, b, h;
		for (c = k + 1; c < n_p && p[c] == p[k]; ++c) {}   /* c calls at p <= p[k] */
		a = log(p[k]);
		b = log(c < n_p ? p[c] : t_max);
		if (b > a) {
			m = (uint32_t)ceil((b - a) * (double)(c + 1));
			h = (b - a) / m;
			for (j = 0; j < m; ++j)
				for (g = 0; g < 4; ++g) {
					double u = a + (b - a) * (double)j / (double)m + h * (gx[g] + 1.0) / 2.0;
					double v = ri_calls_log_lr(c, n, exp(u)) + log(h * gw[g] / 2.0);
					if (v > top) { sum = sum * exp(top - v) + 1.0; top = v; }
					else sum += exp(v - top);
				}
		}
		k = c;
	}
	return top == -INFINITY ? -INFINITY : top + log(sum) - log(log(t_max / t_min));
}

const char *ri_calls_entry(const ri_calls_t *c, int e)
{
	return c->entry[e];
}

static const char *calls_taxon(const ri_idx_t *ri, uint32_t rid)
{
	if (!ri->sig_source_taxon || !ri->sig_source_taxon[rid] || !ri->sig_source_taxon[rid][0]) return "*";
	return ri->sig_source_taxon[rid];
}

static int calls_entry_of(ri_calls_t *c, const char *name)
{
	int e;
	for (e = 0; e < c->n_entry; ++e)
		if (strcmp(c->entry[e], name) == 0) return e;
	c->entry = (char**)realloc(c->entry, (size_t)(c->n_entry + 1) * sizeof(char*));
	c->entry[c->n_entry] = strdup(name);
	return c->n_entry++;
}

static void calls_map(ri_calls_t *c, const char *taxon, int e)
{
	int absent;
	char *key = strdup(taxon);
	khint_t k = kh_put(calls_entry, c->of, key, &absent);
	if (!absent) {
		fprintf(stderr, "[WARNING] --panel lists taxon %s twice; the first row counts\n", taxon);
		free(key);
		return;
	}
	kh_val(c->of, k) = e;
}

/* Columns taxon and name; rows sharing a name are one entry, in order of first appearance. */
static int calls_load_panel(ri_calls_t *c, const char *path)
{
	char line[8192];
	int col_taxon = -1, col_name = -1, rows = 0;
	FILE *fp = fopen(path, "r");
	if (!fp) {
		fprintf(stderr, "[ERROR] cannot open --panel %s\n", path);
		return -1;
	}
	while (fgets(line, sizeof(line), fp)) {
		char *f[64], *s = line, *t;
		int n = 0, i;
		line[strcspn(line, "\r\n")] = '\0';
		if (!line[0]) continue;
		f[n++] = s;
		while (n < 64 && (t = strchr(s, '\t')) != NULL) { *t = '\0'; s = t + 1; f[n++] = s; }
		if (col_taxon < 0) {
			for (i = 0; i < n; ++i) {
				if (strcmp(f[i], "taxon") == 0) col_taxon = i;
				if (strcmp(f[i], "name") == 0) col_name = i;
			}
			if (col_taxon < 0 || col_name < 0) break;
			continue;
		}
		if (col_taxon >= n || col_name >= n) {
			fprintf(stderr, "[ERROR] --panel %s: a row has fewer columns than the header\n", path);
			fclose(fp);
			return -1;
		}
		calls_map(c, f[col_taxon], calls_entry_of(c, f[col_name]));
		++rows;
	}
	fclose(fp);
	if (col_taxon < 0 || col_name < 0 || rows == 0) {
		fprintf(stderr, "[ERROR] --panel %s needs a header with columns taxon and name, and rows\n", path);
		return -1;
	}
	return 0;
}

void ri_calls_destroy(ri_calls_t *c)
{
	int e;
	size_t m;
	if (!c) return;
	if (c->of) {
		khint_t k;
		for (k = kh_begin(c->of); k != kh_end(c->of); ++k)
			if (kh_exist(c->of, k)) free((char*)kh_key(c->of, k));
		kh_destroy(calls_entry, c->of);
	}
	for (e = 0; e < c->n_entry; ++e) free(c->entry[e]);
	free(c->entry); free(c->n_sigs); free(c->calls); free(c->first); free(c->first_s);
	if (c->lst) for (m = 0; m < (size_t)c->n_rule * (size_t)c->n_entry; ++m) free(c->lst[m]);
	free(c->lst); free(c->lst_n); free(c->lst_m); free(c->stat);
	if (c->nc) nc_free(c->nc);
	if (c->read_fp) fclose(c->read_fp);
	if (c->target_fp) fclose(c->target_fp);
	free(c);
}

ri_calls_t *ri_calls_init(const ri_idx_t *ri, const ri_detopt_t *opt)
{
	ri_calls_t *c = (ri_calls_t*)calloc(1, sizeof(ri_calls_t));
	uint32_t rid, n_out = 0;
	int a, e;
	size_t m;
	if (!c) return NULL;
	c->ri = ri;
	c->of = kh_init(calls_entry);
	if (opt->calls_panel) {
		if (calls_load_panel(c, opt->calls_panel) < 0) { ri_calls_destroy(c); return NULL; }
	} else {
		for (rid = 0; rid < ri->n_seq; ++rid) {
			const char *t = calls_taxon(ri, rid);
			if (strcmp(t, "*") != 0 && kh_get(calls_entry, c->of, t) == kh_end(c->of))
				calls_map(c, t, calls_entry_of(c, t));
		}
	}
	if (c->n_entry == 0) {
		fprintf(stderr, "[ERROR] read calls: no panel entries (the index has no taxa)\n");
		ri_calls_destroy(c);
		return NULL;
	}
	c->n_sigs = (uint32_t*)calloc((size_t)c->n_entry, sizeof(uint32_t));
	for (rid = 0; rid < ri->n_seq; ++rid) {
		const char *t = calls_taxon(ri, rid);
		khint_t k = kh_get(calls_entry, c->of, t);
		if (k == kh_end(c->of)) { if (strcmp(t, "*") != 0) ++n_out; continue; }
		++c->n_sigs[kh_val(c->of, k)];
	}
	if (n_out)
		fprintf(stderr, "[WARNING] %u signatures have a taxon missing from --panel; reads hitting them get no call "
		        "for it\n", n_out);

	c->averaged = opt->n_alpha_read == 0;
	c->t_min = opt->t_min;
	c->n_alpha = c->averaged ? 1 : opt->n_alpha_read;
	for (a = 0; a < c->n_alpha; ++a) {
		c->alpha[a] = c->averaged ? CALLS_T_MAX : opt->alpha_read[a];
		if (c->alpha[a] > c->alpha_max) c->alpha_max = c->alpha[a];
	}
	c->alpha_run = opt->alpha_run;
	c->threshold = log((double)c->n_entry / c->alpha_run);
	c->n_rule = 1;
	if (opt->neighbor_regions) {     /* margin 0: the band itself */
		c->nc = nc_load(opt->neighbor_regions, opt->neighbor_pore_model, 0);
		if (!c->nc) { ri_calls_destroy(c); return NULL; }
		c->n_rule = 2;
	}
	m = (size_t)c->n_rule * (size_t)c->n_alpha * (size_t)c->n_entry;
	c->calls    = (uint32_t*)calloc(m, sizeof(uint32_t));
	c->first    = (uint32_t*)calloc(m, sizeof(uint32_t));
	c->first_s  = (float*)calloc(m, sizeof(float));
	if (c->averaged) {
		c->lst   = (double**)calloc(m, sizeof(double*));
		c->lst_n = (uint32_t*)calloc(m, sizeof(uint32_t));
		c->lst_m = (uint32_t*)calloc(m, sizeof(uint32_t));
		c->stat  = (double*)calloc(m, sizeof(double));
	}
	if (opt->read_calls_path) {
		c->read_fp = fopen(opt->read_calls_path, "w");
		if (c->read_fp)
			fputs("read\tposition\tentry\tp\tnonoverlap_hits\thits\tread_events\tsig_len\tn_sigs\tslopes\trid\tsignature\t"
			      "neighbor_check\n", c->read_fp);
	}
	if (opt->target_calls_path) c->target_fp = fopen(opt->target_calls_path, "w");
	if ((opt->read_calls_path && !c->read_fp) || (opt->target_calls_path && !c->target_fp)) {
		fprintf(stderr, "[ERROR] cannot write --read-calls or --target-calls\n");
		ri_calls_destroy(c);
		return NULL;
	}
	fprintf(stderr, "[M::%s] %d panel entries (signatures):", __func__, c->n_entry);
	for (e = 0; e < c->n_entry; ++e) fprintf(stderr, "%s %s (%u)", e ? "," : "", c->entry[e], c->n_sigs[e]);
	if (c->averaged) fprintf(stderr, "; target call averaged over read p-value cutoffs %g to %g", c->t_min, CALLS_T_MAX);
	else for (a = 0, fprintf(stderr, "; target call at alpha_read"); a < c->n_alpha; ++a) fprintf(stderr, " %g", c->alpha[a]);
	fprintf(stderr, "; alpha_run %g; target call threshold ln(K/alpha_run) = %.4f; near-neighbor check %s\n",
	        c->alpha_run, c->threshold, c->nc ? opt->neighbor_regions : "off");
	c->t0 = ri_realtime();
	return c;
}

/* The tests an entry's read p-value covers: its signatures in the index (Bonferroni over signatures). */
static double calls_tests(const ri_calls_t *c, int e)
{
	return (double)c->n_sigs[e];
}

ri_calls_scratch_t *ri_calls_scratch_init(const ri_calls_t *c)
{
	ri_calls_scratch_t *x = (ri_calls_scratch_t*)calloc(1, sizeof(ri_calls_scratch_t));
	if (!x) return NULL;
	x->pick = (int32_t*)calloc((size_t)c->n_entry, sizeof(int32_t));
	x->pick_count = (uint32_t*)calloc((size_t)c->n_entry, sizeof(uint32_t));
	x->nc_buf = c->nc ? malloc(nc_buf_size(c->nc)) : NULL;
	if (!x->pick || !x->pick_count || (c->nc && !x->nc_buf)) {
		free(x->pick); free(x->pick_count); free(x->nc_buf); free(x);
		return NULL;
	}
	return x;
}

void ri_calls_scratch_destroy(ri_calls_t *c, ri_calls_scratch_t *x)
{
	if (!x) return;
	if (c) { c->nc_n += x->nc_n; c->nc_s += x->nc_s; }
	free(x->pick); free(x->pick_count); free(x->nc_buf); free(x);
}

int ri_read_pvalue(const ri_calls_t *c, ri_calls_scratch_t *x, const ri_calls_set_t *sets, uint32_t n_sets,
                uint32_t read_events, uint32_t n_slopes, double *p_out, int *set_out)
{
	int e, best = -1;
	uint32_t i;
	double best_p = 1.0;
	for (e = 0; e < c->n_entry; ++e) { x->pick[e] = -1; x->pick_count[e] = 0; }
	/* each entry's set: the largest count, the later set on ties */
	for (i = 0; i < n_sets; ++i) {
		khint_t k = kh_get(calls_entry, c->of, sets[i].taxon);
		if (k == kh_end(c->of)) continue;
		e = kh_val(c->of, k);
		if (sets[i].nonoverlap_max >= x->pick_count[e]) { x->pick_count[e] = sets[i].nonoverlap_max; x->pick[e] = (int32_t)i; }
	}
	for (e = 0; e < c->n_entry; ++e) {
		const ri_calls_set_t *y;
		uint32_t n;
		double B, p;
		if (x->pick[e] < 0 || sets[x->pick[e]].nonoverlap == 0) continue;
		y = &sets[x->pick[e]];
		n = y->units > y->nonoverlap ? y->units : y->nonoverlap;
		B = ((double)read_events + (double)y->len) / 10.0;
		if (B < 1.0) B = 1.0;
		p = (double)n_slopes * B * exp(ri_log_binom_sf(y->nonoverlap, n, 1.0 / B)) * calls_tests(c, e);
		if (p < best_p) { best_p = p; best = e; }      /* the first entry on ties; p >= 1 is no candidate */
	}
	*p_out = best_p;
	if (set_out) *set_out = best >= 0 ? x->pick[best] : -1;
	return best;
}

/* One more read call of entry e at alpha a under rule r with n reads so far; the target call when log LR reaches the
 * threshold. The ratio only falls between two calls of e, so checking at calls finds the first read where it reaches it. */
static void calls_add(ri_calls_t *c, int r, int a, int e, uint32_t n)
{
	size_t i = ((size_t)r * (size_t)c->n_alpha + (size_t)a) * (size_t)c->n_entry + (size_t)e;
	++c->calls[i];
	if (!c->first[i] && ri_calls_log_lr(c->calls[i], n, c->alpha[a]) >= c->threshold) {
		c->first[i] = n;
		c->first_s[i] = (float)(ri_realtime() - c->t0);
		fprintf(stderr, "[M::ri_calls] target call%s: %s at read %u (alpha_read %g, %u read calls, %.1f s)\n",
		        r ? " after the near-neighbor check" : "", c->entry[e], n, c->alpha[a], c->calls[i], c->first_s[i]);
	}
}

/* Averaged test: one more call of entry e with read p-value p under rule r, n reads so far. The p values stay in a
 * sorted list only until the target call. */
static void calls_add_avg(ri_calls_t *c, int r, int e, double p, uint32_t n)
{
	size_t i = (size_t)r * (size_t)c->n_entry + (size_t)e;
	uint32_t j;
	++c->calls[i];
	if (c->first[i]) return;
	if (c->lst_n[i] == c->lst_m[i]) {
		c->lst_m[i] = c->lst_m[i] ? 2 * c->lst_m[i] : 16;
		c->lst[i] = (double*)realloc(c->lst[i], sizeof(double) * c->lst_m[i]);
	}
	if (p < c->t_min) p = c->t_min;
	for (j = c->lst_n[i]++; j > 0 && c->lst[i][j - 1] > p; --j) c->lst[i][j] = c->lst[i][j - 1];
	c->lst[i][j] = p;
	c->stat[i] = ri_calls_log_mixture(c->lst[i], c->lst_n[i], n, c->t_min, CALLS_T_MAX);
	if (c->stat[i] >= c->threshold) {
		c->first[i] = n;
		c->first_s[i] = (float)(ri_realtime() - c->t0);
		free(c->lst[i]);
		c->lst[i] = NULL;
		c->lst_n[i] = c->lst_m[i] = 0;
		fprintf(stderr, "[M::ri_calls] target call%s (averaged over cutoffs): %s at read %u (%u read calls, %.1f s)\n",
		        r ? " after the near-neighbor check" : "", c->entry[e], n, c->calls[i], c->first_s[i]);
	}
}

int ri_calls_neighbor_wanted(const ri_calls_t *c, double p)
{
	return c && c->nc && p <= c->alpha_max;
}

int ri_calls_neighbor_check(const ri_calls_t *c, ri_calls_scratch_t *x, uint32_t rid, const float *levels,
                            uint32_t n_levels, const nc_band_t *band)
{
	double t = ri_realtime();
	int k = nc_keep(c->nc, rid, levels, n_levels, band, x->nc_buf, NULL);
	x->nc_s += ri_realtime() - t;
	++x->nc_n;
	return k;
}

void ri_calls_record(ri_calls_t *c, uint64_t read_id, const char *read_name, const ri_calls_set_t *set, int e,
                     double p, uint32_t read_events, uint32_t n_slopes, int checked)
{
	int a;
	if (!c || e < 0) return;
	if (c->read_fp) {
		const ri_idx_t *ri = c->ri;
		const char *sig = set->rid < ri->n_seq ? ri->seq[set->rid].name : NULL;
		fprintf(c->read_fp, "%s\t%llu\t%s\t%.17g\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%s\t%s\n", read_name ? read_name : "*",
		        (unsigned long long)read_id, c->entry[e], p, set->nonoverlap, set->units > set->nonoverlap ? set->units : set->nonoverlap,
		        read_events, set->len, c->n_sigs[e], n_slopes, set->rid, sig ? sig : "*",
		        checked == 1 ? "kept" : checked == 0 ? "dropped" : checked == -1 ? "no band" : "-");
	}
	if (c->averaged) {
		if (p > CALLS_T_MAX) return;
		calls_add_avg(c, 0, e, p, (uint32_t)read_id + 1);
		if (c->n_rule > 1 && checked == 1) calls_add_avg(c, 1, e, p, (uint32_t)read_id + 1);
		return;
	}
	for (a = 0; a < c->n_alpha; ++a) {
		if (p > c->alpha[a]) continue;
		calls_add(c, 0, a, e, (uint32_t)read_id + 1);
		if (c->n_rule > 1 && checked == 1) calls_add(c, 1, a, e, (uint32_t)read_id + 1);
	}
}

void ri_calls_finish(ri_calls_t *c, uint64_t n_reads)
{
	int r, a, e;
	if (!c) return;
	if (c->nc) fprintf(stderr, "[M::%s] near-neighbor check: %u calls in %.2f s, summed over threads\n", __func__, c->nc_n, c->nc_s);
	if (c->target_fp) {
		fputs("entry\trule\ttest\talpha_read\tt_min\tcalled\tfirst_call_reads\tfirst_call_time_s\tcalls\treads\tlog_lr\t"
		      "threshold\n", c->target_fp);
		for (r = 0; r < c->n_rule; ++r)
		for (a = 0; a < c->n_alpha; ++a)
			for (e = 0; e < c->n_entry; ++e) {
				size_t i = ((size_t)r * (size_t)c->n_alpha + (size_t)a) * (size_t)c->n_entry + (size_t)e;
				fprintf(c->target_fp, "%s\t%s\t%s\t%g\t", c->entry[e], r ? "read_pvalue+neighbor_check" : "read_pvalue",
				        c->averaged ? "averaged" : "single", c->alpha[a]);
				if (c->averaged) fprintf(c->target_fp, "%.17g\t", c->t_min);
				else fputc('\t', c->target_fp);
				fprintf(c->target_fp, "%d\t", c->first[i] ? 1 : 0);
				if (c->first[i]) fprintf(c->target_fp, "%u\t%.3f\t", c->first[i], c->first_s[i]);
				else fputs("\t\t", c->target_fp);
				fprintf(c->target_fp, "%u\t%llu\t", c->calls[i], (unsigned long long)n_reads);
				if (c->averaged && (c->first[i] || c->lst_n[i]))  /* at the target call, else at the end */
					fprintf(c->target_fp, "%.17g", c->first[i] ? c->stat[i]
					        : ri_calls_log_mixture(c->lst[i], c->lst_n[i], n_reads, c->t_min, CALLS_T_MAX));
				else if (!c->averaged && c->calls[i])
					fprintf(c->target_fp, "%.17g", ri_calls_log_lr(c->calls[i], n_reads, c->alpha[a]));
				else fputs("NA", c->target_fp);
				fprintf(c->target_fp, "\t%.17g\n", c->threshold);
			}
		fclose(c->target_fp);
		c->target_fp = NULL;
	}
	if (c->read_fp) {
		fclose(c->read_fp);
		c->read_fp = NULL;
	}
}
