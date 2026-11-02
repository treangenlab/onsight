/* diagonal_bands.c: a signature's seed hits counted in diagonal bands; see diagonal_bands.h. */
#include <math.h>
#include <pthread.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "diagonal_bands.h"
#include "utils.h"

uint32_t ri_grid_n_slopes(const ri_slope_grid_t *g)
{
	double lo = g->min > 0.0f ? (double)g->min : 1.0, step = (double)g->step;
	return (step > 0.0 && g->max > g->min)
	     ? (uint32_t)floor(((double)g->max - lo) / step + 1e-6) + 1u : 1u;
}

/* Slope 1 is computed in integers, so the default grid (slope 1 only) gives the count of (q - r) / 10. */
uint32_t ri_band_fullest(const ri_sig_hits_t *s, const ri_slope_grid_t *g)
{
	uint32_t i, j, n, k, n_slopes;
	int64_t  diag;
	int64_t *buckets;
	uint32_t *counts;
	uint32_t  best;
	uint32_t  cap;
	double    lo, step;

	if (!s || !s->hits || s->n_hits == 0) return 0;
	n = s->n_hits;
	if (n == 1) return 1;

	lo = g->min > 0.0f ? (double)g->min : 1.0;
	step = (double)g->step;
	n_slopes = (step > 0.0 && g->max > g->min)
	         ? (uint32_t)floor(((double)g->max - lo) / step + 1e-6) + 1u : 1u;

	cap = 1u;
	while (cap < n * 2u) cap <<= 1;
	if (cap < 16) cap = 16;
	buckets = (int64_t*)malloc(sizeof(int64_t) * cap);
	counts  = (uint32_t*)malloc(sizeof(uint32_t) * cap);
	if (!buckets || !counts) {
		free(buckets); free(counts);
		return 0;
	}

	best = 0;
	for (k = 0; k < n_slopes; ++k) {
		double b = lo + step * (double)k;
		int    unit = fabs(b - 1.0) < 1e-6;
		for (j = 0; j < cap; ++j) { buckets[j] = INT64_MIN; counts[j] = 0; }
		for (i = 0; i < n; ++i) {
			uint64_t h;
			uint32_t idx;
			if (unit)
				diag = ((int64_t)s->hits[i].q_pos - (int64_t)s->hits[i].ref_pos) / 10;
			else    /* the cast truncates toward zero, as the integer division above does */
				diag = (int64_t)(((double)s->hits[i].q_pos - b * (double)s->hits[i].ref_pos) / 10.0);
			h = (uint64_t)(diag * 0x9E3779B97F4A7C15ULL);
			idx = (uint32_t)(h & (cap - 1));
			while (buckets[idx] != INT64_MIN && buckets[idx] != diag) {
				idx = (idx + 1) & (cap - 1);
			}
			buckets[idx] = diag;
			counts[idx]++;
			if (counts[idx] > best) best = counts[idx];
		}
	}
	free(buckets);
	free(counts);
	return best;
}

/* Grows one thread's buffers: n hits, a band table of cap slots (kept empty between uses: band_slope resets the slots
 * it fills), n_ev pair events. */
static int band_buf_reserve(ri_band_buf_t *b, uint32_t n, uint32_t cap, uint32_t n_ev)
{
	if (n > b->n) {
		uint64_t *a = (uint64_t*)realloc(b->a, sizeof(uint64_t) * n);
		uint32_t *used = a ? (uint32_t*)realloc(b->used, sizeof(uint32_t) * n) : NULL;
		if (a) b->a = a;
		if (!used) return 0;
		b->used = used;
		b->n = n;
	}
	if (cap > b->cap) {
		int64_t  *keys  = (int64_t*)realloc(b->keys, sizeof(int64_t) * cap);
		uint32_t *hits  = keys ? (uint32_t*)realloc(b->hits, sizeof(uint32_t) * cap) : NULL;
		uint32_t *dec   = hits ? (uint32_t*)realloc(b->dec, sizeof(uint32_t) * cap) : NULL;
		uint32_t *lastq = dec ? (uint32_t*)realloc(b->lastq, sizeof(uint32_t) * cap) : NULL;
		uint32_t j;
		if (keys) b->keys = keys;
		if (hits) b->hits = hits;
		if (dec) b->dec = dec;
		if (!lastq) return 0;
		b->lastq = lastq;
		for (j = 0; j < cap; ++j) b->keys[j] = INT64_MIN;
		b->cap = cap;
	}
	if (n_ev > b->n_ev) {
		uint64_t *ev = (uint64_t*)realloc(b->ev, sizeof(uint64_t) * n_ev);
		if (!ev) return 0;
		b->ev = ev;
		b->n_ev = n_ev;
	}
	return 1;
}

/* One slope b: per band its hits, its non-overlapping hits and the last counted position; the fullest band so far in
 * *best (non-overlapping hits) and *best_x (its hits), the most hits in one band in *xm. Only the table slots it
 * fills are read and reset, so a few hits cost a few steps. */
static inline void band_slope(const uint64_t *a, uint32_t n, double b, uint32_t gap, uint32_t cap, ri_band_buf_t *cb,
                              uint32_t *best, uint32_t *best_x, uint32_t *xm)
{
	int64_t  *keys = cb->keys;
	uint32_t *hits = cb->hits, *dec = cb->dec, *lastq = cb->lastq, *used = cb->used, i, t, nu = 0;
	uint32_t bd = *best, bx = *best_x, m = *xm;
	int unit = fabs(b - 1.0) < 1e-6;
	for (i = 0; i < n; ++i) {
		int64_t  q = (int64_t)(a[i] >> 32), r = (int64_t)(a[i] & 0xffffffffULL);
		int64_t  key = unit ? (q - r) / 10 : (int64_t)(((double)q - b * (double)r) / 10.0);
		uint32_t idx = (uint32_t)((uint64_t)(key * 0x9E3779B97F4A7C15ULL) & (cap - 1));
		while (keys[idx] != INT64_MIN && keys[idx] != key) idx = (idx + 1) & (cap - 1);
		if (keys[idx] == INT64_MIN) {
			keys[idx] = key; hits[idx] = 1; dec[idx] = 1; lastq[idx] = (uint32_t)q; used[nu++] = idx;
		} else {
			++hits[idx];
			if ((uint32_t)q - lastq[idx] >= gap) { ++dec[idx]; lastq[idx] = (uint32_t)q; }
			if (hits[idx] > m) m = hits[idx];
		}
	}
	for (t = 0; t < nu; ++t) {
		uint32_t j = used[t];
		if (hits[j] >= 2 && (dec[j] > bd || (dec[j] == bd && hits[j] > bx))) { bd = dec[j]; bx = hits[j]; }
		keys[j] = INT64_MIN;
	}
	*best = bd; *best_x = bx; *xm = m;
}

/* The most hits h, up to most, that h(h-1)/2 pairs allow */
static inline uint32_t band_bound(int32_t pairs, uint32_t most)
{
	uint32_t h = 1;
	while (h < most && (int64_t)(h + 1) * h / 2 <= (int64_t)pairs) ++h;
	return h;
}

/* Two hits share a band at slope b only if |dq - b dr| < 20 events (truncation toward zero makes band 0 twice as
 * wide); 21 leaves room for rounding and for slope 1, which is computed in integers. */
#define BAND_SHARE 21.0
/* Signatures with more hits are counted at every slope; their pairs would cost more than the slopes. */
#define BAND_PAIRS_MAX 64

/* The hits are sorted by read position once; each slope then needs one pass that keeps, per band, its hits, its
 * non-overlapping hits and the last counted position. Only slopes that can change the counts are counted: a band of
 * h hits needs all h(h-1)/2 pairs of its hits to be able to share a band at that slope, and its non-overlapping hits,
 * which are at least gap apart, as many such pairs at least gap apart. So the slopes are taken in runs over which
 * the same pairs can share a band, and a run is skipped once these bounds cannot beat the counts so far; the counts
 * are those of every slope. */
uint32_t ri_band_nonoverlap(const ri_sig_hits_t *s, const ri_slope_grid_t *g, uint32_t gap, ri_band_buf_t *cb,
                            uint32_t *x, uint32_t *units, uint32_t *xmax)
{
	uint32_t i, j, k, n, n_slopes, cap, best = 1, best_x = 1, nd = 1, xm = 1, n_ev = 0;
	uint64_t *a, *ev, last;
	double lo, step, hi;

	*x = 0; *units = 0; *xmax = 0;
	if (!s || !s->hits || s->n_hits == 0) return 0;
	n = s->n_hits;
	*x = 1; *units = 1; *xmax = 1;
	if (n == 1) return 1;

	lo = g->min > 0.0f ? (double)g->min : 1.0;
	step = (double)g->step;
	n_slopes = (step > 0.0 && g->max > g->min)
	         ? (uint32_t)floor(((double)g->max - lo) / step + 1e-6) + 1u : 1u;
	hi = lo + step * (double)(n_slopes - 1);
	cap = 4;
	while (cap < n * 2u) cap <<= 1;
	if (!band_buf_reserve(cb, n, cap, n <= BAND_PAIRS_MAX ? n * (n - 1) : 0)) return 0;
	a = cb->a; ev = cb->ev;

	for (i = 0; i < n; ++i) a[i] = ((uint64_t)s->hits[i].q_pos << 32) | s->hits[i].ref_pos;
	radix_sort_64(a, a + n);   /* in place, ascending as the earlier qsort */
	for (i = 1, last = a[0] >> 32; i < n; ++i)
		if ((a[i] >> 32) - last >= gap) { ++nd; last = a[i] >> 32; }

	if (n_slopes == 1 || n > BAND_PAIRS_MAX) {
		for (k = 0; k < n_slopes; ++k) band_slope(a, n, lo + step * (double)k, gap, cap, cb, &best, &best_x, &xm);
	} else {
		/* per pair, the slopes kl..kh where it could share a band: events at kl and kh + 1, flagged when the two hits
		 * are at least gap apart in the read (bit 0) and at the start (bit 1); slopes in the top bits for the sort */
		for (i = 0; i < n; ++i)
			for (j = i + 1; j < n; ++j) {
				uint64_t dq_u = (a[j] >> 32) - (a[i] >> 32);
				double dq = (double)dq_u;
				double dr = (double)(a[j] & 0xffffffffULL) - (double)(a[i] & 0xffffffffULL);
				uint64_t far = dq_u >= gap;
				uint32_t kl, kh;
				if (dr == 0.0) {
					if (dq >= BAND_SHARE) continue;
					kl = 0; kh = n_slopes - 1;
				} else if (dr > 0.0 ? dq + BAND_SHARE < (lo - 3.0 * step) * dr || dq - BAND_SHARE > (hi + 3.0 * step) * dr
				                    : dq >= BAND_SHARE || BAND_SHARE - dq < (lo - 3.0 * step) * -dr) {
					continue;   /* the pair cannot share a band anywhere on the grid: found without dividing */
				} else {
					double b1 = (dq - BAND_SHARE) / dr, b2 = (dq + BAND_SHARE) / dr, f0, f1;
					if (b1 > b2) { double t = b1; b1 = b2; b2 = t; }
					f0 = floor((b1 - lo) / step) - 1.0;
					f1 = ceil((b2 - lo) / step) + 1.0;
					if (f1 < 0.0 || f0 > (double)(n_slopes - 1)) continue;
					kl = f0 < 0.0 ? 0 : (uint32_t)f0;
					kh = f1 > (double)(n_slopes - 1) ? n_slopes - 1 : (uint32_t)f1;
				}
				ev[n_ev++] = (uint64_t)kl << 40 | 2u | far;
				ev[n_ev++] = (uint64_t)(kh + 1) << 40 | far;
			}
		radix_sort_64(ev, ev + n_ev);
		{   /* the slopes in runs [k0, k1) of the same pairs */
			int32_t np = 0, ng = 0;
			uint32_t e = 0;
			while (e < n_ev) {
				uint32_t k0 = (uint32_t)(ev[e] >> 40), k1, hb, db;
				for (; e < n_ev && (uint32_t)(ev[e] >> 40) == k0; ++e) {
					int32_t d = (ev[e] & 2) ? 1 : -1;
					np += d;
					if (ev[e] & 1) ng += d;
				}
				k1 = e < n_ev ? (uint32_t)(ev[e] >> 40) : n_slopes;
				if (np <= 0 || k0 >= n_slopes) continue;
				hb = band_bound(np, n);    /* most hits in one band */
				db = band_bound(ng, nd);   /* most non-overlapping hits in one band */
				for (k = k0; k < k1 && k < n_slopes; ++k) {
					if ((db < best || (db == best && hb <= best_x)) && hb <= xm) break;
					band_slope(a, n, lo + step * (double)k, gap, cap, cb, &best, &best_x, &xm);
				}
				if (best == nd && best_x == n && xm == n) break;   /* all hits in one band: nothing can beat it */
			}
		}
	}

	*x = best_x;
	*units = nd;
	*xmax = xm;
	return best;
}

/* lgamma(k + 1) for small k, filled by lgamma itself, so the values are the same */
#define LGF_N 4097
static double lgf[LGF_N];
static pthread_once_t lgf_once = PTHREAD_ONCE_INIT;
static void lgf_fill(void)
{
	uint32_t k;
	for (k = 0; k < LGF_N; ++k) lgf[k] = lgamma((double)k + 1.0);
}
static inline double lgf_at(uint32_t k)
{
	return k < LGF_N ? lgf[k] : lgamma((double)k + 1.0);
}

/* P(Binomial(n, q) >= m), summing terms from m until they stop adding; lq = log(q), l1q = log1p(-q). */
static double binom_tail(uint32_t n, uint32_t m, double q, double lq, double l1q)
{
	double t, s = 0.0;
	uint32_t k;
	if (m == 0 || q >= 1.0) return 1.0;
	if (m > n) return 0.0;
	pthread_once(&lgf_once, lgf_fill);
	t = exp(lgf_at(n) - lgf_at(m) - lgf_at(n - m) + m * lq + (n - m) * l1q);
	for (k = m; k <= n; ++k) {
		s += t;
		if (t < s * 1e-12) break;
		t *= (double)(n - k) / (double)(k + 1) * q / (1.0 - q);
	}
	return s < 1.0 ? s : 1.0;
}

/* The signature's n non-overlapping hits fall uniformly over B = (read + signature length) / 10 bands, so
 * P(fullest band >= m) <= B * P(Binomial(n, 1/B) >= m). */
double ri_band_pvalue(uint32_t m, uint32_t units, uint32_t read_events, uint32_t sig_len)
{
	double B = ((double)read_events + (double)sig_len) / 10.0, p, q;
	if (B < 1.0) B = 1.0;
	q = 1.0 / B;
	p = B * binom_tail(units > m ? units : m, m, q, log(q), log1p(-q));
	return p < 1.0 ? p : 1.0;
}

void ri_band_pv_read(ri_band_buf_t *b, uint32_t read_events)
{
	++b->pv_read;
	b->pv_events = read_events;
}

/* Per signature length: B, q = 1 / B, log(q), log1p(-q) and the p-value of a signature with one hit (NAN until
 * needed), filled once per read; the p-value is the one ri_band_pvalue gives. */
double ri_band_pvalue_buf(ri_band_buf_t *b, uint32_t m, uint32_t units, uint32_t sig_len)
{
	double *v, p;
	if (sig_len >= b->pv_n) {
		uint32_t n = sig_len + 1024, i;
		double *pv_v = (double*)realloc(b->pv_v, sizeof(double) * 5 * n);
		uint32_t *st = pv_v ? (uint32_t*)realloc(b->pv_stamp, sizeof(uint32_t) * n) : NULL;
		if (pv_v) b->pv_v = pv_v;
		if (!st) return ri_band_pvalue(m, units, b->pv_events, sig_len);
		for (i = b->pv_n; i < n; ++i) st[i] = 0;
		b->pv_stamp = st;
		b->pv_n = n;
	}
	v = b->pv_v + 5 * (size_t)sig_len;
	if (b->pv_stamp[sig_len] != b->pv_read) {
		double B = ((double)b->pv_events + (double)sig_len) / 10.0;
		if (B < 1.0) B = 1.0;
		v[0] = B;
		v[1] = 1.0 / B;
		v[2] = log(v[1]);
		v[3] = log1p(-v[1]);
		v[4] = NAN;
		b->pv_stamp[sig_len] = b->pv_read;
	}
	if (m == 1 && units <= 1 && !isnan(v[4])) return v[4];
	p = v[0] * binom_tail(units > m ? units : m, m, v[1], v[2], v[3]);
	p = p < 1.0 ? p : 1.0;
	if (m == 1 && units <= 1) v[4] = p;
	return p;
}

/* Exact ties go to the first slope and then to the band met first in (q, r) order, as the Python near-neighbor check's
 * engine_band. */
int ri_band_box(const ri_sig_hits_t *s, const ri_slope_grid_t *g, uint32_t gap, ri_band_buf_t *x, nc_band_t *box)
{
	uint32_t i, j, k, t, nb, n = s->n_hits, n_slopes = ri_grid_n_slopes(g), cap = 16, best = 1, best_x = 1;
	double lo = g->min > 0.0f ? (double)g->min : 1.0, step = (double)g->step;
	uint64_t *a;
	int64_t *keys;
	uint32_t *slot, *hits, *dec, *lastq, *q0, *q1, *r0, *r1, *nrev;
	memset(box, 0, sizeof(*box));
	if (n == 0) return 1;
	while (cap < n * 2u) cap <<= 1;
	if (n > x->box_n) {
		free(x->box_a); free(x->box_v);
		x->box_a = (uint64_t*)malloc(sizeof(uint64_t) * n);
		x->box_v = (uint32_t*)malloc(sizeof(uint32_t) * 8 * n);
		x->box_n = x->box_a && x->box_v ? n : 0;
	}
	if (cap > x->box_cap) {
		free(x->box_key); free(x->box_slot);
		x->box_key = (int64_t*)malloc(sizeof(int64_t) * cap);
		x->box_slot = (uint32_t*)malloc(sizeof(uint32_t) * cap);
		x->box_cap = x->box_key && x->box_slot ? cap : 0;
	}
	if (x->box_n < n || x->box_cap < cap) return 0;
	a = x->box_a; keys = x->box_key; slot = x->box_slot;
	hits = x->box_v; dec = hits + n; lastq = dec + n; q0 = lastq + n; q1 = q0 + n; r0 = q1 + n; r1 = r0 + n; nrev = r1 + n;
	/* q, r and the strand in one key, sorted in place in (q, r) order (r < 2^31: a position in one signature) */
	for (i = 0; i < n; ++i)
		a[i] = (uint64_t)s->hits[i].q_pos << 32 | (uint64_t)s->hits[i].ref_pos << 1 | (s->hits[i].rev ? 1u : 0u);
	radix_sort_64(a, a + n);
	box->q0 = box->q1 = (uint32_t)(a[0] >> 32); box->r0 = box->r1 = (uint32_t)(a[0] >> 1) & 0x7fffffffu;
	box->strand = (int)(a[0] & 1);
	for (k = 0; k < n_slopes; ++k) {
		double b = lo + step * (double)k;
		int unit = fabs(b - 1.0) < 1e-6;
		for (j = 0; j < cap; ++j) keys[j] = INT64_MIN;
		for (i = 0, nb = 0; i < n; ++i) {     /* bands numbered in the order they are first met */
			uint32_t aq = (uint32_t)(a[i] >> 32), ar = (uint32_t)(a[i] >> 1) & 0x7fffffffu, rev = (uint32_t)(a[i] & 1);
			int64_t q = (int64_t)aq, r = (int64_t)ar;
			int64_t key = unit ? (q - r) / 10 : (int64_t)(((double)q - b * (double)r) / 10.0);
			uint32_t h = (uint32_t)((uint64_t)(key * 0x9E3779B97F4A7C15ULL) & (cap - 1));
			while (keys[h] != INT64_MIN && keys[h] != key) h = (h + 1) & (cap - 1);
			if (keys[h] == INT64_MIN) {
				keys[h] = key; slot[h] = t = nb++;
				hits[t] = dec[t] = 1; lastq[t] = q0[t] = q1[t] = aq; r0[t] = r1[t] = ar; nrev[t] = rev;
			} else {
				t = slot[h];
				++hits[t];
				if (aq - lastq[t] >= gap) { ++dec[t]; lastq[t] = aq; }
				q1[t] = aq;
				if (ar < r0[t]) r0[t] = ar;
				if (ar > r1[t]) r1[t] = ar;
				nrev[t] += rev;
			}
		}
		for (t = 0; t < nb; ++t)
			if (hits[t] >= 2 && (dec[t] > best || (dec[t] == best && hits[t] > best_x))) {
				best = dec[t]; best_x = hits[t];
				box->q0 = q0[t]; box->q1 = q1[t]; box->r0 = r0[t]; box->r1 = r1[t];
				box->strand = 2 * nrev[t] > hits[t] ? 1 : 0;
			}
	}
	return 1;
}

void ri_band_buf_free(ri_band_buf_t *b)
{
	free(b->a); free(b->keys); free(b->hits); free(b->dec); free(b->lastq); free(b->used); free(b->ev);
	free(b->pv_v); free(b->pv_stamp);
	free(b->box_a); free(b->box_key); free(b->box_slot); free(b->box_v);
}
