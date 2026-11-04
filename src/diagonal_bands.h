/* diagonal_bands.h: a signature's seed hits counted in diagonal bands.
 *
 * A seed hit at read event q and signature position r lies on the diagonal q - b r, b the slope (the read's events
 * per position of the signature). Diagonals are grouped into bands 10 events wide, at each slope of a grid. Per
 * signature, the fullest band's hits, and its non-overlapping hits: those at least gap events apart in the read,
 * counted from the left (DALIGNER's diagonal bands; D-SOFT's non-overlapping seed hits). A read p-value multiplies the
 * smallest band p-value by the slopes and signatures tried (calls.h).
 *
 * In: one signature's seed hits, the slope grid, the gap. Out: its band counts, its band p-value, and the box of its
 * best band for the near-neighbor check. State: none; ri_band_buf_t holds one thread's grow-only buffers (the band
 * counters, and for the read being scored the band p-value terms of each signature length).
 */
#ifndef DIAGONAL_BANDS_H
#define DIAGONAL_BANDS_H

#include <math.h>
#include <stdint.h>
#include "seed_hits.h"
#include "neighbor_check.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Slopes min, min + step, ..., max (--stream-diag-slopes); min alone when step is 0. */
typedef struct {
	float min, max, step;
} ri_slope_grid_t;

/* The band of a hit at slope b: floor((q - b r) / 10), so band 0 is 10 events wide like the others (rounding toward
 * zero would put -9 to 9 in it). Slope 1 (unit) is computed in integers. */
static inline int64_t ri_band_key(int64_t q, int64_t r, double b, int unit)
{
	int64_t d = q - r;
	if (unit) return d >= 0 ? d / 10 : -((9 - d) / 10);
	return (int64_t)floor(((double)q - b * (double)r) / 10.0);
}

/* One thread's buffers, grown when a signature has more hits than any before */
typedef struct {
	uint64_t *a;                          /* the band counts */
	int64_t  *keys;
	uint32_t *hits, *dec, *lastq, *used;  /* used: the table slots one slope filled */
	uint32_t  n, cap;
	uint64_t *ev;                         /* per pair of hits, the slopes where they could share a band */
	uint32_t  n_ev;
	double   *pv_v;                       /* per signature length: B, 1 / B, their logs and the p-value of one hit */
	uint32_t *pv_stamp, pv_n, pv_read, pv_events;   /* the read each length was filled for */
	uint64_t *box_a;                      /* the band box */
	int64_t  *box_key;
	uint32_t *box_slot, *box_v, box_n, box_cap;
} ri_band_buf_t;

uint32_t ri_grid_n_slopes(const ri_slope_grid_t *g);

/* The most hits in one band, over the bands and slopes */
uint32_t ri_band_fullest(const ri_sig_hits_t *s, const ri_slope_grid_t *g);

/* The most non-overlapping hits in one band over bands and slopes, ties to the band with more hits; *x the hits in
 * that band, *units the non-overlapping hits over all bands, *xmax ri_band_fullest's count. 0 when out of memory. */
uint32_t ri_band_nonoverlap(const ri_sig_hits_t *s, const ri_slope_grid_t *g, uint32_t gap, ri_band_buf_t *b,
                            uint32_t *x, uint32_t *units, uint32_t *xmax);

/* Chance of m or more of a signature's non-overlapping hits (units of them in all) in one band, for one slope */
double ri_band_pvalue(uint32_t m, uint32_t units, uint32_t read_events, uint32_t sig_len);

/* The same, for the signatures of one read: ri_band_pv_read starts the read; the terms that depend only on the
 * signature length are kept per length for the read. */
void ri_band_pv_read(ri_band_buf_t *b, uint32_t read_events);
double ri_band_pvalue_buf(ri_band_buf_t *b, uint32_t m, uint32_t units, uint32_t sig_len);

/* The box (read events, signature positions, strand) of the band ri_band_nonoverlap counts (at least 2 hits), for the
 * near-neighbor check; without such a band, the first hit. 0 when out of memory. */
int ri_band_box(const ri_sig_hits_t *s, const ri_slope_grid_t *g, uint32_t gap, ri_band_buf_t *b, nc_band_t *box);

void ri_band_buf_free(ri_band_buf_t *b);

#ifdef __cplusplus
}
#endif
#endif // DIAGONAL_BANDS_H
