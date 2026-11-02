/* slope_range.h: the sliding-window slope range (--stream-slope-chains).
 *
 * A read's slope is its events per position of the signature. It is not exactly 1, because event detection splits
 * some bases into two events and merges others (real R9.4.1 reads run at about 1.05-1.1). Reads start on the search
 * grid (--stream-diag-slopes). A read whose call has a read p-value at most alpha on the search grid gives its chain:
 * the slope fitted by least squares through the hits of its best band and the two bands beside it, with the chain's
 * span of signature positions and the residual SD. The last n chains are kept in a ring; after each minibatch with
 * the ring full, the grid becomes the ring's [min, max] slope at step W / S, S the median span and W width events, or
 * 4 times the median residual SD.
 *
 * In: per read, its seed hits per signature and its call by read p-value. Out: the slope grid of the next minibatch;
 * per read, its chain for the read record; the slopes lines of fd 3. State: the ring (ri_slopes_t), changed in read
 * order; ri_slopes_buf_t holds one thread's refit buffers.
 */
#ifndef SLOPE_RANGE_H
#define SLOPE_RANGE_H

#include <stdint.h>
#include <stdio.h>
#include "diagonal_bands.h"
#include "seed_hits.h"
#include "calls.h"

#ifdef __cplusplus
extern "C" {
#endif

/* A read's chain: its read p-value on the search grid, refitted slope, span of signature positions, residual SD */
typedef struct {
	double p, slope, span, sd;
} ri_chain_t;

typedef struct ri_slopes_s ri_slopes_t;

/* One thread's refit buffers, grown when a signature has more hits than any before */
typedef struct {
	uint32_t  cap_n, cap_h;
	uint64_t *a;
	int64_t  *keys;
	uint32_t *dec, *lastq;
	double   *qq, *rr;
	uint8_t  *keep;
} ri_slopes_buf_t;

/* One read's refits, made on any thread (ri_slopes_prepare) and put into the ring in read order (ri_slopes_commit):
 * c1 for the signature of its call, c2 for the one with the most non-overlapping hits when that is another */
typedef struct {
	int        cand, called, other;
	ri_chain_t c1, c2;
	int        e1, e2;                    /* their panel entries */
} ri_slopes_read_t;

/* Prints the first slopes line to emit_fp (NULL: none) */
ri_slopes_t *ri_slopes_init(const ri_slope_grid_t *search, uint32_t n_chains, double alpha, float width, int width_sd,
                            float cover, uint32_t gap, const ri_calls_t *calls, FILE *emit_fp);
void ri_slopes_destroy(ri_slopes_t *r);

void ri_slopes_buf_reserve(ri_slopes_buf_t *x, uint32_t n);
void ri_slopes_buf_free(ri_slopes_buf_t *x);

void ri_slopes_prepare(const ri_slopes_t *r, ri_slopes_buf_t *x, ri_calls_scratch_t *cx, const ri_idx_t *ri,
                       const ri_read_t *rd, int call_e, uint32_t call_rid, uint32_t read_events, ri_slopes_read_t *out);
void ri_slopes_commit(ri_slopes_t *r, int call_e, double call_p, const ri_slopes_read_t *sr);

/* The last committed read's fields of its record: the grid's slopes, its call and its chain */
void ri_slopes_emit(const ri_slopes_t *r, FILE *fp);

/* After a minibatch's last commit: renews *grid from the ring when it is full */
void ri_slopes_end_batch(ri_slopes_t *r, ri_slope_grid_t *grid, uint64_t total_reads);

#ifdef __cplusplus
}
#endif
#endif // SLOPE_RANGE_H
