/* calls.h: read calls and target calls (derivations in the paper's Methods).
 *
 * Read p-value, per read and panel entry s: p_s = n_slopes x B x P[Bin(n, 1/B) >= m] x n_sigs[s], from the signature
 * of s with the smallest band p-value (diagonal_bands.h): m its non-overlapping hit count, n its non-overlapping hits
 * over all bands, B = (read events + signature length) / 10 the bands a chance hit could fall in, n_sigs[s] the
 * index's signatures of s. The read is called for the entry with the smallest p when p <= alpha_read.
 *
 * Target call at one cutoff (--alpha-read): with c read calls in the first n reads, the likelihood ratio of a call
 * rate p against p0 = alpha_read, averaged over a log-uniform p on [p0, 1], is
 *   LR = P[Bin(n, p0) <= c - 1] / (c x P[Bin(n, p0) = c] x ln(1/p0)).
 * Target call averaged over cutoffs (the default): M = (1/L) int LR_t(c(t), n) dt/t over t in [t_min, 1e-3], with
 * L = ln(1e-3 / t_min) and c(t) the calls at p <= t. An entry is called when ln LR (or ln M) >= ln(K / alpha_run), K the
 * panel entries; by Ville's inequality the chance of any false target call in the run is at most alpha_run.
 *
 * Panel entries come from --panel (rows sharing a name are one entry, e.g. the genus Enterobacter from 25 species
 * sets); without it each index taxon is an entry.
 *
 * In: per read, its signature sets and its near-neighbor check. Out: the --read-calls rows and, at the end, the
 * --target-calls table. State: per entry, its read calls and its test, changed in read order.
 */
#ifndef CALLS_H
#define CALLS_H

#include <stdint.h>

#include "neighbor_check.h"
#include "index.h"
#include "options.h"

#ifdef __cplusplus
extern "C" {
#endif

/* One signature set (index taxon) of a read, from the nonoverlap fields of its read event. */
typedef struct {
	const char *taxon;
	uint32_t    nonoverlap_max; /* its largest non-overlapping hit count over signatures */
	uint32_t    rid;     /* its signature with the smallest band p-value, */
	uint32_t    nonoverlap; /* that signature's non-overlapping hit count, */
	uint32_t    units;   /* its non-overlapping hits over all bands */
	uint32_t    len;     /* and its length */
} ri_calls_set_t;

typedef struct ri_calls_s ri_calls_t;

/* One thread's scratch for ri_read_pvalue and ri_calls_neighbor_check, so that the reads of a minibatch can be
 * scored at once. */
typedef struct ri_calls_scratch_s ri_calls_scratch_t;
ri_calls_scratch_t *ri_calls_scratch_init(const ri_calls_t *c);     /* NULL when out of memory */
void ri_calls_scratch_destroy(ri_calls_t *c, ri_calls_scratch_t *x); /* adds its near-neighbor check count and time to c */

/* NULL on error (message on stderr). */
ri_calls_t *ri_calls_init(const ri_idx_t *ri, const ri_detopt_t *opt);
const char *ri_calls_entry(const ri_calls_t *c, int e);     /* entry name */

/* The read p-value of one read: returns the candidate entry (-1: none) and sets *p to its p (1 without one) and *set to the read's
 * set it came from. n_slopes: slopes of the grid the read ran on. */
int ri_read_pvalue(const ri_calls_t *c, ri_calls_scratch_t *x, const ri_calls_set_t *sets, uint32_t n_sets,
                uint32_t read_events, uint32_t n_slopes, double *p, int *set);

/* The near-neighbor check (--neighbor-check) runs on the read: it is on and the read is called at the loosest
 * alpha_read. */
int ri_calls_neighbor_wanted(const ri_calls_t *c, double p);

/* Near-neighbor check of one call on its signature's band: 1 kept, 0 dropped, -1 no check possible (counted as dropped). */
int ri_calls_neighbor_check(const ri_calls_t *c, ri_calls_scratch_t *x, uint32_t rid, const float *levels,
                            uint32_t n_levels, const nc_band_t *band);

/* One scored read, in read order, after ri_read_pvalue: its row in the read table (set: the set the call came from) and
 * the target-call tests, on the read calls and (with the near-neighbor check) on the calls it keeps. checked:
 * ri_calls_neighbor_check's result, or -2 when not checked. */
void ri_calls_record(ri_calls_t *c, uint64_t read_id, const char *read_name, const ri_calls_set_t *set, int entry,
                     double p, uint32_t read_events, uint32_t n_slopes, int checked);

/* End of the run: n_reads every read of the input; writes the target call table. */
void ri_calls_finish(ri_calls_t *c, uint64_t n_reads);
void ri_calls_destroy(ri_calls_t *c);

/* log P[Bin(n, q) >= m] and log P[Bin(n, q) <= m]. */
double ri_log_binom_sf(uint64_t m, uint64_t n, double q);
double ri_log_binom_cdf(uint64_t m, uint64_t n, double q);

/* Target call at one cutoff: log LR of c >= 1 read calls in n reads at p0 = alpha_read. */
double ri_calls_log_lr(uint64_t c, uint64_t n, double p0);

/* Target call averaged over cutoffs: ln M for the read p-values p (sorted, each in [t_min, t_max]) after n reads. */
double ri_calls_log_mixture(const double *p, uint32_t n_p, uint64_t n, double t_min, double t_max);

#ifdef __cplusplus
}
#endif

#endif /* CALLS_H */
