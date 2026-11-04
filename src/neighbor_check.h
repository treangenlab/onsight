#ifndef NEIGHBOR_CHECK_H
#define NEIGHBOR_CHECK_H

/* neighbor_check.h: the near-neighbor check.
 *
 * A read called as a panel entry gets a second p-value, against near neighbors: how well its event levels over the
 * call's band fit the entry's own side (its signature and own strains) against every near neighbor's homologous region.
 * The fit is the subsequence DTW cost (mean squared difference) against the pore model's expected levels.
 *
 * Near-neighbor p-value. With own cost a over n read events, near cost b (the best near-neighbor region) and m own-side
 * regions, Gaussian level noise, the own side's likelihood averaged over its m regions and over a noise variance
 * log-uniform on [V0, V1], and the near side's at its best region and variance:
 *     ln L = (n/2) ln(2e b / (n a)) + lgamma(n/2) - ln m - ln ln(V1 / V0),    p_near = min(1, 1 / L).
 * For a read from any near-neighbor genome L has mean at most 1, so p_near is a p-value (Markov's inequality), and
 * max(p, p_near) one for reads that are unrelated or from a near neighbor. The target call counts the read at cutoff t
 * only if both p and p_near are at most t.
 *
 * The regions come from one file per panel (panel-src/near_neighbors/neighbor_regions.py build). The file is
 * memory-mapped read-only, so all threads share one db; each thread passes its own scratch buffer.
 *
 * In: a read call's signature, the box of its band (diagonal_bands.h) and the read's event levels. Out: the
 * near-neighbor p-value. State: the regions, read only. */

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct nc_db_s nc_db_t;

/* The engine's best diagonal band for the called signature, taken over all of the band's hits. */
typedef struct {
	uint32_t q0, q1;	/* smallest and largest q_pos */
	uint32_t r0, r1;	/* smallest and largest ref_pos (as the index stores it: k-mer index on the hit's strand) */
	int strand;		/* 1 if more than half of the hits are on the reverse strand */
} nc_band_t;

typedef struct {
	double own;		/* best cost of the signature and its own strains */
	double near;		/* best near-neighbor cost; HUGE_VAL if no near-neighbor region overlaps the band */
	uint32_t n, n_own;	/* read events in the band; own-side regions compared */
} nc_costs_t;

/* pore_model: the run's pore model (R9 r94_means_only.txt, R10 uncalled_r1041_model_only_means.txt).
 * pad: margin around the band in events and k-mers (the engine uses 0; at most the file's --max-pad). NULL on error. */
nc_db_t *nc_load(const char *path, const char *pore_model, int pad);
void nc_free(nc_db_t *db);

/* Scratch bytes one thread needs; allocate once per thread and pass to every nc_pvalue call. */
size_t nc_buf_size(const nc_db_t *db);

/* The near-neighbor p-value of a read call, 0 when no near-neighbor region overlaps the band, -1 when there is no check
 * (rid without a sequence, or an empty band; the call counts as dropped). levels: the read's event levels from the
 * engine's segmentation, indexed by q_pos. Only values in [lo, hi] are worked out (lo = max(p, t_min), hi = the
 * loosest cutoff): a region is followed only while it can still raise p_near above lo, and the check stops at the first
 * region that raises it above hi. lo then stands for "at most lo", and the value returned above hi is a lower bound;
 * max(p, p_near) is right for the target call either way. With out, every region is followed and p_near is exact. */
double nc_pvalue(const nc_db_t *db, uint32_t rid, const float *levels, uint32_t n_levels, const nc_band_t *band,
		 void *buf, double lo, double hi, nc_costs_t *out);

#ifdef __cplusplus
}
#endif
#endif
