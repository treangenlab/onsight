#ifndef NEIGHBOR_CHECK_H
#define NEIGHBOR_CHECK_H

/* neighbor_check.h: the near-neighbor check, ported from the Python one (neighbor_check.py, archived 10-25).
 *
 * A read called as a panel entry is kept
 * only if its event levels over the call's band fit the entry's signature, or one of the entry's own strains, better
 * than every near neighbor's homologous region. The fit is the subsequence DTW cost against the pore model's expected
 * levels.
 *
 * The regions come from one file per panel (panel-src/near_neighbors/neighbor_regions.py build). The file is
 * memory-mapped read-only, so all threads share one db; each thread passes its own scratch buffer.
 *
 * In: a read call's signature, the box of its band (diagonal_bands.h) and the read's event levels. Out: kept or
 * dropped. State: the regions, read only. */

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
	int32_t near_genome;	/* that near neighbor, for nc_genome_name; -1 if none */
} nc_costs_t;

/* pore_model: the one the Python check used (R9 r94_means_only.txt, R10 uncalled_r1041_model_only_means.txt).
 * pad: margin around the band in events and k-mers (0 adopted; at most the file's --max-pad). NULL on error. */
nc_db_t *nc_load(const char *path, const char *pore_model, int pad);
void nc_free(nc_db_t *db);

/* Scratch bytes one thread needs; allocate once per thread and pass to every nc_keep call. */
size_t nc_buf_size(const nc_db_t *db);

/* 1 kept, 0 dropped, -1 no check (rid without a sequence, or an empty band); the Python scoring drops a call it
 * cannot check. levels: the read's event levels from the engine's segmentation, indexed by q_pos. With out NULL the
 * check stops at the first near neighbor that ties or beats the own side (about twice as fast); with out the costs
 * are exact. */
int nc_keep(const nc_db_t *db, uint32_t rid, const float *levels, uint32_t n_levels, const nc_band_t *band,
	    void *buf, nc_costs_t *out);

const char *nc_genome_name(const nc_db_t *db, int32_t genome);
uint32_t nc_sig_len(const nc_db_t *db, uint32_t rid);	/* 0 if none; compare with the index to catch a wrong file */

#ifdef __cplusplus
}
#endif
#endif
