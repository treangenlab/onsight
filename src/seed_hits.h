/* seed_hits.h: a read's seed hits, per signature.
 *
 * Each seed of the read found in the index (seed_lookup.h) gives one hit per index position: the signature, the
 * position in it and the strand, at the seed's read event. Per signature the read keeps its distinct hits, up to a cap
 * (--sig-report-anchor-cap) of the heaviest, each weighted by 1 / log2(positions of its hash + 1) with
 * --sig-weight-common-hashes; seeds over --sig-max-hash-occ positions add none. It also counts its distinct hits and
 * distinct seed hashes per signature.
 *
 * In: the read's seeds and the index. Out: per signature with hits, its ri_sig_hits_t, the hits in one block per read.
 * State: none.
 */
#ifndef SEED_HITS_H
#define SEED_HITS_H

#include <stdint.h>
#include "index.h"
#include "options.h"
#include "utils.h"

#ifdef __cplusplus
extern "C" {
#endif

/* One seed hit */
typedef struct ri_seed_hit_s{
	uint32_t q_pos;     // the seed's read event
	uint32_t ref_pos;   // its position in the signature
	float weight;       // 1 / log2(positions of its hash + 1) with --sig-weight-common-hashes, else 1
	uint8_t rev;        // 1: reverse strand
	uint8_t _pad[3];
} ri_seed_hit_t;

/* One read's seed hits on one signature */
typedef struct ri_sig_hits_s{
	uint32_t rid;                     // the signature
	uint32_t min_pos;                 // the span of its positions hit
	uint32_t max_pos;
	uint32_t unique_hits;             // distinct positions hit
	uint32_t distinct_query_hashes;   // distinct seed hashes that hit it
	double weighted_hits;             // the hits' weights summed
	ri_seed_hit_t *hits;              // the hits kept (at most --sig-report-anchor-cap)
	uint32_t n_hits;
	uint32_t m_hits;
} ri_sig_hits_t;

struct ri_tbuf_s;

/* One read: its events and its seed hits per signature */
typedef struct ri_read_s{
	uint32_t read_id;
	const char* read_name;

	uint32_t offset;    // events of the read
	float* levels;      // --neighbor-check: the read's event levels in the output step (in levels_buf from levels_off)
	uint32_t n_levels;
	uint64_t levels_off;
	struct ri_tbuf_s *levels_buf;   // the thread buffer that holds them (pipeline.c)

	// its seed hits, per signature
	ri_sig_hits_t *sig_hits;
	ri_seed_hit_t *hits_block; // the hits of all its signatures, one block (sig_hits[].hits point into it)
	uint32_t n_sig_hits;
	uint32_t m_sig_hits;
#ifdef ONS_PROF
	void *prof; // its ons_prof_read_t
#endif
} ri_read_t;

/* Adds the hits of the read's seeds (positions of each seed's hash, seeds above the repeat cutoff max_occ left out) */
void ri_collect_seed_hits(void *km, const ri_detopt_t *opt, int max_occ, const ri_idx_t *ri, ri_read_t *rd,
                          const mm128_v *seeds);


#ifdef __cplusplus
}
#endif
#endif // SEED_HITS_H
