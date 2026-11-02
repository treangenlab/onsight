#ifndef SEED_LOOKUP_H
#define SEED_LOOKUP_H

/* seed_lookup.h: each seed's positions in the index.
 *
 * The repeat cutoff max_occ, the (1 - --mid_occ_frac) quantile of positions per hash in the index bounded by --q-mid-occ
 * (ri_detopt_update), leaves out the seeds whose hash has more positions.
 *
 * In: a read's seeds, the index. Out: per seed found, its positions. State: none.
 */

#include <string.h> //for memset
#include "index.h"
#include "utils.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	uint32_t n;       // positions of the seed's hash in the index
	uint32_t q_pos;   // the seed's position in the read (as seeds.h stores it)
	uint32_t q_hash;  // 32-bit fingerprint of the hash
	const uint64_t *cr; // the positions
} ri_seed_t;

/**
 * @brief The read's seeds found in the index, without those whose hash has more than max_occ
 * positions (the repeat cutoff, ri_detopt_update).
 *
 * @return Pointer to the collected seeds (*_n_seed_m of them). Needs to be freed.
 */
ri_seed_t *ri_seed_lookup(void *km,
                              int *_n_seed_m,
                              int max_occ,
                              const ri_idx_t *ri,
                              const mm128_v *seeds);

#ifdef __cplusplus
}
#endif
#endif //SEED_LOOKUP_H
