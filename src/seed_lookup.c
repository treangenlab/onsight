/* seed_lookup.c: each seed's positions in the index; see seed_lookup.h. */
#include "seed_lookup.h"
#include "kalloc.h"
#include "seeds.h"
#include "utils.h"
#include "onsprof.h"

static inline uint32_t ri_hash_fingerprint32(uint64_t hash)
{
	return (uint32_t)hash ^ (uint32_t)(hash >> 32);
}

/* The read's seeds found in the index, *n_seed_hits_ of them. */
ri_seed_t *seed_lookup_all(void *km, const ri_idx_t *ri, const mm128_v *seeds, int32_t *n_seed_hits_)
{
	ri_seed_t *seed_hits;
	size_t i;
	int32_t n_seed_hits;

	seed_hits = (ri_seed_t*)ri_kmalloc(km, seeds->n * sizeof(ri_seed_t));
	for (i = n_seed_hits = 0; i < seeds->n; ++i) {
		const uint64_t *cr;
		ri_seed_t *q;
		mm128_t *p = &seeds->a[i];
		uint32_t q_pos = (uint32_t)p->y;
		uint64_t q_hash = p->x >> RI_HASH_SHIFT;
		int t;
		cr = ri_idx_get(ri, q_hash, &t);
		if (t == 0) continue;
		#ifdef ONS_PROF
		if (ons_prof_cur) {
			ons_prof_tls_t *pt = ons_prof_tls();
			ons_prof_cur->n_found++;
			if (t == 1) ons_prof_cur->n_single++;
			ons_prof_cur->post_found += (uint64_t)t;
			if (pt) pt->post[t < ONS_HIST_POST ? t : ONS_HIST_POST - 1]++;
		}
		#endif
		q = &seed_hits[n_seed_hits++];
		q->q_pos = q_pos, q->cr = cr, q->n = t;
		q->q_hash = ri_hash_fingerprint32(q_hash);
	}
	*n_seed_hits_ = n_seed_hits;
	return seed_hits;
}

/* The read's seeds below the repeat cutoff (seed_lookup.h). */
ri_seed_t *ri_seed_lookup(void *km,
                              int *_n_seed_m,
                              int max_occ,
                              const ri_idx_t *ri,
                              const mm128_v *seeds)
{
	int i, n_seed_m, n_seed_m0;
	ri_seed_t *seed_hits = seed_lookup_all(km, ri, seeds, &n_seed_m0);
	for (i = 0, n_seed_m = 0; i < n_seed_m0; ++i)
		if (seed_hits[i].n <= (uint32_t)max_occ) seed_hits[n_seed_m++] = seed_hits[i];
		#ifdef ONS_PROF
		else if (ons_prof_cur) { ons_prof_cur->n_flt_occ++; ons_prof_cur->post_flt_occ += seed_hits[i].n; }
		#endif
	*_n_seed_m = n_seed_m;
	return seed_hits;
}

