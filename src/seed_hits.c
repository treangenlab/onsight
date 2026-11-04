/* seed_hits.c: a read's seed hits, per signature; see seed_hits.h. */
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "seed_hits.h"
#include "seed_lookup.h"
#include "seeds.h"
#include "kalloc.h"
#include "onsprof.h"

/* used only when opt->sig_hits_cap is 0; the option's default is the same */
#define RI_SIG_HITS_CAP_DEFAULT 1024

/* One signature position met by one of the read's seeds, kept per signature in the order the seeds met them */
typedef struct {
	uint32_t ref_pos, q_pos, q_hash;
	float    weight;
	uint8_t  rev;
} ri_occ_t;

/* A signature position's best occurrence so far and its slot in the hits (UINT32_MAX: evicted at the cap) */
typedef struct {
	uint32_t ref_pos, q_pos, hit_idx;
	float    weight;
	uint8_t  rev, used;
} ri_pos_slot_t;

/* Tables over one signature's positions and seed hashes, sized for the read's largest signature */
typedef struct {
	ri_pos_slot_t *pt;
	uint32_t      *hk;
	uint8_t       *hu;
} ri_occ_tabs_t;

static void sig_add_occs(const ri_detopt_t *opt, ri_sig_hits_t *s, const ri_occ_t *o, uint32_t m, ri_occ_tabs_t *t);

/* Adds a read's seed hits, per signature: every index position of every seed found in the index, but for the seeds
 * above the repeat cutoff max_occ (ri_seed_lookup). The positions are first grouped by signature, in the order the
 * seeds met them, then each signature's are added one by one (sig_add_occs): the same hits, order and sums as adding
 * every position in seed order, without hash tables over the whole read. One call per read; km: the thread's
 * memory pool (NULL: malloc). */
void ri_collect_seed_hits(void *km,
								  const ri_detopt_t *opt,
								  int max_occ,
								  const ri_idx_t *ri,
								  ri_read_t* rd,
								  const mm128_v *seeds)
{
	int i, n_seed_m;
	ri_seed_t *seed_hits0;
	uint32_t mask_pos = (1ULL<<31)-1, *cnt, *order, *start, n_order = 0, n_occ = 0, j, k, m;
	uint32_t hits_cap = (opt && opt->sig_hits_cap > 0) ? (uint32_t)opt->sig_hits_cap : (uint32_t)RI_SIG_HITS_CAP_DEFAULT;
	uint64_t n_pos = 0;
	float *weight;
	ri_occ_t *occ;

	#ifdef ONS_PROF
	uint64_t pt0 = ons_prof_cur ? ons_now() : 0;
	#endif
	seed_hits0 = ri_seed_lookup(km, &n_seed_m, max_occ, ri, seeds);
	#ifdef ONS_PROF
	if (ons_prof_cur) { uint64_t pt1 = ons_now(); ons_prof_cur->t_lookup += pt1 - pt0; pt0 = pt1; }
	#endif
	/* each seed's weight; 0: the seed adds no hits */
	weight = (float*)ri_kmalloc(km, (n_seed_m > 0 ? n_seed_m : 1) * sizeof(float));
	for (i = 0; i < n_seed_m; ++i) {
		ri_seed_t *s_match = &seed_hits0[i];
		float q_weight = 1.0f;
		/* --sig-max-hash-occ: without it, very common hashes (thousands of positions) fill a signature's hits before the
		 * rarer hashes land. */
		int sig_skip_admission = 0;
		if (opt->sig_max_hash_occ > 0 && s_match->n > opt->sig_max_hash_occ) {
			sig_skip_admission = 1;
		}
		if (opt->sig_weight_common_hashes) {
			if (sig_skip_admission)
				q_weight = 0.0f;
			else
				q_weight = 1.0f / log2f((float)s_match->n + 1.0f);
		}
		if (sig_skip_admission || !(q_weight > 0.0f)) {
			#ifdef ONS_PROF
			if (ons_prof_cur) { ons_prof_cur->n_skip_adm++; ons_prof_cur->post_skip_adm += s_match->n; }
			#endif
			weight[i] = 0.0f;
			continue;
		}
		#ifdef ONS_PROF
		if (ons_prof_cur) ons_prof_cur->post_visited += s_match->n;
		#endif
		weight[i] = q_weight;
		n_pos += s_match->n;
	}

	/* the positions per signature, signatures in the order the seeds first met them */
	cnt = (uint32_t*)ri_kcalloc(km, ri->n_seq + 1, sizeof(uint32_t));
	order = (uint32_t*)ri_kmalloc(km, (n_pos < ri->n_seq ? n_pos : ri->n_seq) * sizeof(uint32_t) + sizeof(uint32_t));
	for (i = 0; i < n_seed_m; ++i) {
		const uint64_t *hits = seed_hits0[i].cr;
		if (!(weight[i] > 0.0f)) continue;
		for (k = 0; k < seed_hits0[i].n; ++k) {
			uint64_t rid = hits[k] >> RI_ID_SHIFT;
			if (rid >= ri->n_seq) continue;
			if (cnt[rid]++ == 0) order[n_order++] = (uint32_t)rid;
			++n_occ;
		}
	}
	start = (uint32_t*)ri_kmalloc(km, (n_order + 1) * sizeof(uint32_t));
	for (j = 0, m = 0; j < n_order; ++j) {   /* cnt becomes each signature's next free slot */
		start[j] = m;
		m += cnt[order[j]];
		cnt[order[j]] = start[j];
	}
	start[n_order] = m;
	occ = (ri_occ_t*)ri_kmalloc(km, (n_occ > 0 ? n_occ : 1) * sizeof(ri_occ_t));
	for (i = 0; i < n_seed_m; ++i) {
		ri_seed_t *s_match = &seed_hits0[i];
		const uint64_t *hits = s_match->cr;
		uint32_t q_pos = (uint32_t)((s_match->q_pos >> RI_POS_SHIFT) + rd->offset);
		if (!(weight[i] > 0.0f)) continue;
		for (k = 0; k < s_match->n; ++k) {
			uint64_t rid = hits[k] >> RI_ID_SHIFT;
			ri_occ_t *e;
			if (rid >= ri->n_seq) continue;
			e = &occ[cnt[rid]++];
			e->ref_pos = (uint32_t)(hits[k]>>RI_POS_SHIFT)&mask_pos;
			e->q_pos = q_pos;
			e->q_hash = s_match->q_hash;
			e->weight = weight[i];
			e->rev = (hits[k] & 1) ? 1 : 0;
		}
	}

	/* the signatures' hits, in one block: each signature has room for its occurrences, at most the hit cap */
	{
		uint64_t tot = 0;
		for (j = 0; j < n_order; ++j) {
			uint32_t c = start[j + 1] - start[j];
			tot += c < hits_cap ? c : hits_cap;
		}
		rd->hits_block = tot ? (ri_seed_hit_t*)malloc(tot * sizeof(ri_seed_hit_t)) : NULL;
	}
	rd->n_sig_hits = rd->m_sig_hits = n_order;
	if (n_order > 0) rd->sig_hits = (ri_sig_hits_t*)malloc(n_order * sizeof(ri_sig_hits_t));
	{
		ri_occ_tabs_t t;
		uint32_t tab = 4, m_max = 0;
		for (j = 0; j < n_order; ++j)
			if (start[j + 1] - start[j] > m_max) m_max = start[j + 1] - start[j];
		while (tab < 2 * m_max) tab <<= 1;
		t.pt = (ri_pos_slot_t*)ri_kmalloc(km, tab * sizeof(ri_pos_slot_t));
		t.hk = (uint32_t*)ri_kmalloc(km, tab * sizeof(uint32_t));
		t.hu = (uint8_t*)ri_kmalloc(km, tab);
		uint64_t off = 0;
		for (j = 0; j < n_order; ++j) {
			ri_sig_hits_t *s = &rd->sig_hits[j];
			uint32_t c = start[j + 1] - start[j];
			memset(s, 0, sizeof(*s));
			s->rid = order[j];
			s->hits = rd->hits_block + off;
			s->m_hits = c < hits_cap ? c : hits_cap;
			off += s->m_hits;
			sig_add_occs(opt, s, occ + start[j], c, &t);
		}
		ri_kfree(km, t.hu); ri_kfree(km, t.hk); ri_kfree(km, t.pt);
	}
	ri_kfree(km, occ); ri_kfree(km, start); ri_kfree(km, order); ri_kfree(km, cnt); ri_kfree(km, weight);
	if(seed_hits0){ri_kfree(km, seed_hits0); seed_hits0 = NULL;}
	#ifdef ONS_PROF
	if (ons_prof_cur) ons_prof_cur->t_accum += ons_now() - pt0;
	#endif
}

/* The signature's lightest hit: its slot in *out_idx and its weight in *out_w.
 * Caller guarantees s->n_hits > 0. Linear scan; n_hits <= cap
 * (default 1024) so this is cheap. */
static inline void sig_hits_find_min(const ri_sig_hits_t *s, uint32_t *out_idx, float *out_w)
{
	uint32_t min_idx = 0;
	float min_w = s->hits[0].weight;
	uint32_t j;
	for (j = 1; j < s->n_hits; ++j) {
		if (s->hits[j].weight < min_w) {
			min_w = s->hits[j].weight;
			min_idx = j;
		}
	}
	*out_idx = min_idx;
	*out_w = min_w;
}

static inline uint32_t sig_hits_add(ri_sig_hits_t *s, uint32_t q_pos, uint32_t ref_pos, uint8_t rev, float q_weight, uint32_t cap)
{
	uint32_t idx;
	if (!s) return UINT32_MAX;
	if (cap == 0) cap = RI_SIG_HITS_CAP_DEFAULT; /* fallback if opt unset */
	/* Caller (add_seed_hit) handles the cap-full eviction case so it
	 * can also patch the sigpos hash entry that owns the evicted slot.
	 * If we ever reach the cap here it means caller did not pre-check;
	 * fall back to FCFS-reject to stay safe. */
	if (s->n_hits >= cap) return UINT32_MAX;
	if (s->n_hits == s->m_hits) return UINT32_MAX;   /* not reached: the block has room for every occurrence */
	idx = s->n_hits++;
	s->hits[idx].q_pos = q_pos;
	s->hits[idx].ref_pos = ref_pos;
	s->hits[idx].weight = q_weight;
	s->hits[idx].rev = rev;
	return idx;
}

static inline void sig_hits_update(ri_sig_hits_t *s, uint32_t hit_idx, uint32_t q_pos, uint32_t ref_pos, uint8_t rev, float q_weight)
{
	if (!s || !s->hits || hit_idx >= s->n_hits) return;
	s->hits[hit_idx].q_pos = q_pos;
	s->hits[hit_idx].ref_pos = ref_pos;
	s->hits[hit_idx].weight = q_weight;
	s->hits[hit_idx].rev = rev;
}

/* One signature's occurrences o[0..m), in seed order, added one by one: a new position becomes a hit (at the cap it
 * replaces the lightest hit if heavier); a position met again keeps its heaviest occurrence, ties to the earlier read
 * position. Tables over this signature's positions and seed hashes stand in for the read's hash tables. */
static void sig_add_occs(const ri_detopt_t *opt, ri_sig_hits_t *s, const ri_occ_t *o, uint32_t m, ri_occ_tabs_t *t)
{
	uint32_t hits_cap = (opt && opt->sig_hits_cap > 0)
		? (uint32_t)opt->sig_hits_cap
		: (uint32_t)RI_SIG_HITS_CAP_DEFAULT;
	uint32_t i, tab = 4, bits = 2, mask;
	ri_pos_slot_t *pt = t->pt;
	uint32_t *hk = t->hk;
	uint8_t *hu = t->hu;

	while (tab < 2 * m) tab <<= 1, ++bits;
	mask = tab - 1;
	memset(pt, 0, tab * sizeof(ri_pos_slot_t));
	memset(hu, 0, tab);
	for (i = 0; i < m; ++i) {
		const ri_occ_t *e = &o[i];
		uint32_t ref_pos = e->ref_pos, q_pos = e->q_pos, h;
		uint8_t rev = e->rev;
		float q_weight = e->weight;
		ri_pos_slot_t *info;

		h = (ref_pos * 0x9E3779B1u) >> (32 - bits);
		while (pt[h].used && pt[h].ref_pos != ref_pos) h = (h + 1) & mask;
		info = &pt[h];
		#ifdef ONS_PROF
		if (ons_prof_cur && !info->used) ons_prof_cur->n_pos_new++;
		#endif
		if (!info->used) {
			uint32_t hit_idx;
			info->used = 1;
			info->ref_pos = ref_pos;
			info->weight = q_weight;
			info->q_pos = q_pos;
			info->rev = rev;
			hit_idx = sig_hits_add(s, q_pos, ref_pos, rev, q_weight, hits_cap);
			/* At the cap, the incoming hit replaces the lightest one if it is heavier. The replaced hit's entry in the
			 * position table then points nowhere, or it would alias the new hit's slot. The totals (unique_hits,
			 * weighted_hits) count every distinct (rid, ref_pos) hit; the cap only chooses which hits are kept. */
			if (hit_idx == UINT32_MAX && s->n_hits > 0 && s->n_hits >= hits_cap) {
				uint32_t min_idx; float min_w;
				sig_hits_find_min(s, &min_idx, &min_w);
				#ifdef ONS_PROF
				if (ons_prof_cur) { ons_prof_cur->n_cap_scan++; if (q_weight > min_w) ons_prof_cur->n_cap_evict++; }
				#endif
				if (q_weight > min_w) {
					uint32_t old_pos = s->hits[min_idx].ref_pos, g = (old_pos * 0x9E3779B1u) >> (32 - bits);
					while (pt[g].used && pt[g].ref_pos != old_pos) g = (g + 1) & mask;
					if (pt[g].used) pt[g].hit_idx = UINT32_MAX;
					s->hits[min_idx].q_pos = q_pos;
					s->hits[min_idx].ref_pos = ref_pos;
					s->hits[min_idx].weight = q_weight;
					s->hits[min_idx].rev = rev;
					hit_idx = min_idx;
				}
			}
			info->hit_idx = hit_idx;
			if (s->unique_hits == 0) {
				s->min_pos = ref_pos;
				s->max_pos = ref_pos;
			} else {
				if (ref_pos < s->min_pos) s->min_pos = ref_pos;
				if (ref_pos > s->max_pos) s->max_pos = ref_pos;
			}
			s->unique_hits++;
			s->weighted_hits += q_weight;
		} else {
			if (q_weight > info->weight || (q_weight == info->weight && q_pos < info->q_pos)) {
				s->weighted_hits += (double)q_weight - (double)info->weight;
				info->weight = q_weight;
				info->q_pos = q_pos;
				info->rev = rev;
				if (info->hit_idx != UINT32_MAX)
					sig_hits_update(s, info->hit_idx, q_pos, ref_pos, rev, q_weight);
			}
		}

		h = (e->q_hash * 0x9E3779B1u) >> (32 - bits);
		while (hu[h] && hk[h] != e->q_hash) h = (h + 1) & mask;
		if (!hu[h]) { hu[h] = 1; hk[h] = e->q_hash; s->distinct_query_hashes++; }
	}
}
