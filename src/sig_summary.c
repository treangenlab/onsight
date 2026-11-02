/* sig_summary.c: each read's signatures ranked by score (--sig-summary) and the per-taxon table (--sig-profile);
 * see sig_summary.h. */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif
#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 700
#endif

#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "sig_summary.h"
#include "khash.h"

KHASH_MAP_INIT_STR(specprof, uint32_t)
KHASH_SET_INIT_INT(sigseen)

typedef struct {
	char *taxon;
	uint32_t panel_sig_count;
	uint32_t panel_locus_count;
	uint64_t panel_total_emitted_positions;
	uint64_t n_reads_best;
	uint64_t n_reads_top5;
	double sum_best_score;
	double sum_top5_score;
	void *seen_rids;
	void *seen_loci;
	void *panel_loci;
} ri_species_profile_t;

struct ri_profile_s {
	ri_species_profile_t *profiles;
	uint32_t n_profiles;
	uint32_t m_profiles;
	void *profile_index;
	uint32_t *profile_rid_species;  /* per signature: its taxon's profile (UINT32_MAX: none) */
};

static inline uint32_t sig_length(const ri_idx_t *ri, uint32_t rid)
{
	if (rid >= ri->n_seq) return 0;
	return ri->seq[rid].len;
}

static inline float sig_hits_cov(const ri_idx_t *ri, const ri_sig_hits_t *s)
{
	uint32_t slen = sig_length(ri, s->rid);
	if (slen == 0) return 0.0f;
	uint32_t span = s->max_pos - s->min_pos + 1;
	float pct = (float)span / (float)slen;
	return pct > 1.0f ? 1.0f : pct;
}

static inline uint32_t emitted_positions(const ri_idx_t *ri, uint32_t rid)
{
	if (!ri->sig_emitted_positions || rid >= ri->n_seq) return 0;
	return ri->sig_emitted_positions[rid];
}

static inline uint32_t distinct_hash_total(const ri_idx_t *ri, uint32_t rid)
{
	if (!ri->sig_distinct_hashes || rid >= ri->n_seq) return 0;
	return ri->sig_distinct_hashes[rid];
}

static inline float sig_hits_containment(const ri_idx_t *ri, const ri_sig_hits_t *s)
{
	uint32_t denom = emitted_positions(ri, s->rid);
	if (denom == 0) return 0.0f;
	return (float)s->unique_hits / (float)denom;
}

static inline float sig_hits_weighted_containment(const ri_idx_t *ri, const ri_sig_hits_t *s)
{
	uint32_t denom = emitted_positions(ri, s->rid);
	if (denom == 0) return 0.0f;
	return (float)(s->weighted_hits / (double)denom);
}

static inline float sig_hits_hash_containment(const ri_idx_t *ri, const ri_sig_hits_t *s)
{
	uint32_t denom = distinct_hash_total(ri, s->rid);
	if (denom == 0) return 0.0f;
	return (float)s->distinct_query_hashes / (float)denom;
}

static inline float sig_hits_fanout(const ri_sig_hits_t *s)
{
	return s->distinct_query_hashes == 0 ? 0.0f : (float)s->unique_hits / (float)s->distinct_query_hashes;
}

/* A signature's score: its hits over its emitted positions, each hit weighted by 1 / log2(postings + 1) with
 * --sig-weight-common-hashes. */
static inline float sig_hits_score(const ri_idx_t *ri, const ri_detopt_t *opt, const ri_sig_hits_t *s)
{
	if (opt->sig_weight_common_hashes) return sig_hits_weighted_containment(ri, s);
	return sig_hits_containment(ri, s);
}

/* qsort's comparator context, one per thread: the output step sorts the reads of a minibatch at once */
static __thread const ri_idx_t *g_rank_ri = NULL;
static __thread const ri_read_t *g_rank_rd = NULL;
static __thread const ri_detopt_t *g_rank_opt = NULL;
static int cmp_sig_hits_desc(const void *a, const void *b)
{
	uint32_t ia = *(const uint32_t*)a;
	uint32_t ib = *(const uint32_t*)b;
	const ri_sig_hits_t *sa = &g_rank_rd->sig_hits[ia];
	const ri_sig_hits_t *sb = &g_rank_rd->sig_hits[ib];
	float score_a = sig_hits_score(g_rank_ri, g_rank_opt, sa);
	float score_b = sig_hits_score(g_rank_ri, g_rank_opt, sb);
	if (score_a < score_b) return 1;
	if (score_a > score_b) return -1;
	{
		float hash_a = sig_hits_hash_containment(g_rank_ri, sa);
		float hash_b = sig_hits_hash_containment(g_rank_ri, sb);
		if (hash_a < hash_b) return 1;
		if (hash_a > hash_b) return -1;
		{
			float fanout_a = sig_hits_fanout(sa);
			float fanout_b = sig_hits_fanout(sb);
			if (fanout_a > fanout_b) return 1;
			if (fanout_a < fanout_b) return -1;
		}
		if (sa->distinct_query_hashes < sb->distinct_query_hashes) return 1;
		if (sa->distinct_query_hashes > sb->distinct_query_hashes) return -1;
	}
	if (sa->unique_hits < sb->unique_hits) return 1;
	if (sa->unique_hits > sb->unique_hits) return -1;
	{
		float cov_a = sig_hits_cov(g_rank_ri, sa);
		float cov_b = sig_hits_cov(g_rank_ri, sb);
		if (cov_a < cov_b) return 1;
		if (cov_a > cov_b) return -1;
	}
	if (sa->rid > sb->rid) return 1;
	if (sa->rid < sb->rid) return -1;
	return 0;
}

static void format_top_sigs(const ri_idx_t *ri, const ri_detopt_t *opt, const ri_read_t *rd,
	const uint32_t *ord, uint32_t limit, char *buf, size_t buf_size)
{
	uint32_t m;
	size_t pos = 0;
	if (!buf || buf_size == 0) return;
	buf[0] = '\0';
	for (m = 0; m < limit && m < rd->n_sig_hits; ++m) {
		const ri_sig_hits_t *sh = &rd->sig_hits[ord[m]];
		uint32_t rid = sh->rid;
		const char *name = (rid < ri->n_seq) ? ri->seq[rid].name : "*";
		float score = sig_hits_score(ri, opt, sh);
		int n = snprintf(buf + pos, buf_size - pos, "%s%s:%.6f",
			m ? "|" : "", name ? name : "*", score);
		if (n < 0 || (size_t)n >= buf_size - pos) break;
		pos += (size_t)n;
	}
}

/* printf into a new string */
static char *ri_strprintf(const char *fmt, ...)
{
	va_list ap, aq;
	int n;
	char *str;
	va_start(ap, fmt);
	va_copy(aq, ap);
	n = vsnprintf(NULL, 0, fmt, ap);
	va_end(ap);
	str = n >= 0 ? (char*)malloc((size_t)n + 1) : NULL;
	if (str) vsnprintf(str, (size_t)n + 1, fmt, aq);
	va_end(aq);
	return str;
}

void ri_sig_summary_read(const ri_idx_t *ri, const ri_detopt_t *opt, const ri_read_t *rd, int want_row,
                         ri_sig_summary_t *out)
{
	uint32_t *ord = (uint32_t*)malloc(rd->n_sig_hits * sizeof(uint32_t));
	char top5_buf[4096];
	for (uint32_t m = 0; m < rd->n_sig_hits; ++m) ord[m] = m;
	g_rank_ri = ri;
	g_rank_rd = rd;
	g_rank_opt = opt;
	qsort(ord, rd->n_sig_hits, sizeof(uint32_t), cmp_sig_hits_desc);
	const ri_sig_hits_t *best = &rd->sig_hits[ord[0]];
	float best_score = sig_hits_score(ri, opt, best);
	float second_score = rd->n_sig_hits > 1 ? sig_hits_score(ri, opt, &rd->sig_hits[ord[1]]) : 0.0f;
	float best_hash_score = sig_hits_hash_containment(ri, best);
	float best_fanout = sig_hits_fanout(best);
	uint32_t best_rid = best->rid;
	const char *best_name = (best_rid < ri->n_seq) ? ri->seq[best_rid].name : "*";
	const char *best_taxon = ri_sig_taxon(ri, best_rid);
	const char *call_name = best_name;
	const char *call_taxon = best_taxon;
	const char *call_reason = "ok";
	if (best->unique_hits < opt->sig_min_hits) {
		call_name = "*";
		call_taxon = "*";
		call_reason = "insufficient_evidence";
	} else if (best->distinct_query_hashes < opt->sig_min_distinct_hashes) {
		call_name = "*";
		call_taxon = "*";
		call_reason = "low_support_diversity";
	} else if (rd->n_sig_hits > 1 && second_score > 0.0f &&
		best_score / second_score < opt->sig_min_ratio) {
		call_name = "*";
		call_taxon = "*";
		call_reason = "ambiguous";
	}
	out->row = NULL;
	if (want_row) {
		format_top_sigs(ri, opt, rd, ord, 5, top5_buf, sizeof(top5_buf));
		out->row = ri_strprintf("%s\t%s\t%.6f\t%s\t%s\t%u\t%u\t%.6f\t%.6f\t%.6f\n",
			rd->read_name, call_name, best_score, call_reason,
			top5_buf[0] ? top5_buf : "*",
			best->unique_hits, best->distinct_query_hashes,
			best_hash_score, best_fanout, second_score);
	}
	out->ord = ord;
	out->call_taxon = call_taxon;
	out->best_score = best_score;
}

char *ri_sig_summary_empty_row(const char *read_name)
{
	return ri_strprintf("%s\t*\t0.000000\tno_candidates\t*\t0\t0\t0.000000\t0.000000\t0.000000\n", read_name);
}

static uint32_t profile_get_or_add_species(ri_profile_t *p, const char *taxon)
{
	khash_t(specprof) *idxh;
	khiter_t k;
	int absent = 0;
	if (!p || !taxon || taxon[0] == '\0' || strcmp(taxon, "*") == 0) return UINT32_MAX;
	if (!p->profile_index) p->profile_index = kh_init(specprof);
	idxh = (khash_t(specprof)*)p->profile_index;
	k = kh_get(specprof, idxh, taxon);
	if (k != kh_end(idxh)) return kh_val(idxh, k);
	if (p->n_profiles == p->m_profiles) {
		uint32_t new_m = p->m_profiles ? p->m_profiles << 1 : 16;
		p->profiles = (ri_species_profile_t*)realloc(p->profiles, new_m * sizeof(ri_species_profile_t));
		memset(p->profiles + p->m_profiles, 0, (new_m - p->m_profiles) * sizeof(ri_species_profile_t));
		p->m_profiles = new_m;
	}
	p->profiles[p->n_profiles].taxon = strdup(taxon);
	k = kh_put(specprof, idxh, p->profiles[p->n_profiles].taxon, &absent);
	kh_val(idxh, k) = p->n_profiles;
	return p->n_profiles++;
}

static void profile_note_signature(ri_species_profile_t *prof, uint32_t rid)
{
	khash_t(sigseen) *seen;
	int absent = 0;
	if (!prof) return;
	if (!prof->seen_rids) prof->seen_rids = kh_init(sigseen);
	seen = (khash_t(sigseen)*)prof->seen_rids;
	kh_put(sigseen, seen, (khint_t)rid, &absent);
}

static void profile_note_int_set(void **set_ptr, uint32_t value, uint32_t *count_ptr)
{
	khash_t(sigseen) *seen;
	int absent = 0;
	if (!set_ptr || value == UINT32_MAX) return;
	if (!*set_ptr) *set_ptr = kh_init(sigseen);
	seen = (khash_t(sigseen)*)*set_ptr;
	kh_put(sigseen, seen, (khint_t)value, &absent);
	if (absent && count_ptr) ++(*count_ptr);
}

static void profile_seed_from_index(ri_profile_t *p, const ri_idx_t *ri)
{
	uint32_t rid;
	if (!p || !ri) return;
	p->profile_rid_species = (uint32_t*)malloc(sizeof(uint32_t) * (ri->n_seq ? ri->n_seq : 1));
	for (rid = 0; rid < ri->n_seq; ++rid) {
		const char *taxon = ri_sig_taxon(ri, rid);
		uint32_t idx = profile_get_or_add_species(p, taxon);
		if (p->profile_rid_species) p->profile_rid_species[rid] = idx;
		if (idx == UINT32_MAX) continue;
		p->profiles[idx].panel_sig_count++;
		p->profiles[idx].panel_total_emitted_positions += ri->sig_emitted_positions ? ri->sig_emitted_positions[rid] : 0;
		profile_note_int_set(&p->profiles[idx].panel_loci, ri_sig_locus(ri, rid), &p->profiles[idx].panel_locus_count);
	}
}

static int cmp_species_profile_desc(const void *a, const void *b)
{
	const ri_species_profile_t *pa = *(const ri_species_profile_t * const *)a;
	const ri_species_profile_t *pb = *(const ri_species_profile_t * const *)b;
	if (pa->n_reads_best < pb->n_reads_best) return 1;
	if (pa->n_reads_best > pb->n_reads_best) return -1;
	if (pa->sum_best_score < pb->sum_best_score) return 1;
	if (pa->sum_best_score > pb->sum_best_score) return -1;
	if (pa->n_reads_top5 < pb->n_reads_top5) return 1;
	if (pa->n_reads_top5 > pb->n_reads_top5) return -1;
	return strcmp(pa->taxon ? pa->taxon : "*", pb->taxon ? pb->taxon : "*");
}

/* The profile of a signature's taxon: from the table profile_seed_from_index fills (every index taxon is in it) */
static inline uint32_t profile_species_of(ri_profile_t *p, const ri_idx_t *ri, uint32_t rid)
{
	if (p->profile_rid_species && rid < ri->n_seq) return p->profile_rid_species[rid];
	return profile_get_or_add_species(p, ri_sig_taxon(ri, rid));
}

void ri_profile_add_read(ri_profile_t *p, const ri_idx_t *ri, const ri_detopt_t *opt, const ri_read_t *rd,
	const uint32_t *ord, const char *call_name, float best_score)
{
	uint32_t m, top_species_idx[5], top_species_n = 0;
	float top_species_score[5] = {0};
	if (!p || !ri || !rd) return;
	for (m = 0; m < rd->n_sig_hits; ++m) {
		const ri_sig_hits_t *sh = &rd->sig_hits[ord[m]];
		uint32_t idx = profile_species_of(p, ri, sh->rid);
		if (idx == UINT32_MAX) continue;
		profile_note_signature(&p->profiles[idx], sh->rid);
		profile_note_int_set(&p->profiles[idx].seen_loci, ri_sig_locus(ri, sh->rid), NULL);
	}
	if (call_name && strcmp(call_name, "*") != 0) {
		uint32_t best_idx = profile_get_or_add_species(p, call_name);
		if (best_idx != UINT32_MAX) {
			p->profiles[best_idx].n_reads_best++;
			p->profiles[best_idx].sum_best_score += best_score;
		}
	}
	for (m = 0; m < rd->n_sig_hits && m < 5; ++m) {
		const ri_sig_hits_t *sh = &rd->sig_hits[ord[m]];
		uint32_t idx = profile_species_of(p, ri, sh->rid);
		float score = sig_hits_score(ri, opt, sh);
		uint32_t j;
		if (idx == UINT32_MAX) continue;
		for (j = 0; j < top_species_n; ++j) {
			if (top_species_idx[j] == idx) {
				if (score > top_species_score[j]) top_species_score[j] = score;
				break;
			}
		}
		if (j == top_species_n && top_species_n < 5) {
			top_species_idx[top_species_n] = idx;
			top_species_score[top_species_n] = score;
			++top_species_n;
		}
	}
	for (m = 0; m < top_species_n; ++m) {
		p->profiles[top_species_idx[m]].n_reads_top5++;
		p->profiles[top_species_idx[m]].sum_top5_score += top_species_score[m];
	}
}

void ri_profile_write(const ri_profile_t *p, FILE *fp)
{
	ri_species_profile_t **ord;
	uint32_t i;
	if (!p || !fp) return;
	fprintf(fp, "source_taxon\tpanel_sig_count\tpanel_total_emitted_positions\tn_unique_signatures_hit\tn_reads_best\tsum_best_score\tn_reads_top5\tsum_top5_score\tpanel_locus_count\tn_unique_loci_hit\tlocus_breadth\n");
	if (!p->n_profiles) return;
	ord = (ri_species_profile_t**)calloc(p->n_profiles, sizeof(ri_species_profile_t*));
	for (i = 0; i < p->n_profiles; ++i) ord[i] = &p->profiles[i];
	qsort(ord, p->n_profiles, sizeof(ri_species_profile_t*), cmp_species_profile_desc);
	for (i = 0; i < p->n_profiles; ++i) {
		khash_t(sigseen) *seen = ord[i]->seen_rids ? (khash_t(sigseen)*)ord[i]->seen_rids : NULL;
		khash_t(sigseen) *seen_loci = ord[i]->seen_loci ? (khash_t(sigseen)*)ord[i]->seen_loci : NULL;
		uint32_t n_unique_loci_hit = seen_loci ? (uint32_t)kh_size(seen_loci) : 0;
		double locus_breadth = ord[i]->panel_locus_count ? (double)n_unique_loci_hit / (double)ord[i]->panel_locus_count : 0.0;
		fprintf(fp, "%s\t%u\t%llu\t%u\t%llu\t%.6f\t%llu\t%.6f\t%u\t%u\t%.6f\n",
			ord[i]->taxon ? ord[i]->taxon : "*",
			ord[i]->panel_sig_count,
			(unsigned long long)ord[i]->panel_total_emitted_positions,
			seen ? (uint32_t)kh_size(seen) : 0,
			(unsigned long long)ord[i]->n_reads_best,
			ord[i]->sum_best_score,
			(unsigned long long)ord[i]->n_reads_top5,
			ord[i]->sum_top5_score,
			ord[i]->panel_locus_count,
			n_unique_loci_hit,
			locus_breadth);
	}
	free(ord);
}

void ri_profile_destroy(ri_profile_t *p)
{
	uint32_t i;
	if (!p) return;
	if (p->profile_index) kh_destroy(specprof, (khash_t(specprof)*)p->profile_index);
	for (i = 0; i < p->n_profiles; ++i) {
		if (p->profiles[i].taxon) free(p->profiles[i].taxon);
		if (p->profiles[i].seen_rids) kh_destroy(sigseen, (khash_t(sigseen)*)p->profiles[i].seen_rids);
		if (p->profiles[i].seen_loci) kh_destroy(sigseen, (khash_t(sigseen)*)p->profiles[i].seen_loci);
		if (p->profiles[i].panel_loci) kh_destroy(sigseen, (khash_t(sigseen)*)p->profiles[i].panel_loci);
	}
	free(p->profiles);
	free(p->profile_rid_species);
	free(p);
}

ri_profile_t *ri_profile_init(const ri_idx_t *ri)
{
	ri_profile_t *p = (ri_profile_t*)calloc(1, sizeof(ri_profile_t));
	if (p) profile_seed_from_index(p, ri);
	return p;
}
