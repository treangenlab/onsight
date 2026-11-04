/* index.c: the index of the panel's signatures and the signature table; see index.h. */
/* feature-test macros for strdup */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif

#include "index.h"
#include <assert.h>
#include <fcntl.h>
#include "seeds.h"
#include "bseq.h"
#include "khash.h"
#include "rh_kvec.h"
#include "kthread.h"
#include "event_detection.h"
#include "onsprof.h"
#include <string.h>
#include <ctype.h>

#if defined(WIN32) || defined(_WIN32)
#include <io.h> // for open(2)
#else
#include <unistd.h>
#endif

#define idx_hash(a) ((a)>>1)
#define idx_eq(a, b) ((a)>>1 == (b)>>1)
KHASH_INIT(idx, uint64_t, uint64_t, 1, idx_hash, idx_eq)
typedef khash_t(idx) idxhash_t;

KHASH_MAP_INIT_STR(str, uint32_t)
KHASH_SET_INIT_INT64(sigmeta_pos)
KHASH_SET_INIT_INT64(sigmeta_hash)

#define RI_SIGMETA_LOCUS_GAP_BP 1000U

#define kroundup64(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, (x)|=(x)>>32, ++(x))
#define mm_seq4_set(s, i, c) ((s)[(i)>>3] |= (uint32_t)(c) << (((i)&7)<<2))

double ri_realtime0;
int ri_verbose = 1;

typedef struct {
	int mini_batch_size;
	uint64_t batch_size, sum_len;
	mm_bseq_file_t* fp;
	ri_idx_t* ri;
} pipeline_t;

typedef struct {
    int n_seq;
	mm_bseq1_t* seq;
	mm128_v a;
} step_t;

typedef struct {
	uint32_t rid;
	const char *taxon;
	uint32_t assembly_id;
	const char *contig_header;
	uint32_t start_bp;
	uint32_t end_bp;
	int parsed;
} ri_sig_locus_rec_t;

static char *ri_sigmeta_path(const char *path)
{
	size_t n;
	char *meta_path;
	if (!path) return NULL;
	n = strlen(path);
	meta_path = (char *)malloc(n + 13);
	if (!meta_path) return NULL;
	snprintf(meta_path, n + 13, "%s.sigmeta.tsv", path);
	return meta_path;
}

static char *ri_strdup0(const char *s)
{
	size_t n;
	char *out;
	if (!s) return NULL;
	n = strlen(s);
	out = (char *)malloc(n + 1);
	if (!out) return NULL;
	memcpy(out, s, n + 1);
	return out;
}

static char *ri_trim_inplace(char *s)
{
	char *end;
	while (*s && isspace((unsigned char)*s)) ++s;
	end = s + strlen(s);
	while (end > s && isspace((unsigned char)end[-1])) --end;
	*end = '\0';
	return s;
}

static const char *ri_sig_name(const ri_idx_t *ri, uint32_t rid)
{
	if (!ri || rid >= ri->n_seq) return NULL;
	return ri->seq[rid].name;
}

static int ri_csv_col_index(const char *header, const char *target)
{
	int idx = 0;
	const char *p = header;
	const char *field = p;
	while (1) {
		if (*p == ',' || *p == '\0' || *p == '\n' || *p == '\r') {
			size_t len = (size_t)(p - field);
			if (strlen(target) == len && strncmp(field, target, len) == 0) return idx;
			if (*p == '\0' || *p == '\n' || *p == '\r') break;
			field = p + 1;
			++idx;
		}
		++p;
	}
	return -1;
}

static int ri_csv_get_field(const char *line, int target_idx, char *buf, size_t buf_size)
{
	int idx = 0;
	const char *p = line;
	const char *field = p;
	if (!buf || buf_size == 0 || target_idx < 0) return -1;
	buf[0] = '\0';
	while (1) {
		if (*p == ',' || *p == '\0' || *p == '\n' || *p == '\r') {
			if (idx == target_idx) {
				size_t len = (size_t)(p - field);
				if (len >= buf_size) len = buf_size - 1;
				memcpy(buf, field, len);
				buf[len] = '\0';
				return 0;
			}
			if (*p == '\0' || *p == '\n' || *p == '\r') break;
			field = p + 1;
			++idx;
		}
		++p;
	}
	return -1;
}

static void ri_idx_apply_metadata_csv(const char *metadata_path, ri_idx_t *ri)
{
	FILE *fp;
	char line[65536];
	int name_col, tax_col;
	khash_t(str) *name_to_rid;
	uint32_t i;
	if (!metadata_path || !ri || ri->n_seq == 0) return;
	fp = fopen(metadata_path, "r");
	if (!fp) {
		fprintf(stderr, "[WARNING] failed to open metadata CSV '%s'\n", metadata_path);
		return;
	}
	if (!ri->sig_source_taxon)
		ri->sig_source_taxon = (char **)ri_kcalloc(ri->km, ri->n_seq, sizeof(char *));
	name_to_rid = kh_init(str);
	for (i = 0; i < ri->n_seq; ++i) {
		const char *name = ri->seq[i].name;
		int absent = 0;
		if (name) {
			khiter_t k = kh_put(str, name_to_rid, (char*)name, &absent);
			kh_val(name_to_rid, k) = i;
		}
	}
	if (!fgets(line, sizeof(line), fp)) {
		kh_destroy(str, name_to_rid);
		fclose(fp);
		return;
	}
	name_col = ri_csv_col_index(line, "fasta_header");
	tax_col = ri_csv_col_index(line, "source_taxon");
	if (name_col < 0 || tax_col < 0) {
		fprintf(stderr, "[WARNING] metadata CSV '%s' is missing fasta_header/source_taxon columns\n", metadata_path);
		kh_destroy(str, name_to_rid);
		fclose(fp);
		return;
	}
	while (fgets(line, sizeof(line), fp)) {
		char name_buf[4096], tax_buf[4096];
		char *name_trim, *tax_trim;
		if (ri_csv_get_field(line, name_col, name_buf, sizeof(name_buf)) != 0) continue;
		if (ri_csv_get_field(line, tax_col, tax_buf, sizeof(tax_buf)) != 0) continue;
		name_trim = ri_trim_inplace(name_buf);
		tax_trim = ri_trim_inplace(tax_buf);
		if (name_trim[0] == '\0' || tax_trim[0] == '\0') continue;
		{
			khiter_t k = kh_get(str, name_to_rid, name_trim);
			if (k != kh_end(name_to_rid)) {
				uint32_t rid = kh_val(name_to_rid, k);
				if (rid < ri->n_seq) {
					if (ri->sig_source_taxon[rid]) free(ri->sig_source_taxon[rid]);
					ri->sig_source_taxon[rid] = ri_strdup0(tax_trim);
				}
			}
		}
	}
	for (i = 0; i < ri->n_seq; ++i) {
		if (!ri->sig_source_taxon[i]) ri->sig_source_taxon[i] = ri_strdup0("*");
	}
	kh_destroy(str, name_to_rid);
	fclose(fp);
}

static void ri_idx_capture_sigmeta(ri_idx_t *ri, uint32_t rid, const mm128_t *a, size_t start, size_t end)
{
	khash_t(sigmeta_pos) *posh;
	khash_t(sigmeta_hash) *hashh;
	size_t i;
	if (!ri->sig_emitted_positions || !ri->sig_distinct_hashes || start >= end) return;
	posh = kh_init(sigmeta_pos);
	hashh = kh_init(sigmeta_hash);
	for (i = start; i < end; ++i) {
		uint64_t pos_key = ((uint64_t)rid << 32) | (uint64_t)((a[i].y >> RI_POS_SHIFT) & ((1ULL << 31) - 1));
		uint64_t hash_key = a[i].x >> RI_HASH_SHIFT;
		int absent = 0;
		kh_put(sigmeta_pos, posh, (khint64_t)pos_key, &absent);
		kh_put(sigmeta_hash, hashh, (khint64_t)hash_key, &absent);
	}
	ri->sig_emitted_positions[rid] = (uint32_t)kh_size(posh);
	ri->sig_distinct_hashes[rid] = (uint32_t)kh_size(hashh);
	kh_destroy(sigmeta_pos, posh);
	kh_destroy(sigmeta_hash, hashh);
}

static void ri_idx_clear_locus_meta(ri_idx_t *ri)
{
	uint32_t i;
	if (!ri) return;
	if (ri->sig_contig_header) {
		for (i = 0; i < ri->n_seq; ++i) {
			if (ri->sig_contig_header[i]) free(ri->sig_contig_header[i]);
		}
		ri_kfree(ri->km, ri->sig_contig_header);
		ri->sig_contig_header = NULL;
	}
	if (ri->sig_assembly_id) {
		ri_kfree(ri->km, ri->sig_assembly_id);
		ri->sig_assembly_id = NULL;
	}
	if (ri->sig_start_bp) {
		ri_kfree(ri->km, ri->sig_start_bp);
		ri->sig_start_bp = NULL;
	}
	if (ri->sig_end_bp) {
		ri_kfree(ri->km, ri->sig_end_bp);
		ri->sig_end_bp = NULL;
	}
	if (ri->sig_locus_index) {
		ri_kfree(ri->km, ri->sig_locus_index);
		ri->sig_locus_index = NULL;
	}
	if (ri->sig_locus_start_bp) {
		ri_kfree(ri->km, ri->sig_locus_start_bp);
		ri->sig_locus_start_bp = NULL;
	}
	if (ri->sig_locus_end_bp) {
		ri_kfree(ri->km, ri->sig_locus_end_bp);
		ri->sig_locus_end_bp = NULL;
	}
}

static int ri_parse_sig_interval(const char *name, uint32_t *assembly_id, char **contig_header, uint32_t *start_bp, uint32_t *end_bp)
{
	const char *first_dash, *last_dash, *colon;
	char *endptr = NULL;
	unsigned long assembly_ul, start_ul, end_ul;
	size_t contig_len;
	char *contig_dup;
	if (!name || !assembly_id || !contig_header || !start_bp || !end_bp) return 0;
	first_dash = strchr(name, '-');
	last_dash = strrchr(name, '-');
	if (!first_dash || !last_dash || first_dash == last_dash) return 0;
	colon = strchr(last_dash + 1, ':');
	if (!colon) return 0;

	assembly_ul = strtoul(name, &endptr, 10);
	if (endptr != first_dash || assembly_ul > UINT32_MAX) return 0;

	start_ul = strtoul(last_dash + 1, &endptr, 10);
	if (endptr != colon || start_ul > UINT32_MAX) return 0;

	end_ul = strtoul(colon + 1, &endptr, 10);
	if ((*endptr != '\0' && *endptr != '\n' && *endptr != '\r') || end_ul > UINT32_MAX) return 0;

	contig_len = (size_t)(last_dash - first_dash - 1);
	if (contig_len == 0) return 0;
	contig_dup = (char *)malloc(contig_len + 1);
	if (!contig_dup) return 0;
	memcpy(contig_dup, first_dash + 1, contig_len);
	contig_dup[contig_len] = '\0';

	*assembly_id = (uint32_t)assembly_ul;
	*contig_header = contig_dup;
	*start_bp = (uint32_t)start_ul;
	*end_bp = (uint32_t)end_ul;
	if (*end_bp < *start_bp) {
		uint32_t tmp = *start_bp;
		*start_bp = *end_bp;
		*end_bp = tmp;
	}
	return 1;
}

static int ri_sig_locus_rec_cmp(const void *a, const void *b)
{
	const ri_sig_locus_rec_t *ra = (const ri_sig_locus_rec_t *)a;
	const ri_sig_locus_rec_t *rb = (const ri_sig_locus_rec_t *)b;
	const char *taxa = (ra->taxon && ra->taxon[0]) ? ra->taxon : "*";
	const char *taxb = (rb->taxon && rb->taxon[0]) ? rb->taxon : "*";
	int c = strcmp(taxa, taxb);
	if (c != 0) return c;
	if (ra->parsed != rb->parsed) return rb->parsed - ra->parsed;
	if (ra->assembly_id < rb->assembly_id) return -1;
	if (ra->assembly_id > rb->assembly_id) return 1;
	{
		const char *contiga = ra->contig_header ? ra->contig_header : "*";
		const char *contigb = rb->contig_header ? rb->contig_header : "*";
		c = strcmp(contiga, contigb);
		if (c != 0) return c;
	}
	if (ra->start_bp < rb->start_bp) return -1;
	if (ra->start_bp > rb->start_bp) return 1;
	if (ra->end_bp < rb->end_bp) return -1;
	if (ra->end_bp > rb->end_bp) return 1;
	if (ra->rid < rb->rid) return -1;
	if (ra->rid > rb->rid) return 1;
	return 0;
}

static void ri_idx_assign_loci(ri_idx_t *ri)
{
	ri_sig_locus_rec_t *recs;
	uint32_t rid;
	if (!ri || ri->n_seq == 0) return;

	ri_idx_clear_locus_meta(ri);

	ri->sig_assembly_id = (uint32_t *)ri_kcalloc(ri->km, ri->n_seq, sizeof(uint32_t));
	ri->sig_contig_header = (char **)ri_kcalloc(ri->km, ri->n_seq, sizeof(char *));
	ri->sig_start_bp = (uint32_t *)ri_kcalloc(ri->km, ri->n_seq, sizeof(uint32_t));
	ri->sig_end_bp = (uint32_t *)ri_kcalloc(ri->km, ri->n_seq, sizeof(uint32_t));
	ri->sig_locus_index = (uint32_t *)ri_kcalloc(ri->km, ri->n_seq, sizeof(uint32_t));
	ri->sig_locus_start_bp = (uint32_t *)ri_kcalloc(ri->km, ri->n_seq, sizeof(uint32_t));
	ri->sig_locus_end_bp = (uint32_t *)ri_kcalloc(ri->km, ri->n_seq, sizeof(uint32_t));
	if (!ri->sig_assembly_id || !ri->sig_contig_header || !ri->sig_start_bp || !ri->sig_end_bp ||
		!ri->sig_locus_index || !ri->sig_locus_start_bp || !ri->sig_locus_end_bp) {
		ri_idx_clear_locus_meta(ri);
		return;
	}

	for (rid = 0; rid < ri->n_seq; ++rid) {
		ri->sig_assembly_id[rid] = UINT32_MAX;
		ri->sig_locus_index[rid] = UINT32_MAX;
	}

	recs = (ri_sig_locus_rec_t *)calloc(ri->n_seq, sizeof(ri_sig_locus_rec_t));
	if (!recs) {
		ri_idx_clear_locus_meta(ri);
		return;
	}

	for (rid = 0; rid < ri->n_seq; ++rid) {
		const char *name = ri_sig_name(ri, rid);
		const char *taxon = (ri->sig_source_taxon && ri->sig_source_taxon[rid] && ri->sig_source_taxon[rid][0]) ? ri->sig_source_taxon[rid] : "*";
		char *contig_header = NULL;
		uint32_t assembly_id = UINT32_MAX, start_bp = 0, end_bp = 0;
		int parsed = ri_parse_sig_interval(name, &assembly_id, &contig_header, &start_bp, &end_bp);

		recs[rid].rid = rid;
		recs[rid].taxon = taxon;
		recs[rid].parsed = parsed;
		recs[rid].assembly_id = assembly_id;
		recs[rid].contig_header = contig_header;
		recs[rid].start_bp = start_bp;
		recs[rid].end_bp = end_bp;

		ri->sig_assembly_id[rid] = assembly_id;
		ri->sig_contig_header[rid] = contig_header;
		ri->sig_start_bp[rid] = start_bp;
		ri->sig_end_bp[rid] = end_bp;
	}

	qsort(recs, ri->n_seq, sizeof(ri_sig_locus_rec_t), ri_sig_locus_rec_cmp);

	{
		size_t i = 0;
		while (i < ri->n_seq) {
			const char *taxon = (recs[i].taxon && recs[i].taxon[0]) ? recs[i].taxon : "*";
			uint32_t next_locus_index = 1;
			while (i < ri->n_seq) {
				const char *this_taxon = (recs[i].taxon && recs[i].taxon[0]) ? recs[i].taxon : "*";
				uint32_t locus_index;
				if (strcmp(taxon, this_taxon) != 0) break;
				locus_index = next_locus_index++;
				if (!recs[i].parsed) {
					rid = recs[i].rid;
					ri->sig_locus_index[rid] = locus_index;
					ri->sig_locus_start_bp[rid] = recs[i].start_bp;
					ri->sig_locus_end_bp[rid] = recs[i].end_bp;
					++i;
					continue;
				}
				{
					size_t j = i + 1;
					uint32_t locus_start = recs[i].start_bp;
					uint32_t locus_end = recs[i].end_bp;
					while (j < ri->n_seq) {
						const char *next_taxon = (recs[j].taxon && recs[j].taxon[0]) ? recs[j].taxon : "*";
						if (!recs[j].parsed || strcmp(taxon, next_taxon) != 0 ||
							recs[j].assembly_id != recs[i].assembly_id ||
							strcmp(recs[j].contig_header ? recs[j].contig_header : "*",
								recs[i].contig_header ? recs[i].contig_header : "*") != 0 ||
							(uint64_t)recs[j].start_bp > (uint64_t)locus_end + RI_SIGMETA_LOCUS_GAP_BP)
							break;
						if (recs[j].end_bp > locus_end) locus_end = recs[j].end_bp;
						++j;
					}
					while (i < j) {
						rid = recs[i].rid;
						ri->sig_locus_index[rid] = locus_index;
						ri->sig_locus_start_bp[rid] = locus_start;
						ri->sig_locus_end_bp[rid] = locus_end;
						++i;
					}
				}
			}
		}
	}

	free(recs);
}

static void ri_idx_dump_sigmeta_tsv(const char *idx_path, const ri_idx_t *ri)
{
	FILE *fp;
	char *meta_path;
	uint32_t i;
	meta_path = ri_sigmeta_path(idx_path);
	if (!meta_path) return;
	fp = fopen(meta_path, "w");
	if (!fp) {
		free(meta_path);
		return;
	}
	fprintf(fp, "rid\tsignature\tlength_bp\tn_emitted_positions\tn_distinct_hashes\tsource_taxon\tassembly_id\tcontig_header\tsig_start_bp\tsig_end_bp\tlocus_index\tlocus_start_bp\tlocus_end_bp\tlocus_span_bp\n");
	for (i = 0; i < ri->n_seq; ++i) {
		const char *name = ri->seq[i].name;
		uint32_t len = ri->seq[i].len;
		long long assembly_id = (ri->sig_assembly_id && ri->sig_assembly_id[i] != UINT32_MAX) ? (long long)ri->sig_assembly_id[i] : -1LL;
		long long locus_index = (ri->sig_locus_index && ri->sig_locus_index[i] != UINT32_MAX) ? (long long)ri->sig_locus_index[i] : -1LL;
		uint32_t locus_start_bp = ri->sig_locus_start_bp ? ri->sig_locus_start_bp[i] : 0;
		uint32_t locus_end_bp = ri->sig_locus_end_bp ? ri->sig_locus_end_bp[i] : 0;
		uint32_t locus_span_bp =
			(ri->sig_assembly_id && ri->sig_assembly_id[i] != UINT32_MAX && locus_end_bp >= locus_start_bp)
				? (locus_end_bp - locus_start_bp + 1)
				: 0;
		fprintf(fp, "%u\t%s\t%u\t%u\t%u\t%s\t%lld\t%s\t%u\t%u\t%lld\t%u\t%u\t%u\n",
			i, name ? name : "*", len,
			ri->sig_emitted_positions ? ri->sig_emitted_positions[i] : 0,
			ri->sig_distinct_hashes ? ri->sig_distinct_hashes[i] : 0,
			(ri->sig_source_taxon && ri->sig_source_taxon[i]) ? ri->sig_source_taxon[i] : "*",
			assembly_id,
			(ri->sig_contig_header && ri->sig_contig_header[i]) ? ri->sig_contig_header[i] : "*",
			ri->sig_start_bp ? ri->sig_start_bp[i] : 0,
			ri->sig_end_bp ? ri->sig_end_bp[i] : 0,
			locus_index,
			locus_start_bp,
			locus_end_bp,
			locus_span_bp);
	}
	fclose(fp);
	free(meta_path);
}

static void ri_idx_load_sigmeta_tsv(const char *idx_path, ri_idx_t *ri)
{
	FILE *fp;
	char *meta_path, line[8192];
	if (!idx_path || !ri || ri->n_seq == 0) return;
	meta_path = ri_sigmeta_path(idx_path);
	if (!meta_path) return;
	fp = fopen(meta_path, "r");
	if (!fp) {
		free(meta_path);
		return;
	}
	if (!ri->sig_emitted_positions)
		ri->sig_emitted_positions = (uint32_t *)ri_kcalloc(ri->km, ri->n_seq, sizeof(uint32_t));
	if (!ri->sig_distinct_hashes)
		ri->sig_distinct_hashes = (uint32_t *)ri_kcalloc(ri->km, ri->n_seq, sizeof(uint32_t));
	if (!ri->sig_source_taxon)
		ri->sig_source_taxon = (char **)ri_kcalloc(ri->km, ri->n_seq, sizeof(char *));
	if (!fgets(line, sizeof(line), fp)) {
		fclose(fp);
		free(meta_path);
		return;
	}
	while (fgets(line, sizeof(line), fp)) {
		char *fields[6] = {0};
		char *p = line;
		char *field = line;
		int nf = 0;
		while (*p && nf < 6) {
			if (*p == '\t' || *p == '\n' || *p == '\r') {
				char delim = *p;
				*p = '\0';
				fields[nf++] = field;
				if (delim != '\t') break;
				field = p + 1;
			}
			++p;
		}
		if (nf >= 5 && fields[0] && fields[3] && fields[4]) {
			uint32_t rid = (uint32_t)strtoul(fields[0], NULL, 10);
			if (rid < ri->n_seq) {
				ri->sig_emitted_positions[rid] = (uint32_t)strtoul(fields[3], NULL, 10);
				ri->sig_distinct_hashes[rid] = (uint32_t)strtoul(fields[4], NULL, 10);
				if (nf >= 6 && fields[5] && fields[5][0]) {
					if (ri->sig_source_taxon[rid]) free(ri->sig_source_taxon[rid]);
					ri->sig_source_taxon[rid] = ri_strdup0(fields[5]);
				}
			}
		}
	}
	fclose(fp);
	free(meta_path);
}

void ri_idx_stat(const ri_idx_t *ri)
{
	fprintf(stderr, "[M::%s] pore kmer size: %d; concatanated events: %d; quantization bits: %d; w: %d; n: %d; #seq: %d\n", __func__, ri->k, ri->e, ri->q, ri->w, ri->n, ri->n_seq);
}

ri_idx_t* ri_idx_init(int diff, int b, int w, int e, int n, int q, int k, float fine_min, float fine_max, float fine_range, int flag, int strobe_n, int strobe_w_min, int strobe_w_max){
	ri_idx_t* ri;
	ri = (ri_idx_t*)calloc(1, sizeof(ri_idx_t));
  	ri->b = b, ri->w = w; ri->e = e; ri->n = n; ri->q = q; ri->k = k, ri->flag = flag;
	ri->diff = diff;
	ri->fine_min = fine_min, ri->fine_max = fine_max, ri->fine_range = fine_range;
	ri->strobe_n = strobe_n;
	ri->strobe_w_min = strobe_w_min;
	ri->strobe_w_max = strobe_w_max;
	ri->seq = NULL; ri->pore = NULL; ri->h = NULL;
  	ri->B = (ri_idx_bucket_t*)calloc(1<<ri->b, sizeof(ri_idx_bucket_t));
  	ri->km = ri_km_init();

	return ri;
}

void ri_idx_destroy(ri_idx_t *ri){
	uint32_t i;
	if (!ri) return;
	if (ri->h) kh_destroy(str, (khash_t(str)*)ri->h);
	if (ri->B) {
		for (i = 0; i < 1U<<ri->b; ++i) {
			if(ri->B[i].p)free(ri->B[i].p);
			free(ri->B[i].a.a);
			kh_destroy(idx, (idxhash_t*)ri->B[i].h);
		}
	}

	if (ri->sig_source_taxon) {
		for (i = 0; i < ri->n_seq; ++i) {
			if (ri->sig_source_taxon[i]) free(ri->sig_source_taxon[i]);
		}
		ri_kfree(ri->km, ri->sig_source_taxon);
	}
	ri_idx_clear_locus_meta(ri);
	if(ri->sig_emitted_positions) ri_kfree(ri->km, ri->sig_emitted_positions);
	if(ri->sig_distinct_hashes) ri_kfree(ri->km, ri->sig_distinct_hashes);
	if(ri->km) ri_km_destroy(ri->km);
	if (ri->B) free(ri->B);
	free(ri);
}

void ri_idx_add(ri_idx_t *ri, int n, const mm128_t *a){
	int i, mask = (1<<ri->b) - 1;
	for (i = 0; i < n; ++i) {
		mm128_v *p = &ri->B[a[i].x>>RI_HASH_SHIFT&mask].a;
		rh_kv_push(mm128_t, 0, *p, a[i]);
	}
}

static void *worker_pipeline(void *shared, int step, void *in)
{
	int i;
    pipeline_t *p = (pipeline_t*)shared;
    if (step == 0) { // step 0: read sequences
        step_t *s;
		if (p->sum_len > p->batch_size) return 0;
        s = (step_t*)calloc(1, sizeof(step_t));
		s->seq = mm_bseq_read(p->fp, p->mini_batch_size, 0, &s->n_seq); // read a mini-batch
		if (s->seq) {
			uint32_t old_m, m;
			assert((uint64_t)p->ri->n_seq + s->n_seq <= UINT32_MAX); // to prevent integer overflow
			// make room for p->ri->seq
				old_m = p->ri->n_seq, m = p->ri->n_seq + s->n_seq;
				kroundup32(m); kroundup32(old_m);
				if (old_m != m) {
					p->ri->seq = (ri_idx_seq_t*)ri_krealloc(p->ri->km, p->ri->seq, m*sizeof(ri_idx_seq_t));
					p->ri->sig_emitted_positions = (uint32_t*)ri_krealloc(p->ri->km, p->ri->sig_emitted_positions, m*sizeof(uint32_t));
					p->ri->sig_distinct_hashes = (uint32_t*)ri_krealloc(p->ri->km, p->ri->sig_distinct_hashes, m*sizeof(uint32_t));
					memset(p->ri->sig_emitted_positions + old_m, 0, (m - old_m) * sizeof(uint32_t));
					memset(p->ri->sig_distinct_hashes + old_m, 0, (m - old_m) * sizeof(uint32_t));
				}
				for (i = 0; i < s->n_seq; ++i) {
				ri_idx_seq_t *seq = &p->ri->seq[p->ri->n_seq];
				seq->name = (char*)ri_kmalloc(p->ri->km, strlen(s->seq[i].name) + 1);
				strcpy(seq->name, s->seq[i].name);
				seq->len = s->seq[i].l_seq;
				seq->offset = p->sum_len;
				p->sum_len += seq->len;
				s->seq[i].rid = p->ri->n_seq++;
			}
			return s;
		} else free(s);
    } else if (step == 1) { // step 1: the reference signal's seeds
        step_t *s = (step_t*)in;


		for (i = 0; i < s->n_seq; ++i) {
			mm_bseq1_t* t = &s->seq[i];
			if (t->l_seq > 0){
				uint32_t s_len;
				float* s_values = (float*)calloc(t->l_seq, sizeof(float));
				size_t seeds_start = s->a.n;

				ri_reference_signal(t->seq, t->l_seq, p->ri->pore, p->ri->k, 0, &s_len, s_values);
				if (p->ri->flag & RI_I_SEED_MODE_STROBE)
					ri_seeds_strobe(0, s_values, t->rid, 0, s_len, p->ri->diff, p->ri->e, p->ri->strobe_w_min, p->ri->strobe_w_max, p->ri->strobe_n, p->ri->q, p->ri->k, p->ri->fine_min, p->ri->fine_max, p->ri->fine_range, &s->a);
				else
					ri_seeds_contig(0, s_values, t->rid, 0, s_len, p->ri->diff, p->ri->w, p->ri->e, p->ri->n, p->ri->q, p->ri->k, p->ri->fine_min, p->ri->fine_max, p->ri->fine_range, &s->a);

				ri_reference_signal(t->seq, t->l_seq, p->ri->pore, p->ri->k, 1, &s_len, s_values);
				if (p->ri->flag & RI_I_SEED_MODE_STROBE)
					ri_seeds_strobe(0, s_values, t->rid, 1, s_len, p->ri->diff, p->ri->e, p->ri->strobe_w_min, p->ri->strobe_w_max, p->ri->strobe_n, p->ri->q, p->ri->k, p->ri->fine_min, p->ri->fine_max, p->ri->fine_range, &s->a);
				else
					ri_seeds_contig(0, s_values, t->rid, 1, s_len, p->ri->diff, p->ri->w, p->ri->e, p->ri->n, p->ri->q, p->ri->k, p->ri->fine_min, p->ri->fine_max, p->ri->fine_range, &s->a);
				ri_idx_capture_sigmeta(p->ri, t->rid, s->a.a, seeds_start, s->a.n);

				free(s_values);
			}
			free(t->seq); free(t->name);
		}
		free(s->seq); s->seq = 0;
		return s;
    } else if (step == 2) { // seeds to buckets
        step_t *s = (step_t*)in;
		ri_idx_add(p->ri, s->a.n, s->a.a);
		ri_kfree(0, s->a.a); free(s);
	}
    return 0;
}

static void worker_post(void *g, long i, int tid){
	int n, n_keys;
	size_t j, start_a, start_p;
	idxhash_t *h = 0;
	ri_idx_t *ri = (ri_idx_t*)g;
	ri_idx_bucket_t *b = &ri->B[i];
	if (b->a.n == 0) return;

	// sort by seed hash
	radix_sort_128x(b->a.a, b->a.a + b->a.n);

	// count and preallocate
	for (j = 1, n = 1, n_keys = 0, b->n = 0; j <= b->a.n; ++j) {
		//canf: make sure to shift b->a.a[j].x accordingly if we store more data there
		if (j == b->a.n || b->a.a[j].x>>RI_HASH_SHIFT != b->a.a[j-1].x>>RI_HASH_SHIFT) {
			++n_keys;
			if (n > 1) b->n += n;
			n = 1;
		} else ++n;
	}
	h = kh_init(idx);
	kh_resize(idx, h, n_keys);
	b->p = (uint64_t*)calloc(b->n, 8);

	// create the hash table
	for (j = 1, n = 1, start_a = start_p = 0; j <= b->a.n; ++j) {
		if (j == b->a.n || b->a.a[j].x>>RI_HASH_SHIFT != b->a.a[j-1].x>>RI_HASH_SHIFT) {
			khint_t itr;
			int absent;
			mm128_t *p = &b->a.a[j-1];
			itr = kh_put(idx, h, p->x>>RI_HASH_SHIFT>>ri->b<<1, &absent);
			assert(absent && j == start_a + n);
			if (n == 1) {
				kh_key(h, itr) |= 1;
				kh_val(h, itr) = p->y;
			} else {
				int k;
				for (k = 0; k < n; ++k)
					b->p[start_p + k] = b->a.a[start_a + k].y;
				radix_sort_64(&b->p[start_p], &b->p[start_p + n]); // sort by position; needed as in-place radix_sort_128x() is not stable
				kh_val(h, itr) = (uint64_t)start_p<<32 | n;
				start_p += n;
			}
			start_a = j, n = 1;
		} else ++n;
	}
	b->h = h;
	assert(b->n == (int32_t)start_p);

	// deallocate and clear b->a
	ri_kfree(0, b->a.a);
	b->a.n = b->a.m = 0, b->a.a = 0;
}

static void ri_idx_post(ri_idx_t *ri, int n_threads){
	kt_for(n_threads, worker_post, ri, 1<<ri->b);
}

#ifdef ONS_PROF
/* The probe sequence of kh_get(idx, ...), also returning its probe steps. */
static inline khint_t ons_idx_get_steps(const idxhash_t *h, uint64_t key, uint32_t *steps)
{
	*steps = 0;
	if (h->n_buckets) {
		khint_t k, i, last, mask, step = 0;
		mask = h->n_buckets - 1;
		k = idx_hash(key); i = k & mask;
		last = i;
		while (!__ac_isempty(h->flags, i) && (__ac_isdel(h->flags, i) || !idx_eq(h->keys[i], key))) {
			i = (i + (++step)) & mask;
			if (i == last) { *steps = step; return h->n_buckets; }
		}
		*steps = step;
		return __ac_iseither(h->flags, i)? h->n_buckets : i;
	} else return 0;
}
#endif

const uint64_t *ri_idx_get(const ri_idx_t *ri, uint64_t hashval, int *n){

	int mask = (1<<ri->b) - 1;
	khint_t k;
	ri_idx_bucket_t *b = &ri->B[hashval&mask];
	idxhash_t *h = (idxhash_t*)b->h;
	*n = 0;
	#ifdef ONS_PROF
	ons_prof_tls_t *pt = ons_prof_cur ? ons_prof_tls() : NULL;
	if (pt) pt->bucket[hashval&mask]++;
	#endif
	if (h == 0) return 0;
	#ifdef ONS_PROF
	if (pt) {
		uint32_t steps;
		k = ons_idx_get_steps(h, hashval>>ri->b<<1, &steps);
		ons_prof_cur->n_probe += steps;
		pt->probe[steps < ONS_HIST_PROBE ? steps : ONS_HIST_PROBE - 1]++;
		if (ons_prof_keycnt && k != kh_end(h))
			__atomic_fetch_add(&ons_prof_keycnt[ons_prof_slot_off[hashval&mask] + k], 1, __ATOMIC_RELAXED);
	} else
	#endif
	k = kh_get(idx, h, hashval>>ri->b<<1);
	if (k == kh_end(h)) return 0;
	if (kh_key(h, k)&1) { // special casing when there is only one k-mer
		*n = 1;
		return &kh_val(h, k);
	} else {
		*n = (uint32_t)kh_val(h, k);
		return &b->p[kh_val(h, k)>>32];
	}
}

#ifdef ONS_PROF
void ri_idx_prof_dump(const ri_idx_t *ri, FILE *fp)
{
	typedef unsigned long long ull;
	uint64_t n_bkt = 1ULL << ri->b, slots = 0, keys = 0, single = 0, multi = 0, postings = 0, residual_a = 0;
	uint64_t max_post = 0, empty_bkt = 0, i;
	static uint64_t post_hist[ONS_HIST_POST];
	memset(post_hist, 0, sizeof(post_hist));
	for (i = 0; i < n_bkt; ++i) {
		const ri_idx_bucket_t *b = &ri->B[i];
		const idxhash_t *h = (const idxhash_t*)b->h;
		khint_t k;
		residual_a += b->a.m;
		postings += (uint64_t)b->n;
		if (!h) { ++empty_bkt; continue; }
		slots += kh_n_buckets(h);
		keys += kh_size(h);
		for (k = 0; k < kh_end(h); ++k) {
			uint64_t np;
			if (!kh_exist(h, k)) continue;
			if (kh_key(h, k) & 1) { ++single; np = 1; }
			else { ++multi; np = (uint32_t)kh_val(h, k); }
			if (np > max_post) max_post = np;
			post_hist[np < ONS_HIST_POST ? np : ONS_HIST_POST - 1]++;
		}
	}
	fprintf(fp, "key\tvalue\n");
	fprintf(fp, "bucket_bits\t%d\nbuckets\t%llu\nempty_buckets\t%llu\n", ri->b, (ull)n_bkt, (ull)empty_bkt);
	fprintf(fp, "hash_slots\t%llu\nkeys\t%llu\nload_factor\t%.4f\n", (ull)slots, (ull)keys,
		slots ? (double)keys / (double)slots : 0.0);
	fprintf(fp, "keys_single_posting\t%llu\nkeys_multi_posting\t%llu\n", (ull)single, (ull)multi);
	fprintf(fp, "postings_array_entries\t%llu\nmax_postings_per_key\t%llu\n", (ull)postings, (ull)max_post);
	fprintf(fp, "total_postings\t%llu\n", (ull)(postings + single));
	/* bytes of the structures a lookup touches */
	fprintf(fp, "bytes_bucket_array\t%llu\n", (ull)(n_bkt * sizeof(ri_idx_bucket_t)));
	fprintf(fp, "bytes_hash_keys\t%llu\nbytes_hash_vals\t%llu\nbytes_hash_flags\t%llu\n",
		(ull)(slots * 8), (ull)(slots * 8), (ull)(((slots + 15) / 16) * 4));
	fprintf(fp, "bytes_postings_array\t%llu\n", (ull)(postings * 8));
	fprintf(fp, "bytes_unfreed_build_arrays\t%llu\n", (ull)(residual_a * sizeof(mm128_t)));
	fprintf(fp, "signatures\t%u\n", ri->n_seq);
	fprintf(fp, "bits_key_used\t%d\nbits_posting_rid\t32\nbits_posting_pos\t31\nbits_posting_strand\t1\n", 32 - ri->b);
	for (i = 0; i < ONS_HIST_POST; ++i)
		if (post_hist[i]) fprintf(fp, "postings_per_key_hist_%llu\t%llu\n", (ull)i, (ull)post_hist[i]);
}

uint64_t *ri_idx_prof_slot_offsets(const ri_idx_t *ri, uint64_t *n_slots)
{
	uint64_t n_bkt = 1ULL << ri->b, i, off = 0;
	uint64_t *o = (uint64_t*)calloc(n_bkt + 1, sizeof(uint64_t));
	for (i = 0; i < n_bkt; ++i) {
		const idxhash_t *h = (const idxhash_t*)ri->B[i].h;
		o[i] = off;
		if (h) off += kh_n_buckets(h);
	}
	o[n_bkt] = off;
	*n_slots = off;
	return o;
}
#endif

void ri_idx_dump(FILE* idx_file, const ri_idx_t* ri){

	uint32_t pars[7], i;

	pars[0] = ri->w, pars[1] = ri->e, pars[2] = ri->n, pars[3] = ri->q, pars[4] = ri->k, pars[5] = ri->n_seq, pars[6] = ri->flag;
	
	fwrite(RI_IDX_MAGIC, 1, RI_IDX_MAGIC_BYTE, idx_file);
	fwrite(pars, sizeof(uint32_t), 7, idx_file);
	fwrite(&ri->diff, sizeof(int), 1, idx_file);
	fwrite(&ri->fine_min, sizeof(float), 1, idx_file);
	fwrite(&ri->fine_max, sizeof(float), 1, idx_file);
	fwrite(&ri->fine_range, sizeof(float), 1, idx_file);
	/* strobe-events: written only with RI_I_SEED_MODE_STROBE, so contiguous-seed indexes keep the older header */
	if (ri->flag & RI_I_SEED_MODE_STROBE) {
		uint8_t s_n = (uint8_t)ri->strobe_n;
		uint8_t s_wmin = (uint8_t)ri->strobe_w_min;
		uint8_t s_wmax = (uint8_t)ri->strobe_w_max;
		uint8_t s_reserved = 0;
		fwrite(&s_n,        1, 1, idx_file);
		fwrite(&s_wmin,     1, 1, idx_file);
		fwrite(&s_wmax,     1, 1, idx_file);
		fwrite(&s_reserved, 1, 1, idx_file);
	}
	fwrite(ri->pore, sizeof(ri_pore_t), 1, idx_file);
	fwrite(ri->pore->pore_vals, sizeof(float), ri->pore->n_pore_vals, idx_file);
	ri_kfree(ri->km, ri->pore->pore_vals); 
	ri->pore->pore_vals = 0;
	fwrite(ri->pore->pore_inds, sizeof(ri_porei_t), ri->pore->n_pore_vals, idx_file);
	ri_kfree(ri->km, ri->pore->pore_inds); 
	ri->pore->pore_inds = 0;
	ri_kfree(ri->km, ri->pore);

	for (i = 0; i < ri->n_seq; ++i) {
		if (ri->seq[i].name) {
			uint8_t l = strlen(ri->seq[i].name);
			fwrite(&l, 1, 1, idx_file);
			fwrite(ri->seq[i].name, 1, l, idx_file);
		} else {
			uint8_t l = 0;
			fwrite(&l, 1, 1, idx_file);
		}
		fwrite(&ri->seq[i].len, 4, 1, idx_file);
	}

	for (i = 0; i < 1U<<ri->b; ++i) {
		ri_idx_bucket_t *b = &ri->B[i];
		khint_t k;
		idxhash_t *h = (idxhash_t*)b->h;
		uint32_t size = h? h->size : 0;
		fwrite(&b->n, 4, 1, idx_file);
		fwrite(b->p, 8, b->n, idx_file);
		fwrite(&size, 4, 1, idx_file);
		if (size == 0) continue;
		for (k = 0; k < kh_end(h); ++k) {
			uint64_t x[2];
			if (!kh_exist(h, k)) continue;
			x[0] = kh_key(h, k), x[1] = kh_val(h, k);
			fwrite(x, 8, 2, idx_file);
		}
	}

	fflush(idx_file);
}

ri_idx_t* ri_idx_load(FILE* idx_file){

	ri_idx_t* ri;
	uint64_t sum_len = 0;

	char magic[RI_IDX_MAGIC_BYTE];
  	uint32_t i;

	if (fread(magic, 1, RI_IDX_MAGIC_BYTE, idx_file) != RI_IDX_MAGIC_BYTE) return 0;
	if (strncmp(magic, RI_IDX_MAGIC, RI_IDX_MAGIC_BYTE) != 0) return 0;
	int pars[7];
	fread(&pars[0], sizeof(int), 7, idx_file);
	if (pars[6] & (RI_I_STORE_SIG | RI_I_SIG_TARGET)) {
		fprintf(stderr, "[ERROR] the index was built with --store-sig or --sig-target, which onsight no longer supports; rebuild it\n");
		return 0;
	}

	int diff;
	float fine_min, fine_max, fine_range;
	fread(&diff, sizeof(int), 1, idx_file);
	fread(&fine_min, sizeof(float), 1, idx_file);
	fread(&fine_max, sizeof(float), 1, idx_file);
	fread(&fine_range, sizeof(float), 1, idx_file);

	int s_n = 0, s_wmin = 0, s_wmax = 0;
	if (pars[6] & RI_I_SEED_MODE_STROBE) {
		uint8_t b_n = 0, b_wmin = 0, b_wmax = 0, b_reserved = 0;
		fread(&b_n,        1, 1, idx_file);
		fread(&b_wmin,     1, 1, idx_file);
		fread(&b_wmax,     1, 1, idx_file);
		fread(&b_reserved, 1, 1, idx_file);
		s_n = b_n; s_wmin = b_wmin; s_wmax = b_wmax;
	}

	ri = ri_idx_init(diff, 14, pars[0], pars[1], pars[2], pars[3], pars[4], fine_min, fine_max, fine_range, pars[6], s_n, s_wmin, s_wmax);
	ri->n_seq = pars[5];
	ri->seq = (ri_idx_seq_t*)ri_kcalloc(ri->km, ri->n_seq, sizeof(ri_idx_seq_t));

	ri->pore = (ri_pore_t*)ri_kmalloc(ri->km, sizeof(ri_pore_t));
	fread(ri->pore, sizeof(ri_pore_t), 1, idx_file);
	ri->pore->pore_vals = (float*)ri_kmalloc(ri->km, ri->pore->n_pore_vals * sizeof(float));
	fread(ri->pore->pore_vals, sizeof(float), ri->pore->n_pore_vals, idx_file);
	ri->pore->pore_inds = (ri_porei_t*)ri_kmalloc(ri->km, ri->pore->n_pore_vals * sizeof(ri_porei_t));
	fread(ri->pore->pore_inds, sizeof(ri_porei_t), ri->pore->n_pore_vals, idx_file);

	for (i = 0; i < ri->n_seq; ++i) {
		uint8_t l;
		ri_idx_seq_t *s = &ri->seq[i];
		fread(&l, 1, 1, idx_file);
		if (l) {
			s->name = (char*)ri_kmalloc(ri->km, l + 1);
			fread(s->name, 1, l, idx_file);
			s->name[l] = 0;
		}
		fread(&s->len, 4, 1, idx_file);
		s->offset = sum_len;
		sum_len += s->len;
	}
	
	for (i = 0; i < 1U<<ri->b; ++i) {
		ri_idx_bucket_t *b = &ri->B[i];
		uint32_t j, size;
		khint_t k;
		idxhash_t *h;
		fread(&b->n, 4, 1, idx_file);
		b->p = (uint64_t*)malloc(b->n * 8);
		fread(b->p, 8, b->n, idx_file);
		fread(&size, 4, 1, idx_file);
		if (size == 0) continue;
		b->h = h = kh_init(idx);
		kh_resize(idx, h, size);
		for (j = 0; j < size; ++j) {
			uint64_t x[2];
			int absent;
			fread(x, 8, 2, idx_file);
			k = kh_put(idx, h, x[0], &absent);
			assert(absent);
			kh_val(h, k) = x[1];
		}
	}

	return ri;
}

ri_idx_reader_t* ri_idx_reader_open(const char *fn, const ri_idxopt_t *ipt, const char *fn_out)
{
	int64_t is_idx;
	ri_idx_reader_t *r;
	is_idx = ri_idx_is_idx(fn);
	if (is_idx < 0) return 0; // failed to open the index
	r = (ri_idx_reader_t*)calloc(1, sizeof(ri_idx_reader_t));
	r->is_idx = is_idx;
	r->fn = fn ? strdup(fn) : NULL;
	r->fn_out_name = fn_out ? strdup(fn_out) : NULL;
	if (ipt) r->opt = *ipt;
	else ri_idxopt_init(&r->opt);
	if (r->is_idx) {
		r->fp.idx = fopen(fn, "rb");
		r->idx_size = is_idx;
	} else r->fp.seq = mm_bseq_open(fn);
	if (fn_out) r->fp_out = fopen(fn_out, "wb");
	return r;
}

void ri_idx_reader_close(ri_idx_reader_t* r){

	if (r->is_idx) fclose(r->fp.idx);
	else if(r->fp.seq) mm_bseq_close(r->fp.seq);
	if (r->fp_out) fclose(r->fp_out);
	if (r->fn) free(r->fn);
	if (r->fn_out_name) free(r->fn_out_name);
	free(r);
}

ri_idx_t* ri_idx_gen(mm_bseq_file_t* fp, ri_pore_t* pore, int diff, int b, int w, int e, int n, int q, int k, float fine_min, float fine_max, float fine_range, int flag, int strobe_n, int strobe_w_min, int strobe_w_max, int mini_batch_size, int n_threads, int io_n_threads, uint64_t batch_size)
{
	pipeline_t pl;
	if (fp == 0 || mm_bseq_eof(fp) || !pore || !pore->pore_vals) return 0;
	memset(&pl, 0, sizeof(pipeline_t));
	pl.mini_batch_size = (uint64_t)mini_batch_size < batch_size? mini_batch_size : batch_size;
	pl.batch_size = batch_size;
	pl.fp = fp;
	pl.ri = ri_idx_init(diff, b, w, e, n, q, k, fine_min, fine_max, fine_range, flag, strobe_n, strobe_w_min, strobe_w_max);
	
	pl.ri->pore = (ri_pore_t*)ri_kmalloc(pl.ri->km, sizeof(ri_pore_t));
	memcpy(pl.ri->pore, pore, sizeof(ri_pore_t));
	pl.ri->pore->pore_vals = (float*)ri_kmalloc(pl.ri->km, pore->n_pore_vals * sizeof(float));
	memcpy(pl.ri->pore->pore_vals, pore->pore_vals, pore->n_pore_vals * sizeof(float));
	pl.ri->pore->pore_inds = (ri_porei_t*)ri_kmalloc(pl.ri->km, pore->n_pore_vals * sizeof(ri_porei_t));
	memcpy(pl.ri->pore->pore_inds, pore->pore_inds, pore->n_pore_vals * sizeof(ri_porei_t));

	kt_pipeline(n_threads < 3? n_threads : 3, worker_pipeline, &pl, 3);
	ri_idx_post(pl.ri, n_threads);

	return pl.ri;
}

ri_idx_t* ri_idx_reader_read(ri_idx_reader_t* r, ri_pore_t* pore, int n_threads, int io_n_threads){

	ri_idx_t *ri = 0;
	if (r->is_idx) {
		ri = ri_idx_load(r->fp.idx);
		if (ri) ri_idx_load_sigmeta_tsv(r->fn, ri);
	} else{
		ri = ri_idx_gen(r->fp.seq, pore, r->opt.diff, r->opt.b, r->opt.w, r->opt.e, r->opt.n, r->opt.q, r->opt.k, r->opt.fine_min, r->opt.fine_max, r->opt.fine_range, r->opt.flag, r->opt.strobe_n, r->opt.strobe_w_min, r->opt.strobe_w_max, r->opt.mini_batch_size, n_threads, io_n_threads, r->opt.batch_size);
	}

	if (ri) {
		if (!r->is_idx && r->opt.metadata_path)
			ri_idx_apply_metadata_csv(r->opt.metadata_path, ri);
		ri_idx_assign_loci(ri);
		if (r->fp_out) {
			ri_idx_dump(r->fp_out, ri);
			if (r->fn_out_name) ri_idx_dump_sigmeta_tsv(r->fn_out_name, ri);
		}
		ri->index = r->n_parts++;
	}

	return ri;
}

int64_t ri_idx_is_idx(const char* fn){

	int fd, is_idx = 0;
	int64_t ret, off_end;
	char magic[RI_IDX_MAGIC_BYTE];

	if (strcmp(fn, "-") == 0) return 0; // read from pipe; not an index
	fd = open(fn, O_RDONLY);
	if (fd < 0) return -1; // error
#ifdef WIN32
	if ((off_end = _lseeki64(fd, 0, SEEK_END)) >= RI_IDX_MAGIC_BYTE) {
		_lseeki64(fd, 0, SEEK_SET);
#else
	if ((off_end = lseek(fd, 0, SEEK_END)) >= RI_IDX_MAGIC_BYTE) {
		lseek(fd, 0, SEEK_SET);
#endif // WIN32
		ret = read(fd, magic, RI_IDX_MAGIC_BYTE);
		if (ret == RI_IDX_MAGIC_BYTE && strncmp(magic, RI_IDX_MAGIC, RI_IDX_MAGIC_BYTE) == 0)
			is_idx = 1;
	}
	close(fd);
	return is_idx? off_end : 0;
}

int32_t ri_idx_cal_max_occ(const ri_idx_t *ri, float f)
{
	int i;
	size_t n = 0;
	uint32_t thres;
	khint_t *a, k;
	if (f <= 0.) return INT32_MAX;
	for (i = 0; i < 1<<ri->b; ++i)
		if (ri->B[i].h) n += kh_size((idxhash_t*)ri->B[i].h);
	a = (uint32_t*)malloc(n * 4);
	for (i = n = 0; i < 1<<ri->b; ++i) {
		idxhash_t *h = (idxhash_t*)ri->B[i].h;
		if (h == 0) continue;
		for (k = 0; k < kh_end(h); ++k) {
			if (!kh_exist(h, k)) continue;
			a[n++] = kh_key(h, k)&1? 1 : (uint32_t)kh_val(h, k);
		}
	}
	thres = ks_ksmall_uint32_t(n, a, (uint32_t)((1. - f) * n)) + 1;
	free(a);
	return thres;
}

void ri_detopt_update(ri_detopt_t *opt, const ri_idx_t *ri)
{
	if (opt->mid_occ <= 0) {
		opt->mid_occ = ri_idx_cal_max_occ(ri, opt->mid_occ_frac);
		if (opt->mid_occ < opt->min_mid_occ)
			opt->mid_occ = opt->min_mid_occ;
		if (opt->max_mid_occ > opt->min_mid_occ && opt->mid_occ > opt->max_mid_occ)
			opt->mid_occ = opt->max_mid_occ;

		fprintf(stderr, "[M::%s::%.6f] mid_occ = %d, min_mid_occ = %d, max_mid_occ = %d\n", __func__,
				opt->mid_occ_frac, opt->mid_occ, opt->min_mid_occ, opt->max_mid_occ);
	}
}

const char *ri_sig_taxon(const ri_idx_t *ri, uint32_t rid)
{
	if (!ri || rid >= ri->n_seq || !ri->sig_source_taxon || !ri->sig_source_taxon[rid] || ri->sig_source_taxon[rid][0] == '\0')
		return "*";
	return ri->sig_source_taxon[rid];
}

uint32_t ri_sig_locus(const ri_idx_t *ri, uint32_t rid)
{
	if (!ri || rid >= ri->n_seq) return UINT32_MAX;
	if (ri->sig_locus_index && ri->sig_locus_index[rid] != UINT32_MAX)
		return ri->sig_locus_index[rid];
	return rid;
}
