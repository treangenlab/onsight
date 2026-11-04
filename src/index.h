#ifndef INDEX_H
#define INDEX_H

/* index.h: the index of the panel's signatures (index build), and the signature table.
 *
 * Built once per panel: each signature's reference signal (pore_model.h), on both strands, quantized into seeds
 * (seeds.h); the seeds' hashes in a hash table, split into 2^b buckets, with each hash's positions (signature, position,
 * strand). The signature table (--metadata, saved next to the index as .sigmeta.tsv) gives each signature its source
 * taxon and locus, and the emitted positions and distinct hashes that the scores divide by.
 *
 * In: the signatures (FASTA), the pore model, --metadata. Out: the index file, read back at detection. State: the
 * index, read only during detection.
 */

#include "utils.h"
#include "pore_model.h"
#include "options.h"
#include "signal_input.h"

#define RI_IDX_MAGIC   "RI"
#define RI_IDX_MAGIC_BYTE 2

#ifdef __cplusplus
extern "C" {
#endif

typedef struct ri_idx_bucket_s {
	mm128_v a;   // (seed, position) array
	int32_t n;   // size of the _p_ array
	uint64_t *p; // position array for seeds appearing >1 times
	void *h;     // hash table indexing _p_ and seeds appearing once
} ri_idx_bucket_t;

typedef struct ri_idx_seq_s{
	char *name;      // name of the db sequence
	uint64_t offset; // offset in ri_idx_t::S
	uint32_t len;    // length
} ri_idx_seq_t;

typedef struct ri_idx_s{
	int32_t b, w, e, n, q, k, flag;
	int32_t index;
	int diff;
	float fine_min, fine_max, fine_range;
	/* strobe-events (flag & RI_I_SEED_MODE_STROBE): e is the events per strobe, these the rest */
	int32_t strobe_n;
	int32_t strobe_w_min;
	int32_t strobe_w_max;
	struct ri_idx_bucket_s *B; // index (hidden)

	ri_pore_t* pore;

	void *km, *h;

	uint32_t n_seq;
	ri_idx_seq_t *seq;
	uint32_t *sig_emitted_positions;
	uint32_t *sig_distinct_hashes;
	char **sig_source_taxon;
	uint32_t *sig_assembly_id;
	char **sig_contig_header;
	uint32_t *sig_start_bp;
	uint32_t *sig_end_bp;
	uint32_t *sig_locus_index;
	uint32_t *sig_locus_start_bp;
	uint32_t *sig_locus_end_bp;
} ri_idx_t;

// index reader
typedef struct ri_idx_reader_s{
	int is_idx, n_parts;
	int64_t idx_size;
	ri_idxopt_t opt;
	FILE *fp_out;
	union {
		struct mm_bseq_file_s *seq;
		FILE *idx;
	} fp;

	char *fn;
	char *fn_out_name;
} ri_idx_reader_t;

/* The index statistics (k-mer size, events per seed, quantization bits, w, signatures) to stderr. */
void ri_idx_stat(const ri_idx_t *ri);

/* A reader of fn, an index or a FASTA/FASTQ of signatures (the type is tested); with fn_out, an index built from the
 * sequences is also written there. NULL if fn cannot be opened. */
ri_idx_reader_t* ri_idx_reader_open(const char *fn, const ri_idxopt_t *ipt, const char *fn_out);
void ri_idx_reader_close(ri_idx_reader_t* r);

/* The file size if fn is an index file, 0 if not. */
int64_t ri_idx_is_idx(const char* fn);

/* An empty index: 2^b hash buckets, minimizer window w (0: every seed), e events per seed, q quantization bits, k-mer
 * size k. n is unused (RawHash2's BLEND) and kept for the index file. */
ri_idx_t* ri_idx_init(int diff, int b, int w, int e, int n, int q, int k, float fine_min, float fine_max, float fine_range, int flag, int strobe_n, int strobe_w_min, int strobe_w_max);

/* The reader's next index: read from the index file, or built from the sequences with the pore model. */
ri_idx_t *ri_idx_reader_read(ri_idx_reader_t *r, ri_pore_t* pore, int n_threads, int io_n_threads);

/* Adds n seeds; a[i].x >> RI_HASH_SHIFT is the hash (key). */
void ri_idx_add(ri_idx_t* ri, int n, const mm128_t* a);

/* Writes the index to idx_file, or reads one from fp. */
void ri_idx_dump(FILE* idx_file, const ri_idx_t* ri);
ri_idx_t* ri_idx_load(FILE* fp);
void ri_idx_destroy(ri_idx_t* ri);

/* The positions stored under hash hashval, *n of them. */
const uint64_t *ri_idx_get(const ri_idx_t *ri, uint64_t hashval, int *n);

/* The repeat cutoff: the (1 - f) quantile of positions per hash; ri_detopt_update sets opt->mid_occ from it. */
int32_t ri_idx_cal_max_occ(const ri_idx_t *ri, float f);
void ri_detopt_update(ri_detopt_t *opt, const ri_idx_t *ri);

/* The signature table: signature rid's source taxon ("*" if none) and its locus (rid itself if none). */
const char *ri_sig_taxon(const ri_idx_t *ri, uint32_t rid);
uint32_t ri_sig_locus(const ri_idx_t *ri, uint32_t rid);

#ifdef __cplusplus
}
#endif

#endif // INDEX_H
