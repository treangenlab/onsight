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
	// uint32_t is_alt;
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

/**
 * Prints the statistics of the index (pore kmer size, concatanated events, quantization method, w, number of sequences in the reference genome) to stderr.
 *
 * @param ri	index
 * 
 */
void ri_idx_stat(const ri_idx_t *ri);

/**
 * Initialize an index reader
 *
 * @param fn         index or fasta/fastq file name (this function tests the file type)
 * @param ipt        indexing parameters
 * @param fn_out     if not NULL, write built index to this file
 *
 * @return an index reader on success; NULL if fail to open _fn_
 */
ri_idx_reader_t* ri_idx_reader_open(const char *fn, const ri_idxopt_t *ipt, const char *fn_out);

/**
 * Closes, deallocates, and destroys the index reader
 *
 * @param r		index reader
 */
void ri_idx_reader_close(ri_idx_reader_t* r);

/**
 * Checks whether the file is an index
 *
 * @param fn	file name
 *
 * @return the file size if fn is an index file; 0 if fn is not.
 */
int64_t ri_idx_is_idx(const char* fn);

/**
 * Initialize the index with its constant parameters such as window length $w, or number of events to concatanate $n
 *
 * @param b		defines the number of buckets to use in a hash table (2^b buckets). Usually this value is around 14.
 * @param w		window length. if w=0, the minimizer-based seeding is disabled
 * @param e		number of events to concatanate in a single hash value
 * @param n		[Currently disabled] number of items in a seed to generate a hash value for the seed using the BLEND mechanism.
 * @param q		Most significant Q bits to use in the quantization technique. See the RawHash manuscript for details.
 * @param lq    least significant lq bits within the most significant Q bits extracted from the normalized event values.
 * 				If you want to use n as described in the RawHash manuscript, calculate lq as follows: $lq = $Q - 2 - n
 * @param k     k-mer size that a single event represents (usually given in the k-mer model file)
 * 
 * @return		the index
 */
ri_idx_t* ri_idx_init(int diff, int b, int w, int e, int n, int q, int k, float fine_min, float fine_max, float fine_range, int flag, int strobe_n, int strobe_w_min, int strobe_w_max);

/**
 * Reads or constructs the index from file. If the file is not index, it should be a file containing sequences to generate the index for.
 *
 * @param r				opened index reader (see ri_idx_reader_open to open)
 * @param pore_vals		expected event values of each k-mer based on the nanopore model (usually provided as a k-mer model file)
 * @param n_threads		number of threads to use when constructing the index
 * 
 * @return				the index
 */
ri_idx_t *ri_idx_reader_read(ri_idx_reader_t *r, ri_pore_t* pore, int n_threads, int io_n_threads);

/**
 * Adds the values in an array to the index
 *
 * @param ri	index
 * @param n		number of items in the array $a
 * @param a		array of values to add into the index.
 * 				a[i].x>>RI_HASH_SHIFT is used as the hash value (key) in the hash table
 */
void ri_idx_add(ri_idx_t* ri, int n, const mm128_t* a);

/**
 * Writes the index to a file
 *
 * @param idx_file	path to the file to write the index $ri
 * @param ri		index
 * 
 */
void ri_idx_dump(FILE* idx_file, const ri_idx_t* ri);

/**
 * Reads the index from file
 *
 * @param fp	path to the index file
 * 
 * @return		the index
 */
ri_idx_t* ri_idx_load(FILE* fp);

/**
 * Deallocates and destroys the entire index
 *
 * @param ri	index to destroy
 * 
 */
void ri_idx_destroy(ri_idx_t* ri);

/**
 * Queries the hash table to find the values (list of values) of a key (hash value)
 *
 * @param ri		index containing the hash table
 * @param hashval	hash value (key) to query the hash table
 * @param n			number of values that are stored using the same key
 * 
 * @return			pointer to the list of values that share the same key (hash value)
 */
const uint64_t *ri_idx_get(const ri_idx_t *ri, uint64_t hashval, int *n);

int32_t ri_idx_cal_max_occ(const ri_idx_t *ri, float f);

void ri_detopt_update(ri_detopt_t *opt, const ri_idx_t *ri);

/* The signature table: signature rid's source taxon ("*" if none) and its locus (rid itself if none). */
const char *ri_sig_taxon(const ri_idx_t *ri, uint32_t rid);
uint32_t ri_sig_locus(const ri_idx_t *ri, uint32_t rid);

#ifdef __cplusplus
}
#endif

#endif // INDEX_H
