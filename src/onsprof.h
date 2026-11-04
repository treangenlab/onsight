/* onsprof.h: hardware-characterization profiling of the detection path.
 *
 * Compiled only with -DONS_PROF (`make prof` at the root; `make attr` adds -DONS_ATTRIB). The default build never
 * includes this code. Nothing here changes the results: it only reads state, counts and times.
 *
 * Output files use the prefix in env ONS_PROF_OUT (default "ons_prof"):
 *   <prefix>.reads.tsv    one row per read: work counts and ns per stage
 *   <prefix>.batches.tsv  one row per minibatch: pipeline step times, RSS
 *   <prefix>.index.tsv    index layout and size (key/value)
 *   <prefix>.summary.tsv  run totals, wait/work per pipeline thread (key/value)
 *   <prefix>.hist.tsv     histograms (segment length, positions per found seed, probe steps, hits per signature)
 *   <prefix>.buckets.tsv  lookups per index bucket
 *   <prefix>.keys.tsv     lookups absorbed by the hottest keys (only with env ONS_PROF_KEYS=1, a counter per key)
 */
#ifndef ONSPROF_H
#define ONSPROF_H

#include <stdint.h>

/* Attribution build (make attr): hot helpers stay separate symbols, so perf can attribute cycles and cache misses
 * to each stage. */
#if defined(ONS_PROF) && defined(ONS_ATTRIB)
#define ONS_NOINLINE __attribute__((noinline))
#else
#define ONS_NOINLINE
#endif

#ifdef ONS_PROF
#include <stdio.h>
#include <time.h>

#ifdef __cplusplus
extern "C" {
#endif

#define ONS_HIST_SEG   1024  /* segment length (samples), last bin = overflow */
#define ONS_HIST_POST  4096  /* positions per found seed, last bin = overflow */
#define ONS_HIST_PROBE 64    /* hash-table probe steps per lookup */
#define ONS_HIST_HITS  1025  /* hits kept per (read, signature) */
#define ONS_KTP_MAXW   8
#define ONS_KTP_MAXS   4

typedef struct ons_prof_read_s {
	uint64_t read_id;
	uint32_t batch;
	int      tid;

	/* step 0: reading the signal */
	uint32_t l_sig;          /* samples kept after the 30-200 pA filter */
	uint64_t t_read;

	/* step 1: event detection */
	uint32_t n_norm;         /* samples kept after the |z| < 3 filter */
	uint32_t n_peaks;
	uint32_t n_seg_ok, n_seg_rej;  /* segments made into events / rejected by length */
	uint32_t n_events;
	uint32_t below_min_events;
	uint64_t t_norm, t_prefix, t_tstat, t_peaks, t_events;

	/* step 1: seeds (quantization, strobe-events), index lookup, seed hits per signature */
	uint32_t n_seeds;
	uint32_t n_found, n_single;
	uint64_t n_probe;        /* probe steps beyond the first slot */
	uint64_t post_found;     /* index positions of all found seeds */
	uint32_t n_flt_occ;      /* seeds over the repeat cutoff (mid_occ) */
	uint64_t post_flt_occ;
	uint32_t n_skip_adm;     /* seeds over --sig-max-hash-occ */
	uint64_t post_skip_adm;
	uint64_t post_visited;   /* positions added as seed hits */
	uint64_t n_pos_new;      /* new (signature, position) */
	uint64_t n_cap_scan, n_cap_evict;  /* hit cap full: lightest-hit scans, replacements */
	uint32_t n_sigs;         /* signatures with hits */
	uint64_t n_hits_kept;
	uint64_t hits_bytes;     /* the read's signature table and its block of hits */
	uint64_t km_cap;         /* the thread's memory pool after the read (bytes) */
	uint64_t t_seeds, t_lookup, t_accum, t_step1, t_step1_cpu;

	/* step 2, per read on all threads: --sig-summary, then the stream scores */
	uint32_t n_taxa;
	int      call;           /* read p-value <= the top cutoff */
	int      checked;        /* near-neighbor check: -3 skipped, -2 not wanted, -1 no band, 0 dropped, 1 kept */
	uint64_t t_summary, t_bands, t_pvalue, t_slopes, t_neighbor, t_prepare;

	/* step 2, in read order: the rows, the profile, the stream commit (slope ring, read record, read/target calls) */
	uint64_t t_commit;
} ons_prof_read_t;

typedef struct ons_prof_tls_s {
	uint64_t seg_len[ONS_HIST_SEG];
	uint64_t post[ONS_HIST_POST];
	uint64_t probe[ONS_HIST_PROBE];
	uint64_t hits[ONS_HIST_HITS];
	uint64_t *bucket;        /* lookups per index bucket */
} ons_prof_tls_t;

extern int ons_prof_on;
extern __thread ons_prof_read_t *ons_prof_cur;
extern uint32_t *ons_prof_keycnt;          /* per-key lookup counts, or NULL */
extern uint64_t *ons_prof_slot_off;        /* first hash slot of each bucket, for ons_prof_keycnt */
extern uint64_t ons_ktp_wait[ONS_KTP_MAXW][ONS_KTP_MAXS];
extern uint64_t ons_ktp_work[ONS_KTP_MAXW][ONS_KTP_MAXS];

static inline uint64_t ons_now(void)
{
	struct timespec t;
	clock_gettime(CLOCK_MONOTONIC, &t);
	return (uint64_t)t.tv_sec * 1000000000ULL + (uint64_t)t.tv_nsec;
}

static inline uint64_t ons_cpu_now(void)
{
	struct timespec t;
	clock_gettime(CLOCK_THREAD_CPUTIME_ID, &t);
	return (uint64_t)t.tv_sec * 1000000000ULL + (uint64_t)t.tv_nsec;
}

/* Histograms: one slot per step-1 thread (tid). kt_for() starts new threads every minibatch, and one minibatch at a
 * time is in step 1, so a slot is never shared. */
ons_prof_tls_t *ons_prof_tls(void);
void ons_prof_tls_set(int slot);
ons_prof_read_t *ons_prof_read_new(void);
void ons_prof_read_write(const ons_prof_read_t *r, const char *name);
void ons_prof_batch_write(uint32_t batch, uint32_t n_reads, uint64_t n_samples, uint64_t t_start, uint64_t step0,
                          uint64_t t1_start, uint64_t step1, const uint64_t *busy, int n_busy,
                          uint64_t t2_start, uint64_t step2_par, uint64_t step2_ser, uint64_t end_batch);

struct ri_idx_s;
struct ri_detopt_s;
void ons_prof_init(const struct ri_idx_s *ri, const struct ri_detopt_s *opt, int n_threads, int io_threads,
                   int pl_threads);
void ons_prof_finish(void);

/* index.c */
void ri_idx_prof_dump(const struct ri_idx_s *ri, FILE *fp);
uint64_t *ri_idx_prof_slot_offsets(const struct ri_idx_s *ri, uint64_t *n_slots);

#ifdef __cplusplus
}
#endif

#endif /* ONS_PROF */
#endif /* ONSPROF_H */
