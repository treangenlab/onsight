/* onsprof.c: hardware-characterization profiling; see onsprof.h for the output files. */
#ifdef ONS_PROF

#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif

#include "onsprof.h"
#include "index.h"
#include "options.h"
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/resource.h>

int ons_prof_on = 0;
__thread ons_prof_read_t *ons_prof_cur = NULL;
uint32_t *ons_prof_keycnt = NULL;
uint64_t *ons_prof_slot_off = NULL;
uint64_t ons_ktp_wait[ONS_KTP_MAXW][ONS_KTP_MAXS];
uint64_t ons_ktp_work[ONS_KTP_MAXW][ONS_KTP_MAXS];

static __thread ons_prof_tls_t *g_tls = NULL;
static ons_prof_tls_t **g_slots = NULL;
static int g_n_slots_tls = 0;
static FILE *g_reads = NULL, *g_batches = NULL;
static char g_prefix[1024];
static uint32_t g_n_buckets = 0;
static uint64_t g_n_slots = 0;
static uint64_t g_t_start = 0;
static uint64_t g_self_bytes = 0;  /* memory the profiler itself holds */
static int g_n_threads = 0, g_io_threads = 0, g_pl_threads = 0;
static int64_t g_batch_size = 0;
static int32_t g_mid_occ = 0;
static uint32_t g_max_hash_occ = 0;

/* run totals, summed in the output step (one thread, read order) */
static uint64_t g_reads_n = 0, g_samples = 0, g_events = 0, g_seeds = 0;

static FILE *prof_open(const char *suffix)
{
	char fn[1100];
	FILE *fp;
	snprintf(fn, sizeof(fn), "%s.%s", g_prefix, suffix);
	fp = fopen(fn, "w");
	if (!fp) fprintf(stderr, "[ONS_PROF] cannot open %s\n", fn);
	return fp;
}

static uint64_t rss_bytes(void)
{
	long pages = 0, rss = 0;
	FILE *fp = fopen("/proc/self/statm", "r");
	if (!fp) return 0;
	if (fscanf(fp, "%ld %ld", &pages, &rss) != 2) rss = 0;
	fclose(fp);
	return (uint64_t)rss * (uint64_t)sysconf(_SC_PAGESIZE);
}

ons_prof_tls_t *ons_prof_tls(void)
{
	return g_tls;
}

void ons_prof_tls_set(int slot)
{
	g_tls = (ons_prof_on && slot >= 0 && slot < g_n_slots_tls) ? g_slots[slot] : NULL;
}

ons_prof_read_t *ons_prof_read_new(void)
{
	if (!ons_prof_on) return NULL;
	return (ons_prof_read_t*)calloc(1, sizeof(ons_prof_read_t));
}

void ons_prof_init(const ri_idx_t *ri, const ri_detopt_t *opt, int n_threads, int io_threads, int pl_threads)
{
	const char *p = getenv("ONS_PROF_OUT");
	const char *keys = getenv("ONS_PROF_KEYS");
	FILE *fp;

	snprintf(g_prefix, sizeof(g_prefix), "%s", (p && *p) ? p : "ons_prof");
	ons_prof_on = 1;
	g_t_start = ons_now();
	g_n_buckets = 1U << ri->b;
	g_n_threads = n_threads, g_io_threads = io_threads, g_pl_threads = pl_threads;
	g_batch_size = opt->mini_batch_size;
	g_mid_occ = opt->mid_occ;
	g_max_hash_occ = opt->sig_max_hash_occ;
	memset(ons_ktp_wait, 0, sizeof(ons_ktp_wait));
	memset(ons_ktp_work, 0, sizeof(ons_ktp_work));

	g_n_slots_tls = n_threads;
	g_slots = (ons_prof_tls_t**)calloc(g_n_slots_tls, sizeof(ons_prof_tls_t*));
	for (int i = 0; i < g_n_slots_tls; ++i) {
		g_slots[i] = (ons_prof_tls_t*)calloc(1, sizeof(ons_prof_tls_t));
		g_slots[i]->bucket = (uint64_t*)calloc(g_n_buckets, sizeof(uint64_t));
		g_self_bytes += sizeof(ons_prof_tls_t) + (uint64_t)g_n_buckets * sizeof(uint64_t);
	}

	if ((fp = prof_open("index.tsv")) != NULL) {
		ri_idx_prof_dump(ri, fp);
		fclose(fp);
	}

	if (keys && atoi(keys) > 0) {
		ons_prof_slot_off = ri_idx_prof_slot_offsets(ri, &g_n_slots);
		ons_prof_keycnt = (uint32_t*)calloc(g_n_slots ? g_n_slots : 1, sizeof(uint32_t));
		g_self_bytes += (uint64_t)(g_n_buckets + 1) * sizeof(uint64_t) + g_n_slots * sizeof(uint32_t);
	}

	g_reads = prof_open("reads.tsv");
	if (g_reads)
		fputs("read_id\tname\tbatch\ttid\tl_sig\tn_norm\tn_peaks\tn_seg_ok\tn_seg_rej\tn_events\tbelow_min_events\t"
		      "n_seeds\tn_found\tn_single\tn_probe\tpost_found\tn_flt_occ\tpost_flt_occ\tn_skip_adm\tpost_skip_adm\t"
		      "post_visited\tn_pos_new\tn_cap_scan\tn_cap_evict\tn_sigs\tn_hits_kept\thits_bytes\tkm_cap\t"
		      "n_taxa\tcall\tchecked\t"
		      "t_read\tt_norm\tt_prefix\tt_tstat\tt_peaks\tt_events\tt_seeds\tt_lookup\tt_accum\tt_step1\tt_step1_cpu\t"
		      "t_summary\tt_bands\tt_pvalue\tt_slopes\tt_neighbor\tt_prepare\tt_commit\n", g_reads);
	g_batches = prof_open("batches.tsv");
	if (g_batches)
		fputs("batch\tn_reads\tn_samples\tt_start_s\tstep0_ns\tstep1_ns\tbusy_sum_ns\tbusy_max_ns\t"
		      "step2_par_ns\tstep2_ser_ns\tend_batch_ns\tstep2_ns\trss_bytes\tstep1_start_s\tstep2_start_s\n", g_batches);
}

void ons_prof_read_write(const ons_prof_read_t *r, const char *name)
{
	typedef unsigned long long ull;
	if (!g_reads || !r) return;
	g_reads_n++;
	g_samples += r->l_sig;
	g_events += r->n_events;
	g_seeds += r->n_seeds;
	fprintf(g_reads,
		"%llu\t%s\t%u\t%d\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t"
		"%u\t%u\t%u\t%llu\t%llu\t%u\t%llu\t%u\t%llu\t"
		"%llu\t%llu\t%llu\t%llu\t%u\t%llu\t%llu\t%llu\t"
		"%u\t%d\t%d\t"
		"%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t"
		"%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\n",
		(ull)r->read_id, name ? name : "*", r->batch, r->tid, r->l_sig, r->n_norm, r->n_peaks,
		r->n_seg_ok, r->n_seg_rej, r->n_events, r->below_min_events,
		r->n_seeds, r->n_found, r->n_single, (ull)r->n_probe, (ull)r->post_found, r->n_flt_occ,
		(ull)r->post_flt_occ, r->n_skip_adm, (ull)r->post_skip_adm,
		(ull)r->post_visited, (ull)r->n_pos_new, (ull)r->n_cap_scan, (ull)r->n_cap_evict, r->n_sigs,
		(ull)r->n_hits_kept, (ull)r->hits_bytes, (ull)r->km_cap,
		r->n_taxa, r->call, r->checked,
		(ull)r->t_read, (ull)r->t_norm, (ull)r->t_prefix, (ull)r->t_tstat, (ull)r->t_peaks, (ull)r->t_events,
		(ull)r->t_seeds, (ull)r->t_lookup, (ull)r->t_accum, (ull)r->t_step1, (ull)r->t_step1_cpu,
		(ull)r->t_summary, (ull)r->t_bands, (ull)r->t_pvalue, (ull)r->t_slopes, (ull)r->t_neighbor,
		(ull)r->t_prepare, (ull)r->t_commit);
}

void ons_prof_batch_write(uint32_t batch, uint32_t n_reads, uint64_t n_samples, uint64_t t_start, uint64_t step0,
                          uint64_t t1_start, uint64_t step1, const uint64_t *busy, int n_busy,
                          uint64_t t2_start, uint64_t step2_par, uint64_t step2_ser, uint64_t end_batch)
{
	typedef unsigned long long ull;
	uint64_t sum = 0, mx = 0;
	if (!g_batches) return;
	for (int i = 0; i < n_busy; ++i) {
		sum += busy[i];
		if (busy[i] > mx) mx = busy[i];
	}
	fprintf(g_batches, "%u\t%u\t%llu\t%.6f\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%.6f\t%.6f\n",
		batch, n_reads, (ull)n_samples, (double)(t_start - g_t_start) * 1e-9, (ull)step0, (ull)step1, (ull)sum,
		(ull)mx, (ull)step2_par, (ull)step2_ser, (ull)end_batch, (ull)(step2_par + step2_ser + end_batch),
		(ull)rss_bytes(), (double)(t1_start - g_t_start) * 1e-9, (double)(t2_start - g_t_start) * 1e-9);
}

static int cmp_u32_desc(const void *a, const void *b)
{
	uint32_t x = *(const uint32_t*)a, y = *(const uint32_t*)b;
	return (x < y) - (x > y);
}

static void write_hist(FILE *fp, const char *name, const uint64_t *h, int n)
{
	for (int i = 0; i < n; ++i)
		if (h[i]) fprintf(fp, "%s\t%d\t%llu\n", name, i, (unsigned long long)h[i]);
}

void ons_prof_finish(void)
{
	typedef unsigned long long ull;
	static uint64_t seg[ONS_HIST_SEG], post[ONS_HIST_POST], probe[ONS_HIST_PROBE], hits[ONS_HIST_HITS];
	uint64_t *bucket, total_lookups = 0;
	struct rusage ru;
	FILE *fp;
	int i, j;

	if (!ons_prof_on) return;
	bucket = (uint64_t*)calloc(g_n_buckets ? g_n_buckets : 1, sizeof(uint64_t));
	for (i = 0; i < g_n_slots_tls; ++i) {
		ons_prof_tls_t *t = g_slots[i];
		int k;
		for (k = 0; k < ONS_HIST_SEG; ++k) seg[k] += t->seg_len[k];
		for (k = 0; k < ONS_HIST_POST; ++k) post[k] += t->post[k];
		for (k = 0; k < ONS_HIST_PROBE; ++k) probe[k] += t->probe[k];
		for (k = 0; k < ONS_HIST_HITS; ++k) hits[k] += t->hits[k];
		for (k = 0; k < (int)g_n_buckets; ++k) bucket[k] += t->bucket[k];
	}
	for (i = 0; i < (int)g_n_buckets; ++i) total_lookups += bucket[i];

	if ((fp = prof_open("hist.tsv")) != NULL) {
		fputs("hist\tbin\tcount\n", fp);
		write_hist(fp, "segment_len", seg, ONS_HIST_SEG);
		write_hist(fp, "positions_per_found_seed", post, ONS_HIST_POST);
		write_hist(fp, "probe_steps_per_lookup", probe, ONS_HIST_PROBE);
		write_hist(fp, "hits_per_read_sig", hits, ONS_HIST_HITS);
		fclose(fp);
	}
	if ((fp = prof_open("buckets.tsv")) != NULL) {
		fputs("bucket\tlookups\n", fp);
		for (i = 0; i < (int)g_n_buckets; ++i)
			fprintf(fp, "%d\t%llu\n", i, (ull)bucket[i]);
		fclose(fp);
	}

	if (ons_prof_keycnt && (fp = prof_open("keys.tsv")) != NULL) {
		/* Coverage curve: the share of found-seed lookups the N hottest keys absorb. */
		uint64_t n_nz = 0, k, n_hit = 0, cum = 0, next = 1;
		uint32_t *nz;
		for (k = 0; k < g_n_slots; ++k) if (ons_prof_keycnt[k]) ++n_nz;
		nz = (uint32_t*)malloc((n_nz ? n_nz : 1) * sizeof(uint32_t));
		for (k = 0, j = 0; k < g_n_slots; ++k)
			if (ons_prof_keycnt[k]) { nz[j++] = ons_prof_keycnt[k]; n_hit += ons_prof_keycnt[k]; }
		qsort(nz, n_nz, sizeof(uint32_t), cmp_u32_desc);
		fprintf(fp, "# distinct_keys_touched\t%llu\n# found_key_lookups\t%llu\n", (ull)n_nz, (ull)n_hit);
		fputs("top_n_keys\tlookups_covered\tfraction\n", fp);
		for (k = 0; k < n_nz; ++k) {
			cum += nz[k];
			if (k + 1 == next || k + 1 == n_nz) {
				fprintf(fp, "%llu\t%llu\t%.6f\n", (ull)(k + 1), (ull)cum, n_hit ? (double)cum / (double)n_hit : 0.0);
				next <<= 1;
			}
		}
		free(nz);
		fclose(fp);
	}

	getrusage(RUSAGE_SELF, &ru);
	if ((fp = prof_open("summary.tsv")) != NULL) {
		fprintf(fp, "key\tvalue\n");
		fprintf(fp, "wall_s\t%.6f\n", (double)(ons_now() - g_t_start) * 1e-9);
		fprintf(fp, "user_cpu_s\t%.6f\n", ru.ru_utime.tv_sec + ru.ru_utime.tv_usec * 1e-6);
		fprintf(fp, "sys_cpu_s\t%.6f\n", ru.ru_stime.tv_sec + ru.ru_stime.tv_usec * 1e-6);
		fprintf(fp, "peak_rss_bytes\t%llu\n", (ull)ru.ru_maxrss * 1024ULL);
		fprintf(fp, "profiler_self_bytes\t%llu\n", (ull)g_self_bytes);
		fprintf(fp, "workers\t%d\nio_threads\t%d\npipeline_threads\t%d\n", g_n_threads, g_io_threads, g_pl_threads);
		fprintf(fp, "opt_batch_samples\t%lld\n", (long long)g_batch_size);
		fprintf(fp, "opt_mid_occ\t%d\nopt_sig_max_hash_occ\t%u\n", g_mid_occ, g_max_hash_occ);
		fprintf(fp, "reads\t%llu\nsamples\t%llu\nevents\t%llu\nseeds\t%llu\n",
			(ull)g_reads_n, (ull)g_samples, (ull)g_events, (ull)g_seeds);
		fprintf(fp, "index_lookups\t%llu\n", (ull)total_lookups);
		for (i = 0; i < ONS_KTP_MAXW; ++i)
			for (j = 0; j < ONS_KTP_MAXS; ++j)
				if (ons_ktp_wait[i][j] || ons_ktp_work[i][j])
					fprintf(fp, "pipeline_thread%d_step%d_work_s\t%.6f\npipeline_thread%d_step%d_wait_s\t%.6f\n",
						i, j, ons_ktp_work[i][j] * 1e-9, i, j, ons_ktp_wait[i][j] * 1e-9);
		fclose(fp);
	}

	if (g_reads) { fclose(g_reads); g_reads = NULL; }
	if (g_batches) { fclose(g_batches); g_batches = NULL; }
	for (i = 0; i < g_n_slots_tls; ++i) { free(g_slots[i]->bucket); free(g_slots[i]); }
	free(g_slots); g_slots = NULL;
	free(ons_prof_keycnt); ons_prof_keycnt = NULL;
	free(ons_prof_slot_off); ons_prof_slot_off = NULL;
	free(bucket);
	ons_prof_on = 0;
}

#endif /* ONS_PROF */
