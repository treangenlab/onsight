/* pipeline.c: the read pipeline: minibatches of reads, each read's steps on all threads, the output step; see
 * pipeline.h. */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif
#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 700
#endif

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pipeline.h"
#include "kalloc.h"
#include "kthread.h"
#include "rh_kvec.h"
#include "utils.h"
#include "event_detection.h"
#include "seeds.h"
#include "seed_hits.h"
#include "sig_summary.h"
#include "streaming.h"
#include "calls.h"
#include "onsprof.h"

/* One thread's memory: its memory pool, and with --neighbor-check the event levels of its reads of the minibatch,
 * kept for the output step */
typedef struct ri_tbuf_s {
	void *km;
	float *levels;
	uint64_t n_levels, m_levels;
} ri_tbuf_t;

typedef struct {
	int n_processed, n_threads, io_n_threads, n_fp, cur_fp, n_f, cur_f;
	int64_t mini_batch_size;
	const ri_detopt_t *opt;
	char **f;
	ri_signal_file_t *fp;
	const ri_idx_t *ri;
	const char **fn;
	FILE *sig_summary_fp;
	FILE *sig_profile_fp;
	ri_profile_t *profile;      /* --sig-profile */
	/* The stream layer (streaming.h): its state, NULL if off, and one scratch per output-step thread */
	ri_stream_state_t *stream_state;
	ri_stream_scratch_t **stream_scratch;
#ifdef ONS_PROF
	uint32_t prof_batch;
#endif
} pipeline_mt;

/* One read's output, made by the output step's threads and written in read order. */
typedef struct {
	char        synth_name[32];  /* the name of a read without one */
	int         synth;
	ri_sig_summary_t sum;        /* its --sig-summary row and its call for the profile */
	ri_read_scores_t *stream;    /* its stream-layer scores */
} ri_out1_t;

/* One minibatch */
typedef struct {
	const pipeline_mt *p;
	int n_sig;
	ri_signal_t** sig;
	ri_read_t** rd;
	ri_tbuf_t** buf;
	ri_out1_t* out;
#ifdef ONS_PROF
	uint32_t prof_batch;
	uint64_t prof_t_start, prof_step0, prof_t1_start, prof_step1;
	uint64_t *prof_busy; // per step-1 thread
#endif
} step_mt;

static ri_tbuf_t *ri_tbuf_init(void)
{
	ri_tbuf_t *b;
	b = (ri_tbuf_t*)calloc(1, sizeof(ri_tbuf_t));
	b->km = ri_km_init();
	return b;
}

static void ri_tbuf_destroy(ri_tbuf_t *b)
{
	if (b == 0) return;
	ri_km_destroy(b->km);
	free(b->levels);
	free(b);
}

/* One read, as one piece: its events, their seeds (strobe-events or contiguous), and the seeds' index positions added
 * to its hits per signature. A read with fewer than --min-events events gets no seeds. */
static void read_seed_hits(const ri_idx_t *ri,
				const uint32_t s_len,
				const float *sig,
				ri_read_t* rd,
				ri_tbuf_t *b,
				const ri_detopt_t *opt)
{
	uint32_t n_events = 0, n_events_sum = 0;
	double mean_sum = 0, std_dev_sum = 0;

	float* events = detect_events(b->km, s_len, sig,
		opt->window_length1, opt->window_length2, opt->threshold1, opt->threshold2, opt->peak_height,
		opt->min_segment_length, opt->max_segment_length, &mean_sum, &std_dev_sum, &n_events_sum, &n_events);

	#ifdef ONS_PROF
	ons_prof_read_t *pr = ons_prof_cur;
	if (pr) pr->n_events = n_events;
	#endif

	if(n_events < opt->min_events) {
		#ifdef ONS_PROF
		if (pr) pr->below_min_events = 1;
		#endif
		if(events){ri_kfree(b->km, events); events = NULL;}
		return;
	}

	/* --neighbor-check: keep the event levels (q_pos indexes them) in the thread's buffer until the stream step checks
	 * the read's call */
	if (opt->neighbor_regions && rd->n_levels == 0) {
		if (b->n_levels + n_events > b->m_levels) {
			uint64_t m = b->m_levels ? b->m_levels : 65536;
			while (m < b->n_levels + n_events) m <<= 1;
			float *v = (float*)realloc(b->levels, m * sizeof(float));
			if (v) { b->levels = v; b->m_levels = m; }
		}
		if (b->n_levels + n_events <= b->m_levels) {
			memcpy(b->levels + b->n_levels, events, n_events * sizeof(float));
			rd->levels_buf = b;
			rd->levels_off = b->n_levels;
			rd->n_levels = n_events;
			b->n_levels += n_events;
		}
	}

	#ifdef ONS_PROF
	uint64_t pt0 = pr ? ons_now() : 0;
	#endif

	// seeds
	mm128_v seeds = {0,0,0};
	if (ri->flag & RI_I_SEED_MODE_STROBE)
		ri_seeds_strobe(b->km, events, 0, 0, n_events, ri->diff, ri->e, ri->strobe_w_min, ri->strobe_w_max, ri->strobe_n, ri->q, ri->k, ri->fine_min, ri->fine_max, ri->fine_range, &seeds);
	else
		ri_seeds_contig(b->km, events, 0, 0, n_events, ri->diff, ri->w, ri->e, ri->n, ri->q, ri->k, ri->fine_min, ri->fine_max, ri->fine_range, &seeds);

	if(events){ri_kfree(b->km, events); events = NULL;}

	#ifdef ONS_PROF
	if (pr) { pr->t_seeds += ons_now() - pt0; pr->n_seeds = (uint32_t)seeds.n; }
	#endif

	//Seeding
	ri_collect_seed_hits(b->km, opt, opt->mid_occ, ri, rd, &seeds);
	if(seeds.a){ri_kfree(b->km, seeds.a); seeds.a = NULL; seeds.n = seeds.m = 0;}

	rd->offset += n_events;
}

static void hits_worker_for(void *_data,
						   long i,
						   int tid) // kt_for() callback
{
    step_mt *s = (step_mt*)_data; //s->sig and s->n_sig (signals read in this step and num of them)
	ri_tbuf_t* b = s->buf[tid];
	ri_read_t* rd = s->rd[i];
	ri_signal_t* sig = s->sig[i];
	rd->offset = 0;
	rd->levels = NULL, rd->n_levels = 0, rd->levels_off = 0, rd->levels_buf = NULL;
	rd->sig_hits = NULL, rd->n_sig_hits = 0, rd->m_sig_hits = 0, rd->hits_block = NULL;

	#ifdef ONS_PROF
	ons_prof_read_t *pr = (ons_prof_read_t*)rd->prof;
	uint64_t pt0 = 0, pc0 = 0;
	ons_prof_cur = pr;
	ons_prof_tls_set(tid);
	if (pr) {
		pr->batch = s->prof_batch, pr->tid = tid, pr->l_sig = sig->l_sig, pr->t_read = sig->prof_t_read;
		pt0 = ons_now(), pc0 = ons_cpu_now();
	}
	#endif

	if (sig->l_sig > 0) read_seed_hits(s->p->ri, sig->l_sig, sig->sig, rd, b, s->p->opt);
	rd->read_id = sig->rid;
	rd->read_name = sig->name;
	#ifdef ONS_PROF
	if (pr) {
		ons_prof_tls_t *pt = ons_prof_tls();
		ri_km_stat_t ks;
		pr->read_id = sig->rid;
		pr->n_sigs = rd->n_sig_hits;
		pr->hits_bytes = (uint64_t)rd->m_sig_hits * sizeof(ri_sig_hits_t);
		for (uint32_t m = 0; m < rd->n_sig_hits; ++m) {
			uint32_t nh = rd->sig_hits[m].n_hits;
			pr->n_hits_kept += nh;
			pr->hits_bytes += (uint64_t)rd->sig_hits[m].m_hits * sizeof(ri_seed_hit_t);
			if (pt) pt->hits[nh < ONS_HIST_HITS ? nh : ONS_HIST_HITS - 1]++;
		}
		ri_km_stat(b->km, &ks);
		pr->km_cap = ks.capacity;
	}
	#endif
	/* the thread's memory pool is kept for its next read: freeing it per read cost a page fault per 4 kB reused */
	#ifdef ONS_PROF
	if (pr) {
		pr->t_step1 = ons_now() - pt0;
		pr->t_step1_cpu = ons_cpu_now() - pc0;
		s->prof_busy[tid] += pr->t_step1;
	}
	ons_prof_cur = NULL;
	ons_prof_tls_set(-1);
	#endif
}

static ri_signal_t** read_signal_batch(pipeline_mt *pl,
							int64_t chunk_size,
							int *n_)
{	
	*n_ = 0;
	if (pl->n_fp < 1) return 0;

	int64_t size = 0;
	ri_signal_v rsigv = {0,0,0};
	rh_kv_resize(ri_signal_t*, 0, rsigv, 4000);

	while (pl->fp) {
		//Reading data in bulk if the buffer is emptied
		while(pl->fp && pl->fp->cur_read == pl->fp->num_read){
			ri_signal_file_close(pl->fp);
			if(pl->cur_f < pl->n_f){
				if((pl->fp = ri_signal_file_open(pl->f[pl->cur_f++], pl->io_n_threads)) == 0) break;
			}else if(pl->cur_fp < pl->n_fp){
				if(pl->f){
					for(int i = 0; i < pl->n_f; ++i) 
						if(pl->f[i]){free(pl->f[i]); pl->f[i] = NULL;}
					free(pl->f); pl->f = NULL;
				}
				pl->n_f = 0; pl->cur_f = 0;

				ri_char_v fnames = {0,0,0};
				rh_kv_resize(char*, 0, fnames, 256);
				ri_find_signal_files(pl->fn[pl->cur_fp++], &fnames);
				pl->f =  fnames.a;
				if(!fnames.n || ((pl->fp = ri_signal_file_open(pl->f[pl->cur_f++], pl->io_n_threads)) == 0)) break;
				pl->n_f = fnames.n;
			}else {pl->fp = 0; break;}
		}

		if(!pl->fp || pl->fp->cur_read == pl->fp->num_read) break;
		
		ri_signal_t *s = (ri_signal_t*)calloc(1, sizeof(ri_signal_t));
		rh_kv_push(ri_signal_t*, 0, rsigv, s);
		#ifdef ONS_PROF
		uint64_t pt0 = ons_now();
		#endif
		ri_read_signal(pl->fp, s);
		#ifdef ONS_PROF
		s->prof_t_read = ons_now() - pt0;
		#endif
		size += s->l_sig;

		if(size >= chunk_size) break;
	}

	ri_signal_t** a = 0;
	if(rsigv.n) a = rsigv.a;
	else rh_kv_destroy(rsigv);
	*n_ = rsigv.n;

	return a;
}

/* Output step, per read (kt_for callback): what depends only on the read and the slope grid of the minibatch, in
 * any order and thread: its --sig-summary row and call, and its stream scores. */
static void out_worker_for(void *_data, long k, int tid)
{
	step_mt *s = (step_mt*)_data;
	const pipeline_mt *p = s->p;
	const ri_idx_t *ri = p->ri;
	ri_read_t *rd = s->rd[k];
	ri_out1_t *o = &s->out[k];
	if (!rd) return;
	/* A read without a name gets one, so that it still has its rows (FAST5, POD5 and SLOW5 always give one). */
	if (!rd->read_name) {
		snprintf(o->synth_name, sizeof(o->synth_name), "read_%u", rd->read_id);
		rd->read_name = (const char*)o->synth_name;
		o->synth = 1;
	}
	if (rd->n_sig_hits == 0) {
		if (p->sig_summary_fp)
			o->sum.row = ri_sig_summary_empty_row(rd->read_name);
		return;
	}
	#ifdef ONS_PROF
	ons_prof_read_t *pr = (ons_prof_read_t*)rd->prof;
	uint64_t pt0 = pr ? ons_now() : 0;
	#endif
	/* the signatures ranked by score: only --sig-summary and --sig-profile use them */
	if (p->sig_summary_fp || p->profile) ri_sig_summary_read(ri, p->opt, rd, p->sig_summary_fp != NULL, &o->sum);
	#ifdef ONS_PROF
	if (pr) { uint64_t pt1 = ons_now(); pr->t_summary = pt1 - pt0; pt0 = pt1; }
	ons_prof_cur = pr;
	#endif
	if (p->stream_state) {
		if (rd->levels_buf) rd->levels = rd->levels_buf->levels + rd->levels_off;
		o->stream = ri_stream_prepare_read(p->stream_state, ri, rd, p->stream_scratch[tid]);
	}
	#ifdef ONS_PROF
	if (pr) pr->t_prepare = ons_now() - pt0;
	ons_prof_cur = NULL;
	#endif
	/* the hits are not needed after this point: the profile and the stream commit use the counts */
	free(rd->hits_block);
	rd->hits_block = NULL;
	for (uint32_t m = 0; m < rd->n_sig_hits; ++m) rd->sig_hits[m].hits = NULL;
}

static void *read_pipeline(void *shared,
								int step,
								void *in)
{
	int i, k;
    pipeline_mt *p = (pipeline_mt*)shared;
    if (step == 0) { // step 0: a minibatch of signals
        step_mt *s;
        s = (step_mt*)calloc(1, sizeof(step_mt));
		#ifdef ONS_PROF
		s->prof_t_start = ons_now();
		#endif

		s->sig = read_signal_batch(p, p->mini_batch_size, &s->n_sig);

		if (s->n_sig) {
			s->p = p;
			for (i = 0; i < s->n_sig; ++i)
				(*(s->sig[i])).rid = p->n_processed++;
			s->buf = (ri_tbuf_t**)calloc(p->n_threads, sizeof(ri_tbuf_t*));
			for (i = 0; i < p->n_threads; ++i)
				s->buf[i] = ri_tbuf_init();
			s->rd = (ri_read_t**)calloc(s->n_sig, sizeof(ri_read_t*));
			for(i = 0; i < s->n_sig; ++i)
				s->rd[i] = (ri_read_t*)calloc(1, sizeof(ri_read_t));
			#ifdef ONS_PROF
			s->prof_batch = p->prof_batch++;
			s->prof_busy = (uint64_t*)calloc(p->n_threads, sizeof(uint64_t));
			for (i = 0; i < s->n_sig; ++i) {
				ons_prof_read_t *pr = ons_prof_read_new();
				if (pr) pr->checked = -2;
				s->rd[i]->prof = pr;
			}
			s->prof_step0 = ons_now() - s->prof_t_start;
			#endif
			return s;
		} else if(s){
			if(s->sig) {free(s->sig); s->sig = NULL;}
			free(s); s = NULL;
		}
    } else if (step == 1) { // step 1: each read's seed hits
		step_mt *s = (step_mt*)in;
		#ifdef ONS_PROF
		s->prof_t1_start = ons_now();
		#endif
		kt_for(p->n_threads, hits_worker_for, in, s->n_sig);
		#ifdef ONS_PROF
		s->prof_step1 = ons_now() - s->prof_t1_start;
		#endif
		return in;
    } else if (step == 2) { // step 2: output
        step_mt *s = (step_mt*)in;
		const ri_idx_t *ri = p->ri;
		/* Per read, in parallel: the parts that depend only on the read and this minibatch's slope grid. */
		#ifdef ONS_PROF
		uint64_t pt2_start = ons_now(), pt2_par, pt2_ser, pt_end, pt0, n_samples = 0;
		#endif
		s->out = (ri_out1_t*)calloc(s->n_sig, sizeof(ri_out1_t));
		kt_for(p->n_threads, out_worker_for, s, s->n_sig);
		#ifdef ONS_PROF
		pt2_par = ons_now() - pt2_start;
		pt0 = ons_now();
		#endif
		/* In read order: the rows, the profile and the stream state. */
		for (k = 0; k < s->n_sig; ++k) {
			ri_read_t* rd = s->rd[k];
			ri_out1_t *o = &s->out[k];
			if (!rd) continue;
			#ifdef ONS_PROF
			ons_prof_read_t *pr = (ons_prof_read_t*)rd->prof;
			uint64_t pc = pr ? ons_now() : 0;
			#endif
			if (o->sum.row) {
				fputs(o->sum.row, p->sig_summary_fp);
				free(o->sum.row);
			}
			if (rd->n_sig_hits > 0) {
				if (p->profile)
					ri_profile_add_read(p->profile, ri, p->opt, rd, o->sum.ord, o->sum.call_taxon, o->sum.best_score);
				if (p->stream_state)
					ri_stream_commit_read(p->stream_state, rd, o->stream, (uint64_t)rd->read_id);
				free(o->sum.ord);
			}

			free(rd->hits_block);
			rd->hits_block = NULL;
			if(rd->sig_hits){free(rd->sig_hits); rd->sig_hits = NULL;}
			#ifdef ONS_PROF
			if (pr) {
				pr->t_commit = ons_now() - pc;
				n_samples += pr->l_sig;
				ons_prof_read_write(pr, rd->read_name);
				free(pr); rd->prof = NULL;
			}
			#endif
			if (o->synth) rd->read_name = NULL;
			free(rd); s->rd[k] = NULL;
		}
		free(s->out); s->out = NULL;
		#ifdef ONS_PROF
		pt2_ser = ons_now() - pt0;
		pt0 = ons_now();
		#endif
		/* Step 2 runs one minibatch at a time, in read order: the slope grid changes only here. */
		if (p->stream_state) ri_stream_end_batch(p->stream_state);
		#ifdef ONS_PROF
		pt_end = ons_now() - pt0;
		ons_prof_batch_write(s->prof_batch, (uint32_t)s->n_sig, n_samples, s->prof_t_start, s->prof_step0,
		                     s->prof_t1_start, s->prof_step1, s->prof_busy, p->n_threads, pt2_start, pt2_par, pt2_ser,
		                     pt_end);
		free(s->prof_busy); s->prof_busy = NULL;
		#endif

		for (i = 0; i < p->n_threads; ++i) ri_tbuf_destroy(s->buf[i]);
		if(s->buf){free(s->buf); s->buf = NULL;}
		if(s->rd){free(s->rd); s->rd = NULL;}

		for(int i = 0; i < s->n_sig; ++i){
			ri_signal_t *curS = s->sig[i];
			if(!curS) continue;
			if(curS->sig){free(curS->sig); curS->sig = NULL;}
			if(curS->name){free(curS->name); curS->name = NULL;}
			free(curS); s->sig[i] = NULL;
		}
		if(s->sig){free(s->sig); s->sig = NULL;}
		free(s); s = NULL;
	}
    return 0;
}

int ri_detect_files(const ri_idx_t *idx,
					int n_segs,
					const char **fn,
					const ri_detopt_t *opt,
					int n_threads,
					int io_n_threads)
{
	int pl_threads;
	pipeline_mt pl;
	if (n_segs < 1) return -1;
	memset(&pl, 0, sizeof(pipeline_mt));
	pl.n_fp = n_segs;
	pl.n_f = 0; pl.cur_f = 0;
	ri_char_v fnames = {0,0,0};
	rh_kv_resize(char*, 0, fnames, 256);
	ri_find_signal_files(fn[0], &fnames);
	pl.f =  fnames.a;
	pl.io_n_threads = io_n_threads > 1?io_n_threads:1;
	pl.n_threads = (n_threads - pl.io_n_threads > 1)?(n_threads - pl.io_n_threads):1;
	pl.mini_batch_size = opt->mini_batch_size;
	pl_threads = (pl.n_threads > 1)?2:1;
	if(!fnames.n || ((pl.fp = ri_signal_file_open(pl.f[0], pl.io_n_threads)) == 0)){rh_kv_destroy(fnames); return -1;}
	if (pl.fp == 0){rh_kv_destroy(fnames); return -1;}
	pl.fn = fn;
	pl.n_f = fnames.n;
	pl.cur_fp = 1;
	pl.cur_f = 1;
	pl.opt = opt, pl.ri = idx;
	if (opt->sig_summary_path) {
		pl.sig_summary_fp = fopen(opt->sig_summary_path, "w");
		if (!pl.sig_summary_fp) {
			fprintf(stderr, "[ERROR] failed to open signature summary file '%s'\n", opt->sig_summary_path);
			rh_kv_destroy(fnames);
			return -1;
		}
	}
	if (opt->sig_profile_path) {
		pl.sig_profile_fp = fopen(opt->sig_profile_path, "w");
		if (!pl.sig_profile_fp) {
			fprintf(stderr, "[ERROR] failed to open profiling file '%s'\n", opt->sig_profile_path);
			if (pl.sig_summary_fp) fclose(pl.sig_summary_fp);
			rh_kv_destroy(fnames);
			return -1;
		}
		pl.profile = ri_profile_init(idx);
	}
	if (!idx->sig_emitted_positions || !idx->sig_distinct_hashes) {
		fprintf(stderr, "[ERROR] signature scoring needs the emitted-position and distinct-hash metadata of the index (its .sigmeta.tsv)\n");
		if (pl.sig_summary_fp) fclose(pl.sig_summary_fp);
		if (pl.sig_profile_fp) fclose(pl.sig_profile_fp);
		ri_profile_destroy(pl.profile);
		rh_kv_destroy(fnames);
		return -1;
	}

	/* Read calls and target calls (calls.c): made before the stream, finished after the last read. */
	ri_calls_t *calls = NULL;

	{   /* the stream layer */
		ri_stream_config_t scfg;
		ri_stream_config_init(&scfg);
		scfg.emit_fd               = opt->stream_emit_fd;
		scfg.grid.min              = opt->stream_diag_slope_min;
		scfg.grid.max              = opt->stream_diag_slope_max;
		scfg.grid.step             = opt->stream_diag_slope_step;
		scfg.nonoverlap_gap        = opt->stream_nonoverlap_gap;
		scfg.slope_chains          = opt->stream_slope_chains;
		scfg.slope_alpha           = opt->stream_slope_alpha;
		scfg.slope_width           = opt->stream_slope_width;
		scfg.slope_width_sd        = opt->stream_slope_width_sd;
		scfg.slope_cover           = opt->stream_slope_cover;
		if (opt->read_calls_path || opt->target_calls_path || opt->stream_slope_chains > 0) {
			calls = ri_calls_init(idx, opt);
			if (!calls) {
				if (pl.sig_summary_fp) fclose(pl.sig_summary_fp);
				if (pl.sig_profile_fp) fclose(pl.sig_profile_fp);
				ri_profile_destroy(pl.profile);
				rh_kv_destroy(fnames);
				return -1;
			}
			scfg.calls = calls;
		}
		pl.stream_state = ri_stream_init(&scfg, idx);
		if (!pl.stream_state) {
			fprintf(stderr, "[WARN] the stream layer could not start; running without it\n");
		} else {
			pl.stream_scratch = (ri_stream_scratch_t**)calloc(pl.n_threads, sizeof(ri_stream_scratch_t*));
			for (int t = 0; t < pl.n_threads; ++t)
				pl.stream_scratch[t] = ri_stream_scratch_init(pl.stream_state);
		}
	}

	#ifdef ONS_PROF
	ons_prof_init(idx, opt, pl.n_threads, pl.io_n_threads, pl_threads);
	#endif

	kt_pipeline(pl_threads, read_pipeline, &pl, 3);

	if (pl.sig_profile_fp) {
		ri_profile_write(pl.profile, pl.sig_profile_fp);
		fclose(pl.sig_profile_fp);
	}
	if (pl.sig_summary_fp) fclose(pl.sig_summary_fp);
	ri_profile_destroy(pl.profile);
	if (pl.stream_state) {
		for (int t = 0; t < pl.n_threads; ++t)
			ri_stream_scratch_destroy(pl.stream_state, pl.stream_scratch[t]);
		free(pl.stream_scratch); pl.stream_scratch = NULL;
		ri_stream_destroy(pl.stream_state);
		pl.stream_state = NULL;
	}
	if (calls) {
		ri_calls_finish(calls, (uint64_t)pl.n_processed);
		ri_calls_destroy(calls);
	}

	rh_kv_destroy(fnames);

	#ifdef ONS_PROF
	ons_prof_finish();
	#endif

	return 0;
}

