#include <stdint.h>
#include <string.h> //for memset

#ifndef OPTIONS_H
#define OPTIONS_H

/* options.h: the settings of the index build (ri_idxopt_t) and of detection (ri_detopt_t), with their defaults. */

#define RI_CALLS_MAX_ALPHA 8   /* --alpha-read values */

/* RawHash2 index modes that onsight no longer builds; ri_idx_load refuses an index that has either */
#define RI_I_STORE_SIG	0x10
#define RI_I_SIG_TARGET	0x20
/* strobe-events (--seed-mode strobe): the index's e is the events per strobe, strobe_n, strobe_w_min and strobe_w_max
 * the rest; the index header has three more bytes. */
#define RI_I_SEED_MODE_STROBE 0x200

#ifdef __cplusplus
extern "C" {
#endif

// index-build and detection settings
typedef struct ri_idxopt_s{
	short b, w, e, n, q, k, flag, lev_col;
	int64_t mini_batch_size;
	uint64_t batch_size;
	const char *metadata_path;

	int diff;
	float fine_min, fine_max, fine_range;

	/* strobe-events (flag & RI_I_SEED_MODE_STROBE) */
	short strobe_n;       /* strobes per seed; only 2 */
	short strobe_w_min;   /* min gap (events) between strobe-1 last and strobe-2 first */
	short strobe_w_max;   /* max gap (events) between strobe-1 last and strobe-2 first */
} ri_idxopt_t;

typedef struct ri_detopt_s{

	//Seeding parameters: the repeat cutoff mid_occ is the (1 - mid_occ_frac) quantile of positions per hash in the
	//index, bounded to [min_mid_occ, max_mid_occ] (ri_detopt_update); seeds with more positions are ignored
	float mid_occ_frac;
	int32_t min_mid_occ, max_mid_occ;
	int32_t mid_occ;

	uint32_t min_events; // reads with fewer events get no seeds

	//Seed hits per signature (seed_hits.h) and --sig-summary (sig_summary.h)
	uint32_t sig_min_hits;
	uint32_t sig_min_distinct_hashes;
	int sig_weight_common_hashes;
	uint32_t sig_max_hash_occ;
	float sig_min_ratio;
	const char *sig_summary_path;
	const char *sig_profile_path;
	int32_t sig_hits_cap; // --sig-report-anchor-cap: hits kept per read and signature

	int64_t mini_batch_size; // -K: signal samples per minibatch

	//Event detector options
	uint32_t window_length1;
	uint32_t window_length2;
	float threshold1;
	float threshold2;
	float peak_height;

	//Segment length filtering for gen_events()
	uint32_t min_segment_length;  // minimum segment length (skip shorter segments, 0 = no min filter)
	uint32_t max_segment_length;  // maximum segment length (skip longer segments, default 500)

	/* The stream layer (streaming.h): the --stream-* options */
	int          stream_emit_fd;        /* --stream-emit-fd: the read records; -1 (default): none */
	float        stream_diag_slope_min;  /* the search grid of the diagonal bands' slopes, min:max:step */
	float        stream_diag_slope_max;
	float        stream_diag_slope_step;
	uint32_t     stream_nonoverlap_gap;  /* hits at least this many events apart count as non-overlapping */
	uint32_t     stream_slope_chains;    /* > 0: renew the slope grid from the last this many confident chains */
	double       stream_slope_alpha;     /* read p-value at most this makes a chain confident, default 1e-7 */
	float        stream_slope_width;     /* W of the grid step W / S, events, default 10 */
	int          stream_slope_width_sd;  /* 1: W is 4 x the median residual SD of the chains (on with --r10) */
	float        stream_slope_cover;     /* --slope-grid: share of the chains the grid spans, the middle ones; default 0.95 */

	/* Read calls (read p-value) and target calls, calls.h */
	const char  *calls_panel;           /* --panel: TSV with taxon and name; rows sharing a name are one entry */
	const char  *read_calls_path;       /* --read-calls: one row per read with a candidate entry */
	const char  *target_calls_path;     /* --target-calls: one row per entry, rule and test */
	double       alpha_read[RI_CALLS_MAX_ALPHA]; /* --alpha-read: single-cutoff target calls at these cutoffs */
	int          n_alpha_read;          /* 0 (default): the target call averaged over cutoffs t_min to 1e-3 */
	double       t_min;                 /* --t-min: the smallest cutoff of the averaged target call */
	double       alpha_run;             /* --alpha-run: chance of any false target call per run, default 0.01 */
	const char  *neighbor_regions;      /* --neighbor-check: region file (neighbor_regions.py) of the near-neighbor check */
	const char  *neighbor_pore_model;   /* the run's pore model (-p), which the near-neighbor check uses too */
} ri_detopt_t;

/* Fill in the defaults. */
void ri_idxopt_init(ri_idxopt_t *opt);
void ri_detopt_init(ri_detopt_t *opt);

#ifdef __cplusplus
}
#endif
#endif //OPTIONS_H
