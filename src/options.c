/* options.c: default settings; see options.h. */
#include "options.h"

void ri_idxopt_init(ri_idxopt_t *opt)
{
	memset(opt, 0, sizeof(ri_idxopt_t));
	opt->e = 6; opt->w = 2; opt->q = 4; opt->n = 0; opt->k = 6, opt->lev_col = 1;   /* e: 6 events per strobe; contiguous seeds 8 (main.c) */
	opt->b = 14;
	opt->diff = 0; // quantized sig-diff: skip if |quant_diff| <= diff; -1 disables filtering
	opt->mini_batch_size = 50000000;
	opt->batch_size = 4000000000ULL;
	opt->metadata_path = NULL; //--metadata

	opt->fine_min = -2.0f;
	opt->fine_max = 2.0f;
	opt->fine_range = 0.4;

	/* strobe-events (the default; --seed-mode contig for contiguous seeds): two strobes of e events, the second
	 * starting 8-20 events after the first ends */
	opt->flag |= RI_I_SEED_MODE_STROBE;
	opt->strobe_n = 2;
	opt->strobe_w_min = 8;
	opt->strobe_w_max = 20;
}

void ri_detopt_init(ri_detopt_t *opt)
{
	memset(opt, 0, sizeof(ri_detopt_t));

	//seeding
	opt->mid_occ_frac = 1e-2f; //--mid_occ_frac
	opt->min_mid_occ = 50; //--q-mid-occ [I1, I2]
	opt->max_mid_occ = 500000; //--q-mid-occ [I1, I2]

	opt->mini_batch_size = 4000000; //-K: the slope grid changes between minibatches

	opt->min_events = 50; //--min-events

	opt->sig_min_hits = 1; //--sig-min-hits
	opt->sig_min_distinct_hashes = 3; //--sig-min-distinct-hashes
	opt->sig_weight_common_hashes = 1; //--sig-weight-common-hashes
	opt->sig_max_hash_occ = 1000; //--sig-max-hash-occ
	opt->sig_min_ratio = 1.0f; //--sig-min-ratio
	opt->sig_profile_path = NULL; //--sig-profile
	opt->sig_hits_cap = 1024; //--sig-report-anchor-cap

	//Default options for event detection (R9.4.1; --r10 sets the R10.4.1 ones).
	opt->window_length1 = 3; //--seg-window-length1
    opt->window_length2 = 9; //--seg-window-length2
    opt->threshold1 = 4.0f; //--seg-threshold1
    opt->threshold2 = 3.5f; //--seg-threshold2
    opt->peak_height = 0.4f; //--seg-peak-height

	//Segment length filtering for gen_events
	opt->min_segment_length = 0;   // no minimum
	opt->max_segment_length = 500;

	/* The stream layer: the search grid 0.50 to 2.00 at 0.01, the slope grid renewed from the last 199 confident chains,
	 * hits at least 12 events apart counted as non-overlapping. */
	opt->stream_emit_fd            = -1;    /* no read records unless --stream-emit-fd */
	opt->stream_diag_slope_min     = 0.50f;
	opt->stream_diag_slope_max     = 2.00f;
	opt->stream_diag_slope_step    = 0.01f;
	opt->stream_nonoverlap_gap     = 12;    /* one seed's two strobes */
	opt->stream_slope_chains       = 199;
	opt->stream_slope_alpha        = 1e-7;
	opt->stream_slope_width        = 10.0f; /* the engine's band width */
	opt->stream_slope_width_sd     = 0;     /* --r10 sets it */
	opt->stream_slope_cover        = 0.95f; /* 5th smallest to 5th largest of 199: holds the next slope w.p. 190/200 */

	/* Read and target calls: off unless --read-calls or --target-calls. The target call averages over read p-value
	 * cutoffs t_min to 1e-3; alpha_run 0.01 is the chance of any false target call per run. */
	opt->calls_panel                 = NULL;
	opt->read_calls_path             = NULL;
	opt->target_calls_path           = NULL;
	opt->n_alpha_read                = 0;
	opt->t_min                       = 1e-8;  /* about 1 / the unrelated comparisons in the calibration data (100.5 M) */
	opt->alpha_run                   = 0.01;
	opt->neighbor_regions            = NULL;
	opt->neighbor_pore_model         = NULL;
	opt->neighbor_check_all          = 0;
}
