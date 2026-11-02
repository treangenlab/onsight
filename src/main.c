/* main.c: the command line: build an index from signatures, or detect targets in reads with an index. */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <unistd.h>
#include <sys/stat.h>

#ifndef NAN
#define NAN (0.0/0.0)
#endif

#include "utils.h"
#include "options.h"
#include "signal_input.h"
#include "index.h"
#include "pipeline.h"
#include "ketopt.h"
#include "version.h"   /* ONSIGHT_VERSION, written by src/Makefile from ../VERSION */

static ko_longopt_t long_options[] = {
	{ (char*)"level_column",     	ko_required_argument, 	300 },
	{ (char*)"q-mid-occ",     		ko_required_argument, 	301 },
	{ (char*)"mid_occ_frac",     	ko_required_argument, 	302 },
	{ (char*)"min-events",			ko_required_argument, 	303 },
	{ (char*)"seg-window-length1",	ko_required_argument, 	332 },
	{ (char*)"seg-window-length2",	ko_required_argument, 	333 },
	{ (char*)"seg-threshold1",		ko_required_argument, 	334 },
	{ (char*)"seg-threshold2",		ko_required_argument, 	335 },
	{ (char*)"seg-peak-height",		ko_required_argument, 	336 },
	{ (char*)"sig-diff",			ko_required_argument,	347 },
	{ (char*)"r10",					ko_no_argument, 		363 },
	{ (char*)"fine-min",			ko_required_argument, 	364 },
	{ (char*)"fine-max",			ko_required_argument, 	365 },
	{ (char*)"fine-range",			ko_required_argument, 	366 },
	{ (char*)"io-thread",			ko_required_argument, 	369 },
	{ (char*)"version",				ko_no_argument, 	  	371 },
	{ (char*)"help",				ko_no_argument, 	  	'h' },
	{ (char*)"debug",				ko_no_argument, 	  	372 },
	{ (char*)"min-seg-length",		ko_required_argument, 	375 },
	{ (char*)"max-seg-length",		ko_required_argument, 	376 },
	{ (char*)"sig-min-hits",		ko_required_argument, 	392 },
	{ (char*)"sig-min-ratio",		ko_required_argument, 	393 },
	{ (char*)"sig-summary",			ko_required_argument, 	394 },
	{ (char*)"sig-profile",			ko_required_argument, 	395 },
	{ (char*)"metadata",			ko_required_argument, 	396 },
	{ (char*)"sig-min-distinct-hashes", ko_required_argument, 397 },
	{ (char*)"sig-weight-common-hashes", ko_required_argument, 398 },
	{ (char*)"sig-max-hash-occ",	ko_required_argument, 	399 },
	{ (char*)"sig-report-anchor-cap", ko_required_argument, 408 },
	/* strobe-events */
	{ (char*)"seed-mode",           ko_required_argument,   410 },
	{ (char*)"strobe-ks",           ko_required_argument,   411 },
	{ (char*)"strobe-wmin",         ko_required_argument,   412 },
	{ (char*)"strobe-wmax",         ko_required_argument,   413 },
	{ (char*)"strobe-n",            ko_required_argument,   414 },
	/* streaming layer */
	{ (char*)"stream-emit-fd",          ko_required_argument, 428 },
	{ (char*)"stream-diag-slopes",      ko_required_argument, 451 },
	{ (char*)"stream-nonoverlap-gap",   ko_required_argument, 453 },
	{ (char*)"stream-slope-chains",     ko_required_argument, 470 },
	{ (char*)"stream-slope-alpha",      ko_required_argument, 471 },
	{ (char*)"stream-slope-width",      ko_required_argument, 472 },
	{ (char*)"stream-slope-width-sd",   ko_no_argument,       473 },
	{ (char*)"slope-grid",              ko_required_argument, 474 },
	/* Read calls (read p-value) and target calls, calls.c */
	{ (char*)"panel",                   ko_required_argument, 460 },
	{ (char*)"read-calls",              ko_required_argument, 461 },
	{ (char*)"target-calls",            ko_required_argument, 462 },
	{ (char*)"alpha-read",              ko_required_argument, 463 },
	{ (char*)"alpha-run",               ko_required_argument, 464 },
	{ (char*)"neighbor-check",          ko_required_argument, 465 },
	{ (char*)"t-min",                   ko_required_argument, 466 },
	/* chemistry, panel folder and output folder */
	{ (char*)"pore",                    ko_required_argument, 480 },
	{ (char*)"panel-dir",               ko_required_argument, 481 },
	{ (char*)"out-dir",                 ko_required_argument, 'o' },
	{ 0, 0, 0 }
};

/* R10.4.1: 9-mer pore model, the R10.4.1 event segmentation and W of the slope step from the chains' residual SD */
static void r10_preset(ri_idxopt_t *ipt, ri_detopt_t *opt)
{
	ipt->k = 9;
	opt->window_length1 = 3; opt->window_length2 = 6;
	opt->threshold1 = 5.5f; opt->threshold2 = 3.5f;
	opt->peak_height = 0.2f;
	opt->stream_slope_width_sd = 1;   /* W = 4 x the chains' median residual SD, as the R10.4.1 results */
}

static char *path_join(const char *dir, const char *name)
{
	size_t n = strlen(dir) + strlen(name) + 2;
	char *p = (char*)malloc(n);
	snprintf(p, n, "%s/%s", dir, name);
	return p;
}

/* The pore model shipped with the engine: <bin>/../share/onsight/ after make install, <bin>/../extern/local_kmer_models/
 * in the source tree. 0 if neither has it. */
static char *bundled_pore_model(int r10)
{
	const char *name = r10 ? "uncalled_r1041_model_only_means.txt" : "r94_means_only.txt";
	const char *dirs[] = { "../share/onsight", "../extern/local_kmer_models" };
	char exe[4096], *slash, *dir, *p;
	ssize_t n = readlink("/proc/self/exe", exe, sizeof(exe) - 1);
	int i;
	if (n <= 0) return 0;
	exe[n] = 0;
	if (!(slash = strrchr(exe, '/'))) return 0;
	*slash = 0;
	for (i = 0; i < 2; ++i) {
		dir = path_join(exe, dirs[i]);
		p = path_join(dir, name);
		free(dir);
		if (access(p, R_OK) == 0) return p;
		free(p);
	}
	return 0;
}

/* mkdir -p */
static int make_dirs(const char *path)
{
	char *p = strdup(path), *q;
	int r;
	for (q = p + 1; *q; ++q)
		if (*q == '/') { *q = 0; mkdir(p, 0777); *q = '/'; }
	r = mkdir(p, 0777);
	free(p);
	return r == 0 || errno == EEXIST ? 0 : -1;
}

static inline int64_t mm_parse_num(const char *str)
{
	double x;
	char *p;
	x = strtod(str, &p);
	if (*p == 'G' || *p == 'g') x *= 1e9, ++p;
	else if (*p == 'M' || *p == 'm') x *= 1e6, ++p;
	else if (*p == 'K' || *p == 'k') x *= 1e3, ++p;
	return (int64_t)(x + .499);
}

static inline int parse_yes_no_int(const char *arg, int *out)
{
	if (strcmp(arg, "yes") == 0 || strcmp(arg, "y") == 0 || strcmp(arg, "true") == 0 || strcmp(arg, "1") == 0) {
		*out = 1;
		return 0;
	}
	if (strcmp(arg, "no") == 0 || strcmp(arg, "n") == 0 || strcmp(arg, "false") == 0 || strcmp(arg, "0") == 0) {
		*out = 0;
		return 0;
	}
	return -1;
}

int main(int argc, char *argv[])
{
	const char *opt_str = "k:d:p:e:q:w:t:K:o:h";
	ketopt_t o = KETOPT_INIT;
	ri_detopt_t opt;
  	ri_idxopt_t ipt;
	int c, n_threads = 3, io_n_threads = 1, e_given = 0, chem = 0, n_reads, debug = 0;
	char *fnw = 0, *fpore = 0, *s, *idx_path;
	const char *panel_dir = 0, *out_dir = 0, **reads;
	FILE *fp_help = stderr;
	ri_idx_reader_t *idx_rdr;
	ri_idx_t *ri;

	ri_verbose = 3;
	liftrlimit();
	ri_realtime0 = ri_realtime();
	ri_idxopt_init(&ipt);
	ri_detopt_init(&opt);

	// first pass: stop on an unknown option or a missing argument before anything runs
	while ((c = ketopt(&o, argc, argv, 1, opt_str, long_options)) >= 0) {
		if (c == ':') {
			fprintf(stderr, "[ERROR] missing option argument\n");
			return 1;
		} else if (c == '?') {
			fprintf(stderr, "[ERROR] unknown option in \"%s\"\n", argv[o.i - 1]);
			return 1;
		}
	}
	o = KETOPT_INIT;

	while ((c = ketopt(&o, argc, argv, 1, opt_str, long_options)) >= 0) {
		if (c == 'd') fnw = o.arg;
		else if (c == 'p') fpore = o.arg;
		else if (c == 'k') ipt.k = atoi(o.arg);
		else if (c == 'e') { ipt.e = atoi(o.arg); e_given = 1; }
		else if (c == 'q') ipt.q = atoi(o.arg);
		else if (c == 'w') ipt.w = atoi(o.arg);
		else if (c == 't') n_threads = atoi(o.arg);
		else if (c == 'K') {opt.mini_batch_size = mm_parse_num(o.arg);}
		else if (c == 'h') fp_help = stdout;
		else if (c == 300) ipt.lev_col = atoi(o.arg);// --level_column
		else if (c == 301) { //--q-mid-occ
			opt.min_mid_occ = strtol(o.arg, &s, 10); // min
			if (*s == ',') opt.max_mid_occ = strtol(s + 1, &s, 10); //max
			// opt.q_mid_occ = atoi(o.arg);// --q-mid-occ
		}
		else if (c == 302) opt.mid_occ_frac = atof(o.arg);// --occ-frac
		else if (c == 303) opt.min_events = (uint32_t)atoi(o.arg); // --min-events
		else if (c == 332) opt.window_length1 = atoi(o.arg);// --seg-window-length1
		else if (c == 333) opt.window_length2 = atoi(o.arg);// --seg-window-length2
		else if (c == 334) opt.threshold1 = atof(o.arg);// --seg-threshold1
		else if (c == 335) opt.threshold2 = atof(o.arg);// --seg-threshold2
		else if (c == 336) opt.peak_height = atof(o.arg);// --seg-peak-height
		else if (c == 347) {ipt.diff = atoi(o.arg);} // --sig-diff INT
		else if (c == 363) { chem = 10; r10_preset(&ipt, &opt); } // --r10, the same as --pore r10
		else if (c == 480) { // --pore r9|r10
			if (strcmp(o.arg, "r9") == 0) { chem = 9; ipt.k = 6; }
			else if (strcmp(o.arg, "r10") == 0) { chem = 10; r10_preset(&ipt, &opt); }
			else {
				fprintf(stderr, "[ERROR] --pore takes r9 (R9.4.1) or r10 (R10.4.1), got '%s'\n", o.arg);
				return 1;
			}
		}
		else if (c == 481) panel_dir = o.arg; // --panel-dir
		else if (c == 'o') out_dir = o.arg;   // -o, --out-dir
		else if (c == 364) {ipt.fine_min = atof(o.arg);}// --fine-min
		else if (c == 365) {ipt.fine_max = atof(o.arg);}// --fine-max
		else if (c == 366) {ipt.fine_range = atof(o.arg);}// --fine-range
		else if (c == 369) {io_n_threads = atoi(o.arg);}// --io-thread
		else if (c == 371) {printf("onsight %s\n", ONSIGHT_VERSION); return 0;}// --version
		else if (c == 372) debug = 1; // --debug
		else if (c == 375) opt.min_segment_length = (uint32_t)atoi(o.arg); // --min-seg-length
		else if (c == 376) opt.max_segment_length = (uint32_t)atoi(o.arg); // --max-seg-length
		else if (c == 392) opt.sig_min_hits = (uint32_t)atoi(o.arg); // --sig-min-hits
		else if (c == 393) opt.sig_min_ratio = atof(o.arg); // --sig-min-ratio
		else if (c == 394) opt.sig_summary_path = o.arg; // --sig-summary
		else if (c == 395) opt.sig_profile_path = o.arg; // --sig-profile
		else if (c == 396) ipt.metadata_path = o.arg; // --metadata
		else if (c == 397) opt.sig_min_distinct_hashes = (uint32_t)atoi(o.arg); // --sig-min-distinct-hashes
		else if (c == 398) {
			if (parse_yes_no_int(o.arg, &opt.sig_weight_common_hashes) != 0) {
				fprintf(stderr, "[ERROR] --sig-weight-common-hashes only accepts yes/no/true/false/1/0\n");
				return 1;
			}
		}
		else if (c == 399) opt.sig_max_hash_occ = (uint32_t)atoi(o.arg); // --sig-max-hash-occ
		else if (c == 408) { // --sig-report-anchor-cap
			int32_t v = (int32_t)atoi(o.arg);
			if (v < 1) {
				fprintf(stderr, "[ERROR] --sig-report-anchor-cap must be >= 1 (got %s)\n", o.arg);
				return 1;
			}
			opt.sig_hits_cap = v;
		}
		else if (c == 410) { // --seed-mode {contig|strobe}
			if (strcmp(o.arg, "contig") == 0 || strcmp(o.arg, "contiguous") == 0) {
				ipt.flag &= ~RI_I_SEED_MODE_STROBE;
			} else if (strcmp(o.arg, "strobe") == 0 || strcmp(o.arg, "strobemer") == 0 || strcmp(o.arg, "randstrobe") == 0) {
				ipt.flag |= RI_I_SEED_MODE_STROBE;
			} else {
				fprintf(stderr, "[ERROR] --seed-mode expects 'contig' or 'strobe', got '%s'\n", o.arg);
				return 1;
			}
		}
		else if (c == 411) { ipt.e = atoi(o.arg); e_given = 1; } // --strobe-ks (as -e)
		else if (c == 412) { ipt.strobe_w_min = (short)atoi(o.arg); } // --strobe-wmin
		else if (c == 413) { ipt.strobe_w_max = (short)atoi(o.arg); } // --strobe-wmax
		else if (c == 414) { ipt.strobe_n   = (short)atoi(o.arg); } // --strobe-n
		else if (c == 428) { opt.stream_emit_fd = atoi(o.arg); }                       // --stream-emit-fd
		else if (c == 451) {                                                           // --stream-diag-slopes
			float lo = 1.0f, hi = 1.0f, step = 0.0f;
			int n = sscanf(o.arg, "%f:%f:%f", &lo, &hi, &step);
			if (n == 1) { hi = lo; step = 0.0f; }
			if (!(n == 1 || n == 3) || lo <= 0.0f || hi < lo || (hi > lo && step <= 0.0f)) {
				fprintf(stderr, "[ERROR] --stream-diag-slopes takes MIN:MAX:STEP or one slope (got %s)\n", o.arg);
				return 1;
			}
			opt.stream_diag_slope_min = lo; opt.stream_diag_slope_max = hi; opt.stream_diag_slope_step = step;
		}
		else if (c == 453) { opt.stream_nonoverlap_gap = (uint32_t)atoi(o.arg); }      // --stream-nonoverlap-gap
		else if (c == 470) { opt.stream_slope_chains = (uint32_t)atoi(o.arg); }        // --stream-slope-chains
		else if (c == 471) { opt.stream_slope_alpha = atof(o.arg); }                   // --stream-slope-alpha
		else if (c == 472) { opt.stream_slope_width = (float)atof(o.arg); }            // --stream-slope-width
		else if (c == 473) { opt.stream_slope_width_sd = 1; }                          // --stream-slope-width-sd
		else if (c == 474) { opt.stream_slope_cover = (float)atof(o.arg); }            // --slope-grid
		else if (c == 460) { opt.calls_panel = o.arg; }                                // --panel
		else if (c == 461) { opt.read_calls_path = o.arg; }                            // --read-calls
		else if (c == 462) { opt.target_calls_path = o.arg; }                          // --target-calls
		else if (c == 463) {                                                           // --alpha-read A[,B...]
			const char *q = o.arg;
			char *end;
			opt.n_alpha_read = 0;
			while (*q) {
				double a = strtod(q, &end);
				if (end == q || !(a > 0.0 && a < 1.0) || opt.n_alpha_read >= RI_CALLS_MAX_ALPHA || (*end && *end != ',')) {
					fprintf(stderr, "[ERROR] --alpha-read takes up to %d values in (0, 1), comma separated (got %s)\n",
					        RI_CALLS_MAX_ALPHA, o.arg);
					return 1;
				}
				opt.alpha_read[opt.n_alpha_read++] = a;
				q = *end ? end + 1 : end;
			}
		}
		else if (c == 465) { opt.neighbor_regions = o.arg; }                           // --neighbor-check
		else if (c == 466) {                                                           // --t-min
			opt.t_min = atof(o.arg);
			if (!(opt.t_min > 0.0 && opt.t_min < 1e-3)) {
				fprintf(stderr, "[ERROR] --t-min takes a value in (0, 1e-3) (got %s)\n", o.arg);
				return 1;
			}
		}
		else if (c == 464) {                                                           // --alpha-run
			opt.alpha_run = atof(o.arg);
			if (!(opt.alpha_run > 0.0 && opt.alpha_run < 1.0)) {
				fprintf(stderr, "[ERROR] --alpha-run takes a value in (0, 1) (got %s)\n", o.arg);
				return 1;
			}
		}
	}

	/* The index is the first file and the reads the rest, unless --panel-dir gives the index: then every file is reads */
	idx_path = argc > o.ind ? argv[o.ind] : 0;
	reads = (const char**)&argv[o.ind + 1];
	n_reads = argc - o.ind - 1;
	if (argc > o.ind && fp_help != stdout) {
		if (!chem) {
			fprintf(stderr, "[ERROR] give the pore: --pore r9 (R9.4.1) or --pore r10 (R10.4.1)\n");
			return 1;
		}
		if (!fpore && !(fpore = bundled_pore_model(chem == 10))) {
			fprintf(stderr, "[ERROR] no %s pore model found next to the engine; give one with -p\n", chem == 10 ? "R10.4.1" : "R9.4.1");
			return 1;
		}
		if (panel_dir) {   /* <dir>/panel.tsv, <dir>/r9|r10/index.idx, <dir>/regions.lcr; flags given on the command line win */
			idx_path = path_join(panel_dir, chem == 10 ? "r10/index.idx" : "r9/index.idx");
			if (access(idx_path, R_OK) != 0) {
				fprintf(stderr, "[ERROR] the panel folder has no %s index: %s\n", chem == 10 ? "r10" : "r9", idx_path);
				return 1;
			}
			if (!opt.calls_panel) opt.calls_panel = path_join(panel_dir, "panel.tsv");
			if (!opt.neighbor_regions) {
				char *reg = path_join(panel_dir, "regions.lcr");
				if (access(reg, R_OK) == 0) opt.neighbor_regions = reg;
				else fprintf(stderr, "[WARNING] no regions.lcr in %s: read calls are not checked against near neighbors\n", panel_dir);
			}
			reads = (const char**)&argv[o.ind];
			n_reads = argc - o.ind;
		}
		if (out_dir) {
			if (make_dirs(out_dir) != 0) {
				fprintf(stderr, "[ERROR] cannot make the output folder %s: %s\n", out_dir, strerror(errno));
				return 1;
			}
			if (!opt.read_calls_path) opt.read_calls_path = path_join(out_dir, "read_calls.tsv");
			if (!opt.target_calls_path) opt.target_calls_path = path_join(out_dir, "target_calls.tsv");
			if (debug && !opt.sig_summary_path) opt.sig_summary_path = path_join(out_dir, "sig_summary.tsv");
		}
		if (debug && !opt.sig_summary_path) {
			fprintf(stderr, "[ERROR] --debug needs -o DIR (or --sig-summary)\n");
			return 1;
		}
		if (panel_dir && !opt.read_calls_path && !opt.target_calls_path) {
			fprintf(stderr, "[ERROR] --panel-dir needs -o DIR (or --read-calls / --target-calls) for the results\n");
			return 1;
		}
	}

	if (!(ipt.flag & RI_I_SEED_MODE_STROBE) && !e_given) ipt.e = 8;   /* contiguous seeds: 8 events, RawHash2's */

	if (opt.stream_slope_chains > 0 && (opt.stream_nonoverlap_gap == 0 || !(opt.stream_slope_width > 0.0f))) {
		fprintf(stderr, "[ERROR] --stream-slope-chains needs --stream-nonoverlap-gap (the read p-value counts non-overlapping hits) and --stream-slope-width > 0\n");
		return 1;
	}

	/* strobe-events: checked before the index build starts (two strobes, 1 <= w_min <= w_max) */
	if (ipt.flag & RI_I_SEED_MODE_STROBE) {
		if (ipt.strobe_n != 2) {
			fprintf(stderr, "[ERROR] --strobe-n must be 2 (got %d)\n", ipt.strobe_n);
			return 1;
		}
		if (ipt.strobe_w_min < 1 || ipt.strobe_w_max < ipt.strobe_w_min) {
			fprintf(stderr, "[ERROR] strobe gap window invalid: w_min=%d w_max=%d (need 1<=w_min<=w_max)\n", ipt.strobe_w_min, ipt.strobe_w_max);
			return 1;
		}
		if (ipt.e < 1 || (uint32_t)ipt.e * (uint32_t)ipt.q > 32U) {
			fprintf(stderr, "[ERROR] strobe length k_s=%d with q=%d exceeds 32-bit per-strobe entropy budget\n", ipt.e, ipt.q);
			return 1;
		}
	}

	if ((opt.read_calls_path || opt.target_calls_path) && opt.stream_nonoverlap_gap == 0) {
		fprintf(stderr, "[ERROR] --read-calls and --target-calls need --stream-nonoverlap-gap > 0\n");
		return 1;
	}
	if (opt.neighbor_regions && (!(opt.read_calls_path || opt.target_calls_path) || !fpore)) {
		fprintf(stderr, "[ERROR] --neighbor-check needs --read-calls or --target-calls, and the pore model (-p)\n");
		return 1;
	}
	opt.neighbor_pore_model = fpore;

	if (argc == o.ind || fp_help == stdout) {
		/* the name we were invoked as */
		const char *invoked = (argv[0] && *argv[0]) ? argv[0] : "onsight";
		const char *slash = strrchr(invoked, '/');
		if (slash && slash[1]) invoked = slash + 1;
		fprintf(fp_help, "Usage: %s --pore r9|r10 --panel-dir DIR -o OUT [options] <reads.fast5|.pod5|.slow5|.blow5|folder> [...]\n", invoked);
		fprintf(fp_help, "       %s --pore r9|r10 [options] <index.idx> <reads> [...]\n", invoked);
		fprintf(fp_help, "       %s --pore r9|r10 -d <index.idx> --metadata <signatures.csv> <signatures.fa>\n", invoked);
		fprintf(fp_help, "  Detect the panel's targets in the reads, or (-d, no reads) build an index. The defaults are the settings of the\n");
		fprintf(fp_help, "  results (strobe-events, the sliding-window slope range, the target call averaged over read p-value cutoffs).\n");

		fprintf(fp_help, "\n  General:\n");
		fprintf(fp_help, "    --pore STR   the reads' pore, r9 (R9.4.1) or r10 (R10.4.1); required. Sets -k, the event segmentation\n");
		fprintf(fp_help, "                 and the pore model shipped with the engine\n");
		fprintf(fp_help, "    --panel-dir DIR  the panel folder: DIR/panel.tsv, DIR/r9/index.idx and DIR/r10/index.idx (each with its\n");
		fprintf(fp_help, "                 .sigmeta.tsv), DIR/regions.lcr (the near-neighbor check, if there); --panel, --neighbor-check win\n");
		fprintf(fp_help, "    -o DIR       write read_calls.tsv and target_calls.tsv in DIR; --read-calls, --target-calls win\n");
		fprintf(fp_help, "    --debug      also write sig_summary.tsv in DIR: each read's top five signatures. It sorts every read's\n");
		fprintf(fp_help, "                 signatures, so the run is slower; read and target calls do not change\n");
		fprintf(fp_help, "    -h, --help   show this help message\n");
		fprintf(fp_help, "    --version    show version number\n");
		fprintf(fp_help, "    -t INT       number of threads [%d]\n", n_threads);
		fprintf(fp_help, "    --io-thread INT   I/O threads for S/BLOW5 (must be < -t) [%d]\n", io_n_threads);
		fprintf(fp_help, "    -K NUM       minibatch size in signal samples; the slope grid changes between minibatches [4M]\n");

		fprintf(fp_help, "\n  Index build:\n");
		fprintf(fp_help, "    -d FILE      write the index to FILE (and FILE.sigmeta.tsv)\n");
		fprintf(fp_help, "    -p FILE      pore model, also read by --neighbor-check [the --pore one shipped with the engine]\n");
		fprintf(fp_help, "    -k INT       k-mer size of the pore model [6 with --pore r9, 9 with --pore r10]\n");
		fprintf(fp_help, "    --level_column INT   0-based column of the mean levels in the pore model [%d]\n", ipt.lev_col);
		fprintf(fp_help, "    --r10        the same as --pore r10\n");
		fprintf(fp_help, "    --metadata FILE      signature metadata CSV (signature name, source taxon)\n");
		fprintf(fp_help, "    --seed-mode STR      strobe (strobe-events) or contig [%s]\n", (ipt.flag & RI_I_SEED_MODE_STROBE) ? "strobe" : "contig");
		fprintf(fp_help, "    --strobe-ks INT      events per strobe (as -e) [%d]\n", ipt.e);
		fprintf(fp_help, "    --strobe-wmin INT    min gap in events between the strobes [%d]\n", ipt.strobe_w_min);
		fprintf(fp_help, "    --strobe-wmax INT    max gap in events between the strobes [%d]\n", ipt.strobe_w_max);
		fprintf(fp_help, "    --strobe-n INT       strobes per seed (only 2) [%d]\n", ipt.strobe_n);
		fprintf(fp_help, "    -e INT       events per contiguous seed (--seed-mode contig) [8]\n");
		fprintf(fp_help, "    -w INT       minimizer window of contiguous seeds; 0 keeps every seed [%d]\n", ipt.w);
		fprintf(fp_help, "    -q INT       quantization bits [%d]\n", ipt.q);
		fprintf(fp_help, "    --fine-min FLOAT, --fine-max FLOAT, --fine-range FLOAT   quantization range [%g, %g, %g]\n", ipt.fine_min, ipt.fine_max, ipt.fine_range);
		fprintf(fp_help, "    --sig-diff INT       skip an event whose quantized value differs from the last by <= INT; -1 keeps all [%d]\n", ipt.diff);

		fprintf(fp_help, "\n  Event detection:\n");
		fprintf(fp_help, "    (set by --pore; these override it)\n");
		fprintf(fp_help, "    --seg-window-length1 INT  short t-test window [%u]\n", opt.window_length1);
		fprintf(fp_help, "    --seg-window-length2 INT  long t-test window [%u]\n", opt.window_length2);
		fprintf(fp_help, "    --seg-threshold1 FLOAT    peak threshold of the short window [%g]\n", opt.threshold1);
		fprintf(fp_help, "    --seg-threshold2 FLOAT    peak threshold of the long window [%g]\n", opt.threshold2);
		fprintf(fp_help, "    --seg-peak-height FLOAT   min peak height [%g]\n", opt.peak_height);
		fprintf(fp_help, "    --min-seg-length INT      skip events shorter than INT samples [%u]\n", opt.min_segment_length);
		fprintf(fp_help, "    --max-seg-length INT      skip events longer than INT samples [%u]\n", opt.max_segment_length);
		fprintf(fp_help, "    --min-events INT          skip reads with fewer events [%u]\n", opt.min_events);

		fprintf(fp_help, "\n  Seed hits per signature:\n");
		fprintf(fp_help, "    --mid_occ_frac FLOAT      repeat cutoff: ignore seeds whose hash has more positions than this top fraction [%g]\n", opt.mid_occ_frac);
		fprintf(fp_help, "    --q-mid-occ INT1[,INT2]   bounds of the repeat cutoff [%d,%d]\n", opt.min_mid_occ, opt.max_mid_occ);
		fprintf(fp_help, "    --sig-weight-common-hashes STR  weight each hit by 1 / log2(positions + 1): yes/no [%s]\n",
			opt.sig_weight_common_hashes ? "yes" : "no");
		fprintf(fp_help, "    --sig-max-hash-occ INT    also ignore hashes with more positions than INT; 0 = off [%u]\n", opt.sig_max_hash_occ);
		fprintf(fp_help, "    --sig-report-anchor-cap INT   hits kept per read and signature, the heaviest when full [%d]\n",
			opt.sig_hits_cap);
		fprintf(fp_help, "    --sig-min-hits INT        min hits for a signature call [%u]\n", opt.sig_min_hits);
		fprintf(fp_help, "    --sig-min-distinct-hashes INT  min distinct hashes for a signature call [%u]\n", opt.sig_min_distinct_hashes);
		fprintf(fp_help, "    --sig-min-ratio FLOAT     min best/second score ratio for a signature call [%g]\n", opt.sig_min_ratio);
		fprintf(fp_help, "    --sig-summary FILE        one row per read: its best signature, the top five and their scores\n");
		fprintf(fp_help, "    --sig-profile FILE        per-taxon summary at the end of the run\n");

		fprintf(fp_help, "\n  Stream layer (diagonal bands, slope range):\n");
		fprintf(fp_help, "    --stream-emit-fd INT      write the read records (JSON lines: read, slopes) to file descriptor INT [none]\n");
		fprintf(fp_help, "    --stream-diag-slopes MIN:MAX:STEP  search grid of the bands' slopes, or one slope [%.2f:%.2f:%.2f]\n",
			opt.stream_diag_slope_min, opt.stream_diag_slope_max, opt.stream_diag_slope_step);
		fprintf(fp_help, "    --stream-nonoverlap-gap INT  count the hits at least INT events apart in each band; 0 = off [%u]\n", opt.stream_nonoverlap_gap);
		fprintf(fp_help, "    --stream-slope-chains INT  renew the slopes from the last INT confident chains; 0 = off [%u]\n", opt.stream_slope_chains);
		fprintf(fp_help, "    --stream-slope-alpha FLOAT  read p-value at most this makes a chain confident [%g]\n", opt.stream_slope_alpha);
		fprintf(fp_help, "    --stream-slope-width FLOAT  W of the slope step W / S, in events [%g]\n", opt.stream_slope_width);
		fprintf(fp_help, "    --slope-grid FLOAT  share of the chains the grid spans, the middle ones [%g]\n", opt.stream_slope_cover);
		fprintf(fp_help, "    --stream-slope-width-sd   W is 4 x the median residual SD of the chains [on with --r10]\n");

		fprintf(fp_help, "\n  Read calls and target calls:\n");
		fprintf(fp_help, "    --panel FILE          panel targets: TSV with columns taxon and name; taxa sharing a name are one target\n");
		fprintf(fp_help, "                          [each index taxon is a target]\n");
		fprintf(fp_help, "    --read-calls FILE     one row per read with a candidate target: its read p-value and counts\n");
		fprintf(fp_help, "    --target-calls FILE   one row per target and test: target call, first call, test statistic\n");
		fprintf(fp_help, "    --t-min FLOAT         target call averaged over read p-value cutoffs t-min to 1e-3 [%.4g]\n", opt.t_min);
		fprintf(fp_help, "    --alpha-read A[,B]    instead, a target call at each of these single cutoffs\n");
		fprintf(fp_help, "    --alpha-run FLOAT     chance of any false target call in a run [%g]\n", opt.alpha_run);
		fprintf(fp_help, "    --neighbor-check FILE near-neighbor check of the read calls at p <= 1e-3 (or the loosest alpha-read)\n");
		fprintf(fp_help, "                          against the regions in FILE (neighbor_regions.py), with the -p pore model;\n");
		fprintf(fp_help, "                          the target call then also runs on the read calls it keeps\n");

		return fp_help == stdout? 0 : 1;
	}

	if(n_threads < io_n_threads){
		fprintf(stderr, "[ERROR] The overall number of threads (-t [%d]) must NOT be smaller than the number of IO threads (--io-thread [%d).\n", n_threads, io_n_threads);
		return 1;
	}

	idx_rdr = ri_idx_reader_open(idx_path, &ipt, fnw);
	if (idx_rdr == 0) {
		fprintf(stderr, "[ERROR] failed to open file '%s': %s\n", idx_path, strerror(errno));
		return 1;
	}

	if (!idx_rdr->is_idx && fnw == 0 && n_reads < 1) {
		fprintf(stderr, "[ERROR] missing input: give signal files (FAST5, POD5, SLOW5/BLOW5) to detect in, or -d to save the index\n");
		ri_idx_reader_close(idx_rdr);
		return 1;
	}

	
	ri_pore_t pore;
	pore.pore_vals = NULL;
	pore.pore_inds = NULL;
	pore.max_val = -5000.0;
	pore.min_val = 5000.0;
	if(!idx_rdr->is_idx && fpore == 0){
		fprintf(stderr, "[ERROR] missing input: please specify a pore model file with -p when generating the index from a sequence file\n");
		ri_idx_reader_close(idx_rdr);
		return 1;
	}else if(!idx_rdr->is_idx){
		if (ipt.flag & RI_I_SEED_MODE_STROBE)
			fprintf(stderr, "[M::%s] seeds: strobe-events, %d strobes of %d events, the second %d-%d events after the first\n",
			        __func__, ipt.strobe_n, ipt.e, ipt.strobe_w_min, ipt.strobe_w_max);
		else
			fprintf(stderr, "[M::%s] seeds: contiguous, %d events\n", __func__, ipt.e);
		ri_pore_load(fpore, ipt.k, ipt.lev_col, &pore);
		if(!pore.pore_vals){
			fprintf(stderr, "[ERROR] cannot parse the k-mer pore model file. Please see the example k-mer model files in extern/local_kmer_models.\n");
			ri_idx_reader_close(idx_rdr);
			return 1;
		}
	}

	int n_idx = 0;
	while ((ri = ri_idx_reader_read(idx_rdr, &pore, n_threads, io_n_threads)) != 0) {
		int ret;
		++n_idx;
		if (ri_verbose >= 3)
			fprintf(stderr, "[M::%s::%.3f*%.2f] loaded/built the index for %d target sequence(s)\n",
					__func__, ri_realtime() - ri_realtime0, ri_cputime() / (ri_realtime() - ri_realtime0), ri->n_seq);
		if (n_reads > 0) ri_detopt_update(&opt, ri);
		if (ri_verbose >= 3) ri_idx_stat(ri);
		if (n_reads == 0) {
			fprintf(stderr, "[INFO] No signal files: the index is built only.\n");
			ri_idx_destroy(ri);
			continue; // no signal files: index build only
		} else {
			ret = ri_detect_files(ri, n_reads, reads, &opt, n_threads, io_n_threads);
		}
		ri_idx_destroy(ri);
		if (ret < 0) {
			fprintf(stderr, "ERROR: detection failed\n");
			exit(EXIT_FAILURE);
		}
	}
	ri_idx_reader_close(idx_rdr);
	if(pore.pore_vals)free(pore.pore_vals);
	if(pore.pore_inds)free(pore.pore_inds);
	if (n_idx == 0) {
		fprintf(stderr, "[ERROR] could not read or build an index from '%s'\n", idx_path);
		return 1;
	}

	if (fflush(stdout) == EOF) {
		perror("[ERROR] failed to write the results");
		exit(EXIT_FAILURE);
	}

	if (ri_verbose >= 3) {
		fprintf(stderr, "[M::%s] onsight %s\n", __func__, ONSIGHT_VERSION);
		// fprintf(stderr, "[M::%s] CMD:", __func__);
		// for (i = 0; i < argc; ++i) fprintf(stderr, " %s", argv[i]);
		fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, ri_realtime() - ri_realtime0, ri_cputime(), ri_peakrss() / 1024.0 / 1024.0 / 1024.0);
	}
	return 0;
}
