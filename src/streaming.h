/* streaming.h: the stream layer, each read's scores in the output step.
 *
 * Per read: per signature, its diagonal band counts (diagonal_bands.h); per taxon, its best signature; the read's
 * p-value and call (calls.h); the near-neighbor check of a call at the loosest cutoff (1e-3, or the loosest
 * --alpha-read); the slope refits (slope_range.h). This part depends only on the read and its minibatch's slope
 * grid, so the reads of a minibatch are prepared at once, on all threads (ri_stream_prepare_read). The rest runs
 * once per read in read order (ri_stream_commit_read): the slope ring, the read record (--stream-emit-fd), and the
 * read and target calls. The grid changes only in ri_stream_end_batch, after the minibatch's last commit.
 *
 * In: per read, its seed hits per signature (and its event levels for the near-neighbor check). Out: the read
 * records and slopes lines (--stream-emit-fd, off by default); the read and target calls through calls.h. State: the
 * slope grid and ring, the count of reads committed.
 */
#ifndef STREAMING_H
#define STREAMING_H

#include <stdint.h>
#include <stdio.h>

#include "index.h"
#include "seed_hits.h"
#include "diagonal_bands.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Settings, filled in once from the options before the first read. */
typedef struct ri_stream_config_s {
	/* fd of the read records (JSON lines); -1: none */
	int      emit_fd;
	FILE    *emit_fp;            /* fdopen()'d from emit_fd */

	/* The slopes of the diagonal bands; a run starts on the search grid (--stream-diag-slopes). A read whose events run
	 * faster or slower than the signature's positions (real reads: about 1.05-1.1) spreads one true chain over several
	 * slope-1 bands. */
	ri_slope_grid_t grid;

	/* > 0: the read records list, per taxon, the non-overlapping hit counts of its signatures: hits at least this many
	 * events apart in one band (12 = one seed's two strobes; --stream-nonoverlap-gap). Needed by the read p-value. */
	uint32_t nonoverlap_gap;

	/* The sliding-window slope range (slope_range.h): on when slope_chains > 0, which needs nonoverlap_gap > 0 and
	 * calls. */
	uint32_t    slope_chains;
	double      slope_alpha;
	float       slope_width;
	int         slope_width_sd;
	float       slope_cover;

	/* Read calls and target calls (calls.h; --read-calls, --target-calls); NULL: off. Owned by the caller. */
	struct ri_calls_s *calls;
} ri_stream_config_t;

typedef struct ri_stream_state_s ri_stream_state_t;
typedef struct ri_read_scores_s ri_read_scores_t;        /* one read's scores, from prepare to commit */
typedef struct ri_stream_scratch_s ri_stream_scratch_t;  /* one thread's buffers for ri_stream_prepare_read */

ri_stream_state_t *ri_stream_init(const ri_stream_config_t *cfg, const ri_idx_t *ri);
void               ri_stream_destroy(ri_stream_state_t *s);

/* One scratch per thread that prepares reads. Destroy them before ri_calls_finish: each adds its near-neighbor
 * check count and time to the calls. */
ri_stream_scratch_t *ri_stream_scratch_init(const ri_stream_state_t *s);
void                 ri_stream_scratch_destroy(ri_stream_state_t *s, ri_stream_scratch_t *x);

/* Any thread: the scores of a read. rd->levels must be set for the near-neighbor check. */
ri_read_scores_t *ri_stream_prepare_read(const ri_stream_state_t *s, const ri_idx_t *ri, const ri_read_t *rd,
                                         ri_stream_scratch_t *x);

/* Once per read in read order: the slope ring, the read record, the read and target calls. Frees sc. */
void ri_stream_commit_read(ri_stream_state_t *s, const ri_read_t *rd, ri_read_scores_t *sc, uint64_t read_id);

/* Once after each minibatch's last commit: renews the slope grid (--stream-slope-chains). */
void ri_stream_end_batch(ri_stream_state_t *s);

/* Defaults, the same for R9.4.1 and R10.4.1. */
void ri_stream_config_init(ri_stream_config_t *cfg);

#ifdef __cplusplus
}
#endif
#endif // STREAMING_H
