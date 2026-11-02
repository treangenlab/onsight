/* pipeline.h: the read pipeline (detection).
 *
 * Reads come from the signal files in minibatches of -K samples, through three steps that overlap (kt_pipeline):
 * read a minibatch; per read on all threads, its events, seeds and seed hits per signature; then the output step. The
 * output step runs per read on all threads what depends only on the read and the minibatch's slope grid (its
 * --sig-summary row; its stream scores: diagonal bands, read p-value, near-neighbor check, slope refits), then once
 * per read in read order the rows, the profile and the stream state (read records, read calls, target calls, the
 * slope ring). The slope grid changes only between minibatches.
 *
 * In: the index, the signal files, the options. Out: --sig-summary, --sig-profile, --read-calls and --target-calls,
 * and the read records on fd 3 (--stream-emit-fd).
 */
#ifndef PIPELINE_H
#define PIPELINE_H

#include "index.h"
#include "options.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Returns 0, or -1 when a file cannot be opened or the index has no signature metadata. */
int ri_detect_files(const ri_idx_t *idx, int n_segs, const char **fn, const ri_detopt_t *opt, int n_threads, int io_n_threads);

#ifdef __cplusplus
}
#endif
#endif // PIPELINE_H
