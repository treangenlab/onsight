/* sig_summary.h: each read's signatures ranked by score (--sig-summary) and the per-taxon table (--sig-profile).
 *
 * A signature's score is its hits over its emitted positions in the index, each hit weighted by
 * 1 / log2(positions + 1) with --sig-weight-common-hashes. A read's signatures are ranked by score, then hash
 * containment, fanout, distinct hashes, hits, span and rid. Its call is its best signature unless that has fewer hits
 * or distinct hashes than --sig-min-hits and --sig-min-distinct-hashes, or is within --sig-min-ratio of the second.
 *
 * In: a read's seed hits per signature. Out: its --sig-summary row; at the end of the run, the --sig-profile table.
 * State: the profile's counts per taxon (ri_profile_t), updated in read order.
 */
#ifndef SIG_SUMMARY_H
#define SIG_SUMMARY_H

#include <stdint.h>
#include <stdio.h>
#include "index.h"
#include "options.h"
#include "seed_hits.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	uint32_t   *ord;          /* the read's signatures, best first */
	const char *call_taxon;   /* its call, for the profile ("*": none) */
	float       best_score;
	char       *row;          /* its --sig-summary row; NULL without --sig-summary */
} ri_sig_summary_t;

/* For a read with hits; the caller frees out->ord and out->row */
void ri_sig_summary_read(const ri_idx_t *ri, const ri_detopt_t *opt, const ri_read_t *rd, int want_row,
                         ri_sig_summary_t *out);

/* The row of a read without hits */
char *ri_sig_summary_empty_row(const char *read_name);

typedef struct ri_profile_s ri_profile_t;

ri_profile_t *ri_profile_init(const ri_idx_t *ri);
void ri_profile_add_read(ri_profile_t *p, const ri_idx_t *ri, const ri_detopt_t *opt, const ri_read_t *rd,
                         const uint32_t *ord, const char *call_name, float best_score);
void ri_profile_write(const ri_profile_t *p, FILE *fp);
void ri_profile_destroy(ri_profile_t *p);

#ifdef __cplusplus
}
#endif
#endif // SIG_SUMMARY_H
