#ifndef SEEDS_H
#define SEEDS_H

/* seeds.h: quantization and seeds.
 *
 * Each event level is quantized to a bucket number of q bits (dynamic_quantize, RawHash2's adaptive quantization:
 * the levels from --fine-min to --fine-max get a share --fine-range of the buckets). Strobe-events (--seed-mode
 * strobe, the default) join two strobes of e quantized events each, the second starting w_min to w_max events after
 * the first ends, chosen as the window with the smallest h(strobe 1) XOR h(strobe 2) (the randstrobe idea of Sahlin
 * 2021, on events instead of k-mers). Contiguous seeds (--seed-mode contig) join e consecutive quantized events,
 * every one or the minimizers of w.
 *
 * In: a sequence of event levels (a read's, or a signature's reference signal). Out: its seeds (hash, position,
 * strand), each at the position of its first event. State: none.
 */

#include <string.h> //for memset
#include "utils.h"

// a seed's fields in mm128_t: x = hash << RI_HASH_SHIFT | span; y = id << RI_ID_SHIFT | pos << RI_POS_SHIFT | strand
#define RI_HASH_SHIFT 6
#define RI_ID_SHIFT 32
#define RI_POS_SHIFT 1

#ifdef __cplusplus
extern "C" {
#endif

uint32_t dynamic_quantize(float signal,
						  float fine_min,
						  float fine_max,
						  float fine_range,
						  uint32_t n_buckets);

/* Contiguous seeds: the hash of e consecutive quantized event levels, at every position (w = 0) or the minimizers
 * of each w consecutive hashes, appended to p in the layout above. km: the thread's memory pool (NULL: malloc);
 * n: unused (RawHash2's BLEND); quant_bit: q. */
void ri_seeds_contig(void *km,
               const float* s_values,
               uint32_t id,
               int strand,
               uint32_t len,
               int diff,
               int w,
               int e,
               int n,
               uint32_t quant_bit,
               int k,
               float fine_min,
               float fine_max,
               float fine_range,
               mm128_v *p);

/**
 * Strobe-events: two strobes, each a window of e_s quantized events; strobe 2 is the window starting in
 * [strobe 1's end + w_min, strobe 1's end + w_max] with the smallest h(strobe 1) XOR h(window).
 *
 * The seed's position (.y) is that of strobe 1's first event, as for contiguous seeds; its hash mixes both strobes'
 * hashes over 64 bits. With diff >= 0, skipped events do not count toward the strobes, as for contiguous seeds.
 *
 * @param e_s        strobe length k_s in events; e_s*quant_bit must be <= 32
 * @param w_min      minimum gap between strobe-1's last event and strobe-2's first event (events)
 * @param w_max      maximum gap; ring buffer scratch sized as 2*e_s + w_max + 1
 * @param n_strobes  strobes per seed; only 2
 */
void ri_seeds_strobe(void *km,
                      const float* s_values,
                      uint32_t id,
                      int strand,
                      uint32_t len,
                      int diff,
                      int e_s,
                      int w_min,
                      int w_max,
                      int n_strobes,
                      uint32_t quant_bit,
                      int k,
                      float fine_min,
                      float fine_max,
                      float fine_range,
                      mm128_v *p);

#ifdef __cplusplus
}
#endif
#endif //SEEDS_H