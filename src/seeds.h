#ifndef SEEDS_H
#define SEEDS_H

/* seeds.h: quantization and seeds.
 *
 * Each event level is quantized to a bucket number of q bits (dynamic_quantize, RawHash2's adaptive quantization:
 * the levels from --fine-min to --fine-max get a share --fine-range of the buckets). Strobe-events (--seed-mode strobe, used for the results) join two strobes of
 * e quantized events each, the second starting w_min to w_max events after the first ends, chosen as the window with
 * the smallest h(strobe 1) XOR h(strobe 2) (the randstrobe idea of Sahlin 2021, on events instead of k-mers). Contiguous
 * seeds (--seed-mode contig) join e consecutive quantized events, every one or the minimizers of w.
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

/**
 * Contiguous seeds: hash values of e consecutive quantized event levels
 *
 * @param km       thread-local memory pool; using NULL falls back to malloc()
 * @param s_values normalized event values 
 * @param id       ID of the signals (e.g., ref ID or read ID); will be copied to the output $p array
 * @param strand   strand of signal values. Either 1 (forw) or 0 (reverse); will be copied to the output $p array
 * @param len      length of $s_values
 * @param w        if w>0, find a minimizer for every $w consecutive k-mers
 * @param e        number of packed events in a single 32/64-bit hash value
 * @param n        [Currently disabled] if n>0 Enables BLEND. Number of neighbors hash values to use with SimHash.
 * 				   Generates a single hash value from $n hash values; will be stored as a hash value if enabled.
 * 				   Please see the BLEND paper for details of this hash value calculation.
 * @param q        most significant Q bits of normalized event values to use in the quantization step
 * @param lq       least significant lq bits within the most significant Q bits extracted from the normalized event values.
 * 				   If you want to use n as described in the RawHash manuscript, calculate lq as follows: $lq = $Q - 2 - n
 * @param k        k-mer size that a single event represents (usually given in the k-mer model file)
 * @param p        the seeds
 *                 p->a[i].x = hash value<<RI_HASH_SHIFT | span
 *                 p->a[i].y = $id<<RI_ID_SHIFT | lastPos<<RI_POS_SHIFT | $strand
 *                 where lastPos is the position of the i-th signal,
 * 				   span indicates the region that the hash value covers ($k+$e-1),
 *                 and strand indicates whether the seed comes from the top or the bottom $strand.
 *                 Callers may want to set "p->n = 0"; otherwise results are appended to p
 */
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