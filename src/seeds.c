/* seeds.c: quantization and seeds; see seeds.h. */
#include "seeds.h"
#include <assert.h>
#include "rh_kvec.h"
#include <math.h>
#include <float.h>

static inline uint64_t hash64(uint64_t key, uint64_t mask){
    key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
    key = key ^ key >> 24;
    key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
    key = key ^ key >> 14;
    key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
    key = key ^ key >> 28;
    key = (key + (key << 31)) & mask;
    return key;
}

uint32_t dynamic_quantize(float signal,
						  float fine_min,
						  float fine_max,
						  float fine_range,
						  uint32_t n_buckets)
{
    // Total range for normalization
    float minVal = -3.0, maxVal = 3.0;
    float range = maxVal - minVal;
	float coarse_coef1 = (1-fine_range)/2;
	float coarse_coef2 = fine_range + coarse_coef1;

    float normalized = (signal - minVal) / range;

	float a = (fine_min - minVal) / range;
	float b = (fine_max - minVal) / range;

    float quantized = fine_max;
    if (signal >= fine_min && signal <= fine_max) {
        quantized = fine_range * ((normalized - a) / (b - a));
    }
	else {
        if (normalized < 0.5) quantized = fine_range + coarse_coef1 * normalized;
        else quantized = coarse_coef2 + coarse_coef1 * normalized;
    }

    uint32_t quantizedValue = (uint32_t)(quantized * (n_buckets-1));

    return quantizedValue;
}

void seeds_minimizer(void *km,
				   const float* s_values,
				   uint32_t id,
				   int strand,
				   uint32_t len,
				   int diff,
				   int w,
				   int e,
				   uint32_t quant_bit,
				   int k,
				   float fine_min,
				   float fine_max,
				   float fine_range,
				   mm128_v *p)
{
	assert(len > 0 && (w > 0 && w < 256) && e*quant_bit <= (64-RI_HASH_SHIFT));
	
	int j, buf_pos, min_pos;
	mm128_t buf[256], min = { UINT64_MAX, UINT64_MAX };

	uint32_t span = k+e-1;
	
	const uint64_t id_shift = (uint64_t)id<<RI_ID_SHIFT, mask = (1ULL<<32)-1, mask_events = (1ULL<<(quant_bit*e))-1, mask_quant_bit = (1ULL<<quant_bit)-1;

	uint32_t n_buckets = 1UL<<quant_bit;

	memset(buf, 0xff, w * 16);
	rh_kv_resize(mm128_t, km, *p, p->n + len/w);

	int sigBufFull = 0;
	uint32_t f_pos, l, sigBufPos = 0;
	uint32_t f_tmpQuantSignal = 0;
	uint64_t quantVal = 0;

	mm128_t sigBuf[e];
	memset(sigBuf, 0, e*sizeof(mm128_t));

	uint32_t prev_quant_for_diff = UINT32_MAX;

    for (f_pos = l = buf_pos = min_pos = 0; f_pos < len; ++f_pos) {
		f_tmpQuantSignal = dynamic_quantize(s_values[f_pos], fine_min, fine_max, fine_range, n_buckets)&mask_quant_bit;

		/* Quantized sig-diff: skip if |quant_diff| <= diff. diff=-1 disables filtering. */
		if(f_pos > 0 && diff >= 0) {
			int qdiff = (int)f_tmpQuantSignal - (int)prev_quant_for_diff;
			if(qdiff < 0) qdiff = -qdiff;
			if(qdiff <= diff) continue;
		}

		l++;
		mm128_t info = { UINT64_MAX, UINT64_MAX };
		prev_quant_for_diff = f_tmpQuantSignal;

		quantVal = (quantVal<<quant_bit|f_tmpQuantSignal)&mask_events;

		sigBuf[sigBufPos].y = id_shift | (uint32_t)f_pos<<RI_POS_SHIFT | strand;
		if(++sigBufPos == e) {sigBufFull = 1; sigBufPos = 0;}
		sigBuf[sigBufPos].x = hash64(quantVal, mask)<<RI_HASH_SHIFT | span;

		if(!sigBufFull) continue;

		info.x = sigBuf[sigBufPos].x;
		info.y = sigBuf[sigBufPos].y;

		buf[buf_pos] = info; // need to do this here as appropriate buf_pos and buf[buf_pos] are needed below
		if (l == w + e - 1 && min.x != UINT64_MAX) { // special case for the first window - because identical k-mers are not stored yet
			for (j = buf_pos + 1; j < w; ++j)
				if (min.x == buf[j].x && buf[j].y != min.y) rh_kv_push(mm128_t, km, *p, buf[j]);
			for (j = 0; j < buf_pos; ++j)
				if (min.x == buf[j].x && buf[j].y != min.y) rh_kv_push(mm128_t, km, *p, buf[j]);
		}
		if (info.x <= min.x) { // a new minimum; then write the old min
			if (l >= w + e && min.x != UINT64_MAX) rh_kv_push(mm128_t, km, *p, min);
			min = info, min_pos = buf_pos;
		} else if (buf_pos == min_pos) { // old min has moved outside the window
			if (l >= w + e - 1 && min.x != UINT64_MAX) rh_kv_push(mm128_t, km, *p, min);
			for (j = buf_pos + 1, min.x = UINT64_MAX; j < w; ++j) // the two loops are necessary when there are identical k-mers
				if (min.x >= buf[j].x) min = buf[j], min_pos = j; // >= is important s.t. min is always the closest e-mer
			for (j = 0; j <= buf_pos; ++j)
				if (min.x >= buf[j].x) min = buf[j], min_pos = j;
			if (l >= w + e - 1 && min.x != UINT64_MAX) { // write identical k-mers
				for (j = buf_pos + 1; j < w; ++j) // these two loops make sure the output is sorted
					if (min.x == buf[j].x && min.y != buf[j].y) rh_kv_push(mm128_t, km, *p, buf[j]);
				for (j = 0; j <= buf_pos; ++j)
					if (min.x == buf[j].x && min.y != buf[j].y) rh_kv_push(mm128_t, km, *p, buf[j]);
			}
		}
		if (++buf_pos == w) buf_pos = 0;
    }
	if (min.x != UINT64_MAX)
		rh_kv_push(mm128_t, km, *p, min);
}

void seeds_every(void *km,
				   const float* s_values,
				   uint32_t id,
				   int strand,
				   uint32_t len,
				   int diff,
				   int e,
				   uint32_t quant_bit,
				   int k,
				   float fine_min,
				   float fine_max,
				   float fine_range,
				   mm128_v *p){

	assert(len > 0 && (uint32_t)e*quant_bit <= 64);

	uint32_t span = k+e-1;
	
	const uint64_t id_shift = (uint64_t)id<<RI_ID_SHIFT, mask = (1ULL<<32)-1, mask_events = (1ULL<<(quant_bit*e))-1, mask_quant_bit = (1ULL<<quant_bit)-1;

	uint32_t n_buckets = 1UL<<quant_bit;

	int sigBufFull = 0;
	uint32_t f_pos = 0, sigBufPos = 0;
	uint32_t f_tmpQuantSignal = 0;
	uint64_t quantVal = 0;

	rh_kv_resize(mm128_t, km, *p, p->n + len-span-1);

	mm128_t sigBuf[e];
	memset(sigBuf, 0, e*sizeof(mm128_t));

	//First quantization is done here
	f_tmpQuantSignal = dynamic_quantize(s_values[f_pos], fine_min, fine_max, fine_range, n_buckets)&mask_quant_bit;
	uint32_t prev_quant_for_diff = f_tmpQuantSignal;
	sigBuf[sigBufPos].y = id_shift | (uint32_t)f_pos<<RI_POS_SHIFT | strand;
	if(++sigBufPos == e) {sigBufFull = 1; sigBufPos = 0;}
	quantVal = f_tmpQuantSignal&mask_events;
	sigBuf[sigBufPos].x = (hash64(quantVal, mask)<<RI_HASH_SHIFT) | span;
	if(sigBufFull) rh_kv_push(mm128_t, km, *p, sigBuf[sigBufPos]);

    for (f_pos = 1; f_pos < len; ++f_pos) {
		f_tmpQuantSignal = dynamic_quantize(s_values[f_pos], fine_min, fine_max, fine_range, n_buckets)&mask_quant_bit;

		/* Quantized sig-diff: skip if |quant_diff| <= diff. diff=-1 disables filtering. */
		if(diff >= 0) {
			int qdiff = (int)f_tmpQuantSignal - (int)prev_quant_for_diff;
			if(qdiff < 0) qdiff = -qdiff;
			if(qdiff <= diff) continue;
		}

		prev_quant_for_diff = f_tmpQuantSignal;

		sigBuf[sigBufPos].y = id_shift | (uint32_t)f_pos<<RI_POS_SHIFT | strand;
		if(++sigBufPos == e) {sigBufFull = 1; sigBufPos = 0;}

		quantVal = (quantVal<<quant_bit|f_tmpQuantSignal)&mask_events;
		sigBuf[sigBufPos].x = (hash64(quantVal, mask)<<RI_HASH_SHIFT) | span;

		if(!sigBufFull) continue;

		rh_kv_push(mm128_t, km, *p, sigBuf[sigBufPos]);
    }
}

/* On each kept event the hash of the strobe ending there goes into a ring; once 2*e_s + w_max - 1 events are in it,
 * the seed whose strobe 1 is the oldest window is emitted, at strobe 1's first event. */
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
                      mm128_v *p)
{
	assert(len > 0);
	assert(n_strobes == 2);
	assert(e_s > 0 && (uint32_t)e_s * quant_bit <= 32);
	assert(w_min > 0 && w_max >= w_min);

	const int R = 2 * e_s + w_max + 1; /* ring depth covering strobe-1-first..strobe-2-last */
	const uint64_t id_shift = (uint64_t)id << RI_ID_SHIFT;
	const uint64_t mask = (1ULL << 32) - 1;
	const uint64_t mask_strobe = (1ULL << (quant_bit * e_s)) - 1;
	const uint64_t mask_quant_bit = (1ULL << quant_bit) - 1;
	const uint64_t span_mask = (1ULL << RI_HASH_SHIFT) - 1;
	const uint32_t n_buckets = 1UL << quant_bit;

	uint32_t span = (uint32_t)(k + 2 * e_s + w_max - 2);
	if (span > span_mask) span = (uint32_t)span_mask; /* clamp; chain DP only uses .y for ordering */

	rh_kv_resize(mm128_t, km, *p, p->n + (size_t)len);

	uint32_t h_ring[R];
	uint32_t pos_ring[R];
	uint8_t  valid_ring[R];
	memset(h_ring, 0, sizeof(h_ring));
	memset(pos_ring, 0, sizeof(pos_ring));
	memset(valid_ring, 0, sizeof(valid_ring));

	int cursor = 0;
	uint32_t accepted = 0;
	uint64_t quantWin = 0;
	uint32_t prev_quant_for_diff = 0;

	const int delta_s1F = 2 * e_s + w_max - 2; /* ring delta to strobe-1 first event */
	const int delta_s1L = e_s + w_max - 1;     /* ring delta to strobe-1 last event */
	const int max_cand_delta = w_max - w_min;  /* strobe-2 last-event delta range [0, max_cand_delta] */

	for (uint32_t f_pos = 0; f_pos < len; ++f_pos) {
		uint32_t q = dynamic_quantize(s_values[f_pos], fine_min, fine_max, fine_range, n_buckets) & mask_quant_bit;

		if (accepted > 0 && diff >= 0) {
			int qdiff = (int)q - (int)prev_quant_for_diff;
			if (qdiff < 0) qdiff = -qdiff;
			if (qdiff <= diff) continue;
		}
		prev_quant_for_diff = q;
		++accepted;

		quantWin = ((quantWin << quant_bit) | q) & mask_strobe;

		uint8_t valid_now = (accepted >= (uint32_t)e_s) ? 1 : 0;
		uint32_t h_now = valid_now ? (uint32_t)hash64(quantWin, mask) : 0u;

		h_ring[cursor] = h_now;
		pos_ring[cursor] = f_pos;
		valid_ring[cursor] = valid_now;

		if (accepted >= (uint32_t)(2 * e_s + w_max - 1)) {
			int i_s1F = (cursor - delta_s1F + R) % R;
			int i_s1L = (cursor - delta_s1L + R) % R;
			if (valid_ring[i_s1L]) {
				uint32_t h1 = h_ring[i_s1L];
				uint32_t pos_s1F = pos_ring[i_s1F];
				uint32_t best_score = UINT32_MAX;
				uint32_t best_h2 = 0;
				int found = 0;
				int d;
				for (d = 0; d <= max_cand_delta; ++d) {
					int i_cand = (cursor - d + R) % R;
					if (!valid_ring[i_cand]) continue;
					uint32_t score = h1 ^ h_ring[i_cand];
					if (!found || score < best_score) {
						best_score = score;
						best_h2 = h_ring[i_cand];
						found = 1;
					}
				}
				if (found) {
					/* The two strobe hashes in one 64-bit value, mixed by hash64 over all 64 bits, then folded to 32
					 * bits. With a 32-bit mask the high half (strobe 2) would be lost. */
					uint64_t comp = ((uint64_t)h1 << 32) | (uint64_t)best_h2;
					uint64_t mixed = hash64(comp, ~(uint64_t)0);
					uint32_t composite = (uint32_t)(mixed ^ (mixed >> 32)) & (uint32_t)mask;
					mm128_t info;
					info.x = ((uint64_t)composite << RI_HASH_SHIFT) | (uint64_t)(span & span_mask);
					info.y = id_shift | ((uint64_t)pos_s1F << RI_POS_SHIFT) | (uint64_t)(strand & 1);
					rh_kv_push(mm128_t, km, *p, info);
				}
			}
		}

		cursor = (cursor + 1) % R;
	}
}

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
               mm128_v *p)
{
	if(w) seeds_minimizer(km, s_values, id, strand, len, diff, w, e, quant_bit, k, fine_min, fine_max, fine_range, p);
	else seeds_every(km, s_values, id, strand, len, diff, e, quant_bit, k, fine_min, fine_max, fine_range, p);
}
