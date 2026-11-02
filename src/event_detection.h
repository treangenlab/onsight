#ifndef EVENT_DETECTION_H
#define EVENT_DETECTION_H

/* event_detection.h: the raw signal split into events (RawHash2's event detection, after Sigmap).
 *
 * The signal is normalized over the read; two t-tests over sliding windows (--seg-window-length1, -length2) find
 * where the level changes (--seg-threshold1, -threshold2, --seg-peak-height; --r10 sets the R10.4.1 values). An event
 * is the stretch between two changes, its level the mean of its samples; events shorter than --min-seg-length or
 * longer than --max-seg-length samples are skipped.
 *
 * In: one read's raw signal. Out: its event levels (mean current of each event). State: none.
 */

#include <stdint.h>
#include "options.h"

#ifdef __cplusplus
extern "C" {
#endif

float* normalize_signal(void *km,
						const float* sig,
						const uint32_t s_len,
						double* mean_sum,
						double* std_dev_sum,
						uint32_t* n_events_sum,
						uint32_t* n_sig);
/**
 * Detects events from signals
 *
 * @param km	thread-local memory pool; using NULL falls back to malloc()
 * @param s_len	length of $sig
 * @param sig	signal values
 * @param n		number of events
 * 
 * @return		list of event values of length $n
 */
float* detect_events(void *km,
					 const uint32_t s_len,
					 const float* sig,
					 const uint32_t window_length1,
					 const uint32_t window_length2,
					 const float threshold1,
					 const float threshold2,
					 const float peak_height,
					 const uint32_t min_seg_len,
					 const uint32_t max_seg_len,
					 double* mean_sum,
					 double* std_dev_sum,
					 uint32_t* n_events_sum,
					 uint32_t *n_events);

#ifdef __cplusplus
}
#endif
#endif //EVENT_DETECTION_H