#ifndef UDSP_SCOPE_H
#define UDSP_SCOPE_H

#include <stddef.h>

/* samples in an extracted trace -- a caller maps this onto whatever
 * display width it actually has at render time, a fixed-resolution-
 * source/variable-resolution-display convention shared with
 * UDSP_SPECTROGRAM_COLS (spectrogram.h). */
#define UDSP_SCOPE_POINTS 256

/* how far back from "the most recent possible trace" udsp_scope_extract()
 * will search for a trigger -- bounds the search cost regardless of how
 * large the buffer it's searching within is. */
#define UDSP_SCOPE_SEARCH_WINDOW 2048

/* extracts a UDSP_SCOPE_POINTS-sample trace from mono (n samples),
 * starting at the earliest rising zero-crossing (mono[i] <= 0 followed by
 * mono[i+1] > 0) found within the last UDSP_SCOPE_SEARCH_WINDOW samples
 * before "the most recent possible trace" -- a real oscilloscope's
 * trigger, without which the trace would visibly swim/jitter every call
 * since the buffer's start has no relationship to the signal's own
 * phase. Falls back to the earliest sample in that search window,
 * untriggered, if no crossing exists in range (e.g. near-silence, pure
 * DC, or n too small yet).
 *
 * From the trigger point, a short-time autocorrelation estimates the
 * signal's pitch period, and exactly one period is linearly resampled to
 * fill trace_out -- so a periodic signal's trace holds still (one stable
 * cycle) rather than scrolling with phase drift each call. Falls back to
 * a UDSP_SCOPE_POINTS-sample period (i.e. the raw post-trigger samples)
 * if no clear periodicity is found. Reads that run off the end of mono
 * (n too short past the trigger) zero-pad the remainder of trace_out. */
void udsp_scope_extract(const float *mono, size_t n, float *trace_out);

#endif
