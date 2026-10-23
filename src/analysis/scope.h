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

/* how many samples each candidate lag's correlation is measured over --
 * caps the period-estimation cost at O(max_lag * this) per call, which
 * is the dominant cost of udsp_scope_extract(). Long enough to cover
 * several cycles of anything in the lag range, short enough to keep one
 * call comfortably inside a visualizer's frame budget. */
#define UDSP_SCOPE_CORR_WINDOW 256

/* minimum normalized correlation for a lag to count as "this signal is
 * genuinely periodic here" -- a plain correlation coefficient in
 * [-1, 1], so this is a level-independent threshold, unlike a raw
 * correlation sum (which scales with amplitude and would make the same
 * number mean different things at different volumes). Below this for
 * every candidate lag, udsp_scope_extract() falls back to a
 * UDSP_SCOPE_POINTS-sample period rather than locking onto noise. */
#define UDSP_SCOPE_CORR_MIN 0.3f

/* how close to the best score a local peak must come to be accepted as
 * the fundamental -- see udsp_scope_extract()'s own comment on why the
 * first qualifying *peak* is the right choice rather than the global
 * maximum (octave error) or the first threshold crossing (short-lag
 * shoulder). */
#define UDSP_SCOPE_PEAK_FRAC 0.90f

/* period-search bounds, in samples: ~2756Hz down to ~22Hz at 44.1kHz --
 * comfortably outside the fundamental range of anything a music
 * visualizer needs to lock onto at either end. MAX also sizes
 * udsp_scope_extract()'s on-stack score table, so it stays a compile-time
 * constant (this library is allocation-free by construction). */
#define UDSP_SCOPE_MIN_LAG 16
#define UDSP_SCOPE_MAX_LAG 2048

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
