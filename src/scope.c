#include "scope.h"
#include <math.h>

void udsp_scope_extract(const float *mono, size_t n, float *trace_out) {
  if (n < 2) {
    for (size_t i = 0; i < UDSP_SCOPE_POINTS; i++) {
      trace_out[i] = 0.0f;
    }
    return;
  }

  /* rising zero-crossing nearest the tail, same triggering as before */
  size_t offset = 0;
  size_t search_end = (n > UDSP_SCOPE_POINTS) ? n - UDSP_SCOPE_POINTS : n - 1;
  size_t search_start = (search_end > UDSP_SCOPE_SEARCH_WINDOW)
                            ? search_end - UDSP_SCOPE_SEARCH_WINDOW
                            : 0;

  int found_trigger = 0;
  for (size_t i = search_start; i < search_end; i++) {
    if (mono[i] <= 0.0f && mono[i + 1] > 0.0f) {
      offset = i;
      found_trigger = 1;
      break;
    }
  }
  if (!found_trigger) {
    offset = search_start;
  }

  /* short-time autocorrelation over [min_lag, max_lag] to estimate the
   * pitch period, so the trace can be resampled to exactly one cycle
   * instead of a fixed sample count -- lag bounds are ~22Hz-2756Hz at
   * 44.1kHz, comparison window capped at 256 samples purely for
   * performance (this runs once per scope tick). */
  size_t available = n - offset;
  size_t min_lag = 16;
  size_t max_lag = 2048;
  if (max_lag > available / 2) {
    max_lag = available / 2;
  }

  size_t period_samples = UDSP_SCOPE_POINTS; /* fallback: no period found */
  if (max_lag > min_lag) {
    float max_corr = -1.0f;
    size_t best_lag = 0;

    for (size_t lag = min_lag; lag <= max_lag; lag++) {
      size_t cmp_len = available - lag;
      if (cmp_len > 256) {
        cmp_len = 256;
      }

      float corr = 0.0f;
      for (size_t k = 0; k < cmp_len; k++) {
        corr += mono[offset + k] * mono[offset + k + lag];
      }

      if (corr > max_corr) {
        max_corr = corr;
        best_lag = lag;
      }
    }

    if (best_lag > 0 && max_corr > 0.001f) {
      period_samples = best_lag;
    }
  }

  /* resample exactly one period into UDSP_SCOPE_POINTS via linear
   * interpolation, so the displayed waveform holds still instead of
   * scrolling. */
  for (size_t i = 0; i < UDSP_SCOPE_POINTS; i++) {
    float pos = (float)offset +
                (float)i * (float)period_samples / (float)UDSP_SCOPE_POINTS;
    size_t index0 = (size_t)pos;
    size_t index1 = index0 + 1;
    float frac = pos - (float)index0;

    if (index1 < n) {
      trace_out[i] = mono[index0] * (1.0f - frac) + mono[index1] * frac;
    } else if (index0 < n) {
      trace_out[i] = mono[index0];
    } else {
      trace_out[i] = 0.0f;
    }
  }
}
