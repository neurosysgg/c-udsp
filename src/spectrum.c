#include "spectrum.h"

#include "analyze.h"
#include "fft.h"

#include <math.h>

/* single-sided amplitude spectrum normalization for an N-point FFT: 2/N
 * so a full-scale sine lands near (not exactly, since this ignores the
 * Hann window's own ~0.5 coherent-gain loss -- an intentional, documented
 * simplification, see DESIGN.md) 0dBFS rather than scaling with N. */
#define UDSP_SPECTRUM_NORM (2.0f / (float)UDSP_FFT_N)

void udsp_spectrum_bars(const float *mono_tail, int sample_rate,
                        float *scratch_re, float *scratch_im,
                        float *bars_out, int n_bars) {
  udsp_window_hann(mono_tail, scratch_re, UDSP_FFT_N);
  for (int i = 0; i < UDSP_FFT_N; i++) {
    scratch_im[i] = 0.0f;
  }
  udsp_fft_forward(scratch_re, scratch_im, UDSP_FFT_N);

  size_t half = UDSP_FFT_N / 2;
  udsp_fft_magnitude(scratch_re, scratch_im, half, scratch_re);

  float nyquist = (float)sample_rate / 2.0f;
  if (nyquist <= UDSP_SPECTRUM_MIN_HZ) {
    for (int b = 0; b < n_bars; b++) {
      bars_out[b] = UDSP_DBFS_FLOOR;
    }
    return;
  }

  for (int b = 0; b < n_bars; b++) {
    float f0 = UDSP_SPECTRUM_MIN_HZ *
              powf(nyquist / UDSP_SPECTRUM_MIN_HZ, (float)b / (float)n_bars);
    float f1 = UDSP_SPECTRUM_MIN_HZ *
              powf(nyquist / UDSP_SPECTRUM_MIN_HZ, (float)(b + 1) / (float)n_bars);

    int bin0 = (int)(f0 * (float)UDSP_FFT_N / (float)sample_rate);
    int bin1 = (int)(f1 * (float)UDSP_FFT_N / (float)sample_rate);
    if (bin1 <= bin0) {
      bin1 = bin0 + 1;
    }
    if (bin0 < 0) {
      bin0 = 0;
    }
    if (bin1 > (int)half) {
      bin1 = (int)half;
    }

    float peak = 0.0f;
    for (int k = bin0; k < bin1; k++) {
      if (scratch_re[k] > peak) {
        peak = scratch_re[k];
      }
    }
    bars_out[b] = udsp_analyze_dbfs(peak * UDSP_SPECTRUM_NORM);
  }
}
