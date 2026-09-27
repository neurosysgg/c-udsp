#include "spectrum.h"

#include "analyze.h"
#include "fft.h"

#include <math.h>

/* single-sided amplitude spectrum normalization for an N-point FFT: 2/N
 * so a full-scale sine lands near (not exactly, since this ignores the
 * Hann window's own ~0.5 coherent-gain loss -- an intentional, documented
 * simplification, see DESIGN.md) 0dBFS rather than scaling with N. */
#define UDSP_SPECTRUM_NORM (2.0f / (float)UDSP_FFT_N)

void udsp_spectrum_magnitude(const float *mono_tail, float *scratch_re,
                             float *scratch_im, float *mag_out) {
  udsp_window_hann(mono_tail, scratch_re, UDSP_FFT_N);
  for (int i = 0; i < UDSP_FFT_N; i++) {
    scratch_im[i] = 0.0f;
  }
  udsp_fft_forward(scratch_re, scratch_im, UDSP_FFT_N);
  udsp_fft_magnitude(scratch_re, scratch_im, UDSP_FFT_N / 2, mag_out);
}

void udsp_spectrum_bin(const float *mag, int sample_rate, float *bars_out, int n_bars) {
  size_t half = UDSP_FFT_N / 2;
  float nyquist = (float)sample_rate / 2.0f;
  if (nyquist <= UDSP_SPECTRUM_MIN_HZ) {
    for (int b = 0; b < n_bars; b++) {
      bars_out[b] = UDSP_DBFS_FLOOR;
    }
    return;
  }

  /* f0/f1 walk the same log-spaced edges as before, but via one powf()
   * per call (the per-bar step ratio) plus a running multiply instead of
   * two powf() calls per bar -- transcendental calls dominate this loop's
   * cost far more than the multiply they replace. */
  float ratio = powf(nyquist / UDSP_SPECTRUM_MIN_HZ, 1.0f / (float)n_bars);
  float f0 = UDSP_SPECTRUM_MIN_HZ;
  for (int b = 0; b < n_bars; b++) {
    float f1 = f0 * ratio;

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
      if (mag[k] > peak) {
        peak = mag[k];
      }
    }
    bars_out[b] = udsp_analyze_dbfs(peak * UDSP_SPECTRUM_NORM);

    f0 = f1;
  }
}

void udsp_spectrum_bars(const float *mono_tail, int sample_rate,
                        float *scratch_re, float *scratch_im,
                        float *bars_out, int n_bars) {
  udsp_spectrum_magnitude(mono_tail, scratch_re, scratch_im, scratch_re);
  udsp_spectrum_bin(scratch_re, sample_rate, bars_out, n_bars);
}

/* --- the planned path --- */

/* bar edges exactly as udsp_spectrum_bin() walks them, at n points */
static void plan_bands(UDSP_SPECTRUM_PLAN *plan) {
  int n = plan->n, half = n / 2;
  float nyquist = (float)plan->sample_rate / 2.0f;
  if (nyquist <= UDSP_SPECTRUM_MIN_HZ) {
    for (int b = 0; b < plan->n_bars; b++) {
      plan->bin0[b] = plan->bin1[b] = 0; /* an empty band: the floor */
    }
    return;
  }
  float ratio = powf(nyquist / UDSP_SPECTRUM_MIN_HZ, 1.0f / (float)plan->n_bars);
  float f0 = UDSP_SPECTRUM_MIN_HZ;
  for (int b = 0; b < plan->n_bars; b++) {
    float f1 = f0 * ratio;
    int bin0 = (int)(f0 * (float)n / (float)plan->sample_rate);
    int bin1 = (int)(f1 * (float)n / (float)plan->sample_rate);
    if (bin1 <= bin0) {
      bin1 = bin0 + 1;
    }
    if (bin0 < 0) {
      bin0 = 0;
    }
    if (bin1 > half) {
      bin1 = half;
    }
    plan->bin0[b] = bin0;
    plan->bin1[b] = bin1;
    f0 = f1;
  }
}

void udsp_spectrum_plan(UDSP_SPECTRUM_PLAN *plan, int n, int sample_rate,
                        int n_bars) {
  if (n < 16 || n > UDSP_SPECTRUM_MAX_N || (n & (n - 1))) {
    n = UDSP_FFT_N;
  }
  n_bars = n_bars < 1 ? 1 : n_bars;
  n_bars = n_bars > UDSP_SPECTRUM_MAX_BARS ? UDSP_SPECTRUM_MAX_BARS : n_bars;
  plan->n = n;
  plan->sample_rate = sample_rate;
  plan->n_bars = n_bars;
  plan_bands(plan);
  /* the normalisation (2/n) rides in the window: scaled input, scaled
   * power, no pass over the bins for it */
  udsp_window_hann_table(plan->window, (size_t)n);
  for (int i = 0; i < n; i++) {
    plan->window[i] *= 2.0f / (float)n;
  }
  udsp_fft_table(plan->fft, (size_t)n);
}

float *udsp_spectrum_power(const UDSP_SPECTRUM_PLAN *plan,
                           const float *mono_tail, float *scratch) {
  size_t n = (size_t)plan->n;
  float *re = scratch, *im = scratch + n / 2 + 16;
  udsp_fft_real_windowed(mono_tail, plan->window, re, im, n, plan->fft);
  udsp_fft_power(re, im, n / 2 + 1, re);
  return re;
}

void udsp_spectrum_bin_planned(const UDSP_SPECTRUM_PLAN *plan,
                               const float *power, float *bars_out) {
  for (int b = 0; b < plan->n_bars; b++) {
    /* the top bands span hundreds of bins: the vector max (power is never
     * negative, so its peak is its max), where it pays for its setup */
    int k0 = plan->bin0[b], k1 = plan->bin1[b];
    float peak = 0.0f;
    if (k1 - k0 >= 32) {
      peak = udsp_analyze_peak(power + k0, (size_t)(k1 - k0));
    } else {
      for (int k = k0; k < k1; k++) {
        peak = power[k] > peak ? power[k] : peak;
      }
    }
    /* the loudest bin's amplitude: one sqrt per bar, not per bin */
    bars_out[b] = udsp_analyze_dbfs(sqrtf(peak));
  }
}

void udsp_spectrum_bars_planned(const UDSP_SPECTRUM_PLAN *plan,
                                const float *mono_tail, float *scratch,
                                float *bars_out) {
  udsp_spectrum_bin_planned(plan, udsp_spectrum_power(plan, mono_tail, scratch),
                            bars_out);
}
