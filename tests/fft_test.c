#include "udsp_test.h"

#include "fft.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static void test_hann_window(void) {
  float in[8], out[8];
  for (int i = 0; i < 8; i++) {
    in[i] = 1.0f;
  }
  udsp_window_hann(in, out, 8);
  UDSP_TEST_ASSERT(out[0] < 0.01f, "Hann window is ~0 at the first sample");
  UDSP_TEST_ASSERT(out[7] < 0.01f, "Hann window is ~0 at the last sample");
  UDSP_TEST_ASSERT(out[3] > 0.9f || out[4] > 0.9f,
                   "Hann window peaks near the center");
}

static void test_hann_window_in_place(void) {
  float buf[4] = {1.0f, 1.0f, 1.0f, 1.0f};
  udsp_window_hann(buf, buf, 4);
  UDSP_TEST_ASSERT(buf[0] < 0.01f, "Hann window works in-place (in == out)");
}

static void test_fft_dc_signal(void) {
  size_t n = 8;
  float re[8], im[8];
  for (size_t i = 0; i < n; i++) {
    re[i] = 1.0f;
    im[i] = 0.0f;
  }
  udsp_fft_forward(re, im, n);
  float mag[8];
  udsp_fft_magnitude(re, im, n, mag);
  UDSP_TEST_ASSERT(mag[0] > 7.9f,
                   "a constant (DC) signal's entire energy lands in bin 0");
  for (size_t i = 1; i < n; i++) {
    UDSP_TEST_ASSERT(mag[i] < 0.01f,
                     "every non-DC bin is ~0 for a pure DC signal (bin %zu)",
                     i);
  }
}

static void test_fft_nyquist_signal(void) {
  /* alternating +1/-1 is the highest frequency an 8-sample-per-cycle
   * discrete signal can represent -- all its energy belongs at bin n/2. */
  size_t n = 8;
  float re[8], im[8];
  for (size_t i = 0; i < n; i++) {
    re[i] = (i % 2 == 0) ? 1.0f : -1.0f;
    im[i] = 0.0f;
  }
  udsp_fft_forward(re, im, n);
  float mag[8];
  udsp_fft_magnitude(re, im, n, mag);
  UDSP_TEST_ASSERT(
      mag[4] > 7.9f,
      "an alternating +-1 signal's energy lands at the Nyquist bin (n/2)");
}

/* --- the table-driven transforms, against a plain O(n^2) DFT in double --- */

static unsigned int g_seed = 7;
static float noise(void) {
  g_seed = g_seed * 1103515245u + 12345u;
  return (float)((g_seed >> 8) & 0xffff) / 32768.0f - 1.0f;
}

static void dft(const float *re, const float *im, size_t n, double *out_re,
                double *out_im) {
  for (size_t k = 0; k < n; k++) {
    double sr = 0, si = 0;
    for (size_t t = 0; t < n; t++) {
      double a = -2.0 * M_PI * (double)((k * t) % n) / (double)n;
      sr += re[t] * cos(a) - im[t] * sin(a);
      si += re[t] * sin(a) + im[t] * cos(a);
    }
    out_re[k] = sr;
    out_im[k] = si;
  }
}

/* largest error over bins 0..bins-1, relative to the largest magnitude */
static double rel_err(const float *re, const float *im, const double *xr,
                      const double *xi, size_t bins) {
  double err = 0, top = 1e-30;
  for (size_t k = 0; k < bins; k++) {
    double e = hypot(re[k] - xr[k], im[k] - xi[k]);
    double mag = hypot(xr[k], xi[k]);
    err = e > err ? e : err;
    top = mag > top ? mag : top;
  }
  return err / top;
}

static void test_fft_table_matches_dft(void) {
  size_t cap = 4096;
  float *in_re = malloc(cap * sizeof(float)),
        *in_im = malloc(cap * sizeof(float));
  float *re = malloc(cap * sizeof(float)), *im = malloc(cap * sizeof(float));
  float *table = malloc(UDSP_FFT_TABLE_FLOATS(cap) * sizeof(float));
  double *xr = malloc(cap * sizeof(double)), *xi = malloc(cap * sizeof(double));
  for (size_t n = 2; n <= cap; n <<= 1) {
    for (size_t i = 0; i < n; i++) {
      in_re[i] = noise();
      in_im[i] = noise();
    }
    dft(in_re, in_im, n, xr, xi);
    udsp_fft_table(table, n);
    memcpy(re, in_re, n * sizeof(float));
    memcpy(im, in_im, n * sizeof(float));
    udsp_fft_forward_table(re, im, n, table);
    double e_table = rel_err(re, im, xr, xi, n);
    memcpy(re, in_re, n * sizeof(float));
    memcpy(im, in_im, n * sizeof(float));
    udsp_fft_forward(re, im, n);
    double e_old = rel_err(re, im, xr, xi, n);
    UDSP_TEST_ASSERT(e_table < 1e-5,
                     "table FFT n=%zu matches the DFT (rel err %.2g, "
                     "the recurrence FFT: %.2g)",
                     n, e_table, e_old);
  }
  free(in_re), free(in_im), free(re), free(im), free(table), free(xr), free(xi);
}

static void test_fft_real_matches_dft(void) {
  size_t cap = 4096;
  float *in = malloc(cap * sizeof(float)), *zero = calloc(cap, sizeof(float));
  float *re = malloc((cap / 2 + 1) * sizeof(float));
  float *im = malloc((cap / 2 + 1) * sizeof(float));
  float *table = malloc(UDSP_FFT_TABLE_FLOATS(cap) * sizeof(float));
  double *xr = malloc(cap * sizeof(double)), *xi = malloc(cap * sizeof(double));
  for (size_t n = 4; n <= cap; n <<= 1) {
    for (size_t i = 0; i < n; i++) {
      in[i] = noise();
    }
    dft(in, zero, n, xr, xi);
    udsp_fft_table(table, n);
    udsp_fft_real(in, re, im, n, table);
    double e = rel_err(re, im, xr, xi, n / 2 + 1);
    UDSP_TEST_ASSERT(
        e < 1e-5, "real FFT n=%zu matches the DFT's bins 0..n/2 (rel err %.2g)",
        n, e);
  }
  free(in), free(zero), free(re), free(im), free(table), free(xr), free(xi);
}

/* past 2^16 points the bit reversal takes its slow path */
static void test_fft_large(void) {
  size_t n = 1u << 18;
  float *in = malloc(n * sizeof(float));
  float *re = malloc(n * sizeof(float)), *im = malloc(n * sizeof(float));
  float *rre = malloc((n / 2 + 1) * sizeof(float)),
        *rim = malloc((n / 2 + 1) * sizeof(float));
  float *table = malloc(UDSP_FFT_TABLE_FLOATS(n) * sizeof(float));
  for (size_t i = 0; i < n; i++) {
    in[i] = noise();
    re[i] = in[i];
    im[i] = 0.0f;
  }
  udsp_fft_table(table, n);
  udsp_fft_forward_table(re, im, n, table);
  udsp_fft_real(in, rre, rim, n, table);
  double err = 0, top = 0;
  for (size_t k = 0; k <= n / 2; k++) {
    double e = hypot(rre[k] - re[k], rim[k] - im[k]);
    double m = hypot(re[k], im[k]);
    err = e > err ? e : err;
    top = m > top ? m : top;
  }
  UDSP_TEST_ASSERT(
      err / top < 1e-5,
      "2^18 points: real and complex table FFTs agree (rel err %.2g)",
      err / top);
  /* a tone at bin 1000 lands there */
  for (size_t i = 0; i < n; i++) {
    in[i] = (float)cos(2.0 * M_PI * 1000.0 * (double)i / (double)n);
  }
  udsp_fft_real(in, rre, rim, n, table);
  size_t best = 0;
  for (size_t k = 1; k <= n / 2; k++) {
    if (hypot(rre[k], rim[k]) > hypot(rre[best], rim[best])) {
      best = k;
    }
  }
  UDSP_TEST_ASSERT(best == 1000,
                   "2^18 points: a tone at bin 1000 peaks there (got %zu)",
                   best);
  free(in), free(re), free(im), free(rre), free(rim), free(table);
}

static void test_power_and_window_table(void) {
  float re[9] = {3, 0, 1, -2, 5, 0.5f, 0, 1, 2},
        im[9] = {4, 1, 0, 2, 0, 0.5f, 0, -1, 2};
  float p[9];
  udsp_fft_power(re, im, 9, p);
  UDSP_TEST_ASSERT(
      p[0] == 25.0f && p[1] == 1.0f && p[3] == 8.0f && p[8] == 8.0f,
      "power = re^2 + im^2, the tail past the vector width included");
  float in[37], w[37], a[37], b[37];
  for (int i = 0; i < 37; i++) {
    in[i] = noise();
  }
  udsp_window_hann_table(w, 37);
  udsp_window_apply(in, w, a, 37);
  udsp_window_hann(in, b, 37);
  int same = 1;
  for (int i = 0; i < 37; i++) {
    same &= a[i] == b[i];
  }
  UDSP_TEST_ASSERT(
      same, "a Hann table applied gives exactly udsp_window_hann()'s samples");
  udsp_window_apply(in, w, in, 37);
  UDSP_TEST_ASSERT(!memcmp(in, a, sizeof(a)),
                   "udsp_window_apply() works in place");

  /* the fused window: the same bins as windowing first, every path's size */
  for (size_t n = 4; n <= 4096; n <<= 1) {
    static float sig[4096], win[4096], wd[4096], tab[8192];
    static float r1[2049], i1[2049], r2[2049], i2[2049];
    for (size_t i = 0; i < n; i++) {
      sig[i] = noise();
    }
    udsp_window_hann_table(win, n);
    udsp_fft_table(tab, n);
    udsp_window_apply(sig, win, wd, n);
    udsp_fft_real(wd, r1, i1, n, tab);
    udsp_fft_real_windowed(sig, win, r2, i2, n, tab);
    UDSP_TEST_ASSERT(
        !memcmp(r1, r2, (n / 2 + 1) * sizeof(float)) &&
            !memcmp(i1, i2, (n / 2 + 1) * sizeof(float)),
        "n=%zu: the fused window gives windowing first's bins exactly", n);
  }
}

int main(void) {
  test_fft_table_matches_dft();
  test_fft_real_matches_dft();
  test_fft_large();
  test_power_and_window_table();
  test_hann_window();
  test_hann_window_in_place();
  test_fft_dc_signal();
  test_fft_nyquist_signal();
  return udsp_test_summary();
}
