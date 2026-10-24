#include "udsp_test.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "base/fft.h"
#include "offline/irfb.h"

/* a response: a few taps that ring (a decaying, alternating comb) */
static void response(float *h, size_t len, size_t stride) {
  for (size_t i = 0; i < len; i++) {
    h[i * stride] = (float)(pow(0.8, (double)i) * (i % 3 == 0 ? 1 : -0.5));
  }
}

/* the loop simulated: in[t] = x[t] + g * y[t - d], y = h * in, for an
 * impulse x; frames of y */
static void direct(const float *h, size_t len, float g, size_t d, double *y,
                   size_t frames) {
  double *in = calloc(frames, sizeof(double));
  for (size_t t = 0; t < frames; t++) {
    in[t] = (t == 0) + (t >= d ? g * y[t - d] : 0.0);
    double s = 0.0;
    for (size_t k = 0; k < len && k <= t; k++) {
      s += h[k] * in[t - k];
    }
    y[t] = s;
  }
  free(in);
}

static double worst(const float *a, const double *b, size_t n) {
  double w = 0.0;
  for (size_t i = 0; i < n; i++) {
    const double e = fabs(a[i] - b[i]);
    w = e > w ? e : w;
  }
  return w;
}

/* h (len) fed back by g through d: frames of it, over the transform
 * udsp_irfb_size picks */
static float *fed(const float *h, size_t len, size_t stride, float g, size_t d,
                  size_t frames) {
  const size_t n = udsp_irfb_size(frames);
  float *table = malloc(UDSP_FFT_TABLE_FLOATS(n) * sizeof(float));
  float *scratch = malloc(UDSP_IRFB_SCRATCH(n) * sizeof(float));
  float *y = malloc(frames * sizeof(float));
  udsp_fft_table(table, n);
  udsp_irfb(h, len, stride, g, d, y, frames, 1, n, table, scratch);
  free(table);
  free(scratch);
  return y;
}

/* the loop's gain r at its peak frequency (g = r / peak |H|, as a
 * caller keeps it settling), simulated for the frames it asks for */
static void test_against_the_loop(void) {
  enum { LEN = 24, N = 64 };
  float h[LEN], table[UDSP_FFT_TABLE_FLOATS(N)], scratch[UDSP_IRFB_SCRATCH(N)];
  response(h, LEN, 1);
  udsp_fft_table(table, N);
  const float peak = udsp_irfb_peak(h, LEN, 1, N, table, scratch);
  static const struct {
    float r;
    size_t d;
  } CASES[] = {{0.5f, 1}, {-0.8f, 7}, {0.95f, 40}};
  for (size_t c = 0; c < sizeof CASES / sizeof *CASES; c++) {
    const float g = CASES[c].r / peak;
    const size_t frames =
        udsp_irfb_frames(LEN, CASES[c].d, CASES[c].r, 1 << 16);
    double *ref = malloc(frames * sizeof(double));
    direct(h, LEN, g, CASES[c].d, ref, frames);
    float *y = fed(h, LEN, 1, g, CASES[c].d, frames);
    const double w = worst(y, ref, frames);
    UDSP_TEST_ASSERT(w < 2e-5,
                     "%.2f at the peak through %zu frames: the loop's "
                     "response over its %zu frames (worst %.2g)",
                     (double)CASES[c].r, CASES[c].d, frames, w);
    free(y);
    free(ref);
  }
}

static void test_no_feedback(void) {
  enum { LEN = 50 };
  float h[2 * LEN];
  response(h, LEN, 2);
  for (size_t i = 0; i < LEN; i++) {
    h[2 * i + 1] = 99.0f; /* the other channel: not read */
  }
  float *y = fed(h, LEN, 2, 0.0f, 3, LEN + 10);
  double ref[LEN + 10] = {0};
  for (size_t i = 0; i < LEN; i++) {
    ref[i] = h[2 * i];
  }
  UDSP_TEST_ASSERT(worst(y, ref, LEN + 10) < 1e-6,
                   "no feedback: the response itself (one channel of two)");
  free(y);
}

static void test_peak(void) {
  enum { N = 64 };
  float table[UDSP_FFT_TABLE_FLOATS(N)], scratch[UDSP_IRFB_SCRATCH(N)];
  udsp_fft_table(table, N);
  const float two[2] = {1.0f, 1.0f}, dirac[1] = {0.5f};
  UDSP_TEST_ASSERT(
      fabsf(udsp_irfb_peak(two, 2, 1, N, table, scratch) - 2.0f) < 1e-6f &&
          fabsf(udsp_irfb_peak(dirac, 1, 1, N, table, scratch) - 0.5f) < 1e-6f,
      "the peak |H|: two equal taps 2 (at DC), a single one its size");
}

static void test_frames(void) {
  const size_t cap = 48000 * 30;
  /* 0.5^15 is -90.3 dB: 15 trips of 110 frames after the response's 100 */
  UDSP_TEST_ASSERT(udsp_irfb_frames(100, 10, 0.5f, cap) == 100 + 15 * 110 &&
                       udsp_irfb_frames(100, 10, -0.5f, cap) == 100 + 15 * 110,
                   "-90 dB after 15 trips at 0.5 (either sign)");
  UDSP_TEST_ASSERT(udsp_irfb_frames(100, 10, 0.0f, cap) == 100 &&
                       udsp_irfb_frames(100, 10, 1.0f, cap) == cap &&
                       udsp_irfb_frames(100, 10, 0.9999f, cap) == cap &&
                       udsp_irfb_frames(cap + 5, 0, 0.0f, cap) == cap,
                   "none: the response; 1 or near it: the cap; never past it");
  UDSP_TEST_ASSERT(udsp_irfb_size(1) == 4 && udsp_irfb_size(1000) == 2048 &&
                       udsp_irfb_size(1024) == 2048,
                   "the transform: a power of two, twice the frames");
}

/* the loop's delay sets the resonance: a single tap fed back through d
 * rings at rate / d -- its spectrum peaks at bin n / d */
static void test_delay_tunes(void) {
  enum { FRAMES = 4096 };
  const float one = 1.0f;
  for (size_t d = 8; d <= 32; d *= 2) {
    float *y = fed(&one, 1, 1, 0.9f, d, FRAMES);
    int comb = 1;
    for (size_t i = 0; i < FRAMES && i < 6 * d; i++) {
      const size_t trip = i / d;
      const double want = i % d == 0 ? pow(0.9, (double)trip) : 0.0;
      comb &= fabs(y[i] - want) < 1e-5;
    }
    UDSP_TEST_ASSERT(comb, "a tap through %zu frames: a comb, 0.9 a trip", d);
    free(y);
  }
}

int main(void) {
  test_against_the_loop();
  test_no_feedback();
  test_peak();
  test_frames();
  test_delay_tunes();
  return udsp_test_summary();
}
