#include "udsp_test.h"

#include "sinc.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static float *g_table;

static int make(UDSP_SINC *s, uint32_t src, uint32_t dst) {
  size_t n = udsp_sinc_floats(src, dst, UDSP_SINC_ZEROS);
  free(g_table);
  g_table = malloc(n * sizeof *g_table);
  return g_table && udsp_sinc_init(s, src, dst, UDSP_SINC_ZEROS, g_table);
}

static float *sine(double hz, uint32_t rate, size_t n) {
  float *x = malloc(n * sizeof *x);
  for (size_t i = 0; i < n; i++) {
    x[i] = (float)(0.5 * sin(2.0 * M_PI * hz * (double)i / rate));
  }
  return x;
}

static void test_shape(void) {
  UDSP_TEST_ASSERT(udsp_sinc_floats(44100, 48000, 128) == 160u * 256,
                   "44.1 -> 48 kHz: 160 phases of 256 taps (%zu)",
                   udsp_sinc_floats(44100, 48000, 128));
  UDSP_TEST_ASSERT(udsp_sinc_floats(48000, 24000, 128) == 512,
                   "48 -> 24 kHz: one phase, the kernel twice as long");
  UDSP_TEST_ASSERT(udsp_sinc_floats(0, 48000, 128) == 0 &&
                       udsp_sinc_floats(48000, 2000000, 128) == 0 &&
                       udsp_sinc_floats(48000, 44100, 0) == 0,
                   "a rate of 0, past the limit, no zeros: refused");
  UDSP_SINC s;
  UDSP_TEST_ASSERT(make(&s, 44100, 48000) && s.up == 160 && s.down == 147 &&
                       !s.interp,
                   "the ratio reduced");
  UDSP_TEST_ASSERT(udsp_sinc_length(&s, 44100) == 48000 &&
                       udsp_sinc_length(&s, 1) == 2 &&
                       udsp_sinc_length(&s, 147) == 160 &&
                       udsp_sinc_length(&s, 148) == 162,
                   "the output's length: ceil(in * up / down)");
  UDSP_TEST_ASSERT(make(&s, 44100, 47999) && s.interp &&
                       s.rows == UDSP_SINC_PHASES_MAX + 1,
                   "47999 Hz: too many phases, interpolated rows");
}

static void test_dc(void) {
  static const uint32_t rates[][2] = {
      {44100, 48000}, {48000, 44100}, {96000, 48000}, {22050, 48000}};
  for (size_t r = 0; r < sizeof rates / sizeof *rates; r++) {
    UDSP_SINC s;
    make(&s, rates[r][0], rates[r][1]);
    float in[4096], out[512];
    for (size_t i = 0; i < 4096; i++) {
      in[i] = 0.25f;
    }
    /* far enough from both ends for every tap to land in the input */
    udsp_sinc_run(&s, in, 0, 4096, 1000, out, 512);
    double worst = 0.0;
    for (size_t i = 0; i < 512; i++) {
      double e = fabs(out[i] - 0.25);
      worst = e > worst ? e : worst;
    }
    UDSP_TEST_ASSERT(worst < 1e-6, "%u -> %u: a constant passes (off by %.2g)",
                     rates[r][0], rates[r][1], worst);
  }
}

static void test_same_rate(void) {
  UDSP_SINC s;
  make(&s, 48000, 48000);
  float in[64], out[64];
  for (int i = 0; i < 64; i++) {
    in[i] = (float)i - 20.5f;
  }
  udsp_sinc_run(&s, in, 0, 64, 0, out, 64);
  UDSP_TEST_ASSERT(!memcmp(in, out, sizeof in), "the same rate: a copy");
}

/* a tone resampled is the tone at the new rate, in the passband */
static void sine_matches(uint32_t src, uint32_t dst, double hz, double tol) {
  UDSP_SINC s;
  make(&s, src, dst);
  const size_t n = 8192;
  float *in = sine(hz, src, n);
  uint64_t out_n = udsp_sinc_length(&s, n);
  float *out = malloc(out_n * sizeof *out);
  udsp_sinc_run(&s, in, 0, n, 0, out, out_n);
  double worst = 0.0;
  for (uint64_t i = 1000; i + 1000 < out_n; i++) {
    double want = 0.5 * sin(2.0 * M_PI * hz * (double)i / dst);
    double e = fabs(out[i] - want);
    worst = e > worst ? e : worst;
  }
  UDSP_TEST_ASSERT(worst < tol, "%g Hz, %u -> %u: off by %.2g (< %.2g)", hz,
                   src, dst, worst, tol);
  free(in);
  free(out);
}

static void test_sine(void) {
  /* float's own rounding is ~3e-8 on 0.5; the kernel's ripple adds a
   * little near the band's edge */
  sine_matches(44100, 48000, 1000.0, 2e-6);
  sine_matches(44100, 48000, 15000.0, 2e-6);
  sine_matches(48000, 44100, 1000.0, 2e-6);
  sine_matches(96000, 44100, 10000.0, 2e-6);
  sine_matches(44100, 47999, 1000.0, 2e-6);
  sine_matches(44100, 47999, 15000.0, 1e-5);
}

/* past the lower rate's Nyquist: gone (120 dB) */
static void test_stopband(void) {
  UDSP_SINC s;
  make(&s, 48000, 44100);
  const size_t n = 8192;
  float *in = sine(23000.0, 48000, n);
  float out[4096];
  udsp_sinc_run(&s, in, 0, n, 1000, out, 4096);
  double peak = 0.0;
  for (size_t i = 0; i < 4096; i++) {
    peak = fabs(out[i]) > peak ? fabs(out[i]) : peak;
  }
  double db = 20.0 * log10(peak / 0.5 + 1e-30);
  UDSP_TEST_ASSERT(db < -110.0, "23 kHz into 44.1: %.1f dB", db);
  free(in);
}

/* outputs cut into blocks, each from a window of just its span, are
 * the whole run's bits -- the edges (silence outside the input) too */
static void cut(uint32_t src, uint32_t dst) {
  UDSP_SINC s;
  make(&s, src, dst);
  const size_t n = 3000;
  float *in = sine(3000.0, src, n);
  for (size_t i = 0; i < n; i += 7) {
    in[i] += 0.3f; /* not just a tone */
  }
  uint64_t out_n = udsp_sinc_length(&s, n);
  float *whole = malloc(out_n * sizeof *whole);
  float *parts = malloc(out_n * sizeof *parts);
  udsp_sinc_run(&s, in, 0, n, 0, whole, out_n);
  for (uint64_t o = 0; o < out_n;) {
    size_t len = 1 + (size_t)(o * 7 % 97);
    len = o + len > out_n ? (size_t)(out_n - o) : len;
    int64_t lo, hi;
    udsp_sinc_span(&s, o, len, &lo, &hi);
    lo = lo < 0 ? 0 : lo;
    hi = hi > (int64_t)n ? (int64_t)n : hi;
    udsp_sinc_run(&s, in + lo, lo, (size_t)(hi - lo), o, parts + o, len);
    o += len;
  }
  UDSP_TEST_ASSERT(!memcmp(whole, parts, out_n * sizeof *whole),
                   "%u -> %u: blocks from their spans are the whole run", src,
                   dst);
  free(in);
  free(whole);
  free(parts);
}

static void test_cut(void) {
  cut(44100, 48000);
  cut(48000, 44100);
  cut(44100, 47999);
}

/* an impulse's response is the row itself, scaled */
static void test_impulse(void) {
  UDSP_SINC s;
  make(&s, 44100, 48000);
  float in[1] = {1.0f}, out[400];
  udsp_sinc_run(&s, in, 500, 1, 0, out, 400);
  double sum = 0.0, peak = 0.0;
  for (size_t i = 0; i < 400; i++) {
    sum += out[i];
    peak = fabs(out[i]) > peak ? fabs(out[i]) : peak;
  }
  /* an impulse alone at frame 500 of the input: 160/147 outputs a frame,
   * each tap row summing to 1 -> the outputs sum to ~up / down */
  UDSP_TEST_ASSERT(sum == 0.0 && peak == 0.0,
                   "an impulse out of every output's reach: silence");
  udsp_sinc_run(&s, in, 200, 1, 0, out, 400);
  sum = 0.0;
  for (size_t i = 0; i < 400; i++) {
    sum += out[i];
  }
  UDSP_TEST_ASSERT(fabs(sum - 160.0 / 147.0) < 1e-4,
                   "an impulse's response sums to up / down (%.6f)", sum);
}

int main(void) {
  test_shape();
  test_dc();
  test_same_rate();
  test_sine();
  test_stopband();
  test_cut();
  test_impulse();
  free(g_table);
  return udsp_test_summary();
}
