#include "udsp_test.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "conv.h"
#include "fft.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static unsigned g_seed = 1;

static float rnd(void) {
  g_seed = g_seed * 1664525u + 1013904223u;
  return (float)(g_seed >> 8) / 8388608.0f - 1.0f;
}

static void test_inverse_round_trip(void) {
  static const size_t NS[] = {4, 8, 16, 64, 1024, 8192};
  for (size_t t = 0; t < sizeof NS / sizeof *NS; t++) {
    const size_t n = NS[t];
    float *x = malloc(n * sizeof(float)), *y = malloc(n * sizeof(float));
    float *re = malloc((n / 2 + 1) * sizeof(float));
    float *im = malloc((n / 2 + 1) * sizeof(float));
    float *table = malloc(UDSP_FFT_TABLE_FLOATS(n) * sizeof(float));
    udsp_fft_table(table, n);
    for (size_t i = 0; i < n; i++) {
      x[i] = rnd();
    }
    udsp_fft_real(x, re, im, n, table);
    udsp_fft_real_inverse(re, im, y, n, table);
    float worst = 0.0f;
    for (size_t i = 0; i < n; i++) {
      float e = fabsf(y[i] / (float)n - x[i]);
      worst = e > worst ? e : worst;
    }
    UDSP_TEST_ASSERT(worst < 2e-6f,
                     "n=%zu: inverse(forward(x)) = n * x (worst %g)", n,
                     (double)worst);
    free(x), free(y), free(re), free(im), free(table);
  }
}

/* a spectrum built by hand (one cosine, one sine) comes back as them */
static void test_inverse_known(void) {
  enum { N = 32 };
  float re[N / 2 + 1] = {0}, im[N / 2 + 1] = {0}, out[N], table[2 * N];
  udsp_fft_table(table, N);
  /* x[i] = 1 + cos(2 pi 3 i / N) - 2 sin(2 pi 5 i / N): X[0] = N,
   * X[3] = N / 2, X[5] = i N (sin's bin: -i N/2 per unit, times -2) */
  re[0] = N;
  re[3] = N / 2;
  im[5] = N;
  udsp_fft_real_inverse(re, im, out, N, table);
  float worst = 0.0f;
  for (int i = 0; i < N; i++) {
    float want = 1.0f + (float)cos(2 * M_PI * 3 * i / N) -
                 2.0f * (float)sin(2 * M_PI * 5 * i / N);
    float e = fabsf(out[i] / N - want);
    worst = e > worst ? e : worst;
  }
  UDSP_TEST_ASSERT(worst < 1e-5f,
                   "a cosine and a sine from their bins (worst %g)",
                   (double)worst);
}

static void test_parts(void) {
  UDSP_TEST_ASSERT(udsp_conv_parts(0, 64) == 1,
                   "an empty response is one partition");
  UDSP_TEST_ASSERT(udsp_conv_parts(64, 64) == 1,
                   "64 frames in blocks of 64: one");
  UDSP_TEST_ASSERT(udsp_conv_parts(65, 64) == 2, "65 frames: two");
}

/* x (frames at stride) convolved with h directly, in double */
static double direct(const float *x, size_t stride, size_t t, const float *h,
                     size_t hn) {
  double s = 0.0;
  for (size_t j = 0; j < hn && j <= t; j++) {
    s += (double)h[j] * x[(t - j) * stride];
  }
  return s;
}

/* a response of hn frames over blocks of b: every block's output against
 * the direct sum */
static float against_direct(size_t hn, size_t b, size_t n) {
  float *h = malloc(hn * sizeof(float)), *x = malloc(n * sizeof(float));
  float *y = malloc(n * sizeof(float));
  for (size_t i = 0; i < hn; i++) {
    h[i] = rnd() * expf(-3.0f * (float)i / (float)hn);
  }
  for (size_t i = 0; i < n; i++) {
    x[i] = rnd();
  }
  const size_t parts = udsp_conv_parts(hn, b), s = UDSP_CONV_SPECTRUM(b);
  float *spectra = malloc(parts * s * sizeof(float));
  float *line = malloc(parts * s * sizeof(float)),
        *acc = malloc(s * sizeof(float));
  float *frame = malloc(2 * b * sizeof(float)),
        *scratch = malloc(2 * b * sizeof(float));
  float *table = malloc(UDSP_FFT_TABLE_FLOATS(2 * b) * sizeof(float));
  udsp_fft_table(table, 2 * b);
  udsp_conv_ir(h, hn, 1, b, spectra, table, scratch);
  UDSP_CONV_IN in;
  udsp_conv_in_init(&in, b, parts, line, frame);
  for (size_t at = 0; at + b <= n; at += b) {
    udsp_conv_in_push(&in, x + at, 1, table);
    memset(acc, 0, s * sizeof(float));
    udsp_conv_mac(&in, spectra, parts, acc);
    udsp_conv_out(acc, b, y + at, 1, table, scratch);
  }
  double worst = 0.0;
  for (size_t t = 0; t < n / b * b; t++) {
    double e = fabs(y[t] - direct(x, 1, t, h, hn));
    worst = e > worst ? e : worst;
  }
  free(h), free(x), free(y), free(spectra), free(line), free(acc), free(frame);
  free(scratch), free(table);
  return (float)worst;
}

static void test_against_direct(void) {
  float e = against_direct(1000, 64, 8192);
  UDSP_TEST_ASSERT(e < 1e-4f,
                   "1000 frames in blocks of 64 = the direct sum (worst %g)",
                   (double)e);
  e = against_direct(64, 64, 4096);
  UDSP_TEST_ASSERT(e < 1e-4f, "one whole partition (worst %g)", (double)e);
  e = against_direct(5, 256, 4096);
  UDSP_TEST_ASSERT(e < 1e-5f, "a response shorter than a block (worst %g)",
                   (double)e);
  e = against_direct(20000, 512, 48000);
  UDSP_TEST_ASSERT(e < 1e-3f, "20000 frames in blocks of 512 (worst %g)",
                   (double)e);
}

/* a unit impulse at frame d is a delay of d: the input back, d later */
static void test_delay(void) {
  enum { B = 32, D = 45, N = 1024 };
  float h[D + 1] = {0}, x[N], y[N];
  h[D] = 1.0f;
  for (int i = 0; i < N; i++) {
    x[i] = rnd();
  }
  const size_t parts = udsp_conv_parts(D + 1, B), s = UDSP_CONV_SPECTRUM(B);
  float spectra[2 * UDSP_CONV_SPECTRUM(B)], line[2 * UDSP_CONV_SPECTRUM(B)];
  float acc[UDSP_CONV_SPECTRUM(B)], frame[2 * B], scratch[2 * B], table[4 * B];
  udsp_fft_table(table, 2 * B);
  udsp_conv_ir(h, D + 1, 1, B, spectra, table, scratch);
  UDSP_CONV_IN in;
  udsp_conv_in_init(&in, B, parts, line, frame);
  for (int at = 0; at < N; at += B) {
    udsp_conv_in_push(&in, x + at, 1, table);
    memset(acc, 0, s * sizeof(float));
    udsp_conv_mac(&in, spectra, parts, acc);
    udsp_conv_out(acc, B, y + at, 1, table, scratch);
  }
  float worst = 0.0f, lead = 0.0f;
  for (int i = 0; i < N; i++) {
    float e = fabsf(y[i] - (i >= D ? x[i - D] : 0.0f));
    worst = e > worst ? e : worst;
    lead = i < D && fabsf(y[i]) > lead ? fabsf(y[i]) : lead;
  }
  UDSP_TEST_ASSERT(worst < 1e-6f, "an impulse at 45 delays by 45 (worst %g)",
                   (double)worst);
  UDSP_TEST_ASSERT(lead < 1e-6f, "nothing before it (%g)", (double)lead);
}

/* true stereo: two interleaved inputs' lines into one output, each
 * through its own response -- the sum of the two direct convolutions */
static void test_two_inputs(void) {
  enum { B = 64, H = 300, N = 2048 };
  static float hl[H], hr[H], x[2 * N], y[N];
  for (int i = 0; i < H; i++) {
    hl[i] = rnd() * 0.5f;
    hr[i] = rnd() * 0.5f;
  }
  for (int i = 0; i < 2 * N; i++) {
    x[i] = rnd();
  }
  const size_t parts = udsp_conv_parts(H, B), s = UDSP_CONV_SPECTRUM(B);
  static float sl[5 * UDSP_CONV_SPECTRUM(B)], sr[5 * UDSP_CONV_SPECTRUM(B)];
  static float ll[5 * UDSP_CONV_SPECTRUM(B)], lr[5 * UDSP_CONV_SPECTRUM(B)];
  static float acc[UDSP_CONV_SPECTRUM(B)], fl[2 * B], fr[2 * B], scratch[2 * B];
  static float table[4 * B];
  udsp_fft_table(table, 2 * B);
  udsp_conv_ir(hl, H, 1, B, sl, table, scratch);
  udsp_conv_ir(hr, H, 1, B, sr, table, scratch);
  UDSP_CONV_IN l, r;
  udsp_conv_in_init(&l, B, parts, ll, fl);
  udsp_conv_in_init(&r, B, parts, lr, fr);
  for (int at = 0; at < N; at += B) {
    udsp_conv_in_push(&l, x + 2 * at, 2, table);
    udsp_conv_in_push(&r, x + 2 * at + 1, 2, table);
    memset(acc, 0, s * sizeof(float));
    udsp_conv_mac(&l, sl, parts, acc);
    udsp_conv_mac(&r, sr, parts, acc);
    udsp_conv_out(acc, B, y + at, 1, table, scratch);
  }
  double worst = 0.0;
  for (int t = 0; t < N; t++) {
    double want =
        direct(x, 2, (size_t)t, hl, H) + direct(x + 1, 2, (size_t)t, hr, H);
    double e = fabs(y[t] - want);
    worst = e > worst ? e : worst;
  }
  UDSP_TEST_ASSERT(
      worst < 1e-4,
      "left and right through their own responses, summed (worst %g)", worst);
}

int main(void) {
  printf("conv_test\n");
  test_inverse_round_trip();
  test_inverse_known();
  test_parts();
  test_against_direct();
  test_delay();
  test_two_inputs();
  return udsp_test_summary();
}
