#include "udsp_test.h"

#include <math.h>
#include <string.h>

#include "fx/xover.h"

#define RATE 48000
#define N RATE

static float in[2 * N], band[5][2 * N];

/* a sine at hz, both channels */
static void sine(float hz) {
  for (int i = 0; i < N; i++) {
    in[2 * i] = in[2 * i + 1] = (float)sin(6.283185307179586 * hz * i / RATE);
  }
}

/* RMS amplitude of x's left channel over the second half (settled) */
static float amp(const float *x) {
  double sq = 0.0;
  for (int i = N / 2; i < N; i++) {
    sq += (double)x[2 * i] * x[2 * i];
  }
  return (float)sqrt(2.0 * sq / (N / 2));
}

static float db(float g) { return 20.0f * log10f(g); }

/* hz through a crossover of bands at splits: each band's level in
 * dB[b], the bands' sum's */
static float split(int bands, const float *splits, float hz, float *dbs) {
  UDSP_XOVER x;
  udsp_xover_reset(&x);
  udsp_xover_set(&x, bands, splits, RATE);
  sine(hz);
  float *out[5] = {band[0], band[1], band[2], band[3], band[4]};
  udsp_xover_split(&x, in, out, N);
  static float sum[2 * N];
  memset(sum, 0, sizeof sum);
  for (int b = 0; b < bands; b++) {
    dbs[b] = db(amp(band[b]));
    for (int i = 0; i < 2 * N; i++) {
      sum[i] += band[b][i];
    }
  }
  return db(amp(sum));
}

static void test_flat(void) {
  static const float SPLITS[] = {120.0f, 2500.0f, 6000.0f, 12000.0f};
  static const float HZ[] = {30, 120, 400, 2500, 4000, 6000, 9000, 15000};
  float worst = 0.0f, dbs[5];
  for (int bands = 2; bands <= 5; bands++) {
    for (size_t i = 0; i < sizeof HZ / sizeof *HZ; i++) {
      const float s = split(bands, SPLITS, HZ[i], dbs);
      worst = fabsf(s) > worst ? fabsf(s) : worst;
    }
  }
  UDSP_TEST_ASSERT(worst < 0.01f,
                   "2 .. 5 bands summed: flat within %.4f dB (< 0.01)",
                   (double)worst);
}

static void test_bands(void) {
  static const float SPLITS[] = {120.0f, 2500.0f};
  float dbs[3];
  split(3, SPLITS, 40.0f, dbs);
  /* a Linkwitz-Riley low side: 1 / (1 + w^4), w = 40 / 120 */
  const float lr4 = -20.0f * log10f(1.0f + 1.0f / 81.0f);
  UDSP_TEST_ASSERT(fabsf(dbs[0] - lr4) < 0.005f && dbs[1] < -35.0f &&
                       dbs[2] < -60.0f,
                   "40 Hz: the low band, %.3f dB as LR4's %.3f (%.1f, %.1f dB)",
                   (double)dbs[0], (double)lr4, (double)dbs[1], (double)dbs[2]);
  split(3, SPLITS, 600.0f, dbs);
  UDSP_TEST_ASSERT(dbs[1] > -0.3f && dbs[0] < -45.0f && dbs[2] < -35.0f,
                   "600 Hz: the middle (%.1f, %.2f, %.1f dB)", (double)dbs[0],
                   (double)dbs[1], (double)dbs[2]);
  split(3, SPLITS, 120.0f, dbs);
  UDSP_TEST_ASSERT(fabsf(dbs[0] + 6.02f) < 0.1f && fabsf(dbs[1] + 6.02f) < 0.2f,
                   "at a split: both sides -6 dB (%.2f, %.2f)", (double)dbs[0],
                   (double)dbs[1]);
  split(3, SPLITS, 10000.0f, dbs);
  UDSP_TEST_ASSERT(dbs[2] > -0.1f && dbs[1] < -35.0f,
                   "10 kHz: the high band (%.1f, %.2f)", (double)dbs[1],
                   (double)dbs[2]);
}

/* one band: the input as it is; in may be the last output */
static void test_one_and_alias(void) {
  UDSP_XOVER x;
  udsp_xover_reset(&x);
  udsp_xover_set(&x, 1, NULL, RATE);
  sine(1000.0f);
  memcpy(band[0], in, sizeof in);
  float *out[1] = {band[0]};
  udsp_xover_split(&x, band[0], out, N);
  UDSP_TEST_ASSERT(!memcmp(band[0], in, sizeof in), "one band: the input");
  static const float S[] = {500.0f};
  float *o2[2] = {band[1], band[0]};
  udsp_xover_reset(&x);
  udsp_xover_set(&x, 2, S, RATE);
  udsp_xover_split(&x, band[0], o2, N);
  float worst = 0.0f;
  UDSP_XOVER y;
  udsp_xover_reset(&y);
  udsp_xover_set(&y, 2, S, RATE);
  float *o3[2] = {band[2], band[3]};
  udsp_xover_split(&y, in, o3, N);
  for (int i = 0; i < 2 * N; i++) {
    float d = fabsf(band[0][i] - band[3][i]) + fabsf(band[1][i] - band[2][i]);
    worst = d > worst ? d : worst;
  }
  UDSP_TEST_ASSERT(worst == 0.0f, "in as the last output: the same bits");
}

int main(void) {
  printf("xover_test\n");
  test_flat();
  test_bands();
  test_one_and_alias();
  return udsp_test_summary();
}
