#include "udsp_test.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "analysis/truepeak.h"

#include <math.h>
#include <stdlib.h>

#define LEN 4800

/* a sine of amp at rate/4 with its peaks between the samples: every sample
 * sits at amp * 0.707, 3 dB under the wave */
static void between(float *x, float amp) {
  for (int i = 0; i < LEN; i++) {
    const float v = amp * (float)sin(2 * M_PI * i / 4.0 + M_PI / 4);
    x[2 * i] = x[2 * i + 1] = v;
  }
}

static void test_between_samples(void) {
  float *x = malloc(2 * LEN * sizeof *x);
  between(x, 1.0f);
  UDSP_TRUEPEAK t;
  udsp_truepeak_init(&t, 2);
  udsp_truepeak_add(&t, x, LEN);
  /* BS.1770's tolerance is -0.55 .. +0.55 dB round the true value */
  UDSP_TEST_ASSERT(20 * log10(t.peak) > -0.55 && 20 * log10(t.peak) < 0.55,
                   "20 * log10(t.peak) > -0.55 && 20 * log10(t.peak) < 0.55");
  free(x);
}

static void test_a_low_sine_reads_its_level(void) {
  float *x = malloc(2 * LEN * sizeof *x);
  for (int i = 0; i < LEN; i++) {
    x[2 * i] = x[2 * i + 1] = 0.5f * (float)sin(2 * M_PI * 997 * i / 48000.0);
  }
  UDSP_TRUEPEAK t;
  udsp_truepeak_init(&t, 2);
  udsp_truepeak_add(&t, x, LEN);
  UDSP_TEST_ASSERT(fabs(t.peak - 0.5) < 0.01, "fabs(t.peak - 0.5) < 0.01");
  free(x);
}

static void test_any_block_size(void) {
  float *x = malloc(2 * LEN * sizeof *x);
  between(x, 0.8f);
  UDSP_TRUEPEAK a, b;
  udsp_truepeak_init(&a, 2);
  udsp_truepeak_init(&b, 2);
  udsp_truepeak_add(&a, x, LEN);
  for (size_t at = 0; at < LEN;) {
    const size_t n = at % 7 + 1;
    const size_t k = LEN - at < n ? LEN - at : n;
    udsp_truepeak_add(&b, x + 2 * at, k);
    at += k;
  }
  UDSP_TEST_ASSERT(a.peak == b.peak, "a.peak == b.peak");
  free(x);
}

int main(void) {
  test_between_samples();
  test_a_low_sine_reads_its_level();
  test_any_block_size();
  return udsp_test_summary();
}
