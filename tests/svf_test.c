#include "udsp_test.h"

#include <math.h>

#include "svf.h"

static void dc(float *b, int frames) {
  for (int i = 0; i < frames * 2; i++) {
    b[i] = 1.0f;
  }
}

static void test_dc_passes_low_and_not_high(void) {
  static float b[4000 * 2];
  UDSP_SVF s;
  udsp_svf_reset(&s);
  udsp_svf_set(&s, 1000.0f, 0.3f, UDSP_SVF_LOW, 48000);
  dc(b, 4000);
  udsp_svf_stereo(&s, b, 4000);
  UDSP_TEST_ASSERT(fabsf(b[7998] - 1.0f) < 1e-4f &&
                       fabsf(b[7999] - 1.0f) < 1e-4f,
                   "low-pass: DC settles at 1 (%g)", (double)b[7998]);
  udsp_svf_reset(&s);
  udsp_svf_set(&s, 1000.0f, 0.3f, UDSP_SVF_HIGH, 48000);
  dc(b, 4000);
  udsp_svf_stereo(&s, b, 4000);
  UDSP_TEST_ASSERT(fabsf(b[7998]) < 1e-4f, "high-pass: DC settles at 0 (%g)",
                   (double)b[7998]);
}

static void test_split_blocks_identical(void) {
  static float a[300 * 2], b[300 * 2];
  for (int i = 0; i < 600; i++) {
    a[i] = b[i] = (float)((i * 7919) % 200) / 100.0f - 1.0f;
  }
  UDSP_SVF x, y;
  udsp_svf_reset(&x);
  udsp_svf_set(&x, 2500.0f, 0.8f, UDSP_SVF_BAND, 48000);
  y = x;
  udsp_svf_stereo(&x, a, 300);
  udsp_svf_stereo(&y, b, 77);
  udsp_svf_stereo(&y, b + 154, 223);
  int same = 1;
  for (int i = 0; i < 600; i++) {
    same &= a[i] == b[i];
  }
  UDSP_TEST_ASSERT(same, "a block split anywhere filters the same");
}

static void test_cutoff_clamped_stable(void) {
  static float b[2000 * 2];
  UDSP_SVF s;
  udsp_svf_reset(&s);
  udsp_svf_set(&s, 900000.0f, 0.99f, UDSP_SVF_LOW, 48000);
  dc(b, 2000);
  udsp_svf_stereo(&s, b, 2000);
  int bounded = 1;
  for (int i = 0; i < 4000; i++) {
    bounded &= fabsf(b[i]) <= 4.0f;
  }
  UDSP_TEST_ASSERT(bounded && s.f < 2.0f * 0.7f,
                   "a cutoff past fs/4 is clamped, output bounded");
  udsp_svf_reset(&s);
  UDSP_TEST_ASSERT(s.lp[0] == 0.0f && s.bp[1] == 0.0f,
                   "reset clears the state");
}

int main(void) {
  test_dc_passes_low_and_not_high();
  test_split_blocks_identical();
  test_cutoff_clamped_stable();
  return udsp_test_summary();
}
