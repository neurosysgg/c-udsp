#include "udsp_test.h"

#include "fx/gain.h"

#include <math.h>

static void test_unity_gain_is_a_no_op(void) {
  float buf[4] = {0.1f, -0.5f, 1.0f, -1.0f};
  udsp_gain_apply(buf, 4, 1.0f);
  UDSP_TEST_ASSERT(buf[0] == 0.1f && buf[1] == -0.5f && buf[2] == 1.0f &&
                       buf[3] == -1.0f,
                   "gain=1.0 leaves every sample unchanged");
}

static void test_zero_gain_silences(void) {
  float buf[4] = {0.1f, -0.5f, 1.0f, -1.0f};
  udsp_gain_apply(buf, 4, 0.0f);
  UDSP_TEST_ASSERT(buf[0] == 0.0f && buf[1] == 0.0f && buf[2] == 0.0f &&
                       buf[3] == 0.0f,
                   "gain=0.0 silences every sample");
}

static void test_half_gain_halves_every_sample(void) {
  float buf[3] = {1.0f, -0.4f, 0.2f};
  udsp_gain_apply(buf, 3, 0.5f);
  UDSP_TEST_ASSERT(buf[0] == 0.5f, "sample 0 halved");
  UDSP_TEST_ASSERT(buf[1] == -0.2f, "sample 1 halved");
  UDSP_TEST_ASSERT(buf[2] == 0.1f, "sample 2 halved");
}

static void test_channel_agnostic_flat_stream(void) {
  /* same "doesn't need to know frames vs. channels" contract every other
   * function here follows -- a 6-sample buffer standing in for 3 stereo
   * frames gets every sample scaled uniformly, no special-casing. */
  float buf[6] = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
  udsp_gain_apply(buf, 6, 0.25f);
  for (int i = 0; i < 6; i++) {
    UDSP_TEST_ASSERT(buf[i] == 0.25f,
                     "sample %d scaled the same regardless of channel", i);
  }
}

static void test_stereo_gains_per_side(void) {
  float buf[18];
  for (int i = 0; i < 18; i++) {
    buf[i] = 1.0f;
  }
  udsp_gain_stereo(buf, 9, 0.5f, 0.25f);
  int ok = 1;
  for (int i = 0; i < 18; i++) {
    ok &= buf[i] == (i & 1 ? 0.25f : 0.5f);
  }
  UDSP_TEST_ASSERT(ok, "left frames by gl, right by gr, through the tail too");
}

static void test_mix_add(void) {
  float dst[11], src[11];
  for (int i = 0; i < 11; i++) {
    dst[i] = (float)i;
    src[i] = 2.0f;
  }
  udsp_mix_add(dst, src, 11, 0.25f);
  int ok = 1;
  for (int i = 0; i < 11; i++) {
    ok &= dst[i] == (float)i + 0.5f;
  }
  UDSP_TEST_ASSERT(ok, "dst += src * gain");
}

static void test_pan_laws(void) {
  float l, r;
  udsp_pan_gains(0.0f, UDSP_PAN_BALANCE, &l, &r);
  UDSP_TEST_ASSERT(l == 1.0f && r == 1.0f, "balance: unity in the centre");
  udsp_pan_gains(1.0f, UDSP_PAN_BALANCE, &l, &r);
  UDSP_TEST_ASSERT(fabsf(l) < 1e-6f && r == 1.0f,
                   "balance: hard right silences the left");
  udsp_pan_gains(-0.5f, UDSP_PAN_BALANCE, &l, &r);
  UDSP_TEST_ASSERT(l == 1.0f && fabsf(r - 0.70710678f) < 1e-6f,
                   "balance: half left, the right at -3 dB");
  udsp_pan_gains(0.0f, UDSP_PAN_EQUAL_POWER, &l, &r);
  UDSP_TEST_ASSERT(fabsf(l - 0.70710678f) < 1e-6f &&
                       fabsf(r - 0.70710678f) < 1e-6f,
                   "equal power: -3 dB each in the centre");
  udsp_pan_gains(-3.0f, UDSP_PAN_EQUAL_POWER, &l, &r);
  UDSP_TEST_ASSERT(l == 1.0f && fabsf(r) < 1e-6f, "a pan past -1 is -1");
}

int main(void) {
  test_unity_gain_is_a_no_op();
  test_zero_gain_silences();
  test_half_gain_halves_every_sample();
  test_channel_agnostic_flat_stream();
  test_stereo_gains_per_side();
  test_mix_add();
  test_pan_laws();
  return udsp_test_summary();
}
