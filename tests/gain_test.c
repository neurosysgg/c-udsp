#include "udsp_test.h"

#include "gain.h"

static void test_unity_gain_is_a_no_op(void) {
  float buf[4] = {0.1f, -0.5f, 1.0f, -1.0f};
  udsp_gain_apply(buf, 4, 1.0f);
  UDSP_TEST_ASSERT(buf[0] == 0.1f && buf[1] == -0.5f && buf[2] == 1.0f && buf[3] == -1.0f,
                  "gain=1.0 leaves every sample unchanged");
}

static void test_zero_gain_silences(void) {
  float buf[4] = {0.1f, -0.5f, 1.0f, -1.0f};
  udsp_gain_apply(buf, 4, 0.0f);
  UDSP_TEST_ASSERT(buf[0] == 0.0f && buf[1] == 0.0f && buf[2] == 0.0f && buf[3] == 0.0f,
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
    UDSP_TEST_ASSERT(buf[i] == 0.25f, "sample %d scaled the same regardless of channel", i);
  }
}

int main(void) {
  test_unity_gain_is_a_no_op();
  test_zero_gain_silences();
  test_half_gain_halves_every_sample();
  test_channel_agnostic_flat_stream();
  return udsp_test_summary();
}
