#include "udsp_test.h"

#include "synth/ramp.h"

static void ones(float *b, int n) {
  for (int i = 0; i < n; i++) {
    b[i] = 1.0f;
  }
}

static void test_attack_and_release(void) {
  UDSP_RAMP r;
  udsp_ramp_init(&r, 4);
  UDSP_TEST_ASSERT(udsp_ramp_silent(&r), "a new ramp is silent");
  r.target = 1.0f;
  float b[24];
  ones(b, 24);
  udsp_ramp_stereo(&r, b, 12);
  UDSP_TEST_ASSERT(b[0] == 0.25f && b[1] == 0.25f && b[2] == 0.5f &&
                       b[6] == 1.0f && b[23] == 1.0f,
                   "four frames up, then held at the target");
  r.target = 0.5f;
  ones(b, 24);
  udsp_ramp_stereo(&r, b, 12);
  UDSP_TEST_ASSERT(b[0] == 0.75f && b[2] == 0.5f && b[23] == 0.5f,
                   "down to a lower target, from where it was");
  r.target = 0.0f;
  ones(b, 24);
  udsp_ramp_stereo(&r, b, 12);
  UDSP_TEST_ASSERT(b[2] == 0.0f && udsp_ramp_silent(&r),
                   "and released to silence");
}

static void test_overshoot_lands_on_target(void) {
  UDSP_RAMP r;
  udsp_ramp_init(&r, 3);
  r.target = 0.5f;
  float b[8];
  ones(b, 8);
  udsp_ramp_stereo(&r, b, 4);
  UDSP_TEST_ASSERT(b[2] == 0.5f && b[6] == 0.5f && r.value == 0.5f,
                   "a step past the target stops on it");
}

static void test_held_level_scales(void) {
  UDSP_RAMP r = {0.8f, 0.8f, 0.25f};
  float b[19];
  ones(b, 19);
  udsp_ramp_stereo(&r, b, 9);
  int ok = 1;
  for (int i = 0; i < 18; i++) {
    ok &= b[i] == 0.8f;
  }
  UDSP_TEST_ASSERT(ok && b[18] == 1.0f,
                   "a settled ramp scales exactly n frames");
}

int main(void) {
  test_attack_and_release();
  test_overshoot_lands_on_target();
  test_held_level_scales();
  return udsp_test_summary();
}
