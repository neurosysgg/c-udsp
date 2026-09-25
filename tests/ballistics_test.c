#include "udsp_test.h"

#include "ballistics.h"

#include <math.h>

static int near(float a, float b) { return fabsf(a - b) < 0.001f; }

static void test_ballistics_time_constant(void) {
  float level[2] = {0.0f, 1.0f};
  float target[2] = {1.0f, 0.0f};
  udsp_ballistics(level, target, 2, 0.1f, 0.1f, 0.5f);
  UDSP_TEST_ASSERT(near(level[0], 1.0f - expf(-1.0f)),
                   "rising covers 1 - 1/e of the gap in one attack time "
                   "constant (%f)", level[0]);
  UDSP_TEST_ASSERT(near(level[1], expf(-0.2f)),
                   "falling uses the release time constant (%f)", level[1]);
}

static void test_ballistics_frame_rate_independent(void) {
  float a = 0.0f, b = 0.0f, one = 1.0f;
  udsp_ballistics(&a, &one, 1, 0.1f, 0.05f, 0.05f);
  for (int i = 0; i < 4; i++) {
    udsp_ballistics(&b, &one, 1, 0.025f, 0.05f, 0.05f);
  }
  UDSP_TEST_ASSERT(near(a, b),
                   "one 100 ms step lands where four 25 ms steps do "
                   "(%f vs %f)", a, b);
}

static void test_ballistics_instant(void) {
  float level = -60.0f, target = -12.0f;
  udsp_ballistics(&level, &target, 1, 0.01f, 0.0f, 1.0f);
  UDSP_TEST_ASSERT(level == -12.0f, "a zero attack jumps to the target");
  target = -60.0f;
  udsp_ballistics(&level, &target, 1, 0.01f, 1.0f, -1.0f);
  UDSP_TEST_ASSERT(level == -60.0f, "a negative release jumps too");
}

static void test_peak_hold(void) {
  float peak = -60.0f, age = 0.0f, level = -10.0f;
  udsp_peak_hold(&peak, &age, &level, 1, 0.1f, 0.5f, 20.0f);
  UDSP_TEST_ASSERT(peak == -10.0f && age == 0.0f,
                   "the peak jumps up to a louder level");
  level = -40.0f;
  for (int i = 0; i < 5; i++) {
    udsp_peak_hold(&peak, &age, &level, 1, 0.1f, 0.5f, 20.0f);
  }
  UDSP_TEST_ASSERT(near(peak, -10.0f), "it holds for hold_s (%f)", peak);
  udsp_peak_hold(&peak, &age, &level, 1, 0.1f, 0.5f, 20.0f);
  UDSP_TEST_ASSERT(near(peak, -12.0f),
                   "then falls fall_per_s per second (%f)", peak);
  for (int i = 0; i < 50; i++) {
    udsp_peak_hold(&peak, &age, &level, 1, 0.1f, 0.5f, 20.0f);
  }
  UDSP_TEST_ASSERT(peak == -40.0f, "never below the level (%f)", peak);
}

static void test_peak_hold_partial_step(void) {
  float peak = 0.0f, age = 0.4f, level = -60.0f;
  udsp_peak_hold(&peak, &age, &level, 1, 0.2f, 0.5f, 10.0f);
  UDSP_TEST_ASSERT(near(peak, -1.0f),
                   "only the part of a step past the hold falls (%f)", peak);
}

int main(void) {
  test_ballistics_time_constant();
  test_ballistics_frame_rate_independent();
  test_ballistics_instant();
  test_peak_hold();
  test_peak_hold_partial_step();
  return udsp_test_summary();
}
