#include "udsp_test.h"

#include "scope.h"

#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define UDSP_SCOPE_TEST_N 4096

static void test_trigger_finds_rising_crossing(void) {
  static float mono[UDSP_SCOPE_TEST_N];
  /* a slow sine (a handful of cycles across the whole buffer) so there's
   * an unambiguous single rising zero-crossing well before the tail. */
  for (int i = 0; i < UDSP_SCOPE_TEST_N; i++) {
    mono[i] = sinf(2.0f * (float)M_PI * 5.0f * (float)i / (float)UDSP_SCOPE_TEST_N);
  }
  float trace[UDSP_SCOPE_POINTS];
  udsp_scope_extract(mono, UDSP_SCOPE_TEST_N, trace);

  UDSP_TEST_ASSERT(trace[0] >= -0.05f && trace[0] < 0.3f,
                  "trace starts at/near a rising zero-crossing, not an arbitrary phase (got %f)", trace[0]);
  UDSP_TEST_ASSERT(trace[10] > trace[0], "trace is rising just after the trigger point");
}

static void test_untriggered_fallback(void) {
  static float mono[UDSP_SCOPE_TEST_N];
  for (int i = 0; i < UDSP_SCOPE_TEST_N; i++) {
    mono[i] = 0.5f; /* constant DC -- never crosses zero */
  }
  float trace[UDSP_SCOPE_POINTS];
  udsp_scope_extract(mono, UDSP_SCOPE_TEST_N, trace);
  UDSP_TEST_ASSERT(trace[0] == 0.5f && trace[UDSP_SCOPE_POINTS - 1] == 0.5f,
                  "no trigger found anywhere -> falls back to the most recent samples, unmodified");
}

static void test_short_window_zero_pads(void) {
  float mono[10];
  for (int i = 0; i < 10; i++) {
    mono[i] = 1.0f;
  }
  float trace[UDSP_SCOPE_POINTS];
  udsp_scope_extract(mono, 10, trace);
  UDSP_TEST_ASSERT(trace[9] == 1.0f, "a window shorter than UDSP_SCOPE_POINTS keeps its real samples");
  UDSP_TEST_ASSERT(trace[UDSP_SCOPE_POINTS - 1] == 0.0f,
                  "a window shorter than UDSP_SCOPE_POINTS zero-pads the remainder");
}

int main(void) {
  test_trigger_finds_rising_crossing();
  test_untriggered_fallback();
  test_short_window_zero_pads();
  return udsp_test_summary();
}
