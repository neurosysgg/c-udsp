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

/* The period estimator's job is to make one cycle fill trace_out, so a
 * periodic signal holds still instead of scrolling. That's directly
 * observable without reaching into the estimator: extract twice from the
 * same signal at two *different* phases, and if the period lock is
 * right, both traces show the same single cycle -- shape-identical
 * despite the phase difference. A wrong period (an octave multiple, or
 * the short-lag shoulder the first-threshold-crossing approach used to
 * land on) makes the two traces disagree, since they'd be showing
 * different fractions or multiples of a cycle from different starts. */
static void test_period_lock_holds_the_trace_still(void) {
  static float a[UDSP_SCOPE_TEST_N];
  static float b[UDSP_SCOPE_TEST_N];
  /* two harmonics, so a wrong octave is a real possibility rather than
   * something a pure sine would hide */
  const float freq = 440.0f, sr = 44100.0f;
  for (int i = 0; i < UDSP_SCOPE_TEST_N; i++) {
    float t = (float)i / sr;
    float u = (float)(i + 37) / sr; /* same signal, shifted phase */
    a[i] = 0.7f * sinf(2.0f * (float)M_PI * freq * t) +
           0.3f * sinf(2.0f * (float)M_PI * 2.0f * freq * t);
    b[i] = 0.7f * sinf(2.0f * (float)M_PI * freq * u) +
           0.3f * sinf(2.0f * (float)M_PI * 2.0f * freq * u);
  }

  float trace_a[UDSP_SCOPE_POINTS], trace_b[UDSP_SCOPE_POINTS];
  udsp_scope_extract(a, UDSP_SCOPE_TEST_N, trace_a);
  udsp_scope_extract(b, UDSP_SCOPE_TEST_N, trace_b);

  float worst = 0.0f;
  for (int i = 0; i < UDSP_SCOPE_POINTS; i++) {
    float d = trace_a[i] - trace_b[i];
    if (d < 0.0f) {
      d = -d;
    }
    if (d > worst) {
      worst = d;
    }
  }
  UDSP_TEST_ASSERT(worst < 0.15f,
                  "the same tone at two phases yields the same single-cycle "
                  "trace (max deviation %f) -- the trace holds still",
                  worst);

  /* and it really is *one* cycle, not several: a single period of this
   * waveform crosses zero going upward exactly once across the trace. */
  int rising = 0;
  for (int i = 1; i < UDSP_SCOPE_POINTS; i++) {
    if (trace_a[i - 1] <= 0.0f && trace_a[i] > 0.0f) {
      rising++;
    }
  }
  UDSP_TEST_ASSERT(rising == 1,
                  "exactly one rising zero-crossing across the trace, i.e. one "
                  "period fills it (got %d) -- not an octave multiple", rising);
}

/* white-ish noise has no period to find; the estimator must say so
 * (falling back to UDSP_SCOPE_POINTS) rather than locking onto whichever
 * lag rounding noise happened to favour. Asserted via the fallback's own
 * signature: a UDSP_SCOPE_POINTS-sample period resamples 1:1, so the
 * trace is the post-trigger samples verbatim. */
static void test_aperiodic_signal_falls_back(void) {
  static float mono[UDSP_SCOPE_TEST_N];
  unsigned seed = 22222u;
  for (int i = 0; i < UDSP_SCOPE_TEST_N; i++) {
    seed = seed * 1103515245u + 12345u;
    mono[i] = (float)((seed >> 16) & 0x7fff) / 16384.0f - 1.0f;
  }
  float trace[UDSP_SCOPE_POINTS];
  udsp_scope_extract(mono, UDSP_SCOPE_TEST_N, trace);

  /* every trace sample must be an actual sample from the buffer (1:1
   * resampling), not an interpolated blend of two */
  int all_real = 1;
  for (int i = 0; i < UDSP_SCOPE_POINTS && all_real; i++) {
    int found = 0;
    for (int j = 0; j < UDSP_SCOPE_TEST_N && !found; j++) {
      if (mono[j] == trace[i]) {
        found = 1;
      }
    }
    all_real = found;
  }
  UDSP_TEST_ASSERT(all_real,
                  "noise finds no period and falls back to 1:1 resampling of "
                  "the post-trigger samples, rather than locking onto noise");
}

int main(void) {
  test_trigger_finds_rising_crossing();
  test_untriggered_fallback();
  test_short_window_zero_pads();
  test_period_lock_holds_the_trace_still();
  test_aperiodic_signal_falls_back();
  return udsp_test_summary();
}
