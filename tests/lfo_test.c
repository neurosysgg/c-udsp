#include "udsp_test.h"

#include <math.h>

#include "synth/lfo.h"
#include "synth/osc.h"

static void test_shapes(void) {
  static const float WANT[UDSP_LFO_RANDOM][4] = {
      {0.0f, 1.0f, 0.0f, -1.0f},  /* sine (within 1e-5) */
      {0.0f, 1.0f, 0.0f, -1.0f},  /* triangle */
      {-1.0f, -0.5f, 0.0f, 0.5f}, /* saw up */
      {1.0f, 0.5f, 0.0f, -0.5f},  /* saw down */
      {1.0f, 1.0f, -1.0f, -1.0f}, /* square */
  };
  for (int s = 0; s < UDSP_LFO_RANDOM; s++) {
    UDSP_LFO l;
    udsp_lfo_init(&l, (UDSP_LFO_SHAPE)s, 1);
    l.inc = 1u << 30; /* a quarter turn a frame */
    int ok = 1;
    for (int q = 0; q < 4; q++) {
      ok &= fabsf(udsp_lfo_value(&l) - WANT[s][q]) < 1e-5f;
      udsp_lfo_advance(&l, 1);
    }
    UDSP_TEST_ASSERT(ok && l.phase == 0, "shape %d at the quarter turns", s);
  }
}

static void test_random_held_for_a_turn(void) {
  UDSP_LFO l;
  udsp_lfo_init(&l, UDSP_LFO_RANDOM, 7);
  l.inc = udsp_osc_inc_hz(1.0, 1000); /* a turn in 1000 frames */
  float first = udsp_lfo_value(&l);
  udsp_lfo_advance(&l, 999);
  UDSP_TEST_ASSERT(udsp_lfo_value(&l) == first && first >= -1.0f &&
                       first <= 1.0f,
                   "held through the turn");
  udsp_lfo_advance(&l, 2);
  UDSP_TEST_ASSERT(udsp_lfo_value(&l) != first, "a new value past it");
  float now = udsp_lfo_value(&l);
  udsp_lfo_reset(&l, 0);
  UDSP_TEST_ASSERT(l.phase == 0 && udsp_lfo_value(&l) != now,
                   "a reset starts a new one");
}

/* the smooth shape: at each turn's start the value the last turn ended
 * on, still there (a smoothstep), never past the two values it joins,
 * half way between them half way through */
static void test_smooth(void) {
  UDSP_LFO l;
  udsp_lfo_init(&l, UDSP_LFO_SMOOTH, 7);
  l.inc = 1u << 22; /* a turn in 1024 frames */
  float prev = udsp_lfo_value(&l), jump = 0.0f;
  int ok = prev == l.from && l.from != l.held, within = 1, half = 1;
  for (int turn = 0; turn < 50; turn++) {
    const float a = l.from, b = l.held;
    const float lo = a < b ? a : b, hi = a < b ? b : a;
    for (int i = 0; i < 1024; i++) {
      udsp_lfo_advance(&l, 1);
      const float v = udsp_lfo_value(&l);
      jump = fmaxf(jump, fabsf(v - prev));
      prev = v;
      within &= v >= lo - 1e-6f && v <= hi + 1e-6f && v >= -1.0f && v <= 1.0f;
      if (i == 511) {
        half &= fabsf(v - 0.5f * (a + b)) < 1e-3f * fabsf(b - a) + 1e-6f;
      }
    }
    ok &= l.from == b;
  }
  /* the steepest step: 1.5 x the largest gap (2) over 1024 frames */
  UDSP_TEST_ASSERT(ok && within && half && jump < 2.93e-3f,
                   "smooth: turns join (%d), within (%d), half way (%d), "
                   "steepest %.2e",
                   ok, within, half, jump);
  udsp_lfo_init(&l, UDSP_LFO_SMOOTH, 7);
  l.inc = 1u << 22;
  udsp_lfo_advance(&l, 250);
  const float at = udsp_lfo_value(&l);
  udsp_lfo_reset(&l, 0);
  UDSP_TEST_ASSERT(udsp_lfo_value(&l) == at,
                   "a reset goes on from where it is");
  UDSP_LFO r;
  udsp_lfo_init(&r, UDSP_LFO_RANDOM, 7);
  udsp_lfo_init(&l, UDSP_LFO_SMOOTH, 7);
  UDSP_TEST_ASSERT(l.from == r.held, "the same seed's first value");
}

int main(void) {
  test_shapes();
  test_random_held_for_a_turn();
  test_smooth();
  return udsp_test_summary();
}
