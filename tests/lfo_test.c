#include "udsp_test.h"

#include <math.h>

#include "synth/lfo.h"
#include "synth/osc.h"

static void test_shapes(void) {
  static const float WANT[UDSP_LFO_SHAPES - 1][4] = {
      {0.0f, 1.0f, 0.0f, -1.0f},  /* sine (within 1e-5) */
      {0.0f, 1.0f, 0.0f, -1.0f},  /* triangle */
      {-1.0f, -0.5f, 0.0f, 0.5f}, /* saw up */
      {1.0f, 0.5f, 0.0f, -0.5f},  /* saw down */
      {1.0f, 1.0f, -1.0f, -1.0f}, /* square */
  };
  for (int s = 0; s < UDSP_LFO_SHAPES - 1; s++) {
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

int main(void) {
  test_shapes();
  test_random_held_for_a_turn();
  return udsp_test_summary();
}
