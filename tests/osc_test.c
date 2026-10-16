#include "udsp_test.h"

#include <math.h>
#include <string.h>

#include "osc.h"

#define RATE 48000

static UDSP_OSC make(UDSP_OSC_WAVE w, uint32_t inc) {
  UDSP_OSC o;
  udsp_osc_init(&o, w, 1);
  o.inc = inc;
  return o;
}

static void test_increment(void) {
  UDSP_TEST_ASSERT(udsp_osc_inc_millihz(1000000, RATE) == 89478485,
                   "1 kHz at 48 kHz is 2^32/48 per sample, truncated");
  UDSP_TEST_ASSERT(udsp_osc_inc_hz(1000.0, RATE) == 89478485,
                   "the same from hertz, rounded");
  UDSP_TEST_ASSERT(udsp_osc_inc_hz(12000.0, RATE) == 1u << 30,
                   "a quarter of the rate is a quarter turn");
}

static void test_sine_matches_sin(void) {
  UDSP_OSC o = make(UDSP_OSC_SINE, udsp_osc_inc_millihz(1000000, RATE));
  float out[1003];
  udsp_osc_render(&o, out, 1003);
  float worst = 0.0f;
  for (int i = 0; i < 1003; i++) {
    float d = fabsf(out[i] - (float)sin(6.283185307179586 * 1000.0 * i / RATE));
    worst = d > worst ? d : worst;
  }
  UDSP_TEST_ASSERT(worst < 1e-5f,
                   "a 1 kHz sine within 1e-5 of sin() (worst %g)",
                   (double)worst);
}

static void test_triangle_corners_exact(void) {
  UDSP_OSC o = make(UDSP_OSC_TRIANGLE, 1u << 30);
  float out[8];
  udsp_osc_render(&o, out, 8);
  const float want[8] = {-1.0f, 0.0f, 1.0f, 0.0f, -1.0f, 0.0f, 1.0f, 0.0f};
  UDSP_TEST_ASSERT(memcmp(out, want, sizeof(want)) == 0,
                   "a quarter turn per sample walks -1 0 1 0 exactly");
}

static void test_saw_and_square_away_from_edges(void) {
  /* 1/16 turn per sample: the PolyBLEP touches only the samples at t = 0
   * and t = 0.5 (square) */
  UDSP_OSC o = make(UDSP_OSC_SAW, 1u << 28);
  float out[16];
  udsp_osc_render(&o, out, 16);
  UDSP_TEST_ASSERT(out[4] == -0.5f && out[8] == 0.0f && out[12] == 0.5f,
                   "the saw is 2t - 1 between edges");
  UDSP_TEST_ASSERT(out[0] == 0.0f, "and its edge sample is smoothed to 0");
  o = make(UDSP_OSC_SQUARE, 1u << 28);
  udsp_osc_render(&o, out, 16);
  UDSP_TEST_ASSERT(out[3] == 1.0f && out[11] == -1.0f,
                   "the square is +-1 between edges");
  UDSP_TEST_ASSERT(out[0] == 0.0f && out[8] == 0.0f,
                   "its edges are smoothed to the midpoint");
}

static void test_dc_saw_without_increment(void) {
  UDSP_OSC o = make(UDSP_OSC_SAW, 0);
  float out[9];
  udsp_osc_render(&o, out, 9);
  int ok = 1;
  for (int i = 0; i < 9; i++) {
    ok &= out[i] == -1.0f;
  }
  UDSP_TEST_ASSERT(ok, "inc 0: no PolyBLEP (dt 0), the saw sits at -1");
}

/* a block split anywhere gives the same bits: the vector path and the
 * scalar tail are the same arithmetic */
static void test_split_renders_identical(void) {
  for (int w = 0; w < UDSP_OSC_WAVES; w++) {
    UDSP_OSC a = make((UDSP_OSC_WAVE)w, udsp_osc_inc_hz(3277.7, RATE));
    UDSP_OSC b = a;
    float whole[257], parts[257];
    udsp_osc_render(&a, whole, 257);
    udsp_osc_render(&b, parts, 3);
    udsp_osc_render(&b, parts + 3, 101);
    udsp_osc_render(&b, parts + 104, 153);
    UDSP_TEST_ASSERT(memcmp(whole, parts, sizeof(whole)) == 0 &&
                         a.phase == b.phase && a.rng == b.rng,
                     "%s: 257 samples at once == in 3 + 101 + 153",
                     udsp_osc_wave_name((UDSP_OSC_WAVE)w));
  }
}

static void test_phase_advances(void) {
  UDSP_OSC o = make(UDSP_OSC_SINE, 1000);
  o.phase = 0xffffff00u;
  float out[5];
  udsp_osc_render(&o, out, 5);
  UDSP_TEST_ASSERT(o.phase == 0xffffff00u + 5000u,
                   "the phase moves n * inc, wrapping");
}

static void test_noise(void) {
  UDSP_OSC o = make(UDSP_OSC_NOISE, 0);
  float out[4096];
  udsp_osc_render(&o, out, 4096);
  double sum = 0.0;
  int in_range = 1;
  for (int i = 0; i < 4096; i++) {
    sum += out[i];
    in_range &= out[i] >= -1.0f && out[i] < 1.0f;
  }
  UDSP_TEST_ASSERT(in_range, "noise stays in [-1, 1)");
  UDSP_TEST_ASSERT(fabs(sum / 4096) < 0.05, "and averages near 0 (%g)",
                   sum / 4096);
  UDSP_OSC p = make(UDSP_OSC_NOISE, 0);
  p.rng = 2;
  float other[4];
  udsp_osc_render(&p, other, 4);
  UDSP_TEST_ASSERT(memcmp(out, other, sizeof(other)) != 0,
                   "another seed is another sequence");
}

/* phase modulation: none is the plain render; a constant offset is the
 * wave shifted; a modulator's swing splits anywhere the same */
static void test_pm(void) {
  float pm[301], want[301], got[301];
  for (int w = 0; w < UDSP_OSC_WAVES; w++) {
    UDSP_OSC a = make((UDSP_OSC_WAVE)w, udsp_osc_inc_hz(1234.5, RATE));
    UDSP_OSC b = a;
    memset(pm, 0, sizeof(pm));
    udsp_osc_render(&a, want, 301);
    udsp_osc_render_pm(&b, pm, got, 301);
    UDSP_TEST_ASSERT(memcmp(want, got, sizeof(got)) == 0 && a.phase == b.phase,
                     "%s: pm 0 is the plain render, to the bit",
                     udsp_osc_wave_name((UDSP_OSC_WAVE)w));
  }
  UDSP_OSC a = make(UDSP_OSC_SINE, udsp_osc_inc_millihz(1000000, RATE));
  for (int i = 0; i < 301; i++) {
    pm[i] = -1.75f; /* a quarter turn, from below zero */
  }
  udsp_osc_render_pm(&a, pm, got, 301);
  float worst = 0.0f;
  for (int i = 0; i < 301; i++) {
    float d = fabsf(got[i] - (float)cos(6.283185307179586 * 1000.0 * i / RATE));
    worst = d > worst ? d : worst;
  }
  UDSP_TEST_ASSERT(worst < 1e-5f, "a quarter turn on is the cosine (%g)",
                   (double)worst);
  for (int i = 0; i < 301; i++) {
    pm[i] = 1.3f * (float)sin(i * 0.07);
  }
  UDSP_OSC c = make(UDSP_OSC_SAW, udsp_osc_inc_hz(220.0, RATE)), d = c;
  float parts[301];
  udsp_osc_render_pm(&c, pm, got, 301);
  udsp_osc_render_pm(&d, pm, parts, 5);
  udsp_osc_render_pm(&d, pm + 5, parts + 5, 296);
  int bounded = 1;
  for (int i = 0; i < 301; i++) {
    bounded &= fabsf(got[i]) <= 1.1f;
  }
  UDSP_TEST_ASSERT(memcmp(got, parts, sizeof(got)) == 0 && bounded,
                   "modulated past a turn either way: split anywhere, bounded");
  UDSP_TEST_ASSERT(fabsf(udsp_osc_sine(0.25f) - 1.0f) < 1e-5f &&
                       udsp_osc_sine(0.0f) == 0.0f,
                   "the sine at a quarter turn and at zero");
}

static void test_names(void) {
  int ok = 1;
  for (int w = 0; w < UDSP_OSC_WAVES; w++) {
    ok &= udsp_osc_wave_by_name(udsp_osc_wave_name((UDSP_OSC_WAVE)w)) ==
          (UDSP_OSC_WAVE)w;
  }
  UDSP_TEST_ASSERT(ok, "every wave's name reads back");
  UDSP_TEST_ASSERT(udsp_osc_wave_by_name("pulse") == UDSP_OSC_WAVES,
                   "an unknown name is UDSP_OSC_WAVES");
}

int main(void) {
  test_increment();
  test_sine_matches_sin();
  test_triangle_corners_exact();
  test_saw_and_square_away_from_edges();
  test_dc_saw_without_increment();
  test_split_renders_identical();
  test_pm();
  test_phase_advances();
  test_noise();
  test_names();
  return udsp_test_summary();
}
