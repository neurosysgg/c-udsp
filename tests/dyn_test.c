#include "udsp_test.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "fx/dyn.h"

#define RATE 48000

static uint32_t g_rng = 1;
static float noise(void) {
  g_rng ^= g_rng << 13;
  g_rng ^= g_rng >> 17;
  g_rng ^= g_rng << 5;
  return (float)(int32_t)g_rng * (1.0f / 2147483648.0f);
}

static void test_log_exp(void) {
  double worst_l = 0.0, worst_e = 0.0;
  for (int i = 0; i < 200000; i++) {
    float x = (float)ldexp(1.0 + fabs(noise()), (i % 120) - 60);
    double e = fabs(udsp_log2f(x) - log2((double)x));
    worst_l = e > worst_l ? e : worst_l;
    float y = noise() * 100.0f;
    double r = fabs(udsp_exp2f(y) / exp2((double)y) - 1.0);
    worst_e = r > worst_e ? r : worst_e;
  }
  UDSP_TEST_ASSERT(
      worst_l < 4e-6,
      "log2f: %.2g off at worst (on 2^-60..2^60: a float's step at 60)",
      worst_l);
  UDSP_TEST_ASSERT(worst_e < 4e-7, "exp2f: %.2g relative at worst", worst_e);
  UDSP_TEST_ASSERT(udsp_log2f(1.0f) == 0.0f && udsp_exp2f(0.0f) == 1.0f &&
                       udsp_exp2f(3.0f) == 8.0f,
                   "log2 1 = 0, 2^0 = 1, 2^3 = 8 exactly");
}

/* a square wave of amplitude amp as the key: its peak and RMS are amp
 * every frame, so the gain settles where the curve says */
static float comp_settled(UDSP_COMP *c, float amp, size_t frames) {
  static float key[2 * RATE], gain[RATE];
  for (size_t i = 0; i < frames; i++) {
    float v = (i / 24) % 2 ? amp : -amp;
    key[2 * i] = v;
    key[2 * i + 1] = 0.5f * v;
  }
  udsp_comp_gain(c, key, gain, frames);
  return gain[frames - 1];
}

static void test_comp(void) {
  UDSP_COMP c;
  udsp_comp_reset(&c);
  udsp_comp_set(&c, -20.0f, 4.0f, 0.0f, 1.0f, 50.0f, 0, RATE);
  /* -8 dB in: 12 over, 3 out at 4:1, 9 dB of reduction */
  float g = comp_settled(&c, udsp_gain_from_db(-8.0f), RATE);
  UDSP_TEST_ASSERT(fabsf(udsp_db_from_gain(g) + 9.0f) < 0.01f,
                   "4:1, 12 dB over: 9 dB down (%.3f)",
                   (double)udsp_db_from_gain(g));
  UDSP_TEST_ASSERT(fabsf(udsp_comp_curve(&c, -8.0f) + 17.0f) < 1e-4f,
                   "the curve: -8 in, -17 out");
  udsp_comp_reset(&c);
  udsp_comp_set(&c, -20.0f, 4.0f, 0.0f, 1.0f, 50.0f, 1, RATE);
  g = comp_settled(&c, udsp_gain_from_db(-8.0f), RATE);
  UDSP_TEST_ASSERT(fabsf(udsp_db_from_gain(g) + 9.0f) < 0.05f,
                   "RMS of a square: the same (%.3f)",
                   (double)udsp_db_from_gain(g));
  udsp_comp_reset(&c);
  udsp_comp_set(&c, -20.0f, 4.0f, 6.0f, 1.0f, 50.0f, 0, RATE);
  g = comp_settled(&c, udsp_gain_from_db(-20.0f), RATE);
  UDSP_TEST_ASSERT(fabsf(udsp_db_from_gain(g) + 0.5625f) < 0.01f,
                   "a 6 dB knee at the threshold: 3/4 (3)^2 / 12 = 0.56 dB");
  g = comp_settled(&c, udsp_gain_from_db(-24.0f), RATE);
  UDSP_TEST_ASSERT(udsp_db_from_gain(g) > -0.01f,
                   "under the knee: nothing (%.4f)",
                   (double)udsp_db_from_gain(g));
  udsp_comp_set(&c, -20.0f, 1e9f, 0.0f, 1.0f, 50.0f, 0, RATE);
  g = comp_settled(&c, udsp_gain_from_db(-5.0f), RATE);
  UDSP_TEST_ASSERT(fabsf(udsp_db_from_gain(g) + 15.0f) < 0.01f,
                   "ratio infinity: held at the threshold");

  /* attack: from rest, a step 20 dB over at 4:1 (15 dB to take): 63 %
   * of it a time constant later */
  udsp_comp_reset(&c);
  udsp_comp_set(&c, -20.0f, 4.0f, 0.0f, 10.0f, 200.0f, 0, RATE);
  g = comp_settled(&c, udsp_gain_from_db(0.0f), RATE / 100);
  UDSP_TEST_ASSERT(fabsf(udsp_db_from_gain(g) + 15.0f * 0.632f) < 0.1f,
                   "10 ms attack: 63 %% of 15 dB after 10 ms (%.2f)",
                   (double)udsp_db_from_gain(g));
  const float from = c.red;
  g = comp_settled(&c, udsp_gain_from_db(-40.0f), RATE / 5);
  UDSP_TEST_ASSERT(fabsf(udsp_db_from_gain(g) + from * 0.368f) < 0.05f,
                   "200 ms release: 37 %% left after 200 ms (%.2f)",
                   (double)udsp_db_from_gain(g));
}

/* the second knee (Ozone 11's "Lim" stage): past its threshold the slope
 * 1 / ratio2, its term there whatever the first does; the post smoothing
 * a one-pole more on the reduction */
static void test_comp_second(void) {
  UDSP_COMP c;
  udsp_comp_reset(&c);
  udsp_comp_set(&c, -12.0f, 3.0f, 0.0f, 1.0f, 50.0f, 0, RATE);
  udsp_comp_set_second(&c, 0.0f, 10.0f, 0.0f, 1.0f, 50.0f, 0.0f, RATE);
  /* +3 in: 15 over the first at 3:1 (-12 + 4 = -8 at 0 dB in), then 3 over
   * the second at 10:1: -8 + 0.3 = -7.7 out */
  UDSP_TEST_ASSERT(fabsf(udsp_comp_curve(&c, 3.0f) + 7.7f) < 1e-4f &&
                       fabsf(udsp_comp_curve(&c, -6.0f) + 10.0f) < 1e-4f,
                   "two knees: -6 in -10 out, +3 in -7.7 out (%.3f)",
                   (double)udsp_comp_curve(&c, 3.0f));
  float g = comp_settled(&c, udsp_gain_from_db(3.0f), RATE);
  UDSP_TEST_ASSERT(fabsf(udsp_db_from_gain(g) + 10.7f) < 0.01f,
                   "and settled there: 10.7 dB down (%.3f)",
                   (double)udsp_db_from_gain(g));
  /* the second under the first with a gentler ratio: an expansion of the
   * first's line (Ozone's law) -- -20 in, the first not yet acting */
  udsp_comp_set(&c, 0.0f, 3.0f, 0.0f, 1.0f, 50.0f, 0, RATE);
  udsp_comp_set_second(&c, -30.0f, 1.5f, 0.0f, 1.0f, 50.0f, 0.0f, RATE);
  UDSP_TEST_ASSERT(fabsf(udsp_comp_curve(&c, -20.0f) + 16.6667f) < 1e-3f,
                   "a second knee under the first at 1.5:1: -20 in, -16.67 "
                   "out (%.3f)",
                   (double)udsp_comp_curve(&c, -20.0f));
  /* the post one-pole slows a step: 10 ms attack and 10 ms more after */
  udsp_comp_reset(&c);
  udsp_comp_set(&c, -20.0f, 4.0f, 0.0f, 10.0f, 200.0f, 0, RATE);
  udsp_comp_set_smooth(&c, 10.0f, RATE);
  g = comp_settled(&c, udsp_gain_from_db(0.0f), RATE / 100);
  /* two equal poles: 1 - 2 / e of it a time constant later */
  UDSP_TEST_ASSERT(fabsf(udsp_db_from_gain(g) + 15.0f * 0.2642f) < 0.15f,
                   "two 10 ms poles: 26 %% of 15 dB after 10 ms (%.2f)",
                   (double)udsp_db_from_gain(g));
  udsp_comp_set(&c, -20.0f, 4.0f, 0.0f, 10.0f, 200.0f, 0, RATE);
  UDSP_TEST_ASSERT(!c.second && c.post == 0.0f && c.hold == 0.0f &&
                       c.block == 0,
                   "set again: the second knee, the smoothing, the hold, the "
                   "window off");
  /* the hold: a 1 kHz sine at -8 dB, 12 dB over at 4:1 (9 dB to take) --
   * without it the troughs let the reduction sag, with it not */
  static float key[2 * RATE], gain[RATE];
  for (int i = 0; i < RATE; i++) {
    key[2 * i] = key[2 * i + 1] =
        udsp_gain_from_db(-8.0f) *
        (float)sin(6.283185307179586 * 1000 * i / RATE);
  }
  udsp_comp_reset(&c);
  udsp_comp_set(&c, -20.0f, 4.0f, 0.0f, 30.0f, 100.0f, 0, RATE);
  udsp_comp_gain(&c, key, gain, RATE);
  const float sag = udsp_db_from_gain(gain[RATE - 1]);
  udsp_comp_reset(&c);
  udsp_comp_set(&c, -20.0f, 4.0f, 0.0f, 30.0f, 100.0f, 0, RATE);
  udsp_comp_set_hold(&c, 5.0f, RATE);
  udsp_comp_gain(&c, key, gain, RATE);
  const float held = udsp_db_from_gain(gain[RATE - 1]);
  UDSP_TEST_ASSERT(sag > -7.5f && fabsf(held + 9.0f) <
                                      0.3f, /* 5 ms: a little between peaks */
                   "a 1 kHz sine 12 dB over at 4:1: %.2f dB sagging, %.2f held",
                   (double)sag, (double)held);
  /* the window: a click at frame 10, blocks of 1 ms -- its block's peak,
   * then the last block's, then gone (attack and release near instant) */
  memset(key, 0, sizeof key);
  key[20] = key[21] = 1.0f;
  udsp_comp_reset(&c);
  udsp_comp_set(&c, -20.0f, 4.0f, 0.0f, 0.001f, 0.001f, 0, RATE);
  udsp_comp_set_window(&c, 1.0f, RATE);
  udsp_comp_gain(&c, key, gain, 200);
  UDSP_TEST_ASSERT(fabsf(udsp_db_from_gain(gain[95]) + 15.0f) < 0.01f &&
                       fabsf(udsp_db_from_gain(gain[96])) < 0.01f &&
                       fabsf(udsp_db_from_gain(gain[9])) < 0.01f,
                   "a click held through the next block: %.2f at 95, %.2f "
                   "at 96",
                   (double)udsp_db_from_gain(gain[95]),
                   (double)udsp_db_from_gain(gain[96]));
  /* the second's tail: taken at once, let go over its own ms -- a step
   * 6 dB over a 0 dB second knee at 10:1 (the first at 1:1: none) */
  for (int i = 0; i < RATE; i++) {
    key[2 * i] = key[2 * i + 1] = i < RATE / 2 ? 2.0f : 0.0f;
  }
  udsp_comp_reset(&c);
  udsp_comp_set(&c, -60.0f, 1.0f, 0.0f, 1.0f, 1.0f, 0, RATE);
  udsp_comp_set_second(&c, 0.0f, 10.0f, 0.0f, 0.01f, 0.01f, 20.0f, RATE);
  udsp_comp_gain(&c, key, gain, RATE);
  const float on = udsp_db_from_gain(gain[RATE / 2 - 1]);
  const float later = udsp_db_from_gain(gain[RATE / 2 + RATE / 50]);
  UDSP_TEST_ASSERT(fabsf(on + 5.418f) < 0.01f &&
                       fabsf(later + 5.418f / 2.71828f) < 0.05f,
                   "the tail: %.3f dB on, 1 / e of it 20 ms after (%.3f)",
                   (double)on, (double)later);
}

static void test_gate(void) {
  UDSP_GATE g;
  udsp_gate_set(&g, -30.0f, 6.0f, -100.0f, 0.1f, 10.0f, 5.0f, RATE);
  udsp_gate_reset(&g);
  static float key[2 * RATE], gain[RATE];
  /* 50 ms at -20, 50 ms at -33 (inside the return), 50 ms at -40 */
  const float lv[3] = {-20.0f, -33.0f, -40.0f};
  for (int i = 0; i < 7200; i++) {
    float a = udsp_gain_from_db(lv[i / 2400]);
    key[2 * i] = key[2 * i + 1] = (i / 12) % 2 ? a : -a;
  }
  memset(gain, 0, sizeof(gain));
  udsp_gate_gain(&g, key, gain, 7200);
  UDSP_TEST_ASSERT(gain[2399] > 0.999f, "over the threshold: open (%g)",
                   (double)gain[2399]);
  UDSP_TEST_ASSERT(gain[4799] > 0.999f,
                   "under it but inside the return: still open");
  UDSP_TEST_ASSERT(gain[4800 + 480 - 10] > 0.999f,
                   "under the return: held 10 ms");
  UDSP_TEST_ASSERT(gain[7199] == 0.0f, "then shut: silence (%g)",
                   (double)gain[7199]);
  udsp_gate_set(&g, -30.0f, 6.0f, -20.0f, 0.1f, 0.0f, 5.0f, RATE);
  udsp_gate_gain(&g, key + 2 * 4800, gain, 2400);
  UDSP_TEST_ASSERT(fabsf(udsp_db_from_gain(gain[2399]) + 20.0f) < 0.01f,
                   "a -20 dB floor: closed at -20");
}

static void test_limiter(void) {
  enum { CAP = 480 };
  static float mem[UDSP_LIMITER_FLOATS(CAP)];
  UDSP_LIMITER l;
  udsp_limiter_init(&l, mem, CAP);
  udsp_limiter_set(&l, 240, 0.5f, 50.0f, RATE);
  static float x[2 * RATE], y[2 * RATE];
  for (int i = 0; i < 2 * RATE; i++) {
    /* loud noise with bursts to 8x */
    x[i] = noise() * ((i / 4000) % 3 == 1 ? 8.0f : 1.0f);
  }
  memcpy(y, x, sizeof(x));
  for (int i = 0; i < RATE; i += 333) {
    udsp_limiter_stereo(&l, y + 2 * i, i + 333 <= RATE ? 333 : RATE - i);
  }
  float worst = 0.0f;
  for (int i = 0; i < 2 * RATE; i++) {
    worst = fabsf(y[i]) > worst ? fabsf(y[i]) : worst;
  }
  UDSP_TEST_ASSERT(worst <= 0.5f * (1.0f + 1e-6f),
                   "no sample past the ceiling (%.7f of 0.5)", (double)worst);
  UDSP_TEST_ASSERT(worst > 0.49f, "  ... and it gets there");

  /* under the ceiling: the input as it was, lookahead frames late */
  udsp_limiter_init(&l, mem, CAP);
  udsp_limiter_set(&l, 100, 1.0f, 50.0f, RATE);
  for (int i = 0; i < 2000; i++) {
    x[i] = 0.25f * noise();
  }
  memcpy(y, x, 2000 * sizeof(float));
  udsp_limiter_stereo(&l, y, 1000);
  int late = 1;
  for (int i = 0; i < 900; i++) {
    late &=
        y[2 * (i + 100)] == x[2 * i] && y[2 * (i + 100) + 1] == x[2 * i + 1];
  }
  for (int i = 0; i < 100; i++) {
    late &= y[2 * i] == 0.0f;
  }
  UDSP_TEST_ASSERT(late, "under the ceiling: the input, 100 frames late, "
                         "to the bit");

  /* a lone peak: the gain eases down over the lookahead, not a jump */
  udsp_limiter_init(&l, mem, CAP);
  udsp_limiter_set(&l, 48, 1.0f, 100.0f, RATE);
  for (int i = 0; i < 2000; i++) {
    x[2 * i] = x[2 * i + 1] = i == 1000 ? 4.0f : 0.5f;
  }
  udsp_limiter_stereo(&l, x, 2000);
  UDSP_TEST_ASSERT(fabsf(x[2 * 1048] - 1.0f) < 1e-6f,
                   "the peak, 48 late, at the ceiling (%g)",
                   (double)x[2 * 1048]);
  UDSP_TEST_ASSERT(x[2 * 1024] < 0.5f && x[2 * 1024] > 0.125f,
                   "halfway down the ramp before it (%g)", (double)x[2 * 1024]);
  UDSP_TEST_ASSERT(x[2 * 999] == 0.5f, "untouched before the lookahead");
}

/* the up / down compressor's settled gain on a steady key at level dB */
static float ud_settled(UDSP_UDCOMP *c, float level_db) {
  static float key[2 * RATE], gain[RATE];
  const float a = powf(10.0f, level_db / 20.0f);
  for (int i = 0; i < 2 * RATE; i++) {
    key[i] = a;
  }
  udsp_udcomp_reset(c);
  udsp_udcomp_gain(c, key, gain, RATE);
  return 20.0f * log10f(gain[RATE - 1]);
}

static void test_udcomp(void) {
  UDSP_UDCOMP c;
  udsp_udcomp_set(&c, -20.0f, 4.0f, -40.0f, 4.0f, 12.0f, -80.0f, 5.0f, 50.0f, 0,
                  RATE);
  const float above = ud_settled(&c, -10.0f), between = ud_settled(&c, -30.0f);
  const float below = ud_settled(&c, -46.0f), capped = ud_settled(&c, -70.0f);
  const float floor = ud_settled(&c, -90.0f);
  UDSP_TEST_ASSERT(fabsf(above + 7.5f) < 0.01f && fabsf(between) < 1e-4f,
                   "10 dB over at 4:1: -7.5 dB (%.3f); between: 0 (%.4f)",
                   (double)above, (double)between);
  UDSP_TEST_ASSERT(fabsf(below - 4.5f) < 0.01f && fabsf(capped - 12.0f) < 0.01f,
                   "6 dB under at 4:1 up: +4.5 dB (%.3f); 30 under: the 12 dB "
                   "cap (%.3f)",
                   (double)below, (double)capped);
  UDSP_TEST_ASSERT(fabsf(floor) < 1e-4f, "under the floor: no boost (%.4f)",
                   (double)floor);
  UDSP_TEST_ASSERT(fabsf(udsp_udcomp_curve(&c, -10.0f) + 17.5f) < 1e-4f &&
                       fabsf(udsp_udcomp_curve(&c, -46.0f) + 41.5f) < 1e-4f,
                   "the curve: -10 -> -17.5, -46 -> -41.5");
  /* a boost comes in at the release's pace, a cut at the attack's: 5 ms
   * vs 50 ms -- after 5 ms a cut is ~63 % there, a boost ~10 % */
  static float key[2 * 240], gain[240];
  udsp_udcomp_reset(&c);
  for (int i = 0; i < 2 * 240; i++) {
    key[i] = powf(10.0f, -10.0f / 20.0f);
  }
  udsp_udcomp_gain(&c, key, gain, 240);
  const float cut = -20.0f * log10f(gain[239]) / 7.5f;
  udsp_udcomp_reset(&c);
  for (int i = 0; i < 2 * 240; i++) {
    key[i] = powf(10.0f, -46.0f / 20.0f);
  }
  udsp_udcomp_gain(&c, key, gain, 240);
  const float boost = 20.0f * log10f(gain[239]) / 4.5f;
  UDSP_TEST_ASSERT(fabsf(cut - 0.632f) < 0.01f && fabsf(boost - 0.095f) < 0.01f,
                   "5 ms in: the cut %.3f there, the boost %.3f", (double)cut,
                   (double)boost);
}

/* the power detector (Kilohearts' Dynamics'): the two channels' mean
 * square -- one alone reads 3 dB down --, over its window (2 ms here: a
 * tau in, 1 - 1/e of the step's mean square), then the level smoothed
 * (attack 0: at once), the curve's gain as it comes; its floor holding
 * the level, an expansion's cut coming in as the level falls */
static void test_udcomp_power(void) {
  UDSP_UDCOMP c;
  udsp_udcomp_set(&c, -20.0f, 4.0f, -80.0f, 1.0f, 0.0f, -100.0f, 1e-6f, 50.0f,
                  UDSP_UDCOMP_POWER, RATE);
  udsp_udcomp_set_window(&c, 2.0f, RATE);
  static float key[2 * RATE], gain[RATE];
  const float a = powf(10.0f, -10.0f / 20.0f);
  for (int one = 0; one < 2; one++) {
    for (int i = 0; i < RATE; i++) {
      key[2 * i] = a;
      key[2 * i + 1] = one ? 0.0f : a;
    }
    udsp_udcomp_reset(&c);
    udsp_udcomp_gain(&c, key, gain, RATE);
    const float g = 20.0f * log10f(gain[RATE - 1]);
    const float want = one ? -0.75f * (10.0f - 3.0103f) : -7.5f;
    UDSP_TEST_ASSERT(
        fabsf(g - want) < 0.01f, "%s at -10 dB, 4:1 over -20: %.3f dB (%.3f)",
        one ? "the left alone" : "both channels", (double)g, (double)want);
  }
  udsp_udcomp_reset(&c);
  for (int i = 0; i < 2 * 96; i++) {
    key[i] = a;
  }
  udsp_udcomp_gain(&c, key, gain, 96);
  const float at = 20.0f * log10f(gain[95]);
  const float at_want = -0.75f * (10.0f + 10.0f * log10f(1.0f - expf(-1.0f)));
  UDSP_TEST_ASSERT(fabsf(at - at_want) < 0.05f,
                   "its window's tau in: %.3f dB (%.3f)", (double)at,
                   (double)at_want);
  /* the floor holds: 1:2 upward under -40, a floor at -60 -- a key at
   * -80 lifted as one at -60 is, +10 dB (not none) */
  udsp_udcomp_set_slopes(&c, 0.0f, 1.0f, -40.0f, 0.5f, 0.0f, 120.0f, -60.0f,
                         1e-6f, 1e-6f, UDSP_UDCOMP_POWER, RATE);
  udsp_udcomp_set_window(&c, 1.0f, RATE);
  const float held = ud_settled(&c, -80.0f);
  /* an expansion (slope 3 under -30) after a drop from -10 to -50: the
   * level falls at the release's pace (20 ms), so 1 ms on nothing's cut
   * yet; settled, 2 x 20 dB */
  udsp_udcomp_set_slopes(&c, 0.0f, 1.0f, -30.0f, 3.0f, 0.0f, 120.0f, -100.0f,
                         1e-3f, 20.0f, UDSP_UDCOMP_POWER, RATE);
  udsp_udcomp_set_window(&c, 0.1f, RATE);
  const float hi = powf(10.0f, -10.0f / 20.0f),
              lo = powf(10.0f, -50.0f / 20.0f);
  for (int i = 0; i < RATE; i++) {
    key[2 * i] = key[2 * i + 1] = i < RATE / 2 ? hi : lo;
  }
  udsp_udcomp_reset(&c);
  udsp_udcomp_gain(&c, key, gain, RATE);
  const float soon = 20.0f * log10f(gain[RATE / 2 + 48]);
  const float late = 20.0f * log10f(gain[RATE - 1]);
  UDSP_TEST_ASSERT(fabsf(held - 10.0f) < 0.05f && fabsf(soon) < 0.05f &&
                       fabsf(late + 40.0f) < 0.1f,
                   "held at the floor: %+.2f dB; an expansion 1 ms after the "
                   "drop %.2f, settled %.2f",
                   (double)held, (double)soon, (double)late);
}

/* by slopes: compressing above, expanding below, the expansion's depth
 * held, a knee across both, the curve smooth through each */
static void test_udcomp_slopes(void) {
  UDSP_UDCOMP c;
  udsp_udcomp_set_slopes(&c, -20.0f, 0.5f, -40.0f, 4.0f, 0.0f, 12.0f, -80.0f,
                         5.0f, 50.0f, 0, RATE);
  const float above = ud_settled(&c, -10.0f), below = ud_settled(&c, -46.0f);
  const float deep = ud_settled(&c, -95.0f);
  UDSP_TEST_ASSERT(fabsf(above + 5.0f) < 0.01f && fabsf(below + 18.0f) < 0.01f,
                   "10 dB over at a slope of 0.5: -5 dB (%.3f); 6 dB under at "
                   "4: -18 dB (%.3f)",
                   (double)above, (double)below);
  UDSP_TEST_ASSERT(
      fabsf(udsp_udcomp_curve(&c, -95.0f) + 95.0f + 120.0f) < 1e-3f &&
          deep < -100.0f,
      "under the floor, expanding: held at -120 dB (%.1f)", (double)deep);
  /* upward compression to the full: slope 0 pins the level */
  udsp_udcomp_set_slopes(&c, 0.0f, 1.0f, -30.0f, 0.0f, 0.0f, 30.0f, -80.0f,
                         5.0f, 50.0f, 0, RATE);
  UDSP_TEST_ASSERT(fabsf(udsp_udcomp_curve(&c, -50.0f) + 30.0f) < 1e-4f,
                   "slope 0 below: -50 -> -30");
  /* expanding upward above: a boost, capped */
  udsp_udcomp_set_slopes(&c, -20.0f, 2.0f, -40.0f, 1.0f, 0.0f, 6.0f, -80.0f,
                         5.0f, 50.0f, 0, RATE);
  UDSP_TEST_ASSERT(fabsf(udsp_udcomp_curve(&c, -17.0f) + 14.0f) < 1e-4f &&
                       fabsf(udsp_udcomp_curve(&c, -5.0f) + 5.0f - 6.0f) <
                           1e-4f,
                   "slope 2 above: -17 -> -14, -5 at the 6 dB cap");
  /* a 10 dB knee: at a threshold a quarter of the line's 5 dB, the curve
   * continuous with the lines at the knee's ends and its slope too */
  udsp_udcomp_set_slopes(&c, -20.0f, 0.0f, -40.0f, 3.0f, 10.0f, 30.0f, -100.0f,
                         5.0f, 50.0f, 0, RATE);
  const float at = udsp_udcomp_curve(&c, -20.0f) + 20.0f;
  const float at_lo = udsp_udcomp_curve(&c, -40.0f) + 40.0f;
  UDSP_TEST_ASSERT(fabsf(at + 1.25f) < 1e-4f && fabsf(at_lo + 2.5f) < 1e-4f,
                   "the knees at the thresholds: -1.25 (%.4f), -2.5 (%.4f)",
                   (double)at, (double)at_lo);
  float worst = 0.0f;
  for (int i = 0; i < 6000; i++) {
    const float x = -60.0f + 0.01f * (float)i;
    const float d = udsp_udcomp_curve(&c, x + 0.01f) - udsp_udcomp_curve(&c, x);
    const float d2 =
        udsp_udcomp_curve(&c, x + 0.02f) - udsp_udcomp_curve(&c, x + 0.01f);
    worst = fabsf(d2 - d) > worst ? fabsf(d2 - d) : worst;
  }
  UDSP_TEST_ASSERT(worst < 1e-3f, "no corner along the curve (%.5f)",
                   (double)worst);
  /* knees overlapping (thresholds 4 dB apart, a 12 dB knee): both add */
  udsp_udcomp_set_slopes(&c, -20.0f, 0.25f, -24.0f, 0.5f, 12.0f, 30.0f, -100.0f,
                         5.0f, 50.0f, 0, RATE);
  const float mid = udsp_udcomp_curve(&c, -22.0f) + 22.0f;
  const float want = -0.75f * 4.0f * 4.0f / 24.0f + 0.5f * 4.0f * 4.0f / 24.0f;
  UDSP_TEST_ASSERT(fabsf(mid - want) < 1e-4f, "overlapping: %.4f (%.4f)",
                   (double)mid, (double)want);
  /* the settled gain is the curve's, through the vector path too */
  int same = 1;
  for (int i = 0; i < 10; i++) {
    const float lvl = -70.0f + 7.3f * (float)i;
    same &=
        fabsf(ud_settled(&c, lvl) - (udsp_udcomp_curve(&c, lvl) - lvl)) < 0.01f;
  }
  UDSP_TEST_ASSERT(same, "settled gains on the curve");
}

/* side by side, the same bits as one at a time: a peak and two RMS
 * detectors, their own times, over blocks (their state carried) */
static void test_udcomp_gains(void) {
  enum { N = 3, B = 333, BLOCKS = 4 };
  UDSP_UDCOMP one[N], side[N];
  for (int j = 0; j < N; j++) {
    udsp_udcomp_reset(&one[j]);
    udsp_udcomp_set(&one[j], -20.0f + 4.0f * j, 4.0f + j, -40.0f, 3.0f, 12.0f,
                    -80.0f, 1.0f + 7.0f * j, 30.0f + 40.0f * j, j, RATE);
    side[j] = one[j];
  }
  static float key[N][2 * B], a[N][B], b[N][B];
  uint32_t seed = 1;
  int same = 1;
  for (int blk = 0; blk < BLOCKS; blk++) {
    for (int j = 0; j < N; j++) {
      for (int i = 0; i < 2 * B; i++) {
        seed = seed * 1664525u + 1013904223u;
        /* bursts and quiet: every branch of the curve */
        const float amp = (blk + i / 100) % 3 ? 0.9f : 0.002f;
        key[j][i] = amp * ((float)(seed >> 8) / 8388608.0f - 1.0f);
      }
      udsp_udcomp_gain(&one[j], key[j], a[j], B);
    }
    UDSP_UDCOMP *c[N] = {&side[0], &side[1], &side[2]};
    const float *k[N] = {key[0], key[1], key[2]};
    float *g[N] = {b[0], b[1], b[2]};
    udsp_udcomp_gains(c, k, g, N, B);
    same &= !memcmp(a, b, sizeof(a));
    for (int j = 0; j < N; j++) {
      same &= !memcmp(&one[j], &side[j], sizeof(one[j]));
    }
  }
  UDSP_TEST_ASSERT(same, "three side by side (peak, RMS, power): the bits "
                         "of one at a time");
}

/* the vintage compressor: n frames of a sine at hz (its rms level_db) or,
 * hz 0, a square at the half rate (|x| steady); the reduction after them,
 * and the first frame cut (-1: none) */
static float vcomp_run(UDSP_VCOMP *c, UDSP_EQ_BAND *const *key, int keys,
                       float level_db, float hz, size_t n, float *first_cut) {
  const float a = udsp_gain_from_db(level_db);
  float buf[2 * 64];
  size_t done = 0;
  if (first_cut) {
    *first_cut = -1.0f;
  }
  while (done < n) {
    size_t m = n - done < 64 ? n - done : 64;
    for (size_t i = 0; i < m; i++) {
      float x = hz > 0.0f ? sinf(6.2831853f * hz * (float)(done + i) / RATE) *
                                1.41421356f * a
                          : ((done + i) & 1 ? -a : a);
      buf[2 * i] = buf[2 * i + 1] = x;
    }
    udsp_vcomp_stereo(c, key, keys, buf, m);
    for (size_t i = 0; first_cut && hz <= 0.0f && i < m; i++) {
      if (*first_cut < 0.0f && fabsf(buf[2 * i]) < a * 0.9999f) {
        *first_cut = (float)(done + i);
      }
    }
    done += m;
  }
  return c->reduction;
}

static void test_vcomp(void) {
  UDSP_VCOMP c;
  udsp_vcomp_reset(&c);
  udsp_vcomp_set(&c, -20.0f, 4.0f, 0.0f, 0.5f, 100000.0f, RATE);
  udsp_vcomp_set_auto(&c, 0.0f, RATE);
  /* 12 dB over at 4:1: 9 dB down, the static curve's, though the
   * detector reads the output (3 dB over at a slope of 3) */
  float red = vcomp_run(&c, NULL, 0, -8.0f, 0.0f, RATE, NULL);
  UDSP_TEST_ASSERT(fabsf(red - 9.0f) < 0.02f,
                   "feedback, 4:1, 12 dB over: 9 dB down (%.3f)", (double)red);
  UDSP_TEST_ASSERT(fabsf(udsp_vcomp_curve(&c, -17.0f) - 9.0f) < 1e-4f,
                   "the curve on the output's level: -17 reads 9 dB");

  /* the leak: a steady level reads leak / (attack + leak) of itself --
   * attack 10, leak 90: 0.915 dB short, 0.75 of it less reduction */
  udsp_vcomp_reset(&c);
  udsp_vcomp_set(&c, -20.0f, 4.0f, 0.0f, 10.0f, 90.0f, RATE);
  red = vcomp_run(&c, NULL, 0, -8.0f, 0.0f, 2 * RATE, NULL);
  const float want = 0.75f * (12.0f + 20.0f * log10f(0.9f));
  UDSP_TEST_ASSERT(fabsf(red - want) < 0.03f,
                   "the leak's shortfall: %.3f dB down (%.3f)", (double)want,
                   (double)red);

  /* the dead time: nothing until the detector has climbed to the
   * threshold -- from rest toward 0.9 x at 1 / 10 + 1 / 90 per ms, -20
   * of -8 dB (0.251 of it) after -ln(1 - 0.251 / 0.9) x 9 ms = 2.94 ms */
  udsp_vcomp_reset(&c);
  float cut;
  vcomp_run(&c, NULL, 0, -8.0f, 0.0f, RATE / 10, &cut);
  UDSP_TEST_ASSERT(fabsf(cut / RATE * 1000.0f - 2.94f) < 0.1f,
                   "the first frame cut %.2f ms in (2.94)",
                   (double)(cut / RATE * 1000.0f));

  /* auto gain: a steady level comes out as loud as it went in */
  udsp_vcomp_reset(&c);
  udsp_vcomp_set_auto(&c, 500.0f, RATE);
  red = vcomp_run(&c, NULL, 0, -8.0f, 0.0f, 5 * RATE, NULL);
  const float makeup = udsp_db_from_gain(c.in_m / c.out_m);
  UDSP_TEST_ASSERT(red > 8.0f && fabsf(makeup - red) < 0.02f,
                   "auto gain: the makeup the reduction (%.3f vs %.3f)",
                   (double)makeup, (double)red);

  /* idle: the same as running on silence */
  UDSP_VCOMP d = c;
  float sil[2 * 4800] = {0};
  udsp_vcomp_stereo(&c, NULL, 0, sil, 4800);
  udsp_vcomp_idle(&d, 4800);
  UDSP_TEST_ASSERT(fabsf(c.env / d.env - 1.0f) < 1e-3f &&
                       fabsf(c.in_m / d.in_m - 1.0f) < 1e-3f &&
                       fabsf(c.reduction - d.reduction) < 1e-3f,
                   "idle 100 ms: the detector and meters as on silence "
                   "(env %.4g / %.4g)",
                   (double)c.env, (double)d.env);

  /* the key: a high pass at 2 kHz hears a 200 Hz sine 40 dB down: no
   * reduction; the same sine without it, 12 dB over at its peak */
  UDSP_EQ_BAND hp;
  udsp_eq_reset(&hp);
  udsp_eq_set(&hp, UDSP_EQ_LOW_CUT_12, 2000.0f, 0.7071f, 0.0f, RATE);
  UDSP_EQ_BAND *const key[] = {&hp};
  udsp_vcomp_reset(&c);
  udsp_vcomp_set(&c, -20.0f, 4.0f, 0.0f, 0.5f, 100000.0f, RATE);
  udsp_vcomp_set_auto(&c, 0.0f, RATE);
  red = vcomp_run(&c, key, 1, -11.0f, 200.0f, RATE / 2, NULL);
  UDSP_TEST_ASSERT(red == 0.0f, "a key high pass: a low sine unheard (%.3f)",
                   (double)red);
  udsp_vcomp_reset(&c);
  red = vcomp_run(&c, NULL, 0, -11.0f, 200.0f, RATE / 2, NULL);
  UDSP_TEST_ASSERT(red > 8.5f, "without it: compressed (%.3f)", (double)red);
}

int main(void) {
  printf("dyn_test\n");
  test_log_exp();
  test_comp();
  test_comp_second();
  test_gate();
  test_limiter();
  test_udcomp();
  test_udcomp_power();
  test_udcomp_slopes();
  test_udcomp_gains();
  test_vcomp();
  return udsp_test_summary();
}
