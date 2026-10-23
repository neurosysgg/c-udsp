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

/* side by side, the same bits as one at a time: a peak and two RMS
 * detectors, their own times, over blocks (their state carried) */
static void test_udcomp_gains(void) {
  enum { N = 3, B = 333, BLOCKS = 4 };
  UDSP_UDCOMP one[N], side[N];
  for (int j = 0; j < N; j++) {
    udsp_udcomp_reset(&one[j]);
    udsp_udcomp_set(&one[j], -20.0f + 4.0f * j, 4.0f + j, -40.0f, 3.0f, 12.0f,
                    -80.0f, 1.0f + 7.0f * j, 30.0f + 40.0f * j, j > 0, RATE);
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
  UDSP_TEST_ASSERT(same, "three side by side: the bits of one at a time");
}

int main(void) {
  printf("dyn_test\n");
  test_log_exp();
  test_comp();
  test_gate();
  test_limiter();
  test_udcomp();
  test_udcomp_gains();
  return udsp_test_summary();
}
