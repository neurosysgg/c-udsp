#include "udsp_test.h"

#include <math.h>
#include <stdint.h>
#include <string.h>

#include "eq.h"

#define RATE 48000

/* a sine at hz through the band, both channels: its amplitude once it
 * settles (left, from its RMS), and whether right came out the same */
static float gain_at(UDSP_EQ_SHAPE shape, float fc, float q, float db,
                     float hz, int *same) {
  UDSP_EQ_BAND b;
  udsp_eq_reset(&b);
  udsp_eq_set(&b, shape, fc, q, db, RATE);
  static float x[RATE * 2];
  for (int i = 0; i < RATE; i++) {
    x[2 * i] = x[2 * i + 1] = (float)sin(6.283185307179586 * hz * i / RATE);
  }
  udsp_eq_stereo(&b, x, RATE);
  /* RMS, not the peak: a few samples a period miss a sine's top */
  double sq = 0.0;
  *same = 1;
  for (int i = RATE / 2; i < RATE; i++) {
    sq += (double)x[2 * i] * x[2 * i];
    *same &= x[2 * i] == x[2 * i + 1];
  }
  return (float)sqrt(2.0 * sq / (RATE / 2));
}

static float db(float g) { return 20.0f * log10f(g); }

/* the measured gain matches the drawn curve, within a hair */
static void check(UDSP_EQ_SHAPE shape, float fc, float q, float gdb, float hz,
                  float want_db, float tol, const char *what) {
  int same;
  float got = db(gain_at(shape, fc, q, gdb, hz, &same));
  UDSP_EQ_BAND b;
  udsp_eq_reset(&b);
  udsp_eq_set(&b, shape, fc, q, gdb, RATE);
  float drawn = db(udsp_eq_response(&b, hz, RATE));
  UDSP_TEST_ASSERT(fabsf(got - want_db) < tol && fabsf(drawn - got) < 0.05f &&
                       same,
                   "%s: %.2f dB (curve %.2f, want %.2f)", what, (double)got,
                   (double)drawn, (double)want_db);
}

static void test_shapes(void) {
  check(UDSP_EQ_BELL, 1000.0f, 1.0f, 6.0f, 1000.0f, 6.0f, 0.02f,
        "bell +6 at its centre");
  check(UDSP_EQ_BELL, 1000.0f, 1.0f, -12.0f, 1000.0f, -12.0f, 0.02f,
        "bell -12 at its centre");
  check(UDSP_EQ_BELL, 1000.0f, 1.0f, 6.0f, 100.0f, 0.0f, 0.1f,
        "bell: flat far below");
  check(UDSP_EQ_BELL, 1000.0f, 1.0f, 0.0f, 3000.0f, 0.0f, 1e-3f,
        "bell at 0 dB: flat");
  check(UDSP_EQ_LOW_SHELF, 200.0f, 0.71f, 9.0f, 20.0f, 9.0f, 0.1f,
        "low shelf +9 below");
  check(UDSP_EQ_LOW_SHELF, 200.0f, 0.71f, 9.0f, 200.0f, 4.5f, 0.05f,
        "low shelf: half at its corner");
  check(UDSP_EQ_LOW_SHELF, 200.0f, 0.71f, 9.0f, 8000.0f, 0.0f, 0.05f,
        "low shelf: flat above");
  check(UDSP_EQ_HIGH_SHELF, 4000.0f, 0.71f, -6.0f, 18000.0f, -6.0f, 0.15f,
        "high shelf -6 above");
  check(UDSP_EQ_HIGH_SHELF, 4000.0f, 0.71f, -6.0f, 100.0f, 0.0f, 0.05f,
        "high shelf: flat below");
  check(UDSP_EQ_LOW_CUT_12, 1000.0f, 0.7071f, 0.0f, 1000.0f, -3.01f, 0.05f,
        "low cut 12: -3 dB at its corner");
  check(UDSP_EQ_LOW_CUT_12, 1000.0f, 0.7071f, 0.0f, 250.0f, -24.1f, 0.3f,
        "low cut 12: two octaves under, -24");
  check(UDSP_EQ_HIGH_CUT_12, 1000.0f, 0.7071f, 0.0f, 1000.0f, -3.01f, 0.05f,
        "high cut 12: -3 dB at its corner");
  check(UDSP_EQ_LOW_CUT_48, 1000.0f, 0.7071f, 0.0f, 1000.0f, -3.01f, 0.05f,
        "low cut 48: a Butterworth's -3 at its corner");
  check(UDSP_EQ_LOW_CUT_48, 1000.0f, 0.7071f, 0.0f, 500.0f, -48.2f, 0.3f,
        "low cut 48: an octave under, -48");
  check(UDSP_EQ_HIGH_CUT_48, 1000.0f, 0.7071f, 0.0f, 100.0f, 0.0f, 0.01f,
        "high cut 48: flat in its pass band");
  check(UDSP_EQ_HIGH_CUT_48, 1000.0f, 0.7071f, 0.0f, 2000.0f, -48.4f, 0.5f,
        "high cut 48: an octave over, -48");
  check(UDSP_EQ_LOW_CUT_12, 1000.0f, 4.0f, 0.0f, 1000.0f, 12.04f, 0.05f,
        "low cut 12 at Q 4: +12 at its corner");
  int same;
  UDSP_TEST_ASSERT(gain_at(UDSP_EQ_NOTCH, 1000.0f, 2.0f, 0.0f, 1000.0f,
                           &same) < 1e-3f,
                   "the notch takes its frequency out");
}

/* the curve itself, without signals */
static void test_response(void) {
  UDSP_EQ_BAND b;
  udsp_eq_set(&b, UDSP_EQ_BELL, 2000.0f, 2.0f, -9.0f, RATE);
  UDSP_TEST_ASSERT(fabsf(db(udsp_eq_response(&b, 2000.0f, RATE)) + 9.0f) <
                       1e-3f,
                   "the bell's curve at its centre");
  udsp_eq_set(&b, UDSP_EQ_HIGH_CUT_48, 1000.0f, 0.7071f, 0.0f, RATE);
  UDSP_TEST_ASSERT(udsp_eq_response(&b, 23990.0f, RATE) < 1e-6f,
                   "a high cut is closed at Nyquist (the bilinear zero)");
  UDSP_TEST_ASSERT(fabsf(udsp_eq_response(&b, 0.0f, RATE) - 1.0f) < 1e-6f,
                   "and open at DC");
}

/* retuned every 16 frames across its range at full Q: no blow-up */
static void test_swept_stable(void) {
  UDSP_EQ_BAND b;
  udsp_eq_reset(&b);
  uint32_t rng = 7;
  float x[32], worst = 0.0f;
  for (int i = 0; i < 40000; i++) {
    float hz = (i * 37) % 2 ? 23000.0f : 10.0f + (float)(i % 997) * 22.0f;
    udsp_eq_set(&b, (UDSP_EQ_SHAPE)(i / 5000), hz, 40.0f,
                (i % 3) ? 30.0f : -30.0f, RATE);
    for (int k = 0; k < 32; k++) {
      rng ^= rng << 13;
      rng ^= rng >> 17;
      rng ^= rng << 5;
      x[k] = (float)(int32_t)rng * (1.0f / 2147483648.0f);
    }
    udsp_eq_stereo(&b, x, 16);
    for (int k = 0; k < 32; k++) {
      worst = fabsf(x[k]) > worst ? fabsf(x[k]) : worst;
    }
  }
  UDSP_TEST_ASSERT(isfinite(worst) && worst < 1e4f,
                   "swept at Q 40, +-30 dB: bounded (%g)", (double)worst);
}

/* the same bits whatever the block cut */
static void test_blocks(void) {
  static float a[4096], b2[4096];
  uint32_t rng = 3;
  for (int i = 0; i < 4096; i++) {
    rng ^= rng << 13;
    rng ^= rng >> 17;
    rng ^= rng << 5;
    a[i] = b2[i] = (float)(int32_t)rng * (1.0f / 2147483648.0f);
  }
  UDSP_EQ_BAND x, y;
  udsp_eq_reset(&x);
  udsp_eq_reset(&y);
  udsp_eq_set(&x, UDSP_EQ_LOW_CUT_48, 300.0f, 1.0f, 0.0f, RATE);
  udsp_eq_set(&y, UDSP_EQ_LOW_CUT_48, 300.0f, 1.0f, 0.0f, RATE);
  udsp_eq_stereo(&x, a, 2048);
  for (int i = 0; i < 2048; i += 7) {
    udsp_eq_stereo(&y, b2 + 2 * i, i + 7 <= 2048 ? 7 : 2048 - i);
  }
  int same = 1;
  for (int i = 0; i < 4096; i++) {
    same &= a[i] == b2[i];
  }
  UDSP_TEST_ASSERT(same, "a block of 2048 = blocks of 7, to the bit");
}

/* a 12 dB boost ramped in over 480 frames, in uneven calls: no frame
 * steps further than the boosted sine moves on its own, it ends on the
 * target's coefficients exactly, and the same bits as one call */
static void test_ramp(void) {
  static float x[2 * 4800], y[2 * 4800];
  for (int i = 0; i < 4800; i++) {
    x[2 * i] = x[2 * i + 1] = 0.25f * (float)sin(6.283185307179586 * 100.0 *
                                                 i / RATE);
  }
  memcpy(y, x, sizeof(x));
  UDSP_EQ_BAND a, b, to;
  udsp_eq_reset(&a);
  udsp_eq_set(&a, UDSP_EQ_BELL, 100.0f, 1.0f, 0.0f, RATE);
  b = a;
  udsp_eq_set(&to, UDSP_EQ_BELL, 100.0f, 1.0f, 12.0f, RATE);
  udsp_eq_stereo(&a, x, 1000);
  udsp_eq_ramp(&a, &to, 480);
  udsp_eq_stereo(&a, x + 2000, 3800);
  udsp_eq_stereo(&b, y, 1000);
  udsp_eq_ramp(&b, &to, 480);
  for (int i = 1000; i < 4800; i += 61) {
    udsp_eq_stereo(&b, y + 2 * i, i + 61 <= 4800 ? 61 : 4800 - i);
  }
  float worst = 0.0f;
  int same = 1;
  for (int i = 1; i < 4800; i++) {
    float d = fabsf(x[2 * i] - x[2 * i - 2]);
    worst = d > worst ? d : worst;
    same &= x[2 * i] == y[2 * i];
  }
  /* 100 Hz at 0.25 moves 0.0033 a frame; 4x that boosted */
  UDSP_TEST_ASSERT(worst < 0.0033f * 4.0f, "ramped: no steps (%.4f a frame)",
                   (double)worst);
  UDSP_TEST_ASSERT(a.ramp == 0 && a.st[0].m1 == to.st[0].m1 &&
                       a.st[0].a1 == to.st[0].a1 && a.st[0].k == to.st[0].k,
                   "ends on the target's coefficients exactly");
  UDSP_TEST_ASSERT(same, "the same bits in calls of 61 frames");
  udsp_eq_set(&to, UDSP_EQ_LOW_CUT_48, 100.0f, 1.0f, 0.0f, RATE);
  udsp_eq_ramp(&a, &to, 480);
  UDSP_TEST_ASSERT(a.ramp == 0 && a.stages == 4,
                   "a shape with more sections jumps");
}

int main(void) {
  printf("eq_test\n");
  test_shapes();
  test_response();
  test_swept_stable();
  test_blocks();
  test_ramp();
  return udsp_test_summary();
}
