#include "udsp_test.h"

#include <math.h>
#include <stdint.h>
#include <string.h>

#include "fx/eq.h"

#define RATE 48000

/* a sine at hz through the band, both channels: its amplitude once it
 * settles (left, from its RMS), and whether right came out the same */
static float gain_at(UDSP_EQ_SHAPE shape, float fc, float q, float db, float hz,
                     int *same) {
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
  check(UDSP_EQ_BAND_PASS, 1000.0f, 1.0f, 0.0f, 1000.0f, 0.0f, 0.02f,
        "band pass: 0 dB at its centre");
  check(UDSP_EQ_BAND_PASS, 1000.0f, 1.0f, 0.0f, 2000.0f, -5.12f, 0.1f,
        "band pass at Q 1: an octave over, -5.1 (1 / sqrt(1 + 1.5^2))");
  check(UDSP_EQ_BAND_PASS, 1000.0f, 1.0f, 6.0f, 62.5f, -24.07f, 0.1f,
        "band pass: four octaves under, -24.07; the gain ignored");
  int same;
  UDSP_TEST_ASSERT(gain_at(UDSP_EQ_NOTCH, 1000.0f, 2.0f, 0.0f, 1000.0f, &same) <
                       1e-3f,
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
    x[2 * i] = x[2 * i + 1] =
        0.25f * (float)sin(6.283185307179586 * 100.0 * i / RATE);
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

static uint32_t g_rng = 11;
static float noise(void) {
  g_rng ^= g_rng << 13;
  g_rng ^= g_rng >> 17;
  g_rng ^= g_rng << 5;
  return (float)(int32_t)g_rng * (1.0f / 2147483648.0f);
}

/* a chain of n random bands (any shape, cuts of 48 among them: up to 32
 * sections), some ramping, through udsp_eq_chain() in calls of `cut`
 * frames and through udsp_eq_stereo() band by band in one: the same bits,
 * the same state after */
static int chain_same(int n, int cut, int ramps) {
  enum { F = 1500 };
  static float a[2 * F], b[2 * F];
  UDSP_EQ_BAND x[8], y[8], *px[8];
  for (int i = 0; i < n; i++) {
    const UDSP_EQ_SHAPE shape =
        (UDSP_EQ_SHAPE)((uint32_t)((noise() + 1.0f) * 1e6f) % UDSP_EQ_SHAPES);
    udsp_eq_reset(&x[i]);
    udsp_eq_set(&x[i], shape, 300.0f + 4000.0f * (noise() + 1.0f),
                0.5f + 4.0f * (noise() + 1.0f), 12.0f * noise(), RATE);
    if (ramps && i % 2) {
      UDSP_EQ_BAND to;
      udsp_eq_set(&to, shape, 2000.0f, 2.0f, -6.0f, RATE);
      udsp_eq_ramp(&x[i], &to, 100 + 37 * (unsigned)i);
    }
    y[i] = x[i];
    px[i] = &x[i];
  }
  for (int i = 0; i < 2 * F; i++) {
    a[i] = b[i] = noise();
  }
  for (int at = 0; at < F; at += cut) {
    udsp_eq_chain(px, n, a + 2 * at, (size_t)(at + cut <= F ? cut : F - at));
  }
  for (int i = 0; i < n; i++) {
    udsp_eq_stereo(&y[i], b, F);
  }
  int same = !memcmp(a, b, sizeof a);
  for (int i = 0; i < n; i++) {
    same &= !memcmp(x[i].ic1, y[i].ic1, sizeof x[i].ic1) &&
            !memcmp(x[i].ic2, y[i].ic2, sizeof x[i].ic2) &&
            x[i].ramp == y[i].ramp;
  }
  return same;
}

static void test_chain(void) {
  int same = 1;
  static const int CUT[] = {1500, 256, 32, 7, 1};
  for (int n = 1; n <= 8; n++) {
    for (int c = 0; c < 5; c++) {
      same &= chain_same(n, CUT[c], 0) && chain_same(n, CUT[c], 1);
    }
  }
  UDSP_TEST_ASSERT(same, "a chain of 1-8 bands, ramping or not, in calls of "
                         "1500, 256, 32, 7 and 1: the bands' own bits");
  float none[2] = {0.5f, -0.5f};
  udsp_eq_chain(NULL, 0, none, 1);
  UDSP_TEST_ASSERT(none[0] == 0.5f && none[1] == -0.5f, "no bands: as it was");
}

/* the matched bell: its gain at the centre, flat at DC, the cut the
 * boost's inverse; no cramping at the top; the
 * section (Simper's SVF) the biquad it was designed as */
static void test_matched(void) {
  check(UDSP_EQ_BELL_MATCHED, 1000.0f, 2.0f, 6.0f, 1000.0f, 6.0f, 0.02f,
        "matched bell +6 at its centre");
  check(UDSP_EQ_BELL_MATCHED, 50.0f, 2.0f, -12.0f, 50.0f, -12.0f, 0.03f,
        "matched bell -12 at 50 Hz");
  check(UDSP_EQ_BELL_MATCHED, 1000.0f, 2.0f, 6.0f, 20.0f, 0.0f, 0.01f,
        "matched bell: flat far below");
  double b[3], a[3], bc[3], ac[3];
  udsp_eq_bell_matched(12000.0, 2.0, 9.0, RATE, b, a);
  udsp_eq_bell_matched(12000.0, 2.0, -9.0, RATE, bc, ac);
  double worst = 0.0, svf = 0.0;
  UDSP_EQ_BAND f;
  udsp_eq_set(&f, UDSP_EQ_BELL_MATCHED, 12000.0f, 2.0f, 9.0f, RATE);
  for (int n = 0; n < 57; n++) {
    const double hz = 100.0 * pow(1.1, n);
    const double s = sin(3.141592653589793 * hz / RATE), phi = s * s;
    const double p = udsp_eq_biquad_power(b, a, phi);
    worst = fmax(worst, fabs(p * udsp_eq_biquad_power(bc, ac, phi) - 1.0));
    svf = fmax(svf, fabs(10.0 * log10(p) -
                         20.0 * log10(udsp_eq_response(&f, (float)hz, RATE))));
  }
  /* at half the rate, an octave over 12 kHz, the bilinear's bell is
   * back to 0 dB; the matched one keeps its skirt */
  const double top = 10.0 * log10(udsp_eq_biquad_power(b, a, 1.0));
  UDSP_TEST_ASSERT(top > 1.5 && top < 1.8, "half the rate: %.3f dB", top);
  UDSP_TEST_ASSERT(worst < 1e-9, "a cut is the boost's inverse (%.1e)", worst);
  UDSP_TEST_ASSERT(svf < 0.01, "the SVF section is the biquad (%.4f dB)", svf);
}

/* the matched shelves and the first-order shapes: each against its
 * analog prototype across the band, and through the filter itself */
static double proto_db(UDSP_EQ_SHAPE shape, double hz, double q, double g,
                       double f) {
  const double A = pow(10.0, g / 40.0), G = A * A, x = f / hz;
  if (shape == UDSP_EQ_HIGH_SHELF_MATCHED ||
      shape == UDSP_EQ_LOW_SHELF_MATCHED) {
    const double r = sqrt(A) / q;
    const int high = shape == UDSP_EQ_HIGH_SHELF_MATCHED;
    /* A (A s^2 + r s + 1) / (s^2 + r s + A) and its mirror */
    double num_re, num_im, den_re, den_im;
    if (high) {
      num_re = A * (1 - A * x * x), num_im = A * r * x;
      den_re = A - x * x, den_im = r * x;
    } else {
      num_re = A * (A - x * x), num_im = A * r * x;
      den_re = 1 - A * x * x, den_im = r * x;
    }
    return 10 * log10((num_re * num_re + num_im * num_im) /
                      (den_re * den_re + den_im * den_im));
  }
  if (shape == UDSP_EQ_LOW_CUT_6) {
    return 10 * log10(x * x / (1 + x * x));
  }
  if (shape == UDSP_EQ_HIGH_CUT_6) {
    return 10 * log10(1 / (1 + x * x));
  }
  /* first-order shelves: corners hz / sqrt G and hz sqrt G */
  const double sg = sqrt(G), lo = x * sg, hi = x / sg;
  return shape == UDSP_EQ_HIGH_SHELF_6
             ? 10 * log10((1 + lo * lo) / (1 + hi * hi))
             : 20 * log10(G) + 10 * log10((1 + hi * hi) / (1 + lo * lo));
}

static void test_matched_more(void) {
  static const struct {
    UDSP_EQ_SHAPE shape;
    float hz, q, g;
    double tol;
  } C[] = {
      {UDSP_EQ_HIGH_SHELF_MATCHED, 9495, 0.224f, 1.3f, 0.04},
      {UDSP_EQ_HIGH_SHELF_MATCHED, 1000, 0.7f, 6, 0.04},
      {UDSP_EQ_LOW_SHELF_MATCHED, 200, 0.7f, -9, 0.04},
      {UDSP_EQ_LOW_CUT_6, 60, 1, 0, 0.02},
      /* a one-pole matched at two points: within 1 dB at the top */
      {UDSP_EQ_HIGH_CUT_6, 2000, 1, 0, 1.0},
      {UDSP_EQ_HIGH_SHELF_6, 7655, 1, 1.44f, 0.06},
      {UDSP_EQ_LOW_SHELF_6, 300, 1, -6, 0.02},
  };
  for (size_t c = 0; c < sizeof C / sizeof *C; c++) {
    UDSP_EQ_BAND b;
    udsp_eq_reset(&b);
    udsp_eq_set(&b, C[c].shape, C[c].hz, C[c].q, C[c].g, RATE);
    double worst = 0.0;
    for (int k = 0; k < 60; k++) {
      const double f = 20.0 * pow(1000.0, k / 59.0);
      worst =
          fmax(worst, fabs(20 * log10(udsp_eq_response(&b, (float)f, RATE)) -
                           proto_db(C[c].shape, C[c].hz, C[c].q, C[c].g, f)));
    }
    int same;
    const double f = C[c].hz * 1.5;
    const double heard =
        db(gain_at(C[c].shape, C[c].hz, C[c].q, C[c].g, (float)f, &same));
    UDSP_TEST_ASSERT(worst < C[c].tol &&
                         fabs(heard - proto_db(C[c].shape, C[c].hz, C[c].q,
                                               C[c].g, f)) < C[c].tol + 0.02,
                     "shape %d at %.0f Hz: within %.3f dB of its prototype "
                     "to 20 kHz, a sine through it %.3f dB",
                     (int)C[c].shape, (double)C[c].hz, worst, heard);
  }
}

int main(void) {
  printf("eq_test\n");
  test_matched();
  test_matched_more();
  test_shapes();
  test_response();
  test_swept_stable();
  test_blocks();
  test_ramp();
  test_chain();
  return udsp_test_summary();
}
