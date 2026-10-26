#include "udsp_test.h"

#include <math.h>
#include <stdint.h>
#include <string.h>

#include "fx/geq.h"

#define RATE 48000

static UDSP_GEQ geq;

/* Kilohearts' Carve EQ's template in the user's projects (its 31 gains) */
static const float TEMPLATE[UDSP_GEQ_BANDS] = {
    -1.052535f, -0.291816f, 0.449709f,  0.516428f,  0.615663f,  0.567649f,
    0.323983f,  -0.532094f, -1.966427f, -1.836538f, -1.414258f, -1.105063f,
    -0.70698f,  -0.226499f, -0.030785f, 0.105304f,  0.128676f,  0.113437f,
    -0.032287f, -0.154877f, -0.649195f, -1.090956f, -1.34372f,  -1.185577f,
    -0.836268f, -0.534556f, -0.074605f, 0.156441f,  0.257428f,  0.276656f,
    0.238564f};

/* the solved curve passes through what was drawn, at the centres */
static void test_centres(void) {
  float t[UDSP_GEQ_BANDS], g[UDSP_GEQ_BANDS];
  uint32_t rng = 5;
  for (int i = 0; i < UDSP_GEQ_BANDS; i++) {
    rng = rng * 1664525u + 1013904223u;
    t[i] = (float)(rng >> 8) / (float)(1u << 24) * 12.0f - 6.0f;
  }
  udsp_geq_solve(&geq, t, g);
  double worst = 0.0;
  for (int i = 0; i < UDSP_GEQ_BANDS; i++) {
    worst = fmax(worst, fabs(udsp_geq_response(&geq, g, geq.hz[i]) - t[i]));
  }
  UDSP_TEST_ASSERT(worst < 0.01, "+-6 dB at random: centres within %.4f dB",
                   worst);
  memset(t, 0, sizeof t);
  t[0] = 18.0f;
  t[10] = -40.0f;
  udsp_geq_solve(&geq, t, g);
  const double a = udsp_geq_response(&geq, g, geq.hz[0]),
               b = udsp_geq_response(&geq, g, geq.hz[10]);
  UDSP_TEST_ASSERT(fabs(a - 16.0) < 0.05 && fabs(b + 16.0) < 0.05,
                   "targets held to +-16 dB: %.2f, %.2f", a, b);
}

/* against the plugin: one band at +6 (band 15, 630 Hz) a sixth of an
 * octave at a time either side, and the template at five points --
 * Carve EQ's own responses (vst3host, 48 kHz) */
static void test_carve(void) {
  float t[UDSP_GEQ_BANDS] = {0}, g[UDSP_GEQ_BANDS];
  t[15] = 6.0f;
  udsp_geq_solve(&geq, t, g);
  static const double ONE[9] = {-0.0005, -0.4896, 0.0039,  3.2642, 6.0056,
                                3.2642,  0.0039,  -0.4896, -0.0005};
  double worst = 0.0;
  for (int k = 0; k < 9; k++) {
    const double hz = geq.hz[15] * pow(2.0, (k - 4) / 6.0);
    worst = fmax(worst, fabs(udsp_geq_response(&geq, g, hz) - ONE[k]));
  }
  UDSP_TEST_ASSERT(worst < 0.02, "one band +6: within %.4f dB of Carve EQ",
                   worst);
  udsp_geq_solve(&geq, TEMPLATE, g);
  static const double HZ[5] = {50, 200, 1000, 4000, 12000},
                      TPL[5] = {0.6193, -1.4033, 0.1126, -1.1872, 0.2398};
  worst = 0.0;
  for (int k = 0; k < 5; k++) {
    worst = fmax(worst, fabs(udsp_geq_response(&geq, g, HZ[k]) - TPL[k]));
  }
  UDSP_TEST_ASSERT(worst < 0.02, "the template: within %.4f dB of Carve EQ",
                   worst);
}

/* the bells as filters: an impulse through the 31 in series, its
 * spectrum the curve */
static void test_filter(void) {
  float g[UDSP_GEQ_BANDS];
  udsp_geq_solve(&geq, TEMPLATE, g);
  static UDSP_EQ_BAND band[UDSP_GEQ_BANDS];
  UDSP_EQ_BAND *chain[UDSP_GEQ_BANDS];
  for (int i = 0; i < UDSP_GEQ_BANDS; i++) {
    udsp_eq_reset(&band[i]);
    udsp_geq_band(&geq, i, g[i], &band[i]);
    chain[i] = &band[i];
  }
  enum { N = 1 << 15 };
  static float x[2 * N];
  memset(x, 0, sizeof x);
  x[0] = x[1] = 1.0f;
  udsp_eq_chain(chain, UDSP_GEQ_BANDS, x, N);
  double worst = 0.0;
  for (int k = 0; k < 20; k++) {
    const double hz = 40.0 * pow(1.37, k);
    double re = 0.0, im = 0.0;
    for (int n = 0; n < N; n++) {
      re += x[2 * n] * cos(6.283185307179586 * hz * n / RATE);
      im -= x[2 * n] * sin(6.283185307179586 * hz * n / RATE);
    }
    const double got = 10.0 * log10(re * re + im * im);
    worst = fmax(worst, fabs(got - udsp_geq_response(&geq, g, hz)));
  }
  UDSP_TEST_ASSERT(worst < 0.01,
                   "the bells' impulse response is the curve "
                   "(within %.4f dB)",
                   worst);
}

/* at 32 kHz the bands past 0.49 of the rate (16 and 20 kHz) are left
 * out, their gains 0 */
static void test_rates(void) {
  static UDSP_GEQ low;
  udsp_geq_init(&low, 32000);
  float t[UDSP_GEQ_BANDS], g[UDSP_GEQ_BANDS];
  for (int i = 0; i < UDSP_GEQ_BANDS; i++) {
    t[i] = 3.0f;
  }
  udsp_geq_solve(&low, t, g);
  UDSP_TEST_ASSERT(low.n == 29 && g[29] == 0.0f && g[30] == 0.0f &&
                       fabs(udsp_geq_response(&low, g, 1000.0) - 3.0) < 0.01,
                   "32 kHz: %d bands, the top two 0, 1 kHz %.3f dB", low.n,
                   udsp_geq_response(&low, g, 1000.0));
}

int main(void) {
  printf("geq_test\n");
  udsp_geq_init(&geq, RATE);
  test_centres();
  test_carve();
  test_filter();
  test_rates();
  return udsp_test_summary();
}
