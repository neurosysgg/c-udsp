#include "udsp_test.h"

#include <math.h>
#include <stdint.h>
#include <string.h>

#include "fx/zdf.h"

#define RATE 48000

/* a sine at hz through the filter: its amplitude after it settles */
static float gain_at(UDSP_ZDF_MODE mode, float cutoff, float res, float hz) {
  UDSP_ZDF f;
  udsp_zdf_reset(&f);
  udsp_zdf_set(&f, cutoff, res, mode, RATE);
  static float x[RATE];
  for (int i = 0; i < RATE; i++) {
    x[i] = (float)sin(6.283185307179586 * hz * i / RATE);
  }
  udsp_zdf_mono(&f, x, RATE);
  float peak = 0.0f;
  for (int i = RATE / 2; i < RATE; i++) {
    peak = fabsf(x[i]) > peak ? fabsf(x[i]) : peak;
  }
  return peak;
}

static void test_dc(void) {
  float x[4800];
  for (int m = UDSP_ZDF_LOW; m <= UDSP_ZDF_NOTCH; m++) {
    UDSP_ZDF f;
    udsp_zdf_reset(&f);
    udsp_zdf_set(&f, 1000.0f, 0.5f, (UDSP_ZDF_MODE)m, RATE);
    for (int i = 0; i < 4800; i++) {
      x[i] = 1.0f;
    }
    udsp_zdf_mono(&f, x, 4800);
    float want = m == UDSP_ZDF_LOW || m == UDSP_ZDF_NOTCH ? 1.0f : 0.0f;
    UDSP_TEST_ASSERT(fabsf(x[4799] - want) < 1e-4f, "mode %d at DC: %g", m,
                     (double)x[4799]);
  }
}

static void test_cutoff(void) {
  UDSP_TEST_ASSERT(fabsf(gain_at(UDSP_ZDF_LOW, 1000.0f, 0.0f, 1000.0f) - 0.5f) <
                       2e-3f,
                   "no resonance: -6 dB at the cutoff");
  float q = gain_at(UDSP_ZDF_BAND, 1000.0f, 0.75f, 1000.0f);
  UDSP_TEST_ASSERT(fabsf(q - 1.0f / (2.0f - 1.98f * 0.75f)) < 5e-3f,
                   "band at the cutoff: 1 / k (%g)", (double)q);
  UDSP_TEST_ASSERT(gain_at(UDSP_ZDF_NOTCH, 1000.0f, 0.3f, 1000.0f) < 5e-3f,
                   "the notch takes it out");
  UDSP_TEST_ASSERT(gain_at(UDSP_ZDF_LOW, 500.0f, 0.0f, 8000.0f) < 0.005f,
                   "four octaves above a low-pass: -48 dB or so");
  float top = gain_at(UDSP_ZDF_LOW, 23000.0f, 0.0f, 1000.0f);
  UDSP_TEST_ASSERT(fabsf(top - 1.0f) < 1e-3f,
                   "the cutoff clamped under half the rate, still open");
}

/* swept hard every 16 frames at full resonance: no blow-up */
static void test_swept_stable(void) {
  UDSP_ZDF f;
  udsp_zdf_reset(&f);
  uint32_t rng = 1;
  float x[16], worst = 0.0f;
  for (int b = 0; b < 20000; b++) {
    udsp_zdf_set(&f, (b * 37) % 2 ? 23000.0f : 20.0f + (float)(b % 997) * 20.0f,
                 1.0f, (UDSP_ZDF_MODE)(b / 5000), RATE);
    for (int i = 0; i < 16; i++) {
      rng ^= rng << 13;
      rng ^= rng >> 17;
      rng ^= rng << 5;
      x[i] = (float)(int32_t)rng * (1.0f / 2147483648.0f);
    }
    udsp_zdf_mono(&f, x, 16);
    for (int i = 0; i < 16; i++) {
      worst = fabsf(x[i]) > worst ? fabsf(x[i]) : worst;
    }
  }
  UDSP_TEST_ASSERT(
      worst < 200.0f && !isnan(worst),
      "noise through a cutoff jumping 20 Hz <-> 23 kHz: bounded (%g)",
      (double)worst);
}

int main(void) {
  test_dc();
  test_cutoff();
  test_swept_stable();
  return udsp_test_summary();
}
