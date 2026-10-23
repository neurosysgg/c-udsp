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

/* the filter as it was, a frame at a time with the mode's switch inside */
static void ref_mono(UDSP_ZDF *f, float *io, size_t n) {
  float ic1 = f->ic1, ic2 = f->ic2;
  const float k = f->k, a1 = f->a1, a2 = f->a2, a3 = f->a3;
  for (size_t i = 0; i < n; i++) {
    const float v0 = io[i];
    const float v3 = v0 - ic2;
    const float v1 = a1 * ic1 + a2 * v3;
    const float v2 = ic2 + a2 * ic1 + a3 * v3;
    ic1 = 2.0f * v1 - ic1;
    ic2 = 2.0f * v2 - ic2;
    switch (f->mode) {
    case UDSP_ZDF_BAND:
      io[i] = v1;
      break;
    case UDSP_ZDF_HIGH:
      io[i] = v0 - k * v1 - v2;
      break;
    case UDSP_ZDF_NOTCH:
      io[i] = v0 - k * v1;
      break;
    default:
      io[i] = v2;
      break;
    }
  }
  f->ic1 = ic1;
  f->ic2 = ic2;
}

/* alone and side by side (1..6 filters, modes mixed into runs), cutoffs
 * moved between calls: the bits and state of the loop as it was */
static void test_against_ref(void) {
  enum { MAX = 6, LEN = 97 };
  UDSP_ZDF got[MAX] = {0}, want[MAX] = {0};
  static float x[MAX][LEN], y[MAX][LEN];
  uint32_t seed = 3;
  int same = 1;
  for (int trial = 0; trial < 300; trial++) {
    const size_t count = 1 + (size_t)trial % MAX;
    for (size_t j = 0; j < count; j++) {
      seed = seed * 1664525u + 1013904223u;
      const UDSP_ZDF_MODE m =
          (UDSP_ZDF_MODE)((seed >> 9) % 4 == 3 ? (seed >> 13) % 4
                                               : (uint32_t)trial / 75);
      udsp_zdf_reset(&got[j]);
      udsp_zdf_set(&got[j], 40.0f + (float)((seed >> 16) % 12000),
                   (float)((seed >> 4) % 100) / 100.0f, m, RATE);
      want[j] = got[j];
    }
    for (int call = 0; call < 4; call++) {
      for (size_t j = 0; j < count; j++) {
        for (int i = 0; i < LEN; i++) {
          seed = seed * 1664525u + 1013904223u;
          x[j][i] = y[j][i] = (float)(seed >> 8) / 8388608.0f - 1.0f;
        }
        ref_mono(&want[j], y[j], LEN);
      }
      if (count == 1) {
        udsp_zdf_mono(&got[0], x[0], LEN);
      } else {
        UDSP_ZDF *f[MAX];
        float *io[MAX];
        for (size_t j = 0; j < count; j++) {
          f[j] = &got[j];
          io[j] = x[j];
        }
        udsp_zdf_monos(f, io, count, LEN);
      }
      same &= !memcmp(x, y, sizeof(x)) &&
              !memcmp(got, want, count * sizeof(got[0]));
    }
  }
  UDSP_TEST_ASSERT(same, "alone and side by side: the loop's bits and state");
}

int main(void) {
  test_dc();
  test_cutoff();
  test_swept_stable();
  test_against_ref();
  return udsp_test_summary();
}
