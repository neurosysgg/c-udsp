#include "udsp_test.h"

#include "yin.h"

#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define RATE 44100
#define N 8192

static float x[N], scratch[RATE / 20 + 2];

static float detect(float *ap) {
  return udsp_yin(x, N, RATE, 20.0f, 4200.0f, 0.15f, scratch, ap);
}

static void test_sine(void) {
  for (int i = 0; i < N; i++) {
    x[i] = 0.5f * (float)sin(2.0 * M_PI * 440.0 * i / RATE);
  }
  float ap;
  float hz = detect(&ap);
  UDSP_TEST_ASSERT(fabsf(hz - 440.0f) < 0.05f && ap < 0.01f,
                   "a 440 Hz sine: %.3f Hz, aperiodicity %.4f", hz, ap);
}

/* a piano's lowest A, its fundamental weaker than its next partials */
static void test_low_a(void) {
  for (int i = 0; i < N; i++) {
    double t = (double)i / RATE, v = 0;
    for (int k = 1; k <= 8; k++) {
      v += (k == 1 ? 0.2 : 1.0 / k) * sin(2.0 * M_PI * 27.5 * k * t + k);
    }
    x[i] = (float)(0.3 * v);
  }
  float hz = detect(NULL);
  UDSP_TEST_ASSERT(fabsf(hz - 27.5f) < 0.05f,
                   "A0 with a weak fundamental: %.3f Hz, not an octave up", hz);
}

static void test_saw(void) {
  double phase = 0;
  for (int i = 0; i < N; i++) {
    x[i] = (float)(phase * 2.0 - 1.0);
    phase += 146.832 / RATE; /* D3 */
    phase -= floor(phase);
  }
  float hz = detect(NULL);
  float cents;
  int note = udsp_yin_note(hz, &cents);
  UDSP_TEST_ASSERT(note == 50 && fabsf(cents) < 2.0f,
                   "a D3 saw: note %d, %+.2f cents", note, cents);
}

static void test_noise_and_limits(void) {
  unsigned s = 1;
  for (int i = 0; i < N; i++) {
    s = s * 1664525u + 1013904223u;
    x[i] = (float)(s >> 8) / 8388608.0f - 1.0f;
  }
  float ap;
  float hz = detect(&ap);
  UDSP_TEST_ASSERT(hz == 0.0f && ap == 1.0f, "noise: no pitch");
  UDSP_TEST_ASSERT(
      udsp_yin(x, 1000, RATE, 20.0f, 4200.0f, 0.15f, scratch, NULL) == 0.0f,
      "too short for the lowest period: 0");
  float c;
  UDSP_TEST_ASSERT(udsp_yin_note(440.0f, &c) == 69 && fabsf(c) < 1e-3f &&
                       udsp_yin_note(452.0f, &c) == 69 && c > 46.0f &&
                       c < 47.0f && udsp_yin_note(0.0f, &c) == -1,
                   "notes: 440 is A4 (69), 452 its +46.6 cents, 0 none");
}

int main(void) {
  test_sine();
  test_low_a();
  test_saw();
  test_noise_and_limits();
  return udsp_test_summary();
}
