#include "udsp_test.h"

#include <math.h>
#include <string.h>

#include "fx/hilbert.h"

#define RATE 48000

/* a sine's parts 90 degrees apart and as loud, from 30 Hz to 18 kHz (re^2
 * + im^2 steady once the chains settle); shifting by 100 Hz leaves the
 * lower sideband 40 dB down */
static void test_quadrature(void) {
  const double hz[] = {30, 200, 1000, 5000, 18000};
  float worst = 0.0f;
  for (int k = 0; k < 5; k++) {
    UDSP_HILBERT h;
    memset(&h, 0, sizeof h);
    for (int i = 0; i < RATE; i++) {
      float re, im;
      udsp_hilbert_step(
          &h, (float)sin(2 * 3.141592653589793 * hz[k] * i / RATE), &re, &im);
      if (i > RATE / 2) {
        worst = fmaxf(worst, fabsf(sqrtf(re * re + im * im) - 1.0f));
      }
    }
  }
  UDSP_TEST_ASSERT(worst < 0.04f,
                   "30 Hz .. 18 kHz: the parts' magnitude within %.3f of 1",
                   (double)worst);
  UDSP_HILBERT h;
  memset(&h, 0, sizeof h);
  enum { N = RATE };
  static float y[N];
  for (int i = 0; i < N; i++) {
    float re, im;
    udsp_hilbert_step(&h, (float)sin(2 * 3.141592653589793 * 1000 * i / RATE),
                      &re, &im);
    const double ph = 2 * 3.141592653589793 * 100 * i / RATE;
    y[i] = re * (float)cos(ph) + im * (float)sin(ph);
  }
  /* the upper sideband's level and the lower's, by correlation */
  double up[2] = {0, 0}, lo[2] = {0, 0};
  for (int i = N / 2; i < N; i++) {
    const double a = 2 * 3.141592653589793 * 1100 * i / RATE,
                 b = 2 * 3.141592653589793 * 900 * i / RATE;
    up[0] += y[i] * cos(a);
    up[1] += y[i] * sin(a);
    lo[0] += y[i] * cos(b);
    lo[1] += y[i] * sin(b);
  }
  const double u = hypot(up[0], up[1]), l = hypot(lo[0], lo[1]);
  UDSP_TEST_ASSERT(20 * log10(u / l) > 40,
                   "1 kHz shifted 100 Hz: 1100 Hz, 900 Hz %.1f dB under it",
                   20 * log10(u / l));
}

int main(void) {
  printf("hilbert_test\n");
  test_quadrature();
  return udsp_test_summary();
}
