#include "analysis/truepeak.h"

#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* I0, for the Kaiser window */
static double bessel0(double x) {
  double sum = 1.0, term = 1.0;
  for (int k = 1; k < 32; k++) {
    term *= (x / (2.0 * k)) * (x / (2.0 * k));
    sum += term;
  }
  return sum;
}

void udsp_truepeak_init(UDSP_TRUEPEAK *t, int channels) {
  enum { N = UDSP_TRUEPEAK_PHASES * UDSP_TRUEPEAK_TAPS };
  const double beta = 7.0, mid = (N - 1) / 2.0;
  memset(t, 0, sizeof(*t));
  t->channels = channels;
  for (int p = 0; p < UDSP_TRUEPEAK_PHASES; p++) {
    double sum = 0.0;
    double h[UDSP_TRUEPEAK_TAPS];
    for (int k = 0; k < UDSP_TRUEPEAK_TAPS; k++) {
      const int i = k * UDSP_TRUEPEAK_PHASES + p;
      const double x = (i - mid) / UDSP_TRUEPEAK_PHASES;
      const double u = (i - mid) / mid;
      const double sinc = sin(M_PI * x) / (M_PI * x);
      h[k] = sinc * bessel0(beta * sqrt(1.0 - u * u)) / bessel0(beta);
      sum += h[k];
    }
    for (int k = 0; k < UDSP_TRUEPEAK_TAPS; k++) {
      t->h[p][k] = (float)(h[k] / sum);
    }
  }
}

void udsp_truepeak_add(UDSP_TRUEPEAK *t, const float *interleaved,
                       size_t frames) {
  float peak = t->peak;
  int at = t->at;
  for (size_t f = 0; f < frames; f++) {
    for (int c = 0; c < t->channels; c++) {
      const float v = interleaved[f * (size_t)t->channels + (size_t)c];
      t->z[c][at] = t->z[c][at + UDSP_TRUEPEAK_TAPS] = v;
      /* the window, oldest first: taps run the other way round the phase */
      const float *w = &t->z[c][at + 1];
      for (int p = 0; p < UDSP_TRUEPEAK_PHASES; p++) {
        float y = 0.0f;
        for (int k = 0; k < UDSP_TRUEPEAK_TAPS; k++) {
          y += t->h[p][k] * w[k];
        }
        y = fabsf(y);
        peak = y > peak ? y : peak;
      }
    }
    at = at + 1 == UDSP_TRUEPEAK_TAPS ? 0 : at + 1;
  }
  t->at = at;
  t->peak = peak;
}
