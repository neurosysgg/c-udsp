#include "yin.h"

#include "simd.h"

#include <math.h>

/* sum over j < w of (x[j] - x[j + tau])^2: vectors across j */
static double diff_sq(const float *x, const float *y, size_t w) {
  UDSP_VF acc[UDSP_ACC];
  for (int a = 0; a < UDSP_ACC; a++) {
    acc[a] = udsp_vsplat(0.0f);
  }
  size_t i = 0;
  for (; i + UDSP_ACC * UDSP_W <= w; i += UDSP_ACC * UDSP_W) {
    UDSP_UNROLL
    for (int a = 0; a < UDSP_ACC; a++) {
      UDSP_VF d = udsp_vload(x + i + (size_t)a * UDSP_W) -
                  udsp_vload(y + i + (size_t)a * UDSP_W);
      acc[a] += d * d;
    }
  }
  for (int a = 1; a < UDSP_ACC; a++) {
    acc[0] += acc[a];
  }
  double sum = udsp_vsum(acc[0]);
  for (; i < w; i++) {
    double d = (double)x[i] - y[i];
    sum += d * d;
  }
  return sum;
}

float udsp_yin(const float *x, size_t n, uint32_t rate, float fmin, float fmax,
               float threshold, float *scratch, float *aperiodicity) {
  if (aperiodicity) {
    *aperiodicity = 1.0f;
  }
  if (fmin <= 0.0f || fmax <= fmin || rate == 0) {
    return 0.0f;
  }
  size_t tau_max = (size_t)((float)rate / fmin) + 1;
  size_t tau_min = (size_t)((float)rate / fmax);
  tau_min = tau_min < 2 ? 2 : tau_min;
  if (n < 2 * tau_max || tau_max <= tau_min + 1) {
    return 0.0f;
  }
  const size_t w = n - tau_max;
  /* d'(tau): the difference over its mean up to tau (d'(0) = 1) */
  double running = 0.0;
  scratch[0] = 1.0f;
  for (size_t tau = 1; tau <= tau_max; tau++) {
    double d = diff_sq(x, x + tau, w);
    running += d;
    scratch[tau] = running > 0.0 ? (float)(d * (double)tau / running) : 1.0f;
  }
  /* the first dip under the threshold, followed to its bottom; else the
   * lowest dip when even that's too shallow to trust */
  size_t best = 0;
  for (size_t tau = tau_min; tau < tau_max; tau++) {
    if (scratch[tau] < threshold) {
      while (tau + 1 < tau_max && scratch[tau + 1] < scratch[tau]) {
        tau++;
      }
      best = tau;
      break;
    }
  }
  if (!best) {
    return 0.0f;
  }
  /* the parabola through the dip and its neighbours */
  double a = scratch[best - 1], b = scratch[best], c = scratch[best + 1];
  double den = a - 2.0 * b + c;
  double shift = den > 0.0 ? 0.5 * (a - c) / den : 0.0;
  shift = shift > 1.0 ? 1.0 : (shift < -1.0 ? -1.0 : shift);
  if (aperiodicity) {
    *aperiodicity = (float)b;
  }
  return (float)((double)rate / ((double)best + shift));
}

int udsp_yin_note(float hz, float *cents) {
  if (hz <= 0.0f) {
    if (cents) {
      *cents = 0.0f;
    }
    return -1;
  }
  double note = 69.0 + 12.0 * log2((double)hz / 440.0);
  double whole = floor(note + 0.5);
  if (cents) {
    *cents = (float)((note - whole) * 100.0);
  }
  return (int)whole;
}
