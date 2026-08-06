#include "fft.h"

#include <math.h>

#ifndef M_PI
/* M_PI is a POSIX/BSD math.h extension, not standard C11 -- define it
 * ourselves rather than relying on a feature-test macro being set. */
#define M_PI 3.14159265358979323846
#endif

static void swapf(float *a, float *b) {
  float t = *a;
  *a = *b;
  *b = t;
}

void udsp_fft_forward(float *re, float *im, size_t n) {
  /* bit-reversal permutation -- standard iterative in-place reordering so
   * the butterfly passes below can work purely on adjacent/strided pairs */
  for (size_t i = 1, j = 0; i < n; i++) {
    size_t bit = n >> 1;
    for (; j & bit; bit >>= 1) {
      j ^= bit;
    }
    j ^= bit;
    if (i < j) {
      swapf(&re[i], &re[j]);
      swapf(&im[i], &im[j]);
    }
  }

  for (size_t len = 2; len <= n; len <<= 1) {
    double ang = -2.0 * M_PI / (double)len;
    float wr = (float)cos(ang), wi = (float)sin(ang);
    for (size_t i = 0; i < n; i += len) {
      float cur_r = 1.0f, cur_i = 0.0f;
      for (size_t k = 0; k < len / 2; k++) {
        size_t a = i + k, b = i + k + len / 2;
        float ur = re[a], ui = im[a];
        float vr = re[b] * cur_r - im[b] * cur_i;
        float vi = re[b] * cur_i + im[b] * cur_r;

        re[a] = ur + vr;
        im[a] = ui + vi;
        re[b] = ur - vr;
        im[b] = ui - vi;

        float next_r = cur_r * wr - cur_i * wi;
        float next_i = cur_r * wi + cur_i * wr;
        cur_r = next_r;
        cur_i = next_i;
      }
    }
  }
}

void udsp_window_hann(const float *in, float *out, size_t n) {
  if (n < 2) {
    if (n == 1) {
      out[0] = in[0];
    }
    return;
  }
  for (size_t i = 0; i < n; i++) {
    double w = 0.5 - 0.5 * cos(2.0 * M_PI * (double)i / (double)(n - 1));
    out[i] = in[i] * (float)w;
  }
}

void udsp_fft_magnitude(const float *re, const float *im, size_t n, float *out) {
  for (size_t i = 0; i < n; i++) {
    out[i] = sqrtf(re[i] * re[i] + im[i] * im[i]);
  }
}
