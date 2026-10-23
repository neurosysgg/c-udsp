#include "offline/stretch.h"

#include <math.h>
#include <string.h>

#include "base/fft.h"

#define PI 3.14159265358979323846
#define TAU (2.0 * PI)

/* sin and cos of x, reduced to [-pi, pi], by their series (fifteen
 * terms: past 1e-16 there) */
static void sin_cos(double x, double *s, double *c) {
  x -= TAU * floor(x / TAU + 0.5);
  const double x2 = x * x;
  double ts = x, tc = 1.0, sv = x, cv = 1.0;
  for (int i = 1; i <= 15; i++) {
    ts *= -x2 / (double)((2 * i) * (2 * i + 1));
    tc *= -x2 / (double)((2 * i - 1) * (2 * i));
    sv += ts;
    cv += tc;
  }
  *s = sv;
  *c = cv;
}

/* atan(u) for |u| <= tan(pi / 8), by its series */
static double atan_small(double u) {
  const double u2 = u * u;
  double t = u, a = u;
  for (int k = 1; k <= 24; k++) {
    t *= -u2;
    a += t / (double)(2 * k + 1);
  }
  return a;
}

static double atan2_of(double y, double x) {
  const double ax = fabs(x), ay = fabs(y);
  if (ax == 0.0 && ay == 0.0) {
    return 0.0;
  }
  const double t = ay > ax ? ax / ay : ay / ax;
  double a = t > 0.41421356237309503
                 ? PI / 4 + atan_small((t - 1.0) / (t + 1.0))
                 : atan_small(t);
  if (ay > ax) {
    a = PI / 2 - a;
  }
  if (x < 0.0) {
    a = PI - a;
  }
  return y < 0.0 ? -a : a;
}

size_t udsp_stretch_window(uint32_t rate) {
  const double want = 0.0213 * rate;
  size_t w = 64;
  while (w < 16384 && (double)(2 * w) - want < want - (double)w) {
    w *= 2;
  }
  return w;
}

size_t udsp_stretch_frames(size_t n, double factor) {
  const size_t m = (size_t)floor((double)n * factor + 0.5);
  return m ? m : 1;
}

size_t udsp_stretch_scratch(size_t win) {
  return UDSP_FFT_TABLE_FLOATS(win) + 3 * win + 10 * (win / 2 + 1);
}

/* each bin's peak: the nearest (the lower one on a tie), into own as
 * its index; whether there was one */
static int owners(const float *mag, size_t bins, float *own) {
  long last = -1;
  for (size_t k = 0; k < bins; k++) {
    const float m = mag[k];
    const int peak = m > 0.0f && (k < 1 || m > mag[k - 1]) &&
                     (k < 2 || m > mag[k - 2]) &&
                     (k + 1 >= bins || m >= mag[k + 1]) &&
                     (k + 2 >= bins || m >= mag[k + 2]);
    if (peak) {
      last = (long)k;
    }
    own[k] = (float)last;
  }
  if (last < 0) {
    return 0;
  }
  long next = -1;
  for (size_t k = bins; k-- > 0;) {
    if (own[k] == (float)k) {
      next = (long)k;
    }
    const long prev = (long)own[k];
    if (prev < 0 || (next >= 0 && next - (long)k < (long)k - prev)) {
      own[k] = (float)next;
    }
  }
  return 1;
}

void udsp_stretch(const float *x, size_t n, size_t stride, double factor,
                  size_t win, float *y, size_t ystride, float *scratch) {
  const size_t hs = win / 4, half = win / 2, bins = half + 1;
  const size_t nout = udsp_stretch_frames(n, factor);
  float *table = scratch, *w = table + UDSP_FFT_TABLE_FLOATS(win);
  float *frame = w + win, *out = frame + win;
  float *xr = out + win, *xi = xr + bins, *pr = xi + bins, *pi = pr + bins;
  float *zr = pi + bins, *zi = zr + bins, *rr = zi + bins, *ri = rr + bins;
  float *mag = ri + bins, *own = mag + bins;
  udsp_fft_table(table, win);
  for (size_t i = 0; i < win; i++) {
    double s, c;
    sin_cos(TAU * (double)i / (double)win, &s, &c);
    w[i] = (float)(0.5 - 0.5 * c);
  }
  for (size_t i = 0; i < nout; i++) {
    y[i * ystride] = 0.0f;
  }
  /* Hann twice at a quarter's hop overlaps to 1.5 */
  const float scale = (float)(1.0 / (1.5 * (double)win));
  int64_t ta_prev = 0;
  for (int64_t m = -1; m * (int64_t)hs - (int64_t)half < (int64_t)nout; m++) {
    const int64_t ts = m * (int64_t)hs;
    const int64_t ta = (int64_t)floor((double)ts / factor + 0.5);
    for (size_t i = 0; i < win; i++) {
      const int64_t at = ta - (int64_t)half + (int64_t)i;
      frame[i] =
          at >= 0 && at < (int64_t)n ? x[(size_t)at * stride] * w[i] : 0.0f;
    }
    udsp_fft_real(frame, xr, xi, win, table);
    for (size_t k = 0; k < bins; k++) {
      mag[k] = sqrtf(xr[k] * xr[k] + xi[k] * xi[k]);
    }
    const int64_t ha = ta - ta_prev;
    if (m == -1 || !owners(mag, bins, own)) {
      /* the first frame (or one with no peak): as it is */
      for (size_t k = 0; k < bins; k++) {
        zr[k] = mag[k] > 0.0f ? xr[k] / mag[k] : 1.0f;
        zi[k] = mag[k] > 0.0f ? xi[k] / mag[k] : 0.0f;
        pr[k] = xr[k];
        pi[k] = xi[k];
      }
    } else {
      for (size_t k = 0; k < bins; k++) {
        if (own[k] != (float)k) {
          continue;
        }
        /* the peak's phase advance since the last frame, its deviation
         * from the bin's own frequency's, the synthesis hop's share */
        const double om = TAU * (double)k / (double)win;
        double adv = om * (double)hs;
        if (ha > 0) {
          const double cr = (double)xr[k] * pr[k] + (double)xi[k] * pi[k];
          const double ci = (double)xi[k] * pr[k] - (double)xr[k] * pi[k];
          double dev = atan2_of(ci, cr) - om * (double)ha;
          dev -= TAU * floor(dev / TAU + 0.5);
          adv = (om * (double)ha + dev) * (double)hs / (double)ha;
        }
        double s, c;
        sin_cos(adv, &s, &c);
        double nr = zr[k] * c - zi[k] * s, ni = zr[k] * s + zi[k] * c;
        const double nm = sqrt(nr * nr + ni * ni);
        nr /= nm;
        ni /= nm;
        /* the rotation from the analysis phase to the synthesis one */
        rr[k] = (float)((nr * xr[k] + ni * xi[k]) / mag[k]);
        ri[k] = (float)((ni * xr[k] - nr * xi[k]) / mag[k]);
      }
      for (size_t k = 0; k < bins; k++) {
        pr[k] = xr[k];
        pi[k] = xi[k];
        const size_t p = (size_t)own[k];
        const float a = xr[k] * rr[p] - xi[k] * ri[p];
        const float b = xr[k] * ri[p] + xi[k] * rr[p];
        xr[k] = a;
        xi[k] = b;
        if (mag[k] > 0.0f) {
          zr[k] = a / mag[k];
          zi[k] = b / mag[k];
        }
      }
    }
    udsp_fft_real_inverse(xr, xi, out, win, table);
    for (size_t i = 0; i < win; i++) {
      const int64_t at = ts - (int64_t)half + (int64_t)i;
      if (at >= 0 && at < (int64_t)nout) {
        y[(size_t)at * ystride] += out[i] * w[i] * scale;
      }
    }
    ta_prev = ta;
  }
}
