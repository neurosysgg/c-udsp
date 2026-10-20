#include "udsp_test.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "stretch.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define RATE 48000

static float *stretched(const float *x, size_t n, double f, size_t *m) {
  const size_t win = udsp_stretch_window(RATE);
  *m = udsp_stretch_frames(n, f);
  float *y = malloc(*m * sizeof(float));
  float *s = malloc(udsp_stretch_scratch(win) * sizeof(float));
  udsp_stretch(x, n, 1, f, win, y, 1, s);
  free(s);
  return y;
}

static double rms(const float *x, size_t from, size_t to) {
  double s = 0.0;
  for (size_t i = from; i < to; i++) {
    s += (double)x[i] * x[i];
  }
  return sqrt(s / (double)(to - from));
}

/* x's frequency over [from, to): its upward zero crossings' spacing */
static double freq(const float *x, size_t from, size_t to) {
  size_t first = 0, last = 0, count = 0;
  for (size_t i = from + 1; i < to; i++) {
    if (x[i - 1] < 0.0f && x[i] >= 0.0f) {
      first = count ? first : i;
      last = i;
      count++;
    }
  }
  return count > 1 ? (double)(count - 1) * RATE / (double)(last - first) : 0;
}

static void test_window(void) {
  UDSP_TEST_ASSERT(udsp_stretch_window(48000) == 1024 &&
                       udsp_stretch_window(44100) == 1024 &&
                       udsp_stretch_window(96000) == 2048 &&
                       udsp_stretch_window(22050) == 512,
                   "the window: ~21 ms, a power of two");
  UDSP_TEST_ASSERT(udsp_stretch_frames(1000, 2.0) == 2000 &&
                       udsp_stretch_frames(1000, 0.3333) == 333 &&
                       udsp_stretch_frames(1, 0.25) == 1,
                   "the frames: rounded, at least one");
}

/* a sine stretched and squeezed: its pitch and level kept, its length
 * the factor's */
static void test_sine(void) {
  enum { N = RATE };
  float *x = malloc(N * sizeof(float));
  for (int i = 0; i < N; i++) {
    x[i] = 0.5f * (float)sin(2.0 * M_PI * 440.0 * i / RATE);
  }
  static const double F[] = {2.0, 0.5, 3.7, 1.0};
  for (size_t t = 0; t < sizeof F / sizeof *F; t++) {
    size_t m;
    float *y = stretched(x, N, F[t], &m);
    const double hz = freq(y, m / 4, 3 * m / 4);
    const double lv = rms(y, m / 4, 3 * m / 4) / (0.5 / sqrt(2.0));
    UDSP_TEST_ASSERT(m == udsp_stretch_frames(N, F[t]) &&
                         fabs(hz / 440.0 - 1.0) < 0.002 &&
                         fabs(20.0 * log10(lv)) < 0.5,
                     "x%.1f: 440 Hz kept (%.2f), its level (%.2f dB)", F[t], hz,
                     20.0 * log10(lv));
    free(y);
  }
  free(x);
}

/* a decaying noise (a reverb's tail): stretched by two, its level at
 * each point the original's at half the time, its loudness kept */
static void test_tail(void) {
  enum { N = RATE };
  float *x = malloc(N * sizeof(float));
  unsigned seed = 7;
  for (int i = 0; i < N; i++) {
    seed = seed * 1664525u + 1013904223u;
    const float v = (float)(seed >> 8) / 8388608.0f - 1.0f;
    x[i] = v * (float)exp(-6.9 * i / N); /* -60 dB at its end */
  }
  size_t m;
  float *y = stretched(x, N, 2.0, &m);
  int ok = 1;
  double worst = 0.0;
  for (int k = 1; k < 8; k++) {
    const size_t a = (size_t)k * N / 8;
    const double want = rms(x, a - 1000, a + 1000);
    const double got = rms(y, 2 * a - 2000, 2 * a + 2000);
    const double db = fabs(20.0 * log10(got / want));
    worst = db > worst ? db : worst;
    ok = ok && db < 1.5;
  }
  UDSP_TEST_ASSERT(ok && m == 2 * N,
                   "a decaying noise x2: its envelope twice as long (worst "
                   "%.2f dB)",
                   worst);
  free(y);
  free(x);
}

/* silence in, silence out; a stride each way */
static void test_strides(void) {
  enum { N = 4000 };
  static float x[2 * N], y[3 * 2 * N], z[2 * N];
  for (int i = 0; i < N; i++) {
    x[2 * i] = 0.0f;
    x[2 * i + 1] = (float)sin(0.05 * i);
  }
  const size_t win = 1024;
  float *s = malloc(udsp_stretch_scratch(win) * sizeof(float));
  memset(y, 0x7f, sizeof y);
  udsp_stretch(x, N, 2, 2.0, win, y, 3, s);
  udsp_stretch(x + 1, N, 2, 2.0, win, z, 1, s);
  int silent = 1, same = 1;
  for (int i = 0; i < 2 * N; i++) {
    silent = silent && y[3 * i] == 0.0f;
  }
  udsp_stretch(x + 1, N, 2, 2.0, win, y, 3, s);
  for (int i = 0; i < 2 * N; i++) {
    same = same && y[3 * i] == z[i];
  }
  UDSP_TEST_ASSERT(silent && same, "silence stays silent; strides the same");
  free(s);
}

int main(void) {
  printf("stretch_test\n");
  test_window();
  test_sine();
  test_tail();
  test_strides();
  return udsp_test_summary();
}
