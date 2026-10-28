#include "udsp_test.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "fx/onset.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#define RATE 48000.0f
#define LEN 48000

/* kicks: a cosine at 60 Hz from its peak, decaying over 120 ms, at each
 * of the frames in at (n of them), amp high */
static float *kicks(const int *at, int n, float amp) {
  float *b = calloc(2 * LEN, sizeof *b);
  for (int k = 0; k < n; k++) {
    for (int i = at[k]; i < LEN; i++) {
      const double t = (i - at[k]) / (double)RATE;
      const float v = (float)(amp * cos(2 * M_PI * 60 * t) * exp(-t / 0.12));
      b[2 * i] += v;
      b[2 * i + 1] += v;
    }
  }
  return b;
}

/* the onsets in b, block frames at a time */
static size_t onsets(const float *b, size_t block, uint32_t *out, size_t cap) {
  UDSP_ONSET o;
  udsp_onset_init(&o, RATE);
  size_t n = 0;
  for (size_t at = 0; at < LEN; at += block) {
    const size_t k = LEN - at < block ? LEN - at : block;
    uint32_t got[16];
    const size_t m = udsp_onset_run(&o, b + 2 * at, k, got, 16);
    for (size_t i = 0; i < m && n < cap; i++) {
      out[n++] = (uint32_t)(at + got[i]);
    }
  }
  return n;
}

static void test_kicks(void) {
  static const int at[] = {1000, 13000, 25000, 37000};
  float *b = kicks(at, 4, 0.8f);
  uint32_t o[8];
  size_t n = onsets(b, 256, o, 8);
  UDSP_TEST_ASSERT(n == 4 && o[0] == 1000 && o[1] == 13000 && o[2] == 25000 &&
                       o[3] == 37000,
                   "four kicks: an onset on each one's first frame (%zu)", n);
  const size_t blocks[] = {1, 7, 333};
  for (int k = 0; k < 3; k++) {
    uint32_t p[8];
    size_t m = onsets(b, blocks[k], p, 8);
    UDSP_TEST_ASSERT(m == n && !memcmp(p, o, n * sizeof *o),
                     "blocks of %zu: the same frames", blocks[k]);
  }
  free(b);
}

static void test_quiet_and_close(void) {
  float *b = calloc(2 * LEN, sizeof *b);
  uint32_t o[8];
  UDSP_TEST_ASSERT(onsets(b, 256, o, 8) == 0, "silence: none");
  free(b);
  static const int at[] = {2000};
  b = kicks(at, 1, 0.01f); /* -40 dB */
  UDSP_TEST_ASSERT(onsets(b, 256, o, 8) == 0, "under -30 dBFS: none");
  free(b);
  static const int two[] = {2000, 2000 + 480}; /* 10 ms apart */
  b = kicks(two, 2, 0.8f);
  size_t n = onsets(b, 256, o, 8);
  UDSP_TEST_ASSERT(n == 1 && o[0] == 2000, "10 ms apart: one (40 ms guard)");
  free(b);
}

static void test_tone(void) {
  /* a sine held at 0.5 from frame 100: its start only -- at 110 Hz, and
   * at 40 Hz (an 808's tail, its peaks 12.5 ms apart) */
  static const double hz[] = {110, 40};
  for (int k = 0; k < 2; k++) {
    float *b = calloc(2 * LEN, sizeof *b);
    for (int i = 100; i < LEN; i++) {
      const float v = (float)(0.5 * cos(2 * M_PI * hz[k] * (i - 100) / RATE));
      b[2 * i] = b[2 * i + 1] = v;
    }
    uint32_t o[8];
    size_t n = onsets(b, 256, o, 8);
    UDSP_TEST_ASSERT(n == 1 && o[0] == 100,
                     "a tone held at %g Hz: its start, once (%zu)", hz[k], n);
    free(b);
  }
}

int main(void) {
  test_kicks();
  test_quiet_and_close();
  test_tone();
  return udsp_test_summary();
}
