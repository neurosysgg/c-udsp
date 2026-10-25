#include "udsp_test.h"

#include <math.h>
#include <stdint.h>
#include <string.h>

#include "fx/xover.h"

#define RATE 48000
#define N RATE

static float in[2 * N], band[5][2 * N];

/* a sine at hz, both channels */
static void sine(float hz) {
  for (int i = 0; i < N; i++) {
    in[2 * i] = in[2 * i + 1] = (float)sin(6.283185307179586 * hz * i / RATE);
  }
}

/* RMS amplitude of x's left channel over the second half (settled) */
static float amp(const float *x) {
  double sq = 0.0;
  for (int i = N / 2; i < N; i++) {
    sq += (double)x[2 * i] * x[2 * i];
  }
  return (float)sqrt(2.0 * sq / (N / 2));
}

static float db(float g) { return 20.0f * log10f(g); }

static int lr2; /* split() at 12 dB */

/* hz through a crossover of bands at splits: each band's level in
 * dB[b], the bands' sum's */
static float split(int bands, const float *splits, float hz, float *dbs) {
  UDSP_XOVER x;
  udsp_xover_reset(&x);
  if (lr2) {
    udsp_xover_set_lr2(&x, bands, splits, RATE);
  } else {
    udsp_xover_set(&x, bands, splits, RATE);
  }
  sine(hz);
  float *out[5] = {band[0], band[1], band[2], band[3], band[4]};
  udsp_xover_split(&x, in, out, N);
  static float sum[2 * N];
  memset(sum, 0, sizeof sum);
  for (int b = 0; b < bands; b++) {
    dbs[b] = db(amp(band[b]));
    for (int i = 0; i < 2 * N; i++) {
      sum[i] += band[b][i];
    }
  }
  return db(amp(sum));
}

static void test_flat(void) {
  static const float SPLITS[] = {120.0f, 2500.0f, 6000.0f, 12000.0f};
  static const float HZ[] = {30, 120, 400, 2500, 4000, 6000, 9000, 15000};
  float worst = 0.0f, dbs[5];
  for (int bands = 2; bands <= 5; bands++) {
    for (size_t i = 0; i < sizeof HZ / sizeof *HZ; i++) {
      const float s = split(bands, SPLITS, HZ[i], dbs);
      worst = fabsf(s) > worst ? fabsf(s) : worst;
    }
  }
  UDSP_TEST_ASSERT(worst < 0.01f,
                   "2 .. 5 bands summed: flat within %.4f dB (< 0.01)",
                   (double)worst);
}

static void test_bands(void) {
  static const float SPLITS[] = {120.0f, 2500.0f};
  float dbs[3];
  split(3, SPLITS, 40.0f, dbs);
  /* a Linkwitz-Riley low side: 1 / (1 + w^4), w = 40 / 120 */
  const float lr4 = -20.0f * log10f(1.0f + 1.0f / 81.0f);
  UDSP_TEST_ASSERT(fabsf(dbs[0] - lr4) < 0.005f && dbs[1] < -35.0f &&
                       dbs[2] < -60.0f,
                   "40 Hz: the low band, %.3f dB as LR4's %.3f (%.1f, %.1f dB)",
                   (double)dbs[0], (double)lr4, (double)dbs[1], (double)dbs[2]);
  split(3, SPLITS, 600.0f, dbs);
  UDSP_TEST_ASSERT(dbs[1] > -0.3f && dbs[0] < -45.0f && dbs[2] < -35.0f,
                   "600 Hz: the middle (%.1f, %.2f, %.1f dB)", (double)dbs[0],
                   (double)dbs[1], (double)dbs[2]);
  split(3, SPLITS, 120.0f, dbs);
  UDSP_TEST_ASSERT(fabsf(dbs[0] + 6.02f) < 0.1f && fabsf(dbs[1] + 6.02f) < 0.2f,
                   "at a split: both sides -6 dB (%.2f, %.2f)", (double)dbs[0],
                   (double)dbs[1]);
  split(3, SPLITS, 10000.0f, dbs);
  UDSP_TEST_ASSERT(dbs[2] > -0.1f && dbs[1] < -35.0f,
                   "10 kHz: the high band (%.1f, %.2f)", (double)dbs[1],
                   (double)dbs[2]);
}

/* one band: the input as it is; in may be the last output */
static void test_one_and_alias(void) {
  UDSP_XOVER x;
  udsp_xover_reset(&x);
  udsp_xover_set(&x, 1, NULL, RATE);
  sine(1000.0f);
  memcpy(band[0], in, sizeof in);
  float *out[1] = {band[0]};
  udsp_xover_split(&x, band[0], out, N);
  UDSP_TEST_ASSERT(!memcmp(band[0], in, sizeof in), "one band: the input");
  static const float S[] = {500.0f};
  float *o2[2] = {band[1], band[0]};
  udsp_xover_reset(&x);
  udsp_xover_set(&x, 2, S, RATE);
  udsp_xover_split(&x, band[0], o2, N);
  float worst = 0.0f;
  UDSP_XOVER y;
  udsp_xover_reset(&y);
  udsp_xover_set(&y, 2, S, RATE);
  float *o3[2] = {band[2], band[3]};
  udsp_xover_split(&y, in, o3, N);
  for (int i = 0; i < 2 * N; i++) {
    float d = fabsf(band[0][i] - band[3][i]) + fabsf(band[1][i] - band[2][i]);
    worst = d > worst ? d : worst;
  }
  UDSP_TEST_ASSERT(worst == 0.0f, "in as the last output: the same bits");
}

/* the section-at-a-time crossover (before the wavefront), three bands:
 * the reference */
static void ref_tick(const UDSP_XOVER_SPLIT *c, float x, float *ic, float *v1,
                     float *v2) {
  const float v3 = x - ic[1];
  *v1 = c->a1 * ic[0] + c->a2 * v3;
  *v2 = ic[1] + c->a2 * ic[0] + c->a3 * v3;
  ic[0] = 2.0f * *v1 - ic[0];
  ic[1] = 2.0f * *v2 - ic[1];
}

static void ref_split3(UDSP_XOVER *x, const float *from, float *const *out,
                       size_t frames) {
  memmove(out[2], from, frames * 2 * sizeof(float));
  for (int s = 0; s < 2; s++) {
    const UDSP_XOVER_SPLIT *c = &x->c[s];
    for (size_t i = 0; i < frames; i++) {
      for (int ch = 0; ch < 2; ch++) {
        const float v = out[2][2 * i + ch];
        float v1, v2, w1, w2;
        ref_tick(c, v, x->sec[s][0][ch], &v1, &v2);
        const float hp = v - c->k * v1 - v2;
        ref_tick(c, v2, x->sec[s][1][ch], &w1, &w2);
        out[s][2 * i + ch] = w2;
        ref_tick(c, hp, x->sec[s][2][ch], &w1, &w2);
        out[2][2 * i + ch] = hp - c->k * w1 - w2;
      }
    }
  }
  for (size_t i = 0; i < frames; i++) {
    for (int ch = 0; ch < 2; ch++) {
      const float v = out[0][2 * i + ch];
      float v1, v2;
      ref_tick(&x->c[1], v, x->ap[0][1][ch], &v1, &v2);
      out[0][2 * i + ch] = v - 2.0f * x->c[1].k * v1;
    }
  }
}

/* three bands (Multiband Dynamics'), run as a wavefront: the reference's
 * bits at any call length, retuned between calls, the input in place */
static void test_three_as_reference(void) {
  static float got[3][2 * 600], want[3][2 * 600], src[2 * 600];
  uint32_t seed = 7;
  for (int i = 0; i < 2 * 600; i++) {
    seed = seed * 1664525u + 1013904223u;
    src[i] = (float)(seed >> 8) / 16777216.0f - 0.5f;
  }
  UDSP_XOVER a, b;
  udsp_xover_reset(&a);
  udsp_xover_reset(&b);
  static const size_t LENS[] = {1, 2, 3, 4, 5, 7, 64, 256, 1, 3, 254};
  int same = 1, at = 0;
  for (size_t r = 0; r < sizeof LENS / sizeof *LENS; r++) {
    const float hz[2] = {80.0f + 40.0f * (float)r, 2000.0f + 300.0f * (float)r};
    udsp_xover_set(&a, 3, hz, RATE);
    udsp_xover_set(&b, 3, hz, RATE);
    const size_t n = LENS[r];
    float *g[3] = {got[0], got[1], got[2]}, *w[3] = {want[0], want[1], want[2]};
    const float *x = src + 2 * (at % 300);
    if (r % 2) { /* in place: the input in the last band's buffer */
      memcpy(got[2], x, n * 2 * sizeof(float));
      x = got[2];
    }
    udsp_xover_split(&a, x, g, n);
    ref_split3(&b, src + 2 * (at % 300), w, n);
    for (int k = 0; k < 3; k++) {
      same &= memcmp(got[k], want[k], n * 2 * sizeof(float)) == 0;
    }
    same &= memcmp(a.sec, b.sec, sizeof a.sec) == 0 &&
            memcmp(a.ap[0][1], b.ap[0][1], sizeof a.ap[0][1]) == 0;
    at += (int)n;
  }
#if defined(__clang__) && defined(__FMA__)
  /* clang fuses a * b + c by default (DESIGN.md "Bits"): no bits to hold */
  (void)same;
#else
  UDSP_TEST_ASSERT(same, "three bands: the section-at-a-time bits and state, "
                         "calls of 1 .. 256, retuned, in place");
#endif
}

/* 12 dB (LFOTool's split): Linkwitz-Riley's second order -- each side
 * 1 / (1 + w^2) in amplitude, -6 dB at the split; 2 .. 5 bands summed flat
 * (the high side inverted, the lower bands through the allpasses); the
 * same bits whatever the calls' lengths */
static void test_lr2(void) {
  lr2 = 1;
  static const float S[] = {500.0f};
  float dbs[5];
  const float side = -20.0f * log10f(1.25f); /* w = 2 or 1/2 */
  split(2, S, 250.0f, dbs);
  UDSP_TEST_ASSERT(fabsf(dbs[0] - side) < 0.01f &&
                       fabsf(dbs[1] - (side - 12.04f)) < 0.05f,
                   "an octave under: the low side %.3f dB (LR2's %.3f), the "
                   "high %.2f",
                   (double)dbs[0], (double)side, (double)dbs[1]);
  split(2, S, 1000.0f, dbs);
  UDSP_TEST_ASSERT(fabsf(dbs[1] - side) < 0.01f &&
                       fabsf(dbs[0] - (side - 12.04f)) < 0.05f,
                   "an octave over: the high side %.3f dB, the low %.2f",
                   (double)dbs[1], (double)dbs[0]);
  split(2, S, 500.0f, dbs);
  UDSP_TEST_ASSERT(
      fabsf(dbs[0] + 6.02f) < 0.02f && fabsf(dbs[1] + 6.02f) < 0.02f,
      "at the split: both -6 dB (%.2f, %.2f)", (double)dbs[0], (double)dbs[1]);
  static const float SPLITS[] = {120.0f, 2500.0f, 6000.0f, 12000.0f};
  static const float HZ[] = {30, 120, 400, 2500, 4000, 6000, 9000, 15000};
  float worst = 0.0f;
  for (int bands = 2; bands <= 5; bands++) {
    for (size_t i = 0; i < sizeof HZ / sizeof *HZ; i++) {
      const float sum = split(bands, SPLITS, HZ[i], dbs);
      worst = fabsf(sum) > worst ? fabsf(sum) : worst;
    }
  }
  UDSP_TEST_ASSERT(worst < 0.01f,
                   "2 .. 5 bands summed: flat within %.4f dB (< 0.01)",
                   (double)worst);
  lr2 = 0;
  static float whole[3][2 * 600], part[3][2 * 600];
  uint32_t seed = 11;
  for (int i = 0; i < 2 * 600; i++) {
    seed = seed * 1664525u + 1013904223u;
    in[i] = (float)(seed >> 8) / 16777216.0f - 0.5f;
  }
  static const float S3[] = {200.0f, 3000.0f};
  UDSP_XOVER a, b;
  udsp_xover_reset(&a);
  udsp_xover_reset(&b);
  udsp_xover_set_lr2(&a, 3, S3, RATE);
  udsp_xover_set_lr2(&b, 3, S3, RATE);
  float *w[3] = {whole[0], whole[1], whole[2]};
  udsp_xover_split(&a, in, w, 600);
  static const size_t LENS[] = {1, 2, 3, 7, 64, 256, 267};
  size_t at = 0;
  for (size_t r = 0; r < sizeof LENS / sizeof *LENS; r++) {
    float *q[3] = {part[0] + 2 * at, part[1] + 2 * at, part[2] + 2 * at};
    udsp_xover_split(&b, in + 2 * at, q, LENS[r]);
    at += LENS[r];
  }
  int same = at == 600;
  for (int k = 0; k < 3; k++) {
    for (int i = 0; i < 2 * 600; i++) {
      same &= whole[k][i] == part[k][i];
    }
  }
  UDSP_TEST_ASSERT(same, "three bands at 12 dB: the same values in calls of "
                         "1 .. 267");
}

int main(void) {
  test_three_as_reference();
  printf("xover_test\n");
  test_flat();
  test_bands();
  test_one_and_alias();
  test_lr2();
  return udsp_test_summary();
}
