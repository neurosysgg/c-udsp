#include "udsp_test.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "fx/fdn.h"
#include "fx/line.h"

#define RATE 48000

static void test_line(void) {
  float buf[64];
  UDSP_LINE l;
  UDSP_TEST_ASSERT(udsp_line_size(10) == 16 && udsp_line_size(60) == 64,
                   "a ring a power of two, 4 past the longest delay");
  udsp_line_init(&l, buf, 64);
  for (int i = 0; i < 100; i++) {
    udsp_line_push(&l, (float)i);
  }
  UDSP_TEST_ASSERT(udsp_line_tap(&l, 1) == 99.0f &&
                       udsp_line_tap(&l, 40) == 60.0f,
                   "a delay of 1 is the last pushed, 40 back is 60");
  UDSP_TEST_ASSERT(udsp_line_read(&l, 10.0f) == 90.0f &&
                       fabsf(udsp_line_read(&l, 10.25f) - 89.75f) < 1e-4f,
                   "whole delays exact; a ramp read between its frames");
}

/* reads a run ahead, then the run pushed: udsp_line_read's answers
 * frame by frame, its pushes between -- across the ring's wrap, any
 * length, delays down to n + 1 */
static void test_line_reads(void) {
  float buf[128], ref_buf[128], d[40], got[40], x[40];
  float worst = 0.0f;
  int same = 1, ring = 1;
  uint32_t seed = 12345;
  for (int round = 0; round < 60; round++) {
    UDSP_LINE l, r;
    udsp_line_init(&l, buf, 128);
    udsp_line_init(&r, ref_buf, 128);
    const int pre = 90 + round * 7; /* where the run starts in the ring */
    for (int i = 0; i < pre; i++) {
      seed = seed * 1664525u + 1013904223u;
      const float v = (float)(seed >> 8) / 16777216.0f - 0.5f;
      udsp_line_push(&l, v);
      udsp_line_push(&r, v);
    }
    const size_t n = (size_t)(round % 40) + 1;
    for (size_t j = 0; j < n; j++) {
      seed = seed * 1664525u + 1013904223u;
      d[j] = (float)(n + 1) + (float)(seed >> 8) / 16777216.0f * 80.0f;
      x[j] = (float)j - 3.5f;
    }
    d[0] = (float)(n + 1); /* the nearest allowed, whole */
    udsp_line_reads(&l, d, got, n);
    udsp_line_write(&l, x, n);
    for (size_t j = 0; j < n; j++) {
      const float want = udsp_line_read(&r, d[j]);
      udsp_line_push(&r, x[j]);
      const float e = fabsf(got[j] - want);
      worst = e > worst ? e : worst;
      same &= memcmp(&got[j], &want, sizeof want) == 0;
    }
    ring &= l.w == r.w && memcmp(buf, ref_buf, sizeof buf) == 0;
  }
  UDSP_TEST_ASSERT(ring, "the run's pushes: the same ring, the same place");
#if defined(__clang__) && defined(__FMA__)
  /* clang fuses a * b + c by default (DESIGN.md "Bits") */
  UDSP_TEST_ASSERT(worst < 1e-6f, "the reads a frame at a time (%g, fused)",
                   (double)worst);
#else
  (void)worst;
  UDSP_TEST_ASSERT(same, "the reads a frame at a time, bit for bit");
#endif
}

static float *fdn_new(UDSP_FDN *f) {
  float *mem = calloc(udsp_fdn_floats(RATE), sizeof(float));
  udsp_fdn_init(f, mem, RATE);
  return mem;
}

/* an impulse's tail: its level falls 60 dB in the decay time */
static void test_decay(void) {
  UDSP_FDN f;
  float *mem = fdn_new(&f);
  udsp_fdn_set(&f, 1.0f, 1.5f, 20000.0f, 0.6f, 0.0f, 0.0f, 1.0f);
  const int n = RATE * 2;
  float *x = calloc((size_t)n * 2, sizeof(float));
  x[0] = x[1] = 1.0f;
  udsp_fdn_stereo(&f, x, x, (size_t)n);
  /* energy over 100 ms windows at 0.3 s and 1.05 s: 0.75 s apart, 30 dB */
  double e[2] = {0, 0};
  const int at[2] = {RATE * 3 / 10, RATE * 105 / 100};
  for (int w = 0; w < 2; w++) {
    for (int i = at[w]; i < at[w] + RATE / 10; i++) {
      e[w] += (double)x[2 * i] * x[2 * i] + (double)x[2 * i + 1] * x[2 * i + 1];
    }
  }
  double drop = 10.0 * log10(e[0] / e[1]);
  UDSP_TEST_ASSERT(fabs(drop - 30.0) < 3.0,
                   "1.5 s of decay: 30 dB down over 0.75 s (%.1f dB)", drop);
  /* left and right differ: the orthogonal sums */
  double lr = 0, ll = 0, rr = 0;
  for (int i = RATE / 10; i < RATE / 2; i++) {
    lr += (double)x[2 * i] * x[2 * i + 1];
    ll += (double)x[2 * i] * x[2 * i];
    rr += (double)x[2 * i + 1] * x[2 * i + 1];
  }
  double corr = lr / sqrt(ll * rr);
  UDSP_TEST_ASSERT(fabs(corr) < 0.3, "left and right uncorrelated (%.3f)",
                   corr);
  udsp_fdn_clear(&f);
  memset(x, 0, (size_t)n * 2 * sizeof(float));
  x[0] = x[1] = 1.0f;
  udsp_fdn_set(&f, 1.0f, 1.5f, 20000.0f, 0.6f, 0.0f, 0.0f, 0.0f);
  udsp_fdn_stereo(&f, x, x, (size_t)n);
  int mono = 1;
  for (int i = 0; i < n; i++) {
    mono &= x[2 * i] == x[2 * i + 1];
  }
  UDSP_TEST_ASSERT(mono, "width 0: mono");
  free(x);
  free(mem);
}

/* the longest decay, damping open, modulated: 17 s of loud noise, then
 * 43 of silence falling 1 dB a second (no build-up, no blow-up); in
 * blocks of 333, then of 64 for the same bits */
static void test_stable_and_blocks(void) {
  UDSP_FDN f, g;
  float *ma = fdn_new(&f), *mb = fdn_new(&g);
  udsp_fdn_set(&f, 2.0f, 60.0f, 22000.0f, 1.0f, 1.0f, 30.0f, 1.0f);
  udsp_fdn_set(&g, 2.0f, 60.0f, 22000.0f, 1.0f, 1.0f, 30.0f, 1.0f);
  enum { BLOCK = 4800 };
  static float a[2 * BLOCK], b[2 * BLOCK];
  uint32_t rng = 1;
  double e[60] = {0};
  int same = 1;
  for (int k = 0; k < 600; k++) {
    for (int i = 0; i < 2 * BLOCK; i++) {
      rng ^= rng << 13;
      rng ^= rng >> 17;
      rng ^= rng << 5;
      a[i] = b[i] = (k < 170 ? 1.0f : 0.0f) * (float)(int32_t)rng *
                    (1.0f / 2147483648.0f);
    }
    for (int i = 0; i < BLOCK; i += 333) {
      udsp_fdn_stereo(&f, a + 2 * i, a + 2 * i,
                      i + 333 <= BLOCK ? 333 : BLOCK - i);
    }
    for (int i = 0; i < BLOCK; i += 64) {
      udsp_fdn_stereo(&g, b + 2 * i, b + 2 * i,
                      i + 64 <= BLOCK ? 64 : BLOCK - i);
    }
    for (int i = 0; i < 2 * BLOCK; i++) {
      e[k / 10] += (double)a[i] * a[i];
      same &= a[i] == b[i];
    }
  }
  double drop = 10.0 * log10(e[20] / e[59]);
  UDSP_TEST_ASSERT(isfinite(drop) && fabs(drop - 39.0) < 4.0,
                   "a 60 s decay after the noise: 39 dB in 39 s (%.1f)", drop);
  UDSP_TEST_ASSERT(same, "blocks of 333 or 64: the same bits");
  free(ma);
  free(mb);
}

static void test_predelay(void) {
  UDSP_FDN f;
  float *mem = fdn_new(&f);
  udsp_fdn_set(&f, 1.0f, 1.0f, 8000.0f, 0.0f, 0.0f, 100.0f, 1.0f);
  static float x[2 * RATE / 4];
  x[0] = x[1] = 1.0f;
  udsp_fdn_stereo(&f, x, x, RATE / 4);
  int first = -1;
  for (int i = 0; i < RATE / 4 && first < 0; i++) {
    first = x[2 * i] != 0.0f || x[2 * i + 1] != 0.0f ? i : -1;
  }
  /* 100 ms, the diffusers (short delays of their own, even at none, so
   * the timing doesn't jump as the knob moves), the shortest line */
  int want = 4800 + (int)f.len[0];
  for (int k = 0; k < UDSP_FDN_DIFFUSERS; k++) {
    want += (int)f.ap_len[k];
  }
  UDSP_TEST_ASSERT(first >= want - 2 && first <= want + 1,
                   "100 ms of predelay, the diffusers, the first line: "
                   "frame %d (%d)",
                   first, want);
  free(mem);
}

/* idle: the modulation where that much silence through it leaves it */
static void test_idle(void) {
  UDSP_FDN a, b;
  float *ma = fdn_new(&a), *mb = fdn_new(&b);
  enum { N = 3001 };
  static float zero[N * 2], wet[N * 2];
  udsp_fdn_stereo(&a, zero, wet, N);
  udsp_fdn_idle(&b, 1000);
  udsp_fdn_idle(&b, N - 1000);
  UDSP_TEST_ASSERT(a.mc == b.mc && a.ms == b.ms && a.t == b.t,
                   "idle 3001 frames: the modulation silence leaves");
  free(ma);
  free(mb);
}

int main(void) {
  printf("fdn_test\n");
  test_line();
  test_line_reads();
  test_decay();
  test_stable_and_blocks();
  test_predelay();
  test_idle();
  return udsp_test_summary();
}
