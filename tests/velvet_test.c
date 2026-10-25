#include "udsp_test.h"

#include <string.h>

#include "fx/velvet.h"

#define N 4096

static float mem[8192];
static float x[N], got[N], want[N];

/* x through v directly, a frame at a time: the reference */
static void direct(const UDSP_VELVET *v) {
  for (int i = 0; i < N; i++) {
    float y = 0.0f;
    for (uint32_t t = 0; t < v->n; t++) {
      y += i >= (int)v->at[t] ? v->g[t] * x[i - (int)v->at[t]] : 0.0f;
    }
    want[i] = y;
  }
}

/* x written and read in blocks of k */
static void blocks(const UDSP_VELVET *v, size_t k) {
  UDSP_LINE l;
  udsp_line_init(&l, mem, udsp_line_size(1024 + 512));
  for (size_t at = 0; at < N; at += k) {
    const size_t m = N - at < k ? N - at : k;
    udsp_line_write(&l, x + at, m);
    udsp_velvet_read(v, &l, got + at, m);
  }
}

/* a cell one frame wide holds its tap exactly */
static void test_one(void) {
  UDSP_VELVET v;
  const UDSP_VELVET_SEG s = {10, 11, 1, 0.5f};
  udsp_velvet_make(&v, &s, 1, 7);
  UDSP_TEST_ASSERT(
      v.n == 1 && v.at[0] == 10 && (v.g[0] == 0.5f || v.g[0] == -0.5f),
      "one tap at 10, +-0.5 (%u at %u, %g)", v.n, v.at[0], (double)v.g[0]);
}

/* pinned: a seed's taps, each in its cell */
static void test_pinned(void) {
  static const uint32_t AT[] = {19, 36, 58, 82, 129, 299, 342};
  static const float G[] = {-1, -1, 1, 1, 0.5f, -0.5f, 0.5f};
  UDSP_VELVET v, w;
  const UDSP_VELVET_SEG s[2] = {{0, 100, 4, 1.0f}, {100, 400, 3, 0.5f}};
  udsp_velvet_make(&v, s, 2, 1);
  int same = v.n == 7;
  for (uint32_t i = 0; same && i < 7; i++) {
    same = v.at[i] == AT[i] && v.g[i] == G[i];
  }
  UDSP_TEST_ASSERT(same, "seed 1: its seven taps (%u)", v.n);
  udsp_velvet_make(&w, s, 2, 1);
  same = w.n == v.n;
  for (uint32_t i = 0; same && i < v.n; i++) {
    same = w.at[i] == v.at[i] && w.g[i] == v.g[i];
  }
  UDSP_TEST_ASSERT(same, "the same seed, the same taps");
  udsp_velvet_make(&w, s, 2, 2);
  UDSP_TEST_ASSERT(memcmp(v.at, w.at, sizeof v.at) != 0,
                   "another seed, other places");
}

/* every tap in its own cell, ascending; the energy the segments' */
static void test_cells(void) {
  const UDSP_VELVET_SEG s[4] = {{0, 6, 3, 0.2125f},
                                {6, 25, 3, 0.1375f},
                                {25, 92, 3, 0.0875f},
                                {92, 288, 3, 0.05f}};
  int bad = 0;
  for (uint32_t seed = 1; seed < 200; seed++) {
    UDSP_VELVET v;
    udsp_velvet_make(&v, s, 4, seed);
    double e = 0.0;
    for (uint32_t i = 0; i < v.n; i++) {
      const UDSP_VELVET_SEG *g = &s[i / 3];
      const uint32_t c = i % 3, span = g->to - g->from;
      const uint32_t lo = g->from + c * span / 3;
      const uint32_t hi = g->from + (c + 1) * span / 3;
      bad |= v.at[i] < lo || v.at[i] >= (hi > lo ? hi : lo + 1);
      bad |= i && v.at[i] < v.at[i - 1];
      bad |= v.g[i] != g->gain && v.g[i] != -g->gain;
      e += (double)v.g[i] * v.g[i];
    }
    bad |= v.n != 12 || e < 0.22265 || e > 0.22266;
  }
  UDSP_TEST_ASSERT(!bad, "199 seeds: each tap in its cell, ascending, its "
                         "segment's gain; the energy 3 sum g^2");
}

/* more taps than it holds: the first UDSP_VELVET_MAX */
static void test_cap(void) {
  UDSP_VELVET v;
  const UDSP_VELVET_SEG s[2] = {{0, 1000, 50, 1.0f}, {1000, 2000, 50, 0.5f}};
  udsp_velvet_make(&v, s, 2, 3);
  UDSP_TEST_ASSERT(v.n == UDSP_VELVET_MAX && v.g[49] * v.g[49] == 1.0f &&
                       v.g[50] * v.g[50] == 0.25f,
                   "100 asked, %u kept, the second segment's first 14", v.n);
}

/* any block size: the direct convolution's sums */
static void test_blocks(void) {
  UDSP_VELVET v;
  const UDSP_VELVET_SEG s[3] = {
      {0, 4, 2, 0.5f}, {4, 100, 6, 0.25f}, {100, 1000, 10, 0.125f}};
  udsp_velvet_make(&v, s, 3, 11);
  uint32_t r = 5;
  for (int i = 0; i < N; i++) {
    r = r * 1664525u + 1013904223u;
    x[i] = (float)(int32_t)r / 2147483648.0f;
  }
  direct(&v);
  static const size_t K[] = {1, 7, 64, 512};
  float worst = 0.0f;
  for (size_t k = 0; k < 4; k++) {
    memset(got, 0, sizeof got);
    blocks(&v, K[k]);
    for (int i = 0; i < N; i++) {
      const float d = got[i] > want[i] ? got[i] - want[i] : want[i] - got[i];
      worst = d > worst ? d : worst;
    }
  }
  UDSP_TEST_ASSERT(worst < 1e-6f,
                   "blocks of 1, 7, 64, 512: the direct sums (%g)",
                   (double)worst);
}

int main(void) {
  printf("velvet_test\n");
  test_one();
  test_pinned();
  test_cells();
  test_cap();
  test_blocks();
  return udsp_test_summary();
}
