#include "udsp_test.h"

#include <math.h>
#include <stdint.h>
#include <string.h>

#include "fx/shape.h"

static void test_softclip_values(void) {
  float b[11] = {0.0f,  1.0f, -1.0f, 100.0f, -100.0f, 0.5f,
                 0.25f, 2.0f, 4.0f,  8.0f,   0.1f};
  udsp_softclip(b, 11, 4.0f);
  UDSP_TEST_ASSERT(b[0] == 0.0f, "0 stays 0");
  UDSP_TEST_ASSERT(b[1] == 0.8f && b[2] == -0.8f,
                   "drive 4: 1 -> 4/5, odd-symmetric");
  UDSP_TEST_ASSERT(b[3] > 0.99f && b[3] < 1.0f && b[4] < -0.99f,
                   "towards +-1 in the limit");
  UDSP_TEST_ASSERT(b[10] < b[6] && b[6] < b[5] && b[5] < b[1] && b[1] < b[7] &&
                       b[7] < b[8] && b[8] < b[9],
                   "monotonic");
}

static void test_vector_and_tail_agree(void) {
  float a[37], b[37];
  for (int i = 0; i < 37; i++) {
    a[i] = (float)(i - 18) * 0.11f;
    b[36 - i] = a[i];
  }
  udsp_softclip(a, 37, 3.0f);
  udsp_softclip(b, 37, 3.0f);
  int same = 1;
  for (int i = 0; i < 37; i++) {
    same &= a[i] == b[36 - i];
  }
  UDSP_TEST_ASSERT(same, "a sample's result doesn't depend on where it sits");
}

static void test_shapes(void) {
  int small = 1;
  for (int k = 0; k < UDSP_SHAPES; k++) {
    float y = udsp_shape1((UDSP_SHAPE)k, 1e-3f);
    small &= k == UDSP_SHAPE_ROOT || k == UDSP_SHAPE_QUINTIC ||
             fabsf(y - 1e-3f) < 1e-5f;
  }
  UDSP_TEST_ASSERT(small, "every curve but root and quintic passes a small "
                          "signal as it is");
  float worst = 0.0f;
  int odd = 1, capped = 1;
  for (int i = 0; i <= 9000; i++) {
    const float x = (float)i * 1e-3f;
    const float y = udsp_shape1(UDSP_SHAPE_TANH, x);
    const float e = fabsf(y - (float)tanh((double)x));
    worst = e > worst ? e : worst;
    odd &= udsp_shape1(UDSP_SHAPE_TANH, -x) == -y;
    capped &= y <= 1.0f;
  }
  UDSP_TEST_ASSERT(worst < 1.1e-4f && odd && capped,
                   "tanh: within 1e-4 of tanh to 9 (%.2e), odd, never past 1",
                   (double)worst);
  UDSP_TEST_ASSERT(fabsf(udsp_shape1(UDSP_SHAPE_TANH, 1.0f) - 0.7615942f) <
                           2e-6f &&
                       udsp_shape1(UDSP_SHAPE_TANH, 4.97f) == 1.0f,
                   "tanh: 0.76159 at 1 (Kilohearts' Overdrive), 1 from 4.97");
  /* Kilohearts' Saturate, measured through the plugin (DC steps) */
  static const float SAT[][2] = {{0.1f, 0.18625f},
                                 {0.3f, 0.52966f},
                                 {0.5f, 0.79297f},
                                 {0.7f, 0.94678f},
                                 {0.9f, 0.99768f}};
  int sat = 1;
  for (size_t i = 0; i < sizeof SAT / sizeof *SAT; i++) {
    sat &=
        fabsf(udsp_shape1(UDSP_SHAPE_QUINTIC, SAT[i][0]) - SAT[i][1]) < 2e-5f;
  }
  UDSP_TEST_ASSERT(
      sat && udsp_shape1(UDSP_SHAPE_QUINTIC, 1.0f) == 1.0f &&
          udsp_shape1(UDSP_SHAPE_QUINTIC, -7.0f) == -1.0f &&
          fabsf(udsp_shape1(UDSP_SHAPE_QUINTIC, 1e-3f) - 1.875e-3f) < 1e-7f,
      "quintic: Kilohearts' Saturate's values, 1 from 1, "
      "slope 15 / 8 at 0");
  UDSP_TEST_ASSERT(udsp_shape1(UDSP_SHAPE_SOFT, 3.0f) == 1.0f &&
                       udsp_shape1(UDSP_SHAPE_SOFT, 10.0f) == 1.0f &&
                       fabsf(udsp_shape1(UDSP_SHAPE_SOFT, 1.0f) - 0.7778f) <
                           1e-3f,
                   "soft: 1 at 3 and past it, tanh-ish on the way (0.78 at 1)");
  UDSP_TEST_ASSERT(udsp_shape1(UDSP_SHAPE_HARD, -2.0f) == -1.0f &&
                       udsp_shape1(UDSP_SHAPE_HARD, 0.5f) == 0.5f,
                   "hard: clipped at 1, untouched under");
  UDSP_TEST_ASSERT(
      fabsf(udsp_shape1(UDSP_SHAPE_FOLD, 1.5707963f) - 1.0f) < 1e-5f &&
          fabsf(udsp_shape1(UDSP_SHAPE_FOLD, 3.1415927f)) < 1e-5f &&
          fabsf(udsp_shape1(UDSP_SHAPE_FOLD, -1.5707963f) + 1.0f) < 1e-5f,
      "fold: 1 at pi / 2, folded back to 0 at pi, odd");
  UDSP_TEST_ASSERT(
      fabsf(udsp_shape1(UDSP_SHAPE_TUBE, 1.0f) - 0.5f) < 1e-6f &&
          fabsf(udsp_shape1(UDSP_SHAPE_TUBE, -1.0f) + 1.0f / 3.0f) < 1e-6f,
      "tube: the halves flatten unevenly (0.5 and -0.33 at 1)");
  UDSP_TEST_ASSERT(udsp_shape1(UDSP_SHAPE_ROOT, 0.25f) == 0.5f &&
                       udsp_shape1(UDSP_SHAPE_ROOT, -4.0f) == -2.0f &&
                       udsp_shape1(UDSP_SHAPE_ROOT, 1.0f) == 1.0f,
                   "root: sign(x) sqrt(|x|) (0.5 at 0.25, -2 at -4, 1 at 1)");
  const float pz = udsp_shape1(UDSP_SHAPE_ROOT, 0.0f);
  const float nz = udsp_shape1(UDSP_SHAPE_ROOT, -0.0f);
  UDSP_TEST_ASSERT(pz == 0.0f && !signbit(pz) && signbit(nz),
                   "root: +0 and -0 keep their signs");
}

/* udsp_shape_mix the bits of its expression a sample at a time, every
 * curve: the curves' edges (+-1, +-3, 0 and -0, the turns' wrap), large
 * and random values, odd lengths (the vector body and its tail) */
static void test_shape_mix(void) {
  enum { N = 1031 };
  static float in[N], got[N], want[N];
  static const float EDGE[] = {
      0.0f,   -0.0f,  1.0f,    -1.0f,  3.0f,    -3.0f,      2.99f,
      -2.99f, 1e-30f, -1e-30f, 500.0f, -500.0f, 6.2831855f, -6.2831855f,
      0.5f,   4.97f,  -4.97f,  4.96f,  -4.96f,  1.0001f,    -0.9999f};
  uint32_t r = 12345;
  for (int i = 0; i < N; i++) {
    r = r * 1664525u + 1013904223u;
    in[i] = i < 21 ? EDGE[i] : ((float)(r >> 8) / 16777216.0f - 0.5f) * 9.0f;
  }
  static const float DRIVE[] = {1.0f, 3.7f, 0.25f};
  static const size_t LEN[] = {N, 7, 256, 257};
  int same = 1;
  for (int k = 0; k < UDSP_SHAPES; k++) {
    for (int d = 0; d < 3; d++) {
      for (int l = 0; l < 4; l++) {
        const float g = DRIVE[d], m = 0.3f * (float)(d + 1), o = 0.9f;
        memcpy(got, in, sizeof in);
        udsp_shape_mix(got, LEN[l], (UDSP_SHAPE)k, g, m, o);
        for (size_t i = 0; i < LEN[l]; i++) {
          const float x = in[i];
          want[i] =
              (x * (1.0f - m) + udsp_shape1((UDSP_SHAPE)k, x * g) * m) * o;
        }
        same &= !memcmp(got, want, LEN[l] * sizeof(float));
      }
    }
  }
  UDSP_TEST_ASSERT(same, "shape_mix: every curve the bits of a sample at a "
                         "time (edges, -0, odd lengths)");
}

int main(void) {
  test_softclip_values();
  test_vector_and_tail_agree();
  test_shapes();
  test_shape_mix();
  return udsp_test_summary();
}
