#include "udsp_test.h"

#include <math.h>

#include "shape.h"

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
    small &= fabsf(y - 1e-3f) < 1e-5f;
  }
  UDSP_TEST_ASSERT(small, "every curve passes a small signal as it is");
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
  UDSP_TEST_ASSERT(fabsf(udsp_shape1(UDSP_SHAPE_TUBE, 1.0f) - 0.5f) < 1e-6f &&
                       fabsf(udsp_shape1(UDSP_SHAPE_TUBE, -1.0f) +
                             1.0f / 3.0f) < 1e-6f,
                   "tube: the halves flatten unevenly (0.5 and -0.33 at 1)");
}

int main(void) {
  test_softclip_values();
  test_vector_and_tail_agree();
  test_shapes();
  return udsp_test_summary();
}
