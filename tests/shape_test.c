#include "udsp_test.h"

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

int main(void) {
  test_softclip_values();
  test_vector_and_tail_agree();
  return udsp_test_summary();
}
