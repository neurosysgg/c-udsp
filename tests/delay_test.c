#include "udsp_test.h"

#include "delay.h"

static void test_impulse_repeats(void) {
  float ring[UDSP_DELAY_RING_FLOATS(8)];
  UDSP_DELAY d;
  udsp_delay_init(&d, ring, 8);
  float b[10 * 2] = {1.0f, -1.0f};
  udsp_delay_stereo(&d, b, 10, 3, 0.5f, 1.0f);
  UDSP_TEST_ASSERT(b[0] == 1.0f && b[1] == -1.0f, "the dry impulse");
  UDSP_TEST_ASSERT(b[6] == 1.0f && b[7] == -1.0f, "repeated 3 frames later");
  UDSP_TEST_ASSERT(b[12] == 0.5f && b[18] == 0.25f,
                   "each repeat fed back at 0.5");
  UDSP_TEST_ASSERT(b[2] == 0.0f && b[8] == 0.0f, "nothing in between");
}

static void test_split_and_wrap(void) {
  float r1[UDSP_DELAY_RING_FLOATS(5)], r2[UDSP_DELAY_RING_FLOATS(5)];
  UDSP_DELAY x, y;
  udsp_delay_init(&x, r1, 5);
  udsp_delay_init(&y, r2, 5);
  float a[23 * 2], b[23 * 2];
  for (int i = 0; i < 46; i++) {
    a[i] = b[i] = (float)(i % 5) - 2.0f;
  }
  udsp_delay_stereo(&x, a, 23, 4, 0.3f, 0.7f);
  udsp_delay_stereo(&y, b, 9, 4, 0.3f, 0.7f);
  udsp_delay_stereo(&y, b + 18, 14, 4, 0.3f, 0.7f);
  int same = 1;
  for (int i = 0; i < 46; i++) {
    same &= a[i] == b[i];
  }
  UDSP_TEST_ASSERT(same,
                   "split blocks, a ring wrapping 4 times: the same output");
}

static void test_length_clamped(void) {
  float ring[UDSP_DELAY_RING_FLOATS(4)];
  UDSP_DELAY d;
  udsp_delay_init(&d, ring, 4);
  float b[6 * 2] = {1.0f, 1.0f};
  udsp_delay_stereo(&d, b, 6, 100, 0.0f, 1.0f);
  UDSP_TEST_ASSERT(b[8] == 1.0f, "a length past the ring is the ring's");
}

int main(void) {
  test_impulse_repeats();
  test_split_and_wrap();
  test_length_clamped();
  return udsp_test_summary();
}
