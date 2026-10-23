#include "udsp_test.h"

#include "synth/smooth.h"

static void test_glide(void) {
  UDSP_SMOOTH s;
  udsp_smooth_init(&s, 0.0f, 4);
  UDSP_TEST_ASSERT(udsp_smooth_settled(&s) && udsp_smooth_next(&s) == 0.0f,
                   "settled at its start value");
  udsp_smooth_to(&s, 1.0f);
  float v[6];
  udsp_smooth_render(&s, v, 6);
  UDSP_TEST_ASSERT(v[0] == 0.0f && v[1] == 0.25f && v[2] == 0.5f &&
                       v[3] == 0.75f && v[4] == 1.0f && v[5] == 1.0f &&
                       udsp_smooth_settled(&s),
                   "four frames to the target, then held there");
}

static void test_lands_exactly(void) {
  UDSP_SMOOTH s;
  udsp_smooth_init(&s, 0.1f, 7);
  udsp_smooth_to(&s, 0.3f);
  float v[7], last = 0.0f;
  udsp_smooth_render(&s, v, 7);
  int rising = 1;
  for (int i = 1; i < 7; i++) {
    rising &= v[i] > v[i - 1];
  }
  last = udsp_smooth_next(&s);
  UDSP_TEST_ASSERT(rising && last == 0.3f,
                   "a step that doesn't divide evenly still lands on 0.3f");
}

static void test_retarget(void) {
  UDSP_SMOOTH s;
  udsp_smooth_init(&s, 0.0f, 10);
  udsp_smooth_to(&s, 1.0f);
  udsp_smooth_skip(&s, 5);
  udsp_smooth_to(&s, 0.0f);
  float v[11];
  udsp_smooth_render(&s, v, 11);
  UDSP_TEST_ASSERT(v[0] == 0.5f && v[1] < 0.5f && v[10] == 0.0f,
                   "a new target mid-glide: back from 0.5 in 10 frames");
  udsp_smooth_init(&s, 1.0f, 0);
  udsp_smooth_to(&s, 2.0f);
  UDSP_TEST_ASSERT(udsp_smooth_next(&s) == 2.0f, "no length: it jumps");
  udsp_smooth_init(&s, 0.0f, 8);
  udsp_smooth_to(&s, 1.0f);
  udsp_smooth_skip(&s, 3);
  UDSP_SMOOTH t;
  udsp_smooth_init(&t, 0.0f, 8);
  udsp_smooth_to(&t, 1.0f);
  for (int i = 0; i < 3; i++) {
    udsp_smooth_next(&t);
  }
  UDSP_TEST_ASSERT(s.left == t.left && s.value == t.value,
                   "skipping 3 is stepping 3 (%g vs %g)", (double)s.value,
                   (double)t.value);
}

int main(void) {
  test_glide();
  test_lands_exactly();
  test_retarget();
  return udsp_test_summary();
}
