#include "udsp_test.h"

#include <math.h>
#include <string.h>

#include "adsr.h"

#define RATE 48000

static UDSP_ADSR make(float a, float d, float s, float r) {
  UDSP_ADSR e;
  memset(&e, 0, sizeof(e));
  udsp_adsr_set(&e, a, d, s, r, RATE);
  return e;
}

/* the first frame after out[0..n) where the stage is st */
static size_t frames_until(UDSP_ADSR *e, UDSP_ADSR_STAGE st, size_t max) {
  float v;
  size_t n = 0;
  while (e->stage != st && n < max) {
    udsp_adsr_render(e, &v, 1);
    n++;
  }
  return n;
}

static void test_attack_linear(void) {
  UDSP_ADSR e = make(1.0f, 100.0f, 0.5f, 100.0f);
  UDSP_TEST_ASSERT(udsp_adsr_idle(&e), "idle until gated");
  udsp_adsr_gate(&e, 1);
  float v[49];
  udsp_adsr_render(&e, v, 49);
  UDSP_TEST_ASSERT(v[0] == 0.0f && fabsf(v[24] - 0.5f) < 1e-6f && v[48] == 1.0f,
                   "1 ms at 48 kHz: 0, half way at 24 frames, 1 at 48");
  UDSP_TEST_ASSERT(e.stage == UDSP_ADSR_DECAY, "then the decay");
}

static void test_segments_take_their_time(void) {
  UDSP_ADSR e = make(0.0f, 50.0f, 0.0f, 200.0f);
  udsp_adsr_gate(&e, 1);
  size_t a = frames_until(&e, UDSP_ADSR_DECAY, 100000);
  size_t d = frames_until(&e, UDSP_ADSR_SUSTAIN, 100000);
  UDSP_TEST_ASSERT(a == 1, "attack 0 is a frame");
  UDSP_TEST_ASSERT(d >= 2399 && d <= 2401,
                   "a decay from 1 to 0 takes its 50 ms (%zu frames)", d);
  e.value = 1.0f;
  udsp_adsr_gate(&e, 0);
  size_t r = frames_until(&e, UDSP_ADSR_IDLE, 100000);
  UDSP_TEST_ASSERT(r >= 9599 && r <= 9601 && e.value == 0.0f,
                   "a release from 1 takes its 200 ms, ends on 0 (%zu)", r);
}

static void test_sustain_and_release_from_anywhere(void) {
  UDSP_ADSR e = make(10.0f, 10.0f, 0.6f, 10.0f);
  udsp_adsr_gate(&e, 1);
  frames_until(&e, UDSP_ADSR_SUSTAIN, 100000);
  float v[4];
  udsp_adsr_render(&e, v, 4);
  UDSP_TEST_ASSERT(v[3] == 0.6f, "held at the sustain level");
  udsp_adsr_set(&e, 10.0f, 10.0f, 0.3f, 10.0f, RATE);
  udsp_adsr_render(&e, v, 2);
  UDSP_TEST_ASSERT(v[1] == 0.3f, "a sustain moved while held follows");
  /* off during the attack: down from there, no jump */
  UDSP_ADSR f = make(10.0f, 10.0f, 0.6f, 10.0f);
  udsp_adsr_gate(&f, 1);
  float w[200];
  udsp_adsr_render(&f, w, 200);
  udsp_adsr_gate(&f, 0);
  udsp_adsr_render(&f, v, 2);
  UDSP_TEST_ASSERT(v[0] == f.value / f.release / f.release && v[1] < v[0] &&
                       v[0] > w[199],
                   "released mid-attack: falls from where it was");
  /* a retrigger mid-release attacks from there */
  udsp_adsr_render(&f, w, 100);
  float at = f.value;
  udsp_adsr_gate(&f, 1);
  udsp_adsr_render(&f, v, 2);
  UDSP_TEST_ASSERT(v[0] == at && v[1] > at && at > 0.0f,
                   "retriggered mid-release: rises from where it was");
}

int main(void) {
  test_attack_linear();
  test_segments_take_their_time();
  test_sustain_and_release_from_anywhere();
  return udsp_test_summary();
}
