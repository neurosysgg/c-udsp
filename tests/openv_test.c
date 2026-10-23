#include "udsp_test.h"

#include <math.h>
#include <string.h>

#include "synth/openv.h"

#define RATE 48000

static UDSP_OPENV_SHAPE shape(float a, float d, float s, float r) {
  UDSP_OPENV_SHAPE sh;
  memset(&sh, 0, sizeof sh);
  sh.peak = 1.0f;
  sh.sustain = s;
  sh.attack_ms = a;
  sh.decay_ms = d;
  sh.release_ms = r;
  return sh;
}

/* frames rendered until the stage is st */
static size_t until(UDSP_OPENV *e, UDSP_OPENV_STAGE st, size_t max) {
  float v;
  size_t n = 0;
  while (e->stage != st && n < max) {
    udsp_openv_render(e, &v, 1);
    n++;
  }
  return n;
}

static void test_segments(void) {
  UDSP_OPENV_SHAPE sh = shape(1.0f, 50.0f, 0.25f, 200.0f);
  sh.initial = 0.1f;
  sh.end = 0.05f;
  UDSP_OPENV e;
  udsp_openv_init(&e, &sh, RATE);
  float v[49];
  udsp_openv_render(&e, v, 1);
  UDSP_TEST_ASSERT(udsp_openv_idle(&e) && v[0] == 0.05f, "idle at its end");
  udsp_openv_gate(&e, 1);
  udsp_openv_render(&e, v, 49);
  UDSP_TEST_ASSERT(v[0] == 0.1f && fabsf(v[24] - 0.55f) < 1e-6f &&
                       v[48] == 1.0f && e.stage == UDSP_OPENV_DECAY,
                   "from the initial level to the peak in its 48 frames, a "
                   "line");
  size_t d = until(&e, UDSP_OPENV_SUSTAIN, 100000);
  /* its first frame was the 49th rendered: the peak */
  UDSP_TEST_ASSERT(d == 2399 && e.value == 0.25f,
                   "the decay its 50 ms to sustain exactly (%zu)", d);
  udsp_openv_render(&e, v, 3);
  UDSP_TEST_ASSERT(v[2] == 0.25f, "held");
  udsp_openv_gate(&e, 0);
  size_t r = until(&e, UDSP_OPENV_IDLE, 100000);
  UDSP_TEST_ASSERT(r == 9600 && e.value == 0.05f,
                   "the release its 200 ms to the end level (%zu)", r);
}

static void test_slopes(void) {
  float mid[3];
  for (int k = 0; k < 3; k++) {
    UDSP_OPENV_SHAPE sh = shape(10.0f, 10.0f, 1.0f, 10.0f);
    sh.attack_slope = (float)(k - 1);
    UDSP_OPENV e;
    udsp_openv_init(&e, &sh, RATE);
    udsp_openv_gate(&e, 1);
    float v[481];
    udsp_openv_render(&e, v, 481);
    mid[k] = v[240];
    UDSP_TEST_ASSERT(v[0] == 0.0f && v[480] == 1.0f,
                     "slope %d: from 0, on 1 at its end", k - 1);
  }
  /* an exponential of 2^6 between its ends: (2^3 - 1) / (2^6 - 1) */
  UDSP_TEST_ASSERT(mid[1] == 0.5f && fabsf(mid[2] - 7.0f / 63.0f) < 1e-6f &&
                       fabsf(mid[0] - 56.0f / 63.0f) < 1e-6f,
                   "half way: a line's half, 1 late (%g), -1 early (%g)",
                   mid[2], mid[0]);
}

static void test_modes(void) {
  UDSP_OPENV_SHAPE sh = shape(1.0f, 1.0f, 0.5f, 1.0f);
  sh.initial = 0.2f;
  sh.mode = UDSP_OPENV_LOOP;
  sh.loop_ms = 2.0f;
  UDSP_OPENV e;
  udsp_openv_init(&e, &sh, RATE);
  udsp_openv_gate(&e, 1);
  size_t to_return = until(&e, UDSP_OPENV_RETURN, 1000);
  size_t back = until(&e, UDSP_OPENV_ATTACK, 1000);
  UDSP_TEST_ASSERT(to_return == 96 && back == 96 && e.value == 0.2f,
                   "loop: attack and decay, back to the initial level over "
                   "its 2 ms, the attack again");
  udsp_openv_gate(&e, 0);
  UDSP_TEST_ASSERT(e.stage == UDSP_OPENV_RELEASE, "until the gate's off");

  sh.mode = UDSP_OPENV_TRIGGER;
  udsp_openv_init(&e, &sh, RATE);
  udsp_openv_gate(&e, 1);
  udsp_openv_gate(&e, 0);
  size_t all = until(&e, UDSP_OPENV_IDLE, 1000);
  UDSP_TEST_ASSERT(all == 144,
                   "trigger: the gate's end ignored, attack, "
                   "decay and release in their time (%zu)",
                   all);

  sh.mode = UDSP_OPENV_HOLD;
  udsp_openv_init(&e, &sh, RATE);
  udsp_openv_gate(&e, 1);
  until(&e, UDSP_OPENV_SUSTAIN, 1000);
  udsp_openv_restart(&e);
  float v[10];
  udsp_openv_render(&e, v, 2);
  UDSP_TEST_ASSERT(e.stage == UDSP_OPENV_ATTACK && v[0] == 0.5f && v[1] > 0.5f,
                   "restarted while held: the attack from where it is");
  udsp_openv_gate(&e, 0);
  udsp_openv_render(&e, v, 10);
  const float at = e.value;
  udsp_openv_gate(&e, 1);
  udsp_openv_render(&e, v, 1);
  UDSP_TEST_ASSERT(v[0] == at, "retriggered mid release: from where it is");
  udsp_openv_restart(&e);
  udsp_openv_gate(&e, 0);
  udsp_openv_restart(&e);
  UDSP_TEST_ASSERT(e.stage == UDSP_OPENV_RELEASE,
                   "a restart once released: ignored");
}

static void test_blocks(void) {
  UDSP_OPENV_SHAPE sh = shape(3.0f, 7.0f, 0.4f, 5.0f);
  sh.attack_slope = 0.7f;
  sh.decay_slope = -0.3f;
  sh.release_slope = 0.2f;
  sh.mode = UDSP_OPENV_LOOP;
  sh.loop_ms = 4.0f;
  enum { N = 2048 };
  float a[N], b[N];
  UDSP_OPENV x, y;
  udsp_openv_init(&x, &sh, RATE);
  udsp_openv_gate(&x, 1);
  y = x;
  udsp_openv_render(&x, a, N);
  for (size_t at = 0, len = 1; at < N; at += len, len = len * 3 + 1) {
    udsp_openv_render(&y, b + at, at + len > N ? N - at : len);
  }
  UDSP_TEST_ASSERT(memcmp(a, b, sizeof a) == 0,
                   "the same bits however the blocks are cut");
}

int main(void) {
  test_segments();
  test_slopes();
  test_modes();
  test_blocks();
  return udsp_test_summary();
}
