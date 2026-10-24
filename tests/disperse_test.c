#include "udsp_test.h"

#include <math.h>
#include <string.h>

#include "fx/disperse.h"

#define RATE 48000
#define N RATE

static float x[2 * N];

static void sine(float hz) {
  for (int i = 0; i < N; i++) {
    x[2 * i] = x[2 * i + 1] = (float)sin(6.283185307179586 * hz * i / RATE);
  }
}

static float amp(void) {
  double sq = 0.0;
  for (int i = N / 2; i < N; i++) {
    sq += (double)x[2 * i] * x[2 * i];
  }
  return (float)sqrt(2.0 * sq / (N / 2));
}

/* an allpass: a sine's level unchanged at any frequency (whole stages:
 * a fractional one is a mix of a stage's in and out, a passing notch) */
static void test_level(void) {
  static const float HZ[] = {50, 300, 1000, 3000, 12000};
  float worst = 0.0f;
  for (size_t h = 0; h < 5; h++) {
    UDSP_DISPERSE d;
    udsp_disperse_reset(&d);
    udsp_disperse_set(&d, 64, 1000.0f, 2.0f, 2.0f, RATE);
    sine(HZ[h]);
    udsp_disperse_stereo(&d, x, N);
    float e = fabsf(20.0f * log10f(amp()));
    worst = e > worst ? e : worst;
  }
  UDSP_TEST_ASSERT(worst < 0.01f,
                   "64 stages: every sine's level kept (%.4f dB)",
                   (double)worst);
}

/* the smear: an impulse's energy arrives later near the centre than
 * far above it (group delay peaks around the stages' frequency) */
static float delay_at(float hz, float amount) {
  UDSP_DISPERSE d;
  udsp_disperse_reset(&d);
  udsp_disperse_set(&d, amount, 1000.0f, 4.0f, 0.0f, RATE);
  /* a sine burst's envelope: where its energy's centre lands */
  memset(x, 0, sizeof x);
  const int len = (int)(RATE * 20 / hz);
  for (int i = 0; i < len; i++) {
    float w = 0.5f - 0.5f * cosf(6.2831853f * i / len);
    x[2 * i] = x[2 * i + 1] = w * sinf(6.2831853f * hz * i / RATE);
  }
  udsp_disperse_stereo(&d, x, N);
  double e = 0.0, t = 0.0;
  for (int i = 0; i < N; i++) {
    e += (double)x[2 * i] * x[2 * i];
    t += (double)i * x[2 * i] * x[2 * i];
  }
  return (float)(t / e) - len / 2.0f;
}

static void test_smear(void) {
  const float at = delay_at(1000.0f, 32), far = delay_at(8000.0f, 32);
  const float none = delay_at(1000.0f, 0);
  UDSP_TEST_ASSERT(
      at > 10.0f * (far > 1.0f ? far : 1.0f) && fabsf(none) < 0.5f,
      "32 stages at 1 kHz: the centre %.0f frames late, 8 kHz %.1f "
      "(none: %.2f)",
      (double)at, (double)far, (double)none);
  const float more = delay_at(1000.0f, 64);
  UDSP_TEST_ASSERT(fabsf(more - 2.0f * at) < 0.05f * at,
                   "twice the stages: twice the delay (%.0f)", (double)more);
}

/* more stages don't move those already in: 8 stages, then 8 of 16 */
static void test_fan(void) {
  UDSP_DISPERSE a, b;
  udsp_disperse_set(&a, 8, 1000.0f, 1.0f, 3.0f, RATE);
  udsp_disperse_set(&b, 16, 1000.0f, 1.0f, 3.0f, RATE);
  UDSP_TEST_ASSERT(!memcmp(a.st, b.st, 8 * sizeof a.st[0]),
                   "the first 8 stages the same with 16");
}

/* the output the same bits however the blocks are cut: the stages run
 * a wavefront in vector lanes, whose first and last steps of a block
 * leave the lanes without a frame alone (blocks shorter than it too) */
static void test_blocks(void) {
  static float whole[2 * 4096], cut[2 * 4096];
  for (int i = 0; i < 2 * 4096; i++) {
    whole[i] = (float)sin(0.37 * i) * (i % 7 == 0 ? -1.0f : 0.5f);
  }
  memcpy(cut, whole, sizeof whole);
  static const size_t BLOCK[] = {1, 3, 2, 17, 5, 64, 31, 1, 250};
  static const float AMOUNT[] = {64, 40.5f, 17, 7};
  for (int a = 0; a < 4; a++) {
    static UDSP_DISPERSE d1, d2;
    udsp_disperse_reset(&d1);
    udsp_disperse_reset(&d2);
    udsp_disperse_set(&d1, AMOUNT[a], 700.0f, 2.0f, 2.0f, RATE);
    udsp_disperse_set(&d2, AMOUNT[a], 700.0f, 2.0f, 2.0f, RATE);
    udsp_disperse_stereo(&d1, whole, 4096);
    for (size_t at = 0, b = 0; at < 4096; b++) {
      size_t k = BLOCK[b % 9];
      k = k > 4096 - at ? 4096 - at : k;
      udsp_disperse_stereo(&d2, cut + 2 * at, k);
      at += k;
    }
    UDSP_TEST_ASSERT(!memcmp(whole, cut, sizeof whole),
                     "%.1f stages: blocks of 1 .. 250 the bits of one block",
                     (double)AMOUNT[a]);
  }
}

/* every stage count the bits of a stage at a time: the wavefront's lanes
 * past the last stage (a remainder short of a whole wavefront, padded)
 * pass their input on and keep no state. The amount moves down and up
 * between blocks, so stages leave (their state at rest) and come back. */
static float ref_ic[UDSP_DISPERSE_MAX][2][2];

static void ref_run(const UDSP_DISPERSE *d, float *io, size_t frames) {
  const int n = d->whole + (d->frac > 0.0f && d->whole < UDSP_DISPERSE_MAX);
  for (int i = 0; i < n; i++) {
    const UDSP_DISPERSE_STAGE *s = &d->st[i];
    const float w = (i < d->whole ? 1.0f : d->frac) * 2.0f * d->k;
    for (size_t f = 0; f < frames; f++) {
      for (int c = 0; c < 2; c++) {
        float *ic = ref_ic[i][c];
        const float in = io[2 * f + c];
        const float v3 = in - ic[1];
        const float v1 = s->a1 * ic[0] + s->a2 * v3;
        const float v2 = ic[1] + s->a2 * ic[0] + s->a3 * v3;
        ic[0] = 2.0f * v1 - ic[0];
        ic[1] = 2.0f * v2 - ic[1];
        io[2 * f + c] = in - w * v1;
      }
    }
  }
  memset(ref_ic[n], 0, (size_t)(UDSP_DISPERSE_MAX - n) * sizeof ref_ic[0]);
}

static void test_reference(void) {
  enum { F = 1200 };
  static float in[2 * F], got[2 * F], want[2 * F];
  for (int i = 0; i < 2 * F; i++) {
    in[i] = (float)sin(0.21 * i + (i & 1)) * (i % 5 == 0 ? -0.7f : 0.4f);
  }
  static const size_t BLOCK[] = {300, 7, 64, 1, 33, 256, 2, 129};
  float bad = 0.0f; /* the first amount that differs */
  for (int a = 1; a <= 2 * 66 && bad == 0.0f; a++) {
    /* whole and half amounts, and between blocks a few stages less, then
     * more (the last block's back at the amount) */
    const float amount = 0.5f * (float)a;
    static UDSP_DISPERSE d;
    udsp_disperse_reset(&d);
    memset(ref_ic, 0, sizeof ref_ic);
    memcpy(got, in, sizeof in);
    memcpy(want, in, sizeof in);
    for (size_t at = 0, b = 0; at < F; b++) {
      static const float MOVE[] = {0, -5, 3, 0, -1.5f, 9, -12, 0};
      const float m = amount + MOVE[b % 8];
      udsp_disperse_set(&d, m < 0.5f ? 0.5f : m, 700.0f, 1.5f, 2.5f, RATE);
      size_t k = BLOCK[b % 8];
      k = k > F - at ? F - at : k;
      udsp_disperse_stereo(&d, got + 2 * at, k);
      ref_run(&d, want + 2 * at, k);
      at += k;
    }
    if (memcmp(got, want, sizeof got) || memcmp(d.ic, ref_ic, sizeof ref_ic)) {
      bad = amount;
    }
  }
  UDSP_TEST_ASSERT(bad == 0.0f,
                   "0.5 .. 66 stages, moving: a stage at a time's bits "
                   "(first differing: %.1f)",
                   (double)bad);
}

int main(void) {
  printf("disperse_test\n");
  test_blocks();
  test_reference();
  test_level();
  test_smear();
  test_fan();
  return udsp_test_summary();
}
