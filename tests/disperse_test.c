#include "udsp_test.h"

#include <math.h>
#include <string.h>

#include "disperse.h"

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
  UDSP_TEST_ASSERT(worst < 0.01f, "64 stages: every sine's level kept (%.4f dB)",
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
  UDSP_TEST_ASSERT(at > 10.0f * (far > 1.0f ? far : 1.0f) && fabsf(none) < 0.5f,
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

int main(void) {
  printf("disperse_test\n");
  test_level();
  test_smear();
  test_fan();
  return udsp_test_summary();
}
