#include "udsp_test.h"

#include <math.h>
#include <string.h>

#include "fx/bank.h"

#define RATE 48000
#define PI 3.141592653589793

static UDSP_BANK bank;
static float out[RATE * 2 / 4 * UDSP_BANK_MAX];

/* band 0's steady gain at hz: a sine through it, its peak past the
 * first half */
static double gain_at(double hz) {
  enum { F = RATE / 4 };
  static float x[F];
  for (int i = 0; i < F; i++) {
    x[i] = (float)sin(2 * PI * hz * i / RATE);
  }
  udsp_bank_reset(&bank);
  udsp_bank_run(&bank, x, 1, out, F);
  double peak = 0;
  for (int i = F / 2; i < F; i++) {
    peak = fmax(peak, fabs(out[(size_t)i * bank.lanes]));
  }
  return peak;
}

/* every order: unity at the centre, -3 dB at the edges of its width
 * (centre / q, geometric about the centre), far below an octave off */
static void test_bands(void) {
  for (uint32_t order = 1; order <= UDSP_BANK_ORDER; order++) {
    const float hz[3] = {1000, 2000, 4000}, q[3] = {4, 4, 4};
    memset(&bank, 0, sizeof bank);
    udsp_bank_set(&bank, 3, hz, q, NULL, order, RATE);
    const double c = gain_at(1000);
    /* edges f1 f2: f2 - f1 = 250, f1 f2 = 1000^2 */
    const double f2 = 125 + sqrt(125.0 * 125 + 1e6), f1 = f2 - 250;
    const double lo = 20 * log10(gain_at(f1)), hi = 20 * log10(gain_at(f2)),
                 oct = 20 * log10(gain_at(2000));
    UDSP_TEST_ASSERT(fabs(c - 1) < 0.01 && fabs(lo + 3.01) < 0.15 &&
                         fabs(hi + 3.01) < 0.15 &&
                         oct < -15.0 - 6.5 * (order - 1),
                     "order %u: centre %.3f, edges %.2f / %.2f dB, an octave "
                     "up %.1f dB",
                     order, c, lo, hi, oct);
  }
  UDSP_TEST_ASSERT(out[bank.lanes - 1] == 0.0f, "the lanes past n zero");
}

/* a band's gain, and a retune with n and order the same keeps the state
 * (no click: the output carries on from where it was) */
static void test_gain_retune(void) {
  const float hz[1] = {500}, q[1] = {2}, g[1] = {0.5f};
  memset(&bank, 0, sizeof bank);
  udsp_bank_set(&bank, 1, hz, q, g, 2, RATE);
  const double c = gain_at(500);
  UDSP_TEST_ASSERT(fabs(c - 0.5) < 0.005, "gain 0.5: %.4f", c);
  const float hz2[1] = {510};
  udsp_bank_set(&bank, 1, hz2, q, g, 2, RATE);
  float s = 0;
  for (int i = 0; i < UDSP_BANK_ORDER; i++) {
    s += fabsf(bank.ic1[i][0]) + fabsf(bank.ic2[i][0]);
  }
  UDSP_TEST_ASSERT(s > 0, "a retune keeps the state");
  udsp_bank_set(&bank, 1, hz2, q, g, 3, RATE);
  s = 0;
  for (int i = 0; i < UDSP_BANK_ORDER; i++) {
    s += fabsf(bank.ic1[i][0]) + fabsf(bank.ic2[i][0]);
  }
  UDSP_TEST_ASSERT(s == 0, "a new order zeroes it");
}

/* a sine's envelope settles at its amplitude with a fast attack and a
 * slow release (within the ripple that leaves), rises by attack, holds,
 * then falls by release */
static void test_follow(void) {
  enum { F = RATE / 2, L = UDSP_BANK_MAX };
  static float io[(size_t)F * L];
  UDSP_FOLLOW f;
  memset(&f, 0, sizeof f);
  udsp_follow_set(&f, L, 48, 4800, 0);
  for (int i = 0; i < F; i++) {
    io[(size_t)i * L] = (float)(0.5 * sin(2 * PI * 1000 * i / RATE));
  }
  udsp_follow_run(&f, io, F);
  double lo = 1, hi = 0;
  for (int i = F / 2; i < F; i++) {
    lo = fmin(lo, io[(size_t)i * L]);
    hi = fmax(hi, io[(size_t)i * L]);
  }
  UDSP_TEST_ASSERT(lo > 0.48 && hi < 0.501 && io[(size_t)(F - 1) * L + 1] == 0,
                   "a 0.5 sine's envelope %.3f .. %.3f, the other lanes 0", lo,
                   hi);
  /* a step: 1 - 1/e after attack frames; held 100 frames, then 1/e of
   * the way down after release */
  memset(&f, 0, sizeof f);
  udsp_follow_set(&f, L, 100, 1000, 100);
  memset(io, 0, sizeof io);
  for (int i = 0; i < 2000; i++) {
    io[(size_t)i * L] = 1.0f;
  }
  udsp_follow_run(&f, io, 6000);
  const float up = io[(size_t)99 * L], top = io[(size_t)1999 * L],
              hold = io[(size_t)2099 * L], down = io[(size_t)3100 * L];
  UDSP_TEST_ASSERT(fabsf(up - 0.632f) < 0.01f && top > 0.999f && hold == top &&
                       fabsf(down - 0.368f * top) < 0.01f,
                   "attack %.3f at 100 frames, top %.4f, held %.4f, release "
                   "%.3f a 1000 frames on",
                   (double)up, (double)top, (double)hold, (double)down);
}

/* white noise through a band: its power over the input's 2 B / rate, B
 * the noise width, within 5 % at every order */
static void test_noise_width(void) {
  enum { F = RATE * 4 };
  static float x[F];
  static float y[(size_t)F * 8];
  uint32_t seed = 1;
  double vin = 0;
  for (int i = 0; i < F; i++) {
    seed ^= seed << 13;
    seed ^= seed >> 17;
    seed ^= seed << 5;
    x[i] = (float)(int32_t)seed / 2147483648.0f;
    vin += (double)x[i] * x[i];
  }
  vin /= F;
  double worst = 0;
  for (uint32_t order = 1; order <= UDSP_BANK_ORDER; order++) {
    const float hz[1] = {2000}, q[1] = {6};
    memset(&bank, 0, sizeof bank);
    udsp_bank_set(&bank, 1, hz, q, NULL, order, RATE);
    udsp_bank_run(&bank, x, 1, y, F);
    /* past the band's settling, a tenth of a second */
    enum { FROM = RATE / 10 };
    double v = 0;
    for (int i = FROM; i < F; i++) {
      v += (double)y[(size_t)i * bank.lanes] * y[(size_t)i * bank.lanes];
    }
    v /= F - FROM;
    const double want = vin * 2 * udsp_bank_noise_width(2000, 6, order) / RATE;
    worst = fmax(worst, fabs(v / want - 1));
  }
  UDSP_TEST_ASSERT(worst < 0.05, "noise through a band: its width within %.3f",
                   worst);
}

/* the response against the bank's own: its magnitude at the centre, an
 * edge and an octave off, every order (the sampled peak reads a little
 * low: 0.06 dB) */
static void test_response(void) {
  double worst = 0;
  for (uint32_t order = 1; order <= UDSP_BANK_ORDER; order++) {
    const float hz[1] = {1000}, q[1] = {4}, g[1] = {0.7f};
    memset(&bank, 0, sizeof bank);
    udsp_bank_set(&bank, 1, hz, q, g, order, RATE);
    const double f[3] = {1000, 1132, 2000};
    for (int k = 0; k < 3; k++) {
      float re, im;
      udsp_bank_response(&bank, 0, (float)f[k], RATE, &re, &im);
      const double m = hypot(re, im), want = gain_at(f[k]);
      worst = fmax(worst, fabs(20 * log10(m / want)));
    }
  }
  float re, im;
  udsp_bank_response(&bank, 0, 1000, RATE, &re, &im);
  UDSP_TEST_ASSERT(worst < 0.1 && fabsf(im) < 1e-4f,
                   "the response within %.3f dB of the band's, real at the "
                   "centre",
                   worst);
}

int main(void) {
  printf("bank_test\n");
  test_bands();
  test_gain_retune();
  test_follow();
  test_noise_width();
  test_response();
  return udsp_test_summary();
}
