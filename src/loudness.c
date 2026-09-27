#include "loudness.h"

#include <math.h>

/* ITU-R BS.1770's K-weighting cascade, as two second-order sections
 * (Direct Form I): stage 1 is a high-shelf approximating the head's
 * acoustic response, stage 2 is a high-pass approximating the RLB
 * (revised low-frequency B) weighting curve. Coefficients are the spec's
 * published values for a 48kHz design -- see loudness.h's doc comment for
 * why these are applied unmodified regardless of the caller's actual
 * sample rate. */
typedef struct {
  double b0, b1, b2, a1, a2;
} UDSP_BIQUAD;

static const UDSP_BIQUAD UDSP_K_STAGE1 = {
    .b0 = 1.53512485958697,
    .b1 = -2.69169618940638,
    .b2 = 1.19839281085285,
    .a1 = -1.69065929318241,
    .a2 = 0.73248077421585,
};

static const UDSP_BIQUAD UDSP_K_STAGE2 = {
    .b0 = 1.0,
    .b1 = -2.0,
    .b2 = 1.0,
    .a1 = -1.99004745483398,
    .a2 = 0.99007225036621,
};

/* both stages and the mean square in one pass, each stage in transposed
 * Direct Form II. One pass lets stage 2 of a sample overlap stage 1 of
 * the next (as two passes each stage waited on itself); DF-II-T puts one
 * multiply-add, not a chain of them, between a stage's output and its
 * next state. Doubles throughout: stage 2's poles sit at 0.995, where
 * float state drifts. */
static double k_weighted_sumsq(const float *mono, size_t n) {
  const UDSP_BIQUAD *f = &UDSP_K_STAGE1, *g = &UDSP_K_STAGE2;
  double s1 = 0.0, s2 = 0.0, t1 = 0.0, t2 = 0.0, sumsq = 0.0;
  for (size_t i = 0; i < n; i++) {
    double x = mono[i];
    double y = f->b0 * x + s1;
    s1 = f->b1 * x - f->a1 * y + s2;
    s2 = f->b2 * x - f->a2 * y;
    double z = g->b0 * y + t1;
    t1 = g->b1 * y - g->a1 * z + t2;
    t2 = g->b2 * y - g->a2 * z;
    sumsq += z * z;
  }
  return sumsq;
}

float udsp_loudness_lufs_momentary(const float *mono, size_t n, float *scratch) {
  (void)scratch; /* one pass needs none; kept for the callers */
  if (n == 0) {
    return UDSP_LUFS_FLOOR;
  }
  double mean_square = k_weighted_sumsq(mono, n) / (double)n;
  if (mean_square <= 0.0) {
    return UDSP_LUFS_FLOOR;
  }

  float lufs = (float)(-0.691 + 10.0 * log10(mean_square));
  return lufs < UDSP_LUFS_FLOOR ? UDSP_LUFS_FLOOR : lufs;
}
