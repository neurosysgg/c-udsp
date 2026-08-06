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

/* in and out may be the same buffer -- safe because each output sample is
 * written only after the inputs it depends on (x1/x2, held in locals, not
 * re-read from the array) have already been captured. */
static void biquad_apply(const UDSP_BIQUAD *f, const float *in, float *out, size_t n) {
  double x1 = 0.0, x2 = 0.0, y1 = 0.0, y2 = 0.0;
  for (size_t i = 0; i < n; i++) {
    double x0 = in[i];
    double y0 = f->b0 * x0 + f->b1 * x1 + f->b2 * x2 - f->a1 * y1 - f->a2 * y2;
    out[i] = (float)y0;
    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
  }
}

float udsp_loudness_lufs_momentary(const float *mono, size_t n, float *scratch) {
  if (n == 0) {
    return UDSP_LUFS_FLOOR;
  }

  biquad_apply(&UDSP_K_STAGE1, mono, scratch, n);
  biquad_apply(&UDSP_K_STAGE2, scratch, scratch, n);

  double sumsq = 0.0;
  for (size_t i = 0; i < n; i++) {
    double v = scratch[i];
    sumsq += v * v;
  }
  double mean_square = sumsq / (double)n;
  if (mean_square <= 0.0) {
    return UDSP_LUFS_FLOOR;
  }

  float lufs = (float)(-0.691 + 10.0 * log10(mean_square));
  return lufs < UDSP_LUFS_FLOOR ? UDSP_LUFS_FLOOR : lufs;
}
