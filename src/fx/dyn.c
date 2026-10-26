#include "fx/dyn.h"

#include <math.h>
#include <string.h>

#include "base/simd.h"
#include "fx/shape.h"
#include "synth/pitch.h"

typedef union {
  float f;
  uint32_t u;
} BITS;

float udsp_log2f(float x) {
  BITS v = {x};
  int e = (int)((v.u >> 23) & 255u) - 127;
  /* the mantissa in [sqrt(1/2), sqrt(2)): no cancellation just under 1
   * (-1 + 0.99...), and |t| <= 0.172 */
  if ((v.u & 0x007fffffu) > 0x003504f3u) {
    v.u = (v.u & 0x007fffffu) | 0x3f000000u;
    e++;
  } else {
    v.u = (v.u & 0x007fffffu) | 0x3f800000u;
  }
  const float t = (v.f - 1.0f) / (v.f + 1.0f);
  const float t2 = t * t;
  /* ln(m) = 2 atanh(t), over ln 2 */
  const float s =
      1.0f + t2 * (1.0f / 3 + t2 * (1.0f / 5 + t2 * (1.0f / 7 + t2 / 9)));
  return (float)e + 2.8853900817779268f * t * s;
}

float udsp_exp2f(float x) {
  x = x < -126.0f ? -126.0f : (x > 127.0f ? 127.0f : x);
  int n = (int)x;
  n -= x < (float)n; /* floor */
  const float y = (x - (float)n) * 0.69314718f;
  /* e^y's Taylor series to y^8 / 8!, Horner's way, no division */
  static const float INV[9] = {0,        1.0f,     1.0f / 2, 1.0f / 3, 1.0f / 4,
                               1.0f / 5, 1.0f / 6, 1.0f / 7, 1.0f / 8};
  float p = 1.0f;
  UDSP_UNROLL
  for (int k = 8; k >= 1; k--) {
    p = 1.0f + p * y * INV[k];
  }
  BITS s = {.u = (uint32_t)(n + 127) << 23};
  return p * s.f;
}

#if UDSP_W > 1
/* udsp_log2f() and udsp_exp2f() per lane: their operations in their
 * order, so a lane's bits are the scalar call's */
static UDSP_VF vlog2f(UDSP_VF x) {
  UDSP_VU u = (UDSP_VU)x;
  UDSP_VI e = (UDSP_VI)((u >> 23) & 255u) - 127;
  const UDSP_VU m = u & 0x007fffffu;
  const UDSP_VI big = (UDSP_VI)(m > 0x003504f3u);
  u = m | (((UDSP_VU)big & 0x3f000000u) | (~(UDSP_VU)big & 0x3f800000u));
  e -= big;
  const UDSP_VF f = (UDSP_VF)u;
  const UDSP_VF t = (f - 1.0f) / (f + 1.0f);
  const UDSP_VF t2 = t * t;
  const UDSP_VF s =
      1.0f + t2 * (1.0f / 3 + t2 * (1.0f / 5 + t2 * (1.0f / 7 + t2 / 9)));
  return __builtin_convertvector(e, UDSP_VF) + 2.8853900817779268f * t * s;
}

static UDSP_VF vexp2f(UDSP_VF x) {
  x = udsp_vsel(x < -126.0f, udsp_vsplat(-126.0f),
                udsp_vsel(x > 127.0f, udsp_vsplat(127.0f), x));
  UDSP_VI n = __builtin_convertvector(x, UDSP_VI);
  n += x < __builtin_convertvector(n, UDSP_VF); /* floor: -1 where true */
  const UDSP_VF y = (x - __builtin_convertvector(n, UDSP_VF)) * 0.69314718f;
  static const float INV[9] = {0,        1.0f,     1.0f / 2, 1.0f / 3, 1.0f / 4,
                               1.0f / 5, 1.0f / 6, 1.0f / 7, 1.0f / 8};
  UDSP_VF p = udsp_vsplat(1.0f);
  UDSP_UNROLL
  for (int k = 8; k >= 1; k--) {
    p = 1.0f + p * y * INV[k];
  }
  return p * (UDSP_VF)((n + 127) << 23);
}
#endif

/* x[i] = scale * log2(x[i]), x floored at lo (> 0) */
static void log2_scaled(float *x, size_t n, float lo, float scale) {
  size_t i = 0;
#if UDSP_W > 1
  for (; i + UDSP_W <= n; i += UDSP_W) {
    const UDSP_VF v = udsp_vload(x + i);
    udsp_vstore(x + i, scale * vlog2f(udsp_vsel(v < lo, udsp_vsplat(lo), v)));
  }
#endif
  for (; i < n; i++) {
    x[i] = scale * udsp_log2f(x[i] < lo ? lo : x[i]);
  }
}

/* db[i] in dB to a gain; at or under `shut` dB, 0 (the gate's floor) */
static void gains_from_db(float *db, size_t n, float shut) {
  size_t i = 0;
#if UDSP_W > 1
  for (; i + UDSP_W <= n; i += UDSP_W) {
    const UDSP_VF d = udsp_vload(db + i);
    const UDSP_VF g = vexp2f(d * (1.0f / UDSP_DB_PER_LOG2));
    udsp_vstore(db + i, udsp_vsel(d <= shut, udsp_vsplat(0.0f), g));
  }
#endif
  for (; i < n; i++) {
    db[i] = db[i] <= shut ? 0.0f : udsp_gain_from_db(db[i]);
  }
}

/* the key's louder channel a frame into x: its |x|, or its square */
static void key_peak(const float *key, float *x, size_t frames) {
  size_t i = 0;
#if UDSP_W > 1
  for (; i + UDSP_W <= frames; i += UDSP_W) {
    UDSP_VF l, r;
    udsp_vlr(key + 2 * i, &l, &r);
    /* (-0's sign and a NaN's are lost: the floor and log2 never see
     * them) */
    l = udsp_vabs(l);
    r = udsp_vabs(r);
    udsp_vstore(x + i, udsp_vsel(l > r, l, r));
  }
#endif
  for (; i < frames; i++) {
    const float a = key[2 * i] < 0.0f ? -key[2 * i] : key[2 * i];
    const float b = key[2 * i + 1] < 0.0f ? -key[2 * i + 1] : key[2 * i + 1];
    x[i] = a > b ? a : b;
  }
}

static void key_sq(const float *key, float *x, size_t frames) {
  size_t i = 0;
#if UDSP_W > 1
  for (; i + UDSP_W <= frames; i += UDSP_W) {
    UDSP_VF l, r;
    udsp_vlr(key + 2 * i, &l, &r);
    l *= l;
    r *= r;
    udsp_vstore(x + i, udsp_vsel(l > r, l, r));
  }
#endif
  for (; i < frames; i++) {
    const float l = key[2 * i], r = key[2 * i + 1];
    x[i] = l * l > r * r ? l * l : r * r;
  }
}

/* the mean square of the two channels a frame into x */
static void key_mean_sq(const float *key, float *x, size_t frames) {
  size_t i = 0;
#if UDSP_W > 1
  for (; i + UDSP_W <= frames; i += UDSP_W) {
    UDSP_VF l, r;
    udsp_vlr(key + 2 * i, &l, &r);
    udsp_vstore(x + i, 0.5f * (l * l + r * r));
  }
#endif
  for (; i < frames; i++) {
    const float l = key[2 * i], r = key[2 * i + 1];
    x[i] = 0.5f * (l * l + r * r);
  }
}

#define PEAK_FLOOR 1e-15f /* -300 dB */
#define MS_FLOOR 1e-30f

/* s = x + a (s - x) along each of n rows of x, in place, a = att[r]
 * while x falls under s, else rel[r]: the smoother (RMS's mean square
 * too, att = rel). Frame by frame across the rows, so their chains (~12
 * cycles a frame each, the CPU waiting on them) overlap; each row's
 * operations are its own loop's, in order */
static inline __attribute__((always_inline)) void
follow_n(float *const x[], float *s, const float *att, const float *rel,
         size_t frames, size_t n) {
  /* all local: no store into a row can change them (-fno-strict-aliasing
   * builds would reload the pointers every frame) */
  float *row[UDSP_UDCOMP_N];
  float v[UDSP_UDCOMP_N], at[UDSP_UDCOMP_N], re[UDSP_UDCOMP_N];
  for (size_t r = 0; r < n; r++) {
    row[r] = x[r];
    v[r] = s[r];
    at[r] = att[r];
    re[r] = rel[r];
  }
  for (size_t i = 0; i < frames; i++) {
    UDSP_UNROLL
    for (size_t r = 0; r < n; r++) {
      const float want = row[r][i];
      const float a = want < v[r] ? at[r] : re[r];
      v[r] = want + a * (v[r] - want);
      row[r][i] = v[r];
    }
  }
  for (size_t r = 0; r < n; r++) {
    s[r] = v[r];
  }
}

static void follow(float *const x[], float *s, const float *att,
                   const float *rel, size_t n, size_t frames) {
  switch (n) {
  case 1:
    follow_n(x, s, att, rel, frames, 1);
    break;
  case 2:
    follow_n(x, s, att, rel, frames, 2);
    break;
  case 3:
    follow_n(x, s, att, rel, frames, 3);
    break;
  case 4:
    follow_n(x, s, att, rel, frames, 4);
    break;
  default:
    break;
  }
}

/* the key's level a frame into x, dB: its peak (the louder channel's)
 * or its mean square over ms_coef's 10 ms (ms the running one) */
static void key_db(const float *key, float *x, size_t frames, int rms, float mc,
                   float *ms_io) {
  if (rms) {
    key_sq(key, x, frames);
    follow(&x, ms_io, &mc, &mc, 1, frames);
    log2_scaled(x, frames, MS_FLOOR, 0.5f * UDSP_DB_PER_LOG2);
  } else {
    key_peak(key, x, frames);
    log2_scaled(x, frames, PEAK_FLOOR, UDSP_DB_PER_LOG2);
  }
}

/* a one-pole's coefficient for a time constant of ms: exp(-1 / frames) */
static float pole(float ms, unsigned rate) {
  const double frames = (double)ms * rate / 1000.0;
  return frames <= 0.0 ? 0.0f : (float)udsp_exp2(-1.4426950408889634 / frames);
}

/* --- the compressor ---------------------------------------------------- */

void udsp_comp_reset(UDSP_COMP *c) {
  c->red = 0.0f;
  c->ms = 0.0f;
  c->red2 = 0.0f;
  c->held = 0.0f;
  c->r2 = c->t2 = 0.0f;
  c->n = 0;
  c->cur = c->prev = 0.0f;
}

void udsp_comp_set(UDSP_COMP *c, float threshold_db, float ratio, float knee_db,
                   float attack_ms, float release_ms, int rms, unsigned rate) {
  c->thr = threshold_db;
  c->slope =
      (ratio >= 1e9f ? 0.0f : 1.0f / (ratio < 1.0f ? 1.0f : ratio)) - 1.0f;
  c->knee = knee_db < 0.0f ? 0.0f : knee_db;
  c->att = pole(attack_ms, rate);
  c->rel = pole(release_ms, rate);
  c->ms_coef = pole(10.0f, rate);
  c->rms = rms;
  c->second = 0; /* set_second, set_smooth, set_hold, set_window after */
  c->post = 0.0f;
  c->hold = 0.0f;
  c->block = 0;
}

void udsp_comp_set_second(UDSP_COMP *c, float threshold_db, float ratio2,
                          float knee_db, float attack_ms, float release_ms,
                          float tail_ms, unsigned rate) {
  c->second = ratio2 > 0.0f;
  c->thr2 = threshold_db;
  c->slope2 = c->second
                  ? (ratio2 >= 1e9f ? 0.0f : 1.0f / ratio2) - 1.0f - c->slope
                  : 0.0f;
  c->knee2 = knee_db < 0.0f ? 0.0f : knee_db;
  c->att2 = pole(attack_ms, rate);
  c->rel2 = pole(release_ms, rate);
  c->tail = tail_ms > 0.0f ? pole(tail_ms, rate) : 0.0f;
}

void udsp_comp_set_smooth(UDSP_COMP *c, float ms, unsigned rate) {
  c->post = ms > 0.0f ? pole(ms, rate) : 0.0f;
}

void udsp_comp_set_hold(UDSP_COMP *c, float ms, unsigned rate) {
  c->hold = ms > 0.0f ? pole(ms, rate) : 0.0f;
}

void udsp_comp_set_window(UDSP_COMP *c, float ms, unsigned rate) {
  const float b = ms * (float)rate * 0.001f;
  c->block = ms > 0.0f ? (b < 1.0f ? 1u : (uint32_t)(b + 0.5f)) : 0u;
  c->n = 0;
}

/* a knee's reduction (dB) for a slope's change past thr, w wide */
static inline float knee_red(float x, float thr, float slope, float w) {
  const float d = x - thr;
  if (2.0f * d <= -w) {
    return 0.0f;
  }
  if (2.0f * d < w) {
    const float e = d + 0.5f * w;
    return -slope * e * e / (2.0f * w);
  }
  return -slope * d;
}

/* the reduction (dB, >= 0) the curve asks for at level x */
static inline float reduction(const UDSP_COMP *c, float x) {
  return knee_red(x, c->thr, c->slope, c->knee);
}

/* in passes over gain: the level, the curve's reduction (lanes), the
 * smoother (the one serial chain), the gain (lanes) -- per frame the
 * operations a frame at a time would do */
void udsp_comp_gain(UDSP_COMP *c, const float *key, float *gain,
                    size_t frames) {
  if ((c->hold > 0.0f || c->block) && !c->rms) {
    key_peak(key, gain, frames);
    if (c->block) {
      float cur = c->cur, prev = c->prev;
      uint32_t n = c->n;
      for (size_t k = 0; k < frames; k++) {
        cur = gain[k] > cur ? gain[k] : cur;
        gain[k] = cur > prev ? cur : prev;
        if (++n >= c->block) {
          prev = cur;
          cur = 0.0f;
          n = 0;
        }
      }
      c->cur = cur;
      c->prev = prev;
      c->n = n;
    }
    float h = c->held;
    for (size_t k = 0; k < frames; k++) {
      h *= c->hold;
      h = gain[k] > h ? gain[k] : h;
      gain[k] = h;
    }
    c->held = h;
    log2_scaled(gain, frames, PEAK_FLOOR, UDSP_DB_PER_LOG2);
  } else {
    key_db(key, gain, frames, c->rms, c->ms_coef, &c->ms);
  }
  if (c->second) { /* the two smoothers off one level, a frame at a time */
    float red = c->red, red2 = c->red2, r2 = c->r2, t2 = c->t2;
    for (size_t i = 0; i < frames; i++) {
      const float x = gain[i];
      const float w1 = reduction(c, x);
      const float w2 = knee_red(x, c->thr2, c->slope2, c->knee2);
      red = w1 + (w1 > red ? c->att : c->rel) * (red - w1);
      red2 = red + c->post * (red2 - red);
      r2 = w2 + (w2 > r2 ? c->att2 : c->rel2) * (r2 - w2);
      t2 = r2 > t2 ? r2 : r2 + c->tail * (t2 - r2);
      gain[i] = -(red2 + t2);
    }
    c->red = red;
    c->red2 = red2;
    c->r2 = r2;
    c->t2 = t2;
    gains_from_db(gain, frames, -INFINITY); /* finite: never */
    return;
  }
  size_t i = 0;
#if UDSP_W > 1
  const float neg = -c->slope, w = c->knee;
  for (; i + UDSP_W <= frames; i += UDSP_W) {
    const UDSP_VF d = udsp_vload(gain + i) - c->thr;
    const UDSP_VF e = d + 0.5f * w;
    const UDSP_VF knee = neg * e * e / (2.0f * w), line = neg * d;
    const UDSP_VF r = udsp_vsel(2.0f * d < w, knee, line);
    udsp_vstore(gain + i, udsp_vsel(2.0f * d <= -w, udsp_vsplat(0.0f), r));
  }
#endif
  for (; i < frames; i++) {
    gain[i] = reduction(c, gain[i]);
  }
  float red = c->red;
  const float att = c->att, rel = c->rel;
  for (i = 0; i < frames; i++) {
    const float want = gain[i];
    const float a = want > red ? att : rel;
    red = want + a * (red - want);
    gain[i] = -red;
  }
  c->red = red;
  if (c->post > 0.0f) {
    float r2 = c->red2;
    for (i = 0; i < frames; i++) {
      r2 = -gain[i] + c->post * (r2 + gain[i]);
      gain[i] = -r2;
    }
    c->red2 = r2;
  }
  gains_from_db(gain, frames, -INFINITY); /* finite: never */
}

float udsp_comp_curve(const UDSP_COMP *c, float in_db) {
  const float r2 =
      c->second ? knee_red(in_db, c->thr2, c->slope2, c->knee2) : 0.0f;
  return in_db - reduction(c, in_db) - r2;
}

/* --- the up / down compressor ------------------------------------------ */

void udsp_udcomp_reset(UDSP_UDCOMP *c) {
  c->db = 0.0f;
  c->ms = 0.0f;
  c->lvl = -300.0f;
}

static float slope_of(float ratio) {
  return ratio >= 1e9f ? 0.0f : 1.0f / (ratio < 1.0f ? 1.0f : ratio);
}

void udsp_udcomp_set_slopes(UDSP_UDCOMP *c, float above_db, float above_slope,
                            float below_db, float below_slope, float knee_db,
                            float up_max_db, float floor_db, float attack_ms,
                            float release_ms, int rms, unsigned rate) {
  c->a_thr = above_db;
  c->b_thr = below_db < above_db ? below_db : above_db;
  c->a_slope = (above_slope < 0.0f ? 0.0f : above_slope) - 1.0f;
  c->b_slope = 1.0f - (below_slope < 0.0f ? 0.0f : below_slope);
  c->knee = knee_db > 0.0f ? knee_db : 0.0f;
  c->knee_inv = c->knee > 0.0f ? 0.5f / c->knee : 0.0f;
  c->up_max = up_max_db < 0.0f ? 0.0f : up_max_db;
  c->floor = floor_db;
  c->att = pole(attack_ms, rate);
  c->rel = pole(release_ms, rate);
  c->ms_coef = pole(10.0f, rate);
  c->rms = rms;
}

void udsp_udcomp_set(UDSP_UDCOMP *c, float above_db, float above_ratio,
                     float below_db, float below_ratio, float up_max_db,
                     float floor_db, float attack_ms, float release_ms, int rms,
                     unsigned rate) {
  udsp_udcomp_set_slopes(c, above_db, slope_of(above_ratio), below_db,
                         slope_of(below_ratio), 0.0f, up_max_db, floor_db,
                         attack_ms, release_ms, rms, rate);
}

/* the gain (dB) the curve asks for at level x: above's and below's, each
 * a line past its threshold and a parabola across its knee (a knee of 0
 * never takes the parabola: 2 d can't be both under 0 and over -0) */
static inline float ud_gain(const UDSP_UDCOMP *c, float x) {
  const float w = c->knee;
  const float da = x - c->a_thr, db = x - c->b_thr;
  float hi = 0.0f, lo = 0.0f;
  if (2.0f * da >= w) {
    hi = c->a_slope * da;
  } else if (2.0f * da > -w) {
    const float e = da + 0.5f * w;
    hi = c->a_slope * e * e * c->knee_inv;
  }
  if (2.0f * db <= -w) {
    lo = c->b_slope * (c->b_thr - x);
  } else if (2.0f * db < w) {
    const float e = db - 0.5f * w;
    lo = c->b_slope * e * e * c->knee_inv;
  }
  if (lo > 0.0f && x <= c->floor) {
    lo = 0.0f;
  }
  lo = lo < -UDSP_UDCOMP_DEPTH ? -UDSP_UDCOMP_DEPTH : lo;
  const float g = hi + lo;
  return g < c->up_max ? g : c->up_max;
}

/* the curve's gain (dB) for each level in x, in place */
static void ud_want(const UDSP_UDCOMP *c, float *x, size_t frames) {
  size_t i = 0;
#if UDSP_W > 1
  const float w = c->knee, half = 0.5f * w, inv = c->knee_inv;
  const UDSP_VF zero = udsp_vsplat(0.0f);
  /* no knee: the lines alone (the same values: 2 d >= 0 is d >= 0) */
  for (; w == 0.0f && i + UDSP_W <= frames; i += UDSP_W) {
    const UDSP_VF v = udsp_vload(x + i);
    const UDSP_VF hi =
        udsp_vsel(v >= c->a_thr, c->a_slope * (v - c->a_thr), zero);
    UDSP_VF lo = c->b_slope * (c->b_thr - v);
    lo = udsp_vsel((v <= c->b_thr) & ((lo <= 0.0f) | (v > c->floor)), lo, zero);
    lo =
        udsp_vsel(lo < -UDSP_UDCOMP_DEPTH, udsp_vsplat(-UDSP_UDCOMP_DEPTH), lo);
    const UDSP_VF g = hi + lo;
    udsp_vstore(x + i, udsp_vsel(g < c->up_max, g, udsp_vsplat(c->up_max)));
  }
  for (; i + UDSP_W <= frames; i += UDSP_W) {
    const UDSP_VF v = udsp_vload(x + i);
    const UDSP_VF da = v - c->a_thr, db = v - c->b_thr;
    const UDSP_VF ea = da + half, eb = db - half;
    UDSP_VF hi = udsp_vsel(2.0f * da > -w, c->a_slope * ea * ea * inv, zero);
    hi = udsp_vsel(2.0f * da >= w, c->a_slope * da, hi);
    UDSP_VF lo = udsp_vsel(2.0f * db < w, c->b_slope * eb * eb * inv, zero);
    lo = udsp_vsel(2.0f * db <= -w, c->b_slope * (c->b_thr - v), lo);
    lo = udsp_vsel((lo > 0.0f) & (v <= c->floor), zero, lo);
    lo =
        udsp_vsel(lo < -UDSP_UDCOMP_DEPTH, udsp_vsplat(-UDSP_UDCOMP_DEPTH), lo);
    const UDSP_VF g = hi + lo;
    udsp_vstore(x + i, udsp_vsel(g < c->up_max, g, udsp_vsplat(c->up_max)));
  }
#endif
  for (; i < frames; i++) {
    x[i] = ud_gain(c, x[i]);
  }
}

/* comp's passes, the recurrences side by side: RMS's (the louder
 * channel's, or the power's mean) for those detecting it (key_db()'s),
 * then the smoother's */
void udsp_udcomp_gains(UDSP_UDCOMP *const c[], const float *const key[],
                       float *const gain[], size_t n, size_t frames) {
  n = n > UDSP_UDCOMP_N ? UDSP_UDCOMP_N : n;
  /* (zeroed: gcc can't tell only the first m are read) */
  float *rows[UDSP_UDCOMP_N] = {0};
  float ms[UDSP_UDCOMP_N] = {0}, s[UDSP_UDCOMP_N] = {0};
  float att[UDSP_UDCOMP_N] = {0}, rel[UDSP_UDCOMP_N] = {0};
  size_t m = 0;
  for (size_t j = 0; j < n; j++) {
    if (!c[j]->rms) {
      key_db(key[j], gain[j], frames, 0, 0.0f, NULL);
      continue;
    }
    float *x = gain[j];
    rows[m] = x;
    ms[m] = c[j]->ms;
    if (c[j]->rms == UDSP_UDCOMP_POWER) {
      key_mean_sq(key[j], x, frames);
    } else {
      key_sq(key[j], x, frames);
    }
    att[m] = rel[m] = c[j]->ms_coef;
    m++;
  }
  follow(rows, ms, att, rel, m, frames);
  m = 0;
  size_t np = 0;
  float *prow[UDSP_UDCOMP_N] = {0};
  float pl[UDSP_UDCOMP_N] = {0}, pa[UDSP_UDCOMP_N] = {0},
        pr[UDSP_UDCOMP_N] = {0};
  for (size_t j = 0; j < n; j++) {
    float *x = gain[j];
    if (c[j]->rms) {
      c[j]->ms = ms[m++];
      log2_scaled(x, frames, MS_FLOOR, 0.5f * UDSP_DB_PER_LOG2);
    }
    if (c[j]->rms == UDSP_UDCOMP_POWER) {
      /* the level smoothed: follow's att while it falls, the release */
      prow[np] = x;
      pl[np] = c[j]->lvl;
      pa[np] = c[j]->rel;
      pr[np] = c[j]->att;
      np++;
    }
  }
  follow(prow, pl, pa, pr, np, frames);
  np = 0;
  for (size_t j = 0; j < n; j++) {
    float *x = gain[j];
    const int power = c[j]->rms == UDSP_UDCOMP_POWER;
    if (power) {
      c[j]->lvl = pl[np++];
      /* the level held at the floor, the curve then as it is */
      UDSP_UDCOMP held = *c[j];
      held.floor = -INFINITY;
      for (size_t i = 0; i < frames; i++) {
        x[i] = x[i] < c[j]->floor ? c[j]->floor : x[i];
      }
      ud_want(&held, x, frames);
    } else {
      ud_want(c[j], x, frames);
    }
    rows[j] = x;
    s[j] = c[j]->db;
    att[j] = power ? 0.0f : c[j]->att; /* a pole of 0: the gain as it is */
    rel[j] = power ? 0.0f : c[j]->rel;
  }
  follow(rows, s, att, rel, n, frames);
  for (size_t j = 0; j < n; j++) {
    c[j]->db = s[j];
    gains_from_db(gain[j], frames, -INFINITY); /* finite: never */
  }
}

void udsp_udcomp_set_window(UDSP_UDCOMP *c, float ms, unsigned rate) {
  c->ms_coef = pole(ms, rate);
}

void udsp_udcomp_gain(UDSP_UDCOMP *c, const float *key, float *gain,
                      size_t frames) {
  udsp_udcomp_gains(&c, &key, &gain, 1, frames);
}

float udsp_udcomp_curve(const UDSP_UDCOMP *c, float in_db) {
  return in_db + ud_gain(c, in_db);
}

/* --- the gate ---------------------------------------------------------- */

void udsp_gate_reset(UDSP_GATE *g) {
  g->open = 0;
  g->held = 0;
  g->db = g->floor_db;
}

void udsp_gate_set(UDSP_GATE *g, float threshold_db, float return_db,
                   float floor_db, float attack_ms, float hold_ms,
                   float release_ms, unsigned rate) {
  g->open_db = threshold_db;
  g->close_db = threshold_db - (return_db < 0.0f ? 0.0f : return_db);
  g->floor_db = floor_db < -100.0f ? -100.0f : floor_db;
  g->att = pole(attack_ms, rate);
  g->rel = pole(release_ms, rate);
  g->hold = (uint32_t)((double)hold_ms * rate / 1000.0);
  if (g->db < g->floor_db) {
    g->db = g->floor_db;
  }
}

void udsp_gate_gain(UDSP_GATE *g, const float *key, float *gain,
                    size_t frames) {
  key_db(key, gain, frames, 0, 0.0f, NULL);
  float db = g->db;
  for (size_t i = 0; i < frames; i++) {
    const float x = gain[i];
    if (x >= g->open_db) {
      g->open = 1;
      g->held = 0;
    } else if (g->open && x < g->close_db && ++g->held > g->hold) {
      g->open = 0;
    }
    const float want = g->open ? 0.0f : g->floor_db;
    const float a = want > db ? g->att : g->rel;
    db = want + a * (db - want);
    gain[i] = db;
  }
  g->db = db;
  /* at the floor's -100 dB the gate is shut: silence, not -100 dB */
  gains_from_db(gain, frames, -99.9f);
}

/* --- the limiter -------------------------------------------------------
 * mem: the delay (2 cap), the window of held gains (cap), the minimum's
 * queue as values and the frames they came in (cap + 2 each: n + 1 at
 * most, and a ring tells full from empty by one slot more) */

static float *lim_delay(UDSP_LIMITER *l) { return l->mem; }
static float *lim_window(UDSP_LIMITER *l) { return l->mem + 2 * l->cap; }
static float *lim_qv(UDSP_LIMITER *l) { return l->mem + 3 * l->cap; }
static uint32_t *lim_qt(UDSP_LIMITER *l) {
  return (uint32_t *)(void *)(l->mem + 4 * l->cap + 2);
}

static void lim_clear(UDSP_LIMITER *l) {
  memset(l->mem, 0, UDSP_LIMITER_FLOATS(l->cap) * sizeof(float));
  float *w = lim_window(l);
  for (uint32_t i = 0; i < l->n; i++) {
    w[i] = 1.0f;
  }
  l->env = 1.0f;
  l->sum = l->n;
  l->t = 0;
  l->head = l->tail = 0;
  l->reduction = 1.0f;
}

void udsp_limiter_init(UDSP_LIMITER *l, float *mem, uint32_t cap) {
  l->mem = mem;
  l->cap = cap ? cap : 1;
  l->n = l->cap;
  l->ceiling = 1.0f;
  l->rel = 0.0f;
  lim_clear(l);
}

void udsp_limiter_set(UDSP_LIMITER *l, uint32_t lookahead, float ceiling,
                      float release_ms, unsigned rate) {
  lookahead = lookahead < 1 ? 1 : (lookahead > l->cap ? l->cap : lookahead);
  l->ceiling = ceiling;
  /* the recovery's step a frame: 1 - exp(-1 / frames) */
  l->rel = 1.0f - pole(release_ms, rate);
  if (lookahead != l->n) {
    l->n = lookahead;
    lim_clear(l);
  }
}

void udsp_limiter_stereo(UDSP_LIMITER *l, float *io, size_t frames) {
  float *d = lim_delay(l), *w = lim_window(l), *qv = lim_qv(l);
  uint32_t *qt = lim_qt(l);
  const uint32_t n = l->n, qcap = l->cap + 2;
  float g = l->reduction;
  for (size_t i = 0; i < frames; i++) {
    const float x0 = io[2 * i], x1 = io[2 * i + 1];
    const float a = x0 < 0.0f ? -x0 : x0, b = x1 < 0.0f ? -x1 : x1;
    const float p = a > b ? a : b;
    const float r = p > l->ceiling ? l->ceiling / p : 1.0f;
    const uint32_t t = (uint32_t)l->t;
    /* the minimum over the n + 1 frames up to this one: the expired
     * dropped first, so the queue never holds more than n + 1 */
    while (l->tail != l->head && t - qt[l->head] > n) {
      l->head = (l->head + 1) % qcap;
    }
    while (l->tail != l->head && qv[(l->tail + qcap - 1) % qcap] >= r) {
      l->tail = (l->tail + qcap - 1) % qcap;
    }
    qv[l->tail] = r;
    qt[l->tail] = t;
    l->tail = (l->tail + 1) % qcap;
    const float m = qv[l->head];
    /* down at once, back up at the release's rate */
    float e = l->env + (1.0f - l->env) * l->rel;
    e = e < m ? e : m;
    l->env = e;
    const uint32_t slot = t % n;
    l->sum += (double)e - (double)w[slot];
    w[slot] = e;
    g = (float)(l->sum / n);
    g = g > 1.0f ? 1.0f : g;
    /* the frame n back out, this one in */
    io[2 * i] = d[2 * slot] * g;
    io[2 * i + 1] = d[2 * slot + 1] * g;
    d[2 * slot] = x0;
    d[2 * slot + 1] = x1;
    l->t++;
  }
  l->reduction = g;
}

/* --- the vintage compressor -------------------------------------------- */

void udsp_vcomp_reset(UDSP_VCOMP *c) {
  c->env = 0.0f;
  c->g = 1.0f;
  c->in_m = c->out_m = 0.0f;
  c->reduction = 0.0f;
}

void udsp_vcomp_set(UDSP_VCOMP *c, float threshold_db, float ratio,
                    float knee_db, float attack_ms, float leak_ms,
                    unsigned rate) {
  c->thr = threshold_db;
  c->slope = ratio < 1.0f ? 0.0f : ratio - 1.0f;
  c->knee = knee_db < 0.0f ? 0.0f : knee_db;
  c->charge = 1.0f - pole(attack_ms, rate);
  c->leak = pole(leak_ms, rate);
}

void udsp_vcomp_set_auto(UDSP_VCOMP *c, float ms, unsigned rate) {
  c->auto_pole = ms > 0.0f ? pole(ms, rate) : 0.0f;
  c->meter = pole(ms > 0.0f ? ms : 500.0f, rate);
}

/* p^n through log2 / exp2 (no libm), 0 once it underflows */
static float pole_pow(float p, size_t n) {
  if (p <= 0.0f) {
    return 0.0f;
  }
  const float l = udsp_log2f(p) * (float)n;
  return l < -120.0f ? 0.0f : udsp_exp2f(l);
}

void udsp_vcomp_idle(UDSP_VCOMP *c, size_t frames) {
  c->env *= pole_pow(c->leak, frames);
  const float m = pole_pow(c->meter, frames);
  c->in_m *= m;
  c->out_m *= m;
  c->reduction = udsp_vcomp_curve(c, udsp_db_from_gain(c->env));
  c->g = udsp_gain_from_db(-c->reduction);
}

float udsp_vcomp_curve(const UDSP_VCOMP *c, float detector_db) {
  return knee_red(detector_db, c->thr, -c->slope, c->knee);
}

void udsp_vcomp_stereo(UDSP_VCOMP *c, UDSP_EQ_BAND *const *key, int keys,
                       float *io, size_t frames) {
  const float a = c->meter;
  for (size_t i = 0; i < frames; i++) {
    float *f = io + 2 * i;
    const float in = fabsf(f[0]) + fabsf(f[1]);
    f[0] *= c->g;
    f[1] *= c->g;
    float k[2] = {f[0], f[1]};
    for (int j = 0; j < keys; j++) {
      udsp_eq_stereo(key[j], k, 1);
    }
    const float v = fabsf(k[0]) > fabsf(k[1]) ? fabsf(k[0]) : fabsf(k[1]);
    float e = c->env;
    if (v > e) {
      e += (v - e) * c->charge;
    }
    e *= c->leak;
    c->env = e < 1e-20f ? 0.0f : e; /* no denormals in a long silence */
    c->reduction = udsp_vcomp_curve(c, udsp_db_from_gain(c->env));
    c->g = udsp_gain_from_db(-c->reduction);
    /* the meters run with the auto gain off too: turned on, it starts
     * from where the levels are */
    c->in_m = a * c->in_m + (1.0f - a) * in;
    c->out_m = a * c->out_m + (1.0f - a) * (fabsf(f[0]) + fabsf(f[1]));
    if (c->auto_pole > 0.0f) {
      const float makeup = (c->in_m + 1e-12f) / (c->out_m + 1e-12f);
      f[0] *= makeup;
      f[1] *= makeup;
    }
  }
}

/* --- the maximizer ----------------------------------------------------- */

void udsp_maxim_reset(UDSP_MAXIM *m) {
  m->t = 0;
  memset(m->need, 0, sizeof(m->need));
  memset(m->held, 0, sizeof(m->held));
  memset(m->delay, 0, sizeof(m->delay));
  m->env[0] = m->env[1] = 0.0f;
  udsp_eq_reset(&m->lp[0]);
  udsp_eq_reset(&m->lp[1]);
  m->reduction = 0.0f;
}

void udsp_maxim_set(UDSP_MAXIM *m, const UDSP_MAXIM_FIT *fit, float gain_db,
                    float ceiling_db, unsigned rate) {
  uint32_t n = fit->lookahead;
  n = n < 1 ? 1 : (n > UDSP_MAXIM_CAP ? UDSP_MAXIM_CAP : n);
  m->gain = udsp_gain_from_db(gain_db);
  m->ceiling = udsp_gain_from_db(ceiling_db);
  m->aim = udsp_gain_from_db(ceiling_db + fit->over_db);
  m->knee = udsp_gain_from_db(-(fit->knee_db > 0.01f ? fit->knee_db : 0.01f));
  m->link = fit->link < 0.0f ? 0.0f : (fit->link > 1.0f ? 1.0f : fit->link);
  m->rel = pole(fit->release_ms, rate);
  m->wlo = fit->weight_lo;
  m->whi = fit->weight_hi;
  for (int k = 0; k < 2; k++) {
    udsp_eq_set(&m->lp[k], UDSP_EQ_HIGH_CUT_6, fit->crossover_hz, 0.7071f, 0.0f,
                rate);
  }
  if (n != m->n) {
    m->n = n;
    udsp_maxim_reset(m);
  }
}

/* the soft clip: x up to k c, past it c (k + (1 - k) tanh((|x| / c - k) /
 * (1 - k))) -- slope 1 at the knee, c in the limit */
static float maxim_clip(float x, float c, float k) {
  const float a = fabsf(x);
  if (a <= c * k) {
    return x;
  }
  const float y = c * (k + (1.0f - k) * udsp_shape1(UDSP_SHAPE_TANH,
                                                    (a / c - k) / (1.0f - k)));
  return x < 0.0f ? -y : y;
}

void udsp_maxim_stereo(UDSP_MAXIM *m, float *io, size_t frames) {
  const uint32_t n = m->n, d = 2 * n - 1;
  const float aim_db = udsp_db_from_gain(m->aim);
  for (size_t i = 0; i < frames; i++) {
    float *f = io + 2 * i;
    const uint32_t slot = (uint32_t)(m->t % n), back = (uint32_t)(m->t % d);
    /* each channel's need to come down to the aim, mixed toward the
     * louder's */
    float need[2];
    for (int c = 0; c < 2; c++) {
      const float a = fabsf(f[c]) * m->gain;
      need[c] = a > m->aim ? aim_db - udsp_db_from_gain(a) : 0.0f;
    }
    const float lo = need[0] < need[1] ? need[0] : need[1];
    float g[2];
    for (int c = 0; c < 2; c++) {
      m->need[slot][c] = (1.0f - m->link) * need[c] + m->link * lo;
      /* the minimum over the window held, let go at the release's pole */
      float v = 0.0f;
      for (uint32_t j = 0; j < n; j++) {
        v = m->need[j][c] < v ? m->need[j][c] : v;
      }
      float e = m->env[c];
      e = v < e ? v : m->rel * e + (1.0f - m->rel) * v;
      m->env[c] = e;
      m->held[slot][c] = e;
      /* averaged over the window: summed afresh (no drift) */
      float s = 0.0f;
      for (uint32_t j = 0; j < n; j++) {
        s += m->held[j][c];
      }
      g[c] = s / (float)n;
    }
    /* the frame 2 n - 1 back out, split and turned down by band */
    float x[2] = {m->delay[back][0] * m->gain, m->delay[back][1] * m->gain};
    m->delay[back][0] = f[0];
    m->delay[back][1] = f[1];
    float low[2] = {x[0], x[1]};
    udsp_eq_stereo(&m->lp[0], low, 1);
    udsp_eq_stereo(&m->lp[1], low, 1);
    for (int c = 0; c < 2; c++) {
      const float y = low[c] * udsp_gain_from_db(m->wlo * g[c]) +
                      (x[c] - low[c]) * udsp_gain_from_db(m->whi * g[c]);
      f[c] = maxim_clip(y, m->ceiling, m->knee);
    }
    m->reduction = -(g[0] < g[1] ? g[0] : g[1]);
    m->t++;
  }
}

void udsp_maxim_idle(UDSP_MAXIM *m, size_t frames) {
  for (int c = 0; c < 2; c++) {
    m->env[c] *= pole_pow(m->rel, frames);
    for (uint32_t j = 0; j < UDSP_MAXIM_CAP; j++) {
      m->need[j][c] = 0.0f;
      m->held[j][c] = m->env[c];
    }
  }
  memset(m->delay, 0, sizeof(m->delay));
  udsp_eq_reset(&m->lp[0]);
  udsp_eq_reset(&m->lp[1]);
  m->reduction = -m->env[0] > -m->env[1] ? -m->env[0] : -m->env[1];
}
