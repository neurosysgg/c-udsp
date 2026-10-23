#include "fx/dyn.h"

#include <math.h>
#include <string.h>

#include "base/simd.h"
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
  for (int k = 8; k >= 1; k--) {
    p = 1.0f + p * y * INV[k];
  }
  return p * (UDSP_VF)((n + 127) << 23);
}
#endif

/* x[i] = scale * log2(x[i]), x > 0 */
static void log2_scaled(float *x, size_t n, float scale) {
  size_t i = 0;
#if UDSP_W > 1
  for (; i + UDSP_W <= n; i += UDSP_W) {
    udsp_vstore(x + i, scale * vlog2f(udsp_vload(x + i)));
  }
#endif
  for (; i < n; i++) {
    x[i] = scale * udsp_log2f(x[i]);
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

/* the key's level a frame into x, dB: its peak (the louder channel's)
 * or its mean square over ms_coef's 10 ms (ms the running one) */
static void key_db(const float *key, float *x, size_t frames, int rms, float mc,
                   float *ms_io) {
  if (rms) {
    float ms = *ms_io;
    for (size_t i = 0; i < frames; i++) {
      const float l = key[2 * i], r = key[2 * i + 1];
      const float sq = l * l > r * r ? l * l : r * r;
      ms = sq + mc * (ms - sq);
      x[i] = ms < 1e-30f ? 1e-30f : ms;
    }
    *ms_io = ms;
    log2_scaled(x, frames, 0.5f * UDSP_DB_PER_LOG2);
  } else {
    for (size_t i = 0; i < frames; i++) {
      float a = key[2 * i] < 0.0f ? -key[2 * i] : key[2 * i];
      float b = key[2 * i + 1] < 0.0f ? -key[2 * i + 1] : key[2 * i + 1];
      a = a > b ? a : b;
      x[i] = a < 1e-15f ? 1e-15f : a;
    }
    log2_scaled(x, frames, UDSP_DB_PER_LOG2);
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
}

/* the reduction (dB, >= 0) the curve asks for at level x */
static inline float reduction(const UDSP_COMP *c, float x) {
  const float d = x - c->thr;
  const float w = c->knee;
  if (2.0f * d <= -w) {
    return 0.0f;
  }
  if (2.0f * d < w) {
    const float e = d + 0.5f * w;
    return -c->slope * e * e / (2.0f * w);
  }
  return -c->slope * d;
}

/* in passes over gain: the level, the curve's reduction (lanes), the
 * smoother (the one serial chain), the gain (lanes) -- per frame the
 * operations a frame at a time would do */
void udsp_comp_gain(UDSP_COMP *c, const float *key, float *gain,
                    size_t frames) {
  key_db(key, gain, frames, c->rms, c->ms_coef, &c->ms);
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
  gains_from_db(gain, frames, -INFINITY); /* finite: never */
}

float udsp_comp_curve(const UDSP_COMP *c, float in_db) {
  return in_db - reduction(c, in_db);
}

/* --- the up / down compressor ------------------------------------------ */

void udsp_udcomp_reset(UDSP_UDCOMP *c) {
  c->db = 0.0f;
  c->ms = 0.0f;
}

static float slope_of(float ratio) {
  return ratio >= 1e9f ? 0.0f : 1.0f / (ratio < 1.0f ? 1.0f : ratio);
}

void udsp_udcomp_set(UDSP_UDCOMP *c, float above_db, float above_ratio,
                     float below_db, float below_ratio, float up_max_db,
                     float floor_db, float attack_ms, float release_ms, int rms,
                     unsigned rate) {
  c->a_thr = above_db;
  c->b_thr = below_db < above_db ? below_db : above_db;
  c->a_slope = slope_of(above_ratio) - 1.0f;
  c->b_slope = 1.0f - slope_of(below_ratio);
  c->up_max = up_max_db < 0.0f ? 0.0f : up_max_db;
  c->floor = floor_db;
  c->att = pole(attack_ms, rate);
  c->rel = pole(release_ms, rate);
  c->ms_coef = pole(10.0f, rate);
  c->rms = rms;
}

/* the gain (dB) the curve asks for at level x */
static inline float ud_gain(const UDSP_UDCOMP *c, float x) {
  if (x > c->a_thr) {
    return c->a_slope * (x - c->a_thr);
  }
  if (x < c->b_thr && x > c->floor) {
    const float up = c->b_slope * (c->b_thr - x);
    return up < c->up_max ? up : c->up_max;
  }
  return 0.0f;
}

/* comp's passes */
void udsp_udcomp_gain(UDSP_UDCOMP *c, const float *key, float *gain,
                      size_t frames) {
  key_db(key, gain, frames, c->rms, c->ms_coef, &c->ms);
  size_t i = 0;
#if UDSP_W > 1
  for (; i + UDSP_W <= frames; i += UDSP_W) {
    const UDSP_VF x = udsp_vload(gain + i);
    const UDSP_VF down = c->a_slope * (x - c->a_thr);
    UDSP_VF up = c->b_slope * (c->b_thr - x);
    up = udsp_vsel(up < c->up_max, up, udsp_vsplat(c->up_max));
    const UDSP_VF r =
        udsp_vsel((x < c->b_thr) & (x > c->floor), up, udsp_vsplat(0.0f));
    udsp_vstore(gain + i, udsp_vsel(x > c->a_thr, down, r));
  }
#endif
  for (; i < frames; i++) {
    gain[i] = ud_gain(c, gain[i]);
  }
  float db = c->db;
  const float att = c->att, rel = c->rel;
  for (i = 0; i < frames; i++) {
    const float want = gain[i];
    const float a = want < db ? att : rel;
    db = want + a * (db - want);
    gain[i] = db;
  }
  c->db = db;
  gains_from_db(gain, frames, -INFINITY); /* finite: never */
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
