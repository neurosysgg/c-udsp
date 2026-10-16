#include "dyn.h"

#include <string.h>

#include "pitch.h"

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
  static const float INV[9] = {0,         1.0f,      1.0f / 2,
                               1.0f / 3,  1.0f / 4,  1.0f / 5,
                               1.0f / 6,  1.0f / 7,  1.0f / 8};
  float p = 1.0f;
  for (int k = 8; k >= 1; k--) {
    p = 1.0f + p * y * INV[k];
  }
  BITS s = {.u = (uint32_t)(n + 127) << 23};
  return p * s.f;
}

/* a one-pole's coefficient for a time constant of ms: exp(-1 / frames) */
static float pole(float ms, unsigned rate) {
  const double frames = (double)ms * rate / 1000.0;
  return frames <= 0.0 ? 0.0f
                       : (float)udsp_exp2(-1.4426950408889634 / frames);
}

/* the key's level: the louder channel, dB */
static inline float peak_db(const float *k) {
  float a = k[0] < 0.0f ? -k[0] : k[0];
  float b = k[1] < 0.0f ? -k[1] : k[1];
  return udsp_db_from_gain(a > b ? a : b);
}

/* --- the compressor ---------------------------------------------------- */

void udsp_comp_reset(UDSP_COMP *c) {
  c->red = 0.0f;
  c->ms = 0.0f;
}

void udsp_comp_set(UDSP_COMP *c, float threshold_db, float ratio,
                   float knee_db, float attack_ms, float release_ms, int rms,
                   unsigned rate) {
  c->thr = threshold_db;
  c->slope = (ratio >= 1e9f ? 0.0f : 1.0f / (ratio < 1.0f ? 1.0f : ratio)) -
             1.0f;
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

void udsp_comp_gain(UDSP_COMP *c, const float *key, float *gain,
                    size_t frames) {
  float red = c->red, ms = c->ms;
  const float att = c->att, rel = c->rel, mc = c->ms_coef;
  for (size_t i = 0; i < frames; i++) {
    float x;
    if (c->rms) {
      const float l = key[2 * i], r = key[2 * i + 1];
      const float sq = l * l > r * r ? l * l : r * r;
      ms = sq + mc * (ms - sq);
      x = 0.5f * UDSP_DB_PER_LOG2 * udsp_log2f(ms < 1e-30f ? 1e-30f : ms);
    } else {
      x = peak_db(key + 2 * i);
    }
    const float want = reduction(c, x);
    const float a = want > red ? att : rel;
    red = want + a * (red - want);
    gain[i] = udsp_gain_from_db(-red);
  }
  c->red = red;
  c->ms = ms;
}

float udsp_comp_curve(const UDSP_COMP *c, float in_db) {
  return in_db - reduction(c, in_db);
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
  float db = g->db;
  for (size_t i = 0; i < frames; i++) {
    const float x = peak_db(key + 2 * i);
    if (x >= g->open_db) {
      g->open = 1;
      g->held = 0;
    } else if (g->open && x < g->close_db && ++g->held > g->hold) {
      g->open = 0;
    }
    const float want = g->open ? 0.0f : g->floor_db;
    const float a = want > db ? g->att : g->rel;
    db = want + a * (db - want);
    /* at the floor's -100 dB the gate is shut: silence, not -100 dB */
    gain[i] = db <= -99.9f ? 0.0f : udsp_gain_from_db(db);
  }
  g->db = db;
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
