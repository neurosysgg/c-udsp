#include "fx/eq.h"

#include <math.h>

#include "base/simd.h"
#include "synth/osc.h"
#include "synth/pitch.h"

/* an 8th-order Butterworth's sections: 1 / (2 cos((2i + 1) pi / 16)) */
static const float BUTTER8_Q[UDSP_EQ_STAGES] = {0.50979558f, 0.60134489f,
                                                0.89997622f, 2.56291545f};

void udsp_eq_reset(UDSP_EQ_BAND *b) {
  b->ramp = 0;
  for (int s = 0; s < UDSP_EQ_STAGES; s++) {
    b->ic1[s][0] = b->ic1[s][1] = b->ic2[s][0] = b->ic2[s][1] = 0.0f;
  }
}

/* tan(pi t), t a fraction of the rate: a half turn's angle as turns */
static float tan_rate(float t) {
  return udsp_osc_sine(0.5f * t) / udsp_osc_sine(0.5f * t + 0.25f);
}

static void stage(UDSP_EQ_STAGE *s, float g, float k, float m0, float m1,
                  float m2) {
  s->g = g;
  s->k = k;
  s->a1 = 1.0f / (1.0f + g * (g + k));
  s->a2 = g * s->a1;
  s->a3 = g * s->a2;
  s->m0 = m0;
  s->m1 = m1;
  s->m2 = m2;
}

void udsp_eq_set(UDSP_EQ_BAND *b, UDSP_EQ_SHAPE shape, float hz, float q,
                 float gain_db, unsigned rate) {
  float t = hz / (float)rate;
  const float lo = 10.0f / (float)rate;
  t = t < lo ? lo : (t > 0.49f ? 0.49f : t);
  q = q < 0.1f ? 0.1f : (q > 40.0f ? 40.0f : q);
  b->ramp = 0;
  gain_db = gain_db < -30.0f ? -30.0f : (gain_db > 30.0f ? 30.0f : gain_db);
  const float g = tan_rate(t);
  const float k = 1.0f / q;
  /* the amplitude's square root: 10^(dB / 40) */
  const float a = (float)udsp_exp2(gain_db * (3.321928094887362 / 40.0));
  b->stages = 1;
  switch (shape) {
  case UDSP_EQ_LOW_CUT_48:
  case UDSP_EQ_HIGH_CUT_48: {
    const int high = shape == UDSP_EQ_HIGH_CUT_48;
    /* the band's q over a Butterworth's 0.71 sharpens the last section */
    const float res = q * 1.41421356f;
    for (int i = 0; i < UDSP_EQ_STAGES; i++) {
      float sq = BUTTER8_Q[i] * (i == UDSP_EQ_STAGES - 1 ? res : 1.0f);
      float sk = 1.0f / sq;
      if (high) {
        stage(&b->st[i], g, sk, 0.0f, 0.0f, 1.0f);
      } else {
        stage(&b->st[i], g, sk, 1.0f, -sk, -1.0f);
      }
    }
    b->stages = UDSP_EQ_STAGES;
    break;
  }
  case UDSP_EQ_LOW_CUT_12:
    stage(&b->st[0], g, k, 1.0f, -k, -1.0f);
    break;
  case UDSP_EQ_HIGH_CUT_12:
    stage(&b->st[0], g, k, 0.0f, 0.0f, 1.0f);
    break;
  case UDSP_EQ_LOW_SHELF:
    stage(&b->st[0], g / sqrtf(a), k, 1.0f, k * (a - 1.0f), a * a - 1.0f);
    break;
  case UDSP_EQ_HIGH_SHELF:
    stage(&b->st[0], g * sqrtf(a), k, a * a, k * (1.0f - a) * a, 1.0f - a * a);
    break;
  case UDSP_EQ_NOTCH:
    stage(&b->st[0], g, k, 1.0f, -k, 0.0f);
    break;
  default: {
    /* the bell: its bandwidth kept as the gain changes (Q / A) */
    const float kb = 1.0f / (q * a);
    stage(&b->st[0], g, kb, 1.0f, kb * (a * a - 1.0f), 0.0f);
    break;
  }
  }
}

void udsp_eq_ramp(UDSP_EQ_BAND *b, const UDSP_EQ_BAND *to, unsigned frames) {
  if (to->stages != b->stages || frames == 0) {
    for (int k = 0; k < UDSP_EQ_STAGES; k++) {
      b->st[k] = to->st[k];
    }
    b->stages = to->stages;
    b->ramp = 0;
    return;
  }
  const float n = (float)frames;
  for (int k = 0; k < b->stages; k++) {
    const UDSP_EQ_STAGE *x = &b->st[k], *y = &to->st[k];
    b->step[k][0] = (y->g - x->g) / n;
    b->step[k][1] = (y->k - x->k) / n;
    b->step[k][2] = (y->m0 - x->m0) / n;
    b->step[k][3] = (y->m1 - x->m1) / n;
    b->step[k][4] = (y->m2 - x->m2) / n;
    b->end[k] = *y;
  }
  b->ramp = frames;
}

/* one sample through a section's channel: state in, state out */
static inline float step(const UDSP_EQ_STAGE *s, float v0, float *ic1,
                         float *ic2) {
  const float v3 = v0 - *ic2;
  const float v1 = s->a1 * *ic1 + s->a2 * v3;
  const float v2 = *ic2 + s->a2 * *ic1 + s->a3 * v3;
  *ic1 = 2.0f * v1 - *ic1;
  *ic2 = 2.0f * v2 - *ic2;
  return s->m0 * v0 + s->m1 * v1 + s->m2 * v2;
}

/* n frames of a ramp through section k: its coefficients a frame on,
 * then the frame */
UDSP_NO_SLP static void ramp_stage(UDSP_EQ_BAND *b, int k, float *io,
                                   size_t n) {
  UDSP_EQ_STAGE s = b->st[k];
  const float *d = b->step[k];
  float l1 = b->ic1[k][0], l2 = b->ic2[k][0];
  float r1 = b->ic1[k][1], r2 = b->ic2[k][1];
  for (size_t i = 0; i < n; i++) {
    s.g += d[0];
    s.k += d[1];
    s.m0 += d[2];
    s.m1 += d[3];
    s.m2 += d[4];
    s.a1 = 1.0f / (1.0f + s.g * (s.g + s.k));
    s.a2 = s.g * s.a1;
    s.a3 = s.g * s.a2;
    io[i * 2] = step(&s, io[i * 2], &l1, &l2);
    io[i * 2 + 1] = step(&s, io[i * 2 + 1], &r1, &r2);
  }
  b->st[k] = s;
  b->ic1[k][0] = l1;
  b->ic2[k][0] = l2;
  b->ic1[k][1] = r1;
  b->ic2[k][1] = r2;
}

/* a section at a time over the block, both channels in one loop (svf.c's
 * reasoning: two serial chains overlap) */
UDSP_NO_SLP void udsp_eq_stereo(UDSP_EQ_BAND *b, float *io, size_t frames) {
  if (b->ramp) {
    const size_t n = frames < b->ramp ? frames : b->ramp;
    for (int k = 0; k < b->stages; k++) {
      ramp_stage(b, k, io, n);
    }
    b->ramp -= (unsigned)n;
    if (b->ramp == 0) {
      /* the last frame's on the target exactly */
      for (int k = 0; k < b->stages; k++) {
        b->st[k] = b->end[k];
      }
    }
    io += 2 * n;
    frames -= n;
  }
  for (int k = 0; k < b->stages; k++) {
    const UDSP_EQ_STAGE s = b->st[k];
    float l1 = b->ic1[k][0], l2 = b->ic2[k][0];
    float r1 = b->ic1[k][1], r2 = b->ic2[k][1];
    for (size_t i = 0; i < frames; i++) {
      io[i * 2] = step(&s, io[i * 2], &l1, &l2);
      io[i * 2 + 1] = step(&s, io[i * 2 + 1], &r1, &r2);
    }
    b->ic1[k][0] = l1;
    b->ic2[k][0] = l2;
    b->ic1[k][1] = r1;
    b->ic2[k][1] = r2;
  }
}

#if UDSP_W > 1
/* A wavefront over sections (disperse.c's, two lanes a section): lanes
 * 2j and 2j + 1 are section j's left and right, L j frames behind
 * section 0, its input the output section j - 1 made L steps before.
 * Each section runs its frames in order through step()'s operations, so
 * the bits are a section at a time's, but a step's chains are all
 * independent. With L = 1 a step waits on the one before it through the
 * hand-off (a section's output, six operations after its input, shifted
 * one section along: ~19 cycles); L = 2 gives it two steps, so what's
 * left is a section's own recurrence (~15). The channels side by side
 * make the hand-off one shuffle per vector. Lanes past the last section
 * have zero coefficients and state and stay zero; at the block's ends
 * (the first and last L (m - 1) steps) the sections without a frame of
 * their own keep their state. */
#define CHAIN_SECTIONS 16      /* a wavefront's sections at most */
#define CHAIN_PER (UDSP_W / 2) /* sections a vector */
#define CHAIN_MAX_G (CHAIN_SECTIONS / CHAIN_PER)
/* L by the vectors the wavefront takes (G). Two keeps two steps' outputs
 * besides the state: 4 G vectors live, 3 G with one. Past three vectors
 * that spills on SSE's 16 registers and costs more than the wait it
 * saves, so one there; AVX keeps two and, past three, the outputs in
 * memory on purpose (CHAIN_RING_G: a store's forwarding fits in two
 * steps, and the state keeps the registers). bench/README.md "session
 * 262" has the sweep; -DCHAIN_LAG=N forces a lag. */
#ifdef CHAIN_LAG
#define CHAIN_LAG_OF(G) CHAIN_LAG
#elif UDSP_W == 8
#define CHAIN_LAG_OF(G) 2
#else
#define CHAIN_LAG_OF(G) ((G) <= 3 ? 2 : 1)
#endif
#if !defined(CHAIN_RING_G) && UDSP_W == 8
#define CHAIN_RING_G 3
#endif

#if defined(__clang__) || __GNUC__ >= 12
#if UDSP_W == 8
#define CHAIN_SHIFT(a, b)                                                      \
  __builtin_shufflevector(a, b, 6, 7, 8, 9, 10, 11, 12, 13)
#else
#define CHAIN_SHIFT(a, b) __builtin_shufflevector(a, b, 2, 3, 4, 5)
#endif
#else
#if UDSP_W == 8
#define CHAIN_SHIFT(a, b)                                                      \
  __builtin_shuffle(a, b, (UDSP_VI){6, 7, 8, 9, 10, 11, 12, 13})
#else
#define CHAIN_SHIFT(a, b) __builtin_shuffle(a, b, (UDSP_VI){2, 3, 4, 5})
#endif
#endif

typedef struct {
  const UDSP_EQ_STAGE *st;
  float *ic1, *ic2; /* [2]: left, right */
} SECTION;

typedef struct {
  UDSP_VF a1[CHAIN_MAX_G], a2[CHAIN_MAX_G], a3[CHAIN_MAX_G];
  UDSP_VF m0[CHAIN_MAX_G], m1[CHAIN_MAX_G], m2[CHAIN_MAX_G];
  UDSP_VI frame0[CHAIN_MAX_G]; /* each lane's lag: its frame is n - this */
} WAVE;

/* the frame's pair into lanes 0 and 1, the vector's first section's input:
 * x[0] its lanes and the first sections' outputs */
#if defined(__clang__) || __GNUC__ >= 12
#if UDSP_W == 8
#define CHAIN_IN(a, b) __builtin_shufflevector(a, b, 0, 1, 8, 9, 10, 11, 12, 13)
#else
#define CHAIN_IN(a, b) __builtin_shufflevector(a, b, 0, 1, 4, 5)
#endif
#else
#if UDSP_W == 8
#define CHAIN_IN(a, b)                                                         \
  __builtin_shuffle(a, b, (UDSP_VI){0, 1, 8, 9, 10, 11, 12, 13})
#else
#define CHAIN_IN(a, b) __builtin_shuffle(a, b, (UDSP_VI){0, 1, 4, 5})
#endif
#endif

/* a stereo frame as one 8-byte load or store */
typedef float CHAIN_V2 __attribute__((vector_size(8), aligned(4), may_alias));

/* the block's three parts: the head (frames in, none out yet, the later
 * sections idle), the middle (every lane its own frame) and the tail
 * (no frames in, the earlier sections idle) */
enum { CHAIN_MID, CHAIN_HEAD, CHAIN_TAIL };

/* step n: the sections a frame on each, out the outputs made L steps
 * before (overwritten with this step's) */
static inline __attribute__((always_inline)) void
wave_step(const WAVE *w, UDSP_VF *s1, UDSP_VF *s2, UDSP_VF *out, float *io,
          size_t frames, size_t n, int m, int last_g, int last_j,
          const int part, const int G, const int L) {
  UDSP_VF in = {0};
  if (part != CHAIN_TAIL) {
    const CHAIN_V2 p = *(const CHAIN_V2 *)(io + 2 * n);
    in[0] = p[0];
    in[1] = p[1];
  }
  UDSP_VF x[CHAIN_MAX_G];
  x[0] = CHAIN_IN(in, out[0]);
  UDSP_UNROLL
  for (int g = 1; g < G; g++) {
    x[g] = CHAIN_SHIFT(out[g - 1], out[g]);
  }
  UDSP_UNROLL
  for (int g = 0; g < G; g++) {
    const UDSP_VF v3 = x[g] - s2[g];
    const UDSP_VF v1 = w->a1[g] * s1[g] + w->a2[g] * v3;
    const UDSP_VF v2 = s2[g] + w->a2[g] * s1[g] + w->a3[g] * v3;
    const UDSP_VF n1 = 2.0f * v1 - s1[g], n2 = 2.0f * v2 - s2[g];
    out[g] = w->m0[g] * x[g] + w->m1[g] * v1 + w->m2[g] * v2;
    if (part != CHAIN_MID) {
      const UDSP_VI f = (int)n - w->frame0[g];
      const UDSP_VI ok = (f >= 0) & (f < (int)frames);
      s1[g] = udsp_vsel(ok, n1, s1[g]);
      s2[g] = udsp_vsel(ok, n2, s2[g]);
    } else {
      s1[g] = n1;
      s2[g] = n2;
    }
  }
  const size_t behind = (size_t)L * (size_t)(m - 1);
  if (part == CHAIN_MID || (part == CHAIN_TAIL && n >= behind)) {
    const UDSP_VF o = out[last_g];
    *(CHAIN_V2 *)(io + 2 * (n - behind)) =
        (CHAIN_V2){o[2 * last_j], o[2 * last_j + 1]};
  }
}

/* steps [n, end), all in one part: L at a time, so each step's outputs are
 * a fixed buffer (n % L: L is 1 or 2), never an index known only at run
 * time, which keeps the buffers in memory (clang) */
static inline __attribute__((always_inline)) size_t
wave_run(const WAVE *w, UDSP_VF *s1, UDSP_VF *s2, UDSP_VF (*out)[CHAIN_MAX_G],
         float *io, size_t frames, size_t n, size_t end, int m, int last_g,
         int last_j, const int part, const int G, const int L) {
  if (L == 2 && n < end && n % 2) {
    wave_step(w, s1, s2, out[1], io, frames, n, m, last_g, last_j, part, G, L);
    n++;
  }
  for (; n + (size_t)L <= end; n += (size_t)L) {
#ifdef CHAIN_RING_G
    if (G > CHAIN_RING_G) {
      /* out's address hidden: its vectors stored and loaded each step */
      __asm__("" : "+r"(out));
    }
#endif
    UDSP_UNROLL
    for (int l = 0; l < L; l++) {
      wave_step(w, s1, s2, out[l], io, frames, n + (size_t)l, m, last_g, last_j,
                part, G, L);
    }
  }
  if (n < end) {
    /* L = 2's last, n even */
    wave_step(w, s1, s2, out[0], io, frames, n, m, last_g, last_j, part, G, L);
    n++;
  }
  return n;
}

static inline __attribute__((always_inline)) void
wavefront(const SECTION *sec, int m, float *io, size_t frames, const int G) {
  const int L = CHAIN_LAG_OF(G);
  WAVE w;
  UDSP_VF s1[CHAIN_MAX_G], s2[CHAIN_MAX_G], out[2][CHAIN_MAX_G];
  /* the lanes as floats first, then a vector load each: lanes stored one
   * at a time into a vector in memory don't forward to its load (clang
   * kept w there and stalled on each, ~0.3 us a call from G = 2) */
  float f[8][CHAIN_MAX_G * UDSP_W];
  int fr[CHAIN_MAX_G * UDSP_W];
  for (int i = 0; i < G * UDSP_W; i++) {
    const int k = i / 2, c = i % 2;
    const UDSP_EQ_STAGE z = {0}, *s = k < m ? sec[k].st : &z;
    f[0][i] = s->a1;
    f[1][i] = s->a2;
    f[2][i] = s->a3;
    f[3][i] = s->m0;
    f[4][i] = s->m1;
    f[5][i] = s->m2;
    f[6][i] = k < m ? sec[k].ic1[c] : 0.0f;
    f[7][i] = k < m ? sec[k].ic2[c] : 0.0f;
    fr[i] = L * k;
  }
  for (int g = 0; g < G; g++) {
    const int at = g * UDSP_W;
    w.a1[g] = udsp_vload(f[0] + at);
    w.a2[g] = udsp_vload(f[1] + at);
    w.a3[g] = udsp_vload(f[2] + at);
    w.m0[g] = udsp_vload(f[3] + at);
    w.m1[g] = udsp_vload(f[4] + at);
    w.m2[g] = udsp_vload(f[5] + at);
    s1[g] = udsp_vload(f[6] + at);
    s2[g] = udsp_vload(f[7] + at);
    memcpy(&w.frame0[g], fr + at, sizeof w.frame0[g]);
    out[0][g] = out[1][g] = (UDSP_VF){0};
  }
  const int last_g = (m - 1) / CHAIN_PER, last_j = (m - 1) % CHAIN_PER;
  const size_t behind = (size_t)L * (size_t)(m - 1);
  const size_t steps = frames + behind;
  /* the head's steps, the middle's, the tail's */
  const size_t mid = frames > behind ? behind : frames;
  size_t n = wave_run(&w, s1, s2, out, io, frames, 0, mid, m, last_g, last_j,
                      CHAIN_HEAD, G, L);
  n = wave_run(&w, s1, s2, out, io, frames, n, frames, m, last_g, last_j,
               CHAIN_MID, G, L);
  wave_run(&w, s1, s2, out, io, frames, n, steps, m, last_g, last_j, CHAIN_TAIL,
           G, L);
  for (int g = 0; g < G; g++) {
    udsp_vstore(f[6] + g * UDSP_W, s1[g]);
    udsp_vstore(f[7] + g * UDSP_W, s2[g]);
  }
  for (int k = 0; k < m; k++) {
    for (int c = 0; c < 2; c++) {
      sec[k].ic1[c] = f[6][2 * k + c];
      sec[k].ic2[c] = f[7][2 * k + c];
    }
  }
}

/* the sections gathered so far, one wavefront (or a few: CHAIN_SECTIONS
 * at a time) */
static void run_sections(const SECTION *sec, int m, float *io, size_t frames) {
  while (m > 0) {
    const int k = m < CHAIN_SECTIONS ? m : CHAIN_SECTIONS;
    const int g = (k + CHAIN_PER - 1) / CHAIN_PER;
    if (g == 1) {
      wavefront(sec, k, io, frames, 1);
    } else if (g == 2) {
      wavefront(sec, k, io, frames, 2);
    } else if (g == 3) {
      wavefront(sec, k, io, frames, 3);
#if CHAIN_MAX_G > 4 /* SSE: 8; AVX: 4, the default's */
    } else if (g == 4) {
      wavefront(sec, k, io, frames, 4);
    } else if (g <= 6) {
      wavefront(sec, k, io, frames, 6);
#endif
    } else {
      wavefront(sec, k, io, frames, CHAIN_MAX_G);
    }
    sec += k;
    m -= k;
  }
}
#endif

void udsp_eq_chain(UDSP_EQ_BAND *const *bands, int n, float *io,
                   size_t frames) {
#if UDSP_W > 1
  if (frames == 0) {
    return;
  }
  SECTION sec[CHAIN_SECTIONS];
  int m = 0;
  for (int i = 0; i < n; i++) {
    UDSP_EQ_BAND *b = bands[i];
    if (b->ramp) {
      run_sections(sec, m, io, frames);
      m = 0;
      udsp_eq_stereo(b, io, frames);
      continue;
    }
    for (int k = 0; k < b->stages; k++) {
      if (m == CHAIN_SECTIONS) {
        run_sections(sec, m, io, frames);
        m = 0;
      }
      sec[m].st = &b->st[k];
      sec[m].ic1 = b->ic1[k];
      sec[m].ic2 = b->ic2[k];
      m++;
    }
  }
  run_sections(sec, m, io, frames);
#else
  for (int i = 0; i < n; i++) {
    udsp_eq_stereo(bands[i], io, frames);
  }
#endif
}

/* The trapezoidal SVF is its analog prototype at a warped frequency:
 * H = m0 + (m1 s + m2) / (s^2 + k s + 1), s = j tan(pi f / fs) / g */
float udsp_eq_response(const UDSP_EQ_BAND *b, float hz, unsigned rate) {
  float t = hz / (float)rate;
  t = t < 0.0f ? 0.0f : (t > 0.4999f ? 0.4999f : t);
  const double w = tan_rate(t);
  double mag = 1.0;
  for (int k = 0; k < b->stages; k++) {
    const UDSP_EQ_STAGE *s = &b->st[k];
    const double x = w / s->g; /* s = j x */
    /* denominator 1 - x^2 + j k x, numerator m2 + j m1 x */
    const double dr = 1.0 - x * x, di = s->k * x;
    const double nr = s->m2, ni = s->m1 * x;
    const double den = dr * dr + di * di;
    /* m0 + n / d */
    const double hr = s->m0 + (nr * dr + ni * di) / den;
    const double hi = (ni * dr - nr * di) / den;
    mag *= sqrt(hr * hr + hi * hi);
  }
  return (float)mag;
}
