#include "fx/fdn.h"

#include "base/simd.h"
#include "synth/osc.h"
#include "synth/pitch.h"

/* the lines' lengths at size 1, ms: mutually prime-ish, spread over 30-75
 * so their echoes don't line up */
static const float LEN_MS[UDSP_FDN_N] = {29.7f, 37.1f, 41.1f, 43.7f,
                                         53.3f, 59.9f, 67.7f, 73.1f};
/* the diffusers' (Dattorro's input allpasses' at 29.8 kHz, in ms) */
static const float AP_MS[UDSP_FDN_DIFFUSERS] = {4.77f, 3.59f, 12.73f, 9.31f};
/* the most the reads swing, ms */
#define MOD_MS 0.6f
/* the modulation's rate, Hz */
#define MOD_HZ 0.7f

static uint32_t ms_frames(float ms, unsigned rate) {
  return (uint32_t)((double)ms * rate / 1000.0);
}

static uint32_t line_max(int i, unsigned rate) {
  /* the swing reaches twice its depth past the length */
  return ms_frames(LEN_MS[i] * UDSP_FDN_SIZE_MAX + 2.0f * MOD_MS, rate) + 4;
}

static uint32_t tap_max(unsigned rate) {
  return ms_frames(UDSP_FDN_TAP_MAX_MS * UDSP_FDN_SIZE_MAX, rate) + 1;
}

size_t udsp_fdn_floats(unsigned rate) {
  size_t n = udsp_line_size(ms_frames(UDSP_FDN_PREDELAY_MAX_MS, rate) + 1);
  for (int i = 0; i < UDSP_FDN_N; i++) {
    n += udsp_line_size(line_max(i, rate));
  }
  for (int i = 0; i < UDSP_FDN_DIFFUSERS; i++) {
    n += udsp_line_size(ms_frames(AP_MS[i], rate) + 1);
  }
  n += 2 * (size_t)udsp_line_size(tap_max(rate));
  return n;
}

void udsp_fdn_init(UDSP_FDN *f, float *mem, unsigned rate) {
  f->rate = rate;
  uint32_t n = udsp_line_size(ms_frames(UDSP_FDN_PREDELAY_MAX_MS, rate) + 1);
  udsp_line_init(&f->pre, mem, n);
  mem += n;
  for (int i = 0; i < UDSP_FDN_N; i++) {
    n = udsp_line_size(line_max(i, rate));
    udsp_line_init(&f->line[i], mem, n);
    mem += n;
  }
  for (int i = 0; i < UDSP_FDN_DIFFUSERS; i++) {
    n = udsp_line_size(ms_frames(AP_MS[i], rate) + 1);
    udsp_line_init(&f->ap[i], mem, n);
    f->ap_len[i] = ms_frames(AP_MS[i], rate);
    f->ap_len[i] = f->ap_len[i] < 1 ? 1 : f->ap_len[i];
    mem += n;
  }
  n = udsp_line_size(tap_max(rate));
  udsp_line_init(&f->in[0], mem, n);
  udsp_line_init(&f->in[1], mem + n, n);
  f->taps = 0;
  f->early = f->keep = f->feed = 0.0f;
  f->diffuse = 1.0f;
  for (int i = 0; i < UDSP_FDN_N; i++) {
    f->lp[i] = 0.0f;
    f->hf[i] = 0.0f;
    f->lp2[i] = 0.0f;
    f->hf2[i] = 1.0f;
  }
  f->shelf = 0;
  f->shelf_db_s = 0.0f;
  f->shelf2_hz = 1000.0f;
  f->shelf2_db_s = 0.0f;
  f->damp2 = 0.0f;
  f->mc = 1.0f;
  f->ms = 0.0f;
  f->t = 0;
  /* a turn every rate / MOD_HZ frames */
  const float t = MOD_HZ / (float)rate;
  f->rot_c = udsp_osc_sine(t + 0.25f);
  f->rot_s = udsp_osc_sine(t);
  udsp_fdn_set(f, 1.0f, 2.0f, 6000.0f, 0.6f, 0.3f, 0.0f, 1.0f);
}

void udsp_fdn_clear(UDSP_FDN *f) {
  udsp_line_clear(&f->pre);
  for (int i = 0; i < UDSP_FDN_N; i++) {
    udsp_line_clear(&f->line[i]);
    f->lp[i] = 0.0f;
    f->lp2[i] = 0.0f;
  }
  for (int i = 0; i < UDSP_FDN_DIFFUSERS; i++) {
    udsp_line_clear(&f->ap[i]);
  }
  udsp_line_clear(&f->in[0]);
  udsp_line_clear(&f->in[1]);
}

void udsp_fdn_set(UDSP_FDN *f, float size, float decay_s, float damp_hz,
                  float diffusion, float modulation, float predelay_ms,
                  float width) {
  const unsigned rate = f->rate;
  size = size < 0.2f ? 0.2f
                     : (size > UDSP_FDN_SIZE_MAX ? UDSP_FDN_SIZE_MAX : size);
  decay_s = decay_s < 0.05f ? 0.05f : decay_s;
  modulation =
      modulation < 0.0f ? 0.0f : (modulation > 1.0f ? 1.0f : modulation);
  f->depth = modulation * MOD_MS * (float)rate / 1000.0f;
  f->size = size;
  for (int k = 0; k < f->taps; k++) {
    const float ms = f->tap[k].ms * size;
    uint32_t d = (uint32_t)(ms * (float)rate / 1000.0f + 0.5f);
    f->tap_at[k] = d < 1 ? 1 : (d > tap_max(rate) - 1 ? tap_max(rate) - 1 : d);
  }
  for (int i = 0; i < UDSP_FDN_N; i++) {
    f->len[i] = LEN_MS[i] * size * (float)rate / 1000.0f;
    if (f->shelf) {
      /* whole frames: a fraction's cubic would dull the highs a pass */
      f->len[i] = (float)(uint32_t)(f->len[i] + 0.5f);
    }
    /* 60 dB over decay_s: 10^(-3 len / (decay rate)) a pass */
    f->g[i] = (float)udsp_exp2(-3.0 * 3.321928094887362 * f->len[i] /
                               ((double)decay_s * rate));
    /* the shelf's highs: db_s a second, len / rate of one a pass */
    f->hf[i] = f->shelf ? (float)udsp_exp2(-(double)f->shelf_db_s / 20.0 *
                                           3.321928094887362 * f->len[i] / rate)
                        : 0.0f;
    f->hf2[i] = (float)udsp_exp2(-(double)f->shelf2_db_s / 20.0 *
                                 3.321928094887362 * f->len[i] / rate);
  }
  /* a one-pole low-pass's pole for the corner: 1 - e^(-2 pi fc / fs) */
  float fc = damp_hz / (float)rate;
  fc = fc > 0.45f ? 0.45f : (fc < 0.0005f ? 0.0005f : fc);
  f->damp =
      1.0f - (float)udsp_exp2(-6.283185307179586 * 1.4426950408889634 * fc);
  float fc2 = f->shelf2_hz / (float)rate;
  fc2 = fc2 > 0.45f ? 0.45f : (fc2 < 0.0005f ? 0.0005f : fc2);
  f->damp2 =
      1.0f - (float)udsp_exp2(-6.283185307179586 * 1.4426950408889634 * fc2);
  f->ap_g =
      0.75f * (diffusion < 0.0f ? 0.0f : (diffusion > 1.0f ? 1.0f : diffusion));
  uint32_t pd = ms_frames(predelay_ms, rate);
  const uint32_t pd_max = ms_frames(UDSP_FDN_PREDELAY_MAX_MS, rate);
  f->predelay = pd > pd_max ? pd_max : pd;
  f->width = width < 0.0f ? 0.0f : (width > 2.0f ? 2.0f : width);
}

void udsp_fdn_set_diffuse(UDSP_FDN *f, float gain) { f->diffuse = gain; }

void udsp_fdn_set_early(UDSP_FDN *f, const UDSP_FDN_TAP *taps, int n,
                        float early, float keep, float feed) {
  n = n < 0 ? 0 : (n > UDSP_FDN_TAPS ? UDSP_FDN_TAPS : n);
  f->taps = n;
  for (int k = 0; k < n; k++) {
    f->tap[k] = taps[k];
    const float ms = taps[k].ms * f->size;
    uint32_t d = (uint32_t)(ms * (float)f->rate / 1000.0f + 0.5f);
    f->tap_at[k] =
        d < 1 ? 1 : (d > tap_max(f->rate) - 1 ? tap_max(f->rate) - 1 : d);
  }
  f->early = early < 0.0f ? 0.0f : (early > 1.0f ? 1.0f : early);
  f->keep = keep;
  f->feed = feed;
}

void udsp_fdn_set_shelf(UDSP_FDN *f, int on, float db_s, float hz2,
                        float db2_s) {
  f->shelf = on;
  f->shelf_db_s = db_s > 0.0f ? db_s : 0.0f;
  f->shelf2_hz = hz2;
  f->shelf2_db_s = on && db2_s > 0.0f ? db2_s : 0.0f;
}

#if UDSP_W == 1
/* x <- H x / sqrt(8), the fast Walsh-Hadamard transform */
static inline void hadamard(float *x) {
  for (int h = 1; h < UDSP_FDN_N; h <<= 1) {
    for (int i = 0; i < UDSP_FDN_N; i += h << 1) {
      for (int j = i; j < i + h; j++) {
        const float a = x[j], b = x[j + h];
        x[j] = a + b;
        x[j + h] = a - b;
      }
    }
  }
  for (int i = 0; i < UDSP_FDN_N; i++) {
    x[i] *= 0.35355339f;
  }
}
#endif

/* the input's signs into each line, and the two outputs' sums: rows of
 * a Hadamard matrix other than the first, so left and right are
 * uncorrelated */
static const float IN_SIGN[UDSP_FDN_N] = {1, -1, 1, -1, 1, -1, 1, -1};
static const float L_SIGN[UDSP_FDN_N] = {1, 1, -1, -1, 1, 1, -1, -1};
static const float R_SIGN[UDSP_FDN_N] = {1, -1, -1, 1, 1, -1, -1, 1};
/* the taps' row into the lines: another of the Hadamard's */
static const float ALT_SIGN[UDSP_FDN_N] = {1, -1, 1, -1, -1, 1, -1, 1};

/* a frame's taps (er[2]), the frame's input pushed after */
static inline void early_frame(UDSP_FDN *f, float l, float r, float *er) {
  er[0] = er[1] = 0.0f;
  for (int k = 0; k < f->taps; k++) {
    const UDSP_FDN_TAP *t = &f->tap[k];
    er[t->out] += t->gain * udsp_line_tap(&f->in[t->in], f->tap_at[k]);
  }
  udsp_line_push(&f->in[0], l);
  udsp_line_push(&f->in[1], r);
}

/* the outputs' sums l, r and the taps er into a frame of wet: the
 * width's mid and side over both */
static inline void wet_early(float l, float r, const float *er, float early,
                             float width, float diffuse, float *out) {
  const float t = 0.35355339f * 0.5f * diffuse;
  const float mid = t * (l + r) + early * 0.5f * (er[0] + er[1]);
  const float side = (t * (l - r) + early * 0.5f * (er[0] - er[1])) * width;
  out[0] = mid + side;
  out[1] = mid - side;
}

/* each line's swing, a phase of the one modulation (cos, sin weights) */
static const float MOD_C[UDSP_FDN_N] = {1.0f,  0.7071f,  0.0f, -0.7071f,
                                        -1.0f, -0.7071f, 0.0f, 0.7071f};
static const float MOD_S[UDSP_FDN_N] = {0.0f, 0.7071f,  1.0f,  0.7071f,
                                        0.0f, -0.7071f, -1.0f, -0.7071f};

/* the modulation turned a frame on; back onto the circle every 256
 * frames from the start (the rotation's rounding would drift it), so
 * where calls cut the frames doesn't change the bits */
static inline void turn(UDSP_FDN *f, float *mc, float *ms) {
  const float c = *mc * f->rot_c - *ms * f->rot_s;
  *ms = *ms * f->rot_c + *mc * f->rot_s;
  *mc = c;
  if ((++f->t & 255u) == 0) {
    const float k = 1.5f - 0.5f * (*mc * *mc + *ms * *ms);
    *mc *= k;
    *ms *= k;
  }
}

void udsp_fdn_idle(UDSP_FDN *f, size_t frames) {
  float mc = f->mc, ms = f->ms;
  for (size_t n = 0; n < frames; n++) {
    turn(f, &mc, &ms);
  }
  f->mc = mc;
  f->ms = ms;
}

#if UDSP_W > 1
/* the lines in vector lanes (line i in lane i % W of vector i / W): each
 * lane runs the scalar loop's operations in its order -- the read's
 * Hermite, the damping, the Hadamard's butterflies, the push -- and the
 * two outputs' sums stay sums in line order, so the bits are the same
 * as one line at a time */
#define FG (UDSP_FDN_N / UDSP_W)

#if defined(__clang__) || __GNUC__ >= 12
#define PERM(v, ...) __builtin_shufflevector(v, v, __VA_ARGS__)
#else
#define PERM(v, ...) __builtin_shuffle(v, (UDSP_VI){__VA_ARGS__})
#endif

/* v with each lane j and j ^ h swapped */
static inline UDSP_VF swap_lanes(UDSP_VF v, int h) {
#if UDSP_W == 8
  return h == 1   ? PERM(v, 1, 0, 3, 2, 5, 4, 7, 6)
         : h == 2 ? PERM(v, 2, 3, 0, 1, 6, 7, 4, 5)
                  : PERM(v, 4, 5, 6, 7, 0, 1, 2, 3);
#else
  return h == 1 ? PERM(v, 1, 0, 3, 2) : PERM(v, 2, 3, 0, 1);
#endif
}

/* hadamard() over lanes: a butterfly's lower line a + b, its upper a - b */
static inline void vhadamard(UDSP_VF *x) {
  UDSP_VI lane;
  UDSP_UNROLL
  for (int j = 0; j < UDSP_W; j++) {
    lane[j] = j;
  }
  UDSP_UNROLL
  for (int h = 1; h < UDSP_FDN_N; h <<= 1) {
    if (h < UDSP_W) {
      const UDSP_VI upper = (lane & h) != 0;
      UDSP_UNROLL
      for (int g = 0; g < FG; g++) {
        const UDSP_VF sw = swap_lanes(x[g], h);
        x[g] = udsp_vsel(upper, sw - x[g], x[g] + sw);
      }
    } else {
      const int d = h / UDSP_W;
      UDSP_UNROLL
      for (int g = 0; g < FG; g += d << 1) {
        UDSP_UNROLL
        for (int k = g; k < g + d; k++) {
          const UDSP_VF a = x[k], b = x[k + d];
          x[k] = a + b;
          x[k + d] = a - b;
        }
      }
    }
  }
  UDSP_UNROLL
  for (int g = 0; g < FG; g++) {
    x[g] *= 0.35355339f;
  }
}

#endif

#if UDSP_W > 1
/* frames whose pushes wait for the block's end: a read reaches no nearer
 * than its line's length less a rounding of the swing (its circle is
 * kept on 1), so a block shorter than the shortest line reads none of
 * its own pushes. Then each line's run goes in at once and the outputs'
 * sums are taken across frames -- a frame's lanes are never picked
 * apart one by one */
#define FDN_BLOCK 16

static void fdn_lanes(UDSP_FDN *f, const float *in, float *wet, size_t frames) {
  float mc = f->mc, ms = f->ms;
  const float g_ap = f->ap_g;
  /* the lines' rings are one block from udsp_fdn_init(): offsets in it */
  const float *base = f->line[0].buf;
  float *wbase = f->line[0].buf;
  UDSP_VF len[FG], gain[FG], lp[FG], hf[FG], lp2[FG], hf2[FG];
  UDSP_VF mod_c[FG], mod_s[FG];
  UDSP_VF in_sign[FG], alt_sign[FG];
  UDSP_VI off[FG], w[FG], mask[FG];
  float shortest = f->len[0];
  for (int g = 0; g < FG; g++) {
    for (int j = 0; j < UDSP_W; j++) {
      const int i = g * UDSP_W + j;
      len[g][j] = f->len[i];
      gain[g][j] = f->g[i];
      lp[g][j] = f->lp[i];
      hf[g][j] = f->hf[i];
      lp2[g][j] = f->lp2[i];
      hf2[g][j] = f->hf2[i];
      mod_c[g][j] = MOD_C[i];
      mod_s[g][j] = MOD_S[i];
      in_sign[g][j] = IN_SIGN[i];
      alt_sign[g][j] = ALT_SIGN[i];
      off[g][j] = (int)(f->line[i].buf - base);
      w[g][j] = (int)f->line[i].w;
      mask[g][j] = (int)f->line[i].mask;
      shortest = f->len[i] < shortest ? f->len[i] : shortest;
    }
  }
  /* a frame at a time below a few kHz (lines shorter than a block) */
  const size_t room = shortest > FDN_BLOCK + 4 ? FDN_BLOCK : 1;
  const float depth = f->depth, damp = f->damp, width = f->width;
  const float damp2 = f->damp2;
  const int two = f->shelf2_db_s > 0.0f;
  /* a block's frames: each one's lines before the Hadamard (the sums')
   * and its pushes, lines in order */
  UDSP_VF sb[FDN_BLOCK][FG], vb[FDN_BLOCK][FG];
  float erb[FDN_BLOCK][2] = {{0}}; /* the frames' taps */
  const float early = f->early;
  const float *sf = (const float *)sb, *vf = (const float *)vb;
  for (size_t n0 = 0; n0 < frames; n0 += room) {
    const size_t b = frames - n0 < room ? frames - n0 : room;
    for (size_t k = 0; k < b; k++) {
      const size_t n = n0 + k;
      float x = 0.5f * (in[2 * n] + in[2 * n + 1]), xa = 0.0f;
      if (f->taps) {
        early_frame(f, in[2 * n], in[2 * n + 1], erb[k]);
        xa = f->early * f->feed * 0.5f * (erb[k][0] + erb[k][1]);
        x *= 1.0f - f->early + f->early * f->keep;
      }
      if (f->predelay) {
        const float d = udsp_line_tap(&f->pre, f->predelay);
        udsp_line_push(&f->pre, x);
        x = d;
      }
      UDSP_UNROLL
      for (int a = 0; a < UDSP_FDN_DIFFUSERS; a++) {
        const float t = udsp_line_tap(&f->ap[a], f->ap_len[a]);
        const float v = x + g_ap * t;
        udsp_line_push(&f->ap[a], v);
        x = t - g_ap * v;
      }
      turn(f, &mc, &ms);
      UDSP_VF s[FG];
      UDSP_UNROLL
      for (int g = 0; g < FG; g++) {
        const UDSP_VF d =
            len[g] + depth * (1.0f + mod_c[g] * mc + mod_s[g] * ms);
        const UDSP_VI i = __builtin_convertvector(d, UDSP_VI);
        const UDSP_VF fr = d - __builtin_convertvector(i, UDSP_VF);
        UDSP_VF y0, y1, y2;
        /* the next push is k frames on from w: those before it wait */
        const UDSP_VF ym1 =
            udsp_vtaps(base, off[g], (w[g] + (int)k - i) & mask[g], mask[g],
                       &y0, &y1, &y2);
        const UDSP_VF y = udsp_vhermite(ym1, y0, y1, y2, fr);
        lp[g] += damp * (y - lp[g]);
        /* the shelf: the low-pass and hf of the rest (0: the low-pass) */
        s[g] = lp[g] + hf[g] * (y - lp[g]);
        if (two) {
          lp2[g] += damp2 * (s[g] - lp2[g]);
          s[g] = lp2[g] + hf2[g] * (s[g] - lp2[g]);
        }
        sb[k][g] = s[g];
      }
      vhadamard(s);
      UDSP_UNROLL
      for (int g = 0; g < FG; g++) {
        vb[k][g] = in_sign[g] * x + gain[g] * s[g];
        if (f->taps) {
          vb[k][g] += alt_sign[g] * xa;
        }
      }
    }
    /* UDSP_W frames at a time: four lines' values across them (a
     * transpose of the frames' rows) -- the pushes as runs, the sums in
     * line order per frame as one line at a time had them */
    size_t k = 0;
    for (; k + UDSP_W <= b; k += UDSP_W) {
      UDSP_VI at;
      UDSP_UNROLL
      for (int j = 0; j < UDSP_W; j++) {
        at[j] = (int)(k + (size_t)j) * UDSP_FDN_N;
      }
      UDSP_VF c[UDSP_FDN_N], v[UDSP_FDN_N];
      UDSP_UNROLL
      for (int q = 0; q < UDSP_FDN_N; q += 4) {
        udsp_vcols4(sf, at + q, &c[q], &c[q + 1], &c[q + 2], &c[q + 3]);
        udsp_vcols4(vf, at + q, &v[q], &v[q + 1], &v[q + 2], &v[q + 3]);
      }
      UDSP_VF l = udsp_vsplat(0.0f), r = l;
      UDSP_UNROLL
      for (int i = 0; i < UDSP_FDN_N; i++) {
        l += L_SIGN[i] * c[i];
        r += R_SIGN[i] * c[i];
        const int lw = w[i / UDSP_W][i % UDSP_W] + (int)k;
        const int lm = mask[i / UDSP_W][i % UDSP_W];
        float *ring = wbase + off[i / UDSP_W][i % UDSP_W];
        if (lw + UDSP_W - 1 <= lm) {
          udsp_vstore(ring + lw, v[i]);
        } else {
          for (int j = 0; j < UDSP_W; j++) {
            ring[(lw + j) & lm] = v[i][j];
          }
        }
      }
      if (f->taps) {
        UDSP_UNROLL
        for (int j = 0; j < UDSP_W; j++) {
          wet_early(l[j], r[j], erb[k + (size_t)j], early, width, f->diffuse,
                    wet + 2 * (n0 + k + (size_t)j));
        }
      } else {
        l *= f->diffuse;
        r *= f->diffuse;
        const UDSP_VF mid = 0.5f * (l + r), side = 0.5f * (l - r) * width;
        const UDSP_VF wl = 0.35355339f * (mid + side);
        const UDSP_VF wr = 0.35355339f * (mid - side);
        UDSP_UNROLL
        for (int j = 0; j < UDSP_W; j++) {
          wet[2 * (n0 + k + (size_t)j)] = wl[j];
          wet[2 * (n0 + k + (size_t)j) + 1] = wr[j];
        }
      }
    }
    for (; k < b; k++) {
      float l = 0.0f, r = 0.0f;
      for (int i = 0; i < UDSP_FDN_N; i++) {
        const float y = sf[k * UDSP_FDN_N + (size_t)i];
        l += L_SIGN[i] * y;
        r += R_SIGN[i] * y;
        const int lw = w[i / UDSP_W][i % UDSP_W] + (int)k;
        wbase[off[i / UDSP_W][i % UDSP_W] +
              (lw & mask[i / UDSP_W][i % UDSP_W])] =
            vf[k * UDSP_FDN_N + (size_t)i];
      }
      if (f->taps) {
        wet_early(l, r, erb[k], early, width, f->diffuse, wet + 2 * (n0 + k));
      } else {
        l *= f->diffuse;
        r *= f->diffuse;
        const float mid = 0.5f * (l + r), side = 0.5f * (l - r) * width;
        wet[2 * (n0 + k)] = 0.35355339f * (mid + side);
        wet[2 * (n0 + k) + 1] = 0.35355339f * (mid - side);
      }
    }
    UDSP_UNROLL
    for (int g = 0; g < FG; g++) {
      w[g] = (w[g] + (int)b) & mask[g];
    }
  }
  UDSP_UNROLL
  for (int i = 0; i < UDSP_FDN_N; i++) {
    f->lp[i] = lp[i / UDSP_W][i % UDSP_W];
    f->lp2[i] = lp2[i / UDSP_W][i % UDSP_W];
    f->line[i].w = (uint32_t)w[i / UDSP_W][i % UDSP_W];
  }
  f->mc = mc;
  f->ms = ms;
}
#endif

void udsp_fdn_stereo(UDSP_FDN *f, const float *in, float *wet, size_t frames) {
#if UDSP_W > 1
  fdn_lanes(f, in, wet, frames);
#else
  float mc = f->mc, ms = f->ms;
  const float g_ap = f->ap_g;
  for (size_t n = 0; n < frames; n++) {
    float x = 0.5f * (in[2 * n] + in[2 * n + 1]), xa = 0.0f, er[2] = {0};
    if (f->taps) {
      early_frame(f, in[2 * n], in[2 * n + 1], er);
      xa = f->early * f->feed * 0.5f * (er[0] + er[1]);
      x *= 1.0f - f->early + f->early * f->keep;
    }
    if (f->predelay) {
      const float d = udsp_line_tap(&f->pre, f->predelay);
      udsp_line_push(&f->pre, x);
      x = d;
    }
    for (int k = 0; k < UDSP_FDN_DIFFUSERS; k++) {
      const float b = udsp_line_tap(&f->ap[k], f->ap_len[k]);
      const float v = x + g_ap * b;
      udsp_line_push(&f->ap[k], v);
      x = b - g_ap * v;
    }
    turn(f, &mc, &ms);
    float s[UDSP_FDN_N];
    float l = 0.0f, r = 0.0f;
    for (int i = 0; i < UDSP_FDN_N; i++) {
      const float d =
          f->len[i] + f->depth * (1.0f + MOD_C[i] * mc + MOD_S[i] * ms);
      const float y = udsp_line_read(&f->line[i], d);
      f->lp[i] += f->damp * (y - f->lp[i]);
      s[i] = f->lp[i] + f->hf[i] * (y - f->lp[i]);
      if (f->shelf2_db_s > 0.0f) {
        f->lp2[i] += f->damp2 * (s[i] - f->lp2[i]);
        s[i] = f->lp2[i] + f->hf2[i] * (s[i] - f->lp2[i]);
      }
      l += L_SIGN[i] * s[i];
      r += R_SIGN[i] * s[i];
    }
    hadamard(s);
    for (int i = 0; i < UDSP_FDN_N; i++) {
      float v = IN_SIGN[i] * x + f->g[i] * s[i];
      if (f->taps) {
        v += ALT_SIGN[i] * xa;
      }
      udsp_line_push(&f->line[i], v);
    }
    if (f->taps) {
      wet_early(l, r, er, f->early, f->width, f->diffuse, wet + 2 * n);
    } else {
      l *= f->diffuse;
      r *= f->diffuse;
      /* width: mid kept, side scaled */
      const float mid = 0.5f * (l + r), side = 0.5f * (l - r) * f->width;
      wet[2 * n] = 0.35355339f * (mid + side);
      wet[2 * n + 1] = 0.35355339f * (mid - side);
    }
  }
  f->mc = mc;
  f->ms = ms;
#endif
}
