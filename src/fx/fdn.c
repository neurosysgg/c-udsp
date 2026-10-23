#include "fx/fdn.h"

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

size_t udsp_fdn_floats(unsigned rate) {
  size_t n = udsp_line_size(ms_frames(UDSP_FDN_PREDELAY_MAX_MS, rate) + 1);
  for (int i = 0; i < UDSP_FDN_N; i++) {
    n += udsp_line_size(line_max(i, rate));
  }
  for (int i = 0; i < UDSP_FDN_DIFFUSERS; i++) {
    n += udsp_line_size(ms_frames(AP_MS[i], rate) + 1);
  }
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
  for (int i = 0; i < UDSP_FDN_N; i++) {
    f->lp[i] = 0.0f;
  }
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
  }
  for (int i = 0; i < UDSP_FDN_DIFFUSERS; i++) {
    udsp_line_clear(&f->ap[i]);
  }
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
  for (int i = 0; i < UDSP_FDN_N; i++) {
    f->len[i] = LEN_MS[i] * size * (float)rate / 1000.0f;
    /* 60 dB over decay_s: 10^(-3 len / (decay rate)) a pass */
    f->g[i] = (float)udsp_exp2(-3.0 * 3.321928094887362 * f->len[i] /
                               ((double)decay_s * rate));
  }
  /* a one-pole low-pass's pole for the corner: 1 - e^(-2 pi fc / fs) */
  float fc = damp_hz / (float)rate;
  fc = fc > 0.45f ? 0.45f : (fc < 0.0005f ? 0.0005f : fc);
  f->damp =
      1.0f - (float)udsp_exp2(-6.283185307179586 * 1.4426950408889634 * fc);
  f->ap_g =
      0.75f * (diffusion < 0.0f ? 0.0f : (diffusion > 1.0f ? 1.0f : diffusion));
  uint32_t pd = ms_frames(predelay_ms, rate);
  const uint32_t pd_max = ms_frames(UDSP_FDN_PREDELAY_MAX_MS, rate);
  f->predelay = pd > pd_max ? pd_max : pd;
  f->width = width < 0.0f ? 0.0f : (width > 1.0f ? 1.0f : width);
}

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

/* the input's signs into each line, and the two outputs' sums: rows of
 * a Hadamard matrix other than the first, so left and right are
 * uncorrelated */
static const float IN_SIGN[UDSP_FDN_N] = {1, -1, 1, -1, 1, -1, 1, -1};
static const float L_SIGN[UDSP_FDN_N] = {1, 1, -1, -1, 1, 1, -1, -1};
static const float R_SIGN[UDSP_FDN_N] = {1, -1, -1, 1, 1, -1, -1, 1};
/* each line's swing, a phase of the one modulation (cos, sin weights) */
static const float MOD_C[UDSP_FDN_N] = {1.0f,  0.7071f,  0.0f, -0.7071f,
                                        -1.0f, -0.7071f, 0.0f, 0.7071f};
static const float MOD_S[UDSP_FDN_N] = {0.0f, 0.7071f,  1.0f,  0.7071f,
                                        0.0f, -0.7071f, -1.0f, -0.7071f};

void udsp_fdn_stereo(UDSP_FDN *f, const float *in, float *wet, size_t frames) {
  float mc = f->mc, ms = f->ms;
  const float g_ap = f->ap_g;
  for (size_t n = 0; n < frames; n++) {
    float x = 0.5f * (in[2 * n] + in[2 * n + 1]);
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
    /* the modulation turned a frame on; back onto the circle every 256
     * frames from the start (the rotation's rounding would drift it),
     * so where calls cut the frames doesn't change the bits */
    const float c = mc * f->rot_c - ms * f->rot_s;
    ms = ms * f->rot_c + mc * f->rot_s;
    mc = c;
    if ((++f->t & 255u) == 0) {
      const float k = 1.5f - 0.5f * (mc * mc + ms * ms);
      mc *= k;
      ms *= k;
    }
    float s[UDSP_FDN_N];
    float l = 0.0f, r = 0.0f;
    for (int i = 0; i < UDSP_FDN_N; i++) {
      const float d =
          f->len[i] + f->depth * (1.0f + MOD_C[i] * mc + MOD_S[i] * ms);
      const float y = udsp_line_read(&f->line[i], d);
      f->lp[i] += f->damp * (y - f->lp[i]);
      s[i] = f->lp[i];
      l += L_SIGN[i] * s[i];
      r += R_SIGN[i] * s[i];
    }
    hadamard(s);
    for (int i = 0; i < UDSP_FDN_N; i++) {
      udsp_line_push(&f->line[i], IN_SIGN[i] * x + f->g[i] * s[i]);
    }
    /* width: mid kept, side scaled */
    const float mid = 0.5f * (l + r), side = 0.5f * (l - r) * f->width;
    wet[2 * n] = 0.35355339f * (mid + side);
    wet[2 * n + 1] = 0.35355339f * (mid - side);
  }
  f->mc = mc;
  f->ms = ms;
}
