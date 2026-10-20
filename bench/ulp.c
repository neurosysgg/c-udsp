/* ulp: how far each kernel's answers are from the exact ones, on the
 * machine and flags it was built with. Each kernel runs over a corpus of
 * signals (the bench's music, and the edges: silence, DC, an impulse,
 * Nyquist, a full-scale square, noise, quiet, loud, a chirp) and is held
 * to a reference computed in double; the error is counted in ULPs (units
 * in the last place of a float) of a scale the kernel says:
 *
 *   elem   the reference value itself (a correctly rounded answer is 0.5)
 *   input  the input sample's magnitude (a window: its error relative to
 *          what went in, not to a product near zero at the edges)
 *   max    the largest reference output (a transform: its rounding error
 *          spreads over every bin, so a quiet bin is judged against the
 *          loudest, as an FFT's error bound is)
 *   sum    sum |x| (a sum: the error bound of summing, whatever cancels)
 *   exact  no error allowed (peak, min/max, downmix's average)
 *   bits   no reference: only the hash, which matrix.py-style drivers
 *          compare across builds (the synthesis half promises the same
 *          bits everywhere; ctui-daw holds them to µOS's renders)
 *
 * One line per kernel, tab-separated (bench/ulp.py collects the builds):
 *
 *   kernel  metric  max_ulp  worst_signal  tolerance  hash
 *
 * max_ulp is the worst over the corpus, worst_signal where it was; the
 * tolerance is the written-down bound (TOL below). Exits 1 if any kernel
 * is over its tolerance. -k NAME runs only the kernels whose name has
 * NAME in it; -v prints every signal's error too.
 *
 * Subnormal input isn't in the corpus: the kernels don't scale against
 * underflow (a magnitude below ~1e-19 squares to zero), which audio never
 * reaches (-380 dBFS). */
#define _DEFAULT_SOURCE

#include "adsr.h"
#include "analyze.h"
#include "delay.h"
#include "disperse.h"
#include "dyn.h"
#include "eq.h"
#include "fdn.h"
#include "conv.h"
#include "fft.h"
#include "gain.h"
#include "lfo.h"
#include "loudness.h"
#include "openv.h"
#include "osc.h"
#include "pitch.h"
#include "r128.h"
#include "ramp.h"
#include "sampler.h"
#include "scope.h"
#include "shape.h"
#include "sinc.h"
#include "spectrum.h"
#include "svf.h"
#include "wtab.h"
#include "yin.h"
#include "xover.h"
#include "zdf.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXN 8192 /* the longest transform; signals are this long */
#define RMS_BIG (1u << 20)

typedef enum { M_ELEM, M_INPUT, M_MAX, M_SUM, M_EXACT, M_BITS } METRIC;
static const char *const METRIC_NAME[] = {"elem", "input", "max",
                                          "sum",  "exact", "bits"};

/* --- the corpus --- */

#define N_SIG 10
static const char *const SIG_NAME[N_SIG] = {
    "silence", "dc",    "impulse", "nyquist", "square",
    "noise",   "music", "quiet",   "loud",    "chirp"};
static float g_sig[N_SIG][MAXN];
static float *g_big; /* music, RMS_BIG long: rms over many chunks */

static unsigned int g_seed = 12345;
static float noise(void) {
  g_seed = g_seed * 1103515245u + 12345u;
  return (float)((g_seed >> 8) & 0xffff) / 32768.0f - 1.0f;
}

/* bench.c's: a bass tone, a wobbling melody, some noise */
static float music(size_t i) {
  double t = (double)i / 48000.0;
  return (float)(0.4 * sin(2 * M_PI * 55.0 * t) +
                 0.3 * sin(2 * M_PI * (440.0 + 20.0 * sin(t * 3)) * t)) +
         0.05f * noise();
}

static void corpus(void) {
  for (size_t i = 0; i < MAXN; i++) {
    double t = (double)i / 48000.0;
    float m = music(i);
    g_sig[0][i] = 0.0f;
    g_sig[1][i] = 0.5f;
    g_sig[2][i] = i == 0 ? 1.0f : 0.0f;
    g_sig[3][i] = i & 1 ? -1.0f : 1.0f;
    g_sig[4][i] = (i / 50) & 1 ? -1.0f : 1.0f;
    g_sig[5][i] = noise();
    g_sig[6][i] = m;
    g_sig[7][i] = m * 1e-6f;
    g_sig[8][i] = m * 1000.0f;
    /* 20 Hz to 20 kHz over the longest transform, exponentially */
    double k = log(1000.0) / ((double)MAXN / 48000.0);
    g_sig[9][i] = (float)(0.9 * sin(2 * M_PI * 20.0 * (exp(k * t) - 1) / k));
  }
  g_big = malloc(RMS_BIG * sizeof(float));
  for (size_t i = 0; i < RMS_BIG; i++) {
    g_big[i] = music(i);
  }
}

/* --- measuring --- */

/* a float's ULP at |x| (the gap to the next float up; the smallest
 * subnormal at 0) */
static double ulp_at(double x) {
  float f = fabsf((float)x);
  if (isinf(f)) {
    f = FLT_MAX;
  }
  return (double)nextafterf(f, INFINITY) - (double)f;
}

static double err_ulp(float got, double ref, double scale) {
  if (isnan(got) != isnan(ref)) {
    return INFINITY;
  }
  if (isnan(got)) {
    return 0.0;
  }
  return fabs((double)got - ref) / ulp_at(scale);
}

/* the worst error over n outputs, each against the metric's scale;
 * in[] is the input (M_INPUT), s the scale of M_MAX / M_SUM */
static double worst(METRIC m, const float *got, const double *ref,
                    const float *in, double s, size_t n) {
  double e = 0.0;
  for (size_t i = 0; i < n; i++) {
    double scale = m == M_INPUT ? fabs((double)in[i])
                   : m == M_MAX || m == M_SUM ? s
                                              : ref[i];
    double x = m == M_EXACT ? ((double)got[i] == ref[i] ? 0.0 : INFINITY)
                            : err_ulp(got[i], ref[i], scale);
    e = x > e ? x : e;
  }
  return e;
}

static double max_abs(const double *x, size_t n) {
  double m = 0.0;
  for (size_t i = 0; i < n; i++) {
    m = fabs(x[i]) > m ? fabs(x[i]) : m;
  }
  return m;
}

/* FNV-1a over the outputs' bits: equal across builds = the same bits */
static uint64_t g_hash;
static void hash_floats(const float *x, size_t n) {
  const unsigned char *p = (const unsigned char *)x;
  for (size_t i = 0; i < n * sizeof(float); i++) {
    g_hash = (g_hash ^ p[i]) * 0x100000001b3ull;
  }
}

/* --- the references, in double --- */

static double hann_d(size_t i, size_t n) {
  return n < 2 ? 1.0 : 0.5 - 0.5 * cos(2.0 * M_PI * (double)i / (double)(n - 1));
}

/* the DFT of x (n points; xi NULL = real input), bins 0..nout-1, from an
 * exact twiddle table: n^2 work, fine at n <= 8192 */
static void dft_d(const double *xr, const double *xi, size_t n, size_t nout,
                  double *re, double *im) {
  static double c[MAXN], s[MAXN];
  for (size_t j = 0; j < n; j++) {
    c[j] = cos(2.0 * M_PI * (double)j / (double)n);
    s[j] = -sin(2.0 * M_PI * (double)j / (double)n);
  }
  for (size_t k = 0; k < nout; k++) {
    double ar = 0.0, ai = 0.0;
    for (size_t j = 0; j < n; j++) {
      size_t t = (j * k) & (n - 1);
      double yi = xi ? xi[j] : 0.0;
      ar += xr[j] * c[t] - yi * s[t];
      ai += xr[j] * s[t] + yi * c[t];
    }
    re[k] = ar;
    im[k] = ai;
  }
}

/* --- the kernels --- */

typedef struct {
  const char *name;
  METRIC metric;
  double tol;
  /* one signal through the kernel: the worst error (the hash updated) */
  double (*run)(const float *x, size_t n);
  size_t n;
} KERNEL;

static float g_got[2 * MAXN + 32], g_got2[MAXN + 32], g_table[2 * MAXN],
    g_win[MAXN];
static double g_ref[2 * MAXN], g_ref2[MAXN], g_xd[MAXN], g_wd[MAXN];

static void to_d(const float *x, size_t n, double *d) {
  for (size_t i = 0; i < n; i++) {
    d[i] = x[i];
  }
}

static double k_hann(const float *x, size_t n) {
  udsp_window_hann(x, g_got, n);
  hash_floats(g_got, n);
  for (size_t i = 0; i < n; i++) {
    g_ref[i] = (double)x[i] * hann_d(i, n);
  }
  return worst(M_INPUT, g_got, g_ref, x, 0, n);
}

static double k_hann_table(const float *x, size_t n) {
  udsp_window_hann_table(g_win, n);
  udsp_window_apply(x, g_win, g_got, n);
  hash_floats(g_got, n);
  for (size_t i = 0; i < n; i++) {
    g_ref[i] = (double)x[i] * hann_d(i, n);
  }
  return worst(M_INPUT, g_got, g_ref, x, 0, n);
}

/* a complex transform's re and im against the DFT, judged by the largest
 * |bin| component */
static double fft_check(size_t n) {
  hash_floats(g_got, n);
  hash_floats(g_got2, n);
  double s = max_abs(g_ref, n), si = max_abs(g_ref2, n);
  s = si > s ? si : s;
  double a = worst(M_MAX, g_got, g_ref, NULL, s, n);
  double b = worst(M_MAX, g_got2, g_ref2, NULL, s, n);
  return a > b ? a : b;
}

static double k_fft(const float *x, size_t n) {
  memcpy(g_got, x, n * sizeof(float));
  memset(g_got2, 0, n * sizeof(float));
  udsp_fft_forward(g_got, g_got2, n);
  to_d(x, n, g_xd);
  dft_d(g_xd, NULL, n, n, g_ref, g_ref2);
  return fft_check(n);
}

static double k_fft_table(const float *x, size_t n) {
  /* complex input: the signal in re, its reverse halved in im */
  udsp_fft_table(g_table, n);
  memcpy(g_got, x, n * sizeof(float));
  for (size_t i = 0; i < n; i++) {
    g_got2[i] = x[n - 1 - i] * 0.5f;
    g_wd[i] = g_got2[i];
  }
  to_d(x, n, g_xd);
  udsp_fft_forward_table(g_got, g_got2, n, g_table);
  dft_d(g_xd, g_wd, n, n, g_ref, g_ref2);
  return fft_check(n);
}

static double k_fft_real(const float *x, size_t n) {
  udsp_fft_table(g_table, n);
  udsp_fft_real(x, g_got, g_got2, n, g_table);
  to_d(x, n, g_xd);
  dft_d(g_xd, NULL, n, n / 2 + 1, g_ref, g_ref2);
  return fft_check(n / 2 + 1);
}

static double k_fft_real_windowed(const float *x, size_t n) {
  udsp_fft_table(g_table, n);
  udsp_window_hann_table(g_win, n);
  udsp_fft_real_windowed(x, g_win, g_got, g_got2, n, g_table);
  for (size_t i = 0; i < n; i++) {
    g_xd[i] = (double)x[i] * hann_d(i, n);
  }
  dft_d(g_xd, NULL, n, n / 2 + 1, g_ref, g_ref2);
  return fft_check(n / 2 + 1);
}

/* the inverse of udsp_fft_real()'s bins of x, against the exact inverse
 * of those same bins (the forward's own error out of it), judged by the
 * largest |sample| */
static double k_fft_real_inverse(const float *x, size_t n) {
  static float re[MAXN / 2 + 1], im[MAXN / 2 + 1];
  udsp_fft_table(g_table, n);
  udsp_fft_real(x, re, im, n, g_table);
  const size_t m = n / 2;
  for (size_t j = 0; j < n; j++) {
    double v = re[0] + (j & 1 ? -1.0 : 1.0) * re[m];
    for (size_t k = 1; k < m; k++) {
      double a = 2.0 * M_PI * (double)((j * k) & (n - 1)) / (double)n;
      v += 2.0 * ((double)re[k] * cos(a) - (double)im[k] * sin(a));
    }
    g_ref[j] = v;
  }
  udsp_fft_real_inverse(re, im, g_got, n, g_table);
  hash_floats(g_got, n);
  return worst(M_MAX, g_got, g_ref, NULL, max_abs(g_ref, n), n);
}

/* x through a 4-partition response (x's own reverse, decaying) in blocks
 * of 64 */
static double k_conv(const float *x, size_t n) {
  enum { B = 64, P = 4, S = UDSP_CONV_SPECTRUM(B) };
  static float h[P * B], spectra[P * S], line[P * S], acc[S], frame[2 * B],
      scratch[2 * B];
  for (size_t i = 0; i < P * B; i++) {
    h[i] = x[n - 1 - i] * (float)exp(-(double)i / 64.0);
  }
  udsp_fft_table(g_table, 2 * B);
  udsp_conv_ir(h, P * B, 1, B, spectra, g_table, scratch);
  UDSP_CONV_IN in;
  udsp_conv_in_init(&in, B, P, line, frame);
  for (size_t at = 0; at + B <= n; at += B) {
    udsp_conv_in_push(&in, x + at, 1, g_table);
    memset(acc, 0, sizeof acc);
    udsp_conv_mac(&in, spectra, P, acc);
    udsp_conv_out(acc, B, g_got + at, 1, g_table, scratch);
  }
  hash_floats(g_got, n);
  return 0.0;
}

/* the input as one transform's re and im halves */
static double k_power(const float *x, size_t n) {
  size_t h = n / 2;
  udsp_fft_power(x, x + h, h, g_got);
  hash_floats(g_got, h);
  for (size_t i = 0; i < h; i++) {
    g_ref[i] = (double)x[i] * x[i] + (double)x[h + i] * x[h + i];
  }
  return worst(M_ELEM, g_got, g_ref, NULL, 0, h);
}

static double k_magnitude(const float *x, size_t n) {
  size_t h = n / 2;
  udsp_fft_magnitude(x, x + h, h, g_got);
  hash_floats(g_got, h);
  for (size_t i = 0; i < h; i++) {
    g_ref[i] = sqrt((double)x[i] * x[i] + (double)x[h + i] * x[h + i]);
  }
  return worst(M_ELEM, g_got, g_ref, NULL, 0, h);
}

/* the planned spectrum's power (what the bars and ctui-mus's spectrogram
 * read), against |DFT(x * hann * 2/n)|^2 */
static UDSP_SPECTRUM_PLAN g_plan;
static double k_spectrum_power(const float *x, size_t n) {
  udsp_spectrum_plan(&g_plan, (int)n, 48000, 128);
  float *p = udsp_spectrum_power(&g_plan, x, g_got);
  memcpy(g_got2, p, (n / 2 + 1) * sizeof(float));
  hash_floats(g_got2, n / 2 + 1);
  for (size_t i = 0; i < n; i++) {
    g_xd[i] = (double)x[i] * hann_d(i, n) * 2.0 / (double)n;
  }
  dft_d(g_xd, NULL, n, n / 2 + 1, g_ref, g_ref2);
  for (size_t k = 0; k <= n / 2; k++) {
    g_ref[k] = g_ref[k] * g_ref[k] + g_ref2[k] * g_ref2[k];
  }
  return worst(M_MAX, g_got2, g_ref, NULL, max_abs(g_ref, n / 2 + 1),
               n / 2 + 1);
}

static double rms_check(const float *x, size_t n) {
  float got = udsp_analyze_rms(x, n);
  hash_floats(&got, 1);
  double s = 0.0;
  for (size_t i = 0; i < n; i++) {
    s += (double)x[i] * x[i];
  }
  double ref = sqrt(s / (double)n);
  return err_ulp(got, ref, ref);
}

static double k_rms(const float *x, size_t n) { return rms_check(x, n); }

/* the long stretch, once (its signal is music whatever x is) */
static double k_rms_big(const float *x, size_t n) {
  (void)x;
  (void)n;
  return rms_check(g_big, RMS_BIG);
}

static double k_sum(const float *x, size_t n) {
  float got = (float)udsp_analyze_sum(x, n);
  hash_floats(&got, 1);
  double ref = 0.0, abs_sum = 0.0;
  for (size_t i = 0; i < n; i++) {
    ref += x[i];
    abs_sum += fabs((double)x[i]);
  }
  return err_ulp(got, ref, abs_sum);
}

static double k_peak(const float *x, size_t n) {
  float got = udsp_analyze_peak(x, n);
  hash_floats(&got, 1);
  double ref = 0.0;
  for (size_t i = 0; i < n; i++) {
    ref = fabs((double)x[i]) > ref ? fabs((double)x[i]) : ref;
  }
  return worst(M_EXACT, &got, &ref, NULL, 0, 1);
}

static double k_minmax(const float *x, size_t n) {
  float got[2];
  udsp_analyze_minmax(x, n, &got[0], &got[1]);
  hash_floats(got, 2);
  double ref[2] = {x[0], x[0]};
  for (size_t i = 0; i < n; i++) {
    ref[0] = x[i] < ref[0] ? x[i] : ref[0];
    ref[1] = x[i] > ref[1] ? x[i] : ref[1];
  }
  return worst(M_EXACT, got, ref, NULL, 0, 2);
}

/* the signal read as n / 2 stereo frames */
static double k_downmix(const float *x, size_t n) {
  size_t f = n / 2;
  udsp_mono_downmix(x, f, 2, g_got);
  hash_floats(g_got, f);
  for (size_t i = 0; i < f; i++) {
    g_ref[i] = (float)(((double)x[2 * i] + x[2 * i + 1]) * 0.5);
  }
  return worst(M_EXACT, g_got, g_ref, NULL, 0, f);
}

/* BS.1770's K-weighting as two Direct Form I biquads, the textbook form
 * (the kernel's is transposed DF-II, fused) */
static double k_lufs(const float *x, size_t n) {
  float got = udsp_loudness_lufs_momentary(x, n, NULL);
  hash_floats(&got, 1);
  static const double B[2][3] = {
      {1.53512485958697, -2.69169618940638, 1.19839281085285},
      {1.0, -2.0, 1.0}};
  static const double A[2][2] = {{-1.69065929318241, 0.73248077421585},
                                 {-1.99004745483398, 0.99007225036621}};
  double xs[2][2] = {{0}}, ys[2][2] = {{0}}, sumsq = 0.0;
  for (size_t i = 0; i < n; i++) {
    double v = x[i];
    for (int st = 0; st < 2; st++) {
      double y = B[st][0] * v + B[st][1] * xs[st][0] + B[st][2] * xs[st][1] -
                 A[st][0] * ys[st][0] - A[st][1] * ys[st][1];
      xs[st][1] = xs[st][0];
      xs[st][0] = v;
      ys[st][1] = ys[st][0];
      ys[st][0] = y;
      v = y;
    }
    sumsq += v * v;
  }
  double ms = sumsq / (double)n;
  double ref = ms <= 0.0 ? UDSP_LUFS_FLOOR : -0.691 + 10.0 * log10(ms);
  ref = ref < UDSP_LUFS_FLOOR ? UDSP_LUFS_FLOOR : ref;
  return err_ulp(got, ref, ref);
}

static double k_softclip(const float *x, size_t n) {
  const float drive = 3.0f;
  memcpy(g_got, x, n * sizeof(float));
  udsp_softclip(g_got, n, drive);
  hash_floats(g_got, n);
  double inv = 1.0 / (1.0 + (double)drive * 0.25);
  for (size_t i = 0; i < n; i++) {
    double v = (double)x[i] * drive;
    g_ref[i] = v / (1.0 + fabs(v)) * inv * 2.0;
  }
  return worst(M_ELEM, g_got, g_ref, NULL, 0, n);
}

/* dynamics' own log2 and exp2, against libm's in double */
static double k_log2f(const float *x, size_t n) {
  for (size_t i = 0; i < n; i++) {
    float v = fabsf(x[i]) + 0.001f;
    g_got[i] = udsp_log2f(v);
    g_ref[i] = log2((double)v);
  }
  hash_floats(g_got, n);
  return worst(M_ELEM, g_got, g_ref, NULL, 0, n);
}

static double k_exp2f(const float *x, size_t n) {
  for (size_t i = 0; i < n; i++) {
    float v = x[i] * 20.0f;
    v = v > 100.0f ? 100.0f : (v < -100.0f ? -100.0f : v);
    g_got[i] = udsp_exp2f(v);
    g_ref[i] = exp2((double)v);
  }
  hash_floats(g_got, n);
  return worst(M_ELEM, g_got, g_ref, NULL, 0, n);
}

/* the signal as stereo keying a compressor (peak, then RMS) and a gate,
 * and through a limiter: their gains and output, bits only */
static double k_dyn(const float *x, size_t n) {
  for (int rms = 0; rms < 2; rms++) {
    UDSP_COMP c;
    udsp_comp_reset(&c);
    udsp_comp_set(&c, -24.0f, 4.0f, 6.0f, 2.0f, 80.0f, rms, 48000);
    udsp_comp_gain(&c, x, g_got, n / 2);
    hash_floats(g_got, n / 2);
  }
  UDSP_GATE g;
  udsp_gate_set(&g, -30.0f, 4.0f, -40.0f, 0.5f, 5.0f, 20.0f, 48000);
  udsp_gate_reset(&g);
  udsp_gate_gain(&g, x, g_got, n / 2);
  hash_floats(g_got, n / 2);
  static float mem[UDSP_LIMITER_FLOATS(64)];
  UDSP_LIMITER l;
  udsp_limiter_init(&l, mem, 64);
  udsp_limiter_set(&l, 48, 0.25f, 30.0f, 48000);
  memcpy(g_got, x, n * sizeof(float));
  udsp_limiter_stereo(&l, g_got, n / 2);
  hash_floats(g_got, n);
  return 0.0;
}

/* the reverb over the signal as stereo, modulated and diffused, and the
 * saturator's curves at a drive of 4: bits only */
static double k_fdn(const float *x, size_t n) {
  static float *mem;
  if (!mem) {
    mem = calloc(udsp_fdn_floats(48000), sizeof(float));
  }
  UDSP_FDN f;
  udsp_fdn_init(&f, mem, 48000);
  udsp_fdn_set(&f, 0.4f, 1.2f, 5000.0f, 0.7f, 1.0f, 3.0f, 0.8f);
  udsp_fdn_stereo(&f, x, g_got, n / 2);
  hash_floats(g_got, n);
  return 0.0;
}

static double k_waveshape(const float *x, size_t n) {
  for (int k = 0; k < UDSP_SHAPES; k++) {
    memcpy(g_got, x, n * sizeof(float));
    udsp_waveshape(g_got, n, (UDSP_SHAPE)k, 4.0f);
    hash_floats(g_got, n);
  }
  return 0.0;
}

static double k_mix_add(const float *x, size_t n) {
  for (size_t i = 0; i < n; i++) {
    g_got[i] = x[n - 1 - i];
  }
  udsp_mix_add(g_got, x, n, 0.5f);
  hash_floats(g_got, n);
  for (size_t i = 0; i < n; i++) {
    g_ref[i] = (double)x[n - 1 - i] + (double)x[i] * 0.5;
  }
  return worst(M_ELEM, g_got, g_ref, NULL, 0, n);
}

static double k_gain_stereo(const float *x, size_t n) {
  memcpy(g_got, x, n * sizeof(float));
  udsp_gain_stereo(g_got, n / 2, 0.7f, 0.3f);
  hash_floats(g_got, n);
  for (size_t i = 0; i < n; i++) {
    g_ref[i] = (double)x[i] * (double)(i & 1 ? 0.3f : 0.7f);
  }
  return worst(M_ELEM, g_got, g_ref, NULL, 0, n);
}

/* --- the synthesis half's recurrences: bits only --- */

static double k_osc(const float *x, size_t n) {
  (void)x;
  for (int w = 0; w < UDSP_OSC_WAVES; w++) {
    UDSP_OSC o;
    udsp_osc_init(&o, (UDSP_OSC_WAVE)w, 1);
    o.inc = udsp_osc_inc_hz(440.0, 48000);
    udsp_osc_render(&o, g_got, n);
    hash_floats(g_got, n);
    o.inc = udsp_osc_inc_hz(9001.5, 48000);
    udsp_osc_render(&o, g_got, n);
    hash_floats(g_got, n);
  }
  return 0.0;
}

static double k_svf(const float *x, size_t n) {
  for (int mode = 0; mode <= UDSP_SVF_BAND; mode++) {
    UDSP_SVF s;
    udsp_svf_reset(&s);
    udsp_svf_set(&s, 1200.0f, 0.5f, (UDSP_SVF_MODE)mode, 48000);
    memcpy(g_got, x, n * sizeof(float));
    udsp_svf_stereo(&s, g_got, n / 2);
    hash_floats(g_got, n);
  }
  return 0.0;
}

/* every shape of the EQ's band, the signal as stereo, ramping halfway */
static double k_eq(const float *x, size_t n) {
  for (int shape = 0; shape < UDSP_EQ_SHAPES; shape++) {
    UDSP_EQ_BAND b;
    udsp_eq_reset(&b);
    udsp_eq_set(&b, (UDSP_EQ_SHAPE)shape, 1200.0f, 2.0f, 7.5f, 48000);
    memcpy(g_got, x, n * sizeof(float));
    udsp_eq_stereo(&b, g_got, n / 4);
    /* the second half ramping to another setting */
    UDSP_EQ_BAND to;
    udsp_eq_set(&to, (UDSP_EQ_SHAPE)shape, 300.0f, 0.7f, -9.0f, 48000);
    udsp_eq_ramp(&b, &to, (unsigned)(n / 8));
    udsp_eq_stereo(&b, g_got + n / 2, n / 4);
    hash_floats(g_got, n);
  }
  return 0.0;
}

/* the signal as stereo through a 3-band crossover, retuned halfway;
 * through 40 and a half stages of dispersion; keying an up / down
 * compressor (peak, RMS): bits only */
static double k_bands(const float *x, size_t n) {
  static float b0[MAXN], b1[MAXN], b2[MAXN];
  float *out[3] = {b0, b1, b2};
  UDSP_XOVER xo;
  udsp_xover_reset(&xo);
  const float hz[2] = {150.0f, 3000.0f}, hz2[2] = {400.0f, 9000.0f};
  udsp_xover_set(&xo, 3, hz, 48000);
  udsp_xover_split(&xo, x, out, n / 4);
  udsp_xover_set(&xo, 3, hz2, 48000);
  float *half[3] = {b0 + n / 2, b1 + n / 2, b2 + n / 2};
  udsp_xover_split(&xo, x + n / 2, half, n / 4);
  hash_floats(b0, n);
  hash_floats(b1, n);
  hash_floats(b2, n);
  static UDSP_DISPERSE d;
  udsp_disperse_reset(&d);
  udsp_disperse_set(&d, 40.5f, 900.0f, 3.0f, 2.0f, 48000);
  memcpy(g_got, x, n * sizeof(float));
  udsp_disperse_stereo(&d, g_got, n / 2);
  hash_floats(g_got, n);
  for (int rms = 0; rms < 2; rms++) {
    UDSP_UDCOMP c;
    udsp_udcomp_reset(&c);
    udsp_udcomp_set(&c, -20.0f, 8.0f, -40.0f, 4.0f, 24.0f, -80.0f, 2.0f,
                    60.0f, rms, 48000);
    udsp_udcomp_gain(&c, x, g_got, n / 2);
    hash_floats(g_got, n / 2);
  }
  return 0.0;
}

/* every wave's phase moved by the signal (a modulator's swing) */
static double k_osc_pm(const float *x, size_t n) {
  for (int w = 0; w < UDSP_OSC_WAVES; w++) {
    UDSP_OSC o;
    udsp_osc_init(&o, (UDSP_OSC_WAVE)w, 1);
    o.inc = udsp_osc_inc_hz(440.0, 48000);
    udsp_osc_render_pm(&o, x, g_got, n);
    hash_floats(g_got, n);
  }
  return 0.0;
}

/* Operator's (d8): the saw's and the triangle's tables, every recipe's
 * operator under the signal as PM, with and without feedback, the
 * sine's; the envelope through its modes and slopes */
static UDSP_WTAB g_wtab;
static float g_amp[MAXN];
static double k_fmop(const float *x, size_t n) {
  float harm[UDSP_WTAB_HARMONICS];
  for (int s = 0; s <= UDSP_WTAB_TRIANGLE; s++) {
    udsp_wtab_recipe((UDSP_WTAB_SHAPE)s, 64, harm);
    udsp_wtab_build(&g_wtab, harm, 64);
    hash_floats(g_wtab.t[0], UDSP_WTAB_SIZE);
    hash_floats(g_wtab.t[3], UDSP_WTAB_SIZE);
    for (int k = 0; k < 4; k++) {
      UDSP_FMOP o = {k & 1 ? &g_wtab : NULL,
                     0,
                     udsp_osc_inc_hz(k < 2 ? 220.0 : 3100.0, 48000),
                     k & 2 ? 0.45f : 0.0f,
                     0,
                     0};
      for (size_t i = 0; i < n; i++) {
        g_amp[i] = 0.25f + 0.5f * (float)i / (float)n;
      }
      udsp_fmop_render(&o, x, g_amp, g_got, n);
      hash_floats(g_got, n);
    }
  }
  return 0.0;
}

static double k_openv(const float *x, size_t n) {
  (void)x;
  for (int m = 0; m <= UDSP_OPENV_TRIGGER; m++) {
    UDSP_OPENV_SHAPE sh = {0.2f, 1.0f, 0.4f, 0.05f, 2.0f, 6.0f,
                           8.0f, 3.0f, 0.6f, -0.4f, 0.9f, (UDSP_OPENV_MODE)m};
    UDSP_OPENV e;
    udsp_openv_init(&e, &sh, 48000);
    udsp_openv_gate(&e, 1);
    udsp_openv_render(&e, g_got, n / 2);
    udsp_openv_gate(&e, 0);
    udsp_openv_render(&e, g_got + n / 2, n - n / 2);
    hash_floats(g_got, n);
  }
  return 0.0;
}

static double k_zdf(const float *x, size_t n) {
  for (int mode = 0; mode <= UDSP_ZDF_NOTCH; mode++) {
    UDSP_ZDF f;
    udsp_zdf_reset(&f);
    udsp_zdf_set(&f, 1200.0f, 0.8f, (UDSP_ZDF_MODE)mode, 48000);
    memcpy(g_got, x, n * sizeof(float));
    udsp_zdf_mono(&f, g_got, n);
    hash_floats(g_got, n);
  }
  return 0.0;
}

/* an envelope through its stages, and the pitch and LFO math under it */
static double k_adsr(const float *x, size_t n) {
  (void)x;
  UDSP_ADSR e;
  memset(&e, 0, sizeof(e));
  udsp_adsr_set(&e, 3.0f, 10.0f, 0.4f, 15.0f, 48000);
  udsp_adsr_gate(&e, 1);
  udsp_adsr_render(&e, g_got, n / 2);
  udsp_adsr_gate(&e, 0);
  udsp_adsr_render(&e, g_got + n / 2, n - n / 2);
  hash_floats(g_got, n);
  UDSP_LFO l;
  for (int s = 0; s < UDSP_LFO_SHAPES; s++) {
    udsp_lfo_init(&l, (UDSP_LFO_SHAPE)s, 3);
    l.inc = udsp_osc_inc_hz(5.3, 48000);
    for (size_t i = 0; i < n; i++) {
      g_got[i] = udsp_lfo_value(&l);
      udsp_lfo_advance(&l, 32);
    }
    hash_floats(g_got, n);
  }
  for (size_t i = 0; i < n; i++) {
    g_got[i] = (float)udsp_note_inc(20.0 + (double)i * 0.05, 48000);
  }
  hash_floats(g_got, n);
  return 0.0;
}

static float g_ring[UDSP_DELAY_RING_FLOATS(4096)];
static double k_delay(const float *x, size_t n) {
  UDSP_DELAY d;
  udsp_delay_init(&d, g_ring, 4096);
  memcpy(g_got, x, n * sizeof(float));
  udsp_delay_stereo(&d, g_got, n / 2, 1000, 0.4f, 0.35f);
  hash_floats(g_got, n);
  return 0.0;
}

static double k_sampler(const float *x, size_t n) {
  UDSP_SAMPLE smp = {x, (uint32_t)(n / 2), 44100, 2, 0, 0, 60};
  UDSP_SAMPLER_VOICE v;
  memset(&v, 0, sizeof(v));
  udsp_sampler_note_on(&v, &smp, 67, 48000);
  udsp_sampler_render(&v, g_got, n / 2);
  hash_floats(g_got, n);
  return 0.0;
}

/* the cubic resampler at a fractional step, stereo, and YIN's pitch of
 * a tone */
static double k_resample(const float *x, size_t n) {
  uint64_t step = udsp_sampler_step_q32(0.37, 44100, 48000);
  size_t frames = n / 4;
  udsp_resample_cubic(x, 2, 1ull << 32, step, g_got, frames);
  hash_floats(g_got, frames * 2);
  udsp_resample_cubic(x, 1, (1ull << 32) + 777, step, g_got, frames);
  hash_floats(g_got, frames * 2);
  return 0.0;
}

/* the sinc resampler: its table (no libm: the same bits in every libc)
 * and a run across the input's edges, exact and interpolated rows */
static double k_sinc(const float *x, size_t n) {
  static const uint32_t to[2] = {48000, 47999};
  for (int r = 0; r < 2; r++) {
    size_t floats = udsp_sinc_floats(44100, to[r], UDSP_SINC_ZEROS);
    float *table = malloc(floats * sizeof *table);
    UDSP_SINC s;
    udsp_sinc_init(&s, 44100, to[r], UDSP_SINC_ZEROS, table);
    hash_floats(table, floats);
    udsp_sinc_run(&s, x, 0, n, 0, g_got, n);
    hash_floats(g_got, n);
    free(table);
  }
  return 0.0;
}

/* the varispeed sinc: its table and runs across the input's edges at a
 * slower, a faster and a past-the-widest speed (stereo interleaved: the
 * input as two channels) */
static double k_vsinc(const float *x, size_t n) {
  float *table = malloc(udsp_vsinc_floats(UDSP_VSINC_ZEROS) * sizeof *table);
  UDSP_VSINC k;
  udsp_vsinc_init(&k, UDSP_VSINC_ZEROS, table);
  hash_floats(table, udsp_vsinc_floats(UDSP_VSINC_ZEROS));
  static const double speeds[3] = {0.731, 1.618, 5.5};
  for (int s = 0; s < 3; s++) {
    const size_t outs = n / 2;
    udsp_vsinc_run(&k, x, 2, 0, n / 2, -20, 0x9e3779b9u,
                   (uint64_t)(speeds[s] * UDSP_VSINC_ONE), g_got, 2, outs / 2);
    hash_floats(g_got, outs);
  }
  free(table);
  return 0.0;
}

/* YIN in double, the same steps: the reference for the float one (its
 * sums run in vector lanes, so it isn't the same bits across widths:
 * analysis at load time, never in a render) */
static double yin_d(const float *x, size_t n, double rate, double fmin,
                    double fmax, double threshold) {
  size_t tau_max = (size_t)((float)rate / (float)fmin) + 1;
  size_t tau_min = (size_t)((float)rate / (float)fmax);
  size_t w = n - tau_max;
  double running = 0.0;
  g_ref[0] = 1.0;
  for (size_t tau = 1; tau <= tau_max; tau++) {
    double d = 0.0;
    for (size_t j = 0; j < w; j++) {
      double e = (double)x[j] - x[j + tau];
      d += e * e;
    }
    running += d;
    g_ref[tau] = running > 0.0 ? d * (double)tau / running : 1.0;
  }
  size_t best = 0;
  for (size_t tau = tau_min; tau < tau_max; tau++) {
    if (g_ref[tau] < threshold) {
      while (tau + 1 < tau_max && g_ref[tau + 1] < g_ref[tau]) {
        tau++;
      }
      best = tau;
      break;
    }
  }
  if (!best) {
    return 0.0;
  }
  double a = g_ref[best - 1], b = g_ref[best], c = g_ref[best + 1];
  double den = a - 2.0 * b + c;
  double shift = den > 0.0 ? 0.5 * (a - c) / den : 0.0;
  return rate / ((double)best + shift);
}

/* two tones (a sine, a string-like stack of partials) */
static double k_yin(const float *x, size_t n) {
  (void)x;
  double e = 0.0;
  for (int k = 0; k < 2; k++) {
    for (size_t i = 0; i < n; i++) {
      double t = (double)i / 44100.0, v = 0.0;
      for (int h = 1; h <= (k ? 6 : 1); h++) {
        v += sin(2.0 * 3.14159265358979323846 * (k ? 98.0 : 220.0) * h * t + h) / h;
      }
      g_got2[i] = (float)(0.4 * v);
    }
    float hz = udsp_yin(g_got2, n, 44100, 60.0f, 2000.0f, 0.15f, g_got, NULL);
    double ref = yin_d(g_got2, n, 44100.0, 60.0, 2000.0, 0.15);
    hash_floats(&hz, 1);
    double x1 = err_ulp(hz, ref, ref);
    e = x1 > e ? x1 : e;
  }
  return e;
}

static double k_ramp(const float *x, size_t n) {
  UDSP_RAMP r;
  udsp_ramp_init(&r, 300);
  r.target = 0.25f;
  memcpy(g_got, x, n * sizeof(float));
  udsp_ramp_stereo(&r, g_got, n / 2);
  hash_floats(g_got, n);
  return 0.0;
}

/* a program's loudness over the signal as stereo, then as 5.1 (its LFE
 * weighted 0, skipped): the result, the peak and the filter state, bits
 * only */
static UDSP_R128 g_r128;
static double k_r128(const float *x, size_t n) {
  static const int CHANNELS[] = {2, 6};
  for (int t = 0; t < 2; t++) {
    int ch = CHANNELS[t];
    udsp_r128_init(&g_r128, 48000, ch);
    for (int r = 0; r < 8; r++) {
      udsp_r128_add(&g_r128, x, n / (size_t)ch);
    }
    float out[2] = {(float)udsp_r128_integrated(&g_r128.hist),
                    g_r128.hist.peak};
    hash_floats(out, 2);
    for (int c = 0; c < ch; c++) {
      float z[4];
      for (int j = 0; j < 4; j++) {
        z[j] = (float)g_r128.z[c][j];
      }
      hash_floats(z, 4);
    }
    g_hash = (g_hash ^ (uint64_t)g_r128.hops) * 0x100000001b3ull;
  }
  return 0.0;
}

static double k_scope(const float *x, size_t n) {
  udsp_scope_extract(x, n, g_got);
  hash_floats(g_got, UDSP_SCOPE_POINTS);
  return 0.0;
}

/* The tolerances, written down: the arithmetic's own bound where it has a
 * short one (a product of two roundings 1, x*x + y*y 1.5, a sum of
 * squares' square root 1), else a margin over what every build in
 * bench/ulp.py's list measured (bench/README.md "ULPs" has the numbers):
 * a transform's error grows with log n, the power spectrum's doubles it
 * (squared), the old fft's twiddle recurrence drifts (~40 at 2048). A
 * candidate kernel (a sweep's, a rewrite) must stay under them. */
static const KERNEL KERNELS[] = {
    {"hann", M_INPUT, 1.0, k_hann, 2048},
    {"hann_table", M_INPUT, 1.0, k_hann_table, 2048},
    {"fft16", M_MAX, 2.0, k_fft, 16},
    {"fft", M_MAX, 64.0, k_fft, 2048},
    {"fft_table16", M_MAX, 2.0, k_fft_table, 16},
    {"fft_table", M_MAX, 4.0, k_fft_table, 2048},
    {"fft_table8k", M_MAX, 4.0, k_fft_table, 8192},
    {"fft_real16", M_MAX, 2.0, k_fft_real, 16},
    {"fft_real", M_MAX, 4.0, k_fft_real, 2048},
    {"fft_real8k", M_MAX, 4.0, k_fft_real, 8192},
    {"fft_real_windowed", M_MAX, 4.0, k_fft_real_windowed, 2048},
    {"fft_real_inverse", M_MAX, 8.0, k_fft_real_inverse, 2048},
    {"conv", M_BITS, 0.0, k_conv, 2048},
    {"spectrum_power", M_MAX, 4.0, k_spectrum_power, 2048},
    {"spectrum_power8k", M_MAX, 8.0, k_spectrum_power, 8192},
    {"power", M_ELEM, 1.5, k_power, 2048},
    {"magnitude", M_ELEM, 1.5, k_magnitude, 2048},
    {"rms", M_ELEM, 1.0, k_rms, 2048},
    {"rms_big", M_ELEM, 1.0, k_rms_big, RMS_BIG},
    {"sum", M_SUM, 1.0, k_sum, 2048},
    {"peak", M_EXACT, 0.0, k_peak, 2048},
    {"minmax", M_EXACT, 0.0, k_minmax, 2048},
    {"downmix2", M_EXACT, 0.0, k_downmix, 2048},
    {"lufs", M_ELEM, 1.0, k_lufs, 2048},
    {"softclip", M_ELEM, 4.0, k_softclip, 2048},
    {"mix_add", M_ELEM, 1.0, k_mix_add, 2048},
    {"gain_stereo", M_ELEM, 0.5, k_gain_stereo, 2048},
    {"osc", M_BITS, 0.0, k_osc, 2048},
    {"svf", M_BITS, 0.0, k_svf, 2048},
    {"zdf", M_BITS, 0.0, k_zdf, 2048},
    {"eq", M_BITS, 0.0, k_eq, 2048},
    {"dyn", M_BITS, 0.0, k_dyn, 2048},
    {"bands", M_BITS, 0.0, k_bands, 2048},
    {"fdn", M_BITS, 0.0, k_fdn, 2048},
    {"waveshape", M_BITS, 0.0, k_waveshape, 2048},
    {"log2f", M_ELEM, 4.0, k_log2f, 2048},
    {"exp2f", M_ELEM, 4.0, k_exp2f, 2048},
    {"osc_pm", M_BITS, 0.0, k_osc_pm, 2048},
    {"adsr", M_BITS, 0.0, k_adsr, 2048},
    {"fmop", M_BITS, 0.0, k_fmop, 2048},
    {"openv", M_BITS, 0.0, k_openv, 2048},
    {"delay", M_BITS, 0.0, k_delay, 2048},
    {"sampler", M_BITS, 0.0, k_sampler, 2048},
    {"resample", M_BITS, 0.0, k_resample, 2048},
    {"sinc", M_BITS, 0.0, k_sinc, 2048},
    {"vsinc", M_BITS, 0.0, k_vsinc, 2048},
    {"yin", M_ELEM, 4.0, k_yin, 2048},
    {"ramp", M_BITS, 0.0, k_ramp, 2048},
    {"scope", M_BITS, 0.0, k_scope, 2048},
    {"r128", M_BITS, 0.0, k_r128, 8192},
};

int main(int argc, char **argv) {
  const char *only = NULL;
  int verbose = 0;
  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "-k") && i + 1 < argc) {
      only = argv[++i];
    } else if (!strcmp(argv[i], "-v")) {
      verbose = 1;
    } else {
      fprintf(stderr, "usage: %s [-k NAME] [-v]\n", argv[0]);
      return 2;
    }
  }
  corpus();
  int over = 0;
  for (size_t k = 0; k < sizeof(KERNELS) / sizeof(KERNELS[0]); k++) {
    const KERNEL *kn = &KERNELS[k];
    if (only && !strstr(kn->name, only)) {
      continue;
    }
    g_hash = 0xcbf29ce484222325ull;
    double e = 0.0;
    int at = 0;
    for (int s = 0; s < N_SIG; s++) {
      double x = kn->run(g_sig[s], kn->n);
      if (verbose && kn->metric != M_BITS) {
        fprintf(stderr, "  %-18s %-8s %.3f\n", kn->name, SIG_NAME[s], x);
      }
      if (x > e) {
        e = x;
        at = s;
      }
      if (kn->run == k_rms_big) {
        break;
      }
    }
    printf("%s\t%s\t%.3f\t%s\t%.1f\t%016llx\n", kn->name,
           METRIC_NAME[kn->metric], e,
           kn->run == k_rms_big ? "music" : SIG_NAME[at], kn->tol,
           (unsigned long long)g_hash);
    over |= e > kn->tol;
  }
  free(g_big);
  return over;
}
