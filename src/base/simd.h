#ifndef UDSP_SIMD_H
#define UDSP_SIMD_H

/* internal: the vector type every hot loop here is written against
 * (bench/README.md "Reductions" has the measurements behind it).
 *
 * GCC vector extensions (gcc and clang) at the target's native width: 16
 * bytes on baseline x86-64/SSE2 or NEON, 32 with AVX. Explicit vectors
 * because both compilers vectorise these at -O2 already and the same way
 * every time, where auto-vectorisation of the same loops comes and goes
 * with the optimisation level. Never wider than the target: a 32-byte
 * vector on an SSE2 target is split and, for anything but + and *,
 * scalarised (6x slower than a plain loop in the study).
 *
 * UDSP_SCALAR (or a compiler without the extension) makes a vector one
 * float: the same kernels as plain C, the reference the tests hold the
 * vector build to. */

#include <math.h>
#include <string.h>

#if defined(__GNUC__) && !defined(UDSP_SCALAR)
#if defined(__AVX__)
#define UDSP_W 8
#else
#define UDSP_W 4
#endif
typedef float UDSP_VF __attribute__((vector_size(UDSP_W * 4)));
typedef int UDSP_VI __attribute__((vector_size(UDSP_W * 4)));
/* phase accumulators (osc.c): uint32 lanes that wrap like the scalar
 * one, converted with __builtin_convertvector (correctly rounded, as a
 * scalar (float) cast is) */
typedef unsigned UDSP_VU __attribute__((vector_size(UDSP_W * 4)));
/* a where m is set (all ones), else b: the blend vector ?: would be */
static inline UDSP_VF udsp_vsel(UDSP_VI m, UDSP_VF a, UDSP_VF b) {
  return (UDSP_VF)((m & (UDSP_VI)a) | (~m & (UDSP_VI)b));
}
/* four lanes whatever the target: for loops shorter than UDSP_W (an FFT's
 * narrow stages under AVX) */
#define UDSP_HAVE_V4 1
typedef float UDSP_V4 __attribute__((vector_size(16)));

/* the loads and stores below go through these: unaligned, and allowed to
 * alias the floats. Not memcpy: gcc -O2 splits a memcpy into a vector
 * that's then carried round a loop into two 8-byte loads and keeps the
 * vector on the stack (conv_direct's four sums: a store and a load a
 * tap, 1.2x slower than one sum instead of 2x faster) */
typedef float UDSP_V4U __attribute__((vector_size(16), aligned(4), may_alias));

static inline UDSP_V4 udsp_v4load(const float *p) {
  return *(const UDSP_V4U *)p;
}
static inline void udsp_v4store(float *p, UDSP_V4 v) { *(UDSP_V4U *)p = v; }
/* lanes picked from a (0-3) and b (4-7): __builtin_shufflevector (clang,
 * gcc 12+), else gcc's own __builtin_shuffle */
typedef int UDSP_V4I __attribute__((vector_size(16)));
/* udsp_vsel's four-lane form */
static inline UDSP_V4 udsp_v4sel(UDSP_V4I m, UDSP_V4 a, UDSP_V4 b) {
  return (UDSP_V4)((m & (UDSP_V4I)a) | (~m & (UDSP_V4I)b));
}
#if defined(__clang__) || __GNUC__ >= 12
#define UDSP_SHUF4(a, b, i, j, k, l) __builtin_shufflevector(a, b, i, j, k, l)
#else
#define UDSP_SHUF4(a, b, i, j, k, l)                                           \
  __builtin_shuffle(a, b, (UDSP_V4I){i, j, k, l})
#endif

typedef float UDSP_VFU
    __attribute__((vector_size(UDSP_W * 4), aligned(4), may_alias));

static inline UDSP_VF udsp_vload(const float *p) {
  return *(const UDSP_VFU *)p;
}
static inline void udsp_vstore(float *p, UDSP_VF v) { *(UDSP_VFU *)p = v; }
static inline UDSP_VF udsp_vsplat(float x) { return (UDSP_VF){0} + x; }
/* a > b ? a : b per lane, and |a|: the target's own instructions where it
 * has them. The generic form (a compare and a blend: vector ?: is C++-only
 * in gcc) is right everywhere, but gcc doesn't fold it back into maxps and
 * runs it in the integer domain: 3x slower in the peak kernel. maxps is
 * exactly a > b ? a : b, so every build gives the scalar build's answer */
#if defined(__AVX__) || defined(__SSE__)
#include <immintrin.h>
#endif
static inline UDSP_VF udsp_vmax(UDSP_VF a, UDSP_VF b) {
#if defined(__AVX__)
  return (UDSP_VF)_mm256_max_ps((__m256)a, (__m256)b);
#elif defined(__SSE__)
  return (UDSP_VF)_mm_max_ps((__m128)a, (__m128)b);
#else
  UDSP_VI m = a > b;
  return (UDSP_VF)((m & (UDSP_VI)a) | (~m & (UDSP_VI)b));
#endif
}
static inline UDSP_VF udsp_vmin(UDSP_VF a, UDSP_VF b) {
#if defined(__AVX__)
  return (UDSP_VF)_mm256_min_ps((__m256)a, (__m256)b);
#elif defined(__SSE__)
  return (UDSP_VF)_mm_min_ps((__m128)a, (__m128)b);
#else
  UDSP_VI m = a < b;
  return (UDSP_VF)((m & (UDSP_VI)a) | (~m & (UDSP_VI)b));
#endif
}
static inline UDSP_VF udsp_vabs(UDSP_VF a) {
#if defined(__AVX__)
  return (UDSP_VF)_mm256_andnot_ps(_mm256_set1_ps(-0.0f), (__m256)a);
#elif defined(__SSE__)
  return (UDSP_VF)_mm_andnot_ps(_mm_set1_ps(-0.0f), (__m128)a);
#else
  return (UDSP_VF)((UDSP_VI)a & ((UDSP_VI){0} + 0x7fffffff));
#endif
}
/* sqrt per lane: sqrtps (a scalar sqrtf() loop doesn't vectorise without
 * -fno-math-errno: errno for negative input) */
static inline UDSP_VF udsp_vsqrt(UDSP_VF v) {
#if defined(__AVX__)
  return (UDSP_VF)_mm256_sqrt_ps((__m256)v);
#elif defined(__SSE__)
  return (UDSP_VF)_mm_sqrt_ps((__m128)v);
#else
  for (int j = 0; j < UDSP_W; j++) {
    v[j] = sqrtf(v[j]);
  }
  return v;
#endif
}

/* the lanes in reverse order */
static inline UDSP_VF udsp_vrev(UDSP_VF v) {
#if defined(__AVX__)
  __m256 r = _mm256_permute_ps((__m256)v, 0x1b); /* 3 2 1 0 in each half */
  return (UDSP_VF)_mm256_permute2f128_ps(r, r, 0x01);
#elif defined(__SSE__)
  return (UDSP_VF)_mm_shuffle_ps((__m128)v, (__m128)v, 0x1b);
#else
  UDSP_VF r;
  for (int j = 0; j < UDSP_W; j++) {
    r[j] = v[UDSP_W - 1 - j];
  }
  return r;
#endif
}

/* UDSP_W interleaved stereo frames from p into their left and right
 * channels' lanes */
static inline void udsp_vlr(const float *p, UDSP_VF *l, UDSP_VF *r) {
  const UDSP_VF a = udsp_vload(p), b = udsp_vload(p + UDSP_W);
#if defined(__clang__) || __GNUC__ >= 12
#if UDSP_W == 8
  *l = __builtin_shufflevector(a, b, 0, 2, 4, 6, 8, 10, 12, 14);
  *r = __builtin_shufflevector(a, b, 1, 3, 5, 7, 9, 11, 13, 15);
#else
  *l = __builtin_shufflevector(a, b, 0, 2, 4, 6);
  *r = __builtin_shufflevector(a, b, 1, 3, 5, 7);
#endif
#else
#if UDSP_W == 8
  *l = __builtin_shuffle(a, b, (UDSP_VI){0, 2, 4, 6, 8, 10, 12, 14});
  *r = __builtin_shuffle(a, b, (UDSP_VI){1, 3, 5, 7, 9, 11, 13, 15});
#else
  *l = __builtin_shuffle(a, b, (UDSP_VI){0, 2, 4, 6});
  *r = __builtin_shuffle(a, b, (UDSP_VI){1, 3, 5, 7});
#endif
#endif
}

/* the lanes moved up by one, x into lane 0: a wavefront's lane taking
 * the output its neighbour made a step before */
static inline UDSP_VF udsp_vshift_in(UDSP_VF v, float x) {
  const UDSP_VF s = udsp_vsplat(x);
#if defined(__clang__) || __GNUC__ >= 12
#if UDSP_W == 8
  return __builtin_shufflevector(v, s, 8, 0, 1, 2, 3, 4, 5, 6);
#else
  return __builtin_shufflevector(v, s, 4, 0, 1, 2);
#endif
#else
#if UDSP_W == 8
  return __builtin_shuffle(v, s, (UDSP_VI){8, 0, 1, 2, 3, 4, 5, 6});
#else
  return __builtin_shuffle(v, s, (UDSP_VI){4, 0, 1, 2});
#endif
#endif
}

static inline float udsp_vsum(UDSP_VF v) {
  float s = 0.0f;
  for (int j = 0; j < UDSP_W; j++) {
    s += v[j];
  }
  return s;
}
static inline float udsp_vhmax(UDSP_VF v) {
  float m = v[0];
  for (int j = 1; j < UDSP_W; j++) {
    m = v[j] > m ? v[j] : m;
  }
  return m;
}
static inline float udsp_vhmin(UDSP_VF v) {
  float m = v[0];
  for (int j = 1; j < UDSP_W; j++) {
    m = v[j] < m ? v[j] : m;
  }
  return m;
}
#else
#define UDSP_W 1
typedef float UDSP_VF;

static inline UDSP_VF udsp_vload(const float *p) { return *p; }
static inline void udsp_vstore(float *p, UDSP_VF v) { *p = v; }
static inline UDSP_VF udsp_vsplat(float x) { return x; }
static inline UDSP_VF udsp_vmax(UDSP_VF a, UDSP_VF b) { return a > b ? a : b; }
static inline UDSP_VF udsp_vmin(UDSP_VF a, UDSP_VF b) { return a < b ? a : b; }
static inline UDSP_VF udsp_vabs(UDSP_VF a) { return fabsf(a); }
static inline UDSP_VF udsp_vrev(UDSP_VF v) { return v; }
static inline UDSP_VF udsp_vsqrt(UDSP_VF v) { return sqrtf(v); }
static inline float udsp_vsum(UDSP_VF v) { return v; }
static inline float udsp_vhmax(UDSP_VF v) { return v; }
static inline float udsp_vhmin(UDSP_VF v) { return v; }
#endif

/* independent accumulators a reduction keeps: FP add/FMA latency (3-4
 * cycles on Zen 3) x the two pipes that can take one each cycle. Fewer
 * leaves the loop waiting on its own chain -- the one thing unrolling
 * (-funroll-loops) doesn't fix, since without -ffast-math the compiler
 * may not split a float sum itself. Overridable (-DUDSP_ACC=N) for
 * bench/sweep.py, which measures the choice per kernel and build. */
#ifndef UDSP_ACC
#define UDSP_ACC 8
#endif

/* on each loop over the accumulators: gcc's -O2 doesn't fully unroll
 * them by itself, and then the array stays in memory instead of
 * registers (a store-to-load trip per step, 2-3x slower). Written as a
 * _Pragma so a macro can carry it. UDSP_NO_UNROLL drops it (bench/
 * sweep.py's knob; never for a real build) */
#ifdef UDSP_NO_UNROLL
#define UDSP_UNROLL
#else
#define UDSP_UNROLL _Pragma("GCC unroll 16")
#endif

/* on a function running a stereo pair's two serial chains side by side:
 * gcc's SLP vectoriser packs the pair into 2-lane vectors and clears
 * their upper lanes after every operation (a movq on the chain), which
 * doubles the chain's latency (the EQ, the crossover: 1.6-1.8x slower
 * than as two scalar chains). The same operations per lane either way */
#if defined(__GNUC__) && !defined(__clang__)
#define UDSP_NO_SLP __attribute__((optimize("no-tree-slp-vectorize")))
#else
#define UDSP_NO_SLP
#endif

#if UDSP_W > 1
/* base[idx] per lane */
static inline UDSP_VF udsp_vgather(const float *base, UDSP_VI idx) {
#if defined(__AVX2__)
  return (UDSP_VF)_mm256_i32gather_ps(base, (__m256i)idx, 4);
#else
  UDSP_VF v;
  UDSP_UNROLL
  for (int j = 0; j < UDSP_W; j++) {
    v[j] = base[idx[j]];
  }
  return v;
#endif
}

/* fx/line.h's udsp_hermite per lane: the same operations in its order */
static inline UDSP_VF udsp_vhermite(UDSP_VF ym1, UDSP_VF y0, UDSP_VF y1,
                                    UDSP_VF y2, UDSP_VF x) {
  const UDSP_VF c1 = 0.5f * (y1 - ym1);
  const UDSP_VF c2 = ym1 - 2.5f * y0 + 2.0f * y1 - 0.5f * y2;
  const UDSP_VF c3 = 0.5f * (y2 - ym1) + 1.5f * (y0 - y1);
  return ((c3 * x + c2) * x + c1) * x + y0;
}

/* four floats side by side from each lane's base + at[lane], transposed:
 * *c0 holds each lane's first, *c3 its last. A load a lane and a
 * transpose where there are vector loads, else read one by one */
static inline void udsp_vcols4(const float *base, UDSP_VI at, UDSP_VF *c0,
                               UDSP_VF *c1, UDSP_VF *c2, UDSP_VF *c3) {
#if defined(__AVX__)
  __m256 r[4];
  for (int k = 0; k < 4; k++) {
    r[k] = _mm256_set_m128(_mm_loadu_ps(base + at[k + 4]),
                           _mm_loadu_ps(base + at[k]));
  }
  const __m256 t0 = _mm256_unpacklo_ps(r[0], r[1]);
  const __m256 t1 = _mm256_unpackhi_ps(r[0], r[1]);
  const __m256 t2 = _mm256_unpacklo_ps(r[2], r[3]);
  const __m256 t3 = _mm256_unpackhi_ps(r[2], r[3]);
  *c0 = (UDSP_VF)_mm256_shuffle_ps(t0, t2, 0x44);
  *c1 = (UDSP_VF)_mm256_shuffle_ps(t0, t2, 0xee);
  *c2 = (UDSP_VF)_mm256_shuffle_ps(t1, t3, 0x44);
  *c3 = (UDSP_VF)_mm256_shuffle_ps(t1, t3, 0xee);
#elif defined(__SSE__)
  __m128 r0 = _mm_loadu_ps(base + at[0]), r1 = _mm_loadu_ps(base + at[1]);
  __m128 r2 = _mm_loadu_ps(base + at[2]), r3 = _mm_loadu_ps(base + at[3]);
  _MM_TRANSPOSE4_PS(r0, r1, r2, r3);
  *c0 = (UDSP_VF)r0;
  *c1 = (UDSP_VF)r1;
  *c2 = (UDSP_VF)r2;
  *c3 = (UDSP_VF)r3;
#else
  UDSP_UNROLL
  for (int j = 0; j < UDSP_W; j++) {
    (*c0)[j] = base[at[j]];
    (*c1)[j] = base[at[j] + 1];
    (*c2)[j] = base[at[j] + 2];
    (*c3)[j] = base[at[j] + 3];
  }
#endif
}

/* a delay line's four taps per lane (fx/line.h's udsp_hermite's points):
 * the ring at base + off, mask its length - 1, at (masked) the newest
 * tap's neighbour -- at + 1, at, at - 1, at - 2 back as ym1 and *y0..*y2.
 * They lie side by side in the ring: udsp_vcols4 where none of them
 * wraps, else gathered */
static inline UDSP_VF udsp_vtaps(const float *base, UDSP_VI off, UDSP_VI at,
                                 UDSP_VI mask, UDSP_VF *y0, UDSP_VF *y1,
                                 UDSP_VF *y2) {
#if defined(__SSE__)
  const UDSP_VI edge = (at < 2) | (at > mask - 1);
#if defined(__AVX__)
  const int wraps = _mm256_movemask_ps((__m256)edge);
#else
  const int wraps = _mm_movemask_ps((__m128)edge);
#endif
  if (!wraps) {
    UDSP_VF ym1;
    udsp_vcols4(base, off + at - 2, y2, y1, y0, &ym1);
    return ym1;
  }
#endif
  *y0 = udsp_vgather(base, off + at);
  *y1 = udsp_vgather(base, off + ((at - 1) & mask));
  *y2 = udsp_vgather(base, off + ((at - 2) & mask));
  return udsp_vgather(base, off + ((at + 1) & mask));
}
#endif

#endif
