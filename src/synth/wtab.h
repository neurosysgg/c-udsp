#ifndef UDSP_WTAB_H
#define UDSP_WTAB_H

#include <stddef.h>
#include <stdint.h>

/* Band-limited wavetables and an FM operator reading them (ctui-daw's
 * FM synth, d8: Ableton's Operator). A wave is up to 64 harmonics in
 * sine phase (a saw's, a square's, a user's drawn partials); built into
 * one table per octave, each holding only the harmonics a fundamental of
 * that octave keeps below Nyquist (64, 32, ... 1), so a wave played high
 * drops its top partials instead of folding them back. The level is
 * picked from the unmodulated increment: phase modulation pushes the
 * instantaneous frequency past it, and what that adds can still alias
 * (Operator's own simplification). Tables are read linearly between 4096
 * points (harmonic 64 of a table: ~-70 dB off its own level).
 *
 * The harmonics are summed from a table of the oscillator's polynomial
 * sine, never libm's, in a fixed order: the same table on every build.
 * The caller owns it (~115 KB): built once per wave, read-only after. */

#define UDSP_WTAB_BITS 12
#define UDSP_WTAB_SIZE (1u << UDSP_WTAB_BITS)
#define UDSP_WTAB_HARMONICS 64
#define UDSP_WTAB_LEVELS 7 /* 64, 32, 16, 8, 4, 2, 1 harmonics */

typedef struct {
  /* each level's period, its first point again at the end */
  float t[UDSP_WTAB_LEVELS][UDSP_WTAB_SIZE + 1];
} UDSP_WTAB;

typedef enum {
  UDSP_WTAB_SAW,
  UDSP_WTAB_SQUARE,
  UDSP_WTAB_TRIANGLE,
} UDSP_WTAB_SHAPE;

/* harmonic h's amplitude is amp[h - 1], n of them (64 at most); every
 * level scaled by one factor, the full table's peak to 1, so the
 * fundamental is as loud at every level (all zero: silence) */
void udsp_wtab_build(UDSP_WTAB *w, const float *amp, size_t n);

/* a shape's first n harmonics (its Fourier series in sines) into
 * amp[UDSP_WTAB_HARMONICS], zeros past them: Operator's Saw3 .. Saw64 are
 * n = 3 .. 64 (a square's and a triangle's even ones are zero) */
void udsp_wtab_recipe(UDSP_WTAB_SHAPE s, size_t n, float *amp);

/* the level a fundamental advancing inc a sample reads: the most
 * harmonics that stay below Nyquist (the last level past it) */
int udsp_wtab_level(uint32_t inc);

/* level's table at phase t in turns, 0 <= t < 1 */
static inline float udsp_wtab_at(const UDSP_WTAB *w, int level, float t) {
  const float x = t * (float)UDSP_WTAB_SIZE;
  const int i = (int)x;
  const float *p = w->t[level] + i;
  return p[0] + (p[1] - p[0]) * (x - (float)i);
}

/* An FM operator: a phase accumulator reading a table (NULL: the
 * oscillator's sine, exact), its phase moved by a modulator's output in
 * turns and by its own feedback -- the average of its last two outputs
 * (the DX7's, which keeps high feedback from oscillating at Nyquist)
 * times feedback turns. */
typedef struct {
  const UDSP_WTAB *wave;
  uint32_t phase, inc;
  float feedback; /* turns a unit of its output moves its phase */
  float y1, y2;   /* its last two outputs */
} UDSP_FMOP;

/* n samples of amp[i] * wave(phase + pm[i] + feedback): pm NULL is none,
 * amp the operator's level a sample (its envelope times its gain). With
 * feedback each sample waits on the one before (scalar); without it a
 * table's samples don't, and a sine is udsp_osc_render_pm's vectors
 * times amp. The same bits whatever the blocks. */
void udsp_fmop_render(UDSP_FMOP *o, const float *pm, const float *amp,
                      float *out, size_t n);

#endif
