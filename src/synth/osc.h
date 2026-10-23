#ifndef UDSP_OSC_H
#define UDSP_OSC_H

#include <stddef.h>
#include <stdint.h>

/* One oscillator (lifted from µOS's utone, osc.c): a 32-bit phase
 * accumulator, so frequency is exact to rate / 2^32 and never drifts over
 * hours, unlike a rotated phasor. Stateful the way ballistics.c is: the
 * state is this caller-owned struct, nothing is kept in the library.
 *
 * - sine: a quarter-wave fold and a 5-term odd polynomial, ~3.5e-6 off
 *   sin() (about -109 dBFS), no libm in the loop;
 * - saw, square: band-limited with a 2-sample PolyBLEP at each edge;
 * - triangle: naive (its aliases fall as 1/n^2), as utone's is -- a
 *   renderer comparing against µOS's exports depends on that;
 * - noise: xorshift32, the frequency ignored.
 *
 * Every wave but noise is a function of the phase alone, so a block is
 * rendered UDSP_W samples at a time (the phases of a vector are
 * phase + k*inc); noise is a serial recurrence and stays scalar. */

typedef enum {
  UDSP_OSC_SINE,
  UDSP_OSC_SQUARE,
  UDSP_OSC_SAW,
  UDSP_OSC_TRIANGLE,
  UDSP_OSC_NOISE,
  UDSP_OSC_WAVES,
} UDSP_OSC_WAVE;

typedef struct {
  UDSP_OSC_WAVE wave;
  uint32_t phase; /* a turn is 2^32 */
  uint32_t inc;   /* phase advance per sample */
  uint32_t rng;   /* noise state, never 0 */
} UDSP_OSC;

/* phase 0, inc 0 (silent but for noise); seed picks the noise sequence */
void udsp_osc_init(UDSP_OSC *o, UDSP_OSC_WAVE wave, uint32_t seed);

/* the phase increment for a frequency, exact in integers from
 * millihertz (440000 = 440 Hz) or rounded from hertz */
uint32_t udsp_osc_inc_millihz(uint32_t millihz, uint32_t rate);
uint32_t udsp_osc_inc_hz(double hz, uint32_t rate);

/* writes n mono samples in [-1, 1] (the PolyBLEP edges overshoot a
 * little) and advances the phase by n * inc */
void udsp_osc_render(UDSP_OSC *o, float *out, size_t n);

/* the same, each sample's phase moved by pm[i] turns (any sign, any
 * size): phase modulation, the way FM synths do FM -- pm is a
 * modulator's output times its index. The PolyBLEP edges assume the
 * unmodulated rate, so a modulated saw or square aliases somewhat more;
 * noise ignores pm. pm all 0 gives udsp_osc_render's bits. */
void udsp_osc_render_pm(UDSP_OSC *o, const float *pm, float *out, size_t n);

/* the oscillator's sine at a phase in turns, 0 <= turns < 1 (an LFO's,
 * a filter's coefficient; inline: an FM operator's feedback waits on it
 * every sample) */
static inline float udsp_osc_sine(float turns) {
  /* [0,1) -> [-0.5,0.5) (odd about the half turn), then folded to
   * |x| <= 0.25 (symmetric about the quarter turn): five terms suffice */
  float x = turns >= 0.5f ? turns - 1.0f : turns;
  if (x > 0.25f) {
    x = 0.5f - x;
  } else if (x < -0.25f) {
    x = -0.5f - x;
  }
  const float u = 6.283185307179586f * x;
  const float u2 = u * u;
  float p = -1.0f / 362880.0f;
  p = p * u2 + (1.0f / 5040.0f);
  p = p * u2 - (1.0f / 120.0f);
  p = p * u2 + (1.0f / 6.0f);
  p = p * u2 - 1.0f;
  return -u * p;
}

/* udsp_osc_sine() of n turns, each in [0, 1) */
void udsp_osc_sines(const float *turns, float *out, size_t n);

/* "sine", "square", ... <-> the wave; UDSP_OSC_WAVES for an unknown name */
UDSP_OSC_WAVE udsp_osc_wave_by_name(const char *name);
const char *udsp_osc_wave_name(UDSP_OSC_WAVE w);

#endif
