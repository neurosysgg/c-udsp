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

/* "sine", "square", ... <-> the wave; UDSP_OSC_WAVES for an unknown name */
UDSP_OSC_WAVE udsp_osc_wave_by_name(const char *name);
const char *udsp_osc_wave_name(UDSP_OSC_WAVE w);

#endif
