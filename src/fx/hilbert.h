#ifndef UDSP_HILBERT_H
#define UDSP_HILBERT_H

/* An analytic signal: Olli Niemitalo's pair of allpass chains, four
 * second-order sections each (y = a^2 (x + y[-2]) - x[-2]), 90 degrees
 * apart within 0.7 from 0.0005 to 0.4995 of the rate -- the first
 * chain's output a frame late is the real part, the second's the
 * imaginary, leading it. re cos + im sin of a turning phase shifts every
 * partial up by its pace (one sideband: a frequency shifter). Both are
 * allpasses: the level kept, the phase turned. */

typedef struct {
  float x1, x2, y1, y2;
} UDSP_HILBERT_SECTION;

typedef struct {
  UDSP_HILBERT_SECTION s[2][4];
  float late; /* the first chain's last output */
} UDSP_HILBERT;

/* a frame in, its real and imaginary parts out */
void udsp_hilbert_step(UDSP_HILBERT *h, float x, float *re, float *im);

#endif
