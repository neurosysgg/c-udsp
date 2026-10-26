#ifndef UDSP_OS_H
#define UDSP_OS_H

#include <stddef.h>

/* Oversampling by 2 or 4, stereo: each doubling a polyphase IIR halfband
 * (Laurent de Soras' hiir: two chains of first-order allpasses in z^-2,
 * the filter half their sum), its coefficients designed offline with
 * hiir's elliptic recipe and fixed. At 48 kHz: flat to 21 kHz (1e-8 dB);
 * the first doubling (8 allpasses) -92 dB from 27 kHz, the second (4)
 * -82 dB from 75 kHz -- an image or an alias that would land under
 * 21 kHz that far down. Up then down is an allpass at the base rate:
 * 4.6 frames late at 4x near DC, more toward the top (12.8 at 20 kHz).
 * A frame at a time through both doublings: no scratch, any block. */

#define UDSP_OS_STAGES 2
#define UDSP_OS_COEFS 8 /* the most allpasses in a doubling */

typedef struct {
  float x[UDSP_OS_COEFS][2], y[UDSP_OS_COEFS][2]; /* an allpass's, a channel */
} UDSP_OS_CHAIN;

typedef struct {
  UDSP_OS_CHAIN up[UDSP_OS_STAGES], down[UDSP_OS_STAGES];
  float odd[UDSP_OS_STAGES][2]; /* down: the odd frame before */
} UDSP_OS;

void udsp_os_reset(UDSP_OS *o);

/* frames interleaved stereo frames into factor * frames of them (factor 2
 * or 4) */
void udsp_os_up(UDSP_OS *o, const float *in, float *out, size_t frames,
                unsigned factor);

/* factor * frames interleaved stereo frames back into frames */
void udsp_os_down(UDSP_OS *o, const float *in, float *out, size_t frames,
                  unsigned factor);

#endif
