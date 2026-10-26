#ifndef UDSP_SOS_H
#define UDSP_SOS_H

#include <stddef.h>

/* A fixed cascade of biquads, stereo -- a filter fitted offline (a
 * plugin's, measured), its sections the caller's constants: direct form
 * II transposed, a[0] = 1. */

#define UDSP_SOS_MAX 8

typedef struct {
  float b0, b1, b2, a1, a2;
} UDSP_SOS_SECTION;

typedef struct {
  UDSP_SOS_SECTION sec[UDSP_SOS_MAX];
  int n;
  float s[UDSP_SOS_MAX][2][2]; /* a section's two states, a channel */
} UDSP_SOS;

/* n (<= UDSP_SOS_MAX) sections copied in, the state zeroed */
void udsp_sos_init(UDSP_SOS *f, const UDSP_SOS_SECTION *sec, int n);

void udsp_sos_reset(UDSP_SOS *f);

/* frames interleaved stereo frames in place */
void udsp_sos_stereo(UDSP_SOS *f, float *io, size_t frames);

#endif
