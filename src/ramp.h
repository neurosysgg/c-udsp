#ifndef UDSP_RAMP_H
#define UDSP_RAMP_H

#include <stddef.h>

/* A linear ramp toward a target, one fixed step per frame (µOS's
 * utracker's note envelope: a 3 ms attack and release that remove the
 * click a note starting at full level makes). Applied per frame, not per
 * block: a per-block ramp is a staircase, which zippers. The target can
 * change at any time; the ramp goes on from where it is. */

typedef struct {
  float value;
  float target;
  float step; /* > 0 */
} UDSP_RAMP;

/* value and target 0, step 1/frames */
void udsp_ramp_init(UDSP_RAMP *r, unsigned frames);

/* multiplies frames interleaved stereo frames by the ramp, advancing it */
void udsp_ramp_stereo(UDSP_RAMP *r, float *io, size_t frames);

/* settled at 0: the output is silence whatever the source */
static inline int udsp_ramp_silent(const UDSP_RAMP *r) {
  return r->value <= 0.0f && r->target <= 0.0f;
}

#endif
