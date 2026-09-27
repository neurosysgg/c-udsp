#ifndef UDSP_GAIN_H
#define UDSP_GAIN_H

#include <stddef.h>

/* multiplies every one of the n samples in buf by gain, in place --
 * channel-agnostic (a flat float stream, same "doesn't need to know
 * frames vs. channels" convention udsp_analyze_downmix() and every other
 * function here already follows: the caller decides what n spans).
 * gain is not clamped here -- e.g. a caller building a pure-attenuation
 * volume control validates/clamps its own range (see ctui-mus's
 * mus_player_cmd_volume()), this is just the multiply. */
void udsp_gain_apply(float *buf, size_t n, float gain);

#endif
