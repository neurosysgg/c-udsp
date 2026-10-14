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

/* frames interleaved stereo frames: left times gl, right times gr */
void udsp_gain_stereo(float *buf, size_t frames, float gl, float gr);

/* dst[i] += src[i] * gain over n samples: a mixer's sum. The product is
 * rounded before the add (two roundings, never a fused multiply-add
 * unless the compiler is told to contract), so a sum built this way
 * matches one built as "scale, then add" */
void udsp_mix_add(float *dst, const float *src, size_t n, float gain);

/* how a pan position (-1 left .. 0 centre .. 1 right) becomes two gains:
 * - BALANCE: the centre is unity on both sides, a side is attenuated on
 *   an equal-power curve as the pan moves away from it (a track
 *   mixer's pan: a stereo source keeps its level in the middle);
 * - EQUAL_POWER: sin / cos over the whole range, -3 dB each in the
 *   centre (a mono source panned keeps its loudness) */
typedef enum {
  UDSP_PAN_BALANCE,
  UDSP_PAN_EQUAL_POWER,
} UDSP_PAN_LAW;

void udsp_pan_gains(float pan, UDSP_PAN_LAW law, float *gl, float *gr);

#endif
