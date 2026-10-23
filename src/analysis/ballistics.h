#ifndef UDSP_BALLISTICS_H
#define UDSP_BALLISTICS_H

#include <stddef.h>

/* How a displayed value follows a measured one over time: a bar that
 * jumps to every frame's raw value flickers, so a visualizer eases it
 * (attack/release) and marks recent maxima (peak hold). The state
 * (level, peak, age) is caller-owned arrays passed in and updated in
 * place -- the library itself keeps none, so this stays as stateless as
 * spectrum.c: the same inputs always give the same outputs. Values are
 * in whatever unit the caller picked (dBFS bars, linear RMS); only the
 * peak's fall rate is in that unit, per second. dt is the seconds since
 * the previous call, so the result doesn't depend on the frame rate. */

/* level[i] moves toward target[i] as a one-pole filter: with time
 * constant attack_s while rising, release_s while falling (it covers
 * ~63% of the gap in one time constant). A time constant <= 0 jumps
 * straight to the target. */
void udsp_ballistics(float *level, const float *target, size_t n, float dt,
                     float attack_s, float release_s);

/* peak[i] jumps up to level[i] whenever level[i] reaches it (age[i]
 * restarts at 0), stays put for hold_s seconds (age[i] counts them),
 * then falls fall_per_s per second -- never below level[i]. */
void udsp_peak_hold(float *peak, float *age, const float *level, size_t n,
                    float dt, float hold_s, float fall_per_s);

#endif
