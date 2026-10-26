#ifndef UDSP_DYN_H
#define UDSP_DYN_H

#include <stddef.h>
#include <stdint.h>

/* Dynamics: a compressor's and a gate's gain from a key signal, and a
 * lookahead limiter. The level work is in dB, per frame, through a
 * float log2 / exp2 of their own (no libm: the same bits everywhere, as
 * the rest of the synthesis half), both stereo-linked (the louder
 * channel decides, so the image doesn't move). The compressor and the
 * gate only work out a gain from the key: the caller applies it to
 * whatever it likes (the key itself, or the main signal when the key is
 * a sidechain; delayed, for lookahead), so one detector serves both. */

/* log2(x) for x > 0: the exponent's bits and an odd series in (m - 1) /
 * (m + 1), m in [0.71, 1.41); 2^x for -126 < x < 128 -- both within 3 ULPs
 * (bench/ulp.c) */
float udsp_log2f(float x);
float udsp_exp2f(float x);

#define UDSP_DB_PER_LOG2 6.0205999f /* 20 log10(2) */

static inline float udsp_db_from_gain(float g) {
  return UDSP_DB_PER_LOG2 * udsp_log2f(g < 1e-15f ? 1e-15f : g);
}
static inline float udsp_gain_from_db(float db) {
  return udsp_exp2f(db * (1.0f / UDSP_DB_PER_LOG2));
}

/* --- the compressor ----------------------------------------------------
 * Giannoulis, Massberg and Reiss's design (JAES 2012): a gain computer
 * with a soft knee over the key's level, its reduction smoothed in dB by
 * a branching one-pole (attack while it grows, release while it falls).
 * Peak (the louder channel's |x|) or RMS (its mean square over 10 ms).
 */

typedef struct {
  float thr, slope, knee; /* dB; 1 / ratio - 1 (0: no compression) */
  float att, rel;         /* the smoother's poles */
  float ms_coef;          /* RMS's pole */
  int rms;
  float red; /* the smoothed reduction, dB (>= 0) */
  float ms;  /* RMS's mean square */
  /* appended (iZotope Ozone 11's Dynamics, measured): a second knee --
   * its slope's change past its threshold, knee wide, smoothed apart and
   * let go through a tail -- a one-pole after the first's smoother (post
   * 0: none), and the peak's window and hold */
  int second;
  float thr2, slope2, knee2;
  float att2, rel2, tail; /* the second's poles */
  float r2, t2;           /* its smoothed reduction, after the tail */
  float post, red2;
  float hold, held;  /* the peak's hold (its release's pole; 0: none) */
  uint32_t block, n; /* the window's block (frames; 0: none), its fill */
  float cur, prev;   /* this block's peak, the last one's */
} UDSP_COMP;

void udsp_comp_reset(UDSP_COMP *c);

/* threshold and knee in dB, ratio >= 1 (1e9 or more: a limiter's
 * infinity), times in ms (> 0); the second knee and the smoothing off
 * (set them after) */
void udsp_comp_set(UDSP_COMP *c, float threshold_db, float ratio, float knee_db,
                   float attack_ms, float release_ms, int rms, unsigned rate);

/* a second knee, after udsp_comp_set: past threshold_db the level's slope
 * goes from the compressor's 1 / ratio to 1 / ratio2 -- its term (1 /
 * ratio2 - 1 / ratio) (x - threshold) whether the first acts there or not
 * (Ozone 11's "Lim" stage: under the first threshold, or with ratio2
 * under ratio, an expansion); its own quadratic knee, its own attack and
 * release, then a tail: taken at once, let go over tail_ms (0: none).
 * ratio2 <= 0: none */
void udsp_comp_set_second(UDSP_COMP *c, float threshold_db, float ratio2,
                          float knee_db, float attack_ms, float release_ms,
                          float tail_ms, unsigned rate);

/* the smoothed reduction through one more one-pole of ms (0: none) */
void udsp_comp_set_smooth(UDSP_COMP *c, float ms, unsigned rate);

/* the peak held before the curve -- taken at once, let go over ms (a
 * linear one-pole; 0: none) -- so a waveform's troughs don't let the
 * reduction sag (Ozone 11's Dynamics' detector: ripple-free on a sine);
 * the peak detector only */
void udsp_comp_set_hold(UDSP_COMP *c, float ms, unsigned rate);

/* the peak over the last ms to 2 ms (the larger of this block's and the
 * last one's, blocks ms long; 0: none), before the hold -- a window's
 * maximum without a buffer (Ozone 11's: ~5 ms); the peak detector only */
void udsp_comp_set_window(UDSP_COMP *c, float ms, unsigned rate);

/* the gain (linear, <= 1 without a second knee's expansion) for each of
 * frames frames of key (interleaved stereo) into gain */
void udsp_comp_gain(UDSP_COMP *c, const float *key, float *gain, size_t frames);

/* the static curve: the level out (dB) for a steady level in (dB) */
float udsp_comp_curve(const UDSP_COMP *c, float in_db);

/* --- the up / down compressor ------------------------------------------
 * A band of Ableton's Multiband Dynamics (OTT's), Kilohearts' Dynamics
 * too: above one threshold the level's slope is one of its own (under 1
 * pulled down -- compression; over 1 pushed up), below another one too
 * (under 1 pulled up toward it -- upward compression; over 1 pushed down
 * -- expansion), in between left alone; a soft knee of `knee` dB round
 * each threshold (Giannoulis et al.'s, as the compressor's; two knees
 * overlapping add up). The gain (dB, either sign) is smoothed as the
 * compressor's: attack while it falls, release while it rises -- so a
 * boost comes in at the release's pace once the level drops. A boost
 * stops at `up_max` dB, a key under `floor` dB (silence, a tail's last
 * bit) gets none from below -- upward compression would otherwise lift
 * the noise floor without end --, and an expansion's cut stops at
 * UDSP_UDCOMP_DEPTH dB.
 *
 * Its detector (`rms`): 0 the peak (the louder channel's), 1 the RMS
 * over 10 ms (the louder channel's square) -- each then the gain smoothed
 * as above -- or UDSP_UDCOMP_POWER, Kilohearts' Dynamics' (measured): the
 * two channels' mean square over a window of its own
 * (udsp_udcomp_set_window; 10 ms until set; a channel alone reads 3 dB
 * down), then the level smoothed in dB -- attack while it rises, release
 * while it falls -- and the curve's gain taken as it comes. The same as
 * smoothing the gain where the curve compresses; an expansion's cut
 * comes in at the release's pace, as the level falls. Its floor holds
 * the level (a boost from below stays what it is there) rather than
 * ending the boost. */

#define UDSP_UDCOMP_POWER 2

#define UDSP_UDCOMP_DEPTH 120.0f

typedef struct {
  float a_thr, a_slope; /* above: its slope - 1 */
  float b_thr, b_slope; /* below: 1 - its slope */
  float up_max, floor;
  float att, rel, ms_coef;
  int rms;
  float db; /* the smoothed gain, dB */
  float ms;
  float lvl;            /* the power detector's smoothed level, dB */
  float knee, knee_inv; /* dB; 1 / (2 knee), 0 without one */
} UDSP_UDCOMP;

void udsp_udcomp_reset(UDSP_UDCOMP *c);

/* thresholds in dB (below <= above), ratios >= 1 (1e9 or more:
 * infinity), times in ms (> 0), up_max >= 0 and floor in dB */
void udsp_udcomp_set(UDSP_UDCOMP *c, float above_db, float above_ratio,
                     float below_db, float below_ratio, float up_max_db,
                     float floor_db, float attack_ms, float release_ms, int rms,
                     unsigned rate);

/* the same by slopes (dB out a dB in) beyond each threshold, >= 0 --
 * above: 1 / ratio compressing, over 1 expanding upward; below: 1 / ratio
 * compressing upward, over 1 expanding (0: infinity either way) -- and a
 * knee >= 0 dB; udsp_udcomp_set is slopes 1 / ratio and no knee */
void udsp_udcomp_set_slopes(UDSP_UDCOMP *c, float above_db, float above_slope,
                            float below_db, float below_slope, float knee_db,
                            float up_max_db, float floor_db, float attack_ms,
                            float release_ms, int rms, unsigned rate);

/* the power detector's window (ms > 0), after a set (which puts it back
 * at 10 ms) */
void udsp_udcomp_set_window(UDSP_UDCOMP *c, float ms, unsigned rate);

/* the gain (linear) for each of frames frames of key into gain */
void udsp_udcomp_gain(UDSP_UDCOMP *c, const float *key, float *gain,
                      size_t frames);

/* n (<= UDSP_UDCOMP_N) of them side by side, compressor j's gain from
 * key[j] into gain[j]: the same bits as n calls, but their recurrences
 * (RMS's, the smoother's) run frame by frame together, so their chains
 * overlap -- a multiband's bands */
#define UDSP_UDCOMP_N 4
void udsp_udcomp_gains(UDSP_UDCOMP *const c[], const float *const key[],
                       float *const gain[], size_t n, size_t frames);

/* the static curve: the level out (dB) for a steady level in (dB) */
float udsp_udcomp_curve(const UDSP_UDCOMP *c, float in_db);

/* --- the gate ----------------------------------------------------------
 * Opens when the key's peak passes the threshold, closes once it falls
 * `ret` dB under it (hysteresis: no chatter on a level that hovers) and
 * a hold time has passed; the gain moves to 1 or the floor through
 * one-poles in dB, attack opening, release closing. */

typedef struct {
  float open_db, close_db, floor_db;
  float att, rel;
  uint32_t hold, held; /* frames: the hold, and how long it's been */
  int open;
  float db; /* the gain now, dB (floor .. 0) */
} UDSP_GATE;

/* the gate shut, its gain at the floor */
void udsp_gate_reset(UDSP_GATE *g);

/* the floor in dB (-inf is -100 or under: silence); ret >= 0 */
void udsp_gate_set(UDSP_GATE *g, float threshold_db, float return_db,
                   float floor_db, float attack_ms, float hold_ms,
                   float release_ms, unsigned rate);

void udsp_gate_gain(UDSP_GATE *g, const float *key, float *gain, size_t frames);

/* --- the limiter --------------------------------------------------------
 * A brickwall over a lookahead of n frames: the gain each frame needs is
 * held at its minimum over the window to come, recovers at the release
 * rate, and is averaged over the window -- so it reaches every peak's
 * gain smoothly by the time the delayed peak gets there, and no sample
 * out passes the ceiling. Its memory is the caller's:
 * UDSP_LIMITER_FLOATS(cap) floats for a lookahead of up to cap frames.
 * Its delay is the lookahead: a host's latency. */

#define UDSP_LIMITER_FLOATS(cap) ((size_t)(cap) * 5 + 5)

typedef struct {
  float *mem;
  uint32_t cap, n; /* the most frames, and the lookahead now (1..cap) */
  float ceiling, rel;
  float env;           /* the held gain after its release */
  double sum;          /* the window's sum of env */
  uint64_t t;          /* frames in */
  uint32_t head, tail; /* the minimum's queue */
  float reduction;     /* the last frame's gain applied (a meter's) */
} UDSP_LIMITER;

/* takes mem (UDSP_LIMITER_FLOATS(cap)) and clears it */
void udsp_limiter_init(UDSP_LIMITER *l, float *mem, uint32_t cap);

/* the lookahead in frames (1..cap; a change clears the memory -- it is
 * a latency change, a host's to make between plays), the ceiling
 * linear, the release in ms */
void udsp_limiter_set(UDSP_LIMITER *l, uint32_t lookahead, float ceiling,
                      float release_ms, unsigned rate);

/* in place, interleaved stereo: the output lags the input by the
 * lookahead */
void udsp_limiter_stereo(UDSP_LIMITER *l, float *io, size_t frames);

#endif
