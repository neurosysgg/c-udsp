#ifndef UDSP_ANALYZE_H
#define UDSP_ANALYZE_H

#include <stddef.h>

/* floor applied by udsp_analyze_dbfs() below -- silence is -inf dB, which
 * is useless for a meter bar's fill fraction or a spectrum bar's height,
 * so every dB value in this library is clamped to this floor instead.
 * -60dBFS is a common "meter bottom" convention (well below typical noise
 * floors, well above -inf). */
#define UDSP_DBFS_FLOOR -60.0f

/* RMS over a flat float buffer, unscaled (plain ratio-to-full-scale, not
 * gained for visual punch -- that's a caller-side display concern, not a
 * signal-analysis one). Deliberately channel-agnostic: summing squares
 * over every float is the same computation whether or not consecutive
 * samples happen to belong to different channels, so this needs no
 * channel count at all. */
float udsp_analyze_rms(const float *samples, size_t n);

/* max absolute sample value over a flat float buffer, same channel-
 * agnostic reasoning as udsp_analyze_rms(). 0 for an empty buffer. */
float udsp_analyze_peak(const float *samples, size_t n);

/* interleaved (frames * channels floats) -> mono (frames floats), by
 * averaging each frame's channels. channels <= 0 or 1 is treated as
 * already-mono and just copies. This IS channel-aware (unlike the two
 * functions above) because a true per-sample time-domain waveform --
 * unlike a single scalar stat -- has no channel-agnostic equivalent: the
 * oscilloscope/FFT-based analyses need one waveform value per instant,
 * not per raw float. */
void udsp_mono_downmix(const float *interleaved, size_t frames, int channels,
                       float *mono_out);

/* linear amplitude ratio (1.0 = full scale) -> dBFS, clamped to
 * UDSP_DBFS_FLOOR (never -inf, even for exactly 0.0). The one dB
 * conversion every analysis in this library shares -- a meter's RMS/
 * peak, a spectrum bar's bin magnitude -- so the floor/clamp behavior is
 * defined exactly once. */
float udsp_analyze_dbfs(float amplitude);

#endif
