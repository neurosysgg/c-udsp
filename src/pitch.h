#ifndef UDSP_PITCH_H
#define UDSP_PITCH_H

#include <stdint.h>

/* Equal-tempered pitch as exact integers (lifted from µOS's utracker,
 * song.c): a semitone table in Q20 and octaves as shifts, so an octave
 * is exactly 2x and only the pitch class carries a rounding (one table
 * entry; within 0.0006 cents of pow() over the MIDI range). Integers so
 * a renderer's pitch is the same number on every build and target. */

/* the lowest and highest MIDI notes udsp_note_millihz() answers (C0..B8
 * in the tracker's spelling, C-0..B-9) */
#define UDSP_NOTE_MIN 12
#define UDSP_NOTE_MAX 119

/* 2^(semitones/12) in Q32: a playback-rate ratio, any sign; octaves
 * clamp at +20 / -40 */
uint64_t udsp_pitch_ratio_q32(int semitones);

/* a MIDI note's frequency in millihertz (69 = 440000), 0 outside
 * UDSP_NOTE_MIN..UDSP_NOTE_MAX */
uint32_t udsp_note_millihz(int note);

#endif
