#include "pitch.h"

static const uint64_t SEMITONE_Q20[12] = {
    1048576, 1110928, 1176987, 1246974, 1321123, 1399681,
    1482910, 1571089, 1664511, 1763488, 1868350, 1979448,
};

/* MIDI 0 in microhertz (440 / 2^(69/12)): three spare digits keep the
 * top octave's shift (up to 512x) inside a thousandth of a cent */
#define C_MINUS1_UHZ 8175799ull

uint64_t udsp_pitch_ratio_q32(int semitones) {
  /* floor division: -1 semitone is octave -1, pitch class 11 */
  int oct = semitones / 12;
  int pc = semitones % 12;
  if (pc < 0) {
    pc += 12;
    oct -= 1;
  }
  uint64_t ratio = SEMITONE_Q20[pc] << 12;
  if (oct > 0) {
    ratio <<= oct > 20 ? 20 : oct;
  } else if (oct < 0) {
    ratio >>= -oct > 40 ? 40 : -oct;
  }
  return ratio;
}

uint32_t udsp_note_millihz(int note) {
  if (note < UDSP_NOTE_MIN || note > UDSP_NOTE_MAX) {
    return 0;
  }
  uint64_t v = (C_MINUS1_UHZ * SEMITONE_Q20[note % 12]) >> 20;
  v <<= note / 12;
  return (uint32_t)(v / 1000ull);
}
