# c-µdsp design

This library is an extraction: every function here originally lived
under `ctui-mus/src/audiovis/`, built for a terminal music player's four
visualizer widgets (a level meter, an oscilloscope, spectrum bars, a
spectrogram). `ctui-mus`'s `DSP_DESIGN.md` documents that original
context in full (the shared sample window, the consumer thread, the
per-tick pipeline); this doc is the subset of it that's actually about
the numeric library itself, with the app-specific plumbing (threading,
ring buffers, snapshots, a fixed per-tick sample window) stripped out —
none of that belongs here, and nothing in `src/` knows any of it exists.
A caller supplies its own buffers on its own schedule; that's the entire
contract.

## Pipeline layering: transform → domain-specific analysis → stateful history

```
fft.c            generic radix-2 FFT + Hann window. Knows nothing about
                 audio (no Hz, no dB, no channels) -- a pure numerical
                 transform, reusable for anything that needs one.
   │
   ▼
spectrum.c       udsp_spectrum_bars(): windowed FFT -> n log-spaced dBFS
                 bars. Stateless -- one call, one answer. n_bars is a
                 runtime argument, not a constant this file owns, so the
                 same function can back both a "few wide bars" widget and
                 a "many narrow bins" spectrogram column.
   │                                              │
   ▼                                              ▼
a spectrum-bar widget                    spectrogram.c
one call per frame, wider bars           calls udsp_spectrum_bars() with
                                          UDSP_SPECTROGRAM_BINS, then keeps
                                          a scrolling byte-quantized
                                          history -- the one genuinely
                                          stateful module in this library,
                                          and it's its own file *because*
                                          it's stateful, not bundled into
                                          the stateless spectrum.c.

analyze.c        RMS/peak/dBFS/mono-downmix -- basic signal-level
                 primitives with no frequency-domain concept at all.
   │
   ▼
loudness.c       udsp_loudness_lufs_momentary(): K-weighting biquad
                 cascade + mean-square -> LUFS. Built on analyze.c's
                 dBFS floor constant, not on spectrum.c -- loudness and
                 spectral analysis are unrelated axes that happen to
                 often feed the same class of meter widget.

scope.c          udsp_scope_extract(): trigger search + trace copy,
                 built on nothing but a caller-supplied mono buffer.
```

The rule this follows: **a file boundary marks either a change in what
kind of computation is happening (time domain → frequency domain →
perceptual scaling) or a change from stateless to stateful** — not "one
file per visualizer." `spectrum.c` and `spectrogram.c` are separate files
for the second reason (state); `fft.c` and `spectrum.c` are separate for
the first (a numeric transform vs. an audio-domain interpretation of
one). None of `src/*.c` know anything about display sizing beyond their
own defaults (`spectrum.h` owns its own default bar count, `scope.h` its
own trace length, `spectrogram.h` its own grid dimensions) — a caller is
always free to override with a runtime argument where one exists
(`n_bars`), or fork the header's constant locally where one doesn't
(`UDSP_FFT_N`, `UDSP_SCOPE_POINTS`).

**This library only ever computes numbers — it never renders anything.**
Turning a dBFS bar, an intensity byte, or a trace sample into a color or
a terminal cell is entirely the caller's concern.

## Simplifications, and why each one is acceptable here

None of these are secret; each is called out in the relevant header's
doc comment too. Collected here so "is this good enough" has one answer
instead of several scattered ones:

- **LUFS is momentary-only, not BS.1770's integrated/gated loudness.**
  The K-weighting filter cascade and the `-0.691 + 10*log10(meanSquare)`
  formula are the real spec; the multi-block gating algorithm BS.1770
  actually defines for a *program's* loudness is not attempted. Correct
  for "how loud does this sound right now" (a live meter bar); not a
  certified measurement.
- **The K-weighting filter coefficients are the spec's fixed 48kHz
  design, applied unmodified at whatever the caller's real sample rate
  is.** Re-deriving the coefficients per sample rate means implementing
  BS.1770's analog-prototype bilinear-transform annex — real work with
  no payoff for a UI meter reading against a fixed dB scale. At common
  sample rates (44.1/48kHz), the weighting curve's shape is close enough
  that the difference is inaudible-to-the-eye on a bar.
- **The spectrum's dBFS calibration ignores the Hann window's coherent-
  gain loss.** `UDSP_SPECTRUM_NORM` (`spectrum.c`) is the textbook `2/N`
  single-sided-amplitude normalization, not `2/N` divided by the
  window's own ~0.5 gain factor. A full-scale sine reads a few dB below
  true 0dBFS as a result. Fine for a relative, log-scaled bar/waterfall
  display; would need fixing for anything claiming calibrated SPL.
- **The spectrogram quantizes to a byte per cell, not a float.**
  `UDSP_SPECTROGRAM.intensity` is `unsigned char [COLS][BINS]`. A real
  precision loss (256 levels instead of a continuous dB value), but the
  typical consumer is a color lookup that was always going to quantize
  to some number of levels eventually — doing it once here, rather than
  at render time, keeps the struct an order of magnitude smaller for a
  caller that copies it often. A caller needing the pre-quantized float
  should add its own parallel float buffer rather than expect this one
  to carry it.
- **The oscilloscope's trigger is "first rising zero-crossing in the
  search window," not a real scope's configurable edge/level/holdoff.**
  Good enough to stop a trace from visibly swimming frame to frame;
  nothing here attempts more than that.
- **Loudness re-filters its entire input from scratch on every call
  instead of carrying persistent biquad state across calls.** Simpler
  (no cross-call state to get wrong or reset) at the cost of redundant
  work if a caller feeds it overlapping windows repeatedly — an
  acceptable trade for a library with no real-time constraint of its
  own; a caller with a hard per-call budget should size its input
  accordingly.

## What's explicitly out of scope

**No threading, no allocation, no I/O, anywhere in `src/`.** Every
function is a plain, synchronous, single-threaded transform over
caller-owned buffers — deciding how often to call these functions, on
which thread, and how to get audio samples into a `float *` in the first
place is entirely up to the caller. `ctui-mus`'s `audiovis-core` is one
example of that decision (a dedicated consumer thread, a ring-buffer tap,
a fixed-size shared sample window, lock-free snapshot publishing) but
none of it is assumed or required here.

## Adding a new analysis

1. Decide which existing module it's really an alternate *view* of. A
   new frequency-domain analysis almost certainly wants
   `udsp_spectrum_bars()` (possibly with its own bar/bin count) or a new
   sibling to `spectrogram.c` if it needs history; a new time-domain one
   wants the caller's raw buffer directly, the way `scope.c` does.
2. If it needs a genuinely new *kind* of analysis (not a re-parameterized
   call to something that exists), give it its own file under `src/`,
   named for the computation, not for any one caller's use of it.
3. Keep it allocation-free and caller-buffer-based: any working memory a
   new function needs should be a `scratch` parameter the caller owns,
   the same convention `spectrum.c` and `loudness.c` already use.
4. Add a test file under `tests/` following the existing ones' shape —
   known-input/known-output assertions (a DC signal's energy lands in
   FFT bin 0, silence reads the dBFS/LUFS floor, a tone lands in the
   expected log-spaced bar) rather than testing against another
   implementation.
