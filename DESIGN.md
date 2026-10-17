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
                                          history -- a genuinely stateful
                                          module (r128.c is the other),
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
   │
   ▼
r128.c           UDSP_R128: a program's integrated loudness (EBU R128 /
                 BS.1770-4) and sample peak, fed chunk by chunk -- the
                 K-weighting derived per sample rate, 400 ms blocks every
                 100 ms, channel weights, both gates, the blocks in a
                 fixed histogram (UDSP_R128_HIST; two add up to an
                 album). Stateful because a program is too long to hand
                 over whole without the caller allocating it; a separate
                 file from loudness.c for that reason, as spectrogram.c
                 is from spectrum.c.

yin.c            udsp_yin(): a note's fundamental (YIN: the normalized
                 difference function's first dip, a parabola through
                 it) -- for a sampler mapping multi-samples whose roots
                 aren't written down; offline, its sums in vector lanes
                 (so 0.4 ULP of a double reference, not the same bits
                 across widths: it never runs in a render).

scope.c          udsp_scope_extract(): trigger search + pitch-period
                 autocorrelation + resample, built on nothing but a
                 caller-supplied mono buffer.

ballistics.c     udsp_ballistics() / udsp_peak_hold(): how a displayed
                 value follows a measured one over time (attack/release
                 one-pole easing, peak hold then fall). Unit-agnostic
                 (dBFS bars or linear levels alike) and dt-based, so the
                 result doesn't depend on the caller's frame rate. The
                 state (level, peak, age) is caller-owned arrays updated
                 in place: nothing is kept inside the library, which is
                 why this isn't a second stateful module like
                 spectrogram.c.
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

## Tables and plans: precomputation without state

Whatever depends only on sizes and rates -- FFT twiddles, a window, a
spectrum's bar-to-bin ranges -- is computed once into a caller-owned
buffer (`udsp_fft_table()`, `udsp_window_hann_table()`) or struct
(`UDSP_SPECTRUM_PLAN`, `udsp_spectrum_plan()`) and passed back in. That
keeps the library's rules: nothing allocates (the caller sizes the
buffer: `UDSP_FFT_TABLE_FLOATS(n)`, `UDSP_SPECTRUM_SCRATCH`), and nothing
is state -- a table is read-only once filled, so one can serve any number
of calls and threads. The per-call work drops to what actually depends on
the samples. The table-free originals (`udsp_fft_forward()`,
`udsp_window_hann()`, `udsp_spectrum_bars()` and its magnitude/bin halves)
stay, unchanged, for callers that haven't moved.

The FFT size is a plan's, not the library's: `UDSP_FFT_N` (8192) sizes
the old entry points, but consumers want different trade-offs from one
build -- a live meter short windows that follow the music (2048: ~43 ms
at 48 kHz), a detailed display fine bins (8192: ~5.9 Hz). A plan holds
tables for the largest size and uses the first part for smaller ones, so
a small plan touches (and keeps in cache) only what it needs.

The twiddle table stores each stage's twiddles contiguously (stage h's
k-th at `[h + k]`, its imaginary half at `[n + h + k]`), so the butterfly
loop loads plain vectors, and one table built for n serves every smaller
power of two too (`udsp_fft_real()` runs an n/2-point transform plus a
split on the table for n).

## Vectors: `simd.h`

Every loop over samples is written against `simd.h`'s `UDSP_VF`: a GCC
vector extension (gcc, clang) as wide as the target's registers (16 bytes
on SSE2/NEON, 32 with AVX), with the target's own min/max/abs/sqrt where
the generic form compiles badly. Reductions keep `UDSP_ACC` (8)
independent accumulators, since without `-ffast-math` no compiler splits
a float sum's chain itself; the short loops over those accumulator arrays
carry `UDSP_UNROLL`, without which gcc -O2 keeps them in memory. Defining
`UDSP_SCALAR` (or a compiler without the extension) makes a vector one
float: the same kernels as plain C, the reference `make test-variants`
holds the vector builds to. `bench/README.md` has the measurements behind
each of these choices.

## The synthesis half: generators and effects

Added for ctui-daw (session 215 in ctui-wm's log; ctui-wm's
`notes/daw.md`; the synth's blocks -- lfo, adsr, zdf, FM, udsp_exp2 --
in session 223): the blocks a DAW's instruments and effects are built
from, so ctui-daw, ctui-mus, ctui-wm and µOS share one DSP library.
The first ones are lifted from µOS's utone and utracker, where they
were written host-portable (stdint only, caller-owned memory) and
tested against an independent renderer to the bit; their arithmetic is
kept exactly, so a renderer built on them still matches µOS's exports
sample for sample (ctui-daw's `tests/utrk_test.c`).

```
pitch.c     exact equal temperament: a Q20 semitone table, octaves as
            shifts (no pow(), no drift between octaves); udsp_exp2 and a
            fractional note's increment for a synth's pitch (no libm)
osc.c       UDSP_OSC: a 32-bit phase accumulator and five waves; the
            same phase-modulated (2-operator FM)
lfo.c       UDSP_LFO: six shapes read at control rate, a sample and hold
adsr.c      UDSP_ADSR: a linear attack, exponential decay and release
zdf.c       UDSP_ZDF: Simper's trapezoidal SVF, mono: a voice's filter,
            stable to 0.49 of the rate under any modulation
sampler.c   UDSP_SAMPLE + UDSP_SAMPLER_VOICE: Q32 position, linear
            interpolation, loops that keep their fraction; a cubic
            (Hermite) resampler over a caller's window at a fractional
            step (ctui-daw's sampler: a head in memory, the rest
            streamed, copied into the window per block)
ramp.c      UDSP_RAMP: a linear envelope ramp, one step per frame
smooth.c    UDSP_SMOOTH: a parameter's glide to a target, a fixed length
            whatever the distance, the last frame on the target exactly
svf.c       UDSP_SVF: Chamberlin's state-variable filter, stereo
delay.c     UDSP_DELAY: a stereo feedback delay over the caller's ring
shape.c     udsp_softclip: x / (1 + |x|), with makeup; a saturator's
            curves (soft, hard, fold, tube), no makeup
gain.c      + udsp_gain_stereo, udsp_mix_add, udsp_pan_gains
eq.c        UDSP_EQ_BAND: a parametric EQ's band, stereo -- Simper's SVF
            with its outputs mixed into eight shapes (cuts at 12 and 48
            dB an octave, shelves, a bell, a notch), its curve for a
            screen, and a ramp: g, k and the mix moving a frame at a
            time so a glide neither zippers nor leaves the stable region
dyn.c       a float log2 / exp2 of its own (no libm); UDSP_COMP (a soft
            knee, peak or RMS, attack and release in dB), UDSP_GATE
            (hysteresis, hold, a floor) -- both a gain from a key, the
            caller applying it (a sidechain) -- and UDSP_LIMITER, a
            brickwall over a lookahead (the needed gain held at its
            window's minimum, released, averaged over the window: no
            sample past the ceiling)
line.c      UDSP_LINE: a mono delay line a power of two long, read at a
            whole or fractional delay (Catmull-Rom between four)
fdn.c       UDSP_FDN: a reverb -- eight lines through a Hadamard matrix,
            each damped and its gain set from the decay, modulated reads,
            a predelay and four diffusing allpasses ahead
sinc.c      UDSP_SINC: an offline resampler -- a Kaiser-windowed sinc
            (120 dB), polyphase over the rates' reduced ratio (rows
            interpolated past 4096 phases), each row normalised to pass
            DC; input in windows of one channel, the output the same
            however it's cut (ctui-daw's audio clips, converted at import)
dither.c    UDSP_DITHER: TPDF dither ahead of a quantizer (two uniform
            half steps, a xorshift32 the caller seeds: the same bits on
            every build; ctui-daw's export at 16 / 24 bits)
```

The sinc (session 229, ctui-daw's d3d) keeps the same rules where they
matter for an offline kernel: its table comes from a sine and an I0 of
its own in double (libm's sin may round differently between libcs), and
its dot product runs four lanes in a fixed order -- lane j summing taps
j, j + 4, ..., then (0 + 1) + (2 + 3) -- whatever the vector width, so a
file converted on any build is the same file. The caller owns the table
(`udsp_sinc_floats`).

The effects (session 225, ctui-daw's d2f: its EQ Eight, dynamics,
delay, chorus, saturator and reverb are built on these) keep the
synthesis half's rules: the caller's memory (`UDSP_LIMITER_FLOATS`,
`udsp_fdn_floats`, `udsp_line_size`), no libm where it would change
bits between libcs (dB through dyn.c's log2 / exp2 or pitch.c's
`udsp_exp2`, tan as the oscillator's sine over its cosine), and the
same bits however a caller cuts its blocks -- the FDN's modulation
renormalised every 256 frames from its start, not per call.

**State.** A generator or a filter is a recurrence: its state can't be
handed back to the caller between samples the way ballistics.c's
arrays are, so each of these is a struct the caller owns and passes in
-- declared in the header, sized by the compiler, nothing kept in the
library and nothing allocated (utracker's fx.c kept a static pool of
instances and 500 ms delay rings; here the caller sizes the ring:
`UDSP_DELAY_RING_FLOATS`). Two structs never share anything, so any
number of voices run on any threads.

**Vectors.** Where a block's samples don't depend on each other, the
loop is written against `simd.h` like the analysis: the gains, the mix
sum, the soft clip, the ramp's held part, and the oscillators -- every
wave but noise is a function of the phase alone, and a vector's phases
are `phase + k * inc` (`UDSP_VU`, uint32 lanes that wrap as the scalar
does; `__builtin_convertvector` rounds as a scalar cast does). The
vector loop and the scalar tail are the same operations in the same
order, so a block split anywhere gives the same bits (the tests split
renders at odd places and compare with memcmp). Serial recurrences --
the SVF, the delay, noise, the sampler's interpolation -- stay scalar
per voice; their speed-up is several chains at once -- the SVF's two
channels in one loop, r128's channels two at a time (the same
operations per chain, so the same bits), and several voices in the
lanes of one vector, which is a polyphonic synth's layout (d2).

**No libm where it plays.** libm's exp2, tan and pow aren't required to
round the same in every libc (glibc and musl differ in last places), so
what a synth computes while it plays -- a gliding pitch, an envelope's
coefficient, a swept filter's tan -- comes from `udsp_exp2` (a Taylor
series and ldexp, exact scaling) and the oscillator's polynomial sine.

**Bits.** Built with `-ffp-contract=off` (ctui-daw does), every build --
gcc and clang, SSE2 and AVX2, UDSP_SCALAR -- renders the same bits; with
contraction, a compiler may fuse `a * b + c` and a render then differs
from another build in the last place. The tests here hold the kernels to
known values with tolerances where a chain of operations could fuse, so
they pass either way.

## Simplifications, and why each one is acceptable here

None of these are secret; each is called out in the relevant header's
doc comment too. Collected here so "is this good enough" has one answer
instead of several scattered ones:

- **`loudness.c`'s LUFS is momentary-only, not BS.1770's integrated/gated
  loudness.** The K-weighting filter cascade and the `-0.691 +
  10*log10(meanSquare)` formula are the real spec; the multi-block gating
  BS.1770 defines for a *program's* loudness is `r128.c`'s job, not the
  meter's. Correct for "how loud does this sound right now" (a live meter
  bar).
- **The meter's K-weighting coefficients are the spec's fixed 48kHz
  design, applied unmodified at whatever the caller's real sample rate
  is.** `r128.c` derives them per rate (`udsp_r128_kweight()`); the meter
  doesn't take a rate at all, and at common sample rates (44.1/48kHz) the
  weighting curve's shape is close enough that the difference is
  inaudible-to-the-eye on a bar.
- **R128's gated blocks are binned at 0.1 LU, not kept one by one.**
  Each bin sums its blocks' energies exactly, so the only approximation
  is the one bin the relative gate falls inside: it is counted in or out
  whole, by its own blocks' mean loudness. Checked against ffmpeg's
  `ebur128` on real tracks (44.1/48 kHz, stereo and mono): within 0.01 LU.
  A rate whose 100 ms isn't a whole number of samples (11025 Hz) rounds
  the hop; the peak is the sample peak, not BS.1770's true (oversampled)
  peak.
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
- **The tracker blocks keep utracker's shortcuts.** The SVF's
  coefficient uses a 3-term sine (0.01 % off at 10 kHz / 48 kHz) and is
  clamped below fs / 4, the triangle is naive (aliases), the sine is a
  5-term polynomial (~-109 dBFS off), the sampler interpolates linearly:
  each good enough for what it's for, and kept so the renders stay
  comparable with µOS's to the bit. Better ones (a band-limited
  triangle, a cubic or windowed-sinc sampler, a TPT filter) come as new
  functions beside them, not as changes to these (the cubic resampler
  and zdf.c since).
- **Vector sums and maxima combine in a different order than a plain
  loop.** Reductions run in 8 x lanes partial sums (flushed into a double
  every 4096 samples), so a sum can differ from a sequential one in the
  last bits; max/min/peak are exact. RMS over a 20-minute take stays
  within 1e-7 relative of a double reference (the tests pin 1e-6).
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
