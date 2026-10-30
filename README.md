# c-µdsp

_Published on GitHub as `c-udsp` (`µ` → `u`) — GitHub repo names are
ASCII-only. Local checkouts and the ctui-mus submodule path both keep
the literal `c-µdsp` name; only the remote URL differs._

A small, dependency-free C11 DSP library for audio visualizers: FFT,
log-spaced spectrum bars, a spectrogram's scrolling history, momentary
LUFS loudness, integrated (EBU R128) loudness for ReplayGain,
RMS/peak/dBFS, mono downmix, oscilloscope-style triggered trace
extraction, and display ballistics (attack/release smoothing, peak
hold) -- and a synthesis and effects half for
[`ctui-daw`](../ctui-daw): oscillators, a sampler voice, an envelope
ramp, a parameter smoother, a state-variable filter, a feedback delay, a soft clip, pan laws
and a mixer's sum. The analysis was extracted from [`ctui-mus`](../ctui-mus)'s
`audiovis-core`, where this pipeline first took shape as the numeric
backbone behind four terminal visualizers (a level meter, an
oscilloscope, spectrum bars, a spectrogram).

## Status

Every analysis function below is stateless (or, for the two genuinely
stateful modules — the spectrogram's scrolling history and R128's running
measurement — explicit about it); the synthesis half is stateful by
nature (a phase, a filter's memory, a delay line) and keeps all of it in
caller-owned structs. Everything has a passing headless unit test. This is the whole library; there is no
app or example binary here, only `src/`, `tests/`, `bench/` and `tools/`.

## Build

```sh
make -j$(nproc)         # build/O2/: the objects, libudsp.a, the suites
make test [T=NAME]      # the suites (tests/*.c), all cores, a line each
                        # (tools/run_tests.sh; outputs in build/test/)
make SAN=1 test         # the same under ASan+UBSan (build/O2-san/)
make test-valgrind      # the same under memcheck
make test-variants      # gcc/clang, -O2/-O3, SSE2/AVX2, the plain-C build
                        # (UDSP_SCALAR) and ASan+UBSan, each with bench/ulp.c's
                        # tolerances (build/variant/)
make lint               # clang's warnings + analyzer, clang-tidy, cppcheck vs
                        # tools/lint-baseline.txt; the layering (tools/deps.py),
                        # clang-format, ruff, shellcheck
make check              # test-variants, the plain build + test, lint
make COV=1 coverage     # clang coverage per directory (tools/coverage.py;
                        # HTML in build/cov/html/)
make format             # clang-format -i (.clang-format: ctui's)
make compile_commands.json   # clangd's database
make bench              # bench/matrix.py: each kernel's cost per compiler/flags
make clean              # removes build/
```

Requires a C11 compiler and libc's `<math.h>` (linked via `-lm`) —
nothing else. No POSIX, no threads, no allocation anywhere in `src/`.

A consumer compiles `src/*/*.c` into its own build with `-Isrc` (or
links `build/O2/libudsp.a`) and includes headers by their directory.

## Usage

There's no single umbrella header — include only what you need
(`DESIGN.md` "Layout" has what each directory is for):

```c
/* src/base/: the vectors and the transform */
#include "base/simd.h"        /* UDSP_VF, udsp_vload/vstore, ... -- for kernels */
#include "base/fft.h"         /* udsp_fft_forward, udsp_window_hann, udsp_fft_magnitude;
                                 the table-driven ones: udsp_fft_table, udsp_fft_real,
                                 udsp_fft_real_inverse, udsp_fft_forward_table,
                                 udsp_fft_power, udsp_window_hann_table,
                                 udsp_window_apply */

/* src/analysis/: measuring a signal */
#include "analysis/analyze.h"     /* udsp_analyze_rms/peak/minmax/sum/dbfs, udsp_mono_downmix */
#include "analysis/spectrum.h"    /* udsp_spectrum_bars (= udsp_spectrum_magnitude +
                                     udsp_spectrum_bin) -- log-spaced dBFS bars from an
                                     UDSP_FFT_N-point FFT; UDSP_SPECTRUM_PLAN (any FFT
                                     size up to 8192) + udsp_spectrum_bars_planned /
                                     udsp_spectrum_power / udsp_spectrum_bin_planned --
                                     the same, ~18-20x faster */
#include "analysis/spectrogram.h" /* UDSP_SPECTROGRAM + udsp_spectrogram_push(_mag) --
                                     scrolling history */
#include "analysis/loudness.h"    /* udsp_loudness_lufs_momentary -- K-weighted momentary LUFS */
#include "analysis/r128.h"        /* UDSP_R128 + udsp_r128_add/integrated -- a program's
                                     integrated loudness (EBU R128, gated) and sample
                                     peak, fed in chunks; UDSP_R128_HIST adds up an album */
#include "analysis/truepeak.h"    /* UDSP_TRUEPEAK -- BS.1770's true peak */
#include "analysis/chroma.h"      /* peaks, pitch classes, root, chord, key */
#include "analysis/yin.h"         /* udsp_yin -- a note's fundamental */
#include "analysis/scope.h"       /* udsp_scope_extract -- triggered waveform trace */
#include "analysis/ballistics.h"  /* udsp_ballistics, udsp_peak_hold -- attack/release,
                                     peak hold */

/* src/synth/: making one (lifted from µOS's utone / utracker; DESIGN.md) */
#include "synth/pitch.h"     /* udsp_pitch_ratio_q32, udsp_note_millihz, udsp_exp2 --
                                exact equal temperament in integers */
#include "synth/osc.h"       /* UDSP_OSC + udsp_osc_render -- sine, PolyBLEP saw /
                                square, triangle, noise; a 32-bit phase accumulator */
#include "synth/wtab.h"      /* UDSP_WTAB -- band-limited wavetables */
#include "synth/lfo.h"       /* UDSP_LFO -- six shapes, sample and hold */
#include "synth/adsr.h"      /* UDSP_ADSR -- linear attack, exponential decay/release */
#include "synth/openv.h"     /* UDSP_OPENV -- Operator's envelope */
#include "synth/ramp.h"      /* UDSP_RAMP -- a linear per-frame envelope ramp */
#include "synth/smooth.h"    /* UDSP_SMOOTH -- a parameter gliding to a target in
                                a fixed time, landing on it exactly */
#include "synth/sampler.h"   /* UDSP_SAMPLE, UDSP_SAMPLER_VOICE -- a sample at a
                                pitch, Q32 position, linear interpolation, loops */
#include "synth/grain.h"     /* UDSP_GRAIN -- a windowed stretch of a sample at a
                                pitch, either way, cubic: a granular's grain */

/* src/fx/: changing one, in a render */
#include "fx/svf.h"          /* UDSP_SVF -- Chamberlin state-variable filter */
#include "fx/zdf.h"          /* UDSP_ZDF -- Simper's trapezoidal SVF, mono */
#include "fx/eq.h"           /* UDSP_EQ_BAND -- a parametric EQ's band */
#include "fx/xover.h"        /* UDSP_XOVER -- a Linkwitz-Riley crossover */
#include "fx/disperse.h"     /* UDSP_DISPERSE -- allpasses in series */
#include "fx/bank.h"         /* UDSP_BANK, UDSP_FOLLOW -- band-passes side by side
                                and their envelopes (a vocoder's) */
#include "fx/line.h"         /* UDSP_LINE -- a mono delay line, fractional reads */
#include "fx/velvet.h"       /* UDSP_VELVET -- a velvet-noise decorrelator (sparse FIR) */
#include "fx/delay.h"        /* UDSP_DELAY -- stereo feedback delay, caller's ring */
#include "fx/fdn.h"          /* UDSP_FDN -- a feedback delay network reverb */
#include "fx/dyn.h"          /* UDSP_COMP, UDSP_UDCOMP, UDSP_GATE, UDSP_LIMITER --
                                dynamics; udsp_log2f / udsp_exp2f */
#include "fx/shape.h"        /* udsp_softclip -- a rational soft clip */
#include "fx/gain.h"         /* udsp_gain_apply/stereo, udsp_mix_add, udsp_pan_gains */
#include "fx/conv.h"         /* UDSP_CONV_IN + udsp_conv_ir/in_push/mac/out --
                                uniformly partitioned convolution (an IR reverb's) */

/* src/offline/: changing one ahead of a render */
#include "offline/sinc.h"    /* UDSP_SINC -- a Kaiser-windowed sinc resampler */
#include "offline/stretch.h" /* udsp_stretch -- a phase vocoder's time stretch */
#include "offline/dither.h"  /* UDSP_DITHER -- TPDF dither ahead of a quantizer */
#include "offline/irfb.h"    /* udsp_irfb -- an impulse response fed back through itself */
```

Every function takes plain `float *` buffers the caller owns (input,
output, and scratch alike) — nothing in `src/` allocates, so a caller
controls its own memory layout and lifetime entirely. See each header's
doc comment for the exact buffer-size contract; `DESIGN.md` covers how
the modules layer on top of each other and why.

## Speed

The hot loops are written against `src/base/simd.h` (GCC vector extensions at
the target's width, plain C without them) and the transforms take
caller-owned tables filled once, so they run at the machine's speed at
plain `-O2`: a 2048-point spectrum (window, real FFT, power) is 1.8 µs on
a Ryzen 7 5800X, 23x the old path's 41.7 µs; at 8192 points 128 bars take
10 µs planned against 180 µs; RMS 15x, peak 4x (13x with AVX2), the
period-locked scope 6.5x (same traces). `bench/README.md` has the numbers per compiler and flags and what
was learned getting there (accumulators vs unrolling, what -O3, AVX2 and
FMA do or don't buy). The old entry points keep their behaviour and cost
(`udsp_spectrum_bars()` et al.); new code wants the planned/table ones.

## Docs

- `DESIGN.md` — how the modules layer (transform → domain-specific
  analysis → stateful history), and the deliberate simplifications
  each one makes.
- `bench/README.md` — the benchmarks, the measurements, the lessons.
- `CLAUDE.md` — working conventions specific to this repo.

## Relationship to ctui-mus (and ctui-wm)

ctui-wm vendors this repo too, as `c-udsp/` (the GitHub name), for its
audio visualizer and ctui-audio; it builds against a local checkout with
`make UDSP_DIR=../c-µdsp` while a change here isn't pushed yet.


`ctui-mus` pulls this repo in as a git submodule at `vendor/c-µdsp`,
sourced from the local `../c-µdsp` checkout — a real, independent
working copy, not a frozen vendor drop, the same relationship `ctui-mus`
already has with [`ctui`](../ctui) at `vendor/ctui`. `ctui-mus`'s
`audiovis-core` builds its threading, ring-buffer tap, and lock-free
snapshot publishing on top of the pure functions here; none of that
app-specific glue lives in this repo.

## License

MIT.
