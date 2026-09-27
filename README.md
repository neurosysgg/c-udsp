# c-µdsp

_Published on GitHub as `c-udsp` (`µ` → `u`) — GitHub repo names are
ASCII-only. Local checkouts and the ctui-mus submodule path both keep
the literal `c-µdsp` name; only the remote URL differs._

A small, dependency-free C11 DSP library for audio visualizers: FFT,
log-spaced spectrum bars, a spectrogram's scrolling history, momentary
LUFS loudness, RMS/peak/dBFS, mono downmix, oscilloscope-style
triggered trace extraction, and display ballistics (attack/release
smoothing, peak hold). Extracted from [`ctui-mus`](../ctui-mus)'s
`audiovis-core`, where this pipeline first took shape as the numeric
backbone behind four terminal visualizers (a level meter, an
oscilloscope, spectrum bars, a spectrogram).

## Status

Every function below is stateless (or, for the one genuinely stateful
module — the spectrogram's scrolling history — explicit about it) and
has a passing headless unit test. This is the whole library; there is no
app or example binary here, only `src/` and `tests/`.

## Build

```sh
make test            # builds and runs everything under tests/
make test-variants   # the same under gcc/clang, -O2/-O3, SSE2/AVX2, the
                     # plain-C build (UDSP_SCALAR) and ASan+UBSan
make bench           # bench/matrix.py: each kernel's cost per compiler/flags
make clean           # removes built test binaries
```

Requires a C11 compiler and libc's `<math.h>` (linked via `-lm`) —
nothing else. No POSIX, no threads, no allocation anywhere in `src/`.

## Usage

There's no single umbrella header — include only what you need:

```c
#include "fft.h"        /* udsp_fft_forward, udsp_window_hann, udsp_fft_magnitude;
                           the table-driven ones: udsp_fft_table, udsp_fft_real,
                           udsp_fft_forward_table, udsp_fft_power,
                           udsp_window_hann_table, udsp_window_apply */
#include "analyze.h"     /* udsp_analyze_rms/peak/minmax/sum/dbfs, udsp_mono_downmix */
#include "spectrum.h"    /* udsp_spectrum_bars (= udsp_spectrum_magnitude +
                           udsp_spectrum_bin) -- log-spaced dBFS bars from an
                           UDSP_FFT_N-point FFT; UDSP_SPECTRUM_PLAN (any FFT
                           size up to 8192) + udsp_spectrum_bars_planned /
                           udsp_spectrum_power / udsp_spectrum_bin_planned --
                           the same, ~18-20x faster */
#include "gain.h"        /* udsp_gain_apply -- a buffer times a gain */
#include "loudness.h"    /* udsp_loudness_lufs_momentary -- K-weighted momentary LUFS */
#include "scope.h"       /* udsp_scope_extract -- triggered waveform trace */
#include "spectrogram.h" /* UDSP_SPECTROGRAM + udsp_spectrogram_push(_mag) -- scrolling history */
#include "ballistics.h"  /* udsp_ballistics, udsp_peak_hold -- attack/release, peak hold */
```

Every function takes plain `float *` buffers the caller owns (input,
output, and scratch alike) — nothing in `src/` allocates, so a caller
controls its own memory layout and lifetime entirely. See each header's
doc comment for the exact buffer-size contract; `DESIGN.md` covers how
the modules layer on top of each other and why.

## Speed

The hot loops are written against `src/simd.h` (GCC vector extensions at
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
