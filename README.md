# c-µdsp

A small, dependency-free C11 DSP library for audio visualizers: FFT,
log-spaced spectrum bars, a spectrogram's scrolling history, momentary
LUFS loudness, RMS/peak/dBFS, mono downmix, and oscilloscope-style
triggered trace extraction. Extracted from [`ctui-mus`](../ctui-mus)'s
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
make test    # builds and runs everything under tests/
make clean   # removes built test binaries
```

Requires a C11 compiler and libc's `<math.h>` (linked via `-lm`) —
nothing else. No POSIX, no threads, no allocation anywhere in `src/`.

## Usage

There's no single umbrella header — include only what you need:

```c
#include "fft.h"        /* udsp_fft_forward, udsp_window_hann, udsp_fft_magnitude */
#include "analyze.h"     /* udsp_analyze_rms/peak/dbfs, udsp_mono_downmix */
#include "spectrum.h"    /* udsp_spectrum_bars -- log-spaced dBFS bars from an FFT */
#include "loudness.h"    /* udsp_loudness_lufs_momentary -- K-weighted momentary LUFS */
#include "scope.h"       /* udsp_scope_extract -- triggered waveform trace */
#include "spectrogram.h" /* UDSP_SPECTROGRAM + udsp_spectrogram_push -- scrolling history */
```

Every function takes plain `float *` buffers the caller owns (input,
output, and scratch alike) — nothing in `src/` allocates, so a caller
controls its own memory layout and lifetime entirely. See each header's
doc comment for the exact buffer-size contract; `DESIGN.md` covers how
the modules layer on top of each other and why.

## Docs

- `DESIGN.md` — how the modules layer (transform → domain-specific
  analysis → stateful history), and the deliberate simplifications
  each one makes.
- `CLAUDE.md` — working conventions specific to this repo.

## Relationship to ctui-mus

`ctui-mus` pulls this repo in as a git submodule at `vendor/c-µdsp`,
sourced from the local `../c-µdsp` checkout — a real, independent
working copy, not a frozen vendor drop, the same relationship `ctui-mus`
already has with [`ctui`](../ctui) at `vendor/ctui`. `ctui-mus`'s
`audiovis-core` builds its threading, ring-buffer tap, and lock-free
snapshot publishing on top of the pure functions here; none of that
app-specific glue lives in this repo.

## License

MIT.
