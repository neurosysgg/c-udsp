# Working on c-µdsp

`DESIGN.md` has the pipeline layering and the deliberate simplifications
each module makes — read it first. This file is conventions specific to
this repo.

## What this repo is (and isn't)

A small, standalone, libc-only DSP library — FFT, spectrum bars, a
spectrogram's scrolling history, momentary LUFS, RMS/peak/dBFS, mono
downmix, triggered scope extraction. It was extracted out of
`ctui-mus`'s `audiovis-core` because the actual numeric pipeline there
had zero dependency on ctui-mus's own threading/ring-buffer/snapshot
machinery, or on ctui-mus at all — it's a general-purpose transform over
`float *` buffers a caller owns. Nothing in `src/` should ever grow a
dependency on any consumer's app structure; if a change here can only be
justified by "because that's how ctui-mus's audiovis-core happens to
call it," it belongs in the consumer, not here.

## Naming

`UDSP_` for types/macros, `udsp_` for functions, prefixed by module
(`udsp_fft_*`, `udsp_analyze_*`, `udsp_spectrum_*`, `udsp_loudness_*`,
`udsp_scope_*`, `udsp_spectrogram_*`) — the same shape `ctui`'s own
`CTUI_`/`ctui_` and `ctui-mus`'s `MUS_`/`mus_` conventions use, scoped to
this library.

## Code style

- No comments except where the *why* isn't obvious from the code: a
  hidden constraint, a subtle invariant, a documented simplification
  (see `DESIGN.md`). Never restate what a well-named line already says.
- Trust internal callers; validate at boundaries. `n` a power of two for
  `udsp_fft_forward()`, buffer sizes matching a header's documented
  contract — these are the caller's responsibility, not re-checked
  defensively inside the library.
- No allocation, no threading, no I/O, anywhere in `src/`. Every
  function takes its input, output, and scratch buffers as plain
  caller-owned `float *` parameters. This is the one rule every other
  convention here serves — it's what keeps the library trivially
  embeddable in any threading model a consumer picks.
- Stateless by default. `spectrogram.c` is the one exception (a
  spectrogram is inherently "spectrum bars plus memory of recent
  history") and it says so explicitly in its own header — don't add a
  second stateful module without the same kind of justification.

## Testing approach

Plain, dependency-free unit tests under `tests/`, using the header-only
`udsp_test.h` assert/summary pair (`make test` builds and runs all of
them). Prefer known-input/known-output assertions over comparing against
a second implementation — a DC signal's energy lands entirely in FFT bin
0, silence reads the dBFS/LUFS floor exactly, a 1kHz tone at a given
sample rate lands in a specific, pinned log-spaced bar index. Pinning
exact values (not just "greater than the floor") catches a regression in
the math immediately rather than letting it drift.

## Consumers

`ctui-mus` is the first and, so far, only consumer — pulled in via
`git submodule add ../c-µdsp vendor/c-µdsp`, the same pattern it already
uses for `vendor/ctui`. Editing this repo from inside a consumer's
working tree (e.g. `ctui-mus/vendor/c-µdsp/src/...`) is editing this
actual repo directly — commit here, then bump the consumer's pinned
submodule pointer, the same two-step `vendor/ctui` workflow
`ctui-mus/CLAUDE.md` already documents. Never leave a session with an
uncommitted change sitting inside a consumer's `vendor/c-µdsp` checkout.

## Workflow

Only commit when explicitly asked. Stage files explicitly. Never push
without being asked.
