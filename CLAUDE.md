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

## Layout

`src/` is grouped by what a module does (DESIGN.md "Layout"): `base/`
(simd, fft), `analysis/`, `synth/`, `fx/`, `offline/`. Headers are
included by their path under `src/` (`#include "synth/osc.h"`), here and
in the consumers. `tools/deps.py` (in `make lint`) holds the layering --
a new include across directories either fits its ALLOWED table or the
table changes, with DESIGN.md saying why -- and src/'s system headers
(SYSTEM: libc's pure ones and the intrinsics). A new module goes into
the directory it belongs to, its test into `tests/NAME_test.c`.

## Code style

- clang-format (`.clang-format`, ctui's): `make format`; `make lint`
  fails on a file it would change.
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
- Stateless by default in the analysis half. `spectrogram.c` and
  `r128.c` are the exceptions (history, a program's measurement) and say
  so in their headers — don't add another stateful analysis without the
  same kind of justification. The synthesis half (osc, sampler, ramp,
  svf, delay; DESIGN.md) is stateful by nature: its state is a struct
  the caller owns, never a static and never allocated here.

## Speed

- Loops over samples use `src/base/simd.h` (`UDSP_VF`, `udsp_vload`/`vstore`,
  `udsp_vmax`/...), never raw intrinsics outside it and never a vector
  wider than the target. Reductions: `UDSP_ACC` accumulators, and
  `UDSP_UNROLL` on every loop over an accumulator array (gcc -O2 needs
  it; see bench/README.md). Every vector loop keeps a scalar tail.
- Anything that only depends on sizes/rates goes into a caller-owned
  table or plan, filled once, read-only after.
- A speed change is measured, not assumed: `bench/matrix.py --compare
  after` before and after it (add the kernel to bench.c if it's new), and
  `bench/README.md` updated when the picture changes.
- A speed change keeps the answers: `bench/ulp.py` (`make ulp`) before
  calling it done -- every kernel within its written-down tolerance in
  ULPs (`bench/ulp.c`; a new kernel gets a reference there). A knob worth
  trying across values goes into `bench/sweep.py`, which times only the
  candidates that pass.

## Testing approach

`make check` before calling anything done: `make test-variants` (every
build the kernels have a path for, plus ASan/UBSan, each with
bench/ulp.c's tolerances), the plain build's suites, then `make lint` --
lint's only new findings fail (fix them, or add them to
`tools/lint-baseline.txt` under a comment saying why they stay -- never
`--update` alone). `make SAN=1 test` and `make test-valgrind` are the
quick sanitizer runs between; `make COV=1 coverage` shows what the suites
reach. Objects live in `build/FLAVOUR/` (`O2`, `O2-san`, `O2-cov`,
`variant/...`) and don't track flags beyond the flavour: after changing
CFLAGS by hand, `make clean`. Plain, dependency-free unit
tests under `tests/`, using the header-only
`udsp_test.h` assert/summary pair (`make test` builds and runs all of
them). Prefer known-input/known-output assertions over comparing against
a second implementation — a DC signal's energy lands entirely in FFT bin
0, silence reads the dBFS/LUFS floor exactly, a 1kHz tone at a given
sample rate lands in a specific, pinned log-spaced bar index. Pinning
exact values (not just "greater than the floor") catches a regression in
the math immediately rather than letting it drift.

## Consumers

Each one compiles `src/*/*.c` itself, with `-I<this repo>/src`, and
includes headers by their directory; a module moving directories, or a
new directory, is a change in all of them (their Makefiles and includes)
the next time they bump the submodule.

`ctui-daw` (submodule `vendor/c-µdsp`) uses the synthesis half and holds
it to µOS's renders bit for bit (its `tests/utrk_test.c`): a change to
the arithmetic of osc / sampler / ramp / svf / delay / shape / gain is
checked there (`make test` in ../ctui-daw with `UDSP_DIR=../c-µdsp`).

`ctui-wm` is the second (submodule `c-udsp/`, GitHub's ASCII name). The
two lines diverged once (ctui-mus's vendored copy got four commits that
never came back here while this checkout got two others; merged
2026-09-27): a change made in a consumer's vendored copy is pushed back
to this checkout before anything else builds on either. `ctui-mus` is the first consumer — pulled in via
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
