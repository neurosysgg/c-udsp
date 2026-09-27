# bench: what the kernels cost, and why they're written the way they are

```sh
make bench                            # bench/matrix.py: every build, one table
bench/matrix.py -k fft -b gcc-O2,gcc-O3-v3   # some kernels, some builds
bench/matrix.py --save LABEL          # keep a run (bench/results/LABEL.tsv)
bench/matrix.py --compare OLD [NEW]   # old/new per kernel and build
bench/matrix.py -p reduce             # the reduction study (reduce.c)
make test-variants                    # the tests under every build below
```

`bench.c` times each kernel on synthetic audio (a bass tone, a wobbling
melody, some noise): the iteration count calibrated to ~10 ms per run,
then 21 runs, the median reported (ns per call). `matrix.py` builds it
once per compiler and flags and runs each build pinned to one core
(`taskset -c 5`). `reduce.c` is the study behind every reduction here.

The reference machine is an AMD Ryzen 7 5800X (Zen 3: AVX2 + FMA, no
AVX-512, 16 vector registers; ~4.7 GHz boost, `performance` governor),
gcc 16.2 and clang 22.1. `results/baseline.tsv` is the library before the
speed pass, `results/after.tsv` after it, `results/reduce-study.tsv` the
study.

## Builds

| build | flags |
|---|---|
| gcc-O2 | `-O2` (x86-64 baseline: SSE2) |
| gcc-O2-unroll | `-O2 -funroll-loops` |
| gcc-O3 | `-O3` |
| gcc-O2-v3, gcc-O3-v3 | `-march=x86-64-v3` (AVX2 + FMA) |
| gcc-O3-v3-fma | `... -ffp-contract=fast` (see "FMA") |
| gcc-O3-native | `-O3 -march=native` (znver3 tuning) |
| clang-O2, clang-O3-v3 | the same with clang |

## Results

| kernel | before, gcc -O2 | after, gcc -O2 | after, gcc -O3 -march=x86-64-v3 | after, clang -O3 -march=x86-64-v3 |
|---|---|---|---|---|
| Hann window, 2048 | 16.69 µs | 121 ns (137.8x) | 64 ns (259.9x) | 83 ns (200.8x) |
| FFT of 2048 real samples | 23.90 µs | 1.52 µs (15.7x) | 1.17 µs (20.4x) | 1.08 µs (22.1x) |
| window + FFT + power, 2048 (a spectrum's core) | 41.69 µs | 1.73 µs (24.1x) | 1.32 µs (31.6x) | 1.19 µs (35.1x) |
| spectrum bars, 28 of 2048 | 42.47 µs | 2.05 µs (20.7x) | 1.64 µs (25.9x) | 1.49 µs (28.6x) |
| spectrum bars, 128 of 2048 | 44.48 µs | 2.75 µs (16.2x) | 2.34 µs (19.0x) | 2.16 µs (20.6x) |
| spectrum bars, 128 of 8192 (ctui-mus's size; before = the old API now) | 183.01 µs | 10.10 µs (18.1x) | 8.88 µs (20.6x) | 8.31 µs (22.0x) |
| magnitude, 1024 bins | 1.10 µs | 277 ns (4.0x) | 140 ns (7.9x) | 140 ns (7.9x) |
| rms, 2048 | 1.31 µs | 85 ns (15.5x) | 52 ns (25.2x) | 49 ns (26.6x) |
| rms, 1M | 676.74 µs | 49.16 µs (13.8x) | 36.84 µs (18.4x) | 38.14 µs (17.7x) |
| peak, 2048 | 464 ns | 112 ns (4.2x) | 35 ns (13.1x) | 32 ns (14.3x) |
| peak, 1M | 238.76 µs | 55.30 µs (4.3x) | 35.99 µs (6.6x) | 37.95 µs (6.3x) |
| stereo downmix, 2048 frames | 2.04 µs | 233 ns (8.7x) | 236 ns (8.6x) | 347 ns (5.9x) |
| momentary LUFS, 19200 | 86.59 µs | 50.08 µs (1.7x) | 50.01 µs (1.7x) | 57.48 µs (1.5x) |
| period-locked scope, 2048 (before = its scalar version) | 164.93 µs | 25.18 µs (6.6x) | 13.85 µs (11.9x) | 13.07 µs (12.6x) |
| min+max, 2048 (new) | - | 129 ns | 67 ns | 60 ns |
| sum, 1024 (new) | - | 38 ns | 18 ns | 19 ns |

Sizes are 2048 where not said (the baseline's UDSP_FFT_N; 8192 since the
merge with ctui-mus's line, which only the old spectrum API uses now).
"before" of the FFT row is the old complex transform on the same input
(im = 0); of the bars rows `udsp_spectrum_bars()`, which stays as it was
(`udsp_spectrum_bars_planned()` is its fast path); of the spectrum-core
row the old hann + fft + magnitude (now the window fused into the real
FFT, then power). The old functions keep their
behaviour: `udsp_fft_forward()`, `udsp_window_hann()` and
`udsp_spectrum_bars()` cost what they did, `udsp_fft_magnitude()`, the
analysis kernels and the loudness got faster in place. One regression:
clang -O3 -march=x86-64-v3 had auto-vectorised the plain stereo downmix 8
wide (144 ns); the explicit 4-lane shuffle is 345 there (gcc: 233 in every
build).

## What was learned

### Reductions: accumulators, not unrolling

A sum, a max or a min over samples written the plain way
(`for (i) s += x[i] * x[i]`) is one dependency chain: each add waits on
the one before it, 3 cycles on Zen 3. Without `-ffast-math` a compiler may
not reassociate float math, so it can't split the chain itself, and no
flag changes that: the study's plain sum of squares runs at the same ~0.62
ns/sample in every build, gcc or clang, -O2 or -O3, unrolled or not, SSE2
or AVX2. gcc -O3 even reports that loop "vectorized", but as an in-order
("fold-left") reduction: vector loads feeding the same scalar chain, no
faster.

What helps is independent accumulators, so several chains are in flight:
the machine takes two FP adds or FMAs per cycle at 3-4 cycles latency, so
~8 chains keep it busy. `reduce.c` over 2048 floats (L1) and 4 MB (L3):

| sum of squares, 2048 floats | gcc-O2 | gcc-O3 | gcc-O3-v3 |
|---|---|---|---|
| plain loop | 1300 ns | 1283 | 1262 |
| 8 plain-C lanes (`acc[8]`) | 163 | 164 | 155 |
| 16 plain-C lanes | 334 | 96 | 85 |
| 1 vector accumulator (vector extension) | 751 | 752 | 151 |
| 4 vector accumulators | 181 | 232 | **49** |

and the same for a max of |x|. So:

- **`-funroll-loops` doesn't split float accumulators** (gcc's variable
  expansion needs `-fassociative-math`); it only saves loop overhead,
  which on Zen 3's wide front end is noise (±10% either way across the
  kernels). Unrolling isn't the lever; independent chains are.
- **Plain-C lane arrays work but are fragile at -O2**: 8 lanes 163 ns, 16
  lanes 334, 32 lanes 195. gcc's -O2 cost model ("very cheap") takes or
  leaves them unpredictably; -O3 is steadier.
- **GCC vector extensions at the target's width are the reliable form.**
  Wider than the target they are a trap: a 32-byte max on an SSE2 build
  was lowered lane by lane, 2662 ns, 6x *slower* than the plain loop.
- **4 MB runs at the same per-sample speed as 2 KB**: nothing here is
  memory-bound on this machine (the L3 feeds 100+ GB/s).

Hence `src/simd.h`: one vector type at the target's width (16 bytes on
SSE2/NEON, 32 with AVX), `UDSP_ACC` = 8 accumulators in every reduction,
and a scalar build (`UDSP_SCALAR`) as the plain-C reference the tests hold
the vector build to.

### Where unrolling does matter: small loops over vector arrays

The accumulators live in an array (`UDSP_VF acc[UDSP_ACC]`) walked by a
`for (a < UDSP_ACC)` loop. gcc -O3 unrolls that loop completely, the array
becomes eight registers. **gcc -O2 doesn't**, the array stays on the
stack, and every step becomes a store and a reload: rms ran at 189 ns at
-O2 against 83 at -O3 (13 stack accesses in the loop instead of 3). A
`#pragma GCC unroll` on exactly those loops (`UDSP_UNROLL`; clang honours
it too) made -O2 equal to -O3 (85 ns). The same happened in the FFT's
fused first pass (1.95 µs at -O2 vs 1.55 at -O3, equal after the pragma).
Manual unrolling in the source isn't needed; telling gcc -O2 to fully
unroll the few constant-count loops over register arrays is.

### maxps: gcc doesn't fold the blend

A per-lane max in vector extensions is a compare and a blend (vector
`?:` is C++-only in gcc). clang turns it back into `maxps`; gcc emits
`cmpltps` + `pand` + `pandn` + `por`, in the integer domain (a bypass
delay each way on Zen 3): peak ran at 170 ns on SSE2 against 112 with
`_mm_max_ps`. `simd.h` uses the target's instruction where it has one
(SSE/AVX, NEON falls back to the blend); `maxps` is exactly `a > b ? a :
b`, so every build answers like the scalar one.

### The FFT: the scalar parts were the cost

The old transform recomputed its twiddles by a recurrence (serial, and
drifting: 1.9e-5 relative error at 4096 points, where exact twiddles give
1.4e-7) and windowed with a `cos()` per sample (16.7 µs of a 2048-point
spectrum's 42). With tables (caller-owned, filled once, contiguous per
stage so a butterfly loads plain vectors), radix-4 stages (half the loads
and stores of radix-2) and a real-input transform (n real points through
an n/2 complex one), `fft_real` measured, stage by stage, by removing one
part at a time (at -O2):

| part | before fusing | after |
|---|---|---|
| deinterleave copy, scalar | 0.45 µs | fused |
| bit-reversal scatter | 0.35 µs | fused |
| first radix-4 pass, scalar | 0.29 µs | fused: ~0.3 µs for all three |
| radix-4/2 stages, vector | 0.80 µs | 0.80 µs |
| split into the real spectrum, vector | 0.25 µs | 0.25 µs |

A third of the time went into scalar setup before any vector stage ran,
which is also why AVX2 gained almost nothing then. The fix fuses the
three: in bit-reversed order, the four inputs of first-pass group g come
from four contiguous stretches of the input (z[r], z[r + m/2], z[r + m/4],
z[r + 3m/4], r = g reversed), so lanes run over r with plain loads, the
butterfly runs across lanes, and a 4x4 transpose stores each group as
one 16-byte store: 2.33 -> 1.53 µs at -O2, 1.16 with AVX2.

Two effects checked and found small: register spills (the radix-4 body
needs ~20 live vectors; the hot AVX loop has none, the spills are in the
scalar edges) and 4K aliasing between `re` and `im` (5-6% at worst).

### FMA

`-std=c11` puts gcc in `-ffp-contract=off`: no fused multiply-adds unless
asked. Asked (`-ffp-contract=fast`), they help nothing measurably here
and cost the loudness filter 15%: its critical path is a multiply-add
per sample, and a 4-cycle FMA is slower on that chain than the pieces it
replaces when the multiply is off the chain. Not used.

### -O2, -O3, AVX2

With the kernels written against `simd.h`, **-O3 gives nothing over -O2**
for them; it only still helps loops left to the auto-vectoriser (the
scope's early-exit search, ballistics). **AVX2** is worth 1.3x on the
FFT and 1.5-3x on the reductions, but in absolute terms a spectrum is
1.85 vs 1.35 µs, and analysing a 20-minute take ~60 vs ~45 ms: not worth
a run-time dispatch in the library. Build with `-march=x86-64-v3` (or
`native`) where the binary only runs on such machines.

### Loudness

The K-weighting cascade is two IIR stages, serial by nature. Fused into
one pass (stage 2 of a sample overlaps stage 1 of the next) in transposed
Direct Form II, with the mean square in the same loop: 86.6 -> 50 µs per
400 ms window (1.7x), same result to 1e-4 LUFS. Past that it would take
a block-parallel (look-ahead) IIR; a meter at 30 fps spends 1.5 ms/s on
it, not worth it.

### Cache, alignment, copies

The kernels are compute-bound hot: the same work over 4 MB runs at the
per-sample speed of 2 KB (the L3 streams at 100+ GB/s). But a real
caller runs a spectrum now and then, not in a loop, and finds its data
wherever the rest of the machine left it. ctui-wm's `tools/vizbench -c`
measures that (it evicts before each frame): its visualizer's frame
(spectrum, ballistics, paint) was 3.8 µs hot, 4.2 from L3 and 8.6 from
RAM -- with the compute gone, cache warmth decides the frame. What that
led to here:

- **Fewer bytes per spectrum.** `udsp_fft_real_windowed()` applies the
  window while the first pass loads its input: no windowed copy (8 KB)
  and one pass less, bit-exact with windowing first. The planned
  spectrum's scratch went from 16 to 8 KB (`UDSP_SPECTRUM_SCRATCH`).
  Hot: bars 2.25 -> 2.03 µs.
- **Cache lines, not alignment as such.** Unaligned loads cost nothing on
  Zen 3 unless one straddles two 64-byte lines. The scratch's imaginary
  half now starts on a line of its own, and `UDSP_SPECTRUM_PLAN` has its
  tables first: aligned where the caller aligns the plan (`aligned_alloc`,
  `_Alignas(64)`).
- The caller's side (in ctui-wm): the ring of recent samples mirrored
  (one memfd mapped twice) so the spectrum reads it in place, and frames
  whose bars didn't move neither painted nor compared. Together about
  10% hot and 6% cold; the rest of a cold frame is code and page tables,
  not data.

What can't be had: pinning data in a cache (no cache locking on x86;
AMD's L3 QoS is a server feature and needs root). What can: small
working sets, passes fused so data is touched once while hot (a 2048-
sample block is analysed whole while it sits in L1), buffers reused
rather than allocated per call, strides that aren't multiples of 4 KB
(4K aliasing cost the FFT 5-6% when `re` and `im` were 8 KB apart), and
for large output nobody reads back, non-temporal stores.

Copies: the transforms are in place after their first pass, which is
itself the input's copy (split, permuted, windowed). The data movement
left is the algorithm's: every radix-4 pass reads and writes the block
once (a self-sorting FFT would trade the permutation for ping-pong
buffers, no less movement).
