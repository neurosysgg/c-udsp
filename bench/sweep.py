#!/usr/bin/env python3
"""Sweeps one knob of the kernels per build: every value built (bench.c and
ulp.c), each candidate held to ulp.c's tolerances first, the passing ones
timed pinned to one core. One table per build: the kernels the knob
touches, ns per call under each value, the fastest passing value marked.

  bench/sweep.py KNOB [-b BUILD,...] [-k KERNEL] [-r REPS] [--cpu N]
                 [--save LABEL]

KNOB is one of KNOBS below (the default value first, marked "(default)").
A cell is the median ns; "!" = over a tolerance (not a candidate), "*" =
its bits differ from the default's (allowed: within the tolerance), "<"
= the fastest candidate, when it beats the default by more than NOISE.
--save LABEL writes bench/results/sweep-KNOB-LABEL.tsv.

Run on a quiet machine: the timing is matrix.py's (taskset, medians)."""

import argparse
import os
import sys
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matrix  # noqa: E402
import ulp  # noqa: E402

# knob -> (values as extra flags, default first; the kernels it touches:
# bench.c's names, ulp.c's names)
KNOBS = {
    # independent accumulators in every reduction (simd.h)
    "acc": (["-DUDSP_ACC=8", "-DUDSP_ACC=1", "-DUDSP_ACC=2", "-DUDSP_ACC=4",
             "-DUDSP_ACC=12", "-DUDSP_ACC=16"],
            ["rms", "rms_big", "peak", "peak_big", "sum", "minmax",
             "minmax_big", "bars128_planned8k"],
            ["rms", "rms_big", "sum", "peak", "minmax", "spectrum_power"]),
    # the pragma that fully unrolls loops over accumulator arrays (simd.h)
    "unroll": (["", "-DUDSP_NO_UNROLL"],
               ["rms", "peak", "sum", "minmax", "fft_real", "block_fused",
                "scope"],
               ["rms", "sum", "peak", "minmax", "fft_real", "scope"]),
    # floats summed in float lanes before going into the double (analyze.c)
    "chunk": (["-DRMS_CHUNK=4096", "-DRMS_CHUNK=512", "-DRMS_CHUNK=1024",
               "-DRMS_CHUNK=16384", "-DRMS_CHUNK=65536",
               "-DRMS_CHUNK=1048576"],
              ["rms", "rms_big", "sum"],
              ["rms", "rms_big", "sum"]),
    # fused multiply-adds (-std=c11 makes gcc's default off; clang's is on)
    "contract": (["-ffp-contract=off", "-ffp-contract=fast"],
                 ["fft_real", "block_fused", "bars128_planned8k", "lufs",
                  "osc_sine", "osc_saw", "svf", "softclip", "mix_add"],
                 None),
}
DEFAULT = ["gcc-O2", "gcc-O2-v3", "clang-O2-v3"]
# a winner beats the default by more than this, or isn't one; run to run,
# a kernel moves by up to ~4 % (clang's block_fused: 1289-1424 ns), so a
# "<" counts once a second run shows it too
NOISE = 0.05


def label(v):
    return v.split("=")[-1] if v else "(none)"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("knob", choices=KNOBS)
    ap.add_argument("-b", "--builds", default=",".join(DEFAULT))
    ap.add_argument("-k", "--kernel")
    ap.add_argument("-r", "--reps", type=int, default=21)
    ap.add_argument("--cpu", type=int, default=5)
    ap.add_argument("--save")
    a = ap.parse_args()
    values, bkernels, ukernels = KNOBS[a.knob]
    if a.kernel:
        bkernels = [k for k in bkernels if a.kernel in k]
    builds = a.builds.split(",")

    # the ULP check of every candidate, built and run in parallel
    def check(bv):
        b, v = bv
        rows = ulp.run(matrix.build(b, "ulp", v, "-" + a.knob + label(v)))
        if ukernels is not None:
            rows = {k: r for k, r in rows.items() if k in ukernels}
        return bv, rows
    jobs = [(b, v) for b in builds for v in values]
    with ThreadPoolExecutor(os.cpu_count() or 4) as ex:
        checked = dict(ex.map(check, jobs))
        bins = dict(zip(jobs, ex.map(
            lambda bv: matrix.build(bv[0], "bench", bv[1],
                                    "-" + a.knob + label(bv[1])), jobs), strict=True))

    rows_out = []
    for b in builds:
        base = checked[(b, values[0])]
        times = {}
        for v in values:
            ok = not ulp.over(checked[(b, v)])
            times[v] = {}
            if ok:
                for k in bkernels:
                    times[v].update(matrix.run(bins[(b, v)], a.cpu, k,
                                               a.reps))
            print(f"  {b} {label(v)} done", file=sys.stderr)
        print(f"\n{b}: {a.knob}")
        w = max(len(k) for k in bkernels)
        head = f"{'kernel':<{w}}" + "".join(
            f" {label(v) + (' (default)' if v == values[0] else ''):>16}"
            for v in values)
        print(head)
        print("-" * len(head))
        for k in bkernels:
            have = {v: t[k][1] for v, t in times.items() if k in t}
            best = min(have, key=have.get) if have else None
            line = f"{k:<{w}}"
            for v in values:
                rows = checked[(b, v)]
                same_k = [u for u in (k, k.removesuffix("_big"))
                          if u in rows][:1] or list(rows)
                differs = any(rows[u][4] != base[u][4] for u in same_k)
                if v not in have:
                    cell = "!"
                else:
                    t = have[v]
                    cell = f"{t:.0f}" if t >= 100 else f"{t:.1f}"
                    cell += "*" if differs else ""
                    if v == best and values[0] in have and \
                            t < have[values[0]] * (1 - NOISE):
                        cell += "<"
                line += f" {cell:>16}"
                if v in have:
                    rows_out.append((b, k, label(v), have[v],
                                     max((r[1] for r in rows.values()),
                                         default=0.0)))
            print(line)
        for v in values:
            bad = ulp.over(checked[(b, v)])
            if bad:
                print(f"  {label(v)}: over tolerance in {', '.join(bad)}")
    if a.save:
        os.makedirs(matrix.RESULTS, exist_ok=True)
        with open(os.path.join(matrix.RESULTS,
                               f"sweep-{a.knob}-{a.save}.tsv"), "w") as f:
            for b, k, v, t, e in rows_out:
                f.write(f"{b}\t{k}\t{v}\t{t:.1f}\t{e:.3f}\n")


if __name__ == "__main__":
    main()
