#!/usr/bin/env python3
"""Builds bench/bench.c once per compiler x flags, runs each pinned to one
core, and prints one table: the median ns per call of every kernel under
every build.

  bench/matrix.py [-p bench|reduce] [-b BUILD,...] [-k KERNEL] [-r REPS]
                  [--cpu N]
                  [--save LABEL] [--compare LABEL [LABEL2]]

--save LABEL writes bench/results/LABEL.tsv (every build's numbers);
--compare A shows A against the run just made (or A against B), as the
ratio old/new per kernel and build (> 1 = faster now). Builds live in
bench/build/ (gitignored)."""

import argparse
import glob
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RESULTS = os.path.join(HERE, "results")

# name -> (compiler, flags); the baseline build is x86-64 (SSE2 only)
BUILDS = {
    "gcc-O2": ("gcc", "-O2"),
    "gcc-O2-unroll": ("gcc", "-O2 -funroll-loops"),
    "gcc-O3": ("gcc", "-O3"),
    "gcc-O3-unroll": ("gcc", "-O3 -funroll-loops"),
    "gcc-O2-v3": ("gcc", "-O2 -march=x86-64-v3"),
    "gcc-O3-v3": ("gcc", "-O3 -march=x86-64-v3"),
    "gcc-O3-v3-unroll": ("gcc", "-O3 -march=x86-64-v3 -funroll-loops"),
    "gcc-O3-native": ("gcc", "-O3 -march=native"),
    # -std=c11 makes gcc's default -ffp-contract=off: no FMA without this
    "gcc-O2-v3-fma": ("gcc", "-O2 -march=x86-64-v3 -ffp-contract=fast"),
    "gcc-O3-v3-fma": ("gcc", "-O3 -march=x86-64-v3 -ffp-contract=fast"),
    # the plain-C kernels: what the vector paths are measured against
    "gcc-O2-scalar": ("gcc", "-O2 -DUDSP_SCALAR"),
    "clang-O2": ("clang", "-O2"),
    "clang-O3": ("clang", "-O3"),
    "clang-O2-v3": ("clang", "-O2 -march=x86-64-v3"),
    "clang-O3-v3": ("clang", "-O3 -march=x86-64-v3"),
}
DEFAULT = ["gcc-O2", "gcc-O2-unroll", "gcc-O3", "gcc-O3-unroll", "gcc-O2-v3",
           "gcc-O3-v3", "gcc-O3-v3-unroll", "clang-O2", "clang-O3",
           "clang-O3-v3"]


def build(name, prog, extra="", tag=""):
    """prog's binary for build name (bench, ulp: with the library; reduce:
    alone); extra flags on top (a sweep's -D...), tag names that variant"""
    cc, flags = BUILDS[name]
    out = os.path.join(HERE, "build", prog + "-" + name + tag)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    src = sorted(glob.glob(os.path.join(ROOT, "src", "*", "*.c"))) \
        if prog != "reduce" else []
    cmd = [cc, "-std=c11", "-Wall", "-Wextra", "-g", *flags.split(),
           *extra.split(),
           "-I" + os.path.join(ROOT, "src"), "-o", out,
           os.path.join(HERE, prog + ".c"), *src, "-lm"]
    subprocess.run(cmd, check=True)
    return out


def run(binary, cpu, kernel, reps):
    cmd = ["taskset", "-c", str(cpu), binary, "-r", str(reps)]
    if kernel and "/reduce-" not in binary:
        cmd += ["-k", kernel]
    out = subprocess.run(cmd, check=True, capture_output=True, text=True)
    rows = {}
    for line in out.stdout.splitlines():
        k, n, med, mn, per = line.split("\t")
        rows[k] = (int(n), float(med), float(mn))
    return rows


def table(results, builds, ref=None):
    kernels = []
    for b in builds:
        for k in results[b]:
            if k not in kernels:
                kernels.append(k)
    w = max(len(k) for k in kernels)
    head = f"{'kernel':<{w}} {'n':>8}" + "".join(f" {b:>16}" for b in builds)
    print(head)
    print("-" * len(head))
    for k in kernels:
        n = next(results[b][k][0] for b in builds if k in results[b])
        line = f"{k:<{w}} {n:>8}"
        for b in builds:
            cell = ""
            if k in results[b]:
                med = results[b][k][1]
                cell = f"{med:.0f}" if med >= 100 else f"{med:.1f}"
                if ref and b in ref and k in ref[b]:
                    cell += f" x{ref[b][k][1] / med:.2f}"
            line += f" {cell:>16}"
        print(line)


def save(label, results):
    os.makedirs(RESULTS, exist_ok=True)
    with open(os.path.join(RESULTS, label + ".tsv"), "w") as f:
        for b, rows in results.items():
            for k, (n, med, mn) in rows.items():
                f.write(f"{b}\t{k}\t{n}\t{med:.1f}\t{mn:.1f}\n")


def load(label):
    results = {}
    with open(os.path.join(RESULTS, label + ".tsv")) as f:
        for line in f:
            b, k, n, med, mn = line.rstrip("\n").split("\t")
            results.setdefault(b, {})[k] = (int(n), float(med), float(mn))
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-b", "--builds", default=",".join(DEFAULT))
    ap.add_argument("-k", "--kernel")
    ap.add_argument("-r", "--reps", type=int, default=21)
    ap.add_argument("--cpu", type=int, default=5)
    ap.add_argument("-p", "--prog", default="bench",
                    help="bench (the kernels) or reduce (the study)")
    ap.add_argument("--save")
    ap.add_argument("--compare", nargs="+")
    a = ap.parse_args()
    builds = a.builds.split(",")
    for b in builds:
        if b not in BUILDS:
            sys.exit(f"unknown build {b} (know: {', '.join(BUILDS)})")

    if a.compare and len(a.compare) == 2:
        old, new = load(a.compare[0]), load(a.compare[1])
        builds = [b for b in builds if b in new]
        print(f"{a.compare[1]} (ratio: {a.compare[0]} / {a.compare[1]})")
        table(new, builds, old)
        return

    results = {}
    for b in builds:
        results[b] = run(build(b, a.prog), a.cpu, a.kernel, a.reps)
        print(f"  {b} done", file=sys.stderr)
    if a.save:
        save(a.save, results)
    ref = load(a.compare[0]) if a.compare else None
    table(results, builds, ref)


if __name__ == "__main__":
    main()
