#!/usr/bin/env python3
"""Builds bench/ulp.c once per compiler x flags (matrix.py's builds), runs
them, and prints one table: each kernel's worst error in ULPs under every
build (ulp.c says how each is measured and against what).

  bench/ulp.py [-b BUILD,...] [-k KERNEL] [-D FLAGS]

A cell is the worst ULP error over the corpus; "!" after it = over the
kernel's written-down tolerance, "*" = its output bits differ from the
plain-C build's (gcc-O2-scalar, the reference the vector paths are held
to). Kernels measured as "bits" only show "=" or "*". -D adds flags to
every build (a candidate's -DUDSP_ACC=4, ...). Exits 1 if any kernel is
over its tolerance in any build; differing bits are reported, not failed
(contraction, -ffp-contract=fast, is allowed to change them)."""

import argparse
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matrix  # noqa: E402

REF = "gcc-O2-scalar"
DEFAULT = ["gcc-O2-scalar", "gcc-O2", "gcc-O3", "gcc-O2-v3", "gcc-O3-v3",
           "gcc-O3-native", "gcc-O2-v3-fma", "gcc-O3-v3-fma", "clang-O2",
           "clang-O2-v3", "clang-O3-v3"]


def run(binary, kernel=None):
    """ulp's rows: kernel -> (metric, max_ulp, worst_signal, tol, hash)"""
    cmd = [binary] + (["-k", kernel] if kernel else [])
    out = subprocess.run(cmd, capture_output=True, text=True)
    if out.returncode not in (0, 1):
        sys.exit(f"{binary}: exit {out.returncode}\n{out.stderr}")
    rows = {}
    for line in out.stdout.splitlines():
        k, metric, e, at, tol, h = line.split("\t")
        rows[k] = (metric, float(e), at, float(tol), h)
    return rows


def measure(builds, kernel=None, extra="", tag=""):
    """every build's rows, built and run in parallel"""
    def one(b):
        return b, run(matrix.build(b, "ulp", extra, tag), kernel)
    with ThreadPoolExecutor(os.cpu_count() or 4) as ex:
        return dict(ex.map(one, builds))


def over(rows):
    return [k for k, r in rows.items() if r[1] > r[3]]


def table(results, builds, ref):
    kernels = list(next(iter(results.values())))
    w = max(len(k) for k in kernels)
    head = f"{'kernel':<{w}} {'metric':<6} {'tol':>5}" + "".join(
        f" {b:>14}" for b in builds)
    print(head)
    print("-" * len(head))
    for k in kernels:
        metric, _, _, tol, _ = results[builds[0]][k]
        line = f"{k:<{w}} {metric:<6} {tol:>5g}"
        for b in builds:
            _, e, at, tol, h = results[b][k]
            same = ref is None or ref[k][4] == h
            if metric == "bits":
                cell = "=" if same else "*"
            else:
                cell = f"{e:.2f}" + ("!" if e > tol else "") + \
                    ("" if same else "*")
            line += f" {cell:>14}"
        print(line)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-b", "--builds", default=",".join(DEFAULT))
    ap.add_argument("-k", "--kernel")
    ap.add_argument("-D", "--flags", default="")
    a = ap.parse_args()
    builds = a.builds.split(",")
    for b in builds:
        if b not in matrix.BUILDS:
            sys.exit(f"unknown build {b} (know: {', '.join(matrix.BUILDS)})")
    want = builds if REF in builds else [REF] + builds
    results = measure(want, a.kernel, a.flags, "-x" if a.flags else "")
    table(results, builds, results[REF])
    bad = {b: over(results[b]) for b in builds if over(results[b])}
    for b, ks in bad.items():
        print(f"over tolerance in {b}: {', '.join(ks)}")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
