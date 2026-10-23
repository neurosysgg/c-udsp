#!/usr/bin/env python3
"""What the suites ran of src/ (clang's source-based coverage).

  make COV=1 coverage [T=NAME]   build, run the suites, this report
  tools/coverage.py [--html] [--files] [--under PATH] BINARY...

`make COV=1 coverage` builds the suites with -fprofile-instr-generate
-fcoverage-mapping (objects in build/O2-cov/), runs them (each process
writes build/cov/raw/PID-SIG.profraw), then calls this with the
binaries: the raw counts merged (build/cov/merged.profdata) and reported
per directory (src/base, src/fx, ...): lines and functions run, out of
how many. Only src/ counts (not the tests or the bench). Taken from
ctui-wm's tools/coverage.py (by way of ctui-mus's).

--files lists each file instead (the least covered first), --under PATH
only those under PATH (`--under src/fx --files`), --html also writes
llvm-cov's annotated sources to build/cov/html/ (index.html).
"""
import argparse
import collections
import glob
import json
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
COV = os.path.join(ROOT, "build", "cov")
IGNORE = r"(^|/)(tests|tools|bench)/|^/usr/"


def instrumented(path):
    """an object with coverage mapping (COV=1)"""
    p = subprocess.run(["llvm-readelf", "-S", "-W", path],
                       capture_output=True, text=True)
    return "__llvm_covmap" in p.stdout


def module(path):
    parts = path.split("/")
    return "/".join(parts[:2]) if len(parts) > 2 else parts[0]


def pct(a, n):
    return f"{100 * a / n:5.1f} %" if n else "    -  "


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("binaries", nargs="+")
    ap.add_argument("--html", action="store_true")
    ap.add_argument("--files", action="store_true")
    ap.add_argument("--under", default="")
    args = ap.parse_args()

    raw = glob.glob(os.path.join(COV, "raw", "*.profraw"))
    if not raw:
        sys.exit("no counts in build/cov/raw/: run `make COV=1 coverage`")
    prof = os.path.join(COV, "merged.profdata")
    subprocess.run(["llvm-profdata", "merge", "-sparse", "-o", prof] + raw,
                   check=True)
    objs = [b for b in args.binaries if os.path.exists(b) and instrumented(b)]
    if not objs:
        sys.exit("none of the binaries is instrumented (COV=1)")
    base = ["-instr-profile", prof, objs[0]] + \
        [a for o in objs[1:] for a in ("-object", o)] + \
        ["-ignore-filename-regex", IGNORE]
    out = subprocess.run(["llvm-cov", "export", "-summary-only"] + base,
                         capture_output=True, text=True, cwd=ROOT)
    if out.returncode != 0:
        sys.exit(out.stderr)
    files = []
    for f in json.loads(out.stdout)["data"][0]["files"]:
        path = os.path.relpath(f["filename"], ROOT)
        if path.startswith(args.under):
            s = f["summary"]
            files.append((path, s["lines"]["covered"], s["lines"]["count"],
                          s["functions"]["covered"], s["functions"]["count"]))

    if args.files:
        rows = sorted(files, key=lambda r: (r[1] / r[2] if r[2] else 1, -r[2]))
        print(f"{'file':52} {'lines':>16} {'functions':>16}")
        for path, lc, ln, fc, fn in rows:
            print(f"{path:52} {pct(lc, ln)} {ln:>7} {pct(fc, fn)} {fn:>6}")
    else:
        mods = collections.defaultdict(lambda: [0, 0, 0, 0, 0])
        for path, lc, ln, fc, fn in files:
            m = mods[module(path)]
            m[0] += 1
            m[1] += lc
            m[2] += ln
            m[3] += fc
            m[4] += fn
        print(f"{'module':24} {'files':>5} {'lines':>16} {'functions':>16}")
        for name, (nf, lc, ln, fc, fn) in sorted(
                mods.items(), key=lambda kv: kv[1][1] / kv[1][2]
                if kv[1][2] else 1):
            print(f"{name:24} {nf:5} {pct(lc, ln)} {ln:>7} "
                  f"{pct(fc, fn)} {fn:>6}")
    lc = sum(r[1] for r in files)
    ln = sum(r[2] for r in files)
    fc = sum(r[3] for r in files)
    fn = sum(r[4] for r in files)
    print(f"{'all':24} {len(files):5} {pct(lc, ln)} {ln:>7} "
          f"{pct(fc, fn)} {fn:>6}")

    if args.html:
        html = os.path.join(COV, "html")
        subprocess.run(["rm", "-rf", html])
        subprocess.run(["llvm-cov", "show", "-format=html",
                        "-show-line-counts-or-regions",
                        "-output-dir", html] + base, check=True, cwd=ROOT)
        print(f"html: {os.path.relpath(html, ROOT)}/index.html")
    return 0


if __name__ == "__main__":
    sys.exit(main())
