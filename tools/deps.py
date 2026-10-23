#!/usr/bin/env python3
"""The include graph of src/, checked against the layering.

  tools/deps.py            the module graph, then the checks (exit 1 on a
                           finding: `make lint` runs it)
  tools/deps.py --files    every file-level edge, "from -> to"

A module is a directory under src/. tests/ and bench/ include whatever
they test and aren't checked beyond their includes resolving.

Checks:
- layering: each module includes only what ALLOWED names for it (plus
  itself); a new edge either fits the table or the table changes, with
  DESIGN.md saying why
- system headers: src/ includes only SYSTEM's -- no allocation, no I/O,
  no threads (CLAUDE.md: the rule every other one serves)
- unresolved quoted includes, headers nothing includes, a .c whose first
  include isn't its own header (the header proves it stands alone)

An #include "x.h" resolves like the compiler with the Makefile's flags:
the including file's own directory first, then -Isrc (and -Itests).
"""
import argparse
import collections
import glob
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SOURCES = ["src/*/*.[ch]", "tests/*.[ch]", "bench/*.c"]
INCLUDE = re.compile(r'^\s*#\s*include\s*([<"])([^>"]+)[>"]', re.M)

# module -> the modules it may include (itself always). Lower layers first.
ALLOWED = {
    "src/base": [],
    "src/analysis": ["src/base"],
    "src/synth": ["src/base"],
    # a filter's coefficients and an LFO's sweep: synth's sine, exp2 and
    # pitch ratios (DESIGN.md "The synthesis half")
    "src/fx": ["src/base", "src/synth"],
    "src/offline": ["src/base"],
}

# what src/ may take from the system: libc's pure parts and the vector
# intrinsics (simd.h's only)
SYSTEM = {"stddef.h", "stdint.h", "stdbool.h", "string.h", "math.h",
          "float.h", "limits.h", "immintrin.h"}


def module(path):
    return os.path.dirname(path) if path.startswith("src/") else None


def resolve(path, name):
    for d in (os.path.dirname(path), "src", "tests"):
        p = os.path.normpath(os.path.join(d, name))
        if os.path.exists(os.path.join(ROOT, p)):
            return p
    return None


def scan():
    files = sorted({f for pat in SOURCES
                    for f in glob.glob(pat, root_dir=ROOT)})
    edges, problems = [], []
    for f in files:
        with open(os.path.join(ROOT, f)) as fh:
            text = fh.read()
        quoted = []
        for kind, name in INCLUDE.findall(text):
            if kind == "<":
                if f.startswith("src/") and name not in SYSTEM:
                    problems.append(f"{f}: <{name}> isn't one of src/'s "
                                    "system headers (SYSTEM)")
                continue
            to = resolve(f, name)
            if to is None:
                problems.append(f'{f}: "{name}" doesn\'t resolve')
                continue
            quoted.append(to)
            edges.append((f, to))
        if f.startswith("src/") and f.endswith(".c"):
            own = f[:-2] + ".h"
            if os.path.exists(os.path.join(ROOT, own)) and \
                    (not quoted or quoted[0] != own):
                problems.append(f"{f}: its first include isn't {own}")
    return files, edges, problems


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--files", action="store_true")
    args = ap.parse_args()
    files, edges, problems = scan()
    if args.files:
        for a, b in edges:
            print(f"{a} -> {b}")
        return 0
    graph = collections.defaultdict(collections.Counter)
    for a, b in edges:
        ma, mb = module(a), module(b)
        if ma and mb and ma != mb:
            graph[ma][mb] += 1
            if mb not in ALLOWED.get(ma, []):
                problems.append(f"{a} -> {b}: {ma} may not include {mb} "
                                "(ALLOWED)")
    for m in sorted({module(f) for f in files} - {None}):
        if m not in ALLOWED:
            problems.append(f"{m}: not in ALLOWED")
        deps = ", ".join(f"{d} ({n})" for d, n in sorted(graph[m].items()))
        print(f"{m}: {deps or '-'}")
    used = {b for _, b in edges}
    for f in files:
        if f.startswith("src/") and f.endswith(".h") and f not in used:
            problems.append(f"{f}: nothing includes it")
    for p in problems:
        print(p)
    print(f"deps: {len(files)} files, {len(edges)} includes, "
          f"{len(problems)} findings")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
