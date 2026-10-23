#!/usr/bin/env python3
"""Static checks beyond the build's own warnings, and clangd's database.

  tools/lint.py compdb --cc CC < make-dry-run  > build/lint/compile_commands.json
  tools/lint.py [--update] [--gcc] [--ctu] [--all] [PATH...]

Taken from ctui-wm's tools/lint.py (by way of ctui-mus's). compdb
turns `make -nB`'s compile lines into a compilation database (`make lint`
makes it in build/lint/; `make compile_commands.json` copies it to the
top for clangd): every source with the exact flags its target builds it
with, the first line naming a source winning.

The default run (`make lint`) compiles every source (src/, tests/,
bench/) with clang twice: its warnings (a second compiler's view of the
build's, plus STRICT's; --analyze leaves most of them out) and its static
analyzer, then clang-tidy (bugprone-* and cert-*, minus TIDY_QUIET's
policy checks) and cppcheck (--enable=warning,portability, minus
CPPCHECK_QUIET's policy kinds); either is skipped with a note when it
isn't installed. A finding is keyed by file, message and checker, not by
line, so edits elsewhere don't move it; the keys already triaged sit in
tools/lint-baseline.txt, each under a comment saying why it stays. Only
findings beyond the baseline fail the run; baseline keys no longer found
are listed so the file can shrink (--update rewrites it from this run,
keeping the comments of the keys that remain).

PATH... limits the run to sources under those paths (and skips the
baseline's "gone" list, which would be meaningless).

--all lists every finding with its lines, the baseline's too (working
through the baseline).

--gcc also runs gcc -O2 -fanalyzer and prints its findings grouped by
kind, report only: it mostly flags unchecked malloc/realloc (a policy
here, not a bug) and doesn't gate anything.

--ctu also runs the analyzer across files (clang's cross-TU inlining,
the sources parsed on demand) and lists what it finds beyond the one-file run, with
lines; report only. `clang --analyze -Xclang -analyzer-output=text` with the same
-analyzer-config (printed) shows a finding's path.
"""
import argparse
import collections
import concurrent.futures
import json
import os
import re
import shlex
import shutil
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASELINE = os.path.join(ROOT, "tools", "lint-baseline.txt")
FOREIGN = ()
# where `make lint` writes its database
LINTDB = os.path.join(ROOT, "build", "lint")
# flags a compile line keeps for the database (the rest: -o, -MMD, ...)
KEEP = re.compile(r"-(I|D|U|W(?!l,)|std=|O|g|f|m|pthread|include)")
# gcc's own warnings (Makefile WARN): the database is for clang's tools
# (clangd, lint, clang-tidy), which call them unknown options
GCC_ONLY = ("-Wduplicated-cond", "-Wduplicated-branches")
# clang/gcc: path:line:col: warning|error: message [checker or -Wflag]
# (an error counts too: a source that no longer builds, e.g. a fuzzer)
DIAG = re.compile(r"^(.+?):(\d+):(\d+): (?:warning|error): (.*?)(?: \[([^\]]+)\])?$")


def compdb(cc, lines):
    out, seen = [], set()
    for line in lines:
        try:
            argv = shlex.split(line)
        except ValueError:
            continue
        if not argv or argv[0] != cc:
            continue
        flags, srcs, skip = [], [], False
        for a in argv[1:]:
            if skip:
                skip = False
            elif a == "-o":
                skip = True
            elif a.endswith(".c"):
                srcs.append(a)
            elif KEEP.match(a) and a not in GCC_ONLY:
                flags.append(a)
        for src in srcs:
            if src in seen:
                continue
            seen.add(src)
            out.append({"directory": ROOT, "file": src,
                        "arguments": [cc] + flags + ["-c", src]})
    return out


def ours(path):
    rel = os.path.relpath(os.path.join(ROOT, path), ROOT)
    return not rel.startswith("..") and not rel.startswith(FOREIGN)


def run(argv):
    p = subprocess.run(argv, cwd=ROOT, capture_output=True, text=True)
    return p.stdout + p.stderr  # clang-tidy's findings come on stdout


# warnings beside the build's: a local hiding another, a function that's
# neither static nor declared in a header
STRICT = ["-Wshadow", "-Wmissing-prototypes"]


def findings(entries, compiler, extra):
    """{key: count}, key = "file: message [checker]", each (file, line,
    key) once: a header's finding comes back from every source that
    includes it."""
    return diags([[compiler] + e["arguments"][1:-1] + extra +
                  ["-o", "/dev/null", e["file"]] for e in entries])


# clang-tidy's checks that only restate a policy here, or flag an idiom
# the code uses on purpose (session 189 triaged all of bugprone-* and
# cert-*: ~10000 findings, one real -- a conditional with equal arms):
TIDY_QUIET = (
    # return values left unchecked where nothing can be done (fclose of
    # a read, snprintf into a sized buffer, ...); 5000 of them
    "cert-err33-c",
    # flags as bit sets in ints; parameters of one type side by side;
    # `if ((x = f()))` and `!strcmp`; `_NAME` include guards and
    # feature macros; free((void *)p)
    "bugprone-signed-bitwise", "bugprone-easily-swappable-parameters",
    "bugprone-assignment-in-if-condition",
    "bugprone-assignment-in-selection-statement", "cert-exp45-c",
    "bugprone-suspicious-string-compare",
    "bugprone-reserved-identifier", "cert-dcl37-c", "cert-dcl51-cpp",
    "bugprone-multi-level-implicit-pointer-conversion",
    # allocations aren't checked for NULL (cppcheck's memleakOnRealloc
    # likewise); atoi/strtol where a bad number reads as 0 on purpose
    "bugprone-suspicious-realloc-usage",
    "bugprone-unchecked-string-to-number-conversion", "cert-err34-c",
    # int * int for sizes of things already bounded (pictures by their
    # MAX_AREA); int -> float/char narrowing in drawing maths
    "bugprone-implicit-widening-of-multiplication-result",
    "bugprone-narrowing-conversions",
    # the configured shell commands (shot's select/capture, gsettings,
    # timedatectl) run through the shell by design
    "bugprone-command-processor", "cert-env33-c",
    # memcpy into length-counted byte buffers, not strings
    "bugprone-not-null-terminated-result",
    # style: switches without default, enums partly numbered, `1ll`,
    # else-if chains with equal bodies kept apart for reading, macros
    # that paste a field name, `--n` in a condition, signed char as a
    # small number, rand() as a getentropy fallback
    "bugprone-switch-missing-default-case", "bugprone-unhandled-code-paths",
    "cert-int09-c", "cert-dcl16-c", "bugprone-branch-clone",
    "bugprone-macro-parentheses", "bugprone-inc-dec-in-conditions",
    "bugprone-signed-char-misuse", "cert-str34-c",
    "cert-msc30-c", "cert-msc50-cpp",
)


def tidy(entries):
    """clang-tidy over entries: {key: count}, {key: [lines]}"""
    checks = ",".join(["-*", "bugprone-*", "cert-*"] +
                      ["-" + c for c in TIDY_QUIET])
    return diags([["clang-tidy", "-p", LINTDB, "--quiet",
                   f"--checks={checks}", e["file"]] for e in entries])


def diags(jobs):
    """the compiler-style diagnostics of jobs (argv lists, run in
    parallel), keyed as findings() says"""
    seen = set()
    with concurrent.futures.ThreadPoolExecutor(os.cpu_count()) as pool:
        for err in pool.map(run, jobs):
            for line in err.splitlines():
                m = DIAG.match(line)
                if not m or not ours(m[1]):
                    continue
                path = os.path.relpath(os.path.join(ROOT, m[1]), ROOT)
                key = f"{path}: {m[4]} [{m[5] or '?'}]"
                seen.add((path, m[2], key))
    counts = collections.Counter(k for _, _, k in seen)
    lines = collections.defaultdict(list)
    for path, ln, key in seen:
        lines[key].append(int(ln))
    return counts, lines


def ctu(entries, plain):
    """the analyzer with cross-TU inlining over entries; prints the
    findings plain (the one-file run's {key: count}) doesn't have"""
    work = os.path.join(ROOT, "build", "lint", "ctu")
    os.makedirs(work, exist_ok=True)
    with open(os.path.join(LINTDB, "compile_commands.json")) as f:
        db = json.load(f)
    # where each external function is defined; one defined twice (each
    # app's main, a test's stub) is left out: CTU can't pick
    defs = collections.defaultdict(set)
    jobs = [["clang-extdef-mapping", "-p", LINTDB, e["file"]] for e in db]
    with concurrent.futures.ThreadPoolExecutor(os.cpu_count()) as pool:
        for out in pool.map(run, jobs):
            for line in out.splitlines():
                usr, _, path = line.rpartition(" ")
                if re.match(r"\d+:c:", usr):  # "LEN:USR PATH"; not a warning
                    defs[usr].add(path)
    with open(os.path.join(work, "externalDefMap.txt"), "w") as f:
        for usr, paths in sorted(defs.items()):
            if len(paths) == 1:
                f.write(f"{usr} {next(iter(paths))}\n")
    # how to parse each on demand: clang by its full path (its own
    # headers are found from there)
    clang = shutil.which("clang")
    with open(os.path.join(work, "invocations.yaml"), "w") as f:
        for e in db:
            src = os.path.join(ROOT, e["file"])
            argv = [clang] + e["arguments"][1:-2] + ["-c", src]
            f.write(json.dumps(src) + ": [" +
                    ", ".join(json.dumps(a) for a in argv) + "]\n")
    cfg = ("experimental-enable-naive-ctu-analysis=true,ctu-dir=" + work +
           ",ctu-index-name=externalDefMap.txt,ctu-invocation-list=" +
           os.path.join(work, "invocations.yaml"))
    counts, where = findings(entries, "clang", ["--analyze", "-Xclang",
                             "-analyzer-config", "-Xclang", cfg])
    more = counts - plain
    for key in sorted(more):
        path, rest = key.split(": ", 1)
        lns = ",".join(map(str, sorted(where[key])))
        print(f"ctu: {path}:{lns}: {rest}")
    print(f"ctu: {sum(more.values())} beyond the one-file run "
          f"(-analyzer-config {cfg})")


def load_baseline():
    """[(comments, key)] in file order; a key may repeat (one line per
    occurrence)."""
    out, notes, head = [], [], True
    if not os.path.exists(BASELINE):
        return out
    with open(BASELINE) as f:
        for line in f:
            line = line.rstrip("\n")
            if head:  # the file's own comment, up to the first blank line
                head = line.startswith("#")
                if head:
                    continue
            if line.startswith("#") or not line.strip():
                notes.append(line)
            else:
                out.append((notes, line))
                notes = []
    return out


def write_baseline(counts, old):
    notes = {}
    for n, key in old:
        if n and key not in notes:
            notes[key] = n
    head = []
    with open(BASELINE) as f:
        for line in f:
            if not line.startswith("#"):
                break
            head.append(line.rstrip("\n"))
    with open(BASELINE, "w") as f:
        f.write("\n".join(head) + "\n")
        for key in sorted(counts):
            n = notes.get(key)
            if n:
                f.write("\n" + "\n".join(x for x in n if x.strip()) + "\n")
            f.write((key + "\n") * counts[key])


# cppcheck's kinds that only restate a policy here: allocations aren't
# checked for NULL (out of memory ends the process soon enough), and NULL
# ends variadic lists of library calls as their docs show; (void)a, (void)b
CPPCHECK_QUIET = ("nullPointerOutOfMemory", "nullPointerOutOfResources",
                  "nullPointerArithmeticOutOfMemory",
                  "ctunullpointerOutOfMemory", "memleakOnRealloc",
                  "varFuncNullUB", "constStatement")


def cppcheck(entries):
    """cppcheck over entries: {key: count}, {key: [lines]} like
    findings(); its own cache in build/lint/ (only what changed reruns)"""
    work = os.path.join(ROOT, "build", "lint")
    os.makedirs(os.path.join(work, "cppcheck"), exist_ok=True)
    db = os.path.join(work, "cppcheck.json")
    with open(db, "w") as f:
        json.dump(entries, f)
    argv = ["cppcheck", f"--project={db}", f"-j{os.cpu_count()}",
            "--enable=warning,portability", "--inline-suppr", "--quiet",
            f"--cppcheck-build-dir={os.path.join(work, 'cppcheck')}",
            "--template={file}\t{line}\t{message}\t{id}"]
    argv += [f"--suppress={k}" for k in CPPCHECK_QUIET]
    argv += [f"--suppress=*:{p}*" for p in FOREIGN + ("gen/", "/usr/")]
    err = run(argv)
    seen = set()
    for line in err.splitlines():
        parts = line.split("\t")
        if len(parts) != 4 or not ours(parts[0]):
            continue
        path = os.path.relpath(os.path.join(ROOT, parts[0]), ROOT)
        seen.add((path, parts[1], f"{path}: {parts[2]} [cppcheck:{parts[3]}]"))
    counts = collections.Counter(k for _, _, k in seen)
    lines = collections.defaultdict(list)
    for path, ln, key in seen:
        lines[key].append(int(ln))
    return counts, lines


def lint(args):
    db = os.path.join(LINTDB, "compile_commands.json")
    if not os.path.exists(db):
        sys.exit("no build/lint/compile_commands.json: `make lint` makes it")
    with open(db) as f:
        entries = [e for e in json.load(f) if ours(e["file"])]
    if args.paths:
        entries = [e for e in entries
                   if any(e["file"].startswith(p.rstrip("/"))
                          for p in args.paths)]
    counts, where = findings(entries, "clang", ["-fsyntax-only"] + STRICT)
    more, more_where = findings(entries, "clang", ["--analyze"])
    analyzed = more
    counts += more
    for key, lns in more_where.items():
        where[key] += lns
    if shutil.which("clang-tidy"):
        more, more_where = tidy(entries)
        counts += more
        for key, lns in more_where.items():
            where[key] += lns
    else:
        print("clang-tidy isn't installed: its findings left out")
    if shutil.which("cppcheck"):
        more, more_where = cppcheck(entries)
        counts += more
        for key, lns in more_where.items():
            where[key] += lns
    else:
        print("cppcheck isn't installed: its findings left out")
    old = load_baseline()
    base = collections.Counter(k for _, k in old)
    new = counts - base
    if args.all:
        for key in sorted(counts):
            path, rest = key.split(": ", 1)
            lns = ",".join(map(str, sorted(where[key])))
            print(f"{path}:{lns}: {rest}")
    for key in sorted(new):
        path, rest = key.split(": ", 1)
        lns = ",".join(map(str, sorted(where[key])))
        print(f"new: {path}:{lns}: {rest}" +
              (f" (x{new[key]} beyond the baseline)"
               if base[key] else ""))
    if not args.paths:
        for key in sorted(base - counts):
            print(f"gone (drop from the baseline): {key}")
    print(f"clang + clang-tidy + cppcheck: {len(entries)} sources, "
          f"{sum(counts.values())} findings, {sum(new.values())} new")
    if args.update:
        write_baseline(counts, old)
        print(f"baseline rewritten: {os.path.relpath(BASELINE)}")
    if args.ctu:
        ctu(entries, analyzed)
    if args.gcc:
        g, _ = findings(entries, "gcc",
                        ["-O2", "-fanalyzer", "-Wno-format-truncation"])
        kinds = collections.Counter()
        for key, n in g.items():
            kinds[key.rsplit("[", 1)[1].rstrip("]")] += n
        for k, n in kinds.most_common():
            print(f"gcc {k}: {n}")
    return 1 if new and not args.update else 0


def main():
    if sys.argv[1:2] == ["compdb"]:
        ap = argparse.ArgumentParser(prog="lint.py compdb")
        ap.add_argument("--cc", default="cc")
        args = ap.parse_args(sys.argv[2:])
        json.dump(compdb(args.cc, sys.stdin), sys.stdout, indent=1)
        print()
        return 0
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--update", action="store_true")
    ap.add_argument("--gcc", action="store_true")
    ap.add_argument("--ctu", action="store_true")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("paths", nargs="*")
    return lint(ap.parse_args())


if __name__ == "__main__":
    sys.exit(main())
