CC = cc
OPT ?= -O2
# the library's own flags: what a consumer builds src/ with may differ, so
# nothing in src/ may lean on more than -std=c11 and libm
CFLAGS = -Wall -Wextra -std=c11 -Isrc $(OPT) -g \
  -Wstrict-prototypes -Wold-style-definition -Wpointer-arith -Wundef \
  -Walloca -Wnested-externs
CFLAGS += $(VFLAGS)

# `make SAN=1 test`: AddressSanitizer + UBSan, objects and test binaries
# of their own (build/O2-san/). Any undefined behaviour fails the test.
# (No TSAN: nothing here has a thread.)
SAN ?=
ifneq ($(SAN),)
CFLAGS += -fsanitize=address,undefined,float-cast-overflow -fno-omit-frame-pointer \
  -fno-sanitize-recover=all
export UBSAN_OPTIONS = print_stacktrace=1:halt_on_error=1
endif

# `make COV=1 coverage`: the suites built with clang's source-based
# coverage (build/O2-cov/), run, and what they ran reported per directory
# (tools/coverage.py; HTML: build/cov/html/)
COV ?=
ifneq ($(COV),)
ifneq ($(SAN),)
$(error COV can't be combined with SAN)
endif
CC = clang
CFLAGS += -fprofile-instr-generate -fcoverage-mapping -fprofile-continuous
export LLVM_PROFILE_FILE = $(CURDIR)/build/cov/raw/%c%p-%m.profraw
endif

# gcc's -Wduplicated-* only when CC is gcc (clang, as in lint, doesn't
# know them)
ifneq ($(shell $(CC) --version 2>/dev/null | grep -c 'Free Software'),0)
CFLAGS += -Wduplicated-cond -Wduplicated-branches
endif

# objects, the library and the test binaries per flavour: only what
# changed recompiles, and switching flavours rebuilds nothing of another
# (test-variants sets B itself)
B ?= build/$(subst -,,$(OPT))$(if $(SAN),-san)$(if $(COV),-cov)

SRC = $(wildcard src/*/*.c)
OBJ = $(patsubst %.c,$(B)/obj/%.o,$(SRC))
LIB = $(B)/libudsp.a

TEST_SRC = $(wildcard tests/*.c)
TEST_BIN = $(patsubst tests/%.c,$(B)/bin/test-%,$(TEST_SRC))
ULP_BIN = $(B)/bin/ulp

all: $(LIB) $(TEST_BIN)

$(B)/obj/%.o: %.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -MMD -MP -c -o $@ $<

$(LIB): $(OBJ)
	@rm -f $@
	ar rcs $@ $^

$(B)/bin/test-%: $(B)/obj/tests/%.o $(LIB)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -o $@ $< $(LIB) -lm

$(ULP_BIN): $(B)/obj/bench/ulp.o $(LIB)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -o $@ $< $(LIB) -lm

# T=NAME: only the suites whose name has NAME in it. Side by side, a line
# each, a failed one's whole output (tools/run_tests.sh); every suite's
# output in build/test/. TEST_JOBS=1: one after the other
TEST_RUN = $(if $(T),$(foreach t,$(TEST_BIN),$(if $(findstring $(T),$(notdir $(t))),$(t))),$(TEST_BIN))
TEST_JOBS ?= $(shell nproc)
test: $(TEST_RUN)
	@tools/run_tests.sh $(TEST_JOBS) $(TEST_RUN)

# the suites under memcheck: uninitialised reads, which ASan doesn't see
# (a kernel reading a lane of scratch the caller never filled). Arch's
# glibc is stripped: valgrind needs its symbols from debuginfod
# (/etc/debuginfod), as in ctui-wm
test-valgrind: $(TEST_RUN)
	@DEBUGINFOD_URLS="$${DEBUGINFOD_URLS:-$$(cat /etc/debuginfod/*.urls 2>/dev/null | tr '\n' ' ')}" \
	  TEST_WRAP='valgrind -q --error-exitcode=99 --track-origins=yes' \
	  tools/run_tests.sh $(TEST_JOBS) $(TEST_RUN)

# the suites under every build the kernels have a path for: gcc and clang,
# -O2/-O3, baseline (SSE2) and x86-64-v3 (AVX2+FMA), and the plain-C one
# (UDSP_SCALAR, the reference the vector paths are held to), and ASan+UBSan
# with float-cast-overflow (SAN=1's; gcc's -fsanitize=undefined leaves it
# out, and a float cast out of range passed here unseen, session 257);
# bench/ulp.c's tolerances too, in each. Each variant builds in
# build/variant/NAME/, so a second run only recompiles what changed (NAME
# without '=': a .d file's rule naming it would read as an assignment,
# and the headers would stop counting)
VARIANTS = gcc:-O2 gcc:-O3 gcc:-O2,-march=x86-64-v3 gcc:-O3,-march=x86-64-v3 \
	gcc:-O2,-DUDSP_SCALAR clang:-O2 clang:-O3,-march=x86-64-v3 \
	clang:-O2,-DUDSP_SCALAR \
	gcc:-O1,-march=x86-64-v3,-fsanitize=address,-fsanitize=undefined,-fsanitize=float-cast-overflow,-fno-sanitize-recover=all \
	gcc:-O1,-fsanitize=address,-fsanitize=undefined,-fsanitize=float-cast-overflow,-fno-sanitize-recover=all
test-variants:
	@mkdir -p build/variant
	@for v in $(VARIANTS); do \
	  cc=$${v%%:*}; flags=$$(echo $${v#*:} | tr , ' '); \
	  dir=build/variant/$$(echo $$v | tr -c 'a-zA-Z0-9\n' _); \
	  $(MAKE) --no-print-directory -s -j$(TEST_JOBS) CC=$$cc OPT= VFLAGS="$$flags" \
	    SAN= COV= B=$$dir variant-run > $$dir.log 2>&1; rc=$$?; \
	  printf '%-6s %-32s %s\n' $$cc "$$flags" $$([ $$rc = 0 ] && echo ok || echo FAILED); \
	  [ $$rc = 0 ] || { cat $$dir.log; exit 1; }; \
	done

# one variant (test-variants): its suites, then every kernel within its
# ULP tolerance (bench/ulp.c's lines: kernel, ..., ulp, ..., tolerance)
variant-run: $(TEST_BIN) $(ULP_BIN)
	@tools/run_tests.sh $(TEST_JOBS) $(TEST_BIN)
	@./$(ULP_BIN) > $(B)/ulp.out || { \
	  awk -F'\t' '$$3 > $$5 { print "ULP over tolerance:", $$0 }' $(B)/ulp.out; exit 1; }

bench:
	bench/matrix.py

# every kernel's error in ULPs against a double reference, every build
# (bench/ulp.c's tolerances; bench/README.md "ULPs")
ulp:
	bench/ulp.py

# before calling anything done, in one word: every variant (ASan+UBSan
# among them), the plain build and its suites, then lint
check:
	$(MAKE) --no-print-directory test-variants
	$(MAKE) --no-print-directory -j$(TEST_JOBS) SAN= COV= all test
	$(MAKE) --no-print-directory SAN= COV= lint

# the suites' coverage (COV=1 above): T=NAME narrows the run (the report
# then shows what those suites reach); `tools/coverage.py --help` for more
coverage: $(TEST_RUN)
	@[ -n "$(COV)" ] || { echo "coverage builds with COV=1"; exit 1; }
	rm -rf build/cov/raw
	-@tools/run_tests.sh $(TEST_JOBS) $(TEST_RUN) | tail -n 1
	tools/coverage.py --html $(TEST_RUN)

# clang's warnings and static analyzer, clang-tidy and cppcheck over the
# sources, the tests and the bench against tools/lint-baseline.txt (only
# new findings fail; tools/lint.py --help), then the layering between
# src/'s directories (tools/deps.py), clang-format --dry-run, ruff and
# shellcheck -- each skipped when not installed. Its compilation database
# is build/lint/'s, made from this Makefile's own compile lines; `make
# compile_commands.json` copies it to the top for clangd.
FORMAT_SRC = $(wildcard src/*/*.[ch] tests/*.[ch] bench/*.c)
build/lint/compile_commands.json: Makefile $(wildcard src/*/ tests/ bench/)
	@mkdir -p $(dir $@)
	$(MAKE) --no-print-directory -nB SAN= COV= CC=clang all $(B)/obj/bench/bench.o \
		$(B)/obj/bench/ulp.o $(B)/obj/bench/reduce.o | \
		tools/lint.py compdb --cc clang > $@.tmp && mv $@.tmp $@

lint: build/lint/compile_commands.json
	tools/lint.py
	tools/deps.py
	@if command -v clang-format >/dev/null; then \
	  clang-format --dry-run --Werror $(FORMAT_SRC) && echo "clang-format: clean"; \
	else echo "clang-format isn't installed: skipped"; fi
	@if command -v ruff >/dev/null; then ruff check -q tools bench && echo "ruff: clean"; \
	else echo "ruff isn't installed: skipped"; fi
	@if command -v shellcheck >/dev/null; then shellcheck -S warning tools/*.sh && \
	  echo "shellcheck: clean"; else echo "shellcheck isn't installed: skipped"; fi

compile_commands.json: build/lint/compile_commands.json
	cp $< $@

format:
	clang-format -i $(FORMAT_SRC)

clean:
	rm -rf build

-include $(OBJ:.o=.d) $(patsubst %.c,$(B)/obj/%.d,$(TEST_SRC) bench/ulp.c)

.SECONDARY:
.DEFAULT_GOAL := test
.PHONY: all test test-valgrind test-variants variant-run bench ulp check coverage lint \
  format clean
