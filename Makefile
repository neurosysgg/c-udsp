CC = cc
CFLAGS = -Wall -Wextra -std=c11 -Isrc -g

SRC = $(wildcard src/*.c)
HDR = $(wildcard src/*.h)

TEST_SRC = $(wildcard tests/*.c)
TEST_BIN = $(patsubst tests/%.c,test-%,$(TEST_SRC))

test-%: tests/%.c tests/udsp_test.h $(SRC) $(HDR)
	$(CC) $(CFLAGS) -Itests -o $@ tests/$*.c $(SRC) -lm

test: $(TEST_BIN)
	@for t in $(TEST_BIN); do echo "-- $$t --"; ./$$t || exit 1; done

# the suites under every build the kernels have a path for: gcc and clang,
# -O2/-O3, baseline (SSE2) and x86-64-v3 (AVX2+FMA), and the plain-C one
# (UDSP_SCALAR, the reference the vector paths are held to), and ASan+UBSan;
# bench/ulp.c's tolerances too, in each
VARIANTS = gcc:-O2 gcc:-O3 gcc:-O2,-march=x86-64-v3 gcc:-O3,-march=x86-64-v3 \
	gcc:-O2,-DUDSP_SCALAR clang:-O2 clang:-O3,-march=x86-64-v3 \
	clang:-O2,-DUDSP_SCALAR \
	gcc:-O1,-march=x86-64-v3,-fsanitize=address,-fsanitize=undefined,-fno-sanitize-recover=all \
	gcc:-O1,-fsanitize=address,-fsanitize=undefined,-fno-sanitize-recover=all
test-variants:
	@mkdir -p bench/build
	@for v in $(VARIANTS); do \
	  cc=$${v%%:*}; flags=$$(echo $${v#*:} | tr , ' '); fail=0; \
	  for t in $(TEST_SRC); do \
	    $$cc $(CFLAGS) $$flags -Itests -o bench/build/variant $$t $(SRC) -lm \
	      || exit 1; \
	    ./bench/build/variant > bench/build/variant.out || { fail=1; \
	      grep FAIL bench/build/variant.out; }; \
	  done; \
	  $$cc $(CFLAGS) $$flags -o bench/build/variant bench/ulp.c $(SRC) -lm \
	    || exit 1; \
	  ./bench/build/variant > bench/build/variant.out || { fail=1; \
	    awk -F'\t' '$$3 > $$5 { print "ULP over tolerance:", $$0 }' \
	      bench/build/variant.out; }; \
	  printf '%-6s %-32s %s\n' $$cc "$$flags" $$([ $$fail = 0 ] && echo ok || echo FAILED); \
	  [ $$fail = 0 ] || exit 1; \
	done

bench:
	bench/matrix.py

# every kernel's error in ULPs against a double reference, every build
# (bench/ulp.c's tolerances; bench/README.md "ULPs")
ulp:
	bench/ulp.py

clean:
	rm -f $(TEST_BIN)

.DEFAULT_GOAL := test
.PHONY: test test-variants bench ulp clean
