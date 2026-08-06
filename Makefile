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

clean:
	rm -f $(TEST_BIN)

.DEFAULT_GOAL := test
.PHONY: test clean
