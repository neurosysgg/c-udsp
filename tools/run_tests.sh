#!/bin/sh
# Runs test binaries JOBS at a time (`make test`; TEST_JOBS=1 for one at a
# time). Each one's output goes to build/test/NAME.out; a line per suite
# (its last: "N passed, M failed", its seconds), a failed suite's whole
# output, then the totals. Exits 1 when any suite failed or crashed.
# Taken from ctui-wm's tools/run_tests.sh.
#
#   tools/run_tests.sh JOBS BINARY...
#
# TEST_WRAP, when set, goes before each binary (`make test-valgrind`:
# valgrind)
jobs=$1
shift
out=build/test
rm -rf "$out"
mkdir -p "$out"
for t in "$@"; do echo "$t"; done |
	xargs -P"$jobs" -I{} sh -c \
		'n=$(basename {}); s=$(date +%s%N); $TEST_WRAP ./{} > '"$out"'/$n.out 2>&1
		echo $? $((($(date +%s%N) - s) / 1000000)) > '"$out"'/$n.rc'
pass=0 fail=0 bad=
for t in "$@"; do
	n=$(basename "$t")
	read -r rc ms < "$out/$n.rc" 2>/dev/null || rc="?" ms=0
	# udsp_test.h's last line: "-- N passed, M failed --"
	last=$(tail -n 1 "$out/$n.out")
	last=${last#-- }
	last=${last% --}
	case $last in
	*" passed, "*" failed")
		pass=$((pass + ${last%% passed*}))
		f=${last#* passed, }
		fail=$((fail + ${f%% failed*}))
		;;
	esac
	if [ "$rc" = 0 ]; then
		printf '%-24s %-22s %3d.%01d s\n' "${n#test-}" "$last" \
			$((ms / 1000)) $((ms % 1000 / 100))
	else
		bad="$bad ${n#test-}"
		printf -- '-- %s: exit %s --\n' "$n" "$rc"
		cat "$out/$n.out"
	fi
done
echo "$# suites, $pass passed, $fail failed${bad:+; failed:$bad} (outputs in $out/)"
[ -z "$bad" ]
