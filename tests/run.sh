#!/usr/bin/env bash
set -u
cd "$(dirname "$0")"
ROOT=..
DIR=build/cases
BIN=build/driver
FEATURES="conv inv ln exp compose high spsinv spsexp matmul matinv matdet badinv badln badexp badmat prop"

usage() {
    echo "usage: $0 [options] [feature...]"
    echo "  -l, --list       list features"
    echo "  -r, --rebuild    rebuild driver"
    echo "  -q, --quiet      only failures"
    echo "  -a, --all        all features (default)"
    echo "  -g, --gen        regenerate data"
}

list() {
    for f in $FEATURES; do
        printf '%-10s %s cases\n' "$f" "$(find "$DIR/$f" -name '*.in' 2>/dev/null | wc -l)"
    done
}

rebuild=0
quiet=0
gen=0
sel=""
for a in "$@"; do
    case "$a" in
        -l|--list) list; exit 0 ;;
        -h|--help) usage; exit 0 ;;
        -r|--rebuild) rebuild=1 ;;
        -q|--quiet) quiet=1 ;;
        -a|--all) sel="$FEATURES" ;;
        -g|--gen) gen=1 ;;
        -*) echo "unknown option: $a"; usage; exit 2 ;;
        *) sel="$sel $a" ;;
    esac
done
[ -n "$sel" ] || sel="$FEATURES"

if [ "$gen" = 1 ] || [ ! -d "$DIR" ]; then
    python3 gen.py
fi
if [ "$rebuild" = 1 ] || [ ! -x "$BIN" ] || [ "$ROOT/poly.cpp" -nt "$BIN" ] || [ "driver.cpp" -nt "$BIN" ]; then
    mkdir -p build
    if ! g++ -O2 -std=c++17 -o "$BIN" driver.cpp; then
        echo "compile failed"
        exit 1
    fi
fi

pass=0
fail=0
for f in $sel; do
    if [ "$f" = "prop" ]; then
        for op in sqrt trig trans div bit subset interp shift; do
            if printf '%s 64 20261115\n' "$op" | "$BIN" prop >/dev/null; then
                pass=$((pass+1))
                [ "$quiet" = 1 ] || echo "PASS prop/$op"
            else
                fail=$((fail+1))
                echo "FAIL prop/$op"
            fi
        done
        continue
    fi
    [ -d "$DIR/$f" ] || { echo "unknown feature: $f"; fail=$((fail+1)); continue; }
    for input in "$DIR/$f"/*.in; do
        [ -e "$input" ] || continue
        name=$(basename "$input" .in)
        actual="build/$f.$name.out"
        if ! "$BIN" "$f" < "$input" > "$actual"; then
            if [ -e "$DIR/$f/$name.err" ]; then
                pass=$((pass+1))
                [ "$quiet" = 1 ] || echo "PASS $f/$name (expected error)"
            else
                fail=$((fail+1))
                echo "FAIL $f/$name (runtime)"
            fi
            continue
        fi
        if [ -e "$DIR/$f/$name.err" ]; then
            fail=$((fail+1))
            echo "FAIL $f/$name (expected error, got success)"
            continue
        fi
        if [ ! -e "$DIR/$f/$name.out" ]; then
            fail=$((fail+1))
            echo "FAIL $f/$name (missing expected output)"
            continue
        fi
        if diff -q "$DIR/$f/$name.out" "$actual" >/dev/null; then
            pass=$((pass+1))
            [ "$quiet" = 1 ] || echo "PASS $f/$name"
        else
            fail=$((fail+1))
            echo "FAIL $f/$name"
            diff -u "$DIR/$f/$name.out" "$actual" | head -20
        fi
    done
done
echo "passed $pass, failed $fail"
[ "$fail" = 0 ]
