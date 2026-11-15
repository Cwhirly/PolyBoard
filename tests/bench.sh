#!/usr/bin/env bash
set -u
cd "$(dirname "$0")"
BIN=build/driver
n=20000
rep=3
ops="conv inv ln exp compose"
while [ $# -gt 0 ]; do
    case "$1" in
        -n|--size) n="$2"; shift 2 ;;
        -r|--reps) rep="$2"; shift 2 ;;
        -h|--help)
            echo "usage: $0 [-n N] [-r REPS] [conv|inv|ln|exp|compose ...]"
            exit 0 ;;
        -*) echo "unknown option: $1"; exit 2 ;;
        *) ops="$*"; break ;;
    esac
done
if [ ! -x "$BIN" ] || [ "../poly.cpp" -nt "$BIN" ] || [ "driver.cpp" -nt "$BIN" ]; then
    mkdir -p build
    g++ -O2 -std=c++17 -o "$BIN" driver.cpp || exit 1
fi
for op in $ops; do
    printf '%s %d 20261115 %d\n' "$op" "$n" "$rep" | "$BIN" bench
done
