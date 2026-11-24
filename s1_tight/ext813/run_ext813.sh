#!/usr/bin/env bash
# Rebuilds and reruns the outward S1 enclosure on sigma in [8, 13] (ext813) and all its checks (~15 s).
# Needs CAPD 1.8.x (capd-config) and python-flint >= 0.9; reads ../in_br2.txt, ../prof/, ../bracket2.json,
# ../jhat8.json, ../{gen_input2,tail2,laplace}.py, /home/claude/s1_cap/tail.py, /home/claude/s1spec/step3/hstail.py.
# One process at a time, nice 10.  Writes only into this directory (and a temporary directory).
set -euo pipefail
cd "$(dirname "$0")"
CAPDCFG=${CAPDCFG:-$(command -v capd-config || echo /home/claude/ext/CAPD/build/bin/capd-config)}
NICE="nice -n 10"
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
echo "== 0. toolchain: s1prof.cpp + this CAPD reproduce the certified prof_tight.txt byte-for-byte"
$NICE g++ -O2 ../prof/s1prof.cpp $($CAPDCFG --cflags --libs) -o "$TMP/s1prof"
$NICE "$TMP/s1prof" ../in_br2.txt 10.249077651519443 51246 > "$TMP/prof.txt"
cmp "$TMP/prof.txt" ../prof/prof_tight.txt && echo "   prof_tight.txt reproduced"
echo "== 1. outward CAPD integration sigma = 8 -> 13 from the certified sigma = 8 set (../in_br2.txt)"
$NICE g++ -O2 s1ext.cpp $($CAPDCFG --cflags --libs) -o s1ext
$NICE ./s1ext ../in_br2.txt 5 25000 25 ext813_points.txt 20 > ext813_prof.txt
grep -q '^END$' ext813_prof.txt
echo "== 2. checks (seed provenance, tiling, widths, sigma = 8 vs prof_tight, tail envelopes, Laplace representation)"
$NICE python3 -I -B ext813_check.py | tee ext813_check.log
