#!/usr/bin/env bash
# Rebuilds the conditioned tight S1 profile prof_tip.txt (sigma in [sigma_tip, 8]) and its checks (~15 s CPU).
# Needs CAPD 1.8.x (capd-config), python-flint >= 0.9, numpy/scipy (diagnosis and float sanity check only).
# Reads ../in_br2.txt, ../prof/prof_tight.txt, ../bracket2.json.  One process at a time, nice 10.
set -euo pipefail
cd "$(dirname "$0")"
CAPDCFG=${CAPDCFG:-$(command -v capd-config || echo /home/claude/ext/CAPD/build/bin/capd-config)}
NICE="nice -n 10"
echo "== 0. diagnosis: linearised width budget at the tip (float; not part of the certificate)"
$NICE python3 -I -B diag_float.py
echo "== 1. CAPD: C1 Poincare map and C1 solution curve over the certified sigma = 8 set, point runs, box run"
$NICE g++ -O2 s1cond.cpp $($CAPDCFG --cflags --libs) -o s1cond
$NICE ./s1cond ../in_br2.txt 10.249077651519443 51246 25 tip_cond.txt > prof_cond.txt
echo "== 2. hull tightening (arb) and checks"
$NICE python3 -I -B ext_tip_post.py | tee ext_tip_check.log
rm -f prof_cond.txt
