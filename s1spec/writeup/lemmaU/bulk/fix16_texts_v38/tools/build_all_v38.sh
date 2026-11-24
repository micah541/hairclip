#!/bin/bash
# build_all_v38.sh -- regenerates the whole text set v38 and runs the checks (one process at a time, nice 19, python -B).
# Order: numbers -> constants -> interface texts -> step-3 texts -> main texts -> erratum -> balance -> numbers check.
set -e
cd /home/claude/s1spec/writeup/lemmaU/bulk/fix16_texts_v38
export PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1
P="nice -n 19 python3 -B"
$P tools/gen_numbers_A_v38.py . > tools/gen_numbers_A_v38.log
$P tools/gen_ei_constants_v38.py . > tools/gen_ei_constants_v38.log
$P tools/build_interface_v38.py
$P tools/build_texts_v38.py
$P tools/build_main_v38.py
$P tools/build_erratum_v38.py
$P tools/check_balance_v38.py > tools/check_balance_v38.log; grep -c '^PASS' tools/check_balance_v38.log
$P tools/check_numbers_v38.py > tools/check_numbers_v38.log; head -1 tools/check_numbers_v38.log
