#!/bin/bash
# build_all_v42.sh -- (v42: copy of fix21_texts_v41/tools/build_all_v41.sh with v42 names; the macro check also runs on v41 (must PASS) and on v39 (must FAIL))
# regenerates the whole text set v42 and runs the checks (one process at a time, nice 19, python -B).
# Order: records (from disk) -> numbers macros -> EI constants -> texts (build_v42.py from fix21_texts_v41/, read-only) -> balance ->
# macros (used subset of defined, valid names; review fix18 F-2) -> numbers check -> SHA256SUMS.
# Reads fix21_texts_v41/ and the step-4 / step-3 records read-only; writes only into fix22_texts_v42/.
set -e
cd /home/claude/s1spec/writeup/lemmaU/bulk/fix22_texts_v42
export PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1
P="nice -n 19 python3 -B"
$P tools/gen_records_v42.py > tools/gen_records_v42.log
$P tools/gen_numbers_A_v42.py . > tools/gen_numbers_A_v42.log
$P tools/gen_ei_constants_v42.py . > tools/gen_ei_constants_v42.log
$P tools/build_v42.py > tools/build_v42.log
$P tools/check_balance_v42.py > tools/check_balance_v42.log; echo "balance: $(grep -c '^PASS' tools/check_balance_v42.log) PASS, $(grep -c '^FAIL' tools/check_balance_v42.log) FAIL"
$P tools/check_macros_v42.py . v42 tools/macro_check_v42.json > tools/check_macros_v42.log; tail -1 tools/check_macros_v42.log
# the same macro check on the v39 set (read-only there) must FAIL: it shows that the check detects review fix18 F-2
$P tools/check_macros_v42.py ../fix17_texts_v39 v39 tools/macro_check_on_v39.json > tools/macro_check_on_v39.log
grep -q '^MACROS FAIL' tools/macro_check_on_v39.log && echo "macros on text v39: FAIL as expected ($(grep -c '^FAIL' tools/macro_check_on_v39.log) FAIL lines)" || { echo 'macro check does not detect the v39 defects'; exit 1; }
# v42: and it must PASS on the v41 set (read-only there), from which v42 is built
$P tools/check_macros_v42.py ../fix21_texts_v41 v41 tools/macro_check_on_v41.json > tools/macro_check_on_v41.log
grep -q '^MACROS PASS' tools/macro_check_on_v41.log && echo "macros on text v41: PASS" || { echo 'macro check fails on the v41 source set'; exit 1; }
$P tools/check_numbers_v42.py > tools/check_numbers_v42.log; head -1 tools/check_numbers_v42.log | cut -c1-200; grep '^v42:' tools/check_numbers_v42.log
grep -q '^FAIL' tools/check_balance_v42.log && { echo 'BALANCE FAIL'; exit 1; }
grep -q '^MACROS PASS' tools/check_macros_v42.log || { echo 'MACROS FAIL'; exit 1; }
head -1 tools/check_numbers_v42.log | grep -q '^PASS' || { echo 'NUMBERS FAIL'; exit 1; }
find . -maxdepth 1 -type f \( -name '*.tex' -o -name '*.md' -o -name '*.json' \) ! -name 'TEXTS_LOG_v42.md' ! -name 'TEXTS_REPORT_v42.md' | sort | xargs sha256sum > SHA256SUMS
echo "SHA256SUMS: $(wc -l < SHA256SUMS) files"
