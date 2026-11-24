#!/bin/bash
# build_all_v39.sh -- regenerates the whole text set v39 and runs the checks (one process at a time, nice 19, python -B).
# Order: records (from disk) -> numbers macros -> EI constants -> texts (incl. step4_handoff_v39 + patch note) -> balance -> numbers check -> SHA256SUMS.
# Reads fix16_texts_v38/ and the step-4 / step-3 records read-only; writes only into fix17_texts_v39/.
set -e
cd /home/claude/s1spec/writeup/lemmaU/bulk/fix17_texts_v39
export PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1
P="nice -n 19 python3 -B"
$P tools/gen_records_v39.py > tools/gen_records_v39.log
$P tools/gen_numbers_A_v39.py . > tools/gen_numbers_A_v39.log
$P tools/gen_ei_constants_v39.py . > tools/gen_ei_constants_v39.log
$P tools/build_v39.py > tools/build_v39.log
$P tools/check_balance_v39.py > tools/check_balance_v39.log; echo "balance: $(grep -c '^PASS' tools/check_balance_v39.log) PASS, $(grep -c '^FAIL' tools/check_balance_v39.log) FAIL"
$P tools/check_numbers_v39.py > tools/check_numbers_v39.log; head -1 tools/check_numbers_v39.log | cut -c1-200; grep '^v39:' tools/check_numbers_v39.log
grep -q '^FAIL' tools/check_balance_v39.log && { echo 'BALANCE FAIL'; exit 1; }
head -1 tools/check_numbers_v39.log | grep -q '^PASS' || { echo 'NUMBERS FAIL'; exit 1; }
find . -maxdepth 1 -type f \( -name '*.tex' -o -name '*.md' -o -name '*.json' \) ! -name 'TEXTS_LOG_v39.md' ! -name 'TEXTS_REPORT_v39.md' | sort | xargs sha256sum > SHA256SUMS
echo "SHA256SUMS: $(wc -l < SHA256SUMS) files"
