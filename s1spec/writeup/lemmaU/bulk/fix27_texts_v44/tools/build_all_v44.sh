#!/bin/bash
# build_all_v44.sh -- (v44: follows fix24_texts_v43/tools/build_all_v43.sh with v44 names; the macro check (now also seeing the wrapper macros
# \fn, \key, \pth, LaTeX compile fix26 defect 1) must PASS on v44 and FAIL on v43 (\key) and on v39; new title check (fix26 defect 2) must PASS
# on v44 and FAIL on v43.)
# Regenerates the whole text set v44 and runs the checks (one process at a time, nice 19, python -B).
# Order: records (from disk) -> numbers macros -> EI constants -> texts (build_v44.py from fix24_texts_v43/, read-only; copies fig_curves.pdf) ->
# balance -> macros -> titles -> numbers check -> SHA256SUMS.
# Reads fix24_texts_v43/ and the step-4 / step-3 records read-only; writes only into fix27_texts_v44/.
set -e
cd /home/claude/s1spec/writeup/lemmaU/bulk/fix27_texts_v44
export PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1
P="nice -n 19 python3 -B"
$P tools/gen_records_v44.py > tools/gen_records_v44.log
$P tools/gen_numbers_A_v44.py . > tools/gen_numbers_A_v44.log
$P tools/gen_ei_constants_v44.py . > tools/gen_ei_constants_v44.log
$P tools/build_v44.py > tools/build_v44.log
$P tools/check_balance_v44.py > tools/check_balance_v44.log; echo "balance: $(grep -c '^PASS' tools/check_balance_v44.log) PASS, $(grep -c '^FAIL' tools/check_balance_v44.log) FAIL"
$P tools/check_macros_v44.py . v44 tools/macro_check_v44.json > tools/check_macros_v44.log; tail -1 tools/check_macros_v44.log
# the same macro check must FAIL on the v39 set (review fix18 F-2) and on the v43 set (\key, LaTeX compile fix26 defect 1); both read-only
$P tools/check_macros_v44.py ../fix17_texts_v39 v39 tools/macro_check_on_v39.json > tools/macro_check_on_v39.log
grep -q '^MACROS FAIL' tools/macro_check_on_v39.log && echo "macros on text v39: FAIL as expected ($(grep -c '^FAIL' tools/macro_check_on_v39.log) FAIL lines)" || { echo 'macro check does not detect the v39 defects'; exit 1; }
$P tools/check_macros_v44.py ../fix24_texts_v43 v43 tools/macro_check_on_v43.json > tools/macro_check_on_v43.log
[ "$(grep '^FAIL' tools/macro_check_on_v43.log | sort)" = "$(printf '%s\n' "FAIL macros GS_lemma_v43.tex {'undefined_project_macros': ['key']}" "FAIL macros erratum_v43.tex {'undefined_project_macros': ['key']}" "FAIL macros flux2_lemma_v43.tex {'undefined_project_macros': ['key']}" | sort)" ] \
  && echo "macros on text v43: FAIL as expected (\\key in GS_lemma, erratum, flux2_lemma)" || { echo 'macro check does not detect exactly the v43 \key defects'; exit 1; }
$P tools/check_titles_v44.py . tools/title_check_v44.json > tools/check_titles_v44.log; tail -1 tools/check_titles_v44.log
$P tools/check_titles_v44.py ../fix24_texts_v43 tools/title_check_on_v43.json > tools/title_check_on_v43.log
grep -q '^TITLES FAIL' tools/title_check_on_v43.log && grep -q 'step3_v43_main.tex 73 section' tools/title_check_on_v43.log && echo "titles on text v43: FAIL as expected ($(grep -c '^FAIL' tools/title_check_on_v43.log) FAIL lines)" || { echo 'title check does not detect the v43 defect'; exit 1; }
$P tools/check_numbers_v44.py > tools/check_numbers_v44.log; head -1 tools/check_numbers_v44.log | cut -c1-200; grep '^v44:' tools/check_numbers_v44.log
grep -q '^FAIL' tools/check_balance_v44.log && { echo 'BALANCE FAIL'; exit 1; }
grep -q '^MACROS PASS' tools/check_macros_v44.log || { echo 'MACROS FAIL'; exit 1; }
grep -q '^TITLES PASS' tools/check_titles_v44.log || { echo 'TITLES FAIL'; exit 1; }
head -1 tools/check_numbers_v44.log | grep -q '^PASS' || { echo 'NUMBERS FAIL'; exit 1; }
find . -maxdepth 1 -type f \( -name '*.tex' -o -name '*.md' -o -name '*.json' -o -name '*.pdf' \) ! -name 'TEXTS_LOG_v44.md' ! -name 'TEXTS_REPORT_v44.md' ! -name 'LATEX_v44.md' | sort | xargs sha256sum > SHA256SUMS
echo "SHA256SUMS: $(wc -l < SHA256SUMS) files"
