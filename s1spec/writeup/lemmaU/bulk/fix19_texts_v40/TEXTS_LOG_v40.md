# TEXTS_LOG_v40 (resumable log of text set v40)

Dir: `/home/claude/s1spec/writeup/lemmaU/bulk/fix19_texts_v40/`. Source (read only): `fix17_texts_v39/`. Task: apply `fix18_final_review/FINAL_REVIEW_v39.md` F-1..F-7, R-1..R-3.
Marker for the "nothing outside changed" check: scratchpad file `start_marker` touched 2026-11-21 06:34 PST, before the copy.

## Log
- 06:35 copied fix17 -> fix19 (cp -a); then removed the v39 copies (texts and tools) from fix19, so that fix19 holds only v40 files and no v39 tool (which would write into fix17) can be run from here. fix17 is untouched.
- Plan: tools/b40_common.py (retag v39 -> v40 names), gen_records_v40.py, gen_numbers_A_v40.py (F-2 macros), gen_ei_constants_v40.py, build_v40.py (reads fix17 texts, asserted replacements), check_balance_v40.py, check_macros_v40.py (new: used subset of defined, valid names), check_numbers_v40.py, build_all_v40.sh.
- F-1 route: option (i) of the reviewer (B1(a,b) re-derived by first touch in the Open step; no new number). 
- 06:40-07:30 tools written: b40_common.py, gen_records_v40.py (+ distinct-file counts), gen_numbers_A_v40.py (F-2), gen_ei_constants_v40.py, check_macros_v40.py (new), check_balance_v40.py, check_numbers_v40.py, build_v40.py + edits_v40.txt (65 asserted edits as plain text; a first version with Python string literals was discarded because LaTeX '' quotes broke the raw strings -- nothing was written from it except files that were then rebuilt), erratum_part5_v40.tex, records_errata_add_v40.md, build_all_v40.sh.
- Macro check on v39: FAIL (4 invalid names, 8 documents). Found one more instance than the review: the standalone preamble of flux2_lemma used 24 number macros without \input of the numbers file (fixed).
- Number-check fixes: 380.25 (win_arb.py range constant), 0.0466 (MOD_IBP.md [read]), 0.071446 (erratum history) whitelisted with reasons; "(y<=57.75)" dropped from a new sentence (direction check); cubic_E0fine named by its real file.
- Final: build_all_v40.sh: balance 23 PASS / 0 FAIL; macros PASS (v39: FAIL as expected); numbers PASS; SHA256SUMS 24 files; two runs byte-identical (56 files); about 5 s, 61 MB peak. find -newer: nothing outside fix19 changed; fix17 SHA256SUMS verifies.
- DONE. Report: TEXTS_REPORT_v40.md.
