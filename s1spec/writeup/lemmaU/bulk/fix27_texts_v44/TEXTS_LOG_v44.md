# TEXTS_LOG_v44 (resumable log of text set v44)

Task: text set v44 from text v43 (`fix24_texts_v43/`, read-only): apply (A) the check `fix25_check_v43/CHECK_v43.md` (F-1..F-6, Note; Windows kit runs
`kit3v4_run_win1`, `kit4bv4_run_win1`; kit v5 pending, no hash), (B) the LaTeX defects of `fix26_latex_v43/LATEX_v43.md` (\key, moving-argument titles,
fig_curves.pdf) plus a per-document macro check that sees \key and a title check, (C) records/open lists (v43 compiled in fix26, v44 compile pending).
Rules: new files only in fix27_texts_v44/; nothing outside changes; no number changes; step 3 NOT closed; RAM <= 800 MB, nice 19.

## Resume
- Rerun everything: `bash tools/build_all_v44.sh` (from fix27_texts_v44/). Marker for the "nothing outside changed" check: `.marker_start`.
- Tools were made from fix24's by `tools/copy_tools_v44.py` (code lines fix24_texts_v43 -> fix27_texts_v44, v43 -> v44; comments kept);
  b44_common.py, build_v44.py, check_macros_v44.py, check_titles_v44.py written by hand.

## Log
- Nov 24: dir, marker, log created.
- Nov 24: tools copied (copy_tools_v44.py; gen_records_v44 adds the two Windows RUN_NOTES (SUMMARY asserted), CHECK_v43, LATEX_v43, its pdf SHA256SUMS,
  fig_curves.pdf); b44_common.py (D9 = fix24, retag v43 -> v44, TITLE_REPS incl. date-line rule for "\emph{Not compiled.}").
- Nov 24: edits_v44.txt (33 edits in 11 files), erratum_part9_v44.tex (section T), records_errata_add_v44.md, build_v44.py (content assertions).
- Nov 24: check_macros_v44.py (wrapper names \fn,\key,\pth,... counted as used): FAILs exactly on \key in GS_lemma/erratum/flux2_lemma of v43, PASS v44.
  check_titles_v44.py (new): FAIL v43 (6 titles), PASS v44.
- Nov 24: first full build: balance (d) flagged inherited paths lemmaU.tex, make_fig.py in two short-title edits -> REPQ. Then all checks PASS:
  balance 23/0, macros PASS (v39 FAIL, v43 FAIL as expected), titles PASS, numbers PASS (0 STALE/UNMATCHED, 2 NEW-M38 = "1.2" section number,
  0 NEW-UNKNOWN hashes). No decimal number differs between v43 and v44 (per-file multiset diff). Two runs byte-identical (65 files). SHA256SUMS 25 files.
- Nov 24: nothing outside fix27 newer than .marker_start except the concurrent kit-v5 builder's paths (*v5*) and the /home/claude/s1spec dir mtime.
- Nov 24: LaTeX compile of v44 in new scratch C:\Users\micahw\kitruns\latex_v44 (no workaround): 15/15 compile, step3_main 56 pp, 0 errors,
  0 undefined refs/cites in the assembly. LATEX_v44.md, pdf/, latex_logs/. DONE.
