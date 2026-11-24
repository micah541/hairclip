# TEXTS_REPORT_v44: step-3 text set v44 (check fix25 of text v43; LaTeX defects of the v43 compile fix26)

Nov 24, 2026. `D` = `/home/claude/s1spec/writeup/lemmaU/bulk/fix27_texts_v44/`. Source: `fix24_texts_v43/` (read-only, unchanged).
**Step 3 is not closed.** No number computed or changed; no Lean file built or edited; nothing in the proof tree run.

## 1. What changed (`tools/edits_v44.txt`: 33 asserted edits in 11 files; erratum Part 9 = section T; records v44 section; patch-note v44 note)
- **F-1** (records): `RECORDS_ERRATA_v44.md` v44 section corrects TEXTS_REPORT_v43 §1 / TEXTS_LOG_v43 (they describe the first v43 build; the
  final build cites the raw logs; 2 NEW-M38, not 1). The closed v43 files are not edited.
- **F-2** erratum Part 8 and the RECORDS v43 section: "no Lean file was built" corrected in place, marked [text v44, check fix25 F-2]
  (per-file recompile of the three S3 files, sources unchanged, logs rc=0).
- **F-3** U-1 §1.2: axiom checks cover "the dependencies of these theorems", not the import closure. **F-4** C-1 gaps item "Lean" names
  fix23 Part A for (PC), (MS), (GE). **F-5** U-1 Lemma "the Lean velocity is CDF" stated for (a) C^inf X: R -> C with X' != 0 on R and
  (b) the chart lift of an immersion of Circle. **F-6** main_v3 caveat (b): compiled per file, full build and commit pending (status box), plus
  the [hand] steps of U-1 and C-1.
- **Note**: Windows runs of both v4 kits (24 Nov) recorded in step3_main (open list "Re-pin and kits"), main_v3 (§ Independent reproduction,
  status box, open items (a'), Reproduction), step4_handoff §Kits, erratum T, RECORDS (K-W): kit3 reproduces (72 PASS; 6 FAIL CRLF-only;
  3 WARN float last-bit; 2 run-copy fixes), kit4b TOTAL PASS [159/0/1] after 4 run-copy fixes. "kit v5 pending", no hash.
- **LaTeX defects (fix26)**: `\providecommand{\key}` (lemmaU's body) in GS_lemma, erratum, flux2_lemma; short titles for the step3_main
  Lemma U section, the two ei_lemma titles, the G1 caption and two main_v3 captions (all non-starred titles with \fn); `fig_curves.pdf`
  (byte copy of writeup/main/fig_curves.pdf, sha in records_v44.json) beside main_v3_v44.tex.
- **Checkers**: `check_macros_v44.py` counts wrapper macros (\fn, \key, \pth, \texttt, ...) as used, per standalone document — it FAILs on
  v43 exactly on \key in GS_lemma/erratum/flux2_lemma (asserted in build_all), PASS on v44, FAIL on v39. New `check_titles_v44.py`
  (moving arguments): FAIL on v43 (6 titles), PASS on v44.
- **Records/open lists (C)**: v43 compiled (fix26) cited as the v43 compile; defects fixed in v44; "v44 compile pending" in the texts
  (date lines: "Text v43 compiled (MiKTeX, fix26...); text v44 not yet compiled"). Also: the lemmaS_A provenance of numbers_A named a
  nonexistent fix17 generator path (retag artefact) -> fix27 generator; U-1 names check_numbers_v44.py.

## 2. Checks (`bash tools/build_all_v44.sh`)
balance 23 PASS / 0 FAIL; macros PASS (v43 FAIL as expected, v39 FAIL as expected); titles PASS (254 seen; v43 FAIL as expected);
numbers PASS (0 STALE/UNMATCHED/SUSPECT, 56 pinned OK, 126 macros OK; new text: 2 NEW-M38 (the section number 1.2, value coincidence),
0 NEW-MISS, 0 NEW-UNKNOWN hashes; Lean hashes verify). No decimal number differs between v43 and v44 (per-file multiset diff).
Two full runs byte-identical (65 files); SHA256SUMS (25 files) verifies. Nothing outside D newer than `.marker_start`, apart from the
concurrent kit-v5 builder's own paths (`kit_v5_build/`, `kit_step4b_v5_build/`, `kit3v5_run_wsl1/`, `kit_step3_v5.zip`, ...) and the
`/home/claude/s1spec` directory mtime they cause; no Lean repo file changed.

## 3. LaTeX compile of v44 (`LATEX_v44.md`, `pdf/`, `latex_logs/`)
MiKTeX 26.2, new scratch `C:\Users\micahw\kitruns\latex_v44`, no command-line workaround: all 15 standalone documents compile; assembled
step3_v44_main 56 pp, 0 errors, 0 undefined references, 0 undefined citations. The texts still say "text v44 not yet compiled" (built before).

## 4. Open (unchanged otherwise)
Kit v5 (being built separately) and its runs; full `lake build CWFlow.Curve.All` (+ log); commit of `CWFlow/Curve/`; (P3'), Lemma O arc by
arc, single-point statement; step-4 patch (authors' go); final review. The v44 changes are not independently checked; PDFs not proofread.
