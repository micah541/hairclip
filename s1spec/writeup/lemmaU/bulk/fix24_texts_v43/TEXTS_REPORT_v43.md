# TEXTS_REPORT_v43: step-3 text set v43 (delta check fix23 of text v42)

Nov 23, 2026. `D` = `/home/claude/s1spec/writeup/lemmaU/bulk/fix24_texts_v43/`. Source: `fix22_texts_v42/` (read-only, unchanged). **Step 3 is not closed.** No number computed; no certified value changed; no Lean file built or edited; nothing in the proof tree run.

## 1. What changed (`tools/edits_v43.txt`: 41 asserted edits in 9 files; erratum Part 8 = section S; records v43 section; patch-note v43 note)
- **S-1** U-1 §1.2: compile records of the S3 files. The logs `writeup/lean_bridge/raw_{Periodic2Circle,MaxShift,All}.log` did **not** exist at build time, so the text says "no raw compile log was kept; compile logs: see `lean_bridge/IMPL_LOG.md` 'S3 links'; per-file logs pending". `build_v43.py` switches to citing the three logs automatically if all three exist: rerun `bash tools/build_all_v43.sh`. Also: the full `lake build` log is to be kept (records L-2); new records item L-4.
- **S-2** U-1 gaps "Lean": the "All.lean was never compiled" clause removed; (PC), (MS), (GE) added to the compiled list.
- **S-3** erratum R3(b) rewritten (endpoint form proved by `isFlow_glue_endpoint`; `rst_glue` only as (Gl)), marked "[Rewritten in text v43]".
- **S-4** lemmaU standing facts: (P1) in sup form, `sup_{t<T_max} ∫κ²ds = ∞`.
- **W-1** step3_main and main_v3 status boxes, continuation (P1'), C-1 (I4) and label line, U-1 source line and gaps, erratum R4/R6: W7 covers snap_final only; the S3 theorems were reviewed in `fix23_delta_check/DELTA_CHECK_v42.md` Part A (faithful). R6 no longer lists "independent check of the S3 theorems"; text-v42 changes now marked checked (fix23); v41 and v43 wording changes marked not independently checked.
- **W-2** axiom checks for (T4), (T5) in `Criterion.lean`. **W-3** (MS) range `t∈[0,T_max−t₁)`. **W-4** U-1 gluing lemma stated with the existential Ψ of (GE), piecewise step [hand]; C-1 Thm (v) identification T_max = T_*−t₁ written out and [hand]. **W-5** main_v3 caveat: (P1) no longer "literature". **W-6** main_v3 13-Nov box marked superseded. **W-7** CWFlow bib entries (continuation, step3_main), C-1 label, b12_preamble, main_v3 box, continuation l.74, erratum R2 name both `snap_final/SHA256SUMS` and `snap_s3/SHA256SUMS`. **W-8** bootstrap C-1 row: "compiled per file, full build and commit pending". **W-9** (PC) names c(θ)=e^{iθ}; `\tD` removed from U-1.

## 2. Checks (`bash tools/build_all_v43.sh`)
- builder assertions: every edit's old string occurs exactly once; no v42 set-file name left outside history files; no "All.lean never compiled" except as a quote; no limit form of the L² criterion; `\tD` gone; (MS) range present; every non-history file naming snap_final/SHA256SUMS also names snap_s3/SHA256SUMS; fix23 cited in the five status-bearing files; erratum R3(b)/R4/R6 old sentences gone.
- balance 23 PASS / 0 FAIL (path roots + fix22, fix23, snap_s3). macros: PASS on v43, PASS on v42, FAIL on v39 as expected (12).
- numbers: PASS (0 STALE/UNMATCHED; 126 macros; 56 pinned; hashes 0 NEW-UNKNOWN, Lean snapshots re-verified vs repo). 1 "NEW-M38" = the section number 1.2 in erratum Part 8, matched by value coincidence; no new numeric value written.
- two full runs byte-identical (58 files); SHA256SUMS (24 files) verifies.
- nothing outside `D` newer than `.marker_start` (find over /home/claude and the CWFlow/Curve repo dir).

## 3. Open (unchanged)
Full `lake build CWFlow.Curve.All` (+ log); commit of `CWFlow/Curve/`; per-file S3 compile logs; (P3'), Lemma O arc by arc, single-point statement; Windows kit runs; LaTeX compile; step-4 patch (authors' go); final review. The v43 changes are not independently checked.
