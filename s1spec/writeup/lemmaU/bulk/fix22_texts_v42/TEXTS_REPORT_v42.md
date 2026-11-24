# TEXTS_REPORT_v42: step-3 text set v42 (Lean replaces DKS)

Nov 23, 2026. `D` = `/home/claude/s1spec/writeup/lemmaU/bulk/fix22_texts_v42/`. Source: `fix21_texts_v41/` (read-only, unchanged). **Step 3 is not closed.** No number computed; no certified value changed; nothing in the proof tree run.

## 1. What changed (`tools/edits_v42.txt`, 49 asserted edits in 11 files; erratum Part 7 = section R; records section; patch-note sentences)
- **Lean replaces DKS** (REVIEW_W7 §3, items 1-12): U1 Labels ([Lean] = compiled, commit pending, snapshot SHA256s; [DKS] dropped), §1.2 rewritten (source/compile status; setting with Mathlib `Circle` modelled on `EuclideanSpace R (Fin 1)` -- v41's "unit sphere of EuclideanSpace R (Fin 2)" corrected; `bendingEnergy` = ∫κ²|c'|dθ of θ ↦ F(e^{iθ}); theorems (L1) `cdf_exists`, (L2) `cdf_unique`, (T1) `cdf_exists_maximal`, (T2) `cdf_extends_of_bounded_energy`, (T3) `aNormSq_eq_curv_sq`, (T4) `IsFlow.lift_isCDF`, (T5) `IsCDF.curv_sq_le`, (V) `flowVel_eq_neg_kss`, (Sy) `cdf_symmetric`, (Sh) `IsFlow.shift_gen`, (Gl) `rst_glue`, (PC) `exists_isLagImmersion_of_periodic`/`lift_injective`, (MS) `cdf_maximal_shift`, (GE) `isFlow_glue_endpoint`); §1.3 DKS = literature only; §2 Lemmas vel/glue/max/sym/bu now [Lean]; Remark D1 → T2; gaps. C-1: labels, "what is closed", (I4), (I5), lem:max proof, Thm (v) proof, gaps. Status/open-item sentences in continuation, bootstrap, step3_main, step4_handoff, ei_lemma, GS, lemmaU, main_v3, b12_preamble; bibliography `CWFlow` entries ("commit pending", not d72b744; d72b744 kept only as "text v41 cited ..."). DKS/Chou kept as attribution. Author item A-2 removed.
- **S3a/S3b**: first marked [hand]; after the coordinator's update they cite the new Lean theorems (snapshot `writeup/lean_bridge/snap_s3/SHA256SUMS`). One macro per file, `\LeanSthree`, now prints `[Lean]` (provenance kept; switchable). These S3 theorems have **not been independently reviewed** (W7 predates them); I read their statements only.
- **v41 leftovers**: step3_main title now "text v42" (was "text v39"); open list says the v40 check is done (fix20 OK).
- **Open lists** now: full `lake build CWFlow.Curve.All`; commit + citation of `CWFlow/Curve/`; (P3'), Lemma O arc by arc, single-point statement; Windows kit runs; LaTeX compile; step-4 patch (authors' go); final review (text v41/v42 changes not independently checked).
- `main.bib` has no Lean entry; main_v3 says so (cited via \fn).

## 2. Checks (`bash tools/build_all_v42.sh`)
- balance 23 PASS / 0 FAIL (new-text path check now also resolves lean_bridge and the CWFlow repo; toolchain name allowed).
- macros: PASS on v42; FAIL on v39 as expected (12); PASS on v41.
- numbers: PASS (0 STALE/SUSPECT/UNMATCHED; 126 macros; 56 pinned). 16 "NEW-M38" new numbers are all section/theorem numbers (1.2, 3.1, 4.1, 7.1) matched by coincidence of value; no new numeric value was written. Hashes: 0 NEW-UNKNOWN; the registry now includes snap_final + snap_s3, each re-verified against the snapshot and the repository files (0 mismatch; snap_final's All.lean superseded by snap_s3's).
- builder assertions: every file using `\LeanSthree` defines it once before use (step3_main via b12_preamble); "d72b744/read from source/not elaborated" survive only in sentences about text v41 or older.
- Two full runs byte-identical (61 files); SHA256SUMS 24 files verifies.
- Nothing outside D was written by this work. `find -newer` shows outside changes only from other agents during the session (Lean repo `CWFlow/Curve/{MaxShift,Periodic2Circle,All}.lean` + build outputs, `lean_bridge/snap_s3/`, `lean_bridge/IMPL_LOG.md`, `/home/claude/CLAUDE.md`); my tools read those only.

## 3. Tools
Copies of the v41 tools with v42 names (`b42_common.py` new TITLE_REPS; `build_v42.py` new; `check_balance_v42.py` roots; `check_numbers_v42.py` Lean hash registry). Authoring aid: `tools/authoring/` (generator of edits_v42.txt from staged v41 text; not part of the build).
