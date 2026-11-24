# TEXTS_LOG_v42 (resumable log of text set v42)

Task: text set v42 from text v41 (`fix21_texts_v41/`, read-only): (1) Lean replaces DKS per `writeup/lean_bridge/REVIEW_W7.md` section 3 (12 items);
(2) S3a/S3b marked [hand] with one switchable macro `\LeanSthree` per file ("Lean formalization in progress (lean_bridge S3)");
(3) v41 leftovers (step3_main title "text v39"; open list "v40 check pending"); (4) open-item lists synced.
Rules: new files only in fix22_texts_v42/; nothing outside changes; step 3 NOT closed.

## Resume
- Tools are in `tools/`: `b42_common.py`, `build_v42.py` (reads `tools/edits_v42.txt`, `tools/erratum_part7_v42.tex`, `tools/records_errata_add_v42.md`),
  `build_all_v42.sh` (whole pipeline + checks + SHA256SUMS).  Rerun: `bash tools/build_all_v42.sh`.
- Marker for the "nothing outside changed" check: `.marker_start` (created at start).

## Log
- Nov 23: dir + tools copied from fix21 (sed v41->v42); check_balance NAMEDIRS + Lean roots; check_numbers hash registry + Lean snapshot hashes (re-verified vs snapshot and repo).
- Nov 23: verified snap_final/SHA256SUMS == repo CWFlow/CWFlow/Curve/*.lean (SNAP_MATCH); CWFlow HEAD 8718f95 (Curve/ untracked), MM 7b15b01, toolchain v4.33.1.
- Nov 23: TITLE_REPS v42 (incl. step3_main title fix).
- Nov 23: edits authored with tools/authoring/make_edits.py (stage = retagged v41 copies in scratchpad; every old string asserted unique) -> tools/edits_v42.txt (49 edits, 11 files).
- Nov 23: first build: macro def-before-use fix in U1 (macro moved to preamble); toolchain name added to ALLOWED_ABSENT; repository All.lean changed under us (S3 agent).
- Nov 23: coordinator: S3a/S3b PROVED in Lean (snap_s3: Periodic2Circle f7d02857, MaxShift 80b86094, All 91d5788e). Switched \LeanSthree to print [Lean]; new items (PC), (MS), (GE) in U1 1.2; C-1 lem:max and Thm (v) now cite them; U-1 gluing lemma back to endpoint form via isFlow_glue_endpoint; S3a/S3b removed from open lists; checker verifies snap_final + snap_s3 against repo.
- Nov 23: all checks PASS; two runs byte-identical (61 files); SHA256SUMS 24 files OK. DONE.
