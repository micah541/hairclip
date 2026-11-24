# TEXTS_LOG_v43 (resumable log of text set v43)

Task: text set v43 from text v42 (`fix22_texts_v42/`, read-only): apply the delta check `fix23_delta_check/DELTA_CHECK_v42.md`
(should-fix S-1..S-4, wording W-1..W-9). Rules: new files only in fix24_texts_v43/; nothing outside changes; step 3 NOT closed.

## Resume
- Tools in `tools/`: `b43_common.py`, `build_v43.py` (reads `tools/edits_v43.txt`, `tools/erratum_part8_v43.tex`, `tools/records_errata_add_v43.md`),
  `build_all_v43.sh` (whole pipeline + checks + SHA256SUMS). Rerun: `bash tools/build_all_v43.sh`.
- Marker for the "nothing outside changed" check: `.marker_start` (created at start).
- S-1: build_v43.py cites `writeup/lean_bridge/raw_{Periodic2Circle,MaxShift,All}.log` if all three exist at build time, else the IMPL_LOG fallback.
  If the logs appear later, rerun `bash tools/build_all_v43.sh` (texts change; SHA256SUMS regenerated).

## Log
- Nov 23: dir + marker; checker tools copied from fix22 (code lines v42 -> v43; check_balance NAMEDIRS + fix22, fix23, snap_s3); keymap/whitelist file keys renamed.
- Nov 23: b43_common.py (D9 = fix22, retag v42 -> v43, TITLE_REPS v43); edits_v43.txt (S-1..S-4, W-1..W-9); erratum Part 8; records v43 section; build_v43.py with content assertions.
- Nov 23: first build: R6 assertion too broad (Part 8 quotes the phrase) -> narrowed; fallback records line named a brace path -> reworded. 
- Nov 23: all checks PASS (balance 23/0; macros PASS on v43, PASS on v42, FAIL on v39 as expected; numbers PASS). Two runs byte-identical (58 files); SHA256SUMS 24 files OK.
- Nov 23: S-1 raw logs raw_{Periodic2Circle,MaxShift,All}.log did NOT exist at build time -> fallback sentence ("see IMPL_LOG.md 'S3 links'; per-file logs pending"). Rerun build_all_v43.sh once they exist.
- Nov 23: nothing outside fix24_texts_v43/ newer than .marker_start (find over /home/claude and the CWFlow/Curve repo dir). DONE.
