# TEXTS_LOG_v39 (resumable log; text set v39, fix17_texts_v39)
Task: v38 texts -> v39: cite manifest_v38 / driver / REPIN_REPORT / exception records / kits v4 (WSL runs only; NO Windows run of v4 kits: Micah has not run them, list as open); step-4 hand-off wording file + patch note; 0.1475 check; build_all_v39.sh; checks.
Rules: new files only in this dir; never edit fix16 or anything else; nice 19, one process.
Reference timestamp for the outside-change check: scratchpad/ref_start (touched at start).

## Steps
- [x] 1. dir created, start_time.txt
- [x] 2. tools/b39_common.py (retag of v38 file names only; no generic 'v38'->'v39' word retag, because v38's GEN rule 'text v37'->'text v38' damaged historical phrases)
- [x] 3. tools/gen_records_v39.py -> tools/records_v39.json (all hashes/counts/results from disk); tools/v39_values.py (outward rounding)
- [x] 4. tools/make_step4_v39.py -> step4_handoff_v39.tex, STEP4_PATCH_NOTE_v39.md (quotes asserted in step4_proof.tex)
- Findings so far: v38 GEN damage: 'text v38 printed/had' (should be v37) in continuation (3x), U1, ei_lemma (2x), main_v3, step3 main; GS_lemma l.107 points to fix16 continuation_v38 l.41 instead of fix12 continuation_v37 l.41; step3 main 'runs/reg_fin13v4c' (dir is reg_fin13v4b); lemmaU date line self-reference.  GS_boot certified upper bounds round UP to 5.3903e-4 / 0.040936 (caller's 5.3902e-4/0.040935 are truncations).
- [ ] 5. tools/build_v39.py (all text edits)
- [x] 5. tools/build_v39.py (imports make_step4_v39); erratum Part 4 source tools/erratum_part4_v39.md; records addition tools/records_errata_add_v39.md
- [x] 6. tools/gen_numbers_A_v39.py, gen_ei_constants_v39.py (copies): 122 macros and EI table identical to v38
- [x] 7. checks: tools/check_balance_v39.py (+cross-doc refs, inputs, new-text paths, control chars), tools/check_numbers_v39.py (+(8) new numbers vs manifest_v38 data, (9) hash registry, (10) 0.1475; keymap +26 rows keyed to FINAL_manifest_v38.json)
- [x] 8. tools/build_all_v39.sh: balance 23 PASS 0 FAIL; numbers PASS (pinned 56/56, NEW-MISS 0, NEW-UNKNOWN hashes 0); SHA256SUMS
- More findings: v38 lemmaU paired fix10 sha256s with v38 file paths; v38 ei_lemma_v4c.tex has a NUL placeholder in a comment (grep sees it as binary); C1_ke printed '=1.7428e-3' (truncation).
- [ ] 9. outside-change check, TEXTS_REPORT_v39.md
- [x] 9. outside-change check (find -newer start marker on s1spec, s1_*): 0 files outside fix17 (only parent dir mtime); fix16 SHA256SUMS ok; rebuild byte-identical; TEXTS_REPORT_v39.md written.  DONE.
