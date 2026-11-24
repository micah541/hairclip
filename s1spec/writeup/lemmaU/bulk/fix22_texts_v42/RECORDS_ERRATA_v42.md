# RECORDS_ERRATA_v42: RECORDS_ERRATA_v41.md verbatim, then the v42 section at the end

# RECORDS_ERRATA_v41: RECORDS_ERRATA_v40.md verbatim, then the v41 section at the end

# RECORDS_ERRATA_v40: RECORDS_ERRATA_v39.md verbatim, then the v40 section at the end

# RECORDS_ERRATA_v39: RECORDS_ERRATA_v38.md verbatim, then the v39 status section at the end

# RECORDS_ERRATA_v38: corrections to report and record statements (text and records only)

Nov 20, 2026. Directory `fix16_texts_v38/`. Source: `fix15_review/REVIEW_v4c_interface.md` (FREV15) items A-1r, A-2r, A-3r, B-v37-1, B-v37-3, B-v37-4, GS-1; `fix11_review/REVIEW_v4b_v36_proofs.md` (FREV11) items B1-B4, B6, C-7; and the v38 number check. **No report, log, checksum file or run directory was edited, and nothing was rerun.** Each item says what is wrong, what is right, and what to do.

Paths: `V` = `/home/claude/s1spec/step3nl_v4c`, `N` = `/home/claude/s1spec/writeup/lemmaU/bulk/fix13_v37chain`.

## R1. "Same sha256 of every cert input" is false for two run directories (FREV15 A-1r)
- **Wrong statement.**
  - `V/V4C_REPORT.md`, checks table, double-run row: "`runs/P_A` vs `runs/P` ... both pin the same sha256 of every cert input".
  - `N/V37_REPORT.md`, checks table, row "v4c double run (`C/runs/P_A` vs `C/runs/P`) and cert chain": it repeats the claim by reference ("PASS ... record note in V4C_REPORT").
- **Correct statement.** The numbers of `runs/P_A` equal those of `runs/P` (and `runs/reg_fin13v4b` equals `step3nl_v4b/runs/P`), ignoring timing and args [read, FREV15]. But `runs/P_A` and `runs/reg_fin13v4b` pin `v4bcert/out/e2v4b_constants.json` as `37386e13f753...`, which matches no file on disk. The file is now `7d6ada650e18...` in both `out/` and `out_A/`, and that is what `runs/P` pins (checked for this erratum with `sha256sum` and `grep` of the rate tables).
- **Cause** (FREV15 §1.4): `v4bcert/e2_consts.py:24` pins the text of `verifier/constants.py`. That file was edited at 09:24:52, after cert pass A (09:21-09:24:21) and before any verifier run. The change propagates through the pin of `e2_consts.json` into `e2v4b_constants.json`; the numeric content was the same.
- **Consequence.** Only `runs/P` (with `v4bcert/out`, pass B) is a verifiable production object. `V4C_LOG.md` step 4 ("renamed before any run") is true for verifier runs only.

## R2. The third cert pass has no driver record (FREV15 A-2r)
- `V/runs/certC.out` is empty. The cert stages write to a hard-coded `out/`, so the directories must have been swapped by hand to produce `v4bcert/out_A/` (the independent third pass that reproduces `out/` apart from `tc_weighted.json` runtime and its pin in `e0_identities.json`).
- Fix: one line in a new log, or make the pass a driver stage with a `V4BOUT` variable, at the next touch.

## R3. Stale checksum files (FREV15 A-3r, B-v37-1)
- `V/SHA256SUMS_v4c` lists itself, and `V4C_LOG.md` changed after it was written.
- `N/SHA256SUMS_v37` and `N/flux/out/inputs.json` pin `V37_LOG.md`, which was edited after the flux run and after the checksum file. Only the run name is read from that log.
- Fix: regenerate both checksum files (excluding themselves) before the re-pin; at the next flux touch pin a non-log file instead of `V37_LOG.md`.

## R4. Values printed below the certified bound (display roundings; FREV15 B-v37-3 and the v38 number check)
| where | printed | certified (key) | outward value |
|---|---|---|---|
| `N/V37_REPORT.md` (A) table | C0_2 114.763 (N2_ref) | 114.76328... (`A_interface.json` res.tight.r2_export.C0_ref) | 114.764 |
| same | C1_2 0.0698137 | 0.06981372... (C1_ref) | 0.0698138 |
| same | C0_2 116.359 (4 N2_ref) | 116.359231... (C0_ref4) | 116.3593 (7 digits; 116.360 at 6) |
| `V/V4C_REPORT.md` | A-1 margin x1.65695 | 1.6569477... (`e0_A2_budget_v4c.json` cases.2.5e-4.margin), a lower bound | x1.65694 |
| `V/V4C_REPORT.md`, `N/V37_REPORT.md` | whole-tail alpha ratio x1554.46 | 1554.4561... (lhs/rhs_whole_tail; 1554.4560... with the factor 1 + 2C'/L_F, flux C5) | x1554.45 |
| `FREV15` §3.2 | 5.3902e-4, 0.040935 | 5.390215...e-4, 0.0409351908... (`GS_boot_v4c.json`) | 5.3903e-4, 0.040936 |
| `V/e0_lemma/E0_LEMMA_v4c.tex` | d_ab = 0.0131136 | 0.013113620982 | 0.0131137 |
- None is consumed by a script. The v38 texts print the outward values. FREV11 B6 (V36_REPORT 114.782, 0.0698142) is superseded by the v37 values.

## R5. Header comments and labels (FREV15 B-v37-4; FREV11 B7)
- `N/run_chain_v37.sh` ll.2-8: the header comments describe fix6_v36 / fin13v4b. Cosmetic.
- `zone3.py:256`, `u1_ghat.py:96` (in `fix6_v36` and its v37 copies): the strings print c_dU <= 2e-8 as "[hyp]". It is a typed bound that the chain checks against the certified c_dU <= 1.46343e-8; the label should read "[typed, checked]". Stated in `lemmaU_v38.tex`; the scripts are not edited.

## R6. Chain-side record items still open (FREV11 B1-B4, C-7; unchanged in v37, FREV15 §2)
- **B1:** the `base_ibp` outputs (`lemmaS_F10p`, `cap_F10p`, `rigcap`, `weights_exact`, `assembly_arb`) carry no args/sha256 block, although they depend on the verifier. Pin them in the new manifest, or add args at the next chain touch.
- **B2:** `v36cfg.py:67-69` skips the env check silently if the rate table lacks `args`.
- **B3:** `t5_zone3/out/ramp_U.json` `u1_ghat_sha` matches neither the v35 nor the v36 `u1_ghat.json`; the content dependency was verified equal. Rerun `ramp_U` at the re-pin, or record an exception.
- **B4:** "ramp sups <= box" and "no margin < 1" are read off afterwards, not asserted in the chain (minimum margin 1.18446 [read]).
- **C-7:** `weights_t2.py:232` passes the zone-cell value through a nearest-rounded float with no arb assert per cell (deficit <= about 1 ulp against r_2 = 1e-3 [hand]; the same pattern at k = 0, 1).
- These need a chain touch, i.e. a rerun. They are listed in the v38 texts (`lemmaS_A_v38.tex` G-A6, `lemmaU_v38.tex` G-items) and in `TEXTS_REPORT_v38.md`.

## R7. The frame-shift citation (FREV15 GS-1)
- `fix14_interface/GS_lemma.tex` l.57 and `GS_tau1.py:48` cite "true-curve shift <= 3.3e-11 (step4_proof)". No such statement exists in `step4_proof.tex`.
- The computed bound is <= 1.65e-11 at Ybar_box (`fix15_review/c/c_frame_shift.json`, from 5-digit ball strings of the stored matrix), covered x60 by the 1e-9 widening. Corrected in `GS_lemma_v38.tex`. The script comment is not edited.

## R8. Coverage of the v37 number check (`fix12_texts/TEXTS_REPORT.md` §2)
- **What it claimed:** "every printed number with >= 3 significant digits ... is matched". Its tokenizer used the look-behind `(?<![\w.])`, so every number written directly after a control word (`\le116.3782`, `\times4.804`, `\ge0.97623`, ...) was never tokenised. Its record filter also exempted any line that mentions v35, v36, v37, fin13v4b or fix8. The "1077 match" therefore does not cover most bounds of the form `\le x`.
- **The v38 checker** (`tools/check_numbers_v38.py`) removes control words before tokenising, narrows the record filter, records the matched JSON key, flags superseded literals whose successor key differs (SUSPECT), and verifies 30 pinned key numbers with their rounding direction.
- What it found in the v37 texts: three values printed below the certified bound (erratum K10.1-3) and a set of valid but superseded v36 bounds. All are corrected in v38 (`tools/numsub_v38.json`, with keys).

## Recommendation: `runs/P_A` and `runs/reg_fin13v4b`
- **Pin only `V/runs/P`** (with `V/v4bcert/out`, including `strip_cubic_v4c.json` and `e0_A2_budget_v4c.json`). This is what the chain v37 and the flux reassembly consumed (FREV15 part B: all 259 pins in 131 args blocks point to `runs/P`, pass B).
- **Recommended: rerun both at the re-pin** into *new* directories (e.g. `runs/P_A2` with `V4C=1`, and `runs/reg_fin13v4b_2` with `V4C=0 TAG=fin13v4b`) on the current `v4bcert/out`.
  - Cost: about 12 min each, one process, nice 19, PYTHONHASHSEED=0, peak RSS about 80 MB (V4CR).
  - Expected result: equal to `runs/P` (resp. `step3nl_v4b/runs/P`) apart from timing, with every pin verifying.
  - Benefit: a verifiable double run and regression go into the manifest for the cost of two short runs.
  - This is a verifier run, so it was **not** done here and needs Micah's go together with the re-pin.
- **Also, in either case:** write an exception record for the old directories (e.g. `V/EXCEPTION_v4c_PA_reg_pins.json`). It should name the one unverifiable pin (`e2v4b_constants.json` 37386e13... vs 7d6ada65...), the cause (R1) and the numeric equality, and state that the old directories are not pinned.
- **Acceptable alternative** (no compute): the exception record alone, with `runs/P` the only pinned production run. FREV15 rated this adequate; the record then carries a double run whose input pins cannot be verified.


## Status after the step-4 re-pin v38 (text v39, Nov 20, 2026)
Source: `step4/p2/out_v38/REPIN_REPORT.md`, `step4/p2/out_v38/FINAL_manifest_v38.json` (`tools/records_v39.json` has the hashes).
- **R1 / Recommendation (runs/P_A, runs/reg_fin13v4b):** the "acceptable alternative" was taken. Only `runs/P` is pinned (stage 28). The two old directories are accepted only through `step4/p2/v38/EXCEPTION_v4c_PA_reg_e2v4b.json`, with 3 + 3 accepted pins. The recommended reruns into new directories were **not** done. Still open, but optional.
- **R2 (third cert pass has no driver record):** not addressed by the re-pin. Open.
- **R3 (stale checksum files, B-v37-1):** the flux outputs are pinned by stage 21, which reproduces them. The one non-args difference is the `V37_LOG.md` pin in `flux/out/inputs.json` (303647c2... pinned, c36c72c9... on disk). `SHA256SUMS_v37` was not regenerated. Open.
- **R4 (display roundings):** v39 adds O5 of the erratum (C1_ke printed `=1.7428e-3`, now `<= 1.7429e-3`).
- **R6 (chain-side items):**
  - B1: the `t2_A/base_ibp` outputs are now pinned by sha256 as stage-20 inputs of manifest_v38. They still carry no args block. Open.
  - B3: `ramp_U.json` is pinned by sha256 (55317ce5...). Its internal `u1_ghat_sha` mismatch is unchanged and has no exception record. Open.
  - B2, B4, C-7: unchanged. Open.
- **New, v39:**
  - (i) The pinned `GS_boot_v4c` script does not apply review items GS-2 (cell alignment) and GS-6 (the factor (1+2C'/L)^(1/2)). Both are rated negligible. Open.
  - (ii) The v33 manifest has one known F34 difference (`constants_needed_v34.json`), accepted by the v38 driver (`V33_KNOWN`). Earlier statements that manifest_v33 verifies with "673 pins" were false for that path.
  - (iii) The v38 builder's retag damage is listed in erratum O4. It includes a leaked NUL placeholder in a comment of `ei_lemma_v4c.tex`, which made `grep` treat the file as binary. The v39 balance check now rejects control characters.


## Status after the final review fix18 (text v40, Nov 21, 2026)
Source: `fix18_final_review/FINAL_REVIEW_v39.md` (R-1 to R-3). Nothing was rerun or re-pinned for text v40.
- **R-1 (pin counts):** manifest_v38 has 992 pin entries on 646 distinct files: 629 distinct stage paths and 55 scripts, 38 of which are also stage paths (counted from the manifest by `tools/gen_records_v40.py`, `tools/records_v40.json` key `pins`). The texts of v39 printed "629 paths and 55 scripts pinned"; v40 prints the full count. The manifest's own check reports `paths` 629 and `scripts` 55, which is consistent.
- **R-2 (transitive pins):** the verifier code (`rhs_iv.py`, `unc.py`, `bb.py`, `ival.py`) and `step3nl_v4c/cubic_E0fine_0.001_0.01_0.01.json` are pinned only through the args blocks of the pinned rate tables (12 inputs, all matching the disk according to the review); `tail_cert.py`, `mod_sups.py`, `integ.py` are not pinned at all. **Open:** pin them directly at the next manifest touch.
- **R-3 (items missing from the open list), now listed:**
  - fix15 A-4: float cell hulls in `strip_cubic` (immaterial per fix15). Open.
  - fix11 C-8: the count in `PROOFS_REPORT` (records only). Open.
  - fix11 D-6a: optional float check of the beta decomposition with the true functions. Open, optional.
  - fix9 §8: in `closure.py:130` the pad `1e-9*max(1,|U_k|)` is formed in double precision (immaterial per fix9). Open.
- Still open from earlier sections: R2, R3, R6 (B1-B4, C-7), GS-2, GS-6, M6.


## Status after the delta check fix20 (text v41, Nov 21, 2026)
Source: `fix20_delta_check/DELTA_CHECK_v40.md` (wording items W-1 to W-5 and one retag slip). Text v41 changes wording only. Nothing was rerun or re-pinned; no records item was added or closed. The open records items listed above are unchanged.


## Text v42 (Nov 23, 2026; Lean replaces DKS, review lean_bridge/REVIEW_W7.md)
- No step-3 or step-4 record changed; nothing was rerun or re-pinned. The records items listed above are unchanged.
- New citation record: the Lean files of `CWFlow/Curve/` (Micah Warren's repository `agents/lean/CWFlow`, commit pending) are identified by `writeup/lean_bridge/snap_final/SHA256SUMS` (14 files; its `All.lean` superseded) and `writeup/lean_bridge/snap_s3/SHA256SUMS` (`Periodic2Circle.lean`, `MaxShift.lean`, new `All.lean`; the S3 items of review W7). `tools/check_numbers_v42.py` re-verifies every hash against the snapshot and the repository files (read-only) and puts them in its hash registry.
- Records items added: (L-1) cite the commit of `CWFlow/Curve/` once it is made (replace "commit pending"); (L-2) record the exit status and log of the full `lake build CWFlow.Curve.All`; (L-3) `lean_bridge/IMPL_REPORT.md` is stale (review W7 S4); `IMPL_LOG.md` and `snap_final/SHA256SUMS` are the correct record.
