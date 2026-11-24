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
