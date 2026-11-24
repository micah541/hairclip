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
