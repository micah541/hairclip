# Part B: (A) chain v37 and flux reassembly (fork of the fix15 reviewer)

Scope: `N` = `bulk/fix13_v37chain`. Read-only on the proof tree; everything written is in `fix15_review/b/`.
Scripts: `b_check.py` -> `b_check.json` (29 checks, 0 fail); `b_flux.py` -> `b_flux.json` (11 checks, 1 fail = finding B-v37-1); `cmp_reg_vs_v36prod2.txt` (rerun of N/cmp_runs.py, stdout only).
Labels: [cert] my arb recomputation on stored certified inputs; [read] read from the cited JSON/file; [hand] by reading code.

## Verdict: sound modulo items (records/text only; no numerical or wiring defect)

## What was checked
1. **Code diff [read].** All 204 scripts under N classified. 27 changed scripts; each `*_pre_v37` backup is byte-identical to the fix6_v36 (resp. fix8_flux2_impl/code) original. Every hunk is a path repoint (fix6_v36 -> fix13_v37chain, fix8 -> N/flux), the `fin13v4c` branch of `v36cfg.py:31-37`, `KJ_PAD` extended to fin13v4c (`v36cfg.py:49`), the `F2_VTAG` switch (`flux/code/f2common.py:11-15`), tag-parametrised file names, and args records (`V4C`, `F2_VTAG`). No mathematical line changed. 3 new drivers.
2. **Protected trees.** `audit_v4c/writeup_lemmaU_bulk_fix6_v36.sha` (1039 files) and `..._fix8_flux2_impl.sha` pass now; no file in either is newer than the start marker.
3. **Stale outputs.** `run_chain_v37.sh:21-30` deletes `t2_A/out/*`, t4/t5 outs (except ramp_U.json), the n1 chart and the five base_ibp outputs before each run. In-place outputs == `runs/prod2` byte for byte. Only two prod2 JSON equal v36's: `t5_zone3/ramp_U.json` (v35 byte copy, known B3) and `base_ibp/rigor/weights_exact.json` (regenerated; verifier-independent). No prod2 JSON mentions `step3nl_v4b`, `fix6_v36` or `fin13v4b`. Inherited caches consumed by the chain (`tables_exact.pkl`, `ramp_tables.json`, `prop1_forms*.pkl`, dated Nov 13) are byte-equal to fix6_v36's and verifier-independent. So the out_A-type problem did not occur here.
4. **Wiring [cert].** Consumed J0..J2, J4, slopes, peak are >= the raw `step3nl_v4c/runs/P` values (excess 2.6e-7 .. 2.2e-5 rel.); rho 2.679e-3 <= raw (8.4e-5); abar, beta equal mod_sups. The 1.5e-4 row dominates the 1.25e-4 row in all 12 consumed quantities. Rate table env: TAG fin13v4c, V4B=1, FLUX0 1.64e-20, FLUX2 8.0e-20. All 259 sha256 pins in 131 args blocks match the disk; pinned verifier files are `step3nl_v4c/runs/P/...` (pass B) and `step3nl_v4c/verifier/constants.py`.
5. **Numbers.** kappa_J certified 0.04242965327 reproduced, pad 1e-9 present [cert]. c_dU = 1.46342737776e-8, 2e-8/c_dU = x1.366654 [cert]; r0' = 3.75163377348e-9 [cert], equal to stored. Headline margins [read, prod2 JSON]: ramp M1 1.184464; Ghat22 1.916398; late M1 2.353209; H6 far K3 1.860867; zone III box b1 1.238888; C0_2/C1_2 114.76328 / 0.06981372 (N2ref), 116.35923 / 0.26059988 (4 N2ref). All 133 margin rows >= 1.18446; min ratio v37/v36 = 1.000000.
6. **Regressions.** cmp_runs rerun: reg_v36 vs fix6_v36/runs/prod2 = 25 identical + 12 path-only, PASS; prod1 vs prod2 all identical (my filecmp over every JSON too).
7. **Flux.** `flux/out` == `snap/v37_pass1` == `snap/v37_pass2` (9 files, bytes); regression snapshot == fix8 out apart from args/timing. No pin points to v4b/fix6/fix8; chain data are read from `runs/prod2`. Recomputed [cert]: FLUX0/K0bar x4.80971; combined x4.80453 (implied eps_I = 3e-5, mu = 0.08); C1 x9.85e6; C6 x141.32. [read]: C' 2.65189e-26, C3 pad x1547.6 / x2236.4, C5 x1554.46 / x1512.16; C0-C7, C9 pass, C8 = None. Versus fix8: K0bar, K2pt, C', S' all smaller (x0.99984..0.99998); C4 and C5 smaller by 4e-4 / 2e-4 relative (L_F fell), harmless.

## Findings
| # | file:line | severity | finding | fix |
|---|---|---|---|---|
| B-v37-1 | `N/flux/out/inputs.json` (pin of `N/V37_LOG.md`); `N/SHA256SUMS_v37` | minor (records) | V37_LOG.md was edited at 10:41:16, after the flux run (10:39) and after SHA256SUMS_v37 (10:41:07). The flux pin and the checksum line for V37_LOG.md no longer match the disk. The only content the flux stage needs from it is the name of the run (`runs/prod2`), which is unchanged. | Regenerate SHA256SUMS_v37; pin a file that is not a live log (e.g. `runs/prod2/env.txt`) at the next flux touch, or record an exception. |
| B-v37-2 | `N/flux/out/flux2_cert.json` conditions.C8; `flux/texts/*` | text | C8 note still says "OPEN"; texts carry fix-8 wording (`_v36` file names). V37_REPORT says so. | Next text round. |
| B-v37-3 | `V37_REPORT.md` (table "(A) chain") | text | C0_2 "114.763", "116.359", C1_2 "0.0698137" are rounded down (true 114.76328, 116.35923, 0.06981372). Same as B6. Key-table headers "v35"/"R3" mean v36 / reg_v36 (stated). | Round up. |
| B-v37-4 | `N/run_chain_v37.sh:2-8`, script comments | text | Header comments still describe fix6_v36 / fin13v4b. | Cosmetic. |

## Status of fix11 findings in v37
- B1 open: `cap_F10p.json`, `weights_exact.json`, `assembly_arb.json` have no args block; `lemmaS_F10p`, `rigcap` have args without sha pins.
- B2 open (`v36cfg.py:74-78` still `if a:` / `if k in env`); harmless here since the fin13v4c table has the keys and the values were checked.
- B3 open: `ramp_U.json` `u1_ghat_sha` 4b86cc91... vs v37 u1_ghat.json 6f9790f8...
- B4 open (no chain-side assert of "no margin < 1").
- B5 open (`claims_failing` non-empty in asm_*.json).
- B6 superseded by B-v37-3. B7 open (`zone3.py:274` "[hyp]").
- C-7 open (`t2_A/rigor/weights_t2.py:232, 256, 262` unchanged).

## Not checked
The inherited v35/v36 mathematics beyond the diff; n1_chart content; the chain and flux stages were not rerun; flux texts; whether fin13v4c itself is right (Part A).
