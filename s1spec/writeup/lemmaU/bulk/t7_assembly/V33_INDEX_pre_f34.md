# v3.3 step-3 file map (Lemma U and the repaired Lemma S(A))

Task 7b, Nov 13 2026 (Sonnet). Bookkeeping only: no number was changed. Hashes were computed with `sha256sum` when this file was written (Nov 13 2026). Labels: **[cert]** = python-flint arb with outward rounding; **[paper]** = hand proof; **[hyp]** = hypothesis. Status of the whole set: **not yet independently reviewed (task 8), LaTeX not compiled**.

> **SYNC NOTE (Nov 14, 2026, documentation sync).** Sections 1.2 to 1.5, 1.7 and 2 below describe the chain of Nov 13 (t2_A_v33, U_cert_v33, the 10-stage step-4 manifest). They are **superseded by section 4 at the end of this file** (t2_A_v34 with B0 and M14, U_cert_v34, zone3 a381b1f2..., the 19-stage manifest 8f6f7f95...). The old text and hashes are kept for the record; the margins quoted there are the Nov 13 values (some rounded up, e.g. x3.19, x124.8, x1.24; section 4 gives the rounded-down values).

> **SYNC NOTE 2 (Nov 14, 2026, N35).** Section 4 (v34 chain) is in turn superseded in its paths, hashes and counts by **section 4.6 at the end of this file** (chain v35 after the re-review items N1/N2: `fix2/n12/t2_A_v35`, `U_cert_v35`, zone3 `45e8d4dc...`, the 20-stage manifest `a0d348ec...`). The v34 text and hashes of section 4 are kept for the record.

Abbreviations: `B` = `/home/claude/s1spec/writeup/lemmaU/bulk`, `S4` = `/home/claude/s1spec/step4`, `V` = `/home/claude/s1spec/step3nl_v33`.
Closed files (v3.1, v3.2, `step3nl/`, `S4/p2/out/FINAL_manifest.json`) are not edited by any v3.3 component. Environment of every driver: `PYTHONHASHSEED=0`, `nice -n 19 python3 -B`, one process at a time.

## 0. Dependency order (one line)
verifier (V) -> repaired S(A) (t2_A_v33) -> Lemma U certificates (t4_cert: u1..u5, assemble) <-> zone III (t5_zone3) -> proof text (t3_proof) ; step-4 v33 stages (S4/p2/v33) are independent of t2/t4/t5 (they need only the closed step-4 outputs and V's `mod_sups`).

## 1. Components

### 1.1 Verifier, tag fin12f4 (FLUX0 = 1.64e-20, FLUX2 = 8.0e-20, x4 against fin11)
- Path: `V/verifier/` (copy of `step3nl/`; original untouched). Note: `V/V33_VERIFIER.md`.
- Driver: `cd V/verifier && sh run_v33.sh` (tail_cert; closure 1.5e-4; closure 1.25e-4). Then `cd V/verifier && python3 mod_sups.py fin12f4 modsups_fin12f4` (needs both alpha tables; run it in the same environment as run_v33.sh, see V33_VERIFIER.md; the exact command line used was not recorded beyond that).
- Time: closure 1.5e-4 214.8 s, closure 1.25e-4 194.5 s (both recorded in the `.out` files). tail_cert and mod_sups times: not recorded.
- Outputs (sha256 now):

| file | sha256 |
|---|---|
| `V/verifier/rate_table_fin12f4_a0.00015.json` | b20e87153440144b0537afa4cd4b9dff2d788ffd2ce2163626faf973b13d2ad1 |
| `V/verifier/rate_table_fin12f4_a0.000125.json` | 79b60067f0f935eec8807a7641ba0af0b514755f7ac0ae553488f1a0e36fef65 |
| `V/verifier/modsups_fin12f4/final_alpha_points_5e-8.json` | d5b14797e932bb01c3ce2c4cd89c56fd994a42647ac6e4c0133fbb808d3cfde6 |
| `V/verifier/tail_cert_log.txt` | d2237124555e0f64802797ce28a8ad5457427839eb933513d1e7a3aaaaaa945f |

- Result [cert]: tail_cert PASS; closure succeeds at both alpha1; sup|(a,b)| 2.19103e-07, sup|beta| 3.33811e-09, peak sqrt(Lambda) 5.3075e-08; only slopes, sJ4 and F0 move (x4).

### 1.2 Repaired Lemma S(A), v3.3 (t2_A_v33)
- Path: `B/t2_A_v33/`. Reports: `V33_A_REPORT.md`, `A_REPAIR.md`, `A2_REPORT.md`, `T2_REPORT.md`. G1/G3 fix notes: `B/t4_cert/G13_FIXES.md` (A_interface export `BETA_A = '2'`).
- Driver, in order:
  1. `bash B/t2_A_v33/base_ibp/run_chain_v33.sh` (4 steps; reference chain with the fin12f4 slopes; about 2.5 min).
  2. `bash B/t2_A_v33/run_t2.sh` (hjw_rep, ic_rep x4, a3_signs, weights, ramp_rep, outer_L2, 7 assemblies, export_t2; about 40 s per `V33_A_REPORT.md`, log timestamps 20:26:19 to 20:26:58).
- Outputs: `B/t2_A_v33/out/*.json` (19 files) and `out/SHA256SUMS` (verifies: `cd B/t2_A_v33/out && sha256sum -c SHA256SUMS`, all OK now).

| file | sha256 |
|---|---|
| `B/t2_A_v33/out/A_interface.json` (the interface consumed downstream) | 5237a96a1c6ca873e0f6f214f413c1a39996002155b8f128920d9cdacb087aa4 |
| `B/t2_A_v33/out/SHA256SUMS` | 348076020cb7a682aa6106804305b326c3d70ca9810c10db3b2ec1e8df36c371 |

- Result [cert] (tightened (H2) + (H2f) datum, beta_A = 2, (H2f)_0 = 1e-6): E2^A = 0.0566; S0v x2.14; core M1 x4.44; G22 x1.92; FLUX0/K0 x3.51 and FLUX2/K2 x6.20 (thinnest); R5' adds 2.1e-3 at 4 N2_ref. Two consecutive `run_t2.sh` runs were byte-identical (per the report).
- Superseded, kept: `B/t2_A/` (v3.3-pre), `t2_A_v33/out_pre_v33/`, `t2_A_v33/out/A_interface_pre_g13.json`.

### 1.3 Lemma U certificates (t4_cert)
- Path: `B/t4_cert/`. Reports: `T4_REPORT.md`, `U4_REPORT.md`, `G13_FIXES.md`.
- Driver (the v3.3 chain; `run_all.sh` is the pre-v33 default chain and writes `U_cert.json`, which is superseded):
  1. `cd B/t4_cert && R="env PYTHONHASHSEED=0 nice -n 19 python3 -B"`
  2. `$R u1_ghat.py` (needs `t5_zone3/out/ramp_U.json`), `$R u2_K.py --fresh`, `$R u34_sweeps.py --case primary`, `$R u5_enum.py` (the nophi/fallback variants are not used by v3.3).
  3. `$R assemble_U.py --A B/t2_A_v33/out/A_interface.json --out U_cert_v33.json` (under 1 s; writes `out/U_cert_v33.json` and `out/U_cert_for_zone3.json`, see the caveat in section 3).
- Times (recorded in logs): u1 1.2 s, u34 51.9 s, u5 35.1 s, assemble under 1 s; u2_K time not recorded. Whole chain about 5 min per T4_REPORT.

| file | sha256 |
|---|---|
| `B/t4_cert/out/u1_ghat.json` | 4b86cc91ffe63e0fafd2b23133d8f5382047a1f376ca31294042fe1a141b1616 |
| `B/t4_cert/out/u2_K.json` | eddf340e061f3cf59a735a59b2baba58157ef230650f36dfa85e7ca98eb8483a |
| `B/t4_cert/out/u34_primary.json` | cdaf0aba0f4271dccca834bf1d1aa53ef33033313af8cd3f549b7e649df88763 |
| `B/t4_cert/out/u5_enum.json` | 10491e1f4dda0d06d0dc3c0be0867a5250ed2b07e23cfc706be67523235f63f2 |
| `B/t4_cert/out/U_cert_v33.json` (**the Lemma U result**) | 42cf9617941456e36230720674f571e448432d8e0a06a9bab492e75edc6caea5 |
| `B/t4_cert/out/U_cert_for_zone3_v33.json` | 8e58d9d1532e490c822d05a173dd2a764f8a645c1cd3d6afa643aee7646b3446 |
| `B/t4_cert/out/SHA256SUMS` | 23ba5857dbc7a720ae1cd4466165b8edcfcbc9d78c1884446c4433f774531974 |

- Result [cert]: late M1 x2.358, M0 x4.90; present M1 x2.775, M0 x2.225; loop gain 0.0390; rate (R2) x28.2; Ghat_22 = 348341 (with phi); U-fine Bhat = (4.8717e-6, 2.6258e-4, 0.02578); c_dU = 1.7186e-8; M1z = 0.534 (new hypothesis of the bootstrap). Enumerator cross-check OK. `SHA256SUMS` verifies now.
- Superseded, kept: `U_cert.json`, `U_cert_fallback.json`, `U_cert_v33_pre_g13.json`, `U_cert_v33_pre_g3strict.json`.

### 1.4 Zone III (t5_zone3)
- Path: `B/t5_zone3/`. Reports: `Z5_REPORT.md` (final), `T5_REPORT.md` (r1 stop).
- Driver: `cd B/t5_zone3 && env PYTHONHASHSEED=0 nice -n 19 python3 -B ramp_U.py` (once; output unchanged by v33), then `env PYTHONHASHSEED=0 nice -n 19 python3 -B zone3.py --ucert B/t4_cert/out/U_cert_for_zone3_v33.json` (argv as pinned in the output's args block). Times: not recorded in the logs (both are quick; Z5_REPORT says two consecutive runs are byte-identical).

| file | sha256 |
|---|---|
| `B/t5_zone3/out/ramp_U.json` | 07906d61bcf67daf263a18d08553aedef49daf85130fa87eae90c01654a93812 |
| `B/t5_zone3/out/zone3_interface.json` (**the zone III / H6 result**) | 0d6338d0f8580cd0cc3c2b0f829b412f8ca0be37f077173023fe447fe382c754 |

- Result [cert]: box [253, 379.5] (b0, b1, b_z) = (0.1003, 0.4197, 0.0143), thinnest x1.24 on b1; Delta^R_2 = 0.01148; n_III = (0.0714, 0.327, 3.04); cone family t* x2.30; **H6 minimum x1.86**; provisional = false.
- Superseded, kept: `zone3_interface.provisional.json`, `zone3_interface.t2base.json`, `z1_box_prelim.json`.

### 1.5 Proof text (t3_proof)
- Path: `B/t3_proof/lemmaU.tex` (495 lines, **not compiled**), `T3_REPORT.md`, `constants_needed.json` (76 rows), `check_balance.py`.
- Check: `cd B/t3_proof && python3 -I check_balance.py lemmaU.tex` (reported PASS; rerun it before use). Time: under 1 s.

| file | sha256 |
|---|---|
| `B/t3_proof/lemmaU.tex` | 19b0e21dd50472ef09419b6b9461fbf7d1fc1a4f8bbe1fb3d237d892eb052052 |
| `B/t3_proof/constants_needed.json` | 6b075bcc3ca1f3e469b162f7c1b9df1ed9431135d97a8b5564d191596645bb56 |

- Items G5, G6, G7 of T3_REPORT are open for the task-8 review.

### 1.6 Statements (r2, Fable)
- Path: `B/t1_statements/`: `STATEMENTS_r2.tex/.md`, `CONTRACTS_r2.md`, `HYPOTHESIS_O_r2.md`, `R2_DECISIONS.md` (r1 versions kept, superseded).
- No driver (text). `r2_probe/` holds the probe scripts (labelled [probe]/[est]).

| file | sha256 |
|---|---|
| `STATEMENTS_r2.md` | d2dade7390df0e6729c833a3cda9ddd1a4014ac78f031ce7ad00b6940a52f198 |
| `STATEMENTS_r2.tex` | cf3446482e96fa929d6dca8d7dd52c2bf9fe1c013d9da6c777ae7cdedb80fa38 |
| `CONTRACTS_r2.md` | a7b12d8070ce15a9b07fa2793f5812990ab3f2cacdf648d7e15810643f4577d5 |
| `HYPOTHESIS_O_r2.md` | a3f19ae7d167f6b6359269207b0bd843d7b3f08922d36469ef66f3ff9ef39ad7 |
| `R2_DECISIONS.md` | 3ef29f4fcc72ccafdbaf770cdafd483487a54764d7d7bce87718c2ff245bea65 |

### 1.7 Step-4 v33 stages (S4/p2/v33, out_v33)
- Path: `S4/p2/v33/` (scripts, `T6_REPORT.md`), outputs `S4/p2/out_v33/`. The closed `S4/p2/out/FINAL_manifest.json` (sha256 e01cb4da71ba825318b5395f053595d6e980e84c1ac286b9d72ff1e4376a8c94, 34 stages, 523 paths) is unchanged.
- Driver: `cd S4/p2 && nice -n 19 python3 -I -B v33/final_tail_v33.py [--from <stage index>]` (resumable; progress in `out_v33/FINAL_manifest_v33.partial.json`; first it verifies the v3.2 manifest, last it re-verifies it). Total about 280 s (stage times: 0.1 / 252 / 8 / 0 / 0 / 1 / 1 / 19 / 0 / 0 s).
- Stages: 0 h2f_ift, 1 h2f_balls -> h2f_win, 2 t9_h2f, 3 h2f_camp, 4 h2t_v33, 5 t9_cone_v33, 6 h4loc_v33, 7 t9_lemmaO_v33, 8 ratio_v33, 9 notes_v33.

| file | sha256 |
|---|---|
| `out_v33/FINAL_manifest_v33.json` | 52af64e6c5fb923a51429ce7b98bf536502df6338868897cde8f9c1a978e3ca6 |
| `out_v33/h2f_result.json` | 479408f7fdfefc15d6922975e27655c80979bedfd55fd3e23073aa41bc833631 |
| `out_v33/t9_h2f.json` | a069403c733ce788c3eb6d31aa3b0ff8653bf39571cab6486d13985821962abc |
| `out_v33/h2t_v33.json` | 0b6a56a5ce435f9d66cdd31b77aa5ab61e2c4439f3f9ff65e65c0d08e939ef3d |
| `out_v33/t9_cone_v33.json` | f21326e6055c5e3e8ce29825318ca97cb4d2bf0f7249aafa293d1114f596f8eb |
| `out_v33/h4loc_v33.json` | 5bee108abec6ef00bd34d78d71dde93a521648036efb19e6ac3adb3eb847c41c |
| `out_v33/t9_lemmaO_v33.json` | ff0299be6595b844f092a2024794b492e5da82571435c3b73703bbef18a57c9f |
| `out_v33/ratio_v33.json` | 5b5456a9225ed30898ba28d93749077cf3bf799abd175fac34e8c4b8cf8f9d62 |
| `out_v33/notes_v33.json` | c304171ad617d8c56903b24cc71984ef38062bea2aef7be702201d39f18e1917 |

- Result [cert]: (H2f) x3.19 / x4.98 / x4.11; tightened (H2) x124.8 / x2125 / x641; fixed-L cylinders: slope x2.64, h'' x3.04, jets x12.2+, H4 x1119, H1 graph holds on all 8 rows, H7 t* x2.30. 152 paths, 23 scripts, 52 reused pinned v3.2 files (each checked against FINAL_manifest.json).
- The final manifest's relation to FINAL_manifest.json is in `STEP4_RESULT_v33_addendum.md`.

### 1.8 Documents of this directory (t7_assembly)
`V33_INDEX.md` (this file), `STEP4_RESULT_v33_addendum.md`, `main_v3.tex`, `main_v3_CHANGES.md`, `KITS_TODO.md`, `erratum/erratum.tex` (written by another agent in parallel; not checked by this file).

## 2. Rerun everything, in dependency order
All commands from a shell, WSL, tree at `/home/claude`. Expected total: about 20-25 min plus the verifier (about 8 min of closure runs plus tail_cert and mod_sups, not timed).
```
# 1. verifier fin12f4
cd /home/claude/s1spec/step3nl_v33/verifier && sh run_v33.sh
cd /home/claude/s1spec/step3nl_v33/verifier && python3 mod_sups.py fin12f4 modsups_fin12f4   # same environment as run_v33.sh
# 2. repaired Lemma S(A)
bash /home/claude/s1spec/writeup/lemmaU/bulk/t2_A_v33/base_ibp/run_chain_v33.sh
bash /home/claude/s1spec/writeup/lemmaU/bulk/t2_A_v33/run_t2.sh
# 3. zone III ramp (independent of 1 and 2)
cd /home/claude/s1spec/writeup/lemmaU/bulk/t5_zone3 && env PYTHONHASHSEED=0 nice -n 19 python3 -B ramp_U.py
# 4. Lemma U, first pass: u1..u5 and the assembly (this yields the U-fine Bhat)
cd /home/claude/s1spec/writeup/lemmaU/bulk/t4_cert
R="env PYTHONHASHSEED=0 nice -n 19 python3 -B"
$R u1_ghat.py; $R u2_K.py --fresh; $R u34_sweeps.py --case primary; $R u5_enum.py
$R assemble_U.py --A /home/claude/s1spec/writeup/lemmaU/bulk/t2_A_v33/out/A_interface.json --out U_cert_v33.json
# 5. zone III interface, on the U-fine Bhat (U_cert_for_zone3_v33.json is the U-fine Bhat in the key zone3.py reads)
cd /home/claude/s1spec/writeup/lemmaU/bulk/t5_zone3 && env PYTHONHASHSEED=0 nice -n 19 python3 -B zone3.py --ucert /home/claude/s1spec/writeup/lemmaU/bulk/t4_cert/out/U_cert_for_zone3_v33.json
# 6. Lemma U, second pass (reads zone3_interface.json: Delta^R, n_III, b1)
cd /home/claude/s1spec/writeup/lemmaU/bulk/t4_cert && $R assemble_U.py --A /home/claude/s1spec/writeup/lemmaU/bulk/t2_A_v33/out/A_interface.json --out U_cert_v33.json
# 7. proof-text checks
cd /home/claude/s1spec/writeup/lemmaU/bulk/t3_proof && python3 -I check_balance.py lemmaU.tex
# 8. step-4 v33 stages
cd /home/claude/s1spec/step4/p2 && nice -n 19 python3 -I -B v33/final_tail_v33.py
# 9. compare with the hashes of section 1 (modulo `args`/timing fields where an output carries them)
```
Expected agreement: outputs of steps 2, 4, 5, 6 have no timing fields and should be byte-identical to the hashes above. Step-4 outputs are identical apart from `cpu`/`seconds` fields (T6 section 6); `t9_lemmaO_v33` agrees with v3.2 only to about 1e-9 on this machine's python-flint 0.9.0, so a rerun on another flint version may differ in the last digits.

## 3. Caveats found while building this index (not mathematics; for the PI and task 8)
1. **Circular step in the chain.** `assemble_U.py` reads `zone3_interface.json`, and `zone3.py` reads the U-fine Bhat from `U_cert_for_zone3_v33.json`. The latter's args block says `"coordinator: copy ufine.Bhat from U_cert_v33.json"`: the file was **filled by hand (copy), not by a driver**. The assemble script's own `U_cert_for_zone3.json` (non-v33 name) is not what zone3 read. The two-pass order in section 2 reproduces the intended fixed point, but the copy step is manual; check that the rerun reproduces sha 8e58d9d1... before relying on it.
2. **Stale checksum file.** `B/t5_zone3/out/SHA256SUMS` lists `zone3_interface.json` as b643987b..., but the file is now 0d6338d0... (rerun on the v33 U-fine Bhat); the SUMS file also lists four scripts with relative names that are not in `out/`. `sha256sum -c` fails there. The hashes in this index are the current ones. The SUMS of t2_A_v33/out and t4_cert/out verify.
3. **T3_REPORT hashes.** `T3_REPORT.md` quotes A = ec2c22b6... and U = f2543e70... (before the G1/G3 fixes). Current: A = 5237a96a..., U = 42cf9617.... `lemmaU.tex` already cites the current prefixes.
4. `U_cert_v33.json: supersolution.beta_A` was 1.25 before G1; G13_FIXES.md records the relabel to 2.
5. `t4_cert/run_all.sh` is the pre-v33 chain (default A path `t2_A/out/A_interface.json`, writes `U_cert.json`). Use section 2.

## Coordinator fixes (after assembly)
- Stale `t5_zone3/out/SHA256SUMS` regenerated (JSON outputs only; the old file is kept as SHA256SUMS_stale). Scripts are pinned in `t5_zone3/SHA256SUMS_scripts`. Both verify.
- The hand-made U-fine copy is replaced by `t4_cert/make_ucert_for_zone3.py <U_cert_v33.json> <out>`. The Bhat is identical; the old file is kept as U_cert_for_zone3_v33_handmade.json.
- zone3 rerun on the scripted file: all numbers identical; new sha256 0d6338d0f8580cd0cc3c2b0f829b412f8ca0be37f077173023fe447fe382c754.
- **Circular pin, by design:** assemble_U reads zone3_interface.json, and zone3.py reads U-fine from U_cert. The numbers are a fixed point (a second pass changes nothing), but the file hashes cannot all pin each other. The pinned pair of record is U_cert_v33.json 42cf9617… with zone3_interface.json 0d6338d0…; U_cert's args block pins the previous, numerically identical zone3 file.
- Rerun order: assemble_U (pass 1) → make_ucert_for_zone3 → zone3 → assemble_U (pass 2; numbers must match pass 1).
- T3_REPORT and lemmaU.tex now cite the current hash prefixes.


## 4. Current chain (Nov 14 sync): B0, M14, M9. This section supersedes sections 1.2-1.5, 1.7 and 2

Paths: `B` = `/home/claude/s1spec/writeup/lemmaU/bulk`, `S4` = `/home/claude/s1spec/step4`. Labels as above. Margin = threshold / bound, rounded down; bounds rounded up. Sources: `B/fix1/b0/B0_REPORT.md`, `B/fix1/M14_REPORT.md`, `S4/p2/v33/M9_REPORT.md`, `B/fix1/lemmaU_rev/LEMMAU_REV.md`, `B/fix1/b12_text/B12_REPORT.md`, `B/t8_review/REVIEW_v33.md`.

### 4.1 New and changed components
| component | path | what it is |
|---|---|---|
| `t2_A_v34` (B0) | `B/fix1/b0/t2_A_v34/` (copy of `t2_A_v33`, which is untouched) | (A) in the variable f - U with the certified profile difference d = (f-U) - w' inserted at four interfaces; M3 (typed ramp constants read from fin12f4) and G-A1 (integration rule R-sup, with `A0_Rsup`, `Fg` in the export) applied in M14. Driver `run_t2.sh` (38 s), byte-identical on two runs. Includes the 4 x N2_ref assembly `asm_design_N2ref4.json`. |
| B0 step-4 stages | `S4/p2/v33/b0_dprof.py`, `h2u_v33.py` (stages 10, 11) | `b0_dprof`: certified bounds on d on [14, 60] (and ramp/core sups); `h2u_v33`: tightened (H2) for f - U_1 on [14, 60], x27.34 / x562.8 / x600.9 on the M9 arb row. |
| M9 step-4 stages | `S4/p2/v33/` (stages 12-18) | `m9_balls.py` -> `win_arb.py`, `h3_arb.py`, `h2t_v33_arb.py`, `h2u_v33_arb.py`, `h2f_camp_arb.py`, `t9_lemmaO_v33_arb.py`, `h4far_v33.py`. arb replacements of the float cell half-widths and sums, the typed centre map, and the far-family H4 at every centre of hand-off radius >= 490 (min x15.58). The arb rows supersede the float rows. |
| M14 fixes | `B/t4_cert/assemble_U.py`, `u2_K.py`, `B/t5_zone3/ramp_U.py`, `zone3.py` | M1 (K^o_jj with the primary Y_next), M2 (coupling term restored), M3 (typed constants), M4 (u1/u34 regenerated). |
| Lemma U text | `B/fix1/lemmaU_rev/lemmaU_v34.tex`, `make_constants_v34.py`, `check_keys_v34.py`, `constants_needed_v34.json` | revision v34 of `t3_proof/lemmaU.tex` (which is unchanged). Check: `python3 -I check_keys_v34.py` (PASS: all sources pinned, 132 cited keys covered, 184 quoted numbers outward). |
| (A) text | `B/fix1/b12_text/` | `lemmaS_A_v33.tex`, `continuation_v33.tex`, `step3_v33_main.tex`, `numbers_A_v33.tex`, `b12_preamble.tex`, `B12_REPORT.md`. |
| sync report | `B/fix1/SYNC_REPORT.md` | this documentation sync. |

### 4.2 Final certified outputs (sha256)
| file | sha256 |
|---|---|
| `B/fix1/b0/t2_A_v34/out/A_interface.json` | 34832f434b7abd7e2dde14be0ee217d99d64236d07983f03863b4239a3ac4226 |
| `B/fix1/b0/t2_A_v34/out/asm_design_N2ref4.json` | 6ad175c172a65a2ae507e493748808056553d4202c26d7b53db418a6bf3b5cf8 |
| `B/fix1/b0/t2_A_v34/out/SHA256SUMS` | f8cff94a0fbc6febd94ac38453b0a63a4a7f0e8057030494799efbfaef6d9275 |
| `B/t4_cert/out/u1_ghat.json` | 5a06d24cae8591f508f9f7499c7d2b4355e06cd5b9d8d2419a9a0108c5638466 |
| `B/t4_cert/out/u2_K.json` | 10ae3c383f7098284e5ae16f4355db9a2645f64531f62c6565c5824213f73261 |
| `B/t4_cert/out/u34_primary.json` | abaf73724b84a0845b549d74a32050e5efafc8f246360ea2d8c388dd0dd70333 |
| `B/t4_cert/out/U_cert_v34.json` (**the Lemma U result**) | 271430d7ffc0306fdb8b7d4942d1e98262f02ff1abde6706fd46334eb99b2507 |
| `B/t4_cert/out/U_cert_for_zone3_v34.json` | 9feebbb2b672879f332858caf2ac209c29978a3b1df913e7bd780a4127e4c1cd |
| `B/t5_zone3/out/zone3_interface.json` (**zone III / H6**) | a381b1f22248a3adc62fc085086e170d75d1388dd086ea6fc1aff9de7ec704c0 |
| `B/t5_zone3/out/ramp_U.json` | 55317ce5a31e8e3c3b3c4b86e85ef9569bb6d5e99190b68b0f21e4cba1a0c0de |
| `S4/p2/out_v33/FINAL_manifest_v33.json` (**19 stages**) | 8f6f7f95d880f37542482989e76c12a8952392075ed2a573b2892561047d8eb8 |
| `S4/p2/out_v33/b0_dprof.json` | 7fdf378c5f59f098263aec8f39bfff7f365e0d25cbbcd154c3f1010bd4a2cd8e |
| `S4/p2/out_v33/h2u_v33.json` | 009c9d1197f0cd6d7bf79ae1cf6cead0f444d74c5bcb096546f8c6b374d36c2a |
| `S4/p2/out_v33/h2u_v33_arb.json` | 65143fc87c7d5891565a4d69117b0a25961920fb020738e85613655acbc65749 |
| `S4/p2/out_v33/h2f_result_arb.json` | cac17df46c05e54c630b5f42d8c164de32298eda1123c414480249cd9c7260a6 |
| `S4/p2/out_v33/h3_arb.json` | ee1a216e38b23eefcbc969bd19f7debd50a088a6f538fe0e5d5aea372226048a |
| `S4/p2/out_v33/h2t_v33_arb.json` | 59d1adc43668726d22dcb96ea47e588f9b95a4812b2c4f155b0557b1d2f16f48 |
| `S4/p2/out_v33/t9_lemmaO_v33_arb.json` | ad09d6261387c78ded85bb7af15b760139022d58ffbd7f62011066ab05e1772a |
| `S4/p2/out_v33/h4far_v33.json` | c6209a42201cc7c3df98a13e6045af495b8dd8e9e5e4d3cc6f360c23de6577a9 |
| `B/fix1/lemmaU_rev/lemmaU_v34.tex` (after the sync edits) | dc178ba091ebac0a083dddbb79e6ad1309e1ab7f20ecc344d3e29abb623e5a16 |

The verifier (tag fin12f4, `/home/claude/s1spec/step3nl_v33/V33_VERIFIER.md`) is unchanged from section 1.1. The closed `S4/p2/out/FINAL_manifest.json` (e01cb4da...) is unchanged. The previous manifests are kept: `FINAL_manifest_v33_pre_m9.json` (8202d2c6..., 12 stages), and the 10-stage 52af64e6... (not kept as a file).

### 4.3 The 19-stage manifest
`FINAL_manifest_v33.json`: 19 stages, 281 paths, 35 scripts, 85 reused v3.2 pins, all_items_ok = true (`M9_REPORT.md`, Records). Stages: 0 h2f_ift, 1 h2f_balls -> h2f_win, 2 t9_h2f, 3 h2f_camp, 4 h2t_v33, 5 t9_cone_v33, 6 h4loc_v33, 7 t9_lemmaO_v33, 8 ratio_v33, 9 notes_v33, 10 b0_dprof, 11 h2u_v33, 12 m9_balls -> win_arb, 13 h3_arb, 14 h2t_v33_arb, 15 h2u_v33_arb, 16 h2f_camp_arb, 17 t9_lemmaO_v33_arb, 18 h4far_v33. Stage 9 (notes_v33) was rerun in the M9 pass because zone3_interface.json changed (M14); stages 10-11 reran byte-identically. Two M9 passes (each with an emptied `out_v33/m9/`) agree byte for byte (65 m9 JSON files, 9 stage outputs).

Current step-4 margins on the v33 items (arb rows; M9_REPORT):
- (H2f): x3.1885 / x4.9834 / x4.1067 (also x1.59 against the earlier design value 5e-7 at level 0; the canonical (H2f)_0 is 1e-6);
- tightened (H2), w': x124.78 / x2124.96 / x641.53 (the report rounds to x124.79, x2124.96, x641.54; rounded down here: x124.7, x2124, x641);
- tightened (H2), f - U_1 (`h2u_v33_arb`): x27.34 / x562.8 / x600.9;
- (H3): x5.700 / x79.26 / x26.50 / x7.799, binding Y-bar limit 1.575679e-11 (x15.75 against Y-bar = 1e-12);
- far-family Lemma O H4 at every centre of hand-off radius >= 490: minimum x15.58 (stage h4far_v33; the M9 report prints x15.59, which is the value rounded up);
- fixed-L cone and the rest: unchanged (slope x2.641, K2 x3.037, n3 x27.38, H4 x1119, H7 t* x2.302).

### 4.4 Rerun, in dependency order (as run by M14 and M9)
```
# 0. verifier fin12f4 and the step-4 stages 0-8, 10, 11 (b0_dprof, h2u_v33), as in section 2 / 4.3, in index order.
#    (t2_A_v34 reads S4/p2/out_v33/b0_dprof.json.)
# 1. (A), B0 + M3 + G-A1: ~38 s
bash /home/claude/s1spec/writeup/lemmaU/bulk/fix1/b0/t2_A_v34/run_t2.sh
# 2. Lemma U, first pieces
cd B/t4_cert; R="env PYTHONHASHSEED=0 nice -n 19 python3 -B"
$R u1_ghat.py; $R u34_sweeps.py --case primary; $R u2_K.py --fresh      # u1/u34 byte-identical to the previous values (M4)
# (ramp_U.py was changed by M1 and reads u1_ghat.json; its position in the order is not recorded in M14_REPORT: [TODO: source])
# 3. assemble pass 1 (reads the previous zone3 interface), then make_ucert_for_zone3, zone3, assemble pass 2, zone3 again, assemble pass 3
$R assemble_U.py --A B/fix1/b0/t2_A_v34/out/A_interface.json --out U_cert_v34.json     # pass 1
$R make_ucert_for_zone3.py <U_cert_v34.json> U_cert_for_zone3_v34.json
cd B/t5_zone3 && env PYTHONHASHSEED=0 nice -n 19 python3 -B zone3.py --ucert B/t4_cert/out/U_cert_for_zone3_v34.json
$R assemble_U.py ...                                                                      # pass 2 (the assemble pass 2 of the earlier index, now with v34 inputs)
# zone3 again; assemble pass 3. M14: pass 3 == pass 2 apart from args (fixed point).
# 4. Lemma U text checks
cd B/fix1/lemmaU_rev && python3 -I make_constants_v34.py && nice -n 19 python3 -I check_keys_v34.py   # both PASS
python3 -I B/t3_proof/check_balance.py lemmaU_v34.tex
# 5. step-4 v33 stages 9-18 (stage 9 re-pinned to the final zone3 interface; stages 12-18 are the M9 arb stages)
cd S4/p2 && nice -n 19 python3 -I -B v33/final_tail_v33.py --from 9
```
Exact argv of the assemble and zone3 calls are in the `args` blocks of the outputs. The M14 pass files are `U_cert_v34_pass1.json`, `U_cert_v34_pass2.json`; the pre-M14 outputs are kept as `*_pre_m14*`.

### 4.5 Open items carried by this section (not mathematics)
- `t4_cert/run_all.sh` is the pre-v33 chain; `t4lib.T2` still points at `t2_A` (review minor items).
- `zone3_interface.json handoff_thresholds` still lists (H2f)_0 = 5e-7 with a probe label (M9; the certified value passes both).
- `ramp_U.json` pins the pre-M4 `u1_ghat.json` hash 4b86cc91...; the file is 5a06d24c... and only its args block differs.
- `fix1/b12_text/continuation_v33.tex` still cites `U_cert` 4c3f23bc..., zone3 8152c650... and states the far-family coverage item as open (its gap item 3) and (W3) as subject to G-A1; not edited by the sync (outside its list). See `fix1/SYNC_REPORT.md`.
- Review items still open: independent review of the v3.3 text (task 8); compilation.

### 4.6 Current chain (Nov 14, after re-review N1/N2 and the N35 documentation pass). This section supersedes sections 4.1-4.5 where they differ

Paths: `B` = `/home/claude/s1spec/writeup/lemmaU/bulk`, `S4` = `/home/claude/s1spec/step4`. Margin = threshold / bound, rounded down; bounds rounded up. Sources: `B/fix2/n12/N12_REPORT.md` (N1, N2), `B/fix2/N35_REPORT.md` (this pass), `B/t8_review/REREVIEW_v33.md`.

**New and changed components**
| component | path | what it is |
|---|---|---|
| `n1_chart` (step-4 stage 19, new) | `S4/p2/v33/n1_chart.py` (sha256 b42cca46...), output `S4/p2/out_v33/n1_chart.json` (341ed1f4...) | certified chart term between the verifier's S1-chart variable w' and the ray-graph variable of (A) (re-review N1). 8 s. Replaces the typed 2e-7 allowance. |
| driver with explicit environments (N2) | `S4/p2/v33/final_tail_v33.py` (sha256 92308040...); backup `final_tail_v33_pre_n12.py` | every stage runs with `env_child` + `env_read` recorded; a `--from` run asserts the reused stages ran with the same environment. |
| `t2_A_v35` (N1) | `B/fix2/n12/t2_A_v35/` (copy of `fix1/b0/t2_A_v34`, which is untouched) | (A) with `--chart` (reads `n1_chart.json`) in addition to `--dprof`; N35 adds the assertion "verifier ramp sups <= ramp box (0.52, 0.52, 0.05)" in `rigor/assembly_rep.py` (backup `_pre_n35`). Driver `run_t2.sh` (40 s). |
| `assemble_U.py` | `B/t4_cert/assemble_U.py` | N35: asserts `N2_le_4N2ref` (backup `assemble_U_pre_n35.py`). |
| chart lemma | `B/fix2/n12/chart_lemma.tex`, merged as `B/fix1/b12_text/chart_lemma_body.tex` | `\input` by `lemmaS_A_v33.tex` (section "chart identification"). |
| texts | `B/fix1/b12_text/{lemmaS_A_v33,numbers_A_v33,step3_v33_main,continuation_v33}.tex`, `B/fix1/lemmaU_rev/lemmaU_v34.tex`, `B/t7_assembly/main_v3.tex`, `erratum/{ERRATUM.md,erratum.tex}` | v35 numbers and hashes (N35); backups `*_pre_n35.*`. |
| Lemma U text check | `B/fix1/lemmaU_rev/{make_constants_v34.py, check_keys_v34.py, constants_needed_v34.json}` | read the v35 files; PIN of S4 = a0d348ec...; `check_keys_v34.py` PASS (184 quoted numbers, 134 cited keys covered). |

**Final certified outputs (sha256)**
| file | sha256 |
|---|---|
| `B/fix2/n12/t2_A_v35/out/A_interface.json` | 78fe2f8e38a340ba515fbc33cbab60cdd85231095af6ccd192e5f6707f457e8c |
| `B/fix2/n12/t2_A_v35/out/asm_design_N2ref4.json` (the 4 x N2_ref run) | 46bae149800904875e5b510d02de7843dc791f890d53b9a3e45338e63a19aa16 |
| `B/fix2/n12/t2_A_v35/out/SHA256SUMS` | 7791c6afb85a21064cf4649d0d8da3a89d279e37242e3c24cc27799ad69dd3de |
| `B/t4_cert/out/u1_ghat.json` | 5a06d24cae8591f508f9f7499c7d2b4355e06cd5b9d8d2419a9a0108c5638466 |
| `B/t4_cert/out/u2_K.json` | 10ae3c383f7098284e5ae16f4355db9a2645f64531f62c6565c5824213f73261 |
| `B/t4_cert/out/u34_primary.json` | abaf73724b84a0845b549d74a32050e5efafc8f246360ea2d8c388dd0dd70333 |
| `B/t4_cert/out/U_cert_v35.json` (**the Lemma U result**, pass 2; pass 3 equal apart from args) | 960809fdc18dd67c1770fa5b77b2658081d0c81cb24f1b246f48ca4a583d5f54 |
| `B/t4_cert/out/U_cert_for_zone3_v35.json` | 22c13e4403a1b41ab09ab15a9d8852f153832a97742895757fdb1ce3fd0ca19f |
| `B/t5_zone3/out/zone3_interface.json` (**zone III / H6**) | 45e8d4dcaff2bd75fcea725f80e98625af70a40709fc2cf5d3df78b5becb3ea3 |
| `B/t5_zone3/out/ramp_U.json` | 55317ce5a31e8e3c3b3c4b86e85ef9569bb6d5e99190b68b0f21e4cba1a0c0de |
| `S4/p2/out_v33/FINAL_manifest_v33.json` (**20 stages**, 289 paths, 37 scripts, 85 reused pins) | a0d348ec2dfb84bd3c62482517a0aa26f5d29ab3e4230081384054f395f8e357 |
| `S4/p2/out_v33/n1_chart.json` | 341ed1f442ac64114bea6943cf37466a1515361bfcf838b41a7f7d860bf1ece0 |
| `S4/p2/out_v33/b0_dprof.json` | 7fdf378c5f59f098263aec8f39bfff7f365e0d25cbbcd154c3f1010bd4a2cd8e |
| `S4/p2/out_v33/h2u_v33_arb.json` | 4167b5f33d013da845e6170208a8f2c96168dbec7e3a5d5accfbd4e216ca367a |
| `S4/p2/out_v33/h2f_result_arb.json` | c5477282ef01b389d3868d9d1a49340700569b4e5782d4094e7fc9fcc5fbf2c8 |
| `S4/p2/out_v33/h3_arb.json` | 3bb9a58eada8895c4cada985f1bb3d37a8f73bad32aa5fb35a4a99734f99809a |
| `S4/p2/out_v33/h2t_v33_arb.json` | 55045a5ff942b22eaecee1e0354953c93e6d32dbd3be6b869d422ff375267ab8 |
| `S4/p2/out_v33/t9_lemmaO_v33_arb.json` | ad09d6261387c78ded85bb7af15b760139022d58ffbd7f62011066ab05e1772a |
| `S4/p2/out_v33/h4far_v33.json` | eeac024be8320288d93032ee745f39add7ad92ae2452585ecfaec6b5cfd401c8 |
| `B/fix1/lemmaU_rev/lemmaU_v34.tex` (after N35) | 5131776dbf2f003637c70291266f11873e21065ef1ec8ece7f8828eb230cdee6 |
| `B/fix1/lemmaU_rev/constants_needed_v34.json` (after N35) | 0f3c745546a4dbf9ec2dc12e4e3b5571da001ce121039e719e1e2be23c09d37c |

The closed `S4/p2/out/FINAL_manifest.json` (e01cb4da...) is unchanged. The v34 hashes (section 4.2) are superseded: `A_interface` 34832f43..., `U_cert_v34` 271430d7..., zone3 a381b1f2..., manifest 8f6f7f95... (19 stages); the previous manifest is kept as `FINAL_manifest_v33_pre_n12.json` (= 8f6f7f95...).

**The 20-stage manifest.** Stages 0-18 as in section 4.3; **19 n1_chart** (new). Environment block, `env_child`, `env_read` (55 variables): `N12_REPORT.md` section 6. Passes A and B (clean) agree byte for byte apart from `cpu`/`seconds`; pass C (`--from 9`) was needed because `notes_v33` and `h4far_v33` pin `zone3_interface.json`.

**Rerun, in dependency order** (replaces 4.4)
```
# 0. verifier fin12f4; step-4 stages 0-8, 10, 11 and 19 (n1_chart needs hstail.py + hs_uforms.pkl, lintail2.py, /home/claude/s1_cap/tail.py, the fin12f4 rate table)
cd S4/p2 && nice -n 19 python3 -I -B v33/final_tail_v33.py            # or --from <stage>
# 1. (A), B0 + N1 + M3 + G-A1: ~40 s (reads out_v33/b0_dprof.json and out_v33/n1_chart.json)
bash B/fix2/n12/t2_A_v35/run_t2.sh
# 2. Lemma U, zone III: as 4.4 steps 2-3 with A = B/fix2/n12/t2_A_v35/out/A_interface.json and the v35 file names:
cd B/t4_cert; R="env PYTHONHASHSEED=0 nice -n 19 python3 -B"
$R assemble_U.py --A B/fix2/n12/t2_A_v35/out/A_interface.json --zone3 ../t5_zone3/out/zone3_interface.json --out U_cert_v35.json     # pass N
$R make_ucert_for_zone3.py <U_cert_v35.json> U_cert_for_zone3_v35.json ; (cd ../t5_zone3 && ... zone3.py --ucert ../t4_cert/out/U_cert_for_zone3_v35.json)
# passes 1, 2, 3 are kept as U_cert_v35_pass{1,2,3}.json; pass 3 == pass 2 apart from args (pass 2 pins an intermediate zone3 file, pass 3 the final one)
# 3. step-4 stages 9-19 re-pinned to the final zone3 interface
cd S4/p2 && nice -n 19 python3 -I -B v33/final_tail_v33.py --from 9
# 4. text checks
cd B/fix1/lemmaU_rev && python3 -I make_constants_v34.py && nice -n 19 python3 -I check_keys_v34.py          # both PASS
python3 -I B/t3_proof/check_balance.py <file>.tex
```
`t4_cert/run_all.sh` now defaults its `A_interface` argument to the v35 file (backup `run_all_pre_n35.sh`); it is the single-pass chain and does not pass `--zone3`/`--out`. `t4lib.T2` still points at `t2_A`.

**N35 checks.** (1) `t2_A_v35/run_t2.sh` rerun with the new assertion: all 20 JSON outputs identical to `out/` apart from `args`; the rerun is kept in `t2_A_v35/out_n35rerun/` (its `A_interface.json` is 5ee8d3ed..., differing only by the script hashes in `args`), and `out/` (78fe2f8e...) stays the pinned copy. (2) `assemble_U.py` rerun with the assertion on the rerun interface: `t4_cert/out/U_cert_v35_n35rerun.json` (8afbd403...) equals `U_cert_v35.json` and `U_cert_v35_pass3.json` apart from `args`. (3) `check_keys_v34.py` PASS. (4) `check_balance.py` PASS on every edited .tex (expanded for the `\input` files).

**Open items carried by this section** (not mathematics)
- The step-4 reports `step4/p2/v33/T6_REPORT.md` and `M9_REPORT.md` still describe the 19-stage manifest and print x15.59 (rounded up; JSON 15.5878): not edited, `step4/` is closed to this pass. [TODO: authors to add an addendum there if wanted.]
- `ic_rep.py` (t2_A_v35) now cites the `_arb` outputs in the string it writes to `ic_*.json` and `A_interface.json` (key `IC.step4_H2_range`); the pinned outputs were produced by the earlier script text, so a rerun of `run_t2.sh` would change that one string and therefore every downstream hash of this section. Numbers are unaffected.
- `ramp_U.json` pins the pre-M4 `u1_ghat.json` hash 4b86cc91...; the file is 5a06d24c... and only its args block differs (as in 4.5). `zone3_interface.json handoff_thresholds` still lists (H2f)_0 = 5e-7 (as in 4.5).
- Position of `ramp_U.py` in the rerun order: [TODO: source] (not recorded in M14_REPORT).
- (O3) entries of the erratum table on the v35 chain: [TODO: source], not used by v3.3 (`erratum/ERRATUM.md`).
- Review items still open: independent review of the v3.3 text beyond REREVIEW_v33.md; compilation; the literature statements (P1') [to be confirmed by the authors] and (P3').
