# v3.3 step-3 file map (Lemma U and the repaired Lemma S(A))

Task 7b, Nov 13 2026 (Sonnet). Bookkeeping only: no number was changed. Hashes were computed with `sha256sum` when this file was written (Nov 13 2026). Labels: **[cert]** = python-flint arb with outward rounding; **[paper]** = hand proof; **[hyp]** = hypothesis. Status of the whole set: **not yet independently reviewed (task 8), LaTeX not compiled**.

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
