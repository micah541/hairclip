# What the reproduction kits need for v3.3

Task 7b, Nov 13 2026. **Nothing was built and no zip was created.** This is a list of files, sizes and commands for the PI, who decides whether to build and run the kits. Sizes are bytes from `ls -l` today (Nov 13). Paths: `B` = `/home/claude/s1spec/writeup/lemmaU/bulk`, `S4` = `/home/claude/s1spec/step4`, `V` = `/home/claude/s1spec/step3nl_v33`, `P1` = `/home/claude/s1spec/writeup/lemmaU/phase1`, `P2A` = `/home/claude/s1spec/writeup/lemmaU/phase2a`. Hashes of the outputs the kits must reproduce: `V33_INDEX.md`.

> **SYNC NOTE (Nov 14, 2026).** Sections A and B below list the Nov 13 file set. Section D at the end gives the final list: the chain is now t2_A_v34 (B0, M14), U_cert_v34, the 19-stage manifest 8f6f7f95... (including 8 M9 scripts, `b0_dprof.py` and `h2u_v33.py`). Hashes: `V33_INDEX.md` section 4. Nothing was built.

> **SYNC NOTE 2 (Nov 14, 2026, N35).** Section E at the end is the final list: chain v35 (`fix2/n12/t2_A_v35`, `U_cert_v35`, the 20-stage manifest `a0d348ec...`, the new stage `n1_chart`; see section F for the 21-stage manifest of F34). It supersedes sections A, B and D where they differ. Nothing was built.

> **SYNC NOTE 3 (Nov 14, 2026, F34).** Section F at the end: the step-4 v33 manifest is now 21 stages (`8bc47e83...`, 366 paths, 39 scripts); kit B also needs `f12_sliver.py`.

Present state of the kits (not changed): `kit_step3.zip` (1.46 MB; step 3 as of v3.1: parts smoke, verifier tag fin11, ibp, lemmaS, lemmaO; 91 PASS on the Windows run), `kit_step4a.zip` and `kit_step4a_v3.zip` (Gamma_0 and P1; no change needed for v3.3), `kit_step4b.zip` (38.4 MB, built from the 34-stage manifest, **not yet run end to end anywhere**; unpacked build tree `/home/claude/s1spec/kit_step4b`, builder `/home/claude/s1spec/kit_step4b_build/build.py`, `REFRESH.md`). `unzip -l kit_step4b.zip` contains no file with "v33" in its name, so it has none of the v3.3 stages.

## A. kit_step3 (verifier + Lemma S) -> kit_step3 v3 (parts: verifier fin12f4, lemmaS_v33, lemmaU, zone3)
Existing kit parts to change:

| part | change | files |
|---|---|---|
| verifier | environment gets `FLUX0=1.64e-20 FLUX2=8.0e-20`; table tag `fin12f4`; expected values: slopes and F0 x4 (see V33_VERIFIER.md), closure 1.5e-4 about 215 s, 1.25e-4 about 195 s (this machine); `mod_sups.py fin12f4` to regenerate sup\|(a,b)\| = 2.19103e-07, sup\|beta\| = 3.33811e-09 | replace `step3nl/verifier/constants.py`-driven run by the v33 driver `V/verifier/run_v33.sh` (environment of `run_final.sh` + FLUX); ship `rate_table_fin12f4_a0.00015.json` (13,614 B) as the reference table and `modsups_fin12f4/final_alpha_points_5e-8.json` (1,666 B) |
| lemmaS | replace by the repaired Lemma S(A): `B/t2_A_v33/` chain | see list below |
| ibp, lemmaO | unchanged | none |

New files for the repaired Lemma S(A) part (copy of `B/t2_A_v33/`, 9.8 MB as a tree; the minimal subset is):
- `B/t2_A_v33/rigor/` : `a3_signs.py` 2,673; `arlib.py` 3,606; `assembly_rep.py` 20,072; `export_t2.py` 11,989; `hjw_rep.py` 4,862; `ic_rep.py` 2,271; `outer_L2.py` 15,980; `ramp_rep.py` 16,391; `weights_rep.py` 12,692; `weights_t2.py` 24,067 (+ `outer_rep.py`).
- `B/t2_A_v33/base_ibp/` : `run_F10p_ibp.py`, `lemmaS_F1.py` 8,232, `lemmaS_prop1_sym.py` 3,708, `lemmaS_F10p_M0.52.json` 4,831, `lemmaS_cflux.json` 1,235, `rigor/{capture_F2.py, weights_exact.py 19,981, assembly_arb.py 23,807, assembly_arb.json 11,728, cap_F10p.json 6,627, tables_exact.pkl 2,839}`, `run_chain_v33.sh`; `B/t2_A_v33/base_v31H2/rigor/assembly_arb.json` 11,744.
- `/home/claude/s1spec/writeup/lemmaS/rigor/assembly_arb.json` 11,745 (input to the first chain step; also `tables_exact.pkl`).
- Drivers: `base_ibp/run_chain_v33.sh` (about 2.5 min) then `run_t2.sh` (about 40 s). Absolute paths are hard-coded under `/home/claude/` (the kit patches them through its `KIT_ROOT` marker as in kit_step3; **every script of this part needs that patch**, as does `arlib.py` ROOT).
- Reference to compare: `B/t2_A_v33/out/*.json` (19 files, about 1.2 MB) and `out/SHA256SUMS` (348076020c...); the interface `A_interface.json` is 91,639 B.

New part **lemmaU** (Lemma U certificates; one process each, about 5 min in total):
- `B/t4_cert/` : `t4lib.py` 3,008; `u1_ghat.py` 11,244; `u2_K.py` 13,423; `u34_sweeps.py` 10,743; `u5_enum.py` 7,348; `rampphi.py` 2,608; `assemble_U.py` 21,995. Reference outputs: `u1_ghat.json` 26,695; `u2_K.json` 9,127; `u34_primary.json` 32,071; `u5_enum.json`; `U_cert_v33.json`.
- External inputs they read (ship them): `P1/ghat/ghat_arb.py` 18,685; `P1/ghat/weight_consts_arb.py` 9,841; `P1/ramp/ulib.py` 15,760; `P2A/U_ramp/ghat_design_arb.json` 10,161; `/home/claude/s1spec/writeup/lemmaS/rigor/tables_exact.pkl` 2,839; `B/t2_A/rigor/{arlib.py, ramp_rep.py 16,391}` and `B/t2_A/base_ibp/rigor/weights_exact.py` 19,981 (note: `u5_enum` imports these from `t2_A`, not `t2_A_v33`; verify that the copies are identical or repoint, see V33_INDEX caveats).
- Command order: `u1_ghat.py`, `u2_K.py --fresh`, `u34_sweeps.py --case primary`, `u5_enum.py`, `assemble_U.py --A <t2_A_v33/out/A_interface.json> --out U_cert_v33.json`; then the zone3 part; then `assemble_U.py` again (see V33_INDEX section 2). The hand copy `U_cert_for_zone3_v33.json` (600 B) should be generated by a script before shipping (it is a copy of `ufine.Bhat`).

New part **zone3**:
- `B/t5_zone3/` : `ramp_U.py` 14,006; `zone3.py` 22,783; `z3lib.py` 5,963; reference `out/ramp_U.json` 6,454, `out/zone3_interface.json` 18,712.
- Inputs: `/home/claude/s1spec/writeup/lemmaO_v32/lemmaO_v32_consts.json` 9,493; `.../lemmaR_consts.json` 1,982; `B/t1_statements/r2_probe/ramp_out_253_379.5_box0.15_0.52.json` 430; the same `P1/ghat`, `P1/ramp`, `P2A/U_ramp` files as above.
- Commands: `ramp_U.py`; `zone3.py --ucert <t4_cert/out/U_cert_for_zone3_v33.json>`.

Also needed in the kit:
- `check.py` extended with the v3.3 expectations: the margins and sha256 of `V33_INDEX.md`. Compare JSON outputs byte for byte except the `args` blocks (they hold absolute paths and `cwd`), and expect python-flint 0.9.0.
- Pin `python-flint==0.9.0` (the version of all reference runs; `t9_lemmaO_v33` already differs from v3.2 by about 1e-9 across flint versions, T6 section 3).
- README: new parts, time per part (verifier about 8 min + tail_cert; A_v33 about 3.5 min; lemmaU about 5 min; zone3 not timed here; whole new set about 20 min on this machine, measure on the first run).
- `lemmaU.tex` and `main_v3.tex` are text, not run by a kit; the kit may ship `check_balance.py` (`B/t3_proof`) and `check_balance_main.py` (this directory) as the only checks on them.
- Not required by v3.3 but check: the old `lemmaS` part (fin11 table) can stay as a historical v3.1 check or be dropped.

Estimated size of new material for kit_step3 v3: about 1.5 MB of sources and references (the largest items are the `asm_*.json` references, 7 x about 105 kB, and `weights_t2*.json`, 2 x about 49 kB; they can be left out and compared through `SHA256SUMS`).

## B. kit_step4b (step-4 tail + the v33 stages) -> kit_step4b v2
Add stages 17-26 (after the existing 16): the ten v33 stages, run after `tail2` has produced the closed-chain outputs they read.

New files (copy into `step4/p2/v33/`; all under `S4/p2/v33/`): `final_tail_v33.py` 13,314; `h2f_ift.py` 3,047; `h2f_balls.py` 2,594; `h2f_win.py` 5,824; `t9_h2f.py` 11,210; `h2f_camp.py` 4,290; `h2t_v33.py` 2,713; `t9_cone_v33.py` 6,971; `h4loc_v33.py` 16,057; `t9_lemmaO_v33.py` 9,562; `ratio_v33.py` 1,998; `notes_v33.py` 2,812; `verify_v32_manifest.py` 1,399 (13 files, about 84 kB).
Reference outputs to ship as `reference/` (for comparison): `S4/p2/out_v33/` (572 kB: `FINAL_manifest_v33.json` 52,806 B, `h2f_result.json`, `t9_h2f.json` 16,669, `h2t_v33.json`, `t9_cone_v33.json`, `h4loc_v33.json`, `t9_lemmaO_v33.json`, `ratio_v33.json`, `notes_v33.json`, `h2f/` per-ball files, `logs/`) and an extended `reference/outputs.json` (sha256 of the 152 paths).

Inputs the v33 stages read that kit_step4b must already have or must add:
- Already produced by the kit's own chain (verify the names exist in the kit tree): `S4/p2/out/{h1slab/ift_*.json (32 files, about 450 B each), h1slab/win_summary.json, t9_cone_kch.json 10,541, t9_far_v32.json, t9_lemmaO_far.json, t9_h1_clam_s180.json, ycomp_final.json 41,534}`, `S4/p2/t6c/{residuals.json 1,226,642; jets/jets_177..180.npz about 251 kB each}`, `S4/p2/t8se/constants/{amp_A1_w180.json 94,169; n3_b1.json 3,648,471}`, `S4/t7/constants.json` 1,695,200, `S4/t1/polished_v2/slab_180.pkl` 353,880, `S4/p2/{curvjets.py, t9chart.py, t9h/basegeo.py, t9h/win.py}`.
- **Not in the kit** (add): `V/verifier/modsups_fin12f4/final_alpha_points_5e-8.json` 1,666 B (ā = 2.19103e-7 for the H7 ratio), `/home/claude/s1spec/writeup/lemmaO_v32/lemmaO_v32_consts.json` 9,493 B (C_0 sqrt(t_1)). (Both are listed as inputs in the stage args blocks; the kit's `inputs/MANIFEST.sha256` must pin them.)
- The v33 manifest reads the **v3.2 manifest** (`S4/p2/out/FINAL_manifest.json`, 523 paths); in the kit it is `reference/final_manifest.json`: `final_tail_v33.py` and `verify_v32_manifest.py` must be pointed at the kit's own copy of the manifest and at kit-relative paths (the scripts use absolute `/home/claude/...` paths; apply the `[kit patch]` convention, no `-I` flag so that `PYTHONHASHSEED=0` applies, no `os.nice`).
- Run time: stage 1 (`h2f_balls` -> `h2f_win`, 32 balls) is 252 s here and is the only long stage; the whole v33 tail is about 280 s here. On Windows with `--jobs 6` the 32 balls could be parallelised (the reference runs them serially); the builder decides.
- Commands: `run_all.py --jobs 6` (existing stages), then the new stages via `final_tail_v33.py [--from N]` or as stages 17-26 of `run_all.py`; `check.py` then verifies the (H2f), tightened (H2) and fixed-L cone margins of `STEP4_RESULT_v33_addendum.md` section 1 and compares the outputs with the reference.
- Caveat already known to the builder: kit B's H1 points come from the new `h1_points.py`, while the shipped pre-u4 point enclosures still feed the derivative balls (README of kit_step4b, "Differences"). The v33 stages `h2f_ift` read the `h1slab/ift_*.json` files, which are derived from those balls: **watch this stage**; the (H2f) base values in the kit run may differ in low digits from the reference (3.1263e-7, 1.0031e-6, 2.4302e-6 with margins x3.19 / x4.98 / x4.11 have room, but the byte-identity claim will not hold for the `h2f/` files if the balls differ).
- The kit tree currently has size 84 MB (zip 38 MB); the v33 additions are under 1 MB.

## C. Decisions for the PI
1. Build kit_step3 v3 now, or wait for the task-8 review (a review may change numbers and thus the references).
2. Build kit_step4b v2 only after the first end-to-end run of the present kit B (it has never been run), or add the v33 stages first so that one run covers both.
3. The two manual or fragile points to settle before shipping: the hand-copied `U_cert_for_zone3_v33.json` and the `t2_A` versus `t2_A_v33` imports in `u5_enum.py` (V33_INDEX section 3).


## D. Final kit file lists (Nov 14 sync; supersedes the file lists of A and B where they differ)
Paths: `B` = `/home/claude/s1spec/writeup/lemmaU/bulk`, `S4` = `/home/claude/s1spec/step4`. Sizes in bytes from `ls -l` at the sync. Nothing was built and no zip was created.

### D.1 kit_step3 v3: (A) is `t2_A_v34`, not `t2_A_v33`
- Replace the `B/t2_A_v33/` part by the copy `B/fix1/b0/t2_A_v34/` (driver `run_t2.sh`, 38 s; `rigor/` scripts incl. `assembly_rep.py` 91dbd2be..., `export_t2.py` a395f048...; reference `out/A_interface.json` 34832f43..., `out/asm_design_N2ref4.json` 6ad175c1..., `out/SHA256SUMS`). The base chain `base_ibp/`, `base_v31H2/` is as in section A. Its extra input is the step-4 output `S4/p2/out_v33/b0_dprof.json` (7fdf378c...), so the step-4 v33 stages 10 must run first, or the kit must ship that JSON as an input.
- Lemma U part: `B/t4_cert/` as in section A, plus `make_ucert_for_zone3.py` (replaces the hand copy), with the M14 versions of `assemble_U.py` (d14baeea...), `u2_K.py` (059a2734...); reference outputs `u1_ghat.json` 5a06d24c..., `u2_K.json` 10ae3c38..., `u34_primary.json` abaf7372..., `U_cert_v34.json` 271430d7..., `U_cert_for_zone3_v34.json` 9feebbb2.... Run order including **assemble pass 2** and pass 3: `V33_INDEX.md` section 4.4. Still not run by a kit: `u5_enum` (not rerun in M14).
- Zone III part: `B/t5_zone3/` with the M14 `ramp_U.py` (0df7ce69...) and `zone3.py` (b9f23489...); references `ramp_U.json` 55317ce5..., `zone3_interface.json` a381b1f2.... zone3 now also reads the verifier file (a-bar) and `U_cert_v34.json` (c_dU): ship `V/verifier/modsups_fin12f4/final_alpha_points_5e-8.json`.
- Text checks: `B/fix1/lemmaU_rev/{make_constants_v34.py, check_keys_v34.py, constants_needed_v34.json, lemmaU_v34.tex}` (the check reads the JSON files above by absolute path; PASS expected) and `B/t3_proof/check_balance.py`.

### D.2 kit_step4b v2: 19 stages, not 10
Stages 17-35 in `run_all.py` numbering (after the existing 16), or via `final_tail_v33.py [--from N]` (18,875 B now). New scripts compared with section B (all in `S4/p2/v33/`):

| script | bytes | stage |
|---|---|---|
| `b0_dprof.py` | 7,613 | 10 (inputs: `step3nl/agent_int/v3/c4fun.py` 7,680, `step3nl/psi3fun.py`, `step3nl/h2fun.py`, `S4/p2/t9h/s1chart.py` 6,615; the exact input list is in the stage's `args` block) |
| `h2u_v33.py` | 3,001 | 11 |
| `m9_balls.py` | 5,607 | 12 |
| `win_arb.py` | 8,934 | 12 |
| `h3_arb.py` | 5,226 | 13 |
| `h2t_v33_arb.py` | 3,072 | 14 |
| `h2u_v33_arb.py` | 3,169 | 15 |
| `h2f_camp_arb.py` | 4,700 | 16 |
| `t9_lemmaO_v33_arb.py` | 12,547 | 17 |
| `h4far_v33.py` | 9,831 | 18 |

That is the **8 M9 scripts** (`m9_balls`, `win_arb`, `h3_arb`, `h2t_v33_arb`, `h2u_v33_arb`, `h2f_camp_arb`, `t9_lemmaO_v33_arb`, `h4far_v33`) plus `b0_dprof.py` and `h2u_v33.py`. `verify_v32_manifest.py` and the scripts of section B are unchanged. The reference directory is now `S4/p2/out_v33/` (about 1.4 MB: `FINAL_manifest_v33.json` 8f6f7f95..., the stage outputs, `m9/` with 129 files, `h2f/`, `logs/`); `reference/outputs.json` must list the 281 paths of the new manifest, not 152.
- Run time: the M9 stages add a long block (win_arb over 32 balls x 2 ball sets); the M9 report gives no total time here: [TODO: source].
- The M9 stages read, besides the section-B inputs, `S4/p2/out_v33/h2f/ift_M.json` (the widened (H2f) balls), `S4/p2/out/h1slab/ift_M.json`, `S4/p2/out_v33/b0_dprof.json` and the zone III interface (stage 9 pin): `inputs/MANIFEST.sha256` must pin them.
- The H1-points caveat of section B (kit B's `h1_points.py` versus the pre-u4 point enclosures) applies to the arb stages as well: watch stages 0-1 and 12.

### D.3 Documents and decisions
- The decisions of section C stand. Item 1 of C: the task-8 review is still pending.
- Reference values for `check.py` are in `V33_INDEX.md` section 4.3 and `STEP4_RESULT_v33_addendum.md` section 5 (arb margins, rounded down).


## E. Final kit file lists (Nov 14, after N1/N2; N35 documentation pass; supersedes sections A, B, D where they differ)
Paths: `B` = `/home/claude/s1spec/writeup/lemmaU/bulk`, `S4` = `/home/claude/s1spec/step4`. Sizes in bytes from `ls -l`, Nov 14. Nothing was built and no zip was created. Hashes: `V33_INDEX.md` section 4.6.

### E.1 kit_step4b v2: 20 stages
- **New stage 19, `n1_chart`** (after the 19 stages of D.2): `S4/p2/v33/n1_chart.py` (16,601 B, sha256 b42cca46...). Its inputs, which the kit must ship and pin in `inputs/MANIFEST.sha256`:
  - `/home/claude/s1spec/step3/hstail.py` (8,406 B, ab86d9fa...) with `/home/claude/s1spec/step3/hs_uforms.pkl` (4,310 B, 226db25a...);
  - `/home/claude/s1spec/lin/lintail2.py` (20,491 B, a7231c8d...);
  - `/home/claude/s1_cap/tail.py` (4,737 B, 43a8d545...; note the path is `s1_cap/` directly under `/home/claude`, not under `s1spec/`; `hstail.py` puts `/home/claude/s1_cap` and `/home/claude/s1spec/lin` on `sys.path`; the kit copies exist as `kit_step4b/s1_cap/tail.py`, `kit_step4b/lin/lintail2.py`, `kit_step4b/step3/hstail.py`, `hs_uforms.pkl`);
  - the fin12f4 rate table `/home/claude/s1spec/step3nl_v33/verifier/rate_table_fin12f4_a0.00015.json` (13,614 B, b20e8715...);
  - `step3nl/agent_int/v3/c4fun.py`, `step3nl/psi3fun.py`, `step3nl/h2fun.py`, `S4/p2/t9h/s1chart.py` (as for `b0_dprof`, D.2).
- **New driver:** `S4/p2/v33/final_tail_v33.py` (22,551 B, sha256 92308040...; D.2 quoted 18,875 B for the earlier version). It builds each stage's environment from scratch and records it; the kit's `run_all.py` must either call it or reproduce `env_child`/`env_read` (`N12_REPORT.md` section 6; note `-I` makes PYTHONHASHSEED ineffective in the children).
- **Reference outputs:** `S4/p2/out_v33/` with `FINAL_manifest_v33.json` (a0d348ec...; 289 paths, 37 scripts), `n1_chart.json` (5,550 B, 341ed1f4...) and the other stage outputs of V33_INDEX 4.6; `reference/outputs.json` must list the 289 paths (not 281 or 152).
- The H1-points caveat of section B (kit B's `h1_points.py` versus the pre-u4 point enclosures) and the long `win_arb` block of D.2 are unchanged; run time of `n1_chart` is 8 s; run time of the whole chain is not recorded: [TODO: source].
- Kit B (`kit_step4b.zip`) has still not been run end to end anywhere.

### E.2 kit_step3 v3: (A) is `t2_A_v35`
- Replace the `B/fix1/b0/t2_A_v34/` part of D.1 by the copy `B/fix2/n12/t2_A_v35/` (driver `run_t2.sh`, 40 s): `rigor/assembly_rep.py` (4aef3de6..., with the N35 ramp assertion; the pinned outputs were produced with 77d4abe6..., identical numbers), `export_t2.py` (067e8c4d...), `arlib.py` (ROOT patched), `ic_rep.py` (59b66e1a...; its output string now cites the `_arb` files, see V33_INDEX 4.6 open items), reference `out/A_interface.json` 78fe2f8e..., `out/asm_design_N2ref4.json` 46bae149..., `out/SHA256SUMS` 7791c6af.... Extra inputs: `S4/p2/out_v33/b0_dprof.json` (7fdf378c...) **and `S4/p2/out_v33/n1_chart.json` (341ed1f4...)**, so stages 10 and 19 must run first, or the kit ships both JSONs.
- Lemma U part: as D.1 with `assemble_U.py` (6ee50b4d..., N35 assertion `N2_le_4N2ref`; run order with passes 1-3 in V33_INDEX 4.6), references `U_cert_v35.json` 960809fd..., `U_cert_for_zone3_v35.json` 22c13e44..., `u1_ghat.json`, `u2_K.json`, `u34_primary.json` unchanged (5a06d24c..., 10ae3c38..., abaf7372...). `run_all.sh` now defaults to the v35 interface (backup `run_all_pre_n35.sh`).
- Zone III part: `zone3_interface.json` 45e8d4dc... (only the H6 entries differ from D.1's a381b1f2...); `ramp_U.json` 55317ce5... unchanged.
- Text checks: `fix1/lemmaU_rev/{make_constants_v34.py, check_keys_v34.py, constants_needed_v34.json, lemmaU_v34.tex}` read the v35 files (PASS), `B/t3_proof/check_balance.py`.

### F. kit_step4b v3: 21 stages (F34)
- **New stage 20, `f12_sliver`** (after the 20 stages of E.1): `S4/p2/v33/f12_sliver.py` (35,550 B, sha256 124f13be...). Its 67 inputs are existing pins (all read-only); the kit must ship them and pin them in `inputs/MANIFEST.sha256` (the stage file list is in the driver). Output: `S4/p2/out_v33/f12/` (25 JSON files plus `logs/f12_sliver.log`; the logs are not pinned). Run time of the stage: 28 s. [TODO: source for the exact list of kit-side input files beyond those of E.1; take it from the `inputs` block of `out_v33/f12/f12_sliver.json`.]
- **Driver:** `S4/p2/v33/final_tail_v33.py` (sha256 57c20a44...; backup `final_tail_v33_pre_f12.py`); it now has 21 stages.
- **Reference outputs:** `S4/p2/out_v33/` with `FINAL_manifest_v33.json` (**8bc47e83...**; 366 paths, 39 scripts); `reference/outputs.json` must list the 366 paths (not 289, 281 or 152). The 20-stage manifest a0d348ec... is superseded.
- The other points of E.1 are unchanged. Kit B has still not been run end to end anywhere.
