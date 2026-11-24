# STEP 4 RESULT, v3.3 addendum (new hand-off items for step 3 v3.3)

Addendum to `/home/claude/s1spec/step4/STEP4_RESULT.md` (base v2). That file is **not edited**; `STEP4_RESULT.md` stays authoritative for everything not listed here. Nov 13, 2026. Source: `step4/p2/v33/T6_REPORT.md`; hashes in `V33_INDEX.md` section 1.7.

> **SYNC NOTE 2 (Nov 14, 2026, N35).** Section 6 at the end is the final record: 20-stage manifest a0d348ec... with the new stage `n1_chart` and the explicit per-stage environment. It supersedes sections 2 and 5 where they differ.

> **SYNC NOTE 3 (Nov 14, 2026, F34).** Section 7 at the end is the final record: **21-stage manifest 8bc47e83e370fec3e856e9eb61cf0863ec0738320806b294cb2d74772db983ec** (366 paths, 39 scripts) with the new stage `f12_sliver`. It supersedes the 20-stage manifest a0d348ec... of section 6.

> **SYNC NOTE (Nov 14, 2026).** Sections 1-4 below are the Nov 13 text (10-stage manifest 52af64e6...). Section 5 at the end gives the final record: 19-stage manifest 8f6f7f95..., the B0 and M9 stages, and the arb margins. Where section 1 quotes a margin rounded up (x3.19, x4.11, x124.8, x2125), section 5 gives the value rounded down. `STEP4_RESULT.md` itself is still not edited.

Status. Step 3 v3.3 (Lemma U, repaired Lemma S(A); `writeup/lemmaU/bulk/`) needs hand-off data that step 4 v2 did not certify: a datum on the far support [60, 700], a tightened shell threshold, and the fixed-L cone cylinders. All of it is certified below. **Closed step 4 is untouched:** `p2/out/FINAL_manifest.json` (sha256 e01cb4da71ba825318b5395f053595d6e980e84c1ac286b9d72ff1e4376a8c94; 34 stages, 523 paths, 97 upstream, 26 scripts) verified before and after the v33 run, and a size/mtime snapshot of all 10,811 files outside `p2/v33/`, `p2/out_v33/`, `p2/out_v33_pass1/` is identical. Everything is arb with outward rounding, no sampling. Labels: **[cert]**; margin = threshold / bound, rounded down. Not yet independently reviewed (task 8 of the bulk plan).

## 1. New hand-off items (uniform for tau* in slab 180; `p2/v33/`, `p2/out_v33/`)
At Ybar = 1e-12 (certified for Ybar <= 1.7e-12, the level of the n3 boxes).

| item | required | base (arb, slab 180, max over 32 t_c balls) | with error (Ybar = 1e-12) | margin | Ybar limit |
|---|---|---|---|---|---|
| **(H2f)** j = 0, G_0 = int_{ln 60}^{ln 700} (28/y)(f - U_1)^2 ds | 1e-6 | 3.1263e-7 | 3.1363e-7 | **x3.19** | 5.7e-8 |
| (H2f) j = 1 | 5e-6 | 1.0031e-6 | 1.0034e-6 | **x4.98** | 1.7e-8 |
| (H2f) j = 2 | 1e-5 | 2.4302e-6 | 2.4351e-6 | **x4.11** | 1.17e-9 |
| **tightened (H2)** shells [14, 60], j = 0 | 1e-8 | 8.010e-11 | 8.013e-11 | x124.8 | |
| tightened (H2) j = 1 | 1e-6 | 1.335e-10 | 4.71e-10 | x2125 | |
| tightened (H2) j = 2 | 1e-3 | 2.936e-10 | 1.559e-6 | x641 | |
| **fixed-L cone** (H5) slope datum | 0.08 | | 0.03029 (worst row) | **x2.64** | |
| fixed-L cone (H5) sup \|h_0''\| = K2 | 0.03 | | 9.876e-3 | x3.03 | |
| fixed-L cone (H5) \|\|h_0'''\|\|_L2 = N3 | 0.004 | | 1.461e-4 | x27.4 | |
| fixed-L cone H3 jets K2 / K3 / K4 | 0.12 | | 9.88e-3 / 5.96e-5 / 8.36e-4 | x12.2 / x2012 / x143 | |
| fixed-L cone Lemma O H4, int (h_0^(5))^2 | 0.0144 | | 1.2865e-5 (cone r_c 300-310; chart factor 1.0031) | x1119 | |
| fixed-L cone Lemma O H1 (graph) | sup\|h_0\| <= H - C_0 sqrt(t_1) | | graph on all 8 rows: xi <= -179.1 at arc start, >= +170.0 at arc end, XL = 150; sup\|h_0\| <= 0.2162, H <= 0.2185 | holds | |
| fixed-L cone H7 time ratio | t* >= ratio | | ratio <= 1.6000014e-7 (a-bar = 2.19103e-7); t* >= 3.68334e-7 | **x2.30** | |

Notes.
- (H2f) is the new threshold (5e-7, 5e-6, 1e-5) of the r2 statements, with **(H2f)_0 raised to 1e-6** and beta_A = 2 in v3.3 (so (1e-6, 5e-6, 1e-5)). It is the datum of the repaired Lemma S(A) on [60, 700], together with the tightened (H2) = (1e-8, 1e-6, 1e-3). Parts of the level-0 base value: [60, 379.5] 4.84e-8, [379.5, 700] 2.65e-7. The error term changes G by at most 0.3%.
- Method: win.py's cell machinery on x in [45, 725] (fidelity asserted: identical to `win.energies` on [45, 390]; an unwidened ball reproduces the pinned `win_0.json` bit for bit); 32 IFT balls widened by (x_c, lambda) +- 3.28e-11 and alpha +- 1.23e-12; the e-term is summed in L2 over 51 bins of ratio 1.05 (Minkowski), not maximised over sliding windows (T6 section 1).
- Fixed-L cone family: L = 50, X = 3, centres r_c in [300, 380] in 8 rows of width 10, cylinders [r_c - 150, r_c + 150], in L = 50 units. This replaces the unit-L cone family of v3.2 in step 3 v3.3 (H5 of v3.2 is unchanged and still holds for v3.2).
- The thinnest new margins: fixed-L cone H7 t* **x2.30**, slope **x2.64**, (H2f)_0 **x3.19**. The overall thinnest margins of step 4 are unchanged (|alpha| x1.27, ...; STEP4_RESULT section 1).

## 2. The new manifest
- Driver `p2/v33/final_tail_v33.py` -> **`p2/out_v33/FINAL_manifest_v33.json`** (sha256 52af64e6c5fb923a51429ce7b98bf536502df6338868897cde8f9c1a978e3ca6): 10 stages, 152 paths, 23 scripts. Written atomically; resumable with `--from`; progress in `FINAL_manifest_v33.partial.json`; children run `python3 -I -B` under `nice -n 19`.
- Stages: 0 h2f_ift, 1 h2f_balls -> h2f_win, 2 t9_h2f, 3 h2f_camp, 4 h2t_v33, 5 t9_cone_v33, 6 h4loc_v33, 7 t9_lemmaO_v33 (LO_NS = 48), 8 ratio_v33, 9 notes_v33. Total about 280 s.
- **Relation to FINAL_manifest.json.** The new manifest is additional and does not replace it. Its first action verifies FINAL_manifest.json (34 stages, 523 paths, 97 upstream, 26 scripts); its last action re-verifies it. 52 of its inputs are pinned v3.2 files, each checked against the hash that FINAL_manifest.json records for it. No stage reads a v33 output back into a v3.2 stage. A complete v3.3 step-4 record is therefore the pair (FINAL_manifest.json, FINAL_manifest_v33.json).
- Reproduction (T6 section 6): pass 1 kept in `p2/out_v33_pass1/`; pass 2 is the clean run the manifest records. 56 of 74 JSON outputs are byte-identical between passes; the rest differ only in `cpu`/`seconds` or in recorded hashes of per-ball files that differ only in `cpu`. Pass 1's only failure was the H1 end test, see section 4.
- Fidelity of the copies: with v3.2 geometry, `t9_cone_v33.py` and `h4loc_v33.py` reproduce the pinned v3.2 cone rows bit for bit; `t9_lemmaO_v33.py` reproduces the v3.2 sup|h_0|, H and H8 to about 1e-9 relative, not bitwise (likely python-flint 0.9.0 transcendental enclosures).

## 3. Reused v3.2 outputs (not rerun)
| item | pinned output | status |
|---|---|---|
| far family (centres >= 490 in v3.3; a subset of the certified v3.2 family, centres >= 376) | `t9_far_v32.json`, `h4loc_v32.json`, `t9_lemmaO_far.json` (ok) | reused; bases start at r >= 264 >= 150, inside t9_far's region |
| H6 (v3.2) | `t9_h6_v32` | **superseded** by step 3 v3.3's H6 (`zone3_interface.json`, cone_family.H6, minimum x1.86 at the reference Bhat); its pinned output stays in FINAL_manifest.json and is not a v3.3 input |
| (H1) alpha, sqrt(Lambda), (H3) windows and K4, (H4) graph, (H5) far data, t1 interval, all of STEP4_RESULT sections 1-2 | as in STEP4_RESULT | unchanged |
All of these are recorded in `out_v33/notes_v33.json`.

## 4. Arc-end note (H1 of Lemma O; for review)
The v3.2 script `t9_lemmaO.py` computed the arc-end cells but **never tested |xi| > XL there**. The v33 copy asserts it. With v3.2's arc of 1.1 XL the end test fails by a few units (crude enclosures: the width-10 centre row and slab-time motion), so the arc was lengthened to 1.3 XL; cos(theta - theta_c) > 0 is checked on the whole longer arc, so the test is stronger, not weaker. sup|h_0| of about 0.2 reflects the same enclosure widths; Lemma O puts no other condition on H. H8 is informational only (distance >= 0.72 against the old 0.1 rule). **Open question for the review (not decided here): does closed v3.2's cone H1 need an erratum**, because its script did not check the arc ends? The far family's H1 (`t9_lemmaO_far.py`, STEP4_RESULT section 3) is a separate script and does check |xi(+-G)| >= 226.18 > XL.

## 5. Final record after B0 and review M9 (Nov 14, 2026; supersedes section 2 and the margins of section 1)

Sources: `p2/v33/M9_REPORT.md`, `B/fix1/b0/B0_REPORT.md`, `p2/out_v33/FINAL_manifest_v33.json`. Closed step 4 (`p2/out/FINAL_manifest.json`, e01cb4da..., 34 stages, 523 paths) is unchanged and was verified before, at the start and end of every pass, and after. All values arb with outward rounding [cert]; margin = threshold / bound, rounded down.

**Manifest.** `p2/out_v33/FINAL_manifest_v33.json`, sha256 `8f6f7f95d880f37542482989e76c12a8952392075ed2a573b2892561047d8eb8`: **19 stages**, 281 paths, 35 scripts, 85 reused v3.2 pins, all_items_ok = true. Backup of the previous manifest: `FINAL_manifest_v33_pre_m9.json` (8202d2c6..., 12 stages); driver backup `final_tail_v33_pre_m9.py`. Records 0-8 are byte-identical to the earlier run (run with `--from 9`); stage 9 (`notes_v33`) was rerun because `zone3_interface.json` changed (M14); stages 10-11 (B0) reran byte-identically. Two M9 passes agree byte for byte.

Stages: 0-9 as in section 2; **10 b0_dprof** (profile difference d on [14, 60]); **11 h2u_v33** ((H2) for f - U_1); **12 m9_balls -> win_arb** (arb window/shell/(H2f) base sums, 32 t_c balls x 2 ball sets); **13 h3_arb**; **14 h2t_v33_arb**; **15 h2u_v33_arb**; **16 h2f_camp_arb**; **17 t9_lemmaO_v33_arb** (Lemma O H1 with a certified centre map, byte-identical rows); **18 h4far_v33** (far-family H4).

**What M9 changed.** The float cell half-width (`float(mid hh)*1.0000001`; rad(hh)/hh reaches 1.54e-3, so the factor did not cover it) and the float sums were replaced by arb cells and arb sums. The arb values are at most 7.5e-5 (relative) above the float values, so every margin is unchanged to that accuracy. The arb rows supersede the float rows.

| item | threshold | with error [cert] | margin |
|---|---|---|---|
| (H2f) j = 0 | 1e-6 | 3.1363e-7 (3.136242e-7) | x3.1885 (also x1.59 against 5e-7) |
| (H2f) j = 1 | 5e-6 | 1.0034e-6 | x4.9834 |
| (H2f) j = 2 | 1e-5 | 2.4351e-6 (2.435022e-6) | x4.1067 |
| tightened (H2), w', j = 0 | 1e-8 | 8.0135e-11 | x124.78 |
| tightened (H2), w', j = 1 | 1e-6 | 4.7060e-10 | x2124.96 |
| tightened (H2), w', j = 2 | 1e-3 | 1.5588e-6 | x641.53 |
| **tightened (H2), f - U_1 (B0)** j = 0 | 1e-8 | 3.6575e-10 | **x27.34** |
| j = 1 | 1e-6 | 1.7768e-9 | x562.8 |
| j = 2 | 1e-3 | 1.6641e-6 | x600.9 |
| (H3) windows j = 0 / 1 / 2 / 3 | 2e-7, 1e-5, 1e-5, 1e-5 | 3.5087e-8, 1.2616e-7, 3.7724e-7, 1.2822e-6 | x5.700 / x79.26 / x26.50 / x7.799 |
| (H3) Y-bar limit, binding (E_3) | Y-bar = 1e-12 | 1.575679e-11 | x15.75 |
| far-family Lemma O H4, **every** centre at hand-off radius >= 490 (37 rows) | int (h_0^(5))^2 <= M^2 = 0.0225 | worst row (centres 471-584 in base arclength) 1.4435e-3 | **x15.58** |

Notes.
- The margins in the table are rounded down from the JSON (`M9_H2f`, `M9_H2_tightened_fU1`, `M9_H3_H2v32`, `M9_H2_tightened`, `M9_H4_far`). The M9 report prints x3.1885 / x4.983 / x4.107, x124.79, x601.0, x15.59 in places; those are rounded to nearest or up.
- (H2f)_0: the canonical value is **1e-6** (v3.3, beta_A = 2). `STATEMENTS_r2.tex`, `R2_DECISIONS.md` and `zone3_interface.json` still say 5e-7 (stale; the certified 3.1362e-7 passes both).
- **B0 (profile difference).** d = (f - U) - w' = alpha(psi3^' - 1.3056 y^p) + alpha^2 h2^'. At |alpha| <= alpha(t_1) = 1.17387e-4: E_k[d] on [14, 60] = (1.035e-10, 4.186e-10, 1.720e-9, 7.133e-9). (H2) for f - U_1 follows by Minkowski from the w' shells with their error transfer and E[d] (`h2u_v33`). The threshold (1e-8, 1e-6, 1e-3) is unchanged.
- **Far family.** h4far_v33 closes the far-family coverage caveat of v3.2 (certified only at centres 376 and 602; centres in (490, 602) by an unwritten monotonicity argument). It uses h4loc's per-cylinder bound unchanged, slack s = 3 (the excess term reaches 2.157 near the far end) and y_e = min(y_0, y_l/2) in the dissipation term. With h4loc's geometry it reproduces the pinned `h4loc_v32.json` totals for centres 376 and 602 bit for bit.
- The v3.2 H1 script for the cone did not test |xi| > XL at the arc ends (section 4): unchanged; still an open erratum question.
- Not independently reviewed (task 8).

## 6. Final record after N1 and N2 (Nov 14, 2026; supersedes sections 2 and 5)

Sources: `B/fix2/n12/N12_REPORT.md` (`B` = `/home/claude/s1spec/writeup/lemmaU/bulk`), `p2/out_v33/FINAL_manifest_v33.json`. Closed step 4 (`p2/out/FINAL_manifest.json`, e01cb4da..., 34 stages, 523 paths) is unchanged and was verified before the work, at the start and end of every driver pass, and after. Values arb with outward rounding [cert]; margin = threshold / bound, rounded down.

**Manifest.** `p2/out_v33/FINAL_manifest_v33.json`, sha256 `a0d348ec2dfb84bd3c62482517a0aa26f5d29ab3e4230081384054f395f8e357`: **20 stages**, 289 paths, 37 scripts, 85 reused v3.2 pins, all_items_ok = true. Driver `p2/v33/final_tail_v33.py` (sha256 92308040...; backup `final_tail_v33_pre_n12.py`). It supersedes `FINAL_manifest_v33_pre_n12.json` (= 8f6f7f95..., 19 stages); the old outputs are kept in `out_v33/pre_n12/`, the old logs in `out_v33/logs_pre_n12/`.

**New stage 19, `n1_chart`** (`p2/v33/n1_chart.py`, sha256 b42cca46...; output `out_v33/n1_chart.json`, sha256 341ed1f4...; 8 s, about 100 MB). It certifies the chart term between the verifier's S1-chart variable w' and the ray-graph variable of (A) (re-review N1). Inputs: the step-1 tail envelopes (`step3/hstail.py` with `hs_uforms.pkl`, `lin/lintail2.py`, `/home/claude/s1_cap/tail.py`), the fin12f4 rate table, `c4fun.py`, `psi3fun.py`, `h2fun.py`, `t9h/s1chart.py`. Results (N1_chart row; ramp [14 - 2^-19, 28 + 2^-19], core [28, 60]):
- sup |chi^add_k| on the ramp <= (2.087e-8, 4.444e-7, 9.910e-6, 2.307e-4), at sigma = 14; on the core <= (3.95e-19, 2.15e-17, 1.19e-15, 6.63e-14);
- ramp: beta22 <= 2.41e-9, beta32 <= 1.82e-7, beta33 <= 3.20e-9, |D beta22| <= 8.50e-8, |eps| <= 7.83e-10, weight transfer rho_meas <= 1.05e-9; |e| <= 4.39e-10 < 2^-19 (asserted);
- consequence for (A): every margin moves by at most 0.06%; no threshold fails (`N12_REPORT.md` section 4).
The datum is unaffected: (H2) and (H2f) are computed by step 4 directly in the ray chart (`win.py`, `win_arb.py`), and f - U_1 = (f - P) + d exactly. The +-2e-7 ball that `win.py:70` adds at level 0 only widens those enclosures and is harmless; its comment is misleading, and `win.py` is a closed file.

**Environment record (N2).** Every stage runs with `env_child` = {PATH, LC_ALL = C.UTF-8, PYTHONHASHSEED = 0, PYTHONDONTWRITEBYTECODE = 1} plus the stage's own variables; nothing is inherited. `env_read` pins the values of 55 variables found by grepping the stages and their imports (S1TT, S1PROF and FORM_MULT are unset for every stage). A `--from` run asserts that each reused stage ran with the same environment. Children run with `-I`, so PYTHONHASHSEED is recorded but ignored by the interpreter; byte-identity rests on the two clean passes.

**Margins** (arb rows; unchanged by N1 to the quoted digits): (H2f) x3.1885 / x4.9834 / x4.1067; tightened (H2) for w' x124.78 / x2124.96 / x641.53, for f - U_1 x27.34 / x562.8 / x600.9; (H3) x15.75 in Y-bar (binding E_3: 1.575679e-11); far-family Lemma O H4 at every centre >= 490: x15.58 (15.5878); fixed-L cone slope x2.641, K2 x3.037, H4 x1119, H7 t* x2.302. Rounded down from the JSON; the M9 and N12 reports print some of these rounded to nearest or up.

**Reproduction.** Passes A and B (clean: the 64 win_arb and 32 h2f_win sub-runs recomputed) agree on all 289 paths apart from `cpu`, `seconds` and the "(Ns)" timings in the logs. Against the pre-N2 outputs the only differences are `cpu`, the reused-run lists and the input hashes that follow from them.

Not independently reviewed beyond `t8_review/REREVIEW_v33.md`.

## 7. Stage f12_sliver (review items F1, F2; Nov 14, F34 note)

**Manifest.** `p2/out_v33/FINAL_manifest_v33.json`, sha256 `8bc47e83e370fec3e856e9eb61cf0863ec0738320806b294cb2d74772db983ec`: **21 stages**, 366 paths, 39 scripts, 85 reused v3.2 pins, all_items_ok = true. The previous manifest is kept as `FINAL_manifest_v33_pre_f12.json` (= a0d348ec..., 20 stages). New stage 20, `f12_sliver` (`p2/v33/f12_sliver.py`, 124f13be...; output `out_v33/f12/f12_sliver.json`, 1fbc0f80...): certifies F1 (the sliver sups G_i on the extended domain) and F2 (sup|w''| <= sI2' = 0.015456178888449044 on the sliver, M1' = 28 sI2' <= 0.4328) and recomputes every dependent claim from pinned inputs: 326 of 326 checks pass, no quoted number changes. Source: `p2/v33/F12_REPORT.md`. Two clean passes agree byte for byte apart from `seconds`.
