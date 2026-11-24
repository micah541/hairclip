# Task 6 (v3.3): new step-4 hand-off items. All pass [cert]. (H2f) margins ×3.19 / ×4.98 / ×4.11; fixed-L cone thinnest ×2.30 (H7), ×2.64 (slope)

Opus, Nov 13 2026.
- **Closed step 4 untouched.** `p2/out/FINAL_manifest.json` verifies before and after: 34 stages, 523 paths, 97 upstream, 26 scripts, no mismatch. Its sha256 is `e01cb4da…8c94`.
- The step4 tree outside `p2/v33/`, `p2/out_v33/` and `p2/out_v33_pass1/` is unchanged. A size/mtime snapshot of all 10,811 files is identical before and after.
- Driver: `p2/v33/final_tail_v33.py` → `p2/out_v33/FINAL_manifest_v33.json` (sha256 `52af64e6…3ca6`): 10 stages, 152 paths, 23 scripts.
  - 52 inputs are reused pinned v3.2 files. Each is checked against the hash that FINAL_manifest.json records for it.
  - The manifest is written atomically, and a run can resume with `--from` (progress goes to `FINAL_manifest_v33.partial.json`).
  - Every child runs as `python3 -I -B` under `nice -n 19`, one process at a time, peak RSS about 140 MB.
- Labels: everything below is **[cert]** (python-flint arb, outward rounding). margin = threshold / bound, rounded down.

## 1. (H2f): G_j = ∫_{ln 60}^{ln 700} (28/y)(∂_s^j(f − U₁))² ds at t₁, uniformly over slab 180, true curve

| j | base (max over 32 t_c balls) | √(error) at Ȳ = 1e-12 | with error | threshold | margin | Ȳ limit |
|---|---|---|---|---|---|---|
| 0 | 3.1263e-7 | 8.94e-7 | 3.1363e-7 | 1e-6 | **×3.19** | 5.7e-8 |
| 1 | 1.0031e-6 | 1.29e-7 | 1.0034e-6 | 5e-6 | **×4.98** | 1.7e-8 |
| 2 | 2.4302e-6 | 1.56e-6 | 2.4351e-6 | 1e-5 | **×4.11** | 1.17e-9 |

- The base values agree with the probe (3.13e-7, 1.00e-6, 2.43e-6). Parts: [60, 379.5] (4.84e-8, 1.79e-7, 5.44e-7); [379.5, 700] (2.65e-7, 8.25e-7, 1.89e-6).
- The error term changes G by ≤ 0.3%. The stop trigger (×1.2) is far away.
- The Ȳ limits are as camp computes them. The result is certified for Ȳ ≤ YBOX = 1.7e-12, the n3 boxes' level; the actual Ȳ = 1e-12.

How it was computed (stages 0–3):
- **h2f_ift.** The 32 slab-180 IFT balls are widened so that the base value also covers the true curve's hand-off parameters:
  - (x_c, λ) ± 3.28e-11, t9.py's true-curve frame shift;
  - α ± C_α·YBOX = 1.23e-12.
  - (H3) did not include an α shift in U₁; the effect there is ~1e-10 in √E, negligible.
- **h2f_balls → h2f_win.** This is win.py's cell machinery, imported read only, with geometry to x = 730 and cells on x ∈ [45, 725].
  - In-run fidelity: the copied cell routine is asserted identical to `win.energies` on win's own range [45, 390].
  - Scratch check: an unwidened ball 0 reproduces the pinned `win_0.json` windows bit for bit.
  - Coverage asserted: cos φ ≥ 0.992 on every cell, so y increases monotonically, and the cells reach below y = 60 and above y = 700.
  - The sum is in arb, with cell widths rounded up.
- **t9_h2f.** This is t9.py's KCH/e-term block, copied, on y_h ∈ [60, 700] for levels 0–2.
  - In scratch, with t9.py's range it reproduces the pinned KCH = 1.187776680575588 and c0_chain bit for bit.
  - Results: KCH_f = 1.1568; c0 = (8.90e-3, 7.30e-2, 0.161) per √F₀; c0tip = (0.302, 0.311, 0.321) per Ȳ.
  - **The one change:** the e-term is summed in L² over 51 bins of ratio 1.05 (Minkowski), not maximised over sliding windows.
  - Inputs valid on [60, 700]: the n3 boxes are asserted up to y_ℓ ≈ 3.6e5; the |e| core bound holds for every y ≥ 28; (SE) acts only for x ≥ 0.6ℓ ≈ 2.1e5.
- **h2f_camp.** camp's item-3 formula, with A_m and B_m from `amp_A1_w180.json`.

## 2. Tightened (H2) = (1e-8, 1e-6, 1e-3), shells [14, 60] (stage 4, h2t_v33)

This was listed in CONTRACTS_r2 for task 6. It re-assembles pinned outputs only, so it is cheap.

| j | base | with error | margin |
|---|---|---|---|
| 0 | 8.010e-11 | 8.013e-11 | ×124.8 |
| 1 | 1.335e-10 | 4.71e-10 | ×2125 |
| 2 | 2.936e-10 | 1.559e-6 | ×641 |

## 3. Fixed-L cone family

Geometry: L = 50, X = 3, centres r_c ∈ [300, 380] in 8 rows of width 10, cylinders [r_c − 150, r_c + 150], in L = 50 units on the whole cylinder.

| item (stage) | bound (worst row) | threshold | margin |
|---|---|---|---|
| (H5) slope (t9_cone_v33) | 0.03029 | 0.08 | **×2.64** |
| (H5) sup\|h₀″\| = K2 | 9.876e-3 | 0.03 | ×3.04 |
| (H5) ‖h₀‴‖_L² = N3 | 1.461e-4 | 0.004 | ×27.4 |
| H3 jets K2 / K3 / K4 | 9.88e-3 / 5.96e-5 / 8.36e-4 | 0.12 | ×12.2 / ×2012 / ×143 |
| H4 ∫(h₀⁽⁵⁾)² (h4loc_v33) | 1.286e-5 (r_c 300–310; chart factor 1.0031) | 0.0144 | ×1119 |
| H1 graph (t9_lemmaO_v33) | graph on all 8 rows: ξ ≤ −179.1 at arc start, ≥ +170.0 at arc end vs XL = 150. sup\|h₀\| ≤ 0.2162, H = sup\|h₀\| + 0.00227 ≤ 0.2185 | sup\|h₀\| ≤ H − C₀√t₁ (C₀√t₁ = 0.00227) | holds |
| H7 time ratio (ratio_v33) | 50⁻⁴/(1 − 4ā) ≤ 1.6000014e-7, with ā = 2.19103e-7 from the v33 mod_sups | ≤ 1.6001e-7; t* ≥ 3.68334e-7 | t*/ratio **×2.30** |

- **Fidelity of the copies:** in scratch with the v3.2 geometry, `t9_cone_v33.py` and `h4loc_v33.py` reproduce the pinned v3.2 cone rows bit for bit.
  - `t9_lemmaO_v33.py` reproduces the v3.2 sup|h₀|, H and H8 values to about 1e-9 relative, but not bitwise.
  - The likely cause is transcendental arb enclosures under this machine's python-flint 0.9.0.
- **H1 change (please check in review):**
  - t9_lemmaO.py computed the arc-end cells but never tested |ξ| > XL there. The v33 copy asserts it.
  - With the 1.1·XL arc of v3.2 the end test fails by a few units, because of crude enclosures: the width-10 centre row and the slab-time motion.
  - The arc was therefore lengthened to 1.3·XL. cos(θ − θ_c) > 0 is checked on the whole longer arc, so the test is stronger, not weaker.
  - sup|h₀| ≈ 0.2 reflects the same enclosure widths; Lemma O puts no other condition on H. H8 is reported for information only: distance ≥ 0.72, against the old 0.1 rule.
  - The centre rows use t9_lemmaO.py's arclength-to-r convention, arcs [rlo + 2.0, rlo + 12.6].
- C₀√t₁ = 0.00227 is lemmaO_v32_consts.json's cone X = 3 deviation C₀ at the ratio 1.6001e-7, 0.00226524, rounded up.

## 4. Not rerun (stage 9, notes_v33.json)
- **Far family (centres ≥ 490)** is a subset of the certified v3.2 family (centres ≥ 376). Bases start at r ≥ 264 ≥ 150, inside t9_far's region.
  - Pinned and reused: `t9_far_v32.json`, `h4loc_v32.json`, `t9_lemmaO_far.json` (ok).
- **t9_h6_v32** is superseded by task 5's H6 (`zone3_interface.json` cone_family.H6, min margin ×1.86 at B̂_ref). Its pinned output stays in FINAL_manifest.json and is not a v3.3 input.
- (H3), α₁, √Λ, the (H5) far data and the other v3.2 items are unchanged.

## 5. Stage list (FINAL_manifest_v33.json)
| # | stage | s |
|---|---|---|
| 0 | h2f_ift (widened parameter balls) | 0.1 |
| 1 | h2f_balls → h2f_win (32 balls, base sums) | 252 |
| 2 | t9_h2f (chart constants on [60, 700]) | 8 |
| 3 | h2f_camp ((H2f) with error) | 0 |
| 4 | h2t_v33 (tightened H2) | 0 |
| 5 | t9_cone_v33 (fixed-L jets, (H5), H3) | 1 |
| 6 | h4loc_v33 (H4) | 1 |
| 7 | t9_lemmaO_v33 (H1 graph) | 19 |
| 8 | ratio_v33 (H7) | 0 |
| 9 | notes_v33 | 0 |

## 6. Reproduction and records
- Pass 1 is kept as `p2/out_v33_pass1/`. Pass 2 is the clean full run that the manifest records.
- 56 of 74 JSON outputs are byte-identical. The rest differ only in `cpu`/`seconds`, or in recorded sha256s of per-ball files that themselves differ only in `cpu`.
  - Pass 1 had one failure, the H1 end test at 1.1·XL. Section 3 has the fix. Stages 7–9 were rerun with `--from 7`.
- Every output has an args block: argv, env, and input sha256s, including its own script.
- `p2/v33/verify_v32_manifest.py` is a read-only checker for the v3.2 manifest.
- Housekeeping: my early test runs, made without `-B`, wrote bytecode caches. Specifically, 15 `.pyc` files under step4 and 10 under step3/step3nl/s1_tight. All were removed by literal path, and the step4 snapshot confirms it is identical. Every run since uses `-B`.
- Not done: the kits and the STEP4_RESULT/main updates (task 7).


---
**Status note (Nov 14, 2026, F34; appended, body unchanged).** Superseded counts: manifest now 21 stages, 8bc47e83...; see F12_REPORT.md
