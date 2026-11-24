# Review item M9 (v3.3): certified replacements for the step-4 hand-off items. All four items are closed [cert]. No threshold fails and no margin moves by more than 7.5e-5 relative.

Opus, Nov 14 2026. Source: `writeup/lemmaU/bulk/t8_review/REVIEW_v33.md` M9.
- Labels: **[cert]** means python-flint 0.9.0 arb at 128 bits with outward rounding.
- Margin means threshold / bound, rounded down.
- New stages 12–18 in `p2/out_v33/FINAL_manifest_v33.json`.

## Summary

| Item | Stage (script) | Result |
|---|---|---|
| 1. Float cell half-widths and float energy sums | `m9_balls.py` → `win_arb.py`, then `h3_arb.py`, `h2t_v33_arb.py`, `h2u_v33_arb.py`, `h2f_camp_arb.py` | The arb values are at most 7.5e-5 (relative) above the float values. Every margin is unchanged to within 7.5e-5 relative. |
| 2. Typed centre map (`t9_lemmaO_v33.py:98–99`) | `t9_lemmaO_v33_arb.py` | Certified. On all 8 rows, every cell that can hold a centre of hand-off radius in [r_lo, r_lo + 10] is among the typed centre cells. The H1 rows are byte-identical to `t9_lemmaO_v33.json`. |
| 3. Far-family H4 | `h4far_v33.py` | Certified for **every** far centre at hand-off radius ≥ 490, all the way to the far end, so no cap is needed. 37 rows; min margin **×15.59** (centres 471–584 in base arclength). |
| 4. (H2f)₀ | (from item 1) | The canonical threshold is **1e-6** (v3.3, β_A = 2). The certified 3.1362e-7 passes it at ×3.19, and also passes 5e-7 at ×1.59. |

## 1. Cell half-widths and energy sums

### What was wrong
- **The cell half-width.** `win.py:42, 57` and `h2f_win.py:32` use `float(mid hh)·1.0000001`.
  - Measured [cert]: rad(hh)/hh reaches **1.54e-3** on a t_c ball of width 1/16 slab. The element lengths grow with the rescaling.
  - So the 1e-7 factor did **not** cover the element length. The reviewer's concern was real.
  - It enters only the Taylor remainder terms, so the effect on the sums is small, as shown below.
- **Float sums.** `win.py:83, 94, 106` form the per-cell products and the window/shell sums in floats.
- **Two further points**, fixed as well:
  - The float cell centre ± radius did not tile [−1, 1] exactly; rounding could leave gaps of about 1e-16.
  - `float(yr.lower())` rounds to nearest, which is not outward for the window selection.

### Fix (`win_arb.py`; `win.py` is imported read only for the unchanged parts)
- **Cells.** The cell breakpoints are exact floats, with b₀ = −1 and b_n = 1. Each cell is the arb ball arb(b_m) ∪ arb(b_{m+1}), and its centre is the exact midpoint of that ball.
- **Half-width.** |x − x_c| ≤ upper(|hh| · rad(s-ball) / 2), with hh the arb element length on the t_c ball.
- **Per cell.** The bound is upper(28 z²/y²) · (upper y − lower y) in arb. The y-hulls are rounded outward.
- **Sums.** The sums are in arb. The sliding windows have the same cover as `win.sliding`: step 0.5, y ∈ [50, 253], window [y_a, 1.5(y_a + 0.5)].
- **Coverage asserted** for each range: windows [50, 380.25], shells [14, 60], (H2f) [60, 700].
  - The cells are consecutive.
  - cos φ ≥ 0.99203 > 0 on every cell.
  - The first cell lies below the range and the last cell above it.
- **Parameter balls.** All 32 t_c balls were run for two ball sets: v3.2's IFT balls (`out/h1slab/ift_M.json`) and the widened (H2f) balls (`out_v33/h2f/ift_M.json`). The value passed on is the **maximum over both sets**. That value also covers the true-curve frame and α shifts for (H3) and (H2), which win.py did not include (T6 §1 note).

### Base values (slab-wide maximum)

| Quantity | arb [cert] | float (replaced) | relative change |
|---|---|---|---|
| (H3) window E₀ | 3.474466e-8 | 3.474454e-8 | +3.4e-6 |
| (H3) window E₁ | 1.260863e-7 | 1.260860e-7 | +2.5e-6 |
| (H3) window E₂ | 3.753835e-7 | 3.753828e-7 | +2.0e-6 |
| (H3) window E₃ | 9.706184e-7 | 9.706168e-7 | +1.7e-6 |
| (H2) shell E₀ (w′) | 8.010214e-11 | 8.009616e-11 | +7.5e-5 |
| (H2) shell E₁ | 1.335326e-10 | 1.335237e-10 | +6.7e-5 |
| (H2) shell E₂ | 2.935820e-10 | 2.935690e-10 | +4.4e-5 |
| (H2f) G₀ [60, 700] | 3.126231e-7 | 3.126228e-7 | +1.2e-6 |
| (H2f) G₁ | 1.003071e-6 | 1.003070e-6 | +1.0e-6 |
| (H2f) G₂ | 2.430170e-6 | 2.430168e-6 | +9.5e-7 |

- Per ball, the largest relative differences from the pinned float values are: windows 1.8e-6, shells 6.1e-5, (H2f) 2.3e-6.
- The maximising window is at y = 253 on every ball, as before.

### Hand-off items on the arb bases (true curve, Ȳ = 1e-12, slab 180)

| Item | j | with error [cert] | threshold | margin (arb) | margin (float) |
|---|---|---|---|---|---|
| (H3) windows (`h3_arb`) | 0 | 3.5086e-8 | 2e-7 | ×5.700 | ×5.700 |
| | 1 | 1.2615e-7 | 1e-5 | ×79.27 | ×79.27 |
| | 2 | 3.7723e-7 | 1e-5 | ×26.51 | ×26.51 |
| | 3 | 1.2821e-6 | 1e-5 | ×7.800 | ×7.800 |
| (H3) Ȳ limit, binding (item 3, E₃) | | **1.575679e-11** | Ȳ = 1e-12 | ×15.76 | float 1.575680e-11 |
| v3.2 (H2) shells (`h3_arb`) | 0, 1, 2 | 8.01e-11, 4.71e-10, 1.56e-6 | 2e-7, 1e-4, 0.2 | ×2495.8, ×2.1e5, ×1.3e5 | ×2496.0, … |
| tightened (H2), w′ (`h2t_v33_arb`) | 0 | 8.0135e-11 | 1e-8 | **×124.79** | ×124.80 |
| | 1 | 4.7060e-10 | 1e-6 | ×2124.96 | ×2125.04 |
| | 2 | 1.5588e-6 | 1e-3 | ×641.54 | ×641.54 |
| tightened (H2), f − U₁ (`h2u_v33_arb`) | 0 | 3.6575e-10 | 1e-8 | **×27.34** | ×27.34 |
| | 1 | 1.7768e-9 | 1e-6 | ×562.81 | ×562.82 |
| | 2 | 1.6640e-6 | 1e-3 | ×600.96 | ×600.96 |
| (H2f) (`h2f_camp_arb`) | 0 | 3.1362e-7 | 1e-6 | **×3.1885** | ×3.1885 |
| | 1 | 1.0033e-6 | 5e-6 | ×4.9834 | ×4.9834 |
| | 2 | 2.4350e-6 | 1e-5 | ×4.1067 | ×4.1067 |

- The largest relative change in any margin is 7.5e-5 (tightened (H2), j = 0). That is far below the 1% stop rule.
- The (H2f) Ȳ limits are unchanged to 1e-6: 5.737e-8, 1.720e-8 and 1.168e-9.
- **Fidelity.** `h3_arb` recomputes camp's Ȳ limits from the float base and matches the pinned `amp_A1_w180.json` `lim` to 2.2e-16 relative. The formula copy is therefore exact.
- The copies `h2t_v33_arb.py`, `h2u_v33_arb.py` and `h2f_camp_arb.py` differ from their originals only in the input and output paths and a label. The diffs show this.

## 2. Centre map in the Lemma O H1 v33 check (`t9_lemmaO_v33_arb.py`)

- **The hand-off radius** is ρ = |Γ̃ + eN − x_c|/λ in the frame of `t9_cone_kch.json` (λ and x_c, both widened by the true-curve shift 3.3e-11). So ρ ∈ [r_lo, r_lo + 10] implies

  |Γ̃| ∈ [r_lo·λ_lo − |x_c| − |e|, (r_lo + 10)·λ_hi + |x_c| + |e|],

  where |e| is the cell's EB from the same script.
- **The radius of each s-ball cell** is enclosed over the whole slab as |G(s, 0) + t·Γ_t(s)| with t ∈ [−1, 1]. These are the positions that the frames themselves use.
- **Asserted for each row:** the radius cells are consecutive, and each of them is among the typed centre cells (the cells meeting [r_lo + 2.0, r_lo + 12.6]). The centre point is therefore inside the union that defines (Γ_c, θ_c).
- **Result:** this holds on all 8 rows (radius cells / centre cells: 28/28, 24/25, 24/24, 24/24, 24/24, 25/26, 27/27, 28/28). The typed rows are certified as they stand.
  - The output rows are byte-identical to `t9_lemmaO_v33.json`: sup|h₀| ≤ 0.21622, H ≤ 0.21849, graph test passes.
  - Note: the cells' slab-wide x-hulls are about 12 wide, wider than the 10.6 arc, so the guard is membership, not containment.

## 3. Far-family H4 (`h4far_v33.py`)

**What v3.3 needs.**
- The far members (L = 113, X = 2, M = 3/20 from `lemmaO_v32_consts.json` far_v32) have centres at hand-off radius ≥ 490 (`zone3_interface.json` covering; asserted).
- They must cover the whole curve beyond 379.25λ₁, out to the far end and through it on the mirror arm (continuation_v33 §gaps item 3).
- So H4 is needed for every centre from radius 490 to the far end. **No cap applies**, but none is needed: the bound decays.

**Method.**
- h4loc.py's per-cylinder bound is used unchanged (its functions are imported read only). It is applied to rows of centre arclength [a, b] with support [a − 2L − s, b + 2L + s].
  - The cells are the slab-180 jet cells that meet the support by arclength or by radius, as h4loc.py selects them.
  - If the support passes the far end, the mirror cells are added.
  - The sups therefore hold for every centre in the row.
- **First centre.** ρ ≥ 490 implies base arclength x ≥ **471.90** (XSTART). This is the least x_lo among the cells whose radius hull reaches 490·λ_lo − |x_c| − |e|.
- **Rows** have widths max(113, a/5), from 471 up to 369,991 (the largest cell x_hi, which is ≥ ℓ). There are 37 rows, and they are contiguous (asserted).
- **Two conservative changes against h4loc.py:**
  1. **Slack s = 3.** h4loc's excess term X L(1/M_lo − 1) + (1 − λ_lo)x_hi + |x_c| is asserted ≤ s. It reaches **2.157** near the far end, so h4loc's slack of 1 would fail beyond x ≈ 2e5.
  2. **y_e = min(y₀, y_ℓ/2)** in the dissipation term and in the C4 pointwise bounds. This is because of the far-end flattening of the level-2 weight: Y_ℓ ≥ min(y, y_ℓ/2) (LEMMA_F §1.3), with y_ℓ/2 ≥ 178,121.
     - Using the C4 bounds at the smaller weight value is the same convention as `t9_far.py`, which evaluates them at the lowest y.
- **Fidelity.** With h4loc's exact geometry, the code reproduces the pinned `h4loc_v32.json` totals for centres 376 and 602 **bit for bit**.

**Results [cert]:**

| Centres (base arclength) | total ∫(h₀⁽⁵⁾)² | / M² = 0.0225 | margin |
|---|---|---|---|
| [471, 584] (contains radius 490–~600) | 1.443e-3 | 0.0642 | **×15.59** |
| [584, 701] | 6.59e-4 | 0.0293 | ×34.17 |
| [701, 841] | 3.70e-4 | 0.0164 | ×60.8 |
| … | decreasing | | |
| [166353, 199624] … [344950, 369991] (flattened weight, mirror) | ≤ 2.8e-9 | ≤ 1.24e-7 | ≥ ×8.08e6 |

- In every row the error part (D₂ term) dominates. The base and lower-order parts are ≤ 1e-7 of it.
- **This closes M9's far gap:** "(490, 602) by an unwritten monotonicity argument, > 602 uncertified" is replaced by a direct certificate over the whole range.

## 4. (H2f)₀

- **The canonical value is 1e-6** (v3.3, β_A = 2; `h2f_camp.py` THR, lemmaU).
- Certified: G₀ ≤ **3.1362e-7** with error. That passes **1e-6 at ×3.19** and **5e-7 at ×1.59**.
- The documents that say 5e-7 are being fixed separately: `STATEMENTS_r2.tex`, `R2_DECISIONS.md`, `zone3_interface.json`.

## Records

**Closed step 4.**
- `p2/out/FINAL_manifest.json` (sha256 `e01cb4da71ba825318b5395f053595d6e980e84c1ac286b9d72ff1e4376a8c94`) verifies before, at the start and end of every driver pass, and after: 34 stages, 523 paths, 97 upstream, 26 scripts, no mismatch.
- A size/mtime snapshot of all step-4 files shows no change outside `p2/v33/`, `p2/out_v33/` and `p2/out_v33_m9_pass1/`.
- `win.py` and every other v3.2 file is untouched. No `.pyc` files were written.

**v33 manifest.**
- `p2/out_v33/FINAL_manifest_v33.json` → **sha256 `8f6f7f95d880f37542482989e76c12a8952392075ed2a573b2892561047d8eb8`**: 19 stages, 281 paths, 35 scripts, 85 reused v3.2 pins. all_items_ok = true.
- **Backup** of the previous manifest: `FINAL_manifest_v33_pre_m9.json` (8202d2c6…, 12 stages). Driver backup: `final_tail_v33_pre_m9.py` (b58e9b75…).
- **Records 0–8 are kept byte-identical**, including `seconds`. The run used `--from 9`.
- **Stage 9 (`notes_v33`) had to be rerun.** Another session (M14) rewrote `zone3_interface.json` at 08:02 today (now a381b1f2…). That made the old pin stale, and the driver correctly refused `--from 12`.
  - notes_v33.json differs only in that input hash. Its far-family content ("centres ≥ 490") is unchanged.
  - Stages 10–11 (B0) reran with byte-identical outputs.
- **Two passes,** each with a freshly emptied `out_v33/m9/`:
  - all 9 stage outputs and 65 m9 JSON files are byte-identical between the passes;
  - the manifests are equal apart from `seconds` and log hashes.
  - Pass 1 is kept in `p2/out_v33_m9_pass1/`.
- Runs: one process at a time, `nice -n 19`, `python3 -I -B`, peak RSS about 140 MB. The m9 stage took 472 s.

**New files and sha256**

| File | sha256 |
|---|---|
| v33/win_arb.py | aa03d65d2feb3d287b355c7180f41ce133696d5b63b64c9d8bdd52fed249c03b |
| v33/m9_balls.py | 66806a18151141b6dd9a03246afd0cf3452636c69028cfbe2a031f9c90329d15 |
| v33/h3_arb.py | 00dfdc213026e5a96648ef10a87a39bc4ca7107d1020bccfe051801e3bebf6b1 |
| v33/h2t_v33_arb.py | a97b780a6f2e14f6ab4726447aa5442d1afdce82c7ae173b7b8991d9fdbc5395 |
| v33/h2u_v33_arb.py | 02d868f5c84c5342a3511ddbd49decf6cd8eec76111e11666afe4cc8312a31ba |
| v33/h2f_camp_arb.py | 747188eb731068320e6775c0457d1e28a0ab7334f6b5bd9ada9257a994a39052 |
| v33/t9_lemmaO_v33_arb.py | e58ff6c0327f150f3a09faff84b2bc4bf111d39003cbc9eceb7908e051d0e3ca |
| v33/h4far_v33.py | 30fe6a5c04f5b547b6e5a86066eef4459ed6a843f8e49be5da0fbee2b0004bf7 |
| v33/final_tail_v33.py | 6187cd583d28fd8f3db4e11113f6b9eadf78575fa2839e8215ae3f48bbbab780 |
| out_v33/m9/win_arb_summary.json | 3c6cb6d88a4b77ed8b1ddb2a2dd3c70592218d429199b0d386d3416fa6f18abe |
| out_v33/h3_arb.json | ee1a216e38b23eefcbc969bd19f7debd50a088a6f538fe0e5d5aea372226048a |
| out_v33/h2t_v33_arb.json | 59d1adc43668726d22dcb96ea47e588f9b95a4812b2c4f155b0557b1d2f16f48 |
| out_v33/h2u_v33_arb.json | 65143fc87c7d5891565a4d69117b0a25961920fb020738e85613655acbc65749 |
| out_v33/h2f_result_arb.json | cac17df46c05e54c630b5f42d8c164de32298eda1123c414480249cd9c7260a6 |
| out_v33/t9_lemmaO_v33_arb.json | ad09d6261387c78ded85bb7af15b760139022d58ffbd7f62011066ab05e1772a |
| out_v33/h4far_v33.json | c6209a42201cc7c3df98a13e6045af495b8dd8e9e5e4d3cc6f360c23de6577a9 |
| out_v33/notes_v33.json (re-pinned) | 7ec984d16b2676f47ebfbd4df74fcb56605f0e279987ee5252ec62412ac4c0de |

## Open points and follow-ups (stated plainly)

- **Downstream documents** should cite the `_arb` outputs. The numbers change only in the 5th–6th digit: (H2) shell E₀ 8.0102e-11; tightened (H2) ×124.79; (H2f) unchanged to 4 digits.
  - New citations: far H4 ×15.59 for all centres, which replaces "376/602 only".
  - The v33 manifest is now 19 stages, 8f6f7f95….
  - Kit B/v33 must ship the 8 new scripts.
- **The H4 far certificate rests on an assumption**: that h4loc's localized-dissipation and C4 argument applies with the flattened weight Y_ℓ in place of y (LEMMA_F §1.3). This is the same reading `t9_far.py` already uses. It matters only for supports beyond y_ℓ/2 ≈ 1.8e5, where the margin is ×8e6.
- **Not rerun:** the far family's arc-by-arc H1 (`t9_lemmaO_far.py`) and its (H5) data (`t9_far.py`). Both were already uniform over every centre to the far end in v3.2; only H4 was pointwise.
- **`b0_dprof.py`** (fix B0) was not part of M9 and was not audited for float half-widths.


---
**Status note (Nov 14, 2026, F34; appended, body unchanged).** Superseded counts: manifest now 21 stages, 8bc47e83...; see F12_REPORT.md
