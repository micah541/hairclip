# Decision: PRIMARY. Zone-III box (b₀, b₁, b_z) = (0.1003, 0.4197, 0.0143) [cert], margins ×1.50, ×1.24, ×11.2 against (0.15, 0.52, 0.16)

Task 5, revision r2 (r_χ = 253), Opus, Nov 13 2026. This file replaces `T5_REPORT.md`, which was the r1 stop report.
- Only `t5_zone3/` was written. Nothing was run in step4/.
- One process at a time, under 80 MB. Two consecutive runs are byte-identical; checksums are in `out/SHA256SUMS`.

## Files
| file | content |
|---|---|
| `zone3.py` | the assembly: Z1, Z2, ramp norms, Δ^R, n_III, the H6 redo, U-fine, Z6. Writes `out/zone3_interface.json`; `--ucert PATH` reads B̂ from U_cert.json |
| `ramp_U.py` | the outer-ramp coefficients for the **Lemma-U weight**, in arb with branch and bound. Writes `out/ramp_U.json` (see §3 for why `ramp_out.py` was not just rerun) |
| `z3lib.py` | the zone-III box library, copied from the r1 `z1_box.py` |
| `z1_box.py`, `out/z1_box_prelim.json` | the r1 preliminary run, kept |

## 1. Results
| key | value | label |
|---|---|---|
| `box_III_380`, [253, 379.5]λ₁ | b₀ = 0.10032 (at r = 253), b₁ = 0.41973 (at r = 379.25), b_z = 0.01429 | cert |
| `decision` | primary | cert |
| `box_III_700`, [253, 700]λ₁ | \|f\| ≤ 0.2276, \|r f_r\| ≤ 0.5423 (beyond 379.5: far members, curvature chained); `HO1_ok` = true (×2.20, ×1.84) | cert |
| `tau_III` | 4.0708 | cert |
| `n_III` | (0.0714, 0.327, 3.04) | cert, except c_δU ≤ 2e-8 [hyp], whose effect is below 1e-4 relative |
| `DeltaR` | (4.50e-8, 6.80e-6, 0.01148) | cert |
| `DeltaR_b` | (2.7e-12, 3.3e-9, 4.3e-8) | cert |
| cone family t* margin | ×2.30 (L = 50 for every member) | cert |
| R4 at L = 50 | 0.3661 | cert |
| H6, at the fine-energy estimate of B̂ | min margin ×1.86 (far K3) | cert given that B̂ |
| H6, at the provisional B̂ = (5e-6, 1e-3, 0.53) | cone K3 ×0.76: **fails** | cert given that B̂ |

### Why b₁ moved from 0.328 (r1) to 0.420
In r2 the far family is restricted to centres ≥ 490 (§2, item 2). Its conclusion regions then start at 379.25. So on [253, 379.5] b₁ comes from the unit-50 members alone: r/L up to 7.59, times C₂√t = 0.0395, plus the hand-off curvature.

### n_III against the provisional values
n_III is (0.0714, 0.327, 3.04) against the provisional (0.05, 0.12, 1.3): ×1.4, ×2.7, ×2.3. The resulting N̄² is about (0.07, 0.33, 5.9). That is above `N2_ref` = (0.1, 0.2, 5.0) at levels 1 and 2, and well inside 4 × `N2_ref`, task 2's stop test.

### Δ^R in detail
- **Ramp pieces** (asserted in arb, with tilt e = 0.0287):
  - [253, 331.5] lies in member 300's conclusion region [251.44, 348.56];
  - [331.5, 379.5] lies in member 380's region [331.44, 428.56].
  The cut was chosen from 331.5–348.5 to minimise Δ^R₂; any cut in that range is valid.
- **N_k²** = (8.3e-5, 0.101, 1.357, 2.53e4). For k ≤ 2 these are sup-in-time bounds; for k = 3 a time average. The bounds use Lemma R at L = 50 and the Z1 box.
- **Level 2 parts:**
  - z₃ pair: 0.011297;
  - L⁴ foot: 1.6e-4;
  - genuine cubic (2,2,2) and (2,2,3): 9.2e-6;
  - S33 foot: 3.7e-7;
  - other pairs: 1.3e-5.
- **Comparison only:** the contract-mechanical `ramp_out` number (weight 28/y, ×379.5/28) is 0.0167.

## 2. Things the parent must know

1. **Z4 could not be done by rerunning `lemmaS_v32/ramp_out.py`.** That script bounds the (B) strip energy and is not a bound for Lemma U's ramp, for three reasons:
   - its weight is ω = 28/y (weight-derivative factor 1), while Lemma U's is W·D·φ, with |W′/W| up to 4;
   - it uses `run_F1.box` with Umax = 0.05, while U₁ reaches 0.0861 on the ramp;
   - it uses jet-free A-bounds, so the cutoff × jet (cubic) classes are missing.

   `ramp_U.py` redoes it with P2's exact-table machinery mirrored to the decreasing cutoff:
   - box (0.15, 0.52), U = 0.09;
   - centres ρ ≤ 253;
   - g ∈ [3e5, 6e5] and a ∈ [0, 0.3] as branch-and-bound dimensions;
   - the global Y_{j+1} of the widened design. Ĝ₂₂ = 338258.77 is reproduced and asserted.

   The cubic classes are folded with the Z1 sups |z₀| ≤ b_z and |z₁| ≤ b₁ + pU₁. The two classes that remain genuinely cubic, (2,2,2) and (2,2,3), and the S33 foot are bounded by Sobolev on the ramp (s-length ln 1.5) plus Jensen or Cauchy–Schwarz in time. These are [paper] steps; task 3 should write them into Prop. U.
2. **Far family restricted to centres ≥ 490** (a Z5 covering choice; it is a subset of v3.2's centres ≥ 376).
   - With centres from 376, the far cylinders' bases reach down to 150λ₁ and would need U-fine in units of 113. Their K3 then fails: ×1.06 even at the fine estimate.
   - With centres ≥ 490, the bases start at 264, inside the fixed-L coverage (asserted), and the conclusion regions still cover [379.25, ∞).
   - Cost: b₁ margin ×1.585 → ×1.24 (still above the report level of 1.1).
   - Step 4: no new far check is needed (it is a subset).
   - Task 3: B4's far family must be stated with centres ≥ 490.
3. **H6 depends on the U-fine B̂, which task 4 must certify separately from the crude sup.** The fine B̂ is the level for present centres ρ ∈ [150, 168.67]λ₁.
   - At the provisional crude B̂₂ = 0.53, the cone K3 fails (×0.76, from the U-fine strip part).
   - **Requirement:** U-fine B̂₂ ≤ 0.106, with B̂₀, B̂₁ at (5e-6, 1e-3), for every H6 margin ≥ 1.5.
   - At the go/no-go T5 fine estimate (B̂ = w_min × (1.6e-5, 8.22e-4, 0.0446), so B̂₂ ≈ 0.033) every margin is ≥ 1.86 [est input].
   - Rerun with `zone3.py --ucert <U_cert.json>` once task 4 delivers.
4. **New global Y adds for task 4 (small).**
   - The outer L⁴/S33 foot uses Young with δ_o = 1e-3. It adds to the global Y (`ramp_U.json: foot_out`):
     - level 1: Y₂ +0.0044;
     - level 2: Y₂ +0.108 and Y₃ +0.374 (relative 3e-4).
   - It also uses a_min/10 of the level-2 dissipation **on the outer ramp only**. This is admissible because the inner-ramp foot uses its share on [28, 56], a disjoint support. Task 3 should state this.
5. Z1 uses (H3) and α₁ only through the data at t₁ (|f| ≤ 0.0928, y|f′| ≤ 0.113 on the last window), and (H5) s₀ = 0.08 and m₂ = 0.03 in **L = 50 units** on the fixed-L cylinders. (H2f) is not used anywhere in task 5.

## 3. Acceptance
| quantity | threshold | result |
|---|---|---|
| `decision` | primary | **primary** |
| margin of the chosen box | report if < ×1.1 | ×1.24 (b₁); ×1.50 (b₀) |
| `HO1_ok` | true | true (×2.20, ×1.84) |
| Δ^R₂ | ≤ 0.03 | 0.01148 (×2.61) |
| fixed-L t* | ≥ 1.3 | ×2.30 |
| H6 margins | ≥ 1.3 | at the fine estimate: cone slope ×2.21, K2 ×2.55, K3 ×2.27; far slope ×2.29, K2 ×2.29, K3 ×1.86. At the provisional crude B̂: cone K3 ×0.76, so **open until task 4 certifies U-fine B̂₂ ≤ 0.106** |
| cone slope (step-4 data vs s₀ = 0.08 on the fixed-L cylinders) | ≥ 1.3 | task 6 (new geometry) |
| U-fine K₂, K₃ against 0.18 | ≥ 1.5 | fine estimate: ×3.1 (U-fine K₂ 0.057), cylinder K3 ×2.27. Provisional B̂: K₂ ×1.72, K₃ ×0.76 (see item 3) |
| covering assertions (Z4, Z5) | true | true |
| dependencies | B4, thresholds, D1 only | yes. H6 also uses the U-fine B̂ (task 4) by design; never B2′ or B3U |

## 4. Z6: hand-off thresholds and step-4 reruns for task 6
- **(H2f):** y ∈ [60, 700] at t₁, threshold (5e-7, 5e-6, 1e-5), `new: true`. Step-4 margin (1.5, 4.9, 4.0) [probe]. It needs Micah's approval.
- **Fixed-L cone cylinders** (L = 50, X = 3, centres [300, 380], cylinder [r_c − 150, r_c + 150]); the geometry is new, the values are as in v3.2:
  - (H5): slope ≤ 0.08, |h₀″| ≤ 0.03, ‖h₀‴‖_{L²} ≤ 0.004 (L = 50 units, whole cylinder);
  - H3: jets ≤ 0.12;
  - H4: ∫(h₀⁽⁵⁾)² ≤ 0.0144;
  - H1: graph, sup|h₀| ≤ H − 0.00227;
  - H7: ratio 1.6001e-7.
- **Unchanged:** (H3) windows, α₁, the far (H5) data.
- **Step-4 stages to rerun:**
  - `t9_cone_kch` and `h4loc` (cone cylinders), and the Lemma O H1 graph check, all for the fixed-L geometry;
  - `t9_h6_v32` is superseded by this task's H6;
  - the new (H2f) window stage.
