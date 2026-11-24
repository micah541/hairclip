# Go/no-go: uniformly local Lemma S(B) ("Lemma U") — final verdict, Nov 13, 2026

**Verdict: GO-WITH-RISKS. Recommend starting the full proof now, staged, with a kill checkpoint after the first 10–12 h.** No precondition is left. The decision is the PI's.

- Both v1 conditions for a clean go are met: Lemma S(A) rerun at the tightened (H2) gives E₂ᴬ = 0.048 [cert], and the general-weight L⁴ identity is derived and sympy-checked.
- It is not a plain "go" because nothing about Lemma U itself is certified, the thinnest margin is ×1.4–1.9 [est], and this review found items the tests do not model (below). None of them flips a test.

Labels: **[cert]** arb, outward rounding. **[est]** float estimate. **[paper]** hand derivation. **[rev]** this reviewer's scratch estimate, float, no script on disk; redo in the full proof.

## Verdict table

Weight W = 1/(e^(−μ(σ+a)) + 1 + e^(ν(σ−ln1.5−a))), written (μ, ν, a). KS = safety factor on the ramp constants.

| test | question | result | margin | verdict |
|---|---|---|---|---|
| T1 crude family | y\|f′\| ≤ 0.52 and \|f\| ≤ 0.11 on the window? | (4,4,0.3): y\|f′\| ≤ 0.301 / 0.323 / 0.361 at KS = 1 / 2 / 4 [est]; \|f\| ≤ 0.031–0.036 [est] | M1 ×1.72 / 1.61 / 1.44; M0 ×3.0–3.5. (4,4,0): M1 ×1.86–1.55. (4,2,0.5): M1 ×1.53–**1.27** | pass |
| T2 rates | Ĝ₂₂ < 28⁴ with a general weight? | Ĝ₂₂ = 2.398e5, θ_a = 0.390, Gw = 1.639 (ν = 4); 1.965e5 (ν = 2) [est] | ×2.56 (ν = 4), ×3.13 (ν = 2); fine box ×5.6–6.9 | pass |
| T3 outer cutoff of (A) | cost on [350, 700] | Young route adds 2.67e-5 to E₂ᴬ [est, **stale**: not rerun, script now crashes] | ×1800 | pass, Young route only |
| T4 L⁴, general weight | identity and cost | identities I1–I8 hold for any weight [paper + sympy]; L⁴ factor on Ĝ₂₂ ×3.06 (power weight, same box: ×4.08) [est]; ramp foot + S33 = 0.0085 at level 2, 1% of B₂ [est] | inside T2 | pass |
| T5 fine family [150, 250) | K2, K3 ≤ 0.18 | K2 ≤ 0.058–0.059, K3 ≤ 0.076–0.079 [est] | K2 ×3.0–3.1; K3 ×2.3–2.4 | pass |
| T6 outer ramp, unit 50 | level-2 cost | Δ = (9.9e-8, 1.0e-6, 0.0161) [est]; it is 60% of T5's level-2 energy | inside T5 | pass |
| C1 Lemma S(A) rerun, (H2) = (1e-8, 1e-6, 1e-3) | E₂ᴬ < 0.1? | E^A = (1.2110055e-7, 8.1741747e-5, 0.048044818), γ′ = (0.2463, 0.2462, 0.1431) [cert]; floor 0.047045 | ×2.1 below the no-go trigger; core S0v / M0 / M1 ×5.38 / 6.58 / 4.86 [cert] | pass |
| C2 step 4 meets the tightened (H2) | Ȳ limits | 4.97e-8, 9.75e-11, 2.57e-11 against Ȳ = 1e-12 [cert] | ×4.97e4 / ×97.5 / ×25.7 | pass |
| C3 profile mismatch δU | size on the ramp | \|δU\| ≤ 1.06e-6 at y = 56 [paper; inputs cert; needs \|α\| ≤ α₁, a ≤ ā] | negligible | pass |

Flip points of T1 [est], one parameter at a time:
- E₂ᴬ: ×8.0 (KS = 4) to ×13 (KS = 1) for (4,4,0.3); ×4.4–7.5 for (4,2,0.5).
- Ĝ₂₂: ×1.84 (KS = 4) to ×2.07 (KS = 1) for (4,4,0.3); up to ×2.26 for (4,2,0.5).
- KS: 14.8 for (4,4,0.3); 10.2 for (4,2,0.5); 18.5 for (4,4,0).

**The binding tolerance on Ĝ₂₂ is the M1 flip, ×1.84–2.07, not the rate margin ×2.56.** Both are the same parameter; the box fails first, because Gw = 1/(1 − Ĝ₂₂/28⁴) blows up before the rate criterion does.

Joint stress [rev], (4,4,0.3), (H3) data included: KS = 2 with Ĝ₂₂ ×1.25 gives M1 ×1.44; KS = 4 with Ĝ₂₂ ×1.5 gives ×1.17. Still a pass.

## What I checked

- Reran t1, t5 and t7 in a scratch copy: outputs and JSON are byte-identical to those on disk.
- t7's step-4 inputs were rewritten at 17:14 on Nov 12 (final manifest), after t7 ran. The rerun reproduces camp's limits bit for bit, so C2 holds against the final chain.
- CLEAN-GO rows read the right inputs: E^A and γ′ from `assembly_arb.json`, Ĝ₂₂ per weight from `t8_L4_cost.json`, the t7 assertion. t8's fix is the JSON dump only (diffed against `_pre_jsonfix`).
- The C1 log has "claims failing: []" and is identical across the three rechecks.
- Not rerun by me: t8, the (A) rigor scripts, `l4_identity_check.py`, t2, t3.

## What changed since v1

- **T1** no longer depends on an estimated floor: E^A is certified. It now includes δU, the L⁴ ramp foot and the outer ramp. M1 margin ×1.44–1.72 for (4,4,0.3); v1 quoted ×1.6–2.0, and ×1.3–1.7 at KS = 4.
- **T2 corrected: ×3.4 → ×2.56.** v1 used κ = 5; the honest value is κ̃ = ν + 4 = 8. The margin "×14 from y = 40" becomes ×10.7.
- **T4 is no longer "rough".** The identity is derived for an arbitrary weight, so "the identity fails" is excluded. For μ = 4 the damping factor adds nothing to Γ_m, κ̃ or κ̃₂. The cost is still [est].
- **T5** uses Gw = 1.639 from the crude Ĝ and the Prop. B′ δU: K3 ×2.2 → ×2.3–2.4.
- **(H2) tightening** is no longer a condition: C1 and C2 are certified.
- **δU**: v1's 0.0245 was far too pessimistic; 1.06e-6.
- v1's no-go trigger "L⁴ costs more than ×3 in Ĝ₂₂" is superseded by the flip point above.

## Errors and weaknesses found in the evidence

1. **`t3_A_outer.py:40–42`** reads `t1_crude_shell_model.json` as a list; t1 now writes a dict. T3 crashes if rerun, and its output is from the Nov 12 08:24 inputs (Gw 1.41, estimated floor). The effect is negligible (about 3e-5 against 0.048 [rev]), but T3 is unchecked for the clean-go inputs.
2. **`feas_P.py:11,22`** (source of the ramp constants K in `t1:22`, `t5:12`) uses γ = −4 in tables that store |coef| (`lemmaS_prop1_bounds.py:49`). Odd weight-derivative orders then enter with a minus sign. Its Ĝ = (21, 302, 2.7e3), and ε(28)Ĝ₂₂ = 0.0044 in `gap_v32.tex`, are wrong; the honest values are (4.9e4, 2.3e4, 7.8e4) before L⁴. The K_jj entries inherit this through YC (`feas_P.py:33`). The K_{j,j+1} entries, which carry 97% of T1's ramp term, do not. Estimated effect on B₂: ≤ 2% [rev].
3. **`t1:8,36`**: "windows are smaller" is false after the tightening. (H3) is (2e-7, 1e-5, 1e-5), above (H2) = (1e-8, 1e-6) at levels 0 and 1, and T1 omits it. With it: M1 ×1.72 → ×1.67, M0 ×3.51 → ×3.18 [rev].
4. **`t1:4,29` against `t8:4,37` and plan §1**: T1 sweeps the whole weight through the ramp and starts from hand-off data; t8's constants assume the centre is already at y ≥ 56. For a centre launched late, the launch data come from E^A (level 2: 0.18, fifty times T1's data term), and beyond (A)'s outer cutoff from the other centres. Redone that way [rev]: (4,4,0.3) B₂ = 0.52, M1 ×1.98 / 1.80 / 1.56, so T1 is conservative for ν = 4. For **ν = 2** the outer tail adds a neighbour factor of about ×1.5, and M1 falls to ×1.55 / 1.42 / 1.23.
5. **`t1:25`, `t5:55`**: the box is tested only on the window [56, 84] (U = 0.018). The plan's single box |f| ≤ 0.11 out to r = 380λ₁ meets U(380) = 0.086 plus S0 ≈ 0.016: **×1.07** with α₁ = 1.5e-4, ×1.31 with step 4's certified |α| [rev]. v3.1 avoided this with per-window boxes U(1.5y) + 0.15.
6. **Labels.** "Certified E^A" is certified arithmetic for a conditional statement: (A) assumes the box on its whole support and has no outer-end term. That is the gap itself, and Lemma U must supply it. The foot weights are certified for (A)'s power weight on [14, 28]; their use on [28, 56] with a general weight is [est].
7. Minor: `t1:61–62` and `t5:15,33` use typed constants (they match the outputs). T5 models one weight only, (4,4,0.3). The quoted M1 range "×1.44–1.86" is for ν = 4; with (4,2,0.5) it reaches ×1.27.

## Risks the full proof must retire, in priority order

| # | risk | margin now | how it is retired |
|---|---|---|---|
| 1 | Certified Ĝ₂₂ for the general weight (Γ_m, κ̃, κ̃₂ in arb; L⁴ classes with Γ_eff) exceeds the estimate | ×1.84–2.07 to the M1 flip [est]. The same script at Γ = 1 gives 7.8e4 against 6.57e4 certified, so it runs 18% high there | arb enclosure; first task |
| 2 | General-weight ramp constants on [28, 56] exceed the power-weight K (never computed with honest tables; see error 2) | KS flip 14.8; M1 ×1.44 at KS = 4 [est] | honest K_jk in arb, with the L⁴ ramp terms; second task |
| 3 | Launch data and the bootstrap between (A) and Lemma U (error 4): (A)'s outer cutoff [350, 700], the neighbour term, the order of the continuity argument | M1 ×1.56–1.98 [rev]; T3 stale | write the launch lemma first; fix and rerun T3 |
| 4 | Global box at the far end (error 5) | ×1.07 [rev] | radius-dependent box; design choice, no cost in rates |
| 5 | Fine family: outer ramp with unit-50 members (Lemma R ramp recomputed), fixed-L H6, step-4 reruns | K3 ×2.3 [est]; the v3.2 cone slope check is ×1.40 | Lemma R rerun early; then `final_tail.py` stages |
| 6 | Not modelled anywhere: b-terms, coupling between levels, exact Sobolev constants on a non-flat plateau | each a few percent expected, against ×1.4–1.9 | arb assembly |

## Recommended parameters

- **Weight (μ, ν, a) = (4, 4, 0.3)**, with a kept free in [0, 0.3] in the assembly (t8: no constant depends on a; a = 0 gives the best T1 margin but is untested in T5).
  - μ = 4 is forced: it cancels ε on the inner tail and gives κ̃ = ν + 4.
  - ν = 4, not the plan's ν = 2: better M1 margin in every row and no neighbour problem at launch. Its cost is the rate margin, ×2.56 instead of ×3.13.
- **Inner ramp [28, 56], launch y_L = 56.** Fallback [40, 56]: rate margin ×10.7 [est], but a narrower ramp raises the ramp constants.
- **(H2) = (1e-8, 1e-6, 1e-3)** [cert on both sides]. Include (H3) in the hand-off data.
- **Box by radius, not one global box.** Use (A)'s certified core box on y ∈ [28, 60] (M1 ≤ 0.107 [cert] against the 0.52 assumed) for the rate criterion; Ĝ₂₂ on the (0.21, 0.21) box is 1.09e5, so this could lower Gw from 1.64 to about 1.40, the equivalent of a rate margin of about ×3.5 [rev, not run]. Use U(1.5y) + const for r ≥ 150λ₁.
- **(A) outer cutoff:** rescaled ramp [350, 700], Young route.
- **Fine family:** L = 50, X = 3, centres ≥ 300λ₁; damping from the crude Ĝ₂₂.

## Staging

- **Phase 1 (10–12 h):** risks 1–3 in arb or on paper, then T1 again with those inputs.
- **Checkpoint:** continue if the M1 margin is ≥ about ×1.3 with certified Ĝ₂₂ and K. Otherwise stop and report; the fallbacks are the radius-dependent box and the ramp [40, 56].
- **Then** the rest: 40–50 h, plus 3–5 h of step-4 reruns.

## What would still flip the verdict to no-go

- Certified Ĝ₂₂ ≥ about 4.4e5 (×1.84) on the box that governs y = 28, and the radius-dependent box does not rescue it.
- Honest general-weight ramp constants more than about ×10 above the power-weight K.
- The launch argument needs a neighbour factor that does not close for any ν, or (A)'s outer cutoff cannot be fed without circularity.
- E₂ᴬ above about 0.2 once (A) is rerun with its real outer cutoff and the box Lemma U delivers (now 0.048 [cert], conditional).
- The fixed-L family fails Lemma R or H6, or step 4 cannot meet its hand-off data in L = 50 units (not tested here).

## Files

- Certified: `lemmaS_A_rerun/rigor/assembly_arb.json`, `lemmaS_A_rerun/rigor_assembly_arb.log`, `t7_step4_H2_tight.{py,out,json}`.
- Estimates (Nov 13): `t8_L4_cost.*`, `t1_crude_shell_model.*`, `t5_t6_fine.*`. Stale: `t2_rates.*`, `t3_A_outer.*`.
- Derivation: `L4_general_weight.md`, `l4_identity_check.{py,out}`.
- Earlier verdict: `GONOGO_v1.md`. Session state: `HANDOFF_cleango.md`.
