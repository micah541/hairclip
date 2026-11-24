# Go/no-go tests for the uniformly local Lemma S(B)

**Overall verdict: go-with-risks.** Every test passes on float estimates, but the deciding one (the crude family) passes only if the hand-off shell threshold (H2) at level 2 is tightened from 0.2 to about 1e-3. With the v3.1 thresholds it fails. Step 4's certified shell energies are about 3e-10, so the tightening costs nothing there, but Lemma S(A) must be rerun with it.

All numbers are float estimates. Nothing here is a proof. Scripts and outputs are in this directory; `t2`, `t3` run from `writeup/lemmaS_v32`.

## Verdicts

| test | question | result | margin | verdict |
|---|---|---|---|---|
| T1 | Does the crude family close inside the box y\|f′\| ≤ 0.52? | v3.1 (H2): y\|f′\| ≤ 0.51–0.61. Tightened (H2): 0.26–0.33 | ×0.85–1.02 (v3.1); ×1.6–2.0 (tightened), ×1.3–1.7 with ramp constants ×4 | **pass only with tightened (H2)** |
| T2 | Do the rate constants survive general weights? | Ĝ₂₂ = 1.79e5 at Γ = 4 (crude box, with L⁴ extras) against 28⁴ = 6.15e5 | ×3.4 (inner end at y = 28); ×14 if the support starts at y = 40 | pass |
| T3 | Outer cutoff of Lemma S(A) on [350, 700] | Young route adds 1.4e-4 to E₂ᴬ = 0.247 (2.7e-5 to 0.048). Dissipation route adds 0.09–0.20 | ×1700 (Young); fails by the dissipation route | pass, Young route only |
| T4 | General-weight L⁴ cubic terms (rough) | L⁴ extras raise Ĝ₂₂ by ×2.1–2.7; dissipation budget holds for Γ ≤ 10 | inside T2's margin | pass, rough |
| T5 | Fine regime on [150, 250), fixed L = 50 | K2 ≤ 0.060–0.078, K3 ≤ 0.080–0.128 against 0.18 | K2 ×2.3–3.0; K3 ×1.4 (v3.1 (H2)), ×2.2 (tightened) | pass |
| T6 | Outer ramp [250, 379.5] with unit-50 members | adds 0.016 to the level-2 energy (total 0.046–0.19) | included in T5 | pass |

## What each test did

- **T1 (`t1_crude_shell_model.py`).** For any centre, energy ≤ hand-off data + inflow through the rescaled ramp [28, 56] + ramp commutators; no neighbour terms because the weight never vanishes. Inflow is 1.70·∫W dτ·E^A; z₃ comes from the level-1 dissipation of (A), z₂ at level 1 from the level-0 dissipation. Sobolev on the window, then compare with the box.
  - The level-2 energy is dominated by inflow, i.e. by E₂ᴬ. With v3.1's E₂ᴬ = 0.247 (of which 0.2 is the hand-off threshold) the bound is 2.0–4.4; with E₂ᴬ = 0.048 it is 0.30–0.81.
  - Best weight tried: μ = ν = 4, plateau margin 0–0.3.
- **T2 (`t2_rates.py`).** The tables already store absolute values of the weight-derivative coefficients, so "loss of sign cancellations" is already paid in v3.1; a general weight only replaces γ^m by Γ_m. For W = 1/(e^(−μσ) + 1 + e^(νσ)), Γ_m = max(μ, ν)^m exactly. The damping factor (1 + g y⁻⁴) needs Ĝ < y_a⁴ at the inner end of the support.
  - Sanity check: at Γ = 1 the script gives Ĝ₂₂ = 7.8e4 against v3.1's certified 6.58e4.
- **T3 (`t3_A_outer.py`).** Coefficients of a decreasing ramp on rescaled [350, 700] with the core box. Because ε(350) is tiny, absorbing z₃ in (A)'s own dissipation is cheap. Feeding z₃ from the uniformly local dissipation is not.
- **T4.** Run inside T2 with `lemmaS_L4terms.l4_extra_YC` at γ = Γ. Rough: the L⁴ identity for a general weight has not been derived; its weight-derivative constant is taken as Γ + 4.
- **T5/T6 (`t5_t6_fine.py`).** Centres 150 and 225 λ1; data from (H2) with the tail weight and from the (H3) windows; inflow and inner ramp as in T1 with the tail weight; the profile mismatch on the inner ramp (U(α(τ)) − U₁, bounded by 0.0245) included; outer ramp recomputed with unit-50 members and Lemma R (R4 = 0.366).

## What would change the verdict

- **To no-go:**
  - (H2) level 2 cannot be tightened, or the floor of E₂ᴬ (0.047, from the certified inflow and ramp parts of `lemmaS_F10p_M0.52.json`) rises above about 0.1 when (A) is rerun. Then T1 fails.
  - The true general-weight ramp constants exceed the power-weight ones used here by more than about ×8 (T1 at ×4 still has ×1.3–1.7).
  - The general-weight L⁴ identity fails or costs more than ×3 in Ĝ₂₂ (T2 margin ×3.4; ×14 with the support starting at y = 40, which narrows the inner ramp).
- **To a clean go:** a rerun of Lemma S(A) with (H2) = (1e-8, 1e-6, 1e-3) confirming E₂ᴬ ≈ 0.05, and a paper derivation of the L⁴ identity with a general weight. About 3–5 agent-hours together.
- **Not modelled:** the b-terms, the exact Sobolev constants on a non-flat plateau, the higher-derivative structure of the damping factor, and the coupling between levels. Each is expected to be a few percent, but T1's margin is only ×1.3–2.0.

## Changes to the plan

- Tighten (H2) to about (1e-8, 1e-6, 1e-3) and rerun Lemma S(A): new, required.
- Lemma S(A) outer cutoff: use the Young (order ≤ j) form, not the dissipation of the new lemma.
- Weight: μ = ν = 4. Inner ramp [28, 56] is adequate; [40, 56] is a fallback for the rate criterion.
- Effort estimate unchanged: about 50–65 agent-hours after these tests.

## Corrections from the referee, applied here

- `extra_outer` in `lemmaS_F1.assemble` is the inner-ramp foot/S33 term, not an outer-end term. The outer end of (A) is still unquantified; T3 addresses it.
- Whether a ramp needs z_k for k ≤ j or k ≤ j+1 depends on the method. T3 uses the k ≤ j (Young) form; T1 and T5 use k ≤ j+1 with dissipation.
- The crude far-field bounds are in W4, not W3.
