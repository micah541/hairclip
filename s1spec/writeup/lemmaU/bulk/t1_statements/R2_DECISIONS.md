# Bulk task 1, revision r2: the two design decisions

Fable, Nov 13 2026. Resolves the stops of task 2 (`t2_A/T2_REPORT.md`) and task 5 (`t5_zone3/T5_REPORT.md`).
Files of this revision (r1 files untouched): `STATEMENTS_r2.tex` (balance script passes; not compiled), `STATEMENTS_r2.md`, `CONTRACTS_r2.md`, `HYPOTHESIS_O_r2.md`, `r2_probe/` (scripts and outputs of the checks below, with `SHA256SUMS`).
Labels: **[cert]** arb, project script. **[est]** float model or fit on certified inputs. **[probe]** arb cells on step-4 data, base curve only, scratch run, not a pinned stage; treat as [est].
Nothing outside `t1_statements/` was written. Step-4 scripts were imported read-only (`-I -B`, outputs to scratch).

## Summary

| | decision |
|---|---|
| 1 | Option (a), in a sharper form: a **new hand-off threshold (H2f)** on the whole far support of (A), [60, 700]: G_j ≤ (5e-7, 5e-6, 1e-5). Certified on step-4 data by task 6. The datum of (A) is H2 + H2f. |
| 1b | **The binding quantity is the flux constant K0, not S0v.** B2′ is restated at the level β_A·Ê with **β_A = 5/4** (r1: 2). |
| 2 | Option A with **r_χ = 253 = 379.5/1.5**: χ ramp on [253, 379.5]λ₁, zone III = [253, 379.5]λ₁, B3U(b, c) on {y ≥ 56, r < 253λ₁}. Box (0.15, 0.52) unchanged. |

Three items need Micah (last section). None blocks the resume of tasks 2 and 5.

## Decision 1: the initial datum of (A)

### What step 4 has (read-only check)
- Pinned outputs: (H2) shells on [14, 60]; (H3) sliding windows [y, 1.5y], y ≤ 253, so up to 379.7 (`t9h/win.py`); (H5) far-cylinder jets on x ≥ 150.
- **No existing output bounds the datum on [379.5, 700] at the needed level.**
  - From the thresholds (task 5's chain): |f| ≤ 0.2276 there, so G₀ ~ 1e-3.
  - From step 4's certified far curvature, |h″| ≤ 1.80e-4: G₀ ≲ 1.1e-4 [est]. That bound is the curvature of the profile itself (U₁′ ≈ 1.8e-4), so it cannot see f − U₁.
- **A direct computation with step 4's own window code does.** `r2_probe/far_probe2.py` imports `t9h/win.py` and sums its cells on [60, 700] for all 32 t_c-balls of slab 180 (3.5 s per ball) [probe]:

| range | G₀ | G₁ | G₂ |
|---|---|---|---|
| [60, 379.5] | 4.84e-8 | 1.79e-7 | 5.44e-7 |
| [379.5, 700] | 2.65e-7 | 8.25e-7 | 1.89e-6 |
| **[60, 700]** | **3.13e-7** | **1.00e-6** | **2.43e-6** |

  - Check of the probe: its last window [253, 379.5] gives 3.455e-8; step 4's certified sliding maximum is 3.474e-8.
  - The window energies grow ×3.3 per window, so the far part carries 85% of the total.
  - On [379.5, 700] the base curve has |f| ≤ 0.1234, |y f_y| ≤ 0.1248, |z_k| ≤ (2.9, 4.9, 6.6, 6.1)e-3.
  - Not included: the error term e (1% on the last (H3) window; a few per cent expected here).

### What task 2's report missed
T2_REPORT names S0v as the binding quantity. Its own probe JSONs show that the flux constant is tighter. K0 and K2 are linear in E^A (`assembly_arb.py`, `f4b`), and the verifier consumes them at the bootstrap level of E^A.

| datum (tightened H2) | Ê₀ | FLUX0/K0 at Ê | at 2Ê | S0v vs 0.0153 |
|---|---|---|---|---|
| Phase 2a (H2 only) | 1.96e-7 | ×9.9 | ×4.9 | ×3.53 |
| + 5·H3, H3 = (2e-7, …), nothing beyond 379.5 | 1.20e-6 | ×1.72 | **×0.86** | ×1.97 |
| same + 2·H3 on [379.5, 700] (T2's "passes" row) | 1.60e-6 | ×1.30 | **×0.65** | ×1.77 |
| + 5·H3, H3 = 1.5 × step 4, nothing beyond | 4.66e-7 | ×4.36 | ×2.18 | ×2.78 |

All rows [cert] (`t2_A/probe/asm_*.json`, κ_J applied, outer term at placeholders).
- So the "5·H3" part of (C7) is not harmless at the v3.2 (H3) thresholds: it alone breaks the flux input at the level 2Ê.
- T2's S0v margins are against 0.02; against the B2′ value 0.0153 they are ×0.765 of that.
- Fit to these rows [est]: Ê₀ = datum₀ + 1.93e-7; K0 = 3.5e-23 + 1.94e-15·Ê₀; K2 = 6.7e-22 + 4.9e-15·Ê₀.
- With β = 2 and margin ×2, the datum would have to be ≤ 3.2e-7 at level 0. Step 4's actual value is 3.1e-7 before the error term. So the factor 2 cannot be kept with the present FLUX0.

### Decision
1. **(H2f), new hand-off threshold.** At t₁, in the hand-off chart, for j = 0, 1, 2:

   G_j := ∫ from ln 60 to ln 700 of (28/y)·(∂_s^j (f − U₁))² ds ≤ H2f_j, **H2f = (5e-7, 5e-6, 1e-5)**.

   Same integrand, weight and U as the (H3) windows. Then E^A_j(τ₁) ≤ H2_j + H2f_j (tightened H2 = (1e-8, 1e-6, 1e-3)).
2. **β_A = 5/4 in B2′.** Only the verifier uses E^A at the hypothesis level, and only through K0, K2. Lemma U, Lemma D and Lemma N use the derived Ê. The open–closed argument needs Ê < β_A·Ê and nothing more. B3U keeps its factor 2.
3. (H3) and (H5) stay as they are. Z3 no longer delivers n_III1.

### Why not the other options
- **(b) Z3 at τ₁ from task 5's machinery:** impossible. Thresholds give G₀ ~ 1e-3; the need is ≤ 5e-7.
- **(c) move (A)'s outer cutoff to [253, 379.5]:** the datum would be covered by (H3), but the level-0 outer source rises ×5 to ×40 and eats the same Ê₀ budget, and the launch lemma's neighbour constant q goes 0.0185 → about 0.068 (loop gain 0.04 → 0.14, late M1 about −5%) [est]. It also reopens tasks 3 and 4. Rejected.
- **5·H3 with H3 tightened:** counts the largest window five times; needs H3₀ ≤ 6e-8, which step 4 meets by only ×1.7, and still needs a bound on [379.5, 700].

### Expected numbers at the thresholds [est]; task 2 certifies them
| quantity | value | margin |
|---|---|---|
| Ê | (7.0e-7, 1.0e-4, 0.0726) | E₂ rule 0.1 |
| K0(Ê) | 1.40e-21 | FLUX0: ×2.9 at Ê, **×2.3 at (5/4)Ê**, ×1.47 at 2Ê |
| K2 at (5/4)Ê | 5.2e-21 | FLUX2: ×3.9 |
| S0v | 6.2e-3 | ×2.5 vs 0.0153 |
| core M1 | — | about ×4 |
| step 4 vs H2f | (3.13e-7, 1.00e-6, 2.43e-6) [probe] + e | **about ×1.5**, ×4.9, ×4.0 |

- The two thin margins trade through H2f₀: at 4e-7 they are ×2.7 (flux) and ×1.2 (step 4); at 6e-7, ×2.05 and ×1.8.
- Late M1 of Lemma U: E₁ rises 8%, not ×1.5 as r1 feared; about ×2.2 [est].

## Decision 2: the χ ramp

**r_χ = 253.** Task 5's option A, with 253 in place of "≈ 252".
- The fixed-L family covers r ≥ 251.44λ₁ [cert, task 5]; 253 leaves 1.56λ₁ of slack.
- 253 = 379.5/1.5: the ramp is exactly the last hand-off window, and the last zone-I/II window [168.67, 253] ends where the ramp begins.
- Not chosen: 250 with centres from 298.5 (cylinders would start at 148.5, below the U-fine strip); lower centres (break H6); v3.2 members (bring back (F)).

### Consequences, checked
| item | r1 | r2 | check |
|---|---|---|---|
| B3U(b, c) region | r < 190λ₁ | r < 253λ₁ | — |
| U_max = 0.09 | U₁ ≤ 0.0861 on the support | same | the box and U_max were already global on r ≤ 379.5λ₁ |
| rate margin ×29.07, Ĝ, g, K | — | **unchanged** [cert] | depend on the box only |
| present centres (window Sobolev) | y_c ≤ 126.7 | y_c ≤ 168.67 | float model, K_S = 1.33: \|f\| ≤ 0.0673 at 168.67, **M0 ×2.23** (was ×2.84); M1 ≥ ×2.66 unchanged [est] |
| late centre | M1 ×2.281 | ×2.279 | same model [est] |
| Lemma N, τ_III | 5.22 | 4.071 | ln(700/253)/c⁻ |
| b-term coefficient | 2.5e-8 | 3.21e-8 | 8b̄·(1.875/ln 1.5)/253; negligible against Δ^R |
| D on the ramp | D(190) = 1 + 2.7e-4 | D(253) = 1 + 8.6e-5 | smaller |
| Δ^R | (9.86e-8, 1.04e-6, 0.0161), T6 | **(1.0e-7, 1.09e-6, 0.0175)** [est] | `ramp_out.py` rerun in scratch on [253, 379.5], box (0.15, 0.52), T6 member data; limit 0.03 |
| U-fine strip | [150, 250) | [150, 253) | windows ρ ∈ [150, 168.67]; T6's window [225, 337.5] crossed the ramp |
| zone-III box | — | primary [cert, prelim] | task 5's covering A, [251.44, 379.5] ⊃ zone III; thinnest ×1.495 |

- Model: `r2_probe/rchi_model.py` (`kc2_model.py` with the ramp moved), output `rchi_model.out`. Δ^R: `r2_probe/t6re.py`, `t6re.out`.
- **Two errors in the go/no-go T6 estimate**, both small: its ramp starts at 250 and its pieces are cut at 330, but member 300 covers [251.44, 348.56] and member 380 covers [331.44, 428.56]. The cut must lie in [331.44, 348.56]. With the cut at 340 the T6 value itself becomes 0.0170.
- `ramp_out.py` used the v3.2 box (0.236, 0.24); zone III's certified b₁ is 0.328. Task 5 runs it with (0.15, 0.52). Effect on Δ^R₂: none at 3 digits.
- The mixed lines (`U_ramp.tex` l.28, l.60, l.71, l.114; `U_model2.py` l.175) are read with r_χ in place of 190. Those files are not edited; task 3's text supersedes them.

## Verifier error (v3.3 erratum list, item E-V1)
- **Error.** `closure.unconditional()` marches with `set_mode(delta, True)`: `rhs_iv.py:198` multiplies the absolute flux term by 0, and `unc.phi_val` (`unc.py:17–27`) has no constant term. The sum at `closure.py:182` therefore omits F0·D_k on the march cells, although the docstring's inequality has +2F0. F0 = 3.28e-22 per unit τ; it enters only the tail slope (`closure.py:197`). Confirmed by reading the code.
- **Size.** F0·τ_march/X_c = 4.8e-7 relative, on every J_m. v3.1's J_m (windows from τ₁) carry it too.
- **Handling.** No verifier file is changed and nothing is rerun. `STATEMENTS_r2.tex`, Lemma 4, now has F0 in c_k and gives
  κ_J = (½·max_k(L_up − Ly₁ − 2P_k)⁺ + F0·τ_march)/X_c = 0.042347 + 4.8e-7 = **0.042348** [cert, `t2_A/out/hjw_rep.json`], limit 0.13.
  v3.3 quotes the verifier's J_m, s_m and states H-J-W with J_m(1 + κ_J). All of (A) runs with these.
- **v3.1 erratum:** one line. The script omits the term; no v3.1 margin changes at the quoted digits.

## Hypothesis (O)
(O1)–(O3) are unchanged. (O1) stays at 0.11: that is what the certified segments use. It may be relaxed to 0.15 once task 2 (A5) certifies the segments at (0.15, 0.52).
What changes is the sentence around it (`HYPOTHESIS_O_r2.md`):
- "Certified given (O)" also needs the datum: the tightened (H2) and (H2f). (H2f) is a hand-off condition that step 4 has not yet certified.
- The r1 facts at the v3.1 thresholds are withdrawn. With the honest datum at the v3.1/v3.2 thresholds, the flux input fails at the level 2Ê even with nothing beyond 379.5 (×0.73) [cert, `asm_0_h3_v32.json`, setting v31].

## Needs Micah
1. **Step 4 gets a new hand-off item, (H2f).** It is the existing window stage with its range extended from 379.7 to 700 (geometry to x = 730), plus the error constants on the new bins. No existing step-4 result changes. This is more than a rerun, so task 6 must not start it without approval. If it is refused, the fallback is option (c) above, with its costs.
2. **β_A = 5/4, or a verifier rerun.** r2 adopts β_A = 5/4, which touches no closed step; the flux margin is then ×2.3 and step 4's margin on H2f₀ ×1.5 (both [est]). The roomier alternative reruns the closed verifier for v3.3 with FLUX0, FLUX2 raised ×4 (new rate table; 10 minutes; the tail floor has room ×10⁴). Then β_A = 2 and H2f₀ = 1e-6 (step 4 ×3) would both be possible. Recommended if the certified margins come out below the estimates.
3. **main_v2 wording.** "(F2) follows from (F1) and (O)" is true only with the tightened (H2) and (H2f) added as hand-off conditions. Text in `HYPOTHESIS_O_r2.md`.
4. For information: the v3.1 erratum gains two items (E-V1 above; Lemma S(A)'s datum covers [14, 60] of a support that reaches 700).
