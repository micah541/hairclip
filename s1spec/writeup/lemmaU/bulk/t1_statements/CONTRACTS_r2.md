# Task contracts, revision r2 (changes to `CONTRACTS.md`)

Bulk task 1 (Fable), Nov 13 2026. **`CONTRACTS.md` stands except where this file says otherwise.** Specification: `STATEMENTS_r2.tex` (supersedes `STATEMENTS.tex`). Reasons and checks: `R2_DECISIONS.md`. Labels as there.

## New fixed values (all tasks)

| symbol | value | replaces |
|---|---|---|
| r_χ | 253 (= 379.5/1.5) | 190 everywhere (χ ramp, zone III, B3U region, Lemma N, D(·), b-term) |
| τ_III | ln(700/253)/c⁻ = 4.071 | 5.22 |
| present centres with a window in zones I/II | y_c ∈ [56, 168.67] | [56, 126.7] |
| U-fine strip | [150, 253)λ₁, windows ρ ∈ [150, 168.67]λ₁ | [150, 250) |
| (H2f) | G_j = ∫_{ln 60}^{ln 700} (28/y)(∂_s^j(f − U₁))² ds ≤ (5e-7, 5e-6, 1e-5) at t₁ | "5·H3 + (28/379.5)·n_III1" |
| datum of (A) | E^A_j(τ₁) ≤ H2_j + H2f_j, H2 = (1e-8, 1e-6, 1e-3) | Remark 15 of r1 |
| β_A | 5/4: B2′ is E^A_j ≤ β_A·Ê_j | 2 |
| κ_J | 0.042348 [cert]; every J_m enters as J_m(1 + κ_J) | "≤ 0.13 [est]" |
| provisional Δ^R | (1.0e-7, 1.09e-6, 0.0175) [est] | (9.9e-8, 1.0e-6, 0.0161) |
| provisional n_III | (0.05, 0.12, 1.3) [est], unchanged; **n_III1 is deleted** | — |

Rule added for all tasks: a "margin" is always threshold/bound, and the threshold is named in the JSON next to it.

---

## Task 2: resume

**Keep:** A1 (done), `hjw_rep.py`, `ic_rep.py`, the `--ic/--h3key/--hjw/--out` options. The probe runs stay as records.

### Changed items

| # | change |
|---|---|
| A1 | Closed. Export in `HJW`: `kappa_J`, `kappa_J_W1`, `kappa_J_W3`, `F0`, `tau_march`, `X_const`, verifier sha256. Statement: Lemma 4 of `STATEMENTS_r2.tex` (F0 is in c_k). |
| A4 | **Replaced.** `ic_rep.py` takes `--h2f a,b,c` and sets E^A_j(τ₁) = H2_j + H2f_j. No H3 term, no zone-III term, no `--zone3` input. Settings, all at the tightened (H2): `design` (5e-7, 5e-6, 1e-5); `lo` (4e-7, 5e-6, 1e-5); `hi` (6e-7, 5e-6, 1e-5); `probe` (3.5e-7, 1.1e-6, 2.7e-6). One more run: `design` at the v3.1 (H2), for the (O) note. The adopted setting is `design`. |
| A4b (new) | **Flux at the bootstrap level.** For each setting report K0, K2 at Ê and at β_A·Ê for β_A = 5/4 and 2, against FLUX0 = 4.1e-21 and FLUX2 = 2.0e-20; and `FLUX0_needed` = 4·K0(Ê), `FLUX2_needed` = 4·K2(Ê) (the inputs that would give ×2 at β_A = 2). Report the parts of Ê₀ (hand-off, inflow, ramp, TOP, outer). |
| A7 (new) | **Audit of two claims.** (i) State the y-range on which `R['S0']`, `R['S1']` bound sup\|v_A\|, sup\|y v_A′\| (expected [28, 60], from `conv = (60/28)^γ`). (ii) `S0_B (y_E=50) ≤ S0V` fails in the probe rows and is close to 1 at the design datum. State whether any statement of `STATEMENTS_r2.tex` consumes it (expected: none; it is the (B)-exit of v3.1, replaced by Lemma U). Keep it in the JSON under `v31_only`. If something does consume it, stop and report. |
| A2 | As r1, with Lemma N at r_χ. In addition report the level-0 outer source ô₀/γ′₀ at `N2_ref` and at 4 × `N2_ref`: Ê₀ is now budgeted by the flux. |
| A3, A5 | As r1. A5's segment redo at (0.15, 0.52) also decides whether (O1) can be relaxed from 0.11 to 0.15; say so in the report. |
| A6 | `IC` becomes {`H2`, `H2f`, `E_tau1`, setting name}. `core` gains `beta_A`, `K0_at_E`, `K0_at_betaE`, `K2_at_E`, `K2_at_betaE`, `FLUX0_needed`, `FLUX2_needed`, and `S0v_margin_vs_0.0153`. |

### Acceptance (replaces the rows of the same name)

| quantity | threshold | expected [est] |
|---|---|---|
| FLUX0/(β_A·K0(Ê)), FLUX2/(β_A·K2(Ê)), β_A = 5/4, setting `design` | ≥ 2.0; report if < 2.2; **stop trigger** if < 1.5 | ×2.3, ×3.9 |
| S0v against 0.0153 | ≥ 2.0 | ×2.5 |
| core M1 | ≥ 2.0 | about ×4 |
| E₂ᴬ at `N2_ref` | ≤ 0.1; report if > 0.08 | 0.0726 |
| ô₀/γ′₀ at 4 × `N2_ref` | ≤ 1e-7 (report if larger) | 2.5e-8 at the placeholders |
| H-J-W | closed: κ_J = 0.042348 ≤ 0.13 | — |

All other rows of r1 stand. Report file name: `A2_REPORT.md`, as before.

---

## Task 5: resume

**Decision taken: option A with r_χ = 253.** The preliminary Z1/Z2 run already covers it (covering A, [251.44, 379.5]).

### Changed items

| # | change |
|---|---|
| Z1 | Region [253, 379.5]λ₁. Rerun `z1_box.py` on that region for the final JSON (the result can only improve on the preliminary one). |
| Z2 | Region [253, 700]λ₁. |
| Z3 | n_III_k ≥ sup_τ ∫ over {350 ≤ y ≤ 700, r ≥ 253λ₁} of (z^A_k)² ds, k = 0, 1, 2; nonempty only for τ − τ₁ ≤ τ_III = 4.071. **n_III1 is deleted** ((H2f) is the datum at τ₁). (H2f) may be used as a hand-off input here; list it in `inputs_used` if it is. |
| Z4 | Ramp [253, 379.5]λ₁, quintic in ln r over ln 1.5. Copy `lemmaS_v32/ramp_out.py`, set its box to (0.15, 0.52) (it has the v3.2 box (0.236, 0.24)), then `ramp_cert` with unit-50 members and Lemma R at L = 50. **Every piece of the ramp must lie inside the conclusion region of its member**, asserted in arb (member 300: [251.44, 348.56]; member 380: [331.44, 428.56]; T6's cut at 330 is outside). Weight Ψ_ρ ≤ D(253)·φ·W. b-term: coefficient 8b̄·sup\|χ′\|/(r_χ(1 − 4ā)) with sup\|χ′\| = 1.875/ln 1.5, times sup_τ ∫_outer W D φ z_j²; export as `DeltaR_b`. |
| Z5 | Family as in r1 (L = 50λ₁, X = 3, centres ≥ 300λ₁). Covering to assert: cone members cover [253, 379.5]λ₁; far members cover [379.5, 700]λ₁. H6 for a cylinder [r_c − 150, r_c + 150]: the part below 253λ₁ from U-fine, the rest from the family (overlap [251.44, 253)). U-fine constants K₂, K₃ on [150, 253)λ₁, from windows [ρ, 1.5ρ] with ρ ∈ [150, 168.67]λ₁ (two windows, ρ = 150 and 168.67, suffice). Do not use a window that enters the ramp. |
| Z6 | Add (H2f): region y ∈ [60, 700] at t₁, values (5e-7, 5e-6, 1e-5), `new: true`, `step4_margin_if_known`: (1.5, 4.9, 4.0) [probe, base curve only; `t1_statements/r2_probe/far_summary.json`]. Keep the fixed-L per-cylinder items of your report. |

### `zone3_interface.json`
- add `r_chi` (253) and `tau_III`;
- `box_III_380` is over [r_χ, 379.5]λ₁;
- delete `n_III1`;
- `DeltaR` carries the piece list and the member of each piece.

### Acceptance
- Δ^R₂ ≤ 0.03 stands (expected 0.0175 [est]).
- K₂, K₃ of U-fine against 0.18: ≥ 1.5 stands (expected about ×3.0 and ×2.0 [est]; the second window is now at 168.67, not 225).
- New: the covering assertions of Z4 and Z5 are true.
- All other rows stand. Report file: `Z5_REPORT.md` (or keep `T5_REPORT.md` and say so in its first line).

---

## Task 3 (not started): text changes
- χ ramp, zone III, B3U and Lemma N with r_χ = 253; Lemma T(iv) with 1/r_χ and sup|χ′| = 1.875/ln 1.5.
- B2′ at β_A·Ê; the paragraph "The level β_A" of `STATEMENTS_r2.tex` §6.
- Datum of (A): (H2), (H2f). Present centres of Lemma U still use (H2), (H3).
- Lemma 4 with F0; H-J-W with J_m(1 + κ_J).
- Window Sobolev for windows with 1.5ρ ≤ r_χλ₁. U-fine on [150, r_χ)λ₁.

## Task 4 (not started): changes
- U3: cover y_c ∈ [56, 253], as r1.
- U6: present centres on a cover of y_c ∈ [56, 168.67]; outer term Δ^R·D(253)·φ; read `beta_A` and the flux rows from `A_interface.json`; (K1) is "fluxes at β_A·Ê".
- Acceptance: present M0 ≥ 2.0 stands; report if < 2.1 (expected ×2.23 at y_c = 168.67 [est]).
- If the certified late M1 is below 1.8, report the row with H2f at `lo` (task 2), not "H3 tightened".

## Task 6 (after 2–5): additions
- **(H2f) on the step-4 data. Do not start before Micah approves** (new step-4 item; `R2_DECISIONS.md`, last section).
  - Base part: the window stage with geometry to x = 730 and the cell sum over y ∈ [60, 700], all 32 balls (the probe script shows the call; 3.5 s per ball).
  - Error part: the chart constants of `t9.py` on the new bins and a `camp` item with its Ȳ limit; check first that T7's far bound on |e| and the N3′ boxes are stated for x up to 730.
  - Back up changed scripts as `*_pre_h2f.py`; add the stage to `final_tail.py`; existing outputs must stay byte-identical.
- The fixed-L per-cylinder checks listed by task 5.
- The tightened (H2), as planned.

## Stop triggers: one added

| trigger | test | task |
|---|---|---|
| flux input not met at the bootstrap level | `A_interface.json`: FLUX0/K0_at_betaE < 1.5 or FLUX2/K2_at_betaE < 1.5 at setting `design` | 2 |
| step 4 does not meet (H2f) by ×1.2 | task 6's certified G_j against (5e-7, 5e-6, 1e-5) | 6 |

The four triggers of `BULK.md` stand.
