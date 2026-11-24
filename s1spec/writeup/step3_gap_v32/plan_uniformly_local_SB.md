# Plan: a uniformly local Lemma S(B) ("Lemma U") to close the far-field gap

Status: design only. Nothing below is proved. Numbers marked (est.) are float estimates from `gap_v32.tex` / `feas_P.py`; the rest are certified values quoted from v3.1/v3.2.

## 1. What Lemma U must say

**Setting.** s = ln y, f the slope over the S1 ray, levels j = 0, 1, 2. For a centre ρ (a physical radius, 0 < ρ ≤ 380 λ1):
- **Launch time** τ_ρ = τ1 if ρ ≥ y_L λ1, otherwise the time at which ρ/λ(τ) = y_L. Proposed y_L = 56.
- **Deviation** v_ρ = f − U_ρ, with U_ρ(r) = 1.3056 α(τ_ρ) (r/λ(τ_ρ))^p the profile frozen at launch (Prop. B′: forcing is only the b-term and the y⁻⁴ residual). z_j = (y∂_y)^j v_ρ.
- **Weight** Θ_ρ = W(ln(r/ρ)) · (1 + g y⁻⁴), where
  - W is fixed in physical space, equals 1 on the window [ρ, 1.5ρ], decays like (r/ρ)^μ inward and (r/ρ)^(−ν) outward, and satisfies |W^(m)| ≤ Γ_m W for m ≤ 4. Proposed μ = 4, ν = 2.
  - the factor (1 + g y⁻⁴) supplies damping g y⁻⁴ that cancels the rate ε(y)Ĝ; g ≈ Ĝ.
- **Cutoffs.** Inner: a ramp fixed in rescaled coordinates inside the core of Lemma S(A), y ∈ [28, 56]. Outer: the physical ramp [190, 379.5] λ1 of Lemma R. No other cutoff.
- **Energy** 𝔈_j^ρ(τ) = ∫ η_in² χ_out² Θ_ρ z_j² ds.

**Statement.** Under the bootstrap box (|f| ≤ 0.11, y|f′| ≤ 0.52 on 28 ≤ y, r ≤ 380 λ1) and the inputs of §3, for all ρ and all τ ≥ τ_ρ:

    𝔈_j^ρ(τ) ≤ B_j(ρ),   j = 0, 1, 2,

with two regimes:
- **Crude (ρ < y_L λ1, launched late):** B_j is a fixed multiple of the Lemma S(A) level, (60/28)·E^A ≈ (6.7e-7, 3.9e-4, 0.53) (est.). Consequence: on every window, sup|v| and sup y|v′| inside the box with a stated factor. This replaces B3/W3 (E_B, |f| ≤ ½, curvature ≤ S1ᶠ/r).
- **Fine (ρ ∈ [150, 250] λ1, present at τ1):** B_j = hand-off data plus explicit small terms; target unweighted level-2 energy ≈ 0.03 (est.) against a tolerance ≈ 0.3. Consequence: the a-priori bounds K2, K3 on [150, 250) for the fixed-L family (§4).

It also returns the level-1 space–time dissipation ∬ εΘ_ρ z_3², which feeds the outer cutoff of Lemma S(A) (§3).

## 2. Tables and identities to derive or regenerate

| item | what exists | what is new |
|---|---|---|
| Divergence-form identities, levels 0–2 | `lemmaS_F1.div_form`, `prop1_forms.pkl` (sympy-verified) | none; the weight enters linearly |
| Rate tables | `lemmaS_prop1_bounds.table`: monomials (m, i, r) with factor γ^m for ω^(m)/ω | replace γ^m by Γ_m with absolute values of coefficients; re-verify that m counts weight derivatives only (read `lemmaS_prop1_sym.py`). Loses sign cancellations |
| Interpolation step (z_{j+1} against dissipation) | uses Θ″ − κ²Θ with κ = γ + 4 | redo with \|(εΘ)″\| ≤ Γ̃_2 εΘ |
| Weight sups Γ_m | none | arb sups of W^(m)/W and of the (1 + g y⁻⁴) factor |
| Inner ramp [28, 56], increasing cutoff | `Ramp` class handles η²·power | extend `rem`, `rho` to η²·Θ_ρ; inflow term c∫(η²)′Θ z² |
| L⁴ / cubic classes | `lemmaS_L4terms.py`, used as "extras" in `lemmaS_Bbox.py` for power weights | general-weight version of the L⁴ identity (weight derivative terms) |
| Outer ramp | `ramp_out.py`, `ramp_cert.py`, Lemma R (certified) | rerun with Θ_ρ ≤ 1 there (only smaller) |
| Forcing | Prop. B′, `cres.py` | add δU = U(α(τ)) − U_ρ on the inner ramp (\|δU\| ≤ 2·1.3056 α1 y^p ≈ 0.02 at y = 56) |
| Assembly | none | new script: crude fixed point over ρ, then the fine family; arb |

## 3. Bootstrap with Lemma S(A) and the verifier

- **Verifier:** unchanged. It still needs only S0ᵛ, S1ᵛ on [34.25, 60] and the flux constants from Lemma S(A).
- **Lemma S(A) → Lemma U:** on the inner ramp [28, 56], z_0, z_1, z_2 from the shell energies E^A; z_3 from the level-1 dissipation of (A), about 1.3e-4 per unit τ in the εω norm (est.; the affine-in-T budget of `lemmaS_F4B.py` is the certified source).
- **Lemma U → Lemma S(A):** (A) needs (a) the crude box on its support and (b) a quantified outer cutoff. Proposal: cut (A) on a ramp fixed in rescaled coordinates, y ∈ [350, 700] (its forcing integral already stops at 350), with z_{0,1,2} from the crude sups and z_3 from the level-1 dissipation of Lemma U.
- **Order (no circularity beyond the usual continuity argument):** B1 verifier → (A) on [14, 60] given the box and the outer-cutoff inputs → U-crude (box strictly improved) → Lemma O conclusions D5a and Lemma R → U-fine → H6.
- **Lemma O region:** unchanged role; supplies the outer ramp.

## 4. Fixed-L cone family and [150, 250)

- **Family:** L = 50 λ1, X = 3, centres r_c ∈ [300, c_max], c_max ≈ 380. Cylinders [r_c − 150, r_c + 150], conclusion regions [250, c_max + 50]. M = 0.12, K = 0.18, t* margin 2.30 (existing certificate; all members have the same time ratio 1.6e-7).
- **Covering:** r ∈ [150, 250): Lemma U fine. r ∈ [250, c_max + 50]: fixed-L members, unit ratio 1. Beyond: far family (L = 113, conclusion regions from 263), ratio 50/113. Far cylinders' inner part [150, 250): Lemma U fine, ratio 113/50 on a curvature that is small there.
- **Why it helps:** K2 = (50/y)·y|f′| and K3 ∝ (50/y)^(3/2) at y ≥ 150, so the tolerated unweighted level-2 window energy is ≈ 0.3 (est.) instead of ≈ 3e-4 at y = 50.
- **H6:** redo `lemmaO_v32_h6.py` for this geometry (simpler: one unit, no (50/L_m)² scaling).
- **Lemma R / outer ramp:** members of unit 50 on [250, 379.5]: R4 ≤ 0.366 (certified) and (y/L)^5 up to 7.6^5, against 5^5 now; the ramp contribution rises by roughly 10× (est.), still inside the ≈ 0.3 tolerance but to be recomputed. Alternative: keep cone members r/4 only for the ramp.

## 5. Step-4 reruns

- Fixed-L cylinders: H1 (graph, sup|h0|), H2 slope, H3 jets j = 2, 3, 4, H5 data thresholds in L = 50 units, H4 per cylinder (error part ≈ 1e-3 M², est.), Lemma O H1 with C0.
- `t9_h6` closure with the new covering; far-family inner part.
- (H2), (H3) hand-off energies: unchanged thresholds expected; add the weighted hand-off energy ∫Θ_ρ z_j² at τ1 for ρ ∈ [150, 250] as a derived quantity (from the shell and window energies).
- Not affected: the chain, T7, T8, the far family at (0.15, 0.3).

## 6. Risks, cheapest tests, effort

| # | risk | cheapest test that retires it | hours |
|---|---|---|---|
| R1 | The crude family does not close: two-sided tails add ≈ 2.9× to each energy (est.) against a tolerance ≈ 5.7× in √(E1·E2) | float "shell model": steady-state fixed point of B_j(ρ) with tails 1/μ, 1/ν, e^C, inflow; scan μ, ν | 2–3 |
| R2 | Ĝ with Γ_m and absolute values is too large (loss of cancellations; L⁴ extras) | evaluate the existing tables with \|co\|·Γ^m for μ = 4, ν = 2, with and without extras; need ε(28)Ĝ ≲ 0.2 | 2–3 |
| R3 | (A)'s outer cutoff: z_3 from the crude level-1 dissipation gives ≈ 0.07 against E2^A = 0.25 (est.) | compute K for the ramp [350, 700] with `Ramp(increasing=False)` and the dissipation budget | 2 |
| R4 | General-weight L⁴ identity does not go through | derive on paper for the {2,3,3} class; check with sympy | 4–6 |
| R5 | Inner ramp [28, 56]: commutators plus δU exceed the fine budget | `feas_P.py` extended with Γ_m and δU; current estimate 0.017 of 0.3 | 1–2 |
| R6 | Fixed-L outer ramp too costly | rerun `ramp_cert.py` with L = 50 members | 1 |

- **Go/no-go:** R1–R3 and R5–R6, about 8–11 agent-hours. Stop if R1 or R2 fails.
- **Full work, if go:** tables and weights 8–12; inner ramp and forcing 6–8; L⁴ classes 8–12; crude and fine assemblies in arb 8–10; (A) outer cutoff 4–6; fixed-L H6 and Lemma R rerun 4; write-up and integration into v3.2 6–8; step-4 reruns 3–5 (step-4 agent). Total about 50–65 agent-hours beyond the go/no-go, with R4 the largest uncertainty.
