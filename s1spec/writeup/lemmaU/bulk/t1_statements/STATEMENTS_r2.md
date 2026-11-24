# Bulk task 1: statements for Lemma U, revision r2 (summary)

Fable, Nov 13 2026. Supersedes `STATEMENTS.md`. Full text: `STATEMENTS_r2.tex` (not compiled; braces, environments, `$`, `\[ \]`, labels checked by script). Contracts: `CONTRACTS.md` as amended by `CONTRACTS_r2.md`. Decisions and checks: `R2_DECISIONS.md`. Hypothesis (O): `HYPOTHESIS_O_r2.md`.
Labels: **[cert]** arb project script; **[est]** float or fit; **[rev]** kill-check 2 scratch; **[paper]** by hand; **[probe]** arb cells on step-4 data, base curve only, not pinned.

## What r2 changes (§2 of the tex)

| # | change |
|---|---|
| R1 | χ ramp on **[r_χ, 379.5]λ₁, r_χ = 253**. Zone III = [253, 379.5]λ₁. B3U(b, c) on {y ≥ 56, r < 253λ₁}. |
| R2 | Datum of (A): new hand-off threshold **(H2f)**, G_j = ∫_{ln 60}^{ln 700} (28/y)(∂_s^j(f − U₁))² ds ≤ (5e-7, 5e-6, 1e-5). E^A_j(τ₁) ≤ H2_j + H2f_j. |
| R3 | The flux constant K0 is the binding quantity of (A). B2′ is E^A_j ≤ **β_A·Ê_j, β_A = 5/4**. |
| R4 | Lemma 4 carries the flux floor F0 (verifier erratum). **κ_J = 0.042348** [cert]. |
| R5 | U-fine on [150, 253)λ₁, windows ρ ∈ [150, 168.67]λ₁. |

Everything else is as r1. The sections below repeat r1 with the changes worked in.

## Decision on the box (unchanged)

**Global box (0.15, 0.52)**, U_max 0.09, sup|v| ≤ 0.16, design data56: Ĝ₂₂ = 338258.8, g = 3.503e5, rate margin ×29.07 [cert]. Fallback (0.20, 0.60). None of these depends on r_χ. Task 5, preliminary [cert]: zone III is primary, b₀ = 0.1003, b₁ = 0.3281, b_z = 0.0143 (×1.495, ×1.585, ×11.2).

## 1. Lemma D in window form (§4)

- Statement, hypotheses, layer-cake corollary, S_X, Σ_X and the data term 𝓓_k: as r1.
- **H-J-W** now holds with J_m(1 + κ_J), κ_J = (½·max_k(L_up − Ly₁ − 2P_k)⁺ + F0·τ_march)/X_c = 0.042347 + 4.8e-7.
  - (W1) fails on every march cell when checked with P_k ≥ 0. (W2) holds. (W3) holds only with F0 added to c_k: `closure.py:182` omits it.
  - The verifier files are not changed. v3.3 quotes its J_m, s_m and applies Lemma 4.

## 2. R5 in L² form (§5)

As r1. Two additions:
- Lemma N splits O at r = 253λ₁; τ_III = 4.071.
- **Remark (datum of (A))** is rewritten around (H2f). Probe on step 4's base curve: G = (3.13e-7, 1.00e-6, 2.43e-6); the part on [379.5, 700] is (2.65e-7, 8.25e-7, 1.89e-6). Expected step-4 margins ×1.5, ×4.9, ×4.0 [est].

## 3. The bootstrap (§6)

- **Set.** B1; B2′ (core box, E^A ≤ (5/4)·Ê on [14, 700]); B3U: (a) 𝔈^ρ_j ≤ 2B̂_j, (b) box (0.15, 0.52) and (c) |v| ≤ 0.16 on {y ≥ 56, r < 253λ₁}; B4.
- **Why 5/4.** Only the verifier uses E^A at the hypothesis level, through K0 and K2, which are linear in E^A. Everything else uses the derived Ê. Open–closed needs Ê < β_A·Ê only.
- **Order.** As r1, with: fluxes at β_A·Ê in D1; (H2), (H2f) in D3′; r_χ in the box rows; U-fine on [150, 253)λ₁.
- **Closure** (Proposition 16): (K1) is "fluxes at β_A·Ê below FLUX0, FLUX2"; (K3) is for windows [ρ, 1.5ρ] with 1.5ρ ≤ 253λ₁.
- **Loop gain:** as r1.

## 4. Box and zone III (§7)

Zone III delivers: Z1 on [253, 379.5]λ₁; Z2 out to 700λ₁; Z3 n_III on O ∩ {r ≥ 253λ₁} (n_III1 is dropped); Z4 Δ^R for the ramp [253, 379.5] with the box (0.15, 0.52) and every piece inside its member's conclusion region, plus the b-term (coefficient 3.21e-8). Expected Δ^R = (1.0e-7, 1.09e-6, 0.0175) [est].

## 5. Hypothesis (O) (§8, `HYPOTHESIS_O_r2.md`)

(O1)–(O3) unchanged. "Certified given (O)" needs in addition the tightened (H2) and (H2f). r1's numbers at the v3.1 thresholds are withdrawn.

## Expected margins after r2 [est]

| quantity | value | acceptance |
|---|---|---|
| FLUX0/(β_A·K0(Ê)) | ×2.3 | ≥ 2.0 |
| FLUX2/(β_A·K2(Ê)) | ×3.9 | ≥ 2.0 |
| step 4 against H2f₀ | ×1.5 | ≥ 1.2 (trigger) |
| S0v against 0.0153 | ×2.5 | ≥ 2.0 |
| E₂ᴬ | 0.0726 | ≤ 0.1 |
| late M1 (Lemma U) | about ×2.2 | ≥ 1.5 (trigger) |
| present M0, worst at y_c = 168.67 | ×2.23 | ≥ 2.0 |
| zone III b₀ | ×1.495 [cert, prelim] | report if < 1.1 |
| Δ^R₂ | 0.0175 | ≤ 0.03 |

## Findings of r1 and r2

| # | finding | status |
|---|---|---|
| C5 | one frozen profile U₁ for all centres | as r1 |
| C6 | point terms of the C² cutoffs | as r1; task 2 (A3) |
| C7 | (A)'s datum must cover [14, 700] | **settled by (H2f)**; r1's "not a stop trigger" was wrong: with 5·H3 the flux input fails at 2Ê (×0.86 [cert]) |
| C8, C9 | localisation at y = 28; segments at the far-field box | as r1; task 2 (A5) |
| C10 (new) | K0, not S0v, binds (A); the bootstrap factor on E^A is free | β_A = 5/4 |
| C11 (new) | verifier omits F0 on march cells | absorbed in κ_J; erratum item E-V1 |
| C12 (new) | T6's ramp [250, 379.5] and its cut at 330 lie outside the members' conclusion regions | ramp [253, 379.5], cut in [331.44, 348.56] |
