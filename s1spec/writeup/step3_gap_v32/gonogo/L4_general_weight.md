# The L⁴ identity and the cubic terms for a general weight (Lemma U)

Status: derivation on paper (this note) plus a sympy check of every identity used (`l4_identity_check.py`: identities I1–I8 verified for an *arbitrary* weight). The numbers in §4 are **float estimates** (`t8_L4_cost.py`, `t8_L4_cost.out`/`.json`). The one exception is the power-weight foot weights, which are certified (`lemmaS_A_rerun/rigor/weights_exact`).

## 1. Setting and notation

- s = ln y, and z_k = ∂_s^k v (= (y∂_y)^k v), so z_{k+1} = z_k′.
- Energy weight Θ(s) > 0 (C³). Dissipation weight Ξ = εΘ with ε = e^{−4s}. On a cutoff ramp, the weights are Θη² and Ξη², with 0 ≤ η ≤ 1 and |η′| ≤ η′_max.
- N_k = ‖z_k‖_Ξ = (∫Ξ z_k² ds)^{1/2}. In the time-integrated (variation-of-constants) bounds, the same with ∬ ds dτ. X = ‖z₂²‖_Ξ = (∫Ξ z₂⁴)^{1/2}.
- Box: |z₁| ≤ M₁ on supp Ξ. The crude box gives M₁ = 0.52, as in v3.1's core; the bootstrap closes it with S₁ + 1.05U ≤ 0.52, so |z₁| ≤ S₁ ≤ M₁.
- Constants of the weight:
  - L = Ξ′/Ξ, the log-derivative;
  - κ̃ = sup|L|;
  - κ̃₂ = sup(Ξ″/Ξ)₊ = sup(L² + L′)₊;
  - Γ_m = sup|Θ^{(m)}/Θ|.

  For the power weight of Lemma S(A) (Θ = ω = e^{−γs}): L ≡ −κ with κ = γ + 4, Ξ″/Ξ = κ², and Γ_m = γ^m.

## 2. Lemma (general-weight L⁴ bounds)

Assume κ̃ < ∞ and that the boundary terms of the integrations by parts vanish. This holds because Ξ is cut off: η_in at the inner ramp and χ_out at Lemma R's outer ramp.

**(a) L⁴ bound.** X ≤ M₁ (3N₃ + κ̃ N₂).

**(b) Class {2,3,3}.** ∫Ξ|z₂|z₃² ≤ ½ X (κ̃N₃ + N₄) ≤ ½ M₁ (3N₃ + κ̃N₂)(κ̃N₃ + N₄).

**(c) On a cutoff ramp,** with Ξ → Ξη² and the cutoff-free ramp datum F₂ = ∫_ramp Ξ z₂² (space–time ∬ in the time-integrated form):

  X_η ≤ M₁ (3N₃ + κ̃N₂ + 2η′_max F₂^{1/2}),

  ∫Ξη²|z₂|z₃² ≤ ½ X_η (κ̃N₃ + N₄) + ⅓ ∫(Ξηη′)′ |z₂|³,  with |(Ξηη′)′| ≤ C_η^Ξ Ξ and C_η^Ξ = sup_ramp |L ηη′ + η′² + ηη″|.

**(d) Interpolation.** N_{k+1}² ≤ N_k N_{k+2} + ½ κ̃₂ N_k².

**Proofs.** All are integrations of exact pointwise identities, followed by Cauchy–Schwarz.

- **(a).** Integrate (I1): Ξ z₂⁴ = (Ξ z₁ z₂³)′ − 3Ξ z₁ z₂² z₃ − Ξ′ z₁ z₂³. This gives
  X² = −3∫Ξ z₁z₂²z₃ − ∫Ξ′ z₁z₂³.
  - |first term| ≤ 3M₁ (∫Ξz₂⁴)^{1/2}(∫Ξz₃²)^{1/2} = 3M₁ X N₃.
  - |second term| ≤ κ̃ M₁ ∫Ξ|z₂|³ ≤ κ̃ M₁ X N₂.

  Divide by X.
- **(b).** Put φ = ½ z₂|z₂|. It is C¹ with φ′ = |z₂| z₃ and |φ| = ½ z₂².
  - By (I2), ∫Ξ|z₂|z₃² = ∫Ξ z₃ φ′ = −∫(Ξ′ z₃ + Ξ z₄) φ.
  - This is at most ½(κ̃ ∫Ξ z₂²|z₃| + ∫Ξ z₂²|z₄|) ≤ ½ X (κ̃N₃ + N₄).
  - Then apply (a).
- **(c).** By (I3), (Ξη²)′ = Ξ′η² + 2Ξηη′. Two extra terms arise.
  - In (a), the extra term satisfies |2∫Ξηη′ z₁z₂³| ≤ 2M₁ (∫Ξη²z₂⁴)^{1/2} (∫Ξη′²z₂²)^{1/2} ≤ 2M₁ η′_max X_η F₂^{1/2}.
  - In (b), the extra term is −∫2Ξηη′ z₃φ. By (I4) it equals ⅓∫(Ξηη′)′|z₂|³ exactly, using z₃φ = ⅙(|z₂|³)′ (checked on both signs of z₂). Expand (Ξηη′)′ = Ξ(Lηη′ + η′² + ηη″).
- **(d).** Integrate (I5): Ξ z_{k+1}² = (Ξ z_k z_{k+1} − ½Ξ′z_k²)′ − Ξ z_k z_{k+2} + ½ Ξ″ z_k².

None of the steps uses the form of Θ: the identities are exact for every weight (sympy, `l4_identity_check.py`). So the general-weight L⁴ identity cannot "fail"; only its constants change.

## 3. What changes against the power-weight tables

| ingredient | power weight (Lemma S(A)) | general weight (Lemma U) |
|---|---|---|
| divergence-form tables (weight derivatives, m ≤ 3 occur) | \|coef\|·γ^m | \|coef\|·Γ_m ≤ \|coef\|·Γ_eff^m, with Γ_eff = max_m Γ_m^{1/m}. The tables already store \|coef\|, so no further loss of cancellations. |
| L⁴ class `('L4', r, w)`: w·X·N_r | w M₁ (3N₃ + κN₂) N_r | the same with κ → κ̃ (linear) |
| class {2,3,3} (`s33='ibp'`) | ½ w M₁ (3N₃+κN₂)(κN₃+N₄) | κ → κ̃ (quadratic in κ̃) |
| interpolation in Ĝ_jj = Y_j + Y_{j+1}(Y_{j+1}/(4B_c) + κ²/2) | κ² | κ̃₂ |
| exact ramp foot term (C_η, (A): 15.84 certified) | C_η with L = −κ | C_η^Ξ: L ranges over its values on the ramp |
| cutoff-free ramp data | J₂^free (verifier) | (A)'s energies E^A on the inner ramp [28,56] (inside (A)'s core, η_A = 1); Lemma R on the outer ramp |

No new class of terms appears and no new data are needed. In `lemmaS_L4terms.l4_extra_YC`, κ enters as `kap`. In `ghat_from_YC`, κ² enters as `kap**2`, so pass √κ̃₂ there.

## 4. The constants for the Lemma-U weight

The weight is Θ = W(σ) D(s), with
- σ = s − ln(ρ/λ(τ)),
- W = 1/(e^{−μ(σ+a)} + 1 + e^{ν(σ−ln1.5−a)}),
- D = 1 + g e^{−4s} (the damping factor),
- θ = g e^{−4s}/D, and θ_a = θ(ln 28) = g/(28⁴ + g).

Exact formulas (sympy, I6–I8):
- L = −4 + W′/W − 4θ, θ′ = −4θ(1−θ), and Ξ″/Ξ = L² + L′.
- W′/W = μp − νq ∈ (−ν, μ), where p, q are the normalised weights of e^{−μ(σ+a)} and e^{ν(σ−ln1.5−a)}.
- (ln W)″ = −[p(1−p)μ² + q(1−q)ν² + 2pqμν] ≤ 0, so W is log-concave.
- Θ^{(m)}/Θ = Σ_k C(m,k)(W^{(k)}/W)(−4)^{m−k}θ for m − k ≥ 1, plus the k = m term W^{(m)}/W.

**Bounds without using the geometry:**
- κ̃ ≤ max(4 + ν + 4θ_a, μ − 4).
- |Θ^{(m)}/Θ| ≤ (1−θ_a)Γ_W^m + θ_a(Γ_W + 4)^m.

At θ_a ≈ 0.39 the second bound would give Γ_eff ≈ 6.5. That is too pessimistic, because the damping factor is large only near y = 28, where W is on its inner tail.

**With the geometry:**
- The support starts at y = 28 (the foot of the inner ramp).
- ρ/λ(τ) ≥ 56 for the crude family (from launch on) and ≥ 150 for the fine family.
- On the inner tail with μ = 4, the factor e^{4σ}(1 + g y^{−4}) ∝ y⁴ + g, so the y^{−4} of the damping factor is cancelled: Θ^{(m)}/Θ = 4^m y⁴/(y⁴+g) ≤ 4^m, and L ∈ (−4, 0).
- On the outer tail, θ ≤ g/(84⁴e^{4a} + g) ≈ 4·10^{−3}.

Hence (t8, grid plus the exact σ → ±∞ limits; damping g fixed by the rate criterion g = Ĝ₂₂/(1 − Ĝ₂₂/28⁴), i.e. θ_a = Ĝ₂₂/28⁴, iterated to a fixed point):

| case | θ_a | g | Γ_m (m=1,2,3) | κ̃ (ν+4) | κ̃₂ | Ĝ_jj, no L⁴ | Ĝ_jj, with L⁴ | L⁴ factor (22) [power, same box] | L⁴ extras general/power | rate margin 28⁴/Ĝ₂₂ | Gw = 1/(1−θ_a) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| crude, mu=4,nu=4,a=0.3 | 0.3901 | 3.931e+05 | 4, 16, 64 | 8.000 | 64.00 | 4.93e+04, 2.303e+04, 7.845e+04 | 4.932e+04, 2.856e+04, 2.398e+05 | ×3.06 [×4.08] | ×2.75 | **2.564** | 1.639 |
| fine, mu=4,nu=4,a=0.3 | 0.1776 | 1.328e+05 | 4, 16, 64 | 8.000 | 64.00 | 5.053e+04, 2.196e+04, 5.466e+04 | 5.054e+04, 2.441e+04, 1.092e+05 | ×2.00 [×2.45] | ×3.64 | **5.629** | 1.216 |
| crude, mu=4,nu=2,a=0.5 | 0.3197 | 2.888e+05 | 4, 16, 64 | 6.000 | 36.00 | 4.3e+04, 1.907e+04, 7.024e+04 | 4.301e+04, 2.384e+04, 1.965e+05 | ×2.80 [×4.08] | ×2.15 | **3.128** | 1.470 |
| fine, mu=4,nu=2,a=0.5 | 0.1458 | 1.05e+05 | 4, 16, 64 | 6.000 | 36.00 | 4.43e+04, 1.82e+04, 4.817e+04 | 4.43e+04, 2.027e+04, 8.965e+04 | ×1.86 [×2.45] | ×2.77 | **6.856** | 1.171 |
| crude, mu=4,nu=4,a=0.0 | 0.3901 | 3.931e+05 | 4, 16, 64 | 8.000 | 64.00 | 4.93e+04, 2.303e+04, 7.845e+04 | 4.932e+04, 2.856e+04, 2.398e+05 | ×3.06 [×4.08] | ×2.75 | **2.564** | 1.639 |
| fine, mu=4,nu=4,a=0.0 | 0.1776 | 1.328e+05 | 4, 16, 64 | 8.000 | 64.00 | 5.053e+04, 2.196e+04, 5.466e+04 | 5.054e+04, 2.441e+04, 1.092e+05 | ×2.00 [×2.45] | ×3.64 | **5.629** | 1.216 |

References (crude box, no damping factor): power weight (γ=1, κ=5): Ĝ₂₂ = 1.906e4 → 7.772e4 with L⁴ (×4.08); previous T2 v1 (Γ=4, κ=5): 6.702e4 → 1.786e5. All [est] (float grid + exact σ → ±∞ limits). Every case has a fixed point and `budget_ok` (L⁴ dissipation budget) holds.

## 5. Cost of the general-weight L⁴ terms

- **The identity costs only constants** (§2). For the Lemma-U weight with μ = 4, the damping factor adds nothing to Γ_m, κ̃ or κ̃₂: the grid sups equal the σ → −∞ limits Γ_m = 4^m, κ̃ = ν + 4, κ̃₂ = κ̃². So the previous agent's shortcut (Γ = 4, κ = ν + 4) is exact, and the plateau offset a does not move any constant (rows a = 0 and a = 0.3 agree).
- **Rates [est].** Crude family, μ = ν = 4: Ĝ₂₂ = 2.398e5 (L⁴ factor ×3.06, *smaller* than the power weight's ×4.08 on the same box), rate margin 28⁴/Ĝ₂₂ = **×2.56**, θ_a = 0.390, Gw = 1.639. With ν = 2: ×3.13. Fine family (box 0.21): ×5.6–6.9; but the fine family's damping must use the crude Ĝ (its weight reaches y = 28), so downstream uses Gw = 1.639.
- The GONOGO.md v1 figure "T2 ×3.4" used κ = 5 and is superseded by ×2.56.
- **Inner-ramp L⁴ foot terms [est]** (crude, μ = ν = 4, a = 0.3): cutoff-free ramp data F_j = 2.67e-12, 1.8e-09, 1.06e-06; foot₁ = 4.58e-06, foot₂ = 0.00834, S33 = 0.000106; level-2 total **0.00845**, ≈1% of the level-2 bound B₂ ≈ 0.69 (t1). These use (A)'s certified power-weight foot weights scaled by the extras ratio ×2.75.
- **Effect downstream:** see t1/t5 CLEAN-GO rows and GONOGO.md.

## 6. What remains for the full proof (not done here)

- arb enclosures of Γ_m, κ̃, κ̃₂ and C_η^Ξ over the geometry (σ, s). The float values above are suprema approached only at σ → ±∞ or on grids.
- The ramp versions (c) in space–time form, with F₂ and J_G taken from (A): E^A on [28,56] and (A)'s level-1 dissipation. These follow (A)'s certified ramp code (`lemmaS_F2.py` foot/S33 blocks) with κ → κ̃ and C_η → C_η^Ξ.
- The class weights in arb with Γ_eff (as `rigor/weights_exact.py` does for γ).
