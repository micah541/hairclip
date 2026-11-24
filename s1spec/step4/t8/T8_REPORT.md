# T8 report: Lemma F — certified constants, amplification, Ȳ requirement

**Certified run (arb, polished base, slabs 5–187).** All runs used ≤ 2 processes at nice 10, with relative paths and subprocess workers (Windows-safe).
- `python cforms.py`: identity (I_j) with generic β, b₂, c₀, Euler-checked for j ≤ 3 (10 s).
- `python cconst.py --slabs 5-187 --out constants/slabs --jobs 2`: per-slab constants, restartable. 186 CPU-s in all (1.0 s per slab); 4.7 min wall on the shared machine.
- `python merge.py`: builds `constants/constants.json`, with the T7 inputs and their sources.
- `python camp.py`: Lemma J and Corollary C in arb, giving `constants/amp.json` (1 s). Option `--Cg0` for sensitivity, `--F0` for T4's initial data.

Statement and proof are in `LEMMA_F.md`, with `LEMMA_F.tex`/`.pdf` built by `mktex.sh`. The float study (`fb.py`, `lfnum.py`, `jcas.py`, `amp.py`) is kept as a cross-check (§6).

## 1. Lemma F: what changed for the certified version
- **Half-curve energies, far-end boundary term (§3.1).** F_j' gains 𝔅^D_j + 𝔅^S_j.
  - 𝔅^D_j = D·T_j(ℓ̃)w_j(ℓ̃)². Its coefficient is the transport defect **D = b₁(ℓ̃) + ℓ̃_τ = A(ℓ̃) + ε_ℓ, |D| ≤ 6.9e-13.** T6's D_end = ϖ(ℓ̃) + ℓ̃_τ (≤ 1.2e-9) is not this coefficient: its extra terms cq̃(ℓ̃) + b̃t̃₁(ℓ̃) cancel exactly in b₁ = −5κ̃κ̃' − cq̃ − b̃t̃₁ + ϖ.
  - **Lemma B** (trace inequality on the flattened plateau) gives |𝔅^D_j| ≤ ε_j[(2/y_ℓ + 5/4)F_j + D_j + ½∫(T'')⁺w_j²], with ε_j ≤ 6.5e-12. This is included.
  - 𝔅^S_j (the other fluxes) vanishes under the symmetric-end hypothesis (SE); see §4, flag 1.
- **Exact transport:** β := Q − cX0 − 6κ̃κ̃' + A (T6 §5), with A^(m) from the hulls of κ̃r, r, r', r'', r'''.
- **Two fixes to the float study.**
  - θ₃ was −0.01 for opt4n (δ_P = 1.5 together with the ½D_j of (S6)). The certified run uses δ_s = ¼ at j = 3, which gives θ₃ = 0.240.
  - k_jj := sup_R(q̂ + (λ_j + ΔF_j)T_j)⁺/W. The ramp share ΔF_jT_j of the Young steps was missing (+0.1%).
- **Small additions:** the e-self-coupling at j = 0 now lowers λ₀ by 2A₂ (≤ 1e-8); the residual terms of (I0) and (I2) are in Lemma J (item 6), assuming Ȳ ≥ 1e-12.

## 2. Method (no sampling)
- **Inputs.**
  - T6 `jets/`: hulls per element sub-cell over the whole slab, for κ̃⁽⁰⁻⁵⁾, Q, p, q, n, t₁, V*, r, r'–r''', A, κ̃r.
  - T6 `residuals.json`: 𝓡_j, the hat terms RFa/RFb (opt4n weights), R1, Rᴵs, ‖χr‖_ρ, A(ℓ̃), ε_ℓ, and the ranges of c, b̃, ℓ̃, Q(ℓ̃).
  - `t7/constants.json` (§4, flag 2).
- **Pieces.**
  - Fixed elements (x ≤ 40) are split exactly at y_a and y_b, and into 64 pieces per sub-cell on the ramp.
  - Moving elements use their hull. Lengths are bounded by h_e(τ₁)/4, which is valid because the knots increase in τ.
- **Core (η = 1).** Everything is divided by T analytically. The transport/∂_τT cancellation, ∂_τ ln T + cyτ₁ + (2−2j)c = (2j−2)(D − β(ℓ̃))(Φ − sΦ')/Y, is done exactly. This identity and the weight jets were checked against `weights.py` to 1e-16.
- **Ramp.** The foot ratios are written as bounded rationals in t: S^(l) = t^{4−l}σ_l, with σ₀ ≥ 1 handled by monotonicity.
- **Norm lower bounds.** ‖χp̃‖, ‖χñ‖ and ‖χV*‖_ρ come from cell infima, minus ‖χr‖_ρ. Mine are 3.5% below T7's certified nV_min; the larger of the two is used.
- **Integration.** Corollary C bounds are valid for all τ in each slab: decay is taken from slab ends, sources are slab sups, and inflow uses running sups and λ* = running min.
- **Ȳ/residual split.** √F_j ≤ Ȳa_j + b_j, using (u+v)² ≤ (1+10⁻³)u² + 1001v².
- **Hand-off.** The good-time window is slabs 184–187 (k ∈ [2.03e4, 2.36e4]). All items use sups over the window. Item 4 is solved by bisection, and the final value is re-checked in arb.

## 3. Results
**Constants** (certified, against the float study over slabs 5–187):

| j | λ_j min (float) | k_jj max (Δ vs float) | θ_j | γ_j max (× float) | max c_ji |
|---|---|---|---|---|---|
| 0 | 0.2387 (0.2389) | 9.70e4 (+2.0%) | 1.245 | 1.0e2 (×1.26) | – |
| 1 | 0.2488 (0.2488) | 3.075e7 (+1.5%) | 1.245 | 5.6e2 (×1.18) | 7.0e-3 |
| 2 | 0.2485 (0.2485) | 1.027e10 (+1.6%) | 1.245 | 2.0e4 (×1.79) | 0.33 |
| 3 | 0.2452 (0.2457) | 1.135e12 (+3.9%) | 0.240 | 7.9e5 (×2.45) | 10.4 (bend, k ≈ 24) |

**Amplification and Ȳ requirement** (window sup; √F_j ≤ A_jȲ + B_j):

| T7 inputs | A₀ | A₁ | A₂ | A₃ | **Ȳ item 3 (E₃)** | Ȳ item 4 (K4, good time) | item 2 |
|---|---|---|---|---|---|---|---|
| **t7/constants.json as is** (C_γ0 max 1.72, provisional) | 6.7e3 | 7.6e4 | 1.6e6 | **1.46e8** | **1.60e-11** | 3.4e-11 | ≥ 6.4e-8 |
| same, C_γ0 = 0.096 (T7's component value at slab 100, expected from c2r3) | 6.7e3 | 6.1e4 | 1.1e6 | 1.09e8 | 2.15e-11 | 4.8e-11 | ≥ 6.4e-8 |
| float study with T7 prototype constants (§6) | 5.8e3 | 5.5e4 | 7.7e5 | 8.9e7 | 2.67e-11 | 6.3e-11 | – |

- **Residual parts.** B₃ = 1.6e-5 (0.65% of the item-3 allowance) and B₀ = 1.4e-6; both are already subtracted in the limits.
- **A₃ split** (C_γ0 = 0.096):
  - ramp inflow √(k₃₃J₃) = 8.4e7, with J₃ = 6.2e3;
  - modulation G₃ = 2.5e7;
  - Duhamel part 5.6e4.
- **J₀..J₄ at the end:** 112, 59, 21, 6.2e3, 1.81e6. The routes are Λ / (I0) / (I2). 74% of J₄ is the c₃√J_{D4} (1/√ε) term.
  - With C_γ0 = 1.72, J₂ = 90 (the (I0) modulation term) and J₃ = 1.28e4.
- **Requirement.** **Ȳ ≤ 1.6e-11 certified with T7's current file** (2.15e-11 once C_γ0 is final). The budget is Ȳ ≈ 2.3e-12 (T7), so **the margin is ×7.0 (×9.4)**, and no fallback is needed.
- **Initial data.** The window bound gains e^{−Λ_j/2}√F_j(τ₀) ≈ 0.036√F_j(τ₀) (F₀ = 0 in the run; T4 to supply, via `--F0`).

## 4. Flags
1. **Far end, hypothesis (SE).** The base is exactly symmetric about its own far-end normal line N; the true curve is symmetric about the axis. N is off the axis by |q̃(ℓ̃)| ≤ 4.6e-9 (≈ 1.3e-14·ℓ̃) and tilted by ≤ 2.4e-13; the turning defect is ≤ 1.65e-14. This looks like a fixed rigid rotation of about 1.6e-14 rad. It enters only through the error's parity at ℓ̃ (𝔅^S_j), and Lemma F assumes (SE).
   - To restore (SE), P2 should compare with the rotated base. Rotations are symmetries up to the axial drift, |δr| ≤ 1.1e-15 (for T6); the tip moves by ≈ 4e-14, which is T4's start tolerance. This is for T5/T9.
   - A direct bound on 𝔅^S would need far-end traces of w_{j+3}. The effect is nil (mode energies ≤ 1e-20 relative), but this is an input, not a proof.
2. **T7 inputs.**
   - From `t7/constants.json`: μ₁ = 0.15, μ₂ ≥ 0.6415, a_f = 0.0696, b_f = 0.0394, c_I = (0.425, 0.449, 0.539, 0.448), Ψ₁ = 0.3885, c₁W = 0.0943, C_γ2 ≤ 0.029.
   - C_x = 7.90 is the max of the 4 slabs that have C3; 179 slabs are pending c2r3.
   - C_γ0 is provisional (C_x at τ_c) and costs ×1.34 in Ȳ.
   - μ_e = 0.036, ϱ_c = 0.043 and K_e = 1.4e4 are not in the file: they are prototype placeholders.
   - When T7 updates: `python merge.py && python camp.py` (about 2 s; cconst only reads nV_min from T7).
3. **Conservatism of γ_j.** On the moving outer cells the slab hulls contain the knot motion, which makes γ₃ ×2.45 the float value; G₃ is about 23% of A₃. Jets resolved in τ (T6's engine with τ-subintervals, about 1 CPU-h) would recover about 10–15% of Ȳ. This is not needed.
4. **Lemma N3′ (T7 flag 4; STEP4_RESULT gap 8; review A1, B2, B4, B6).** Replaced by LEMMA_F §6.1. The exact operator's nonlinear remainder has no x-growing coefficient (T7's x²n₂ growth is an interval-hull artifact), and its weak form has no trilinear classes. Levels 0–3 and 𝒬 are certified per slab by `p2/t8se/n3far.py` → `camp.py --n3`. The boxes assumed there are checked by `n3boxes.py`; the margins at Ȳ = 1e-12 are ×2.76, ×3.97, ×1.69 and ×1.70. The modulation bound now includes the residual part of (M). C5 is corrected: Lemma O's H4 needs the good time chosen in step 4. With W′ (|W′| = 1.46e-4), H4 fails by ×16–279; with one slab it holds (≤ 0.31·M², `h4loc.py`). Requirement: 4.99e-12 on W′, and 1.79e-11 (item 3) on a one-slab window. Refresh after T7's flags: the chain now runs on the outward-rounded SE data from se_corr.py, with x_r taken from the 2-norm of c_φ (x_r ≤ 2.6e-16). Every reported number is unchanged to 3 digits. **Review B1:** the residual of the (SE)-corrected base is now computed directly by T6's engine (`p2/t6c/`, LEMMA_F §6.1). On these data the |e′| margin is ×6.43 (was ×2.76), the W′ requirement is 5.15e-12 and the (Y) level is 3.2e-14. se_corr's 𝓡₀ increment was not valid: it left out the shift of S₁(ℓ).
5. **Other inputs outside T8:**
   - T1's base window energies and T9's chart factor 1.02 (as in the float study);
   - the bootstrap hypothesis Ȳ ≥ 1e-12 (Lemma J item 6);
   - slabs 0–4 stay outside P2 (T6 starts at slab 5).

## 5. Lemma N3′, Lemma J, consistency
- **Level 4 is not used anywhere.** Item 4 uses D₂ at a good time (C5), and the nonlinear level 3 uses N3′. The certified run has j ≤ 3.
- **λ < 2μ₁.** This holds with margin: λ* ≤ 0.2527 against 2(μ₁ − 10⁻³) = 0.298.
- **Agreement with the float study.** The λ_j agree to 2e-4 and k_jj to within 4%. The T7-prototype requirement (2.67e-11) is reproduced by the float code, and the certified value differs from it only through the three items above (C_γ0, the γ hulls, C_x 7.9 against 7.5).

## 6. Float study (superseded, kept for reference)
The float study used fb.py (the float T1 base, whose far-field κ⁽⁴⁾ and κ⁽⁵⁾ are noisy: T6) and 26 points per element.
- **T5 design [14,28]:** A₃ = 1.8e8; Ȳ = 1.35e-11 (good time).
- **Recommended opt4n:** A₃ = 1.2e8; Ȳ = 2.0e-11 with placeholders, and 2.67e-11 with T7's prototype constants.
- **Levers.**
  - The (1+x²)⁻¹ weight lever is withdrawn (T7: ×1.17 at most).
  - Remaining levers: T5 fallback 1 restricted to D₄ (×2.6); E₃ ≤ 3e-5 (×1.9); T7's ζ-ramp moved in (×1.4).
- **What dominates levels 3–4.** Ramp inflow √k₃₃ (the price of measuring the flux in W = (1+x²)⁻² at y ≈ 20), then the 1/√ε loss in J₄. The bend cross terms c_ji are negligible.

## 7. Rerun on the v2 base (t1/polished_v2, slabs 0–180), done by T7's `run_all.py --t8`
- **Code changes, slabs and τ only, plus T7 inputs:**
  - `merge.py` takes its slab range from `t6/residuals.json`.
  - `camp.py` takes the window as the trailing slabs with k ≥ 2.0e4 at their start: slabs 177–180, k 2.04e4 → 2.37e4 (184–187 on v1). It checks that slabs are contiguous with length 0.15, and F₀ is at the start of slab 0.
  - `cconst.py` reads μ_e, ϱ_c and K_e from `t7/constants.json` when present (K_e = 890.6 replaces the 1.4e4 placeholder).
  - On v1 data, merge and camp reproduce Ȳ = 1.964e-11 exactly.
- The v1 outputs are kept as `constants/slabs_v1`, `constants/constants_v1.json` and `constants/amp_v1.json`.
- **Result: Ȳ requirement 2.092e-11** (binding item 3, E₃).
  - T7 inputs: μ_e 0.036, ϱ_c 0.0471, K_e 890.6, C_γ0 0.207, C_γ2 0.0159, μ1 0.15, C_x 7.82.
  - λ* ≤ 0.2520 (certified) against 2μ1 = 0.30.
  - ϱ_c does not affect Ȳ: it enters only through (ϱ_c − λ*/2)₊ = 0.
- BASE_T1 and KCH are still the v1 values.
