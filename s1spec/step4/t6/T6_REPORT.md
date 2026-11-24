# T6 report: rigorous residual enclosures for P2 (slabs 5–187)

Code in `t6/` (arb, 128 bits, python-flint 0.9; portable, no fork, T1's `polstore.load_slab`):
- `tp.py`: Taylor models in τ; bivariate (s, τ) ball polynomials.
- `base.py`: the base as an exact object.
- `engine.py`: the per-element residual.
- `weights.py`, `slab.py`.
- `run.py`, `merge.py`, `crosscheck.py`: the command-line steps.
- `reader.py`: reading the outputs.

`python run.py --slabs 5-187 --out out --jobs 2` (restartable) → `python merge.py --in out` → `python crosscheck.py --cert out`.

## 1. Method (no sampling in the certified path)
- **Base.** On slab i each unknown is the degree-8 polynomial in τ through the stored stage balls. The nodes are the Radau nodes, taken as exact rationals; c₈ = 1, so slabs join exactly. Per element it is T1's polynomial in s = 2u − 1. The outer knots XT(1 + Cℓ̃³)^{1/3} are Taylor models in τ. All objects are polynomials in (s, τ) with ball coefficients; products are exact, and every truncation folds into the radii.
- **Residual.**
  - ρ_κ = κ̃_τ|ₓ − Φₓ and ρ_θ = ∫ρ_κ are exact. They equal T1's arb residual, and the Leg₀ moments are below 1e-37 at the nodes.
  - Re-derived identity: r + iA = e^{iθ}∫₀ˣe^{−iθ}ρ_θ, so that Ṽ = V* + r and ϖ = ϖ* + A, where V* = −κ̃'' + cp̃ + b̃ñ and ϖ* is T1's identity.
  - e^{iφ} per element: power series plus a Cauchy tail.
  - r′, r″, r‴ from (r, A)′ = (ρ_θ − κA, κr); S₁ and S₂ by exact integration.
- **Norms.** Per cell (≤ 0.25 on [0, 40], split at the weight breakpoints; 2 per outer element) I form the exact polynomial P_c(τ) = ∫_c f² dx.
  - Weights are bounded on cell × ¼-slab: monotone end values; interval bounds for T′²/T and |T″|. There are no violations against T8's weight jets.
  - sup_τ Σ w̄_c P_c is bounded by a 32-piece Taylor-shift subdivision.
  - `norms_int` = (∫_slab N² dτ)^{1/2} uses the same bound.
- **Definitions.**
  - T5 (R) for R2, Rs, R1 and Rᴵs.
  - Lemma F 𝓡_j: T8's flattened T_j, ∂⁻¹r = S₁ − S₁(ℓ̃), designs `opt4n` and T5 ramps.
  - `RFa`/`RFb`: residual parts of 𝓡̂_j.
  - `RFlit`: T5/resid.py literal form.

## 2. Results: max certified sup over slabs 5–187; **no slab fails**

| norm | target | max (slab) | margin |
|---|---|---|---|
| R2 | 2e-14 | 5.70e-16 (5) | 35 |
| Rs | 5e-11 | 4.21e-14 (5) | 1.2e3 |
| R1 | 5e-13 | 6.75e-15 (6) | 74 |
| Rᴵs | 1e-10 | 3.06e-18 (5) | 3e7 |
| 𝓡₀ opt4n (T5 ramps 5.43e-9) | 1e-7 | 5.82e-9 (182) | 17 |
| 𝓡₁ / 𝓡₂ | 1e-7 | 5.80e-11 / 3.69e-10 (182) | 1.7e3 / 270 |
| 𝓡₃ | 1e-7 | 5.81e-9 (182) | 17 |
| 𝓡₄ | 1e-5 | 1.43e-7 (182) | 70 |
| T5-literal RF0–RF4 | same | 2.4e-11 … 6.9e-9, 1.9e-7 | ≥ 15 |

- Hat terms (max): RFa/RFb₀ 1.8e-9; RFa₄ 1.0e-9; RFb₄ 1.9e-9.
- Slab 182 has two knots inserted at its start. It carries every Lemma F maximum. Elsewhere: 𝓡₃ ≤ 2.5e-9, 𝓡₄ ≤ 6.0e-8, 𝓡₀ ≤ 1.2e-10.

## 3. Cross-check
Against resid.py on T1's 9 sample slabs plus slab 182, with 41 equispaced times plus the nodes: **certified ≥ estimate in all 90 cases.** The ratio lies between 1.001 and 1.091, with the largest values at RF4 and Rᴵs; `crosscheck.json` has every ratio. The jets contain T1's arb values of the polished κ⁽⁰⁻⁵⁾ and Q at 5 times × 17 points per element.

## 4. Findings
- **Time structure.** Residuals peak at τ₀, where the stage-polynomial derivative jumps. Between the nodes the far-field ρ_κ is ~10³× its nodal value (τ-degrees 1–8), so the sup must use the exact τ-polynomial. The integrals are smaller: (∫R2²)^{1/2} = 4.2e-17 against √Δτ·sup = 2.2e-16 at slab 5.
- **Level 0.** S₁(ℓ̃) = ∫₀^ℓ̃ r (≤ 6.5e-9; ≤ 1.3e-10 off slab 182) dominates 𝓡₀.
- **Far end (T8, T4/T9).** The defects are:
  - closure |q̃(ℓ̃)| ≤ 4.6e-9;
  - |t̃₁(ℓ̃)| ≤ 2.4e-13;
  - scalar ODE defects ≤ 2.3e-14;
  - D_end = b₁(ℓ̃) + ℓ̃_τ ≤ 1.2e-9.

  Lemma F's closed-curve reduction assumes these are 0, so it needs a boundary term of relative size ≲ 1e-13.
- **T8.** `fb.py` samples the *float* base. In the far field its κ⁽⁴⁾ and κ⁽⁵⁾ carry element-scale noise up to ~50× the polished values. Use the polished base.

## 5. T7's open corrections: clock direction and b₁
T7 used Ṽ_s = cp̃ + b̃ñ − κ̃'' and b₁ = −6κ̃κ̃' + Q − cx. Exactly, Ṽ = Ṽ_s + r and b₁ → b₁ + A, so b₁′ gains κr. The certified maxima are:
- **g₃ = Ṽ/n.**
  - |n − n_s| ≤ ‖χr‖_ρ ≤ 4.2e-14 (`Rs_chi`).
  - ‖χ(g₃ − g₃,s)‖_ρ ≤ 2‖χr‖_ρ/(n_s − ‖χr‖_ρ).
  - γ₀ gains Rᴵs/n ≤ 3.1e-18/n, and γ₂ gains ‖ζr‴‖_W/n ≤ 2.7e-15/n.
- **(C1).**
  - δ⟨𝓛u, u⟩_ρ = −½∫(ρA)′u², with sup over [0, 34.25] of |A′| ≤ 1.2e-10 and |Aρ′/ρ| ≤ 9.9e-13.
  - ‖δℓ_j‖_ρ ≤ 1.3e-10‖φ_j‖_ρ + 2.7e-12‖φ_j′‖_ρ.
- **(C4),** for x ≥ 10.25: |A|/x ≤ 2.0e-13 and |A′| ≤ 2.3e-12.

All of these are ≥ 1e7 below T7's 1e-3 margin. For T8, β → β + A, with A⁽ᵐ⁾ from the hulls of κr, r, r1–r3 and κ⁽ᵐ⁾.

## 6. Runtime, storage, how to read
- **Runtime.** 3054 CPU-s: 13–21 s per slab, mean 16.7 s. Wall time was 1 h 45 min with 2 workers at nice 10, sharing the 2 cores with 3 other jobs. The same command (`--jobs N`) runs on Windows.
- **Storage.** `residuals.json` 1.0 MB; `jets/` 39 MB (183 npz files plus `index.json`); `out/` 1.5 MB.
- **`residuals.json`.** `per_slab[k]` holds `norms` (certified sup, rounded up), `norms_int`, `scalars` and `t7`. `summary` and `definitions` cover every name. `reader.py` provides `norm`, `interval_max(name, s0, s1)`, `jets` and `hull`.
- **P2 driver.** Take `interval_max` over each certificate interval: R2, Rs, R1 and Rᴵs for (Y), `RF{j}_opt4n` for Lemma F. The hat terms satisfy 𝓡̂_j² ≤ (8/a)(RFa_j + |x₃|γ̂ᵃ_j)² + ν⁻¹(RFb_j + |x₃|γ̂ᵇ_j)². In time integrals, `norms_int`² may replace Δτ·sup².
- **Jets.** `jets_XXX.npz` stores (N_el, 4, 2) arrays of [lo, hi] per element sub-cell over the slab, for x, κ⁽⁰⁻⁵⁾, Q, θ, p, q, n, t₁, Ṽ, ϖ, b₁, V*, ϖ*, r, A, r1–r3 and κr. For tighter constants, `base.SlabBase` + `engine.element_a` rebuild the exact polynomials of the same base.
