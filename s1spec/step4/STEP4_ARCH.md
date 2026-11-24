# Step 4 architecture (recommendation, Nov 9)

Scope: design only. Inputs read: step3_proof §1–2 (H1–H6), lemmaS p.1, lemmaO_v3 Revision, STEP4_DESIGN, FEAS, FEAS2, hls.py.
One new float check was run for this note (`arch/jaccheck.py`, output `arch/jaccheck.out`, §0). Everything numerical below is float unless marked.

## Decision in one paragraph
Use **two engines, one switch**.
(I) From t = 0 to k_max ≈ 10–14 (t ≈ 4.598): a **Fourier validated integrator in the HLS variables** (physical frame, no modulation, no S1 weights), with a small Lohner block for the low modes and an analytic tail, in arb.
(II) From there to the hand-off: **energy shadowing of step-3 type** around a tip-adapted numerical trajectory in the rescaled frame, 3 distinguished directions, about 5–8 coercivity certificates.
The 42 macroscopic certificates of FEAS2 (N up to 76, margins down to 0.001, frame-dependent counts 4–16) are dropped. The initial layer disappears by choosing the initial curve smooth. The hand-off is a base switch by the triangle inequality, not a glued bootstrap.

## 0. New float check: the macroscopic phase in the HLS Fourier frame
HLS in the time τ' (dτ' = (2π/L)⁴ dt) is θ_τ' = −(∂⁴+∂²)θ + N(θ), with N cubic, nonlocal, at most 3 derivatives and **L-free** (L and t are slaved quadratures; this is already how hls.py steps). I linearised it by finite differences on the lowest 260 modes of the feas2 snapshots. Rates are per unit t.

| snapshot | k_max | L | modes for 1e-13 | eig(J): n₊, max | sym L²: n₊, max | sym H¹: n₊, max |
|---|---|---|---|---|---|---|
| t=0.002 | 2.8 | 103.6 | 450 | 9, 25 | 26, 112 | 7, 49 |
| t=0.01 | 3.4 | 103.1 | 397 | 10, 9.4 | 24, 71 | 6, 32 |
| t=0.1 | 3.1 | 99.6 | 366 | 10, 2.1 | 21, 50 | 6, 23 |
| t=0.5 | 2.7 | 90.2 | 289 | 7, 1.3 | 17, 33 | 6, 12 |
| t=2 | 2.2 | 70.4 | 188 | 5, 1.1 | 12, 23 | 6, 5.5 |
| t=3 | 2.0 | 61.3 | 152 | 3, 1.2 | 11, 21 | 6, 3.7 |
| t=4.5 | 2.8 | 52.1 | 165 | 3, 7.1 | 11, 31 | 4, 14 |
| t=4.59 | 5.0 | 50.0 | 278 | 3, 102 | 14, 159 | 3, 147 |
| k=10 | 10 | 48.9 | 473 | 3, 1.8e3 | 22, 2.3e3 | 3, 2.2e3 |
| k=20 | 20 | 48.4 | 928 | 3, 2.9e4 | 36, 3.5e4 | 3, 3.4e4 |

Reading:
- The macroscopic phase is a small problem in this frame: 150–450 Fourier modes (half that in the symmetric sine class), 3–10 frozen unstable directions with rates O(1) per unit t, and ∫(max rate) dt ≈ 6–7 over [0.01, 4.5]. So the worst-case transverse growth is ≲ 1e3 (the measured growth modulo symmetries is ≤ 0.4, FEAS2).
- The late growth, 0.18 k⁴ = 1/(T−t), is the T-shift mode. It contributes a further factor ≈ 180 from t = 4.5 to k = 10, in one direction only.
- Already at t = 0.002 the curve has k_max = 2.8 and needs 450 modes. The k = 42 initial layer is gone by then.
- Caveats: frozen-time eigenvalues, finite differences, full (not symmetric) class, 260-mode truncation.

## 1. Phases and methods

| phase | range | frame / clock | norm | error propagation | modulation |
|---|---|---|---|---|---|
| P0 initial curve | t = 0 | physical | — | exact data (tiny interval family) | — |
| P1 macroscopic + approach | t ∈ [0, t_a], k_max ≤ k_a ≈ 10–14 | physical HLS θ(α), sine class; clock τ' | weighted ℓ¹ of Fourier coefficients (C³ control) | Lohner block (M ≈ 100–200 low modes) + analytic tail, per Chebyshev-in-time slab | none |
| P2 rapid | k_a → k₁ (≥ 6600) | rescaled, tip-centred; clock τ (dτ = λ⁻⁴ dt) | step-3 energies of the error: E (ρ, H²), E_I (W), then far-field energies | energy shadowing, coercive tail, 3 directions | time shift, axial translation; third direction as a scalar |
| P3 hand-off | at t₁ | step-3 chart over S1 | H1–H6 as stated | triangle inequality: enclosure of ũ + error ball | IFT for (x₁, λ₁, α₁) |

**P0. Initial curve.** Define Γ₀ by a finite sine series for θ₀(α) (about 250 coefficients, copied from the float state at t ≈ 0.001–0.002, or from a clip with rounded corners) and a length L₀. The two closure conditions are solved by interval Newton on two coefficients, so Γ₀ is a member of a tiny interval family. Embeddedness is certified in arb: local injectivity from |κ|·ℓ < π on arcs, and positive distance between non-adjacent arcs by subdivision. Take the state well before the crossing at t ≈ 0.0028 so that the slot gap has a comfortable margin. This removes the uncovered layer t < 0.01 (k = 42, coercivity −5.5e5 in FEAS2).

**P1. Macroscopic + approach: Fourier, Lohner block + tail.** Why:
- The equation is semilinear with a constant-coefficient diagonal principal part and polynomial nonlinearity. This is the standard setting for rigorous parabolic integrators (self-consistent bounds with a Lohner block, or Chebyshev-in-time Newton–Kantorovich with a finite propagator block). No frame, weight or coercivity design is needed.
- A propagator block captures non-normal transients that energy estimates cannot. The sym-L² count (11–26) and the frame-dependent counts in FEAS2 (4–16) show that energy shadowing would pay for non-normality here.
- The tail contracts once (2πn/L)² ≫ κ_max². That gives M ≈ 100 at L = 100, and M ≈ 200 at k = 10.
- Keep the block error as a parallelepiped (QR). The long T-shift axis is then absorbed exactly at the switch.
- About 300–400 slabs (≈ 30 for t < 0.01, 250 macroscopic, 50 approach).
- End at k_a = 10–14. Choose k_a so that the first P2 certificate has margin ≥ 0.05. At k = 10 FEAS2 has N = 29 with margin 0.02–0.08, which is thin. Raising k_a to 20–30 costs only more modes (N ≈ 47 k).

**P2. Rapid phase: energy shadowing around ũ.** Why: FEAS F3 gives exactly 3 non-coercive directions for all k ≥ 10, and F2 gives no amplification. FEAS2 gives certificates N = 29 → 6 → 3 and one certificate for all k ≥ 1600. Fourier is impossible here (N ≈ 47 k).
- Base: half curve, exactly symmetric, rescaled by λ(τ) with a normalisation (κ_tip λ = const). Tip-adapted geometric spectral elements over the whole half curve (300–500 DOF, +50% for the residual).
- Error e: normal graph over ũ.
- Directions 1 and 2 (time shift, axial translation) are exact symmetries. Compare the true solution with a shifted and translated ũ. Their modulation equations create no residual, and their size is unconstrained (H1 only needs the IFT).
- Direction 3 is not unstable (eigenvalue −0.256; only the symmetric part is +0.02 to +0.07). Carry it as a scalar a₃ = ⟨e, f₃⟩_ρ in a 2×2 comparison system with ‖q‖, as step 3 does for α. Closure needs c₁₂c₂₁ < μ₃c.
- Energies: the step-3 E and E_I applied to e, with the same ρ, W, ε, ε_I, so the step-3 coercivity code and constants are reused with L_ũ = L_S1 + certified perturbation. Between certificate times use the P2 spacing bound of FEAS2.
- Far field (beyond the W-zone): physical local energies on dyadic annuli. An annulus of radius R that leaves the core at time τ_R has remaining time ≤ (R/Y)⁴ in its own units, so its Gronwall factor is O(1). No spectral information is needed.

**P3. Hand-off.** See §3 and §4.

## 2. Validated numerics engine

**Arithmetic.** arb (python-flint 0.9.0 is installed) for every rigorous quantity. Floats only for approximate inverses and for producing ũ. No outward-rounded doubles and no double-double: the problem sizes (≤ 700 DOF per slab) make arb cheap, and it removes the rounding question. Step 3's own remaining "IEEE double without outward rounding" assembly (step3_proof §6.1) is a separate open item.

**P1 engine.**
- Unknown: sine coefficients of θ − α on a slab, Chebyshev degree 8–12 in τ'.
- Base: float predictor from hls.py, then Newton refinement of the slab collocation system with a float Jacobian inverse and an arb residual, at 128 bits, down to a residual ≲ 1e-22.
- Convolutions by arb_poly products (no FFT needed at this size).
- Bounds: Banach-algebra estimates for the cubic nonlocal term (it has the form θ_α ∂⁻¹(θ_α θ_ααα − mean), integrated by parts once), semigroup bounds e^{−(n⁴−n²)τ'} for the tail, interval fundamental matrix for the block.
- L, t and the base point X₀ are enclosed as quadratures afterwards.

**P2 engine.**
- Space: unknown κ(ξ, τ) on geometric elements, element size ∝ 1 + σ, degree p ≤ 12–16 (prefer more elements over higher degree). Use an integrated basis so that κ is exactly C¹ across elements. θ, the position and (p, q) = (Γ·ν, Γ·t) are then enclosed by validated quadrature. Exact symmetry at the tip is built into the basis.
- Time: Chebyshev degree 6–10 per slab, slab length Δτ ≈ 0.1–0.2, about 190 slabs for k = 10 → 6600. Each further doubling of k costs about 2.8 units of τ.
- Base: a native float collocation solver in this representation, seeded by HLS where HLS reaches, then iterative refinement with arb residuals. Do **not** interpolate HLS snapshots. Their floors are too high: T varies by 1e-5 with the step sequence, the odd asymmetry of 4.3e-9 (physical) grows like k in λ units, and filter and time-stepping errors are 1e-11 to 1e-12.
- Residual: r = (normal velocity defect of ũ) on each element × slab cell, enclosed by arb evaluation with subdivision and mean-value forms. Use the **weak form**: the error energy at L² level pairs r with e″ after two integrations by parts, so it needs only κ̃ (one derivative of θ̃). The ε-level (H²) needs r in L², i.e. θ̃‴.
- Clock: τ with d ln λ/dτ = −(¼ − a(τ)); ln λ is an unknown. t never appears, so the k ≈ 2500 clock floor is avoided. t₁ = t_a + ∫λ⁴ dτ is an enclosure computed once. s = −ln(T − t) agrees with τ up to O(a) but would need T, which is not known a priori.

**Certificates.** Port feas2/p1 (Galerkin + tail relative to S1) to arb: 5–8 certificates (N = 29, about 15, 6, 3, plus spares).

## 3. Error budget (λ₁ units at t₁; k₁ = 6600, λ₁ ≈ 2.3e-4, |α₁| ≈ 4.5–5e-4)

**Tip (H1).** √Λ₁ = 5e-8, truth ≈ 45|α|³ ≈ 4–6e-9. Allow 2e-8 for the validated error, split equally:

| quantity | allowed |
|---|---|
| ‖e‖_ρ | 8e-9 |
| ‖e″‖_ρ (ε = 1e-5) | 2.5e-6 |
| ‖ζe′‖_W (μ = 0.08) | 2.8e-8 |
| ‖(ζe′)″‖_W (ε_I = 3e-5) | 5e-6 |
| α error | 2e-5 (choose k₁ so the truth is 4.5e-4) |

Sources:
- **Incoming ball at k_a.** The rigorous contraction from the certified margins (0.02–0.08 near k = 10, 0.14–0.18 beyond k = 100, over Δτ ≈ 26) is only about 0.06–0.1. The float truth is 1e-3 to 1e-4, but it cannot be used. So require ≤ 1e-9 in rescaled units at k_a, i.e. ≈ 1.5e-10 physical in C².
- **P1 delivers** residual 1e-22 per slab × 400 slabs × growth ≤ 1e3 × 180 ≈ 1e-14. Margin 4 orders.
- **P2 residual.** Steady state ‖e‖ ≤ ‖r‖/c. Targets: weak residual ≤ 1e-11, strong L²_ρ residual ≤ 1e-8.
- **Nonlinear terms.** O(e²) ≈ 1e-16, negligible. a₃ must stay ≲ 1e-7 so that its coupling does not dominate.

**Far field.**

| condition | threshold | needed error | comment |
|---|---|---|---|
| H2, y ∈ [14, 60] | E₀,₁,₂ ≤ 2e-7, 1e-4, 0.2 | slope error ≲ 1e-4 | 3–4 orders looser than the tip |
| H3, strip and windows to 350 | E₁,₂,₃ ≤ 3e-5 (norm 5.5e-3) | ≲ 1e-3 in levels 1–3 of the slope | error is scale-invariant on exit from the core (≈ 1e-7), so derivative constants up to 1e4 are affordable |
| H4 | \|f\| ≤ ½, y\|f′\| ≤ 0.12 | 1e-3 | enclosure of ũ + error |
| H5 far jets, r ≥ 226 λ₁ | 8.9e-3, 7.0e-5, 6.2e-7 = physical κ ≤ 38, κ_s ≤ 1.3e3, κ_ss ≤ 4.9e4 | O(1) relative, C⁴ | expected truth: κ ≈ 2–4 with margins ≳ 10 (not measured) |

The far field is limited by the truth versus the thresholds, not by the validation error (risk R1).

**Is double precision enough?**
- Tip, weak form: yes. The floor is p²·1e-16 ≈ 1e-13 against 1e-11.
- Tip, strong L² residual: marginal at p = 20 (p⁶·1e-16 ≈ 6e-9 against 1e-8), fine at p ≤ 12 (3e-10).
- Far-field levels 2–3 by strong-form energies: no. They need θ̃⁽⁵⁾, and the floor p¹⁰·1e-13 is O(1) at p = 20. FEAS F1 already shows θ⁗ wrong by 2.4e-6 to 9e-5.
- Macroscopic phase: no, by 2–3 orders if worst-case growth is assumed (1e-13 × 400 × 1e3 ≈ 4e-8 against 1.5e-10).
- Recommendation: the float solvers stay in double; the base is polished by iterative refinement against arb residuals at 128 bits. This is cheap and makes every item above attainable. For the far-field derivative levels use interior parabolic smoothing of the error as a second line (Lemma-O-type cylinders, which have 4 orders of slack).

## 4. Cheaper alternatives

- **Different initial curve.** The crossing cannot be avoided (STEP4_DESIGN §2), but only P0 needs embeddedness.
  - Adopt: the smoothed clip (P0). It removes the initial layer at no cost.
  - Do not chase a shorter macroscopic phase now. With the Fourier engine it costs compute (about 250 slabs of the same code), not design effort.
  - Optional float experiments, one day, off the critical path: (a) arm 12 instead of 24, which halves L and the mode count if it still singularises (design logs exist); (b) the state at t ≈ 4.3 with the two crossing strands pulled apart and smoothed, which would cut the macroscopic time by about 90% if the flow re-crosses and still converges to S1. Untested.
- **Gluing steps 3 and 4.** No full glue. Step 3 is closed with fixed constants around u₀(α); re-deriving it around ũ would reopen it. Partial glue only:
  - P2 uses the step-3 energies, weights and coercivity code for the error.
  - At t₁ write w = (ũ − u₀(α₁)) + e. The first term is an explicit arb enclosure of the step-3 functionals of ũ; the second is the error ball.
  - k₁ is a free parameter. Later hand-off costs 15–30 slabs per doubling and no new certificate, so pick k₁ where every H has margin ≥ 3 in float.
- **Rejected.**
  - Energy shadowing through the macroscopic phase: 42 fragile certificates, weight redesign.
  - Uniform Fourier through the rapid phase: N ≈ 3e5.
  - Lohner block in the rapid phase: unnecessary with 3 directions.

## 5. Task list (S ≈ 1–2 days, M ≈ 1 week, L ≈ 2–3 weeks)

| # | task | size |
|---|---|---|
| T0 | Float hand-off audit: all H1–H5 quantities and margins at k = 3000–6600 (and beyond) from an exactly symmetric trajectory; fix k₁ | S |
| T1 | Native float solver: rescaled frame, tip-adapted elements, clock τ, symmetric, arb-refined; k = 10 → k₁; residual report (weak and strong) | M |
| T2 | Initial curve Γ₀: series, closure by interval Newton, embeddedness certificate | S |
| T3 | Fourier validated integrator (sine class, cubic nonlocal term, slabs, Lohner block + tail, arb), with its bounds note; run t = 0 → t_a | L |
| T4 | Frame-change lemma at t_a: Fourier ball → energy ball over ũ; extraction of time shift and translation | M |
| T5 | Shadowing theorem for P2 (paper proof): error equation over a moving base, 3 directions, nonlinear terms, weak residual, perturbation between certificates | L |
| T6 | Residual enclosure engine (arb, space-time piecewise polynomials, ρ and W norms) | M |
| T7 | Coercivity certificates in arb (port of feas2/p1), spacing bound | M |
| T8 | Far field: dyadic physical energies, derivative levels, H2–H5 error bounds, jets | L |
| T9 | Hand-off assembly: enclosure of the step-3 functionals of ũ, IFT for (x₁, λ₁, α₁), H1–H6 check, parametric in the step-3 thresholds | M |
| T10 | Production runs, logs, write-up | M |

**Total.** 3 L + 6 M + 2 S ≈ 15–16 weeks nominal. Honest estimate: 5–6 months single lane, 3–4 months in two lanes (T3/T2/T4 alongside T1/T5–T9). T5 and T8 are mathematics, not code, and carry most of the uncertainty.

**Risks.**
- R1 (highest). The truth may not meet the far-field thresholds at k ≈ 6600. The hand-off tables show the S1 chart ending at x ≈ 16 (k = 800) and 28 (k = 1600), growing like k^0.75. That extrapolates to about 80 at k = 6600, against the 350 required by H3/H4. STEP4_DESIGN §3 instead estimates the core at 0.25/|α| ≈ 500. If the first is right, k₁ must rise to several 1e4 (cheap for step 4: +8 units of τ; fine for step 3, since α₁ gets smaller), or the windows must shrink. T0 decides.
- R2. Thin certificate margins at k ≈ 10–30 (0.004–0.08). Mitigation: raise k_a.
- R3. Explicit constants for the far-field derivative levels (T8).
- R4. Wrapping over 400 slabs in P1. Mitigation: larger block, more bits; both cheap.
- R5. The step-3 thresholds for H3/H4 are still moving. Keep T9 parametric.
- R6. The float "no amplification" and "3 directions" facts rest on one clip and one trajectory. The H²-level count was never verified (FEAS caveat).
- The theorem is in the symmetric class throughout; odd perturbations are not covered and need not be.

**First milestone (about 1.5 weeks): T0 + T1.** An exactly symmetric float trajectory in the clock τ from k = 10 to k ≥ 6600 in the tip-adapted representation, with measured weak residual ≤ 1e-11 and strong residual ≤ 1e-8, and a table of every H1–H5 quantity against its threshold at candidate k₁. Go/no-go on k₁ and on the far-field thresholds before any rigorous code is written.

**Second milestone: one validated Fourier slab at t ≈ 2 (T3 prototype).** It proves the P1 technique and gives its real cost.
