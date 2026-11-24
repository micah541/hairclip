# Step 3 assembly blueprint (Nov 7, evening)

This file lists the full closed system of differential inequalities for the trapping argument, with the source of every constant. The float implementation is `inner_closure_float2.py`; a rigorous verifier must reproduce it with explicit splits.

Status tags:
- **[C]** certified (interval arithmetic);
- **[F]** float or modelled, to be certified;
- **[A]** analytic step still to be written out.

## 0. Variables and norms
- **Coordinates.** x = σ − σ_tip, on the half line (even class).
- **Weights.**
  - ρ = e^{−φ}, φ = (1+x²)^{2/3}/4.
  - W = e^{−ψ}, ψ′ = (1−s)·4x/(1+x²) + s·φ′, s = S((x−30.25)/8). So W = (1+x²)^{−2} for x ≤ 30.25.
- **Decomposition.** u = αψ̂3 + α²ĥ2 + w.
  - w̃ = χw, χ = χ_32 (C³ smoothstep on σ ∈ [29, 32]).
  - v = ζw′, ζ = 0 → 1 on x ∈ [10.25, 16.25] and 1 → 0 on [57, 60].
- **Shape variables (inner).** A = ‖w̃‖_ρ, W1 = ‖w̃′‖_ρ, B = ‖w̃″‖_ρ, W3 = ‖w̃‴‖_ρ, D4 = ‖w̃⁗‖_ρ.
  - Interpolation **[C]**: W1² ≤ AB + c1A², W3² ≤ B·D4 + c1B², with c1 = 0.5804.
- **Shape variables (intermediate).** V_j = ‖v^{(j)}‖_W.
  - Interpolation **[C]**: V1² ≤ V0V2 + c1W·V0², with c1W = 0.8477 (R = 60).
- **Energies.**
  - E = A² + εB², with ε = 1e−5.
  - E_I = V0² + ε_I·V2², with ε_I = 3e−5.
  - Λ = E + μE_I, with μ = 0.08 (admissible window [2.906e−3, 0.1135] **[C]**).
- **Pointwise bounds** **[C]**: sup_reg g² ≤ K_reg(Φ1‖g‖² + 2‖g‖‖g′‖).
  - ρ-regions [0,1], [1,2], [2,4], [4,8], [8,12], [12,16.25]: K = 1.488, 2.078, 5.222, 56.92, 992.4, 3.0208e4; Φ1 = 1.083.
  - W-regions [16.25,20], [20,30], [30,34.25]: K_W = 1.609e5, 8.119e5, 2.323e6; Ψ1 = 1.305.
  - sup_{[16.25,34.25]} |w| ≤ |w(16.25)| + Φ(34.25)^{1/2}·V0, with Φ(34.25) = 9.8547e6.

## 1. Modulation [A + C-parts]
The constraints ⟨w̃, f_j⟩_ρ = 0 (j = 1, 2, 3) give, for x = (a, b, β):

  Σ_i x_i (G_ij + E_ij(α, w)) = ⟨w̃, L*f_j⟩ − ⟨[L,χ]w, f_j⟩ + ⟨χΦ_F, f_j⟩ + ⟨χ𝒩(w), f_j⟩,

with E_1j = −⟨χA1(u), f_j⟩, E_2j = −⟨χB1(u), f_j⟩, E_3j = 2α⟨χĥ2, f_j⟩.

Used bounds:
- |(a, b)| ≤ 31.2 A (Minv 2.0005, ‖(L*f1, L*f2)‖ 15.59) **[C]**, plus corrections from E, Φ_F, 𝒩 and the strip **[F]**.
- |β| ≤ (ε1 A + 0.63 ‖χΦ_F‖ + ⟨χ𝒩, f3⟩-part)/G33, with ε1 = 6.645e−5, G33 ≥ 0.031105 **[C]**.
- TODO **[A]**: a rigorous perturbation bound for (G + E)^{−1}, and bounds on the 𝒩- and commutator contributions to (a, b, β).

## 2. Slow mode
- α_τ = λ̃α + β, with λ̃ = −0.2561532.
- Comparison: |α|′ ≤ λ̃|α| + |β|.

## 3. Inner E0 level
(1/2) d/dτ A² ≤
- −(1/32 − Π3) A² − (1/9) B², where Π3 = ε1·n3/G33 = 0.01607 and n3 = ‖χψ̂3‖ = 7.52 **[C]**;
- + |α| Q_α(A, W1, B), the α-linear form. The regional spectral version is used: 59.3(A² + W1²) + 35.9AB + 7.9W1B + 6.2B² **[C]**;
- + R1|α|A², the exact modulation α-form, R1 = 13.549 **[C]**;
- + (1 + n3‖f3‖/G33)·‖χΦ_F‖·A = 5.72‖χΦ_F‖A, with ‖χΦ_F‖ ≤ 1.3e−6|α| + 1.6e−6α² + 12.34|α|³ **[C]**;
- + Σ_reg (1/2) Σ_mono C_mono,reg·(two highest-derivative factors in L²)·(the others in sup_reg) **[C coefficients, cubic_E0fine]**.
  - In the region [16.25, 31.25], the sups come from E_I (§0);
- + strip ≤ sB·V0² + sA·V0² + 2·IA·K_ρ(16.25)(A² + 2A·W1) **[C]**: sB 3.5955e−5, sA 1.3372e−4, IA 9.60e−12.
  - This comes from −⟨[L,χ]w, χw⟩ = ∫A w² + B w′². The drift part has a good sign and is kept inside A.
- TODO **[A]**: the β-loop contribution of ⟨χ𝒩, f3⟩ (cubic size), and the (a, b)-corrections times α‖(gA, gB)‖.

## 4. Inner E2 level (weight ε)
(ε/2) d/dτ B² ≤ ε × [
- −(1/2) D4² + 14.34 B² + 16.3 W1·B + 10.7 A·B (Gårding) **[C]**;
- + ‖χD‖·(D4 + 2.2 W3 + 1.2 B), where D = F(u0+w) − F(u0) − Lw:
  - tip coefficients 6.9e−3 D4, 1.2e−2 W3, 4.2e−2 B, 0.54 A, 0.62 W1, and 4.75, 0.136, 3.05 × sup|w″| terms **[C]**;
  - far [12, 34.25] coefficients 4.3e−3, 6.9e−4, 2.9e−5, 0.014, 0.46, 3.06; A0/A1 8.6e−6 (q/4 fix) **[C]**;
- + 43.33 A², the exact E2 modulation form **[C]**;
- + ‖χΦ_F‖(D4 + 2.2 W3 + 1.2 B);
- + strip, level 2: (d3 V2 + d2 V1 + d1 V0 + d0 ‖w‖_{W,strip})·‖D_ρ w̃‖, with d = 6.84e−3, 1.54e−2, 2.36e−2, 4.56e−2 **[C]**. ]

## 5. Intermediate level 0 [derivation written out Nov 7; constants C]
v_τ = L1 v − [L1, ζ]g + ζc0′w + ζΦ_F′ + ζ𝒩(w)′ + ζM′, where g = w′ and M = −βψ̂3 − 2αβĥ2 − aψ1 − bψ2 + aA1(u) + bB1(u).

(1/2) d/dτ V0² ≤
- −μ1 V0² − (1/2) V2², with μ1 = 0.1552 (0.3357 on the blend, which contains the strip) **[C, L-int]**;
- + a_f W1² + b_f B², with a_f 6.956e−2, b_f 3.936e−2 (inner ramp) **[C]**;
- + outer ramp ≤ 3(sup|A_zo|·G1² + sup|B_zo|·G2²), with sup|A_zo| 3.2e−18 (G = crude bounds from Lemma O) **[C]**;
- + ‖ζc0′‖_W·sup|w|·V0, where ‖ζc0′‖_W = 9.9e−12 **[C]**;
- + ‖Φ_F‖_W·(V1 + (sup ζ′ + Ψ1) V0), with ‖Φ_F‖_W ≤ 1.09e−7|α| + 6.5e−9α² + 0.1823|α|³ (after one integration by parts) **[C]**;
- + 0.0185 V2² + 0.025 V1V2 (top-order 1 − J⁻⁴) + ‖ζ D_low‖_W (V1 + 1.67 V0).
  - D_low has d_010 1.5e−3 (w‴), d_100 1.1e−4 (w″), d_200 0.015 (w″²), d_110 0.98 (w″w‴), d_300 3.26 (w″³), A0 1.7e−5, A1 7.1e−5 **[C, agent_int/v2]**;
- + modulation: |β|·0.26 V0 + 2|α||β|·1.8 V0 + |a|·1.7e−7 V0 + |b|·1.7e−8 V0 + |a||α|·0.273 V0 + |b||α|·0.016 V0 **[C norms]**;
  - plus |a|·[(1/2) sup|1 − qψ′|·V0² + sup_ramp(q ζ′ W/ρ)·W1²] + |b|·0.66 V0², where sup|1 − qψ′| ≤ 77 on [10.25, 60] **[F constant]**;
- + inner-ramp couplings of 𝒩 to ρ-norms of w̃ (coefficients ≤ T3 [10.25, 16.25] × (W/ρ)^{1/2} ≤ 0.66) **[A, small]**.

## 6. Intermediate level 2 (weight ε_I)
(ε_I/2) d/dτ V2² ≤ ε_I × [
- −(1/2) V4² − 0.6555 V2² + 5.7e−6(V3V2 + V2²) (Gårding) **[C]**;
- + (c3 D4 + c2 W3 + c1 B + c0 W1)·‖D2 v‖_W, with c = 0.442, 0.526, 0.424, 0.425 **[C]**, and ‖D2 v‖ ≤ V4 + 2Ψ1 V3 + 2c1W V2;
- + 0.0185 V4² + ‖ζ𝒩′_low‖·‖D2 v‖. The latter uses sup|w″| and sup|w‴| on the far region **[C coefficients]**;
- + ‖ζΦ_F′‖_W·‖D2 v‖ **[F: derivative-level forcing not yet certified; needs fifth jets of ψ̂3, ĥ2]**;
- + modulation: |a|(38.5 + 0.66) V2² + 2|a| V2V3 + |a||α|·‖ζA1(u0)‴‖·‖D2 v‖ **[F constants]**. ]

## 7. Outer region (Lemma O, to be written) [A]
Statement (scale-invariant short-time stability for CDF):
- **Hypothesis.** At t0, γ ∩ B(p, ρ) is one graph arc over a line, with |h′| ≤ 1/8 and |h^{(j)}| ≤ Mρ^{1−j} for 2 ≤ j ≤ 6.
- **Conclusion.** On [t0, t0 + δρ⁴], with δ ≤ δ0(M), the curve stays a graph on B(p, ρ/2) with constants 2M, and moves by at most C(M)δρ.
- **Use.**
  - Every physical point at distance r from the centre is handed over at y_h = 30, i.e. at r = 30λ(t_h). The remaining-time ratio is (T − t_h)/(r/4)⁴ ≈ 3.2e−4.
  - This gives:
    1. the crude bounds G1, G2 on [57, 60] for §5;
    2. the far-field nonlinear coefficients of E_I on [30, 60];
    3. the absence of other curve pieces in B(x, 60λ);
    4. the chart for y ≤ 60.
- The far curve outside 60λ1 at hand-off is handled the same way, with data from step 4.
- **Proof plan.** Use the graph PDE for the slope f = h_x: f_t = −(a(f)f_xx + b(f)f_x²)_xx, with a = (1+f²)^{−2} ≥ 0.886 and |b| ≤ 0.63 for |f| ≤ 1/4.
  - Hand-over regularity: E_I controls only (w′, w″, w‴) pointwise in time, but its dissipation gives ∫‖v⁗‖²dτ < ∞. So pick a *good* hand-over time per point in the window y_h ∈ [30, 34] (Chebyshev), where ‖w^{(5)}‖_{L²} is bounded. Then h0 ∈ H⁶ locally, with h0^{(j)} pointwise for j ≤ 5.
  - A priori bounds (continuity argument, variable-scale weights ρ(q) = |q − x_c|/4): |f| ≤ 1/4 and ρ^j|f^{(j)}| ≤ 2M for j ≤ 3.
  - Energies E_k = ∫η_k²(f^{(k)})² for k = 4, 5, with nested cutoffs: dE_k/dt ≤ −1.77∫η²(f^{(k+2)})² + C(K)(E_k + 1), so E_k(t) ≤ (E_k(0) + Ct)e^{Ct}.
  - Improvement step:
    - ‖f‴(t) − f0‴‖_{L²} ≤ √t·(∫‖∂⁵G‖²)^{1/2}, using the k = 5 dissipation;
    - ‖(f‴(t) − f0‴)′‖_{L²} ≤ √E_4 + M_4;
    - hence sup|f‴(t) − f0‴| ≲ t^{1/4}(CM_4)^{1/2}; lower derivatives get smaller powers of t;
    - with t/ρ⁴ ≤ 3.2e−4, t^{1/4} = 0.13.
  - Displacement: |h(t) − h0| ≤ ∫|G_x| ≤ C(K)t. Foreign pieces are excluded globally, since every point moves at most Cδρ(p).

## 8. Closure
Comparison system for (α, Λ) on τ ≥ τ1:
- d|α|/dτ ≤ λ̃|α| + |β|max(α, Λ);
- (1/2) dΛ/dτ ≤ r(α, Λ)·Λ, where r is the rigorous maximum over shapes of the right-hand sides of §3 + §4 + μ(§5 + §6) at fixed Λ.

Float results (all [C] constants):
- α1 = 1e−3: √Λ1 = 1.3e−6 → max 1.92e−6;
- α1 = 2e−3: max 4.43e−6;
- cap ≈ 7.5e−6 to 1e−5.

Rigorous verifier, TODO [A]:
- explicit Young splits, or interval branch-and-bound over the shapes (λ, θ, D4/B, θ_I, V4/V2);
- a monotone upper envelope of r on an (|α|, log Λ) grid;
- interval integration of the comparison ODE.

## 9. Hand-off requirements (interface with step 4)
At τ1, all of the following must hold:
- (a) the constraints hold, i.e. (x1, λ1, α1) are chosen by the implicit function theorem;
- (b) |α1| ≤ 1e−3 (or 2e−3);
- (c) ‖w̃1‖_ρ² + ε‖w̃1″‖² ≤ (1e−6)² (scaled with α1);
- (d) ‖v1‖_W² + ε_I‖v1″‖² ≤ (3e−6)²;
- (e) a graph over S1 for y ≤ 60 with |u′| ≤ 0.1;
- (f) the Lemma O hypotheses for all curve points at distance ≥ 30λ1, with M ≤ M0.

Numerics (agent_handoff2, interim) give dist ≈ 170|α|³, i.e. 1.7e−7 at α = 1e−3, reached at k_max ≈ 3300.
