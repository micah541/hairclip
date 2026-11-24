# Step 3: current argument (Nov 9, evening). This supersedes ASSEMBLY.md §7 (old Lemma O plan) and the "Plan B" entries of CERT_PLAN.md.

Tags: [C] certified (interval/ball arithmetic); [A] argued, still to be written or machine-checked; [F] float only; [O] open.

## 0. Claim
- **Hypotheses:** step 4 delivers the hand-off conditions of §8 at some time t1.
- **Conclusion:** the CDF solution exists on [t1, T). Its rescaled tip profile converges to S1 as τ → ∞ (modulo symmetries). The rest of the curve stays smooth up to T. So the singularity at T is modelled on S1.

## 1. Frame and decomposition
- **Frame.** Rescaled CDF in a modulated frame:
  - centre x_c(τ) on the symmetry axis;
  - scale λ(τ) ~ (T − t)^{1/4};
  - τ ~ −ln(T − t).
- **Chart.** Even class. The curve is a normal graph u over S1, with x = arclength from the tip.
- **Decomposition.** u = αψ̂3 + α²ĥ2 + w, with C⁴ versions of the slow-manifold functions (agent_int/v3/c4fun.py).
- **Modulation.** (a, b, β) are fixed by ⟨w̃, f_j⟩_ρ = 0 for j = 1, 2, 3. Then α_τ = λ̃α + β with λ̃ = −0.2561532.
- **Weight and spectrum.** ρ = exp(−(1+x²)^{2/3}/4). Even spectrum of L at S1: 1, 1/4, −0.2562, −0.6999, … [C, step 2].

## 2. Linear coercivity
- **Statement.** Re⟨Lv, v⟩_ρ ≤ −(3/64)‖v‖²_ρ for all even v ⊥ f1, f2, f3 [C].
- **Method** (writeup/step3_linear_certificate): Morse index of A = −H_sym − ω equal to 3, A invertible, Gram matrix negative definite.
- **Use in the energy.** The θ-form gives coercivity 1/9 with ω = 1/32 (the form used in the energy) [C].

## 3. Energies
- **Inner energy.** E = ‖w̃‖²_ρ + ε‖w̃″‖²_ρ, with w̃ = χ_32 w and ε = 1e-5.
- **Intermediate energy.** E_I = ‖v‖²_W + ε_I‖v″‖²_W, with v = ζw′ and ε_I = 3e-5.
  - ζ: 0 → 1 on x ∈ [10.25, 16.25] and 1 → 0 on [57, 60].
  - W = (1+x²)^{−2} for x ≤ 30.25, blended into a ρ-type weight beyond (ASSEMBLY.md §0).
  - Constants [C]: μ1 = 0.1552 (L-int), μ2 = 0.6555 (Gårding), c1W = 0.8477, K_W.
- **Lyapunov functional.** Λ = E + μE_I with μ = 0.08. The admissible window is μ ∈ [2.906e-3, 0.1135] [C, agent_coup].

## 4. Forcing, modulation, nonlinearity (all [C] unless marked)
- **Forcing.**
  - Inner: Fn(α) = 7.93e-7|α| + 2.769e-6α² + 12.104|α|³.
  - Far: FnW = 1.09e-7|α| + 6.5e-9α² + 0.1823|α|³.
- **Modulation.**
  - α-level form ≤ 6.19|α|, including the β-loop. E2-level form ≤ 43.33 (agent_mod).
  - |(a, b)| ≤ 3.708A linearly, and ≤ 7.57A with nonlinear feedback (jet box 1e-4).
  - |β| ≤ (0.00244 + 63.5|α|)A + 0.626‖χΦ_F‖.
- **Inner nonlinear terms.** KREG pointwise constants with cubic coefficients (cubic_E0fine), assuming a jet box of 1e-2.
- **Inner-ramp linear deviation** ("term (a)", [10.25, 16.25]): 0.19|α| + 578α² plus a w-part (agent_int/planB/planB_ramp.json). It needs S_ramp = max(sup|w|, sup|w′|) on [10.25, 16.25] ≤ 1e-2.
- **Far T3 coefficients** on [16.25, 60], for |w′| ≤ 0.02 (r60_T3_boxw1_0.02.json). These are exact polynomials in w″, w‴, w⁗, and the |w|-dependence is nil.
- **Other closure terms.**
  - (b) μ·9.9e-12·S_all·V0.
  - (c) E_I outer ramp: CFLUX·(levels). CFLUX is an INPUT, safe up to 3e-2, being certified (expected ≤ 0) [O].
  - (d) strip cubic terms with E_I sup bounds.

## 5. Closure (rigorous verifier, step3nl/verifier/)
- **Comparison system.**
  - d|α|/dτ ≤ λ̃|α| + |β|max(α, Λ);
  - (1/2)dΛ/dτ ≤ r(α, Λ)Λ, where r is the interval branch-and-bound maximum over shape variables (λ, θ, D4/B, θ_I, V4/V2).
- **Result at α1 = 1e-3, √Λ1 = 3.5e-7:**
  - peak √Λ = 6.39e-7, margin 11.7× to the cap 7.5e-6;
  - Λ → 0, with certified decay rate ≈ 0.011;
  - max S_ramp = 9.0e-3 (box 1e-2 OK);
  - E_I-based sup|w′| on x ≥ 16.25: 1.40e-2 (box 0.02 OK);
  - J0, J1, J2 = ∫_{τ1}^∞ V_j² dτ = 7.21e-11, 1.013e-10, 2.24e-11.
- **Backup** √Λ1 = 2.5e-7: margin 15.1×.
- **Caveats being fixed** [A]:
  - rounding is by 5e-16 inflation, not directed rounding;
  - the tail lemma (radial scaling, degree ≥ 2) is argued, not proved;
  - the meaning of the branch-and-bound leaf acceptance must be stated.

## 6. Far field
- **Key fact** [C, symbolic]: in shell variables ξ = y/y_k the slope obeys V_τ = −(ξ/4)V_ξ − y_k^{−4}∂²(aV″ + bV′²) + (S1 tail ~1e-14) + (modulation transport terms).
- **Lemma S(A): deviation shells, v = w′ on [24, 60]** (writeup/lemmaS):
  - inflow only through J_j, via K_j = 1.72e5, 8.56e7, 4.29e10;
  - the β-forcing accumulates only over transit time 3.7 (≈ 1e-5);
  - gives S0_v = sup|w′| ≤ 0.013 and S1_v = sup y|w″| ≤ 0.22 at √Λ1 = 3.5e-7;
  - the energy-inequality rate constants C_j are hand-derived, and Proposition 1 is being written and machine-checked [A].
- **Lemma S(B): full slope f = u′ on [48, R(τ)].**
  - R(τ) = R1·e^{(τ−τ1)/4} is a physically fixed edge inside the Lemma O region; it uses moving windows of fixed s-length.
  - No forcing. The modulation costs ∫|a|, ∫|b|.
  - Gives a priori S0_f, S1_f = O(1): parabolicity and graph property [A, in progress].
- **Lemma O, used ONCE at t1** for physical points with y(t1) ≥ y_O ≈ 150–200, with step-4 data (writeup/lemmaO/lemmaO_v3) [C].
  - Near the cone: X = 3 gives t* ≥ 3.70e-7 against a needed (3/y_O)⁴ ≤ 1.6e-7 at y_O = 150.
  - Far curve: M ≤ 3.69 at X = 2.
  - Hypotheses H1–H9 to be checked on the data.
  - It also gives: no other curve pieces in the balls, and smoothness of the far curve up to T.
- **Region bookkeeping.**
  - Inner energy: x ≲ 32.
  - E_I: [10.25, 60].
  - S(A): [24, 60]. S(B): [48, R(τ)].
  - Lemma O: physical r ≥ y_Oλ1 (for all t ≥ t1).
  - The S(B) edge lies inside the Lemma O region.

## 7. Bootstrap (continuity argument) [A, not written]
- **Assume** on [τ1, τ*]: the jet boxes (inner ramp 1e-2, far |w′| ≤ 0.02, inner cubic boxes) and the parabolicity/graph bounds hold.
- **Derive** the closure bounds (§5), Lemma S, and Lemma O. Then show the boxes hold with margin at τ*.
- **Conclude** τ* = ∞.
- **Also needed:**
  - local well-posedness/continuation for CDF (cite);
  - modulation by the implicit function theorem all along (|α|, A small);
  - the curve stays a graph over S1 for x ≤ 60.

## 8. Hand-off requirements (to be supplied by step 4)
- (a) (x1, λ1, α1) given by the implicit function theorem (constraints satisfied).
- (b) |α1| ≤ 1e-3.
- (c) Λ1 = E + μE_I at τ1 with √Λ1 ≤ 3.5e-7.
- (d) Lemma S inputs: shell energies at τ1 on [24, ∞) (sliding windows), and far slope ≤ 0.5.
- (e) Graph over S1 out to the S(B)/Lemma O overlap.
- (f) Lemma O H1–H9 for y ≥ y_O, including the L² bound on h‴ and the distance to the other arm.
- (g) The ‖χΦ_F‖ bound (forcing norm) at the hand-off α.

## 9. Conclusion step [A, not written]
- Λ → 0 and α → 0, so the rescaled tip profile converges to S1 and λ(τ) → 0 (finite-time singularity at T).
- The far curve stays smooth (Lemmas O and S), so the singular set is the single tip point.
- The blow-up rate is self-similar.
