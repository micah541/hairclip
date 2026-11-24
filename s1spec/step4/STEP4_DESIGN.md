# Step 4 design note (Nov 9): validated integration from the clip into the step-3 trapping set

## 1. Target: the hand-off conditions (ASSEMBLY.md §9 with Plan B)
At some time t1 (k_max ≈ 2900 for the clip with arm 24, chamber 0.4, W 0.85, at |α1| = 1e-3):
- (a) modulation (x1, λ1, α1), obtained from the implicit function theorem;
- (b) |α1| ≤ 1e-3;
- (c) ‖w̃1‖²_ρ + ε‖w̃1″‖²_ρ ≤ (1e-6)²;
- (d) ‖v1‖²_W + ε_I‖(1 + (x/60)⁴)^{1/2} v1″‖²_W ≤ (3e-6)² on [10.25, 260];
- (e) the curve is a graph over S1 for y ≤ 260, with the slope bounds;
- (f) the Lemma O hypotheses hold beyond 200λ1 and on the far curve.

Error budget. The true remainder there is about 45|α|³ ≈ 4.5e-8. So the validated error must be ≲ 1e-6 in the ρ-norm at the tip. It may be much looser further out: W(260)/W(10) ≈ 9e-9, and the far field needs only O(1) derivative bounds.

## 2. Why step 4 cannot be short-cut
S1 is not embedded: its loop has a transversal self-crossing at y = O(1). Condition (c) controls w in C⁰ near the crossing (H² ⊂ C¹ locally). A C⁰-small normal perturbation of two transversally crossing strands still crosses. So every curve in the trapping set self-intersects. An embedded initial curve must therefore be carried through the crossing by the flow itself, with validated error bounds.

"S1 with its crossing resolved" is exactly the slot-plus-chamber design. The shape study shows that the long arms are also needed: a disk with a slot, or arms ≤ 5, do not singularize.

## 3. Structure of the trajectory (float runs)
- **Time span.** In s = −ln(T − t), the run covers Δs = 4 ln(k*/k_start) ≈ 20–25.
- **Spectrum at S1** (s-units): 1 (T-shift), 1/4 (axial translation), −0.256 (α), −0.70, …
- **Shape of the solution.**
  - From k ≈ 100 on, α·k_max ≈ constant (2.9–3.3). So the S1-like core extends to y ≈ 0.25/|α|, i.e. physical r ≈ 0.13.
  - Beyond that the curve is quasi-static: the arms and the far bend move very little during the remaining time.
- **Errors near S1.** Modulo (T, translation, α), errors contract like e^{−0.7 Δs} spectrally. Step 3's energies already give a contraction rate of about 0.1 near S1.
- **Uniform-arclength cost.** HLS needs N ≈ 80 k_max (2.6e5 at k*). In float a run to k = 1600 takes about 200 s on this 2-core box. This is fine for float, but not for validation: in a uniform Fourier basis the tip's slow dynamics is spread over about 1e5 modes. So uniform Fourier with self-consistent bounds (CAPD/Zgliczyński style) is ruled out.

## 4. Proposed method (to be confirmed by §5)
1. **Approximate solution ũ.**
   - Double precision (double-double near the tip if needed).
   - Half curve (symmetry), in a modulated rescaled frame: centre on the axis, scale λ(t) fixed by a normalization.
   - Space: spectral elements (Chebyshev) on a geometric partition of the rescaled arclength σ, with element size ∝ 1 + |σ|.
   - Time: piecewise polynomial in s on each step, so the residual R = ∂_s ũ − F(ũ) is an explicit function that can be enclosed with interval arithmetic.
2. **Error equation** in the normal-graph chart over ũ: e_s = L_ũ e + τ_ũ e_σ + R + N_ũ(e) + modulation terms.
3. **Error propagation.**
   - *Core* (ρ-weighted, y ≲ 15–30): split e = P_M e + Q_M e, with M low modes in a fixed ρ-orthonormal basis.
     - Low modes: rigorous linear propagator (numerical fundamental matrix plus an a-posteriori enclosure, i.e. Lohner/QR). This captures contraction and expansion exactly.
     - Tail: certified coercivity ⟨sym(L_ũ) q, q⟩_ρ ≤ −c‖q‖²_ρ on ran Q_M (Gårding plus interval bounds on ũ's coefficients), as in step 3 but with ũ in place of S1.
   - *Intermediate* (scale-adapted weights as in E_I, out to where the curve is quasi-static): energy estimate whose coercivity comes from the weight design, so no spectral information is needed.
   - *Far field* (arms, far bend): crude local energies; only modest accuracy is needed.
   - *Symmetry directions* (T, translation, α): absorbed by modulation. Their errors are allowed to grow; only the hand-off implicit-function step needs them bounded.
4. **Precision.** Final error ≈ Σ_j ‖Φ(s_end, s_j)‖·‖R_j‖. If the amplification modulo symmetries is ≤ 1e4, a residual floor of about 1e-12 suffices.

## 5. Feasibility study (float only; run before any rigorous code)
- **F1 (resolution).** Find the DOF of the geometric spectral-element representation needed for a residual ≲ 1e-12 at k = 100, 400, 1600, and extrapolate to 2900. Expected: a few hundred to about 1e3, against 2.6e5 uniform.
- **F2 (sensitivity).** Measure the amplification modulo (T, translation, α) of perturbations introduced at k = 10, 30, 100, 400, observed at k = 1600, in ‖·‖_ρ (core) and in the W-norm.
- **F3 (tail dimension).** Find the smallest M such that sym(L_ũ) ≤ −c (c = 0.05, 0.1) on the ρ-orthogonal complement of the first M basis functions, at k = 10, 30, 100, 400, 1600. Check: for S1 itself the count of non-coercive directions must come out small (3 plus a few).
- **Decision rule.**
  - If F2 ≲ 1e4 and F3 ≲ 60 throughout: build the hybrid (Lohner core, coercive tail, weighted energies).
  - If F3 is large only early: give the early phase (short Δs) a larger Lohner block.
  - If F2 is huge: revisit the shape and choose the variant with the smallest early amplification.

## 6. Feasibility results (Nov 9, feas/FEAS_REPORT.md) and revised plan
- **F1.** Geometric Chebyshev representation, θ coefficients ≤ 1e-13: 462 / 351 / 295 DOF at k = 1600 for r = 0.25 / 0.5 / 1, and 484 / 368 / 305 at k = 2900. Growth is ~ln k, against 2.6e5 uniform.
  - θ⁗ is not resolved to 1e-12 at that truncation (errors 2.4e-6 to 9e-5); roughly +160–410 DOF are needed for the residual.
  - The snapshot floor is ~1e-13 to 2e-12.
- **F2.** No amplification modulo (T, translation, α):
  - ≤ 0.33 at the perturbation, ≤ 7e-3 by k = 1600 for k_j ≥ 100;
  - linear in δ; Δα/δ ≤ 0.25;
  - for k_j = 10 and 30 the response falls to a noise floor of ~1e-10.
- **F3.** The linearization was validated by finite differences and against the S1 spectrum (1, 0.25, −0.25615, −0.69989). The number of eigenvalues of sym_ρ(L_Γ + transport) above −0.1 is exactly 3 at EVERY k from 10 to 1600, the same as for S1. The complement of the first 3 S1 eigenvectors already works (M = 3).
- **Consequence: no Lohner block is needed.** Step 4 becomes a shadowing argument of step-3 type around the numerical trajectory ũ(τ):
  - the same energies;
  - 3 modulation parameters (T-shift, axial translation, α-shift);
  - coercivity ≳ 0.1 certified at a grid of times plus perturbation bounds;
  - forcing = the residual of ũ.

  Steps 3 and 4 can be glued: base ũ(τ) for τ ≤ τ1, base u0(α(τ)) after. The hand-off then becomes continuity of the energy bounds, not a separate (a)–(f) check.
- **Numerical base.** HLS float runs (autonomous in θ, L, so the time bookkeeping does not limit precision), converted to tip-adapted Chebyshev in rescaled σ and Chebyshev in s per step; the residual is enclosed with interval arithmetic.
- **Shared critical path with step 3.** The far-field framework (E_I to 60, Lemma S shells, Lemma O once at the end) is used by both steps.
- **Caveats.** One clip, one trajectory, float only. The H²-level (ε) variant of F3 was not verified. The early phase k < 10 is still to be checked.
- **(Nov 9) Macroscopic phase found.** In the target run (design/run3_clip_24_0.4_0.05_0.85.log) k_max = 42 at t = 0. It smooths to about 2 by t ≈ 0.45, and the curve SELF-CROSSES by then (ncross = 4; 2 from t ≈ 3 on, the loop). T ≈ 4.5988639.
  - So step 4 must also cover a long macroscopic phase, t ∈ [0, ~4.59] with k_max 2–42, before the rapid phase studied in F2/F3 (k ≥ 10).
  - Feasibility part 2 (feas2/) measures, for this phase: sensitivity, the tail count, the cost of a Galerkin+tail coercivity certificate per profile, and the spacing between certificates.

## 7. Feasibility part 2 (Nov 10, feas2/FEAS2_REPORT.md)
- **Crossings.** 0 → 2 at t ≈ 0.003, 2 → 4 at t ≈ 0.016, back to 2 at t ≈ 2.64. Two persist; T ≈ 4.5989.
- **F2b.** No amplification from the macroscopic phase. Perturbations injected at t ∈ [0, 4.5] shrink about 100× between k = 100 and 1600 (worst response 0.40 per unit δ).
- **F3b.** Non-coercive count above −0.05 in the actual rescaled frame: 3–5 in the macroscopic phase (4 typically, 5 near t = 3), and 3 at t = 4.59. A frame-consistent count gives 10–16 for t ≤ 3.5.
- **P1.** The certificate closes at Galerkin size N = 3 (k = 1600), 6 (k = 100), 29 (k = 10), 43–47 (t = 4.59), 54–57 (t = 4.5) and 76 (t = 0.5, 2.0). Margins are 0.001–0.29.
- **P2.** About 42 certificates for t ∈ [0.01, 4.55], 2 for the approach, 3 for the rapid phase out to k = 6600. The initial layer t < 0.01 is not covered.
- **Floors.**
  - Beyond k ≈ 2500, double-precision t no longer advances, so step 4 must use T − t or s as the clock.
  - T varies by about 1e-5 with the step sequence.
- **Verdict.**
  - The 3-direction shadowing holds for t ≥ 4.55 (k ≥ 5).
  - The macroscopic phase needs 4–5 modulated or low-mode directions in the rescaled frame, or a different frame or norm.
  - Idea: choose the initial curve as the smoothed clip at t ≈ 0.002, still embedded, given by an exact finite series, to skip the stiff initial layer (k_max = 42 at t = 0).

## 8. Architecture (Fable, Nov 10): STEP4_ARCH.md
- **Engine 1, macroscopic and approach phases** (t = 0 to k ≈ 10–14):
  - validated Fourier integrator in HLS variables (semilinear, cubic, L-free in HLS time), physical frame, no modulation;
  - Lohner block of ~100–200 low modes plus an analytic tail, on Chebyshev-in-time slabs, in arb;
  - float check: 150–450 modes for 1e-13; 3–10 frozen unstable directions; growth ~1–2 per unit t (worst transverse ≲ 1e3).
- **Engine 2, rapid phase** (k ≈ 10–14 to hand-off):
  - step-3-type energy shadowing around a tip-adapted spectral-element trajectory in the rescaled frame;
  - T-shift and translation modulated exactly; the third direction carried as a scalar;
  - about 5–8 certificates; clock τ with ln λ as an unknown.
- **Initial curve:** a smoothed clip as a finite sine series. Closure by interval Newton; embeddedness certified in arb. k_max = 2.8 already at t = 0.002.
- **Hand-off:** base switch by the triangle inequality, so step 3 is not reopened.
- **Error budget:**
  - tip: allowed ≈ 2e-8 in √Λ; the rapid-phase rigorous contraction is only ≈ 0.06–0.1, so the ball entering it must be ≤ 1e-9;
  - double precision is OK at the tip with a weak-form residual;
  - elsewhere use iterative refinement against arb residuals at 128 bits.
- **Size:** 11 tasks (3 large, 6 medium, 2 small); about 15–16 weeks nominal (5–6 months in one lane, 3–4 in two).
- **Main risk:** the far-field truth at the hand-off k. The chart-end extrapolation is ~80 (reviewer) or ~450–500 (chart_end ∝ k, from the uncapped 27.6 at k = 400) against the 350 needed.
- **Milestone 1** (go/no-go, float): an exactly symmetric trajectory from k = 10 to k ≥ 6600 in clock τ, with every H1–H5 quantity tabulated against its threshold.

## 9. Milestone 1 (float go/no-go, m1/STEP4_M1_REPORT.md): NO-GO at k <= 1.01e4
- **Pass:**
  - tip: alpha <= 5e-4 from k ~ 5.6e3; sqrt Lyap <= 5e-8 from k ~ 4.7e3 (2.2x at 6081, 6.5x at 1.01e4);
  - inner shells [14,60] from k ~ 2020;
  - S1 chart extends to ~0.117 k (past 380 for k >~ 3250);
  - Lemma O cone/far jets from k ~ 5.2e3 / 4.2e3.
- **Fail:** outer windows at levels 3-4. At k = 1.01e4: E3[169,253] = 5.1e-6 (1.7x over), E4[112,169] = 6.2e-6 (2.1x), E4[169,253] = 4.4e-5 (15x). At k = 6600: 21x and 183x over.
- **Cause:** the deviation from U is a clean cubic, v ~ -c3 (y/k)^3 with c3 ~ 38-40.6 (missing alpha^3 term in U), so the window energies scale like k^-5.9.
- **Extrapolation:** all items hold from k ~ 1.6e4 (2.2-2.3e4 including the sliding windows to 380).
- **Floor:** N = 2^19 insufficient at k >= 8000, so 2^20 was used.
- **DECISION: raise the hand-off to k ~ 2.3e4** (alpha1 ~ 1.25e-4).
  - Cost: about +1.8 doublings of k, i.e. ~30-55 more step-4 slabs.
  - Step 3 should only get easier (smaller alpha means smaller forcing, which is ~69% of J0).
  - Verifier and Lemma S reruns follow.

## 10. T3 prototype (GO) and T5 design (Fable), Nov 10-11
- **T3 (t3/T3_REPORT.md): validated Fourier slab works.**
  - Method: Duhamel bounds, Lohner block, n^-6 tail, Chebyshev base.
  - t ~ 2: 30 s per slab; 5 chained slabs are stable.
  - Extrapolation: ~80-91 slabs to k = 10-14 (to k = 20: 5-8 CPU-h).
  - Final C^2 error ~3e-14.
  - **Conditions:**
    - consistent AND embedded Gamma0. The initial layer from the corner clip is present at t = 0.001. The crossing at t ~ 0.0028 rules out long relaxation; use graded slabs or a later embedded start.
    - a tight b2 interval;
    - block resizing.
  - The agent was resumed for these, plus the approach test, the X0/k_max quadratures and the Galerkin-limit lemma.
- **T5 design (t5/T5_DESIGN.md): the P2 shadowing theorem.**
  - **Headline risk:** the far field at levels 3-4. Hand-off item 3 (E3 windows) and item 4 (K4) require a tip error Ybar <= 2.5e-11, so the target is 1e-11.
    - Estimated amplifications: 8e7 at level 3 and 3e9 at level 4.
    - This forces weak residuals of 2e-14, strong of 5e-11, an arb-polished C^4 base, and a new linear 5-level far-field Lemma F (T8).
    - Fallbacks: time-differentiated energies (10-60x); an interior smoothing lemma; asking step 3 for less (level 3 from levels <= 2).
  - **FEAS2 certificates certify the wrong form.** The modulation term breaks coercivity on the complement of the S1-sym eigenvectors.
    - Fix: constraints {f1, f2, chi psi3}, with the base's own modes as directions.
    - Coercivity is 0.061 (k=10) to 0.083 (k >= 40). Expect N ~ 40-80 at small k.
  - **Third direction:** handled by clock modulation, since the tau-shift along the base orbit is an exact symmetry. The 2x2 system is triangular, so the closure condition is empty.
  - **Recommended start:** k_a ~ 20.
  - **P2 cost:** ~190 slabs, ~45-50 certificates, 10-30 CPU-h. The human cost is in T8 and T7.
- **Next:** T1 (native rescaled-frame base solver, C^4, arb residual feasibility, hand-off check at 2e4) launched; T3 continuing.

## 11. T3 complete: GO for P1 production (Nov 10)
- **Gamma0 final** (t2/gamma0_final.json):
  - modes 186-1000 slaved to the slow manifold, changing coefficients by <= 1.6e-10;
  - b2 re-closed at 128 bits (half-width 8e-24);
  - embedded at t = 0.001 (slot half-gap >= 0.0248199).
- **Graded initial slabs work:** residual 6e-12, then 8.6e-14; error set ~1e-20.
- **Block resizing works across the swing:** 25 slabs, growth 1.06-1.15 per slab (true linearised growth).
- **Approach test at k ~ 13-16:**
  - needs M = 1300 at k = 12, 1500 at 14, 2000 at 16;
  - 24 time coefficients;
  - 720-800 s per slab.
- **Quadratures:** t, L, X0 and k_max done. **Galerkin-limit lemma** written (T3_DESIGN §6).
- **Production estimate:** 103 slabs to k ~ 12, 3.8 CPU-h (~6 h with allowances).
- **Open:** worst-case final C^2 error is 3e-13, against T5's P1 budget of C^2 <= 3e-14 and C^5 <= 1e-9 (P2 target Ybar <= 1e-11). Asked T3 to localise the growth, find the cheapest fix and write a resumable driver run_p1.py. A cheap runner agent will then run the chain.

## 12. T1 complete (Nov 10): native P2 base, k = 20 -> 2.36e4 (t1/T1_REPORT.md)
- **Representation:** degree-11 Hermite-C^4 elements, exactly symmetric, kappa(0) = K0, ln lambda unknown. Radau IIA (s = 8, dtau = 0.15).
- **Run:** 188 slabs, 840 -> 1218 DOF, 950 s float; arb-128 polish of the whole base 40 min.
- **Deviations from T5, each forced by a measured failure:**
  - hybrid Legendre/Galerkin test space (pure moments have condition 8e12);
  - Lagrangian outer knots (a lagging mesh is unstable).
- **Signs:** T5's signs are confirmed symbolically.
- **Residuals** after the full polish, k >~ 60: R2 1.5e-17, R1 4e-16, Rs 3e-14, RF3 1.2e-9, RF4 3.3e-8. All T6 targets are met with margin >= 80.
  - The first 4 slabs miss, so start the base about 0.6 in tau before tau_a.
  - The whole base must be polished from the seed.
- **Hand-off at k = 2e4 / 2.3e4:**
  - |alpha| 1.345e-4 / 1.166e-4 (x1.12 / x1.29) — binding;
  - sqrt(Lambda) 6.5e-10 / 4.3e-10;
  - window E0 x8 / x14, E3 x15 / x22;
  - K4 x250 / x290;
  - shells >= 1400x under.
  - **k1 = 2.3e4 recommended.**
- **M1 cross-check:** the base reproduces M1 to ~1e-3. M1's "excess A" at k = 1e4 was an HLS floor.
- **Next:** T8 (far-field Lemma F, levels 0-4; the headline risk) using the T1 base; then T7 (coercivity certificates, projected form) and T6 (rigorous residual enclosures, mostly built on T1's arbops).

## 13. P1 error budget met; production started (Nov 10)
- **Cause of the earlier 10x overshoot:** a 5.5e-18 mismatch at mode 250 between Gamma0 and slab 1's start. It was absorbed into the edge of the slow block (x1560 jump in slab 1).
- **Fix:** a second slaving pass on Gamma0 (modes >= 250), then b2 re-closed and embeddedness re-certified (t2/gamma0_final.json).
- **Expected error at the switch** (float estimate x3):
  - C^2 2.6e-16 at k = 12, 3.3e-15 at k = 20; budget 3e-14.
  - C^5 8e-14 and 4.8e-12; budget 1e-9.
- **Optional:** about 20x more margin from a tighter closure quadrature.
- **Production:** t3/run_p1.py, resumable, one checkpoint per slab. Started with --target-k 20 (expected 5.1 CPU-h, x1.5 for retries), with a cheap runner agent babysitting it. The validated C^2/C^5 are reported per slab in t3/p1/status.json.

## 14. T8 complete (Nov 10): Lemma F proved (linear), far-field accuracy requirement (t8/T8_REPORT.md)
- **Amplification at tau1** (T5 design, 188 slabs): A0..A4 = (3.4e3, 5.7e4, 5.6e5, 1.8e8, 8.2e10). T5's guess was 2.2x optimistic at level 3 and 27x at level 4.
- **Requirement on Ybar:**
  - 8e-12 as designed (K4 via F4 binds);
  - 1.35e-11 with a "good-time" hand-off (level 4 drops out);
  - **2.0e-11 recommended:** level-3 ramp [14,22] and level 4 from level-2 dissipation. Item 3 binds; the rest are >= 2.9x clear.
  - Expected budget (T5): Ybar ~ 7e-12, so the margin is ~3x, resting on S1 reference values for the T7 constants.
- **Dominant term at level 3:** the inflow constant 1.6e6, from the weight conversion (1+x^2)^-2 -> far-field weight at the ramp.
- **Lever:** an E_I weight decaying like (1+x^2)^-1 in P2 (a T7 decision) gives Ybar ~ 2e-10, about 10x.
- **Corrections to T5:**
  - the far-end weight must be even at l~ (flattened weight);
  - T7 must reserve half the dissipation surplus of (Y), so Xi^2 roughly doubles.
- **Open items:**
  - the nonlinear class w3^2 w4 at level 3 is not yet written; level 4 stays in the bootstrap until it is;
  - the T7 inputs are placeholders;
  - the arb plan is 3-4 CPU-h; T1 must re-polish with storage.
- **Next:** T7. First a float study of coercivity with W = (1+x^2)^-1 vs (1+x^2)^-2 along the base; then the arb certificate prototype.

## 15. T7 prototype (Nov 10): coercivity certified at k = 20 (t7/T7_REPORT.md)
- **Weight:** keep W = (1+x^2)^-2.
  - T8's x10 lever held the J's fixed. Recomputed per weight, the best slow variant gains only x1.17.
  - Pure (1+x^2)^-1 collapses at large k.
- **Requirement and budget:** Ybar <= 2.7e-11 (2.0e-11 with T8 placeholders). The budget is ~2.3e-12, dominated by P1's start value. With T3's improved P1 numbers it is lower still. Margin >= ~x11.
- **Certificate in arb:**
  - c(0.03) = 0.0676 at k = 20.4, so m_c = 0.038; also certified at k = 40 (0.0708) and 2e4 (0.0679).
  - Method: the block is the whole Hermite-cubic FEM space (N = 743), with Hermite-interpolation tail estimates and Haynsworth inertia (congruence + Gershgorin).
  - Loss vs float 0.0013, scaling like h^2. 30 s per certificate.
- **Plan:**
  - (C2) with element-wise margins (T5's global s_c is 15x too crude);
  - about 29 certificates in total;
  - T7 total 4-7 CPU-h, mostly the (C4) far-curve constants, shared with T8.
- **Open:**
  - the prototype used the float state, not the polished one;
  - T6 corrections to the clock direction and to b1;
  - the Garding eps-level constants are placeholders;
  - phi1, phi2 come from the late base, and T9 must use the same.
- **P1 production:** slab 24, t = 0.375, error radius 4.5e-22.
- **Next:** T1 re-polish with lossless storage (running). T8 is closing the w3^2 w4 class. Then T6, T7 production and the T8 arb run on the stored base.

## 16. Nov 10 (late): base stored; T8 closed; T6 and T7 production running
- **Base stored:** t1/polished/ holds all 188 slabs, arb-128 lossless (55 MB). Slabs >= 5 meet all T6 targets, so P2 starts at slab 5 (tau_a = 0.75, k ~ 24). The P1 target is extended to k = 25.
- **T8 Lemma N3:** the level-3 classes w3^2 w4 and w3^3 close without level 4, so level 4 is out of the bootstrap. Ybar <= 2.7e-11 against a budget of 2.3e-12 (x11). T7 and T8 agree exactly.
- **T8 flags:**
  - lambda < 2 mu1 holds by only 2.7% at k ~ 20; T7 is to certify it from k ~ 24, else move k_a;
  - K_e and the N3 coefficients a0, a1 need T7's n2-n4.
- **Running:**
  - T6: rigorous residual enclosures and jets for slabs 5-187;
  - T7 production: 29 certificates, (C2) element-wise, (C3)/(C4) constants, t7/constants.json;
  - P1 production to k = 25.

## 17. P1 production: error diagnosis and M-rule fix (Nov 10, evening)
- **Status:** slab 90 at t = 4.59775, k = 15.0; C2 4.25e-15 at k = 13.45, about 9x the T3 estimate.
- **Cause of the C2 excess:** wrapping in slabs 10-68, from slab 1's absorption spread into near-neutral directions; x1.1e3 from wrapping alone. It is frozen in; only a restart from Gamma0 could cut it (x5-20, untested, 3.8 h).
- **Extrapolation:** C2 = 7.6e-14 at k = 24, against T5's 3e-14. The overall Ybar margin stays about x5. Accepted.
- **C5:** use the "lead" figure for Fbar_j(tau0). With the old M rule (+80 per slab), C5 would reach about 1e-8 at k = 24 because the box takes over.
  - Fix (a): M step max(80, 0.15 M), Mmax 3600. Applied at the slab-90 boundary (run_p1.py; backup .bak_mstep).
  - Expected C5 lead 2.2e-10 (3.5e-10 with far-field terms) at k = 24; budget 1e-9.
- **Restart command:** `run_p1.py --target-k 25 --Mmax 3600`.

## 18. T6 complete (Nov 10, evening): residuals certified, slabs 5-187 (t6/T6_REPORT.md)
- **All norms certified** (max over slabs, margin vs target):
  - R2 5.7e-16 (35x), Rs 4.2e-14 (1200x), R1 6.8e-15 (74x), R^I_s 3.1e-18.
  - Lemma F R0..R4: 5.8e-9, 5.8e-11, 3.7e-10, 5.8e-9, 1.4e-7. The tightest margin is 17x, at slab 182 (knot insertion).
  - 51 CPU-min, run here.
- **Method notes:**
  - the sup over tau uses the exact time polynomial, because between Radau nodes the far residual is ~1e3x larger;
  - certified >= resid.py estimates, by at most 1.09x.
- **T7 corrections provided:** Vtilde = V* + r; b1 gains the tangential residual A. The bounds are >= 1e7 below T7's margin.
- **Flags:**
  - the far end is off-axis by <= 4.6e-9, and the defect b1(l~) + l~_tau <= 1.2e-9. Lemma F needs a boundary term (T8 asked).
  - T8's fb.py must use the polished base (float far-field noise is up to 50x).
- **Running:** T7 production, T8 certified evaluation, P1 production (k -> 25).

## 19. T8 certified (Nov 10, 16:20 PT): Lemma F constants in arb, slabs 5-187 (t8/T8_REPORT.md, t8/constants/)
- **Requirement:** Ybar <= 1.60e-11 (item 3 binds) against a budget of ~2.3e-12, so the margin is x7.0. It becomes x9.4 once T7's provisional C_gamma0 is replaced.
- **Run:** 186 CPU-s.
- **Certified amplification:** A0..A3 = 6.7e3, 7.6e4, 1.6e6, 1.46e8. The good-time item-4 limit is 3.4e-11.
- **Boundary terms:**
  - the far-end transport coefficient is <= 6.9e-13, handled by Lemma B;
  - the 4.6e-9 off-axis offset of the base's far end enters only through the error's symmetry at l~. Lemma F assumes hypothesis (SE). P2 must compare with the base rotated by ~1.6e-14 rad; to be settled in the P2 driver/T9/T4.
- **Open:**
  - Lemma N3 needs a sharper far box (kappa e ~ 0.19 on the bend; T7 flag 4);
  - T4 must supply F_j(tau0) (damped by ~0.036);
  - T7's mu_e, rho_c and K_e are still placeholders. Rerun merge.py and camp.py when t7/constants.json lands (~2 s).

## 20. P1 production complete (Nov 10, 19:35 PT)
- **Run:** 96 slabs (0-95), 5.49 CPU-h. The final k_max is certified at 26.2537786005.
- **Final enclosures:** t = 4.597846334999 +- 6.9e-19 (time since Gamma0; clip time = t + 0.001), L = 48.32562154110 +- 4e-15, X0 = 0.84121520310 +- 1.5e-15.
- **Errors at the end:** C2 1.13e-13 (T5 assumed 3e-14), C5 lead 3.52e-10 (budget 1e-9). At k ~ 12: C2 2.5e-15.
- **Switch:** at base slab ~7 (k 26.25). The assembly agent computes Y(tau0) exactly. The expected margin against Ybar <= 1.6e-11 is ~x2-3, rather than x7.

## 21. Assembly stopped at the switch (Nov 10 night): base seeded from the wrong trajectory
- **Finding:** the T1 base was seeded from the HLS float snapshot p2_k20 (the original clip's float run), not from P1's validated trajectory (Gamma0 with slaved modes).
  - Mismatch at the switch (P1 slab 95, k = 26.25 vs base slab 7): Y(tau0) ~ 2.8e-9 (float, minimised over scale/shift/clock). That is x173 over Ybar <= 1.6e-11.
  - Clock offset 7e-6. The P1 ball itself contributes <= 1e-11.
- **Consequence:** P2 cannot close on this base. The true mismatch decays at ~0.53 per unit tau, but the comparison argument only certifies ~0.05.
- **(SE) also inherited from the seed:** the far tangent is tilted by 1.65e-14.
- **Decision:** re-seed the base from P1's exact 128-bit centre (slab 95, or slab 93 with warm-up). Store it in t1/polished_v2, then rerun T6 (~1 CPU-h), T7 (re-certify; heavy batch paused), T8 (minutes), then the assembly. Expected Y(tau0) ~ 5-10e-12, margin ~x2-3.
- **Lesson:** the P2 base must be a continuation of the validated P1 trajectory, not an independent float run of "the same" curve.

## 22. Base rebuilt from P1 (Nov 11 morning; t1/POLISH_V2.md)
- **Downtime:** a server restart at ~20:50 PT killed the rebuild. It was idle ~10 h until it was resumed at 06:45.
- **polished_v2:** slabs 0-180, tau 0 -> 27.15, k 26.2538 -> 23718.7.
  - Slab 0 = P1 slab 95's centre exactly (192-bit conversion, 128-bit storage); it matches P1's curvature to 4.7e-17.
  - It meets the T6 targets from slab 0, so no warm-up is needed.
- **Float switch mismatch:** Y(tau0) ~ 1e-14 (float floor), against 2.8e-9 for the old base. The certified Y(tau0) will be dominated by the P1 ball (C2 1.13e-13).
- **Residuals:** no misses on 24 sampled slabs. Tightest: slab 0 R2 x8; slab 160 RF3 x10, RF4 x38.
- **(SE):**
  - exact at slab ends (tilt <= 4e-28);
  - inside slabs, tilt <= 3.2e-18 and offset <= 1e-15 (slab 0; <= 2.8e-20 / 2.7e-16 for k >= 100), because integral(kappa) is conserved only at the collocation nodes;
  - the assembly needs the STEP4_RESULT §2(ii) correction.
- **Running:**
  - T6 v2 (t6/run_v2_chain.sh -> t6/T6V2_DONE);
  - then T7 v2 + T8 rerun (T7 agent, re-planned for slabs 0-180);
  - then the assembly.
- **Old-base T7 results (superseded, indicative):** min c_I 0.0528; lambda < 2 mu1 by 18.7%; rho_c >= 0.0381; K_e <= 891; C_gamma0 <= 0.389; T8 Ybar <= 1.96e-11.

## 23. T7 v2 complete (Nov 11, ~11:30 PT)
- **Certificates:** all 26 certified on tau in [0, 27.15]. Min c_I 0.0528; min rho_c 0.0471, after a finer C_gamma0 recomputation on 34 slabs before element insertions.
- **Constants:** lambda < 2 mu1 holds by 19.0% (lambda_max 0.2520 vs 0.30). K_e <= 890.6; C_gamma0 <= 0.207.
- **T8 on v2:** Ybar <= 2.092e-11 (item 3 binds; good-time window slabs 177-180). rho_c does not enter Ybar.
- **Open (to the assembly):**
  - BASE_T1 / KCH for v2;
  - Lemma N3 per x-range far box (x^2 n2 <= 6.2e14).
- **Assembly resumed on v2:** switch (arb), P2 comparison run, T9, STEP4_RESULT.md.

## 24. FINAL ASSEMBLY: end-to-end chain certified modulo 4 gaps (Nov 11, ~14:00 PT; STEP4_RESULT.md)
- **Switch:** P1 was run one more slab (96, inside base slab 2).
  - Y(tau_a) <= 6.20e-13 (raw 2.2e-14).
  - The IFT contracts (q 4.2e-9); clock shift <= 1.2e-12.
  - F_j(tau_a) <= 3.7e-23 ... 3.2e-15.
- **(SE):** holds exactly with a curvature bump (|a2| <= 1.6e-14).
- **P2:** Ybar = 1e-12 (residual level 1.33e-13). The requirement is 1.948e-11 (item 3), so the margin is x19.5.
- **Boxes:** tip x950; far |e'| x18, x|e''| x133. Binding: far level 3, x^2|e'''| <= 5.5e-3 vs 1e-2 (x1.8).
- **Hand-off:**
  - item 5: base certified (T6 hulls; worst x44);
  - items 1-4: error parts certified, base values float (gap 7);
  - thinnest margin |alpha| x1.14.
- **Theorem 4 (draft):** t1 in [4.5978580262, 4.5978580263] since Gamma0.
- **Gaps:**
  - (7) arb base values for items 1-4 and KCH — assembly;
  - (8) Lemma N3 per x-range — T8 (the numerical 4% theta3 loss goes into the assembly);
  - (9) T7 tip box |e''| <= 1e-2 not implied pointwise by Lambda_e — T7;
  - (10) T7 far constants on the corrected base (dkappa <= 3e-17) — T7.

## 25. Assembly update (Nov 11, afternoon)
- **theta3 loss included:** theta3 0.240 -> 0.220.
- **KCH in arb = 1.152** (not 1.02). The requirement becomes Ybar <= 1.725e-11 (x17.3 at Ybar = 1e-12). The binding certified margin is still the far box at level 3 (x1.8).
- **Cone jets in arb:** H3 K4 4.5e-3 vs 0.06 (x13). H5 cone K2/K3/K4 <= 1.9e-2 / 4.8e-4 / 9.8e-3 vs 0.12.
- **Blocked:** alpha by IFT and sqrt(Lambda), and with them the windows and shells (which need alpha to +-0.3%), need step 3's f3, psi3 and h2 in certified form.
  - Pointed the assembly to step 3's interval representations (agent_mod/ibp boxes, verifier rhs_iv, step3/ linear certificate).
  - The base values for these items are float until then.
- **Gap 10:** waiting on T7. **Gaps 8 and 9:** with T8 and T7.

## 26. Hand-off H1 needs a tighter S1 on sigma in [8,13] (Nov 11, late afternoon)
- **What step 3 relies on:** interval hull boxes for f1, f2, f3, psi3, h2 and rho (ibp_boxes.py, from the wide CAPD run). S1 itself is tight only for sigma <= 8 (s1_tight); beyond that there are only tail-lemma envelopes (~5e-5).
- **Effect:** alpha can be pinned to +-0.4% (f3 is usable as is). sqrt(Lambda) only to ~3.5e-6, against the needed 5e-8 (x70 short).
- **Commissioned:** a tight S1 enclosure on [8,13] (width <= 1e-8), into /home/claude/s1_tight/ext813/. Routes: outward integration from the 4.4e-13 sigma = 8 state, or sharper tail bounds.
- **Meanwhile:** the assembly builds the windows and shells now (alpha +-0.4% suffices; Theta certified to +-2e-7).

## 27. Assembly near-final (Nov 11, evening)
- **alpha certified:** -1.12964e-4 +- 5.2e-9 (|alpha| <= 1.13e-4 vs 1.5e-4). THINF 1.8617364190 +- 1.7e-9. The hand-off window is slab 180, t_c in [1-1/512, 1].
- **Windows and shells in arb:**
  - windows E0-E3 <= 1.44e-8, 5.19e-8, 1.70e-7, 4.70e-7 (vs 2e-7 and 1e-5);
  - shells <= 2.8e-12, 3.7e-12, 8.9e-11;
  - base K4 3.2e-4.
- **Gaps 9/10 folded in:** T7 Lemma 4a; 3e-9 of eps D4^2 reserved; T7's republished constants used.
- **Rounding fix:** se_corr uses nextafter, with outward rounding throughout p2/. Y(tau_a) <= 6.20e-13 unchanged.
- **Margins:**
  - item 4 binds T8 at Ybar <= 5.11e-12 (x5.1);
  - cone K4 x1.67;
  - far box level 3 x1.8 (Ybar <= 1.9e-12), binding overall;
  - Ybar = 1e-12.
- **Open:**
  - G1 sqrt(Lambda): certified 1.09e-7 vs 5e-8 (float 6e-10). Limited by the step-1 tip-end theta widths (~8.8e-10). Fixes: a scalar tip offset in the projection (assembly) and/or a tighter tip-end profile (S1-ext agent).
  - G4: N3 lower-order terms and Q (T8).
  - G5: far-cylinder item 5 composed by hand; to be scripted (assembly).

## 28. STEP 4 CHAIN CERTIFIED, no open gaps (Nov 11, ~16:30 PT)
- **H1:** certified over all of W' with 32 t_c-balls, using the tight tip profile (s1_tight/ext_tip) and the [8,13] extension.
  - sqrt(Lambda) <= 2.30e-8 vs 5e-8 (x2.18); at the W' midpoint, 6.72e-9.
  - alpha = -1.1296366e-4 +- 3.0e-9.
  - Error constants: C_Lambda <= 24.1, C_alpha <= 0.72 (t9h/clam.py).
- **G4 closed** (N3' levels 0-3 and Q). The requirement is Ybar <= 5.048e-12; the far box gives <= 1.85e-12 (binding); the actual Ybar is 1e-12. TOPLOSS is kept, so the top term is counted twice (conservative).
- **G5 closed** (t9_far.py).
- **Thinnest margins:** far box x1.85 (in Ybar), cone K4 x1.61, H1 x2.18, |alpha| x1.33.
- **Fix:** t9.py now writes the cone +0.01 c0 allowance itself.
- **Next:**
  - independent Fable review (STEP4_REVIEW.md);
  - kit B after the review;
  - the step-4 proof document;
  - the main-theorem write-up.
- **Kit A fix:** a Windows fsync on a read-only handle failed (EBADF). Patched in the kit and in t3/run_p1.py; v2 zip sent to the user.

## 29. Independent review of step 4 (Fable, Nov 11 ~16:15 PT; STEP4_REVIEW.md)
- **Verdict:** Theorem 4 is not proved as stated, but nothing found makes it false. Findings: 3 (A), 9 (B), 4 (C). The "no open gap" status is withdrawn.
- **(A) gaps:**
  - **A1 — Lemma O H4** (integral of (h0^(5))^2 <= M^2 per cylinder, required by step-3 H5) is unchecked. It fails on the certified D2(tau*) <= 1.22e-8 over the 1.46e-4-long W': cone x10.8-252, innermost far cylinder x4.8. Fix: a longer good-time window and/or a localized D2 (assembly + T8).
  - **A2 — the nonlinear terms are missing from the certified (Y).** closure_cert is linear and ycomp uses a typed DN = 1e-6. The far top-coefficient perturbation (sup A+ <= 0.0119) and the tip D4-terms need dissipation slices (T7).
  - **A3 — the tip pointwise box is attributed in a circle** and is false at the outer edge (1.35e4 Ybar vs 1.054e4 Ybar). Fix: a short lemma from Lambda_e (T7).
- **(B) fixes:**
  - B1: se_corr majorants underived (dominate the far residuals x309);
  - B2: true bootstrap margin |e'| x1.29, not x1.85;
  - B3: the H1 IFT lacks an existence test;
  - B4: camp typed values rounded the wrong way;
  - B5: typed allowances in T9;
  - B6: C_x omits residual/nonlinear parts of (M);
  - B7: unscripted steps;
  - B8: S1 non-uniqueness survives, but needs a hull-containment script and a theorem phrasing;
  - B9: rounding spots.
- **Checked and found in order:** Gamma0, the P1 start, the Galerkin-limit lemma, the switch, the step-3 definitions in T9.
- **Not audited:** slab.py internals, the T6/T7 engines, cconst/n3forms, the CAPD programs.
- **Assignments:**
  - assembly: A1 lead, B1-B3, B5, B7, B9;
  - T7: A2, A3;
  - T8: localized D2 for A1, plus B2, B4, B6;
  - S1 agent: B8 script;
  - writer: TODO markers and C items.
- **Kit A on the user's machine** (granted access to C:\Kit4a only): slab 2 matches the reference t enclosure to 22 digits; running.
- **B8 closed (Nov 11):** s1_tight/consistency/check_outA5.py (12 s). prof_tip, prof_tight and ext813 lie inside the step-3 outA5 hulls (ibp_boxes) wherever they overlap: 1025 tip rows, 501 far rows. Every downstream object therefore refers to the same certified set of shrinkers.
