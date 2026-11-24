# Step 3 fix plan (Nov 9, late) — for findings 1–7 of STEP3_REVIEW.md

Status of the evidence: everything below is a FLOAT/design estimate unless tagged [C]. My scripts and logs are in `step3nl/fixplan/` (patched copies; nothing else was modified). The closure runs there use the verifier's code with its "conditional" J integrals (see V1).

## 0. Decisions

**Hand-off: α1 = 5e-4, √Λ1 ≤ 2e-7** (fallback 3.5e-7). One change of parameters removes or shrinks four findings:

| quantity (verifier run, α1 = 5e-4) | √Λ1 = 2e-7 | √Λ1 = 3.5e-7 | limit |
|---|---|---|---|
| peak √Λ | 2.49e-7 | 4.18e-7 | cap 7.5e-6 |
| S_ramp | 3.5e-3 | 5.9e-3 | 1e-2 |
| sup\|w′\| on [16.25, 34.25] (E_I) | 5.4e-3 | 9.1e-3 | **1e-2 (old boxes hold: finding 7 disappears)** |
| tail ρ (for Λ) | 2.61e-3 | 2.08e-3 | > 0 |
| J0, J1, J2 (conditional) | 1.37e-11, 1.93e-11, 4.27e-12 | 4.0e-11, 5.6e-11, 1.24e-11 | — |

(`fixplan/run_base_L2e-7.out`, `run_base_L3.5e-7.out`.)

**Theorem version proved by this plan: W (trapping with a floor).** For all τ ≥ τ1: |α| ≤ α1, Λ ≤ Λ_peak, modulation valid, far curve a graph with bounded slope and curvature ≤ C/r. This gives λ → 0, T < ∞, a single singular point, λ = (T−t)^{1/4}(1 + O(sup|a|)), and the rescaled tip within ~1e-6 (H²_ρ) of the S1 family for all time, eventually within the floor √Λ_F ≈ 2e-8. The only absolute forcing is term (c) (F4).

**Upgrade C (convergence to S1)** costs one extra analytic lemma and no new certified constants: (c) is bounded by W(57)·(shell energies of (A)), and those obey Ė_j ≤ −γ′E_j + cK_jV_j² + forcing(β, a, b)². The joint system (α, Λ, E^A) is homogeneous of degree ≥ 2 with loop gain ~K_j·W(57) ≤ 1e10·1.4e-18. Size: small-to-medium, do it last. Rate would be ≈ ρ/2 ≈ 1.3e-3 for √Λ (not 0.011).

**Far field: option (ii)**, modified (F1). Option (i) needs a re-certification of L-int/Gårding with O(0.1) coefficient deviations (large) and still needs level-2 inflow; option (iii) only moves the blend region outward. All three need sharp time integrals at level 2, so V1 is unavoidable; given V1, (ii) is the cheapest.

## V1. Verifier upgrade: unconditional sharp J0, J1, J2, J4 (critical path)

- **Argument.** For every state with s = √Ly ≤ s_up (no lower bound needed): H − (1−δ)NEG_θ ≤ r1·s², with r1 the B&B sup at s = s_up (valid for smaller s because H/s² is nondecreasing and NEG/s² is s-independent, both already certified by `tail_cert.py` mode s). L has degree exactly 1, so L − δ·NEG_θ ≤ Φ(α) := sup_shapes l̂²/(4δ n̂), an absolute α-dependent constant. Hence Λ′ ≤ 2r1Λ + 2Φ(α) − θ·NEGV on the true trajectory using only Ly ≤ barrier, and θ∫NEGV ≤ Λ1 + Σ2r1⁺L_upD + ∫2Φ. Add the V4 dissipation μ(ε_I/2)V4² to NEGV to get J4.
- **Float estimate.** Φ ≈ 2·540·Fn(α)²/δ (the 5.72·Fn·A term dominates). At α1 = 5e-4: Fn = 1.9e-9, ∫2Φ ≈ 5e-15/δ = 2.6e-14 at δ = 0.2, against Λ1 = 4e-14 and the conditional total 5.1e-14. So unconditional J ≈ 2× the conditional values. At α1 = 1e-3 the same term is 1.2e-12, 10× Λ1: another reason for 5e-4. **Design values used below: 3× conditional**, J = (4.1e-11, 5.8e-11, 1.3e-11), J4 = 4.3e-7.
- **Tasks.** One extra B&B pass per step on H − (1−δ)NEG_θ; one sup for Φ per α-level; V4 in NEGV. Size: small-to-medium.

## F1. Lemma S(A) redesigned (finding 1)

- **Design.** Deviation shells v = w′, levels 0–2, core y ∈ [28, 60], inner ramp y ∈ [14, 28] (x ∈ [16.25, 30.25]: ζ ≡ 1 and W = (1+x²)⁻² exactly; s-length ln 2), γ_w = 1.
- **Argument.** Write every level-j estimate in divergence form: test ∂_s^j ṽ, move two derivatives to the test function. Then
  - the pointwise coefficients are (V, V_s) only, i.e. the bootstrap bounds |V| ≤ 1/2 and y|V′| ≤ S1 box, plus jets of U up to U⁗ (u0 ∈ C⁴, piecewise C⁵). **No V-jets beyond first order; the 7th-order box of Prop. 1 is removed.** Cubic terms use ‖ṽ_ss‖²_{L⁴} ≤ 3 sup|ṽ_s|·‖ṽ_sss‖;
  - with Θ = ωη²e^{−4s}: −∫Θ w ∂⁴w = −∫Θ(w″)² + 2∫Θ″(w′)² − ½∫Θ⁗w². So the ramp commutators involve only orders j and j+1 on the ramp, with coefficients linear in η-derivatives. Bound them by Cauchy–Schwarz against J_j, J_{j+1} (J3 ≤ √(J2J4) + c1W·J2). No Young against the shell dissipation, hence no ϑ′⁴/ϑ³ constants (the 2.66e7 of the review).
- **Float estimate** (`fixplan/shellA_est.py`; linear model, commutators inflated ×3; hand-off energies not included):

| inputs | E0, E1, E2 (core) | S0_v = sup\|w′\| | S1_v = sup y\|w″\| |
|---|---|---|---|
| (5e-4, 2e-7), conditional J | 8.9e-7, 4.5e-4, 34 | 0.0064 | 0.50 |
| **(5e-4, 2e-7), J×3 (design)** | 2.7e-6, 1.4e-3, 103 | **0.011** | **0.87** |
| (5e-4, 3.5e-7), J×3 | 7.8e-6, 3.9e-3, 295 | 0.019 | 1.47 |
| weak rigorous J of VERIFIER.md | — | 3.4 | 1.7e3 (fails) |

  Requirements: S0_v ≤ 0.02 (existing T3 box) and S1_v ≤ 3 (F2). Margins 1.8× and 3.5× at the design point. E2 is dominated by the level-2 commutator through J3 (i.e. J4); if needed, widen the box to 0.04 (F7).
- **Forcing.** As in lemmaS.tex §3; negligible (≤ 1e-9 in E0) with the β of F3.
- **Hand-off (d), now explicit:** E0(τ1) ≤ 5e-7, E1(τ1) ≤ 3e-4, E2(τ1) ≤ 20 on y ∈ [14, 60].
- **Tasks.** (1) Prop. 1 in divergence form, symbolic derivation + interval sups of the η-coefficients — medium (largest single item). (2) K-type constants for the new ramp — small (reuse `lemmaS_interval.py`). (3) Assembly script with rigorous J's — small.
- **Fallback if V1 slips.** Levels 0–1 with pointwise inflow (V0, V1 bounded by Λ) give S0_v ≈ 0.011–0.026 at peak 2.5e-7 (level-1 commutator through pointwise V2 is the weak point); level 2 has no J-free version. Not recommended.

## F2. Pointwise data on x ∈ [34.25, 60] (finding 2)

- **Argument.** In `rhs_iv.py` use sup|w′| ≤ max(sI1, S0_v) and sup|w″| ≤ max(sI2, S1_v/32) as coefficients; S0_v, S1_v are bootstrap constants supplied by F1 (the (A) core starts at x = 30.25, so there is overlap with the K_W region). Remove sup|w‴| altogether by the exact identity (integrate ∫W v″⁴ by parts once): ‖v″²‖_W ≤ sup|v′|·(3V3 + sup|ψ′|·V2), sup|ψ′| ≤ 1.31.
- **Float estimate** (`fixplan/cell_test.py`, tail-anchor cell α ≤ 7.5e-5, Λ = 1.74e-13; baseline U = −0.00199):
  - GN replacement of supv2: U unchanged;
  - S1_v = 3 (coefficient 0.094 on the whole far region, conservative): U unchanged;
  - S1_v = 6: U = +0.49, fails. So the requirement is **S1_v ≤ 3**.
  - Full march at (5e-4, 3.5e-7) with S1_v = 3, GN and the F3 terms (q = 50): closes, tail ρ = 2.3e-3 (`run_mod.out`).
- **Tasks.** Two-line change in `rhs_iv.py` plus the identity in the write-up; ramp versions of the identity (ζ′ terms) — small.

## F3. Modulation (finding 3)

- **Answer: yes, integration by parts is enough.** 𝒩 is quasilinear: the w⁗ coefficient depends on (w, w′) only, the w‴ coefficient on jets ≤ 2. Moving derivatives onto χf_jρ (smooth; u0 ∈ C⁴ suffices) leaves
  - a linear part, exact in α: ⟨w̃, L*f_j⟩ + α⟨w, (∂_αL)*χf_j⟩ ≤ (c_j + m_j|α|)A. Its β-loop effect in E0 is already the certified 6.19|α| form;
  - a quadratic part |Q_i(w,w)| ≤ q_i(A + W1 + B)² (forms in w, w′, w″ only, sup-free);
  - a cubic remainder ≤ (box)·quadratic, plus sup|w″|·B² from w″³ (allowed: shape-dependent coefficient, sublinear in D4).
  No D4 or W3 terms; a, b, β are pointwise bounded by Λ-quantities. **Jet box needed: only |w|, |w′| ≤ 0.02 on supp χ, and only for the cubic remainder** (a factor 1 + O(box)); 7.57 is replaced by 3.708 + m|α|.
- **Float estimate.** Expected q_β ≈ 31·sup|f3|·(d200 ≈ 4.8) ≈ 3–10, q_(a,b) ≈ 10–30. Sensitivity of the tail rate to the new E0 term n3·q_β(A+W1+B)²A: +0.0034 at q = 300, +0.0118 at q = 1000 (√Λ = 4.2e-7); i.e. ≈ 23·q_β·A. At the 2e-7 hand-off (plateau 2.5e-7, ρ = 2.6e-3), q_β ≤ 200 keeps half the tail margin: 20× room. The β ODE gets + q_β·1.13e5·Λ; the tail α-condition has margin 39 with q = 50. Effects of q_(a,b) in I0/I2 are ≤ 1e-3 of the dissipation.
- **Downstream.** ∫|a|dτ is no longer needed (F5: windows fixed in physical space); only sup|a|, sup|b| ≲ 1e-5 enter.
- **Tasks.** Symbolic IBP of the three pairings and interval bounds for c_j, m_j, q_i — medium (reuse `nlin_E2.py` forms and `lem_integrals.py`). Insert in `rhs_iv.py`/`closure.py` as in `fixplan/rhs_iv.py`, `closure_fix.py` — small.

## F4. Term (c) (finding 4)

- **Argument.** Treat (c) as absolute forcing (the verifier's non-QUAD path for (c); (b) stays quadratic). I0 flux: needs ‖w′, w″, w‴‖_{L²(57,60)} pointwise in time, which are the (A) energies E0, E1, E2 on y ∈ [54.75, 57.75] (inside the core). I2 flux: needs ‖w⁗‖_{L²(ramp)}, available only time-integrated from the level-1 dissipation of (A); as a forcing term an L¹_loc-uniform bound suffices (one-line extension of the tail lemma).
- **Float estimate.** With lemmaS_cflux.json: μ·6.1e-18 in I0, μ ε_I·1.5e-17 in I2; floor Λ_F ≈ μ·6.1e-18/ρ = 1.9e-16, √Λ_F ≈ 1.4e-8 (peak 2.5e-7). Even 100× larger ramp norms give 1.4e-7. The J's become affine in T with slope ~1e-16; (A) is damped (γ′ ≈ 1/4), so sup E_j changes by < 1e-9 relative.
- **Tasks.** Recompute the two absolute constants from the F1 energy bounds instead of an assumed box — small. Tail lemma with L¹_loc forcing — small.

## F5. Lemma O geometry and the outer far field (finding 5)

- **Design at α1 = 5e-4:** y_O = 200, cylinder half-width 0.75·r (L = 50λ1 at the innermost centre, cylinder y ∈ [50, 350]), X = 3, M = 0.12. S(B) on [48, R] with R1 = 175, outer ramp [150, 175] inside the Lemma O conclusion region [150, 250].
- **Float check.** Time ratio (3/150)⁴ = 1.6e-7 against t* ≥ 3.70e-7 [C]: margin 2.3. Tangent-frame slope of the bent cone across the cylinder ≤ 0.120 against the hard limit 0.2437: room 0.12 for the data. Curvature in cylinder units 0.042 against M = 0.12. S1-ray chart slope at R1: 0.14 (≤ 1/2). Alternative y_O = 200, half-width 0.675·r: slope 0.108, time margin 1.5.
- **Argument for H6 and (B).** Take every (B) window fixed in physical space (y_edge′ = (¼ − a)y_edge). Then the a-flux vanishes and the b-flux costs ∫|b|e^{−s}dτ ≤ 4 sup|b|/y_entry, so (B) needs only sup|a|, sup|b| (floor-compatible). Window energies change by a factor e^{C y_entry⁻⁴}. The Lemma O inner strip y ∈ [50, 150]·λ1 consists of windows present at τ1, so **its K2, K3 bounds for all t are the hand-off data times e^{C·50⁻⁴}**; they do not depend on (A) or on J. Windows entering later through y = 48 stay at r < 48λ1 and need only finite bounds (|f| ≤ 1/2 from levels 0–1; level 2 from E2^A ≤ ~100, any finite value gives curvature ≤ S1_f/r).
- **Hand-off additions.** H2 with the true limit; on y ∈ [48, 350]: y|f′| ≤ 0.12 and ∫(y²f″)²dy/y small enough for K3 (step 4 data); T − t1 ≤ λ1⁴/(1 − 4 sup|a|).
- **Tasks.** (B) Prop. for physically fixed windows, levels 0–2, same divergence-form computation as F1 without ramp inflow — small once F1 exists. Region table and the two strict-improvement statements (O ⇄ B) — small (writing). Lemma O constants unchanged [C].

## F6. Decay (finding 6)

Use ρ from the run: Λ decays at ρ = 2.6e-3 (2e-7 hand-off), √Λ and α at ρ/2 = 1.3e-3. In version W nothing downstream needs a rate (F5 removed ∫|a|). In version C the rate is min(ρ/2, γ′/2) = ρ/2. Delete "0.011" from the blueprint. No task.

## F7. Jet boxes (finding 7)

At (5e-4, 2e-7): sup|w′| ≤ 5.4e-3 on [16.25, 34.25] and S_ramp ≤ 3.5e-3, inside the existing 1e-2 boxes with factor 1.8–2.9; S0_v ≤ 0.011 on [34.25, 60] inside the T3 box 0.02. No re-certification needed. Optional safety (small, ~minutes each): rerun `cubic_E0fine`, `nlin_E2`, T3 with α ≤ 5e-4 and boxes 0.02 (inner far) / 0.04 (T3); with α halved the coefficients do not grow.

## Bootstrap ring (to be written, not certified)

Assume on [τ1, τ*]: Λ ≤ Λ_peak, |α| ≤ α1; S0_v ≤ 0.02, S1_v ≤ 3 on y ∈ [28, 60]; |f| ≤ 1/2 and window bounds on [48, R]; Lemma O K-bounds on r ≥ 150λ1. Then: closure (needs S0_v, S1_v, (c) constants) ⇒ Λ, α, J strictly inside; (A) (needs J, |V| ≤ 1/2 on its outer ramp from (B)) ⇒ S0_v ≤ 0.011, S1_v ≤ 0.87; (B) (needs (A) on [48, 56], Lemma O at [150, 175]) ⇒ strict; Lemma O (needs (B) on [50, 150], hand-off) ⇒ strict by its own slack. No arrow uses a time-integrated quantity as a coefficient.

## Order of work

1. V1 (unconditional J, J4) — everything in the far field feeds on it.
2. F1 Prop. 1 in divergence form + constants; in parallel F3 (modulation IBP lemma).
3. F2 + F3 + F4 edits in the verifier, then one production run at (5e-4, 2e-7) and one at (5e-4, 3.5e-7).
4. F5: (B) proposition, region table, updated hand-off list for step 4.
5. Bootstrap/assembly write-up; then upgrade C if wanted.

## How much of step 3 remains

Certified and reusable: linear coercivity, all inner and E_I constants, forcing, the verifier with its rigor pass, Lemma O constants. Not yet existing in certifiable form: Lemma S (Prop. 1 for (A) and (B)), the modulation lemma in usable form, unconditional J's, the (c) constants, and the whole bootstrap write-up. My estimate: roughly 60% done. Remaining: 2 medium tasks (F1, F3), 1 small-to-medium (V1), about 6 small ones, and the assembly write-up. Main risk: the level-2 shell energy (E2 ≈ 100, dominated by J4). If the certified Prop. 1 constants exceed my linear-model estimate by more than ~10× in E2, S1_v passes 3; the remedies then are a lower Λ1, or dropping sup|w″| from E_I by a weighted Gagliardo–Nirenberg bound in sup|w′| (medium).

## Caveats on this plan

- `shellA_est.py` uses the constant-coefficient model of the fourth-order operator with a ×3 allowance; the quasilinear terms are not in it.
- q_β, q_(a,b) are order-of-magnitude guesses; only the tolerance (≈ 200) is measured.
- The patched runs are not rigorous (float-level edits of a rigorous code) and used conditional J's; the stress run (q = 300, S1_v = 6) was stopped after its first cells showed rate +0.49, which the cell tests attribute to S1_v = 6.
- The (5e-4, 2e-7) hand-off assumes step 4 can deliver it; the intrinsic distance at α = 5e-4 (≈ 170|α|³ ≈ 2e-8) is well below it.
