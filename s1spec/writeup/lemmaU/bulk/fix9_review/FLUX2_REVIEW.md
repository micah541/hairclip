# FLUX2_REVIEW: independent review of flux lemma v2 (fix round 9)

Nov 20, 2026. Fable, independent reviewer. The proof tree was read-only; everything I wrote is in this directory.

- Paths: `B` = `/home/claude/s1spec/writeup/lemmaU/bulk`; `V` = `/home/claude/s1spec/step3nl_v4b/verifier`; `D8` = `B/fix8_flux2_impl`; `N36` = `B/fix6_v36`.
- Labels: **[cert]** = arb with outward rounding (mine unless a file is named); **[sym]** = exact sympy; **[hand]** = my derivation by hand, not machine-checked; **[read]** = verified by reading the cited file:line; **[est]** = estimate.
- My scripts: `recompute.py` → `recompute.json` (arb, 160 bits, one process, nice 19, PYTHONHASHSEED=0, `python3 -B`, under 100 MB); `chi4_check.py` → `chi4_check.json` (exact sympy, 79 MB, 185 s).

## 0. Verdict

**Sound modulo listed items.** I found no error in the mathematics or in the certified numbers. The lemma is not yet a proof, because of the open items in §7, all of which were already known except the text items T1–T4.

- The verifier-side readings behind C2–C5 are correct. I re-derived them from `closure.py`, `bb.py`, `rhs_iv.py` and `TAIL_LEMMA.tex` (§1).
- The reflected-integral lemma is correct, and adding 2C′ the way `flux2_assemble.py` does is valid (§2).
- The level-4 chart structure holds: my own construction confirms that (u‴)² is the only quadratic high-jet term (§3).
- The k = 2 export is a valid certified computation inside (A)'s framework; its written proof is still owed (§4).
- All four departures from the design are valid (§5).
- My recomputation reproduces K2^pt, C′, S′, W(57), N̄ and every margin. K̄0 is confirmed as a valid bound and is conservative by ×3.4 (§6).

**Step 3 is not closed.** §7 lists what still stands in the way.

## 1. Verifier-side readings (C2–C5)

### 1.1 What a march cell certifies [read]
- `closure.py:129` calls `orc.rate(aq, fdn(Ldn), fup(Lup))`, which is `bb.bound_cell(0, aq, Ly_lo, Ly_hi, floors=True)` (`closure.py:27-32`).
- `bb.py:142-185` returns U ≥ sup of tot/Λ over all shapes, α ∈ [0, aq] and Λ ∈ [Ly_lo, Ly_hi]. The Λ range is handled by evaluating H − NEG at s_hi and L at s_lo (`bb.py:5-6`, `rhs_iv.py:391-398`; TAIL_LEMMA.tex, Lemma "s-monotone split").
- The flux constants sit in L as one degree-0 constant, μ·FLUX0 + μ·ε_I·FLUX2 (`rhs_iv.py:229`), so the cell bounds (tot_s + F0)/Λ, with F0/Λ taken at Ly_lo. The leaf argument (TAIL_LEMMA.tex, "Branch-and-bound leaf certificate") makes U an upper bound whatever the tolerances.
- So the reading "a cell bounds tot/Λ for every state with Λ ∈ [Ly_lo, Ly_hi] and α ≤ aq" is right.
- Cell ends: Ly_hi = max(L_k, g)·(1 + 1e-3) and Ly_lo = min(L_k, g)·(1 − 1e-3) (`closure.py:112-115`); the step is accepted only if L_{k+1} ∈ [Ly_lo, Ly_hi] (`closure.py:131-132`). The cushion above the barrier is therefore not automatically 1e-3·L when the barrier rises, so C2 has to be checked step by step against max(L_k, L_{k+1}). `flux2_assemble.py:93-107` does exactly that.
- Barrier: L′ = 2·rr·L with rr = U_k + 1e-9·max(1, |U_k|) (`closure.py:130`).

### 1.2 First touch on the march with Y = Λ − Z [hand]
Let G = Y − L on step k, and suppose G > 0 somewhere. Take t0 = the last time before it with G ≤ 0, so G(t0) = 0 and G > 0 just after.
- Just after t0: Λ ≥ Y > L ≥ Ly_lo, and Λ ≤ Y + 2C′ ≤ L + 2C′ + G < Ly_hi for G small, by C2 (strict, room ×5e7). The α-barrier is unchanged, because β is evaluated at Ly_hi (`closure.py:119`) and Λ ≤ Ly_hi still holds. So the state is in the cell and ½Y′ ≤ tot ≤ U_k·Λ.
- If U_k ≥ 0: Y′ ≤ 2U_k(Y + 2C′), so G′ ≤ 2U_k·G + (4U_k·C′ − 2·pad_k·L). The bracket is negative by C3, and Gronwall from G(t0) = 0 gives G ≤ 0. Contradiction.
- If U_k < 0: Λ ≥ Y gives Y′ ≤ 2U_k·Y, so G′ ≤ 2U_k·G − 2·pad_k·L < 0 at G = 0. No condition.
- Start: Z(τ_1) = 0, so Y(τ_1) = Λ(τ_1) ≤ Ly1 < L(0).
- C3 as coded uses Ly_lo_k as the lower bound of L on the step (`flux2_assemble.py:101`), which is conservative.
- The implementation's reading and conditions C2, C3 are correct.

### 1.3 Tail [hand, from `closure.py:157-207` and TAIL_LEMMA.tex "Tail contraction" and "Floor"]
- Anchor cell: rate_s ≤ −ρ on Λ ∈ [L*, (1 + η)L*], α ≤ a*, floors off, η = ETA = 0.01 (`closure.py:165`; env ETA=0.01 in `run_v4b.sh:7`).
- With G = Y − ℓ > 0 small: Λ ∈ (ℓ, ℓ + 2C′ + G]. Scale by λ² = ℓ/L*: the scaled state has Λ/λ² ∈ [L*, (1 + η)L*] as soon as 2C′ < η·L_F (C4), and α/λ ≤ a*. Joint scaling gives tot_s ≤ −ρΛ, so Y′ ≤ −2ρΛ + 2F0 ≤ −2ρY + 2F0.
- This is the differential inequality of the floor paragraph with Y in place of Λ. With ℓ′ = −ρℓ and ρℓ ≥ 2F0: G′ ≤ −2ρG. On the flat part: G′ ≤ −2ρG − ρL_F. So G ≤ 0, with the same L_F = 2F0/ρ.
- α-barrier: β is nondecreasing in Λ and each Λ-dependent term has degree ≤ 1 in Λ (`closure.py:22`: √Λ and QBL·Λ). With Λ ≤ (1 + 2C′/L_F)·ℓ, the right side of the whole-tail condition (`closure.py:185`) grows by at most the factor (1 + 2C′/L_F). C5 is correct.
- The tail anchor boxes are asserted at L* without cushion (`closure.py:162`). They are covered because L* + 2C′ ≤ Ly_hi of the last march cell, whose boxes are asserted at Ly_hi (`closure.py:142-143`), and a* ≤ alpha_up of that cell. I checked both: margins ×1.97e8 and ×5.1e7, a* = 3.68e-5 and 3.58e-5 against 5e-5 [cert, `recompute.json`].

### 1.4 Unconditional integrals (C6) [hand, from `closure.py:209-256`]
- The march uses r1 at the point cell Ly_hi with L dropped (`closure.py:222`), valid for every Λ ≤ Ly_hi (TAIL_LEMMA.tex §V1). With C2, Λ ≤ Ly_hi still holds.
- Then 2θ∫_a^b NEGV ≤ Y(a) − Y(b) + 2∫(r1⁺Λ + Φ + F0), with Y(a) ≤ L(a) and −Y(b) ≤ 2C′. The constant gains C′ and the slope nothing.
- κ_J enters as J_m·(1 + κ_J) with κ_J = (…)/X_const (`N36/t2_A/rigor/hjw_rep.py:38-39`), so an absolute pad of 1e-9 on κ_J covers C′ ≤ 1e-9·X_const. Margin ×141 [cert, mine; same as reported].
- If someone bounds the march and the tail pieces separately instead of telescoping, the gain is 2C′; the margin is then ×70 [cert]. Either way it passes.

### 1.5 Other consumers (C7) [read]
- `mod_sups.py:33-37`: march rows at Ly_hi; one tail row at L*. The scripted reserve check is the right one.
- `N36/v36cfg.py:66-84`: J's, slopes and peak rounded up to 5 digits, ρ down.
- `assembly_arb.py:170-183` (base_ibp): raw sqrt_LF enters PHI_S, which feeds only I3c/I3s. `assembly_rep.py:295-309` uses them for the v3.1 page claims and the honest I3 only. Flux lemma v2 does not use I3. Confirmed.
- No march step of either production table used the 1e-6 "slack" path of `closure.py:126-128` (n_slack_steps = 0 [cert, `recompute.json`]), so that pre-existing [est] device is not in play.

## 2. The absorption lemma

- **Lemma (reflected integral): correct** [hand]. Z(τ) = 2·sup_{σ≤τ} ∫_σ^τ (φ − s) ≤ 2C′ + 2(S′ − s)(τ − σ) ≤ 2C′ when S′ ≤ s. The running minimum m of the absolutely continuous P is nonincreasing and satisfies |m(t) − m(r)| ≤ ∫_r^t |P′|, so it is absolutely continuous and Z′ = P′ − m′ ≥ P′ = 2(φ − s). Subtracting from (★) gives ½Y′ ≤ tot.
- **Why an integrated term may enter a pointwise verifier.** The certificates are statements about all states of a cell, not about a trajectory. Y obeys the pointwise inequality with the constant fluxes, evaluated at a state whose energy lies in [Y, Y + 2C′]. Each certificate tolerates that shift: the cell cushion (C2), the strictness pad (C3), the tail's monotonicity in Λ and the anchor width η (C4, C5), and one additive constant in the integrals (C6). I agree with the design's paragraph; I found no consumer it misses.
- **Bootstrap order (C9): valid** [hand]. Under the bootstrap set with (BΦ), the flux lemma and the absorption give the verifier's conclusions; (A) then returns X_k ≤ C0_k + C1_k·T, a strict improvement. (BΦ) is closed (monotone convergence) and open (the extra piece ∫_{τ′}^{τ″} is independent of the window start and tends to 0, given local integrability past τ′).
- **The typed ā = 1e-6 must be a bootstrap hypothesis, not a remark.** It is improved strictly by the verifier's conclusion (2.14041e-7; `mod_sups.py`), so this works, but no continuation text states it yet (T3).

## 3. Level-4 chart identity

- **Chart lemma (d) at k = 4: valid** [hand]. The proof in `B/fix1/b12_text/chart_lemma_body.tex` (d) uses only D_y = T along the curve, D_y^k y^p = p^k y^p and f∘Y = b + Φ. Nothing in it depends on k ≤ 3; the restriction came from `n1_chart.py` building k ≤ 3 only.
- **Structure: confirmed independently** [sym, `chi4_check.py`]. I rebuilt χ_k with jet symbols and a total-derivative operator, with sin δ and cos δ as algebraic symbols, so the construction shares no code with `D8/code/chart4_shell.py`. Results, all exact:
  - χ_4 has no b_5, b_6 and no k^(5);
  - ∂²χ_4/∂b_i∂b_j = 0 for (2,3), (2,4), (3,3), (3,4), (4,4), and ∂³χ_4/∂b_2³ = 0;
  - ∂²χ_4/∂b_2² is not identically zero (leading term 10·σ⁴·sin δ), so the quadratic term is real;
  - the coefficients depend on (σ, e, δ, k…k⁗, u, b, b_1) only;
  - χ_k = 0 on a ray for k ≤ 4; k ≤ 3 are affine in (b_2, b_3) without b_4.
- The assertion in `chart4_shell.py:98-118` tests the right things. Its envelope |k⁗| ≤ (m_b(σ + q_s) + k_b(1 + k_b·p_b))/4 follows from k‴ = −kq/4 [hand].
- Size: all chart-4 coefficients are ≤ 2.2e-35 on the shell [cert, `D8/out/chart4_shell.json`], consistent in size with my leading term 10·σ⁴·sin δ at |δ| of order 1e-44 [est]. Their effect on any number is nil.

## 4. The k = 2 export

**The export is valid as a certified computation** [read].
- Retention: `N36/t2_A/rigor/weights_t2.py:177-181` lowers the core B_c by r_2 plus a 1e-6 cushion; `:241` imposes the absolute amount r_2·a_c on the zones [40, 60], the segments and the transition up to 350; `:284-285` asserts r_2 ≤ the certified core retained fraction. Cost: Ĝ_22 ×1.00053.
- Window form: `assembly_rep.py:230-240` computes C0_2 = (Ê_2 + D_2)/(r_2·a_c) and C1_2 = (2F_2²/γ′_2 + coup_2 + o_2 + D_2^slope)/(r_2·a_c), and asserts D^c ≤ D. Every level-2 source is a `Tr` (`:100-117`) built by the same three rules as at k = 0, 1.
- Integration: E_2′ + γ′E_2 + r_2·a_c·X_2′ ≤ sources + 2F√E_2, and 2F√E ≤ γ′E + F²/γ′; the code charges 2F²/γ′, which is conservative.
- The flux lemma reads the 4·N2_ref setting (116.378, 0.260600), which is the bootstrap level of the outer data. Correct choice.
- Weight reading: E^A_j = ∫ω·η_A²·(z^A_j)² ds with ω = 28/y (`lemmaS_A_v33.tex:39, 45`), so ∫_{Y(S)} z_k² ds ≤ (y_max/28)·E_k on the core. The implementer's [hand] reading is right.

**`lemmaD_k2_sketch.tex` states the right lemma; it is a sketch.** Missing for a written proof:
1. The level-2 energy inequality with the retained term. `lemmaS_A_v33.tex:106, 199, 269` still say r_2 = 0, and `:302` says "Level 2 exports nothing".
2. Condition (D2) at level 2, source by source. The level-2 source list differs from k = 0, 1 in one way: the top term goes by IBP with Young users against the ramp slack (`assembly_rep.py:164-166, 190-199`, TOP_MODE[2] = 'IBP'), and there are `classes_IBP` and the `ibp223` handler. Each needs its window form written: the Young costs are linear in V_m² after conversion (kind (i) of `lemmaS_A_v33.tex:300`); the absorbed part is pointwise in τ; `ibp223` is Φ times a box constant.
3. The window form of the cubic Φ with slope ρ_J·Φ (`trc`, `assembly_rep.py:111-114`): Hölder in τ for J_F^{3/4}·J_G^{1/4}, with a pointwise F_max.
4. The statement that E_2(τ_a) ≤ Ê_2 at every τ_a, and that the couplings use E_k ≤ Ê_k (the conclusion level, legitimate a posteriori on [τ_1, τ′]).
5. New numbers: Ĝ_22, the Ĝ22 margin 1.91640, B_c^core = 1688999/10⁶.
6. An independent review of the v36 chain as a whole (R3 patches plus r_2); none exists yet.

## 5. Departures from the design

| # | Departure | Finding |
|---|---|---|
| 1 | W(57) = 1.3403e-18 | **Valid.** The verifier's constants come from `agent_sob/sob_constants_P2.json` (`v4cert/c5_constants.py:14`): W = (1 + x²)⁻² on x ≤ 30.25, S7 blend on [30.25, 38.25] (`agent_sob/sob_constants_P2.py:172, 192`). My own validated integral gives W(57) ≤ 1.340304797627754e-18 [cert], identical to `flux2_coef.py:322` to 16 digits. The v3.1 value 1.43716e-18 belonged to a quintic blend that the verifier never used. |
| 2 | OR sups read from `v4cert/out/c2_t3d.json`, sub-region [57, 60] | **Valid.** The sub-region is produced by the pinned C2 run itself (`c2_t3d.py:149`) with boxes \|w\| ≤ 0.5 and \|w′\| ≤ 0.02 ≥ S0V. Each value is asserted ≤ the ZB value (`flux2_coef.py:71-72`). |
| 3 | Minkowski in the shell conversions | **Valid and tighter.** I checked each: M_l² = ∫σ^{1−2l}(Σ s(l,m)D^m g)² ds; the coefficients (1), (1,1), (1,3,2), (1,6,11,6); the factor σ_min⁻¹ in ‖g″²‖; the change of variable with 1 + ρ; the absorption of B̄_22, B̄_33, B̄_44; and `square()` (`flux2_data.py:101-110`), which is Cauchy–Schwarz with Σθ_i ≤ 1. My direct, non-affine evaluation at seven values of T lies below the affine bounds every time (§6). |
| 4 | Typed ā = 1e-6 | **Valid**, conservative by ×4.6, sensitivity 1e-6 relative. It must appear in the bootstrap set (T3). |

## 6. Spot recomputation [cert, `recompute.py`]

Inputs taken as given: Ê and the exports from `N36/runs/prod2/t2_A/asm_design_N2ref4.json`; the chart and d data from `D8/out/chart4_shell.json` (sizes 1e-36 to 1e-4); the OR coefficient sups from `c2_t3d.json`; κ⁽⁰⁾ only for the comparison row. Recomputed with my own code: W(57), the sups of ζ^(i), the D2v coefficients, ϰ_2, the level-0 pair coefficients, N̄, the X2 table, K̄0, K2^pt, the window norms, C′, S′, and C0–C6 on the rate tables.

| Quantity | Mine | Reported | Note |
|---|---|---|---|
| W(57) | 1.340304797627754e-18 | same | identical |
| N̄_0, N̄_1, N̄_2 | 1.99251e-11, 3.25769e-12, 1.46196e-12 | same | agree to 15 digits |
| z_1…z_4 | 0.729168, 0.834803, 1.944448, 10.370373 | 0.729171, 0.834815, 1.944460, 10.370379 | reported ≥ mine; exact values 35/48, 280/27 |
| 𝔞_0…𝔞_3 | 10.3839, 10.4120, 8.9701, 4.4323 | 10.3974, 10.4181, 8.9734, 4.4332 | reported ≥ mine |
| K2^pt | 1.38803e-19 | 1.39001e-19 | ratio 0.9986 |
| C′ | 2.65224e-26 | 2.65232e-26 | from the reported affine window constants |
| S′ | 1.05496e-28 | 1.05510e-28 | |
| K̄0, κ form with my N̄ | 3.409853e-21 | 3.409853e-21 | reproduces |
| **K̄0, direct pair form** | **1.00289e-21** | — | valid bound by a different reduction; see below |
| Combined margin, κ form | ×4.8044 | ×4.8044 | reproduces |
| Combined margin, direct form | ×16.29 | — | |
| C2 | ×5.34e7 / ×5.00e7 | same | α = 1.5e-4 / 1.25e-4 |
| C3 | ×1547.25 / ×2235.83 | ×1547 / ×2236 | I allowed for the float rounding of rr |
| C4 | ×1.848e5 / ×1.736e5 | same | |
| C5 | ×1554.8 / ×1512.4 | same | |
| C6 | ×141.3 | ×141 | |

- **K̄0 is conservative by ×3.4.** I reduced level 0 by the design's hand form (pairs (0,0), (1,0), (2,0), (1,1); §4 of the design) and summed Γ_ab·N̄_a·N̄_b directly: 1.003e-21. The implementation's own Γ (F1 form) summed directly gives about 9.9e-22 [est from its four Γ values], so two different reductions agree within 1%. The reported 3.41e-21 comes from the fixed AM–GM weights (3/4, 3/4, 1/2), which suit badly because N̄_0 ≫ N̄_1, N̄_2. This is not an error. The slack would be s ≥ 1.2316e-21 instead of 1.0391e-21.
- **Windows.** At T = 0, 1e-3, 1, 100, 1e4, 1e6, 1e9 my direct Minkowski values of ∫M_3², ∫M_4², ∫s_2² and of ∫φ are all below the affine bounds; at T = 0 the ratios are 0.983 (m3) and 0.99993 (m4), so the affine constants are tight where it matters.
- Not recomputed by me: the chart and d integrals of X3 beyond the structure check, the profile envelopes, the 35 [sym] checks of X1 as such (covered indirectly by the agreement of the two level-0 reductions), and the OR coefficient sups.

## 7. What still stands between the current state and "step 3 closed"

**A. Flux lemma v2 itself**
1. **U-1, the regularity lemma (C8).** Needed: Λ absolutely continuous; the level-0 and level-2 identities for v = ζg with a C³ cutoff whose fourth derivative jumps at 16.25 and 57, and with the piecewise-C³ weight; g⁗ ∈ L²(OR × I); Agmon for g″ a.e.; local integrability of εω(z^A_{k+2})² past τ′ (used for (BΦ) being open); the same for (A)'s level-2 identity. All of it follows from smoothness of the flow for τ < τ_max, but it has to be written.
2. **G1, the written Lemma D for k = 2** (§4, items 1–5), and the corresponding edits of the (A) text.
3. **Text items** (none changes a number):
   - T1. `flux2_lemma_v2.tex:182` still says "Status. Design. Open before this can be called proved: the certified constants …"; `:70` still says "To be asserted in sympy"; `:3` describes placeholders; the data lemma's proof shows the squared-sum conversions, while the implemented Minkowski form is only in §impl.
   - T2. The march first-touch argument exists only as the docstring of `closure.py:2-8` and the sketch in Proposition "abs". TAIL_LEMMA.tex proves the tail only. The Y version of §1.2 has to be written as a proof.
   - T3. The bootstrap set must list (BΦ), ā = 1e-6 and the boxes B-a to B-e with their strict improvements; `BPhi_continuation_v36.tex` is one paragraph and no continuation text carries it.
   - T4. The hypothesis of the E_I lemma must be the combined (I3′). K2^pt alone is 1.74 × FLUX2, so any sentence "K2 ≤ FLUX2" is false for the new K^OR. V36_REPORT's "FLUX2/K2 6.0186" and "FLUX0/K0 3.4005" come from the old flux block (old κ, old W(57)) and should be marked as records, superseded by `flux2_cert.json`.
   - B1 → Λ ≤ L + 2C′; tail floor quotes with (1 + 2C′/L_F)^{1/2}, relative 2.7e-8.
4. **Pinning.** `D8/out/*` is in no manifest and no kit.

**B. Verifier fin13v4b** (from `V4B_REPORT.md`, unchanged by this work)
5. `ei_lemma_v4.tex`, the written E_I lemma with hypothesis (I3′), does not exist.
6. No independent review of the v4b delta (E0 F-B, E2 fix, R-2, R-3) is on record; the v4 review predates it.
7. Not audited: the θ-certificate, `cubic_E0fine` and ADD_STRIP_CUBIC, the Q2 Gram, the Gårding identity behind μ_2, the level-2 commutator constants; the float check of the β decomposition.

**C. (A) chain v36, Lemma U, zone III**
8. Independent review of v36 (item 6 of §4).
9. Texts: `lemmaS_A_v33.tex` step 6 (R3, r_2), `numbers_A_v33.tex`, `lemmaU_v34.tex` Lemma dU (c_dU 1.46343e-8), `continuation_v33.tex`, `main_v3.tex`, erratum, index.

**D. Step 4 and kits** (needs Micah's go)
10. Re-pin stages 8, 9, 18, 19, 20; a new final_tail driver that reads the EXCEPTION record; a new manifest.
11. Rebuild and run kit_step3 and kit_step4b with the v4b verifier, the v36 chain and the flux stages; independent Windows runs.

**E. Author items**
12. (P1′) a reference for local existence with uniqueness (the Lean CWFlow check, or Giga–Ito / Escher–Mayer–Simonett); (P2′) the sentence on modulation not breaking down before T_max; (P3′); Lemma O arc by arc; the single-point statement.

**F. Whole document**
13. LaTeX compile of every text (none has been compiled).
14. A final independent review of the assembled step 3 once A–C are done.

## 8. Minor remarks (no action required for soundness)
- C3: `closure.py:130` forms rr in double precision, so the pad can be short by one ulp of U (relative 2e-7 of the pad). Margin ×1547; immaterial.
- The (A)-chain flux block still carries the old W(57) (records only), as the implementer says.
- The fixed AM–GM weights at level 0 cost ×3.4 of margin (§6). If a larger pointwise margin is ever wanted, quote the pair form; nothing needs it now.

## Files
- `FLUX2_REVIEW.md`: this file.
- `recompute.py`, `recompute.json`: the arb recomputation (inputs pinned by sha256 in the `args` block).
- `chi4_check.py`, `chi4_check.json`: the exact structure check of χ_4.
