# STEP 4 REVIEW (independent referee, Nov 11)

Scope: Theorem 4 of `STEP4_RESULT.md` (status of Nov 11, 15:37 PT) and the chain Γ0 → P1 → switch → P2 → T9 → step-3 hypotheses (v3.1).
No existing file was changed. My scripts are in `/tmp/rev/` (`chk_u4.py`, `h4.py`).

**Verdict in one line.** Theorem 4 is not proved as stated. One step-3 hypothesis is not delivered and does not follow from the certified numbers (A1). Two parts of the P2 bootstrap are cited but not certified (A2, A3). I found nothing that makes the conclusion false, and I expect all three to close, A1 only after a change to the hand-off window or a new local estimate.

Severity: **(A)** the claim is unproven at this point (each entry says how likely it is to change the conclusion); **(B)** a needed fix that is fine in principle; **(C)** presentation.

---

## (A) Gaps

### A1. Lemma O hypothesis H4 (fifth derivative) is not checked, and fails on the certified numbers
- **What step 3 asks.** (H5) of `step3_proof.tex` l.97 requires "hypotheses H1–H7, H9 of Lemma O". Line 246 lists, among the items "to be checked on the step-4 data", H4: ∫(h₀⁽⁵⁾)² ≤ M² per cylinder (cylinder units). It is hypothesis (I) of `lemmaO.tex` l.26 and is used in its Step 3 (l.128, ‖ζ²a h₀⁽⁵⁾‖ ≤ M). `lemmaO_v3.tex` l.50 repeats it ("good-time selection").
- **What step 4 delivers.** Nothing. `t9.py` and `t9_far.py` bound K2, K3, K4, ‖h‴‖_L² and the slope only. The (H5) rows of STEP4_RESULT §2 and Theorem 4(b) do not mention it. `LEMMA_F.md` C5 (l.200) says "Lemma O (H4) already uses such a selection", i.e. T8 assumed the selection happens inside Lemma O. It does not: it is a hypothesis on the data at t₁.
- **It does not follow from what is certified.** The only control of e⁽⁵⁾ is D₂(τ*) = ∫T₂w₄², with T₂ = 28Y_ℓ². In cylinder units ∫(h₀⁽⁵⁾)²dX = L⁷∫(h⁽⁵⁾)²dx, so the error part is L⁷D₂/(28y₀²). With the certified D₂(τ*) ≤ 1.2236e-8 (`amp_n3_fin.json`, Ȳ = 1e-12), recomputed by `/tmp/rev/h4.py`:

| cylinder | L | error part | M² | over by |
|---|---|---|---|---|
| cone, r_c = 200 | 50 | 0.156 | 0.0144 | ×10.8 |
| cone, r_c = 290 | 72.5 | 0.959 | 0.0144 | ×67 |
| cone, r_c = 380 | 95 | 3.63 | 0.0144 | ×252 |
| far, centre 376 | 113 | 4.77 | 1 | ×4.8 |
| far, centre 602 | 113 | 0.74 | 1 | passes |

  The base part is negligible (≤ 6.3e-6, from the slab-180 hulls of κ̃‴).
- **Cause.** D₂(τ*) is a mean over W′, and W′ is only 1.46e-4 long. The bound is F₂/(θ₂|W′|) ≈ (1.36e-6)²/1.46e-4 ≈ 1.27e-8. Lowering Ȳ does not help: the residual part alone is D₂(0) = 2.9e-9.
- **Fix (one of).**
  1. Lengthen the good-time window to a whole slab or more. I estimate D₂ ≈ 1e-11, which would give ≈ 3e-3 at r_c = 380 (about ×5 under); this is an estimate, not a run. H1–H4 must then be certified over the longer window: α is fine (≈ 1.17e-4 at the slab start), but √Λ needs many more t_c-balls or a Lipschitz bound in τ.
  2. Prove a localised D₂ on the window y ∈ [L, 7L] instead of the global one.
  3. Ask step 3 whether Lemma O can carry a separate, larger constant for ‖h₀⁽⁵⁾‖.
- A smaller D₂ would also widen the K4 margins (×1.61 and ×1.81), which use √(F₃D₂).
- Also unscripted, but I expect them to hold: Lemma O H1 (graph height against H − C₀√t₁) and H6 (neighbourhood of centres).

### A2. Nonlinear terms are missing from the certified (Y) comparison
- **Evidence.**
  - `t7/closure_cert.py` certifies ϱ_c "in linear mode" (docstring l.1). Its term list has no N-terms.
  - `p2/ycomp.py` l.23 has a single allowance, `DN = 1e-6`, subtracted from the rate. Its comment points to a "P2 section" of STEP4_RESULT that no longer exists, and justifies the value by "actual box |e|, |e′| ≤ 1e-8", which is the tip box only.
  - T7's certified n₂–n₄ (tip and far, `t7/constants.json` → `nonlinear`) are read by nothing in the (Y) path. `merge.tip_lemma` uses only m₂₃.
  - T5 §2.5 is a sketch (unspecified constants C, "…" in the ε-level bound; N_I0, N_I2 never written out).
- **Why it matters.**
  - *Far part of (I0), (I2).* T5's E_I runs over the whole far curve [10.25, ℓ̃]. The top coefficient of the graph operator is perturbed there by Ã ≈ 4κ̃e, and T8's own certified number is sup Ã⁺ ≤ 0.0119 (`n3_fin.json`, `sup_Atp`). T8 charges this to every θ_j of the F_j levels (TOPLOSS), but nothing charges it to the ½V₂² and ½V₄² dissipation in (Y). DN = 1e-6 is four orders too small for it.
  - *Tip.* Terms such as η₁n₄D₄A in (E0) and εη₁n₄D₄² in (E2) contain D₄, which Λ_e does not bound. They cannot be absorbed by a rate loss; they need a slice of dissipation.
  - *No slice is left.* ycomp's `shares` (l.31) spend the whole half-reserve on residual Young terms, the other half belongs to Lemma J (κ_D = ½), and the only extra slice is `D4CUB = 3e-9`, taken by the cubic term.
- **Expected effect.** Small but not nil. A loss of about 1e-2 in ϱ_c raises the (Y) level from 1.33e-13 to about 1.7e-13, far below Ȳ = 1e-12. Taking the top loss out of Lemma J's share raises the J's by roughly 20–25%, so the level-3 far box (now ×1.76) and item 4 (×5.0) shrink by about 10%. These are my estimates.
- **Fix.** Write the tip and E_I nonlinear lemma (the analogue of N3′ for (E0), (E2), (I0), (I2), with 𝒬 and E_ij of (M)), with T7's n_j per x-range and T8's cancellation for the x-growing terms. Then rerun `closure_cert.py` with reduced dissipation coefficients and Gårding constants, and reallocate the shares.

### A3. The tip pointwise box is attributed in a circle and is not certified
- **Evidence.** The box |e|, |e′| ≤ 2e3·Ȳ/√μ_e on x ≤ 34.25 is used by:
  - T7 Lemma 4a ("verified box", `T7_REPORT.md` l.134, `merge.py` l.172);
  - T8 N3′ (`LEMMA_F.md` l.272: "T5's tip box, certified by the assembly"; `n3far.py` l.380);
  - the assembly (`clam.py` l.5: "T8's certified tip box"; STEP4_RESULT §1: "T8's certified box").

  Its only source is T5 §2.5, "≲ 2·10³Ȳ/√μ_e", tagged [D]. No script computes it (`boxes.py` l.32 only prints it).
- **The constant is wrong at the outer edge.** Λ_e ≤ Ȳ² alone allows:
  - |e′(x)| up to (4.777/(μ_eW(x)))^{1/2}Ȳ, the sharp constant for ∫v² + ε_I∫v″² on the line. That is 1.056e4·Ȳ at x = 30.25 (so "2e3/√μ_e" = 1.054e4 is exactly this value there) and **1.35e4·Ȳ at x = 34.25**, above even 1.25 × 1.054e4 = 1.32e4.
  - |e(30.25)| ≥ 1.16e4·Ȳ (take e′ ∝ 1/W on [16.25, 30.25]), against the claimed |e(28)| ≤ 1.054e4·Ȳ; and up to about 1.69e4·Ȳ at x = 34.25.
- **Consequences.** Numerically none at Ȳ = 1e-12: n3far was run with Y = 1.85e-12, so its η₁ = 2.44e-8 exceeds 1.69e-8. But:
  - the statement "the boxes close for every Ȳ ≤ 1.85e-12" is false for the tip part;
  - Lemma 4a's β scales with η₁, so the reserved slice D4CUB = 3.0e-9 is short (about 4.8e-9 needed);
  - `clam.py`'s tail term uses the false constant (its weight is ≈ 1e-6, so no effect).
- **Fix.** A short certified lemma giving sup|e|, sup|e′| on [0, 34.25] from (A, B, V₀, V₂) and K_e, with arb constants. Pass the result to n3far (`--eta1`), `merge.tip_lemma`, `clam.py` and `boxes.py`, and assert it against T7's η₁ = 1e-5.

---

## (B) Needed fixes, fine in principle

### B1. `se_corr.py`: underived majorants now dominate the far residual norms
- `out/se_corr_v2.log`: the (SE) increments are up to ×124 (R1), ×157 (RF1), ×190 (RF2), ×309 (RF3) of the T6 values. So R1 in (Y) and the residual parts B_j of √F_j are set by this script, not by T6.
- The formulas contain unexplained constants: `quad = 20·dsum²(1+Ks)³` (l.102), `4^j·3·Ks²` (l.105), `8(…)²` (l.131), `10·Lh²(n+4)·Yn_rad_max` (l.67). There is no derivation: STEP4_RESULT cites "(SE) §2(ii)", which exists only in `out/STEP4_RESULT_v1_stopped.md` as a different proposal (static shift in physical arclength, s ∈ [6, 12]); the implemented correction is a comoving bump on x/ℓ̃ ∈ [0.6, 0.8]. Existence of the root a₂ is also not argued.
- It matters because the thinnest bootstrap margin (B2) is almost entirely this residual part.
- **Fix.** Write the correction as a lemma with its derivation, or compute the corrected residual directly with T6's engine on the corrected base.

### B2. The far-box margins are reported against the wrong boxes
- N3′ was certified on |e′| ≤ 7e-4, x|e″| ≤ 1.65e-4, x²|e‴| ≤ 1e-2 (`n3_fin.json` args: b0, b1, b2). `boxes.py` l.27–30 compares all three with η_F = 1e-2.
- Rerun at Ȳ = 1e-12: 5.421e-4, 7.532e-5, 5.689e-3. Real margins: **×1.29, ×2.19, ×1.76**, not ×18, ×132, ×1.85.
- The thinnest margin of the whole bootstrap is therefore |e′| at ×1.29. It hardly depends on Ȳ (5.47e-4 at 1.85e-12): it is B₀ = 1.79e-6, the residual part (B1).
- No script asserts the closed values against the boxes n3far assumed, nor that √F₀ with N3′ stays inside the factor 1.25 of the |e| box.
- **Fix.** Let boxes.py read the boxes from `n3_fin.json` and assert all four. Consider rerunning n3far with a larger b0.

### B3. H1: the step-3 IFT is a Newton step with a typed remainder
- `h1.py` l.126–129: `solve` adds a ball Q = 1e-11 / 1e-13 to G for the quadratic remainder. Nothing bounds the second derivative, and no Krawczyk or Brouwer test proves that the three constraints have an exact zero in the reported box. (H1) asks for the constraints to hold exactly.
- By my estimate the remainder is ≈ 1e-17, so the result will stand.
- `fU4_safe` (l.122–124) returns a typed bound 50 on any exception. I counted on ball K = 31: 0 exceptions in 6166 calls, max |u⁗| = 0.107. Make it a hard error.
- Typed bounds in the rotation-mode norm (l.197: 1.01, 1e-3, 1.6·60^1.03, factor 50) and the first-order factor 1.01 (l.204). They enter as 1.7e-12 × 17.6, so they are harmless.
- `clam.py`'s factor 1.01 is justified by "|u′|, |ku| ≤ 1e-6", which is false (u′ ≈ α ψ̂₃′ ≈ 1.5e-4). The factor is still ample.
- H1 itself looks sound: √Λ ≤ 2.291e-8 over 32 balls, C_ΛȲ = 2.4e-11.

### B4. Typed inputs on the wrong side (the step-3 v3.1 error class)
- `camp.py` l.34, `BASE_T1_CERT`, against `out/t9_win_tight.json`:
  - E3 4.696e-7 < 4.69602e-7;
  - shell E0 2.832e-12 < 2.83221e-12;
  - shell E1 3.659e-12 < 3.65948e-12;
  - K4 3.163e-4 against ≈ 3.1633e-4 (my recomputation from the r_c ∈ [200, 210] row).
- Effect ≤ 1e-4 relative, against margins ≥ 17.
- **Fix.** Read the JSON outputs instead of typing them.

### B5. Typed allowances in T9
- `t9.py` l.61–62 (2e-9, 5e-11, 1e-13, 1.001), l.66 and l.114 (1.0001), l.116 (c0: 0.002 and "+0.01 per level, as used since v2").
- `t9_far.py` l.32–34: (1.003)^{j+1}. The e″ term carries (1 − κ̃e)⁻² ≤ 1.0048, which is above 1.003; harmless, since e″ ≪ κ̃ there.
- I checked the three additive ones by hand; they hold with room. The c0 allowance and the factors have no derivation. Margins: E3 ×17.9, far jets ≥ ×14.7, cone K4 ×1.61.
- **Fix.** Replace by computed arb expressions, or state each as a short lemma.

### B6. T8's hypothesis |x| ≤ C_xȲ omits the residual and nonlinear parts of (M)
- `camp.py` l.124 and n3far use C_x = 7.82, the linear constant. (M) adds x_r = m̄c_φℛ₂ (≤ 3.3e-14 at slab 2, about 0.4% of C_xȲ) and nonlinear terms.
- **Fix.** Use the constant majorant Y ≤ 6.3e-13, which ycomp's numbers already support (level 1.33e-13, Y(τ_a) ≤ 6.20e-13), or add x_r explicitly.

### B7. Small unscripted steps
- ycomp is run without `Y0` (`ycomp_fin.log`: "Y0 0.000e+00"), so Ȳ ≥ Y(τ_a) is not machine-checked. It holds: 6.20e-13 < 1e-12.
- camp places F_j(τ₀) at the start of slab 2; the switch is 94% through that slab. The decay of the initial data is over-counted by about 1.7% in √F. The effect is negligible, since the F₀ terms are tiny.
- Theorem 4(a): t₁ ∈ [4.5978580262, 4.5978580263] and k_max ∈ [2.36e4, 2.38e4] are tagged certified, but no script produces them. They need ∫|x₁| and the clock shift ∫ϑ; both are tiny (ϑ̄ ≤ 6.2e-8).
- Step-3 (H6) and the choice of S1 are not stated in Theorem 4.

### B8. S1 non-uniqueness: the argument survives, with two additions
- Every object I traced is an interval enclosure valid for every zero of M in Y: prof_tip (conditioned, mean-value form over the convex Y), ext813 (unconditioned box), the 3425 IBP boxes (hulls from `outA5.txt`), and `c4fun`'s ψ̂₃, ĥ₂ (fixed functions). T9 multiplies hulls, so its bounds hold for each zero. Distinct zeros differ by ≤ 3e-12, far below the step-3 floor 4.9e-10, so there is no contradiction.
- Additions:
  1. State Theorem 4 for "a shrinker S1 in the certified set", or prove uniqueness (ext_tip's enclosure A₀ = −4850 over all of Y is most of it).
  2. Add a scripted check that the step-2/3 hulls (`outA5.txt`) contain prof_tip. The READMEs only prove nesting in prof_tight. I did not verify this inclusion.

### B9. Rounding spots
- `ycomp.py` l.12: `up()` is round-to-nearest, despite its comment (output only).
- `switch.py` l.295, l.338: `arb(float(a.lower()))` is not outward (≤ 1 ulp in a weight argument).
- `h1.py` l.166: cells are tiled by float accumulation `x += h`, so ulp-size gaps or overlaps are possible.
- T6 `slab.py` l.304–306: scalars and ranges are stored round-to-nearest.
- All are far below every margin, but they are not outward.

---

## (C) Presentation
- STEP4_RESULT cites "(SE) §2(ii)" and (through ycomp) a "P2 section"; neither is in the file. The pairing of the residual terms in ycomp is documented only in code comments. I rechecked the constants: sup|ρ′/ρ| = 1.082 ≤ 1.083; sup|ρ″/ρ| ≈ 1.16 ≤ 1.3; sup|χ′| = 0.729; sup|χ″| = 0.835.
- Theorem T5 as written is no longer the theorem used: φ₁, φ₂ are T7's late-base functions, level 4 is dropped, the box differs, and the base is (SE)-corrected. A consolidated statement is needed.
- Stale text: `T2_REPORT.md` describes the unslaved curve (closure width 4e-16; the final value 8.1e-24 is only in `final.log`); `T8_REPORT.md` §1–5 and `LEMMA_F.md` §6.4 still speak of slabs 5–187 and placeholders; the `h1.py` docstring says K = 0..7 and prof_tight; `T7_REPORT.md` gives ℓ̃ = 6.7e5 at slab 180 where the jets give 3.7e5.
- STEP4_RESULT §1 lists the thinnest margins as ×1.85 and ×1.61; see B2.

---

## Checked and found in order
- **Γ0 and P1 start.** M0 = 1000 and N_c = 4000 ≥ 3762, so every Γ0 mode is in the centre or the box. The b₂ interval enters w and the slaved rows correctly. L0 matches. `final.log` certifies closure and embeddedness for `gamma0_final.json`.
- **Galerkin-limit lemma** (T3_DESIGN §6): first exit, compactness, mild limit and uniqueness are complete as written.
- **Switch.** The map α = π − s·x, κ = −sθ_α, λ_P = K0/κ_tip, the ball norms with the n⁻⁶ tail, the length factor, and the a-priori ball feeding the IFT are consistent. F_j(τ_a) are rounded up into `--F0`.
- **Step-3 definitions in T9.** w = u − αψ̂₃ − α²ĥ₂; χ₃₂; ζ on [10.25, 16.25] and [57, 60]; ε, ε_I, μ = 0.08. W_step3 ≤ (1+x²)⁻² (I checked the blend: φ′ ≥ 4x/(1+x²) for x ≥ 6.4). Windows: ω = 28/y, U = 1.3056αy^1.0246, sliding cover of y ∈ [50, 253]. The 32 balls cover W′ exactly, and SUBW = 0.15/1024.
- **Lemma O jets.** The h‴ and h⁗ formulas are right; cylinder scalings K_j = L^{j−1}sup|h⁽ʲ⁾|; far thresholds use 8.8e-3 (v3.1).

## Not audited
`t3/slab.py` internals, T7's certificate engine (`cert7`, `c2r3`, `far7`), T6's engine, `cconst.py` and the generated `n3forms.py`, and the CAPD programs. For these I relied on the reports.

## Overall verdict
- **Γ0, P1 and the switch:** I found no defect.
- **P2:** the comparison argument is certified only in its linear part plus the far-field levels. A2 and A3 must be closed before "bootstrap closed" is true. I expect both to close, with the level-3 far box and the |e′| box (×1.29, B2) getting thinner.
- **Hand-off:** (H1)–(H4) are delivered, modulo B3–B5. **(H5) is not:** Lemma O's H4 is missing and fails by ×5 to ×250 on the current certified D₂ (A1). This is the one item that needs a design change.
- The status line "All hand-off items are certified on W′, and no step-4 gap is open" should be withdrawn until A1–A3 are closed.
