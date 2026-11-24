# Step 3, final gap review (of writeup/step3_proof/step3_proof.tex)

Read: step3_proof.tex §2–§6 in full; TAIL_LEMMA.tex (tail, floor, L¹ lemma); closure.py tail; tail_cert_log.txt; constants.py (MOD block); CERT_PLAN.md θ-form entry; lemmaS.tex (resolutions 7, 8, strip). Not re-read: INCONSISTENCIES.md, MOD_IBP.md, lemmaO_v3 revision, the Lemma S assembly scripts. No computations beyond hand arithmetic.

## Verdict

**Not closed as written, but close.** The bootstrap logic (§4) is sound and non-circular. Three statement-level defects must be fixed before the theorem is true as stated (findings 1–3); each is cheap and none touches the certified closure. After those, step 3 is closed modulo the final-rigor items (4–6) and polish.

## Ranked findings

### 1. Embeddedness is assumed and concluded, but S1 is non-embedded — SERIOUS (statement)
- **Location.** Theorem hypothesis (H6) "embedded"; (H5) "every other arm at distance ≥ 0.1"; (W4) "does not meet any other arm"; §5.3 "the curve stays embedded"; (P1) "smooth closed embedded data"; §6 item 8 "embeddedness of S1".
- **Why.** The project's premise is that S1 is a non-embedded shrinker. A C¹-small even normal graph over S1 (|w|, |w′| ≤ 2e-3 by W2) inherits its transversal self-crossing. So (H1) and (H6) cannot hold together and the theorem is vacuous as stated; the conclusion "stays embedded" is false. If instead S1 is embedded on the chart region and the self-crossing lies elsewhere, that must be stated; I could not confirm it from the proof text.
- **Fix.** Replace "embedded" by "immersed" throughout. (P1) holds for immersed closed curves (Dziuk–Kuwert–Schätzle; the continuation criterion is ‖κ‖_{L²}, which E, E_I, the shells and Lemma O do bound). Read Lemma O per arc of the parametrisation; drop H8 and the "other arm" clauses, or keep them only as a separate, optional embeddedness statement away from the crossing.

### 2. Hand-off (H3)/(H4) are inconsistent with the bent-cone profile — SERIOUS (interface with step 4)
- **Location.** (H4) "y|f′| ≤ 0.12 on y ∈ [48, 350]"; (H3) "E1, E2, E3 ≤ 3e-5 on the strip [50, 150] and windows [60, 350]"; §3.7 "full slope f = U + v".
- **Why.** With U = 1.3056·α1·y^{1.0246}: y|U′| = 0.113 at 150, 0.120 at 158, 0.133 at 175, 0.270 at 350. So a curve on the slow manifold violates (H4) for y ≳ 158. And E1 of U alone on [50, 150] is 5.6e-3, 190× above 3e-5, so (H3) is unattainable if it means the full slope; if it means the deviation v, then K2 ≤ 0.080 must add the U-part (0.037 in cylinder units) and the proof must say so. The strip value "y|f′| ≤ 0.168" suggests deviation + U, but the text does not define it.
- **Fix.** State (H3) for v explicitly; state (H4) as y|v′| ≤ c, or as y|f′| ≤ 0.12 only on [48, 150] plus the Lemma O tangent-frame slope beyond; recompute K2, K3 as U-part + deviation and re-tabulate the margin (0.037 + 0.08 = 0.117 against M = 0.12 would be a 1.03× margin, so this needs an actual number).

### 3. The E0 α-form constant is an assumption, not a certificate — SERIOUS-to-MINOR
- **Location.** constants.py line 72: `MOD_ALPHA = 6.19 * 66.5 / 63.5  # ... (assumption)`; §3.4 caveat (i); §6 item 3.
- **Why.** 6.19 is the result of a cancellation between Q1 ([−13.47, 11.02]) and the β-loop form ([−11.16, 9.15]) computed with the *old* modulation functional. Rescaling by a ratio of upper bounds of β1 does not preserve a cancellation; if the IBP lemma changed the α-linear functional, the true constant could be anywhere up to ≈ 25. Not "a conservative reading".
- **Size.** At α1 = 5e-4: 6.5α = 3.2e-3 against 1/32 − Π3 = 0.0152. With 25α = 0.0123 the early rate grows but stays below coercivity; closure would probably survive with a higher peak.
- **Fix.** Either recompute the rank form with the MOD_IBP α-linear functional (small), or rerun the verifier once with MOD_ALPHA = 25 and quote that run (one run; removes the question).

### 4. Floor: the L¹_loc lemma is not what closure.tail implements — MINOR (real)
- **Location.** TAIL_LEMMA.tex, Lemma "L¹_loc forcing"; §3.5 "L_F ≥ Z_max/η"; §6 item 2(a).
- **Why.** The constant-forcing floor L_F = 2F0/ρ (implemented) is correct for K0, K2, which are pointwise bounds. The time-integrated piece I3 needs the L¹ lemma, whose anchor cell must be Λ ∈ [L*, (1+η)L*] with L_F ≥ Z_max/η; closure.tail certifies the anchor on [L*, (1+1e-9)L*], i.e. η = 1e-9, for which the lemma gives nothing. So the reported 5.49e-10 does not contain the (1+η) factor.
- **Size.** Z_max ~ I3-scale/ρ ~ 1e-23 or below; with η = 0.01 the extra floor is ~1e-21 against L_F = 3e-19. Negligible.
- **Fix.** Certify the anchor cell on [L*, 1.01L*] (H/s² is monotone, so the rate moves in the 5th digit), add Z_max/η, restate the floor as (1+η)(L_F + Z_max/η).

### 5. Ramp M1 allowance 0.01 for α²ĥ2′ is not derived — MINOR
- **Location.** D2; lemmaS.tex resolution (7); §6 item 4.
- **Why.** Needs sup y·α²|ĥ2″| ≤ 0.01 on y ∈ [14, 28]. Heuristically ĥ2″ ~ c2·y^{1.05}, so the term is ≈ 2.3e-4·c2 and the box 0.8 (0.67 + this) tolerates c2 up to ~560. Almost certainly fine, but it is a number in a bootstrap input.
- **Fix.** One interval sup of y|ĥ2″| from the existing C⁴ jets (minutes).

### 6. Lemma S arithmetic and Lemma O statement — MINOR (final-rigor items)
- Lemma S assembly in doubles without outward rounding, sympy-float class weights, the neglected sliver: genuine rigor items, as listed in §6 item 1. With minimum margin 1.83 rounding cannot change the outcome, but the tag [C] is not yet earned.
- lemmaO_v3 is not self-contained; B4 quotes K = 5.535 for the far cylinders while §3.8 and (H5) use M = 1, K = 1.5. One of the two is stale.
- (B): the derived value of |f| ≤ 1/2 is not tabulated (U(175) = 0.13 plus S0, so it is ≈ 0.15; state it).

### 7. Cosmetic
- **"3/64 vs ω = 1/32" (§6 item 7): not a gap.** The θ-form is its own certified theorem: Re⟨Lv,v⟩ + (1/9)‖v″‖² ≤ −(1/32)‖v‖² on the same constraint set (CERT_PLAN.md line 154, `certificate_step3_theta_T1.json`, pipeline with FORM_MULT = 9/8, OMEGA_Q = 1/32). It is not derived from the 3/64 statement and need not be. Cite it in §3.1 and delete the caveat.
- S0V passed by environment variable and not echoed in logs: make it a logged constant.
- (P1)–(P3) need references; k(0) ≠ 0 is a certified property of S1 from step 1 and should be cited as such.

## Focus questions

1. **Bootstrap (§4).** Closes as written. Hypotheses are B1–B4 only; D1 uses B1, B2; D2 uses D1 only (s_I2 depends on V0, V2, so it is pointwise in Λ); D3 uses D1, D2, B2, B3; D4 uses D3, hand-off, B4; D5 uses hand-off, D4; D6 uses D1. Every hypothesis is re-derived strictly (smallest factors 1.74 and 1.83; the 1.18 on the ramp is a derived-vs-box check, not a hypothesis). No time integral is used as a coefficient. First touch is correct for the march (strict via r + 1e-9, cell contains the touching state) and for the constant-forcing floor. The only defect is finding 4.
2. **§6 caveats.** Real: item 3 (E0 α-term, finding 3), item 2(a) (finding 4), item 4's allowance (finding 5), item 1 (rigor), items 5–6 (hand-off, sharpened in finding 2). Cosmetic: item 7, item 2(b), item 8 apart from embeddedness, which is finding 1.
3. **tail_cert.py.** The current log ends "TAIL_CERT PASS". The split is the right one: scale invariance is required of tot_s only (floors off: H and L both of joint valuation ≥ 2, certified), and the degree < 2 terms sit in L with hi(L) ≤ 2, which is what the Λ-cell rule needs, and enter the tail only through the floor. The earlier "FAIL" tested joint homogeneity of the full right-hand side, which is false by design in version W. The tail logic is sound, subject to finding 4.
4. **Hand-off H1–H6.** Not complete or consistent yet: (H3)/(H4) as in finding 2; (H6) as in finding 1. Missing: whether (H2)/(H3) energies are of v or f and with which weight normalisation on the windows; the far-cylinder family (centres, radii) as an explicit list step 4 can test; an explicit statement that |w|, |w′| boxes at τ1 follow from (H1). Nothing else is required by a component and absent.

## What looks solid
- The dependency order and strictness table of §4.
- The verifier: outward rounding, sound leaves, unconditional J's, certified structure, floor for pointwise forcing.
- Lemma S design (order-j commutators, two-zone M1, M0/M1 only) with margins ≥ 1.83.
- Modulation by IBP without D4/W3, and the resulting sup|a|, sup|b| ≤ 3.4e-7, which makes §5 (λ → 0, T finite, single point, rate) immediate.
- The θ-form coercivity certificate.
