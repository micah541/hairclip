# Step 3 review (Nov 9, late) — gaps in STEP3_ARGUMENT.md

Scope: STEP3_ARGUMENT.md (blueprint), ASSEMBLY.md, CERT_PLAN.md (Nov 9 entries), verifier/ (VERIFIER.md, constants.py, rhs_iv.py, closure.py, logs), lemmaO_v3.tex, lemmaS.tex + lemmaS_assemble.py, MOD_LEMMA.md. step3_linear_certificate.tex was not re-read (§2 taken as given). No heavy runs; only hand arithmetic, one float scan, and one run of the existing lemmaS_assemble.py.

Caveat on freshness: files in writeup/lemmaS/ and verifier/ were being modified during the review (lemmaS_assemble.py and lemmaS_prop1_bounds.py at 15:29; rhs_iv.py and ival.py at 15:24, after the last complete logs at 15:16; verify_log_final_* contain 4 lines only). Numbers below are from the state I saw.

## Ranked findings

### 1. Lemma S(A) does not close with the numbers the blueprint quotes — FATAL as stated
- **Location.** Blueprint §6 "gives S0_v ≤ 0.013 and S1_v ≤ 0.22"; CERT_PLAN line 429; lemmaS.tex §3–4; lemmaS_assemble.py.
- **Why.**
  - The quoted values are reproduced exactly by E_j = (1/4)K_j·J_j alone: E0 = 3.1e-6, E1 = 2.17e-3, E2 = (1/4)(K1J1 + K2J2) = 0.242, giving S0 = 0.0129 and S1 = 0.2195. That omits (i) the factor κ = (60/24)(24/19)e^{C ε_in} = 3.2–4.3 of lemmaS.tex §3, (ii) the factor 2 in In_2 ≤ c(2K1J1 + 2K2J2), (iii) the hand-off data E_j(τ1), (iv) the ramp commutators.
  - With κ and the factor 2 only: S0_v ≈ 0.023–0.027 (above the 0.02 box for which T3 was certified) and S1_v ≈ 0.46–0.54.
  - With the machine-derived Prop. 1 constants (lemmaS_assemble.py as of 15:29, γ_w = 1, ramp [19,24], √Λ1 = 3.5e-7): E = (3.7e-4, 0.39, 1.4e3), **S0_v = 0.247, S1_v = 10.9**. At 2.5e-7: 0.19 and 8.3. The script was mid-edit (a second call raised KeyError), so treat these as indicative, but the mechanism is checkable by hand: the ramp-commutator constant is GhatR_0 = 2.66e7, so ε(19)·GhatR_0 ≈ 200. The y⁻⁴ smallness is eaten by the cutoff (ϑ = S⁴ on an s-interval of length 0.234 gives ϑ'⁴/ϑ³ ~ 1e6), and these terms multiply ramp norms (K_c·J_j), not E_j. The statement "insensitive to C_j ≤ 1e5" in lemmaS.tex does not cover them.
  - lemmaS.tex itself still says (Finding 4) that S1_v does not close from (I2) and asks for unweighted ramp integrals; the blueprint says it closes. The two documents disagree.
- **Consequence.** The far box |w′| ≤ 0.02 on x ∈ [34.25, 60] (see finding 2), the |V| ≤ 1/2 hypothesis feeding (B), and the level-2 entry data of (B) are all unsupplied.
- **Cheapest fix.** (a) Treat the ramp commutators at lower order directly with the E_I dissipation (they live on x ∈ [21.25, 26.25] where E_I has all derivatives up to v⁗ time-integrated) instead of Young against the shell dissipation; (b) widen the ramp (s-length 0.234 → ~0.7 gains ~80 in the η-derivative constants; K_j also drop at smaller y); (c) recertify T3 for a larger slope box (coefficients depend on w′ only through (1+p²) factors) rather than force S0_v ≤ 0.02; (d) state explicit thresholds for E_j(τ1) in hand-off (d).

### 2. Pointwise bounds on x ∈ [34.25, 60] used by the E_I estimate are never supplied — SERIOUS
- **Location.** verifier/rhs_iv.py lines 92–94, 129, 158 (sI1, sI2 = supv1, supv2 with KW = 2.323e6); closure.py `farsup`; VERIFIER.md "it is for x ≥ 16.25"; blueprint §5 "E_I-based sup|w′| on x ≥ 16.25: 1.40e-2".
- **Why.** KW = 2.323e6 is the Sobolev constant of the W-region [30, 34.25] (ASSEMBLY §0). Beyond 38.25 the weight is ρ-type: W(57) ≤ 1.44e-18 (lemmaS_cflux.json). So E_I gives no pointwise control on [34.25, 60], yet ζ = 1 up to 57 and the T3 terms d200, d110, d300 are applied there with supv1, supv2, and the 0.02 box is "checked" there with sI1. Same for S_all ("sup|w| over every region (0..60)").
- **What is needed and from where.** sup|w′|, sup|w″|, sup|w‴| on [34.25, 60], pointwise in time. Lemma S(A) can give the first two (S0_v, S1_v/y) if finding 1 is repaired; nothing gives sup|w‴| (levels 0–2 give only y²w‴ in L²). This is the Plan-B failure mode again in a smaller place.
- **Cheapest fix.** Split the far T3 terms at 34.25; beyond, use S0_v and S1_v/y as *coefficients* (keeps degree 2, so the tail lemma survives). Rewrite the I2 term d110·supv2·V2 (from w‴²) by a weighted Gagliardo–Nirenberg bound in V2, V3, V4 so no sup|w‴| is needed, or add level 3 to (A). For reference, the verifier already tolerates sup|w″| ≈ 0.19 worst case on the far region (my float scan of sI2 at √Λ = 6.39e-7), so S1_v up to ~6 would be tolerable there.

### 3. Modulation bounds in the verifier are not the ones MOD_LEMMA proves under the bootstrap boxes — SERIOUS
- **Location.** Blueprint §4 "|(a,b)| ≤ 7.57A (jet box 1e-4)", "|β| ≤ (0.00244 + 63.5|α|)A + 0.626‖χΦ_F‖"; constants.py CAB = 7.57, BETA0/1/F; closure.py `beta`; MOD_LEMMA.md statement and caveats.
- **Why.**
  - 7.57 needs |w|, |w′| ≤ 1e-4 on supp χ (x ≤ 32). The bootstrap only gives S_ramp ≤ 9.0e-3 and |w′| ≤ 1.4e-2 on [16.25, 32]. Under η = 1e-2 the lemma gives c_A = 34.7 (β: 2.70), and η = 1e-2 is itself violated on [16.25, 32] (see finding 7).
  - The β bound is the "linear corollary" with the 13.4·W1·A term and every nonlinear term dropped. The full statement has c_W1·W1, c_B·B, c_W3·W3, c_D4·D4 and S2-terms. Worst case along the trajectory W1 ~ 1e-5 and B ≤ √(Λ/ε) ~ 2e-4, so 33.2·W1 + 2.03·B ~ 8e-4 against 7.57A ~ 5e-6.
  - D4 and W3 are controlled only in L²(dτ). So a, b, β are **not pointwise bounded** by Λ as the comparison system assumes. Inserted as stated, β ⊃ 0.0471·D4 gives in E0 the term n3·0.0471·D4·A = 0.354·D4·A, which (ε/2)D4² with ε = 1e-5 cannot absorb (Young ratio ~4e5); in I0 it gives μ·0.0122·D4·V0, the same size as the "term (a)" that produced rate +5.6.
  - ASSEMBLY §3 TODO "β-loop contribution of ⟨χ𝒩, f3⟩" is still open; constants.py is the "BLOOP=0" configuration with MOD_ALPHA = 6.19 patched in.
- **Cheapest fix.** Integrate ⟨χ𝒩(w), f_j⟩ by parts (f_j smooth, ρ-weighted) so no D4, W3 appear; recompute A0, A1 as α-proportional (MOD_LEMMA caveat 2) at the real box (0.02 on [16.25, 32]); then put the resulting W1- and B-terms into rhs_iv.py as shape-dependent coefficients. Also needed downstream: ∫|a|, ∫|b| for Lemma S(B) and for λ(τ)e^{τ/4} → const must use the corrected bound (with the B-term, ∫|a| can reach ~0.3 at the certified tail rate, so exp(c_j∫|a|) is no longer ≈ 1).

### 4. Term (c) (E_I outer-ramp flux) cannot have the form the verifier assumes; Λ → 0 is not yet proved — SERIOUS
- **Location.** Blueprint §4(c) "CFLUX … expected ≤ 0" [O], §5 "Λ → 0", §9; lemmaS_cflux.json; VERIFIER.md Follow-up 2.
- **Why.** The verifier needs CFLUX·(V0² + V1² + V2²) with CFLUX ≤ 3e-2 in the ζ-cut norms. lemmaS_cflux.json gives 10.9 (I0) and 113.7 (I2), and only in *uncut* ramp norms; the alternative is an absolute constant 6.1e-18 / 1.5e-17 (8e-21 / 2.4e-20 "sharpened"). A cutoff flux is never bounded by the cut norms. With an absolute constant F0 the tail lemma returns a floor Λ_F ~ μF0/ρ (√Λ_F ~ 1e-8 to 1e-9), so no convergence to S1, and J_j become affine in T. Also the I2 flux needs v‴ = w⁗ on [57, 60], which Lemma S(A) (levels ≤ 2) does not give pointwise in time.
- **Cheapest fix.** Treat (c) as a forcing g(τ) = W(57)·C·(ramp quantities from Lemma S(A))², not as part of r. Show the (A) energies decay (they obey Ė_j ≤ −γ′E_j + cK_jV_j² + forcing, and V_j, β → 0), close the pair (Λ, E^A) as a 2×2 system (loop gain ~ K_j·W(57) ≪ 1), and take the level-3 ramp quantity time-integrated (enough for a forcing term). The tail lemma must be restated to allow an integrable forcing.

### 5. Lemma O near the cone: the quoted time margin ignores the slope hypothesis; overlap with S(B) is unquantified — SERIOUS
- **Location.** Blueprint §6 "X = 3 gives t* ≥ 3.70e-7 against (3/y_O)⁴ ≤ 1.6e-7 at y_O = 150", "y_O ≈ 150–200"; lemmaO_v3.tex H2, H6, H7.
- **Why.**
  - The ratio (3/y_O)⁴ means cylinder half-width X·L = y_Oλ1 = r: the cylinder centred at y_O reaches the tip, where H1/H3 (graph, |h″| ≤ 0.12 in cylinder units) fail.
  - H2 "|h_0,x| ≤ 0.12 suffices" is not met by the bent cone. With u0′ = α·1.4·y^{1.0246}, α1 = 1e-3: slopes 0, 0.238, 0.483 at y = 0, 150, 300, i.e. |h_x| up to 0.238 in the tangent frame at 150, against the hard limit 1/4 − C1√t = 0.2437. At y_O = 200 with half-width r: 0.274 > 1/4, fails.
  - In general, with half-width c·r: time needs c·y_O ≥ 122 (X = 3), slope needs c·y_O ≲ 150 before any w′ allowance. The two margins trade off (e.g. c·y_O = 135: time margin 1.5, slope 0.21 vs 0.244). "Margin 2.3" exists only at zero slope margin. Also lemmaO_v3.tex says needed 1.2e-7 (margin 3.08), the blueprint 1.6e-7 (2.3).
  - H6: the inner end strip of the innermost cylinder needs K2 (pointwise h_xx) and K3 (L² h_xxx) for all t ∈ [t1, T). The tex still points to "the intermediate energy up to y = 260" (Plan B). In the blueprint this must come from S(B), which needs S1_f and E2^{f,w}, i.e. the level-2 inflow of (A) (finding 1), and (1−c)·y_O ≥ 48. Conversely S(B) needs Lemma O at its outer ramp (I3). The mutual dependence is fine in a continuity argument only if each side's conclusion region covers the other's hypothesis strip; no such region bookkeeping exists yet.
- **Cheapest fix.** Fix (y_O, c, R1) explicitly, e.g. y_O ≈ 180, c ≈ 0.73, R1 ≈ y_O; put the true H2 requirement (≤ 0.2437 − allowance) in the hand-off; state the K2 bound S(B) must deliver as "y|f′| ≤ 0.54 on the strip" (K2 = 0.18 at X = 3). Lowering α1 relaxes the slope side linearly.

### 6. Decay rate misquoted — MINOR for closure, matters for §9 and Lemma S(B)
- **Location.** Blueprint §5 "certified decay rate ≈ 0.011".
- **Why.** −0.0110 is the certified *lower* bound of the max rate at α = 0 (VERIFIER.md caveat 1, "lb only informational"), i.e. the decay cannot be faster than that in this comparison system. The certified tail gives ρ = 2.04e-3 for Λ (verify_log_abcdQ_L3.5e-7), so √Λ decays at ≈ 1.0e-3, and α is slaved (α ≤ κ√Λ, κ ≈ 120).
- **Consequence.** ∫A dτ ≲ 6.4e-7/1.0e-3 ≈ 6e-4, not 6e-5. Still fine for 7.57A, not with finding 3's terms.
- **Fix.** Quote ρ/2 = 1.0e-3 everywhere and recompute ∫|a|, ∫|b|.

### 7. Jet box 1e-2 is violated on [16.25, 32] for the inner coefficients — MINOR (cheap)
- **Location.** cubic_E0fine_0.001_0.01_0.01.json and nlin_E2_0.001_0.01_0.01.json (eta0 = eta1 = 0.01), region "16.25 < x ≤ 31.25" and far [12, 34.25]; closure log: far sup|w′| ≤ 1.398e-2 (1.083e-2 for the backup 2.5e-7).
- **Why.** Only the E_I T3 coefficients were recertified at 0.02. The inner E0 cubic, inner E2 far and MOD_LEMMA coefficients assume |w′| ≤ 1e-2 on that region; both hand-off levels exceed it. The blueprint §7 lists "inner cubic boxes" without saying which box holds where.
- **Fix.** Recertify those three at 0.02 on [16.25, 34.25], or require √Λ peak ≤ 4.5e-7.

### 8. Bootstrap and conclusion: items required but not listed — MINOR
- §7 does not list, among bootstrap hypotheses, the Lemma S outputs used as inputs of the closure (finding 2) nor |V| ≤ 1/2 on the (A) outer ramp (y > 60), which comes from (B). Order of derivation in §7 ("derive closure, Lemma S, Lemma O") hides that closure ⇄ S(A) ⇄ S(B) ⇄ O is a ring; each arrow needs a strict-improvement statement at τ*.
- Hand-off §8 lacks: numerical thresholds for (d); an upper bound for T − t1 in terms of λ1 (T − t1 ≤ λ1⁴·e^{4∫|a|}, needed for H7); the choice of (y_O, c, R1); the S(B) entry data for windows on [48, R1].
- §9 is correct in outline given Λ → 0 and ∫(|a| + |b|)dτ < ∞: H²_ρ-convergence on x ≤ 32 gives C¹_loc convergence of the tip, λe^{τ/4} → const gives T < ∞ and the self-similar rate, |ẋ_c| ≲ λ|b| gives a single limit point, and for fixed physical r > 0 the bound y|f′| ≤ S1_f gives curvature ≤ S1_f/r up to T. It therefore depends on finding 4 (Λ → 0) and on S1_f (level 2 of (B)).
- I2 forcing in rhs_iv.py is 1·α³·D2v only; the α- and α²-parts of ‖ζΦ_F′‖_W (2.7e-8, 1.16e-7, CERT_PLAN line 321) are missing. Numerically negligible, but ASSEMBLY §6 tags this [F] while the blueprint says [C].
- H8 (distance to the other arm) is unnecessary for immersed CDF (each arc evolves locally) and should not be imposed where the limit S1 self-crosses.

### 9. Presentational
- Region bookkeeping mixes x and y = x − 2.249: S(A) "[24, 60]" and S(B) "[48, R]" are in y; E_I, χ, ζ are in x. The overlaps do hold (S(A) ramp x ∈ [21.25, 26.25] inside the polynomial-W zone; E_I outer ramp y ∈ [54.75, 57.75] inside the (A) core).
- "d|α|/dτ ≤ λ̃|α| + |β|max(α, Λ)" should read β(α_up, Λ_up).
- R(τ) = R1e^{(τ−τ1)/4} is not exactly physically fixed (differs by e^{∫a}); define the edge as r = const and the a-flux disappears from (B).
- lemmaO_v3.tex is titled "for Plan B"; lemmaS.tex §5–7 predates the verifier's J's.
- Rounding by 5e-16 inflation, tail lemma text (TAIL_LEMMA.tex.in), B&B leaf semantics: already listed as caveats; tail_cert.py does certify the degree structure, but for the QUAD=1 right-hand side that finding 4 invalidates.

## What looks solid
- The transport structure itself: slope equation V_τ = −(¼ − a)yV′ + … with all fourth-order terms carrying y⁻⁴ (symbolic residual 0); (B) having no forcing; the diagnosis that the modulation forcing in (A) is not y⁻⁴-suppressed but acts over a transit time 3.7 (size ~1e-5, negligible even with β's α-part).
- The time-integrated → pointwise step in (A) is legitimate: inflow enters as ∫In_j dτ ≤ cK_jJ_j into a damped inequality, so sup_τ E_j follows from J_j. This is not a repeat of the Plan B error.
- Inner closure architecture: every sup is taken at fixed time from (A, W1, B, V0, V2) with D4, V4 as free dissipated variables; monotone P-class arithmetic; radial/degree lemmas machine-checked by tail_cert.py; J0, J1, J2 from the θ-split.
- Margins of the comparison system itself: peak √Λ 6.39e-7 against cap 7.5e-6 (11.7×); S_ramp 9.0e-3 against 1e-2; term (a) with certified ramp coefficients has 5–10× room.
- Lemma O v3 constants (ball arithmetic, monotone in M, K) and its hypothesis list as a self-contained statement; using it once at t1 removes the per-point hand-over problem.
- Modulation algebra: Gram/Schur inversion, Neumann factor ν = 0.08, leading 3.708A, β-loop cancellation to 6.19|α|.
- §9 logic, conditional on findings 3–5.
