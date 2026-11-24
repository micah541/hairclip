# U2 trace: does step 3 use Theorem B (D4)? Audit of the inherited shell-lemma items (R3); Evans multiplicity (R4)

Unit U2 of `OUTLINE.md` (§6 D4, §7.1 R3/R4, §7.2). Nov 16, 2026, Opus. Read-only: nothing in the proof tree was edited or rerun. The only computations were two small arb/float checks (nice 19, under 100 MB, reported below as [est] or [cert-check]).
Labels: **[read]** = verified by reading the cited file:line; **[cert-check]** = my own arb evaluation from step-1 data; **[est]** = my estimate (float model on certified inputs); **[not verified]** = I did not check it.

## Verdicts

- **A (D4): CONFIRMED. Step 3 v3.3 does not use Theorem B**, either directly or through any cited constant. Its linear input is its own constrained-coercivity certificate: a Morse count by a Plücker/Sturm zero count, a Gram matrix, and positivity at infinity, computed with step-1 profile data only. Its slow rate λ̃ = −0.2561532224577463 is a floating-point parameter fixed by Newton's method. Step 3 certifies residuals against that λ̃ but never certifies λ̃ as an eigenvalue, and step 2 does not certify it either: Theorem B covers only Re λ > −1/8. The verifier reads no step-2 output. Three things are wrong labels or code reuse, not dependencies:
  1. v3.2 text line 52 and `STEP3_ARGUMENT.md:17` tag the even spectrum "(step 2)" as certified.
  2. `main_v3.tex:118` says step 3 "rests on the linear theory of step 2", and the figure has an arrow from step 2 to step 3.
  3. Step 3 and step 4 import the step-2 module `lin/lintail2.py`, but only for the step-1 tail constants and the formula μ(s) = (s/4)^{1/3}.

  The one place step-2 constants appear is the generation lemma of the planning note `step3_blueprint.tex`, which uses step 2's numerical-range bound Λ₀ = 1.7909. That lemma feeds the semigroup decay corollary (‖e^{τ𝓛}‖ ≤ 14.8e^{−0.0213τ}), which v3.2 displays (line 131) but which nothing in v3.3 consumes.
- **B (R3): none of the three items breaks a margin, but two of them are not justified as written.**
  1. **cos θ constant along the ray: justified, and in fact exact.** In (A)'s ray chart, θ is the fixed angle between the translation direction (the symmetry axis) and the S₁ ray. Even if it is read as the varying S₁ tangent angle, the change to c_s is ≤ 7e-15 [cert-check].
  2. **Forcing constant r₀^F = 1.4(β̄ + 1.02α₁ā): unjustified as written but harmless.** It is derived under the identification λ̃ = −p/4, which is off by 3.22e-6 because p = 1.0246 exactly. That leaves out α(λ̃ + p/4)·1.3056y^p, which is 14.5% of β̄ at α₁. It also leaves out the fourth-order U term O(αy^{p−4}). The honest F_j are ×1.069, ×1.074, ×1.089 the certified ones [est]. The forcing carries at most 0.57% of Ê, so Ê moves by less than about 0.1% [est].
     - **Same root cause outside R3:** Lemma U's δU lemma assumes α_τ = −(p/4)α + β with the D1 bound β̄. Its c_δU is therefore about ×1.143 too small [est]. This is very likely harmless, but no script has checked it.
  3. **Flux lemma (`lemmaS_cflux.json`): a real rigor gap with huge numerical room.** There are four problems:
     - **(a)** The κ tables are computed for a **quintic** cutoff that is marked "ASSUMED". The verifier's E_I cutoff on [57, 60] is the **septic** C³ smoothstep S7.
     - **(b)** They rest on an **a-priori higher-jet box** (M₂…M₅ = 1.4, 5.8, 29.4, 150). Nothing in v3.3 establishes it.
     - **(c)** The input table `cflux_raw.json` is not on disk, so the κ values cannot be reproduced from the records.
     - **(d)** I did not find where the level-3 part I₃ is consumed by the verifier.

     The "linear with nonnegative coefficients" claim, by contrast, is verified by reading `f4b`. Every flux term carries W(57) ≤ 1.44e-18. The margins are ×3.37 and ×5.96 against the verifier inputs, and the verifier is insensitive to the flux level: the ×4 rerun moved only the flux-driven quantities. So this is a documentation-and-rerun item, not a threat to the result. It still has to be repaired or listed in the paper.
- **C (R4): moot for the Main Theorem, because A holds.** Theorem B as a standalone result needs the hypotheses listed in §4. `step2_proof.tex` verifies the ODE-level ones: an analytic admissible plane, a 2-dimensional decaying subspace, and zeros of D equal to eigenvalues. It does **not** verify the operator-level ones: L as a closed operator on L²(ρ) with isolated eigenvalues of finite algebraic multiplicity and index-0 Fredholm, and the equality of zero order with Jordan-chain length in this half-line weighted setting. `step2_proof.tex:149-151` says so itself. Only the sentence "algebraically simple" depends on that citation. "Simple zero of D_even, D_odd" is certified as it stands.

## 1. Trace A: every use of spectral information in step 3

Abbreviations: `M3` = `writeup/lemmaU/bulk/fix1/b12_text/step3_v33_main.tex`; `P32` = `writeup/step3_proof_v32/step3_proof_v32.tex`; `LC` = `writeup/step3_linear_certificate.tex` (the "fourth supplement", the linear certificate); `CP` = `step3nl/CERT_PLAN.md`; `S2` = `writeup/step2_proof.tex`. All paths are under `/home/claude/s1spec/`.

| # | Spectral item used by step 3 | Where it is used | Where it comes from | Source class |
|---|---|---|---|---|
| 1 | Coercivity Re⟨Lv,v⟩_ρ ≤ −(3/64)‖v‖² on even v ⊥ f₁,f₂,f₃ | v3.3 cites v3.2 §3.1 unchanged: M3:45; P32:120-131 | LC:60-66 (theorem), proof LC:132 from Lemma "positivity at infinity" (LC:76-84), constrained inertia (LC:86-88), Morse count (LC:96-108), Prop. "the count" (LC:112-117), Gram matrix (LC:119-129) | step-3 certificate (`step3/`); profile data from step 1 (LC:117, LC:153) |
| 2 | θ-form Re⟨Lv,v⟩ + (1/9)‖v″‖² ≤ −(1/32)‖v‖² (what the energy actually uses; constant PI3) | P32:133; verifier `step3nl_v33/verifier/constants.py:14` (`PI3 = 0.01607`) | `step3/certificate_step3_theta_T1.json`; CP:155-161 (positivity, tail theorem at λ = −1/32, Morse count 3, Gram Schur complement in [−4.93e-4, −1.73e-4]) | step-3 certificate |
| 3 | "Morse index of A = −H_sym − ω is 3" | LC:132, P32:129 | Plücker zero count of W₀₁ on [σ_tip, 32] (CAPD) plus the tail theorem `step3/hstail.py` on σ ≥ 32 (LC:116); Morse index theorem cited from Reid (LC:103) | step-3 certificate + **classical citation (Reid)**; not Evans |
| 4 | Count of unstable or slow modes (three constraints) | constraint set f₁,f₂,f₃ (P32:124) | The inertia lemma gives n(A) − n(G) = 3 − 3 = 0 negative directions (LC:86-88, 132). No statement that "L has exactly two unstable even eigenvalues" is used anywhere | step-3 certificate |
| 5 | Slow rate λ̃ = −0.2561532224577463 | decomposition α_τ = λ̃α + β (M3:17; P32:59; `STEP3_ARGUMENT.md:16`); verifier `LAMT = -0.2561532` (`verifier/constants.py:14`, `closure.py:3-6`) | **A chosen parameter**: "fixed in floating point once and for all … tuned by Newton's method so that F′(0) = F‴(0) = 0" (LC:59). Residuals against it are certified: ‖(L*−λ̃)f₃‖, ‖r₃‖ (LC:198-217); ‖(L−λ̃)ψ̂₃‖_ρ ≤ 1.2772e-6 (CP:211, `step3nl/psi3res3.py`). The verifier's −0.2561532 has smaller modulus than λ̃, so it is conservative for d\|α\|/dτ ≤ LAMT\|α\| + \|β\| [read] | step-3 certificate. **Step 2 does not certify λ̃**: S2:40-47 covers only Re λ > −1/8, and S2:50-52 gives −0.2561532 as a numerical value |
| 6 | Symmetry modes ψ₁ = p (λ = 1), ψ₂ = ±⟨e₁,ν⟩ (λ = 1/4); L maps span{ψ₁,ψ₂} into itself | constraints f₁,f₂ (LC:44-45); triangular lemma (LC:209: "(I−Π)LΠ = 0") | Exact algebraic identities Lp = p, L⟨a,ν⟩ = ¼⟨a,ν⟩. Checked symbolically by `s1spec/lincheck.py:15-17` (a file of the step-2 bundle, but the identities are elementary and stated in S2:35 as "one verifies directly") | elementary identity; not Theorem B. The paper should prove it in §2 (rescaled frame) |
| 7 | Modulation solvability (Gram ⟨ψ_i,f_j⟩ invertible, ‖M⁻¹‖ ≤ 2.0005), projection Π, triangular constants | P32:129; LC:198-217 | `step3/run_lyap.py` integrals (LC:234) | step-3 certificate |
| 8 | Slow-manifold functions ψ̂₃, ĥ₂ (exponents A = 1−4λ̃, 1−8λ̃) and far exponent p = 1.0246 | P32:59; (A) and U: `U = 1.3056αy^p` | `step3nl/psi3fun.py:16`, `h2fun.py:18` (built from λ̃); certified residuals (CP:136-143, 211-225) | step-3 certificate |
| 9 | Spectrum "1, ¼, −0.2562, −0.6999 … [C, step 2]" | P32:52 (tag \tC), `step3nl/STEP3_ARGUMENT.md:17`, `step3nl_v33/STEP3_ARGUMENT.md:17` | Descriptive only. −0.6999 is used nowhere, and the gap between the slow mode and −0.6999 is never used; coercivity (item 1) replaces it | **mislabelled.** Theorem B certifies only 1 and ¼ here; the rest is the numerical remark S2:50-52 |
| 10 | Semigroup decay ‖e^{τ𝓛}u‖ ≤ 14.8e^{−0.0213τ} | displayed in P32:131 and LC:228-231 | Needs a semigroup. Generation is the planning note `writeup/step3_blueprint.tex:63-69` (Lemma "generation"), which uses **step 2's numerical-range constants Λ₀ = 1.7909, θ₀ = 0.48** | **the only step-2 constant in the step-3 lineage.** Not consumed: `grep` finds 14.8 and 0.0213 in no v3.3 file, not in the verifier and not elsewhere in P32. The energy argument uses item 2 (P32:133) |
| 11 | Nonlinear flow and continuation | `fix1/b12_text/continuation_v33.tex` | Literature (P1′)–(P3′), not the linear semigroup | literature (L1–L3) |

Verifier inputs (`step3nl_v33/verifier/constants.py`, `closure.py`) [read]:
- The verifier reads only step-3 JSON files: `cubic_E0fine_0.001_0.01_0.01.json` (line 18) and `agent_int/planB/planB_ramp.json` (line 46). Its other inputs are typed constants with step-3 provenance: LAMT, PI3, the mod_ibp constants (lines 63-74) and FLUX0/FLUX2 (line 65).
- No step-2 file is read (`cap/`, `apr/`, `contour_result`, the Evans enclosures).

Code-level contacts with the step-2 bundle (code reuse, not logical dependence) [read]:
- **`lin/lintail2.py`.** Imported by `step3/hstail.py:18`, which feeds the Morse count and the step-1 `profile_bounds`. Also imported by `step4/p2/v33/n1_chart.py:34`, which uses only `LT.mu`, at `n1_chart.py:109`. Pinned in `final_tail_v33.py:76,146`.
  - The objects used are step-1 data: `RBAR` from `s1_cap/tail.radii` at `lintail2.py:34-40`, plus the formulas `mu` and `Efun` at `lintail2.py:44-50`.
  - The paper should attribute these to Theorem A. The archive must ship `lintail2.py` with the step-3 and step-4 code.
- **`evans.py`** (float profile `PROF`, `STIP`). Imported by `step3nl/psi3md.py:4`, `h2md.py:9`, `psi3hat.py:7`, `adjcut.py:6`, `nlflow.py:10`, `agent_mod/float_gram.py:6`, and by `step4/measure.py:7`. These are untrusted float fitting and diagnostics. The certified ψ̂₃ and ĥ₂ are checked by `psi3res3.py` and `h2fun.py` against step-1 hulls.

## 2. Trace B: the R3 items of `lemmaS_A_v33.tex` (G-A5, line 320)

`A` = `writeup/lemmaU/bulk/fix1/b12_text/lemmaS_A_v33.tex`; `LS` = `writeup/lemmaS/lemmaS.tex`; `ASM` = `writeup/lemmaU/bulk/fix2/n12/t2_A_v35/base_ibp/rigor/assembly_arb.py`; `LU` = `writeup/lemmaU/bulk/fix1/lemmaU_rev/lemmaU_v34.tex`.

| Item | Where used | Finding | Class | Effect |
|---|---|---|---|---|
| **R3.1** cos θ constant in c_s | A:76 (c = ¼ − a − b cosθ e^{−s}); A:190 (loss (\|1−2j\|+1)b̄/14); LU:40, LU:112 (uses β_s′ = −β_s) | (A) works in the ray chart y = Re(e^{−iθ∞}Γ), h = Im(…) (`chart_lemma_body.tex:9`). A translation by b·e₁ in this chart gives h_τ ∋ −b cos θ₀ h_y − b sin θ₀ with θ₀ the **fixed** angle between e₁ and the ray, so f_τ ∋ −b cos θ₀ f_y. The term is exact, and only \|cos θ₀\| ≤ 1 is used. The ray passes through the dilation centre because p = ⟨Γ,ν⟩ → 0 on the S₁ tail (step-1 tail bounds). The modulation terms of the slope equation are written by hand (LS:94): the symbolic check `lemmaS/lemmaS_symbolic.py` covers only a = b = 0 | **justified (exact)** once the derivation is written | 0. Under the alternative reading (θ = S₁ tangent angle, y = σ), the extra term in c_s is ≤ b̄·sup_{σ≥14} σ\|k\|/14 ≤ 2.19e-7 · 4.43e-7/14 ≈ 7e-15 [cert-check: `step3/hstail.profile_bounds`, σk ≤ 4.43e-7 at σ = 14, ≤ 2.1e-17 at 28] |
| **R3.1′** cos θ in R | LS:96 (R ∋ b_m cos θ U′) | same argument: exact in the ray chart | **justified** | 0 |
| **R3.2** forcing constant r₀^F = 1.4(β̄ + 1.02α₁ā), \|∂_s^jR\| ≤ r₀^F p^j y^p | A:266; ASM:105, 113, 166; B0_REPORT §1 ("the forcing of f − U(α(τ),·)") | For z^A = f − U(α(τ),y), U = 1.3056αy^p, the exact source is −U_τ − c·pU − ε(∂²−5∂+6)𝒢_U + (b-transport of U) = −1.3056y^p[β + (λ̃ + p/4)α − paα] − (b cosθ/y)pU − ε-term. The scripts use λ̃ = −p/4, explicitly at P32:235-236. With p = 1.0246 exactly, λ̃ + p/4 = −3.2225e-6, so three things are missing: **(i)** 1.3056·3.2225e-6·α₁ = 6.31e-10 (the 1.4-vs-1.3056 slack on β covers only 3.15e-10); **(ii)** the ε-part 0.0634·αy^{p−4}, whose s-derivatives carry \|p−4\|^j; **(iii)** for y < 14.8 the a/b coefficient 1.3377(1+1/y) exceeds 1.4·1.02 = 1.428 (by ≤ 0.4% of a term that is 1% of r₀) | **unjustified as written; harmless** | Honest ‖∂^jR‖ over the certified one: ×1.069, ×1.074, ×1.089 for j = 0, 1, 2 [est]. My reconstruction of the certified F² reproduces `numbers_A_v33.tex:40` (5.6506e-13, 5.9320e-13). The forcing share of Ê is ≤ 0.57%, 0.063%, 0.0055% (the G-A1 numbers, A:316), so Ê rises by ≲ 0.1% at level 0 [est]. No margin is affected (the smallest is Ĝ₂₂ ×1.917, which does not depend on F) |
| **R3.2′** (new, same root cause) Lemma U δU | LU:55-60 (hypothesis "α_τ = −(p/4)α + β", \|β\| ≤ β̄ = 3.33811e-9 from D1) | The D1 β̄ bounds α_τ − λ̃α, not α_τ + (p/4)α. The correct bound is β̄ + 3.2225e-6·α₁ = 3.8215e-9. Then c_δU = 1.3056(β̄′ + pāα₁)/(p(¼ − ā)) ≈ 1.965e-8 instead of 1.7187e-8 (×1.143) [est] | **unjustified hypothesis; very likely harmless** | c_δU enters U's launch data, N̄ (LU:372, term p^k c_δU 700^p √ℓ₁ ≈ 1.4e-5 against √(βB̂/w_min)), and zone III n^III. Not quantified by a script. Rerun needed |
| **R3.3a** flux lemma: cutoff shape | A:309; LS:140-148; `lemmaS/lemmaS_cflux.py:65` (`zeta_expr = 1 - Sx((x-57)/3)  # quintic … (ASSUMED shape)`) | The verifier's E_I cutoff is ζ = S7((60−x)/3), with S7 = u⁴(35−84u+70u²−20u³) (`step3nl/agent_int/v3/int_f23.py:14-18`; `ASSEMBLY.md:17`). The κ tables were built for the quintic, which has ζ⁽⁶⁾ ≡ 0 and lower maxima of the derivatives. The IBP reduction (`reduce()` in `lemmaS_cflux.py`) moves derivatives onto coefficients and treats ζ as smooth across x = 57 and 60, so endpoint jump terms of ζ⁽⁴⁾/ζ⁽⁵⁾ are not visibly accounted for [read; not re-derived] | **real gap (wrong object)** | Unknown factor on κ; plausibly ×2–×50 at level 2 [est, from the ratio of S7 to quintic derivative maxima]. Each K carries W(57) ≤ 1.44e-18. FLUX margins at β_A = 2 are ×3.365 and ×5.955 (`numbers_A_v33.tex:50-51`). A ×4 FLUX rerun moved only the flux-driven verifier quantities (`step3nl_v33/V33_VERIFIER.md:10`). So it is fixable by a rerun, and if κ grows past the margins, by raising FLUX0/FLUX2 again |
| **R3.3b** flux lemma: higher-jet box | `lemmaS_cflux_num.py:30` (`MBOX = M0 0.25, M1 0.4, M2 1.4, M3 5.8, M4 29.4, M5 150`); LS:154 lists "the jet box (H) on V" as an input | v3.3 supplies only levels 0 and 1 pointwise on the core (B2′: \|f\| ≤ 0.11, \|yf_y\| ≤ 0.52). I found no v3.3 statement that bounds V^{(2..5)} pointwise on [57, 60] (searched M3, A, `continuation_v33.tex`, LU) | **real gap (unproved hypothesis)** | Numerically irrelevant at W ~ 1e-18 unless the jets are astronomically large, but it is an unproved hypothesis of a certified constant |
| **R3.3c** flux lemma: records | `lemmaS_cflux_final.py:3` reads `cflux_raw.json`, written by `lemmaS_cflux_num.py:102` | `cflux_raw.json` exists nowhere under `/home/claude/s1spec` (find). The κ in `lemmaS_cflux.json` cannot be traced to a stored arb output | **records gap** | — |
| **R3.3d** flux lemma: linearity and nonnegativity ("asserted", A:59, A:309) | ASM:159-173 (`f4b`) | K₀ = Σk⁽⁰⁾_l N_l, K₂ = Σ_{l≤2} k⁽²⁾_l N_l with k ≥ 0 (sums of sup\|c_ab\|/W; `lemmaS_cflux_final.py:8-13`). N₀ = W·57.75·E^u₀, N₁ = W·E^u₁/54.75, N₂ = 2W(E^u₂+E^u₁)/54.75³, E^u = (60/28)E. The conversions are correct for ω = 28/y, dx = y ds, y ∈ [54.75, 57.75] [read] | **justified** (by code reading; upgrade "asserted" to "by construction") | — |
| **R3.3e** level-3 part of the flux (I₃ = k⁽²⁾₃N₃) | ASM:173 (`I3c`, `I3s` computed); verifier F0 = μ(FLUX0 + ε_I·FLUX2) (`closure.py:144`) | `f4b`'s K₂ excludes the N₃ term. I did not find I₃ consumed in the verifier (`grep I3` gives nothing in `step3nl_v33/verifier`) | **[not verified]**: either wired in somewhere I did not find, or missing | v3.1 size: I₃ ≤ 1.84e-21 + 1.39e-27(T−τ₁) (LS:11). As an impulse in Λ that is ~μ ε_I · 2e-21 ≈ 5e-27 against Λ ~ 1e-15 [est] |
| (W2) verifier tail lemma | A:320 | outside the three R3 items; not audited here | **[not verified]** | — |

## 3. Recommendations

**Step 2 in the paper (D4).**
- Present Theorem B as an independent result that explains the spectral picture: the four symmetry modes and a certified gap of 1/8. Say plainly that the Main Theorem does not use it.
- The proof chain is A + C + D + (L1)–(L3). Step 3's linear input is its own coercivity certificate (Morse count 3, Gram Schur complement, positivity at infinity), built from Theorem A data.
- This removes the Evans multiplicity citation from the trust base of the Main Theorem.
- Concrete edits:
  - Fig. `dep`: drop the arrow from step 2 to step 3, or mark it "context only".
  - Delete "it rests on the linear theory of step 2" (`main_v3.tex:118`).
  - In Remark (c) (`main_v3.tex:122`), say that the multiplicity citation concerns Theorem B only.
  - In the consolidated step-3 text (D7), replace P32:52's "(step 2) \tC" with: λ̃ is a parameter; Theorem B places no other eigenvalue in Re λ > −1/8 (numerically the next even one is −0.6999); none of this is used.
  - Prove Lp = p and L⟨a,ν⟩ = ¼⟨a,ν⟩ in §2 rather than citing `lincheck.py`.
  - Do not state the semigroup-decay corollary (14.8e^{−0.0213τ}), or state it with its own dependence on Theorem B's numerical range. Nothing uses it.
  - List among the cited classical facts of step 3 the Morse index theorem (Reid) and the constrained-inertia lemma (LC:86, 103). They replace the Evans citation in the cited class.
  - In the archive notes, say that `lin/lintail2.py` is a step-1 constants module used by steps 3 and 4.

**R3 items.**
1. **cos θ.** Remove the item from the "inherited" list. In the paper, derive the modulation terms of the slope equation in the ray chart in three lines (translation along e₁ gives −b cos θ₀ ∂_y; dilation about the origin gives −(¼−a)y∂_y) and define θ₀ as the fixed angle between the axis and the ray.
2. **Forcing.**
   - Replace r₀^F by the honest bound r₀′ = 1.3056(β̄ + \|λ̃ + p/4\|α₁) + 1.3377α₁ā(1 + 1/14), and add the ε-term c_ε α₁ \|p−4\|^j y^{p−4} with c_ε = 1.3056(p²−5p+6)(p²−p) + (cubic U terms). Or keep 1.4 and add the two missing terms.
   - In Lemma U, replace the hypothesis α_τ = −(p/4)α + β by α_τ = λ̃α + β, and use β̄ + \|λ̃ + p/4\|α₁ in c_δU.
   - Then rerun t2 → t4 → zone3 (minutes, but it rewrites v35 outputs, so ask first). Expected change: Ê below 0.1%, c_δU ×1.143, all margins essentially unchanged [est].
   - Until the rerun, list this in the paper as a known, quantified correction.
3. **Flux lemma.**
   - Recompute the κ tables with ζ = S7((60−x)/3). Do the IBP without moving derivatives across the endpoints, or keep the endpoint jump terms. Store `cflux_raw.json` with an args block.
   - Replace the a-priori jet box by a bound the bootstrap supplies, or restructure. One option: bound ⟨v⁗, H⟩ with H containing A_i, A_i′, A_i″ only, which needs jets up to V‴, paired with the verifier's V₄ (time-integrated J₄).
   - Find or wire in the consumption of I₃.
   - Until then, the paper must list the flux lemma as "certified for a model cutoff under an a-priori jet hypothesis; numerically ≥ 10¹⁰ below any closure quantity", not as certified.

## 4. Trace C: Evans multiplicity hypotheses (Theorem B only)

Theorem B certifies two things: the eigenvalues in Re λ > −1/8 are exactly the four symmetry values (S2:40-47), and each is a simple zero of D_even/D_odd (winding number 2 against two known zeros). Both are certified as stated. "Algebraically simple" needs the Gardner–Jones theorem: order of zero = algebraic multiplicity. Its hypotheses in this weighted half-line setting:

| Hypothesis | Status in `step2_proof.tex` |
|---|---|
| H1. Eigenvalue problem as a 4-dimensional first-order system on [0,∞), analytic in λ, with a λ-independent 2-dimensional boundary subspace at x = 0 (even u′ = u‴ = 0; odd u = u″ = 0) | verified (S2:124-130) |
| H2. Admissible (decaying) subspace at +∞ is 2-dimensional and analytic in λ on a neighbourhood of the region | verified: Lyapunov–Perron at σ₀ = 32, uniqueness and contraction; interior cover by 174 boxes (S2:114-115, 143) |
| H3. Admissible solutions are exactly the solutions with u^{(j)} ∈ L²(ρ); the other two grow faster than ρ^{−1/2} | stated: "every other solution grows like exp((a/2 − o(1))σ^{4/3})" (S2:115). The comparison with the weight exponent ¼ is not written out [not verified in detail] |
| H4. Zeros of D = eigenvalues | verified (S2:129-130) |
| H5. L is a closed operator on L²(ρ) (or on its even/odd parts), the λ₀ in question are isolated eigenvalues of finite algebraic multiplicity, and L − λ is Fredholm of index 0 nearby. Needed so that "algebraic multiplicity" is the rank of a Riesz projection | **not verified.** S2 defines eigenvalues through the ODE only (S2:37). The planning note `step3_blueprint.tex:63-80` sketches m-sectoriality and compact resolvent ("sketch", "modulo standard ODE asymptotics"). It is not part of the step-2 proof |
| H6. Jordan chains of the operator correspond to λ-derivatives of the analytic admissible basis that stay in the domain uniformly in λ, and the order of the zero of D equals the chain length (adjoint system and a constant bilinear concomitant, with boundary terms vanishing at 0 by parity and at ∞ by the weight) | **not verified.** S2:149-151 says so explicitly |

If simplicity is wanted without the citation, there is a direct route [not done]: for λ₀ ∈ {1, ¼}, certify that the adjoint eigenfunction ψ* (step 3 already computes adjoint data for λ = 1, ¼ in `step3/wstar_all.json`, CP:10-13) has ⟨ψ, ψ*⟩_ρ ≠ 0. That excludes a Jordan chain once H5 holds. Since A holds, none of this is needed for the Main Theorem.
