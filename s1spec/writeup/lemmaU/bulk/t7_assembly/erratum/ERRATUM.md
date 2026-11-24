# Erratum to step 3 v3.1 / v3.2, the step-3 verifier (fin11) and step 4 (task 7a)

Nov 13, 2026. Bulk task 7a (`../../BULK.md`). Nothing closed was edited. The tex version of this document is `erratum.tex` (balanced, not compiled: no LaTeX on this machine).

**Labels.**
- **[cert]**: arb with outward rounding, computed by a project script whose output is cited.
- **[est]**: float model or estimate.
- **[rev]**: a reviewer's scratch estimate.
- **[paper]**: a hand argument.
- **[hyp]**: a hypothesis.

**Line references.** Every file:line below was re-read for this erratum. Where a source report's citation is wrong, the corrected line is given and marked *(cit. corr.)*.

**Closed files covered.**
- v3.1 step 3:
  - `writeup/step3_proof/step3_proof.tex` and `writeup/lemmaS/lemmaS.tex`;
  - the scripts `writeup/lemmaS/*.py` and `writeup/lemmaS/rigor/*.py`.
- v3.2: `writeup/step3_proof_v32/step3_proof_v32.tex`, `writeup/step3_amend_v32/amend_v32.tex`, `writeup/lemmaS_v32/`.
- The verifier: `step3nl/verifier/` (fin11).
- Step 4: `writeup/step4_proof/step4_proof.tex`, `step4/p2/t9_lemmaO.py`, `step4/p2/out/FINAL_manifest.json`.
- Informational (working notes): `writeup/step3_gap_v32/gap_v32.tex`, `feas_P.py`.

## Summary

| group | items | effect certified? |
|---|---|---|
| A. Lemma S(A) | 10 (A1–A10) | yes, except A5 (the I₃ consumer was not identified) and the (B)-window part of A7 (not re-audited) |
| B. Verifier | 2 (B1 = E-V1, B2 = H-J window uniformity) + the v3.3 rerun note | yes |
| C. Go/no-go and gap notes (informational) | 2 (C1, C2) | no: corrected values are [est]; superseded by certified Phase-1 numbers |
| D. Step 4 | 1 (D1) | far family: yes (script). Cone family: [paper] argument on certified output; no script checks the implication |
| E. Status | — | see §E |

**The status needs one correction to the plain wording "v3.1/v3.2 holds with corrected constants under (F1) and (O)".**
- That holds only if the step-3 hand-off list is **extended by (H2f)** (A9).
- With the closed hand-off list alone, the repaired Lemma S(A) **fails** its flux input at the bootstrap level 2Ê: ×0.7317 [est from cert components]. This is hand arithmetic from a probe JSON (2K₀ = 5.60e-21 against FLUX0 = 4.1e-21), not a certified output, so it carries no [cert] label (review M8).
- (H2f) is certified by step 4, but in the v3.3 stage set (`FINAL_manifest_v33.json`, 21 stages, sha256 8bc47e83...), not in the closed manifest. The final arb row (review M9, stage `h2f_camp_arb`) gives G = (3.1363e-7, 1.0034e-6, 2.4351e-6), margins ×3.188 / ×4.983 / ×4.106 against (1e-6, 5e-6, 1e-5), and ×1.59 at level 0 against the design value 5e-7 used below.

---

## A. Lemma S(A) (v3.1; v3.2 inherits it unchanged)

v3.2 states that Lemma S(A) is unchanged (`step3_proof_v32.tex:211`, `amend_v32.tex:19`, `:74`). Every item in this group therefore applies to v3.2 too.

**Repaired chains.**
- Phase 2a: `writeup/lemmaU/phase2a/A_repair/` (`run_chain.sh`; `out/*.json`).
- Bulk task 2 (fin11 verifier): `writeup/lemmaU/bulk/t2_A/` (`run_t2.sh`; `out/A_interface.json`).
- **v3.3 chain** (fin12f4 verifier): `writeup/lemmaU/bulk/t2_A_v33/` (`out/A_interface.json`).

The v3.3 proof text is `writeup/lemmaU/bulk/t3_proof/lemmaU.tex` together with `t1_statements/STATEMENTS_r2.tex`. Assembly into a v3.3 step-3 file is task 7.

### A1. Core class {2,3,3}: certified with an assumed sup|z₂| ≤ 3, while the tex claims integration by parts
- **Closed code:**
  - `lemmaS/rigor/weights_exact.py:24` (`M2STAR = 3`);
  - `:172` (default `s33='assumed'`);
  - `:246` (`core(..., s33='assumed')`).
  - The certified run is `lemmaS_F10p_M0.52.json` (args: `'s33': 'assumed'`).
- **Closed tex:**
  - `lemmaS.tex:19` (H-L4: "closed by the L⁴ identity and one integration by parts") and `:349`;
  - `lemmaS.tex:89` calls M₂* = 3 superseded.
  - The quotes: γ′₂ = 0.143 at `lemmaS.tex:11`, `step3_proof.tex:211`, `:300`; Ĝ₂₂ = 6.58e4 (×2.33) at `lemmaS.tex:27`, `step3_proof.tex:224`.
- **Wrong:** the script did not do what the text claims.
- **Effect** [cert], with `--s33 ibp` at the v3.1 (H2) (`writeup/step3_gap_v32/gonogo/lemmaS_A_ibp_v31H2/rigor/assembly_arb.json`):
  - E₂ = 0.24704515 (v3.1: 0.24704507);
  - γ′₂ 0.14306 → **0.12355**;
  - core Ĝ₂₂ 6.58e4 → **77720.9** (margin ×2.338 → ×1.977);
  - every box margin is unchanged.
  - In the edge-free (B) lines: E_B (y_E = 50) 0.250 → 0.25014; γ′_B 0.2394 → 0.2376.
  - The claim holds; only quotes change.
- **Repair:** every repaired chain uses `s33 = ibp` (`t2_A_v33/base_ibp/`).

### A2. Top pair (j+1, j+1) on the ramp [14, 28] silently zeroed
- **Closed code:**
  - `lemmaS/lemmaS_F2.py:81` (`out[j + 1] = {m: 0.0 ...}`);
  - `lemmaS/rigor/ramp_arb.py:72` (tables for k ≤ j only) *(cit. corr.: A_REPAIR's `ramp_arb.py:105` does not exist; the file has 80 lines)*;
  - `lemmaS/rigor/assembly_arb.py:111` (comment "k=j+1 table is 0").
- **Closed tex:** `lemmaS.tex:209` ("no Young step, only J_k, k ≤ j"). This sits in a section that `lemmaS.tex:89` labels history, but the page-1 statement (`lemmaS.tex:11`) and `step3_proof.tex:211`, `:290` (D3) rest on that computation.
- **Wrong:** the cutoff-derivative part of the top pair cannot be reduced by the stated integration by parts. It needs a Young step against the ramp dissipation.
- **Effect** [cert] (`phase2a/A_repair/out/ramp_rep.json`, `assembly_rep.json`):

  | level | route | added to E |
  |---|---|---|
  | 0 | J | +4.96e-8 |
  | 1 | J | +1.04e-5 |
  | 2 | IBP + Young, δ = 0.3178 a_min, inside the certified slack 0.3200 a_min | **+0.00984** |

  - At level 2 the parts are: Young 0.00317, c″ jet-free 0.00129, jets of c″ 0.00539.
  - Route J at level 2 would cost 3.03.
  - The proof as written is incomplete; the claim holds with corrected constants (combined in §E).
- **Repair:** `t2_A_v33/rigor/ramp_rep.py`, Lemma P (`STATEMENTS_r2.tex`). One caution: phase 2a's sentence "c, c′ vanish at t = 0, 1" is false as written (A2_REPORT). Level 2 is still correct by Lemma P(i).

### A3. Cutoff × jet classes on (A)'s ramp, including the jets created by the IBP derivatives
- **Closed code:**
  - `lemmaS/rigor/ramp_arb.py:45`: every Leibniz derivative is bounded with the jet-free `Abr[(i, r + u)]`;
  - `ramp_arb.py:30–35` and `:50` (`expand_ab`, the full reduction);
  - `lemmaS_F2.py:49`, `:24–30`, `:56`, and `:138` (`Abr = A_bounds_s(...)`, jet-free).
  - *(cit. corr.: A_REPAIR cites `ramp_arb.py:78, 82–83`, which do not exist; the lines are 45 and 30–35 / 50.)*
- **Closed tex:** `lemmaS.tex:20` (H-L4r "closed by the space–time identity, J₂,₃,₄, pointwise Λ").
- **Wrong:**
  - The jet monomials of A_i^{(r+u)} are omitted.
  - The full reduction of the pairs (0,3) and (1,3) at level 2 needs A^{(r+u)} up to order 7. With honest jets it creates z₅ and z₄ against non-vanishing weights, and no datum controls those.
- **Effect** [cert] (`phase2a/A_repair/out/ramp_rep.json`):
  - level 0: 1.9e-13; level 1: 4.6e-7;
  - **level 2: +0.01176** in total:
    - (2,2,2) 0.00777, via Φ = 3.18e-5 with sup weight 245, of which the IBP-derived jets contribute 190;
    - (2,3) 0.00316;
    - (2,2,3) 7.1e-4;
    - (2,2) 1.05e-4;
    - (2,4) 1.9e-5.
  - The earlier "≲ 1.5e-3" [rev] was low by about ×8.
  - The claim holds with corrected constants.
- **Repair:** pairs (a, j+1) with a < j get one integration by parts only, and the class handlers are asserted (`t2_A_v33/rigor/ramp_rep.py`).

### A4. Profile jets U_j (j ≥ 2) dropped from the quadratic coefficients
- **Closed code:** `lemmaS/rigor/weights_exact.py:81–83` (`ab_exact`: V = [M0, M1, 0, …]) and `:88`.
- **Effect** [cert] (`phase2a/A_repair/out/weights_rep.json`):
  - core Ĝ_jj = 2280.85 / 2351.49 / 77720.88 → **2281.53 / 2372.05 / 79270.71** (Ĝ₂₂ margin against 153664: ×1.977 → ×1.936);
  - ramp Ĝ₂₂: 235592 → 239176;
  - γ′₂: −0.0025.
  - The claim holds.
- **Repair:** the honest Ab (`t2_A_v33/rigor/arlib.py`).

### A5. I₃ assumes the full 2·a_min of level-1 dissipation
- **Closed code:** `lemmaS/rigor/assembly_arb.py:175` (the factor `fac`, with the 2 a_min). v3.1's core budget slack is 0.3 / 0.2 a_min at levels 0 / 1 (`lemmaS/rigor/weights_exact.out:2–3`). I₃ uses level 1, so it has 0.2 against the 2 assumed, a factor 10. *(cit. corr.: CHECKPOINT finding 4 quotes 0 at level 2; line 4 of the file shows 0.4.)*
- **Closed tex:** `lemmaS.tex:11` (I₃ ≤ 1.84e-21 + 1.39e-27(T − τ₁)); `lemmaS.tex:21` and `step3_proof.tex:221` (H-F4 "closed").
- **Effect** [cert], honest I₃ with the retained r₁ = 1.249 a_min after the re-split:
  - at the v3.1 (H2) with (H2f) design: **3.370e-21 + 1.603e-25(T − τ₁)** (`t2_A/out/A_interface.json`, `res.v31.v31_only`);
  - v3.3 setting: slope 1.617e-25 (`t2_A_v33/out/A_interface.json`).
- **Consumer not identified.** No verifier input named I₃ exists (`step3nl/verifier/constants.py`), and `step3_proof.tex:200` lists only K₀ and K₂ as flux inputs. `lemmaS.tex:27` says the F4 constants are "≥ 10³ below the verifier placeholders"; CHECKPOINT gives ×10³ room [rev].
- **Not certified here:** the effect on any downstream inequality.

### A6. UQ = 1/20 used beyond its range
- **Closed code:** `lemmaS/rigor/weights_exact.py:24` (`UQ = 1/20`). It is used:
  - in `window()` (`:277–287`) for windows up to y = 380, where U(380) = 0.086;
  - in the core for y ∈ [220, 350], where U reaches 0.079.
- **Effect** [cert]:
  - segment check on [220, 350]: ε(220)·Ĝ(U = 0.08) ≤ ε(28)·Ĝ_core for every jk (`phase2a/A_repair/out/weights_rep.json`);
  - it also passes with the B3U boxes (0.161, 0.61) (`t2_A/out/weights_t2.json`);
  - windows with the honest U(1.5 y_E): Ĝ₂₂ of window 2025/8, 33467 → 37268; window margins unchanged.
  - The claim holds; it was unproved in v3.1.

### A7. Per-region convention at y = 28
- **Closed code:**
  - `lemmaS/rigor/weights_exact.py:246`, `:254`: `core()` and `ramp()` build separate L⁴/S33 constants and separate YC[j+1];
  - `lemmaS/rigor/ramp_arb.py:29` (ρ) and `:57` (`YCr[j + 1]*rho(t)`, ramp-only constants).
  - The boundary terms at y = 28 (the interpolation jump and [Θ z₁ z₂³]) are dropped.
- **Wrong:** the identities are taken separately on the ramp (box (0.05, 0.9)) and on the core (box (0.11, 0.52)). No single μ exists: Y_ramp[3] = 1156.7 against Y_core[3] = 644.2.
- **Effect** [cert] (`t2_A/out/weights_t2.json`; `A_interface.json` `interface28`). The repair:
  1. sets the ramp coefficient box to the core box (0.11, 0.52), which is valid because the verifier certifies ramp |y f_y| ≤ 0.4396;
  2. uses one C^{1,1} Ỹ on [14, 700];
  3. localises L⁴/S33.
- **Results:**
  - core Ĝ₂₂ = **80141.6** (×1.917);
  - γ′ = (0.242926, 0.241064, 0.119615);
  - **ramp M1 margin ×2.05 (box 0.9) → ×1.185 (box 0.52).** This replaces the ramp-box statement of `lemmaS.tex:9`, `:17`.
  - The tightened ramp box is a design choice made by task 2's agent, to be checked in the task-8 review.
- **Not certified:** v3.1's (B) windows (`window()`, `weights_exact.py:277`) use the same convention and were not re-audited. They feed the (B) boxes, K2/K3 and the (B) exit. Those conclusions belong to Lemma S(B), which v3.1/v3.2 already take conditionally through (F).

### A8. The outer cutoff of (A) was never treated: new Hypothesis (O)
- **Closed code:**
  - `lemmaS/rigor/assembly_arb.py:102–104` (Rout = 350; the forcing integral stops there; no outer-cutoff term);
  - the `lemmaS_F1.py` docstring ("outer cutoff far away").
- **Closed tex:**
  - `step3_proof.tex:81` ("shell (A) core y ∈ [28, 60]");
  - `:272` (D3 uses B3, |V| ≤ ½ on y > 60, as the coefficient bound);
  - `step3_proof_v32.tex:261` and `amend_v32.tex:74` ("nothing in Lemma S(A) … fixed rescaled y ∈ [14, 60]"); `amend_v32.tex:19`.
- **Wrong:** (A)'s energies live beyond y = 60: to 350 in the code, to 700 with an honest cutoff. Nothing bounds the cutoff terms.
- **Effect:** a new **Hypothesis (O)**, (O1)–(O3) (`t1_statements/HYPOTHESIS_O_r2.md`):
  - (O1): |f| ≤ 0.11, |y f_y| ≤ 0.52 on (60, 350];
  - (O2): |f| ≤ ½, |y f_y| ≤ 1 on [350, 700];
  - (O3): |z_k| ≤ 1 on [350, 700], k ≤ 3.
- **Given (O)**, the outer cutoff is bounded through L² data (R5′, `t2_A/out/outer_L2.json`). Since (O3) gives N̄_k² = sup_τ ∫_{[350,700]} z_k² ds ≤ ln 2 [paper], I evaluated the certified affine export at N̄² = (ln 2, ln 2, ln 2) (`erratum/o3_outer_check.py` → `o3_outer_check.json`, two runs byte-identical) [cert]:

  | setting (verifier) | Ê at (O3) | FLUX0/(2K₀) | FLUX2/(2K₂) |
  |---|---|---|---|
  | v3.1 (H2) + (H2f) design (fin11) | (8.9538e-7, 2.02941e-4, 0.255138) | **×1.1245** | ×1.6036 |
  | tightened (H2) + (H2f) design (fin11) | (7.0538e-7, 1.03941e-4, 0.0561375) | ×1.4598 | ×2.4468 |
  | v3.3: tightened (H2) + (H2f) 1e-6 (fin12f4) | (1.20538e-6, 1.03941e-4, 0.0561379) | ×3.4546 | ×6.1153 |

  - **Scope of ×1.1245 (review M8).** It certifies the **flux line only** (FLUX0/(2K₀)), computed with the A5 ramp box (0.52, not the closed 0.9), the C^{1,1} Ỹ, and (H2f) at the design values (5e-7, 5e-6, 1e-5). It is not a certificate of the whole repaired Lemma S(A) at the closed settings. The closed-version claim is: "flux line certified under (F1) + (O) + (H2f), with the ramp box of A5". The step-4 value of (H2f)₀ is below 5e-7 (3.1363e-7, margin ×1.59, `M9_H2f` of `FINAL_manifest_v33.json`), so the design value dominates it; no script asserts G₀ ≤ 5e-7 itself.
  - Against the reference outer data N2_ref [est], Ê₀ rises 2.0%, Ê₁ 0.03%, and Ê₂ falls 0.2% [cert].
  - S0v, M1 and Ĝ₂₂ were not recomputed at (O3). Given those shifts they move by under 1% [est].
- **Repair:** in v3.3 the outer data come from Lemma U (zone III, Lemma N) and (O) is removed, subject to the task-8 review.

### A9. (A)'s initial datum covers [14, 60] only
- **Closed tex:**
  - `lemmaS.tex:25` (hand-off (2): energies "on [14, 60]");
  - `step3_proof.tex:265` and `step3_proof_v32.tex:301` ("(H2) gives E_j(τ₁) ≤ (2e-7, 1e-4, 0.2) < 2Ê").
- **Wrong:** (A)'s energy at τ₁ integrates over its whole support [14, 700]. The closed hand-off list covers [60, 379.7] only through the (H3) windows, which are pitched at the v3.1/v3.2 levels, and says nothing on [379.5, 700].
- **Effect** [cert] (`t2_A/probe/asm_0_h3_v32.json`, setting v31):
  - with the honest datum from (H2) plus the (H3) windows, and nothing counted beyond 379.5: Ê = (1.393e-6, 2.455e-4, 0.2716) and K₀ = 2.802e-21;
  - at the bootstrap level 2Ê of `step3_proof.tex:261` (B2), 2K₀ = 5.60e-21 exceeds the verifier input FLUX0 = 4.1e-21: **×0.7317** [est from cert components: hand arithmetic from the probe JSON, not a script output].
  - **The closed claim fails at the closed hand-off thresholds.**
- **Repair:** a new hand-off condition **(H2f)**: ∫_{60}^{700} (28/y) z_j² dy/y ≤ H2f_j at t₁.
  - Step 4 certifies G = (3.1363e-7, 1.0034e-6, 2.4351e-6) against (1e-6, 5e-6, 1e-5), i.e. ×3.188 / ×4.983 / ×4.106 [cert] (arb row of review M9, stage `h2f_camp_arb`; the float row gave the same values to 7.5e-5 relative) (`step4/p2/out_v33/h2f_result.json`; stage h2f_camp of `FINAL_manifest_v33.json`; not part of the closed `FINAL_manifest.json`).
  - With (H2f) at the design values (5e-7, 5e-6, 1e-5), which dominate step 4's certified G, the repaired (A) closes at the v3.1 (H2) with the fin11 verifier and factor 2: FLUX0/(2K₀) ×1.146 at N2_ref [cert, given the est outer data], ×1.1245 at (O3) [cert; flux line only, see A8].

### A10. Point terms of the C² cutoffs at levels 0–2
- **Closed code:** the reduction z_a z_{a+1} → ½ c′ z_a² drops the boundary term at y = 28 (`lemmaS/rigor/ramp_arb.py:32`, `lemmaS_F2.py:27`).
  - The cutoff S² (S the quintic smoothstep) has (η²)^{(m)}(1) = 0 for m = 1, 2 but not for m = 3. So the m = 3 monomial leaves a point term.
  - (Location [paper], from A2_REPORT A3 and the code.)
- **Effect** [cert] (`t2_A/out/a3_signs.json`; `A_interface.json` `pointterm.dE`):
  - on the stored forms the (2,2,0) coefficient on (j+1, j+1) is +4 and the (3,2,0) coefficient on (j, j+1) is +2;
  - g_c(28) = 5.862e-4;
  - dE = (+6.41e-9, +3.00e-6, +9.7e-9), i.e. 0.9% of E₀ᴬ and 2.9% of E₁ᴬ (tightened (H2)). At level 2, Lemma P(i) gives −g_c z₂² ≤ 0.
  - The claim holds.
- **Repair:** `t2_A_v33/rigor/a3_signs.py`, `ramp_rep.py`.

### Also flagged (assembly check, not a stated claim)
The (B)-exit check S0_B (y_E = 50) ≤ S0V = 0.0153 (`phase2a/A_repair/rigor/assembly_rep.py:151`, from `lemmaS_F4B.py:44–52`) fails with the repaired constants at the v3.1 (H2): **0.01854 > 0.0153** [cert] (`t2_A/out/A_interface.json`, `res.v31.v31_only`).
- No consumer was found. The verifier input S₀ᵛ is the core sup on [34.25, 60] (`step3_proof.tex:200`, `:282`), which holds at ×1.981.
- The v3.1 tex does not quote S0_B.
- For the task-8 review to confirm.

---

## B. Verifier (`step3nl/verifier/`, fin11)

### B1. E-V1: the flux floor F0 is missing on the march cells
- **Closed code:**
  - `closure.py:182`: the march sum Σ (r1⁺·Ly_hi + Φ)·D_k has no F0·D_k;
  - its docstring `closure.py:170` has +2F0;
  - `closure.py:175`: `set_mode(delta, True)`;
  - `rhs_iv.py:172`: the absolute flux term is added to L;
  - `rhs_iv.py:191`: `if K.DROPL: L = L * cc(0.0)` *(cit. corr.: R2_DECISIONS cites `rhs_iv.py:198`, which is inside `circ()`)*;
  - `unc.py:16–27`: `lk` has no constant term;
  - F0 enters only the tail slope, `closure.py:197`.
- **Effect** [cert] (`bulk/t2_A/out/hjw_rep.json`):
  - F0 = 3.28048e-22 per unit τ;
  - F0·τ_march/X_c = 4.81e-7 relative on every J_m (τ_march = 5.4813, X_c = 3.7417e-15).
  - No quoted v3.1 digit changes.
- **Handling:** no verifier file is changed. v3.3 uses J_m(1 + κ_J) (Lemma 4 of `STATEMENTS_r2.tex`).

### B2. H-J used window-uniformly without a check
- **Closed tex:** `lemmaS.tex:256` (layer-cake: ∫^τ e^{−γ′(τ−σ)}V_m² ≤ J_m + s_m/γ′). This needs J_m(T) ≤ J_m + s_m T on every window, not only on windows from τ₁ (KILLCHECK2 §2).
- **Effect** [cert] (`t2_A/out/hjw_rep.json`):
  - (W1) fails at every cell with P_k ≥ 0; (W2) holds (it relies on the verifier's tail lemma, not re-run);
  - **κ_J = 0.0423470 (W1) + 4.8e-7 (W3) = 0.0423475**, limit 0.13.
  - Every repaired (A) number in this erratum uses J_m(1 + κ_J). The claim holds with J ×1.0424.

### The v3.3 verifier rerun (fin12f4)
- `step3nl_v33/` (copy; `V33_VERIFIER.md`) with FLUX0 = 1.64e-20 and FLUX2 = 8.0e-20 (×4); `rate_table_fin12f4_a0.00015.json`, sha256 b20e8715….
- [cert]:
  - tail_cert PASS, log byte-identical to fin11;
  - J₀…J₄ unchanged to 1e-6;
  - slopes, sJ4 and F0 ×4 (F0 = 1.312192e-21); √L_F and α_F ×2;
  - α-condition margin a*/α_F 107.3 → 53.7.
- **The E-V1 code is unchanged in the copy.** It is absorbed again: κ_J = 0.0423485 (W1) + 1.92e-6 (W3) = **0.0423504** [cert] (`t2_A_v33/out/hjw_rep.json`).

---

## C. Go/no-go notes and gap note (informational; not closed steps)

### C1. `feas_P.py` sign bug
- `writeup/step3_gap_v32/feas_P.py:11` sets γ = −μ, and `:20–22` (`co*gam**m`) applies it to tables that store |coef| (`lemmaS/lemmaS_prop1_bounds.py:49`). Odd weight-derivative orders enter with a minus sign; the K_jj inherit it through YC (`feas_P.py:36`). *(cit. corr.: GONOGO cites `:22` and `:33`.)*
- **Effect:**
  - Ĝ = (21, 302, 2.7e3) is wrong; the honest values are about (4.9e4, 2.3e4, 7.8e4) before L⁴ [est] (GONOGO error 2);
  - the K_{j,j+1} entries are not affected;
  - B₂ moves by ≤ 2% [rev].
- **Superseded by:** Phase 1, certified Ĝ₂₂ ≤ 2.77e5 and the ramp K table [cert] (`lemmaU/phase1/ghat/`, `ramp/`); then `t4_cert`.

### C2. The `gap_v32.tex` quote
- `writeup/step3_gap_v32/gap_v32.tex:23` quotes Ĝ_jj = (21, 302, 2.7e3), ε(28)Ĝ₂₂ = 0.0044 and K_{2·}, all from C1.
- **Corrected:** ε(28)Ĝ₂₂ ≈ 7.8e4/28⁴ ≈ 0.13 [est].
- The gap statement itself (`gap_v32.tex:11`) is unaffected. The feasibility discussion is superseded by GONOGO.md and Lemma U.

(GONOGO.md lists further defects of the go/no-go scripts: the `t3_A_outer.py` crash, `t1` omitting (H3), the box test on [56, 84]. They are working-note issues and are not erratum items.)

---

## D. Step 4

### D1. `t9_lemmaO.py`: the arc-end condition |ξ| > XL is never tested
- **Closed code:**
  - `step4/p2/t9_lemmaO.py:63` computes `ends`, and nothing uses it;
  - the docstring `:6` claims the test.
- **Closed tex:** `step4_proof.tex:715` ("the graph arc … leaves [−XL, XL] at both ends").
- **Pinned outputs:** `p2/out/t9_lemmaO_v32.json` (sha256 07ea3d75…, pinned) and `t9_lemmaO.json` (v3.1, C₀ = 0.0021).

**Far family (v3.2 H1): not affected.**
- H1 for the far family comes from `t9_lemmaO_far.py` (stage 28), not from these rows. Its far rows fail H8 by design and are not used.
- That script checks coverage explicitly: ξ(±G) ≥ 226.18 > X L = 226 (`graph_covers_X`), turning ≤ 0.04095 rad, `ok = true` [cert] (`p2/out/t9_lemmaO_far.json`, sha256 8915a8a4…, pinned).

**Cone family (v3.2: 18 cylinders, L = r/4, X = 3): the claim holds by an argument on certified output; no script checks the implication.** [paper]
- The cells are consecutive closed s-intervals (`t9_lemmaO.py:36–39`). So the first cell past each arc end shares the arc's boundary point P with the last arc cell (both enclosures contain P(t) at each t).
- Suppose |ξ(P)| ≤ XL at some t.
  - Then the arc cell has Lo|ξ| ≤ XL, so it enters `hmax` (`:67`), and |η(P)| ≤ hmax < HL = hmax + C₀L.
  - The neighbouring non-arc cell then has distance 0 to the box. But H8 requires ≥ 0.1L for every non-arc cell (`:70–76`).
- **All 18 cone rows have `ok = true`, min H8 distance 0.1360 L** [cert]. Hence |ξ(P)| > XL at both arc ends, for every t.
- With cos(θ − θ_c) > 0 on the arc (`:62`, asserted) and the centre inside, ξ runs from < −XL to > XL. [−X, X] is covered once.
- Preconditions:
  - the arc start is a₀ ≥ 28.75 (rlo = 200), and cells exist below it (the anchor cell at x ≈ 4.8, `t9_lemmaO_far.json`);
  - the cone arcs end at a₁ ≤ 705, far below ℓ̃, so no mirror cells are involved.
- **So the closed v3.2 cone H1 (sup|h₀| ≤ 0.1872, H ≤ 0.1895) stands.** The same argument covers v3.1's `t9_lemmaO.json`.
- **Erratum:** `step4_proof.tex:715` should say that the end condition follows from H8, not that it is checked.

**v3.3.**
- The fixed-L family (L = 50, X = 3, centres 300–380) asserts the end condition on an arc lengthened to 1.3·XL: ξ ≤ −179.1 at the start and ≥ +170.0 at the end, against XL = 150; sup|h₀| ≤ 0.2162 [cert] (`step4/p2/out_v33/t9_lemmaO_v33.json`).
- The v3.2 r/4 cone family is superseded in v3.3 and is not a v3.3 input.
- The v33 report notes that a 1.1·XL arc fails *its* test. That test is stricter: it applies to whole end cells, in the v33 geometry.
- `t9_lemmaO_v33.py` reproduces v3.2 only to about 1e-9 relative (python-flint 0.9.0), not bitwise.

---

## E. Status statement

**v3.1 / v3.2 step 3.** Lemma S(A) holds with corrected constants (A1–A10, B1–B2) under three conditions:
1. (F1), as before;
2. **(O)** (new, A8);
3. the hand-off list **extended by (H2f)** (A9).

Notes on these conditions:
- (H2f) is certified by step 4 in the v3.3 stage set (`FINAL_manifest_v33.json`), not in the closed manifest.
- Without (H2f), the flux input fails at 2Ê, ×0.7317 [est from cert components; hand arithmetic from a probe].
- With it, the binding margin of the **flux line** at the closed settings (v3.1 (H2), fin11, factor 2, A5 ramp box, (H2f) at 5e-7) is **FLUX0/(2K₀) = ×1.1245** [cert, (O3) data; flux line only]. It is thin but above 1. The other quantities (S0v, M1, Ĝ₂₂) at (O3) are "under 1% [est]".
- (F2) is then a consequence (the repaired Lemma S(A)), so the accurate conditionality is "(F1) and (O)", plus (H2f) as a certified hand-off item.

Also open for the closed versions:
- the A7 ramp-box design (task-8 review);
- the (B)-window convention, not re-audited (A7);
- the I₃ consumer, not identified (A5).

**v3.3 (Lemma U).** It removes (F1) and (O), **subject to the task-8 independent review (Fable)**. Its datum is the tightened (H2) and (H2f), both certified by step 4 (v33 stages).

### Certified repaired numbers

Settings:
- **v3.1 column:** `t2_A`, v3.1 (H2) = (2e-7, 1e-4, 0.2) + (H2f) design (5e-7, 5e-6, 1e-5), fin11 verifier, κ_J applied, outer data at N2_ref [est] unless "(O3)".
- **v3.3 column:** `t2_A_v33`, tightened (H2) = (1e-8, 1e-6, 1e-3) + (H2f) (1e-6, 5e-6, 1e-5), fin12f4, β_A = 2.

All [cert].

| quantity | v3.1 quote (`lemmaS.tex:11,27`; `step3_proof.tex:211,224,294–302`) | repaired, v3.1 setting | v3.3 setting |
|---|---|---|---|
| Ê₀, Ê₁, Ê₂ | 3.12e-7, 1.81e-4, 0.248 | 8.777e-7, 2.0288e-4, 0.25558 ((O3): 8.954e-7, 2.0294e-4, 0.25514) | 1.18766e-6, 1.03883e-4, 0.056585 |
| γ′₀, γ′₁, γ′₂ | 0.246, 0.246, 0.143 | 0.242926, 0.241064, 0.119615 | same |
| S₀ᵛ (vs 0.0153) | 0.00575 (×2.66) | 0.007724 (×1.981) | 0.007137 (×2.144) |
| S₁ᵛ | 0.171 | 0.1773 | 0.1034 |
| core M0 / M1 margin | ×5.86 / ×2.82 | ×5.309 / ×2.723 | ×5.464 / ×4.443 |
| core Ĝ₂₂ (vs 153664) | 6.58e4 (×2.33) | 80141.6 (×1.917) | same |
| ramp M1 box | 0.9 (×2.05) | 0.52 (×1.185) | same |
| K₀, K₂ | 6.83e-22, 3.26e-21 | 1.789e-21, 6.151e-21 ((O3): 1.823e-21, 6.236e-21) | 2.339e-21, 6.456e-21 |
| FLUX0/(2K₀), FLUX2/(2K₂) | ×3.0, ×3.0 | ×1.146, ×1.626 ((O3): ×1.1245, ×1.6036) | ×3.505, ×6.196 ((O3): ×3.455, ×6.115) |
| I₃ | 1.84e-21 + 1.39e-27 T | 3.370e-21 + 1.603e-25 T | slope 1.617e-25 |
| E_B (y_E = 50), γ′_B | (3.15e-7, 1.83e-4, 0.250), 0.2394 | (8.890e-7, 2.0550e-4, 0.25888), 0.23718 | (v3.1 only) |
| (B) boxes M1 / M0 / z₂, \|f\| ≤ ½, K2 / K3 | ×2.05 / 2.54 / 4.26, ×5.38, ≥ ×4.8 | ×2.052 / 2.544 / 4.268, ×5.388, ×5.080 / 4.825 (A7 caveat) | same |
| κ_J | — | 0.0423475 | 0.0423504 |
| C0 (z₂, z₃) export | — | 8.623e-7, 2.4663e-4 | 1.11659e-6, 1.65440e-4 |

v3.3 only [cert]:
- Lemma U (`t4_cert/out/U_cert_v35.json`, sha256 960809fd…; v3.3 final, after B0, M14 and N1):
  - late M1 ×2.353, M0 ×4.892;
  - present M1 ×2.770, M0 ×2.225;
  - loop gain ≤ 0.03904.
- Zone III: H6 ×1.860 (`t5_zone3/zone3_interface.json`, sha256 45e8d4dc…).
- Step 4, v33 stages:
  - (H2f) ×3.188 / ×4.983 / ×4.106 (M9 arb row);
  - tightened (H2), w' ×124.7 / ×2124 / ×641; tightened (H2) in the variable f − U₁ (B0, `h2u_v33_arb`) ×27.34 / ×562.8 / ×600.9;
  - fixed-L cone slope ×2.64, H7 t*/ratio ×2.30, H4 ×1119;
  - far-family Lemma O H4 at every centre of hand-off radius ≥ 490 (M9, `h4far_v33`): minimum ×15.58.

**v3.3 final (B0, M14, N1; Nov 14).** The v3.3 column of the table above comes from `t2_A_v33` (before B0, M3, G-A1 and N1). **Open item (O3-v35), explicit and not blocking:** its (O3) entries have not been recomputed on the new chain; the table above supplies no numbers for them and no source exists [TODO: source for (O3) values on t2_A_v35; `o3_outer_check.py` would need a `t2_A_v35` entry in its TREES list]. They are not used by v3.3, where Lemma N derives the outer data instead of assuming (O3). Final values from `fix2/n12/t2_A_v35/out/A_interface.json` (78fe2f8e38a340ba515fbc33cbab60cdd85231095af6ccd192e5f6707f457e8c) at N2_ref, and `fix2/n12/t2_A_v35/out/asm_design_N2ref4.json` (46bae149…) at 4·N2_ref, all [cert] (rounded outward: bounds up, margins down):
  - Ê at N2_ref (EA_ref) ≤ (1.2016e-6, 1.0439e-4, 0.056690); Ê supersolution of Lemma U (`U_cert_v35.json`) ≤ (1.2018e-6, 1.0441e-4, 0.056694);
  - core, at N2_ref: S0v vs 0.0153 ×2.122 (the verifier input is S0v + sup_core|d|), core M0 ×5.455, core M1 ×4.435, Ĝ₂₂ ×1.917, ramp M1 ×1.184, flux at β_A = 2 ×3.396 / ×6.012;
  - the same at 4·N2_ref (the margins quoted in the texts): S0v ×2.116, core M0 ×5.450, core M1 ×4.409, Ĝ₂₂ ×1.917, ramp M1 ×1.184, flux ×3.365 / ×5.955;
  - B0: B2′ now reads sup|v_A| + sup_core|d| ≤ 0.0153 (margin ×2.116 for the derived value), and (A) carries the profile difference d at the datum, the ramp J conversion, the fluxes and S0V/S1V (`B0_REPORT.md`);
  - N1 (re-review): v3.1 identified the verifier's S1-chart variable w′ with the ray-graph quantity without proof, and v3.3 (B0) covered it only by a typed level-0 allowance 2e-7. The identification is now an exact lemma with certified bounds (`fix2/n12/chart_lemma.tex`, merged into `lemmaS_A_v33.tex` §chart; step-4 stage `n1_chart`, 341ed1f4…); every margin moved by at most 0.06%;
  - G-A1 (integration rule R-sup) is applied: Ê rose by (6.7e-9, 6.5e-8, 3.0e-6) [est, measured in M14].

## Files
- `ERRATUM.md` (this file) and `erratum.tex` (same content, balanced, uncompiled).
- `o3_outer_check.py` → `o3_outer_check.json` (item A8; args block with input sha256; reruns byte-identical).
