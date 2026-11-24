# ERRATUM_v41 (text v41): ERRATUM_v40.md verbatim, then a pointer to Part 6

Part 6 (text v41, delta check fix20) is section Q of `erratum_v41.tex`; its source is `tools/erratum_part6_v41.tex`.

---

# ERRATUM_v40 (text v40): ERRATUM_v39.md verbatim, then a pointer to Part 5

Part 5 (text v40, review fix18) is section P of `erratum_v40.tex`; its source is `tools/erratum_part5_v40.tex`.

---

# ERRATUM_v39 (text v39): ERRATUM_v38.md verbatim, then Part 4

Part 4 (text v39) is section O of `erratum_v39.tex`; its source is `tools/erratum_part4_v39.md`.

---

# Erratum v38 to step 3 v3.1 / v3.2 / v3.3, the step-3 verifiers (fin11, fin12f4, fin13v4, fin13v4b) and step 4

**v38, Nov 20, 2026** (`writeup/lemmaU/bulk/fix16_texts_v38/`, built by `tools/build_erratum_v38.py`). Parts 1 and 2 below are `fix12_texts/ERRATUM_v37.md` verbatim (its first line, the v37 title, is kept as a heading); Part 3 (sections K-N) is new. Where Part 3 settles or supersedes a statement of Part 2 (sections G, I, J), Part 3 governs.

## (v37 title) Erratum v37 to step 3 v3.1 / v3.2 / v3.3, the step-3 verifiers (fin11, fin12f4, fin13v4) and step 4

**v37, Nov 20, 2026** (`writeup/lemmaU/bulk/fix12_texts/`). The earlier erratum (`t7_assembly/erratum/ERRATUM.md`, title "Erratum to step 3 v3.1 / v3.2, the step-3 verifier (fin11) and step 4 (task 7a)", Nov 13-14) is kept verbatim below as Part 1 (sections A-E and its Files list; its relative paths are relative to `t7_assembly/erratum/`); the new material is Part 2, sections F-J, at the end.


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


---
---

# Part 2 (v37, Nov 20, 2026): sections F to J

Written for fix round 12 (`fix12_texts/`). Nothing closed or existing was edited; this file and `erratum_v37.tex` are new. Everything above this line is the Nov 13-14 erratum, kept verbatim.

**Paths in Part 2.** Relative to `/home/claude/s1spec/writeup/lemmaU/bulk/` unless absolute. Short names:
- `V33/` = `/home/claude/s1spec/step3nl_v33/` (verifier tag fin12f4);
- `V4/` = `/home/claude/s1spec/step3nl_v4/` (tag fin13v4);
- `V4B/` = `/home/claude/s1spec/step3nl_v4b/` (tag fin13v4b; production run `V4B/runs/P`).

**Labels in Part 2.** As in Part 1, plus:
- **[hand]**: a hand derivation, not machine-checked;
- **[read]**: verified by reading the cited file:line, or a certified input taken as given;
- **[sym]**: exact sympy identity;
- **[probe]**: an interval branch-and-bound evaluation of an inequality whose new constants are themselves [est] (so the result is [est] as a statement about the true inequality);
- **[float of cert]**: a float quotient of certified values.
Every number in Part 2 is copied from the named report or JSON. Where a report prints a ball midpoint or rounds a bound down, the outward value is given here (with the digits taken from the JSON) and the fact is noted.

**Source reports (short names used below).**
| short | file |
|---|---|
| AUDIT | `fix3/audit/AUDIT_RAMP.md` (G0, Opus, Nov 16) |
| GI | `fix5/g0_indep/G0_INDEP.md` (independent check of G0, Fable, Nov 19) |
| SPEC | `fix5/ei_spec/EI_SPEC.md` (E_I derivation; §7 old-to-new code table; §11 SF/U items) |
| V4R | `V4B/V4_REPORT.md` (fin13v4) |
| REV4 | `fix5/v4_review/V4_REVIEW.md` (review of fin13v4; items R-1..R-10) |
| E0A | `fix5/e0_audit/E0_AUDIT.md` |
| E2A | `fix5/e2_audit/E2_AUDIT.md` |
| V4BR | `V4B/V4B_REPORT.md` (fin13v4b) |
| R3R | `fix4_r3/R3_REPORT.md` (tag v36pre, on fin12f4) |
| V36R | `fix6_v36/V36_REPORT.md` and `fix6_v36/r3/r3_forcing.json` (tag v36, on fin13v4b) |
| FLX1 | `fix3/flux/FLUX_REPORT.md` (flux lemma v1, fix round 3) |
| ROUTE | `fix3/ROUTE_DECISION.md` |
| EF | `fix8_flux2_impl/texts/erratum_flux2.tex` (items E-F1..E-F5) and `fix8_flux2_impl/out/flux2_cert.json` |
| FREV | `fix9_review/FLUX2_REVIEW.md` |
| PRF | `fix10_proofs/PROOFS_REPORT.md` |

---

## F. Defects found Nov 16-20, 2026 (fix rounds 3-10)

Summary:
| item | subject | status |
|---|---|---|
| F1 | E_I inequality of the step-3 tip-zone verifier (fin12f4) is not an upper bound (G0) | coded terms fixed in fin13v4, then fin13v4b; written E_I lemma `ei_lemma_v4.tex` **open** |
| F2 | KREG[12,16.25] typed below its certificate | fixed in fin13v4 |
| F3 | E0 energy block (beta coupling, missing terms, rounding) | fixed in fin13v4b (lemma `V4B/e0_lemma/E0_LEMMA_v4b.tex`, partly [hand]/[read]) |
| F4 | E2 energy block (modulation couplings, forcing factor, nonlinear constants) | fixed in fin13v4b (Cauchy-Schwarz forms [hand]) |
| F5 | Strip hypothesis sup\|w'''\| <= 0.02 and the tail alpha-condition (REV4 R-2, R-3) | fixed in fin13v4b; one text conflict (F5c) |
| F6 | Forcing constant r0 of Lemma S(A) assumed lambda~ = -p/4 (R3); c_dU | fixed in v36 (`fix6_v36/`); texts not edited |
| F7 | Outer-ramp flux lemma: v3.1/v1 defects and the form of hypothesis (I3) | flux lemma v2 certified (`fix8_flux2_impl/`); text items **open** |
| F8 | V36_REPORT rounds two design-setting quotes down | record only; nothing consumes them |
| F9 | Text defects in the v3.3 texts (and one closed v3.2 line) | **open** (texts not edited) |

### F1. G0: the coded E_I inequality of fin12f4 is not an upper bound

**What was wrong.** The verifier's E_I inequality (level 0 and level 2, `V33/verifier/rhs_iv.py`, `constants.py`, `unc.py`, `closure.py`) has no written derivation, and the code is not an upper bound for (1/2) d/dtau (V0^2 + eps_I V2^2):
- AUDIT verdict KILL (process verdict: negative rounding reserve plus uncovered terms; a rerun is needed);
- GI verdict FIX-AND-RERUN: 18 terms of the v-equation absent (11) or under-covered (7), 17 typed constants rounded the wrong way; ASSEMBLY.md and the code disagree in at least 7 places [read].
- SPEC (Fable) wrote the derivation and the old-to-new code table (SPEC §7). It was implemented as fin13v4 (V4R), then extended as fin13v4b (V4BR).

**Where; what; fix (SPEC §7 rows marked corrected / new / deleted).** All rows below are fixed in fin13v4 (`V4/verifier/rhs_iv.py`, `_ei_v4`, NEW branch) unless stated otherwise; fin13v4b inherits them.

F1a. **Typed constants below their certified values** (AUDIT §4, GI §1.3; SPEC §7 rows 2, 15, 35, 39) [cert, exact decimal comparison]:

| constant (`V33/verifier/constants.py` line) | typed | certified raw | rel. |
|---|---|---|---|
| AF (23) | 6.9555e-2 | 0.069555037634 | -5.4e-7 |
| CF0 = c3 (24) | 0.4421 | 0.442142686616 | -9.65e-5 |
| CF2 = c1 (24) | 0.4235 | 0.423527657147 | -6.53e-5 |
| CF3 = c0 (24) | 0.4250 | 0.425017231881 | -4.05e-5 |
| DS1, DS2 (25) | 1.5418e-2, 2.3648e-2 | 1.54183748e-2, 2.36481789e-2 | -2.4e-5, -7.6e-6 |
| SB, SA, IA (22) | 3.5955e-5, 1.3372e-4, 9.5997e-12 | 3.59551663e-5, 1.33720281e-4, 9.59970542e-12 | -4.6e-6, -2.1e-6, -5.7e-7 |
| TQ d010, d110, d300 (28) | 1.67e-3, 1.086, 3.314 | 1.67114e-3, 1.086422, 3.3140088 | -6.8e-4, -3.9e-4, -2.7e-6 |
| TQ d100 (28) | 1.1e-4 | 1.129e-4 (ramp region only) | -2.6e-2 (far value 2.5e-5) |
| FN[0], FN[1] (30); C3 (16) | 7.93e-7, 2.769e-6; 12.104 | 7.93250e-7, 2.769321e-6; 12.104385 | -3.2e-4, -1.2e-4; -3.2e-5 |
| FNW[0] (31) | 1.09e-7 | 1.092e-7 | -1.8e-3 |

(Source: GI §1.3 table; the CF/AF raw values from AUDIT §4, `agent_coup/coupling_constants.json`.) The verifier's +-1 ulp on typed literals makes the double contain the decimal, not the certified value (AUDIT §4, `rhs_iv.py:14-16`).
*Effect:* no single item is large; but no typed number was certified. The TQ sups were computed for alpha <= 1e-3, so the true sups on the run are smaller (GI §1.3).
*Fix:* `V4/v4cert/out/constants_v4.json` (stage C5), 269 entries with (raw, value, dir, source, key) and an assert that every value is outward. The far TQ are replaced by the zone x level tables of stage C2 (lever L1, alpha-level 2.5e-4: e.g. CBAR 3.365e-3 against typed 2.27e-2) [cert, V4R]. Status: **fixed in fin13v4**.

F1b. **MU1 provenance** (SPEC row 40; GI 0-1). `constants.py:19` MU1 = 0.1552 is above the L-int certificate it cites (0.155167, `agent_lint/runs/cert_PX2_XA10.25_A230.25.json`) and valid only through `agent_int/v2/v2_T4_P2_0.25.json` (0.15552) [read]. *Fix:* fin13v4 reads MU1 = 0.1555226 (certified) from `v2_T4_P2_0.25.json`, rounded down (V4R). J0, J1 fall slightly because of it (V4R). **Fixed in fin13v4.**

F1c. **(N1) top-order cross terms of 1 - J^-4 missing** (AUDIT (N1); GI 0-3, 2-3; SPEC rows 3, 16).
- Level 0: `rhs_iv.py:147` codes d001 V2^2 only; needed in addition sup\|c' - 2psi'c\| V1V2 + sup\|psi''c + psi'c' - psi'^2 c\| V0V2 [hand; identity machine-checked by GI `check_identity.py`]. ASSEMBLY.md:78 lists "0.025 V1V2"; the code has nothing.
- Level 2: `rhs_iv.py:178` codes d001 V4^2 where the pairing gives d001 V4 \|\|D2v\|\|.
- *Size:* at the typed box d001 = 0.0227, ratio 0.239 of the level-0 E_I dissipation and 0.110 at level 2 [cert ratio of hand coefficients, GI §2.1]; 0.030 / 0.011 at the boxes reached [est input]. AUDIT estimated 10-20% [est] (see Conflicts, item 1).
- *Fix:* SPEC §4.3: level 0 CBAR V2^2 + (C1 + 2 PSI1 CBAR) V1V2 + (PSI1 C1 + PSI2 CBAR) V0V2; level 2 CBAR V4 DV + C1 V3 DV; the quasilinear top order is kept globally as (c v''')', c = 1 - J^-4. **Fixed in fin13v4.**

F1d. **(N2) beta-, psi1-, psi2-couplings absent from I0 and I2** (AUDIT (N2); GI 0-10, 0-11, 2-12; SPEC row 33). ASSEMBLY.md:80 lists them; `rhs_iv.py` has no beta [read]. *Size:* cross term kappa/2 up to 17% of rho in rate units at alpha = 1.5e-4; effect on the top eigenvalue 1.5e-6 [cert, GI §2.1]; the degree-0 forcing part 2.73e-11 at alpha = 1.5e-4 doubles the coded I0 forcing 2.84e-11 [cert, GI §2.1]. *Fix:* SPEC §4.9 T11 at both levels, with the strip part of beta in state form. **Fixed in fin13v4.**

F1e. **Inner ramp [10.25, 16.25].**
- Level-0 quadratic and cubic monomials d200 w''^2, d110 w''w''', d300 w''^3 not in term (a) (AUDIT 0.5b; GI 0-6; SPEC row 12, new).
- Term (a) W3, B, W1 coefficients low by 3.0x, 2.7x, 2.1x if the top order is done globally [hand, GI 0-5]; term (a) D4 entry replaced by V0 RQ (SPEC row 10); the commutator c K4 + c' K3 added (SPEC §4.4, row 24, new).
- The S-dependence used a secant (chord) between eta = 0 and 1e-2 (`constants.py:39-48`); replaced by the row eta = 1e-3 plus an assert (SPEC row 36; GI 0-7: the secant happens to hold there with 0.01% to spare, still not a proof).
- Level 2: quasilinear commutator coefficients (7.74, 7.74, 7.74, 7.73)e-5 of CF [cert, AUDIT §2]; the level-2 nonlinear block "not derivable as the code is organised" (GI 2-10; radial degree 2.25); SPEC §4.3-4.6 reformulates it (one factor converted to v).
- Level-2 modulation commutator -(aq + b cos theta) zeta' g: <= 1.2709 \|(a,b)\| W1 \|\|D2v\|\| [cert, AUDIT 2.6b]; coded as TC2 = 1.27089 (SPEC row 32, new; V4R C1). See Conflicts, item 4.
- **Fixed in fin13v4.**

F1f. **Level-2 far block** (AUDIT 2.5b, 2.5c; GI 2-4..2-9; SPEC rows 19-23):
- `rhs_iv.py:179` 2 d200 supv1 V3 -> 2 d200 sZ V2 (index; row 19, corrected);
- d300 supv1^2 V2 -> 3 d300 sZ^2 V2 (factor 3; row 21);
- d110 supv1 (3V3 + 1.31V2) -> d110 sG (3V3 + PSI1 V2) with the global sup sG (row 20);
- the else-branch d110 supv2 V2 (GN = 0) deleted (row 22, not valid beyond 34.25);
- rows with d^a, d^b, E1c, E0c new (row 23); Dd001 and c' terms (GI 2-4).
- **Fixed in fin13v4.**

F1g. **Modulation constants.**
- `c62` (62 Am V2V3, `rhs_iv.py:180`) **deleted** (SPEC row 26): no derivation; it was silently covering the deficit of the level-2 a-transport constant: certified need 39.1842 (Euclidean with b) against the coded 39.16 (AUDIT 2.6a). Now T2 = 39.18420 and Q2B 2.10e-4 typed honestly (V4R C1; SPEC row 25, corrected).
- `rhs_iv.py:144` am_phi FnW(al) -> PHIAB Fn(al) in the (a,b) bound (SPEC row 38; GI 0-14: FnW is 7.3x too small at small alpha, 66x in the alpha^3 term; `mod_sups.py:25` already used Fn).
- **Fixed in fin13v4.**

F1h. **Forcing.** Level 2 coded as 1 x \|alpha\|^3 \|\|D2v\|\| (`rhs_iv.py:181`), not an upper bound for \|alpha\| < 1.64e-4, i.e. on the whole run (needed/coded 1.20, 19.9, 5.7e4 at alpha = 1.5e-4, 3.68e-5, 6.9e-7 [cert, GI 2-11]); replaced by FNW2(al) DV (SPEC row 28). Level 0: FNW typed from the C^3 profiles, valid only as a sum; retyped with the C^4 triple (SPEC row 9; GI 0-9). **Fixed in fin13v4.**

F1i. **Small omitted terms** (SPEC rows 4, 13, 30, 31, 34): Garding e-terms at level 2 (\|e1\| <= 5.66e-6, \|e0\| <= 5.11e-6 [AUDIT 2.1]); c0'w at level 2; range of S_all in the c0'w term plus C0Hb (1.9e-42 [cert, V4R C4]); A0 w; the O(k) terms T12 (SPEC U-5: constants [est]). **Fixed in fin13v4** (T12 with \|cos\|, \|sin\| <= 1, V4R "remains unjustified" item 6).

F1j. **SF-3: `unc.lk()` was not a majorant of L** (SPEC §11 SF-3; row 37). `V33/verifier/unc.py:19-28` omits CABQ Amc al V0 (`rhs_iv.py:151`) and CAB03 Amc al D2v (`:182`); size tiny [est], but the J constants are certified numbers consumed by (A). *Fix:* `unc.lk_v4` linearises `rhs_iv.L_list`, the same list `rhs_pn` sums. **Fixed in fin13v4.**

F1k. **SF-4 / U-4: box hypotheses logged, not asserted** (`closure.py`). *Fix:* B-b, B-c, B-d (including \|w(34.25)\| + 25.75 S0V <= 0.5), strip S0..S2 <= 0.02 and the modulation box are asserted in `closure.box_assert` at every step and at the tail anchor. B-e (S0V, S1V from (A)) remains an input (V4R). **Fixed in fin13v4**; worst ratios in fin13v4b: 0.789, 0.773, 0.765, 0.747 (V4BR).

F1l. **SF-1: E0 had no term for the strip part of beta** (SPEC §11). *Fix:* + 7.52 BTs A in E0, BTs in state form (V4R). **Fixed in fin13v4.**

F1m. **U-3: the strip part of the modulation bounds as a linear form** (`V4/v4cert/u3_strip_form.py`, [cert/read]). Two corrections were found:
1. `linear_alpha_S0` multiplies \|alpha\| S0 (`ibp_run.py`), so the old STRIP_BETA counted it with alpha = 1;
2. the lemma omitted the Neumann feedback on the strip sources (`ibp_final.py:75-81`); v4 adds it.
*Effect:* STRIP_BETA 1.06e-9 -> 8.04e-11, STRIP_AB 5.52e-9 -> 4.37e-10 (`V33/verifier/constants.py:76-77` -> fin13v4), used in the alpha-equation and `mod_sups`; sup\|beta\| falls 29% (3.33811e-9 -> 2.35741e-9 at alpha1 = 1.5e-4). The fallback (constant kept, F0 x7.4) was not needed. **Fixed in fin13v4.** (Further change in F5a.)

F1n. **U-1 regularity lemma** for the energy identities (SPEC §11): not written in fin13v4/fin13v4b. Now `fix10_proofs/U1_regularity.tex`, proved [hand] given (H6) smoothness at t_1, the Lean theorems (read, not compiled, n = 1 not elaborated) and the DKS smooth-extension step (PRF). **Open** (independent review; inputs).

**Effect on the verifier outputs** (alpha1 = 1.5e-4, sqrt(Lambda_1) = 5e-8). [cert]; values are the stored doubles as printed in V4R and V4BR (7 significant digits; not outward-rounded by the reports; the consumers round them, V36R §3):

| quantity | fin12f4 | fin13v4 | fin13v4b |
|---|---|---|---|
| steps to the tail | 10 | 10 | 10 |
| peak sqrt(Lambda) | 5.307458e-8 | 5.307698e-8 | 5.308562e-8 |
| tail rho | 2.693042e-3 | 2.684366e-3 | 2.677666e-3 |
| sqrt(L_F) | 9.871704e-10 | 9.887644e-10 | 9.900007e-10 |
| alpha_F | 6.866519e-7 | 6.876746e-7 | 6.884239e-7 |
| J0 const / slope | 6.027237e-12 / 3.454743e-18 | 6.020438e-12 / 3.452361e-18 | 6.025378e-12 / 3.455292e-18 |
| J4 const / slope | 6.236181e-8 / 3.574508e-14 | 6.242096e-8 / 3.579469e-14 | 6.247218e-8 / 3.582508e-14 |
| sup\|(a,b)\| (6 digits up) | 2.19103e-7 | 2.14024e-7 | 2.14041e-7 |
| sup\|beta\| (6 digits up) | 3.33811e-9 | 2.35741e-9 | 2.35488e-9 |
| first cell rate U | +3.034406e-2 | +3.035649e-2 | +3.04299e-2 (printed +0.0304299) |

At alpha1 = 1.25e-4: rho 2.860963e-3 (fin12f4) -> 2.855715e-3 (fin13v4) -> 2.849387e-3 (fin13v4b) (V4R, V4BR).
Checks: regression with the new terms off reproduces the previous tag bit for bit; double runs byte-identical apart from timing; tail_cert PASS (V4R, V4BR).
**Status.** The coded terms of SPEC §5 are in fin13v4b. Still open: the written E_I lemma `ei_lemma_v4.tex` (not written; V4R "W"; FREV §7 item 5), an independent review of the v4b delta (FREV §7 item 6), and the items in §J.
Sources: AUDIT, GI, SPEC §7/§11, V4R, V4BR.

### F2. KREG[12,16.25] typed below its certificate

- **Where:** `V33/verifier/constants.py:17`, `KREG["12<x<=16.25"] = 3.0208e4` (Sobolev constant K_rho on (12, 16.25]).
- **What was wrong:** the raw certified value is 30208.006, so the typed 3.0208e4 is below it (V4R, C5 row).
- **Effect:** relative -2e-7 [est from the two printed values]; KREG enters the ramp sup S_ramp (`closure.py:35`).
- **Fix:** fin13v4 uses the arb value 1/rho(16.25) (V4R, C5: "(new) KREG[12,16.25] = 3.0208e4 < raw 30208.006, now the arb 1/rho(16.25)"). **Fixed in fin13v4.**
- Source: V4R (C5 row of "Spec items: status").

### F3. E0 energy block (fin12f4 / fin13v4)

**What was wrong** (E0A §0-§3; file:line in `V4/verifier/`):
- E0-1/E0-2: `constants.py:14` PI3 = 0.01607 and `rhs_iv.py:40` f572 = 5.72 are valid numbers, but their derivation in ASSEMBLY.md:49 is wrong (n3 = 7.52 is not \|\|chi psi3^\|\|; the certified norm is 7.6970 [cert, `gram7_g3.json`]; beta_12 dropped).
- **E0-3: beta_Lin missing** (`rhs_iv.py:221`, MODA = 25, `run_v4.sh:7`): the alpha-form needs Qfull 6.194 + n3_123 m_new = 29.527 [cert].
- E0-4: second-order modulation pieces, +1.48 alpha A^2 [cert]; E0-5: nonlinear drift absent, alpha-form 0.549 [est]; E0-6: Phi1 = 1.083 missing in the strip term `rhs_iv.py:237` (+4.8e-8 A^2 [cert]).
- **E0-7:** `constants.py:16` ALIN sq = 59.3 < 59.32715 and BB = 6.2 < 6.21990 (`alin_E0_spec_0.001.json`); effect first cell +8.4e-6, anchor +1.9e-6 [cert probe].
- E0-8: the alin form covers x <= 31.25 only; the strip [31.25, 34.25] is not derived (W1B about 8.7 > 7.9 if continued [est]).
- E0-9: MOD_ALPHA literal 6.19 < 6.194179 (dead; superseded by MODA = 25). E0-10: the other 7.52 literals are valid as bounds of n3_123 only. E0-11: H31, H32 parts of beta_E missing (negligible).

**The V4 probe value.** V4R ("What remains unjustified", item 1) replaced PI3, f572 by 0.01835, 5.937 (n3 = 7.52 times the lemma constants) and found the tail anchor rate -2.684e-3 -> **-2.85e-4** and the first cell +3.340e-2 [probe] (`V4/probe_E0/e0probe.json`), "above the -0.002 the tail lemma needs at this anchor".
How it was resolved:
- E0A §2, §4: the probe is mechanically right (reproduced exactly: -2.850049e-4, +3.340046e-2), but its premise is wrong: 0.01835 / 5.937 are neither the documented nor the correct constants. With w~ orthogonal to f1, f2, f3 only the projected norm enters: **n3_123 <= 6.5882** [cert] (S1); the rank-one symmetrisation S2 halves the beta_3 term. Pi3 needed = **0.010385** [cert + read] against 0.01607 coded; forcing factor needed 5.334 against 5.72 (E0A §2).
- E0A §4: -0.002 is the design threshold rho_min of `closure.py:153` (the march declines the anchor and goes on), not a lemma requirement.
- Fix F-B (E0A §6) was implemented in fin13v4b: PI3 = 0.01607, MODA = 25, f572 = 5.72 **kept**, justified by the E0 modulation lemma `V4B/e0_lemma/E0_LEMMA_v4b.tex` (V4BR).

**Fix and new numbers (fin13v4b; V4BR "E0 part")** [cert given reads unless stated]:
- A^2 budget: Pi3_true (S1 + S2) = 0.0103854, slack 5.685e-3. Alpha-form need at alpha <= 2.5e-4: 38.541 = Qfull 6.194 + n3_123 m_new 29.527 + E0-4 pieces 2.27 + **drift 0.549308 [cert]** (`e0_strip_drift.py`; was 0.549 [est]). Excess 3.385e-3, **margin x1.679** (x2.786 at alpha <= 1.6e-4). The lemma is stated for alpha <= AMAX_E0 = 2.5e-4, asserted.
- Forcing budget: 5.3439 <= 5.72 (x1.070) at alpha <= 2.5e-4 [cert], every Fn piece included.
- **n3_123 = 6.588186** [cert] at `rhs_iv.py:184` (beta_S coupling); `qb752` (beta_Q) stays 7.52 >= n3_123.
- Strip alpha-linear form (E0-8) [cert], coded in state form; W1B is now the raw core maximum 7.8967. Its weighted divergence form T_c^{rho chi^2} [sym].
- ALIN (E0-7) from raw `alin_E0_spec_0.001.json`, rounded up: sq 59.327153, AB 35.868736, W1B 7.896699, BB 6.219904.
- Phi1 in the strip commutator term (E0-6): 2 IA K_rho(16.25)(PHX A^2 + 2 A W1).
- **New terms, all coded with certified constants** (V4BR "New finding"): H3 (a,b) w-parts of A1(u), B1(u) (CA1 32.0, CB1 1.032); H2 cubic drift (CD1 0.607); H4 QAB Q and strip parts of (a,b) in Q1 and the loop (GABN 8.03, LOOPC 5.65); H5 alpha^2 remainders of A1(u0), B1(u0) (RA2 26.9, RB2 5.50); H6 beta beyond its lead in the h2/E33 terms (HOB 33.37); H7 H31/H32 part of beta_E. (Printed values as in V4BR; the outward values are in `V4B/v4bcert/out/e0v4b_constants.json` and the E0 lemma's table.)
- **New hypothesis:** \|w\| <= 0.02 on [0, 34.25], asserted (worst ratio 0.030).
- Machine checks `e0_identities.py` I1-I8 [sym]; the beta decomposition and the Cauchy-Schwarz forms of H3-H7 stay [hand]; the float check of the beta decomposition suggested by E0A was not done (V4BR "Open" 2).
- Effect (E0 part alone): costs 0.6% of rho and +0.18% of J0 (V4BR "Reading the numbers").
- **Status: fixed in fin13v4b.** E0_LEMMA_v4b.tex not compiled. Sources: E0A, V4R item 1, V4BR.

### F4. E2 energy block (fin12f4 / fin13v4)

**What was wrong** (E2A §0, §2; `V4/verifier/rhs_iv.py:243-251`, byte-identical in fin12f4):
- **E2-1:** `rhs_iv.py:250` `modE2 AA` = 43.33 A^2 bounds only the (a,b)-lead coupling (certified 43.3271 [read `agent_mod/mod_forms.json`]). Missing: the beta coupling, (a,b) beyond the 2x2 lead (f3 row 0.0131 A; alpha-linear 622.4 alpha A; quadratic 51.34 Q; forcing 3.84 Fn; strip), and the couplings through A1(u), B1(u). Added A Dq coefficient 0.1131 at the tail anchor, 0.3380 at alpha = 1.6e-4, 0.5025 at 2.5e-4 [cert].
- **E2-2:** E2 forcing `eps Fn Dq` (`rhs_iv.py:106`): needed factor **13.984** at alpha <= 2.5e-4 [cert], coded 1.
- **E2-3:** `rhs_iv.py:43` r143 = 14.34 < 14.343002687670673 (`garding_cert.json`); effect 5.2e-10 on the rate [cert probe].
- **E2-4:** 8 of nd1..nd8, nf1..nf7 typed below their own certificate (e.g. nd4 0.54 < 0.543880, nd8 3.05 < 3.052044).
- **E2-5, E2-6:** the certificate's jet box (\|w\|, \|w'\| <= 1e-4 on [0, 34.25], alpha <= 1e-3) asserted nowhere and false in the worst case (1.355e-4 on [0,12] [cert]); at alpha level 1e-3 four constants fail (up to x1.45).
- **E2-7:** strip part of chi N(w) treated as if chi w^(j) = w~^(j); needs SRW = 4.4292e-3 [cert] W-norm terms.
- E2 audit size (probe, Cauchy-Schwarz form): first cell +0.03035649 -> +0.03036965, anchor -0.00268437 -> -0.00268121 [cert given the constants].

**Fix (fin13v4b; V4BR "E2 part"):**
1. RM + NS terms times Dq exactly as E2A §5 F2; modE2 stays 43.33 (>= 43.3271).
2. **R143 = 14.343002687670673**, read from JSON and rounded up.
3. FE2 = 1 + BETAF (n3 + 2 alpha \|\|h2\|\|) + PHIAB (NF12 + u0pr) multiplies the three E2 entries of L_list; about 13.98 at alpha = 2.5e-4 (V4BR prints "about 13.98"; E2A gives 13.984).
4. Assertions: E2_alpha <= 2.5e-4 (ratio 0.64); E2_tip_jet sup(\|w\|, \|w'\|) <= 1e-4 on [0, 8] with K = KREG[4,8] (ratio 0.324).
5. nd = max over [0,8] (box 1e-4) and [8,12]; nf = max over [12,16.25] and [16.25,34.25], from the recertified JSON, rounded up; every typed value covers its new value (e.g. nd4 0.54 -> 0.205).
- Effect: the E2 part gives back +0.37% of rho and -0.10% of J (V4BR).
- **Tip-jet hypothesis** (H-E2b) is new and asserted; the Cauchy-Schwarz forms of RM stay [hand]; `e2_identity_check.py` gives the D_rho expansion [sym] and a float check.
- LAMT: no change; the verifier's -0.2561532 is 2.2458e-8 less negative than the residuals' -0.2561532224577463, the safe direction (E2A §6).
- **Status: fixed in fin13v4b.** Sources: E2A, V4BR.

### F5. Strip hypothesis on w''' and the tail alpha-condition (REV4 R-2, R-3)

F5a. **R-2: strip hypothesis sup\|w'''\| <= 0.02 on sigma in [29, 32]** (SPEC §11 SF-2; REV4 R-2).
- Where: used in STRIP_BETA, STRIP_AB (`V4/verifier/closure.py:22, 157`, `mod_sups.py:24-25`). Not implied by Lambda (needs V3, i.e. V4). REV4: S3 share 4.92% of STRIP_BETA [cert], 0.17% of sup\|beta\|.
- Fix (fin13v4b; `V4B/v4bcert/r2_strip_ibp.py` [cert]): one integration by parts of the 4 chi' w''' commutator term (chi' vanishes at both strip ends); the S3 coefficient is then 0. **STRIP_BETA 8.040e-11 -> 7.787e-11, STRIP_AB 4.368e-10 -> 4.227e-10**; these now need only S0, S1, S2 <= 0.02 (asserted). Used in the alpha equation, the tail and mod_sups.
- **Status: fixed in fin13v4b.**

F5b. **R-3: the tail alpha-condition was checked at the anchor only** (`V4/verifier/closure.py:156-158`), while STRIP_BETA/sqrt(l) grows as l falls (REV4: 0.0813 at l = L_F against margin 177 [hand]; inherited from fin12f4).
- Fix: `closure.tail` checks rhs_all along the whole tail (monotone pieces). **Margins 177.08 (alpha1 = 1.5e-4) and 175.10 (1.25e-4)**, stored as `tail.alpha_condition_whole_tail`; cross-check `V4B/verifier/t_tail_alpha.py`, min margins 177.11 / 175.13 [num-arb]. **Fixed in fin13v4b.**

F5c. **Text conflict, not resolved here.** `V4B/e0_lemma/E0_LEMMA_v4b.tex`, last paragraph ("Not covered here", lines 99-101), still lists "the strip hypothesis sup\|w'''\| <= 0.02 in STRIP_BETA/STRIP_AB (V4_REVIEW R-2)" as not covered. V4BR says R-2 removed that hypothesis, and the same file's hypothesis list (lines 28-29) names only S0, S1, S2 <= 0.02. The tex and the report disagree; the tex needs a new version. **Open (text).**

Sources: REV4 R-2/R-3, V4BR "V4_REVIEW items".

### F6. R3: honest forcing constant of Lemma S(A); c_dU of Lemma U

- **What was wrong:** the forcing constant of (A) and the delta-U lemma assumed lambda~ = -p/4, i.e. alpha_tau = -(p/4) alpha + beta. The true lambda~ = -0.2561532224577463 differs from -p/4 (p = 1.0246) by DL = max(\|lambda~ + p/4\|, \|LAMT + p/4\|) = 3.2224577463e-6 (exact) (R3R §2).
- **Where:**
  - `t2_A_v35/rigor/assembly_rep.py` (E_j forcing, export C1, `f4b_honest` c_F), `base_ibp/rigor/assembly_arb.py` (`window_loop`, `f4b` c_F and (B)-exit) with r0 = 1.4(betabar + 1.02 alpha_1 abar) (R3R §3);
  - `t4_cert/u34_sweeps.py` c_dU with betabar;
  - texts: `fix1/lemmaU_rev/lemmaU_v34.tex:55-57` (Lemma dU hypothesis "alpha_tau = -(p/4) alpha + beta", value "<= 1.7187e-8"), `lemmaS_A_v33.tex` step 6 (r0^F), `numbers_A_v33.tex`;
  - closed text: `writeup/step3_proof_v32/step3_proof_v32.tex:235` (R3R cites "P32:235-236"): "d/dtau[alpha lambda^-p] = (beta - p a alpha) lambda^-p (lambda~ = -p/4)". Not edited (closed).
- **Fix** (R3R §2): r0' = c_U (betabar + DL alpha_1) + c_U p alpha_1 (abar + bbar/14); F_j <= r0' p^j sqrt(I_int) + sqrt(Feps_j^2) with the eps part bounded cellwise (no series); c_dU = c_U (betabar + DL alpha_1 + p abar alpha_1) / (p (1/4 - abar)) with hypothesis alpha_tau = lambda~ alpha + beta.
- **Effect on fin12f4 (v36pre, `fix4_r3/`)** [cert]: r0' = 5.03643e-9 (old 4.72029e-9, x1.066975); F_j (A, [14,700]) x1.06961, x1.07464, x1.08923 (j = 0, 1, 2); Ehat_j up by at most x1.000394; **c_dU 1.71860e-8 -> 1.96498e-8** (x1.14336), room against the typed bound 2e-8 x1.16373 -> x1.01782 [float of cert; `c_dU < 2e-8` asserted in arb]. Largest relative margin loss x0.99944. No margin fails.
- **Effect in v36 (fin13v4b + R3 + r_2; `fix6_v36/r3/r3_forcing.json`)** [cert; JSON balls, digits rounded outward here]:
  - r0_new <= 3.7516343e-9 (ball [3.7516342034570612758e-9 +/- 1.85e-29]; r0_old 3.3426796e-9, ratio 1.122343351);
  - F_j^2 (A, [14,700]) j = 0, 1, 2: <= 3.5931248e-13, 3.8197758e-13, 4.1572077e-13 (balls 3.59312473631829e-13, 3.81977579632679e-13, 4.15720760970261e-13, radii < 4e-28);
  - **c_dU <= 1.4634276e-8** (ball [1.4634275344469975636e-8 +/- 3.00e-28]), betabar' = 2.83824866194500e-9; typed bound 2e-8 / c_dU = x1.36665 [float of cert] (V36R §4).
  - (V36R prints r0 3.75163e-9 in the work order and F_0^2 3.59312e-13: truncated midpoints, see Conflicts, item 12.)
- **Status:** fixed in v36 (`fix6_v36/`, double run byte-identical). Texts not edited: lemmaU_v34 Lemma dU (hypothesis and value; v36 value 1.46343e-8), lemmaS_A_v33 step 6, numbers_A_v33 (R3R §7 item 2; V36R §7; FREV §7 item 9). **Text open.**
- Sources: R3R §2-§3, §5, §7; V36R §3-§4; `fix6_v36/r3/r3_forcing.json`.

### F7. Outer-ramp flux lemma ([57, 60]); hypothesis (I3)

F7a. **Flux lemma v1 (fix round 3) not closed at level 2** (FLX1):
- Level 0 closed: K0 <= 3.9038e-21 at beta_A = 2 against FLUX0 = 1.64e-20 (x4.201) [cert].
- Level 2, pointwise part: P2 <= 5.8864e-20 against FLUX2 = 8.0e-20 (x1.359; was quoted x5.956) [cert].
- Not bounded: T3 (g'''^2, coefficient nonzero for every cutoff; 12.488 N3^2, N3 time-integrated only) and O4 (fourth derivative of the septic cutoff, zeta''''(57+) = -280/27, needs z4 data). No time-integrated consumer existed in the verifier (`grep I3` empty) [read]. Unbounded pieces [est] about 1e-21 and 1e-25 against 8e-20.
- Route decision (ROUTE): time-integrated absorption lemma (a) plus a level-2 export r_2 = 1/1000 from (A); no verifier or step-4 rerun.

F7b. **Defects of the v3.1 / v1 flux block** (EF, items E-F1..E-F5):
- **E-F1 cutoff shape:** the v3.1 flux tables used a quintic cutoff; the verifier's cutoff is the septic S7((60 - x)/3), C^3, with zeta^(4)(57+) = -280/27. v2 tables are generated from S7 (`flux2_coef.py`).
- **E-F2 jet box:** v3.1 modelled the commutator with quasilinear coefficients and an a-priori jet box (M2..M5); v2 uses the commutator of L1 with j = J^-4 inside it, j in [1 - cbar, 1 + cbar], \|j'\| <= cbar_1, no jet box.
- **E-F3 drift double use:** F1 (v1) used the drift sign for its one-sided diagonal and also to absorb the outer-ramp transport term (ROUTE §4). v2 uses the sign once, with the modulation transport inside the same ball C_eff in C +- abar(\|Q_u\| + \|C_u\|).
- **E-F4 trace route:** the endpoint route through g^(3)(57) does not close K2^OR (residue (5/2) g''(g^(3))^2, `flux2_forms.py` (f)); v2 uses the level-2 export of (A) as a space-time L^2 bound of g^(4) on [57, 60] (this departs from ROUTE §1 item 3).
- **E-F5 W(57):** 1.43716e-18 (`lemmaS_cflux.json`) belongs to the quintic blend; for the verifier's S7 blend W(57) <= 1.3404e-18 (FREV: 1.340304797627754e-18 [cert], `flux2_coef.py:322`). The (A)-chain flux block (records only, also in v36) still uses the old value.

F7c. **Level-2 modulation bound had no ramp term** (ROUTE §4; FLX1 §4 item 4): `V33/verifier/rhs_iv.py:178-180` has no ramp term for either ramp at level 2 [read]. Inner ramp: fixed in fin13v4 (F1e, TC2). Outer ramp: inside flux lemma v2 (E-F3).

F7d. **Flux lemma v2 numbers** (`fix8_flux2_impl/out/flux2_cert.json`) [cert, outward digits from the JSON]:
- Kbar0 <= 3.4098535e-21 (up 3.40985343319526e-21); K2^pt <= 1.3900119e-19 (up 1.3900118221562857e-19);
- combined margin (FLUX0 + eps_I FLUX2)/(Kbar0 + eps_I K2pt) >= x4.8044192 (dn 4.804419228679062);
- K0 vs FLUX0 >= x4.8095908; slack s >= 1.0390701e-21; C' <= 2.6523207e-26; S' <= 1.0551032e-28;
- conditions C0-C7 and C9 pass; C8 (regularity U-1) "OPEN, not checkable by script".
- FREV reproduces all of these; Kbar0 is conservative by x3.4 (direct pair form 1.00289e-21, combined margin x16.29).

F7e. **The hypothesis must be the combined (I3')** (FREV §7 T4). `K2pt_vs_FLUX2_record` = 1.737514777695357 (flux2_cert.json; FREV rounds it to 1.74): K2^pt alone exceeds FLUX2. So **any sentence "K2 <= FLUX2" is false for the new K^OR**. The E_I lemma must assume mu(K0 + eps_I K2) <= mu(FLUX0 + eps_I FLUX2), pointwise plus the time-integrated part through the absorption lemma (C' channel). The verifier consumes only the sum F0 = mu(FLUX0 + eps_I FLUX2) (`closure.py:144`; ROUTE §2.1), so the coded runs are unaffected.
- **Status:** v2 certified, reviewed "sound modulo listed items" (FREV §0). Open: U-1 review, text items T1-T4 (FREV §7), the written E_I lemma with (I3'), pinning of `fix8_flux2_impl/out/*` (in no manifest or kit).
- Sources: FLX1, ROUTE, EF, FREV §5-§7.

### F8. V36_REPORT rounds two design-setting quotes down

- **Where:** `fix6_v36/V36_REPORT.md` §6 table (and its verdict line: "C0_2 = 114.78, C1_2 = 0.0698").
- **What:** C0_2 = 114.782 and C1_2 = 0.0698142 are quoted; the certified values are 114.78222... and 0.069814238... (PRF item 7), so the quotes are below the bounds.
- **Effect:** none; nothing consumes the design setting. A_interface's `C0_ref`/`C1_ref` are correct upper bounds, and flux lemma v2 uses the 4 N2_ref values C0_2 <= 116.3782, C1_2 <= 0.2606004 [cert, PRF].
- **Status:** record only; quote 114.7823 / 0.06981424 (or the A_interface values) in any text. Source: PRF "Gaps and findings" item 7.

### F9. Text defects in the v3.3 texts found in fix rounds 9-10 (and earlier reports)

None of these texts was edited. All **open** in the v3.3 texts. (While this erratum was written, other agents were writing v37 texts in `fix12_texts/` (`lemmaS_A_v37.tex`, `lemmaU_v37.tex`, `continuation_v37.tex`, `numbers_A_v37.tex`, `main_v3_v37.tex`, ...). This erratum did not read them; whether they fix these items is not checked here.)
1. `fix1/b12_text/lemmaS_A_v33.tex:106` "r_2 = 0"; `:302` "Level 2 exports nothing (r_2 = 0)"; `:199` core B_c = (11/20, 9/20, 169/100), retained (1.2500, 1.2500, 0). v36: B_c^(2) = 1688999/10^6, retained 0.001001 >= r_2 = 1/1000 [cert, V36R §6; PRF item 6].
2. `fix1/lemmaU_rev/lemmaU_v34.tex:150, 190` Lemma D stated for k in {0, 1} (r_k = 1.249). It holds verbatim for k = 2 with r_2 = 1/1000; proof `fix10_proofs/G1_export_k2.tex` [hand + cert] (PRF item 6).
3. `lemmaU_v34.tex:55-57` (and `:478`) Lemma dU: hypothesis alpha_tau = -(p/4) alpha + beta and c_dU <= 1.7187e-8. Corrected hypothesis alpha_tau = lambda~ alpha + beta; v36 value 1.46343e-8 (F6).
4. **(GS) missing** from the hand-off list and from the continuation bootstrap: (GS) = the tip arc is a normal graph over S1 on \|x\| <= 60 at tau_1. It is implicit in step 4's certification of Lambda(tau_1). The hand-off list should state it, and `fix1/b12_text/continuation_v33.tex` should add (GS) next to (G); the open/closed argument is in U1_regularity.tex (PRF item 4).
5. `V4B/e0_lemma/E0_LEMMA_v4b.tex:99-101`: stale "not covered" sentence on sup\|w'''\| (F5c).
6. `V4B/e0_lemma/E0_LEMMA_v4b.tex` "E2 couplings" paragraph prints d_ab = 0.0131136; E2A gives DAB raw 0.013113620982 (and "dab = 0.0131137" in its §1), so the printed value is rounded down in the last digit [read]. Display only; the code uses the JSON value (V4BR).
7. Flux lemma v2 text `fix8_flux2_impl` `flux2_lemma_v2.tex:182` still says "Status. Design. Open ..."; `:70` "To be asserted in sympy"; `:3` describes placeholders (FREV §7 T1); march first-touch argument only as a docstring (T2); bootstrap set must list (BPhi), abar = 1e-6 and boxes B-a..B-e (T3); combined hypothesis (I3') (T4, F7e).
8. Closed: `writeup/step3_proof_v32/step3_proof_v32.tex:235` identity with lambda~ = -p/4 (F6). Not edited (closed); this erratum is the correction.

---

## G. Status (v37, Nov 20, 2026)

**Step 3 is not closed.** Nothing in Part 2 is called closed.

- **v3.1 / v3.2:** as in §E. In addition, the fin12f4 verifier certificate is invalid as written (F1-F4); fin12f4 is a copy of fin11 with changed flux inputs (§B), so the same presumably holds for fin11 [not re-checked here]. Therefore every closed-version statement that rests on them (including §E's "×1.1245" flux line and the fin12f4-based v3.3 column of §E) is a record, not a certificate.
- **Verifier:** fin13v4b (`V4B/runs/P`) closes at alpha1 = 1.5e-4 and 1.25e-4 with every term of SPEC §5 and the E0/E2/R-2/R-3 fixes [cert, V4BR]. Conditional on: the written E_I lemma (missing), U-1 (proved [hand], not reviewed), the [read] inputs and [hand] steps of V4BR "Open", hypothesis (I3') (F7e), B-e (S0V, S1V from (A)).
- **(A) chain, Lemma U, zone III:** v36 (`fix6_v36/`, on fin13v4b + R3 + r_2) closes; no margin fails; thinnest margins ramp M1 x1.18446, ramp z1 x1.20115, zone III b1 x1.2389, H6 far K3 x1.8609, Ghat22 x1.9164; c_dU x1.367 under 2e-8 [cert / float of cert, V36R]. No independent review of v36 on record (FREV §7 item 8).
- **Flux lemma v2:** certified, combined margin x4.8044 [cert]; reviewed sound modulo listed items (FREV).
- **G1 (Lemma D, k = 2):** proved [hand + cert] (PRF); not reviewed.
- **Step 4 / kits / manifests:** not re-pinned. Stages 8, 9, 18, 19, 20 of manifest_v33, a new final_tail driver that reads `EXCEPTION_v32_*.json`, a new manifest, and both kits (V36R §7). The closed FINAL_manifest is unaffected.
- **Texts:** F9 items, `ei_lemma_v4.tex`, main_v3, numbers_A_v33, continuation_v33, V33_INDEX successor; LaTeX compile of every text.

---

## H. Superseded records

These outputs or quotes are kept as records only; do not cite them as certificates.

| record | where | superseded by | why |
|---|---|---|---|
| FLUX0/K0 at beta = 2 (tight, N2ref) 3.40048; FLUX2/K2 6.01860 (and the v35/R3 columns 3.39666/3.39535, 6.01283/6.01070) | V36R §4 table | `fix8_flux2_impl/out/flux2_cert.json`: Kbar0 <= 3.4098535e-21, K2pt <= 1.3900119e-19, combined margin >= x4.8044 | old flux block: old kappa (quasilinear commutator, v3.1 form), old W(57) = 1.43716e-18; the separate "FLUX2/K2" line is meaningless for the new K^OR (F7e) (FREV §7 T4) |
| K0 <= 3.9038e-21 (x4.201), P2 <= 5.8864e-20 (x1.359) | FLX1 (`fix3/flux/flux_cert.json`) | flux lemma v2 (F7d) | v1 not closed at level 2; v1 inherits E-F3 |
| flux line "FLUX0/(2K0) = x1.1245" and the §E v3.3 flux columns | §A8, §E (Part 1) | v36 + flux lemma v2 | computed with the old flux block and fin11/fin12f4 J's |
| fin12f4 verifier outputs: rate tables, modsups, tail logs (`V33/verifier/rate_table_fin12f4_*`, `modsups_fin12f4`) | `step3nl_v33/` | fin13v4b, `V4B/runs/P` | E_I inequality not an upper bound (F1); E0/E2 defects (F3, F4) |
| fin13v4 outputs (`V4/runs/P`) | `step3nl_v4/` | fin13v4b, `V4B/runs/P` | E0 and E2 blocks not certified bounds as coded (F3, F4); R-2, R-3 (F5) |
| (A) chain v35 (A_interface 78fe2f8e..., U_cert_v35 960809fd..., zone3 45e8d4dc...) | `fix2/n12/t2_A_v35/`, `t4_cert/`, `t5_zone3/` | v36 (A_interface 799b087c..., U_cert_v36_pass2 d512596a..., zone3 d6813ddb...) | consumes fin12f4; r0 with lambda~ = -p/4 (F6) |
| v36pre (R3 on fin12f4) | `fix4_r3/` | v36 | consumes fin12f4 (R3R §7 item 1) |
| c_dU 1.7187e-8 / 1.71860e-8 / 1.96498e-8 | lemmaU_v34, v35, R3R | 1.46343e-8 (v36) | R3 and fin13v4b |
| V4 probe -2.85e-4 (anchor), +3.340e-2 (first cell) | V4R item 1 | E0A §2, §4; fin13v4b | premise wrong (F3) |
| STRIP_BETA 1.06e-9 / 8.04e-11; STRIP_AB 5.52e-9 / 4.37e-10 | fin12f4 / fin13v4 | 7.787e-11 / 4.227e-10 (fin13v4b) | F1m, F5a |
| W(57) = 1.43716e-18 (also "<= 1.44e-18" in AUDIT, FLX1, ROUTE) | `lemmaS_cflux.json`; (A)-chain flux block | W(57) <= 1.3404e-18 | E-F5 |
| C0_2 = 114.782, C1_2 = 0.0698142 | V36R §6 | 114.78222..., 0.069814238... (design; unused); 116.3782, 0.2606004 (used) | F8 |

---

## I. Conflicts between reports

Listed, not resolved, unless a later report resolves it (then the resolution is named; nothing was chosen silently here).

1. **Size of (N1).** AUDIT: about 10-20% of the E_I dissipation [est]. GI: 0.239 (level 0) and 0.110 (level 2) of the dissipation at the typed box [cert ratio of hand coefficients]; 0.030 / 0.011 at the boxes reached. GI §4.1 notes the disagreement. Both are pre-fix sizes; fin13v4 codes the term.
2. **Verdict on G0.** AUDIT: KILL. GI: FIX-AND-RERUN, disputing "KILL" if read as "plausibly fatal" (GI §4.2 item 1). AUDIT's KILL is a process verdict ("the verifier may well close after the fixes"). Resolved in practice by fin13v4/fin13v4b closing [cert].
3. **Inner-ramp monomials.** AUDIT: "~100x the coded term-(a) W3 coefficient" [est, s_ramp ~ 1e-2]. GI: "about 50x" [est, sup\|w''\| = 5e-3] (GI §4.2 item 3). Different assumed sup\|w''\|; both [est].
4. **Level-2 transport commutator constant.** AUDIT: <= 1.2709 \|(a,b)\| W1 \|\|D2v\|\| (table 2.6b) and "1.271" in the text; V4R C1: TC2 = 1.27089. The printed values differ in the 5th digit; whether 1.27089 is a rounded-up bound or the AUDIT's 1.2709 is (V4R states every C5 value is outward-asserted) is not visible from the reports.
5. **kappa on the outer ramp, sup(\|zeta'\| + psi' zeta) on [57, 60].** AUDIT: 1.5883 [cert]; V4R C1: KAP OR 1.58832 [cert]; GI §4.2 item 5 and `check_kappa_outer.py` docstring: 1.5833 [cert], while calling the audit's 1.5883 "right". 1.5833 < 1.58832, so the two cannot both be upper bounds of the same sup; probably a digit slip in GI, not checked here. No consequence (1.67 covers both).
6. **"First march cell" of fin12f4.** GI §2.2: U = +0.0271016 (alpha <= 1.5e-4, Ly about 2.5e-15). SPEC §0 and V4R: first cell U = +3.034406e-2 (V4R: "first cell U = 0.03034405538749633"). Both are called the first march cell of fin12f4; the cells probably differ (GI's cell definition is in `fix5/g0_indep/cells.jsonl`). Not reconciled by any report.
7. **E0: constants V4 proposed.** V4R item 1: "the right values are 0.01835 and 5.937". E0A §2: these are "neither the documented nor the correct values"; correct: Pi3 0.010385, forcing 5.334 (5.325 + extras). Resolved by E0A and V4BR (constants kept, lemma added).
8. **E0 A^2 budget and forcing.** E0A §3: need 38.30 (37.204 [cert] + 1.1 [est]); ratio x2.67 at alpha <= 1.6e-4, "x1.61 at 2.5e-4"; forcing 5.334 (5.3331 + extras, alpha <= 1.6e-4). V4BR: need 38.541 at alpha <= 2.5e-4 (drift certified 0.549308; E0-4 pieces 2.27 instead of 1.48); margin x1.679 at 2.5e-4, x2.786 at 1.6e-4; forcing 5.3439 at alpha <= 2.5e-4. V4BR is later and certified; E0A's numbers are superseded. Note that E0A's own "x1.61 at 2.5e-4" does not follow from its slack 5.685e-3 and need 38.30 (that gives about x1.71 [est, my arithmetic]); not explained in E0A.
9. **R-2 cost.** REV4: the S3 term is 4.92% of STRIP_BETA [cert]. V4BR: after the IBP, STRIP_BETA 8.040e-11 -> 7.787e-11 (-3.1%). Not a contradiction (the IBP adds a w'' coefficient), but the reports do not reconcile the two figures.
10. **E0_LEMMA_v4b.tex vs V4BR on the w''' hypothesis** (F5c): the tex's last paragraph says not covered; V4BR says removed. Text conflict, open.
11. **G1 status.** V36R §6: "G1 ... has not been done"; FREV §7 item 2: G1 open. PRF (later): G1 proved [hand + cert]. Time order, not a disagreement on substance; the proof is not yet reviewed.
12. **Printed values below the certified bounds (display rounding).** V36R §4: F_0^2 = 3.59312e-13 (ball midpoint 3.59312473...e-13); the work order's r0_new 3.75163e-9 (ball 3.75163420...e-9); FREV §5 heading "W(57) = 1.3403e-18" (certified <= 1.340304797627754e-18); FREV §6 / work order K2^pt 1.39001e-19 and Kbar0 3.40985e-21 (JSON up values 1.3900118e-19, 3.4098534e-21); E0_LEMMA_v4b d_ab 0.0131136 (raw 0.013113620982); V36R §6 C0_2, C1_2 (F8). None is consumed by a script; this erratum quotes the outward values.
13. **ROUTE trace route vs implementation.** ROUTE §1 item 3 chose the endpoint (4,0) term via the trace g'''(57) and Agmon; EF E-F4 says the trace route does not close K2^OR and v2 uses a space-time L^2 bound of g^(4) instead. FREV §5 accepts the departure.
14. **P2 margin.** FLX1 treats P2 <= FLUX2 (x1.359) as the level-2 requirement; ROUTE §2.1 says the real condition is the combined one (x4.2) and "x1.359 is not a real constraint"; FREV T4 makes this the hypothesis (I3'). The FLX1 framing is superseded.
15. **W(57).** AUDIT §5, FLX1, ROUTE: W(57) <= 1.44e-18; EF/FREV: 1.43716e-18 belongs to a blend the verifier never used; the true value is <= 1.3404e-18. The old figure is still an upper bound, so the old statements stay true but use the wrong source.

---

## J. Open items (Nov 20, 2026)

From V4BR "Open":
1. [read] inputs of the E0 and E2 lemmas: gram7/gram8 Grams, pair_f3, bloop.json (Qfull), final_5e-4.json (m_new, BETA*, QB, CAB0, MAB, QAB, PHIAB), lyap_consts2, mod_forms Q2 (43.3271), garding_cert, coupling constants, profile hulls.
2. [hand] steps: the beta decomposition and its identification with coded terms; the Cauchy-Schwarz forms of H3-H7 and of E2's RM/NS; the float check of the beta decomposition on random w (not done).
3. Not audited: the theta-certificate (1/32, 1/9), `cubic_E0fine` and ADD_STRIP_CUBIC = 2.6e-3, the Q2 Gram, the Garding identity behind mu2/e1/e0, the level-2 commutator constants DS.
4. Outer-ramp FLUX0/FLUX2 for the new K^OR: now certified by flux lemma v2 (F7d) but not wired as a written hypothesis (I3'); U-1 regularity lemma (now [hand], F1n); B-e (S0V, S1V inputs from (A)); `ei_lemma_v4.tex` not written; legacy typed constants MOD_E2, t22/t12, qb752, LAMT (valid per the audits).
5. `E0_LEMMA_v4b.tex` not compiled (and F5c, F9 item 6).
6. Downstream re-pinning on fin13v4b: done for the (A) chain (v36); step 4, kits, manifests not done.

From FREV §7:
7. U-1 (C8), G1 written proof: now in `fix10_proofs/` [hand], not independently reviewed (PRF).
8. Text items T1-T4 of flux lemma v2 (F9 item 7); B1 -> Lambda <= L + 2C'; tail floor quotes with (1 + 2C'/L_F)^(1/2).
9. Pinning: `fix8_flux2_impl/out/*` in no manifest and no kit.
10. No independent review of the v4b delta (E0 F-B, E2 fix, R-2, R-3).
11. Independent review of v36.
12. Texts: lemmaS_A_v33 step 6 (R3, r_2), numbers_A_v33, lemmaU_v34 Lemma dU, continuation_v33, main_v3, the erratum (this file), index.
13. Step 4: re-pin stages 8, 9, 18, 19, 20; new final_tail driver reading `EXCEPTION_v32_*.json`; new manifest. Kits rebuilt and rerun; independent Windows runs.
14. Author items: (P1') local existence with uniqueness (Lean CWFlow check, or Giga-Ito / Escher-Mayer-Simonett); (P2') (proved [hand] in U1_regularity.tex per PRF); (P3'); Lemma O arc by arc; the single-point statement.
15. LaTeX compile of every text; final independent review of the assembled step 3.

From PRF "Gaps and findings":
16. (H6) smoothness at t_1 is step 4's [A] item (`step4_proof.tex` l.672), not proved in U-1.
17. Lean: build not run here; the n = 1 instance not elaborated; gluing and closure of `IsFlow` are written mathematics, not Lean.
18. DKS (D1) "bounded int kappa^2 => smooth extension" is a step inside the proof of Thm 3.1, not its statement.
19. (GS) at tau_1 missing from the hand-off list and the continuation bootstrap (F9 item 4).
20. Read, not re-derived: graph equation (`graph_eq.py`, `graph_struct.pkl`); modulation matrix = G + E and the Neumann bound nu <= 0.0466 (MOD_IBP.md); the C^4 matching of psi3^, h2^ (`c4fun.py`).
21. G1 (s5): F_max from the verifier's peak sqrt(Lambda); that Lambda <= L + 2C' stays inside it is flux-lemma check C7, not redone in PRF.

Also open from earlier parts of this erratum: (O3-v35) (§E), the A7 (B)-window convention and the A5 I_3 consumer for the closed versions.

## Files (v37)
- `fix12_texts/ERRATUM_v37.md` (this file) and `fix12_texts/erratum_v37.tex` (same entries and numbers, slightly more compact; brace and environment balance checked by script; not compiled, no LaTeX on this machine).

---

# Part 3 (v38, November 20, 2026): sections K-N

Sources of this part: `fix11_review/REVIEW_v4b_v36_proofs.md` (FREV11), `fix15_review/REVIEW_v4c_interface.md` (FREV15), `step3nl_v4c/V4C_REPORT.md` (V4CR), `fix13_v37chain/V37_REPORT.md` (V37R), `fix14_interface/INTERFACE_REPORT.md` (IFR), and the number check of this text set (`fix16_texts_v38/tools/numbers_check_v38.json`, NC38). Nothing closed was edited; texts v37 and earlier are unchanged.

## K. Defects and items found Nov 20, 2026 (fix rounds 11-15) and by the v38 number check

### K1. ADD_STRIP_CUBIC was underived (FREV11 D-3, blocking at the time)
- **Defect.** The strip part of the cubic E-block term was coded as `ADD_STRIP_CUBIC` = 2.6e-3 times the far-region cubic term; this scaling was not derived and is not a bound.
- **Fix.** Verifier fin13v4c: term H8 with certified coefficients C3S_abc (`step3nl_v4c/v4bcert/out/strip_cubic_v4c.json`, arb, keys asserted) [cert]. FREV15 recomputed the coefficients independently: coded/recomputed 2.8 to 17, all covering [cert].
- **Effect.** Favourable: tail rho 2.677666e-3 -> 2.679225e-3, J x0.99984 (V4CR) [cert].

### K2. E0 A^2 budget omitted two pieces (FREV11 A-1)
- The (a,b) lead mismatch, 0.179383 per alpha A^2, and the third-order loop piece, 2.641e-3 at alpha_E, were missing.
- Certified with them (`v4bcert/out/e0_A2_budget_v4c.json`): need 38.7232, margin x1.65694 at |alpha| <= 2.5e-4 (x2.747 at 1.6e-4). V4BR and text v37 printed x1.679 (without A-1). The verifier constants are unchanged [cert].

### K3. The cubic_E0fine box was not asserted (FREV11 D-2a)
- The box |w|, |w'| <= 0.01 on [16.25, 31.25], |alpha| <= 1e-3 behind `cubic_E0fine_0.001_0.01_0.01.json` was a hypothesis of the coded constants but not asserted. fin13v4c asserts it at every march step and at the tail anchor; worst ratios 0.0591, 0.116, 0.16 [cert, read].

### K4. The whole-tail alpha-condition "margin" was a difference (FREV11 A-2; texts)
- v37 texts printed "margins x177.07 and x175.09" (`main_v3_v37.tex`, `lemmaU_v37.tex` Lemma H-J-W and Table, and the macro `\tailAlphaMarV` in `continuation_v37.tex` (W1) and `step3_v37_main.tex`). These numbers are lhs - rhs, not lhs/rhs.
- The ratio lhs/rhs_whole_tail is x1554.45 (alpha_1 = 1.5e-4) and x1512.16 (1.25e-4), from `rate_table_fin13v4c_a*.json` key `tail.alpha_condition_whole_tail`, rounded down; with the factor 1 + 2C'/L_F it is condition C5 of the flux lemma (same values). In v38 `\tailAlphaMarV` is the ratio; the difference is `\tailAlphaDiffV`.

### K5. Step 4's two [A] sentences and the identification of the step-4 solution (FREV11 C-1)
- Step 4 asserted smoothness for t > 0 "by parabolic regularity" and continuation "by the literature fact (P1)" (`step4_proof.tex` l.130, 336, 370, 672; [A]); no text identified the validated object (P1 Galerkin solution, P2 normal graph) with a smooth flow.
- Closed by `C1_identification.tex` (fix 14) and, for the a-priori form of Theorem P2, by Lemma (P2-a) with proof in `fix16_texts_v38/C1_identification_v38.tex` (FREV15 C1-1) [hand]. Inherits: Lean read, not built; DKS not re-read; P2's analytic steps [D]; Remark C1:rem:C1 (regularity in tau at slab ends, flagged).

### K6. The v37 strictness sketch for (GS) is invalid (IFR; FREV15 part C)
- `fix12_texts/continuation_v37.tex` Lemma C:lem:GS: the global-sup argument fails (|u| is about 1.1 near x = 60, so sup|k| sup|u| > 1/2); trajectory values and asserted boxes were mixed; [8, 10.25] had no w' bound; X_S = 60 1/4 against boxes ending at 60; openness and closedness were attributed to the wrong lemmas; "at tau_1, implicit in step 4" is not a proof.
- Replaced by `fix14_interface/GS_lemma.tex` (GS-1 at tau_1 [cert] `GS_tau1.json`; GS-2 strictness [cert] `GS_boot.json`, rerun on fin13v4c as `fix15_review/c/GS_boot_v4c.json`; GS-3, GS-4 [hand]), in text v38 `GS_lemma_v38.tex`, `continuation_v38.tex`.

### K7. One bootstrap set (FREV15 D-1, D-2, D-5)
- The bootstrap existed in two versions (`lemmaU_v37.tex` BS_U without BF and GS; `step3_v37_main.tex` with them by reference only); U-1 posed (GS) on |x| <= 60 1/4, which no bootstrap member improves; B1 was stated as sqrt(Lambda) <= sqrt(Lambda_bar) instead of Lambda <= L + 2C'.
- Fixed in `bootstrap_v38.tex` (Definition BS, order table, Proposition BS:prop:close), `U1_regularity_v38.tex` (X_S = 60, O through GS-circ, B1 form), `lemmaU_v38.tex`, `continuation_v38.tex`, `flux2_lemma_v38.tex`.

### K8. Level-2 Garding step and regularity (FREV15 D-4; FREV11 C-6)
- `ei_lemma_v4.tex` item (i), level 2, uses six derivatives of v while U-1 gives v'''' bounded with jumps. The inequality holds by continuity on H^4 with vanishing boundary values and density of smooth functions; the sentence and a list of the integrations by parts with their regularity are in `ei_lemma_v4c.tex` [hand].

### K9. Records (FREV15 A-1r, A-2r, A-3r, B-v37-1, B-v37-3)
- See `fix16_texts_v38/RECORDS_ERRATA_v38.md`: the "same sha256" statement of V4CR / V37R is false for `runs/P_A` and `runs/reg_fin13v4b` (pin 37386e13... of an `e2v4b_constants.json` that no longer exists); the third cert pass has no driver record; `SHA256SUMS_v4c` and `SHA256SUMS_v37` are stale; V37R rounds C0_2, C1_2 down.

### K10. Text slips found by the v38 number check (NC38)
1. `lemmaS_A_v33.tex` .. `lemmaS_A_v37.tex` (ramp box line): "|y f_y| <= ... <= 0.4390", below the asserted value 0.43902 = 0.52/x1.184464. v38: 0.4391.
2. `continuation_v33.tex` .. `continuation_v37.tex` (Lemma G): "the verifier's ramp value <= 7.12e-3", below the asserted 7.1221e-3 = 0.05/x7.020454. v38: 7.13e-3.
3. `step3_v37_main.tex` (K3): "S_1(2 Bhat) <= 0.28578", below the certified 0.285945 (v36) and 0.285924 (v37) (`U_cert_v36_pass2.json` key `supersolution.M1z_hyp_check.zones_I_II.S1_2Bhat`); the conclusion "<= M_1^z = 0.534" is unaffected. v38: 0.28593.
4. The v37 number checker (`fix12_texts/tools/check_numbers_v37.py`) never saw a number written directly after a control word (`\le116.3782`, `\times4.804`, ...) because of its look-behind `(?<![\w.])`, and it treated every line that mentions v35, v36, v37, fin13v4b or fix8 as a record. Its "1077 match" therefore excluded every number written directly after `\le`, `\ge`, `\times` and similar. The v38 checker removes control words first and uses a narrower record filter; the v36 literals it then found in the v37 texts (e.g. C0_2 116.3782, C1_2 0.2606004, abar 2.14041e-7, Ehat_2 0.056798, D_2 0.0552611, kappa_J 0.042432, K0bar 3.4099e-21, C' 2.6524e-26) are replaced in v38 by the v37 values, each with its JSON key (`tools/numsub_v38.json`). These were valid but superseded bounds (v36 data), not wrong bounds, except items 1-3.
5. KAP_OR: AUDIT 1.5883 and V4R 1.58832 are display roundings below the certified 1.5883238 (`constants_v4.json` key kappa/KAP_OR); GI's 1.5833 is a digit slip; the texts use 1.5884 (up). No consequence.

### K11. Other text items (FREV15)
- GS-1: the citation "true-curve shift <= 3.3e-11 (step4_proof)" is found nowhere; the computed bound is <= 1.65e-11 [cert, `fix15_review/c/c_frame_shift.json`], covered x60 by the 1e-9 widening.
- C1-2: Y4's inputs are typed, so the drift formula is [hand] + [sym]. C1-3: Lemma kP2 mixed x and y; conclusion right. A-6: `E0_LEMMA_v4c.tex` quoted fin13v4b R-3 ratios. GS-3: GS-2 uses the B2' input S0V. GS-5: closure.py lines are l.145, l.164 in v4c. All corrected in the v38 copies.
- FREV11 B5 (claims_failing are v3.1 record claims) and B7 ("[hyp]" strings in zone3.py:256, u1_ghat.py:96) are stated in `lemmaU_v38.tex`; B1-B4 and C-7 are chain-side code items and stay open.

## L. Status (v38, November 20, 2026)
**Step 3 is not closed.** Nothing here is called closed.
- **Reviewed, sound modulo listed items, no blocking mathematical item:** fin13v4c and its fin13v4b delta; (A) chains v36 and v37; flux lemma v2 and its v37 reassembly; U-1, G1; C-1 and (GS); the coherence of the bootstrap (FREV11, FREV15, FLUX2_REVIEW).
- **Texts:** v38 set in `fix16_texts_v38/` (`TEXTS_REPORT_v38.md`): one bootstrap set, (P2-a) written, numbers generated from or keyed to the v37 / fin13v4c JSON. Not compiled; not yet reviewed as a set.
- **Open:** step-4 re-pin and kits (the authors' go); the records items K9; the Lean build and an n = 1 instance; a re-reading of DKS for (D1) and the uniform bounds used in C-1; the [read] and not-audited verifier inputs (DS0-DS3, K_rho, K_W, W monotone on [31.25, 34.25], the tail monotonicity statement of VERIFIER.md, the D-1a pins); chain-side B1-B4 and C-7; G1_check rerun on the v37 data; GS_boot_v4c as a pinned stage (with GS-2, GS-6); the regularity of the base in tau at slab ends (C1:rem:C1); LaTeX compile; the final independent review.

## M. Conflicts: what the reports now settle
Settled (the resolution is named):
1. I.4 (TC2): both 1.2709 (AUDIT) and 1.27089 (V4R) are outward roundings of 1.2708883534 (`constants_v4.json`, transport/TC2). Settled.
2. I.5 (KAP_OR): settled for the texts (K10.5); GI's 1.5833 stays a digit slip in a record.
3. I.8 (E0 budget): superseded by K2 (need 38.7232, x1.65694).
4. I.10 and TEXTS_REPORT v37 item 3.1 (E0_LEMMA_v4b w''' note): settled by `E0_LEMMA_v4c.tex` (FREV11 A-3: "Removed since fin13v4b").
5. I.11 (G1 status): reviewed, FREV11 part C (sound modulo items).
6. I.12 and TEXTS_REPORT 3.2, 3.4 (display roundings): v38 prints up-rounded v37 values (C0_2 <= 116.3593, C1_2 <= 0.2605999 at 4 N2_ref; 114.764, 0.0698138 at N2_ref; K2pt <= 1.3900e-19; C' <= 2.6519e-26; D_2 <= 0.0552518). New display items, same kind: V37R C0_2/C1_2 (B-v37-3); FREV15 "5.3902e-4 / 0.040935" (up: 5.3903e-4 / 0.040936); V4CR "x1.65695" (down: x1.65694); E0_LEMMA_v4c d_ab 0.0131136 (raw 0.013113620982; v38 copy 0.0131137).
7. TEXTS_REPORT 3.3 (Delta^R_1): v37 6.79715e-6 as v35/v36; 6.7972e-6 kept.
8. LEMMAU_MAIN_v37_NOTES C3 (D^c_2 only in V36R): the key exists, `A_interface.json` res.tight.r2_export.window_data.D_window_const, v37 <= 0.0552015; used in v38.
9. Notes C4 (kappa_J pad): the v37 key `kappa_J` includes the pad, 0.0424297; text 0.042430.
10. Notes C6 (n1_chart): v38 cites `fix13_v37chain/n1/out/n1_chart_fin13v4c.json` (459fa058...), the file chain v37 consumes; "moves every margin by at most 0.06%" is a v35 measurement, marked [est].
11. Notes T5 (k = 2 range and weight): `G1_export_k2_v38.tex` states them (weight y^-4 (28/y) on [ln 28, ln 350]); consistent with lemmaU. Settled.

Remaining (records or process; no consequence for a quoted bound): I.1, I.3 (pre-fix sizes, [est]); I.2 (process verdict); I.6 (first march cell of fin12f4); I.9 (R-2 cost figures); I.13-I.15 (resolved by FREV in substance); Notes C5 (old "complete/unconditional" sentences in main_v3, covered by the status box), C7 and T1 (`n35_asserts.py`, `check_keys_v34.py`, `constants_needed_v34.json` not regenerated for v36/v37), T2 (two v33-era quotes in Lemma D inputs not recomputed), T4 (main_v3 source for 0.2686, x4.9, x5.3; erratum (O3) column).

## N. Open items (v38)
1. Step-4 re-pin (stages 8, 9, 18, 19, 20 on `step3nl_v4c/runs/P`, `v4bcert/out` incl. `strip_cubic_v4c.json` and `e0_A2_budget_v4c.json`, `fix13_v37chain/runs/prod2`, `fix13_v37chain/flux/out`; add `C1_ke.json`, `GS_tau1.json`, `GS_boot_v4c.json`), a new final_tail driver reading `EXCEPTION_v32_*.json`, a new manifest; kits; Windows runs. Needs the authors' go.
2. Records (K9, `RECORDS_ERRATA_v38.md`).
3. Lean: `lake build` at d72b744 and an n = 1 instantiation file.
4. DKS pp. 1234-35 re-read (D1 and the uniform bounds on d_s^m kappa used in C-1).
5. (P3'), Lemma O arc by arc, the single-point statement (author items).
6. Verifier inputs [read] / not audited (L above); chain-side B1-B4, C-7; G1_check on v37 data.
7. LaTeX compile of every v38 file; final independent review of the v38 text set.
