# REVIEW_v4b_v36_proofs: independent review of the verifier delta fin13v4 → fin13v4b, the (A) chain v36, the proofs U-1 and G1, and the missing audits (fix round 11)

Nov 20, 2026. Fable, independent reviewer. The proof tree was read-only; everything written is in this directory (`fix11_review/`).

- Paths: `V` = `/home/claude/s1spec/step3nl_v4b`; `B` = `/home/claude/s1spec/writeup/lemmaU/bulk`; `S3` = `/home/claude/s1spec/step3`.
- Labels: **[cert]** arb with outward rounding on certified inputs taken as given (mine unless a file is named); **[sym]** exact sympy; **[num]** float or mpmath check; **[hand]** my derivation, not machine-checked; **[read]** verified by reading the cited file:line; **[est]** estimate.
- Parts B and C were done by two forks of this reviewer (same model, same context) and their findings are relayed here; I spot-checked the cited lines of C-1, C-2, C-7. Parts A and D are mine.
- All runs: one process at a time (shared lock), nice 19, PYTHONHASHSEED=0, `python3 -B`, peak RSS 144 MB.
- **Not reviewed:** `B/fix12_texts/` (v37 texts, `ei_lemma_v4.tex`) was being written by another agent while this review ran. Those texts must take the findings below into account, D-3 in particular.

## 0. Verdicts

| Part | Verdict | One line |
|---|---|---|
| **A** verifier delta fin13v4 → fin13v4b | **sound modulo items** | Code equals the E0/E2 specs term by term; every constant I recomputed is covered; the A² budget omits two small pieces and still holds (×1.657 instead of ×1.679 [cert]). No verifier change needed for part A. |
| **B** (A) chain v36 | **sound modulo items** | R3 patches, wiring, κ_J pad and regressions all hold; only record/text items. |
| **C** U-1 with (P2′) | **sound modulo items** | Internal arguments hold; Lean hypotheses taken correctly. Two inputs are closed nowhere: smoothness of γ(t₁) / identification with step 4's solution, and the (GS) graph hypothesis in the bootstrap. |
| **C** G1 (Lemma D, k = 2) | **sound modulo minor items** | Proof and code agree; C0_2, C1_2 are what the proof yields. |
| **D** missing audits | **needs fixes (one)** | θ-certificate, Gårding identity, Q2 Gram, `cubic_E0fine`: sound. **`ADD_STRIP_CUBIC = 2.6e-3` is not a bound** of the strip part of the E0 cubic form; it needs a certified replacement in the verifier, hence one more verifier pass and chain rerun. Size is small: ≤ 1.3e-4 of rate at the peak [est, over-estimate]. |

**Step 3 is not closed.** The one new blocking item is D-3. I recommend holding the step-4 re-pin and the kits until it is decided, because fixing it changes the verifier tag and (slightly) the numbers v36 consumes.

## 1. Part A: fin13v4 → fin13v4b

### 1.1 What I did
- Diffed `verifier/{constants,rhs_iv,closure,unc,mod_sups}.py` against `step3nl_v4/verifier/`. The `_pre_v4b.py` backups are byte-identical to fin13v4. Every hunk is gated by `V4B` and is one of: E0 F-B, E2 fix, R-2, R-3, args records. Nothing else.
- Read `v4bcert/{tc_weighted, e0_strip_drift, e0_lemma_consts, e0_identities, r2_strip_ibp, e2_v4b_constants}.py` and `E0_LEMMA_v4b.tex` line by line against `rhs_iv.py:192-219, 274-294`, `closure.py:70-95, 180-190`.
- Recomputed in arb with my own code: `a_recompute.py` → `a_recompute.json`, 23 checks, 0 fail.
- `sha256sum -c V/SHA256SUMS` passes (651 files); `step3nl_v4/SHA256SUMS` passes; `audit_v4b/v4_before.sha` = `v4_after.sha`.

### 1.2 E0 fix
**Projection [cert].** From `gram7_g3.json` with my own formula: ‖g0‖ ≤ 7.6969814466, n3_123 ≤ 6.5881854016. Coded N3P = 6.5881854035 covers it. The literal 7.52 (`rhs_iv.py:50`, q_β term) is valid only as a bound of n3_123.

**Rank-one symmetrisation [hand + cert].** The identity sup|⟨w,u⟩⟨w,v⟩| = ½(‖u‖‖v‖ + |⟨u,v⟩|) is right. It applies because β_3 = H33⟨w̃,u_3⟩ with u_3 = (L* − λ̃)f3, and u_3 is supported on the cutoff strip and the tip bumps: `S3/lyap_consts2.py:226` has r3 = Q2 + |a|Q4 + |b|Q5 with no core term, i.e. F is treated as an exact solution of the adjoint equation [read]. I confirm r3 = r3_far + r3_bump and r3 ≤ ε1 [cert]. My values: X ≤ 1.02383e-5, Π3_true ≤ 0.01038536 ≤ 0.01607 [cert].

**A² absorption [cert given reads]: holds, but the written budget omits two pieces.**
- I reproduce the reported budget: need 38.5411 at α = 2.5e-4, excess 3.3853e-3, slack 5.6846e-3, margin ×1.6792.
- **A-1.** `Qfull` is built with the 2×2 lead (a,b) = M⁻¹(⟨w,L*f1⟩,⟨w,L*f2⟩) and H33 = 1/G33 (`agent_mod/bloop/bloop_forms.py:5, 50-52`; BLOOP.md caveats). The modulation lemma's lead is the 3×3 H form. The difference, bounded by d_ab·A with d_ab = 0.0131136 [read, `e2v4b_constants.json` DAB] plus dlead 3.8e-9, couples in E0 through α(g_A,g_B) and through the loop. It is in neither the budget (`e0_lemma_consts.py:144-151`) nor H4 (`rhs_iv.py:205`). Size: (d_ab + dlead)(‖(g_A,g_B)‖ + L_C) = **0.17938 per αA²** [cert]. The E2 block has the same piece and codes it (DAB in `rhs_iv.py:213`); E0 does not.
- **A-1b.** The third-order piece α²|r|·M_ab·α in the loop is also omitted: 0.00264 per αA² at α = 2.5e-4 [cert].
- With both added: need 38.7232, excess 3.4308e-3, **margin ×1.6569** [cert]. The inequality of Lemma 1 is therefore true as stated; its proof (iv) is incomplete.
- Fix: two lines in the lemma's item (iv) and in the budget script. The constants the verifier consumes do not change, so no rerun for this.

**Forcing budget [cert].** 5.343946 ≤ 5.72 (×1.0704) at α ≤ 2.5e-4, all seven pieces recomputed. FE2 = 13.98386 [cert].

**Higher-order terms H1–H7 [hand, read].** Each coded line matches the lemma and I re-derived each bound:
- H1 (strip α-linear form): weighted divergence form with X = χ²; pair bounds via IR, IRW, RW are right; the drift coefficient is added to m00. Taylor form b(α)/α = b′(0) + α·b″(ξ)/2 is coded correctly (`e0_strip_drift.py:90-95`).
- H2 (cubic drift), H3 ((a,b) w-parts), H5 (α² remainders): algebra as in I4–I6; ‖χw′‖ ≤ W1 + ‖χ′w‖ and ‖χ′w‖ ≤ (35/48)·IR^½·sup|w| are right.
- H4, H6, H7: the split of (a,b) and β into lead and rest matches the modulation lemma's bound, except for A-1.
- Drift identity c_dr: re-derived by hand, equals `tc_weighted.py:105`.
- LOOPC, GABN, HOB recomputed, coded values cover mine [cert].

**Assertions.** `E0_alpha`, `E0_w_all34`, `E2_alpha`, `E2_tip_jet` are in `closure.box_assert` and are evaluated at every march step and at the tail anchor. All box ratios ≤ 1 in both rate tables [cert read].

### 1.3 E2 fix
- `rhs_iv.py:111-114, 211-219` is E2_AUDIT §5 F2/F3 literally: RABx, n3e, wprt, u0pr, RM, NS, FE2. R143 is read from JSON. `closure.py:77-79, 95` is F4.
- nd/nf: every coded value ≥ the region maxima of `nl_regions_a2.5e-4.json`, with the column map of `rhs_iv.py:282-284`; the L² coefficients nd+nf cover the maximum over all four regions [cert].
- The far box 0.02 on [16.25, 34.25] is asserted (`closure.py:89`), the IR box 1e-3 on [0, 16.25] by `Bb_IR_jet`, the tip box 1e-4 by the new assert.
- D_ρ constants: sup 2|φ′| ≤ 2.16483, sup|φ′² − φ″| ≤ 1.16109 on [0, 34.25] [cert], so 2.2 and 1.2 are valid.
- The E2 derivation (E2_AUDIT §1) is Fable's own and remains [hand]; I checked T4 strip, T5b, T5c, T5d by hand and found no error.

### 1.4 R-3 (whole-tail α-condition): sound
- The monotone split in `closure.py:182-186` is right [hand]: with α = κ√ℓ, the terms BETA1·κ√ℓ, BETAF·Fn(κ√ℓ)/√ℓ, QBL·√ℓ are nondecreasing in ℓ; STRIP_BETA/√ℓ is decreasing. On the flat part ℓ = L_F the condition is weaker.
- Recomputed [cert]: lhs 177.1921, rhs_all 0.1139682 (1.5e-4); lhs 175.2133, rhs_all 0.1158487 (1.25e-4). Equal to the stored values. Direct evaluation of β(κ√ℓ, ℓ)/√ℓ on 4001 points gives max 0.0821 / 0.0845, below the monotone bound.
- **A-2 (wording).** "Margin 177" in V4B_REPORT and CLAUDE.md is the difference lhs − rhs. The ratio is **×1554.8 / ×1512.4** [cert].

### 1.5 R-2 (w‴ strip hypothesis removed): sound
- The integration by parts is right [hand]: ∫ρf·4χ′w‴ = −∫(4ρfχ′)′w″, new w″ coefficient ρ(2fχ″ − 4f′χ′ + 4φ′fχ′). χ′ vanishes at both strip ends (Sm′ = 140t³(1−t)³).
- f_j′ conventions agree with the closed lemma code: F′ = ws1 + φ′ws0 as in `S3/lyap_consts2.py:140`; χ derivatives carry (−1/3)^j.
- The old κ are reproduced exactly; the S3 coefficient is zero (stored as 5e-324); the verifier reads `STRIP_BETA_R2`, `STRIP_AB_R2` (`constants.py:143-144`).
- The state form BTs/ABs in `rhs_iv.py:149-150` still uses κ_{j,3}·sup|w‴| from E_I norms. That is not a hypothesis and stays valid.

### 1.6 Findings, part A
| # | file:line | severity | finding | fix |
|---|---|---|---|---|
| A-1 | `V/v4bcert/e0_lemma_consts.py:144-151`; `E0_LEMMA_v4b.tex` item (iv); `rhs_iv.py:205` | minor (proof incomplete, inequality true) | Lead mismatch of (a,b) (0.1794 αA²) and third-order loop piece (0.0026 αA²) missing from the A² budget. Margin ×1.679 → ×1.657 [cert]. | Add both to (iv) and to the budget script; cite `a_recompute.json` until then. |
| A-2 | `V4B_REPORT.md` (R-3 rows), CLAUDE.md | text | "margin 177" is a difference; the ratio is ×1554.8 / ×1512.4. | Reword. |
| A-3 | `E0_LEMMA_v4b.tex:101` | text | "Not covered here" still lists the sup\|w‴\| ≤ 0.02 strip hypothesis, which R-2 removed. | Delete; add a line on R-2 and R-3. |
| A-4 | `E0_LEMMA_v4b.tex` item (ii) | text | The support claim for u_3 needs its reason: f3 = χF − bumps with (L* − λ̃)F = 0, so u_3 = [L*,χ]F + bump terms (`lyap_consts2.py:226`). | One sentence. |
| A-5 | `V/v4bcert/e0_strip_drift.py:41, 175` | minor ([read] input) | RW, IRW use W ≥ 1/KW on the strip, KW = 1/W(34.25) for the blended weight; this needs W decreasing on [31.25, 34.25] (`agent_sob/sob_constants_P2.py:251-258`). Not re-derived by me. | List among the [read] inputs. |
| A-6 | `rhs_iv.py:40`, `unc.py` | info | 5.72 and 7.52 remain typed literals; valid by the lemma. | None. |

## 2. Part B: (A) chain v36 (fork; scripts `b_check.py` → `b_check.json`, 23 checks, 0 fail)

**Sound modulo items.** No blocking or major finding.

- **Diff audit [read].** Every hunk of the 15 backed-up scripts is R3 (gated `V36_R3`), r_2 (gated `V36_R2`), the κ_J pad, wiring, or args records. `patches/*.diff` are byte-equal to fix4_r3's.
- **R3 [hand + cert].** r0′ = c_U(β̄ + DL·α₁) + c_U·p·α₁(ā + b̄/14), DL = 3.2224577463e-6; r0′ = 3.75163420346e-9 and c_dU = 1.46342753445e-8 recomputed, equal to stored; 2e-8/c_dU = 1.36665; asserted in `u34_sweeps.py:112`, `zone3.py:135`. c_dU uses the v4b sups 2.14041e-7 and 2.35488e-9.
- **Wiring [cert].** Consumed J's, slopes and peak ≥ raw `V/runs/P` values (excess 1.5e-6 to 8.4e-5 relative); ρ = 2.677e-3 ≤ raw. Only the 1.5e-4 row is consumed and it dominates the 1.25e-4 row in every consumed quantity. All 259 sha256 pins in 131 args blocks match the disk.
- **κ_J pad [cert].** Certified κ_J = 0.04243186976; the pad 1e-9 is added after the certified sum (`hjw_rep.py:39`), conservative. C′/X_const = 7.1e-12, covered ×141.3.
- **r_2.** Touches level 2 only; retained 0.001001 ≥ r_2 asserted; the only margin moved is Ĝ22, 1.91741 → 1.91640.
- **Regressions.** `cmp_runs.py` rerun: v36 code with v35 inputs = v35, and R3-only = fix4_r3 (path strings only differ). prod1 = prod2 byte for byte in every JSON. `audit/tree_before.sha` (5391 files) and `step4_before.sha` (1812) pass now.
- **Headline numbers [read].** Ramp M1 1.184464; late/present M1 2.353047 / 2.770166; H6 1.860867; zone III 1.238888; C0_2 114.78222, C1_2 0.0698142389 (design) and 116.37817 / 0.26060039 (4·N2_ref).

| # | file:line | severity | finding | fix |
|---|---|---|---|---|
| B1 | `fix6_v36/t2_A/base_ibp/` outputs (lemmaS_F10p, cap_F10p, rigcap, weights_exact, assembly_arb .json); `assembly_arb.py:164` | minor (records, inherited) | No args/sha256 block although they now depend on the verifier. | Add at the re-pin, or pin in the new manifest. |
| B2 | `v36cfg.py:67-69` | minor | Env check silently skipped if the rate table lacks `args`. | Require the keys for fin13v4b. |
| B3 | `t5_zone3/out/ramp_U.json` `u1_ghat_sha` | minor (known) | Hash matches neither v35 nor v36 u1_ghat.json; content dependency verified equal. | Rerun ramp_U at re-pin or record an exception. |
| B4 | `t2_A/rigor/assembly_rep.py` | minor | "Ramp sups ≤ box" and "no margin < 1" are read off afterwards, not asserted in the chain. Minimum margin 1.18446 [read]. | Add the asserts. |
| B5 | `asm_*.json` `claims_failing` | text | Non-empty in every run (v3.1 record claims); not consumed by Lemma U or zone III. | Say so in the (A) text. |
| B6 | `V36_REPORT.md` §4, §6 | text | FLUX ratios superseded by `flux2_cert.json`; C0_2 "114.782", C1_2 "0.0698142" are rounded down. | Mark and correct. |
| B7 | `zone3.py:256`, `u1_ghat.py:96` | text | Strings still say c_δU ≤ 2e-8 "[hyp]". | Wording. |

Not checked in B: the inherited v35 mathematics beyond the diff; the numerical content of n1_chart on the v4b sups; the chain was not rerun.

## 3. Part C: the two proofs (fork; reruns in `c_rerun/`: G1_check 35/35, U1_sym 32/32, both equal to the originals)

### 3.1 U-1 with (P2′): sound modulo items
- **Lean hypotheses are taken correctly [read].** The source is on disk (`/mnt/c/Users/micahw/agents/lean/HSLflow-Whitney`, HEAD d72b7447; the five cited files hash to the prefixes quoted in `U1_regularity.tex:32`). In `Flow/Main.lean:29-36`, `Flow/UniqueMain.lean:29-34`, `Flow/GraphUnique.lean:42-46`, `Weinstein/Global.lean:53-57`, `Scalar/Pointwise.lean:60-90`:
  - the only dimension hypotheses are `NeZero (finrank E)` and `finrank V = 2·finrank E`; no dim ≥ 2 hypothesis;
  - the Lagrangian clause is vacuous for curves;
  - `IsFlow` is jointly C^∞ on [0,T]×M with purely normal velocity, so uniqueness is literal equality of the maps, with no reparametrisation quotient. This is the class DKS uses.
  - The repository was not built here, and the n = 1 instance is not elaborated.
- **Gluing (`U1:63-67`): sound [hand].** Every ∂_t^m∂_x^k Φ is a smooth function of spatial jets, so the one-sided values agree at the junction.
- **Reflection (`U1:74-78`): sound [hand + sym].** σJ = −Jσ and the velocity has two J's; S3 covers both orientations. It rests on uniqueness (L2) and on step 4's exact symmetry of γ(t₁).
- **Blow-up (`U1:80-85`):** logic sound. DKS2002.pdf could not be opened here; (D1) rests on the quote in P1_CONFIRMATION §1.
- **(P2′) (`U1:107-123`): sound [hand].** IFT, conserved constraints and the compactness step, including T_max = ∞.
- **Theorem U-1 (a)–(d):** differentiation under the integral with the C³ cutoff, the integration by parts with Lipschitz data, and the Agmon step are correct.

### 3.2 G1: sound modulo minor items
- The Young step, the integration and `assembly_rep.py:235-236` agree; `trq`, `trc`, `trp` (`:108-118`) are valid window forms; the cubic bound is degree 1 in J and Lemma S holds on any interval; no level-2 source escapes the window form.
- C0_2 ≤ 116.3782 and C1_2 ≤ 0.2606004 are what the proof yields and what flux lemma v2 used.

### 3.3 Findings, part C
| # | where | severity | finding | fix |
|---|---|---|---|---|
| **C-1** | `U1:23, 191`; `step4_proof.tex:672` | **major (interface)** | Smoothness of γ(t₁) is step 4's unproved "[A] parabolic regularity" sentence. No text identifies step 4's validated C⁰H⁵ solution (its own parametrisation) with a smooth normal-velocity flow. | Apply (L1)/(L2) from Γ0 (a trigonometric polynomial) and prove T_max > t₁ and equality with the validated solution up to reparametrisation; or prove smoothing. |
| **C-2** | `U1:92-95, 108, 130, 194` | **major (lemma + script owed)** | (GS), the tip arc as a normal graph over S1, is assumed and is in no bootstrap set; nothing shows it improves strictly. It is posed on \|x\| ≤ 60.25 while the verifier boxes stop at x = 60. Room is large [est]. | Add (GS) to the continuation bootstrap with an arb check from the boxes and (G)/chart on [60, 60.25]; or pose it on \|x\| ≤ 60. |
| C-3 | `U1:32, 45, 192` | minor | n = 1 instance (E = ℝ¹, M = S¹, V = ℂ) not elaborated; repository not built. | A short Lean file instantiating both theorems, plus `lake build` (Micah). |
| C-4 | `U1:84` | minor | "Φ(T_max) is an immersion" is attributed to DKS. | Add: ∂_t ln\|Φ_x\| = k·k_ss is bounded, so \|Φ_x\| ≥ c > 0. |
| C-5 | `U1:99-101` vs `108` | minor/text | Openness is stated for strict (GS) but 𝔅 is the closed version. | Define the open set with relaxed bounds. |
| C-6 | `U1:136, 162` | minor | The claim that all E_I integrations by parts need only "v‴ Lipschitz, v⁗ bounded" cannot be checked until `ei_lemma_v4.tex` exists. | That text must list each spatial IBP with the regularity used. |
| C-7 | `fix6_v36/t2_A/rigor/weights_t2.py:232` | minor rigor | Zone-cell `bc = floor(float(avail.lower())·1e6)/1e6` passes through a nearest-rounded float with no arb assert per cell; `G1_check` B3 only reads flags. Deficit ≤ about 1 ulp against r₂ = 1e-3 [hand]. Same pattern at k = 0, 1. | Add the arb assert at the next chain touch. |
| C-8 | `PROOFS_REPORT.md` | text | Says 36 G1 checks; there are 35. | Edit. |

Inputs left unproved by the two proofs: (H6) and the symmetry of γ(t₁) (step 4); the Lean build and n = 1 instance; the DKS extension step; (GS) at τ₁ and its persistence; the graph equation, the identification of the modulation matrix with G + E and the Neumann bound 0.0466; the C⁴ matching in `c4fun.py`; for G1, Steps 1–7 of Lemma S(A) at j = 2 and H-J-W on windows.

## 4. Part D: the missing audits

### 4.1 θ-certificate (−1/32, 1/9): sound
- **Statement and scaling [read + hand].** Re⟨Lv,v⟩ + θ‖v″‖² + ω‖v‖² = (1−θ)[−‖v″‖² + m∫ρ(A v′² + (V+ω)v²)], m = 9/8. The generators multiply both A and V − λ by FORM_MULT (`S3/gen_hs3.py:58`, `gen_hs2.py:45`; `positivity.py:24`; `hstail.py:42-47`), so the shift is m·ω as it must be. The theorem really is ≤ −(1/32)‖v‖², not −(1/36).
- **Constraint set.** f1 = χp, f2 = χψ2, f3 = χF − bumps (`gram_post8.py:2-6`), the same f_j as the modulation lemma.
- **Post-processing recomputed [cert].** `gram_post8.certify` on the 34 `outT1*.txt` reproduces the stored G3 and S ∈ [−4.9266e-4, −1.7298e-4]. With plain intervals from the stored mid/rad I get the same enclosure: g11 < 0, det > 0, S < 0, so G3 is negative definite and n(G) = 3 = the Morse count. All 34 CAPD outputs and the Morse run say RESULT OK; `positivity_cert_T1.json` carries the same (9/8, 1/32).
- **Not done:** the CAPD runs themselves were not rerun; the pipeline is the closed linear certificate's with two environment switches. The certificate JSON has no args/sha block (pins of the 37 files are in `d_checks.json`).

### 4.2 `cubic_E0fine`: sound modulo one unasserted box
- **Derivation [hand].** R3 = ∫₀¹(1−t)·D²(T_c/ρ)(u0+tw)[w,w;w]dt; the ½, the symmetric double-sum bookkeeping (`cubic_E0fine.py:25-34`), the exact expansion in U2 (degree ≤ 2) and the region/L²/sup assembly in `rhs_iv.py:29-37, 262-272` are right. Coverage x ≤ 31.25 is contiguous.
- **D-2a.** The coefficients are certified for |w|, |w′| ≤ 0.01 and α ≤ 1e-3. On [0, 16.25] the asserted box is 1e-3. On [16.25, 31.25] the asserted box is 0.02 (`closure.py:88-89`), which does not imply 0.01. In the runs the worst-case sups are |w| ≤ 5.91e-4 and |w′| ≤ 1.161e-3 [cert values of the rate tables], ×8.6 inside. So the hypothesis holds; it is not asserted.

### 4.3 `ADD_STRIP_CUBIC = 2.6e-3`: **not a bound (D-3)**
- `constants.py:38`, `rhs_iv.py:273`: the strip part of the E0 cubic form is coded as 2.6e-3 × (the [16.25, 31.25] region term). That term is a product of ρ-norms of w̃ = χw. No derivation exists (VERIFIER.md:7 only states it).
- It cannot be a bound [hand]. With weight X = χ² the strip form has terms proportional to X₁ = 2χχ′ and X₂ = 2χ′² + 2χχ″, in which w is not multiplied by χ. Near σ = 32, χ′² ≫ χ², so a state can have a strip cubic term far above anything expressed in ‖χw^(j)‖. And for the X-part the factor would have to be about 1, not 2.6e-3, with a strip coefficient ≥ the far one.
- **Certified replacement data [cert, `d_cubic_strip.py`].** Sups over the strip of ½·Σ|∂²C_l/∂U_i∂U_j| for the weighted form of `tc_weighted.pkl`, box |α| ≤ 2.5e-4, |w|,|w′|,|w″| ≤ 0.02 (all asserted). Monomial (a,b,c) = |w|^a|w′|^b|w″|^c:

  | (0,1,2) | (0,2,1) | (0,3,0) | (1,0,2) | (1,1,1) | (1,2,0) | (2,0,1) | (2,1,0) | (3,0,0) |
  |---|---|---|---|---|---|---|---|---|
  | 0.56703 | 2.57194 | 0.76568 | 0.28319 | 3.54183 | 1.06797 | 2.5e-17 | 5.6e-18 | 4.1e-37 |

- **Proposed term H8** (state form, same pattern as H1; sup on the lowest-order factor):
  RW·[(c030·sI1 + c120·S0s)·V0² + (c021·sI1 + c111·S0s)·V0V1 + (c012·sI1 + c102·S0s)·V1²] + IRW·S0s²·(c201·V1 + c210·V0) + c300·IR·S0s³, replacing `astrip·cubfar`. IR, IRW, RW as in `e0v4b_constants.json`. All terms have x-degree 3, so the tail lemma's degree conditions are unaffected.
- **Size [est].** At the production peak, with every norm at its own maximum (mutually exclusive, so an over-estimate): 3.67e-19, i.e. **1.3e-4 per unit Λ**. Compare first-cell rate +0.0304 and tail ρ 2.678e-3. In the shape where V1 is maximal the E_I dissipation is about 1.6e4 per unit Λ, so the realistic effect on ρ is far below 1%. Closure is not at risk; the coded inequality is simply not an upper bound until this is fixed.
- No existing slack absorbs it as a statement one can write: the PI3 slack is A²-type, H1 has no reserve, and the coded `astrip` term is in the wrong variables.

### 4.4 Q2 Gram (43.33): sound
- The form is Q2(w) = −a⟨D_ρf1,w⟩ − b⟨D_ρf2,w⟩ with (a,b) = M⁻¹ℓ; `mod_forms.py:51-70` builds it correctly, and the Loewner-majorant argument (G ≤ mid + diag(row radii)) is valid [hand].
- My own code on `gram8_g1s4.json`: float λ_max 43.32424; certified λ_max ≤ 43.3271 by the majorant route with an arb Cholesky [cert]. Coded 43.33 covers it. A direct congruence test (cG − GKG ⪰ 0) is inconclusive because the projected Gram is nearly singular; recorded as information only.
- Not re-derived: the Gram integrals themselves (`mod_gram.py`).

### 4.5 Gårding identity (r143, r163, r107): sound
- Re-derived by hand, including every boundary term at x = 0 (even v; c1 = −(5km + q/4) and c1t = −(9km + q/4) vanish at the tip): ⟨Mv,v⟩_ρ = ∫ρ(−v″² + Av′² + Vv²) and ⟨Mv,v⟩ ≤ −½‖v″‖² − ∫ρμv².
- [num] mpmath on generic even data: identity to 0.0 at 30 digits; inequality holds.
- [sym] (Lw)″ = Mt w″ + e1 w′ + e0 w with the stated Mt; residual 0.
- [cert] my own expressions on the same hulls with finer grouping: sup(−μ2) ≤ 14.34288, sup|e1| ≤ 16.24913, sup|e0| ≤ 10.66858. Coded 14.343003, 16.3, 10.7 cover them.

### 4.6 Float check of the β decomposition: passes on a model
- `d_beta_float.py` [num], 40 seeds: a finite-dimensional model with the same algebraic structure (arbitrary L, D, cutoff, nonlinear F with DF(0) = L, A1, B1 with α_g). Checked to rounding (≤ 1.2e-13 relative): x = HR − HEᵀx with the lemma's E and R_j; β = β_3 + β_12 + β_Lin + β_Q + β_F + β_S + β_E; the E0 pairing identity with 𝓜 as in the lemma (no ψ1, ψ2 terms); ⟨g0,w̃⟩ = ⟨Πg0,w̃⟩; the oblique identity with ⟨g0,h⟩ = 1; E33 = 2α⟨χh2,f3⟩ exactly; E13, E23 first order with the signs of BLOOP.md; the rank-one bound (worst ratio 0.246); the self-adjoint step of the Q2 lead.
- This confirms the hand identification of the pieces and their signs. It does **not** test the actual profile functions or any certified constant; a check with the true ψ̂3, f_j was not done.

### 4.7 Findings, part D
| # | file:line | severity | finding | fix |
|---|---|---|---|---|
| **D-3** | `V/verifier/constants.py:38`, `rhs_iv.py:273` | **major (blocking for "closed"; numerically small)** | `ADD_STRIP_CUBIC` is an underived scaling, not a bound of the strip cubic form. ≤ 1.3e-4 of rate [est]. | New verifier tag with term H8 (§4.3), coefficients recertified in the chain; then rerun the (A) chain and re-read the flux assembly. |
| D-2a | `V/verifier/closure.py:88-89` vs `cubic_E0fine_0.001_0.01_0.01.json` | minor | Box 0.01 on [16.25, 31.25] not asserted (0.02 is). Holds ×8.6 in both runs. | Assert `far_w`, `far_w1` ≤ 0.01 at the same verifier touch; until then `d_cubic_strip.json` is the record. |
| D-1a | `S3/certificate_step3_theta_T1.json` | minor (records) | No args/sha256 block; CAPD runs not rerun by anyone since Nov 6. | Pin the 37 files in the next manifest. |
| D-6a | — | minor | β float check done on a model only. | Optional: repeat with the true functions when the E_I lemma text is final. |

Still not audited by anyone: the level-2 commutator constants DS0..DS3, the K_ρ/K_W Sobolev constants, and the [read] inputs of the E0/E2 lemmas (Gram integrals, `pair_f3`, `bloop.json`, `final_5e-4.json`, `lyap_consts2`).

## 5. Updated closing list for "step 3 closed"

Numbering follows the 14 items of `fix9_review/FLUX2_REVIEW.md` §7. New items carry N.

**A. Flux lemma v2**
1. ~~U-1 written~~ **done and reviewed (sound modulo items).** Remaining: (1a) the (GS) bootstrap lemma with an arb script (C-2); (1b) text fixes C-4, C-5.
2. ~~G1 written~~ **done and reviewed.** Remaining: text edits of `lemmaS_A` (r_2, level-2 export) and `lemmaU` Lemma D (k = 2); per-cell arb assert (C-7); PROOFS_REPORT count (C-8).
3. Flux text items T1–T4: unchanged (not in this review's scope; `fix12_texts` unreviewed).
4. Pin `fix8_flux2_impl/out/*` in a manifest and kit: unchanged.

**B. Verifier**
5. `ei_lemma_v4.tex`: still owed (a draft is being written in `fix12_texts`, unreviewed). It must (i) list each spatial integration by parts with the regularity used (C-6); (ii) not present 2.6e-3 as derived (D-3); (iii) carry the corrected A² budget (A-1).
6. ~~Independent review of the v4b delta~~ **done: sound modulo A-1..A-5.** Remaining: lemma text and budget script fixes (A-1, A-3, A-4), wording (A-2).
7. Audits: θ-certificate, Gårding, Q2 Gram, `cubic_E0fine`, β float check **done** (§4). Remaining: D-2a assert; D-1a pins; DS constants and Sobolev constants unaudited.
- **N1 (new, blocking). Verifier pass with a certified strip cubic term (D-3)**: new tag, term H8, assert D-2a, optional asserts B4. Compute: about 2 × 12 min for the verifier, about 7 min per (A) chain pass, then `flux2_assemble` and `mod_sups`. Needs Micah's go, since it supersedes fin13v4b and v36 as the objects to pin.

**C. (A) chain, Lemma U, zone III**
8. ~~Independent review of v36~~ **done: sound modulo B1–B7.** If N1 is done, the chain is rerun on the new tag (same code; regressions as in v36).
9. Texts: as before, plus B5, B6, B7.

**D. Step 4 and kits (Micah's go)**
10. Re-pin stages 8, 9, 18, 19, 20, new final_tail driver reading the EXCEPTION record, new manifest. **Hold until N1 is decided.** Add B1 (args blocks or pins for the five base_ibp outputs), B2, B3.
11. Kits rebuilt and run; Windows runs.

**E. Author items**
12. (P1′) local existence with uniqueness = the Lean CWFlow theorems; hypotheses confirmed by reading the source; **owed: `lake build` and an n = 1 instantiation file (C-3)**. (P2′): proved (Proposition in U-1). (P3′), Lemma O arc by arc, single-point statement: unchanged. DKS extension step: authors re-read pp. 1234–35 (C-4).
- **N2 (new). Step-4 to step-3 interface (C-1):** smoothness of γ(t₁) and identification of step 4's validated solution with the smooth flow; plus (GS) at τ₁ in the hand-off list.

**F. Whole document**
13. LaTeX compile: unchanged.
14. Final independent review of the assembled step 3: unchanged; it should include `fix12_texts` and the N1 delta.

## 6. Files in this directory
- `REVIEW_v4b_v36_proofs.md`: this file.
- `a_recompute.py` → `a_recompute.json`: part A, arb, 23 checks.
- `b_check.py` → `b_check.json`: part B, arb, 23 checks.
- `c_rerun/`: copies and reruns of `G1_check.py`, `U1_sym.py` with outputs.
- `d_checks.py` → `d_checks.json`: θ-certificate post-processing, Gårding, D_ρ constants, Q2 Gram (11 checks).
- `d_cubic_strip.py` → `d_cubic_strip.json`: cubic box check and certified strip cubic coefficients.
- `d_beta_float.py` → `d_beta_float.json`: float check of the β decomposition (model).
- `SHA256SUMS`.

Proof-tree check after my runs: no file under `step3`, `step3nl`, `step3nl_v33`, `step3nl_v4`, `step3nl_v4b`, `s1_tight`, `step4` is newer than the start of this review; no `__pycache__` was written (`S3/__pycache__` holds two pre-existing files of Nov 13).
