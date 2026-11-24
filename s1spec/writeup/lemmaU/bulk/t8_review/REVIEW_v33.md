**Verdict: v3.3 does not yet close step 3.** No inequality was found to fail and every certified margin I re-derived holds, but four items block "closed" and "unconditional": the tightened (H2) is certified for a different deviation than the one v3.3 uses, a main lemma has no written proof, the passage from the bootstrap to the theorem is not written for v3.3, and main_v3 quotes verifier numbers that v3.3 cannot use.

Task 8, independent review (Fable), Nov 13 2026. Paths: `B` = `/home/claude/s1spec/writeup/lemmaU/bulk`, `V` = `/home/claude/s1spec/step3nl_v33`, `S4` = `/home/claude/s1spec/step4/p2`. Labels: **[checked]** = I read the code or re-ran it; **[agent]** = found by a sub-reviewer, spot-checked by me where marked.

## Findings

### Blocking

**B0. The shell datum (H2) and the variable of (A) do not match the deviation v3.3 uses.** [agent, definitions confirmed by me]
- v3.3 takes z^A_j = ∂_s^j(f − U), U = 1.3056·α·y^p, on all of [14, 700]: `lemmaU.tex:46`, the datum E^A(τ₁) ≤ H2 + H2f (l.403, (K6)), launch lemma (a) c₂·H₂ on [28, 60] (l.116–120), and l.60 says "(H2) … already stated for f − U₁".
- Step 4 certifies (H2) for the v3.2 **shell variable** w′ = f − α·ψ̂₃′(x) − α²·ĥ₂′(x): `S4/t9h/win.py:3, 64–72`; `S4/v33/h2t_v33.py` only re-assembles those shells. The verifier's J_m and (A)'s ramp/core data are also for w′ (`writeup/lemmaS/lemmaS.tex:69, 323`; `step3_proof_v32.tex:71`).
- So E^A as used in v3.3 joins two different deviations at y = 60, and nothing bounds the difference d = α(ψ̂₃′ − 1.3056·y^p) + α²ĥ₂′ on [14, 60], at τ₁ or later.
- Failure scenario: the certified 8.0e-11 bounds ∫ω w′²; ∫ω (f − U₁)² on [14, 60] could be near 1e-7 [est, from the certified 4.84e-8 on [60, 379.5]], against the tightened 1e-8. With 1e-7 the datum is about 1.17e-6 against Ê₀ = 1.188e-6, and Ê₀, the flux margin and the late-centre launch all move. The same d enters Lemma U's ramp conversions (Step 6) and launch (b) at every τ, next to δU.
- Fix: either certify (H2) for f − U₁ on [14, 60] (win-mode cells; the error part carries over), or keep w′ in (A) and add a certified bound on ∂_s^j d (j ≤ 3) on [14, 60] to Lemma dU, the launch lemma and (K6). Then rerun t2 → t4.

**B1. The repaired Lemma S(A) has no proof text.** [checked]
- Lemma U takes its level-k inequality as hypothesis (D1) (`B/t3_proof/lemmaU.tex:153–157`) and (K1) rests on it. main_v3 cites "sources `t2_A_v33`, erratum" for it (`B/t7_assembly/main_v3.tex:302`).
- What exists: agent reports (`B/t2_A_v33/A_REPAIR.md`, `A2_REPORT.md`), scripts, and `B/t1_statements/STATEMENTS_r2.tex`, where R5′ and the closure have a "Proof sketch" (l.233, l.296).
- Not written anywhere as a proof: the ramp repairs R1–R3, the point-term routes (Lemma P case (i) at level 2), the y = 28 interface repair A5 (C^{1,1} Ỹ, partition-of-unity localisation of the L⁴/S33 identities, `rigor/weights_t2.py:14–20`), the tools T1–T7 class by class on O, and why (D2) holds on every window (all data sources of degree ≤ 1 in the J's).
- Missing step: a `lemmaS_v33.tex` that states the inequality (D1) with its constants and proves it term by term, as lemmaU.tex does for Prop. U. Until then step 3 is "certified by scripts", not proved.

**B2. No v3.3 argument leads from the bootstrap to the theorem.** [agent, confirmed by me]
- `lemmaU.tex:397` proves BS_U on [τ₁, τ_max) only. Global existence and (W4)/(W5) are left to v3.2.
- The v3.2 continuation proof (`writeup/step3_proof_v32/step3_proof_v32.tex:358`) uses "the shell energies up to level 2 (with the strip to level 3)" and "(B) and Lemma O beyond". The single-point argument (l.382) uses "the windows (B) (curvature ≤ S₁^f/r for r ≤ 150λ₁)".
- v3.3 drops B3 and (B), and Lemma U has levels 0–2 only: between y = 60 and r = 150λ₁ the curvature is bounded pointwise (|y f_y| ≤ 0.52) but its derivative only in L².
- Missing step: restate the continuation criterion (P1) and §5 of v3.2 with B2′, B3U, Z1–Z3, U-fine and B4, and say which norm (P1) needs. I expect it to go through; it is not written.
- There is also no assembled v3.3 step-3 document (bootstrap set, order table and D1–D6 in one place).

**B3. main_v3 states fin11 tail numbers; v3.3 runs on fin12f4.** [agent, confirmed by me]
- `main_v3.tex:111` (Main Theorem (iv)) and l.299: floor √Λ_F = 4.94e-10, eventually |α| ≤ 3.43e-7.
- fin12f4 (`V/V33_VERIFIER.md` table): sqrt_LF = 9.8717e-10, alpha_F = 6.8665e-7, both ×2. fin11 cannot be kept: 2·K0(Ê) = 4.68e-21 exceeds its FLUX0 = 4.1e-21.
- So the Main Theorem as drafted asserts a decay floor the v3.3 chain does not prove. Fix: 9.88e-10 and 6.87e-7. "Only slopes, sJ4 and F0 move" in `V33_INDEX.md` §1.1 and `main_v3_CHANGES.md:48` is wrong for the same reason.

### Must-fix (each is small; none moved a margin in my reruns)

**M1. Outer-ramp K° uses the wrong Y_{j+1}.** [checked] `B/t5_zone3/ramp_U.py:59–75, 189` builds the interpolation remainder with the contract-convention Y (M1z = M1, no φ, no outer-foot adds): (457.075, 354.489, 1299.674). The Ĝ that Lemma U uses has the primary Y: (457.092, 356.940, 1320.219) (`u1_ghat.json: cases.primary.Y_next`). Step 8 of the proof needs the same Y in both. Effect: K°₂₂ low by 1.6% of its ρ-part; Δ^R₂ low by about 2e-7 of 0.01148 [est]. Fix: two-pass (ramp_U foot adds → u1 → ramp_U K with u1's Y).

**M2. `assemble_U.py:113` drops the coupling term of the rule `EA_of_N2`.** [checked, re-run] With the term restored (scratch `assemble_U_coup.py`), Ê₂ goes 0.0565883051 → 0.0565883053 and (K1)–(K3) pass with the same margins. The exported "certified" Ê is low by 3.5e-9 relative.

**M3. Typed constants in the certified chain.** [checked]
- `B/t2_A_v33/rigor/assembly_rep.py:142–144`: ramp values `0.4328`, `0.00116`, `2e-5` are typed, not read from the verifier. `0.00116` is **below** the verifier's sI1 = 1.160338e-3 (fin12f4; lemmaS.tex:17 says 1.161e-3). Harmless against the box 0.11, but it is a wrong typed number inherited from v3.1 (`writeup/lemmaS/rigor/assembly_arb.py:270`). The check "verifier ramp sups ≤ ramp box" is only a reported margin, never asserted.
- `B/t5_zone3/zone3.py:24` types ā = 2.19103e-7; l.132 types c_δU = 2e-8 as [hyp] (u34 certifies 1.7187e-8 and asserts < 2e-8, so the value is fine).
- `B/t4_cert/u2_K.py:40` (`A_W = 0.3`) and `ramp_U.py:181` (`(0.0, 0.3)`) use the double 0.3, which is 1.1e-17 below 3/10: the certified cover does not contain the weight parameter of the proof.

**M4. Stale provenance in the pinned records.** [checked]
- `u1_ghat.json` and `u34_primary.json` were produced with the pre-g13 `t4lib.py` (hash e05cadd6…) and pin the **fin11** file `step3nl/verifier/final_alpha_points_5e-8.json`. The strings are identical in fin12f4 (2.19103e-07, 3.33811e-09), so no number changes; rerun for the record.
- `U_cert_v33.json` pins `zone3_interface.json` 9611888c…, a file no longer on disk. I re-ran the assembly against the current file (0d6338d0…): output identical apart from the args block.
- `u5_enum.py` imports `t2_A/rigor/ramp_rep.py`: same sha as the v33 copy, and `base_ibp/rigor/weights_exact.py` is also identical. It is a cross-check, not a proof input. It does not matter.

**M5. Proof-text corrections in lemmaU.tex.** [checked]
- Step 1 (l.272), dissipation reserve: "the reserve is used locally … disjoint supports" is wrong as written. The foot terms multiply the **global** ‖z₄‖_Ξ, and the level-2 budget has exactly zero slack (`u1_ghat.json: budget_slack[2] = 0.0`), so two separate Young steps would overdraw it by a_min/10. The conclusion is still right: X_in² + X_out² ≤ X² gives one Young step on √(a² + b²)·N₄, which costs exactly the two data terms the scripts charge. Write that argument.
- Step 0 (l.258): the split should read η²χ²(Θφ)^(m), not η²(χ²Θφ)^(m).
- l.142, l.466, l.479: "(W2) … not rerun [hyp]" is out of date; tail_cert and closure.tail were rerun for fin12f4.
- §"What remains a hypothesis" and (K6) still say (H2f) and the fixed-L data are pending.
- (A)'s core uses sup|z₁ᴬ| ≤ 0.52 (`weights_t2.py:22`, "B2′"). B2′ gives only |y f_y| ≤ 0.52, hence 0.5333. The bound does follow from E^A ≤ 2Ê by (A)'s own Sobolev step (√2·S1v = 0.146), but no text says so. This is the same defect task 4 found in Lemma U (M1z).

**M6. Margins quoted as [cert] at the estimated N2_ref.** [agent] Flux ×3.505/×6.196, S0v ×2.143, core M1 ×4.442 are evaluated at N2_ref = (0.1, 0.2, 5.0) [est]; Lemma N's certified N̄² = (0.0715, 0.330, 5.034) exceeds it at k = 1, 2. The dominating certified run exists (`asm_design_N2ref4.json`: flux ×3.472, S0v ×2.137, M1 ×4.415). (K1) holds; quote those numbers and have assemble_U assert them.

**M7. main_v3 wording against what v3.3 delivers.** [agent; I confirmed (W3)/(W4) are verbatim from v3.2]
- (W4) "|f| ≤ ½ on y ∈ [48, R(τ)]" and "graph over Lemma O cylinders on r ≥ 150λ₁": v3.3 gives |f| ≤ 0.11 / 0.15 / 0.1003 / 0.2276 on [28,60] / zones I–II / [253,379.5]λ₁ / [253,700]λ₁ (the last only for τ − τ₁ ≤ 4.07), cone conclusion regions from 251.44λ₁, and U-fine Sobolev bounds on [150,253)λ₁. R(τ) and S₁^f are not defined in v3.3.
- (W3) "B̂ at all centres" means centres ρ ≤ 253λ₁.
- (H4) still says graph on y ∈ [48, 350]; (H2f) needs f on [60, 700].
- l.316 quotes phase-2a numbers (flux ×4.9/×5.3, "all margins at least ×2") that contradict the erratum (ramp M1 ×1.185; flux ×0.73 without (H2f)).
- l.130: Lemma N derives H-O1 and the L² form H-O2′, not pointwise (O3).
- Several roundings are not outward (e.g. loop gain "≤ 0.0390" for 0.039034; c_δU 1.7186e-8 for 1.718604e-8).

**M8. Erratum: scope of "×1.1245".** [agent] It is the flux line only, computed with the A5 ramp box (0.52, not the closed 0.9), the C^{1,1} Ỹ, and (H2f) at (5e-7, 5e-6, 1e-5). S0v, M1 and Ĝ at (O3) are "<1% [est]". ×0.7317 is hand arithmetic from a probe JSON, so "[cert]" overstates. (H2f) for the closed chain is certified only in `FINAL_manifest_v33.json`, and against 1e-6; no script asserts G₀ ≤ 5e-7 (the stored value 3.14e-7 implies it). State the closed-version claim as "flux line certified under (F1) + (O) + (H2f), with the ramp box of A5".

**M9. Step-4 hand-off items that are labelled certified but are not fully so.** [agent; not re-checked by me]
- `S4/t9h/win.py:42, 57` and `S4/v33/h2f_win.py:99`: the cell Taylor half-width uses `float(mid h)·1.0000001`; element lengths depend on t_c, and nothing asserts rad(h)/h < 1e-7. Affects (H2f) and the v3.2 (H2), (H3) base values.
- `win.py:83, 94, 106`: the v3.2 (H2)/(H3) base energies are float sums at the unwidened base frame; `h2t_v33.py:22` adds one nextafter. Harmless against ×5.7 and more, but not outward. (H2f) itself is done correctly (widened balls, arb sums): ×3.19 / ×4.98 / ×4.11.
- `S4/v33/t9_lemmaO_v33.py:98–99`: the map from centre r_c to base arclength ("x ≈ r + 2.26") is typed, not asserted (about 0.3 slack [est]).
- Far-family H4 exists for centres 376 and 602 only (`S4/t8se/h4loc.py:117`); centres in (490, 602) follow by a monotonicity argument that is not written, and centres above 602 have no certificate. Inherited from v3.2.
- (H2f)₀ is 5e-7 in `STATEMENTS_r2.tex`, `R2_DECISIONS.md` and `zone3_interface.json`, and 1e-6 in lemmaU and the code. The certified 3.14e-7 passes both.

### Minor
- `lemmaU.tex` is not compiled (balance check passes: I re-ran it).
- `zone3.py:189–191`: the far-cylinder K₃ uses a fixed list of 4 pieces; a pure far cylinder needs 3 far pieces (conclusion regions are 226(1 − e_far) long). The 4-piece value 0.161 dominates every case (my arithmetic: at most √(0.104² + 3·0.0458²) = 0.131). Enumerate the cases in code.
- `t4_cert/run_all.sh` is the pre-v33 chain; `t4lib.T2` still points at `t2_A`.
- `export_t2.py:41` computes `rule_reproduces_EA_ref` but nothing asserts it (it is True).
- Erratum A2 cites lemmaS.tex:209 for a sentence that is not there; `hjw_rep.py` cites rhs_iv.py:198 for line 191. [agent]

## The questions asked

1. **Lemma U.**
   - **G6 (point terms): passes.** [checked] The tables have weight-derivative order m ≤ 3, and m = 3 occurs only on pairs (j, j+1) (my probe of `tables_exact.pkl`). `u2_K.py` and `ramp_U.py` bound every rem_m, including (η²)‴ and (χ²)‴, directly (`coef_arrays`, `coef_Z`; no m = 4, no integration by parts). The global identities (L⁴, {2,3,3}, interpolation) are applied to the positive weight Ξ after the pointwise bound, and differentiate the cutoffs at most twice. The outer cubic classes go by Sobolev sup, not by parts.
   - **G5 (N₃²): passes.** [checked] `zone3.py:54, 68` uses R4 of `lemmaR_consts.json`, which bounds the time average over [0, t] and is increasing in t (`lemmaR.tex:31, 47`), so horizon × N₃² bounds the integral over any shorter interval.
   - **Open–closed and new centres: sound.** The energy is jointly continuous in (τ, ρ) and a newly launched centre is handled by uniform continuity next to (τ*, 56λ(τ*)).
   - **Lemma D, window form: the lemma is correct.** Its input (D1)/(D2) is B1 above.
   - **M1z = 0.534: holds in all three regions** (core 0.5333; zones I/II 0.285; zone III 0.508), asserted in arb.
   - Every term of Prop. U that I traced into `assemble_U.py` is there with the stated constant (launch, inflow, K columns, X-terms, S33/T223, couplings, b-term, outer ramp, forcing).
2. **Design choices.**
   - **A5 (ramp box = core box): legitimate and not circular.** The verifier (D1) consumes B1 and B2′ at the hypothesis level and outputs the ramp sups; (A) (D3′) consumes them. The ×1.185 is a derived inequality, not a bootstrap margin, so any factor above 1 suffices. Conditions: M3 (assert it against the verifier output) and B1 (write it).
   - **Far family ≥ 490: legitimate.** The removed conclusion regions [263, 377] are covered by the fixed-L members ([251.44, 428.56]); far coverage starts at 379.25 ≤ 379.5; H6 is recomputed for the new neighbours.
   - **H1 arc 1.3·XL: legitimate.** [agent] The v33 check requires ξ < −XL and ξ > XL at the arc ends (actual −179.1, +170.0 against 150) and takes sup|h₀| over |ξ| ≤ XL only. Erratum D1 is right: v3.2 never tested the arc ends, but its certified H8 distance (0.136 L) implies the end condition; `step4_proof.tex:715` needs the wording fix.
3. **Circular dependency: not circular.** `zone3.py` uses the U-fine B̂ only in H6 (l.204–207); `assemble_U.py` reads only ΔR, ΔR_b, n_III and b₁, none of which depend on it. In the bootstrap, D5a → U → U-fine → D5b is the usual hypothesis-level order. Only the file hashes loop (M4).
4. **u5_enum reads t2_A: it does not matter** (M4).
5. **Erratum:** M8.
6. **Certified versus estimated:** M1–M4, M6, M9. Thresholds against the manifests: (H2f) and the fixed-L cone data are in `FINAL_manifest_v33.json` and match what lemmaU uses (pins verified on 10 files [agent]); (H3) is in the closed manifest and its sliding windows contain all five windows of the launch lemma, including [253.125, 379.69]; the tightened (H2) is B0; α₁ was not traced.
7. **main_v3:** B0, B2, B3, M7. "Unconditional, subject to review" is not yet justified; it should stay "conditional on the items of this review".

## What I verified and how
- Read in full: `lemmaU.tex`, `assemble_U.py`, `u1_ghat.py`, `u2_K.py`, `u34_sweeps.py`, `ramp_U.py`, `zone3.py`, `z3lib.py`, `phase1/ghat/ghat_arb.py`, `phase1/ramp/ulib.py` (incl. `bb_sup`), `hjw_rep.py`, `export_t2.py`; the reports of tasks 2–7; `STATEMENTS_r2.tex` §4–§7.
- Re-ran `assemble_U.py` in `t8_review/scratch/t4copy/` against the current inputs: byte-identical to `U_cert_v33.json` apart from the args block (B̂, Ê, late M1 ×2.358, M0 ×4.900, present ×2.775 / ×2.225, loop gain 0.0390). Re-ran it with the coupling term (M2).
- Checked by hand: Lemma T(i), Lemma dU, the launch lemma (both parts), Lemma N's covering, the layer-cake corollary, the comparison lemma, Lemma S and its use for S33/T223, the forcing integral, the L⁴/S33 foot algebra and the dissipation budget, Lemma P's signs against `a3_signs.json`, the H-J-W lemma.
- Two sub-reviewers (Opus, read-only) audited the step-4 v33 stages and the erratum/main_v3; I confirmed their blocking claims against the cited lines (B0, B2, B3) and mark the rest [agent].
- Checksums: `t2_A_v33/out`, `t4_cert/out`, `t5_zone3/out` and `SHA256SUMS_scripts` all verify.

## What I could not check
- The inside of (A): `weights_t2.py`, `ramp_rep.py`, `outer_L2.py` and the base chain line by line (see B1). I checked their interfaces, boxes and the sign logic only.
- The verifier itself (closed, double-interval code); I checked only that fin12f4 differs from fin11 by the flux inputs and that its outputs are what the chain reads.
- `piece_norms` in `zone3.py:56–70` (conversion from cylinder graphs to ray derivatives) against `lemmaR.tex` §4, beyond the normalisation.
- The weight constants Γ_m (`weight_consts_arb.py`) beyond their stated domain.
- Whether (P1) holds with the v3.3 norms (B2).
- The size of ψ̂₃′ − 1.3056·y^p on [14, 60] (B0): not computed.
- Step-4 internals beyond the sub-reviewer's read: the n3 boxes up to x = 730, α₁ in the manifest, the cause of the 1e-9 disagreement of `t9_lemmaO_v33` with v3.2, whether Lemma O may be applied arc by arc.
