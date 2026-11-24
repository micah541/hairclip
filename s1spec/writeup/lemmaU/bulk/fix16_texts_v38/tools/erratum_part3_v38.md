
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
