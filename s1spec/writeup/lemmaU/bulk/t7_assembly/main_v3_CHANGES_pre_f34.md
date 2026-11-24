# main_v3.tex: changes from main_v2.tex

Task 7b, Nov 13 2026. Line numbers are in `main_v3.tex` (v2 line numbers in brackets). `main_v2.tex`, `main.tex` and `main.bib` untouched; the bibliography was not edited and no citation was added. No number was altered: every new number is quoted from certified outputs (see `V33_INDEX.md`, `STEP4_RESULT_v33_addendum.md`) and carries its [cert] label or is the T3/T6 report headline. **Not compiled** (no LaTeX here). Balance check: `python3 -I check_balance_main.py main_v3.tex` gives braces net 0, environments matched, no undefined `\ref`, `$` even per paragraph: PASS. (The same script gives PASS on main_v2.tex.)

Decision applied: the Main Theorem is **unconditional**, with a boxed status note "unconditional subject to the pending independent review (task 8) and compilation; see V33_INDEX.md". The (F)/(O) material is kept as a historical remark.

| lines v3 [v2] | change |
|---|---|
| 12 [12] | Main Theorem title no longer "conditional on (F1) and (O)". |
| 15-16 [15-16] | The `hypF`/`hypO` environments are retitled "Former Hypothesis (F)/(O) ... superseded, no longer assumed" (kept for the historical remark). |
| 33-34 [33-34] | Comment on the `\GAPSTATEMENT` macro updated (historical; theorem not conditional). The macro text itself is unchanged. |
| 56 [56] | Date: "Draft v3 of November 13, 2026 (v3.3 of step 3)". |
| 60 [60] | Abstract: "conditionally on two stated hypotheses (F1) and (O)" replaced by "subject to the pending independent review and compilation recorded in the status box". |
| 64 [64] | Status box rewritten: steps 1, 2, 3 (v3.3), 4 complete as written; theorem unconditional subject to the pending review (task 8) and compilation, pointing to V33_INDEX.md; earlier drafts were conditional on (F1) and (O); Lemma U, repaired Lemma S(A), new hand-off item (H2f) and tightened (H2); the uncompiled text and the unseen new scripts are stated; erratum at `erratum/erratum.tex`; kits pending. |
| 86-100 [86-100] | The Hypothesis (F) and (O) environments and the connecting paragraph are removed from section 1.3 (blank lines remain at 86-89, 91-99); the boxed note "Status: unconditional subject to ..." is at 90. |
| 103 [103] | Main Theorem: "Assume Hypotheses (F1) and (O)" deleted. |
| 118 [118] | Proof structure: step 3's bootstrap "closes by Lemma U and the repaired Lemma S(A)". |
| 121-122 [121-122] | Remark retitled "What the theorem uses beyond the certified computations"; item (a) no longer mentions (F1)/(O); new item (a2): step 3 v3.3 uses (H2f), tightened (H2) and the fixed-L cone cylinders from the v33 stages. |
| 127-142 [new] | New Remark "Earlier drafts were conditional on (F1) and (O)" (label `rem:hist`): what (F1) and (O) were, that both are replaced by proofs (Lemma U; repaired S(A) with (O1)-(O3) derived in Lemma N), the price (two hand-off items and M1z = 0.534 [cert]); the (F) macro text and the (O) wording are retained inside it. The "wording provisional" remark of v2 L99 is dropped. |
| 227-229 [211-213] | Table of regions: windows (B) and strip now Lemma U; cylinders: fixed-L cone family (L = 50, centres [300, 380]) and far family (centres >= 490). |
| 243, 256 [227, 240] | Dependency figure: node "Hypotheses (F1), (O)" (dashed) replaced by an ordinary node "Lemma U, repaired S(A); uniformly local far field; (H2f)"; caption says earlier drafts had the dashed hypothesis box. |
| 282 [266] | Step-3 documents line: adds Lemma U (`\dref{LU}`), statements (`\dref{ST}`), the file map (`\dref{IDX}`), the erratum (`\dref{ERR}`). |
| 284 [268] | Theorem W title: "as amended by v3.2 and v3.3". |
| 288 [272] | (H2) tightened to (1e-8, 1e-6, 1e-3); new (H2f) = (1e-6, 5e-6, 1e-5) on [60, 700]. |
| 291 [275] | (H5): fixed-L cone family (L = 50, centres [300, 380]) and far family centres >= 490. |
| 294 [278] | "If the bootstrap closes (under (F1) and (O))" replaced by "The bootstrap closes (by Lemma U and the repaired Lemma S(A); unconditional subject to the pending review)". |
| 294 [278] | (W3): v3.3 values E-hat^A = (1.1879e-6, 1.0391e-4, 0.056589) and B-hat = (4.8732e-6, 5.5667e-4, 0.36561) [cert]; the v3.1 values and |f| <= 0.093 are marked as having omitted the edge terms. |
| 302-305 [286-289] | Components list: the Lemma S(A) item rewritten as the repaired Lemma S(A) (beta_A = 2, tightened (H2) + (H2f), fin12f4; E2 0.0566, core M1 x4.44, S0v x2.14, G22 x1.92, flux x3.51 / x6.20 [cert]); the Lemma S(B) item replaced by a **Lemma U** item (Lemma T, Lemma D window form, Prop. U, closure; margins late M1 x2.358, M0 x4.90, present M1 x2.775, M0 x2.225, loop gain <= 0.0390, rate x28.2; c_dU = 1.7186e-8; M1z = 0.534; zone III box (0.1003, 0.4197, 0.0143), Delta^R_2 = 0.01148 [cert]; proof text not compiled; open points G5-G7); Lemma O/R item: fixed-L cone family and the ramp [253, 379.5]; H6 item: v3.3 factors (x1.86 minimum, cone slope x2.46, K2 x2.54, K3 x2.51, t* x2.30 [cert]; v3.2 smallest 1.69 conditional on (F1), (O)). |
| 309 [293] | Bootstrap: B2 with beta_A = 2 (B2'), B3 as B3U, B4 far centres >= 490. |
| 311-314 [295-296] | Subsection retitled "The far-field estimate: Lemma U, and the former Hypotheses (F1) and (O)"; new paragraph "Status in v3.3" (replacement of (F1), (O); new hand-off items with margins x3.19 / x4.98 / x4.11 and x124.8+ [cert]; not yet done: review, compilation, kits); a lead-in sentence that the present tense below refers to v3.1/v3.2; opens `\begin{remark}[Historical: the gap of v3.1 and v3.2]`. |
| 315-326 [297-306] | The v2 text of the gap subsection is kept as the content of that remark. Edited only where it asserted the present status: "Hypothesis (F) was then unchanged as a text ... was a consequence of the repaired Lemma S(A) under (O)" (v2 L296); "erratum \dref{ERR}" (v2 L296); "in v3.2, under (F1) and (O), the bootstrap would close as written" (v2 L303). |
| 326-327 [307] | "Proposed route" paragraph: added "This route is the one carried out in v3.3 as Lemma U" and `\end{remark}`. |
| 333 [313] | Step-4 documents line: adds the addendum `\dref{R4v33}` and states that the step-4 document is closed and unchanged. |
| 344 [new] | New item in "The chain": "v3.3 stages" (10 stages by `final_tail_v33.py`; manifest `FINAL_manifest_v33.json`, 152 paths, 23 scripts, 52 reused v3.2 files; FINAL_manifest.json unchanged). |
| 354-355 [333] | Hand-off table: (H2) row replaced by tightened (H2) (x124.8 base; x2125 and x641 with error); new (H2f) row (x3.19 / x4.98 / x4.11). |
| 359-362 [337-339] | Hand-off table: cone slope and sup\|h0''\| rows replaced by the fixed-L values (0.03029, x2.64; 0.009876, x3.04); Lemma O H4 row: 1.286e-5 against 0.0144 (x1119), far unchanged; new fixed-L H7 row (x2.30). |
| 365 [342] | Review paragraph: thinnest margins mention both the v3.2 and v3.3 cone values; the last sentence (H6 "conditional on (F1) and (O)") replaced by: superseded by the step-3 v3.3 H6 (x1.86), and the open arc-end question (v3.2 H1 script did not test \|xi\| > XL at the arc ends; v33 copy does, arc 1.1 XL -> 1.3 XL; erratum undecided). |
| 379 [356] | "What is computed": adds the repaired S(A) constants, Lemma U certificates, zone III box, and the v3.3 stages. |
| 382 [359] | Manifests and kits: adds FINAL_manifest_v33.json; states that no kit yet contains the v3.3 additions (see KITS_TODO.md). |
| 387 [364] | Step-3 kit note: "predates v3.2 and v3.3". |
| 401 [378] | Related work: "conditional on Hypotheses (F1) and (O)" replaced by "unconditional only subject to the pending independent review of the new step-3 text and its compilation". |
| 405-416 [383-389] | Open items rewritten: 1 independent review (G5-G7); 2 compilation; 3 hypothesis status (F1)/(O) removed, hand-off items assumed by step 3 v3.3 and certified in step 4; 4 erratum; 5 Lemma O H1 arc ends; 6 step-3 caveats; 7 step 4; 8 reproduction (kit B not run end to end; no v3.3 kit; the authors decide); 9 literature. |
| 424 [398] | AI-use statement: erratum cited; the far-field gap "formerly (F1) and (O)"; new sentence that the v3.3 text (statements, proof, arb certificates, zone III, step-4 stages, this assembly) was written by separate Claude agent runs, and that this v3.3 material has not yet been independently reviewed. (Also fixes the v2 grammar "is stated".) |
| 448-453 [422] | Appendix table: new rows for `\dref{LU}`, `\dref{ST}`, `\dref{R4v33}`, `\dref{IDX}`, `\dref{ERR}`. |

## Not changed (deliberately)
- The \GAPSTATEMENT macro text, the bibliography (`main.bib`), sections 2, 3 (strategy), steps 1 and 2, the Main Theorem items (i)-(vi) and their numbers, the "What is not claimed" remark, section 4.3 caveats, the "Independent reproduction" bullets other than the step-3 kit sentence.
- The Main Theorem numbers (lambda_1 interval, 2.1e-17, 5.31e-8 ...) were not rechecked against v3.3; Theorem W's (W1), (W2), (W4), (W5) keep their v3.2 values. The v3.3 verifier moved only slopes, sJ4 and F0 (x4); sup|(a,b)| = 2.19103e-07 and sup|beta| = 3.33811e-09 are unchanged to six digits.

## For the PI / review
- `\dref{ERR}` points to `erratum/erratum.tex`, written in parallel by another agent; not read here.
- Theorem W (W3) quotes B-hat from the T3 report headline (4.8732e-6, 5.5667e-4, 0.36561); the (W4) box and K_2 = K_3 = 0.18 / 0.3 are the v3.2 values, and the v3.3 box is derived inside the bootstrap. Check in review that the wording of (W3), (W4) matches `lemmaU.tex` Prop. closure.
- Statement of the hand-off item (H2f) in Theorem W uses the integral form of `HYPOTHESIS_O_r2.md` with the v3.3 value (H2f)_0 = 1e-6 (the r2 statement had 5e-7; raised by decision, `BULK.md`).

## Documentation sync, Nov 14, 2026 (`fix1/SYNC_REPORT.md`)
Backup `main_v3_pre_sync.tex`. Every number is from the final certified sources named in `fix1/SYNC_REPORT.md`; bounds rounded up, margins down. The boxed status note ("unconditional subject to the pending independent review ... and compilation", l.64 and l.90) is kept. Line numbers refer to the file after the edits.

| change | review item |
|---|---|
| (iv) of the Main Theorem: floor sqrt(Lam_F) = 9.872e-10 and eventual |alpha| <= 6.867e-7 (verifier fin12f4, `V33_VERIFIER.md`: 9.871704e-10, 6.866519e-7), replacing the fin11 values 4.94e-10, 3.43e-7. The tail rate 2.69e-3 and peak 5.31e-8 stand. | B3 |
| (H4): graph over the S_1 ray on y in [48, 700] (was [48, 350]). | M7 |
| (W1): sup\|a\|, sup\|b\| <= 2.192e-7 (was 2.19e-7, not outward); (W2): \|w'\| <= 1.161e-3 (was 1.16e-3, not outward). | outward roundings |
| (W3) rewritten: E-hat = (1.2021e-6, 1.0440e-4, 0.056686), core sup\|v_A\| <= 7.1881e-3, B-hat = (4.8735e-6, 5.5927e-4, 0.36670) at every launched centre rho <= 253 lambda_1. | M7; lemmaU_v34 |
| (W4) rewritten to the delivery of `lemmaU_v34.tex` and `continuation_v33.tex`: \|f\| <= 0.15, \|y f_y\| <= 0.52, \|f - U_1\| <= 0.16 on {y >= 56, r < 253 lambda_1}; zone III \|f\| <= 0.10032, \|r f_r\| <= 0.41974; graphs over the Lemma O cylinders with K_2 = K_3 = 0.18 (cone), 0.3 (far); the v3.2 box \|f\| <= 1/2 on [48, R(tau)] and S_1^f are no longer used. | M7 |
| Components: forcing bullet (peak <= 5.308e-8, floor <= 9.872e-10; was 5.307e-8 and 4.936e-10); Lemma S(A) bullet (final margins at 4 x N2_ref: E2 <= 0.056686, core M1 x4.409, S0v + sup\|d\| x2.116, G22 x1.917, flux x3.363 / x5.953, replacing x4.44, x2.14, x1.92, x3.51 / x6.20); Lemma U bullet (late M1 x2.353, M0 x4.892; present M1 x2.770, M0 x2.225; loop gain <= 0.03904; (K6) x1.190; c_dU <= 1.7187e-8; zone III box (0.10032, 0.41974, 0.014291), thinnest x1.238); H6 bullet (U-fine B-hat (4.8723e-6, 2.6266e-4, 0.025840)). | phase-2a / v3.3 margins superseded |
| Status of v3.3 and step-4 stage bullet: 19 stages, 281 paths, 35 scripts, 85 reused pins (was 10 / 152 / 23 / 52), sha256 8f6f7f95...; the M9 arb rows. | M9 |
| Hand-off table: (H2) base 8.02e-11 (was 8.01e-11, not outward), x124.7 / x2124 / x641 (were x124.8 / x2125); new row (H2) in f - U_1: (3.658e-10, 1.777e-9, 1.665e-6), x27.34 / x562.8 / x600.9; (H2f) base (3.127e-7, 1.004e-6, 2.431e-6), with error (3.137e-7, 1.0034e-6, 2.436e-6), margins x3.188 / x4.983 / x4.106 (were x3.19 / x4.98 / x4.11); cone sup\|h_0''\| x3.03 (was x3.04); new far-family H4 row, minimum x15.58. | M9; outward roundings |
| Open item 3 (hypothesis status): (H4) on [48, 700]; far-family coverage closed by `h4far_v33`. | M7, M9 |

## N35 documentation pass, Nov 14, 2026 (`fix2/N35_REPORT.md`)
Backup `main_v3_pre_n35.tex`. Every number from the v35 sources (A 78fe2f8e..., U 960809fd..., zone3 45e8d4dc..., manifest a0d348ec..., 20 stages); bounds rounded up, margins down.

| change | review item |
|---|---|
| (W3): E-hat = (1.2018e-6, 1.0441e-4, 0.056694), core sup\|v_A\| <= 7.1878e-3, B-hat = (4.8735e-6, 5.5933e-4, 0.36675). | N3 (v35 numbers) |
| Lemma S(A) bullet: sources `t2_A_v35`; final margins at 4 x N2_ref: E2 <= 0.056694, core M1 x4.409, S0v + sup\|d\| x2.116, G22 x1.917, flux x3.365 / x5.955; chart term (stage `n1_chart`) mentioned. | N3, N1 |
| Lemma U bullet: (K6) datum margin x1.189 (was x1.190); U-fine B-hat_2 <= 0.025842. | N3 |
| Bootstrap B3: `R(tau)` defined as v3.2's window-region end, not used in v3.3; the v3.3 form (B3U) stated. | M7 (l.309) |
| "Lemma N derives (O1)-(O3)" replaced by: Lemma N derives (H-O1) and the L^2 form (H-O2') of the outer data, not the pointwise (O3) (two places). | M7 (l.130, l.312) |
| "all repaired margins at least x2" replaced by the erratum's repaired-v3.1-setting values (thinnest: flux x1.146 at beta = 2, ramp M1 x1.185); the unsourced E_2 = 0.2686, x2.73, x4.9, x5.3 are withdrawn with a visible [TODO: source]. | M7 (l.315) |
| Step-4 stage bullets and manifest: 20 stages, 289 paths, 37 scripts, a0d348ec...; stage `n1_chart`; explicit environments (N2). Hand-off table: H4 1.2865e-5 (was 1.286e-5, not outward). | N3, minor |
| (P1) cited as Dziuk-Kuwert-Schatzle 2002 Thm 3.1 and/or Chou 2003, marked [to be confirmed by the authors] (l.122 and the open-items list). | N5 |
| Status note: re-review found no blocking item; N1 closed; remaining dependences listed. Open item 3: chart identification closed. | N1 |
