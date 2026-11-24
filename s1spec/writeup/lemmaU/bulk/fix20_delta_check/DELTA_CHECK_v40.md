# DELTA_CHECK_v40: independent check of the v39 -> v40 text delta (step 3)

Nov 21, 2026. Fable, independent checker. Scope: the delta only (`fix19_texts_v40/`, the 65 edits of `tools/edits_v40.txt`) against the findings F-1..F-7, R-1..R-3 of `fix18_final_review/FINAL_REVIEW_v39.md`. Nothing outside this directory was written; nothing in the proof tree was run. Checks were read-only Python (nice 19, < 100 MB, a few seconds).

- `T` = `/home/claude/s1spec/writeup/lemmaU/bulk/fix19_texts_v40`; `V` = `/home/claude/s1spec/step3nl_v4c`; `N` = `/home/claude/s1spec/writeup/lemmaU/bulk/fix13_v37chain`.
- Labels: **[read]** checked by reading the cited file:line; **[hand]** my argument, not machine-checked; **[cert-read]** a certified value read from a pinned JSON and compared.

## 0. Verdict

**OK. Every item F-1..F-7, R-1..R-3 is FIXED. No edit introduced a mathematical error.** Five wording remarks (W-1..W-5, one of them introduced by v40) and one cosmetic retag slip are listed in §4; none affects the argument and none needs a computation.

**The F-1 continuation argument is now correct and non-circular [hand/read].** Step 3 is still not closed (external and author items unchanged); the texts say so everywhere.

| item | verdict | note |
|---|---|---|
| F-1 closure / Open step | **FIXED** | §1; wording W-1..W-4 |
| F-2 macros | **FIXED** | §2; values and rounding correct |
| F-3 symmetry in (H6) | **FIXED** | §3 |
| F-4 T_max convention | **FIXED** | §3 |
| F-5 DKS statement | **FIXED** | §3; one precision remark |
| F-6 n^III | **FIXED** | 0.071444 >= 0.0714437 |
| F-7 (H4) range | **FIXED** | matches `win_arb.py` header (e) |
| R-1, R-2, R-3 | **FIXED** | §3 |
| new errors | **none (math)**; W-2 is a new wording over-statement | §4 |

## 1. F-1: is the continuation argument correct and non-circular?

Read: `T/bootstrap_v40.tex:36-38, 43-48, 64-67, 85-99`; `T/flux2_lemma_v40.tex:143-201`; `T/lemmaU_v40.tex:420, 449, 455`; `T/U1_regularity_v40.tex:120-174`; `V/verifier/TAIL_LEMMA.tex:67-83`; `V/verifier/closure.py:15-22, 96-185`; `T/ei_lemma_v4c_v40.tex` §Hypotheses; `T/GS_lemma_v40.tex:98-107`; `T/continuation_v40.tex` Lemma C:lem:G; `N/flux/code/{chart4_shell,flux2_data,flux2_coef}.py`.

### 1.1 Step 1 (members with a gap on [tau*, tau*+delta])

What each extension actually uses beyond tau* [read]:

| member | extended by | uses beyond tau* |
|---|---|---|
| GS | `GS_lemma_v40.tex:98` (openness with endpoint) | continuity of P ((P2')(b), needs only tau < tau_mod), smooth flow, the GS-2 gap at tau* |
| G | Lemma C:lem:G "Open" | smooth flow, continuity of lambda, slope gap at tau* |
| U-1 on [tau_1, tau*+delta] | `U1_regularity_v40.tex:143, 165-173` | GS, G, tau < tau_mod only (see 1.3) |
| B1(c), (F2), B2', scalar parts of B3U, B4 | continuity + gap | U-1 continuity |
| (F1) | extra window piece -> 0, C0_k > 0 | U-1(d) (needs G) |
| B3U(a) | `lemmaU_v40.tex:455` | joint continuity of E^rho; G, lambda |

None of these uses B1(a) or B1(b) beyond tau*. The gap of B1(c) *at* tau* comes from B1(b) on [tau_1, tau*], which is available. So the sentence "No step of Step 1 uses B1(a,b) beyond tau*" (`bootstrap_v40.tex:91`) is true.

### 1.2 Step 2 (first touch on [tau_1, tau*+delta]): do the hypotheses hold from Step 1 alone?

Inputs of the first-touch argument and where each comes from:

1. **Regularity** (Lambda in C^1, alpha in C^1, identities at every tau): U-1 under GS, G. OK.
2. **E_I inequality**: hypotheses (B-a)-(B-f), (E0), (E2), (H8) are all box inequalities of B1(c), plus S^V_0, S^V_1 of B2', plus (I3'). No hypothesis of the form "Lambda <= barrier" or "|alpha| <= alpha-barrier" [read, `ei_lemma_v4c_v40.tex` §Hypotheses]. OK.
3. **(I3') = flux lemma (d)**, hence (F:eq:star) with the certified (s, C', S'): needs (B Phi) = (F1), (F2) (typed 1e-6), E^A <= 2 Ehat (B2'), and chart/coefficient boxes. I checked the code, because the text of Lemma F:lem:data says "|alpha| <= alpha_1", which would be a B1(a) consequence: the certified computation uses the **box level** `AL = 2.5e-4` (`N/flux/code/chart4_shell.py:47, 249`: "standing boxes B-a, B-d, B-e"), `ABAR_HYP = 1e-6` (`flux2_coef.py:39`), and beta_A Ehat with the doubled exports (`flux2_data.py` header). So (F:eq:star) rests only on BF, B2', B1(c). **This is the point on which non-circularity depends, and it holds** [read].
4. **alpha-equation** alpha' <= lambda~ alpha + beta(alpha, Lambda), beta nondecreasing in Lambda: modulation lemma on the validity box (part of B1(c)). OK.
5. **Cell certificates**: static; valid for every state with |alpha| <= alpha_k^up, Lambda in [Ly_lo, Ly_hi], under the boxes. OK.
6. **(C0)-(C7)**: static, certified. OK.
7. **Start**: (H1) at tau_1. OK.

The pair (|alpha|, Y) is then handled jointly without assuming either beyond the first touch [hand, from `closure.py:113-120, 134` and Prop. F:prop:abs]:

- March, alpha: on step k, a_up is a fixed point with beta(a_up, Ly_hi)(1+1e-9)/|lambda~| <= a_up, so alpha cannot cross a_up while Lambda <= Ly_hi_k (strict by the factor 1e-9), and a_next is the linear comparison. This needs only **Lambda <= Ly_hi_k**, not Lambda <= L + 2C'.
- March, Y: at a first touch t0 of Y = L, Lambda <= L + G + 2C' < Ly_hi_k on a right-neighbourhood by (C2) (margin >= 5e7) and continuity; so the alpha-barrier persists there, the state is in cell k, and the Gronwall step with the pad (C3) gives G <= 0.
- Tail: the whole-tail alpha-condition has room x1512 (C5), so alpha' < (kappa sqrt(ell))' strictly while Lambda <= (1 + 2C'/L_F) ell; (C4) puts the scaled state in the anchor cell just after a touch of Y = ell.

So B1(a,b) on [tau_1, tau*+delta] follow from the Step-1 members only. The contradiction with the definition of tau*, and the "End of the interval" step through (P2')(c), are correct.

Two robustness facts [read/hand]:
- Step 2 re-derives B1(a,b) on the **whole** interval from tau_1, so the closedness of B1(a,b) at tau* is not even needed.
- "The right side is the continuous barrier function" (`bootstrap_v40.tex:87`) is true for both production runs: a* equals the last march alpha (3.684592528e-5, 3.580561555e-5) and L* the last L [cert-read, `V/runs/P/rate_table_fin13v4c_a0.00015.json`, `..._a0.000125.json`]. Remark (2) (`bootstrap_v40.tex:99`) also checks: L* + 2C' <= Ly_hi of the last cell (2.8076e-15 vs 2.8180e-15; 2.7096e-15 vs 2.7124e-15), a* <= alpha_up = 5e-5.

### 1.3 U-1 under (GS), (G) only

The proof (`U1_regularity_v40.tex:165-173`) uses: P in C^1 and A in C^infinity on the open set O ((P2')(a,b)); the graph function U(P; x) on O_GS (GS); far-field graph (G) for (d). It uses neither B1(b) nor the validity box [read]. B_mod enters only (P2')(c) (compactness of K). The restatement is correct for all *regularity* claims and identities. See W-2 for the two *inequalities* quoted inside (c), (d).

### 1.4 The other F-1 edits

- Table rows B1(a), B1(b) ("no gap; reproduced, not improved") and the "Members without a gap" paragraph: correct.
- Order table rows (P2'), U-1, D1: consistent with 1.1-1.3.
- Flux lemma (F6) (`flux2_lemma_v40.tex:196`) and Remark F:rem:boot (`:201`): correct.
- `lemmaU_v40.tex:420` (D6 row), `:449` (proof sentence), `:455` (Open): correct; one leftover sentence, W-5.

**F-1: FIXED.**

## 2. F-2: macros

- Values against the pinned JSON (both files' sha256 occur in `FINAL_manifest_v38.json`: A_interface c48791b2..., flux2_data f54a7774...) [cert-read]:

| macro | printed (`T/numbers_A_v40.tex:150-158`) | source | direction |
|---|---|---|---|
| `\HtwoV` | (1e-8, 1e-6, 1e-3) | `res.tight.IC.H2` | exact |
| `\HtwofV` | (1e-6, 5e-6, 1e-5) | `res.tight.IC.H2f` | exact |
| `\EtauOneV` | (1.01e-6, 6e-6, 1.01e-3) | exact sum; equals the midpoints of the balls `IC.E_tau1` (bound "up") | exact |
| `\BPhiCzerooneV` | 1.6635e-4 | `res.export.1.C0` = 1.66340772...e-4 | up, correct |
| `\BPhiConeoneV` | 7.6840e-8 | `res.export.1.C1` = 7.68394205...e-8 | up, correct |
| `\BPhiCzerotwoV` | 116.3593 | `res.export.2.C0` = 116.35923142... | up, correct |
| `\BPhiConetwoV` | 0.2605999 | `res.export.2.C1` = 0.26059987543... | up, correct |
| `\GAoneV` | (6.7e-9, 6.5e-8, 3.0e-6) | verbatim from `fix1/b12_text/numbers_A_v33.tex:74`; no pinned JSON; used only as [est] record (`lemmaS_A_v40.tex`, G-A1) | record |

- All four BPhi macros are used as upper bounds ("C0_1 <= ..."), so rounding up is the right direction. (K6): Ehat_0/1.01e-6 = 1.2004/1.01 = 1.1885, consistent with the printed x1.188.
- My own scan (independent of `check_macros_v40.py`): every `\...V` macro used in every v40 tex is defined; all 126 names in `numbers_A_v40.tex` are letters only; `\input{numbers_A_v40}` is in the preamble (before `\begin{document}`) of `lemmaU_v40.tex:25`, `main_v3_v40.tex:61`, `flux2_lemma_v40.tex:22`. The numbers file uses `\providecommand` + `\renewcommand`, so double input is harmless. `erratum_v40.tex` names the macros only inside `\verb`.
- Not done: a LaTeX compile (none installed). The static check does not replace it.

**F-2: FIXED.**

## 3. F-3 .. F-7, R-1 .. R-3

- **F-3 FIXED.** `step3_v40_main.tex` (H6) now states symmetry about the x-axis, with source C-1 (i), (iv). C-1 (i) does state the symmetry and (iv) gives (S4.2) [read, `C1_identification_v40.tex` Theorem C-1]. `main_v3_v40.tex` hand-off row changed too.
- **F-4 FIXED.** T_max is a length in every remaining occurrence; I scanned all v40 tex files for `T_{\max}` and found no absolute-time use left. `continuation_v40.tex`: "[t_1, t_1+T_max)", "T_max = int lambda^4 ((P2')(c))", convention paragraph at `:64`. `GS_lemma_v40.tex:107`: "t(tau_b) < t_1 + T_max".
- **F-5 FIXED.** `U1_regularity_v40.tex:61` quotes the statement of DKS Theorem 3.1 exactly as recorded in `writeup/lit_P1/P1_CONFIRMATION.md:20`; Lemma U:lem:bu (`:88-92`) now uses it; (D1) is kept as Remark U:rem:D1 and for C-1. Precision remark (not an error): the proof says "Phi has no smooth extension to [0, T_max]"; what Lemma U:lem:max gives is "no extension as a smooth *solution* (a family of immersions, an IsFlow) on [0, T_max]". That is the sense needed for "maximal solution", and it is covered by the listed author item A-2.
- **F-6 FIXED.** 0.071444 in `step3_v40_main.tex`; rounds 0.0714437 up.
- **F-7 FIXED.** The cited ranges [14,60], [50,380.25], [60,700] and "cos(phi) > 0 on every cell, cells consecutive" are item (e) of the header of `step4/p2/v33/win_arb.py:16-18` [read]. Texts now say "certified"; only step-4's wording stays open.
- **R-1 FIXED.** 629 + 55 - 38 = 646; `tools/records_v40.json` `pinned` has 646 entries. The builder's reading of the 38 (scripts that are also stage paths) differs from the reviewer's ("two spellings"); the counts agree either way.
- **R-2, R-3 FIXED.** `RECORDS_ERRATA_v40.md:89-96` lists the transitive pins and the four missing open items.

## 4. Remarks on the delta (none blocking, no computation)

| # | file:line | remark | suggested wording |
|---|---|---|---|
| W-1 | `T/bootstrap_v40.tex:91` vs `:93` | Order inside Step 1 is implicit. Continuity of the scalar members, (F1) and B3U(a) on [tau*, tau*+delta] needs U-1 there, which needs GS and G there first. U-1 is invoked only in Step 2. The order GS, G -> U-1 -> the rest is available (GS-3 and Lemma G "Open" need only continuity of P and the smooth flow), so this is presentation, not a gap. Related: `:87` calls Q "continuous on [tau_1, tau_mod)", but Q is defined only where GS and G hold. (Inherited from v39.) | Start Step 1 with GS and G, then "U-1 now applies on [tau_1, tau*+delta]", then the scalar members. |
| W-2 | `T/U1_regularity_v40.tex:143` with `:148` (c) and `:149` (d) | **Introduced by v40 (wording).** The theorem now says no member other than (GS), (G) is needed, but (c) still asserts "the inequality (F:eq:star) holds at every tau" and (d) "the level-j inequality (D1) holds at every tau". Those inequalities need the hypotheses of the E_I and flux lemmas (B1(c), BF, B2'), respectively of Lemma S(A). U-1 only upgrades "a.e." to "every tau". `bootstrap_v40.tex:93` states the dependence correctly, so the argument is unaffected. | "...; wherever the hypotheses of the E_I and flux lemmas (resp. of Lemma S(A)) hold, (F:eq:star) (resp. (D1)) holds at every such tau, not only a.e." |
| W-3 | `T/bootstrap_v40.tex:93`; `T/flux2_lemma_v40.tex:176` | For alpha, Step 2 cites only TAIL_LEMMA "Tail contraction". On the march the alpha-barrier is the per-cell fixed point of `closure.py:113-120, 134`; Prop. F:prop:abs covers it with one clause ("the alpha-barrier is unchanged"). Also "Lambda <= Ly_hi_k up to the first touch" should include a right-neighbourhood (that is what (C2) gives). The mathematics is as in §1.2. | One sentence: "on a march step the alpha-barrier needs only Lambda <= Ly_hi_k, which holds up to and just beyond a first touch of Y = L by (C2)". |
| W-4 | `T/flux2_lemma_v40.tex:77` | Lemma F:lem:data's sketch says "D^k d is a rate per unit tau (|alpha| <= alpha_1)". The certified code uses the box level 2.5e-4 (§1.2 item 3). With "alpha_1" the text would make (F:eq:star) depend on B1(a), contrary to `bootstrap_v40.tex:93`. (Inherited; it matters only now that F-1 leans on it.) | "(|alpha| <= 2.5e-4, box (B-a))". |
| W-5 | `T/lemmaU_v40.tex:449` (last sentence of the implication step); `T/main_v3_v40.tex:320` | Leftovers of the old claim: "every inequality of BS_U holds ... with a strictly smaller right-hand side"; "D6 (barrier) derive each bound strictly better ... a continuity argument". | Add "except B1(a,b), which are re-derived by first touch (bootstrap, Open step)". |
| cosmetic | `T/lemmaS_A_v40.tex:16-17` | Retag slip: title says "text v39", date line says "text v40; text v40: chain v37". | "text v40" / "text v40; text v39: chain v37". |

## 5. Integrity of the delta

- `fix17_texts_v39/SHA256SUMS` verifies (v39 untouched); `T/SHA256SUMS` verifies.
- Diff v39 -> v40 per file (after normalising the version tags): besides the 65 listed edits only headers, dates and "text v39" history retags changed. I read every changed line of C1_identification, E0_LEMMA, chart_lemma_body, ei_lemma, b12_preamble, G1_export, GS_lemma, U1_regularity, continuation, lemmaS_A, and the listed edits of bootstrap, flux2_lemma, lemmaU, step3 main, main_v3, step4_handoff. No number changed except the macros of §2 and F-6.
- The set's own logs: balance 23 PASS, macros PASS, numbers PASS (0 stale / unmatched).

## 6. Not checked

- Prop. F:prop:abs and TAIL_LEMMA were checked for their **dependency structure** and re-read once; their internal estimates were reviewed before (fix9, fix15, fix18) and not re-derived here.
- No LaTeX compile. Erratum Part 5 and the .md companions were read once, not cross-checked line by line.
- External and author items are unchanged: Lean build and n = 1 instance, DKS confirmations, (P3'), Lemma O arc by arc, Windows kit runs, step-4 text.
