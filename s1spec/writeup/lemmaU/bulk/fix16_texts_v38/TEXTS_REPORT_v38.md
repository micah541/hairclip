# TEXTS_REPORT_v38: the step-3 text set v38 (verifier fin13v4c, (A) chain v37, flux v37)

Nov 20, 2026. Directory `D` = `/home/claude/s1spec/writeup/lemmaU/bulk/fix16_texts_v38/`.
- **Work order:** the closing list §5 and the findings of `fix15_review/REVIEW_v4c_interface.md` (FREV15), plus the open fix-11 items B1-B5, B7 and C-7 of `fix11_review/REVIEW_v4b_v36_proofs.md` (FREV11).
- **What was touched.** Nothing outside `D` was written, no existing file was edited, and nothing in the proof tree was run.
- **Resources.** One process at a time, `nice -n 19`, `python3 -B`, PYTHONHASHSEED=0. The full rebuild with all checks (`tools/build_all_v38.sh`) takes 4 s at 50 MB peak RSS.
- **Not compiled.** LaTeX is not installed, so nothing is compiled.

**Step 3 is not closed.** The texts now state one bootstrap set with one closure proof and write out (P2-a). Every current number points to fin13v4c / chain v37 / flux v37. What remains is listed in §5.

## 1. File map
| new file (D) | built from (unchanged) | how | review items |
|---|---|---|---|
| `bootstrap_v38.tex` | (new) | hand-written. Definition BS (GS, B1 with its box part, BF, B2', B3U, B4, G); member-by-member strict-improvement table with certifying file and key; order table; Proposition BS:prop:close (nonempty, implication, closed, open, end of interval, joint continuation with (P2')) | FREV15 D-1, D-2, D-5; §5 items 1b, 9 |
| `C1_identification_v38.tex` | `fix14_interface/C1_identification.tex` | `tools/build_interface_v38.py`. **New Lemma (P2-a) with proof** (§C1:sec:P2a); Remark C1:rem:C1; Y4 label; x/y wording of Lemma kP2 | C1-1, C1-2, C1-3 |
| `GS_lemma_v38.tex` | `fix14_interface/GS_lemma.tex` | same builder. fin13v4c, `GS_boot_v4c.json`, closure.py l.145/164, frame-shift citation, B2' input, sliver and GS-6 notes, GS as a member of BS | GS-1, GS-2, GS-3, GS-5, GS-6 |
| `U1_regularity_v38.tex` | `fix10_proofs/U1_regularity.tex` | same builder. X_S = 60; GS-circ and the open set O; B1 as Lambda <= L + 2C'; (S4.1), (S4.2) by C-1; (S4.3) GS by GS-1; gap items 1 and 4 closed | D-1, D-5, C-5 |
| `continuation_v38.tex` | `fix12_texts/continuation_v37.tex` | `tools/build_texts_v38.py`. GS definition and lemma replaced by the fix-14 text (the v37 sketch is withdrawn); BS and Proposition BS:prop:close; (S4.1) by C-1; (W1) ratio; v37/v4c numbers | D-1, D-3, GS §6, A-2 |
| `lemmaU_v38.tex` | `fix12_texts/lemmaU_v37.tex` | same builder. Chain v37 hashes and numbers; BS_U = BS (BF, GS added; order rows); n1_chart of chain v37; D^c_2 from its key; B5 and B7 wording | D-2, D-3, B5, B7 |
| `lemmaS_A_v38.tex` | `fix12_texts/lemmaS_A_v37.tex` | same builder. Chain v37 sources; B1, B4, C-7 notes; claims_failing (B5); rounding statements | D-3, B1, B4, B5, C-7, B-v37-3 |
| `ei_lemma_v4c.tex` | `fix12_texts/ei_lemma_v4.tex` | same builder. fin13v4c (H8 hypotheses, D-2a box, A-1 budget x1.65694); **D-4 density sentence** at the level-2 Garding step; list of integrations by parts with their regularity; audits done/not done; settled E0 conflict | D-3, D-4, C-6 |
| `E0_LEMMA_v4c_v38.tex` | `step3nl_v4c/e0_lemma/E0_LEMMA_v4c.tex` | same builder. R-3 ratios of fin13v4c; d_ab rounded up | A-6, I.12 |
| `flux2_lemma_v38.tex` | `fix12_texts/flux2_lemma_v37.tex` | same builder. v37 reassembly numbers (`fix13_v37chain/flux/out/flux2_cert.json`); status (C8 = U-1, reviewed); bootstrap items as members of BS; C3/C5/C6 values | B-v37-2, D-3, §5 item 3 |
| `G1_export_k2_v38.tex` | `fix10_proofs/G1_export_k2.tex` | same builder. Chain v37 paths and numbers (table at 4 N2_ref from `asm_design_N2ref4.json`) | §5 item 2 |
| `main_v3_v38.tex` | `fix12_texts/main_v3_v37.tex` | `tools/build_main_v38.py`. v38 status box; ratio of the tail alpha-condition; chain v37 hashes; open list | D-3, A-2 |
| `step3_v38_main.tex` | `fix12_texts/step3_v37_main.tex` | same builder. `\input`s `bootstrap_v38` in place of the old bootstrap section; (H6) by C-1 and (GS) in the hypothesis list; verifier section on fin13v4c; "What remains" | D-2, D-3 |
| `erratum_v38.tex`, `ERRATUM_v38.md` | `fix12_texts/erratum_v37.tex`, `ERRATUM_v37.md` | `tools/build_erratum_v38.py`. Parts 1-2 verbatim; **Part 3**: K (D-3, A-1, D-2a, A-2, C-1, the GS sketch, bootstrap, D-4, records, text slips, other), L (status), M (conflicts settled and remaining), N (open items) | work-order item 3 |
| `numbers_A_v38.tex/.json` | (generated) | `tools/gen_numbers_A_v38.py`, successor of `gen_numbers_A_v37.py`, sources repointed. `tailAlphaMarV` is now the **ratio** lhs/rhs; new macros for A-1, D-2a, GS, C-1, flux C0/C3/C5/C6, BPhi exports, D^c_2. Each macro carries file, key and rounding direction in the JSON | A-2 |
| `ei_constants_v38.tex/.json` | (generated) | `tools/gen_ei_constants_v38.py` from `step3nl_v4c/v4cert/out/constants_v4.json` (byte-identical content to v4b; table unchanged) | |
| `b12_preamble_v38.tex`, `chart_lemma_body_v38.tex` | v37 copies | copy | |
| `RECORDS_ERRATA_v38.md` | (new) | records corrections R1-R8 and the recommendation for `runs/P_A`, `runs/reg_fin13v4b` | A-1r, A-2r, A-3r, B-v37-1/3/4, GS-1, B1-B4, B6, B7, C-7 |

### Tools (`D/tools/`)
- `b38_common.py`: the builder class with asserted-unique replacements, a generic retag map (fin13v4b -> fin13v4c, fix6_v36 -> fix13_v37chain, v36 -> v37 hashes, file names), protected comment lines, and markers for intentional historical references.
- `build_interface_v38.py`, `build_texts_v38.py`, `build_main_v38.py`, `build_erratum_v38.py` (with `erratum_part3_v38.md`).
- `numsub_v38.json`: 80 literal-number replacements. Each row names its JSON key and direction; produced by `make_numsub_v38.py` from the checker's suggestions and reviewed by hand (wrong key associations removed, manual rows added).
- `build_all_v38.sh`: the whole rebuild with checks; it is reproducible.

## 2. The two pieces of mathematics owed
### (P2-a), FREV15 C1-1 (`C1_identification_v38.tex` §C1:sec:P2a, Lemma C1:lem:P2a)
- **Statement.** For a smooth symmetric CDF flow on [t96, T') starting at the P1 curve, there are tau_e in (tau_a, 27.15] and continuous (lambda, c, t), C^1 on each slab, with t_tau = lambda^4(1 - theta), |theta| <= 3.37e-8, the global normal-graph representation, the constraints, sqrt(Lambda_e) <= Ybar, the far bounds and the boxes on [tau_a, tau_e). Either tau_e = 27.15 with t(27.15) < T', or t -> T'.
- **Proof** (about 30 lines): (a) local structure over the closed base (foot-point map, implicit function theorem, D_pF = M + E); (b) the set S and uniqueness of the path; (c) the conclusions of Steps 2-6; (d) limits of (lambda, c, t) at tau_e from the Lipschitz bounds of ln lambda and c; (e) closedness of the representation by Arzela-Ascoli with J >= 1/2; (f) openness by the IFT at the limit with m-bar ||E|| <= 1.04e-4.
- **Labels.** Steps (a)-(f) are **[hand]**. Steps 2-6 of Theorem P2 and the switch lemma are **[cited]** with step 4's [C]/[D] labels. sup|kappa-tilde e| <= 1.743e-3 is **[cert]** (`C1_ke.json`).
- **Flag** (Remark C1:rem:C1, [hand], for the final review): the base is "continuous across slabs" in tau. If it is not C^1 there, (lambda, c, t) are C^1 only slab-wise. Nothing downstream needs more.

### One bootstrap set and one closure proof (`bootstrap_v38.tex`)
- **Members** (Definition BS). Each comes with its strict improvement and the file that certifies it (table §BS:sec:strict):

| member | improved to | file |
|---|---|---|
| GS | 5.3903e-4 / 0.040936 | `fix15_review/c/GS_boot_v4c.json` |
| B1(a) | tail alpha ratio x1554.45 (x1512.16), C5 | rate tables; `flux2_cert.json` |
| B1(b) | Y <= L with pad (C3 x1547.6) and cushion | `flux2_cert.json` |
| B1(c) box part: alpha <= 2.5e-4, [0,8] 1e-4, IR 1e-3, 0.02 boxes, 1/2 box, **D-2a 0.01 box**, validity box | asserted ratios <= 0.79 | `box_hypotheses.max_ratio` |
| BF (F1) (BPhi) | factor 2 | Lemma D, G1; `flux2_data.json` res.export |
| BF (F2) a_bar <= 1e-6 | 2.14039e-7, x4.67 | `modsups_fin13v4c` |
| B2' | S^v_0 + sup\|d\| x2.116, S^v_1 + sup\|Dd\| x5.745, core M1, M0 | `asm_design_N2ref4.json` |
| B3U | U_cert margins | U_cert (chain v37) |
| B4 | H6 x1.860 | `zone3_interface.json` |
| G | slopes < 1/2 | Lemma G |

- **Order table.** C-1, GS-1, (P2'), U-1, D1, flux, D5a, (A), U-crude/fine, D5b, G, GS-2.
- **Proof.** Proposition BS:prop:close covers every member, including the joint continuation GS / (P2') and the end of the interval. The other texts refer to it: lemmaU_v38 (BS_U = BS), continuation_v38, flux2_lemma_v38 (Definition F:def:boot as members of BS), U1_regularity_v38 (the closed set of (P2') is that of Definition BS), ei_lemma_v4c (the hypotheses are the box part B1(c)), GS_lemma_v38.

## 3. Checks
### Balance (`tools/check_balance_v38.py` -> `tools/balance_v38.json`)
- Braces, `\begin`/`\end`, `$` parity and `\ref`/`\label`: **PASS for all 18 files** and for the `\input`-expanded `step3_v38_main.tex`.
- In the expanded document all references resolve and there are no duplicate labels.
- Standalone files list cross-document references (BS:, C:, F:, EI:), as in v37.

### Numbers (`tools/check_numbers_v38.py` -> `tools/numbers_check_v38.json`): **PASS**
Improvements over v37:
1. Every ledger value carries file and key path. A match is **KEYED** when the line names that key, **FILE** when it names the file, and **VALUE** otherwise.
2. Control words are removed before tokenising. **The v37 checker never saw numbers written after `\le`, `\ge`, `\times`** (RECORDS_ERRATA R8).
3. The record filter is narrower.
4. **SUSPECT**: a value-only match that also matches a superseded value whose successor key (v36 -> v37, v4b -> v4c) differs.
5. The 122 macros are checked against their sources' sha256.
6. **30 pinned rows** (`tools/keymap_v38.json`) map a printed literal to file + key + rounding direction and are verified (30/30 pass).
7. STALE suggestions are written out with the successor key (`stale_suggestions_v38.json`).

Totals:

| status | count |
|---|---|
| KEYED | 365 |
| FILE | 109 |
| VALUE | 1043 |
| WHITELIST | 137 (135 entries with reasons; 98 INHERITED-UNVERIFIED step-1/2/4 or v3.1/v3.2 quotes, mostly in main_v3) |
| RECORD | 26 |
| HISTORICAL | 654 (erratum, by design) |
| STALE | 0 |
| SUSPECT | 0 |
| UNMATCHED | 0 |

Limits:
- VALUE matches prove only that a value of that size exists.
- Numbers with fewer than 3 significant digits (8638 tokens) are not checked.
- The erratum is not checked.
- The CURRENT ledger includes the old audits fix3/fix5 and fix4_r3, whose constants are still cited.

Found and corrected by the check (erratum K10):
- 0.4390 -> 0.4391 (lemmaS_A);
- 7.12e-3 -> 7.13e-3 (continuation);
- S_1(2 Bhat) 0.28578 -> 0.28593 (step3 main): v37 printed a value below the certified one;
- about 75 valid but superseded v36 literals updated to v37 (keys in `numsub_v38.json`).

## 4. Conflicts of `fix12_texts/TEXTS_REPORT.md` §3 and `ERRATUM_v37.md` §I
Settled by the reports, with the resolution named in `ERRATUM_v38.md` §M:
- I.4 TC2: both quotes are outward.
- I.5 KAP_OR: the text uses 1.5884.
- I.8: superseded by A-1.
- I.10 / v37 3.1: settled by E0_LEMMA_v4c A-3.
- I.11: G1 reviewed.
- I.12 / v37 3.2, 3.4: up-rounded v37 values.
- v37 3.3: Delta^R_1 6.7972e-6.
- Notes C3: the D^c_2 key exists, 0.0552015.
- Notes C4: kappa_J with the pad.
- Notes C6: n1_chart of chain v37.
- Notes T5: k = 2 range and weight.

**Remaining** (records or process; no consequence for a quoted bound):
- I.1, I.2, I.3, I.6, I.9, I.13-I.15;
- Notes C5: old "complete/unconditional" sentences in main_v3, covered by the status box;
- Notes C7 and T1: `n35_asserts.py`, `check_keys_v34.py`, `constants_needed_v34.json` not regenerated;
- Notes T2: two v33-era quotes in the Lemma D inputs;
- Notes T4: main_v3 source of 0.2686, x4.9, x5.3; erratum (O3) column.

## 5. Still open (step 3 is not closed)
1. **Re-pin and kits** (Micah's go):
   - step-4 stages 8, 9, 18, 19, 20 on `step3nl_v4c/runs/P` + `v4bcert/out` (incl. `strip_cubic_v4c.json`, `e0_A2_budget_v4c.json`), `fix13_v37chain/runs/prod2`, `fix13_v37chain/flux/out`;
   - add `C1_ke.json`, `GS_tau1.json`, `GS_boot_v4c.json` (with GS-2, GS-6);
   - a new final_tail driver reading `EXCEPTION_v32_*.json`, and a new manifest;
   - kits `kit_step3`, `kit_step4b`; Windows runs.
   - Note: files under `step4/p2/v38/`, `step4/p2/out_v38/` (incl. `FINAL_manifest_v38.json`) and `kit_*_v4*/` appeared while this task ran. They were written by another process, not by this task, and the texts do not cite them; they must be updated once that re-pin is reported.
2. **Records**:
   - `RECORDS_ERRATA_v38.md` R1-R3, R5, R6;
   - recommended: rerun `runs/P_A` and the regression into new directories at the re-pin, plus an exception record for the old ones; pin only `runs/P` in any case.
   - Chain-side FREV11 B1-B4 and C-7 need a chain touch.
   - `G1_check.py` has not been rerun on the v37 data. It is cited as the v36 arithmetic; the code is unchanged.
3. **Lean**: `lake build` at commit d72b744 and an n = 1 instantiation file.
4. **DKS re-read**: pp. 1234-35, (D1) and the uniform bounds on d_s^m kappa used in C-1.
5. **Step-4 wording**: (H6) by C-1, GS in the hand-off list (Theorem 4(c)).
6. **[read] and not-audited verifier inputs**:
   - DS0-DS3, K_rho, K_W;
   - W monotone on [31.25, 34.25] (now also under H8);
   - the tail monotonicity statement of `verifier/VERIFIER.md`;
   - the D-1a pins;
   - the regularity of the base in tau at slab ends (C1:rem:C1);
   - 98 INHERITED-UNVERIFIED numbers in the whitelist (step-1/2/4 and v3.1/v3.2 quotes, not produced by these chains).
7. **Author items**: (P3'), Lemma O arc by arc, the single-point statement.
8. **LaTeX compile** of every v38 file (none compiled).
9. **Final independent review** of the v38 text set, including a line-by-line check of `ei_lemma_v4c.tex` against `rhs_iv.py`, of Lemma (P2-a), and of `bootstrap_v38.tex`.

## 6. Checksums
`D/SHA256SUMS` (all files in `D` and `D/tools`, excluding itself).
