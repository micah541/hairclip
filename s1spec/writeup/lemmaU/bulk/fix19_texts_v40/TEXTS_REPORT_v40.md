# TEXTS_REPORT_v40: step-3 text set v40 (fixes of the final review fix18)

Nov 21, 2026. `D` = `/home/claude/s1spec/writeup/lemmaU/bulk/fix19_texts_v40/`. Source: `fix17_texts_v39/` (read only; its SHA256SUMS still verifies). Review applied: `fix18_final_review/FINAL_REVIEW_v39.md`.

**Step 3 is not closed.** Every status box says so. No number was computed for v40 and no certified value changed; nothing in the proof tree was run.

## 1. Per finding
| item | what changed | files |
|---|---|---|
| F-1 closure | Reviewer's option (i), no compute. B1(a,b) are marked as members *without a gap*: the implication step returns them in the same form. The Open step is rewritten. Step 1: every other member is extended to [tau*, tau*+delta] by continuity and its gap. Step 2: B1(a,b) on [tau_1, tau*+delta] by first touch (Prop. F:prop:abs; TAIL_LEMMA "Tail contraction (first touch)"), whose inputs (U-1, (F:eq:star) under BF, B2', B1(c), the cell certificates, (C0)-(C7)) hold there by Step 1. Also fixed: the table rows B1(a), B1(b) ("no gap; reproduced, not improved"), the order-table rows (P2'), U-1, D1, the implication step, and Remark (2) (GS-4 box part, as closed by the reviewer). Flux lemma: (F6) and Remark F:rem:boot. lemmaU: row D6 (l.417), proof sentence (l.446), Open step. U-1 Theorem is stated under (GS) and (G) only; B_mod is used only by (P2')(c). | bootstrap, flux2_lemma, lemmaU, U1_regularity |
| F-2 macros | `\HtwoV`, `\HtwofV` from A_interface.json `res.tight.IC` (exact thresholds); `\EtauOneV` = their exact sum, asserted inside the certified balls `IC.E_tau1`; `\GAoneV` = the record [est] value, verbatim from `fix1/b12_text/numbers_A_v33.tex` (no pinned JSON holds it; it is labelled [est] in the text). BPhi macros renamed `\BPhiCzerooneV` etc. (flux2_data.json `res.export`, rounded up as before: 1.6635e-4, 7.6840e-8, 116.3593, 0.2605999). The generator asserts letters-only names. `\input{numbers_A_v40}` added to lemmaU, main_v3 and the flux lemma's standalone preamble. The last one is a further instance the review did not list: 24 macros were undefined there. | gen_numbers_A_v40.py, numbers_A_v40.tex, 3 texts |
| F-2 check | New `tools/check_macros_v40.py`, run by the build. Per standalone document (15) it checks used ⊂ defined, with \input expansion and the stepmain blocks handled. It also checks that every defined name is valid. Undefined names fail if they are project macros or not on a fixed list of LaTeX/package commands. On v39: **FAIL** (4 invalid names, 8 documents; `tools/macro_check_on_v39.log`); the build asserts this. On v40: **PASS**. | tools |
| F-3 | (H6) now reads "smooth closed immersed, symmetric about the x-axis". It names its source (C-1 (i), (iv); step-4 symmetric Γ0 + uniqueness; step4_handoff item (3)) and its uses (U-1 (S4.2), the even verifier class). The main_v3 hand-off table row was changed too. | step3_v40_main, main_v3 |
| F-4 | T_max is a length everywhere. In continuation, [t_1,T_max) became [t_1,t_1+T_max), and T_max = t_1+∫λ⁴ became ∫λ⁴; a convention paragraph was added. In GS_lemma, t(τ_b) < T_max became < t_1+T_max. | continuation, GS_lemma |
| F-5 | The criterion is now the *statement* of DKS Thm 3.1 (quoted from P1_CONFIRMATION §1). The (D1) route is kept as Remark U:rem:D1. An author-confirmation note says that "maximal solution" is read as maximal among smooth solutions. (D1) is still used by C-1 and is still listed for re-reading. | U1_regularity, continuation, step3 main, main_v3 |
| F-6 | n^III level 0 is printed 0.071444 in step3 main (certified 0.0714437…). | step3_v40_main |
| F-7 | (H4) on [14,700] is now cited as certified by stage 12 (win_arb.py item (e): cos φ>0 on every cell of [14,60], [50,380.25], [60,700]). Only step-4's wording stays open. | continuation, step3 main, step4_handoff (+ stage-12 row), patch note |
| R-1 | "629 paths and 55 scripts" now reads "992 pin entries on 646 distinct files: 629 stage paths and 55 scripts, 38 of them also stage paths". The counts are computed from the manifest by gen_records_v40.py. The 38 are scripts that are also stage inputs; there are no relative/absolute double spellings (646 either way). | 9 texts |
| R-2, R-3 | Records items for the transitive verifier-code pins; the four missing open items (fix15 A-4, fix11 C-8, fix11 D-6a, fix9 §8 closure.py:130). | RECORDS_ERRATA_v40.md, step4_handoff §1 |
| erratum | Part 5 (section P) lists all of the above. | erratum_v40.tex, ERRATUM_v40.md |

## 2. Checks (`tools/build_all_v40.sh`, about 5 s, 61 MB peak, nice 19, one process)
- balance: 23 PASS, 0 FAIL (braces, environments, $, refs/labels across the set, control characters, inputs, paths in new text).
- macros: PASS on v40; FAIL on v39 as expected (12 FAIL lines).
- numbers: PASS (0 STALE/SUSPECT/UNMATCHED; 126 macros, 56 pinned keys). New-text numbers: 16 matched to manifest_v38, 4 to other current certified files, 6 whitelisted with reasons (380.25 is a win_arb.py range constant; 0.0466 is from MOD_IBP [read]; 0.071446 is history), 0 missing; hashes: 0 new unknown.
- Two full runs are byte-identical (56 files). SHA256SUMS lists 24 files.
- Nothing outside `D` changed (`find /home/claude -newer <marker>`: 0 files; only the parent directory's mtime changed, from creating `D`).

## 3. Needs no new compute; nothing deferred
F-1 needed no certified number: the first-touch argument uses only the existing (C0)-(C7) margins and the box ratios.

## 4. Limits
- The macro check is static. It does not check that environments are defined, and lowercase non-project names are accepted only if they are on its fixed command list. It does not replace a LaTeX compile (none is installed).
- The F-1 rewrite is [hand]; nobody has re-reviewed it.

## 5. Open (unchanged in kind)
- Lean `lake build` at d72b744 and the n=1 instance.
- DKS: confirm the statement of Thm 3.1 and "maximal solution"; re-read (D1) and the uniform bounds for C-1.
- (P3'); Lemma O arc by arc; the single-point statement.
- Windows runs of kit_step3_v4 / kit_step4b_v4.
- LaTeX compile.
- A step-4 text version applying STEP4_PATCH_NOTE_v40.
- An independent check of the v40 changes (mainly F-1).
- Unaudited: DS0-DS3, the [read] inputs of E0/E2, the D-1a pins, the tail_cert degree table, M6.
- Records: R2, R3, R6, GS-2, GS-6, the transitive pins (R-2), and the four R-3 items.
