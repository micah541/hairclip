# TEXTS_REPORT_v41: step-3 text set v41 (wording items of the delta check fix20)

Nov 21, 2026. `D` = `/home/claude/s1spec/writeup/lemmaU/bulk/fix21_texts_v41/`. Source: `fix19_texts_v40/` (read only; its SHA256SUMS still verifies, and so does `fix17_texts_v39/`). Items applied: `fix20_delta_check/DELTA_CHECK_v40.md` §4, W-1 to W-5 and the cosmetic retag slip. Nothing else was changed.

**Step 3 is not closed.** No number was computed for v41 and no certified value changed. `numbers_A_v41.tex` and `ei_constants_v41.tex` are identical to v40 apart from file names. Nothing in the proof tree was run.

## 1. Edits (`tools/edits_v41.txt`: 16 edits, each asserted to occur exactly once)
| item | file | what changed |
|---|---|---|
| W-1 | bootstrap_v41.tex | *Closed* step: GS and G are first taken at tau*, then U-1 is applied on [tau_1, tau*], then the other members. Q is now "defined and continuous wherever GS and G hold" (v40: "on [tau_1, tau_mod)"). *Open* Step 1 states its order: GS (GS-3) and G (Lemma G "Open"), then U-1 on [tau_1, tau*+delta], then the scalar members, (F1) and B3U(a). Step 2 refers back to Step 1 for U-1. (4 edits) |
| W-2 | U1_regularity_v41.tex | The theorem's parenthesis says that the "GS, G only" statement is about regularity. (c): (F:eq:star) holds at every tau *wherever the hypotheses of the E_I and flux lemmas (B1(c), BF, B2') hold*. (d): (D1) holds at every tau *wherever the hypotheses of Lemma S(A) hold*. (3 edits) |
| W-3 | bootstrap_v41.tex Step 2; flux2_lemma_v41.tex, march of Prop. F:prop:abs | Tail: the TAIL_LEMMA comparison. March step: the alpha-barrier is the per-cell fixed point (closure.py:113-120, 134) and needs only Lambda <= Ly_hi_k, which holds up to and just beyond a first touch by (C2). (2 edits) |
| W-4 | flux2_lemma_v41.tex, Lemma F:lem:data sketch | "(\|alpha\| <= alpha_1)" now reads "(\|alpha\| <= 2.5e-4, box (B-a) of the E_I lemma, part of B1(c))". A note says that (F:eq:star) therefore does not depend on B1(a). |
| W-5 | lemmaU_v41.tex implication step; main_v3_v41.tex bootstrap paragraph | B1(a,b) are excepted from "strictly better"; they are re-derived by first touch (bootstrap, Open step). |
| cosmetic | lemmaS_A_v41.tex | The title says "text v41" (v40: "text v39"). The date line says "text v41; text v40; text v39: chain v37" (v40: "text v40; text v40: chain v37"). (2 REPQ edits) |
| headers | bootstrap, lemmaU | The date lines name the items applied (2 REPQ edits). Every other file is only retagged and retitled, as in v40. |
| erratum | erratum_v41.tex | Part 6 (section Q) is new (`tools/erratum_part6_v41.tex`); Parts 1-5 are verbatim. ERRATUM_v41.md and RECORDS_ERRATA_v41.md are the v40 files verbatim plus a v41 pointer or section; STEP4_PATCH_NOTE_v41.md has a one-line v41 section (no step-4 edit). |

## 2. Checks (`tools/build_all_v41.sh`: 5 s, 61 MB peak, nice 19, one process)
- balance: 23 PASS, 0 FAIL.
- macros: PASS on v41. The same check FAILs on v39, as expected (12 FAIL lines); this is the negative control kept from v40. New in v41: it must also PASS on the v40 source set, and it does.
- numbers: PASS (0 STALE/SUSPECT/UNMATCHED; 126 macros, 56 pinned keys). The only new-text number is `2.5\cdot10^{-4}` (W-4, and the erratum Q4 that repeats it), status NEW-M38. That status is a match by value (to `chart4_shell.json res.Dd_sup[1]`). The actual source is the code constant `AL = arb('2.5e-4')` in `fix13_v37chain/flux/code/chart4_shell.py:47`, which the delta check cites. Hashes: 0 new unknown. `tools/whitelist_v41_new.json` is empty.
- Two full runs are byte-identical: 57 files, all of D except `start_time.txt`. SHA256SUMS lists 24 files and verifies.
- Diff against v40 after normalising the version tags: only the listed edits, the two header comment lines per file, and the date and title retags.
- Nothing outside D changed: `find /home/claude -newer <marker>` finds only the parent directory `bulk/`, whose mtime changed when D was created.

## 3. Tools (copies of the v40 tools, v41 names)
- `b41_common.py` has new TITLE_REPS for v41.
- `build_v41.py` follows `build_v40.py` and asserts that exactly the six edited files changed.
- `build_all_v41.sh` adds the macro-check PASS on v40.
- In `check_balance_v41.py`, NAMEDIRS also contain fix19_texts_v40/ and fix20_delta_check/.
- `gen_*_v41.py`, `check_numbers_v41.py` and `check_macros_v41.py` change only in names and paths.
- In `keymap_v41.json` and `whitelist_v41.json`, the file names were retagged.

## 4. Seen but not changed (outside the requested items)
- `step3_v41_main.tex` \title still says "text v39". This is the same kind of retag slip as the cosmetic item, inherited from v40, and the delta check did not list it.
- The open list of `step3_v41_main.tex` and erratum P7 (history) still list "an independent check of the text-v40 changes" as open. Erratum Q6 says it is done (fix20). The v41 wording changes themselves have not been checked independently.

## 5. Open (unchanged)
Lean build at d72b744 and the n=1 instance; DKS (Theorem 3.1 statement, "maximal solution", (D1) and the uniform bounds of C-1); (P3'), Lemma O arc by arc, the single-point statement; Windows runs of kit_step3_v4 and kit_step4b_v4; LaTeX compile; a step-4 text version; the unaudited items and records items listed in RECORDS_ERRATA_v41.md.
