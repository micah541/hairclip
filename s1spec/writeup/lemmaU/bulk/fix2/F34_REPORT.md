# F34_REPORT: final documentation pass for v3.3 (Nov 14, 2026)

No mathematics and no certified output changed. Every edited file has a backup `*_pre_f34.<ext>` next to it (reports in `step4/p2/v33/` too). Line numbers are after editing. Paths relative to `/home/claude/s1spec/writeup/lemmaU/bulk/` (`B`).

## Checks
- `fix1/lemmaU_rev/make_constants_v34.py` regenerated `constants_needed_v34.json` (104 rows, 184 numbers). `check_keys_v34.py`: **RESULT: PASS** (K6 1.18983, B2' 2.11619). Only change in the JSON: the S4 pin (backup `constants_needed_v34_pre_f34.json` reconstructed by swapping the pin back; diff is the 4 pin lines).
- `t3_proof/check_balance.py`: PASS on `continuation_v33.tex`, `numbers_A_v33.tex`, `main_v3.tex`, `lemmaU_v34.tex`, `chart_lemma.tex`, `chart_lemma_body.tex`, `erratum.tex`, and on the expanded `step3_v33_main.tex` (with its `\input` files, so it covers `lemmaS_A_v33.tex`). `lemmaS_A_v33.tex` and `step3_v33_main.tex` alone report only undefined refs, identical to their `_pre_f34` backups (the refs live in the `\input` files). Nothing compiled.
- No remaining "20 stages", "a0d348ec" or "19 stages 8f6f7f95" as a current statement in the edited .tex files and ERRATUM; the older sections of V33_INDEX, KITS_TODO, STEP4 addendum and main_v3_CHANGES keep them as history, with a new sync note on top and a new final section.

## Edits
1. Manifest pin (21 stages, 366 paths, 39 scripts, 8bc47e83...):
   - `make_constants_v34.py:189` PIN.
   - `step3_v33_main.tex:12, 41`; `continuation_v33.tex:109`; `lemmaS_A_v33.tex:19, 286`; `numbers_A_v33.tex:10` (header comment); `lemmaU_v34.tex:27, 535`.
   - `main_v3.tex:344` ("Twenty-one further stages", manifest line, plus a new sentence on stage `f12_sliver`) and `:384`.
   - `V33_INDEX.md`: sync note 3 (l.9), PIN in the text-check row (l.268), new section 4.7 at the end. `STEP4_RESULT_v33_addendum.md`: sync note 3 (l.7), new section 7. `KITS_TODO.md`: sync note 3 (l.9), new section F (kit B needs `f12_sliver.py`, driver 57c20a44..., 366 reference paths). `main_v3_CHANGES.md`: new F34 section at the end.
   - `ERRATUM.md:36` (was 19 stages, 8f6f7f95; now 21, 8bc47e83). Also `erratum.tex:152`, the same sentence (not in your list; same content, so I kept them in step).
2. F1/F2 in the text:
   - `fix2/n12/chart_lemma.tex:69` and `fix1/b12_text/chart_lemma_body.tex:58`: Remark (4) now cites the certified sI2' = 0.015456178888449044 (M1' = 28 sI2' <= 0.4328) from stage `f12_sliver`; the frame conventions stay marked inherited. `chart_lemma.tex:42` and `chart_lemma_body.tex:31`: "same convention on the sliver" now refers to sI2' of Remark (4).
   - `lemmaS_A_v33.tex:164`: one sentence that the sups G_i are certified on the extended domain [14 - 2^-19, 16.25 - X0] by `f12_sliver`, no quoted number changed. This is the only place where the G_i are described.
3. F3 (asserts): reworded to cite `fix2/n35_asserts.py` (check on the pinned outputs, PASS) and to say the `*_n35assert.py` copies with the assertion inside give identical numbers: `lemmaS_A_v33.tex:19, 87, 317`; `lemmaU_v34.tex:530`; `V33_INDEX.md:318` (the N35-checks paragraph).
4. F4:
   - `main_v3.tex:315`: the erratum table does **not** supply 0.2686, x4.9, x5.3 (checked), so it is now an explicit "Open item" sentence, marker kept.
   - `ERRATUM.md:361`, `erratum.tex:240`: the (O3)-v35 entries: no source in the table; relabelled "explicit and not blocking", marker kept.
   - Appended status notes (bodies unchanged): `step4/p2/v33/T6_REPORT.md:107`, `M9_REPORT.md:189`; closing note in `t8_review/REREVIEW_v33.md:106`.

## Remaining TODOs (markers left)
- `main_v3.tex:315` and `main_v3_CHANGES.md` (older row): [TODO: source for 0.2686, x4.9, x5.3] (withdrawn figures, unused).
- `ERRATUM.md:361`, `erratum.tex:240`: [TODO: source for (O3) values on t2_A_v35] (`o3_outer_check.py` would need a `t2_A_v35` entry in TREES; not used by v3.3).
- `KITS_TODO.md` section F: [TODO: source] for the exact kit-side input file list of `f12_sliver` (take from the `inputs` block of `out_v33/f12/f12_sliver.json`).
- Pre-existing, untouched: the open items at the end of `V33_INDEX.md` section 4.6 (`ramp_U.py` rerun position; (O3) entries), `ic_rep.py` label string, `BULK.md:120` (still 19 stages / x15.59; not in this task's list), the authors' literature confirmations (P1', P2', P3'), LaTeX compile.
