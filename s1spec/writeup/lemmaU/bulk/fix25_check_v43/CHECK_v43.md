# CHECK_v43: independent check of text set v43 (24 Nov 2026)

Read-only check of `fix24_texts_v43/` against `fix22_texts_v42/`, `fix23_delta_check/DELTA_CHECK_v42.md`, the Lean
sources in `agents/lean/CWFlow/CWFlow/Curve/`, `lean_bridge/REVIEW_W7.md` and `lean_bridge/raw_*.log`. Nothing was
built in Lean; no file outside this directory was written (the checkers were run on a scratch copy).

## Verdict: OK WITH FIXES

No blocking finding. One should-fix (records), five wording items, one note. Step 3 is described as **not closed**
everywhere; that stays correct.

## 1. DELTA_CHECK_v42 findings

Method: v42 retagged to v43 and diffed against v43, file by file; every changed hunk maps to one of the 41 edits of
`tools/edits_v43.txt`, the header lines, erratum Part 8, or the v43 sections of RECORDS / PATCH_NOTE / ERRATUM.md.
No unlisted change.

| Item | Status | Where |
|---|---|---|
| S-1 raw logs of the S3 files | fixed; the texts cite the three logs (rebuild of 24 Nov 00:00), which exist with rc=0 | `U1_regularity_v43.tex:46`, `RECORDS_ERRATA_v43.md:118`; but see F-1, F-2 |
| S-2 "All.lean never compiled" | fixed; (PC), (MS), (GE) added | `U1_regularity_v43.tex:198` |
| S-3 erratum R3(b) | fixed (endpoint form by `isFlow_glue_endpoint`) | `erratum_v43.tex:680` |
| S-4 limit form of (P1) | fixed (sup form); no other limit form of the L² criterion in the set | `lemmaU_v43.tex:386` |
| W-1 review coverage of S3 | fixed in the status boxes, continuation, C-1 labels and (I4), U-1, erratum R4/R6; one residue, F-4 | |
| W-2 axiom checks for (T4), (T5) | fixed; matches `Criterion.lean:112-115` | `U1_regularity_v43.tex:46` |
| W-3 (MS) range | fixed | `U1_regularity_v43.tex:70` |
| W-4 [hand] steps | fixed; the C-1 (v) argument T_max = T_* − t₁ is correct | `U1_regularity_v43.tex:91-93`, `C1_identification_v43.tex:167` |
| W-5, W-6 main_v3 | fixed | `main_v3_v43.tex:349`, `:83` |
| W-7 both SHA256SUMS | fixed in all seven places | |
| W-8 bootstrap row | fixed | `bootstrap_v43.tex:69` |
| W-9 (PC) inhabited; `\tD` | fixed | `U1_regularity_v43.tex:69`, `:22` |

## 2. Lean citations against the source

- **Files.** All 16 files of `CWFlow/Curve/` hash to `snap_final/SHA256SUMS` (13 files) and `snap_s3/SHA256SUMS`
  (3 files); the ten prefixes in U-1 §1.2 match. No `sorry`, `axiom`, `native_decide` in `Curve/` (one comment in
  `N1.lean:6`). HEAD 8718f95, MM 7b15b01, toolchain v4.33.1 as cited. `Curve/` is still untracked, no remote. The 84
  modified tracked files differ from HEAD by line endings only (`git diff --ignore-cr-at-eol` is empty).
- **Statements.** Read and compared with U-1 §1.2: (L1), (L2), `cdf_unique_up_to_reparam`, (T1)–(T5), (V), (Sy),
  (Sh), (PC) with `lift_injective`, (MS), (GE), `cdf_exists_maximal_A`, `IsCDF`, `Defs.lean`. All paraphrases agree:
  hypotheses, `Icc`/`Ico`/`Ioo`, `0 < T`, `0 ≤ t₁`, ENNReal `Tmax` with `toReal` only under `≠ ⊤`, sup (not limit)
  form of the criterion, `Circle` on `EuclideanSpace ℝ (Fin 1)`, class = C^∞ on `[0,T] × Circle`, interior times only
  in (T4). No hypothesis is dropped. Not read in full: the library `rst_glue`, `isFlow_restart`, `flowVel_comp`
  beyond its statement header.
- **Compile status.** The texts say "compiled per file against the repository's existing build outputs; full
  `lake build CWFlow.Curve.All` pending; commit pending" in every file with a substantive [Lean] claim. That is
  what the records support: `raw_{LiftV,Lift3,N1,Lift,Criterion,Periodic2Circle,MaxShift,All}.log` all end rc=0, and
  each `.olean` is newer than its source. No text says "built" or cites a commit.
- **DKS.** Attribution and marked history only.
- Over-claims found: F-3 and F-5 (both inherited from v42, both wording).

## 3. New inconsistencies; checkers

Run on a scratch copy, one process, nice 19, peak 61 MB: balance 23 PASS / 0 FAIL (log identical to the set's);
macros PASS; numbers PASS (0 STALE, 0 UNMATCHED, 56 pinned, 0 NEW-UNKNOWN hashes, Lean hashes verify).
`SHA256SUMS` verifies (24 files). No number changed between v42 and v43. No new cross-reference or macro problem.

## 4. "Step 3 closed"

Searched every file for closed / closes / complete / unconditional near step 3 or the theorem. Every status box and
open list says not closed. `main_v3_v43.tex:137` ("whose bootstrap closes by Lemma U ...") is the proof summary and
sits under the status box of `:81`, which governs it. Nothing to fix.

## Findings

**Blocking:** none.

**Should-fix**

- **F-1 (records).** `TEXTS_REPORT_v43.md:6` and `TEXTS_LOG_v43.md:18` say the three S3 raw logs did not exist and
  that the text carries the fallback sentence. The set was rebuilt on 24 Nov 00:00 (`tools/build_v43.log:1`:
  "compile logs on disk: True") and the texts now cite the logs. The report (mtime 23:19) describes the earlier
  build. It also says "1 NEW-M38"; the current log has 2 (the section number 1.2, in `erratum_v43.tex` and
  `RECORDS_ERRATA_v43.md`; harmless). Update both files.

**Wording**

- **F-2.** `erratum_v43.tex:689` and `RECORDS_ERRATA_v43.md:117` say "no Lean file was built or edited". The logs
  cited two lines later (L-4) come from a per-file recompile of the three S3 files on 23 Nov 23:58 – 24 Nov 00:00
  (sources unchanged, hashes equal to `snap_s3`). Say so, and say the logs end with rc=0, as U-1 does for the other
  five.
- **F-3 (inherited from v42).** `U1_regularity_v43.tex:46`: "the compile fails if any other axiom, in particular
  `sorryAx`, is used anywhere in the import closure". `#print axioms` covers the dependencies of the named theorem,
  not every imported declaration. Write "anywhere in the dependencies of these theorems".
- **F-4 (W-1 residue).** `C1_identification_v43.tex:176`: "the theorems used are compiled (...; review
  `REVIEW_W7.md`)" is followed by (PC), (MS) with no pointer to fix23 Part A. Add it, as at `:31` and `:52`.
- **F-5 (inherited from v42).** `U1_regularity_v43.tex:75`, Lemma "the Lean velocity is CDF", is stated for an
  immersion on an open `U ⊂ ℝ`. (V) `flowVel_eq_neg_kss` is for `X : ℝ → ℂ`, C^∞ with `X' ≠ 0` on all of ℝ, and the
  chart form `flowVel_chart_eq_lift` is for lifts of immersions of `Circle`. State the lemma for these two cases;
  nothing else is used.
- **F-6.** `main_v3_v43.tex:141`: "(P1) ... is proved from compiled theorems (..., commit pending; ...)". Add "full
  build pending", or point to the box of `:81`; the proof also uses the [hand] steps of U-1 and C-1.

**Note**

- The open lists (`step3_v43_main.tex:101`, erratum S3, main_v3) still list both Windows kit runs as open. The
  `kit_step3_v4` Windows run was done on 24 Nov (reproduces; two kit bugs, kit v5 needed). The texts are
  conservative here; update with the next text round.

## Unchanged open items

Full `lake build CWFlow.Curve.All` and its log; commit of `CWFlow/Curve/`; (P3′), Lemma O arc by arc, the
single-point statement; Windows run of kit 4b and kit v5; LaTeX compile; step-4 patch (authors' go); final review.
