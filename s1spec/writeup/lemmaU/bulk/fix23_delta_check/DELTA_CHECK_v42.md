# DELTA_CHECK_v42: the S3 Lean theorems and the text set v42 (23 Nov 2026)

Read-only review. Nothing was built, edited or run in Lean; this file is the only output.

## Verdict

**OK with fixes; no blocking finding.**

- **Part A (Lean, `snap_s3/`).** The four new theorems state what C-1 Theorem (v), the U-1 gluing lemma and the
  Γ₀-feeding step need. No vacuity, no ENNReal or interval-endpoint problem.
- **Part B (text v42).** The Lean citations match the files, the 12 items of REVIEW_W7 §3 are applied, DKS is
  attribution only, and the rewritten arguments (U-1 §1.2–§2; C-1 (I4), (I5), `C1:lem:max`, Thm (v)) are correct.
  Four should-fix items and nine wording items, all text or records.
- Step 3 is described as **not closed** in every status box and open list. That stays correct.

What I checked and what I did not: I read statements and definitions and compared them with the texts. Proofs are the
compiler's job, and I did not run it. The compile evidence for the S3 files is the `IMPL_LOG.md` entry and the
`.olean` files on disk; no raw log was kept (S-1).

---

## Part A. `snap_s3/` (Periodic2Circle.lean, MaxShift.lean, All.lean)

### A.0 Records

| Check | Result |
|---|---|
| `snap_s3/*.lean` vs `snap_s3/SHA256SUMS` | match (f7d02857…, 80b86094…, 91d5788e…) |
| `snap_s3/*.lean` vs repository `CWFlow/Curve/` | byte-identical |
| `snap_final/` (13 files other than `All.lean`) vs repository | byte-identical; `All.lean` differs only by the two new imports |
| `.olean` newer than source | Periodic2Circle 22:49:34 > 22:49:04; MaxShift 22:58:07 > 22:57:28; All 22:58:49 > 22:58:18 |
| `sorry`, `axiom`, `native_decide` in the two new files | none (text search) |
| `#guard_msgs` axiom checks | present for all four theorems, expected `[propext, Classical.choice, Quot.sound]` |
| Raw compile logs for the three S3 files | **none kept** (`scratch_curve/last_out.txt` is overwritten per run; it is empty, which fits a clean `All.lean` compile) |
| Repository state | HEAD 8718f95, `CWFlow/Curve/` untracked, no remote; MM at 7b15b01; toolchain v4.33.1 |

### A.1 `exists_isLagImmersion_of_periodic`, `lift_injective` — faithful

Statement: `ContDiff ℝ ∞ c`, `Periodic c (2π)`, `∀ θ, deriv c θ ≠ 0` (c : ℝ → ℂ) give
`∃ F : Circle → ℂ, lift F = c ∧ IsLagImmersion E1 Jc F`. `lift_injective`: `lift` is injective.

- The conclusion has exactly the type of the hypothesis of T1 (`cdf_exists_maximal`).
- `lift F = c` is equality of functions, so `F(e^{iθ}) = c(θ)` for all θ.
- The hypotheses are what C-1 Lemma `C1:lem:G0` supplies for Γ₀: X⁰ real-analytic, 2π-periodic in the normalised
  arclength parameter, |X⁰_α| = L₀/2π > 0. `deriv` is the real derivative of a map ℝ → ℂ, as wanted.
- No embeddedness or orientation hypothesis. Not vacuous: the hypotheses are satisfiable (any regular closed curve).
- This also answers note N1 of REVIEW_W7 in substance (T1's hypothesis is inhabited as soon as one such c exists);
  there is still no explicit Lean `example`, which is harmless.

### A.2 `cdf_maximal_shift` — faithful

Hypotheses: `h1`, `h2` are clauses (i), (ii) of T1 with F₀ replaced by `Φ 0` (T1 gives `Φ 0 = F₀`, so they apply after
one rewrite); `0 ≤ t₁`; `ofReal t₁ < Tmax`. Conclusions: (1) `0 < Tmax − ofReal t₁`; (2) `Φ(t₁+·)` is an `IsFlow` on
`[0,T]` for `T > 0`, `ofReal T < Tmax − ofReal t₁`; (3) every `IsFlow Ψ` on `[0,T]`, `T > 0`, with `Ψ 0 = Φ t₁` has
`ofReal T < Tmax − ofReal t₁` and `Ψ = Φ(t₁+·)` on `Icc 0 T`; (4) if `Tmax − ofReal t₁ ≠ ⊤`, then for every K there is
`t ∈ Ico 0 (Tmax − ofReal t₁).toReal` with `K < bendingEnergy (lift (Φ (t₁ + t)))`.

- **ENNReal.** Truncated subtraction cannot bite, because `ofReal t₁ < Tmax`. `ofReal` of a negative number is 0, and
  the hypothesis `0 ≤ t₁` excludes that case. `Tmax = ⊤` gives `Tmax − ofReal t₁ = ⊤`: clauses (2), (3) hold for every
  T and (4) is void, as it should be. `toReal` appears only under `≠ ⊤`, so the junk value `⊤.toReal = 0` is not used.
  Under `≠ ⊤` the bound is `Tmax.toReal − t₁`, so (4) speaks of absolute times in `[t₁, Tmax)`.
- **Icc / Ico.** Agreement on the closed `[0,T]`; blow-up on the half-open `[0, T_max − t₁)`. Both correct. `0 < T` is
  required wherever `IsFlow T` appears, as in T1. `t₁ = 0` is allowed and gives T1 back.
- **Convention.** The new maximal time is a *length* (`Tmax − t₁`), the convention fixed in text v40 (F-4) for
  U-1's `T_max`.
- **Not vacuous.** The hypotheses are exactly T1's output; nothing unproved is carried. Clause (4) is not assumed for
  Φ: the proof gets it from T1 applied to `Φ t₁`.
- **What it does not state** (note, not a defect): the equation `Tm' = Tmax − ofReal t₁` for the T1-flow from `Φ t₁` is
  internal to the proof (`heqT`). C-1 Thm (v) recovers it by a two-line hand step from clauses (2), (3) and T1(ii)
  (checked in B.4). A one-line Lean corollary would remove that hand step.

### A.3 `isFlow_glue_endpoint` — faithful

`0 < a`, `0 < b`, `IsFlow a Φ`, `IsFlow b Θ`, `Θ 0 = Φ a` give
`∃ Ψ, IsFlow (a+b) Ψ ∧ (∀ t ∈ Icc 0 a, Ψ t = Φ t) ∧ ∀ t ∈ Icc 0 b, Ψ (a+t) = Θ t`.
This is the endpoint form of U-1 Lemma `U:lem:glue` (same hypotheses a, b > 0, seam condition only). The two agreement
clauses are consistent at t = a. See W-4 for the one hand step the text adds.

### A.4 Library statements used

`rst_glue`, `isFlow_restart` (`Flow/Restart.lean`), `IsFlow.shift_gen`, `cdf_unique`, `cdf_exists_maximal` were read:
statements as quoted in U-1 §1.2 ((Gl), (Sh), (L2), (T1)). `rst_glue` is used in the proof of T2 (Criterion.lean l.69),
as the text says. `isFlow_restart` carries only `hJi`, `hJJ`, `hdim`, discharged by `Jc_isom`, `Jc_JJ`, `hdim1`.

---

## Part B. Text set v42 against v41

Method: per-file diff of `fix22_texts_v42/` against `fix21_texts_v41/`; `tools/edits_v42.txt` (49 edits, 11 files),
`tools/erratum_part7_v42.tex`, `tools/records_errata_add_v42.md`; a sentence-level search of all v42 files for DKS,
(D1), (D2), A-2, d72b744, "read from source", "not elaborated", "re-read". `SHA256SUMS` of the set verifies (24 files).

### B.1 Lean citations vs the Lean files — match

- Every cited name exists where the text says: `cdf_exists`, `cdf_unique`, `cdf_unique_up_to_reparam`,
  `cdf_symmetric`, `IsFlow.shift_gen`, `IsFlow.conj_reflect`, `isLagImmersion_iff`, `IsCurveImmersion` (N1.lean);
  `flowVel_eq_neg_kss`, `flowVel_conj_any`, `Jc_isom`, `Jc_JJ` (N1Basic.lean); `cdf_exists_maximal`,
  `cdf_extends_of_bounded_energy` (Criterion.lean); `aNormSq_eq_curv_sq` (Lift3.lean); `IsFlow.lift_isCDF` (Lift.lean);
  `IsCDF.curv_sq_le` (Energy.lean); `flowVel_chart_eq_lift` (LiftV.lean); `contDiff_lift_of_chart`,
  `deriv_lift_ne_zero_of_chart` (CircleChart.lean); `flowVel_comp` (library, Scalar/Reparam.lean); `rst_glue`
  (library, namespace `CWFlow`); the four S3 theorems.
- Paraphrases (L1), (L2), (T1)–(T5), (V), (Sy), (Sh), (Gl), (PC), (MS), (GE) of U-1 §1.2 agree with the Lean
  statements, including `Icc`/`Ico`, `0 < T`, and "T_max real value used only when finite". One omission: W-3.
- Definitions (`curv`, `dS`, `lift`, `bendingEnergy`, `Circle` modelled on `EuclideanSpace ℝ (Fin 1)`) agree with
  Defs.lean; v41's "unit sphere of `EuclideanSpace R (Fin 2)`" is corrected.
- Hashes: all ten 8-hex prefixes in U-1 §1.2 match `snap_final/SHA256SUMS` and `snap_s3/SHA256SUMS`. HEAD 8718f95,
  MM 7b15b01, toolchain v4.33.1 and "16 files" match the repository.

### B.2 The 12 items of REVIEW_W7 §3 — applied

| Item | Status |
|---|---|
| 1 U-1 label line | applied |
| 2 U-1 §1.2 source line; "Not checked here" deleted | applied |
| 3 `Circle`, not the sphere in ℝ² | applied |
| 4 `U:lem:vel` proof [Lean], S1–S3 kept as cross-check | applied |
| 5 §1.3 DKS as literature only | applied |
| 6 §2: glue, max, sym, bu, Remark D1 | applied; glue is in the endpoint form through (GE), which is fine (see S-3 for the erratum sentence that says otherwise) |
| 7 U-1 gaps list | applied, with one stale clause (S-2) |
| 8 C-1 label line and (I4) | applied |
| 9 C-1 (I5), l.28, `C1:lem:max`, Thm (v), gaps | applied |
| 10 dependent status sentences | applied in continuation, bootstrap, step3_main, step4_handoff, ei_lemma, GS_lemma, lemmaU, main_v3. In `erratum_v42.tex` the old lines are kept verbatim and Part 7 (section R) says to read them subject to it; acceptable for an erratum. Two residues: W-5, W-6 |
| 11 bibliography | DKS, Chou kept as attribution; `CWFlow` entry "commit pending" (W-7) |
| 12 CLAUDE.md / status lists | outside the text set; not done here |

### B.3 DKS as input; over-claims

- DKS is not an input anywhere in the v42 text. Remaining mentions are attribution, marked history ("text v41: …"),
  bibliography, or dated older status boxes that the v42 box declares superseded. `\tD` is still defined in U-1 but
  never used.
- The caveats "compiled per file", "commit pending", "full `lake build CWFlow.Curve.All` pending" appear in every file
  that makes a substantive [Lean] claim: U-1, C-1, continuation, step3_main, step4_handoff, main_v3, lemmaU, erratum
  R2. `bootstrap_v42` and `ei_lemma_v4c_v42` mention the Lean theorems in one clause each; ei_lemma names the open
  build and commit, bootstrap's table row does not (W-8).
- One statement claims more than Lean proves: S-4.

### B.4 Mathematical check of the rewritten arguments — correct

- **U-1 `U:lem:max`.** F = γ(t₁) is the value at t₁ of an `IsFlow` on [0,T] ∋ t₁ (C-1 (i), (iv)), so it is an
  `IsLagImmersion` by the `imm` clause; T1 applies. T1(ii) with T = T_max gives "no flow on [0,T_max]" and uniqueness
  on [0,T_max). No circularity: C-1 does not use U-1's lemmas any more.
- **U-1 `U:lem:sym`.** (S4.2) says conj F(p) = F(p⁻¹); with p → p⁻¹ this is the hypothesis of `cdf_symmetric`. Applied
  on each [0,T], T < T_max; the union is [0,T_max). Correct.
- **U-1 `U:lem:bu`.** T1(iii), in the sup form. Correct.
- **U-1 `U:lem:glue`.** (GE) plus the remark that the `IsFlow` clauses on [0,a+b] see only values on [0,a+b]. That
  remark is true (`ContMDiffOn` on `Icc 0 T ×ˢ univ`, `imm` on `Icc`, `HasDerivWithinAt … (Icc 0 T)`), and it is a hand
  step (W-4). The lemma is no longer used by any other lemma of the set.
- **C-1 `C1:lem:max`.** (PC) gives F with `lift F = X⁰`; F₀ is defined by F₀(e^{iα}) = X⁰(α), i.e. `lift F₀ = X⁰`;
  `lift_injective` gives F = F₀. Then T1 and (Sy). The symmetry hypothesis is `C1:lem:G0` (X⁰(−α) = conj X⁰(α)).
  Correct.
- **C-1 (I5).** Only T1(iii) is used (Props P1, P2 argue by contradiction with a bound on ∫k² ds on [0,T_*), with
  T_* < ∞). Correct. The C-4 completion is indeed not needed.
- **C-1 Thm (v).** (MS) needs 0 ≤ t₁ < T_*; the theorem has T_* > t(27.15) ≥ t₁. With Φ, T_max from T1 for
  F = γ(t₁): for T < T_max, Φ|[0,T] is a flow from γ(t₁), so (MS)(3) gives T < T_* − t₁ and Φ = γ(t₁+·); for
  T < T_* − t₁, γ(t₁+·) is a flow from γ(t₁) by (MS)(2), so T1(ii) gives T < T_max. Hence T_max = T_* − t₁ and the
  flows agree. Correct.

### B.5 Open-item lists

Consistent across U-1 gaps, C-1 gaps, continuation, step3_main "What remains", main_v3 box and §Open, lemmaU,
erratum R6: full `lake build`; commit and citation; (P3′), Lemma O arc by arc, single-point statement; Windows runs of
the v4 kits; LaTeX compile; step-4 patch (authors' go); final review. Differences: S-2 and W-1.

---

## Findings by severity

**Blocking:** none.

**Should-fix**

- **S-1 (records).** No raw compile log exists for `Periodic2Circle.lean`, `MaxShift.lean`, `All.lean`; U-1 §1.2 lists
  the logs kept and the missing `N1Basic` log, but not this. Add the sentence, and keep the log of the full
  `lake build CWFlow.Curve.All` (records item L-2 already asks for it). Until that build, the status "compiled per
  file" is right.
- **S-2 (contradiction inside U-1).** `U1_regularity_v42.tex` l.196 (gaps, item Lean (a)) says "the aggregator
  `All.lean` was never compiled". §1.2 of the same file and `IMPL_LOG.md` say it was compiled after the S3 additions.
  Delete the clause. The same item lists the compiled theorems without (PC), (MS), (GE) in its first sentence.
- **S-3 (erratum contradicts the text).** `erratum_v42.tex` R3(b) (l.679) says text v42 states the gluing lemma "as
  `rst_glue` (gluing with an overlap)". The text states it in the endpoint form, proved by `isFlow_glue_endpoint`
  (R4 says so). R3(b) is a leftover of the first draft; rewrite it.
- **S-4 (more than Lean proves).** `lemmaU_v42.tex` l.384 states (P1) as "if the maximal time is finite then
  ∫κ² ds → ∞". T1 gives sup = ∞, not a limit. No argument uses the limit (every use is a contradiction with a uniform
  bound). Write sup_{t} ∫κ² ds = ∞, as continuation (P1′) and U-1 do.

**Wording**

- **W-1** Review status of the S3 theorems. Erratum R4/R6 and TEXTS_REPORT say they are unreviewed; the status boxes of
  step3_main, main_v3, continuation l.134 and C-1 (I4) put "review W7: sound" next to lists that include (PC), (MS),
  (GE), which W7 did not see. Cite this file for the S3 theorems and drop "an independent check of the S3 theorems"
  from R6.
- **W-2** U-1 §1.2: "each theorem … is followed in its file by an axiom check". The checks for (T4) and (T5) are in
  `Criterion.lean`, not in `Lift.lean` / `Energy.lean`. Say "in its file or in `Criterion.lean`".
- **W-3** U-1 (MS): the last clause omits the range. Write "unbounded for t ∈ [0, T_max − t₁)".
- **W-4** U-1 `U:lem:glue` proof is labelled [Lean], but the step from the existential Ψ of (GE) to the piecewise Φ of
  the statement is a hand congruence. Either state the lemma with "there is an `IsFlow` Ψ with …", or mark that step
  [hand]. Likewise the last two lines of C-1 Thm (v) (T_max = T_* − t₁) are [hand] on top of (MS).
- **W-5** `main_v3_v42.tex` l.347 still lists "the literature facts (P1)" among the step-3 caveats (unchanged from
  v41); ll.432, 435 say the opposite.
- **W-6** `main_v3_v42.tex` l.81: the 13 November status box still reads "(P1′) … Dziuk–Kuwert–Schätzle Theorem 3.1
  and/or Chou; [to be confirmed by the authors]". It is dated and superseded by the v42 box, but a reader of the
  compiled page sees both; add "(superseded, see the v42 box)" in that box.
- **W-7** The `CWFlow` bibliography entries (continuation, step3_main) and the labels in C-1 and `b12_preamble` name
  only `snap_final/SHA256SUMS`. `All.lean` there is superseded and the S3 files are in `snap_s3/SHA256SUMS`; name both.
- **W-8** `bootstrap_v42.tex` (order table, row C-1): add "compiled per file; build and commit pending", as elsewhere.
- **W-9** U-1 (PC): "it also shows that the hypothesis of (T1) is inhabited" needs one regular periodic curve
  (e.g. θ ↦ e^{iθ}); say so, or drop the remark. `\tD` in the U-1 preamble can be removed.

**Notes (no action in the texts)**

- `lean_bridge/IMPL_REPORT.md` is still stale (W7 S4; records item L-3 notes it).
- `ERRATUM_v42.md` is a pointer to Part 7 plus the older text verbatim; its old Lean/DKS items (ll.777–778, 846,
  867–868) are history and are read subject to section R.
- Item 12 of W7 §3 (CLAUDE.md status) is not part of the text set.
- External items unchanged: full `lake build CWFlow.Curve.All`; commit of `CWFlow/Curve/`; (P3′), Lemma O, the
  single-point statement; Windows kit runs; LaTeX compile; step-4 patch; final review.
