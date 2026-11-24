# REVIEW_W7: independent review of the Lean bridge (CWFlow/Curve), 23 Nov 2026

Read-only review. Nothing was built, run in Lean, or edited; this file is the only output.
Reviewed: `snap_final/*.lean` (14 files, byte-identical to `/mnt/c/Users/micahw/agents/lean/CWFlow/CWFlow/Curve/`, and
matching `snap_final/SHA256SUMS`), the library definitions they rest on, `SCOPE.md`, `IMPL_REPORT.md`, `IMPL_LOG.md`,
the raw compile logs, and `U1_regularity_v41.tex` §1–2, `C1_identification_v41.tex`.

## Verdict

**Sound; no blocking finding.** T1–T5 and the n = 1 port state what the paper needs. The Lean results replace DKS
completely in U-1 and C-1 (both (D1) with the C-4 completion and (D2) with author item A-2). Four should-fix items,
none mathematical: (S1) settle the build with one `lake build`, (S2) commit and push so the paper has something to cite,
(S3) two small hand bridges that are not Lean theorems (periodic curve → `IsLagImmersion` on `Circle`; the restart
identity C-1(v)), (S4) stale `IMPL_REPORT.md`. Text changes are listed in §3.

What I could and could not check: I read statements and definitions; proofs are the compiler's job. I did not re-run the
compiler. The evidence that the files compiled is the raw logs (rc=0 for LiftV, Lift3, N1, Lift, Criterion), the `.olean`
files on disk, and the `#guard_msgs` axiom checks inside `Criterion.lean`, `N1.lean`, `Lift3.lean`, which fail the
compile if any axiom other than `propext`, `Classical.choice`, `Quot.sound` (in particular `sorryAx`) is used.

---

## 1. Statement fidelity

### Definitions (`Defs.lean`, `N1Basic.lean`, library)

| Object | Lean | Reading | OK |
|---|---|---|---|
| `E1` | `EuclideanSpace ℝ (Fin 1)` | model space of Mathlib's charts on `Circle` | yes |
| `Jc` | `ContinuousLinearMap.mul ℝ ℂ I` | multiplication by i; `Jc_isom`, `Jc_JJ` proved | yes |
| `lift F θ` | `F (Circle.exp θ)` | `Circle.exp θ = e^{iθ}`, period 2π (`Circle.periodic_exp` is used for `IsCDF.periodic`) | yes |
| `curv c θ` | `⟪c'', i c'⟫ / ‖c'‖³` | signed curvature, positive for the counterclockwise circle, with N = iT | yes |
| `dS c f θ` | `f'(θ) / ‖c'(θ)‖` | ∂_s | yes |
| `bendingEnergy c` | `∫₀^{2π} curv² ‖c'‖ dθ` | ∫κ² ds over one period, i.e. over the parametrised closed curve once | yes |
| `IsLagImmersion E1 Jc F` | chartwise C^∞, injective derivative, Lagrangian | see below | yes |
| `IsFlow E1 Jc T Φ` | C^∞ on `Icc 0 T × Circle`, immersion at each t ∈ [0,T], `HasDerivWithinAt … (Icc 0 T) t` with value `flowVel` in every chart | smooth up to both ends; one-sided at 0 and T | yes |
| `IsCDF I γ` | I open, γ C^∞ on I×ℝ, 2π-periodic, regular, `∂_t γ = −κ_ss · (i γ_θ/|γ_θ|)` | classical CDF in the θ-parametrisation, purely normal | yes |

- **Sign and normal.** `flowVel_eq_neg_kss`: `flowVel Jc X θ = −(κ_ss) · (i X'/|X'|)` with the same `curv`, `dS`. So the
  curvature vector is κN, N = iT, and the flow is ∂_tγ = −κ_ss N. This is the paper's CDF (normal velocity −κ_ss) and
  DKS's ∂_t f = −∇_s²κ; it is invariant under orientation reversal (κ and N both flip). Length decreases: L' = −∫κ_s².
- **`IsLagImmersion E1 Jc F₀` is exactly "C^∞ immersion S¹ → ℂ".** `isLagImmersion_iff : IsLagImmersion E1 Jc F ↔
  IsCurveImmersion F` (N1.lean) with the Lagrangian clause automatic (`lag_dim_one`). `IsCurveImmersion` is chartwise
  smoothness plus injective derivative in Mathlib's (stereographic) charts of `Circle`. No embeddedness, no orientation,
  no turning-number hypothesis. In terms the paper uses: `contDiff_lift_of_chart` and `deriv_lift_ne_zero_of_chart` show
  that then `lift F` is C^∞ on ℝ with nonvanishing derivative. The converse direction is not in Lean (S3a below).
- **Uniqueness class.** `cdf_unique` is equality of maps among `IsFlow`s (purely normal velocity, C^∞ on the closed
  time interval including t = 0, parametrised by `Circle`). Flows with a tangential term are covered by
  `cdf_unique_up_to_reparam` (equality up to a family of diffeomorphisms χ_t, χ_0 = id). Same class as (L2) in the texts.

### T1 `cdf_exists_maximal`  — faithful

`∃ Tmax : ℝ≥0∞, 0 < Tmax ∧ ∃ Φ, Φ 0 = F₀ ∧ (∀ T>0, ofReal T < Tmax → IsFlow T Φ) ∧ (every IsFlow Ψ on [0,T], T>0, from F₀
has ofReal T < Tmax and Ψ = Φ on [0,T]) ∧ (Tmax ≠ ⊤ → ∀ K, ∃ t ∈ [0, Tmax.toReal), K < bendingEnergy (lift (Φ t)))`.

- This is Lemma "maximal smooth solution" (existence, uniqueness, every flow is a restriction, no flow on [0,Tmax])
  together with Lemma "blow-up criterion": sup over [0,Tmax) of ∫κ² ds is +∞.
- `Tmax.toReal` occurs only under `Tmax ≠ ⊤`, so the junk value `⊤.toReal = 0` is never used.
- `Φ : ℝ → Circle → ℂ` is total; its values for t ∉ [0,Tmax) are arbitrary and never referenced.
- Not vacuous: `0 < T` is required wherever `IsFlow T` appears (for T < 0 `IsFlow` is trivially true, for T = 0 the
  velocity clause is trivial; neither case can enter).
- Weaker than DKS in one respect only: unboundedness (sup = ∞), no rate (T−t)^{-1/4} and no statement that the limit is
  ∞. Every use in the v41 texts is in the sup form (`U:lem:bu`, C-1 Props P1/P2, `continuation_v41` l.70,
  `bootstrap_v41` l.99), so nothing is lost.

### T2 `cdf_extends_of_bounded_energy` — faithful

Hypotheses: `0 < T`; `IsFlow T' Φ` for all T' ∈ (0,T); t₀ ∈ (0,T); `bendingEnergy (lift (Φ t)) ≤ K` on `[t₀, T)`.
Conclusion: `∃ ε>0, Ψ, IsFlow (T+ε) Ψ ∧ Ψ = Φ on [0,T)`.

- This is (D1) in the strong form: the extension is an `IsFlow` on the closed interval [0, T+ε], so Ψ(T) is an immersion
  and the velocity identity holds at T. The C-4 completion in C-1 (I5) is therefore covered.
- `Ico` vs `Icc`: the bound is asked only on the half-open [t₀,T) and agreement is given on [0,T); correct.
- Note: `hT` is unused (implied by `ht₀`); `t₀ > 0` is required, harmless (a bound near T is all that is ever available).

### T3 `aNormSq_eq_curv_sq` — faithful, internal

For an immersion F and a chart point y, there is θ with `Circle.exp θ = chartInv ι1 z y` and `aNormSq = curv(lift F)(θ)²`.
It is used only inside T2 to feed the library's `flow_extends_of_bounded_A`. Neither T1 nor T2 mentions `aNormSq`, so the
meaning of the library's `aNormSq` is not an input to the paper; I did not audit it.

### T4 `IsFlow.lift_isCDF` — faithful, with one limitation (note)

`IsFlow E1 Jc T Φ → IsCDF (Ioo 0 T) (fun q => Φ q.1 (Circle.exp q.2))`. Interior times only: smoothness on (0,T)×ℝ and
the two-sided time derivative for t ∈ (0,T). Smoothness up to t = 0 and t = T is in `IsFlow.smooth` itself, not in T4.
The paper needs T4 only as the Lean form of Lemma `U:lem:vel` (the lift moves by −κ_ss N); at t = 0 the identity follows
by continuity from `IsFlow.smooth`, or by citing `IsFlow.vel` with `flowVel_eq_neg_kss`.

### T5 `IsCDF.curv_sq_le` — faithful

A CDF on (a,b) with ∫κ² ds ≤ K on [t₀,b) has κ² ≤ C on [t₀,b)×ℝ (b real, so finite). Route: Ė₁ ≤ (14/3)⁶K⁷ for
E₁ = ∫κ_s² ds, and κ² ≤ K²/2 + K + E₁ (uses L·∫κ² ≥ 2). Scaling of the rate bound is consistent (L⁻⁷ both sides).

### Port in `N1.lean`

`cdf_exists`, `cdf_unique`, `cdf_unique_up_to_reparam`, `cdf_symmetric` (ι(p) = p⁻¹ = p̄ on `Circle`, σ = conj: exactly
(S4.2)'s symmetry), `IsFlow.shift_gen` (time shift, any a ≥ 0), `IsFlow.conj_reflect`, `cdf_exists_maximal_A`. All are
plain instantiations; statements match the texts' (L1), (L2) and Lemma "symmetry". The five `sorry`s of the 21 Nov draft
are gone (S1–S4 proved; S5 gluing dropped in favour of the library's `rst_glue`). The one text hit for "sorry" in
`CWFlow/` is a comment in `N1.lean` l.6.

### Vacuity / degeneracy checks

- No unproved `Prop` hypothesis is carried by T1/T2 or by `exists_maximal_flow`, `flow_extends_of_bounded_A`,
  `exists_flow_of_immersion`, `flow_unique_of_immersion` (only `hJi`, `hJJ`, `hdim`, discharged by `Jc_isom`, `Jc_JJ`,
  `Complex.finrank_real_complex`).
- `IsFlow` is satisfiable given an immersion (`cdf_exists`). **Note (N1):** there is no Lean example that
  `IsLagImmersion E1 Jc` is inhabited (e.g. the inclusion `Circle → ℂ`). Mathematically it obviously is; a two-line
  `example` would make non-vacuity machine-checked, and S3a below would do it for Γ₀-type data.

## 2. Instances and junk values

- **New global instance** `instIsManifoldCircleInfty : IsManifold 𝓘(ℝ, E1) ∞ Circle` (CircleChart.lean l.26), from
  Mathlib's `ω` structure. `IsManifold` is `Prop`-valued, so it cannot create a diamond or change the meaning of any
  term; it only lets the library's `[IsManifold 𝓘(ℝ,E) ∞ M]` variable block resolve at `Circle`. Harmless. **Note.**
- `instNeZeroFinrankE1 : NeZero (finrank ℝ E1)` (N1Basic l.30): `Prop`-valued, true, harmless.
- No `axiom`, `native_decide`, `implemented_by`, `unsafe`, `maxHeartbeats`, local notation or attribute changes in the
  Curve files (text search). Only `set_option linter.unusedSectionVars false`.
- **Junk values.** `deriv` of a non-differentiable function is 0 and a non-integrable interval integral is 0 in Mathlib.
  Every occurrence of `curv`, `dS`, `bendingEnergy` in T1, T2, T4, T5 (statements and hypotheses) is at a time where the
  curve is an immersion: in T1, t ∈ [0,Tmax) and `IsFlow` holds on some [0,T] ∋ t; in T2, t ∈ [t₀,T) ⊂ (0,T); in T5,
  t ∈ [t₀,b) ⊂ I. There `lift (Φ t)` is C^∞ with nonvanishing derivative (`contDiff_lift_of_chart`,
  `deriv_lift_ne_zero_of_chart`), so `curv` is the true curvature and the integrand of `bendingEnergy` is continuous.
  No statement relies on a junk value. `IsCDF.periodic` is asked for all t (not only t ∈ I); stronger than needed,
  satisfied by every lift.

## 3. Do the Lean results replace DKS? Text changes

**Yes, fully**, for U-1 and C-1, modulo the two hand bridges S3a, S3b. Mapping:

| Paper item | Now | Remaining hand content |
|---|---|---|
| U1 (L1), (L2), instance paragraph | `cdf_exists`, `cdf_unique`, compiled at n = 1 | none |
| U1 Lemma `U:lem:vel` (n = 1 velocity) | `flowVel_eq_neg_kss`, T4; conj: `flowVel_conj_any`; reparam: library `flowVel_comp` | none ([sym] S1–S3 become cross-checks) |
| U1 Lemma `U:lem:glue` (gluing) | not needed for `U:lem:max` any more (T1 has maximality built in). Lean has `rst_glue` (needs agreement on an overlap [0,s], s>0) and `isFlow_restart`; the endpoint-matching form stated in the text is **not** a Lean theorem | see S3b |
| U1 Lemma `U:lem:max` | T1 clauses 1–2 | S3a (data on `Circle`) |
| U1 Lemma `U:lem:sym` | `cdf_symmetric` on every [0,T], T < Tmax | none |
| U1 Lemma `U:lem:bu`, (D2), author item A-2 | T1 last clause | none; "maximal solution in the sense of DKS" disappears |
| U1 Remark `U:rem:D1`, (D1) | T2 | none |
| C-1 (I5) (D1) + C-4 completion (immersion at the extension time) | T2 (extension is an `IsFlow` on [0,T+ε]); in fact C-1 uses it only through `C1:lem:max`, i.e. T1 | none |
| C-1 `C1:lem:max` | T1 + `cdf_symmetric` with F₀ from Γ₀ | S3a |
| C-1 Thm (v) (flow from γ(t₁) = shifted flow from Γ₀, T_max = T_* − t₁) | `IsFlow.shift_gen` + T1 uniqueness clause + `rst_glue` + T1 maximality | S3b: the 10-line combination is not a Lean theorem |

Still hand and unaffected by this work: the identification of step 4's objects (P1 Galerkin in its own gauge, P2 normal
graph) with the Lean flow (C-1 Props P1, P2, Lemma (P2-a)); and that ∫κ² ds of the step-4 curve equals
`bendingEnergy (lift (Φ t))` (reparametrisation invariance, one sentence).

**Text changes needed (not applied):**

1. `U1_regularity` label line: redefine [Lean] as "compiled theorem of `CWFlow/Curve` (commit …), axioms `propext`,
   `Classical.choice`, `Quot.sound`"; drop "statement read from source".
2. `U1_regularity` §1.2: new source line (repo, commit of CWFlow **after** `Curve/` is committed, MM commit, toolchain
   v4.33.1); list T1, T2, T4, `cdf_exists`, `cdf_unique`, `cdf_symmetric`, `flowVel_eq_neg_kss`, `IsFlow.shift_gen`,
   `rst_glue`. Delete "Not checked here: that this instance elaborates … that the repository builds".
3. `U1_regularity` §1.2 instance: M is Mathlib's `Circle` (unit complex numbers) with model `EuclideanSpace ℝ (Fin 1)`,
   **not** "the unit sphere of `EuclideanSpace ℝ (Fin 2)`". Parametrisation θ ↦ e^{iθ}; ∫κ² ds is `bendingEnergy ∘ lift`.
4. `U1_regularity` Lemma `U:lem:vel`: proof becomes [Lean] (`flowVel_eq_neg_kss`, T4, `flowVel_conj_any`, `flowVel_comp`);
   keep S1–S3 as an independent check.
5. `U1_regularity` §1.3 (DKS): delete (D1), (D2) as inputs and the author-confirmation sentence; keep DKS (and Chou) as
   a literature remark only ("the criterion goes back to …; DKS also give the rate (T−t)^{-1/4}, not used").
6. `U1_regularity` §2: Lemma `U:lem:glue`: delete, or restate as `rst_glue` (overlap form) [Lean]. Lemma `U:lem:max`:
   proof = T1 [Lean]. Lemma `U:lem:sym`: proof = `cdf_symmetric` [Lean]; delete "it is not a Lean lemma". Lemma
   `U:lem:bu`: proof = T1 last clause [Lean]; statement unchanged (sup form). Remark `U:rem:D1`: replace by one line
   citing T2. Closing sentence "(P1′) with …": existence, uniqueness and criterion all [Lean].
7. `U1_regularity` gaps list (ll.180–181, 184): close the Lean and DKS items; add the two hand bridges S3a, S3b if they
   are not formalised.
8. `C1_identification` label line (l.26) and (I4) (l.47): compiled, n = 1 elaborated; delete "read from source, not
   compiled … (review C-3)"; "gluing, maximal, symmetry, blow-up" are [Lean] for any immersion F₀ of `Circle`.
9. `C1_identification` (I5) (l.49): delete (D1) and the C-4 completion; if kept, replace by T2. l.28 "the criterion from
   DKS (D1)" → "from the Lean theorem T1". `C1:lem:max` proof (l.56): T1 + `cdf_symmetric`, plus the S3a sentence for
   F₀ := Γ₀ on `Circle`. Thm (v) proof (l.162): cite `IsFlow.shift_gen`, T1, `rst_glue` (S3b). Gaps ll.171–172: close.
10. Dependent status sentences that say "Lean read, not built", "n = 1 not elaborated", "DKS (D1) / Theorem 3.1 / re-read
    DKS / author item A-2": `continuation_v41` ll.24, 70–73, 131; `bootstrap_v41` l.64 (table row "DKS (D1)");
    `step3_v41_main` l.99; `step4_handoff_v41` l.42; `ei_lemma_v4c_v41` l.159; `GS_lemma_v41` l.137; `erratum_v41`
    ll.353, 497, 538, 584, 610, 636, 654; `main_v3_v41` status boxes and (P1) paragraph; `lemmaU_v41` l.381 ("quoted
    from the literature" → proved in Lean). `P1_CONFIRMATION.md` becomes historical.
11. Bibliography: keep DKS2002 and Chou2003 as attribution; cite the Lean repository at the new commit.
12. CLAUDE.md / status lists: "Lean build n = 1" and "DKS re-read" leave the list of external items once S1, S2 are done.

## 4. The compile caveat

Facts found (all read-only):

- The 84 "modified" tracked files in the CWFlow repo differ from HEAD (8718f95) **only by line endings**:
  `git diff --ignore-cr-at-eol` is empty. `../MM` (7b15b01) is clean. So the sources the `.olean`s were built from and
  the committed sources have the same content; the only new content is the untracked `CWFlow/Curve/` and `scratch_curve/`.
- Over the whole import closure of `CWFlow.Curve.Criterion` (372 local modules in CWFlow and MM), every `.olean` is newer
  than its source and newer than every `.olean` it imports (mtime check; 0 violations). The Curve `.olean`s were written
  in dependency order (22:27 N1Basic, 22:28 LiftV, 22:29 Lift3, 22:30 N1, 22:31 Lift, 22:32 Criterion).
- The axiom guards in `Criterion.lean` would have failed on any `sorryAx` anywhere in the closure, including inside
  the pre-existing `.olean`s.
- `CWFlow/Curve/All.lean` (the aggregator) was never compiled (no `All.olean`).
- There is no raw log for `N1Basic` in `lean_bridge/` (its `.olean` exists and everything downstream compiled against it).

Risk: **low**. What remains is (i) mtimes are not content hashes, and plain `lean` does not check that an imported
`.olean` corresponds to its source; (ii) the per-file route was never confirmed by Lake. A mismatch would most likely
show as a compile error, not as a wrong theorem.

Single command to settle it (from `C:\Users\micahw\agents\lean\CWFlow`, one build at a time, per the repo's rules):

    lake build CWFlow.Curve.All

Expect it either to finish quickly (only `All` to build) or, if Lake's content hashes see the CRLF change, to rebuild
the whole closure (hours; the `IsFlow` closure needs about 10 GB of commit). Either outcome settles it; exit code 0
with no `error:` lines is the pass condition.

## Findings by severity

**Blocking:** none.

**Should-fix**

- **S1** Run `lake build CWFlow.Curve.All` once (§4). Until then the texts should say "compiled per file" rather than
  "built".
- **S2** `CWFlow/Curve/` is untracked and the repo has no remote configured; the texts cite GitHub commit `d72b744`,
  which does not contain these files. Commit, push, and cite the new commit (and MM's 7b15b01) before any text says [Lean].
- **S3a** Interface for the data: T1 takes `F₀ : Circle → ℂ` with `IsLagImmersion`. Γ₀ and γ(t₁) are given as periodic
  parametrised curves. "c smooth, 2π-periodic, c′ ≠ 0 ⇒ F with `lift F = c` is an `IsLagImmersion`" is not in Lean
  (only the other direction is). Either add this lemma (short, uses `contDiffOn_liftTau`-type facts) or state it as a
  [hand] sentence in `C1:lem:max` / `C1:lem:G0`.
- **S3b** C-1 Theorem (v) (the maximal flow from γ(t₁) is the shifted maximal flow from Γ₀, T_max = T_* − t₁) follows
  from `IsFlow.shift_gen`, T1 and `rst_glue`, but is not a Lean theorem, and the endpoint gluing lemma `U:lem:glue` as
  written in U-1 is not one either. Either add a Lean corollary or keep a [hand] proof that cites exactly these three.
- **S4** `IMPL_REPORT.md` is stale: it still says T1–T4 are not compiled and lists the pre-fix hash for `LiftV.lean`
  (e9bc…; final 540e2ce9…). `IMPL_LOG.md` (last entry) and `snap_final/SHA256SUMS` are the correct record.

**Notes**

- N1 No machine-checked inhabitant of `IsLagImmersion E1 Jc` (non-vacuity of T1's hypothesis); S3a would give one.
- N2 T4 covers interior times only; t = 0 via `IsFlow.smooth`/`IsFlow.vel` directly.
- N3 T1 gives sup = ∞ without a rate; sufficient for every use in the v41 texts.
- N4 Global `Prop` instance `IsManifold 𝓘(ℝ,E1) ∞ Circle`: harmless.
- N5 Unused hypothesis `hT` in T2; `push_neg` deprecation warning in `Criterion.lean`; cosmetic.
- N6 `Curve.IsFlow.lift_isCDF` lives in `CWFlow.Curve`, so dot notation on `h : CWFlow.IsFlow …` does not find it.
- N7 I did not audit the library below the statements used (`exists_maximal_flow`, `flow_extends_of_bounded_A`,
  `rst_glue`, `aNormSq`): their proofs are covered by the compiler and the axiom guards, and their statements, as read,
  carry no extra hypotheses.
