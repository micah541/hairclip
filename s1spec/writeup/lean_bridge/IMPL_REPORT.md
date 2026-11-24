# IMPL_REPORT: Lean bridge for the L² blow-up criterion (SCOPE W0–W6)

Date: 2026-11-23 (evening). Repo `/mnt/c/Users/micahw/agents/lean/CWFlow`, new files only in `CWFlow/Curve/`
(nothing else in the repo touched, nothing committed, no `lake update/clean/cache`). Log: `IMPL_LOG.md`.

## Bottom line

- **Compiled, sorry-free, axiom-clean (real leanjob builds):** W0 (sympy), W2, W3, W4 (including **T5**), and the Mathlib/`Extension` half
  of W5 (circle charts read through `θ ↦ e^{iθ}`; `|A|²` in a chart = `|A|²` of the lift; joint smoothness of the lifted
  family from chart readings, `LiftC`). K0 and K2 passed. K3 (chart
  friction) was resolved without the fallback.
- **Written but NOT compiled:** `N1Basic`, `LiftV`, `Lift3` (rest of W5: T3, T4 core), `N1` (W1 port, K1), `Lift` (T4),
  `Criterion` (W6: T2, T1). **The reason is resources, not mathematics.** Every module that imports the `Weinstein`/`Scalar` chain
  (`flowVel`, `IsLagImmersion`) or `IsFlow` needs about 6.4 GB (Scalar/Weinstein closure, measured: killed at import, peak commit
  6.40 GB) or about 9.8 GB (IsFlow closure, leanjob estimate) of Windows commit. Only about 7.0 GB was free all evening
  (commit limit 24.9 GB; vmmemWSL, Zoom, Edge and others running). I did not bypass leanjob's memory guard.
- So **K1 could not be decided**. Porting risk was reduced by checking the dependencies: every repo declaration that the N1 port uses is
  byte-identical between d72b744 and the current repo. One real drift was found: the current repo has its own
  `CWFlow.IsFlow.comp_diffeo`, so the port's lemma is renamed `IsFlow.comp_diffeo_n1`.
- **To finish:** with ≥ 8 GB free commit, run stage A; with ≥ 11.5 GB, run stage B. The script is `tools/watch_build.sh`, or by hand
  (from WSL, one at a time):
  `LJ_THREADS=1 tools/lj.sh CWFlow/Curve/N1Basic.lean CWFlow/Curve/LiftV.lean CWFlow/Curve/Lift3.lean`, then
  `LJ_THREADS=1 tools/lj.sh CWFlow/Curve/N1.lean CWFlow/Curve/Lift.lean CWFlow/Curve/Criterion.lean`
  (`tools/lj.sh` wraps the repo's `../tools/leanjob.py`, i.e. `scripts/compile1.sh`). `Criterion.lean` ends with
  `#guard_msgs` axiom checks for T1–T5. Expect a few rounds of small fixes: these files were written without compiler feedback.

## Final target statements (verbatim from the files)

(T1) `Criterion.lean`, **not compiled**:
```lean
theorem cdf_exists_maximal {F₀ : Circle → ℂ} (hF₀ : IsLagImmersion E1 Jc F₀) :
    ∃ Tmax : ℝ≥0∞, 0 < Tmax ∧ ∃ Φ : ℝ → Circle → ℂ, Φ 0 = F₀ ∧
      (∀ T : ℝ, 0 < T → ENNReal.ofReal T < Tmax → IsFlow E1 Jc T Φ) ∧
      (∀ (T : ℝ) (Ψ : ℝ → Circle → ℂ), 0 < T → IsFlow E1 Jc T Ψ → Ψ 0 = F₀ →
        ENNReal.ofReal T < Tmax ∧ ∀ t ∈ Icc (0 : ℝ) T, Ψ t = Φ t) ∧
      (Tmax ≠ ⊤ → ∀ K : ℝ, ∃ t ∈ Ico (0 : ℝ) Tmax.toReal, K < bendingEnergy (lift (Φ t)))
```
(T2) `Criterion.lean`, **not compiled**:
```lean
theorem cdf_extends_of_bounded_energy {T : ℝ} (hT : 0 < T) {Φ : ℝ → Circle → ℂ}
    (hΦ : ∀ T' ∈ Ioo 0 T, IsFlow E1 Jc T' Φ) {K t₀ : ℝ} (ht₀ : t₀ ∈ Ioo 0 T)
    (hK : ∀ t ∈ Ico t₀ T, bendingEnergy (lift (Φ t)) ≤ K) :
    ∃ ε > 0, ∃ Ψ : ℝ → Circle → ℂ, IsFlow E1 Jc (T + ε) Ψ ∧ ∀ t ∈ Ico 0 T, Ψ t = Φ t
```
(T3) `Lift3.lean`, **not compiled**. Its mathematical content is `aNormSq_chart_eq_lift` + `aNormSq_lift_real` in `LiftA.lean`, which compiled:
```lean
theorem aNormSq_eq_curv_sq {F : Circle → ℂ} (hF : IsLagImmersion E1 Jc F) (z : Circle)
    {y : E1} (hy : y ∈ Extension.chartDomE ι1 z) :
    ∃ θ : ℝ, Circle.exp θ = Extension.chartInv ι1 z y ∧
      Extension.aNormSq Jc (F ∘ Extension.chartInv ι1 z) y = curv (lift F) θ ^ 2
```
(T4) `Lift.lean`, **not compiled**. The core is `isCDF_of_chartData` in `Lift3.lean`, also not compiled:
```lean
theorem IsFlow.lift_isCDF (h : IsFlow E1 Jc T Φ) :
    IsCDF (Ioo 0 T) (fun q : ℝ × ℝ => Φ q.1 (Circle.exp q.2))
```
(T5) `Energy.lean`, **compiled**:
```lean
theorem IsCDF.curv_sq_le {a b : ℝ} {γ : ℝ × ℝ → ℂ} (h : IsCDF (Ioo a b) γ) {t₀ K : ℝ}
    (ht₀ : t₀ ∈ Ioo a b) (hK : ∀ t ∈ Ico t₀ b, bendingEnergy (fun θ => γ (t, θ)) ≤ K) :
    ∃ C : ℝ, ∀ t ∈ Ico t₀ b, ∀ θ, curv (fun θ' => γ (t, θ')) θ ^ 2 ≤ C
```
`IsCDF` (`Arclength.lean`, compiled):
```lean
structure IsCDF (I : Set ℝ) (γ : ℝ × ℝ → ℂ) : Prop where
  isOpen : IsOpen I
  smooth : ContDiffOn ℝ ∞ γ (I ×ˢ univ)
  periodic : ∀ t θ, γ (t, θ + 2 * π) = γ (t, θ)
  reg : ∀ t ∈ I, ∀ θ, deriv (fun θ' => γ (t, θ')) θ ≠ 0
  vel : ∀ t ∈ I, ∀ θ, HasDerivAt (fun s => γ (s, θ))
    (-((dS (fun θ' => γ (t, θ')) (dS (fun θ' => γ (t, θ')) (curv (fun θ' => γ (t, θ')))) θ : ℝ) : ℂ) *
      (Complex.I * (deriv (fun θ' => γ (t, θ')) θ / (‖deriv (fun θ' => γ (t, θ')) θ‖ : ℂ)))) t
```
`curv`, `dS`, `lift`, `bendingEnergy`, `E1` (`Defs.lean`) are as in SCOPE §1.3.

## `#print axioms` from a real build

Source: `axioms_light.log`, a leanjob run of `scratch_curve/AxLight.lean`, which imports `Energy` and `LiftA`. Every result below printed
`[propext, Classical.choice, Quot.sound]`:

- **T5 and the energy estimate:** `IsCDF.curv_sq_le` (T5), `IsCDF.hasDerivAt_energy1`, `integral_rate` (identity (I)),
  `integral_rate_le`, `curv_sq_le`, `two_le_length_mul_bendingEnergy` (P5), `rate_le` (C);
- **evolution equations:** `pt_cKap` (E4), `pt_cDs` (E2), `pt_cSpd` (E1);
- **circle charts:** `aNormSq_chart_eq_lift`, `aNormSq_lift_real`, `contDiff_lift_of_chart`, `deriv_lift_ne_zero_of_chart`,
  `exists_equiv_fderiv_liftTau`.

There are no `#print axioms` results yet for T1–T4, because those files have not been built.

## Files (all in `CWFlow/Curve/`)

| File | Status | sha256 |
|---|---|---|
| Defs.lean | compiled | e759837587536f9948312f22c75be8843cf5c764cb8067cafba27c748c7c6d9a |
| Periodic.lean (W2) | compiled | 220836752dca634bc4337e2b2696a85ce03d7ba35551c6448867c94e13e3411a |
| Arclength.lean (W3) | compiled | a2f6172228708d1acf9a4ba5cec260dc87f50c145160d97ee1aa8e206531bc8d |
| Energy.lean (W4, T5) | compiled | dbd9441fb455141d51e350d25a847029980d024bba1d70b1e7390322707b9c63 |
| CircleChart.lean (W5) | compiled | 831f495dd37d9bb6549be5b9cfeb8e14a9336ff539893f9ea02ba113d29fb3a7 |
| LiftA.lean (W5) | compiled | c7c5c75e3a9a8f01d10d89a675d09b22ffa148e8a0e2386b873cdbd128fb8b2d |
| LiftC.lean (W5, T4 smoothness part) | compiled | cc397d2ea4d8176acb9cdf33a3b6eae4f0398ea33b4bbc963e0a7ac36a9a9694 |
| N1Basic.lean (W1 part) | not compiled (memory) | 357cf86acc5b4d481a190e76f57ef5a2dd1b9680e738dc9c4fc22023ecc481c1 |
| LiftV.lean (W5) | not compiled (memory) | e9bc336ec8ece97cb6d0c4e2ffb254a019f378f513ca080580112c560099535f |
| Lift3.lean (W5, T3) | not compiled (memory) | 95d693ca9c29606886633437edaa271457cb10f84c7e261595f9034ec1ab8171 |
| N1.lean (W1) | not compiled (memory) | f1d3d641263df03ca09f4301cdf71e986d10f5a5f321102365cd225777ed8320 |
| Lift.lean (W5, T4) | not compiled (memory) | 0833e4dfa0c35f1d4df3e833b2be68b83bd3ce9592710ab183d38624dc6bed77 |
| Criterion.lean (W6, T1, T2) | not compiled (memory) | fe5677083c1d1fc3f5454bc90bfb91568371c38394e88f129fb9860605d952f1 |
| All.lean | aggregator, not compiled | f5c734e53bcbc3d5e839aecfed968d4733f31a6268a71e95ee82b802879304e4 |

- **Snapshots:** `snap_W2/`, `snap_W3/`, `snap_W4/` (compiled states) and `snap_W5W6_draft/` (all files, draft state).
- **W0:** `w0/check_E_I.py` and its output `w0/check_E_I.out`. All checks PASS: (E1)–(E4) for a general normal speed V, the identity (I), and
  an explicit antiderivative for the integration by parts.

## Differences from SCOPE

1. **`IsCDF.vel`:** stated in the exact form of `flowVel_eq_neg_kss`, `−(κ_ss : ℂ) * (I * (c'/|c'|))`, rather than SCOPE's schematic
   `−(κ_ss) • N`. The meaning is the same; SCOPE already says "as in `flowVel_eq_neg_kss`".
2. **N1 port split in two:** `N1Basic.lean` (Jc, the `n = 1` velocity identity, conj algebra) and `N1.lean` (the IsFlow parts and
   `cdf_exists_maximal_A`). Also:
   - `curv`, `dS`, `E1` are moved to `Defs.lean`;
   - `ι1` is moved to `CircleChart.lean`;
   - `IsFlow.comp_diffeo` is renamed `IsFlow.comp_diffeo_n1`;
   - the glue section (and the `sorry`) is dropped, as SCOPE says.

   No statement was changed.
3. **T4 core as a separate theorem:** `isCDF_of_chartData` takes the IsFlow fields in chart form. This lets more of the code compile
   without the IsFlow import closure.
4. **Names live in `CWFlow.Curve`:** `IsFlow.lift_isCDF` is `CWFlow.Curve.IsFlow.lift_isCDF`. Call it explicitly; dot notation on
   `h : CWFlow.IsFlow …` will not find it.
5. **Route as SCOPE §2.2, with simpler constants:**
   - rate bound `∫ rate ≤ (14/3)⁶ K⁷`;
   - sup bound `κ² ≤ K²/2 + K + E₁`, with `A₀a₁ ≤ (A₀² + a₁²)/2`, so no square root in the final constant;
   - (P3)/(P4) are one lemma, `sq_le_sq_add_of_periodic`.

## Remaining gaps

- **Compilation of the 7 unbuilt files.** This needs ≥ 8 GB, then ≥ 11.5 GB, of free Windows commit. Until then, T1–T4 are drafts.
- **Global instance in `CircleChart.lean`:** `instIsManifoldCircleInfty : IsManifold 𝓘(ℝ, E1) ∞ Circle`, from Mathlib's `ω`. It is harmless
  (Prop-valued), but it is new and global.
- **Not done (not mine):** W7 review; the text patches (U1_regularity, C-1).

---
## UPDATE Nov 23 night (supersedes the "not compiled" status above)
All 14 files of `CWFlow/Curve/` compile (per-file `lean`, tools/raw.sh), no sorry; Criterion's #guard_msgs axiom checks pass for T1-T5 ([propext, Classical.choice, Quot.sound]). LiftV changed (one proof step); current hashes in `snap_final/SHA256SUMS`. `All.lean` not yet compiled. Review: `REVIEW_W7.md` (sound, no blocking).
