# N1.lean sorry-closing report (2026-11-23)

**Result:** S1, S2, S3 and S4 are closed, with no statement changes. S5 (gluing) is narrowed to a single seam-smoothness lemma `glue_seam`, which is still `sorry`. Everything was checked by a real `lake build CWFlow.N1` (Lean v4.33.1, Mathlib, `C:\hsl\CWFlow`): "Build completed successfully (10199 jobs)".

- **Installed file:** `/mnt/c/hsl/CWFlow/CWFlow/N1.lean` = `CWFlowN1_v5.lean`, sha256 `207bcb037decf663d54910e326761021e3aad94b112aef37a6a8e7bcbeb6a8ee`.
- **Version history** (each one built successfully before it was installed):

  | Version | Change | sha256 |
  |---|---|---|
  | v2 | baseline, 5 sorries | d2dc6b96… |
  | v3 | S2 + S4 closed | a75ca2a3… |
  | v4 | S1 + S3 closed as well | e81ba0ed… |
  | v5 | S5 narrowed to `glue_seam` | 207bcb03… |

- **Backups:** `N1_v2.lean.bak`, `N1_v3.lean.bak`, `N1_v4.lean.bak` sit next to N1.lean. N1_v1.lean.bak was already there.
- **Build logs:** `build_v3.log`, `build_v4.log`, `build_v5.log` (this directory).
- **Library files:** no library file was edited.
  - Three repo modules that were already built are now imported: `Scalar.Smooth`, `Flow.DiffeoChart`, `Flow.TanAlgebra`.
  - Scratch files are in `/mnt/c/hsl/CWFlow/scratch/` (T1 = S1, T2 = S2, T3 = S3, T4 = S4, T5 = S5 reduction).

## Final `#print axioms` (from build_v5.log)
```
cdf_exists                  [propext, Classical.choice, Quot.sound]
cdf_unique                  [propext, Classical.choice, Quot.sound]
cdf_unique_up_to_reparam    [propext, Classical.choice, Quot.sound]
cdf_symmetric               [propext, Classical.choice, Quot.sound]
flowVel_eq_neg_kss          [propext, Classical.choice, Quot.sound]
flowVel_conj                [propext, Classical.choice, Quot.sound]
IsFlow.conj_reflect         [propext, Classical.choice, Quot.sound]
IsFlow.shift                [propext, Classical.choice, Quot.sound]
IsFlow.glue_of_seam         [propext, Classical.choice, Quot.sound]
IsFlow.glue                 [propext, sorryAx, Classical.choice, Quot.sound]
```
- The file contains no `axiom`, `native_decide` or `admit`.
- The only `set_option` lines are `linter.unusedSectionVars false` in three sections. That option silences a style linter and does not affect sorry detection.

## Per item

**S4 `IsFlow.shift`: closed.**
- Proved generically as `IsFlow.shift_gen`, for any E, M, V and J. It uses `ContMDiffOn.comp` with (t,x) ↦ (t+a,x), and `HasDerivWithinAt.scomp`.
- It needs `0 ≤ a` (for the MapsTo). The hypothesis `haT : a ≤ T` is unused; the statement is unchanged.

**S2 `flowVel_conj`: closed.**
- It holds with **no hypotheses at all**: `flowVel_conj_any X y`. The reason is that conj is an ℝ-linear isometry equivalence, so `LinearIsometryEquiv.comp_fderiv` and `fderiv_fun_neg` apply unconditionally.
- Intermediate results:
  - g, g⁻¹, tr and perp are unchanged;
  - `meanCurv_conj`: H(conj∘X) = conj H;
  - `flowPot_conj`: flowPot(conj∘X) = −flowPot X.
- The drafted statement is kept verbatim. Its hypotheses hU, hX, hinj and hy are now unused (linter warnings only).

**S1 `flowVel_eq_neg_kss`: closed, statement unchanged** (E = ℝ, X : ℝ → ℂ, `ContDiff ℝ ∞ X`, X' ≠ 0 everywhere).
- **Basis issue avoided:** `Module.finBasis ℝ ℝ` is an unknown nonzero vector, so the proof never uses it directly:
  - `ginvOf_dim1`: g⁻¹ξ = ξ(1)/‖L 1‖²;
  - `trOf_dim1`: tr_g B = B(1,1)/‖L 1‖² for any B with B a c = (ac)·B 1 1. It is proved through the repo's `sum_dualBas_smul`, so no basis vector appears.
- **Proof chain:**
  1. `meanCurv_real`: H = (X'' − ⟪X'',X'⟫X'/|X'|²)/|X'|²;
  2. `JmeanCurv_real`: JH = −(k/|X'|)·X';
  3. `flowPot_real`: flowPot = −k'/|X'| = −k_s;
  4. `flowVel_real`.
- **Sign conventions** (as in the draft):
  - J = i·;
  - k = ⟪X'', iX'⟫/|X'|³;
  - T = X'/|X'|, N = iT;
  - result: flowVel = −k_ss·N.
  - For the counterclockwise unit circle, k = +1 and N points inward.
- **Why the statement stays on ℝ:**
  - An E1 = EuclideanSpace ℝ (Fin 1) corollary through the repo's `flowVel_comp_linearEquiv` (Scalar/LinearChange.lean) was written and attempted.
  - `LinearChange.olean` is not built (that module is not in the import closure), so the corollary was dropped rather than building an extra library module.
  - Going from ℝ to E1 charts on Circle is one application of that lemma together with `flowVel_congr_of_eventuallyEq` (locality).

**S3 `IsFlow.conj_reflect`: closed, statement unchanged.** `cdf_symmetric` is now sorry-free. It is built from two generic lemmas and one diffeomorphism:
- `IsFlow.conj` (any E and M, with V = ℂ and J = Jc):
  - the immersion clause comes from the new `isLagImmersion_conj`;
  - the velocity clause comes from the CLM chain rule plus `flowVel_conj_any`.
- `IsFlow.comp_diffeo` (any E, M, V and J, with φ a fixed diffeomorphism of M): it uses the repo's `isLagImmersion_comp_diffeo`, `chartRep_comp_diffeo`, `flowVel_comp`, `flowVel_congr_of_eventuallyEq` and `contDiffOn_flowPot`.
- `invDiffeo : Circle ≃ₘ Circle`, p ↦ p⁻¹, comes from Mathlib's `LieGroup (𝓡 1) ω Circle`.
- With these, `conj_reflect := IsFlow.conj (IsFlow.comp_diffeo h invDiffeo)`.

**S5 `IsFlow.glue`: narrowed, not closed.**
- **Proved:** `IsFlow.glue_of_seam` (generic, axiom-clean) derives the full `IsFlow` for the glued family from one hypothesis, joint smoothness at the seam points (a, p). It covers:
  - the immersion clause;
  - the velocity clause, including one-sided derivatives at t = a (`HasDerivWithinAt.union`);
  - smoothness at every point with t ≠ a.
- **Still sorry:** the remaining gap is the new lemma
  ```lean
  theorem glue_seam ... (p : Circle) :
      ContMDiffWithinAt (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E1)) 𝓘(ℝ, ℂ) ∞
        (fun q : ℝ × Circle => (if q.1 ≤ a then Φ₁ q.1 else Φ₂ (q.1 - a)) q.2)
        (Icc 0 (a + b) ×ˢ univ) (a, p) := by sorry -- SORRY (S5)
  ```
  The statement of `IsFlow.glue` is unchanged.
- **Blocker** (not attempted; estimated at several days of Lean work):
  1. Mathlib has no lemma for gluing `ContDiffOn` across a closed seam; it only has `ContDiffOn.union_of_isOpen`. One would have to build a `HasFTaylorSeriesUpToOn` on the union from matching one-sided `iteratedFDerivWithin`s.
  2. That all one-sided time and mixed derivatives match at t = a has to come from the PDE, by induction on m: ∂ₜ^m f is a function of the spatial jets. This needs one-sided t-derivatives to commute with spatial derivatives.
  3. The result then has to be moved between `ContMDiffWithinAt` and chart `ContDiffWithinAt`.
- **What this means for U1:** gluing stays a hand argument, but its only unformalized step is now the seam regularity.

## Practical notes
- `lake env lean scratch\X.lean` with narrow imports takes about 35 s. With `import CWFlow.N1` it takes about 7 min, and a full build of N1 takes about 7.5 min.
