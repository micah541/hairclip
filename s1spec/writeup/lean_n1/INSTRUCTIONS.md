# CWFlowN1.lean: how to build it

The file is a draft that has never been elaborated (there is no Lean toolchain in this WSL). It targets HSLflow-Whitney commit `d72b7447`, package `CWFlow/` (toolchain `leanprover/lean4:v4.33.1`, Mathlib `v4.33.1`, ShortTime `7b15b01b`, DifferentialGeometry `v0.1.3`).

## Where to put it
Use your own checkout of HSLflow-Whitney at d72b744. The copy in `src_d72b744/` is untrusted download data; do not build it.

Option A (recommended): make it a module of the CWFlow library.
```
cp CWFlowN1.lean  <repo>/CWFlow/CWFlow/N1.lean      # module name CWFlow.N1
cd <repo>/CWFlow
lake exe cache get                                   # compiled Mathlib
lake build CWFlow.N1
```
This needs no change to `lakefile.toml` or `CWFlow.lean`. Lake treats any module under the lib root `CWFlow` as local, and `lake build CWFlow.N1` builds it and its import cone: the whole of CWFlow plus ShortTime and DifferentialGeometry, compiled from source on the first build. That takes several GB of RAM and a while, as the repo README warns.

Option B, if `lake build CWFlow.N1` reports "unknown target":
```
cd <repo>/CWFlow && lake build                       # builds the library
lake env lean /path/to/CWFlowN1.lean
```

If you already have a built CWFlow (for example from building WhitneyFlow), option B needs only the second line.

## What success looks like
- The build ends with no `error:`.
- There are 5 warnings `declaration uses 'sorry'`, for `flowVel_eq_neg_kss`, `flowVel_conj`, `IsFlow.conj_reflect`, `IsFlow.shift` and `IsFlow.glue`, and no others.
- The four `#print axioms` lines print:
  - `cdf_exists`, `cdf_unique`, `cdf_unique_up_to_reparam`: `[propext, Classical.choice, Quot.sound]`.
  - `cdf_symmetric`: the same plus `sorryAx`, which comes through S3.

  If the first three show `sorryAx`, something is wrong.
- Optional: once the file compiles, replace the first three `#print axioms` with `#guard_msgs in`, as in Main.lean l.41-42.

## Spots most likely to fail, and why

1. **Mathlib lemma names (all marked `[GUESS]`).** Mathlib v4.33.1 is newer than the drafting model's knowledge, so names may have moved.
   - `real_inner_eq_norm_add_mul_self_sub_norm_mul_self_sub_norm_mul_self_div_two` in `Jc_isom`. A fallback is in a comment, or use `Complex.inner`/`re (conj x * y)` with `Complex.conj_I`.
   - `finrank_eq_one_iff_of_nonzero'` in `lag_dim_one`. An alternative: write `u = (u 0) • EuclideanSpace.single 0 1`.
   - `Complex.finrank_real_complex` and `finrank_euclideanSpace_fin`, reached through `simp [E1, …]`.
   - `ContinuousLinearMap.mul` / `mul_apply'` in `Jc_apply`. Any CLM `v ↦ I * v` works, for example `(LinearMap.mulLeft ℝ Complex.I).toContinuousLinearMap`.
   - Module path `Mathlib.Geometry.Manifold.Instances.Sphere`. If that import fails, CWFlow may already import it transitively; delete the line.
2. **Instances on `Circle`** (the probe `example`s near the top).
   - Mathlib may give `IsManifold (𝓡 1) ω Circle` (analytic) instead of `∞`. If so, add `instance : IsManifold 𝓘(ℝ, E1) ∞ Circle := IsManifold.of_le (n := ω) le_top` (name `[GUESS]`).
   - `CompactSpace Circle` might need `inferInstanceAs (CompactSpace (Metric.sphere (0:ℂ) 1))`.
   - Alternative M: `Metric.sphere (0 : EuclideanSpace ℝ (Fin 2)) 1` (the choice in U1) with `[Fact (finrank ℝ (EuclideanSpace ℝ (Fin 2)) = 1 + 1)]`. Then ι is `fun p => ⟨reflect p, _⟩`, which is more work than `p⁻¹` on `Circle`.
3. **`NeZero (Module.finrank ℝ E1)`.** `simp [E1]` should reduce the goal to `1 ≠ 0`. If the instance is not found at the call sites, give it explicitly: `(E := E1)` is already there in `cdf_exists`.
4. **Universe.** The main theorems put `E` and `V` in the same universe `u`. `EuclideanSpace ℝ (Fin 1)` and `ℂ` are both in `Type`, so this should be fine. If it is not, the error mentions `Type u`.
5. **The `ContDiff` `∞` notation.** The file opens `ContDiff` and `Manifold` as the repo does. If `∞` elaborates as `ω` somewhere, the probe `IsManifold … ∞ Circle` will show it.
6. **`cdf_unique_up_to_reparam`.** It depends on the exact elaboration of `dtW`/`flowRep` and on the diffeomorphism notation `≃ₘ⟮𝓘(ℝ, E1), 𝓘(ℝ, E1)⟯`; both are copied from NormalFlow.lean l.282-291. The step `show G 0 (χ 0 x) = Φ 0 x` relies on beta reduction plus the coercion of `χ 0` to a function. If it fails, try `simp only [hχ0]; exact congrFun h0 x`.
7. **`lag_dim_one`, `inner_Jc_self`.** The `rw` steps assume which occurrence is rewritten; the comments say which. `nlinarith`/`linarith` closes the variants.
8. **`IsFlow.glue`.** The `if t ≤ a` needs `Decidable (t ≤ a)` on ℝ. That is classical; the file is `noncomputable section`, so it should be fine.

## The sorries (see LEAN_N1_REPORT.md for what each one is worth)
- **S1** `flowVel_eq_neg_kss`: at n = 1, flowVel = −k_ss N.
- **S2** `flowVel_conj`: conj-equivariance of flowVel.
- **S3** `IsFlow.conj_reflect`: closure of IsFlow under conj∘·∘ι.
- **S4** `IsFlow.shift`.
- **S5** `IsFlow.glue`.

The main targets (i), (ii) and (ii′) contain no sorry. (iv) `cdf_symmetric` is proved modulo S3.
