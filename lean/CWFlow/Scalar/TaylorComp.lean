import CWFlow.Scalar.Multilinear

/-!
# The Taylor composition of Faà di Bruno's formula

Mathlib's `iteratedFDeriv_comp` writes `Dⁿ(g ∘ f)(x)` as `(ftaylorSeries g (f x)).taylorComp (ftaylorSeries f x) n`.
Here:
* `taylorComp_congr`: the `n`-th coefficient only uses `q k` for `k ≤ n` and `p k` for `1 ≤ k ≤ n`;
* `contDiffAt_taylorComp`: it is smooth in the coefficients it uses.
-/

open scoped ContDiff

noncomputable section

namespace CWFlow

universe u

variable {E F G : Type u} [NormedAddCommGroup E] [NormedSpace ℝ E] [NormedAddCommGroup F]
  [NormedSpace ℝ F] [NormedAddCommGroup G] [NormedSpace ℝ G]

theorem taylorComp_congr {q q' : FormalMultilinearSeries ℝ F G} {p p' : FormalMultilinearSeries ℝ E F}
    {n : ℕ} (hq : ∀ k ≤ n, q k = q' k) (hp : ∀ k, 0 < k → k ≤ n → p k = p' k) :
    q.taylorComp p n = q'.taylorComp p' n := by
  simp only [FormalMultilinearSeries.taylorComp]
  refine Finset.sum_congr rfl fun c _ => ?_
  simp only [FormalMultilinearSeries.compAlongOrderedFinpartition]
  rw [hq _ c.length_le]
  congr 1
  funext m
  exact hp _ (c.partSize_pos m) (c.partSize_le m)

theorem contDiffAt_taylorComp {P : Type*} [NormedAddCommGroup P] [NormedSpace ℝ P]
    {Q : P → FormalMultilinearSeries ℝ F G} {Pp : P → FormalMultilinearSeries ℝ E F} {x₀ : P} {n : ℕ}
    (hQ : ∀ k ≤ n, ContDiffAt ℝ ∞ (fun x => Q x k) x₀)
    (hP : ∀ k ≤ n, ContDiffAt ℝ ∞ (fun x => Pp x k) x₀) :
    ContDiffAt ℝ ∞ (fun x => (Q x).taylorComp (Pp x) n) x₀ := by
  simp only [FormalMultilinearSeries.taylorComp]
  refine ContDiffAt.sum fun c _ => ?_
  have h1 :=
    (ContinuousLinearMap.contDiff (𝕜 := ℝ) (E := F [×c.length]→L[ℝ] G)
      (F := ContinuousMultilinearMap ℝ (fun m : Fin c.length => E [×c.partSize m]→L[ℝ] F)
        (E [×n]→L[ℝ] G)) (n := ∞)
      (c.compAlongOrderedFinpartitionL ℝ E F G)).contDiffAt.comp x₀ (hQ _ c.length_le)
  have h2 : ContDiffAt ℝ ∞ (fun x => fun m : Fin c.length => Pp x (c.partSize m)) x₀ :=
    contDiffAt_pi.2 fun m => hP _ (c.partSize_le m)
  have h := (contDiff_cml_apply c.length (Ei := fun m => E [×c.partSize m]→L[ℝ] F)
    (G := E [×n]→L[ℝ] G)).contDiffAt.comp x₀ (h1.prodMk h2)
  refine h.congr_of_eventuallyEq (Filter.Eventually.of_forall fun x => ?_)
  rfl

end CWFlow
