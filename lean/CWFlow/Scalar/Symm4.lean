import Mathlib.Analysis.Calculus.FDeriv.Symmetric
import Mathlib.Analysis.Calculus.ContDiff.Basic
import Mathlib.Analysis.Calculus.ContDiff.FTaylorSeries

open scoped ContDiff
open Filter

namespace CWFlow

section
variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
  {F : Type*} [NormedAddCommGroup F] [NormedSpace ℝ F]

/-- Derivative of a multilinear-map-valued function, evaluated at a fixed point. -/
private lemma fderiv_eval {n : ℕ} {T : E → ContinuousMultilinearMap ℝ (fun _ : Fin n => E) F}
    {y : E} (hT : DifferentiableAt ℝ T y) (p : Fin n → E) (a : E) :
    fderiv ℝ T y a p = fderiv ℝ (fun z => T z p) y a := by
  have h := ((ContinuousMultilinearMap.apply ℝ (fun _ : Fin n => E) F p).hasFDerivAt.comp y
    hT.hasFDerivAt)
  rw [show (fun z => T z p) = (ContinuousMultilinearMap.apply ℝ (fun _ : Fin n => E) F p) ∘ T
    from rfl, h.fderiv]
  rfl

private lemma iter_succ_eval {w : E → F} {y : E} (k : ℕ) (hw : ContDiffAt ℝ ((k + 1 : ℕ) : ℕ∞ω) w y)
    (m : Fin (k + 1) → E) :
    iteratedFDeriv ℝ (k + 1) w y m
      = fderiv ℝ (fun z => iteratedFDeriv ℝ k w z (Fin.tail m)) y (m 0) := by
  rw [iteratedFDeriv_succ_apply_left, fderiv_eval]
  exact hw.differentiableAt_iteratedFDeriv (by exact_mod_cast Nat.lt_succ_self k)

private lemma two_symm {w : E → F} {z : E} (hw : ContDiffAt ℝ 2 w z) (c d : E) :
    iteratedFDeriv ℝ 2 w z ![c, d] = iteratedFDeriv ℝ 2 w z ![d, c] :=
  IsSymmSndFDerivAt.iteratedFDeriv_cons (hf := hw.isSymmSndFDerivAt (by simp))

private lemma three_symm23 {w : E → F} {z : E} (hw : ContDiffAt ℝ 3 w z) (a b c : E) :
    iteratedFDeriv ℝ 3 w z ![a, b, c] = iteratedFDeriv ℝ 3 w z ![a, c, b] := by
  rw [iter_succ_eval 2 hw, iter_succ_eval 2 hw]
  have hev : (fun x => iteratedFDeriv ℝ 2 w x (Fin.tail ![a, b, c]))
      =ᶠ[nhds z] (fun x => iteratedFDeriv ℝ 2 w x (Fin.tail ![a, c, b])) := by
    filter_upwards [hw.eventually (by simp)] with x hx
    exact two_symm (hx.of_le (by exact_mod_cast (by norm_num : (2:ℕ) ≤ 3))) b c
  simp only [Matrix.cons_val_zero]
  rw [hev.fderiv_eq]

private lemma three_symm12 {w : E → F} {z : E} (hw : ContDiffAt ℝ 3 w z) (a b c : E) :
    iteratedFDeriv ℝ 3 w z ![a, b, c] = iteratedFDeriv ℝ 3 w z ![b, a, c] := by
  have hG : ContDiffAt ℝ 2 (fderiv ℝ w) z := hw.fderiv_right (m := 2) (by norm_num)
  have key : ∀ p : Fin 3 → E, iteratedFDeriv ℝ 3 w z p
      = iteratedFDeriv ℝ 2 (fun y => fderiv ℝ w y) z (Fin.init p) (p (Fin.last 2)) :=
    fun p => iteratedFDeriv_succ_apply_right p
  rw [key, key]
  have h2 := two_symm hG a b
  have e1 : Fin.init ![a, b, c] = ![a, b] := by ext i; fin_cases i <;> rfl
  have e2 : Fin.init ![b, a, c] = ![b, a] := by ext i; fin_cases i <;> rfl
  rw [e1, e2]
  have : (fun y => fderiv ℝ w y) = fderiv ℝ w := rfl
  rw [this, h2]
  rfl

private lemma four_step {w : E → F} {y : E} (hw : ContDiffAt ℝ 4 w y) (a b c d b' c' d' : E)
    (h : ∀ z, ContDiffAt ℝ 3 w z →
      iteratedFDeriv ℝ 3 w z ![b, c, d] = iteratedFDeriv ℝ 3 w z ![b', c', d']) :
    iteratedFDeriv ℝ 4 w y ![a, b, c, d] = iteratedFDeriv ℝ 4 w y ![a, b', c', d'] := by
  have e : ∀ n : Fin 4 → E, iteratedFDeriv ℝ 4 w y n
      = fderiv ℝ (fun z => iteratedFDeriv ℝ 3 w z (Fin.tail n)) y (n 0) :=
    fun n => iter_succ_eval 3 hw n
  rw [e, e]
  have t1 : Fin.tail ![a, b, c, d] = ![b, c, d] := by ext i; fin_cases i <;> rfl
  have t2 : Fin.tail ![a, b', c', d'] = ![b', c', d'] := by ext i; fin_cases i <;> rfl
  rw [t1, t2]
  have hev : (fun z => iteratedFDeriv ℝ 3 w z ![b, c, d])
      =ᶠ[nhds y] (fun z => iteratedFDeriv ℝ 3 w z ![b', c', d']) := by
    filter_upwards [hw.eventually (by simp)] with x hx
    exact h x (hx.of_le (by exact_mod_cast (by norm_num : (3:ℕ) ≤ 4)))
  simp only [Matrix.cons_val_zero]
  rw [hev.fderiv_eq]

end

/-- The fourth derivative of a smooth function is invariant under the cyclic permutation of its last three
arguments. -/
theorem iteratedFDeriv_four_perm {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    {w : E → ℝ} {y : E} (hw : ContDiffAt ℝ ∞ w y) (m : Fin 4 → E) :
    iteratedFDeriv ℝ 4 w y ![m 0, m 3, m 1, m 2] = iteratedFDeriv ℝ 4 w y m := by
  have hw4 : ContDiffAt ℝ 4 w y := hw.of_le (by first | exact WithTop.coe_le_coe.2 le_top | exact le_of_lt (by exact_mod_cast ENat.coe_lt_top 4) | (norm_cast; exact le_top))
  have h12 := four_step hw4 (m 0) (m 3) (m 1) (m 2) (m 1) (m 3) (m 2)
    (fun z hz => three_symm12 hz _ _ _)
  have h23 := four_step hw4 (m 0) (m 1) (m 3) (m 2) (m 1) (m 2) (m 3)
    (fun z hz => three_symm23 hz _ _ _)
  have hm : ![m 0, m 1, m 2, m 3] = m := by ext i; fin_cases i <;> rfl
  rw [h12, h23, hm]

end CWFlow
