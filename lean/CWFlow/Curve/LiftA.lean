import CWFlow.Curve.CircleChart
import CWFlow.Extension.ImmReparam

/-!
# `|A|²` of a closed curve read in a chart (SCOPE W5, part (T3))

At `n = 1` the squared second fundamental form `aNormSq J X y` of the repository (`Extension/ImmCoord.lean`)
is `g⁻³ ⟪X'', J X'⟫²` with `g = |X'|²`. For the lift `w ↦ c (w₀)` of a curve `c : ℝ → ℂ` this is
`⟪c'', J c'⟫² / |c'|⁶` (`aNormSq_lift_real`); by the repository's reparametrisation invariance
`aNormSq_comp` the value in any chart of the circle equals the value of the lift (`aNormSq_chart_eq_lift`).
-/

open Set Function Filter
open scoped Manifold Topology ContDiff RealInnerProductSpace

noncomputable section

namespace CWFlow.Curve

open CWFlow.Extension

theorem eb0_coord : Aθ.symm (eb (0 : Fin 1)) = 1 := by
  simp [eb]

/-- **`|A|²` of the lift**, in the coordinate of `E1`. -/
theorem aNormSq_lift_real (J : ℂ →L[ℝ] ℂ) {c : ℝ → ℂ} (hc : ContDiff ℝ ∞ c) (θ : ℝ) :
    aNormSq J (fun w : E1 => c (Aθ.symm w)) (Aθ θ) =
      ((‖deriv c θ‖ ^ 2)⁻¹) ^ 3 * ⟪deriv (deriv c) θ, J (deriv c θ)⟫ ^ 2 := by
  have hcd : Differentiable ℝ c := hc.differentiable (by simp)
  have hd1 : ContDiff ℝ ∞ (deriv c) := (contDiff_infty_iff_deriv.mp hc).2
  have hd1d : Differentiable ℝ (deriv c) := hd1.differentiable (by simp)
  set S : ℂ →L[ℝ] (E1 →L[ℝ] ℂ) := ContinuousLinearMap.smulRightL ℝ E1 ℂ (Aθ.symm : E1 →L[ℝ] ℝ)
    with hS
  have hfd : fderiv ℝ (fun w : E1 => c (Aθ.symm w)) = fun w => S (deriv c (Aθ.symm w)) :=
    funext fun w => (fderiv_lift_real hcd w).fderiv
  have h2 : HasFDerivAt (fun w : E1 => S (deriv c (Aθ.symm w)))
      (S.comp (ContinuousLinearMap.smulRightL ℝ E1 ℂ (Aθ.symm : E1 →L[ℝ] ℝ)
        (deriv (deriv c) (Aθ.symm (Aθ θ))))) (Aθ θ) :=
    S.hasFDerivAt.comp (Aθ θ) (fderiv_lift_real hd1d (Aθ θ))
  have hdX : dX (fun w : E1 => c (Aθ.symm w)) (Aθ θ) 0 = deriv c θ := by
    simp only [dX, hfd, hS, ContinuousLinearMap.smulRightL_apply_apply, ContinuousLinearMap.smulRight_apply, ContinuousLinearEquiv.coe_coe,
      eb0_coord, one_smul, Aθ_symm_Aθ]
  have hddX : ddX (fun w : E1 => c (Aθ.symm w)) (Aθ θ) 0 0 = deriv (deriv c) θ := by
    rw [ddX, hfd, h2.fderiv]
    simp only [hS, ContinuousLinearMap.comp_apply,
      ContinuousLinearMap.smulRightL_apply_apply, ContinuousLinearMap.smulRight_apply, ContinuousLinearEquiv.coe_coe, eb0_coord, one_smul,
      Aθ_symm_Aθ]
  have hg : gX (fun w : E1 => c (Aθ.symm w)) (Aθ θ) 0 0 = ‖deriv c θ‖ ^ 2 := by
    simp only [gX, hdX, real_inner_self_eq_norm_sq]
  have hginv : (gX (fun w : E1 => c (Aθ.symm w)) (Aθ θ))⁻¹ 0 0 = (‖deriv c θ‖ ^ 2)⁻¹ := by
    rw [Matrix.inv_def, Matrix.adjugate_fin_one, Matrix.det_fin_one, hg]
    simp [Ring.inverse_eq_inv']
  unfold aNormSq gNormSq
  rw [Fintype.sum_unique, Fintype.sum_unique]
  simp only [ginvZ, cX, PiLp.toLp_apply, Fin.prod_univ_three]
  have e0 : ∀ i : Fin 3, (default : Fin 3 → Fin 1) i = 0 := fun i => Subsingleton.elim _ _
  simp only [e0, Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.head_cons, hginv, hdX, hddX]
  ring

theorem det_jacM_liftTau_ne_zero {z : Circle} {w : E1} (hw : w ∈ liftDom z) :
    (jacM (liftTau z) w).det ≠ 0 := by
  rw [Matrix.det_fin_one]
  intro h0
  have hv : fderiv ℝ (liftTau z) w (eb 0) = 0 := by
    ext i
    have hi : i = 0 := Subsingleton.elim _ _
    subst hi
    exact h0
  have := injective_fderiv_liftTau hw (hv.trans (map_zero _).symm)
  have h1 : (eb (0 : Fin 1) : E1) 0 = 1 := by simp [eb]
  rw [this] at h1
  simp at h1

/-- **`|A|²` in a chart of the circle equals `|A|²` of the lift** (any `J`; chartwise smooth, immersed,
Lagrangian). -/
theorem aNormSq_chart_eq_lift (J : ℂ →L[ℝ] ℂ) {F : Circle → ℂ}
    (hsm : ∀ z : Circle,
      ContDiffOn ℝ ∞ (F ∘ (extChartAt 𝓘(ℝ, E1) z).symm) (extChartAt 𝓘(ℝ, E1) z).target)
    (hinj : ∀ z : Circle, ∀ y ∈ (extChartAt 𝓘(ℝ, E1) z).target,
      Injective (fderiv ℝ (F ∘ (extChartAt 𝓘(ℝ, E1) z).symm) y))
    (hlag : ∀ z : Circle, ∀ y ∈ (extChartAt 𝓘(ℝ, E1) z).target, ∀ u v,
      ⟪J (fderiv ℝ (F ∘ (extChartAt 𝓘(ℝ, E1) z).symm) y u),
        fderiv ℝ (F ∘ (extChartAt 𝓘(ℝ, E1) z).symm) y v⟫ = 0)
    {z : Circle} {y : E1} (hy : y ∈ (extChartAt 𝓘(ℝ, E1) z).target) :
    ∃ θ : ℝ, Circle.exp θ = (extChartAt 𝓘(ℝ, E1) z).symm y ∧
      aNormSq J (F ∘ (extChartAt 𝓘(ℝ, E1) z).symm) y =
        ((‖deriv (lift F) θ‖ ^ 2)⁻¹) ^ 3 * ⟪deriv (deriv (lift F)) θ, J (deriv (lift F) θ)⟫ ^ 2 := by
  obtain ⟨θ, hθ⟩ := Circle.exp_surjective ((extChartAt 𝓘(ℝ, E1) z).symm y)
  refine ⟨θ, hθ, ?_⟩
  have hw : Aθ θ ∈ liftDom z := by
    simp only [liftDom, mem_ofPred_eq, Aθ_symm_Aθ, hθ]
    exact (extChartAt 𝓘(ℝ, E1) z).map_target hy
  have hτ : liftTau z (Aθ θ) = y := by
    simp only [liftTau, Aθ_symm_Aθ, hθ]
    exact (extChartAt 𝓘(ℝ, E1) z).right_inv hy
  have hc : ContDiff ℝ ∞ (lift F) := contDiff_lift_of_chart hsm
  have h1 := aNormSq_comp (J := J) (isOpen_extChartAt_target z) (hsm z) (hinj z) (hlag z)
    (isOpen_liftDom z) (contDiffOn_liftTau z) (fun w hw' => liftTau_mem_target hw')
    (fun w hw' => det_jacM_liftTau_ne_zero hw') hw
  rw [hτ] at h1
  rw [← h1, ← aNormSq_congr (lift_eventuallyEq F hw)]
  exact aNormSq_lift_real J hc θ

end CWFlow.Curve
