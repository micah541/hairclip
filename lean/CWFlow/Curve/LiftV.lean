import CWFlow.Curve.CircleChart
import CWFlow.Curve.N1Basic
import CWFlow.Scalar.LinearChange
import CWFlow.Flow.TanAlgebra

/-!
# The velocity `J∇div(JH)` read in a chart of the circle is `−κ_ss N` of the lift (SCOPE W5, part (T4))

For `F : Circle → ℂ` smooth with injective derivative in every chart, `flowVel Jc` of the chart representation
`F ∘ φ_z⁻¹` at `φ_z (e^{iθ})` equals `−κ_ss N` of the lift `θ ↦ F (e^{iθ})` at `θ`. Proof: the repository's
invariance lemmas `flowVel_comp'` (reparametrisation by `liftTau z`), `flowVel_congr_of_eventuallyEq`
(locality) and `flowVel_comp_linearEquiv` (`Aθ : ℝ ≃ E1`), then `flowVel_eq_neg_kss` for the lift.
-/

open Set Function Filter
open scoped Manifold Topology ContDiff RealInnerProductSpace

noncomputable section

namespace CWFlow.Curve

open CWFlow

theorem injective_fderiv_lift_real {c : ℝ → ℂ} (hc : Differentiable ℝ c)
    (hreg : ∀ θ, deriv c θ ≠ 0) (w : E1) :
    Injective (fderiv ℝ (fun w : E1 => c (Aθ.symm w)) w) := by
  rw [(fderiv_lift_real hc w).fderiv]
  intro u v huv
  simp only [ContinuousLinearMap.smulRightL_apply_apply, ContinuousLinearMap.smulRight_apply,
    ContinuousLinearEquiv.coe_coe] at huv
  have h1 : Aθ.symm u = Aθ.symm v := smul_left_injective ℝ (hreg _) huv
  exact Aθ.symm.injective h1

variable {F : Circle → ℂ}
  (hsm : ∀ z : Circle,
    ContDiffOn ℝ ∞ (F ∘ (extChartAt 𝓘(ℝ, E1) z).symm) (extChartAt 𝓘(ℝ, E1) z).target)
  (hinj : ∀ z : Circle, ∀ y ∈ (extChartAt 𝓘(ℝ, E1) z).target,
    Injective (fderiv ℝ (F ∘ (extChartAt 𝓘(ℝ, E1) z).symm) y))
include hsm hinj

/-- **(T4), pointwise form**: the chart velocity is `−κ_ss N` of the lift. -/
theorem flowVel_chart_eq_lift {z : Circle} {θ : ℝ}
    (hθ : Circle.exp θ ∈ (extChartAt 𝓘(ℝ, E1) z).source) :
    flowVel Jc (F ∘ (extChartAt 𝓘(ℝ, E1) z).symm) (extChartAt 𝓘(ℝ, E1) z (Circle.exp θ)) =
      -((dS (lift F) (dS (lift F) (curv (lift F))) θ : ℝ) : ℂ) *
        (Complex.I * (deriv (lift F) θ / (‖deriv (lift F) θ‖ : ℂ))) := by
  have hc : ContDiff ℝ ∞ (lift F) := contDiff_lift_of_chart hsm
  have hreg : ∀ t, deriv (lift F) t ≠ 0 := deriv_lift_ne_zero_of_chart hsm hinj
  have hcd : Differentiable ℝ (lift F) := hc.differentiable (by simp)
  set Y : E1 → ℂ := fun w => lift F (Aθ.symm w) with hY
  have hYs : ContDiff ℝ ∞ Y := hc.comp Aθ.symm.contDiff
  have hYi : ∀ w, Injective (fderiv ℝ Y w) := injective_fderiv_lift_real hcd hreg
  have hw : Aθ θ ∈ liftDom z := by
    simp only [liftDom, mem_ofPred_eq, Aθ_symm_Aθ]; exact hθ
  -- step 1: linear change `Aθ`
  have e1 : flowVel Jc (lift F) θ = flowVel Jc Y (Aθ θ) := by
    have h := flowVel_comp_linearEquiv Aθ hYs hYi Jc θ
    have hYA : Y ∘ Aθ = lift F := funext fun t => by
      show lift F (Aθ.symm (Aθ t)) = lift F t
      rw [Aθ_symm_Aθ]
    rw [hYA] at h
    exact h
  -- step 2: locality, `Y = X ∘ liftTau z` near `Aθ θ`
  have e2 : flowVel Jc Y (Aθ θ) =
      flowVel Jc ((F ∘ (extChartAt 𝓘(ℝ, E1) z).symm) ∘ liftTau z) (Aθ θ) :=
    flowVel_congr_of_eventuallyEq (lift_eventuallyEq F hw)
  -- step 3: reparametrisation by `liftTau z`
  have e3 := flowVel_comp' (isOpen_extChartAt_target z) (hsm z) (hinj z) Jc (isOpen_liftDom z)
    (contDiffOn_liftTau z) (fun w hw' => exists_equiv_fderiv_liftTau hw')
    (by intro w hw'; exact liftTau_mem_target hw') hw
  have e4 : liftTau z (Aθ θ) = extChartAt 𝓘(ℝ, E1) z (Circle.exp θ) := by
    simp [liftTau]
  rw [← e4, ← e3, ← e2, ← e1]
  exact flowVel_eq_neg_kss hc hreg θ

end CWFlow.Curve
