import CWFlow.Curve.Defs
import Mathlib.Geometry.Manifold.Instances.Sphere
import Mathlib.Geometry.Manifold.ContMDiff.NormedSpace
import Mathlib.Geometry.Manifold.ContMDiff.Atlas
import Mathlib.Analysis.SpecialFunctions.Complex.Circle
import Mathlib.Analysis.SpecialFunctions.ExpDeriv
import Mathlib.LinearAlgebra.FiniteDimensional.Basic

/-!
# The circle's charts read through `θ ↦ e^{iθ}` (SCOPE W5, Mathlib part)

For `z : Circle` the chart `φ_z = extChartAt 𝓘(ℝ, E1) z` composed with `θ ↦ e^{iθ}` gives a local
diffeomorphism `liftTau z : E1 → E1`, `w ↦ φ_z (e^{i w₀})`, on the open set `liftDom z`. For any map `F` on
the circle, `w ↦ F (e^{i w₀})` agrees with `F ∘ φ_z⁻¹ ∘ liftTau z` on `liftDom z`. From chartwise
smoothness and immersion of `F` this gives smoothness and regularity of the lift `θ ↦ F (e^{iθ})`.

`Aθ : ℝ ≃L[ℝ] E1` is the identification `a ↦ (a)`.
-/

open Set Function Filter
open scoped Manifold Topology ContDiff

noncomputable section

/-- `Circle` is a smooth (`C^∞`) manifold; Mathlib states it for `ω`. -/
instance instIsManifoldCircleInfty : IsManifold 𝓘(ℝ, EuclideanSpace ℝ (Fin 1)) ∞ Circle :=
  IsManifold.of_le (n := ω) le_top

namespace CWFlow.Curve

/-- The linear isomorphism `ℝ¹ ≃ E1` used for the charts `chartInv ι1 z` of the repository: the identity. -/
abbrev ι1 : EuclideanSpace ℝ (Fin 1) ≃L[ℝ] E1 := ContinuousLinearEquiv.refl ℝ E1

/-- `ℝ ≃ E1 = ℝ¹`. -/
def Aθ : ℝ ≃L[ℝ] E1 :=
  ((EuclideanSpace.equiv (Fin 1) ℝ).trans (ContinuousLinearEquiv.funUnique (Fin 1) ℝ ℝ)).symm

@[simp] theorem Aθ_symm_apply (w : E1) : Aθ.symm w = w 0 := rfl

@[simp] theorem Aθ_symm_Aθ (a : ℝ) : Aθ.symm (Aθ a) = a := Aθ.symm_apply_apply a

theorem Aθ_one_ne_zero : Aθ 1 ≠ 0 := by
  intro h
  have := congrArg Aθ.symm h
  simp at this

/-- The domain of `liftTau z`. -/
def liftDom (z : Circle) : Set E1 :=
  {w | Circle.exp (Aθ.symm w) ∈ (extChartAt 𝓘(ℝ, E1) z).source}

/-- The chart at `z` read through `θ ↦ e^{iθ}`. -/
def liftTau (z : Circle) (w : E1) : E1 := extChartAt 𝓘(ℝ, E1) z (Circle.exp (Aθ.symm w))

theorem continuous_expA : Continuous fun w : E1 => Circle.exp (Aθ.symm w) :=
  Circle.exp.continuous.comp Aθ.symm.continuous

theorem isOpen_liftDom (z : Circle) : IsOpen (liftDom z) :=
  (isOpen_extChartAt_source z).preimage continuous_expA

theorem contMDiff_expA : ContMDiff 𝓘(ℝ, E1) 𝓘(ℝ, E1) ∞ fun w : E1 => Circle.exp (Aθ.symm w) :=
  contMDiff_circleExp.comp Aθ.symm.contDiff.contMDiff

theorem contDiffOn_liftTau (z : Circle) : ContDiffOn ℝ ∞ (liftTau z) (liftDom z) := by
  have h : ContMDiffOn 𝓘(ℝ, E1) 𝓘(ℝ, E1) ∞ (liftTau z) (liftDom z) := by
    refine ContMDiffOn.comp (t := (chartAt E1 z).source) ?_ contMDiff_expA.contMDiffOn ?_
    · exact contMDiffOn_extChartAt
    · intro w hw
      have : Circle.exp (Aθ.symm w) ∈ (extChartAt 𝓘(ℝ, E1) z).source := hw
      rwa [extChartAt_source] at this
  exact contMDiffOn_iff_contDiffOn.1 h

theorem liftTau_mem_target {z : Circle} {w : E1} (hw : w ∈ liftDom z) :
    liftTau z w ∈ (extChartAt 𝓘(ℝ, E1) z).target :=
  (extChartAt 𝓘(ℝ, E1) z).map_source hw

theorem chartSymm_liftTau {z : Circle} {w : E1} (hw : w ∈ liftDom z) :
    (extChartAt 𝓘(ℝ, E1) z).symm (liftTau z w) = Circle.exp (Aθ.symm w) :=
  (extChartAt 𝓘(ℝ, E1) z).left_inv hw

/-- `w ↦ F (e^{i w₀})` is `F ∘ φ_z⁻¹ ∘ liftTau z` on `liftDom z`. -/
theorem lift_eq_comp {V : Type*} (F : Circle → V) {z : Circle} {w : E1} (hw : w ∈ liftDom z) :
    F (Circle.exp (Aθ.symm w)) = (F ∘ (extChartAt 𝓘(ℝ, E1) z).symm) (liftTau z w) := by
  simp only [Function.comp_apply, chartSymm_liftTau hw]

theorem lift_eventuallyEq {V : Type*} (F : Circle → V) {z : Circle} {w : E1} (hw : w ∈ liftDom z) :
    (fun w => F (Circle.exp (Aθ.symm w))) =ᶠ[𝓝 w] (F ∘ (extChartAt 𝓘(ℝ, E1) z).symm) ∘ liftTau z := by
  filter_upwards [(isOpen_liftDom z).mem_nhds hw] with w' hw' using lift_eq_comp F hw'

/-- `θ ↦ e^{iθ}` as a map into `ℂ`, read on `E1`. -/
theorem hasFDerivAt_coe_expA (w : E1) :
    HasFDerivAt (fun w : E1 => ((Circle.exp (Aθ.symm w) : Circle) : ℂ))
      ((ContinuousLinearMap.smulRight (1 : ℝ →L[ℝ] ℝ)
        (Complex.exp (Aθ.symm w * Complex.I) * Complex.I)).comp (Aθ.symm : E1 →L[ℝ] ℝ)) w := by
  have h1 : HasDerivAt (fun t : ℝ => Complex.exp (t * Complex.I))
      (Complex.exp (Aθ.symm w * Complex.I) * Complex.I) (Aθ.symm w) := by
    have := ((Complex.ofRealCLM.hasDerivAt (x := Aθ.symm w)).mul_const Complex.I).cexp
    simpa using this
  exact h1.hasFDerivAt.comp w (Aθ.symm : E1 →L[ℝ] ℝ).hasFDerivAt

/-- The derivative of `liftTau z` is injective. -/
theorem injective_fderiv_liftTau {z : Circle} {w : E1} (hw : w ∈ liftDom z) :
    Injective (fderiv ℝ (liftTau z) w) := by
  -- `ψ = coe ∘ φ_z⁻¹` is smooth on the target, and `ψ ∘ liftTau z = coe ∘ exp ∘ Aθ⁻¹` near `w`
  have hψ : ContDiffOn ℝ ∞ (fun y : E1 => (((extChartAt 𝓘(ℝ, E1) z).symm y : Circle) : ℂ))
      (extChartAt 𝓘(ℝ, E1) z).target := by
    have h : ContMDiffOn 𝓘(ℝ, E1) 𝓘(ℝ, ℂ) ∞
        (fun y : E1 => (((extChartAt 𝓘(ℝ, E1) z).symm y : Circle) : ℂ))
        (extChartAt 𝓘(ℝ, E1) z).target := by
      have : Fact (Module.finrank ℝ ℂ = 1 + 1) := ⟨by simp⟩
      have hc : ContMDiff (𝓡 1) 𝓘(ℝ, ℂ) ∞ (fun p : Circle => (p : ℂ)) := contMDiff_coe_sphere
      exact hc.comp_contMDiffOn (contMDiffOn_extChartAt_symm (I := 𝓘(ℝ, E1)) (n := ∞) z)
    exact contMDiffOn_iff_contDiffOn.1 h
  have hT := liftTau_mem_target hw
  have hψd : DifferentiableAt ℝ (fun y : E1 => (((extChartAt 𝓘(ℝ, E1) z).symm y : Circle) : ℂ))
      (liftTau z w) :=
    (hψ.contDiffAt ((isOpen_extChartAt_target z).mem_nhds hT)).differentiableAt (by simp)
  have hτd : DifferentiableAt ℝ (liftTau z) w :=
    ((contDiffOn_liftTau z).contDiffAt ((isOpen_liftDom z).mem_nhds hw)).differentiableAt (by simp)
  have hev : (fun w : E1 => ((Circle.exp (Aθ.symm w) : Circle) : ℂ)) =ᶠ[𝓝 w]
      (fun y : E1 => (((extChartAt 𝓘(ℝ, E1) z).symm y : Circle) : ℂ)) ∘ liftTau z := by
    filter_upwards [(isOpen_liftDom z).mem_nhds hw] with w' hw'
    simp only [Function.comp_apply, chartSymm_liftTau hw']
  have hchain := (hψd.hasFDerivAt.comp w hτd.hasFDerivAt).congr_of_eventuallyEq hev
  have huniq := (hasFDerivAt_coe_expA w).unique hchain
  intro u v huv
  have hu : (fderiv ℝ (liftTau z) w) (u - v) = 0 := by rw [map_sub, huv, sub_self]
  have h0 := congrArg (fun L : E1 →L[ℝ] ℂ => L (u - v)) huniq
  simp only [ContinuousLinearMap.smulRight_apply, ContinuousLinearMap.comp_apply, hu, map_zero,
    ContinuousLinearMap.one_apply, ContinuousLinearEquiv.coe_coe] at h0
  have hne : Complex.exp (Aθ.symm w * Complex.I) * Complex.I ≠ 0 :=
    mul_ne_zero (Complex.exp_ne_zero _) Complex.I_ne_zero
  have h1 : Aθ.symm (u - v) = 0 := by
    rcases smul_eq_zero.1 h0 with h | h
    · exact h
    · exact absurd h hne
  have : u - v = 0 := Aθ.symm.injective (h1.trans (map_zero _).symm)
  exact sub_eq_zero.1 this

/-- The derivative of `liftTau z` is an isomorphism. -/
theorem exists_equiv_fderiv_liftTau {z : Circle} {w : E1} (hw : w ∈ liftDom z) :
    ∃ e : E1 ≃L[ℝ] E1, (e : E1 →L[ℝ] E1) = fderiv ℝ (liftTau z) w := by
  set L := fderiv ℝ (liftTau z) w
  have hinj : Injective (L : E1 →ₗ[ℝ] E1) := injective_fderiv_liftTau hw
  let e : E1 ≃ₗ[ℝ] E1 := LinearEquiv.ofInjectiveEndo (L : E1 →ₗ[ℝ] E1) hinj
  refine ⟨e.toContinuousLinearEquiv, ?_⟩
  ext v
  rfl

theorem fderiv_lift_real {c : ℝ → ℂ} (hc : Differentiable ℝ c) (w : E1) :
    HasFDerivAt (fun w : E1 => c (Aθ.symm w))
      (ContinuousLinearMap.smulRightL ℝ E1 ℂ (Aθ.symm : E1 →L[ℝ] ℝ) (deriv c (Aθ.symm w))) w := by
  have h := ((hc (Aθ.symm w)).hasDerivAt.hasFDerivAt).comp w (Aθ.symm : E1 →L[ℝ] ℝ).hasFDerivAt
  refine h.congr_fderiv (ContinuousLinearMap.ext fun v => ?_)
  simp [ContinuousLinearMap.smulRight_apply]

/-! ### Smoothness and regularity of the lift -/

section Lift

variable {V : Type*} [NormedAddCommGroup V] [NormedSpace ℝ V] {F : Circle → V}

theorem exists_liftDom (θ : ℝ) : Aθ θ ∈ liftDom (Circle.exp θ) := by
  simp only [liftDom, mem_ofPred_eq, Aθ_symm_Aθ]
  exact mem_extChartAt_source _

/-- The lift `θ ↦ F (e^{iθ})` near `θ`, as `F ∘ φ⁻¹ ∘ liftTau ∘ Aθ`. -/
theorem lift_eventuallyEq_real (F : Circle → V) (θ : ℝ) :
    (fun t => F (Circle.exp t)) =ᶠ[𝓝 θ]
      ((F ∘ (extChartAt 𝓘(ℝ, E1) (Circle.exp θ)).symm) ∘ liftTau (Circle.exp θ)) ∘ Aθ := by
  have h := lift_eventuallyEq F (exists_liftDom θ)
  have hc : Tendsto Aθ (𝓝 θ) (𝓝 (Aθ θ)) := Aθ.continuous.tendsto θ
  filter_upwards [hc.eventually h] with t ht
  simp only [Function.comp_apply] at ht ⊢
  rw [← ht, Aθ_symm_Aθ]

variable (hsm : ∀ z : Circle,
    ContDiffOn ℝ ∞ (F ∘ (extChartAt 𝓘(ℝ, E1) z).symm) (extChartAt 𝓘(ℝ, E1) z).target)
include hsm

theorem contDiffAt_chart_liftTau (θ : ℝ) :
    ContDiffAt ℝ ∞ (((F ∘ (extChartAt 𝓘(ℝ, E1) (Circle.exp θ)).symm) ∘
      liftTau (Circle.exp θ)) ∘ Aθ) θ := by
  set z := Circle.exp θ
  have hw : Aθ θ ∈ liftDom z := exists_liftDom θ
  have h1 : ContDiffAt ℝ ∞ (F ∘ (extChartAt 𝓘(ℝ, E1) z).symm) (liftTau z (Aθ θ)) :=
    (hsm z).contDiffAt ((isOpen_extChartAt_target z).mem_nhds (liftTau_mem_target hw))
  have h2 : ContDiffAt ℝ ∞ (liftTau z) (Aθ θ) :=
    (contDiffOn_liftTau z).contDiffAt ((isOpen_liftDom z).mem_nhds hw)
  exact (h1.comp _ h2).comp θ Aθ.contDiff.contDiffAt

/-- The lift of a chartwise smooth map is smooth. -/
theorem contDiff_lift_of_chart : ContDiff ℝ ∞ (fun t => F (Circle.exp t)) :=
  contDiff_iff_contDiffAt.2 fun θ =>
    (contDiffAt_chart_liftTau hsm θ).congr_of_eventuallyEq (lift_eventuallyEq_real F θ)

/-- The lift of a chartwise immersion is regular. -/
theorem deriv_lift_ne_zero_of_chart
    (hinj : ∀ z : Circle, ∀ y ∈ (extChartAt 𝓘(ℝ, E1) z).target,
      Injective (fderiv ℝ (F ∘ (extChartAt 𝓘(ℝ, E1) z).symm) y)) (θ : ℝ) :
    deriv (fun t => F (Circle.exp t)) θ ≠ 0 := by
  set z := Circle.exp θ
  have hw : Aθ θ ∈ liftDom z := exists_liftDom θ
  have hT := liftTau_mem_target hw
  have h1 : DifferentiableAt ℝ (F ∘ (extChartAt 𝓘(ℝ, E1) z).symm) (liftTau z (Aθ θ)) :=
    ((hsm z).contDiffAt ((isOpen_extChartAt_target z).mem_nhds hT)).differentiableAt (by simp)
  have h2 : DifferentiableAt ℝ (liftTau z) (Aθ θ) :=
    ((contDiffOn_liftTau z).contDiffAt ((isOpen_liftDom z).mem_nhds hw)).differentiableAt
      (by simp)
  have hd : HasFDerivAt (((F ∘ (extChartAt 𝓘(ℝ, E1) z).symm) ∘ liftTau z) ∘ Aθ)
      (((fderiv ℝ (F ∘ (extChartAt 𝓘(ℝ, E1) z).symm) (liftTau z (Aθ θ))).comp
        (fderiv ℝ (liftTau z) (Aθ θ))).comp (Aθ : ℝ →L[ℝ] E1)) θ :=
    (h1.hasFDerivAt.comp _ h2.hasFDerivAt).comp θ Aθ.hasFDerivAt
  have hd' := (hd.congr_of_eventuallyEq (lift_eventuallyEq_real F θ)).hasDerivAt
  rw [hd'.deriv]
  simp only [ContinuousLinearMap.comp_apply, ContinuousLinearEquiv.coe_coe]
  intro h0
  have := (hinj z _ hT) (h0.trans (map_zero _).symm)
  exact Aθ_one_ne_zero ((injective_fderiv_liftTau hw) (this.trans (map_zero _).symm))

end Lift

end CWFlow.Curve
