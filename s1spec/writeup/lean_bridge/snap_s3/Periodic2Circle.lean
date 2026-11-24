import CWFlow.Curve.CircleChart
import CWFlow.Curve.N1Basic
import CWFlow.Weinstein.Global

/-!
# From a `2π`-periodic regular curve to an immersion of `Circle` (review W7, item S3a)

Converse of `contDiff_lift_of_chart` / `deriv_lift_ne_zero_of_chart` (CircleChart.lean): a `2π`-periodic
`C^∞` curve `c : ℝ → ℂ` with `c' ≠ 0` everywhere is `lift F` for a unique `F : Circle → ℂ`
(`circleOfPeriodic`), and this `F` is an `IsLagImmersion E1 Jc` (smooth and immersive in every chart of
`Circle`; the Lagrangian clause is automatic at `n = 1`, `lag_dim_one`).

Chart argument: near `y` in the target of the chart at `z`, `F ∘ φ_z⁻¹ = c ∘ Aθ⁻¹ ∘ g`, where `g` is a
smooth local inverse (inverse function theorem) of `liftTau z = φ_z ∘ exp ∘ Aθ⁻¹`.
-/

open Set Function Filter
open scoped Manifold Topology ContDiff

noncomputable section

namespace CWFlow.Curve

open CWFlow

/-! ### Chart readings of a map whose lift is a smooth regular curve -/

section Chart

set_option linter.unusedSectionVars false

variable {V : Type*} [NormedAddCommGroup V] [NormedSpace ℝ V]

/-- Near any point of the target of the chart at `z`, `F ∘ φ_z⁻¹` is `c ∘ Aθ⁻¹ ∘ g` for a smooth `g`
with invertible derivative, provided `F (e^{iθ}) = c θ`. -/
theorem exists_chart_local_form (F : Circle → V) (c : ℝ → V) (hF : ∀ θ, F (Circle.exp θ) = c θ)
    (z : Circle) {y : E1} (hy : y ∈ (extChartAt 𝓘(ℝ, E1) z).target) :
    ∃ g : E1 → E1, ∃ e : E1 ≃L[ℝ] E1, ContDiffAt ℝ ∞ g y ∧ HasFDerivAt g (e : E1 →L[ℝ] E1) y ∧
      (F ∘ (extChartAt 𝓘(ℝ, E1) z).symm) =ᶠ[𝓝 y] (fun w => c (Aθ.symm w)) ∘ g := by
  set p := (extChartAt 𝓘(ℝ, E1) z).symm y with hp_def
  have hp : p ∈ (extChartAt 𝓘(ℝ, E1) z).source := (extChartAt 𝓘(ℝ, E1) z).map_target hy
  obtain ⟨θ₀, hθ₀⟩ := Circle.exp_surjective p
  set w₀ : E1 := Aθ θ₀ with hw₀_def
  have hw₀ : w₀ ∈ liftDom z := by
    show Circle.exp (Aθ.symm (Aθ θ₀)) ∈ (extChartAt 𝓘(ℝ, E1) z).source
    rw [Aθ_symm_Aθ, hθ₀]
    exact hp
  have hyw : liftTau z w₀ = y := by
    show extChartAt 𝓘(ℝ, E1) z (Circle.exp (Aθ.symm (Aθ θ₀))) = y
    rw [Aθ_symm_Aθ, hθ₀, hp_def]
    exact (extChartAt 𝓘(ℝ, E1) z).right_inv hy
  subst hyw
  have hτ : ContDiffAt ℝ ∞ (liftTau z) w₀ :=
    (contDiffOn_liftTau z).contDiffAt ((isOpen_liftDom z).mem_nhds hw₀)
  obtain ⟨e, he⟩ := exists_equiv_fderiv_liftTau hw₀
  have hτ' : HasFDerivAt (liftTau z) (e : E1 →L[ℝ] E1) w₀ := by
    rw [he]; exact (hτ.differentiableAt (by simp)).hasFDerivAt
  have hn : (∞ : WithTop ℕ∞) ≠ 0 := by simp
  have hs := hτ.hasStrictFDerivAt' hτ' hn
  set g := hs.localInverse (liftTau z) e w₀ with hg_def
  have hg_sm : ContDiffAt ℝ ∞ g (liftTau z w₀) := hτ.to_localInverse hτ' hn
  have hg_d : HasFDerivAt g (e.symm : E1 →L[ℝ] E1) (liftTau z w₀) :=
    hs.to_localInverse.hasFDerivAt
  have hg_c : Tendsto g (𝓝 (liftTau z w₀)) (𝓝 w₀) := hs.localInverse_tendsto
  refine ⟨g, e.symm, hg_sm, hg_d, ?_⟩
  filter_upwards [hs.eventually_right_inverse,
    hg_c.eventually ((isOpen_liftDom z).mem_nhds hw₀)] with y' h1 h2
  simp only [Function.comp_apply]
  rw [← hF, ← chartSymm_liftTau h2, h1]

variable {c : ℝ → V} {F : Circle → V}

/-- Chartwise smoothness of `F` when `F (e^{iθ}) = c θ` with `c` smooth. -/
theorem contDiffOn_chart_of_lift (hF : ∀ θ, F (Circle.exp θ) = c θ) (hc : ContDiff ℝ ∞ c)
    (z : Circle) :
    ContDiffOn ℝ ∞ (F ∘ (extChartAt 𝓘(ℝ, E1) z).symm) (extChartAt 𝓘(ℝ, E1) z).target := by
  intro y hy
  obtain ⟨g, e, hg, -, hev⟩ := exists_chart_local_form F c hF z hy
  have h1 : ContDiffAt ℝ ∞ (fun w : E1 => c (Aθ.symm w)) (g y) :=
    hc.contDiffAt.comp _ Aθ.symm.contDiff.contDiffAt
  exact ((h1.comp y hg).congr_of_eventuallyEq hev).contDiffWithinAt

/-- Chartwise immersion of `F` when `F (e^{iθ}) = c θ` with `c` regular. -/
theorem injective_fderiv_chart_of_lift {c : ℝ → ℂ} {F : Circle → ℂ}
    (hF : ∀ θ, F (Circle.exp θ) = c θ)
    (hc : ContDiff ℝ ∞ c) (hreg : ∀ θ, deriv c θ ≠ 0) (z : Circle)
    (y : E1) (hy : y ∈ (extChartAt 𝓘(ℝ, E1) z).target) :
    Injective (fderiv ℝ (F ∘ (extChartAt 𝓘(ℝ, E1) z).symm) y) := by
  obtain ⟨g, e, -, hg', hev⟩ := exists_chart_local_form F c hF z hy
  have hcd : Differentiable ℝ c := hc.differentiable (by simp)
  have hd := (fderiv_lift_real hcd (g y)).comp y hg'
  rw [hev.fderiv_eq, hd.fderiv]
  intro u v huv
  have h0 : ((ContinuousLinearMap.smulRightL ℝ E1 ℂ (Aθ.symm : E1 →L[ℝ] ℝ)
      (deriv c (Aθ.symm (g y)))).comp (e : E1 →L[ℝ] E1)) (u - v) = 0 := by
    rw [map_sub, huv, sub_self]
  simp only [ContinuousLinearMap.comp_apply, ContinuousLinearEquiv.coe_coe] at h0
  change (Aθ.symm (e (u - v))) • deriv c (Aθ.symm (g y)) = 0 at h0
  rcases smul_eq_zero.1 h0 with h | h
  · have h' : e (u - v) = 0 := Aθ.symm.injective (h.trans (map_zero _).symm)
    exact sub_eq_zero.1 (e.injective (h'.trans (map_zero _).symm))
  · exact absurd h (hreg _)

end Chart

/-! ### The map on `Circle` defined by a periodic curve -/

/-- The map `F : Circle → V` with `F (e^{iθ}) = c θ`, for `c` `2π`-periodic
(via `AddCircle (2π) ≃ₜ Circle`, `AddCircle.homeomorphCircle'`). -/
def circleOfPeriodic {V : Type*} {c : ℝ → V} (hc : Periodic c (2 * Real.pi)) : Circle → V :=
  fun z => hc.lift (AddCircle.homeomorphCircle'.symm z)

theorem circleOfPeriodic_exp {V : Type*} {c : ℝ → V} (hc : Periodic c (2 * Real.pi)) (θ : ℝ) :
    circleOfPeriodic hc (Circle.exp θ) = c θ := by
  have h : AddCircle.homeomorphCircle'.symm (Circle.exp θ) = (θ : AddCircle (2 * Real.pi)) := by
    rw [Homeomorph.symm_apply_eq, AddCircle.homeomorphCircle'_apply_mk]
  simp only [circleOfPeriodic, h]
  exact Periodic.lift_coe hc θ

/-- **S3a.** A `2π`-periodic `C^∞` curve `c : ℝ → ℂ` with `c' ≠ 0` everywhere is the lift of an immersion
`F : Circle → ℂ` (in the sense of `IsLagImmersion E1 Jc`, the hypothesis of `cdf_exists_maximal`). -/
theorem exists_isLagImmersion_of_periodic {c : ℝ → ℂ} (hc : ContDiff ℝ ∞ c)
    (hper : Periodic c (2 * Real.pi)) (hreg : ∀ θ, deriv c θ ≠ 0) :
    ∃ F : Circle → ℂ, lift F = c ∧ IsLagImmersion E1 Jc F := by
  have hF : ∀ θ, circleOfPeriodic hper (Circle.exp θ) = c θ := circleOfPeriodic_exp hper
  refine ⟨circleOfPeriodic hper, funext hF, ⟨fun z => ?_, fun z y hy => ?_, fun _ _ _ u w => ?_⟩⟩
  · exact contDiffOn_chart_of_lift hF hc z
  · exact injective_fderiv_chart_of_lift hF hc hreg z y hy
  · exact lag_dim_one _ u w

/-- Uniqueness of the map in S3a: `lift` is injective (`Circle.exp` is surjective). -/
theorem lift_injective : Injective (lift : (Circle → ℂ) → ℝ → ℂ) := by
  intro F G h
  funext p
  obtain ⟨θ, rfl⟩ := Circle.exp_surjective p
  exact congrFun h θ

/-! ## Axiom audit -/

/-- info: 'CWFlow.Curve.exists_isLagImmersion_of_periodic' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_isLagImmersion_of_periodic
/-- info: 'CWFlow.Curve.lift_injective' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms lift_injective

end CWFlow.Curve
