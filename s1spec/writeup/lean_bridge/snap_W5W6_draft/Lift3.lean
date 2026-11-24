import CWFlow.Curve.LiftA
import CWFlow.Curve.LiftC
import CWFlow.Curve.LiftV
import CWFlow.Curve.Arclength
import CWFlow.Weinstein.Global
import CWFlow.Extension.IntM

/-!
# `|A|²` in charts and the lifted flow, from chart data (SCOPE W5)

* `aNormSq_eq_curv_sq` (T3): `|A|²` read in the chart `chartInv ι1 z` is `κ²` of the lift `θ ↦ F (e^{iθ})`.
* `isCDF_of_chartData`: the core of (T4), with the three fields of an `IsFlow` stated in chart form (joint
  smoothness of `(t, y) ↦ Φ t (chartInv ι1 z y)`, immersion, velocity `flowVel` at interior times); the lift
  `(t, θ) ↦ Φ t (e^{iθ})` is then `IsCDF` on `(0, T)`. `CWFlow.Curve.Lift` feeds an `IsFlow` into it.

This file avoids the import closure of `IsFlow` (about 4.5 GB of `.olean`), so it can be compiled separately.
-/

open Set Function Filter
open scoped Manifold Topology ContDiff RealInnerProductSpace

noncomputable section

namespace CWFlow.Curve

open CWFlow

/-- **(T3)** `|A|²` in an `ι1`-chart is `κ²` of the lift. -/
theorem aNormSq_eq_curv_sq {F : Circle → ℂ} (hF : IsLagImmersion E1 Jc F) (z : Circle)
    {y : E1} (hy : y ∈ Extension.chartDomE ι1 z) :
    ∃ θ : ℝ, Circle.exp θ = Extension.chartInv ι1 z y ∧
      Extension.aNormSq Jc (F ∘ Extension.chartInv ι1 z) y = curv (lift F) θ ^ 2 := by
  have hy' : y ∈ (extChartAt 𝓘(ℝ, E1) z).target := hy
  obtain ⟨θ, hθ, h⟩ := aNormSq_chart_eq_lift Jc (fun z => hF.smooth z) (fun z => hF.inj z)
    (fun z => hF.lag z) hy'
  refine ⟨θ, hθ, ?_⟩
  have e : F ∘ Extension.chartInv ι1 z = F ∘ (extChartAt 𝓘(ℝ, E1) z).symm := rfl
  rw [e, h]
  have hreg : deriv (lift F) θ ≠ 0 :=
    deriv_lift_ne_zero_of_chart (fun z => hF.smooth z) (fun z => hF.inj z) θ
  have hn : ‖deriv (lift F) θ‖ ≠ 0 := norm_ne_zero_iff.2 hreg
  simp only [curv, Jc_apply]
  field_simp
  try ring

section Flow

variable {T : ℝ} {Φ : ℝ → Circle → ℂ}

/-- **(T4), chart-data form**: the lift of a family that is, on every `[0, T']` with `T' < T`, jointly smooth in
charts, immersed, and moving by `flowVel` at interior times, is a curve diffusion flow on `(0, T)`. -/
theorem isCDF_of_chartData
    (hsm : ∀ T' ∈ Ioo 0 T, ∀ z : Circle,
      ContDiffOn ℝ ∞ (fun q : ℝ × E1 => Φ q.1 (Extension.chartInv ι1 z q.2))
        (Icc 0 T' ×ˢ Extension.chartDomE ι1 z))
    (himm : ∀ t ∈ Ioo 0 T, IsLagImmersion E1 Jc (Φ t))
    (hvel : ∀ t ∈ Ioo 0 T, ∀ x z : Circle, x ∈ (chart E1 z).source →
      HasDerivAt (fun s => Φ s x) (flowVel Jc (chartRep E1 (Φ t) z) (chart E1 z x)) t) :
    IsCDF (Ioo 0 T) (fun q : ℝ × ℝ => Φ q.1 (Circle.exp q.2)) where
  isOpen := isOpen_Ioo
  smooth := by
    intro q hq
    have hm : (q.1 + T) / 2 ∈ Ioo 0 T := ⟨by linarith [hq.1.1, hq.1.2], by linarith [hq.1.2]⟩
    exact (contDiffAt_lift_of_chart (hsm _ hm) ⟨hq.1.1, by linarith [hq.1.2]⟩ q.2).contDiffWithinAt
  periodic := fun t θ => by
    show Φ t (Circle.exp (θ + 2 * Real.pi)) = Φ t (Circle.exp θ)
    rw [Circle.periodic_exp θ]
  reg := fun t ht θ =>
    deriv_lift_ne_zero_of_chart (fun z => (himm t ht).smooth z) (fun z => (himm t ht).inj z) θ
  vel := fun t ht θ => by
    set z := Circle.exp θ
    have hx : Circle.exp θ ∈ (chart E1 z).source := mem_extChartAt_source _
    have h2 := hvel t ht (Circle.exp θ) z hx
    have e : flowVel Jc (chartRep E1 (Φ t) z) (chart E1 z (Circle.exp θ)) =
        -((dS (lift (Φ t)) (dS (lift (Φ t)) (curv (lift (Φ t)))) θ : ℝ) : ℂ) *
          (Complex.I * (deriv (lift (Φ t)) θ / (‖deriv (lift (Φ t)) θ‖ : ℂ))) :=
      flowVel_chart_eq_lift (fun z => (himm t ht).smooth z) (fun z => (himm t ht).inj z) hx
    rw [e] at h2
    exact h2

end Flow

/-- info: 'CWFlow.Curve.aNormSq_eq_curv_sq' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms aNormSq_eq_curv_sq
/-- info: 'CWFlow.Curve.isCDF_of_chartData' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms isCDF_of_chartData

end CWFlow.Curve
