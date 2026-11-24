import CWFlow.Curve.N1
import CWFlow.Curve.Lift3
import CWFlow.Curve.Energy
import CWFlow.Extension.FlowClosed

/-!
# From `IsFlow` on the circle to the lifted curve diffusion flow (SCOPE W5, (T4))

`isCDF_of_forall_isFlow`, `IsFlow.lift_isCDF`: the lift `(t, θ) ↦ Φ t (e^{iθ})` of a flow is a curve
diffusion flow (`IsCDF`) at interior times. The work is in `CWFlow.Curve.Lift3` (`isCDF_of_chartData`); here the
three fields of `IsFlow` are put in chart form (`Extension.fc_contDiffOn_chart_of_isFlow` for smoothness).
-/

open Set Function Filter
open scoped Manifold Topology ContDiff RealInnerProductSpace

noncomputable section

namespace CWFlow.Curve

open CWFlow

variable {T : ℝ} {Φ : ℝ → Circle → ℂ}

/-- **(T4), family form**: if `Φ` is a flow on every `[0, T']`, `T' < T`, its lift is a curve diffusion
flow on `(0, T)`. -/
theorem isCDF_of_forall_isFlow (hΦ : ∀ T' ∈ Ioo 0 T, IsFlow E1 Jc T' Φ) :
    IsCDF (Ioo 0 T) (fun q : ℝ × ℝ => Φ q.1 (Circle.exp q.2)) := by
  refine isCDF_of_chartData (fun T' hT' z => Extension.fc_contDiffOn_chart_of_isFlow ι1 (hΦ T' hT') z)
    (fun t ht => ?_) (fun t ht x z hx => ?_)
  · obtain ⟨hm, hlt⟩ := Extension.fc_mid ht
    exact (hΦ _ hm).imm t ⟨ht.1.le, hlt.le⟩
  · obtain ⟨hm, hlt⟩ := Extension.fc_mid ht
    exact ((hΦ _ hm).vel t ⟨ht.1.le, hlt.le⟩ x z hx).hasDerivAt (Icc_mem_nhds ht.1 hlt)

/-- **(T4)** The lifted flow is a curve diffusion flow at interior times. -/
theorem IsFlow.lift_isCDF (h : IsFlow E1 Jc T Φ) :
    IsCDF (Ioo 0 T) (fun q : ℝ × ℝ => Φ q.1 (Circle.exp q.2)) :=
  isCDF_of_forall_isFlow fun _ hT' => h.mono hT'.2.le

end CWFlow.Curve
