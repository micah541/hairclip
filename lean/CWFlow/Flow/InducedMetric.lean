import CWFlow.Flow.GraphImm
import DifferentialGeometry.Geometry.Metric.Construction.SmoothMetricFromCoefficients

/-!
# The metric induced by an immersion

An immersion `F : M → V` into an inner product space induces a smooth Riemannian metric on `M`,
`g_x(v, w) = ⟨DF v, DF w⟩` (`exists_inducedMetric`), with `DF` read in the chart at `x` itself. In particular a
manifold carrying an immersion carries a smooth Riemannian metric, which is all the Mantegazza–Martinazzi theorem
needs as background.
-/

open Set Filter Bundle
open DifferentialGeometry DifferentialGeometry.Geometry
open scoped InnerProductSpace RealInnerProductSpace Topology ContDiff Manifold

noncomputable section

set_option linter.unusedSectionVars false

namespace CWFlow

universe u

variable {E : Type u} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
variable {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
variable {V : Type u} [NormedAddCommGroup V] [InnerProductSpace ℝ V]

variable (E) (F : M → V)

/-- The first fundamental form of `F` at `x`, read in the chart at `x`. -/
def indForm (x : M) : TangentSpace 𝓘(ℝ, E) x →L[ℝ] TangentSpace 𝓘(ℝ, E) x →L[ℝ] ℝ :=
  (show E →L[ℝ] E →L[ℝ] ℝ from
    ((innerSL ℝ : V →L[ℝ] V →L[ℝ] ℝ).bilinearComp (fderiv ℝ (chartRep E F x) (chart E x x))
      (fderiv ℝ (chartRep E F x) (chart E x x))))

variable {E F}

theorem indForm_apply (x : M) (v w : E) :
    (show E →L[ℝ] E →L[ℝ] ℝ from indForm E F x) v w =
      ⟪fderiv ℝ (chartRep E F x) (chart E x x) v, fderiv ℝ (chartRep E F x) (chart E x x) w⟫ := rfl

variable {J : V →L[ℝ] V} (hF : IsLagImmersion E J F)
include hF

/-- The derivative of `F` in the chart at `x`, on a frame vector of the chart at `x₀`, is the derivative of `F` in
the chart at `x₀`. -/
theorem fderiv_chartRep_symmL {x₀ x : M} (hx : x ∈ (chart E x₀).source) (e : E) :
    fderiv ℝ (chartRep E F x) (chart E x x)
        ((trivializationAt E (TangentSpace 𝓘(ℝ, E)) x₀).symmL ℝ x e) =
      fderiv ℝ (chartRep E F x₀) (chart E x₀ x) e := by
  have hs : x ∈ (chartAt E x₀).source := by rwa [extChartAt_source] at hx
  have hxx : x ∈ (chart E x).source := mem_extChartAt_source (I := 𝓘(ℝ, E)) x
  have hy : chart E x₀ x ∈ overlap E x₀ x := chart_mem_overlap hx hxx
  have hτ : transition E x₀ x (chart E x₀ x) = chart E x x := transition_chart hx
  have heq : chartRep E F x₀ =ᶠ[𝓝 (chart E x₀ x)] chartRep E F x ∘ transition E x₀ x := by
    filter_upwards [(isOpen_overlap (E := E) x₀ x).mem_nhds hy] with w hw
    exact chartRep_transition F hw
  have hd : DifferentiableAt ℝ (chartRep E F x) (transition E x₀ x (chart E x₀ x)) :=
    ((hF.smooth x).contDiffAt ((isOpen_chart_target (E := E) x).mem_nhds
      (mapsTo_transition (E := E) x₀ x hy))).differentiableAt (by simp)
  have hdτ : DifferentiableAt ℝ (transition E x₀ x) (chart E x₀ x) :=
    ((contDiffOn_transition (E := E) x₀ x).contDiffAt
      ((isOpen_overlap (E := E) x₀ x).mem_nhds hy)).differentiableAt (by simp)
  rw [heq.fderiv_eq, fderiv_comp _ hd hdτ, hτ, TangentBundle.symmL_trivializationAt_eq_core hs,
    tangentBundleCore_coordChange_achart, ModelWithCorners.range_eq_univ, fderivWithin_univ]
  rfl

/-- **The metric induced by an immersion.** -/
theorem exists_inducedMetric :
    ∃ g : SmoothRiemannianMetric 𝓘(ℝ, E) M, ∀ (x : M) (v w : TangentSpace 𝓘(ℝ, E) x),
      g.inner x v w = indForm E F x v w := by
  refine smoothMetric_of_localCoeff (I := 𝓘(ℝ, E)) (indForm E F) (fun x v w => ?_) (fun x v hv => ?_)
    (fun x₀ i j => ?_)
  · exact real_inner_comm _ _
  · have hne : fderiv ℝ (chartRep E F x) (chart E x x) v ≠ 0 := by
      intro h0
      have hz : fderiv ℝ (chartRep E F x) (chart E x x) (0 : E) = 0 := map_zero _
      exact hv (hF.inj x _ (mem_extChartAt_target (I := 𝓘(ℝ, E)) x) (a₁ := (v : E)) (a₂ := (0 : E))
        (h0.trans hz.symm))
    exact real_inner_self_pos.2 hne
  · have key : ∀ e e' : E, ContMDiffOn 𝓘(ℝ, E) 𝓘(ℝ) ∞
        (fun x => indForm E F x ((trivializationAt E (TangentSpace 𝓘(ℝ, E)) x₀).symmL ℝ x e)
          ((trivializationAt E (TangentSpace 𝓘(ℝ, E)) x₀).symmL ℝ x e'))
        (trivializationAt E (TangentSpace 𝓘(ℝ, E)) x₀).baseSet := by
      intro e e'
      have hD := (hF.smooth x₀).fderiv_of_isOpen (isOpen_chart_target (E := E) x₀) (m := ∞) le_rfl
      have h1 : ContDiffOn ℝ ∞ (fun y => ⟪fderiv ℝ (chartRep E F x₀) y e,
          fderiv ℝ (chartRep E F x₀) y e'⟫) (chart E x₀).target :=
        (hD.clm_apply contDiffOn_const).inner ℝ (hD.clm_apply contDiffOn_const)
      have hbase : (trivializationAt E (TangentSpace 𝓘(ℝ, E)) x₀).baseSet = (chartAt E x₀).source :=
        TangentBundle.trivializationAt_baseSet x₀
      rw [hbase]
      have h2 : ContMDiffOn 𝓘(ℝ, E) 𝓘(ℝ) ∞
          ((fun y => ⟪fderiv ℝ (chartRep E F x₀) y e, fderiv ℝ (chartRep E F x₀) y e'⟫) ∘
            chart E x₀) (chartAt E x₀).source :=
        h1.contMDiffOn.comp (contMDiffOn_extChartAt (I := 𝓘(ℝ, E)) (n := ∞) (x := x₀))
          fun x hx => (chart E x₀).map_source (by rwa [extChartAt_source])
      refine h2.congr fun x hx => ?_
      have hx' : x ∈ (chart E x₀).source := by rwa [extChartAt_source]
      show ⟪fderiv ℝ (chartRep E F x) (chart E x x) _, fderiv ℝ (chartRep E F x) (chart E x x) _⟫ = _
      rw [fderiv_chartRep_symmL hF hx', fderiv_chartRep_symmL hF hx']
      rfl
    exact key _ _

end CWFlow
