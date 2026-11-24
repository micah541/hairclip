import CWFlow.Flow.Velocity

/-!
# The family of immersions of a scalar solution is jointly smooth

`famRep g J F u α (t, y)` is the family `t ↦ graphImm (fnOf g (u t))` in the chart at `α`. For a solution `u` of
the scalar equation it is smooth on `Icc 0 T ×ˢ (chart E α).target` (`contDiffOn_famRep`): near each point it is
the chart graph of a good patch composed with a transition.
-/

open Set Filter Bundle
open DifferentialGeometry DifferentialGeometry.Tensor0SBundle ShortTime ShortTime.HFlow CWFlow.JetRec
open DifferentialGeometry.Geometry.Curvature DifferentialGeometry.Geometry.Connection
open DifferentialGeometry.Integral.L2
open DifferentialGeometry.Analysis.Sobolev DifferentialGeometry.Analysis.Spectral
open DifferentialGeometry.Analysis.Parabolic.TensorSpectral DifferentialGeometry.Analysis.Elliptic
open scoped InnerProductSpace RealInnerProductSpace Topology ContDiff Manifold

noncomputable section

set_option linter.unusedSectionVars false

namespace CWFlow.Bridge

universe u

variable {E : Type u} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
  [NeZero (Module.finrank ℝ E)] [CompleteSpace E]
variable {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
  [CompactSpace M] [T2Space M] [SigmaCompactSpace M]
variable {V : Type u} [NormedAddCommGroup V] [InnerProductSpace ℝ V]

variable (g : SmoothRiemannianMetric 𝓘(ℝ, E) M) (J : V →L[ℝ] V) (F : M → V)

/-- The family of immersions `t ↦ graphImm (fnOf g (u t))` in the chart at `α`. -/
def famRep (u : ℝ → SmoothCcTensor g 0 0) (α : M) (q : ℝ × E) : V :=
  chartRep E (graphImm E J F (fnOf g (u q.1))) α q.2

variable {J F}
variable (hJi : ∀ u v, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v, J (J v) = -v) (hF : IsLagImmersion E J F)
variable {u₀ : SmoothCcTensor g 0 0} {T : ℝ} {u : ℝ → SmoothCcTensor g 0 0}
variable (hu : (cwData g hJi hJJ hF).toEquationData.IsSolution g u₀ T u)
include hu

/-- **The family of immersions of a solution is jointly smooth in each chart.** -/
theorem contDiffOn_famRep (α : M) :
    ContDiffOn ℝ ∞ (famRep g J F u α) (Icc 0 T ×ˢ (chart E α).target) := by
  rintro ⟨t, y⟩ ⟨ht, hy⟩
  obtain ⟨i, hi, hgi⟩ := exists_good g hJi hJJ hF hu ht ((chart E α).symm y)
  have hyβ := mem_graphDom_of_good g J F (u t) hi hgi
  have hov : y ∈ overlap E α i.1 := ⟨hy, hi.1⟩
  have hwβ := contDiffOn_wOf g hu.smooth i.1
  have hΨ : ContDiffAt ℝ ∞ (chartPsi E J F i.1)
      (graphMap (fun z => wOf g u i.1 (t, z)) (transition E α i.1 y)) :=
    (contDiffOn_chartPsi hJi hJJ hF i.1).contDiffAt ((isOpen_chartDom hJi hJJ hF i.1).mem_nhds hyβ.2)
  have hK0 := contDiffWithinAt_graphFu_time hu.pos (isOpen_chart_target (E := E) i.1) hwβ ht hyβ.1 hΨ
  have hφ : ContDiffOn ℝ ∞ (fun p : ℝ × E => (p.1, transition E α i.1 p.2))
      (Icc 0 T ×ˢ overlap E α i.1) :=
    contDiffOn_fst.prodMk ((contDiffOn_transition (E := E) α i.1).comp contDiffOn_snd fun p hp => hp.2)
  have hmaps : MapsTo (fun p : ℝ × E => (p.1, transition E α i.1 p.2))
      (Icc 0 T ×ˢ overlap E α i.1) (Icc 0 T ×ˢ (chart E i.1).target) :=
    fun p hp => ⟨hp.1, mapsTo_transition (E := E) α i.1 hp.2⟩
  have hK' := ContDiffWithinAt.comp (g := fun p : ℝ × E =>
      graphFu (chartPsi E J F i.1) (fun z => wOf g u i.1 (p.1, z)) p.2)
    (f := fun p : ℝ × E => (p.1, transition E α i.1 p.2)) (t, y) hK0 (hφ (t, y) ⟨ht, hov⟩) hmaps
  have hK : ContDiffWithinAt ℝ ∞
      (fun p : ℝ × E => graphFu (chartPsi E J F i.1) (fun z => wOf g u i.1 (p.1, z))
        (transition E α i.1 p.2)) (Icc 0 T ×ˢ overlap E α i.1) (t, y) := hK'
  have hc1 : ContinuousWithinAt
      (fun p : ℝ × E => fderiv ℝ (fun z => wOf g u i.1 (p.1, z)) (transition E α i.1 p.2))
      (Icc 0 T ×ˢ overlap E α i.1) (t, y) :=
    ContinuousWithinAt.comp (g := fun p : ℝ × E => fderiv ℝ (fun z => wOf g u i.1 (p.1, z)) p.2)
      (f := fun p : ℝ × E => (p.1, transition E α i.1 p.2))
      ((contDiffOn_slice_fderiv hu.pos (isOpen_chart_target (E := E) i.1) hwβ).continuousOn
        (t, transition E α i.1 y) ⟨ht, hyβ.1⟩) (hφ.continuousOn (t, y) ⟨ht, hov⟩) hmaps
  have hc2 : ContinuousWithinAt (fun p : ℝ × E => transition E α i.1 p.2)
      (Icc 0 T ×ˢ overlap E α i.1) (t, y) :=
    ((contDiffOn_transition (E := E) α i.1).continuousOn.comp continuousOn_snd fun p hp => hp.2)
      (t, y) ⟨ht, hov⟩
  have hgm : ContinuousWithinAt
      (fun p : ℝ × E => graphMap (fun z => wOf g u i.1 (p.1, z)) (transition E α i.1 p.2))
      (Icc 0 T ×ˢ overlap E α i.1) (t, y) := hc2.prodMk hc1
  have hev : famRep g J F u α =ᶠ[𝓝[Icc 0 T ×ˢ overlap E α i.1] (t, y)]
      fun p : ℝ × E => graphFu (chartPsi E J F i.1) (fun z => wOf g u i.1 (p.1, z))
        (transition E α i.1 p.2) := by
    filter_upwards [hgm.preimage_mem_nhdsWithin ((isOpen_chartDom hJi hJJ hF i.1).mem_nhds hyβ.2),
      self_mem_nhdsWithin] with p hp hpS
    exact graphImm_eq (contDiffOn_chartRep_fnOf g (u p.1)) hJi hJJ hF
      ⟨hpS.2.2, mapsTo_transition (E := E) α i.1 hpS.2, hp⟩
  refine (hK.congr_of_eventuallyEq hev (hev.self_of_nhdsWithin ⟨ht, hov⟩)).mono_of_mem_nhdsWithin ?_
  exact mem_nhdsWithin.2 ⟨univ ×ˢ overlap E α i.1, isOpen_univ.prod (isOpen_overlap (E := E) α i.1),
    ⟨mem_univ _, hov⟩, fun p hp => ⟨hp.2.1, hp.1.2⟩⟩

end CWFlow.Bridge
