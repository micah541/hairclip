import CWFlow.Flow.FlowExistence
import CWFlow.Flow.InducedMetric

/-!
# Short-time existence for the Hamiltonian gradient flow of volume

`exists_flow_of_immersion`: for a Lagrangian immersion `F : M → V` of a closed manifold `M` (`dim V = 2 dim M`),
there are `T > 0` and a family `Φ t`, `t ∈ [0, T]`, of Lagrangian immersions, jointly smooth in `(t, x)`, with
`Φ 0 = F` and `∂ₜΦ = J∇div(JH)[Φ_t]`. No background metric is assumed: the Mantegazza–Martinazzi theorem is run
with the metric induced by `F`.
-/

open Set Filter
open scoped InnerProductSpace RealInnerProductSpace Topology ContDiff Manifold

noncomputable section

namespace CWFlow

universe u

variable {E : Type u} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
  [NeZero (Module.finrank ℝ E)]
variable {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
  [CompactSpace M] [T2Space M] [SigmaCompactSpace M]
variable {V : Type u} [NormedAddCommGroup V] [InnerProductSpace ℝ V] [FiniteDimensional ℝ V]

/-- **Short-time existence for the flow `∂ₜF = J∇div(JH)`** of a Lagrangian immersion of a closed manifold. -/
theorem exists_flow_of_immersion {J : V →L[ℝ] V} {F : M → V}
    (hJi : ∀ u v, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v, J (J v) = -v) (hF : IsLagImmersion E J F)
    (hdim : Module.finrank ℝ V = 2 * Module.finrank ℝ E) :
    ∃ T > 0, ∃ Φ : ℝ → M → V, Φ 0 = F ∧
      ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, V) ∞ (fun q : ℝ × M => Φ q.1 q.2) (Icc 0 T ×ˢ univ) ∧
      (∀ t ∈ Icc (0 : ℝ) T, IsLagImmersion E J (Φ t)) ∧
      ∀ t ∈ Icc (0 : ℝ) T, ∀ x z : M, x ∈ (chart E z).source →
        HasDerivWithinAt (fun s => Φ s x) (flowVel J (chartRep E (Φ t) z) (chart E z x)) (Icc 0 T) t := by
  haveI : CompleteSpace E := FiniteDimensional.complete ℝ E
  obtain ⟨g, -⟩ := exists_inducedMetric hF
  exact Bridge.exists_flow g hJi hJJ hF hdim

/-- info: 'CWFlow.exists_flow_of_immersion' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_flow_of_immersion

end CWFlow
