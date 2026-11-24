import CWFlow.Bridge.Smooth

/-!
# The scalar equation as Mantegazza–Martinazzi data

`cwData g hJi hJJ hF : ShortTime.MMData 𝓘(ℝ, E) M 2` is the scalar equation `∂ₜu = −div(JH)[F_u]` of a Lagrangian
immersion `F : M → V` of a closed manifold, written in the form of Mantegazza–Martinazzi:
`∂ₜu = −E^{ij}E^{kl}∇⁴_{ijkl}u + b`, with `E` the glued cometric of the graph and `b` of lower order. The fields are
proved in `Bridge/Fields.lean` and `Bridge/Smooth.lean`.
-/

open Set Filter Bundle
open DifferentialGeometry DifferentialGeometry.Tensor0SBundle ShortTime ShortTime.HFlow CWFlow.JetRec
open DifferentialGeometry.Geometry.Curvature DifferentialGeometry.Geometry.Connection
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

variable (g : SmoothRiemannianMetric 𝓘(ℝ, E) M) {J : V →L[ℝ] V} {F : M → V}
variable (hJi : ∀ u v, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v, J (J v) = -v) (hF : IsLagImmersion E J F)

/-- **The scalar equation of the flow as Mantegazza–Martinazzi data.** -/
def cwData : MMData 𝓘(ℝ, E) M 2 where
  T₀ := 1
  T₀_pos := one_pos
  dom := mmDom g J F
  Ef := fun _ => mmEf g J F
  b := mmb g J F
  dom_relOpen := mmDom_relOpen hJi hJJ hF 1
  Ef_lower := fun _ => mmEf_lower g J F
  b_lower := mmb_lower g J F
  Ef_symm := fun _ => mmEf_symm g J F
  Ef_pos := fun _ _ _ _ _ hζ _ hξ => mmEf_pos hJi hJJ hF hζ hξ
  Ef_smooth := fun x₀ _ i j => contMDiffOn_mmEf hJi hJJ hF 1 x₀ i j
  b_smooth := contMDiffOn_mmb hJi hJJ hF 1

/-- info: 'CWFlow.Bridge.cwData' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms cwData

end CWFlow.Bridge
