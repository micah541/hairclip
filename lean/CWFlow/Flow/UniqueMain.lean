import CWFlow.Flow.Main
import CWFlow.Flow.FlowUnique

/-!
# Uniqueness of the flow of Lagrangian immersions: final statement

`flow_unique_of_immersion`: on a compact manifold, two solutions (`IsFlow`) of `∂ₜF = J∇div(JH)` on `[0, T]`
with the same initial immersion coincide. No auxiliary Riemannian metric is assumed (it is the metric induced
by the initial immersion, `exists_inducedMetric`).
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

/-- **Uniqueness of the flow `∂ₜF = J∇div(JH)` of Lagrangian immersions of a compact manifold.** -/
theorem flow_unique_of_immersion {J : V →L[ℝ] V}
    (hJi : ∀ u v, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v, J (J v) = -v)
    (hdim : Module.finrank ℝ V = 2 * Module.finrank ℝ E) {T : ℝ} {Φ₁ Φ₂ : ℝ → M → V}
    (hΦ₁ : IsFlow E J T Φ₁) (hΦ₂ : IsFlow E J T Φ₂) (h0 : Φ₁ 0 = Φ₂ 0) :
    ∀ t ∈ Icc (0 : ℝ) T, Φ₁ t = Φ₂ t := by
  intro t ht
  haveI : CompleteSpace E := FiniteDimensional.complete ℝ E
  obtain ⟨g, -⟩ := exists_inducedMetric (hΦ₁.imm 0 ⟨le_rfl, ht.1.trans ht.2⟩)
  exact Bridge.flow_unique g hJi hJJ hdim hΦ₁ hΦ₂ h0 t ht

/-- info: 'CWFlow.flow_unique_of_immersion' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms flow_unique_of_immersion

end CWFlow
