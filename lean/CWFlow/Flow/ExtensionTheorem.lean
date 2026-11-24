import CWFlow.Extension.FlowClosed
import CWFlow.Flow.Restart

/-!
# The extension theorem in `ℂⁿ` (Chen–Warren, Thm. 4.1)

A flow `∂ₜF = J∇div(JH)` of Lagrangian immersions of a closed manifold into a Hermitian vector space, defined on
`[0, T)` with bounded second fundamental form, extends to a flow on `[0, T + ε]` for some `ε > 0`.

The proof combines blueprint step 4 (`Extension.isFlow_closed_of_bound`: the flow extends smoothly to `[0, T]`,
from the energy estimates of Chen–Warren §4) with step 5 (`isFlow_restart`: a flow on `[0, T]` can be continued,
through the gauge theorem for the reversed path, continuous dependence for the scalar equation and uniqueness).

The bound on the second fundamental form is stated in the charts `chartInv ι z` determined by a fixed linear
isomorphism `ι : ℝⁿ ≃ E`; `aNormSq` is the squared norm `|A|²_g` read in a chart, and does not depend on the
chart (`aNormSq_comp`).
-/

open Set
open scoped ContDiff Manifold Topology RealInnerProductSpace

namespace CWFlow

universe u

variable {E : Type u} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
  [NeZero (Module.finrank ℝ E)]
variable {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
  [CompactSpace M] [T2Space M] [SigmaCompactSpace M]
variable {V : Type u} [NormedAddCommGroup V] [InnerProductSpace ℝ V] [FiniteDimensional ℝ V]

/-- **The extension theorem** (Chen–Warren, Thm. 4.1, in `ℂⁿ`): a flow on `[0, T)` (an `IsFlow` on every
`[0, T']`, `T' < T`) whose second fundamental form is bounded extends to a flow on `[0, T + ε]`. -/
theorem flow_extends_of_bounded_A {n : ℕ} (ι : EuclideanSpace ℝ (Fin n) ≃L[ℝ] E) {J : V →L[ℝ] V}
    (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    (hdim : Module.finrank ℝ V = 2 * n)
    {T : ℝ} (hT : 0 < T) {Φ : ℝ → M → V} (hΦ : ∀ T' ∈ Ioo 0 T, IsFlow E J T' Φ) {K : ℝ}
    (hK : ∀ t ∈ Ico 0 T, ∀ z : M, ∀ y ∈ Extension.chartDomE ι z,
      Extension.aNormSq J (Φ t ∘ Extension.chartInv ι z) y ≤ K ^ 2) :
    ∃ ε > 0, ∃ Ψ : ℝ → M → V, IsFlow E J (T + ε) Ψ ∧ ∀ t ∈ Ico 0 T, Ψ t = Φ t := by
  obtain ⟨Φbar, hΦbar, hagree⟩ := Extension.isFlow_closed_of_bound ι hJi hJJ hdim hT hΦ hK
  have hdim' : Module.finrank ℝ V = 2 * Module.finrank ℝ E := by
    rw [hdim, ← ι.toLinearEquiv.finrank_eq, finrank_euclideanSpace_fin]
  obtain ⟨ε, hε, Ψ, hΨ, hΨeq⟩ := isFlow_restart hJi hJJ hdim' hT hΦbar
  exact ⟨ε, hε, Ψ, hΨ, fun t ht => by rw [hΨeq t (Ico_subset_Icc_self ht), hagree t ht]⟩

/-- info: 'CWFlow.flow_extends_of_bounded_A' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms flow_extends_of_bounded_A

end CWFlow
