import CWFlow.Flow.GraphUnique

/-!
# Flows with velocity `κ • J∇div(JH)`

`IsFlowS E J κ T Φ` is `IsFlow E J T Φ` with velocity `κ • flowVel`. For `κ = 1` it is `IsFlow`; for `κ = -1`
it describes the time reversal of a flow. The gauge theorem is proved for this generality.
-/

open Set Filter
open scoped InnerProductSpace RealInnerProductSpace Topology ContDiff Manifold

noncomputable section

set_option linter.unusedSectionVars false

namespace CWFlow

section Defs

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
variable {M : Type*} [TopologicalSpace M] [ChartedSpace E M]
variable {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]

variable (E) in
/-- **A solution of `∂ₜF = κ J∇div(JH)` on `[0, T]`**: as `IsFlow`, with velocity `κ • flowVel`. -/
structure IsFlowS (J : V →L[ℝ] V) (κ : ℝ) (T : ℝ) (Φ : ℝ → M → V) : Prop where
  smooth : ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, V) ∞ (fun q : ℝ × M => Φ q.1 q.2) (Icc 0 T ×ˢ univ)
  imm : ∀ t ∈ Icc (0 : ℝ) T, IsLagImmersion E J (Φ t)
  vel : ∀ t ∈ Icc (0 : ℝ) T, ∀ x z : M, x ∈ (chart E z).source →
    HasDerivWithinAt (fun s => Φ s x) (κ • flowVel J (chartRep E (Φ t) z) (chart E z x)) (Icc 0 T) t

/-- A flow is an `IsFlowS` with `κ = 1`. -/
theorem IsFlow.toS {J : V →L[ℝ] V} {T : ℝ} {Φ : ℝ → M → V} (h : IsFlow E J T Φ) :
    IsFlowS E J 1 T Φ where
  smooth := h.smooth
  imm := h.imm
  vel := fun t ht x z hx => by simpa only [one_smul] using h.vel t ht x z hx

/-- An `IsFlowS` with `κ = 1` is a flow. -/
theorem IsFlowS.toIsFlow {J : V →L[ℝ] V} {T : ℝ} {Φ : ℝ → M → V} (h : IsFlowS E J 1 T Φ) :
    IsFlow E J T Φ where
  smooth := h.smooth
  imm := h.imm
  vel := fun t ht x z hx => by simpa only [one_smul] using h.vel t ht x z hx

/-- An `IsFlowS` on `[0, T]` is one on `[0, δ]`, `δ ≤ T`. -/
theorem IsFlowS.mono {J : V →L[ℝ] V} {κ T δ : ℝ} {Φ : ℝ → M → V} (h : IsFlowS E J κ T Φ)
    (hδ : δ ≤ T) : IsFlowS E J κ δ Φ where
  smooth := h.smooth.mono (prod_mono (Icc_subset_Icc_right hδ) subset_rfl)
  imm := fun t ht => h.imm t ⟨ht.1, ht.2.trans hδ⟩
  vel := fun t ht x z hx => (h.vel t ⟨ht.1, ht.2.trans hδ⟩ x z hx).mono (Icc_subset_Icc_right hδ)

end Defs

end CWFlow
