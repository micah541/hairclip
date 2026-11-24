import CWFlow.Extension.Cover

/-!
# Graph charts

`IsGraphChart ι J F c₀ p U σ φ` packages the conclusion of `graph_at`: `σ` parametrizes an open piece of `M`
around `p` by the ball of radius 6, on which `F ∘ σ` is the graph `x ↦ F p + U x + J U ∇φ(x)` over the Lagrangian
plane of `U`, with `|D²φ| ≤ c₀`, `|∇φ(x)| ≤ c₀|x|`, and `σ(B_ρ)` is the connected component of `p` in
`{|horC F| < ρ, |verC F| < 1}`. `exists_graphChart_cover` restates `exists_cover` in these terms.
-/

open Metric Set MeasureTheory
open CWFlow.Interp
open scoped ContDiff Manifold Topology RealInnerProductSpace

namespace CWFlow.Extension

variable {n : ℕ} {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
  {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
  {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
  (ι : EuclideanSpace ℝ (Fin n) ≃L[ℝ] E)

/-- A graph chart of `F` of slope `c₀` centred at `p`. -/
structure IsGraphChart (J : V →L[ℝ] V) (F : M → V) (c₀ : ℝ) (p : M)
    (U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V) (σ : EuclideanSpace ℝ (Fin n) → M)
    (φ : EuclideanSpace ℝ (Fin n) → ℝ) : Prop where
  frame : IsLagFrame J U
  param : IsParam ι (ball 0 6) σ
  center : σ 0 = p
  smooth : ContDiffOn ℝ ∞ φ (ball 0 6)
  zero : φ 0 = 0
  graph : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6,
    F (σ x) = graphX J U (F p) φ x ∧ horC U (F p) (F (σ x)) = x ∧
    ‖jet 2 φ x‖ ≤ c₀ ∧ ‖gradE φ x‖ ≤ c₀ * ‖x‖
  comp : ∀ ρ, 0 < ρ → ρ ≤ 6 → σ '' ball 0 ρ = connectedComponentIn
    {q | ‖horC U (F p) (F q)‖ < ρ ∧ ‖verC J U (F p) (F q)‖ < 1} p

/-- `exists_cover` in terms of graph charts. -/
theorem exists_graphChart_cover [T2Space M] [CompactSpace M] (A : IntAtlas (E := E) M) {J : V →L[ℝ] V}
    (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    [FiniteDimensional ℝ V] (hdim : Module.finrank ℝ V = 2 * n) {F : M → V} (hF : IsLagImm ι J F)
    {c₀ c₁ : ℝ} (hc₀ : 0 < c₀) (hc₀1 : 6 * c₀ < 1) (hc₁ : 0 ≤ c₁)
    (hsmall : 6 * (c₁ * (1 + c₀ ^ 2) ^ ((3 : ℝ) / 2)) < c₀)
    (hA : ∀ z : M, ∀ y ∈ chartDomE ι z, aNormSq J (F ∘ chartInv ι z) y ≤ c₁ ^ 2) :
    ∃ (N : ℕ) (p : Fin N → M) (U : Fin N → (EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V))
      (σ : Fin N → EuclideanSpace ℝ (Fin n) → M) (φ : Fin N → EuclideanSpace ℝ (Fin n) → ℝ),
      (N : ℝ) * (volume (ball (0 : EuclideanSpace ℝ (Fin n)) (1 / 16))).toReal ≤ volM ι A F ∧
      (∀ j, IsGraphChart ι J F c₀ (p j) (U j) (σ j) (φ j)) ∧
      ∀ q : M, ∃ j, q ∈ σ j '' ball 0 1 := by
  obtain ⟨N, p, U, σ, φ, hvol, hch, hcov⟩ :=
    exists_cover ι A hJi hJJ hdim hF hc₀ hc₀1 hc₁ hsmall hA
  refine ⟨N, p, U, σ, φ, hvol, fun j => ?_, hcov⟩
  obtain ⟨h1, h2, h3, h4, h5, h6, h7⟩ := hch j
  exact ⟨h1, h2, h3, h4, h5, h6, h7⟩

end CWFlow.Extension
