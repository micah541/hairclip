import CWFlow.Extension.IntM

/-!
# Flows on `M`, read in charts (light interface)

`IsFlowM ι J I Φ`: on the open time set `I`, `Φ` is jointly smooth in every chart, every `Φ t` is a Lagrangian
immersion (`IsLagImm`), and the velocity in every chart is the coordinate velocity `velX` (`= J∇div(JH)`).
This is the interface used by the extension theorem; the existing `CWFlow.IsFlow` (closed intervals, basis-free
`flowVel`) implies it on the interior through `flowVel_eq_velX` (heavy file `ImmPot`).

`volM` and `energyM` are the volume and the energy `E_k = ∫ |∇^k C|²_g dV_g` over `M`, defined with a fixed
`IntAtlas`.
-/

open scoped ContDiff Manifold Topology

namespace CWFlow.Extension

variable {n : ℕ} {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
  {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
  {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
  (ι : EuclideanSpace ℝ (Fin n) ≃L[ℝ] E)

/-- A flow of Lagrangian immersions on the open time set `I`, read in the charts of `M`. -/
structure IsFlowM (J : V →L[ℝ] V) (I : Set ℝ) (Φ : ℝ → M → V) : Prop where
  isOpen : IsOpen I
  smooth : ∀ z : M, ContDiffOn ℝ ∞ (fun q : ℝ × EuclideanSpace ℝ (Fin n) => Φ q.1 (chartInv ι z q.2))
    (I ×ˢ chartDomE ι z)
  imm : ∀ t ∈ I, IsLagImm ι J (Φ t)
  vel : ∀ t ∈ I, ∀ z : M, ∀ y ∈ chartDomE ι z,
    HasDerivAt (fun s => Φ s (chartInv ι z y)) (velX J (Φ t ∘ chartInv ι z) y) t

/-- The volume of `F`. -/
noncomputable def volM (A : IntAtlas (E := E) M) (F : M → V) : ℝ :=
  intM ι A F (fun _ => 1) (densOf volX)

/-- The energy `E_k(F) = ∫_M |∇^k C|²_g dV_g`. -/
noncomputable def energyM (A : IntAtlas (E := E) M) (J : V →L[ℝ] V) (k : ℕ) (F : M → V) : ℝ :=
  intM ι A F (fun _ => 1) (densOf fun X => eDens J X k)

end CWFlow.Extension
