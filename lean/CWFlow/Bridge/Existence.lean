import CWFlow.Bridge.Realizable
import ShortTime.MM.Main

/-!
# Short-time existence and uniqueness for the scalar equation

For a Lagrangian immersion `F : M → V` of a closed manifold, the scalar equation `∂ₜu = −div(JH)[F_u]` with
initial datum `u = 0` (the immersion `F` itself) has a smooth solution for a short time, and solutions are unique
(`exists_scalar_solution`, `scalar_solution_unique`): the Mantegazza–Martinazzi theorem applied to `cwData`.

`hasDerivWithinAt_fnOf` says what a solution of the Mantegazza–Martinazzi equation is, in terms of the flow:
`∂ₜ u(t, x) = −div(JH)[F_{u(t)}](x)`, computed in any chart around `x` in whose Weinstein domain the jet lies.
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

variable (g : SmoothRiemannianMetric 𝓘(ℝ, E) M)

theorem fnOf_zero : fnOf g (0 : SmoothCcTensor g 0 0) = 0 := by
  funext x
  unfold fnOf
  rw [toSection_zero_apply, map_zero]
  rfl

/-- The scalar value of a `(0,0)`-tensor, through its model. -/
theorem scal_toModel (x : M) (T : TensorRSSpace 0 0 𝓘(ℝ, E) x) :
    scal E (TensorRSSpace.toModel T) = lam (I := 𝓘(ℝ, E)) x T := by
  rw [scal_apply, ← toModel_unitZeroSec x, TensorRSSpace.toModel_apply_toModel]
  rfl

variable {J : V →L[ℝ] V} {F : M → V}
variable (hJi : ∀ u v, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v, J (J v) = -v) (hF : IsLagImmersion E J F)
include hJi hJJ hF

/-- The zero function (the immersion `F` itself) lies in the domain of the equation. -/
theorem covJet_zero_mem_mmDom (t : ℝ) (x : M) :
    covJet g 0 0 (2 * 2) (0 : SmoothCcTensor g 0 0) x ∈ mmDom g J F t x := by
  intro i hx
  have hi := tsupport_patchρ_subset i hx
  unfold good
  rw [cJet_covJet g 0 hi]
  have h0 : (jetOf (chartRep E (fnOf g (0 : SmoothCcTensor g 0 0)) i.1) (chart E i.1 x)).1 = 0 := by
    show fderiv ℝ (chartRep E (fnOf g (0 : SmoothCcTensor g 0 0)) i.1) (chart E i.1 x) = 0
    rw [fnOf_zero]
    exact fderiv_const_apply (0 : ℝ)
  rw [h0]
  exact ((weinstein_global hJi hJJ hF).1 i.1).2 _ ((chart E i.1).map_source hi.1)

/-- **Short-time existence for the scalar equation**, starting from `u = 0`. -/
theorem exists_scalar_solution :
    ∃ T > 0, ∃ u : ℝ → SmoothCcTensor g 0 0,
      (cwData g hJi hJJ hF).toEquationData.IsSolution g 0 T u :=
  (mantegazza_martinazzi (cwData g hJi hJJ hF) (by norm_num) g).1 0
    fun x => covJet_zero_mem_mmDom g hJi hJJ hF 0 x

/-- **Uniqueness for the scalar equation.** -/
theorem scalar_solution_unique (u₀ : SmoothCcTensor g 0 0) (T : ℝ) (u v : ℝ → SmoothCcTensor g 0 0)
    (hu : (cwData g hJi hJJ hF).toEquationData.IsSolution g u₀ T u)
    (hv : (cwData g hJi hJJ hF).toEquationData.IsSolution g u₀ T v) :
    ∀ t ∈ Icc (0 : ℝ) T, u t = v t :=
  (mantegazza_martinazzi (cwData g hJi hJJ hF) (by norm_num) g).2.1 u₀ T u v hu hv

/-- **A solution of the Mantegazza–Martinazzi equation solves the scalar equation of the flow:**
`∂ₜ u(t, x) = −div(JH)[F_{u(t)}](x)`, in any chart `α` around `x` in whose Weinstein domain the jet lies. -/
theorem hasDerivWithinAt_fnOf {u₀ : SmoothCcTensor g 0 0} {T : ℝ} {u : ℝ → SmoothCcTensor g 0 0}
    (hu : (cwData g hJi hJJ hF).toEquationData.IsSolution g u₀ T u) {t : ℝ} (ht : t ∈ Icc (0 : ℝ) T)
    {α x : M} (hα : x ∈ patchU E α) (hgα : good g J F α x (covJet g 0 0 (2 * 2) (u t) x)) :
    HasDerivWithinAt (fun τ => fnOf g (u τ) x)
      (-flowPot J (chartGraph E J F (fnOf g (u t)) α) (chart E α x)) (Icc 0 T) t := by
  have h1 := (scal E).hasFDerivAt.comp_hasDerivWithinAt t (hu.eqn t ht x)
  have h2 : (scal E ∘ fun τ => TensorRSSpace.toModel ((u τ).toSection x)) =
      fun τ => fnOf g (u τ) x := funext fun τ => scal_toModel x _
  rw [h2] at h1
  have h3 := lam_cwData_F hJi hJJ hF (u t) t (hu.jets_mem t ht x) hα hgα
  rw [← h3, ← scal_toModel]
  exact h1

/-- info: 'CWFlow.Bridge.exists_scalar_solution' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_scalar_solution

/-- info: 'CWFlow.Bridge.hasDerivWithinAt_fnOf' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms hasDerivWithinAt_fnOf

end CWFlow.Bridge
