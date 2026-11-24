import Mathlib.Analysis.Calculus.InverseFunctionTheorem.ContDiff
import Mathlib.Analysis.Calculus.ContDiff.Operations
import Mathlib.Analysis.Normed.Module.FiniteDimension

/-!
# Local inverses of a smooth map with injective derivative between spaces of equal dimension

`exists_openPartialHomeomorph_of_injective_fderiv`: around every point of an open set `Q` on which `f` is smooth
with injective derivative, `f` restricts to an open partial homeomorphism with source inside `Q` and smooth
inverse. Used for the Weinstein maps `chartPsi` (lifting a flow to cotangent coordinates).
-/

open Set Filter

open scoped Topology ContDiff

namespace CWFlow

variable {W V : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W] [FiniteDimensional ℝ W]
  [NormedAddCommGroup V] [NormedSpace ℝ V] [FiniteDimensional ℝ V]

/-- The derivative of `f` at a point where it is injective, as a continuous linear equivalence. -/
theorem exists_equiv_fderiv_of_injective (hdim : Module.finrank ℝ W = Module.finrank ℝ V) {f : W → V}
    {q : W} (hinj : Function.Injective (fderiv ℝ f q)) :
    ∃ L : W ≃L[ℝ] V, (L : W →L[ℝ] V) = fderiv ℝ f q :=
  ⟨(LinearMap.linearEquivOfInjective (fderiv ℝ f q : W →ₗ[ℝ] V) hinj hdim).toContinuousLinearEquiv,
    ContinuousLinearMap.ext fun _ => rfl⟩

theorem exists_openPartialHomeomorph_of_injective_fderiv
    (hdim : Module.finrank ℝ W = Module.finrank ℝ V) {f : W → V} {Q : Set W} (hQ : IsOpen Q)
    (hf : ContDiffOn ℝ ∞ f Q) (hinj : ∀ q ∈ Q, Function.Injective (fderiv ℝ f q)) {q₀ : W}
    (hq₀ : q₀ ∈ Q) :
    ∃ e : OpenPartialHomeomorph W V, q₀ ∈ e.source ∧ e.source ⊆ Q ∧ (∀ q, e q = f q) ∧
      ContDiffOn ℝ ∞ e.symm e.target := by
  have hW : CompleteSpace W := FiniteDimensional.complete ℝ W
  have hcd : ∀ q ∈ Q, ContDiffAt ℝ ∞ f q := fun q hq => hf.contDiffAt (hQ.mem_nhds hq)
  have hder : ∀ q ∈ Q, ∃ L : W ≃L[ℝ] V, HasFDerivAt f (L : W →L[ℝ] V) q := by
    intro q hq
    obtain ⟨L, hL⟩ := exists_equiv_fderiv_of_injective hdim (hinj q hq)
    refine ⟨L, ?_⟩
    rw [hL]
    exact ((hcd q hq).differentiableAt (by simp)).hasFDerivAt
  obtain ⟨L, hL⟩ := hder q₀ hq₀
  refine ⟨((hcd q₀ hq₀).toOpenPartialHomeomorph f hL (by simp)).restrOpen Q hQ,
    ⟨(hcd q₀ hq₀).mem_toOpenPartialHomeomorph_source hL (by simp), hq₀⟩, fun q hq => hq.2,
    fun q => rfl, ?_⟩
  intro b hb
  set e := ((hcd q₀ hq₀).toOpenPartialHomeomorph f hL (by simp)).restrOpen Q hQ with he
  have ha : e.symm b ∈ Q := (e.map_target hb).2
  obtain ⟨La, hLa⟩ := hder _ ha
  exact (e.contDiffAt_symm hb hLa (hcd _ ha)).contDiffWithinAt

end CWFlow
