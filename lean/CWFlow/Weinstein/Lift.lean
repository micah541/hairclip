import CWFlow.Weinstein.ZeroSection

/-!
# The metric lift `T*L → TL → ℂⁿ`

For an immersion `F : E → V` (in a chart, `E` finite-dimensional), `g_x(u, w) = ⟨DF u, DF w⟩` and
`metricMap F x : E →L E* `, `u ↦ g_x(u, ·)`. Where `DF x` is injective, `g_x` is an isomorphism, and the lift
`liftA F x = DF x ∘ g_x⁻¹ : E* →L V` satisfies `⟨liftA ξ, DF u⟩ = ξ u` (`inner_liftA`), so that `λ_A = λ_can`
(`lamA_liftA`); on a Lagrangian, `⟨J liftA ξ, liftA η⟩ = 0` (`inner_J_liftA`). `liftA` is smooth wherever `F` is
smooth and `DF` injective (`contDiffAt_liftA`).
-/

open Filter

open scoped InnerProductSpace RealInnerProductSpace Topology ContDiff

noncomputable section

namespace CWFlow

variable {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
variable {F : E → V}

/-- `u ↦ g(u, ·)` with `g(u, w) = ⟨DF u, DF w⟩`. -/
def metricMap (F : E → V) (x : E) : E →L[ℝ] (E →L[ℝ] ℝ) :=
  ((ContinuousLinearMap.compL ℝ E V ℝ).flip (fderiv ℝ F x)).comp
    ((innerSL ℝ : V →L[ℝ] V →L[ℝ] ℝ).comp (fderiv ℝ F x))

@[simp] theorem metricMap_apply (x u w : E) : metricMap F x u w = ⟪fderiv ℝ F x u, fderiv ℝ F x w⟫ :=
  rfl

/-- The metric lift `ξ ↦ DF(g⁻¹ ξ)`. -/
def liftA (F : E → V) (x : E) : (E →L[ℝ] ℝ) →L[ℝ] V :=
  (fderiv ℝ F x).comp (ContinuousLinearMap.inverse (metricMap F x))

theorem metricMap_injective {x : E} (hinj : Function.Injective (fderiv ℝ F x)) :
    Function.Injective (metricMap F x) := by
  rw [injective_iff_map_eq_zero]
  intro u hu
  have h := congrArg (fun f => f u) hu
  simp only [metricMap_apply, ContinuousLinearMap.zero_apply, inner_self_eq_zero] at h
  exact hinj (by rw [h, map_zero])

/-- Where `DF x` is injective, `g_x : E ≃L E*`. -/
def metricEquiv (F : E → V) (x : E) (hinj : Function.Injective (fderiv ℝ F x)) :
    E ≃L[ℝ] (E →L[ℝ] ℝ) :=
  (LinearMap.linearEquivOfInjective (metricMap F x : E →ₗ[ℝ] (E →L[ℝ] ℝ)) (metricMap_injective hinj)
    (by
      rw [← (LinearMap.toContinuousLinearMap : Module.Dual ℝ E ≃ₗ[ℝ] (E →L[ℝ] ℝ)).finrank_eq,
        Subspace.dual_finrank_eq])).toContinuousLinearEquiv

theorem coe_metricEquiv {x : E} (hinj : Function.Injective (fderiv ℝ F x)) :
    (metricEquiv F x hinj : E →L[ℝ] (E →L[ℝ] ℝ)) = metricMap F x :=
  ContinuousLinearMap.ext fun _ => rfl

theorem inner_liftA {x : E} (hinj : Function.Injective (fderiv ℝ F x)) (ξ : E →L[ℝ] ℝ) (u : E) :
    ⟪liftA F x ξ, fderiv ℝ F x u⟫ = ξ u := by
  have h : metricEquiv F x hinj ((metricEquiv F x hinj).symm ξ) u = ξ u := by
    rw [ContinuousLinearEquiv.apply_symm_apply]
  rw [liftA, ContinuousLinearMap.comp_apply, ← coe_metricEquiv hinj,
    ContinuousLinearMap.inverse_equiv]
  rw [← h]
  rfl

theorem inner_J_liftA {J : V →L[ℝ] V} {x : E}
    (hL : ∀ u w, ⟪J (fderiv ℝ F x u), fderiv ℝ F x w⟫ = 0) (ξ η : E →L[ℝ] ℝ) :
    ⟪J (liftA F x ξ), liftA F x η⟫ = 0 :=
  hL _ _

/-- In the cotangent model, `λ_A = λ_can`: `λ_A(x, ξ)(a) = ξ a₁`. -/
theorem lamA_liftA {x : E} (hinj : Function.Injective (fderiv ℝ F x)) (ξ : E →L[ℝ] ℝ)
    (a : E × (E →L[ℝ] ℝ)) : lamA F (liftA F) (x, ξ) a = ξ a.1 := by
  simp only [lamA, ContinuousLinearMap.coe_comp', Function.comp_apply, innerSL_apply_apply,
    ContinuousLinearMap.coe_fst']
  exact inner_liftA hinj ξ a.1

/-- `liftA` is smooth wherever `F` is smooth and `DF x` is injective. -/
theorem contDiffAt_liftA {x : E} (hF : ContDiffAt ℝ ∞ F x)
    (hinj : Function.Injective (fderiv ℝ F x)) : ContDiffAt ℝ ∞ (liftA F) x := by
  have hD := hF.fderiv_right (m := ∞) le_rfl
  have hg : ContDiffAt ℝ ∞ (metricMap F) x := by
    have h1 := (ContinuousLinearMap.contDiff (𝕜 := ℝ) (E := E →L[ℝ] V)
      (F := (V →L[ℝ] ℝ) →L[ℝ] E →L[ℝ] ℝ) (n := ∞)
      ((ContinuousLinearMap.compL ℝ E V ℝ).flip)).contDiffAt.comp x hD
    have h2 := (ContinuousLinearMap.contDiff (𝕜 := ℝ) (E := E →L[ℝ] V)
      (F := E →L[ℝ] V →L[ℝ] ℝ) (n := ∞) ((ContinuousLinearMap.compL ℝ E V (V →L[ℝ] ℝ))
      (innerSL ℝ : V →L[ℝ] V →L[ℝ] ℝ))).contDiffAt.comp x hD
    exact h1.clm_comp h2
  have hI : ContDiffAt ℝ ∞ (fun y => ContinuousLinearMap.inverse (metricMap F y)) x := by
    have h0 := contDiffAt_map_inverse (𝕜 := ℝ) (n := ∞) (metricEquiv F x hinj)
    rw [coe_metricEquiv hinj] at h0
    exact h0.comp x hg
  exact hD.clm_comp hI

end CWFlow
