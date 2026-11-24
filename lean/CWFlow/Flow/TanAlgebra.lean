import CWFlow.Flow.Tangential

/-!
# Algebra of the tangential coordinates `tanPart`

Uniqueness, linearity, behaviour under `L ↦ L ∘ A`, smoothness in `(L, v)` and locality of `flowVel`.
-/

open Filter Set

open scoped InnerProductSpace RealInnerProductSpace Topology ContDiff

noncomputable section

namespace CWFlow

variable {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]

theorem tanPart_unique {L : E →L[ℝ] V} (hL : Function.Injective L) {v : V} {e : E}
    (h : ∀ a, ⟪v - L e, L a⟫ = 0) : tanPart L v = e := by
  have h1 : ⟪L (tanPart L v - e), L (tanPart L v - e)⟫ = 0 := by
    have := congrArg₂ (· - ·) (h (tanPart L v - e)) (inner_sub_tanPart hL v (tanPart L v - e))
    simp only [← inner_sub_left, sub_zero] at this
    rw [← this]
    congr 1
    simp
  have h2 : L (tanPart L v - e) = 0 := inner_self_eq_zero.mp h1
  have h3 : tanPart L v - e = 0 := hL (by rw [h2, map_zero])
  exact sub_eq_zero.mp h3

theorem tanPart_apply_self {L : E →L[ℝ] V} (hL : Function.Injective L) (c : E) :
    tanPart L (L c) = c :=
  tanPart_unique hL (fun a => by simp)

theorem tanPart_add (L : E →L[ℝ] V) (v w : V) :
    tanPart L (v + w) = tanPart L v + tanPart L w := by
  simp [tanPart]

theorem tanPart_sub (L : E →L[ℝ] V) (v w : V) :
    tanPart L (v - w) = tanPart L v - tanPart L w := by
  simp [tanPart]

theorem tanPart_eq_zero {L : E →L[ℝ] V} (hL : Function.Injective L) {v : V}
    (h : ∀ a, ⟪v, L a⟫ = 0) : tanPart L v = 0 :=
  tanPart_unique hL (fun a => by simp [h a])

theorem tanPart_J_apply {L : E →L[ℝ] V} (hL : Function.Injective L) {J : V →L[ℝ] V}
    (hLag : ∀ a b, ⟪J (L a), L b⟫ = 0) (c : E) : tanPart L (J (L c)) = 0 :=
  tanPart_eq_zero hL (fun a => hLag c a)

theorem tanPart_comp_equiv {L : E →L[ℝ] V} (hL : Function.Injective L) (A : E ≃L[ℝ] E) (v : V) :
    tanPart (L.comp (A : E →L[ℝ] E)) v = A.symm (tanPart L v) := by
  have hinj : Function.Injective (L.comp (A : E →L[ℝ] E)) := hL.comp A.injective
  refine tanPart_unique hinj (fun a => ?_)
  have := inner_sub_tanPart hL v (A a)
  simpa using this

theorem sub_apply_tanPart_of_eq {L : E →L[ℝ] V} (hL : Function.Injective L) {J : V →L[ℝ] V}
    (hLag : ∀ a b, ⟪J (L a), L b⟫ = 0) {v : V} {b c : E} (h : v = J (L b) + L c) :
    v - L (tanPart L v) = J (L b) := by
  have : tanPart L v = c := by
    rw [h, tanPart_add, tanPart_J_apply hL hLag, tanPart_apply_self hL, zero_add]
  rw [this, h]; abel

theorem contDiffAt_tanPart {L : E →L[ℝ] V} (hL : Function.Injective L) (v : V) :
    ContDiffAt ℝ ∞ (fun p : (E →L[ℝ] V) × V => tanPart p.1 p.2) (L, v) := by
  have h1 : ContDiffAt ℝ ∞ (fun p : (E →L[ℝ] V) × V => ginvOf p.1) (L, v) :=
    ContDiffAt.comp (g := (ginvOf : (E →L[ℝ] V) → (E →L[ℝ] ℝ) →L[ℝ] E))
      (f := fun p : (E →L[ℝ] V) × V => p.1) (L, v) (contDiffAt_ginvOf hL) contDiffAt_fst
  have h2 : ContDiffAt ℝ ∞ (fun p : (E →L[ℝ] V) × V =>
      ((innerSL ℝ : V →L[ℝ] V →L[ℝ] ℝ) p.2).comp p.1) (L, v) :=
    ContDiffAt.clm_comp ((innerSL ℝ : V →L[ℝ] V →L[ℝ] ℝ).contDiff.contDiffAt.comp (L, v)
      contDiffAt_snd) contDiffAt_fst
  exact h1.clm_apply h2

theorem flowVel_congr_of_eventuallyEq {J : V →L[ℝ] V} {X₁ X₂ : E → V} {y : E}
    (h : X₁ =ᶠ[𝓝 y] X₂) : flowVel J X₁ y = flowVel J X₂ y := by
  have h1 : fderiv ℝ X₁ =ᶠ[𝓝 y] fderiv ℝ X₂ := h.fderiv
  have h2 : fderiv ℝ (fderiv ℝ X₁) =ᶠ[𝓝 y] fderiv ℝ (fderiv ℝ X₂) := h1.fderiv
  have h3 : (fun z => J (meanCurv X₁ z)) =ᶠ[𝓝 y] fun z => J (meanCurv X₂ z) := by
    filter_upwards [h1, h2] with z hz1 hz2
    simp only [meanCurv, hz1, hz2]
  have h4 : flowPot J X₁ =ᶠ[𝓝 y] flowPot J X₂ := by
    filter_upwards [h1.eventuallyEq_nhds, h3.eventuallyEq_nhds] with z hz1 hz3
    simp only [flowPot, hz1.self_of_nhds, hz3.fderiv_eq]
  simp only [flowVel, h1.self_of_nhds, h4.fderiv_eq]

end CWFlow
