import CWFlow.Weinstein.Existence

/-!
# The Weinstein neighbourhood theorem for immersed Lagrangians in `ℂⁿ` (chart version)

The `ℂⁿ` instance of `isLocalWeinstein_weinsteinMap`. Take `Φ₀ = Phi0 J F A` and `β = beta J F A` with `A = liftA F`
the metric lift. Then `λ_A = λ_can` (`lamA_liftA`), `dβ = Φ₀*ω₀ − ω_can` (`dOne_beta`), and `Φ₀*ω₀ = ω_can` on
the zero section (`pb2_Phi0_zero`, `inner_liftA`).

`isLocalWeinstein_Cn`: let `F` be a smooth Lagrangian immersion on an open `S ⊆ E`. Then the canonical map is a
local Weinstein map for `ω₀` on an open set containing the zero section over `S`. No embeddedness is assumed: `F` need
not be injective.
-/

open ContinuousAlternatingMap Filter Set

open scoped InnerProductSpace RealInnerProductSpace Topology ContDiff

noncomputable section

namespace CWFlow

variable {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
variable {J : V →L[ℝ] V} {F : E → V}

omit [FiniteDimensional ℝ E] in
theorem contDiffAt_Phi0 {G : Type*} [NormedAddCommGroup G] [NormedSpace ℝ G] {A : E → G →L[ℝ] V}
    {x : E} {ξ : G} (hF : ContDiffAt ℝ ∞ F x) (hA : ContDiffAt ℝ ∞ A x) :
    ContDiffAt ℝ ∞ (Phi0 J F A) (x, ξ) :=
  (hF.comp (x, ξ) contDiffAt_fst).sub
    (J.contDiff.contDiffAt.comp (x, ξ) ((hA.comp (x, ξ) contDiffAt_fst).clm_apply contDiffAt_snd))

/-- Over the set where `DF` is injective, `λ_A = λ_can` for the metric lift. -/
theorem lamA_liftA_eventuallyEq {S : Set E} (hS : IsOpen S)
    (hinj : ∀ x ∈ S, Function.Injective (fderiv ℝ F x)) {q : E × (E →L[ℝ] ℝ)} (hq : q.1 ∈ S) :
    lamA F (liftA F) =ᶠ[𝓝 q] lamCan (E := E) := by
  filter_upwards [continuous_fst.continuousAt.preimage_mem_nhds (hS.mem_nhds hq)] with ⟨x, ξ⟩ hx
  refine ContinuousLinearMap.ext fun a => ?_
  rw [lamCan_apply]
  exact lamA_liftA (hinj x hx) ξ a

/-- **Weinstein neighbourhoods for immersed Lagrangians in `ℂⁿ`, in a chart.** The canonical map
`weinsteinMap (Phi0 J F A) (beta J F A) O` (`A = liftA F`, `O = S × E*`) is a local Weinstein map on an open set
containing the zero section over `S`. -/
theorem isLocalWeinstein_Cn (hJi : ∀ u v, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v, J (J v) = -v)
    {S : Set E} (hS : IsOpen S) (hF : ContDiffOn ℝ ∞ F S)
    (hinj : ∀ x ∈ S, Function.Injective (fderiv ℝ F x))
    (hL : ∀ x ∈ S, ∀ u w, ⟪J (fderiv ℝ F x u), fderiv ℝ F x w⟫ = 0) :
    IsLocalWeinstein (omega0 J) F (moserDom lamCan (beta J F (liftA F)) (Prod.fst ⁻¹' S))
        (weinsteinMap (Phi0 J F (liftA F)) (beta J F (liftA F)) (Prod.fst ⁻¹' S)) ∧
      ∀ x ∈ S, (x, (0 : E →L[ℝ] ℝ)) ∈ moserDom lamCan (beta J F (liftA F)) (Prod.fst ⁻¹' S) := by
  have hFA : ∀ x ∈ S, ContDiffAt ℝ ∞ F x := fun x hx => hF.contDiffAt (hS.mem_nhds hx)
  have hAA : ∀ x ∈ S, ContDiffAt ℝ ∞ (liftA F) x := fun x hx =>
    contDiffAt_liftA (hFA x hx) (hinj x hx)
  have hO : IsOpen (Prod.fst ⁻¹' S : Set (E × (E →L[ℝ] ℝ))) := hS.preimage continuous_fst
  have hdlam : ∀ q ∈ (Prod.fst ⁻¹' S : Set (E × (E →L[ℝ] ℝ))), ∀ a b,
      dOne (lamA F (liftA F)) q a b = omegaCanVal a b := by
    intro q hq a b
    rw [dOne, (lamA_liftA_eventuallyEq hS hinj hq).fderiv_eq, ← dOne, dOne_lamCan]
  have h := isLocalWeinstein_weinsteinMap (Ω := omega0 J) (F := F) (Φ₀ := Phi0 J F (liftA F))
    (β := beta J F (liftA F)) hO ?_ ?_ ?_ ?_ ?_ ?_
  · exact ⟨h.1, fun x hx => h.2 x hx⟩
  · rintro ⟨x, ξ⟩ hx
    exact (contDiffAt_Phi0 (hFA x hx) (hAA x hx)).contDiffWithinAt
  · rintro ⟨x, ξ⟩ hx
    exact (contDiffAt_beta (hFA x hx) (hAA x hx)).contDiffWithinAt
  · intro x _
    simp [Phi0]
  · intro x hx
    refine ContinuousLinearMap.ext fun v => ?_
    rw [ContinuousLinearMap.zero_apply]
    exact beta_zero hJi ((hFA x hx).differentiableAt (by simp)) ((hAA x hx).differentiableAt (by simp)) v
  · rintro ⟨x, ξ⟩ hx a b
    rw [pullback2_omega0 hJi hJJ, dOne_beta hJi hJJ ((hFA x hx).of_le two_le_infty)
      ((hAA x hx).of_le two_le_infty) (hL x hx), hdlam _ hx]
    ring
  · intro x hx a b
    have hxS : x ∈ S := hx
    rw [pullback2_omega0 hJi hJJ, pb2_Phi0_zero hJi hJJ ((hFA x hxS).differentiableAt (by simp))
      ((hAA x hxS).differentiableAt (by simp)) (hL x hxS) (fun ξ η => inner_J_liftA (hL x hxS) ξ η),
      inner_liftA (hinj x hxS), inner_liftA (hinj x hxS)]
    rfl

end CWFlow
