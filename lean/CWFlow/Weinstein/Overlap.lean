import CWFlow.Weinstein.Glue

/-!
# The canonical Weinstein maps of two charts agree

`weinsteinMap_natural`: let `F₁` on `S₁` and `F₂` on `S₂` be the same immersion in two charts, with
`F₁ = F₂ ∘ τ` on the overlap `S₁₂`. Then the canonical Weinstein maps satisfy
`Ψ₂ (cotLift τ q) = Ψ₁ q` whenever the Moser curve from `q` stays over `S₁₂`. This is the gluing statement for
the Weinstein map on `T*L`.
-/

open Set Function

open scoped ContDiff

noncomputable section

namespace CWFlow

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
variable {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
variable {J : V →L[ℝ] V} {F₁ F₂ : E → V} {τ : E → E}

theorem contDiffOn_beta_liftA {F : E → V} {S : Set E} (hS : IsOpen S) (hF : ContDiffOn ℝ ∞ F S)
    (hinj : ∀ x ∈ S, Function.Injective (fderiv ℝ F x)) :
    ContDiffOn ℝ ∞ (beta J F (liftA F)) (Prod.fst ⁻¹' S) := by
  rintro ⟨x, ξ⟩ hx
  have hFx := hF.contDiffAt (hS.mem_nhds hx)
  exact (contDiffAt_beta hFx (contDiffAt_liftA hFx (hinj x hx))).contDiffWithinAt

/-- **Gluing.** The canonical Weinstein maps of two charts agree under the cotangent lift of the change of
chart, along Moser curves staying over the overlap. -/
theorem weinsteinMap_natural {S₁ S₂ S₁₂ : Set E} (hS₁ : IsOpen S₁) (hS₂ : IsOpen S₂)
    (hS₁₂ : IsOpen S₁₂) (hsub : S₁₂ ⊆ S₁) (hmaps : MapsTo τ S₁₂ S₂)
    (hF₁ : ContDiffOn ℝ ∞ F₁ S₁) (hF₂ : ContDiffOn ℝ ∞ F₂ S₂)
    (hinj₁ : ∀ x ∈ S₁, Function.Injective (fderiv ℝ F₁ x))
    (hinj₂ : ∀ x ∈ S₂, Function.Injective (fderiv ℝ F₂ x))
    (hτ : ContDiffOn ℝ ∞ τ S₁₂) (hτinv : ∀ y ∈ S₁₂, ∃ e : E ≃L[ℝ] E, (e : E →L[ℝ] E) = fderiv ℝ τ y)
    (hFτ : ∀ y ∈ S₁₂, F₁ y = F₂ (τ y)) {q : E × (E →L[ℝ] ℝ)} {γ : ℝ → E × (E →L[ℝ] ℝ)}
    (hγ : IsMoserCurve lamCan (beta J F₁ (liftA F₁)) (Prod.fst ⁻¹' S₁) q γ)
    (hγN : ∀ t ∈ Icc (0 : ℝ) 1, (γ t).1 ∈ S₁₂) :
    weinsteinMap (Phi0 J F₂ (liftA F₂)) (beta J F₂ (liftA F₂)) (Prod.fst ⁻¹' S₂) (cotLift τ q) =
      weinsteinMap (Phi0 J F₁ (liftA F₁)) (beta J F₁ (liftA F₁)) (Prod.fst ⁻¹' S₁) q := by
  have hτA : ∀ y ∈ S₁₂, ContDiffAt ℝ ∞ τ y := fun y hy => hτ.contDiffAt (hS₁₂.mem_nhds hy)
  have hF₂A : ∀ y ∈ S₁₂, ContDiffAt ℝ ∞ F₂ (τ y) := fun y hy =>
    hF₂.contDiffAt (hS₂.mem_nhds (hmaps hy))
  have hinj₁' : ∀ y ∈ S₁₂, Function.Injective (fderiv ℝ F₁ y) := fun y hy => hinj₁ y (hsub hy)
  have hinj₂' : ∀ y ∈ S₁₂, Function.Injective (fderiv ℝ F₂ (τ y)) := fun y hy =>
    hinj₂ _ (hmaps hy)
  have hO₁ : IsOpen (Prod.fst ⁻¹' S₁ : Set (E × (E →L[ℝ] ℝ))) := hS₁.preimage continuous_fst
  have hO₂ : IsOpen (Prod.fst ⁻¹' S₂ : Set (E × (E →L[ℝ] ℝ))) := hS₂.preimage continuous_fst
  have hN : IsOpen (Prod.fst ⁻¹' S₁₂ : Set (E × (E →L[ℝ] ℝ))) := hS₁₂.preimage continuous_fst
  have hTN : MapsTo (cotLift τ) (Prod.fst ⁻¹' S₁₂) (Prod.fst ⁻¹' S₂) := fun p hp => hmaps hp
  have hT : ContDiffOn ℝ ∞ (cotLift τ) (Prod.fst ⁻¹' S₁₂) := fun p hp => by
    obtain ⟨e, he⟩ := hτinv p.1 hp
    exact (contDiffAt_cotLift (hτA p.1 hp) he).contDiffWithinAt
  have hTinv : ∀ p ∈ (Prod.fst ⁻¹' S₁₂ : Set (E × (E →L[ℝ] ℝ))),
      ∃ e : (E × (E →L[ℝ] ℝ)) ≃L[ℝ] (E × (E →L[ℝ] ℝ)),
        (e : _ →L[ℝ] _) = fderiv ℝ (cotLift τ) p := fun p hp => by
    obtain ⟨e, he⟩ := hτinv p.1 hp
    exact exists_equiv_fderiv_cotLift (hτA p.1 hp) he
  have hlamT : ∀ p ∈ (Prod.fst ⁻¹' S₁₂ : Set (E × (E →L[ℝ] ℝ))),
      lamCan p = pb1form (cotLift τ) lamCan p := fun p hp => by
    obtain ⟨e, he⟩ := hτinv p.1 hp
    exact lamCan_natural (hτA p.1 hp) he
  have hβT : ∀ p ∈ (Prod.fst ⁻¹' S₁₂ : Set (E × (E →L[ℝ] ℝ))),
      beta J F₁ (liftA F₁) p = pb1form (cotLift τ) (beta J F₂ (liftA F₂)) p := fun p hp =>
    beta_natural hS₁₂ hτA hτinv hFτ hF₂A hinj₁' hinj₂' hp
  have hend := moserEnd_natural hO₂ contDiff_lamCan.contDiffOn (contDiffOn_beta_liftA hS₂ hF₂ hinj₂)
    hN hTN hT hTinv hlamT hβT hO₁ contDiff_lamCan.contDiffOn (contDiffOn_beta_liftA hS₁ hF₁ hinj₁)
    hγ hγN
  have h1 : (moserEnd lamCan (beta J F₁ (liftA F₁)) (Prod.fst ⁻¹' S₁) q).1 ∈ S₁₂ := by
    rw [moserEnd_eq hO₁ contDiff_lamCan.contDiffOn (contDiffOn_beta_liftA hS₁ hF₁ hinj₁) hγ]
    exact hγN 1 ⟨zero_le_one, le_rfl⟩
  simp only [weinsteinMap, Function.comp_apply]
  rw [hend, Phi0_natural hS₁₂ hτA hτinv hFτ hF₂A hinj₁' hinj₂' h1]

end CWFlow
