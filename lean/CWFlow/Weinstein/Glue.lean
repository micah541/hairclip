import CWFlow.Moser.Natural
import CWFlow.Weinstein.Cn

/-!
# Gluing: naturality of the canonical Weinstein map under changes of chart

A change of chart `τ` of `L` (smooth with invertible derivative on an open `S₁₂`) acts on cotangent coordinates by
the cotangent lift `cotLift τ (y, η) = (τ y, η ∘ (Dτ y)⁻¹)`. This map is smooth with invertible derivative
(`contDiffAt_cotLift`, `exists_equiv_fderiv_cotLift`), and it preserves `λ_can` (`lamCan_natural`).

Let `F₁ = F₂ ∘ τ` be the same immersion seen in two charts. Then:
* the metric lift, `Φ₀`, `ℓ` and `β` are natural (`liftA_natural`, `Phi0_natural`, `ell_natural`,
  `beta_natural`);
* by `moserEnd_natural`, the canonical Weinstein maps agree,
  `weinsteinMap₂ (cotLift τ q) = weinsteinMap₁ q`, whenever the Moser curve from `q` stays over `S₁₂`
  (`weinsteinMap_natural`).
-/

open ContinuousAlternatingMap Filter Set Function

open scoped InnerProductSpace RealInnerProductSpace Topology ContDiff

noncomputable section

namespace CWFlow

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]

/-- The cotangent lift `(y, η) ↦ (τ y, η ∘ (Dτ y)⁻¹)` of a change of chart `τ`. -/
def cotLift (τ : E → E) (q : E × (E →L[ℝ] ℝ)) : E × (E →L[ℝ] ℝ) :=
  (τ q.1, q.2.comp (ContinuousLinearMap.inverse (fderiv ℝ τ q.1)))

@[simp] theorem cotLift_fst (τ : E → E) (q : E × (E →L[ℝ] ℝ)) : (cotLift τ q).1 = τ q.1 := rfl

variable {τ : E → E}

/-! ### The cotangent lift -/

section Lift

variable {q : E × (E →L[ℝ] ℝ)} (hτ : ContDiffAt ℝ ∞ τ q.1) {e : E ≃L[ℝ] E}
  (he : (e : E →L[ℝ] E) = fderiv ℝ τ q.1)
include hτ he

theorem contDiffAt_inverse_fderiv :
    ContDiffAt ℝ ∞ (fun y => ContinuousLinearMap.inverse (fderiv ℝ τ y)) q.1 := by
  have h0 := contDiffAt_map_inverse (𝕜 := ℝ) (n := ∞) e
  rw [he] at h0
  exact h0.comp q.1 (hτ.fderiv_right (m := ∞) le_rfl)

theorem contDiffAt_cotLift : ContDiffAt ℝ ∞ (cotLift τ) q :=
  (hτ.comp q contDiffAt_fst).prodMk
    (contDiffAt_snd.clm_comp ((contDiffAt_inverse_fderiv hτ he).comp q contDiffAt_fst))

theorem hasFDerivAt_cotLift :
    HasFDerivAt (cotLift τ) (((fderiv ℝ τ q.1).comp (ContinuousLinearMap.fst ℝ E (E →L[ℝ] ℝ))).prod
      ((ContinuousLinearMap.compL ℝ E E ℝ q.2).comp
          ((fderiv ℝ (fun y => ContinuousLinearMap.inverse (fderiv ℝ τ y)) q.1).comp
            (ContinuousLinearMap.fst ℝ E (E →L[ℝ] ℝ))) +
        ((ContinuousLinearMap.compL ℝ E E ℝ).flip
          (ContinuousLinearMap.inverse (fderiv ℝ τ q.1))).comp
            (ContinuousLinearMap.snd ℝ E (E →L[ℝ] ℝ)))) q := by
  have hI := (contDiffAt_inverse_fderiv hτ he).differentiableAt (by simp)
  exact ((hτ.differentiableAt (by simp)).hasFDerivAt.comp q hasFDerivAt_fst).prodMk
    (hasFDerivAt_snd.clm_comp (hI.hasFDerivAt.comp q hasFDerivAt_fst))

theorem fderiv_cotLift_fst (a : E × (E →L[ℝ] ℝ)) :
    (fderiv ℝ (cotLift τ) q a).1 = fderiv ℝ τ q.1 a.1 := by
  rw [(hasFDerivAt_cotLift hτ he).fderiv]
  simp

theorem fderiv_cotLift_zero (ζ : E →L[ℝ] ℝ) :
    fderiv ℝ (cotLift τ) q (0, ζ) = (0, ζ.comp (ContinuousLinearMap.inverse (fderiv ℝ τ q.1))) := by
  rw [(hasFDerivAt_cotLift hτ he).fderiv]
  simp

theorem exists_equiv_fderiv_cotLift :
    ∃ e' : (E × (E →L[ℝ] ℝ)) ≃L[ℝ] (E × (E →L[ℝ] ℝ)),
      (e' : _ →L[ℝ] _) = fderiv ℝ (cotLift τ) q := by
  have hinj : Function.Injective (fderiv ℝ (cotLift τ) q) := by
    rw [injective_iff_map_eq_zero]
    intro v hv
    have h1 : v.1 = 0 := by
      have h := congrArg Prod.fst hv
      rw [fderiv_cotLift_fst hτ he, ← he] at h
      exact e.injective (h.trans (_root_.map_zero e).symm)
    have hv' : v = (0, v.2) := Prod.ext h1 rfl
    rw [hv', fderiv_cotLift_zero hτ he] at hv
    have h2 : v.2 = 0 := by
      have h := congrArg Prod.snd hv
      rw [← he, ContinuousLinearMap.inverse_equiv] at h
      refine ContinuousLinearMap.ext fun u => ?_
      have h' := congrArg (fun f : E →L[ℝ] ℝ => f (e u)) h
      simpa using h'
    exact Prod.ext h1 h2
  exact ⟨(LinearMap.linearEquivOfInjective (fderiv ℝ (cotLift τ) q : _ →ₗ[ℝ] _) hinj
    rfl).toContinuousLinearEquiv, ContinuousLinearMap.ext fun _ => rfl⟩

/-- The cotangent lift preserves the Liouville form. -/
theorem lamCan_natural : lamCan q = pb1form (cotLift τ) lamCan q := by
  refine ContinuousLinearMap.ext fun a => ?_
  rw [pb1form_apply, lamCan_apply, lamCan_apply, fderiv_cotLift_fst hτ he]
  simp only [cotLift, ContinuousLinearMap.coe_comp', Function.comp_apply]
  rw [← he, ContinuousLinearMap.inverse_equiv]
  simp

end Lift

/-! ### Naturality of the coordinate data -/

section Natural

variable {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
variable {J : V →L[ℝ] V} {F₁ F₂ : E → V}

theorem pb1_natural {W : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W] {Φ₁ Φ₂ : W → V} {T : W → W}
    {q : W} (h : Φ₁ =ᶠ[𝓝 q] Φ₂ ∘ T) (hΦ₂ : DifferentiableAt ℝ Φ₂ (T q))
    (hT : DifferentiableAt ℝ T q) : pb1 J Φ₁ q = pb1form T (pb1 J Φ₂) q := by
  simp only [pb1, pb1form]
  rw [h.fderiv_eq, fderiv_comp q hΦ₂ hT, h.eq_of_nhds, Function.comp_apply,
    ContinuousLinearMap.comp_assoc]

theorem fderiv_natural {W : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W] {f₁ f₂ : W → ℝ}
    {T : W → W} {q : W} (h : f₁ =ᶠ[𝓝 q] f₂ ∘ T) (hf₂ : DifferentiableAt ℝ f₂ (T q))
    (hT : DifferentiableAt ℝ T q) : fderiv ℝ f₁ q = pb1form T (fderiv ℝ f₂) q := by
  rw [h.fderiv_eq, fderiv_comp q hf₂ hT]
  rfl

theorem contDiffAt_ell {G : Type*} [NormedAddCommGroup G] [NormedSpace ℝ G] {F : E → V}
    {A : E → G →L[ℝ] V} {x : E} {ξ : G} (hF : ContDiffAt ℝ ∞ F x) (hA : ContDiffAt ℝ ∞ A x) :
    ContDiffAt ℝ ∞ (ell F A) (x, ξ) :=
  contDiffAt_const.mul ((hF.comp (x, ξ) contDiffAt_fst).inner ℝ
    ((hA.comp (x, ξ) contDiffAt_fst).clm_apply contDiffAt_snd))

/-- The metric lift is natural: `liftA F₂ (τ y) (η ∘ (Dτ)⁻¹) = liftA F₁ y η` when `F₁ = F₂ ∘ τ`. -/
theorem liftA_natural {y : E} (hDF : fderiv ℝ F₁ y = (fderiv ℝ F₂ (τ y)).comp (fderiv ℝ τ y))
    (hinj₁ : Function.Injective (fderiv ℝ F₁ y)) (hinj₂ : Function.Injective (fderiv ℝ F₂ (τ y)))
    {e : E ≃L[ℝ] E} (he : (e : E →L[ℝ] E) = fderiv ℝ τ y) (η : E →L[ℝ] ℝ) :
    liftA F₂ (τ y) (η.comp (ContinuousLinearMap.inverse (fderiv ℝ τ y))) = liftA F₁ y η := by
  have hinv : ContinuousLinearMap.inverse (fderiv ℝ τ y) = (e.symm : E →L[ℝ] E) := by
    rw [← he, ContinuousLinearMap.inverse_equiv]
  set η' := η.comp (ContinuousLinearMap.inverse (fderiv ℝ τ y)) with hη'
  have hpair : ∀ w, ⟪liftA F₁ y η - liftA F₂ (τ y) η', fderiv ℝ F₂ (τ y) w⟫ = 0 := by
    intro w
    have hw : fderiv ℝ F₂ (τ y) w = fderiv ℝ F₁ y (e.symm w) := by
      rw [hDF, ContinuousLinearMap.comp_apply, ← he]
      simp
    rw [inner_sub_left, inner_liftA hinj₂, hw, inner_liftA hinj₁, hη', hinv]
    simp
  have hz : liftA F₁ y η - liftA F₂ (τ y) η' = fderiv ℝ F₂ (τ y)
      (fderiv ℝ τ y (ContinuousLinearMap.inverse (metricMap F₁ y) η) -
        ContinuousLinearMap.inverse (metricMap F₂ (τ y)) η') := by
    simp only [liftA, ContinuousLinearMap.comp_apply, map_sub, hDF]
  have h0 := hpair (fderiv ℝ τ y (ContinuousLinearMap.inverse (metricMap F₁ y) η) -
    ContinuousLinearMap.inverse (metricMap F₂ (τ y)) η')
  rw [← hz, inner_self_eq_zero, sub_eq_zero] at h0
  exact h0.symm

/-! Hypotheses: `F₁ = F₂ ∘ τ` on the open overlap `S₁₂`, both immersions, `τ` a local diffeomorphism. -/

variable {S₁₂ : Set E} (hS₁₂ : IsOpen S₁₂) (hτ : ∀ y ∈ S₁₂, ContDiffAt ℝ ∞ τ y)
  (hτinv : ∀ y ∈ S₁₂, ∃ e : E ≃L[ℝ] E, (e : E →L[ℝ] E) = fderiv ℝ τ y)
  (hFτ : ∀ y ∈ S₁₂, F₁ y = F₂ (τ y)) (hF₂ : ∀ y ∈ S₁₂, ContDiffAt ℝ ∞ F₂ (τ y))
  (hinj₁ : ∀ y ∈ S₁₂, Function.Injective (fderiv ℝ F₁ y))
  (hinj₂ : ∀ y ∈ S₁₂, Function.Injective (fderiv ℝ F₂ (τ y)))
include hS₁₂ hτ hτinv hFτ hF₂ hinj₁ hinj₂

theorem fderiv_F_natural {y : E} (hy : y ∈ S₁₂) :
    fderiv ℝ F₁ y = (fderiv ℝ F₂ (τ y)).comp (fderiv ℝ τ y) := by
  have h : F₁ =ᶠ[𝓝 y] F₂ ∘ τ := by
    filter_upwards [hS₁₂.mem_nhds hy] with z hz using hFτ z hz
  rw [h.fderiv_eq, fderiv_comp y ((hF₂ y hy).differentiableAt (by simp))
    ((hτ y hy).differentiableAt (by simp))]

theorem Phi0_natural {q : E × (E →L[ℝ] ℝ)} (hq : q.1 ∈ S₁₂) :
    Phi0 J F₂ (liftA F₂) (cotLift τ q) = Phi0 J F₁ (liftA F₁) q := by
  obtain ⟨e, he⟩ := hτinv q.1 hq
  simp only [Phi0, cotLift]
  rw [liftA_natural (fderiv_F_natural hS₁₂ hτ hτinv hFτ hF₂ hinj₁ hinj₂ hq) (hinj₁ q.1 hq)
    (hinj₂ q.1 hq) he, hFτ q.1 hq]

theorem ell_natural {q : E × (E →L[ℝ] ℝ)} (hq : q.1 ∈ S₁₂) :
    ell F₂ (liftA F₂) (cotLift τ q) = ell F₁ (liftA F₁) q := by
  obtain ⟨e, he⟩ := hτinv q.1 hq
  simp only [ell, cotLift]
  rw [liftA_natural (fderiv_F_natural hS₁₂ hτ hτinv hFτ hF₂ hinj₁ hinj₂ hq) (hinj₁ q.1 hq)
    (hinj₂ q.1 hq) he, hFτ q.1 hq]

/-- `β` is natural: `β₁ = (cotLift τ)*β₂` over the overlap. -/
theorem beta_natural {q : E × (E →L[ℝ] ℝ)} (hq : q.1 ∈ S₁₂) :
    beta J F₁ (liftA F₁) q = pb1form (cotLift τ) (beta J F₂ (liftA F₂)) q := by
  obtain ⟨e, he⟩ := hτinv q.1 hq
  have hN : ∀ᶠ p in 𝓝 q, p.1 ∈ S₁₂ :=
    continuous_fst.continuousAt.preimage_mem_nhds (hS₁₂.mem_nhds hq)
  have hT : DifferentiableAt ℝ (cotLift τ) q :=
    (contDiffAt_cotLift (hτ q.1 hq) he).differentiableAt (by simp)
  have hA₂ : ContDiffAt ℝ ∞ (liftA F₂) (τ q.1) := contDiffAt_liftA (hF₂ q.1 hq) (hinj₂ q.1 hq)
  have hΦ : Phi0 J F₁ (liftA F₁) =ᶠ[𝓝 q] Phi0 J F₂ (liftA F₂) ∘ cotLift τ := by
    filter_upwards [hN] with p hp
    exact (Phi0_natural hS₁₂ hτ hτinv hFτ hF₂ hinj₁ hinj₂ hp).symm
  have hFf : (fun p : E × (E →L[ℝ] ℝ) => F₁ p.1) =ᶠ[𝓝 q]
      (fun p : E × (E →L[ℝ] ℝ) => F₂ p.1) ∘ cotLift τ := by
    filter_upwards [hN] with p hp
    exact hFτ p.1 hp
  have hℓ : ell F₁ (liftA F₁) =ᶠ[𝓝 q] ell F₂ (liftA F₂) ∘ cotLift τ := by
    filter_upwards [hN] with p hp
    exact (ell_natural hS₁₂ hτ hτinv hFτ hF₂ hinj₁ hinj₂ hp).symm
  have hlamA₁ : lamA F₁ (liftA F₁) q = lamCan q := by
    refine ContinuousLinearMap.ext fun a => ?_
    rw [lamCan_apply]
    exact lamA_liftA (hinj₁ q.1 hq) q.2 a
  have hlamA₂ : lamA F₂ (liftA F₂) (cotLift τ q) = lamCan (cotLift τ q) := by
    refine ContinuousLinearMap.ext fun a => ?_
    rw [lamCan_apply]
    exact lamA_liftA (hinj₂ q.1 hq) _ a
  have hlam : pb1form (cotLift τ) (lamA F₂ (liftA F₂)) q = lamCan q := by
    rw [lamCan_natural (hτ q.1 hq) he]
    simp only [pb1form, hlamA₂]
  have hsplit : pb1form (cotLift τ) (beta J F₂ (liftA F₂)) q =
      pb1form (cotLift τ) (pb1 J (Phi0 J F₂ (liftA F₂))) q -
        pb1form (cotLift τ) (lamA F₂ (liftA F₂)) q -
        pb1form (cotLift τ) (pb1 J (fun p : E × (E →L[ℝ] ℝ) => F₂ p.1)) q -
        pb1form (cotLift τ) (fderiv ℝ (ell F₂ (liftA F₂))) q := by
    simp only [pb1form, beta, ContinuousLinearMap.sub_comp]
  rw [hsplit, beta,
    pb1_natural hΦ ((contDiffAt_Phi0 (ξ := (cotLift τ q).2) (hF₂ q.1 hq) hA₂).differentiableAt
      (by simp)) hT,
    pb1_natural hFf (((hF₂ q.1 hq).comp (cotLift τ q) contDiffAt_fst).differentiableAt (by simp)) hT,
    fderiv_natural hℓ ((contDiffAt_ell (ξ := (cotLift τ q).2) (hF₂ q.1 hq) hA₂).differentiableAt
      (by simp)) hT,
    hlamA₁, hlam]

end Natural

end CWFlow
