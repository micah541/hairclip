import CWFlow.Extension.FlowM
import CWFlow.Extension.VolEvol

/-!
# Scaling of the coordinate geometry, and the parabolic rescaling of a flow

For a real `c ≠ 0` and a map `X : ℝⁿ → V`, the coordinate quantities of `c • X` are:
`g ↦ c² g`, `Γ` invariant, `C ↦ c² C`, `∇^k C ↦ c² ∇^k C`, `|C|² ↦ c⁻² |C|²`,
`potX ↦ c⁻² potX`, `velX ↦ c⁻⁵ velX`.
Consequently, if `Φ` is a flow (`IsFlowM`) then so is `Φ̃(s) = c • Φ(c⁻⁴ s)`.
-/

open Set
open scoped ContDiff Matrix Topology RealInnerProductSpace Manifold

set_option linter.unusedSectionVars false

namespace CWFlow.Extension

variable {n : ℕ} {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V] {J : V →L[ℝ] V}

/-! ### Derivatives -/

theorem rs_fderiv_smul {F G : Type*} [NormedAddCommGroup F] [NormedSpace ℝ F]
    [NormedAddCommGroup G] [NormedSpace ℝ G] (c : ℝ) (X : F → G) (x : F) :
    fderiv ℝ (fun y => c • X y) x = c • fderiv ℝ X x := by
  rw [show (fun y => c • X y) = c • X from rfl, fderiv_const_smul_field, Pi.smul_apply]

theorem dX_smul (c : ℝ) (X : EuclideanSpace ℝ (Fin n) → V) (x a) :
    dX (fun y => c • X y) x a = c • dX X x a := by
  simp [dX, rs_fderiv_smul]

theorem ddX_smul (c : ℝ) (X : EuclideanSpace ℝ (Fin n) → V) (x a b) :
    ddX (fun y => c • X y) x a b = c • ddX X x a b := by
  have h : fderiv ℝ (fun y => c • X y) = fun y => c • fderiv ℝ X y := by
    funext y; exact rs_fderiv_smul c X y
  simp only [ddX]
  rw [h, rs_fderiv_smul c (fderiv ℝ X)]
  simp

theorem gX_smul (c : ℝ) (X : EuclideanSpace ℝ (Fin n) → V) (x) :
    gX (fun y => c • X y) x = c ^ 2 • gX X x := by
  ext a b
  simp only [gX, dX_smul, real_inner_smul_left, real_inner_smul_right, Matrix.smul_apply,
    smul_eq_mul]
  ring

theorem rs_inv_smul {c : ℝ} (hc : c ≠ 0) (A : Matrix (Fin n) (Fin n) ℝ) :
    (c ^ 2 • A)⁻¹ = (c ^ 2)⁻¹ • A⁻¹ := by
  have hc2 : c ^ 2 ≠ 0 := pow_ne_zero 2 hc
  by_cases h : IsUnit A.det
  · have := Matrix.inv_smul' (A := A) (Units.mk0 (c ^ 2) hc2) h
    simpa [Units.smul_def] using this
  · have h2 : ¬ IsUnit (c ^ 2 • A).det := by
      rw [Matrix.det_smul, isUnit_iff_ne_zero]
      rw [isUnit_iff_ne_zero] at h
      simpa [hc2] using h
    rw [Matrix.nonsing_inv_apply_not_isUnit _ h2, Matrix.nonsing_inv_apply_not_isUnit _ h]
    simp

theorem rs_gX_inv_smul {c : ℝ} (hc : c ≠ 0) (X : EuclideanSpace ℝ (Fin n) → V) (x a b) :
    (gX (fun y => c • X y) x)⁻¹ a b = (c ^ 2)⁻¹ * (gX X x)⁻¹ a b := by
  rw [gX_smul, rs_inv_smul hc]
  simp

theorem rs_ginvZ_smul {c : ℝ} (hc : c ≠ 0) (X : EuclideanSpace ℝ (Fin n) → V) (x) :
    ginvZ (fun y => c • X y) x = (c ^ 2)⁻¹ • ginvZ X x := by
  ext I
  simp [ginvZ, rs_gX_inv_smul hc]

theorem cX_smul (c : ℝ) (X : EuclideanSpace ℝ (Fin n) → V) (x) :
    cX J (fun y => c • X y) x = c ^ 2 • cX J X x := by
  ext I
  simp only [cX, ddX_smul, dX_smul, map_smul, real_inner_smul_left, real_inner_smul_right,
    PiLp.smul_apply, smul_eq_mul]
  ring

theorem gamX_smul {c : ℝ} (hc : c ≠ 0) (X : EuclideanSpace ℝ (Fin n) → V) (x)
    (hinj : Function.Injective (fderiv ℝ X x)) :
    gamX (fun y => c • X y) x = gamX X x := by
  ext I
  simp only [gamX, rs_gX_inv_smul hc, ddX_smul, dX_smul, real_inner_smul_left,
    real_inner_smul_right]
  refine Finset.sum_congr rfl fun d _ => ?_
  field_simp

/-! ### `SX` -/

theorem SX_smul {c : ℝ} (hc : c ≠ 0) {X : EuclideanSpace ℝ (Fin n) → V}
    {Ω : Set (EuclideanSpace ℝ (Fin n))} (hΩ : IsOpen Ω) (hX : ContDiffOn ℝ ∞ X Ω)
    (hinj : ∀ x ∈ Ω, Function.Injective (fderiv ℝ X x)) (k : ℕ) {x} (hx : x ∈ Ω) :
    SX J (fun y => c • X y) k x = c ^ 2 • SX J X k x := by
  induction k generalizing x with
  | zero => exact cX_smul c X x
  | succ k ih =>
    have hev : SX J (fun y => c • X y) k =ᶠ[𝓝 x] fun y => c ^ 2 • SX J X k y := by
      filter_upwards [hΩ.mem_nhds hx] with y hy using ih hy
    rw [SX_succ, SX_succ, hev.fderiv_eq, rs_fderiv_smul, gamX_smul hc X x (hinj x hx), ih hx]
    simp [smul_sub]

/-! ### `aNormSq`, `potX`, `velX` -/

theorem rs_gNormSq_smul_left {m : ℕ} (a : ℝ) (h : Z n 2) (T : Z n m) :
    gNormSq (a • h) T = a ^ m * gNormSq h T := by
  unfold gNormSq
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl fun I _ => ?_
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl fun J _ => ?_
  simp only [PiLp.smul_apply, smul_eq_mul, Finset.prod_mul_distrib, Finset.prod_const,
    Finset.card_univ, Fintype.card_fin]
  ring

theorem rs_gNormSq_smul_right {m : ℕ} (b : ℝ) (h : Z n 2) (T : Z n m) :
    gNormSq h (b • T) = b ^ 2 * gNormSq h T := by
  unfold gNormSq
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl fun I _ => ?_
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl fun J _ => ?_
  simp only [PiLp.smul_apply, smul_eq_mul]
  ring

theorem aNormSq_smul {c : ℝ} (hc : c ≠ 0) (X : EuclideanSpace ℝ (Fin n) → V) (x)
    (hinj : Function.Injective (fderiv ℝ X x)) :
    aNormSq J (fun y => c • X y) x = (c ^ 2)⁻¹ * aNormSq J X x := by
  unfold aNormSq
  rw [rs_ginvZ_smul hc, cX_smul, rs_gNormSq_smul_left, rs_gNormSq_smul_right]
  field_simp

theorem potX_smul {c : ℝ} (hc : c ≠ 0) {X : EuclideanSpace ℝ (Fin n) → V}
    {Ω : Set (EuclideanSpace ℝ (Fin n))} (hΩ : IsOpen Ω) (hX : ContDiffOn ℝ ∞ X Ω)
    (hinj : ∀ x ∈ Ω, Function.Injective (fderiv ℝ X x)) {x} (hx : x ∈ Ω) :
    potX J (fun y => c • X y) x = (c ^ 2)⁻¹ * potX J X x := by
  unfold potX
  rw [SX_smul hc hΩ hX hinj 1 hx, mul_neg]
  congr 1
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl fun a _ => ?_
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl fun b _ => ?_
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl fun d _ => ?_
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl fun e _ => ?_
  simp only [rs_ginvZ_smul hc, PiLp.smul_apply, smul_eq_mul]
  field_simp

theorem velX_smul {c : ℝ} (hc : c ≠ 0) {X : EuclideanSpace ℝ (Fin n) → V}
    {Ω : Set (EuclideanSpace ℝ (Fin n))} (hΩ : IsOpen Ω) (hX : ContDiffOn ℝ ∞ X Ω)
    (hinj : ∀ x ∈ Ω, Function.Injective (fderiv ℝ X x)) {x} (hx : x ∈ Ω) :
    velX J (fun y => c • X y) x = ((c ^ 2)⁻¹ * (c ^ 2)⁻¹ * c) • velX J X x := by
  have hev : potX J (fun y => c • X y) =ᶠ[𝓝 x] fun y => (c ^ 2)⁻¹ • potX J X y := by
    filter_upwards [hΩ.mem_nhds hx] with y hy
    simpa using potX_smul (J := J) hc hΩ hX hinj hy
  have hd : fderiv ℝ (potX J (fun y => c • X y)) x = (c ^ 2)⁻¹ • fderiv ℝ (potX J X) x := by
    rw [hev.fderiv_eq, rs_fderiv_smul]
  unfold velX
  rw [Finset.smul_sum]
  refine Finset.sum_congr rfl fun a _ => ?_
  rw [Finset.smul_sum]
  refine Finset.sum_congr rfl fun b _ => ?_
  rw [hd, rs_gX_inv_smul hc, dX_smul, map_smul, smul_smul]
  simp only [smul_apply, smul_eq_mul, smul_smul]
  congr 1
  ring

/-! ### Lagrangian immersions and flows -/

section Flow

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
  {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
  (ι : EuclideanSpace ℝ (Fin n) ≃L[ℝ] E)

theorem IsLagImm.smul {F : M → V} (hF : IsLagImm ι J F) {c : ℝ} (hc : c ≠ 0) :
    IsLagImm ι J (fun p => c • F p) where
  smooth z := (hF.smooth z).const_smul c
  inj z y hy := by
    have h : fderiv ℝ ((fun p => c • F p) ∘ chartInv ι z) y = c • fderiv ℝ (F ∘ chartInv ι z) y :=
      rs_fderiv_smul c (F ∘ chartInv ι z) y
    rw [h]
    intro u w huw
    apply hF.inj z y hy
    have : c • fderiv ℝ (F ∘ chartInv ι z) y u = c • fderiv ℝ (F ∘ chartInv ι z) y w := huw
    exact smul_right_injective V hc this
  lag z y hy u w := by
    have h : fderiv ℝ ((fun p => c • F p) ∘ chartInv ι z) y = c • fderiv ℝ (F ∘ chartInv ι z) y :=
      rs_fderiv_smul c (F ∘ chartInv ι z) y
    rw [h]
    simp only [smul_apply, map_smul, real_inner_smul_left,
      real_inner_smul_right, hF.lag z y hy u w, mul_zero]

theorem rs_image_eq_preimage {c : ℝ} (hc : 0 < c) (I : Set ℝ) :
    (fun s => c ^ 4 * s) '' I = (fun s => (c ^ 4)⁻¹ * s) ⁻¹' I := by
  have h4 : c ^ 4 ≠ 0 := (pow_pos hc 4).ne'
  ext s
  constructor
  · rintro ⟨t, ht, rfl⟩
    simpa [h4] using ht
  · intro hs
    exact ⟨_, hs, by field_simp⟩

theorem IsFlowM.rescale {I : Set ℝ} {Φ : ℝ → M → V} (hΦ : IsFlowM ι J I Φ) {c : ℝ} (hc : 0 < c) :
    IsFlowM ι J ((fun s => c ^ 4 * s) '' I) (fun s p => c • Φ ((c ^ 4)⁻¹ * s) p) := by
  have hc0 : c ≠ 0 := hc.ne'
  have h4 : c ^ 4 ≠ 0 := (pow_pos hc 4).ne'
  rw [rs_image_eq_preimage hc]
  refine ⟨hΦ.isOpen.preimage (continuous_const.mul continuous_id), ?_, ?_, ?_⟩
  · intro z
    have hA : ContDiff ℝ ∞ (fun q : ℝ × EuclideanSpace ℝ (Fin n) => ((c ^ 4)⁻¹ * q.1, q.2)) :=
      (contDiff_const.mul contDiff_fst).prodMk contDiff_snd
    have h1 := (hΦ.smooth z).comp hA.contDiffOn (s := ((fun s => (c ^ 4)⁻¹ * s) ⁻¹' I) ×ˢ chartDomE ι z)
      (by
        rintro ⟨s, y⟩ ⟨hs, hy⟩
        exact ⟨hs, hy⟩)
    exact h1.const_smul c
  · intro t ht
    exact (hΦ.imm _ ht).smul ι hc0
  · intro t ht z y hy
    have hv := hΦ.vel _ ht z y hy
    have hlin : HasDerivAt (fun s : ℝ => (c ^ 4)⁻¹ * s) ((c ^ 4)⁻¹) t := by
      simpa using (hasDerivAt_id t).const_mul ((c ^ 4)⁻¹)
    have hcomp := (hv.scomp t hlin).const_smul c
    have hX := hΦ.imm _ ht
    have hvel := velX_smul (J := J) hc0 (im_isOpen_chartDomE ι z) (hX.smooth z) (hX.inj z) hy
    beta_reduce at hvel hcomp
    have hcomp' : HasDerivAt (fun s => c • Φ ((c ^ 4)⁻¹ * s) (chartInv ι z y))
        (c • (c ^ 4)⁻¹ • velX J (Φ ((c ^ 4)⁻¹ * t) ∘ chartInv ι z) y) t := hcomp
    have e : (fun p => c • Φ ((c ^ 4)⁻¹ * t) p) ∘ chartInv ι z =
        fun y => c • (Φ ((c ^ 4)⁻¹ * t) ∘ chartInv ι z) y := rfl
    rw [e, hvel]
    convert hcomp' using 1
    rw [smul_smul]
    congr 1
    ring

end Flow

/-- info: 'CWFlow.Extension.potX_smul' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms potX_smul

/-- info: 'CWFlow.Extension.velX_smul' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms velX_smul

/-- info: 'CWFlow.Extension.aNormSq_smul' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms aNormSq_smul

/-- info: 'CWFlow.Extension.IsFlowM.rescale' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms IsFlowM.rescale

end CWFlow.Extension
