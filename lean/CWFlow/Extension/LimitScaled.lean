import CWFlow.Extension.LimitM
import CWFlow.Extension.EnergyBound
import CWFlow.Extension.Rescale

/-!
# The extension theorem without smallness hypotheses

`exists_limit_flow` needs the second fundamental form to be small (`≤ c₁ = 1/1000`).  A flow with
`|A|² ≤ K²` becomes small after the parabolic rescaling `Φ̃(s) = λ Φ(λ⁻⁴ s)`, `λ ≥ |K| / c₁`
(`aNormSq ↦ λ⁻² aNormSq`).  The energy bounds on `[t̃₁, λ⁴ T)` come from `energyM_bounded`, and the
extension is scaled back.  Result: `exists_limit_flow_of_bound`.
-/

open Metric Set Filter
open scoped ContDiff Manifold Topology RealInnerProductSpace

namespace CWFlow.Extension

theorem ls_isFlowM_mono {n : ℕ} {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    [FiniteDimensional ℝ E] {M : Type*} [TopologicalSpace M] [ChartedSpace E M]
    [IsManifold 𝓘(ℝ, E) ∞ M]
    {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
    (ι : EuclideanSpace ℝ (Fin n) ≃L[ℝ] E) {J : V →L[ℝ] V} {I I' : Set ℝ} {Φ : ℝ → M → V}
    (hΦ : IsFlowM ι J I Φ) (hI' : IsOpen I') (hsub : I' ⊆ I) : IsFlowM ι J I' Φ where
  isOpen := hI'
  smooth z := (hΦ.smooth z).mono (prod_mono hsub subset_rfl)
  imm t ht := hΦ.imm t (hsub ht)
  vel t ht z y hy := hΦ.vel t (hsub ht) z y hy

theorem ls_image_Ioo {L : ℝ} (hL : 0 < L) (a T : ℝ) :
    (fun s => L ^ 4 * s) '' Ioo a T = Ioo (L ^ 4 * a) (L ^ 4 * T) := by
  have h4 : 0 < L ^ 4 := pow_pos hL 4
  ext s
  rw [rs_image_eq_preimage hL]
  simp only [mem_preimage, mem_Ioo]
  rw [lt_inv_mul_iff₀ h4, inv_mul_lt_iff₀ h4]

theorem ls_consts :
    6 * ((1 / 1000 : ℝ) * (1 + (1 / 100 : ℝ) ^ 2) ^ ((3 : ℝ) / 2)) < 1 / 100 ∧
    (1 + (1 / 100 : ℝ) ^ 2) ^ ((3 : ℝ) / 2) * (1 / 1000) ≤ 1 := by
  have h1 : (1 : ℝ) ≤ 1 + (1 / 100 : ℝ) ^ 2 := by nlinarith [sq_nonneg (1 / 100 : ℝ)]
  have hr : (1 + (1 / 100 : ℝ) ^ 2) ^ ((3 : ℝ) / 2) ≤ 11 / 10 := by
    calc (1 + (1 / 100 : ℝ) ^ 2) ^ ((3 : ℝ) / 2) ≤ (1 + (1 / 100 : ℝ) ^ 2) ^ (2 : ℝ) :=
          Real.rpow_le_rpow_of_exponent_le h1 (by norm_num)
      _ = (1 + (1 / 100 : ℝ) ^ 2) ^ 2 := Real.rpow_two _
      _ ≤ 11 / 10 := by norm_num
  constructor <;> nlinarith

theorem exists_limit_flow_of_bound {n : ℕ} {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    [FiniteDimensional ℝ E] {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
    [T2Space M] [CompactSpace M]
    {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V] [FiniteDimensional ℝ V]
    (ι : EuclideanSpace ℝ (Fin n) ≃L[ℝ] E) {J : V →L[ℝ] V}
    (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    (hdim : Module.finrank ℝ V = 2 * n)
    {a T : ℝ} (haT : a < T) {Φ : ℝ → M → V} (hΦ : IsFlowM ι J (Ioo a T) Φ) {K : ℝ}
    (hK : ∀ t ∈ Ioo a T, ∀ z : M, ∀ y ∈ chartDomE ι z, aNormSq J (Φ t ∘ chartInv ι z) y ≤ K ^ 2) :
    ∃ t₁ ∈ Ioo a T, ∃ Φbar : ℝ → M → V,
      (∀ t ∈ Ioo a T, Φbar t = Φ t) ∧
      (∀ z : M, ContDiffOn ℝ ∞ (fun q : ℝ × EuclideanSpace ℝ (Fin n) => Φbar q.1 (chartInv ι z q.2))
        (Icc t₁ T ×ˢ chartDomE ι z)) ∧
      (∀ t ∈ Icc t₁ T, IsLagImm ι J (Φbar t)) ∧
      (∀ t ∈ Icc t₁ T, ∀ z : M, ∀ y ∈ chartDomE ι z,
        HasDerivWithinAt (fun s => Φbar s (chartInv ι z y)) (velX J (Φbar t ∘ chartInv ι z) y)
          (Icc t₁ T) t) := by
  classical
  obtain ⟨hsmall, hK₃⟩ := ls_consts
  have hc₀ : (0 : ℝ) < 1 / 100 := by norm_num
  have hc₁ : (0 : ℝ) ≤ 1 / 1000 := by norm_num
  obtain ⟨A⟩ := exists_intAtlas (E := E) (M := M)
  set L : ℝ := max 1 (|K| / (1 / 1000)) with hLdef
  have hL1 : 1 ≤ L := le_max_left _ _
  have hL : 0 < L := lt_of_lt_of_le one_pos hL1
  have hL0 : L ≠ 0 := hL.ne'
  have h4 : 0 < L ^ 4 := pow_pos hL 4
  have hKL : |K| ≤ (1 / 1000) * L := by
    have := le_max_right 1 (|K| / (1 / 1000 : ℝ))
    rw [← hLdef, div_le_iff₀ (by norm_num)] at this
    linarith
  have himg := ls_image_Ioo hL a T
  have hΨ : IsFlowM ι J (Ioo (L ^ 4 * a) (L ^ 4 * T)) (fun s p => L • Φ ((L ^ 4)⁻¹ * s) p) := by
    have := hΦ.rescale ι hL
    rwa [himg] at this
  have hAΨ : ∀ s ∈ Ioo (L ^ 4 * a) (L ^ 4 * T), ∀ z : M, ∀ y ∈ chartDomE ι z,
      aNormSq J ((fun s p => L • Φ ((L ^ 4)⁻¹ * s) p) s ∘ chartInv ι z) y ≤ (1 / 1000 : ℝ) ^ 2 := by
    intro s hs z y hy
    have hs' : (L ^ 4)⁻¹ * s ∈ Ioo a T := by
      rw [← himg, rs_image_eq_preimage hL] at hs
      exact hs
    have hX := hΦ.imm _ hs'
    have e : ((fun s p => L • Φ ((L ^ 4)⁻¹ * s) p) s ∘ chartInv ι z) =
        fun y => L • (Φ ((L ^ 4)⁻¹ * s) ∘ chartInv ι z) y := rfl
    rw [e, aNormSq_smul hL0 _ y (hX.inj z y hy)]
    have h1 := hK _ hs' z y hy
    have hK2 : K ^ 2 ≤ ((1 / 1000 : ℝ) * L) ^ 2 := by
      rw [← sq_abs K]
      exact pow_le_pow_left₀ (abs_nonneg _) hKL 2
    rw [inv_mul_le_iff₀ (pow_pos hL 2)]
    nlinarith
  have hlt : L ^ 4 * a < L ^ 4 * T := mul_lt_mul_of_pos_left haT h4
  set t₁t : ℝ := (L ^ 4 * a + L ^ 4 * T) / 2 with ht₁t
  have ht₁tm : t₁t ∈ Ioo (L ^ 4 * a) (L ^ 4 * T) := ⟨by rw [ht₁t]; linarith, by rw [ht₁t]; linarith⟩
  have hB := energyM_bounded ι A hJi hJJ hdim hc₀ (le_refl _) hc₁ hsmall hK₃ hΨ hAΨ ht₁tm
  choose B hB using hB
  have hsub : Ioo t₁t (L ^ 4 * T) ⊆ Ioo (L ^ 4 * a) (L ^ 4 * T) := Ioo_subset_Ioo ht₁tm.1.le le_rfl
  have hΨ' := ls_isFlowM_mono ι hΨ isOpen_Ioo hsub
  obtain ⟨t₁', ht₁', Ψbar, hag, hsm, hlag, hvel⟩ := exists_limit_flow ι A hJi hJJ hdim hc₀ (le_refl _)
    hc₁ hsmall hK₃ B ht₁tm.2 hΨ' (fun s hs => hAΨ s (hsub hs))
    (fun s hs l => hB l s ⟨hs.1.le, hs.2⟩)
  -- scale back
  have ht₁'a : L ^ 4 * a < t₁' := lt_trans ht₁tm.1 ht₁'.1
  set t₁ : ℝ := (L ^ 4)⁻¹ * t₁' with ht₁def
  have hLt₁ : L ^ 4 * t₁ = t₁' := mul_inv_cancel_left₀ h4.ne' _
  have ht₁m : t₁ ∈ Ioo a T := by
    rw [ht₁def]
    exact ⟨(lt_inv_mul_iff₀ h4).2 ht₁'a, (inv_mul_lt_iff₀ h4).2 ht₁'.2⟩
  have hmaps : MapsTo (fun s => L ^ 4 * s) (Icc t₁ T) (Icc t₁' (L ^ 4 * T)) := by
    intro s hs
    refine ⟨?_, mul_le_mul_of_nonneg_left hs.2 h4.le⟩
    rw [← hLt₁]
    exact mul_le_mul_of_nonneg_left hs.1 h4.le
  let Φbar : ℝ → M → V := fun t p => if t < T then Φ t p else L⁻¹ • Ψbar (L ^ 4 * T) p
  have hEq : ∀ t ∈ Icc t₁ T, Φbar t = fun p => L⁻¹ • Ψbar (L ^ 4 * t) p := by
    intro t ht
    by_cases htT : t < T
    · have hmem : L ^ 4 * t ∈ Ioo t₁t (L ^ 4 * T) := by
        refine ⟨lt_of_lt_of_le ht₁'.1 ?_, mul_lt_mul_of_pos_left htT h4⟩
        rw [← hLt₁]
        exact mul_le_mul_of_nonneg_left ht.1 h4.le
      funext p
      simp only [Φbar, if_pos htT]
      rw [hag _ hmem]
      show Φ t p = L⁻¹ • L • Φ ((L ^ 4)⁻¹ * (L ^ 4 * t)) p
      rw [inv_mul_cancel_left₀ h4.ne', smul_smul, inv_mul_cancel₀ hL0, one_smul]
    · have : t = T := le_antisymm ht.2 (not_lt.1 htT)
      subst this
      funext p
      simp only [Φbar, if_neg htT]
  refine ⟨t₁, ht₁m, Φbar, ?_, ?_, ?_, ?_⟩
  · intro t ht
    funext p
    simp only [Φbar, if_pos ht.2]
  · intro z
    have hA' : ContDiff ℝ ∞ (fun q : ℝ × EuclideanSpace ℝ (Fin n) => (L ^ 4 * q.1, q.2)) :=
      (contDiff_const.mul contDiff_fst).prodMk contDiff_snd
    have h1 := ((hsm z).comp hA'.contDiffOn (s := Icc t₁ T ×ˢ chartDomE ι z)
      (by rintro ⟨s, y⟩ ⟨hs, hy⟩; exact ⟨hmaps hs, hy⟩)).const_smul L⁻¹
    refine h1.congr ?_
    rintro ⟨t, y⟩ ⟨ht, hy⟩
    simpa using congrFun (hEq t ht) (chartInv ι z y)
  · intro t ht
    rw [hEq t ht]
    exact (hlag _ (hmaps ht)).smul ι (inv_ne_zero hL0)
  · intro t ht z y hy
    have hv := hvel (L ^ 4 * t) (hmaps ht) z y hy
    have hlin : HasDerivWithinAt (fun s : ℝ => L ^ 4 * s) (L ^ 4) (Icc t₁ T) t := by
      simpa using ((hasDerivAt_id t).const_mul (L ^ 4)).hasDerivWithinAt
    have hcomp := (hv.scomp t hlin hmaps).const_smul L⁻¹
    have hX := hlag (L ^ 4 * t) (hmaps ht)
    have hvelX := velX_smul (J := J) (inv_ne_zero hL0) (im_isOpen_chartDomE ι z) (hX.smooth z)
      (hX.inj z) hy
    have e : Φbar t ∘ chartInv ι z = fun y => L⁻¹ • (Ψbar (L ^ 4 * t) ∘ chartInv ι z) y := by
      rw [hEq t ht]; rfl
    rw [e, hvelX]
    refine HasDerivWithinAt.congr_of_mem (f := fun s => L⁻¹ • Ψbar (L ^ 4 * s) (chartInv ι z y))
      ?_ ?_ ht
    · have hval : ∀ w : V, ((L⁻¹ ^ 2)⁻¹ * (L⁻¹ ^ 2)⁻¹ * L⁻¹) • w = L⁻¹ • L ^ 4 • w := by
        intro w
        rw [smul_smul]
        congr 1
        field_simp
      rw [hval]
      exact hcomp
    · intro s hs
      simpa using congrFun (hEq s hs) (chartInv ι z y)

/-- info: 'CWFlow.Extension.exists_limit_flow_of_bound' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_limit_flow_of_bound

end CWFlow.Extension
