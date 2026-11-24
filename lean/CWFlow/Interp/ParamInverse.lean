import Mathlib.Analysis.Calculus.InverseFunctionTheorem.ContDiff
import Mathlib.Analysis.Calculus.InverseFunctionTheorem.ApproximatesLinearOn
import Mathlib.Analysis.Calculus.MeanValue
import Mathlib.Analysis.InnerProductSpace.PiL2
import Mathlib.Topology.Compactness.Compact

/-!
# A quantitative inverse function theorem with a time parameter

For a smooth family `χ_t : ℝⁿ → ℝⁿ` on `(a, b) × B(0, R)` with `χ_{t₀} = id`, there is `δ > 0` such
that for `|t - t₀| < δ` the map `χ_t` is injective on `B(0, R')`, has derivative within `1/4` of the
identity, and has an inverse `ξ(t, ·)` on `B(0, r)` (`r < R' < R`) which is jointly smooth in `(t, x)`.
-/

namespace CWFlow.Interp

open Set Metric Filter
open scoped ContDiff Topology NNReal

section Helpers

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]

/-- The slice `y ↦ χ (t, y)` has the partial derivative as derivative. -/
theorem pi_slice_hasFDerivAt {χ : ℝ × E → E} {t : ℝ} {y : E} (h : DifferentiableAt ℝ χ (t, y)) :
    HasFDerivAt (fun y' => χ (t, y'))
      ((fderiv ℝ χ (t, y)).comp (ContinuousLinearMap.inr ℝ ℝ E)) y :=
  HasFDerivAt.comp y (g := χ) (g' := fderiv ℝ χ (t, y)) h.hasFDerivAt
    (hasFDerivAt_prodMk_right t y)

/-- The derivative of a slice `y ↦ χ (t, y)` is the partial derivative. -/
theorem pi_fderiv_slice {χ : ℝ × E → E} {t : ℝ} {y : E} (h : DifferentiableAt ℝ χ (t, y)) :
    fderiv ℝ (fun y' => χ (t, y')) y =
      (fderiv ℝ χ (t, y)).comp (ContinuousLinearMap.inr ℝ ℝ E) :=
  (pi_slice_hasFDerivAt h).fderiv

/-- A linear map within `1/4` of the identity is injective. -/
theorem pi_inj_of_le {D : E →L[ℝ] E} (h : ‖D - ContinuousLinearMap.id ℝ E‖ ≤ 1 / 4) :
    Function.Injective D := by
  rw [injective_iff_map_eq_zero]
  intro v hv
  have h1 : ‖(D - ContinuousLinearMap.id ℝ E) v‖ ≤ 1 / 4 * ‖v‖ :=
    ((D - ContinuousLinearMap.id ℝ E).le_opNorm v).trans (by gcongr)
  have h2 : (D - ContinuousLinearMap.id ℝ E) v = -v := by simp [hv]
  rw [h2, norm_neg] at h1
  have : ‖v‖ = 0 := by linarith [norm_nonneg v]
  exact norm_eq_zero.1 this

/-- Mean value inequality: `g` approximates the identity on a convex set. -/
theorem pi_approx {g : E → E} {s : Set E} (hs : Convex ℝ s)
    (hd : ∀ y ∈ s, DifferentiableAt ℝ g y)
    (hb : ∀ y ∈ s, ‖fderiv ℝ g y - ContinuousLinearMap.id ℝ E‖ ≤ 1 / 4) :
    ∀ x ∈ s, ∀ y ∈ s, ‖g x - g y - (x - y)‖ ≤ 1 / 4 * ‖x - y‖ := by
  intro x hx y hy
  have h := Convex.norm_image_sub_le_of_norm_fderiv_le (𝕜 := ℝ) (f := fun z => g z - z) (C := 1 / 4)
    (s := s) (x := y) (y := x) (fun z hz => (hd z hz).sub differentiableAt_id)
    (fun z hz => by
      have hh : HasFDerivAt (fun z => g z - z) (fderiv ℝ g z - ContinuousLinearMap.id ℝ E) z :=
        (hd z hz).hasFDerivAt.sub (hasFDerivAt_id z)
      rw [hh.fderiv]
      exact hb z hz) hs hy hx
  have e : g x - g y - (x - y) = (g x - x) - (g y - y) := by abel
  rw [e]
  exact h

/-- A map approximating the identity on an open set covers the balls in the image. -/
theorem pi_surj [CompleteSpace E] {f : E → E} {s : Set E} (hs : IsOpen s)
    (hf : ∀ x ∈ s, ∀ y ∈ s, ‖f x - f y - (x - y)‖ ≤ 1 / 4 * ‖x - y‖) {b : E} {ρ : ℝ}
    (hρ : 0 ≤ ρ) (hsub : closedBall b ρ ⊆ s) :
    closedBall (f b) (3 / 4 * ρ) ⊆ f '' s := by
  by_cases hE : Subsingleton E
  · intro y _
    exact ⟨b, hsub (mem_closedBall_self hρ), Subsingleton.elim _ _⟩
  · have hnt : Nontrivial E := not_subsingleton_iff_nontrivial.1 hE
    have hN : ‖((ContinuousLinearEquiv.refl ℝ E).symm : E →L[ℝ] E)‖₊ = 1 := by
      simp
    have hA : ApproximatesLinearOn f ((ContinuousLinearEquiv.refl ℝ E : E ≃L[ℝ] E) : E →L[ℝ] E) s
        (1 / 4 : ℝ≥0) := by
      intro x hx y hy
      have := hf x hx y hy
      simpa using this
    have hc : (1 / 4 : ℝ≥0) < ‖((ContinuousLinearEquiv.refl ℝ E).symm : E →L[ℝ] E)‖₊⁻¹ := by
      rw [hN]
      rw [inv_one]
      exact_mod_cast (by norm_num : (1 / 4 : ℝ) < 1)
    have := hA.closedBall_subset_target (Or.inr hc) hs hρ hsub
    simp only [ApproximatesLinearOn.toOpenPartialHomeomorph_target] at this
    have h34 : (((‖((ContinuousLinearEquiv.refl ℝ E).symm : E →L[ℝ] E)‖₊)⁻¹ : ℝ≥0) : ℝ)
        - ((1 / 4 : ℝ≥0) : ℝ) = 3 / 4 := by
      rw [hN, inv_one]
      push_cast
      norm_num
    rw [h34] at this
    exact this

end Helpers

section Local

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]

/-- Smooth local inverse of `(s, y) ↦ (s, χ (s, y))`, giving a smooth local inverse in `y`. -/
theorem pi_local_inverse {χ : ℝ × E → E} {t : ℝ} {y : E} (hχ : ContDiffAt ℝ ∞ χ (t, y))
    (hD : Function.Injective (fderiv ℝ (fun y' => χ (t, y')) y)) :
    ∃ g : ℝ × E → E, ContDiffAt ℝ ∞ g (t, χ (t, y)) ∧ g (t, χ (t, y)) = y ∧
      ∀ᶠ z in 𝓝 (t, χ (t, y)), χ (z.1, g z) = z.2 := by
  have hdiff : DifferentiableAt ℝ χ (t, y) := hχ.differentiableAt (by simp)
  have hsl := pi_fderiv_slice hdiff
  set Ψ : ℝ × E → ℝ × E := fun z => (z.1, χ z) with hΨ
  have hΨc : ContDiffAt ℝ ∞ Ψ (t, y) := contDiffAt_fst.prodMk hχ
  set Ψ' : ℝ × E →L[ℝ] ℝ × E :=
    (ContinuousLinearMap.fst ℝ ℝ E).prod (fderiv ℝ χ (t, y)) with hΨ'
  have hΨd : HasFDerivAt Ψ Ψ' (t, y) :=
    (hasFDerivAt_fst).prodMk hdiff.hasFDerivAt
  have hinj : Function.Injective Ψ'.toLinearMap := by
    rw [← LinearMap.ker_eq_bot, LinearMap.ker_eq_bot']
    rintro ⟨σ, v⟩ h
    have h1 : σ = 0 := by
      have := congrArg Prod.fst h
      simpa [hΨ'] using this
    subst h1
    have h2 : fderiv ℝ χ (t, y) (0, v) = 0 := by
      have := congrArg Prod.snd h
      simpa [hΨ'] using this
    have h3 : fderiv ℝ (fun y' => χ (t, y')) y v = 0 := by
      rw [hsl]; simpa using h2
    have := hD (h3.trans (map_zero _).symm)
    simp [this]
  set e : (ℝ × E) ≃L[ℝ] (ℝ × E) :=
    (LinearEquiv.ofInjectiveEndo Ψ'.toLinearMap hinj).toContinuousLinearEquiv with he
  have hecoe : (e : ℝ × E →L[ℝ] ℝ × E) = Ψ' := by
    ext z <;> simp [he, LinearEquiv.coe_ofInjectiveEndo]
  have hΨe : HasFDerivAt Ψ (e : ℝ × E →L[ℝ] ℝ × E) (t, y) := by rw [hecoe]; exact hΨd
  have hn : (∞ : WithTop ℕ∞) ≠ 0 := by simp
  have hinv := hΨc.to_localInverse hΨe hn
  have hstrict := hΨc.hasStrictFDerivAt' hΨe hn
  refine ⟨fun z => (hΨc.localInverse hΨe hn z).2, ?_, ?_, ?_⟩
  · exact contDiffAt_snd.comp _ hinv
  · have := hΨc.localInverse_apply_image hΨe hn
    simp only [hΨ] at this ⊢
    rw [this]
  · have h := hstrict.eventually_right_inverse
    filter_upwards [h] with z hz
    have h1 := congrArg Prod.fst hz
    have h2 := congrArg Prod.snd hz
    simp only [hΨ] at h1 h2
    have h3 : ((z.1, (hΨc.localInverse hΨe hn z).2) : ℝ × E) = hΨc.localInverse hΨe hn z :=
      Prod.ext h1.symm rfl
    rw [h3]
    exact h2

end Local

section Main

variable {n : ℕ}

local notation "𝓔" => EuclideanSpace ℝ (Fin n)

/-- Quantitative inverse function theorem with a time parameter (cosmetic deviation: `n` is a
section variable and `𝓔` abbreviates `EuclideanSpace ℝ (Fin n)`). -/
theorem exists_param_inverse {a b t₀ : ℝ} (ht₀ : t₀ ∈ Ioo a b) {R : ℝ}
    {χ : ℝ × 𝓔 → 𝓔}
    (hχ : ContDiffOn ℝ ∞ χ (Ioo a b ×ˢ ball (0 : 𝓔) R))
    (hid : ∀ y ∈ ball (0 : 𝓔) R, χ (t₀, y) = y)
    {r R' : ℝ} (hr : 0 < r) (hrR' : r < R') (hR'R : R' < R) :
    ∃ δ > 0, Ioo (t₀ - δ) (t₀ + δ) ⊆ Ioo a b ∧
      ∃ ξ : ℝ × 𝓔 → 𝓔,
        ContDiffOn ℝ ∞ ξ (Ioo (t₀ - δ) (t₀ + δ) ×ˢ ball (0 : 𝓔) r) ∧
        ∀ t ∈ Ioo (t₀ - δ) (t₀ + δ),
          InjOn (fun y => χ (t, y)) (ball (0 : 𝓔) R') ∧
          (∀ y ∈ ball (0 : 𝓔) R',
            Function.Injective (fderiv ℝ (fun y' => χ (t, y')) y) ∧
            ‖fderiv ℝ (fun y' => χ (t, y')) y - ContinuousLinearMap.id ℝ _‖ ≤ 1 / 4) ∧
          (∀ x ∈ ball (0 : 𝓔) r,
            ξ (t, x) ∈ ball (0 : 𝓔) R' ∧ χ (t, ξ (t, x)) = x) ∧
          (∀ y ∈ ball (0 : 𝓔) R', χ (t, y) ∈ ball (0 : 𝓔) r → ξ (t, χ (t, y)) = y) := by
  classical
  have hU : IsOpen (Ioo a b ×ˢ ball (0 : 𝓔) R) := isOpen_Ioo.prod isOpen_ball
  have hχc : ∀ p ∈ Ioo a b ×ˢ ball (0 : 𝓔) R, ContDiffAt ℝ ∞ χ p :=
    fun p hp => hχ.contDiffAt (hU.mem_nhds hp)
  have hdiff : ∀ p ∈ Ioo a b ×ˢ ball (0 : 𝓔) R, DifferentiableAt ℝ χ p :=
    fun p hp => (hχc p hp).differentiableAt (by simp)
  have hGc : ContinuousOn
      (fun z : ℝ × 𝓔 => (fderiv ℝ χ z).comp (ContinuousLinearMap.inr ℝ ℝ 𝓔))
      (Ioo a b ×ˢ ball (0 : 𝓔) R) :=
    (hχ.continuousOn_fderiv_of_isOpen hU (by simp)).clm_comp continuousOn_const
  have hD0 : ∀ y ∈ ball (0 : 𝓔) R,
      (fderiv ℝ χ (t₀, y)).comp (ContinuousLinearMap.inr ℝ ℝ 𝓔) = ContinuousLinearMap.id ℝ 𝓔 := by
    intro y hy
    rw [← pi_fderiv_slice (hdiff _ ⟨ht₀, hy⟩)]
    have : (fun y' => χ (t₀, y')) =ᶠ[𝓝 y] id := by
      filter_upwards [isOpen_ball.mem_nhds hy] with y' hy' using hid y' hy'
    rw [this.fderiv_eq]
    simp
  -- uniform smallness
  have hmem : ∀ᶠ t in 𝓝 t₀, t ∈ Ioo a b := isOpen_Ioo.mem_nhds ht₀
  have hev : ∀ᶠ t in 𝓝 t₀, t ∈ Ioo a b ∧ ∀ y ∈ closedBall (0 : 𝓔) R',
      ‖(fderiv ℝ χ (t, y)).comp (ContinuousLinearMap.inr ℝ ℝ 𝓔) - ContinuousLinearMap.id ℝ 𝓔‖
        < 1 / 4 ∧ ‖χ (t, y) - y‖ < (R' - r) / 4 := by
    refine hmem.and ((isCompact_closedBall (0 : 𝓔) R').eventually_forall_of_forall_eventually ?_)
    intro y hy
    have hyR : y ∈ ball (0 : 𝓔) R := by
      rw [mem_closedBall_zero_iff] at hy
      rw [mem_ball_zero_iff]; linarith
    have hnh : Ioo a b ×ˢ ball (0 : 𝓔) R ∈ 𝓝 (t₀, y) := hU.mem_nhds ⟨ht₀, hyR⟩
    have hc1 : ContinuousAt (fun z : ℝ × 𝓔 =>
        (fderiv ℝ χ z).comp (ContinuousLinearMap.inr ℝ ℝ 𝓔) - ContinuousLinearMap.id ℝ 𝓔)
        (t₀, y) := (hGc.continuousAt hnh).sub continuousAt_const
    have hc2 : ContinuousAt (fun z : ℝ × 𝓔 => χ z - z.2) (t₀, y) :=
      (hχ.continuousOn.continuousAt hnh).sub continuousAt_snd
    have h1 : ‖(fderiv ℝ χ (t₀, y)).comp (ContinuousLinearMap.inr ℝ ℝ 𝓔)
        - ContinuousLinearMap.id ℝ 𝓔‖ < 1 / 4 := by
      rw [hD0 y hyR]; simp
    have h2 : ‖χ (t₀, y) - y‖ < (R' - r) / 4 := by
      rw [hid y hyR]; simp; linarith
    have e1 := hc1.norm.tendsto.eventually (gt_mem_nhds h1)
    have e2 := hc2.norm.tendsto.eventually (gt_mem_nhds h2)
    filter_upwards [e1, e2] with z hz1 hz2
    exact ⟨hz1, hz2⟩
  obtain ⟨ε, hε, hεP⟩ := Metric.eventually_nhds_iff.1 hev
  have hI : ∀ t ∈ Ioo (t₀ - ε) (t₀ + ε), t ∈ Ioo a b ∧ ∀ y ∈ closedBall (0 : 𝓔) R',
      ‖(fderiv ℝ χ (t, y)).comp (ContinuousLinearMap.inr ℝ ℝ 𝓔) - ContinuousLinearMap.id ℝ 𝓔‖
        < 1 / 4 ∧ ‖χ (t, y) - y‖ < (R' - r) / 4 := by
    intro t ht
    apply hεP
    rw [Real.dist_eq, abs_lt]
    constructor <;> linarith [ht.1, ht.2]
  -- the per-time facts
  have per : ∀ t ∈ Ioo (t₀ - ε) (t₀ + ε),
      InjOn (fun y => χ (t, y)) (ball (0 : 𝓔) R') ∧
      (∀ y ∈ ball (0 : 𝓔) R',
        ‖fderiv ℝ (fun y' => χ (t, y')) y - ContinuousLinearMap.id ℝ 𝓔‖ ≤ 1 / 4) ∧
      (∀ x ∈ ball (0 : 𝓔) r, ∃ y ∈ ball (0 : 𝓔) R', χ (t, y) = x) := by
    intro t ht
    obtain ⟨htab, hsm⟩ := hI t ht
    have hyR : ∀ y ∈ ball (0 : 𝓔) R', y ∈ ball (0 : 𝓔) R :=
      fun y hy => ball_subset_ball hR'R.le hy
    have hdy : ∀ y ∈ ball (0 : 𝓔) R', HasFDerivAt (fun y' => χ (t, y'))
        ((fderiv ℝ χ (t, y)).comp (ContinuousLinearMap.inr ℝ ℝ 𝓔)) y :=
      fun y hy => pi_slice_hasFDerivAt (hdiff _ ⟨htab, hyR y hy⟩)
    have hDle : ∀ y ∈ ball (0 : 𝓔) R',
        ‖fderiv ℝ (fun y' => χ (t, y')) y - ContinuousLinearMap.id ℝ 𝓔‖ ≤ 1 / 4 := by
      intro y hy
      rw [(hdy y hy).fderiv]
      exact (hsm y (ball_subset_closedBall hy)).1.le
    have happ := pi_approx (g := fun y' => χ (t, y')) (convex_ball (0 : 𝓔) R')
      (fun y hy => (hdy y hy).differentiableAt) hDle
    refine ⟨?_, hDle, ?_⟩
    · intro x hx y hy hxy
      have h := happ x hx y hy
      have h' : ‖x - y‖ ≤ 1 / 4 * ‖x - y‖ := by
        simp only at hxy
        rw [hxy, sub_self, zero_sub, norm_neg] at h
        exact h
      have : ‖x - y‖ = 0 := by linarith [norm_nonneg (x - y)]
      exact sub_eq_zero.1 (norm_eq_zero.1 this)
    · intro x hx
      have hxR' : x ∈ ball (0 : 𝓔) R' := ball_subset_ball hrR'.le hx
      have hsub : closedBall x ((R' - r) / 2) ⊆ ball (0 : 𝓔) R' := by
        intro z hz
        rw [mem_closedBall, dist_eq_norm] at hz
        rw [mem_ball_zero_iff]
        rw [mem_ball_zero_iff] at hx
        have : ‖z‖ ≤ ‖x‖ + ‖z - x‖ := by
          calc ‖z‖ = ‖x + (z - x)‖ := by rw [add_sub_cancel]
            _ ≤ ‖x‖ + ‖z - x‖ := norm_add_le _ _
        linarith
      have hs := pi_surj (f := fun y' => χ (t, y')) isOpen_ball happ (b := x)
        (ρ := (R' - r) / 2) (by linarith) hsub
      have hxm : x ∈ closedBall (χ (t, x)) (3 / 4 * ((R' - r) / 2)) := by
        rw [mem_closedBall, dist_comm, dist_eq_norm]
        have := (hsm x (ball_subset_closedBall hxR')).2
        linarith
      exact hs hxm |> fun ⟨y, hy, hyx⟩ => ⟨y, hy, hyx⟩
  refine ⟨ε, hε, fun t ht => (hI t ht).1, ?_⟩
  set ξ : ℝ × 𝓔 → 𝓔 := fun p => Function.invFunOn (fun y => χ (p.1, y)) (ball (0 : 𝓔) R') p.2
    with hξ
  have hpos : ∀ t ∈ Ioo (t₀ - ε) (t₀ + ε), ∀ x ∈ ball (0 : 𝓔) r,
      ξ (t, x) ∈ ball (0 : 𝓔) R' ∧ χ (t, ξ (t, x)) = x := fun t ht x hx =>
    Function.invFunOn_pos ((per t ht).2.2 x hx)
  have hleft : ∀ t ∈ Ioo (t₀ - ε) (t₀ + ε), ∀ y ∈ ball (0 : 𝓔) R', ξ (t, χ (t, y)) = y :=
    fun t ht y hy => (per t ht).1.leftInvOn_invFunOn hy
  refine ⟨ξ, ?_, fun t ht => ⟨(per t ht).1, fun y hy => ⟨?_, (per t ht).2.1 y hy⟩,
    hpos t ht, fun y hy _ => hleft t ht y hy⟩⟩
  · -- smoothness
    intro p hp
    obtain ⟨t, x⟩ := p
    obtain ⟨ht, hx⟩ := hp
    obtain ⟨y, hy, hχy⟩ := (per t ht).2.2 x hx
    have htab : t ∈ Ioo a b := (hI t ht).1
    obtain ⟨g, hg, hgy, hginv⟩ := pi_local_inverse
      (hχc (t, y) ⟨htab, ball_subset_ball hR'R.le hy⟩)
      (pi_inj_of_le ((per t ht).2.1 y hy))
    rw [hχy] at hg hgy hginv
    have hnb : ∀ᶠ z in 𝓝 (t, x), z ∈ Ioo (t₀ - ε) (t₀ + ε) ×ˢ ball (0 : 𝓔) r :=
      (isOpen_Ioo.prod isOpen_ball).mem_nhds ⟨ht, hx⟩
    have hgb : ∀ᶠ z in 𝓝 (t, x), g z ∈ ball (0 : 𝓔) R' :=
      hg.continuousAt.eventually_mem (isOpen_ball.mem_nhds (by rw [hgy]; exact hy))
    have hxi : ξ (t, x) = y := by
      rw [← hχy]; exact hleft t ht y hy
    have hat : ContDiffAt ℝ ∞ ξ (t, x) := by
      refine hg.congr_of_eventuallyEq ?_
      filter_upwards [hnb, hginv, hgb] with z hz hz2 hz3
      have h1 := hleft z.1 hz.1 (g z) hz3
      rw [hz2] at h1
      exact h1
    exact hat.contDiffWithinAt
  · -- unused slot: injectivity of the derivative
    exact pi_inj_of_le ((per t ht).2.1 y hy)

end Main

/-- info: 'CWFlow.Interp.exists_param_inverse' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms exists_param_inverse

end CWFlow.Interp
