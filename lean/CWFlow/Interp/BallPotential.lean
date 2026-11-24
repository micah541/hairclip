import Mathlib.Analysis.InnerProductSpace.Calculus
import Mathlib.Analysis.InnerProductSpace.PiL2
import Mathlib.Analysis.Calculus.ParametricIntervalIntegral
import Mathlib.Analysis.Calculus.MeanValue
import Mathlib.Analysis.Calculus.ContDiff.Operations
import Mathlib.Analysis.Normed.Module.FiniteDimension
import Mathlib.MeasureTheory.Integral.IntervalIntegral.FundThmCalculus

/-!
# The Poincaré lemma for 1-forms on a ball

A smooth vector field with symmetric derivative on a ball of a finite-dimensional real inner product
space is a gradient (`exists_potential_ball`), the potential `φ x = ∫₀¹ ⟪s (x₀ + t (x - x₀)), x - x₀⟫ dt`
being unique up to the value at the centre (`eqOn_ball_of_fderiv_eq`). Also a derivative bound along
segments from the centre (`norm_lt_of_fderiv_le_ball`).
-/

open Metric Set MeasureTheory

open scoped ContDiff RealInnerProductSpace Topology

noncomputable section

namespace CWFlow.Interp

variable {F : Type*} [NormedAddCommGroup F] [InnerProductSpace ℝ F]

/-- The `z`-derivative of the integrand `⟪s (x₀ + t (z - x₀)), z - x₀⟫`. -/
def bp_D (s : F → F) (x₀ z : F) (t : ℝ) : F →L[ℝ] ℝ :=
  t • ((innerSL ℝ (z - x₀)).comp (fderiv ℝ s (x₀ + t • (z - x₀)))) +
    innerSL ℝ (s (x₀ + t • (z - x₀)))

lemma bp_mem_ball {x₀ z : F} {R t : ℝ} (hz : z ∈ ball x₀ R) (ht : t ∈ Icc (0 : ℝ) 1) :
    x₀ + t • (z - x₀) ∈ ball x₀ R := by
  rw [mem_ball_iff_norm] at hz ⊢
  rw [add_sub_cancel_left, norm_smul, Real.norm_of_nonneg ht.1]
  calc t * ‖z - x₀‖ ≤ 1 * ‖z - x₀‖ := by gcongr; exact ht.2
    _ < R := by rwa [one_mul]

lemma bp_hasFDerivAt_integrand {x₀ z : F} {R t : ℝ} {s : F → F}
    (hs : DifferentiableOn ℝ s (ball x₀ R)) (hz : z ∈ ball x₀ R) (ht : t ∈ Icc (0 : ℝ) 1) :
    HasFDerivAt (fun w => inner ℝ (s (x₀ + t • (w - x₀))) (w - x₀)) (bp_D s x₀ z t) z := by
  have hp := bp_mem_ball hz ht
  have hsd : HasFDerivAt s (fderiv ℝ s (x₀ + t • (z - x₀))) (x₀ + t • (z - x₀)) :=
    (hs.differentiableAt (isOpen_ball.mem_nhds hp)).hasFDerivAt
  have hg : HasFDerivAt (fun w : F => x₀ + t • (w - x₀)) (t • ContinuousLinearMap.id ℝ F) z :=
    (((hasFDerivAt_id z).sub_const x₀).const_smul t).const_add x₀
  have h1 := hsd.comp z hg
  have h2 : HasFDerivAt (fun w : F => w - x₀) (ContinuousLinearMap.id ℝ F) z :=
    (hasFDerivAt_id z).sub_const x₀
  have h := h1.inner ℝ h2
  refine h.congr_fderiv ?_
  ext a
  simp [bp_D, fderivInnerCLM_apply, real_inner_comm, inner_smul_right, inner_sub_left]
  ring

lemma bp_path_cont (x₀ z : F) : Continuous (fun t : ℝ => x₀ + t • (z - x₀)) := by fun_prop

lemma bp_slice_cont {x₀ z : F} {R : ℝ} {s : F → F} (hs : ContDiffOn ℝ ∞ s (ball x₀ R))
    (hz : z ∈ ball x₀ R) :
    ContinuousOn (fun t : ℝ => inner ℝ (s (x₀ + t • (z - x₀))) (z - x₀)) (Icc 0 1) ∧
    ContinuousOn (fun t : ℝ => bp_D s x₀ z t) (Icc 0 1) := by
  have hmaps : MapsTo (fun t : ℝ => x₀ + t • (z - x₀)) (Icc 0 1) (ball x₀ R) :=
    fun t ht => bp_mem_ball hz ht
  have hc : ContinuousOn (fun t : ℝ => s (x₀ + t • (z - x₀))) (Icc 0 1) :=
    hs.continuousOn.comp (bp_path_cont x₀ z).continuousOn hmaps
  have hd : ContinuousOn (fun t : ℝ => fderiv ℝ s (x₀ + t • (z - x₀))) (Icc 0 1) :=
    (hs.continuousOn_fderiv_of_isOpen isOpen_ball (by simp)).comp
      (bp_path_cont x₀ z).continuousOn hmaps
  refine ⟨hc.inner continuousOn_const, ?_⟩
  unfold bp_D
  exact ((continuousOn_id.smul (continuousOn_const.clm_comp hd))).add
    ((innerSL ℝ).continuous.comp_continuousOn hc)

lemma bp_ftc {x₀ x : F} {R : ℝ} {s : F → F} (hs : ContDiffOn ℝ ∞ s (ball x₀ R))
    (hsym : ∀ x ∈ ball x₀ R, ∀ a b : F,
      inner ℝ (fderiv ℝ s x a) b = inner ℝ a (fderiv ℝ s x b))
    (hx : x ∈ ball x₀ R) (a : F) :
    ∫ t in (0 : ℝ)..1, bp_D s x₀ x t a = inner ℝ (s x) a := by
  have hsd : DifferentiableOn ℝ s (ball x₀ R) := hs.differentiableOn (by simp)
  have hder : ∀ t ∈ uIcc (0 : ℝ) 1, HasDerivAt (fun t : ℝ => t * inner ℝ (s (x₀ + t • (x - x₀))) a)
      (bp_D s x₀ x t a) t := by
    intro t ht
    rw [uIcc_of_le zero_le_one] at ht
    have hp := bp_mem_ball hx ht
    have hsf : HasFDerivAt s (fderiv ℝ s (x₀ + t • (x - x₀))) (x₀ + t • (x - x₀)) :=
      (hsd.differentiableAt (isOpen_ball.mem_nhds hp)).hasFDerivAt
    have hpath : HasDerivAt (fun t : ℝ => x₀ + t • (x - x₀)) (x - x₀) t := by
      simpa using ((hasDerivAt_id t).smul_const (x - x₀)).const_add x₀
    have h1 := hsf.comp_hasDerivAt t hpath
    have h2 := h1.inner ℝ (hasDerivAt_const t a)
    have h3 := (hasDerivAt_id t).mul h2
    refine h3.congr_deriv ?_
    have hs' := hsym _ hp (x - x₀) a
    simp only [bp_D, add_apply, smul_apply,
      ContinuousLinearMap.comp_apply, innerSL_apply_apply, id]
    simp only [Function.comp, inner_zero_right]
    rw [hs']
    simp only [smul_eq_mul]
    ring
  have hcont : ContinuousOn (fun t : ℝ => bp_D s x₀ x t a) (uIcc 0 1) := by
    rw [uIcc_of_le zero_le_one]
    exact (bp_slice_cont hs hx).2.clm_apply continuousOn_const
  have := intervalIntegral.integral_eq_sub_of_hasDerivAt hder (hcont.intervalIntegrable)
  rw [this]
  simp

lemma bp_hasFDerivAt_pot [FiniteDimensional ℝ F] {x₀ x : F} {R : ℝ} {s : F → F}
    (hs : ContDiffOn ℝ ∞ s (ball x₀ R))
    (hsym : ∀ x ∈ ball x₀ R, ∀ a b : F,
      inner ℝ (fderiv ℝ s x a) b = inner ℝ a (fderiv ℝ s x b))
    (hx : x ∈ ball x₀ R) :
    HasFDerivAt (fun y => ∫ t in (0 : ℝ)..1, inner ℝ (s (x₀ + t • (y - x₀))) (y - x₀))
      (innerSL ℝ (s x)) x := by
  have hsd : DifferentiableOn ℝ s (ball x₀ R) := hs.differentiableOn (by simp)
  have hxn : ‖x - x₀‖ < R := mem_ball_iff_norm.1 hx
  set r : ℝ := (R - ‖x - x₀‖) / 2 with hr
  have hr0 : 0 < r := by rw [hr]; linarith
  set R' : ℝ := ‖x - x₀‖ + r with hR'
  have hR'R : R' < R := by rw [hR', hr]; linarith
  have hR'0 : 0 ≤ R' := by rw [hR']; positivity
  have hKsub : closedBall x₀ R' ⊆ ball x₀ R := closedBall_subset_ball hR'R
  obtain ⟨Ms, hMs⟩ := (isCompact_closedBall x₀ R').exists_bound_of_continuousOn
    (hs.continuousOn.mono hKsub)
  obtain ⟨Md, hMd⟩ := (isCompact_closedBall x₀ R').exists_bound_of_continuousOn
    ((hs.continuousOn_fderiv_of_isOpen isOpen_ball (by simp)).mono hKsub)
  have hMd0 : 0 ≤ Md := (norm_nonneg _).trans (hMd x₀ (mem_closedBall_self hR'0))
  have hzR : ∀ z ∈ ball x r, ‖z - x₀‖ < R' := by
    intro z hz
    rw [mem_ball_iff_norm] at hz
    calc ‖z - x₀‖ = ‖(z - x) + (x - x₀)‖ := by rw [sub_add_sub_cancel]
      _ ≤ ‖z - x‖ + ‖x - x₀‖ := norm_add_le _ _
      _ < R' := by rw [hR']; linarith
  have hzball : ∀ z ∈ ball x r, z ∈ ball x₀ R := fun z hz => by
    rw [mem_ball_iff_norm]; exact (hzR z hz).trans hR'R
  have hbound : ∀ z ∈ ball x r, ∀ t ∈ Icc (0 : ℝ) 1, ‖bp_D s x₀ z t‖ ≤ R' * Md + Ms := by
    intro z hz t ht
    have hpK : x₀ + t • (z - x₀) ∈ closedBall x₀ R' := by
      rw [mem_closedBall_iff_norm, add_sub_cancel_left, norm_smul, Real.norm_of_nonneg ht.1]
      calc t * ‖z - x₀‖ ≤ 1 * ‖z - x₀‖ := by gcongr; exact ht.2
        _ ≤ R' := by rw [one_mul]; exact (hzR z hz).le
    have h1 : ‖t • ((innerSL ℝ (z - x₀)).comp (fderiv ℝ s (x₀ + t • (z - x₀))))‖ ≤ R' * Md := by
      rw [norm_smul, Real.norm_of_nonneg ht.1]
      calc t * ‖(innerSL ℝ (z - x₀)).comp (fderiv ℝ s (x₀ + t • (z - x₀)))‖
          ≤ 1 * (‖innerSL ℝ (z - x₀)‖ * ‖fderiv ℝ s (x₀ + t • (z - x₀))‖) := by
            gcongr
            · exact ht.2
            · exact ContinuousLinearMap.opNorm_comp_le _ _
        _ ≤ R' * Md := by
            rw [one_mul, innerSL_apply_norm]
            exact mul_le_mul (hzR z hz).le (hMd _ hpK) (norm_nonneg _) hR'0
    have h2 : ‖innerSL ℝ (s (x₀ + t • (z - x₀)))‖ ≤ Ms := by
      rw [innerSL_apply_norm]; exact hMs _ hpK
    exact (norm_add_le _ _).trans (add_le_add h1 h2)
  have key := intervalIntegral.hasFDerivAt_integral_of_dominated_of_fderiv_le (μ := volume)
    (a := 0) (b := 1) (s := ball x r)
    (F := fun z t => inner ℝ (s (x₀ + t • (z - x₀))) (z - x₀))
    (F' := fun z t => bp_D s x₀ z t) (bound := fun _ => R' * Md + Ms) (x₀ := x)
    (ball_mem_nhds x hr0)
    (by
      filter_upwards [ball_mem_nhds x hr0] with z hz
      rw [uIoc_of_le zero_le_one]
      exact ((bp_slice_cont hs (hzball z hz)).1.mono Ioc_subset_Icc_self).aestronglyMeasurable
        measurableSet_Ioc)
    (by
      have := (bp_slice_cont hs hx).1
      rw [← uIcc_of_le (zero_le_one' ℝ)] at this
      simpa using this.intervalIntegrable)
    (by
      rw [uIoc_of_le zero_le_one]
      exact ((bp_slice_cont hs hx).2.mono Ioc_subset_Icc_self).aestronglyMeasurable
        measurableSet_Ioc)
    (Filter.Eventually.of_forall fun t ht z hz => by
      rw [uIoc_of_le zero_le_one] at ht
      exact hbound z hz t ⟨ht.1.le, ht.2⟩)
    intervalIntegrable_const
    (Filter.Eventually.of_forall fun t ht z hz => by
      rw [uIoc_of_le zero_le_one] at ht
      exact bp_hasFDerivAt_integrand hsd (hzball z hz) ⟨ht.1.le, ht.2⟩)
  have hint : ∫ t in (0 : ℝ)..1, bp_D s x₀ x t = innerSL ℝ (s x) := by
    ext a
    rw [ContinuousLinearMap.intervalIntegral_apply, bp_ftc hs hsym hx a]
    · simp
    · have := (bp_slice_cont hs hx).2
      rw [← uIcc_of_le (zero_le_one' ℝ)] at this
      exact this.intervalIntegrable
  rw [hint] at key
  exact key

/-- **Poincaré lemma on a ball.** A smooth vector field on a ball whose derivative is symmetric is the
gradient of a smooth function vanishing at the centre. -/
theorem exists_potential_ball [FiniteDimensional ℝ F] {x₀ : F} {R : ℝ} {s : F → F}
    (hs : ContDiffOn ℝ ∞ s (ball x₀ R))
    (hsym : ∀ x ∈ ball x₀ R, ∀ a b : F, inner ℝ (fderiv ℝ s x a) b = inner ℝ a (fderiv ℝ s x b)) :
    ∃ φ : F → ℝ, ContDiffOn ℝ ∞ φ (ball x₀ R) ∧ φ x₀ = 0 ∧
      ∀ x ∈ ball x₀ R, ∀ a : F, fderiv ℝ φ x a = inner ℝ (s x) a := by
  refine ⟨fun y => ∫ t in (0 : ℝ)..1, inner ℝ (s (x₀ + t • (y - x₀))) (y - x₀), ?_, by simp, ?_⟩
  · rw [contDiffOn_infty_iff_fderiv_of_isOpen isOpen_ball]
    refine ⟨fun x hx => (bp_hasFDerivAt_pot hs hsym hx).differentiableAt.differentiableWithinAt, ?_⟩
    have : ContDiffOn ℝ ∞ (fun x => innerSL ℝ (s x)) (ball x₀ R) :=
      (innerSL ℝ : F →L⋆[ℝ] F →L[ℝ] ℝ).contDiff.comp_contDiffOn hs
    exact this.congr fun x hx => (bp_hasFDerivAt_pot hs hsym hx).fderiv
  · intro x hx a
    rw [(bp_hasFDerivAt_pot hs hsym hx).fderiv]
    simp

/-- Two functions on a ball with the same derivative and the same value at the centre agree. -/
theorem eqOn_ball_of_fderiv_eq {x₀ : F} {R : ℝ} {φ ψ : F → ℝ}
    (hφ : DifferentiableOn ℝ φ (ball x₀ R)) (hψ : DifferentiableOn ℝ ψ (ball x₀ R))
    (hd : ∀ x ∈ ball x₀ R, fderiv ℝ φ x = fderiv ℝ ψ x) (h0 : φ x₀ = ψ x₀) :
    Set.EqOn φ ψ (ball x₀ R) := by
  rcases le_or_gt R 0 with hR | hR
  · intro x hx
    exact absurd hx (by rw [ball_eq_empty.2 hR]; exact notMem_empty x)
  · refine (convex_ball x₀ R).eqOn_of_fderivWithin_eq hφ hψ isOpen_ball.uniqueDiffOn
      (fun x hx => ?_) (mem_ball_self hR) h0
    rw [fderivWithin_of_isOpen isOpen_ball hx, fderivWithin_of_isOpen isOpen_ball hx]
    exact hd x hx

/-- A bound on the derivative along the segment from the centre: if `u` is differentiable on the ball,
vanishes at the centre and `‖fderiv u x‖ ≤ a * (1 + ‖u x‖ ^ 2) ^ (3/2 : ℝ)` wherever `‖u x‖ ≤ c`, and
`R * a * (1 + c ^ 2) ^ (3/2 : ℝ) < c`, then `‖u x‖ < c` on the whole ball. -/
theorem norm_lt_of_fderiv_le_ball {W : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W]
    {x₀ : F} {R a c : ℝ} (ha : 0 ≤ a) (hc : 0 < c) {u : F → W}
    (hu : DifferentiableOn ℝ u (ball x₀ R)) (h0 : u x₀ = 0)
    (hbound : ∀ x ∈ ball x₀ R, ‖u x‖ ≤ c → ‖fderiv ℝ u x‖ ≤ a * (1 + ‖u x‖ ^ 2) ^ ((3 : ℝ) / 2))
    (hsmall : R * (a * (1 + c ^ 2) ^ ((3 : ℝ) / 2)) < c) :
    ∀ x ∈ ball x₀ R, ‖u x‖ < c := by
  intro x hx
  by_contra hcon
  push Not at hcon
  have he : ‖x - x₀‖ < R := mem_ball_iff_norm.1 hx
  set e := x - x₀ with he_def
  set f : ℝ → W := fun t => u (x₀ + t • e) with hf
  have hfd : ∀ t ∈ Icc (0 : ℝ) 1, HasDerivAt f (fderiv ℝ u (x₀ + t • e) e) t := by
    intro t ht
    have hp := bp_mem_ball hx ht
    have hsf : HasFDerivAt u (fderiv ℝ u (x₀ + t • e)) (x₀ + t • e) :=
      (hu.differentiableAt (isOpen_ball.mem_nhds hp)).hasFDerivAt
    have hpath : HasDerivAt (fun t : ℝ => x₀ + t • e) e t := by
      simpa using ((hasDerivAt_id t).smul_const e).const_add x₀
    exact hsf.comp_hasDerivAt t hpath
  have hfc : ContinuousOn f (Icc 0 1) := fun t ht => (hfd t ht).continuousAt.continuousWithinAt
  set S : Set ℝ := Icc (0 : ℝ) 1 ∩ (fun t => ‖f t‖) ⁻¹' Ici c with hS
  have hSc : IsClosed S := hfc.norm.preimage_isClosed_of_isClosed isClosed_Icc isClosed_Ici
  have h1S : (1 : ℝ) ∈ S := ⟨⟨zero_le_one, le_rfl⟩, by simpa [hf, he_def] using hcon⟩
  have hSb : BddBelow S := ⟨0, fun t ht => ht.1.1⟩
  have hts : sInf S ∈ S := hSc.csInf_mem ⟨1, h1S⟩ hSb
  set ts := sInf S with hts_def
  have hts01 : ts ∈ Icc (0 : ℝ) 1 := hts.1
  have hts_c : c ≤ ‖f ts‖ := hts.2
  have hlt : ∀ t ∈ Ico 0 ts, ‖f t‖ < c := by
    intro t ht
    by_contra h
    push Not at h
    have : ts ≤ t := csInf_le hSb ⟨⟨ht.1, ht.2.le.trans hts01.2⟩, h⟩
    exact absurd ht.2 (not_lt.2 this)
  have hf0 : f 0 = 0 := by simp [hf, h0]
  set K : ℝ := a * (1 + c ^ 2) ^ ((3 : ℝ) / 2) with hK_def
  have hK : 0 ≤ K := by positivity
  have hderb : ∀ t ∈ Ico 0 ts, ‖fderiv ℝ u (x₀ + t • e) e‖ ≤ K * ‖e‖ := by
    intro t ht
    have ht1 : t ∈ Icc (0 : ℝ) 1 := ⟨ht.1, ht.2.le.trans hts01.2⟩
    have hp := bp_mem_ball hx ht1
    have hfc' : ‖u (x₀ + t • e)‖ < c := hlt t ht
    have h1 := hbound _ hp hfc'.le
    have h2 : a * (1 + ‖u (x₀ + t • e)‖ ^ 2) ^ ((3 : ℝ) / 2) ≤ K := by
      rw [hK_def]
      gcongr
    exact (ContinuousLinearMap.le_opNorm _ _).trans
      (mul_le_mul_of_nonneg_right (h1.trans h2) (norm_nonneg _))
  have hmv := norm_image_sub_le_of_norm_deriv_le_segment' (f := f) (a := 0) (b := ts)
    (f' := fun t => fderiv ℝ u (x₀ + t • e) e) (C := K * ‖e‖)
    (fun t ht => (hfd t ⟨ht.1, ht.2.trans hts01.2⟩).hasDerivWithinAt) hderb ts
    ⟨hts01.1, le_rfl⟩
  rw [hf0, sub_zero, sub_zero] at hmv
  have hKe : 0 ≤ K * ‖e‖ := by positivity
  have : c ≤ R * K := by
    calc c ≤ ‖f ts‖ := hts_c
      _ ≤ K * ‖e‖ * ts := hmv
      _ ≤ K * ‖e‖ * 1 := by gcongr; exact hts01.2
      _ = K * ‖e‖ := mul_one _
      _ ≤ K * R := by gcongr
      _ = R * K := mul_comm _ _
  linarith

/-- info: 'CWFlow.Interp.exists_potential_ball' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms exists_potential_ball

/-- info: 'CWFlow.Interp.eqOn_ball_of_fderiv_eq' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms eqOn_ball_of_fderiv_eq

/-- info: 'CWFlow.Interp.norm_lt_of_fderiv_le_ball' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms norm_lt_of_fderiv_le_ball

end CWFlow.Interp
