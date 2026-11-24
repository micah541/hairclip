import Mathlib.Analysis.Calculus.MeanValue
import Mathlib.Analysis.Calculus.ParametricIntervalIntegral
import Mathlib.Analysis.Calculus.ContDiff.Operations
import Mathlib.Analysis.Normed.Module.FiniteDimension
import Mathlib.MeasureTheory.Integral.IntervalIntegral.FundThmCalculus

/-!
# The time-integral potential

For `θ(s, y)` on `[0, T] × B` put `timePot θ (t, y) = −∫₀ᵗ θ(s, y) ds`. If a family of covectors `ξ(t, y)`
satisfies `ξ(0, ·) = 0` and `∂ₜξ = −D_yθ`, then `D_y (timePot θ) = ξ` (`hasFDerivAt_timePot`), and
`∂ₜ (timePot θ) = −θ` (`hasDerivWithinAt_timePot`). A function on `[0, T] × B` whose two partial derivatives
are smooth is smooth (`contDiffOn_of_partial`, from `hasFDerivWithinAt_of_partial`).
-/

open Set Filter MeasureTheory

open scoped Topology ContDiff Interval

set_option linter.unusedSectionVars false

noncomputable section

namespace CWFlow

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
variable {W : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W]

/-- A function on `I × B` (`B` open) with a spatial derivative `A` at the points of `I × B`, continuous at
`(t, y)`, and a time derivative at `(t, y)`, is differentiable at `(t, y)` within `I × B`. -/
theorem hasFDerivWithinAt_of_partial {I : Set ℝ} {B : Set E} (hB : IsOpen B) {f : ℝ × E → W}
    {A : ℝ × E → (E →L[ℝ] W)}
    (hA : ∀ p ∈ I ×ˢ B, HasFDerivAt (fun z => f (p.1, z)) (A p) p.2)
    {t : ℝ} {y : E} (ht : t ∈ I) (hy : y ∈ B)
    (hAc : ContinuousWithinAt A (I ×ˢ B) (t, y)) {b : W}
    (hb : HasDerivWithinAt (fun s => f (s, y)) b I t) :
    HasFDerivWithinAt f ((ContinuousLinearMap.fst ℝ ℝ E).smulRight b +
      (A (t, y)).comp (ContinuousLinearMap.snd ℝ ℝ E)) (I ×ˢ B) (t, y) := by
  have h2 : HasFDerivWithinAt (fun q : ℝ × E => f (q.1, y))
      ((ContinuousLinearMap.fst ℝ ℝ E).smulRight b) (I ×ˢ B) (t, y) := by
    have h := HasFDerivWithinAt.comp (g := fun s => f (s, y)) (f := fun q : ℝ × E => q.1) (t, y)
      hb.hasFDerivWithinAt hasFDerivWithinAt_fst (fun q (hq : q ∈ I ×ˢ B) => hq.1)
    refine h.congr_fderiv ?_
    ext <;> simp
  have h1 : HasFDerivWithinAt (fun q : ℝ × E => f q - f (q.1, y))
      ((A (t, y)).comp (ContinuousLinearMap.snd ℝ ℝ E)) (I ×ˢ B) (t, y) := by
    rw [hasFDerivWithinAt_iff_isLittleO, Asymptotics.isLittleO_iff]
    intro c hc
    obtain ⟨r, hr, hrA⟩ := Metric.continuousWithinAt_iff.mp hAc c hc
    obtain ⟨r', hr', hrB⟩ := Metric.isOpen_iff.mp hB y hy
    have hρ : 0 < min r r' := lt_min hr hr'
    filter_upwards [self_mem_nhdsWithin, nhdsWithin_le_nhds (Metric.ball_mem_nhds (t, y) hρ)]
      with q hqS hqb
    obtain ⟨s, z⟩ := q
    rw [Metric.mem_ball, Prod.dist_eq, max_lt_iff] at hqb
    have hsI : s ∈ I := hqS.1
    have hzb : z ∈ Metric.ball y (min r r') := hqb.2
    have hder : ∀ w ∈ Metric.ball y (min r r'),
        HasFDerivWithinAt (fun w => f (s, w) - A (t, y) w) (A (s, w) - A (t, y))
          (Metric.ball y (min r r')) w := fun w hw =>
      ((hA (s, w) ⟨hsI, hrB (Metric.ball_subset_ball (min_le_right _ _) hw)⟩).sub
        (A (t, y)).hasFDerivAt).hasFDerivWithinAt
    have hbd : ∀ w ∈ Metric.ball y (min r r'), ‖A (s, w) - A (t, y)‖ ≤ c := fun w hw => by
      have hw' : dist w y < min r r' := hw
      have hd : dist ((s, w) : ℝ × E) (t, y) < r := by
        rw [Prod.dist_eq, max_lt_iff]
        exact ⟨hqb.1.trans_le (min_le_left _ _), hw'.trans_le (min_le_left _ _)⟩
      have := hrA (x := (s, w)) ⟨hsI, hrB (Metric.ball_subset_ball (min_le_right _ _) hw)⟩ hd
      rw [dist_eq_norm] at this
      exact this.le
    have hmv := (convex_ball y (min r r')).norm_image_sub_le_of_norm_hasFDerivWithin_le hder hbd
      (Metric.mem_ball_self hρ) hzb
    have he : f (s, z) - f (s, y) - (f (t, y) - f (t, y)) - (A (t, y)) (z - y) =
        f (s, z) - A (t, y) z - (f (s, y) - A (t, y) y) := by
      rw [map_sub]
      abel
    show ‖f (s, z) - f (s, y) - (f (t, y) - f (t, y)) - (A (t, y)) (z - y)‖ ≤
      c * ‖((s, z) : ℝ × E) - (t, y)‖
    rw [he]
    refine hmv.trans (mul_le_mul_of_nonneg_left ?_ hc.le)
    exact norm_snd_le (((s, z) : ℝ × E) - (t, y))
  have h : HasFDerivWithinAt (fun q : ℝ × E => f q - f (q.1, y) + f (q.1, y))
      ((A (t, y)).comp (ContinuousLinearMap.snd ℝ ℝ E) +
        (ContinuousLinearMap.fst ℝ ℝ E).smulRight b) (I ×ˢ B) (t, y) := h1.add h2
  have hfun : (fun q : ℝ × E => f q - f (q.1, y) + f (q.1, y)) = f := funext fun q => by simp
  rw [hfun] at h
  refine h.congr_fderiv ?_
  rw [add_comm]

/-- A function on `[0, T] × B` with smooth partial derivatives `ξ` (space) and `−θ` (time) is smooth. -/
theorem contDiffOn_of_partial {T : ℝ} (hT : 0 < T) {B : Set E} (hB : IsOpen B) {f θ : ℝ × E → ℝ}
    {ξ : ℝ × E → (E →L[ℝ] ℝ)} (hθ : ContDiffOn ℝ ∞ θ (Icc 0 T ×ˢ B))
    (hξ : ContDiffOn ℝ ∞ ξ (Icc 0 T ×ˢ B))
    (h1 : ∀ p ∈ Icc (0 : ℝ) T ×ˢ B, HasFDerivAt (fun z => f (p.1, z)) (ξ p) p.2)
    (h2 : ∀ p ∈ Icc (0 : ℝ) T ×ˢ B, HasDerivWithinAt (fun s => f (s, p.2)) (-θ p) (Icc 0 T) p.1) :
    ContDiffOn ℝ ∞ f (Icc 0 T ×ˢ B) := by
  have hfd : ∀ p ∈ Icc (0 : ℝ) T ×ˢ B, HasFDerivWithinAt f
      ((-θ p) • ContinuousLinearMap.fst ℝ ℝ E + (ξ p).comp (ContinuousLinearMap.snd ℝ ℝ E))
      (Icc 0 T ×ˢ B) p := by
    rintro ⟨t, y⟩ hp
    have h := hasFDerivWithinAt_of_partial hB h1 hp.1 hp.2 (hξ.continuousOn (t, y) hp) (h2 (t, y) hp)
    refine h.congr_fderiv ?_
    ext <;> simp [mul_comm]
  have hUD : UniqueDiffOn ℝ (Icc (0 : ℝ) T ×ˢ B) := (uniqueDiffOn_Icc hT).prod hB.uniqueDiffOn
  rw [contDiffOn_infty_iff_fderivWithin hUD]
  refine ⟨fun p hp => (hfd p hp).differentiableWithinAt, ?_⟩
  have hsm : ContDiffOn ℝ ∞ (fun p : ℝ × E =>
      (-θ p) • ContinuousLinearMap.fst ℝ ℝ E + (ξ p).comp (ContinuousLinearMap.snd ℝ ℝ E))
      (Icc 0 T ×ˢ B) := (hθ.neg.smul contDiffOn_const).add (hξ.clm_comp contDiffOn_const)
  exact hsm.congr fun p hp => (hfd p hp).fderivWithin (hUD p hp)

/-- The time-integral potential `−∫₀ᵗ θ(s, y) ds`. -/
def timePot (θ : ℝ × E → ℝ) (p : ℝ × E) : ℝ := -∫ s in (0 : ℝ)..p.1, θ (s, p.2)

theorem timePot_zero (θ : ℝ × E → ℝ) (y : E) : timePot θ (0, y) = 0 := by
  simp [timePot]

variable {T : ℝ} {B : Set E} {θ : ℝ × E → ℝ}

/-- `∂ₜ (timePot θ) = −θ`. -/
theorem hasDerivWithinAt_timePot (hθ : ContinuousOn θ (Icc 0 T ×ˢ B)) {t : ℝ}
    (ht : t ∈ Icc (0 : ℝ) T) {y : E} (hy : y ∈ B) :
    HasDerivWithinAt (fun s => timePot θ (s, y)) (-θ (t, y)) (Icc 0 T) t := by
  have hcy : ContinuousOn (fun s => θ (s, y)) (Icc 0 T) :=
    hθ.comp (continuous_id.prodMk continuous_const).continuousOn (fun s hs => ⟨hs, hy⟩)
  have : Fact (t ∈ Icc (0 : ℝ) T) := ⟨ht⟩
  have hint : IntervalIntegrable (fun s => θ (s, y)) volume 0 t :=
    (hcy.mono (Icc_subset_Icc_right ht.2)).intervalIntegrable_of_Icc ht.1
  have h := intervalIntegral.integral_hasDerivWithinAt_right hint (s := Icc 0 T) (t := Icc 0 T)
    (hcy.stronglyMeasurableAtFilter_nhdsWithin measurableSet_Icc t) (hcy t ht)
  exact h.neg

variable [FiniteDimensional ℝ E]

/-- `D_y (timePot θ) = ξ` when `ξ(0, ·) = 0` and `∂ₜξ = −D_yθ`. -/
theorem hasFDerivAt_timePot (hB : IsOpen B) (hθ : ContinuousOn θ (Icc 0 T ×ˢ B))
    (hθd : ∀ p ∈ Icc (0 : ℝ) T ×ˢ B, DifferentiableAt ℝ (fun z => θ (p.1, z)) p.2)
    (hD : ContinuousOn (fun p : ℝ × E => fderiv ℝ (fun z => θ (p.1, z)) p.2) (Icc 0 T ×ˢ B))
    {ξ : ℝ × E → (E →L[ℝ] ℝ)} (hξ0 : ∀ y ∈ B, ξ (0, y) = 0)
    (hξd : ∀ t ∈ Icc (0 : ℝ) T, ∀ y ∈ B,
      HasDerivWithinAt (fun s => ξ (s, y)) (-fderiv ℝ (fun z => θ (t, z)) y) (Icc 0 T) t)
    {t : ℝ} (ht : t ∈ Icc (0 : ℝ) T) {y : E} (hy : y ∈ B) :
    HasFDerivAt (fun z => timePot θ (t, z)) (ξ (t, y)) y := by
  obtain ⟨ε, hε, hεB⟩ := Metric.isOpen_iff.mp hB y hy
  have hcb : Metric.closedBall y (ε / 2) ⊆ B :=
    (Metric.closedBall_subset_ball (half_lt_self hε)).trans hεB
  have hKS : Icc (0 : ℝ) t ×ˢ Metric.closedBall y (ε / 2) ⊆ Icc (0 : ℝ) T ×ˢ B :=
    prod_mono (Icc_subset_Icc_right ht.2) hcb
  have hK : IsCompact (Icc (0 : ℝ) t ×ˢ Metric.closedBall y (ε / 2)) :=
    isCompact_Icc.prod (isCompact_closedBall y (ε / 2))
  obtain ⟨C, hC⟩ := hK.exists_bound_of_continuousOn (hD.mono hKS)
  have hslice : ∀ z ∈ B, ContinuousOn (fun s => θ (s, z)) (Icc 0 t) := fun z hz =>
    hθ.comp (continuous_id.prodMk continuous_const).continuousOn
      (fun s hs => ⟨⟨hs.1, hs.2.trans ht.2⟩, hz⟩)
  have hDslice : ContinuousOn (fun s => fderiv ℝ (fun z => θ (s, z)) y) (Icc 0 t) :=
    ContinuousOn.comp (g := fun p : ℝ × E => fderiv ℝ (fun z => θ (p.1, z)) p.2)
      (f := fun s : ℝ => ((s, y) : ℝ × E)) hD
      (continuous_id.prodMk continuous_const).continuousOn
      (fun s hs => ⟨⟨hs.1, hs.2.trans ht.2⟩, hy⟩)
  have hIoc : Ι (0 : ℝ) t = Ioc 0 t := uIoc_of_le ht.1
  have key := intervalIntegral.hasFDerivAt_integral_of_dominated_of_fderiv_le (μ := volume)
    (a := 0) (b := t) (s := Metric.ball y (ε / 2)) (F := fun z s => θ (s, z))
    (F' := fun z s => fderiv ℝ (fun w => θ (s, w)) z) (bound := fun _ => C) (x₀ := y)
    (Metric.ball_mem_nhds y (half_pos hε))
    (by
      filter_upwards [Metric.ball_mem_nhds y (half_pos hε)] with z hz
      rw [hIoc]
      exact ((hslice z (hcb (Metric.ball_subset_closedBall hz))).mono
        Ioc_subset_Icc_self).aestronglyMeasurable measurableSet_Ioc)
    ((hslice y hy).intervalIntegrable_of_Icc ht.1)
    (by
      rw [hIoc]
      exact (hDslice.mono Ioc_subset_Icc_self).aestronglyMeasurable measurableSet_Ioc)
    (Eventually.of_forall fun s hs z hz => by
      rw [hIoc] at hs
      exact hC (s, z) ⟨⟨hs.1.le, hs.2⟩, Metric.ball_subset_closedBall hz⟩)
    intervalIntegrable_const
    (Eventually.of_forall fun s hs z hz => by
      rw [hIoc] at hs
      exact (hθd (s, z) ⟨⟨hs.1.le, hs.2.trans ht.2⟩,
        hcb (Metric.ball_subset_closedBall hz)⟩).hasFDerivAt)
  have hftc : ∫ s in (0 : ℝ)..t, fderiv ℝ (fun z => θ (s, z)) y = -ξ (t, y) := by
    have hcont : ContinuousOn (fun s => -ξ (s, y)) (Icc 0 t) := fun s hs =>
      ((hξd s ⟨hs.1, hs.2.trans ht.2⟩ y hy).continuousWithinAt.mono
        (Icc_subset_Icc_right ht.2)).neg
    have hder : ∀ s ∈ Ioo (0 : ℝ) t,
        HasDerivAt (fun s => -ξ (s, y)) (fderiv ℝ (fun z => θ (s, z)) y) s := fun s hs => by
      have h := ((hξd s ⟨hs.1.le, hs.2.le.trans ht.2⟩ y hy).hasDerivAt
        (Icc_mem_nhds hs.1 (hs.2.trans_le ht.2))).neg
      rwa [neg_neg] at h
    have h := intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le ht.1 hcont hder
      (hDslice.intervalIntegrable_of_Icc ht.1)
    rw [h, hξ0 y hy]
    simp
  rw [hftc] at key
  have h := key.neg
  rwa [neg_neg] at h

end CWFlow
