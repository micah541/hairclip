import CWFlow.Extension.LimitT
import CWFlow.Extension.GraphLocal
import CWFlow.Extension.FlowM
import CWFlow.Extension.ImmReparam
import CWFlow.Flow.LocalInverse

/-!
# Helpers for the extension of the flow to the final time

* `lm_local_rep`: if `σ` parametrizes a piece of `M` and `σ x₀ = chartInv ι z y`, the transition
  `τ` from the chart at `z` to the parameter is a smooth map on an open neighbourhood of `y` with
  invertible Jacobian and `chartInv ι z = σ ∘ τ` there (inverse function theorem for
  `chartFwd ι z ∘ σ`).
* `lm_local_T`: transfer of the facts at the final time from the graph parameter to the chart at `z`.
* `lm_chart_T`: the three hypotheses of `limit_chart_at_T` for the extended graph chart.
-/

open Metric Set Filter
open scoped ContDiff Manifold Topology RealInnerProductSpace

namespace CWFlow.Extension

theorem lm_local_rep {n : ℕ} {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    [FiniteDimensional ℝ E] {M : Type*} [TopologicalSpace M] [ChartedSpace E M]
    [IsManifold 𝓘(ℝ, E) ∞ M]
    (ι : EuclideanSpace ℝ (Fin n) ≃L[ℝ] E) {σ : EuclideanSpace ℝ (Fin n) → M}
    (hσ : IsParam ι (ball (0 : EuclideanSpace ℝ (Fin n)) 6) σ) {z : M}
    {y x₀ : EuclideanSpace ℝ (Fin n)} (hy : y ∈ chartDomE ι z)
    (hx₀ : x₀ ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 2) (hσx : σ x₀ = chartInv ι z y) :
    ∃ (τ : EuclideanSpace ℝ (Fin n) → EuclideanSpace ℝ (Fin n))
      (N : Set (EuclideanSpace ℝ (Fin n))), IsOpen N ∧ y ∈ N ∧ N ⊆ chartDomE ι z ∧
      ContDiffOn ℝ ∞ τ N ∧ (∀ w ∈ N, τ w ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 2) ∧
      (∀ w ∈ N, (jacM τ w).det ≠ 0) ∧ τ y = x₀ ∧ ∀ w ∈ N, chartInv ι z w = σ (τ w) := by
  classical
  have hQ := im_isOpen_W ι hσ z
  have hx₀Q : x₀ ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6 ∩ σ ⁻¹' (extChartAt 𝓘(ℝ, E) z).source :=
    ⟨ball_subset_ball (by norm_num) hx₀, by rw [mem_preimage, hσx]; exact im_chartInv_mem_source ι hy⟩
  have hinj : ∀ q ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6 ∩ σ ⁻¹' (extChartAt 𝓘(ℝ, E) z).source,
      Function.Injective (fderiv ℝ (chartFwd ι z ∘ σ) q) :=
    fun q hq => ir_inj_tau (hσ.trans_det z q hq)
  obtain ⟨e, hxe, heQ, hecoe, hesm⟩ := exists_openPartialHomeomorph_of_injective_fderiv
    (W := EuclideanSpace ℝ (Fin n)) (V := EuclideanSpace ℝ (Fin n)) rfl hQ (hσ.trans_smooth z) hinj hx₀Q
  have hNopen : IsOpen (e.target ∩ e.symm ⁻¹' ball (0 : EuclideanSpace ℝ (Fin n)) 2) :=
    e.continuousOn_symm.isOpen_inter_preimage e.open_target isOpen_ball
  have hye : e x₀ = y := by
    rw [hecoe x₀]
    show chartFwd ι z (σ x₀) = y
    rw [hσx]; exact im_chartFwd_chartInv ι hy
  have hyt : y ∈ e.target := by rw [← hye]; exact e.map_source hxe
  have hτy : e.symm y = x₀ := by rw [← hye]; exact e.left_inv hxe
  have hwQ : ∀ w ∈ e.target, e.symm w ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6 ∩
      σ ⁻¹' (extChartAt 𝓘(ℝ, E) z).source := fun w hw => heQ (e.map_target hw)
  have hg : ∀ w ∈ e.target, chartFwd ι z (σ (e.symm w)) = w := by
    intro w hw
    have := e.right_inv hw
    rw [hecoe] at this
    exact this
  refine ⟨e.symm, e.target ∩ e.symm ⁻¹' ball 0 2, hNopen, ⟨hyt, by rw [mem_preimage, hτy]; exact hx₀⟩,
    ?_, hesm.mono inter_subset_left, fun w hw => hw.2, ?_, hτy, ?_⟩
  · intro w hw
    rw [← hg w hw.1]
    exact im_chartFwd_mem ι (hwQ w hw.1).2
  · intro w hw
    refine im_det_ne_zero_of_left_inv (τ := e.symm) (τ' := chartFwd ι z ∘ σ) ?_ ?_ ?_
    · exact (hesm.contDiffAt (e.open_target.mem_nhds hw.1)).differentiableAt (by simp)
    · exact ((hσ.trans_smooth z).contDiffAt (hQ.mem_nhds (hwQ w hw.1))).differentiableAt (by simp)
    · filter_upwards [e.open_target.mem_nhds hw.1] with w' hw'
      exact hg w' hw'
  · intro w hw
    have h1 := im_chartInv_chartFwd ι (hwQ w hw.1).2
    rw [hg w hw.1] at h1
    exact h1

theorem lm_local_T {n : ℕ} {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    [FiniteDimensional ℝ E] {M : Type*} [TopologicalSpace M] [ChartedSpace E M]
    [IsManifold 𝓘(ℝ, E) ∞ M] {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
    (ι : EuclideanSpace ℝ (Fin n) ≃L[ℝ] E) {J : V →L[ℝ] V} {F : ℝ → M → V}
    {σ : EuclideanSpace ℝ (Fin n) → M}
    (hσ : IsParam ι (ball (0 : EuclideanSpace ℝ (Fin n)) 6) σ)
    {t₁ T : ℝ} (hT : t₁ < T) {Y : ℝ × EuclideanSpace ℝ (Fin n) → V}
    (hY : ContDiffOn ℝ ∞ Y (Icc t₁ T ×ˢ closedBall (0 : EuclideanSpace ℝ (Fin n)) 2))
    (hYF : ∀ s ∈ Icc t₁ T, ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 2, F s (σ y) = Y (s, y))
    (L : V →L[ℝ] EuclideanSpace ℝ (Fin n))
    (hbd : ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 2,
      ‖L ∘L fderiv ℝ (fun y' => Y (T, y')) y - ContinuousLinearMap.id ℝ _‖ ≤ 1 / 8)
    (hlg : ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 2, ∀ u w,
      ⟪J (fderiv ℝ (fun y' => Y (T, y')) y u), fderiv ℝ (fun y' => Y (T, y')) y w⟫ = 0)
    (hvl : ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 2,
      HasDerivWithinAt (fun s => Y (s, y)) (velX J (fun y' => Y (T, y')) y) (Icc t₁ T) T)
    {z : M} {y x₀ : EuclideanSpace ℝ (Fin n)} (hy : y ∈ chartDomE ι z)
    (hx₀ : x₀ ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 2) (hσx : σ x₀ = chartInv ι z y) :
    ContDiffWithinAt ℝ ∞ (fun q : ℝ × EuclideanSpace ℝ (Fin n) => F q.1 (chartInv ι z q.2))
      (Icc t₁ T ×ˢ chartDomE ι z) (T, y) ∧
    ContDiffAt ℝ ∞ (F T ∘ chartInv ι z) y ∧
    Function.Injective (fderiv ℝ (F T ∘ chartInv ι z) y) ∧
    (∀ u w, ⟪J (fderiv ℝ (F T ∘ chartInv ι z) y u), fderiv ℝ (F T ∘ chartInv ι z) y w⟫ = 0) ∧
    HasDerivWithinAt (fun s => F s (chartInv ι z y)) (velX J (F T ∘ chartInv ι z) y)
      (Icc t₁ T) T := by
  classical
  have hTt : T ∈ Icc t₁ T := right_mem_Icc.mpr hT.le
  obtain ⟨τ, N, hN, hyN, hNdom, hτ, hmaps, hdet, hτy, hinv⟩ := lm_local_rep ι hσ hy hx₀ hσx
  set X : EuclideanSpace ℝ (Fin n) → V := fun u => Y (T, u) with hX
  have hXsm : ContDiffOn ℝ ∞ X (ball (0 : EuclideanSpace ℝ (Fin n)) 2) := by
    refine hY.comp (contDiff_const.prodMk contDiff_id).contDiffOn ?_
    intro u hu
    exact ⟨hTt, ball_subset_closedBall hu⟩
  have hXinj : ∀ u ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 2, Function.Injective (fderiv ℝ X u) :=
    fun u hu => lt_inj_of_bound L _ (lt_of_le_of_lt (hbd u hu) (by norm_num))
  have hrep : ∀ s ∈ Icc t₁ T, ∀ w ∈ N, F s (chartInv ι z w) = Y (s, τ w) := by
    intro s hs w hw
    rw [hinv w hw]
    exact hYF s hs (τ w) (hmaps w hw)
  have hev : (F T ∘ chartInv ι z) =ᶠ[nhds y] X ∘ τ := by
    filter_upwards [hN.mem_nhds hyN] with w hw
    exact hrep T hTt w hw
  have hcd : ContDiffAt ℝ ∞ (F T ∘ chartInv ι z) y := by
    refine ContDiffAt.congr_of_eventuallyEq ?_ hev
    exact ContDiffAt.comp y (hXsm.contDiffAt (isOpen_ball.mem_nhds (hmaps y hyN)))
      (hτ.contDiffAt (hN.mem_nhds hyN))
  have hfd : fderiv ℝ (F T ∘ chartInv ι z) y = fderiv ℝ (X ∘ τ) y := hev.fderiv_eq
  refine ⟨?_, hcd, ?_, ?_, ?_⟩
  · have hf : ContDiffOn ℝ ∞ (fun q : ℝ × EuclideanSpace ℝ (Fin n) => (q.1, τ q.2))
        ((univ : Set ℝ) ×ˢ N) :=
      contDiffOn_fst.prodMk (hτ.comp contDiffOn_snd fun q hq => hq.2)
    have hg : ContDiffWithinAt ℝ ∞ Y (Icc t₁ T ×ˢ closedBall (0 : EuclideanSpace ℝ (Fin n)) 2)
        (T, τ y) := hY _ ⟨hTt, ball_subset_closedBall (hmaps y hyN)⟩
    have hcdw : ContDiffWithinAt ℝ ∞ (fun q : ℝ × EuclideanSpace ℝ (Fin n) => Y (q.1, τ q.2))
        (Icc t₁ T ×ˢ chartDomE ι z ∩ (univ : Set ℝ) ×ˢ N) (T, y) := by
      refine ContDiffWithinAt.comp (g := Y) (f := fun q : ℝ × EuclideanSpace ℝ (Fin n) => (q.1, τ q.2))
        (T, y) hg ((hf.contDiffWithinAt ⟨mem_univ _, hyN⟩).mono inter_subset_right) ?_
      rintro ⟨s, w⟩ ⟨⟨hs, -⟩, -, hw⟩
      exact ⟨hs, ball_subset_closedBall (hmaps w hw)⟩
    have h2 : ContDiffWithinAt ℝ ∞ (fun q : ℝ × EuclideanSpace ℝ (Fin n) => F q.1 (chartInv ι z q.2))
        (Icc t₁ T ×ˢ chartDomE ι z ∩ (univ : Set ℝ) ×ˢ N) (T, y) :=
      hcdw.congr (fun q hq => hrep q.1 hq.1.1 q.2 hq.2.2) (hrep T hTt y hyN)
    exact (contDiffWithinAt_inter ((isOpen_univ.prod hN).mem_nhds ⟨mem_univ _, hyN⟩)).mp h2
  · rw [hfd]
    exact injective_fderiv_comp isOpen_ball hXsm hXinj hN hτ hmaps hdet hyN
  · intro u w
    rw [hfd]
    exact lagrangian_comp isOpen_ball hXsm hlg hN hτ hmaps hyN u w
  · have hv : velX J (F T ∘ chartInv ι z) y = velX J (fun y' => Y (T, y')) x₀ := by
      rw [cge_velX_congr hev, velX_comp isOpen_ball hXsm hXinj hlg hN hτ hmaps hdet hyN, hτy]
    rw [hv]
    refine (hvl x₀ hx₀).congr (fun s hs => ?_) ?_
    · show F s (chartInv ι z y) = Y (s, x₀)
      rw [← hσx]
      exact hYF s hs x₀ hx₀
    · show F T (chartInv ι z y) = Y (T, x₀)
      rw [← hσx]
      exact hYF T hTt x₀ hx₀

theorem lm_chart_T {n : ℕ} {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    [FiniteDimensional ℝ E] {M : Type*} [TopologicalSpace M] [ChartedSpace E M]
    [IsManifold 𝓘(ℝ, E) ∞ M] {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
    [FiniteDimensional ℝ V]
    (ι : EuclideanSpace ℝ (Fin n) ≃L[ℝ] E) {J : V →L[ℝ] V} {a T t₀ t₁ : ℝ} {Φ : ℝ → M → V}
    (hΦ : IsFlowM ι J (Ioo a T) Φ) (ht₀ : a < t₀) (ht₁ : t₀ < t₁) (hT : t₁ < T)
    {σ : EuclideanSpace ℝ (Fin n) → M}
    (hσ : IsParam ι (ball (0 : EuclideanSpace ℝ (Fin n)) 6) σ)
    {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} {y₀ : V}
    (hq : ∀ t ∈ Ioo t₀ T, ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) (11 / 2),
      ‖fderiv ℝ (fun y' => horC U y₀ (Φ t (σ y'))) y - ContinuousLinearMap.id ℝ _‖ ≤ 1 / 8)
    {Y : ℝ × EuclideanSpace ℝ (Fin n) → V}
    (hY : ContDiffOn ℝ ∞ Y (Icc t₁ T ×ˢ closedBall (0 : EuclideanSpace ℝ (Fin n)) 2))
    (hYeq : EqOn Y (fun q => Φ q.1 (σ q.2))
      (Ioo t₁ T ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) 2)) :
    ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 2,
      ‖horCL U ∘L fderiv ℝ (fun y' => Y (T, y')) y - ContinuousLinearMap.id ℝ _‖ ≤ 1 / 8 ∧
      (∀ u w, ⟪J (fderiv ℝ (fun y' => Y (T, y')) y u), fderiv ℝ (fun y' => Y (T, y')) y w⟫ = 0) ∧
      HasDerivWithinAt (fun s => Y (s, y)) (velX J (fun y' => Y (T, y')) y) (Icc t₁ T) T := by
  have hIoo : ∀ t ∈ Ioo t₁ T, t ∈ Ioo a T := fun t ht => ⟨by linarith [ht.1], ht.2⟩
  have hslice : ∀ t ∈ Ioo t₁ T, ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 2,
      (fun y' => Y (t, y')) =ᶠ[nhds y] fun y' => Φ t (σ y') := by
    intro t ht y hy
    filter_upwards [isOpen_ball.mem_nhds hy] with y' hy'
    exact hYeq (show ((t, y') : ℝ × EuclideanSpace ℝ (Fin n)) ∈ _ from ⟨ht, hy'⟩)
  have hdiff : ∀ t ∈ Ioo a T, ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6,
      DifferentiableAt ℝ (fun y' => Φ t (σ y')) y := fun t ht y hy =>
    (((hΦ.imm t ht).comp_param ι hσ).1.contDiffAt (isOpen_ball.mem_nhds hy)).differentiableAt (by simp)
  have hfd : ∀ t ∈ Ioo t₁ T, ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 2,
      fderiv ℝ (fun y' => Y (t, y')) y = fderiv ℝ (fun y' => Φ t (σ y')) y :=
    fun t ht y hy => (hslice t ht y hy).fderiv_eq
  have hq' : ∀ t ∈ Ioo t₁ T, ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 2,
      ‖horCL U ∘L fderiv ℝ (fun y' => Y (t, y')) y - ContinuousLinearMap.id ℝ _‖ ≤ 1 / 8 := by
    intro t ht y hy
    have hy6 : y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) (11 / 2) :=
      ball_subset_ball (by norm_num) hy
    have h3 : HasFDerivAt (fun y' => horC U y₀ (Φ t (σ y')))
        ((horCL U).comp (fderiv ℝ (fun y' => Φ t (σ y')) y)) y :=
      HasFDerivAt.comp (g := horC U y₀) (f := fun y' => Φ t (σ y')) y
        (hasFDerivAt_horC (U := U) (y₀ := y₀) (Φ t (σ y)))
        (hdiff t (hIoo t ht) y (ball_subset_ball (by norm_num) hy)).hasFDerivAt
    have h2 : fderiv ℝ (fun y' => horC U y₀ (Φ t (σ y'))) y =
        horCL U ∘L fderiv ℝ (fun y' => Φ t (σ y')) y := h3.fderiv
    rw [hfd t ht y hy, ← h2]
    exact hq t ⟨lt_trans ht₁ ht.1, ht.2⟩ y hy6
  have hlag' : ∀ t ∈ Ioo t₁ T, ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 2, ∀ u w,
      ⟪J (fderiv ℝ (fun y' => Y (t, y')) y u), fderiv ℝ (fun y' => Y (t, y')) y w⟫ = 0 := by
    intro t ht y hy u w
    rw [hfd t ht y hy]
    exact ((hΦ.imm t (hIoo t ht)).comp_param ι hσ).2.2 y (ball_subset_ball (by norm_num) hy) u w
  have hvel' : ∀ t ∈ Ioo t₁ T, ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 2,
      HasDerivAt (fun s => Y (s, y)) (velX J (fun y' => Y (t, y')) y) t := by
    intro t ht y hy
    have h1 := cge_hasDerivAt_family ι hΦ hσ (hIoo t ht) (y := y) (ball_subset_ball (by norm_num) hy)
    rw [cge_velX_congr (hslice t ht y hy)]
    refine h1.congr_of_eventuallyEq ?_
    filter_upwards [Ioo_mem_nhds ht.1 ht.2] with s hs
    exact hYeq (show ((s, y) : ℝ × EuclideanSpace ℝ (Fin n)) ∈ _ from ⟨hs, hy⟩)
  exact limit_chart_at_T hT hY (horCL U) hq' hlag' hvel'

end CWFlow.Extension
