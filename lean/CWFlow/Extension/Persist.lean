import CWFlow.Extension.VelBound
import CWFlow.Extension.GraphChart
import CWFlow.Extension.ChartGraphEvol
import CWFlow.Extension.Slope
import CWFlow.Extension.GraphBridge
import Mathlib.Analysis.ODE.Gronwall
import Mathlib.Analysis.Calculus.FDeriv.Symmetric

/-!
# Uniform persistence of graph charts along a flow

Let `Φ` be a flow of Lagrangian immersions of a compact manifold with bounded second fundamental form and
bounded energies (`IsFlowM`).  In a graph chart of `Φ t₀` the horizontal map `χ_t = horC ∘ Φ_t ∘ σ` stays
`C¹`-close to the identity, and the slope and the height of `Φ_t ∘ σ` stay small, for `t ∈ [t₀, t₀ + τ₀)` with
`τ₀` independent of `t₀` and of the chart (`exists_persist`).  The proof is a Gronwall argument using the
compiled velocity bound `exists_const_vel_bound`: the position moves by at most `C (t - t₀)`, and
`∂ₜ (D_y Y) v = D(velX) v` is bounded by `C |D_y Y v|`.
-/

open Metric Set
open CWFlow.Interp
open scoped ContDiff Manifold Topology RealInnerProductSpace

namespace CWFlow.Extension

section Abstract

variable {n : ℕ} {V : Type*} [NormedAddCommGroup V] [NormedSpace ℝ V]

local notation "𝔼" => EuclideanSpace ℝ (Fin n)

theorem ps_hasFDerivAt_slice {Y : ℝ × 𝔼 → V} {a T : ℝ}
    (hY : ContDiffOn ℝ ∞ Y (Ioo a T ×ˢ ball (0 : 𝔼) 6)) {s : ℝ} (hs : s ∈ Ioo a T) {y : 𝔼}
    (hy : y ∈ ball (0 : 𝔼) 6) :
    HasFDerivAt (fun y' => Y (s, y')) ((fderiv ℝ Y (s, y)).comp (ContinuousLinearMap.inr ℝ ℝ 𝔼)) y := by
  have hO : IsOpen (Ioo a T ×ˢ ball (0 : 𝔼) 6) := isOpen_Ioo.prod isOpen_ball
  have hat : ContDiffAt ℝ ∞ Y (s, y) := hY.contDiffAt (hO.mem_nhds ⟨hs, hy⟩)
  exact (hat.differentiableAt (by simp)).hasFDerivAt.comp y (hasFDerivAt_prodMk_right s y)

theorem ps_hasDerivAt_dY {Y : ℝ × 𝔼 → V} {W : ℝ → 𝔼 → V} {a T : ℝ}
    (hY : ContDiffOn ℝ ∞ Y (Ioo a T ×ˢ ball (0 : 𝔼) 6))
    (hW1 : ∀ t ∈ Ioo a T, ∀ y ∈ ball (0 : 𝔼) 6, HasDerivAt (fun s => Y (s, y)) (W t y) t)
    {t : ℝ} (ht : t ∈ Ioo a T) {y : 𝔼} (hy : y ∈ ball (0 : 𝔼) 6) (v : 𝔼) :
    HasDerivAt (fun s => fderiv ℝ (fun y' => Y (s, y')) y v) (fderiv ℝ (W t) y v) t := by
  have hO : IsOpen (Ioo a T ×ˢ ball (0 : 𝔼) 6) := isOpen_Ioo.prod isOpen_ball
  have hY2 : ∀ p ∈ Ioo a T ×ˢ ball (0 : 𝔼) 6, ContDiffAt ℝ ∞ Y p := fun p hp =>
    hY.contDiffAt (hO.mem_nhds hp)
  have hg : ∀ p ∈ Ioo a T ×ˢ ball (0 : 𝔼) 6, ContDiffAt ℝ ∞ (fderiv ℝ Y) p := fun p hp =>
    (hY2 p hp).fderiv_right (m := ∞) (by simp)
  have hgd : HasFDerivAt (fderiv ℝ Y) (fderiv ℝ (fderiv ℝ Y) (t, y)) (t, y) :=
    ((hg (t, y) ⟨ht, hy⟩).differentiableAt (by simp)).hasFDerivAt
  have hsym := (hY2 (t, y) ⟨ht, hy⟩).isSymmSndFDerivAt (by simp; exact WithTop.coe_le_coe.2 le_top)
  have hline : HasDerivAt (fun s : ℝ => ((s, y) : ℝ × 𝔼)) (((1 : ℝ), (0 : 𝔼)) : ℝ × 𝔼) t :=
    (hasDerivAt_id' t).prodMk (hasDerivAt_const t y)
  have h1 : HasDerivAt (fun s => fderiv ℝ Y (s, y))
      (fderiv ℝ (fderiv ℝ Y) (t, y) ((1 : ℝ), (0 : 𝔼))) t :=
    hgd.comp_hasDerivAt t hline
  have h2 : HasDerivAt (fun s => fderiv ℝ Y (s, y) ((0 : ℝ), v))
      (fderiv ℝ (fderiv ℝ Y) (t, y) ((1 : ℝ), (0 : 𝔼)) ((0 : ℝ), v)) t :=
    (ContinuousLinearMap.apply ℝ V ((0 : ℝ), v)).hasFDerivAt.comp_hasDerivAt t h1
  have h3 : HasDerivAt (fun s => fderiv ℝ (fun y' => Y (s, y')) y v)
      (fderiv ℝ (fderiv ℝ Y) (t, y) ((1 : ℝ), (0 : 𝔼)) ((0 : ℝ), v)) t := by
    refine h2.congr_of_eventuallyEq ?_
    filter_upwards [isOpen_Ioo.mem_nhds ht] with s hs
    rw [(ps_hasFDerivAt_slice hY hs hy).fderiv]
    simp
  have hev : W t =ᶠ[𝓝 y] fun y' => fderiv ℝ Y (t, y') ((1 : ℝ), (0 : 𝔼)) := by
    filter_upwards [isOpen_ball.mem_nhds hy] with y' hy'
    have hl' : HasDerivAt (fun s : ℝ => ((s, y') : ℝ × 𝔼)) (((1 : ℝ), (0 : 𝔼)) : ℝ × 𝔼) t :=
      (hasDerivAt_id' t).prodMk (hasDerivAt_const t y')
    have hd2 : HasDerivAt (fun s => Y (s, y')) (fderiv ℝ Y (t, y') ((1 : ℝ), (0 : 𝔼))) t :=
      ((hY2 (t, y') ⟨ht, hy'⟩).differentiableAt (by simp)).hasFDerivAt.comp_hasDerivAt t hl'
    exact (hW1 t ht y' hy').unique hd2
  have hk : HasFDerivAt (fun y' => fderiv ℝ Y (t, y'))
      ((fderiv ℝ (fderiv ℝ Y) (t, y)).comp (ContinuousLinearMap.inr ℝ ℝ 𝔼)) y :=
    hgd.comp y (hasFDerivAt_prodMk_right t y)
  have hk2 := (ContinuousLinearMap.apply ℝ V ((1 : ℝ), (0 : 𝔼))).hasFDerivAt.comp y hk
  have hfd : fderiv ℝ (W t) y v =
      fderiv ℝ (fderiv ℝ Y) (t, y) ((0 : ℝ), v) ((1 : ℝ), (0 : 𝔼)) := by
    rw [hev.fderiv_eq]
    rw [show (fun y' => fderiv ℝ Y (t, y') ((1 : ℝ), (0 : 𝔼))) =
      (ContinuousLinearMap.apply ℝ V ((1 : ℝ), (0 : 𝔼))) ∘ (fun y' => fderiv ℝ Y (t, y')) from rfl]
    rw [hk2.fderiv]
    simp
  rw [hfd, ← hsym ((1 : ℝ), (0 : 𝔼)) ((0 : ℝ), v)]
  exact h3

theorem ps_flow_estimates {Y : ℝ × 𝔼 → V} {W : ℝ → 𝔼 → V} {a T C : ℝ} (hC : 0 < C)
    (hY : ContDiffOn ℝ ∞ Y (Ioo a T ×ˢ ball (0 : 𝔼) 6))
    (hW1 : ∀ t ∈ Ioo a T, ∀ y ∈ ball (0 : 𝔼) 6, HasDerivAt (fun s => Y (s, y)) (W t y) t)
    (hW2 : ∀ t ∈ Ioo a T, ∀ y ∈ ball (0 : 𝔼) 6, ‖W t y‖ ≤ C)
    (hW3 : ∀ t ∈ Ioo a T, ∀ y ∈ ball (0 : 𝔼) 6, ∀ v,
      ‖fderiv ℝ (W t) y v‖ ≤ C * ‖fderiv ℝ (fun y' => Y (t, y')) y v‖)
    {t₀ : ℝ} (ht₀ : t₀ ∈ Ioo a T) {t : ℝ} (ht : t ∈ Ico t₀ T) {y : 𝔼} (hy : y ∈ ball (0 : 𝔼) 6) :
    ‖Y (t, y) - Y (t₀, y)‖ ≤ C * (t - t₀) ∧
      ∀ v, ‖fderiv ℝ (fun y' => Y (t, y')) y v - fderiv ℝ (fun y' => Y (t₀, y')) y v‖ ≤
        (Real.exp (C * (t - t₀)) - 1) * ‖fderiv ℝ (fun y' => Y (t₀, y')) y v‖ := by
  have hsub : Icc t₀ t ⊆ Ioo a T := fun s hs => ⟨ht₀.1.trans_le hs.1, lt_of_le_of_lt hs.2 ht.2⟩
  constructor
  · exact norm_image_sub_le_of_norm_deriv_le_segment' (f := fun s => Y (s, y))
      (f' := fun s => W s y)
      (fun s hs => (hW1 s (hsub hs) y hy).hasDerivWithinAt)
      (fun s hs => hW2 s (hsub (Ico_subset_Icc_self hs)) y hy) t ⟨ht.1, le_rfl⟩
  · intro v
    set D : ℝ → V := fun s => fderiv ℝ (fun y' => Y (s, y')) y v with hD
    have hd : ∀ s ∈ Icc t₀ t, HasDerivAt D (fderiv ℝ (W s) y v) s := fun s hs =>
      ps_hasDerivAt_dY hY hW1 (hsub hs) hy v
    have key := norm_le_gronwallBound_of_norm_deriv_right_le (f := fun s => D s - D t₀)
      (f' := fun s => fderiv ℝ (W s) y v) (δ := 0) (K := C) (ε := C * ‖D t₀‖) (a := t₀) (b := t)
      (fun s hs => ((hd s hs).sub_const _).continuousAt.continuousWithinAt)
      (fun s hs => ((hd s (Ico_subset_Icc_self hs)).sub_const _).hasDerivWithinAt)
      (by simp)
      (fun s hs => by
        have h1 := hW3 s (hsub (Ico_subset_Icc_self hs)) y hy v
        have h2 := norm_le_norm_add_norm_sub' (D s) (D t₀)
        have h3 := mul_le_mul_of_nonneg_left h2 hC.le
        show ‖fderiv ℝ (W s) y v‖ ≤ C * ‖D s - D t₀‖ + C * ‖D t₀‖
        change ‖fderiv ℝ (W s) y v‖ ≤ C * ‖D s‖ at h1
        nlinarith)
    have := key t ⟨ht.1, le_rfl⟩
    rw [gronwallBound_of_K_ne_0 hC.ne'] at this
    refine this.trans (le_of_eq ?_)
    field_simp
    ring

end Abstract

section Local

variable {n : ℕ} {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
  {J : V →L[ℝ] V} {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V}

local notation "𝔼" => EuclideanSpace ℝ (Fin n)

theorem ps_horCL_le (u : V) : ‖horCL U u‖ ≤ ‖u‖ := by
  have := cv_horC_lip U 0 u 0
  simpa [horC_eq] using this

theorem ps_verCL_le (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (u : V) : ‖verCL J U u‖ ≤ ‖u‖ := by
  have := cv_verC_lip hJi U 0 u 0
  simpa [verC_eq] using this

theorem ps_norm_J (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (u : V) : ‖J u‖ = ‖u‖ := by
  have h : ‖J u‖ ^ 2 = ‖u‖ ^ 2 := by
    rw [← real_inner_self_eq_norm_sq, ← real_inner_self_eq_norm_sq, hJi]
  exact (sq_eq_sq₀ (norm_nonneg _) (norm_nonneg _)).mp h

theorem ps_horver (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hU : IsLagFrame J U) (v w : 𝔼) :
    horCL U (U v + J (U w)) = v ∧ verCL J U (U v + J (U w)) = w := by
  constructor
  · ext a
    rw [horCL_apply]
    simp only [inner_add_right, gb_UU, gb_UJU hU, add_zero]
    simp [eb, EuclideanSpace.inner_single_left]
  · ext a
    rw [verCL_apply]
    simp only [inner_add_right, gb_JUU hU, gb_JUJU hJi, zero_add]
    simp [eb, EuclideanSpace.inner_single_left]

theorem ps_global_ext {φ : 𝔼 → ℝ} (hφ : ContDiffOn ℝ ∞ φ (ball (0 : 𝔼) 6)) {y : 𝔼}
    (hy : y ∈ ball (0 : 𝔼) 6) : ∃ φ' : 𝔼 → ℝ, ContDiff ℝ ∞ φ' ∧ φ' =ᶠ[𝓝 y] φ := by
  classical
  obtain ⟨r, hr, hsub⟩ : ∃ r > 0, closedBall y (2 * r) ⊆ ball (0 : 𝔼) 6 := by
    obtain ⟨ε, hε, hb⟩ := Metric.isOpen_iff.1 isOpen_ball y hy
    exact ⟨ε / 4, by positivity,
      (closedBall_subset_ball (show 2 * (ε / 4) < ε by linarith)).trans hb⟩
  let b : ContDiffBump y := ⟨r, 2 * r, hr, by linarith⟩
  refine ⟨fun x => if x ∈ ball (0 : 𝔼) 6 then b x * φ x else 0, ?_, ?_⟩
  · rw [contDiff_iff_contDiffAt]
    intro x
    by_cases hx : x ∈ ball (0 : 𝔼) 6
    · have h1 : ContDiffAt ℝ ∞ (fun x => b x * φ x) x :=
        b.contDiff.contDiffAt.mul (hφ.contDiffAt (isOpen_ball.mem_nhds hx))
      refine h1.congr_of_eventuallyEq ?_
      filter_upwards [isOpen_ball.mem_nhds hx] with x' hx'
      simp [hx']
    · have hnot : x ∉ closedBall y (2 * r) := fun h => hx (hsub h)
      have hev : (fun x => if x ∈ ball (0 : 𝔼) 6 then b x * φ x else 0) =ᶠ[𝓝 x]
          fun _ => (0 : ℝ) := by
        filter_upwards [isClosed_closedBall.isOpen_compl.mem_nhds hnot] with x' hx'
        have hb0 : b x' = 0 := image_eq_zero_of_notMem_tsupport (by
          rw [b.tsupport_eq]; exact hx')
        split_ifs <;> simp [hb0]
      exact (contDiffAt_const (c := (0 : ℝ))).congr_of_eventuallyEq hev
  · filter_upwards [closedBall_mem_nhds y hr] with x hx
    have hx6 : x ∈ ball (0 : 𝔼) 6 := hsub (closedBall_subset_closedBall (by linarith) hx)
    simp [hx6, b.one_of_mem_closedBall hx]

theorem ps_graph_data (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    (hU : IsLagFrame J U) {φ : 𝔼 → ℝ} (hφ : ContDiffOn ℝ ∞ φ (ball (0 : 𝔼) 6)) {y : 𝔼}
    (hy : y ∈ ball (0 : 𝔼) 6) (y₀ : V) :
    (∀ v, horCL U (fderiv ℝ (graphX J U y₀ φ) y v) = v) ∧
    (∀ v, ‖verCL J U (fderiv ℝ (graphX J U y₀ φ) y v)‖ ≤ ‖jet 2 φ y‖ * ‖v‖) ∧
    (∀ v, ‖fderiv ℝ (graphX J U y₀ φ) y v‖ ≤ (1 + ‖jet 2 φ y‖) * ‖v‖) ∧
    slopeZ J U (fderiv ℝ (graphX J U y₀ φ) y) = jet 2 φ y := by
  obtain ⟨φ', hφ', hev⟩ := ps_global_ext hφ hy
  have hg : fderiv ℝ (graphX J U y₀ φ) y = fderiv ℝ (graphX J U y₀ φ') y := by
    refine Filter.EventuallyEq.fderiv_eq ?_
    filter_upwards [hev.eventually_nhds] with x hx
    exact (cge_graphX_congr hx).symm
  have hj : jet 2 φ' y = jet 2 φ y := by
    ext I
    simp only [jet_apply]
    rw [← iteratedFDerivWithin_univ, ← iteratedFDerivWithin_univ]
    rw [Filter.EventuallyEq.iteratedFDerivWithin_eq (by rwa [nhdsWithin_univ]) hev.self_of_nhds]
  rw [hg, ← hj]
  have hv : ∀ v, fderiv ℝ (graphX J U y₀ φ') y v =
      U v + J (U (fderiv ℝ (gradE φ') y v)) := fun v => gb_fderiv_graphX_apply hφ' y v
  refine ⟨fun v => by rw [hv]; exact (ps_horver hJi hU _ _).1, fun v => ?_, fun v => ?_,
    slopeZ_graphX hJi hJJ hU hφ' y₀ y⟩
  · rw [hv, (ps_horver hJi hU _ _).2]
    exact (ga_gradE_at hφ'.contDiffAt).2 v
  · rw [hv]
    refine (norm_add_le _ _).trans ?_
    rw [ps_norm_J hJi, U.norm_map, U.norm_map]
    have := (ga_gradE_at (x := y) hφ'.contDiffAt).2 v
    linarith

theorem ps_hor_diff {L L₀ : 𝔼 →L[ℝ] V} {ε : ℝ} (hP₀ : ∀ v, horCL U (L₀ v) = v)
    (hd : ∀ v, ‖L v - L₀ v‖ ≤ ε * ‖v‖) (u : 𝔼) : ‖horCL U (L u) - u‖ ≤ ε * ‖u‖ := by
  have : horCL U (L u) - u = horCL U (L u - L₀ u) := by rw [map_sub, hP₀]
  rw [this]
  exact (ps_horCL_le _).trans (hd u)

theorem ps_slope_frob (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) {L L₀ : 𝔼 →L[ℝ] V} {c ε : ℝ}
    (hε0 : 0 ≤ ε) (hε1 : ε ≤ 1 / 8) (hc0 : 0 ≤ c) (hc1 : c ≤ 1)
    (hP₀ : ∀ v, horCL U (L₀ v) = v) (hQ₀ : ∀ v, ‖verCL J U (L₀ v)‖ ≤ c * ‖v‖)
    (hd : ∀ v, ‖L v - L₀ v‖ ≤ ε * ‖v‖) :
    ‖slopeZ J U L - slopeZ J U L₀‖ ≤ ((n : ℝ) + 1) * (3 * ε) := by
  have hlow : ∀ u, (7 / 8) * ‖u‖ ≤ ‖horCL U (L u)‖ := by
    intro u
    have h := ps_hor_diff hP₀ hd u
    have h2 := norm_le_norm_add_norm_sub' u (horCL U (L u))
    rw [norm_sub_rev] at h2
    nlinarith [norm_nonneg u]
  have hinj : Function.Injective ((horCL U).comp L) := by
    intro x x' h
    have h0 : horCL U (L (x - x')) = 0 := by
      simp only [ContinuousLinearMap.comp_apply] at h
      rw [map_sub, map_sub, h, sub_self]
    have := hlow (x - x')
    rw [h0, norm_zero] at this
    have : ‖x - x'‖ = 0 := by nlinarith [norm_nonneg (x - x')]
    exact sub_eq_zero.1 (norm_eq_zero.1 this)
  obtain ⟨E, hEx, h2, hinv⟩ := sl_inverse hinj
  have hid : (horCL U).comp L₀ = ContinuousLinearMap.id ℝ 𝔼 := ContinuousLinearMap.ext hP₀
  have e1 : ∀ a b : Fin n, slopeZ J U L ![a, b] = verCL J U (L (E.symm (eb a))) b := by
    intro a b
    rw [sl_slope_apply, hinv]
    simp
  have e0 : ∀ a b : Fin n, slopeZ J U L₀ ![a, b] = verCL J U (L₀ (eb a)) b := by
    intro a b
    rw [sl_slope_apply, hid, ContinuousLinearMap.inverse_id]
    simp
  have col : ∀ a : Fin n,
      ‖verCL J U (L (E.symm (eb a))) - verCL J U (L₀ (eb a))‖ ≤ 3 * ε := by
    intro a
    set w : 𝔼 := eb a with hw
    set u : 𝔼 := E.symm w with hu
    have hEu : horCL U (L u) = w := by rw [← hEx]; simp [hu]
    have hnw : ‖w‖ = 1 := by simp [hw, eb]
    have hu8 : ‖u‖ ≤ 8 / 7 := by
      have := hlow u
      rw [hEu, hnw] at this
      linarith
    have huw : ‖u - w‖ ≤ ε * ‖u‖ := by
      have := ps_hor_diff hP₀ hd u
      rw [hEu, norm_sub_rev] at this
      exact this
    have hdec : verCL J U (L u) - verCL J U (L₀ w) =
        verCL J U (L₀ (u - w)) + verCL J U (L u - L₀ u) := by
      simp only [map_sub]; abel
    rw [hdec]
    have t1 : ‖verCL J U (L₀ (u - w))‖ ≤ c * (ε * ‖u‖) :=
      (hQ₀ _).trans (mul_le_mul_of_nonneg_left huw hc0)
    have t2 : ‖verCL J U (L u - L₀ u)‖ ≤ ε * ‖u‖ := (ps_verCL_le hJi _).trans (hd u)
    have t3 : c * (ε * ‖u‖) ≤ ε * ‖u‖ :=
      mul_le_of_le_one_left (mul_nonneg hε0 (norm_nonneg _)) hc1
    have t4 : ε * ‖u‖ ≤ ε * (8 / 7) := mul_le_mul_of_nonneg_left hu8 hε0
    refine (norm_add_le _ _).trans ?_
    linarith
  have hsq : ‖slopeZ J U L - slopeZ J U L₀‖ ^ 2 ≤ n * (3 * ε) ^ 2 := by
    rw [sl_normsq2]
    calc ∑ a, ∑ b, ((slopeZ J U L - slopeZ J U L₀) ![a, b]) ^ 2
        = ∑ a, ‖verCL J U (L (E.symm (eb a))) - verCL J U (L₀ (eb a))‖ ^ 2 := by
          refine Finset.sum_congr rfl fun a _ => ?_
          rw [EuclideanSpace.norm_sq_eq]
          refine Finset.sum_congr rfl fun b _ => ?_
          simp only [PiLp.sub_apply, e1, e0, Real.norm_eq_abs, sq_abs]
      _ ≤ ∑ _a : Fin n, (3 * ε) ^ 2 :=
          Finset.sum_le_sum fun a _ => pow_le_pow_left₀ (norm_nonneg _) (col a) 2
      _ = n * (3 * ε) ^ 2 := by simp
  refine (sq_le_sq₀ (norm_nonneg _) (by positivity)).1 (hsq.trans ?_)
  have hn : (0 : ℝ) ≤ n := Nat.cast_nonneg n
  nlinarith [mul_nonneg (add_nonneg (add_nonneg (sq_nonneg (n : ℝ)) hn) zero_le_one)
    (sq_nonneg (3 * ε))]

end Local

section Main


theorem exists_persist {n : ℕ} {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    [FiniteDimensional ℝ E] {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
    [T2Space M] [CompactSpace M]
    {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V] [FiniteDimensional ℝ V]
    (ι : EuclideanSpace ℝ (Fin n) ≃L[ℝ] E) (A : IntAtlas (E := E) M) {J : V →L[ℝ] V}
    (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    (hdim : Module.finrank ℝ V = 2 * n)
    {c₀ c₁ : ℝ} (hc₀ : 0 < c₀) (hc₀1 : c₀ ≤ 1 / 100) (hc₁ : 0 ≤ c₁)
    (hsmall : 6 * (c₁ * (1 + c₀ ^ 2) ^ ((3 : ℝ) / 2)) < c₀)
    (hK₃ : (1 + c₀ ^ 2) ^ ((3 : ℝ) / 2) * c₁ ≤ 1) (B : ℕ → ℝ) :
    ∃ τ₀ > 0, ∀ {a T : ℝ} {Φ : ℝ → M → V}, IsFlowM ι J (Ioo a T) Φ →
      (∀ t ∈ Ioo a T, ∀ z : M, ∀ y ∈ chartDomE ι z, aNormSq J (Φ t ∘ chartInv ι z) y ≤ c₁ ^ 2) →
      (∀ t ∈ Ioo a T, ∀ l, energyM ι A J l (Φ t) ≤ B l) →
      ∀ {t₀ : ℝ}, t₀ ∈ Ioo a T → ∀ {p : M} {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V}
        {σ : EuclideanSpace ℝ (Fin n) → M} {φ₀ : EuclideanSpace ℝ (Fin n) → ℝ},
      IsGraphChart ι J (Φ t₀) c₀ p U σ φ₀ →
      ∀ t ∈ Ico t₀ (t₀ + τ₀), t < T → ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) (11 / 2),
        ‖fderiv ℝ (fun y' => horC U (Φ t₀ p) (Φ t (σ y'))) y - ContinuousLinearMap.id ℝ _‖ ≤ 1 / 8 ∧
        ‖horC U (Φ t₀ p) (Φ t (σ y)) - y‖ ≤ 1 / 8 ∧
        ‖verC J U (Φ t₀ p) (Φ t (σ y))‖ ≤ 1 / 2 ∧
        slopeLE J U (2 * c₀) (fderiv ℝ (fun y' => Φ t (σ y')) y) ∧
        ‖slopeZ J U (fderiv ℝ (fun y' => Φ t (σ y')) y)‖ ≤ 2 * c₀ := by
  obtain ⟨C₀, hC₀⟩ := exists_const_vel_bound ι A hJi hJJ hdim hc₀ hc₀1 hc₁ hsmall hK₃ B
  set C : ℝ := max C₀ 1 with hCdef
  have hC1 : 1 ≤ C := le_max_right _ _
  have hCpos : 0 < C := by linarith
  have hn0 : (0 : ℝ) ≤ n := Nat.cast_nonneg n
  set ε₀ : ℝ := min (1 / 8) (c₀ / (4 * ((n : ℝ) + 1))) with hε
  have hε0 : 0 < ε₀ := lt_min (by norm_num) (by positivity)
  have hε1 : ε₀ ≤ 1 / 8 := min_le_left _ _
  have hε2 : ε₀ * (4 * ((n : ℝ) + 1)) ≤ c₀ := by
    have h := min_le_right (1 / 8 : ℝ) (c₀ / (4 * ((n : ℝ) + 1)))
    have hp : 0 < 4 * ((n : ℝ) + 1) := by positivity
    calc ε₀ * (4 * ((n : ℝ) + 1)) ≤ c₀ / (4 * ((n : ℝ) + 1)) * (4 * ((n : ℝ) + 1)) :=
          mul_le_mul_of_nonneg_right h hp.le
      _ = c₀ := by field_simp
  obtain ⟨τ₀, hτ₀, hτ1, hτ2⟩ : ∃ τ₀ : ℝ, 0 < τ₀ ∧ C * τ₀ ≤ 1 / 8 ∧
      Real.exp (C * τ₀) - 1 ≤ ε₀ / 2 := by
    refine ⟨min (1 / (8 * C)) (Real.log (1 + ε₀ / 2) / C), ?_, ?_, ?_⟩
    · exact lt_min (by positivity) (div_pos (Real.log_pos (by linarith)) hCpos)
    · calc C * min (1 / (8 * C)) (Real.log (1 + ε₀ / 2) / C) ≤ C * (1 / (8 * C)) :=
            mul_le_mul_of_nonneg_left (min_le_left _ _) hCpos.le
        _ = 1 / 8 := by field_simp
    · have h1 : C * min (1 / (8 * C)) (Real.log (1 + ε₀ / 2) / C) ≤ Real.log (1 + ε₀ / 2) :=
        calc C * min (1 / (8 * C)) (Real.log (1 + ε₀ / 2) / C)
            ≤ C * (Real.log (1 + ε₀ / 2) / C) :=
              mul_le_mul_of_nonneg_left (min_le_right _ _) hCpos.le
          _ = Real.log (1 + ε₀ / 2) := by field_simp
      have h2 := Real.exp_le_exp.2 h1
      rw [Real.exp_log (by linarith)] at h2
      linarith
  refine ⟨τ₀, hτ₀, ?_⟩
  intro a T Φ hΦ hA hE t₀ ht₀ p U σ φ₀ hch t ht htT y hy
  have hy6 : y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6 := ball_subset_ball (by norm_num) hy
  have hyn : ‖y‖ < 11 / 2 := by simpa using hy
  have htm : t ∈ Ioo a T := ⟨ht₀.1.trans_le ht.1, htT⟩
  have hs0 : 0 ≤ t - t₀ := by linarith [ht.1]
  have hsτ : t - t₀ < τ₀ := by linarith [ht.2]
  have hCs : C * (t - t₀) ≤ 1 / 8 := by
    have := mul_le_mul_of_nonneg_left hsτ.le hCpos.le
    linarith
  have hexp : Real.exp (C * (t - t₀)) - 1 ≤ ε₀ / 2 := by
    have := Real.exp_le_exp.2 (mul_le_mul_of_nonneg_left hsτ.le hCpos.le)
    linarith
  -- the family
  have hY : ContDiffOn ℝ ∞ (fun q : ℝ × EuclideanSpace ℝ (Fin n) => Φ q.1 (σ q.2)) (Ioo a T ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) 6) :=
    cge_contDiffOn_family ι hΦ.smooth hch.param hΦ.isOpen
  have hW1 : ∀ s ∈ Ioo a T, ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6,
      HasDerivAt (fun r => Φ r (σ y)) (velX J (fun y' => Φ s (σ y')) y) s :=
    fun s hs y hy => cge_hasDerivAt_family ι hΦ hch.param hs hy
  have hbd := fun s (hs : s ∈ Ioo a T) y (hy : y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6) =>
    hC₀ (hΦ.imm s hs) (hA s hs) (hE s hs) hch.param y hy
  have hW2 : ∀ s ∈ Ioo a T, ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6, ‖velX J (fun y' => Φ s (σ y')) y‖ ≤ C :=
    fun s hs y hy => ((hbd s hs y hy).1).trans (le_max_left _ _)
  have hW3 : ∀ s ∈ Ioo a T, ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6, ∀ v,
      ‖fderiv ℝ (velX J (fun y' => Φ s (σ y'))) y v‖ ≤
        C * ‖fderiv ℝ (fun y' => (fun q : ℝ × EuclideanSpace ℝ (Fin n) => Φ q.1 (σ q.2)) (s, y')) y v‖ :=
    fun s hs y hy v => ((hbd s hs y hy).2 v).trans
      (mul_le_mul_of_nonneg_right (le_max_left _ _) (norm_nonneg _))
  obtain ⟨h1, h2⟩ := ps_flow_estimates hCpos hY hW1 hW2 hW3 ht₀ (⟨ht.1, htT⟩ : t ∈ Ico t₀ T) hy6
  have h1' : ‖Φ t (σ y) - Φ t₀ (σ y)‖ ≤ C * (t - t₀) := h1
  -- graph data at t₀
  obtain ⟨hgraph, hhor, hjet, hgrad⟩ := hch.graph y hy6
  have hev : (fun y' => Φ t₀ (σ y')) =ᶠ[𝓝 y] graphX J U (Φ t₀ p) φ₀ := by
    filter_upwards [isOpen_ball.mem_nhds hy6] with x hx
    exact (hch.graph x hx).1
  obtain ⟨hP₀', hQ₀', hN₀', hS₀⟩ := ps_graph_data hJi hJJ hch.frame hch.smooth hy6 (Φ t₀ p)
  set L : EuclideanSpace ℝ (Fin n) →L[ℝ] V := fderiv ℝ (fun y' => Φ t (σ y')) y with hL
  set L₀ : EuclideanSpace ℝ (Fin n) →L[ℝ] V := fderiv ℝ (fun y' => Φ t₀ (σ y')) y with hL₀
  have hL₀g : L₀ = fderiv ℝ (graphX J U (Φ t₀ p) φ₀) y := hev.fderiv_eq
  rw [← hL₀g] at hP₀' hQ₀' hN₀' hS₀
  have hjet' : ‖jet 2 φ₀ y‖ ≤ c₀ := hjet
  have hQ₀ : ∀ v, ‖verCL J U (L₀ v)‖ ≤ c₀ * ‖v‖ := fun v =>
    (hQ₀' v).trans (mul_le_mul_of_nonneg_right hjet' (norm_nonneg _))
  have hN₀ : ∀ v, ‖L₀ v‖ ≤ 2 * ‖v‖ := fun v =>
    (hN₀' v).trans (mul_le_mul_of_nonneg_right (by linarith) (norm_nonneg _))
  have hd : ∀ v, ‖L v - L₀ v‖ ≤ ε₀ * ‖v‖ := by
    intro v
    have h := h2 v
    have h3 : (Real.exp (C * (t - t₀)) - 1) * ‖L₀ v‖ ≤ ε₀ / 2 * (2 * ‖v‖) :=
      mul_le_mul hexp (hN₀ v) (norm_nonneg _) (by linarith)
    change ‖L v - L₀ v‖ ≤ (Real.exp (C * (t - t₀)) - 1) * ‖L₀ v‖ at h
    linarith
  have hLd : HasFDerivAt (fun y' => Φ t (σ y')) L y :=
    (ps_hasFDerivAt_slice hY htm hy6).differentiableAt.hasFDerivAt
  refine ⟨?_, ?_, ?_, ?_, ?_⟩
  · have hc := (hasFDerivAt_horC (U := U) (y₀ := Φ t₀ p) (Φ t (σ y))).comp y hLd
    have hfd : fderiv ℝ (fun y' => horC U (Φ t₀ p) (Φ t (σ y'))) y = (horCL U).comp L :=
      hc.fderiv
    rw [hfd]
    refine ContinuousLinearMap.opNorm_le_bound _ (by norm_num) fun v => ?_
    have := ps_hor_diff hP₀' hd v
    simp only [sub_apply, ContinuousLinearMap.comp_apply,
      ContinuousLinearMap.id_apply]
    nlinarith [norm_nonneg v]
  · calc ‖horC U (Φ t₀ p) (Φ t (σ y)) - y‖
        = ‖horC U (Φ t₀ p) (Φ t (σ y)) - horC U (Φ t₀ p) (Φ t₀ (σ y))‖ := by rw [hhor]
      _ ≤ ‖Φ t (σ y) - Φ t₀ (σ y)‖ := cv_horC_lip U _ _ _
      _ ≤ 1 / 8 := h1'.trans hCs
  · have hv0 : verC J U (Φ t₀ p) (Φ t₀ (σ y)) = gradE φ₀ y := by
      rw [hgraph]; exact verC_graphX hJi hJJ hch.frame φ₀ y
    have h3 := cv_verC_lip hJi U (Φ t₀ p) (Φ t (σ y)) (Φ t₀ (σ y))
    have h4 := norm_le_norm_add_norm_sub' (verC J U (Φ t₀ p) (Φ t (σ y)))
      (verC J U (Φ t₀ p) (Φ t₀ (σ y)))
    rw [hv0] at h4 h3
    have h5 : ‖gradE φ₀ y‖ ≤ c₀ * (11 / 2) := hgrad.trans (by nlinarith [norm_nonneg y])
    linarith
  · intro v
    have h3 : ‖verCL J U (L v)‖ ≤ ‖verCL J U (L₀ v)‖ + ε₀ * ‖v‖ := by
      have h4 := norm_le_norm_add_norm_sub' (verCL J U (L v)) (verCL J U (L₀ v))
      have h5 : ‖verCL J U (L v) - verCL J U (L₀ v)‖ ≤ ε₀ * ‖v‖ := by
        rw [← map_sub]; exact (ps_verCL_le hJi _).trans (hd v)
      linarith
    have h6 := ps_hor_diff hP₀' hd v
    have h7 := norm_le_norm_add_norm_sub' v (horCL U (L v))
    rw [norm_sub_rev] at h7
    have h8 := hQ₀ v
    have hnv := norm_nonneg v
    have hn' : 0 ≤ ‖horCL U (L v)‖ := norm_nonneg _
    nlinarith [mul_nonneg hc₀.le hnv, mul_nonneg hε0.le hnv, mul_nonneg hc₀.le hn',
      mul_nonneg hc₀.le (mul_nonneg hε0.le hnv)]
  · have hfr := ps_slope_frob hJi hε0.le hε1 hc₀.le (by linarith) hP₀' hQ₀ hd
    rw [hS₀] at hfr
    have h4 := norm_le_norm_add_norm_sub' (slopeZ J U L) (jet 2 φ₀ y)
    have h5 : ((n : ℝ) + 1) * (3 * ε₀) ≤ 3 / 4 * c₀ := by nlinarith
    linarith

end Main

/-- info: 'CWFlow.Extension.exists_persist' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_persist

end CWFlow.Extension
