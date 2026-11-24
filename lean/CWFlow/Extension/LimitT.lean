import CWFlow.Extension.TimeGraphFlow
import CWFlow.Extension.ChartGraphEvol
import DifferentialGeometry.Analysis.Calculus.SmoothExtension.BorelHalfLine.Parametric

/-!
# A chart extension at the final time

If `Y` is smooth on `[t₁, T] × B̄₂`, with the graph-type bound `‖L D Y_t − 1‖ ≤ 1/8`, the Lagrangian
condition and the flow equation on `(t₁, T)`, then the same three statements hold at `t = T`
(`limit_chart_at_T`).  The proof extends `Y` smoothly to all times near `y₀` (a parametric Borel
extension) and passes to the limit `t → T⁻` of the continuous quantities `D Y_t`, `velX J Y_t`.
-/

open Metric Set
open scoped ContDiff Topology RealInnerProductSpace

namespace CWFlow.Extension

/-- A closed bound passes to the endpoint. -/
theorem lt_le_at_end {t₁ T : ℝ} (hT : t₁ < T) {f : ℝ → ℝ} (hf : ContinuousAt f T) {c : ℝ}
    (h : ∀ t ∈ Ioo t₁ T, f t ≤ c) : f T ≤ c :=
  ContinuousWithinAt.closure_le (s := Ioo t₁ T) (by rw [closure_Ioo hT.ne]; exact right_mem_Icc.mpr hT.le)
    hf.continuousWithinAt continuousWithinAt_const h

theorem lt_inj_of_bound {n : ℕ} {V : Type*} [NormedAddCommGroup V] [NormedSpace ℝ V]
    (L : V →L[ℝ] EuclideanSpace ℝ (Fin n)) (D : EuclideanSpace ℝ (Fin n) →L[ℝ] V)
    (h : ‖L ∘L D - ContinuousLinearMap.id ℝ _‖ < 1) : Function.Injective D := by
  rw [injective_iff_map_eq_zero]
  intro x hx
  have h1 : (L ∘L D) x = 0 := by simp [hx]
  have h2 : ‖(L ∘L D - ContinuousLinearMap.id ℝ (EuclideanSpace ℝ (Fin n))) x‖ ≤ ‖L ∘L D - ContinuousLinearMap.id ℝ _‖ * ‖x‖ :=
    ContinuousLinearMap.le_opNorm _ _
  have h3 : (L ∘L D - ContinuousLinearMap.id ℝ (EuclideanSpace ℝ (Fin n))) x = -x := by simp [h1]
  rw [h3, norm_neg] at h2
  by_contra hx0
  have : 0 < ‖x‖ := norm_pos_iff.mpr hx0
  nlinarith

theorem limit_chart_at_T {n : ℕ} {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
    [FiniteDimensional ℝ V] {J : V →L[ℝ] V}
    {t₁ T : ℝ} (hT : t₁ < T) {Y : ℝ × EuclideanSpace ℝ (Fin n) → V}
    (hY : ContDiffOn ℝ ∞ Y (Icc t₁ T ×ˢ closedBall (0 : EuclideanSpace ℝ (Fin n)) 2))
    (L : V →L[ℝ] EuclideanSpace ℝ (Fin n))
    (hq : ∀ t ∈ Ioo t₁ T, ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 2,
      ‖L ∘L fderiv ℝ (fun y' => Y (t, y')) y - ContinuousLinearMap.id ℝ _‖ ≤ 1 / 8)
    (hlag : ∀ t ∈ Ioo t₁ T, ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 2, ∀ u w,
      ⟪J (fderiv ℝ (fun y' => Y (t, y')) y u), fderiv ℝ (fun y' => Y (t, y')) y w⟫ = 0)
    (hvel : ∀ t ∈ Ioo t₁ T, ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 2,
      HasDerivAt (fun s => Y (s, y)) (velX J (fun y' => Y (t, y')) y) t) :
    ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 2,
      ‖L ∘L fderiv ℝ (fun y' => Y (T, y')) y - ContinuousLinearMap.id ℝ _‖ ≤ 1 / 8 ∧
      (∀ u w, ⟪J (fderiv ℝ (fun y' => Y (T, y')) y u), fderiv ℝ (fun y' => Y (T, y')) y w⟫ = 0) ∧
      HasDerivWithinAt (fun s => Y (s, y)) (velX J (fun y' => Y (T, y')) y) (Icc t₁ T) T := by
  intro y₀ hy₀
  have hTt : T ∈ Icc t₁ T := right_mem_Icc.mpr hT.le
  -- the extension
  have hmem : y₀ ∈ interior (closedBall (0 : EuclideanSpace ℝ (Fin n)) 2) := by
    rw [interior_closedBall _ (by norm_num)]; exact hy₀
  have hg : ContDiffOn ℝ ∞ (Function.uncurry (fun (t : ℝ) (z : EuclideanSpace ℝ (Fin n)) => Y (t + t₁, z)))
      (Icc (0 : ℝ) (T - t₁) ×ˢ closedBall (0 : EuclideanSpace ℝ (Fin n)) 2) := by
    have hmap : ContDiff ℝ ∞ (fun p : ℝ × EuclideanSpace ℝ (Fin n) => ((p.1 + t₁ : ℝ), p.2)) :=
      (contDiff_fst.add contDiff_const).prodMk contDiff_snd
    refine hY.comp hmap.contDiffOn ?_
    rintro ⟨s, z⟩ ⟨hs, hz⟩
    refine ⟨⟨by linarith [hs.1], by linarith [hs.2]⟩, hz⟩
  obtain ⟨gext, V₀, hV₀, hgsm, hgeq⟩ := DifferentialGeometry.Analysis.borel_interval_extend_param
    (fun (t : ℝ) (z : EuclideanSpace ℝ (Fin n)) => Y (t + t₁, z)) (T - t₁) (by linarith)
    (closedBall (0 : EuclideanSpace ℝ (Fin n)) 2) y₀ hmem hg
  set Yt : ℝ × EuclideanSpace ℝ (Fin n) → V := fun q => gext (q.1 - t₁) q.2 with hYt
  set O : Set (ℝ × EuclideanSpace ℝ (Fin n)) := (univ : Set ℝ) ×ˢ interior V₀ with hO
  have hOopen : IsOpen O := isOpen_univ.prod isOpen_interior
  have hYtsm : ContDiffOn ℝ ∞ Yt O := by
    have hmap : ContDiff ℝ ∞ (fun p : ℝ × EuclideanSpace ℝ (Fin n) => ((p.1 - t₁ : ℝ), p.2)) :=
      (contDiff_fst.sub contDiff_const).prodMk contDiff_snd
    refine hgsm.comp hmap.contDiffOn ?_
    rintro ⟨s, z⟩ ⟨-, hz⟩
    exact ⟨mem_univ _, interior_subset hz⟩
  have hagree : ∀ t ∈ Icc t₁ T, ∀ z ∈ V₀, Yt (t, z) = Y (t, z) := by
    intro t ht z hz
    have := hgeq (t - t₁) ⟨by linarith [ht.1], by linarith [ht.2]⟩ z hz
    simp only [hYt]
    rw [this]
    simp
  have hy₀V : y₀ ∈ interior V₀ := mem_interior_iff_mem_nhds.mpr hV₀
  have hslice : ∀ t ∈ Icc t₁ T, (fun y => Yt (t, y)) =ᶠ[nhds y₀] fun y => Y (t, y) := by
    intro t ht
    filter_upwards [hV₀] with z hz using hagree t ht z hz
  have hDeq : ∀ t ∈ Icc t₁ T, tf_sd Yt (t, y₀) = fderiv ℝ (fun y' => Y (t, y')) y₀ := by
    intro t ht
    exact (hslice t ht).fderiv_eq
  have hveleq : ∀ t ∈ Icc t₁ T, velX J (fun y => Yt (t, y)) y₀ = velX J (fun y' => Y (t, y')) y₀ :=
    fun t ht => cge_velX_congr (hslice t ht)
  have hOnhds : ∀ t : ℝ, O ∈ nhds (t, y₀) := fun t => hOopen.mem_nhds ⟨mem_univ _, hy₀V⟩
  -- D is continuous
  have hDsm : ContDiffOn ℝ ∞ (tf_sd Yt) O := tf_sd_smooth hOopen hYtsm
  have hDcont : ContinuousAt (fun t : ℝ => tf_sd Yt (t, y₀)) T :=
    (hDsm.continuousOn.continuousAt (hOnhds T)).comp (f := fun t : ℝ => (t, y₀))
      (continuous_id.prodMk continuous_const).continuousAt
  have hLD : ContinuousAt (fun t : ℝ => L ∘L tf_sd Yt (t, y₀)) T :=
    ContinuousAt.clm_comp continuousAt_const hDcont
  -- Part 1: the bound at T
  have hbound : ‖L ∘L fderiv ℝ (fun y' => Y (T, y')) y₀ - ContinuousLinearMap.id ℝ _‖ ≤ 1 / 8 := by
    rw [← hDeq T hTt]
    refine lt_le_at_end hT (f := fun t => ‖L ∘L tf_sd Yt (t, y₀) - ContinuousLinearMap.id ℝ _‖) ?_ ?_
    · exact (hLD.sub continuousAt_const).norm
    · intro t ht
      rw [hDeq t (Ioo_subset_Icc_self ht)]
      exact hq t ht y₀ hy₀
  -- Part 2: Lagrangian condition at T
  have hlagT : ∀ u w, ⟪J (fderiv ℝ (fun y' => Y (T, y')) y₀ u), fderiv ℝ (fun y' => Y (T, y')) y₀ w⟫ = 0 := by
    intro u w
    rw [← hDeq T hTt]
    have hc : ContinuousAt (fun t : ℝ => ⟪J (tf_sd Yt (t, y₀) u), tf_sd Yt (t, y₀) w⟫) T := by
      have hu : ContinuousAt (fun t : ℝ => tf_sd Yt (t, y₀) u) T := hDcont.clm_apply continuousAt_const
      have hw : ContinuousAt (fun t : ℝ => tf_sd Yt (t, y₀) w) T := hDcont.clm_apply continuousAt_const
      exact (J.continuous.continuousAt.comp hu).inner hw
    refine le_antisymm (lt_le_at_end hT hc ?_) ?_
    · intro t ht
      rw [hDeq t (Ioo_subset_Icc_self ht), hlag t ht y₀ hy₀ u w]
    · have hc' : ContinuousAt (fun t : ℝ => -⟪J (tf_sd Yt (t, y₀) u), tf_sd Yt (t, y₀) w⟫) T := hc.neg
      have := lt_le_at_end hT hc' (c := 0) (by
        intro t ht
        rw [hDeq t (Ioo_subset_Icc_self ht), hlag t ht y₀ hy₀ u w]; simp)
      simpa using this
  -- Part 3: injectivity near (T, y₀) and continuity of velX
  set O' : Set (ℝ × EuclideanSpace ℝ (Fin n)) :=
    O ∩ (fun q => ‖L ∘L tf_sd Yt q - ContinuousLinearMap.id ℝ _‖) ⁻¹' Iio 1 with hO'
  have hcontN : ContinuousOn (fun q => ‖L ∘L tf_sd Yt q - ContinuousLinearMap.id ℝ _‖) O :=
    ((continuousOn_const.clm_comp hDsm.continuousOn).sub continuousOn_const).norm
  have hO'open : IsOpen O' := hcontN.isOpen_inter_preimage hOopen isOpen_Iio
  have hTO' : (T, y₀) ∈ O' := by
    refine ⟨⟨mem_univ _, hy₀V⟩, ?_⟩
    show ‖L ∘L tf_sd Yt (T, y₀) - ContinuousLinearMap.id ℝ _‖ < 1
    rw [hDeq T hTt]; linarith
  have hinj : ∀ q ∈ O', Function.Injective (fderiv ℝ (fun y => Yt (q.1, y)) q.2) :=
    fun q hq => lt_inj_of_bound L _ hq.2
  have hvsm := contDiffOn_velX_family (J := J) hO'open (hYtsm.mono inter_subset_left) hinj
  have hvcont : ContinuousAt (fun t : ℝ => velX J (fun y => Yt (t, y)) y₀) T :=
    (hvsm.continuousOn.continuousAt (hO'open.mem_nhds hTO')).comp (f := fun t : ℝ => (t, y₀))
      (continuous_id.prodMk continuous_const).continuousAt
  -- Part 4: the velocity
  set h : ℝ → V := fun s => Yt (s, y₀) with hh
  have hhsm : ContDiff ℝ ∞ h := by
    rw [contDiff_iff_contDiffAt]
    intro s
    have h1 : ContDiffAt ℝ ∞ Yt (s, y₀) := hYtsm.contDiffAt (hOnhds s)
    exact h1.comp s ((contDiff_id.prodMk contDiff_const).contDiffAt)
  have hhY : ∀ s ∈ Icc t₁ T, h s = Y (s, y₀) := fun s hs => hagree s hs y₀ (mem_of_mem_nhds hV₀)
  have hhd : ∀ t ∈ Ioo t₁ T, HasDerivAt h (velX J (fun y => Yt (t, y)) y₀) t := by
    intro t ht
    rw [hveleq t (Ioo_subset_Icc_self ht)]
    refine (hvel t ht y₀ hy₀).congr_of_eventuallyEq ?_
    filter_upwards [Ioo_mem_nhds ht.1 ht.2] with s hs using hhY s (Ioo_subset_Icc_self hs)
  have hne : (nhdsWithin T (Ioo t₁ T)).NeBot := by
    rw [← mem_closure_iff_nhdsWithin_neBot, closure_Ioo hT.ne]; exact hTt
  have hlim : deriv h T = velX J (fun y => Yt (T, y)) y₀ := by
    refine tendsto_nhds_unique_of_eventuallyEq (l := nhdsWithin T (Ioo t₁ T))
      ((hhsm.continuous_deriv (by simp)).continuousAt.tendsto.mono_left nhdsWithin_le_nhds)
      (hvcont.tendsto.mono_left nhdsWithin_le_nhds) ?_
    filter_upwards [self_mem_nhdsWithin] with t ht using (hhd t ht).deriv
  have hderivT : HasDerivAt h (velX J (fun y' => Y (T, y')) y₀) T := by
    rw [← hveleq T hTt, ← hlim]
    exact (hhsm.differentiable (by simp) T).hasDerivAt
  refine ⟨hbound, hlagT, ?_⟩
  exact hderivT.hasDerivWithinAt.congr (fun s hs => (hhY s hs).symm) (hhY T hTt).symm

/-- info: 'CWFlow.Extension.limit_chart_at_T' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms limit_chart_at_T

end CWFlow.Extension
