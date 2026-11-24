import CWFlow.Extension.GraphRadiusLocal
import CWFlow.Interp.RayLift

/-!
# Radial lifts of a Lagrangian immersion with small second fundamental form

Let `F : M → V` be a Lagrangian immersion of a compact manifold with `|A| ≤ c₁` (in the charts), `p ∈ M` with
`F p = y₀` and tangent plane at `p` equal to the horizontal plane (`verCL ∘ dFq p = 0`).  Writing
`f = horC U y₀ ∘ F` and `G = goodSet … 6`:

* `grf_apriori`: along a partial lift of a ray `s ↦ s • v` (`‖v‖ < 6`, `s < b ≤ 1`) the slope stays below `c₀` and
  `‖verC‖ ≤ 6 c₀ s`;
* `exists_radial_section`: the lifts of all rays assemble into a continuous section `σ` of `f` over the ball of
  radius `6`, with the slope bounds, whose image is the connected component of `p` in `G`.
-/

open Metric Set Filter
open CWFlow.Interp
open scoped ContDiff Manifold Topology RealInnerProductSpace

set_option linter.unusedSectionVars false

namespace CWFlow.Extension

variable {n : ℕ} {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
  {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
  {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
  (ι : EuclideanSpace ℝ (Fin n) ≃L[ℝ] E)
  {J : V →L[ℝ] V} {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} {y₀ : V} {F : M → V}

/-! ### Part 1: derivative bounds for the jets -/

theorem grf_clm_norm_sq_le {W : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W]
    (L : EuclideanSpace ℝ (Fin n) →L[ℝ] W) (v : EuclideanSpace ℝ (Fin n)) :
    ‖L v‖ ^ 2 ≤ ‖v‖ ^ 2 * ∑ i, ‖L (eb i)‖ ^ 2 := by
  have h1 : ‖L v‖ ≤ ∑ i, |v i| * ‖L (eb i)‖ := by
    rw [ir_lin_expand L v]
    refine (norm_sum_le _ _).trans (le_of_eq ?_)
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [norm_smul, Real.norm_eq_abs]
  have h2 : ‖v‖ ^ 2 = ∑ i, |v i| ^ 2 := by
    rw [EuclideanSpace.norm_sq_eq]
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [Real.norm_eq_abs]
  calc ‖L v‖ ^ 2 ≤ (∑ i, |v i| * ‖L (eb i)‖) ^ 2 := pow_le_pow_left₀ (norm_nonneg _) h1 2
    _ ≤ (∑ i, |v i| ^ 2) * ∑ i, ‖L (eb i)‖ ^ 2 := Finset.sum_mul_sq_le_sq_mul_sq _ _ _
    _ = _ := by rw [h2]

theorem grf_fderiv_jet2 {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ)
    (x v : EuclideanSpace ℝ (Fin n)) : ‖fderiv ℝ (jet 2 φ) x v‖ ≤ ‖jet 3 φ x‖ * ‖v‖ := by
  have h := grf_clm_norm_sq_le (fderiv ℝ (jet 2 φ) x) v
  have h3 : ∑ i, ‖fderiv ℝ (jet 2 φ) x (eb i)‖ ^ 2 = ‖jet 3 φ x‖ ^ 2 :=
    sum_norm_fderiv_jet_sq hφ 2 x
  rw [h3] at h
  refine (sq_le_sq₀ (norm_nonneg _) (by positivity)).1 ?_
  calc _ ≤ _ := h
    _ = _ := by ring

theorem grf_fderiv_gradE {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ)
    (x v : EuclideanSpace ℝ (Fin n)) : ‖fderiv ℝ (gradE φ) x v‖ ≤ ‖jet 2 φ x‖ * ‖v‖ := by
  have h1 : fderiv ℝ (gradE φ) x v =
      WithLp.toLp 2 (fun b => ∑ a, v a * jet 2 φ x ![a, b]) := by
    rw [ir_lin_expand (fderiv ℝ (gradE φ) x) v]
    ext b
    simp [gb_gradE_comp hφ]
  rw [h1]
  exact sl_frob (jet 2 φ x) v

theorem grf_real_le {A B c : ℝ} (hA : 0 ≤ A) (hc : 0 ≤ c) (h : A ^ 2 ≤ (1 + B ^ 2) ^ 3 * c ^ 2) :
    A ≤ c * (1 + B ^ 2) ^ ((3 : ℝ) / 2) := by
  have e : ((1 + B ^ 2) ^ ((3 : ℝ) / 2)) ^ 2 = (1 + B ^ 2) ^ 3 := by
    rw [← Real.rpow_natCast, ← Real.rpow_mul (by positivity)]
    rw [show ((3 : ℝ) / 2 * ((2 : ℕ) : ℝ)) = ((3 : ℕ) : ℝ) by norm_num, Real.rpow_natCast]
  refine (sq_le_sq₀ hA (by positivity)).1 ?_
  rw [mul_pow, e]
  linarith

/-! ### Part 2: a mean value inequality with a singular left endpoint -/

theorem grf_mvt {W : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W] {w : ℝ → W} {b C : ℝ}
    (hb : 0 < b) (hcont : ContinuousOn w (Ico 0 b))
    (hder : ∀ s ∈ Ioo 0 b, ∃ w' : W, HasDerivAt w w' s ∧ ‖w'‖ ≤ C) (h0 : w 0 = 0) :
    ∀ s ∈ Ico 0 b, ‖w s‖ ≤ C * s := by
  intro s hs
  rcases hs.1.eq_or_lt with h | hpos
  · rw [← h, h0]; simp
  choose! w' hw' using hder
  have hlim : Tendsto (fun ε => ‖w ε‖ + C * s) (𝓝[>] 0) (𝓝 (C * s)) := by
    have h1 : ContinuousWithinAt w (Ioi 0) 0 :=
      (hcont 0 ⟨le_rfl, hb⟩).mono_of_mem_nhdsWithin
        (mem_of_superset (Ioo_mem_nhdsGT hb) fun t ht => ⟨ht.1.le, ht.2⟩)
    have h2 : Tendsto (fun ε => ‖w ε‖) (𝓝[>] 0) (𝓝 ‖w 0‖) :=
      (continuous_norm.tendsto _).comp h1
    rw [h0, norm_zero] at h2
    simpa using h2.add_const (C * s)
  refine ge_of_tendsto hlim ?_
  filter_upwards [Ioo_mem_nhdsGT hpos] with ε hε
  have hI : Icc ε s ⊆ Ico 0 b := fun t ht => ⟨hε.1.le.trans ht.1, lt_of_le_of_lt ht.2 hs.2⟩
  have hmv := norm_image_sub_le_of_norm_deriv_right_le_segment (f := w) (a := ε) (b := s)
    (f' := w') (C := C) (hcont.mono hI)
    (fun t ht => (hw' t ⟨hε.1.trans_le ht.1, lt_of_lt_of_le ht.2 hs.2.le⟩).1.hasDerivWithinAt)
    (fun t ht => (hw' t ⟨hε.1.trans_le ht.1, lt_of_lt_of_le ht.2 hs.2.le⟩).2) s
    ⟨hε.2.le, le_rfl⟩
  have h3 : ‖w s‖ ≤ ‖w ε‖ + ‖w s - w ε‖ := by
    have := norm_add_le (w ε) (w s - w ε)
    simpa using this
  have hC : 0 ≤ C := by
    have := (hw' (s / 2) ⟨by linarith, by linarith [hs.2]⟩).2
    exact (norm_nonneg _).trans this
  nlinarith [mul_nonneg hC hε.1.le]

/-! ### Part 3: the point `p` and the injectivity of `dFq` -/

theorem grf_dFq_inj (hF : IsLagImm ι J F) (q : M) : Function.Injective (dFq ι F q) :=
  hF.inj q _ (im_chartFwd_mem ι (mem_extChartAt_source q))

theorem grf_slopeZ_zero {L : EuclideanSpace ℝ (Fin n) →L[ℝ] V} (hL : verCL J U ∘L L = 0) :
    slopeZ J U L = 0 := by
  ext I
  rw [sl_slope_apply]
  have : verCL J U (L ((horCL U ∘L L).inverse (eb (I 0)))) = 0 := by
    have := DFunLike.congr_fun hL ((horCL U ∘L L).inverse (eb (I 0)))
    simpa using this
  rw [this]
  simp

theorem grf_mem_good (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    (hU : IsLagFrame J U) [FiniteDimensional ℝ V] (hdim : Module.finrank ℝ V = 2 * n)
    (hF : IsLagImm ι J F) {p : M} (hp0 : F p = y₀) (hpv : verCL J U ∘L dFq ι F p = 0) :
    p ∈ goodSet ι J U y₀ F 6 := by
  refine ⟨?_, ?_, ?_⟩
  · rw [hp0, horC_eq, sub_self, map_zero, norm_zero]; norm_num
  · rw [hp0, verC_eq, sub_self, map_zero, norm_zero]; norm_num
  · refine injective_horCL_of_slopeLE hJi hJJ hU hdim (c := 0) (grf_dFq_inj ι hF p) ?_
    intro a
    have := DFunLike.congr_fun hpv a
    simp only [ContinuousLinearMap.comp_apply, zero_apply] at this
    rw [this]
    simp

/-! ### Part 4: the local form along a lift -/

theorem grf_local (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    (hU : IsLagFrame J U) [FiniteDimensional ℝ V] (hdim : Module.finrank ℝ V = 2 * n)
    (hF : IsLagImm ι J F) {c₁ : ℝ}
    (hA : ∀ z : M, ∀ y ∈ chartDomE ι z, aNormSq J (F ∘ chartInv ι z) y ≤ c₁ ^ 2)
    {v : EuclideanSpace ℝ (Fin n)} {S : Set ℝ} {γ : ℝ → M} {s₁ : ℝ}
    (hγc : ContinuousWithinAt γ S s₁) (hγG : γ s₁ ∈ goodSet ι J U y₀ F 6)
    (hγf : ∀ s ∈ S, horC U y₀ (F (γ s)) = s • v) :
    ∃ φ : EuclideanSpace ℝ (Fin n) → ℝ, ContDiff ℝ ∞ φ ∧ ∀ᶠ s in 𝓝[S] s₁,
      slopeZ J U (dFq ι F (γ s)) = jet 2 φ (s • v) ∧ verC J U y₀ (F (γ s)) = gradE φ (s • v) ∧
      ‖jet 3 φ (s • v)‖ ^ 2 ≤ (1 + ‖jet 2 φ (s • v)‖ ^ 2) ^ 3 * c₁ ^ 2 := by
  obtain ⟨r, hr, sl, φ, hφ, hs, hmaps, hx, hev⟩ := exists_local_graph_M ι hJi hJJ hU hdim hF hγG
  refine ⟨φ, hφ, ?_⟩
  have e1 := Filter.Tendsto.eventually hγc hev
  filter_upwards [e1, self_mem_nhdsWithin] with s h hs'
  obtain ⟨-, hball, hinv⟩ := h
  rw [hγf s hs'] at hball hinv
  obtain ⟨-, h2, -, h4, h5⟩ := hx _ hball
  rw [hinv] at h2 h4
  refine ⟨h4, ?_, ?_⟩
  · rw [h2, verC_graphX hJi hJJ hU]
  · have h6 := norm_jet3_sq_le φ (s • v)
    have hA' := hA (γ s₁) _ (hmaps hball)
    rw [h5] at hA'
    exact h6.trans (mul_le_mul_of_nonneg_left hA' (by positivity))

/-! ### Part 5: the a priori estimate -/

theorem grf_apriori (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    (hU : IsLagFrame J U) [FiniteDimensional ℝ V] (hdim : Module.finrank ℝ V = 2 * n)
    (hF : IsLagImm ι J F) {c₀ c₁ : ℝ} (hc₀ : 0 < c₀) (hc₁ : 0 ≤ c₁)
    (hsmall : 6 * (c₁ * (1 + c₀ ^ 2) ^ ((3 : ℝ) / 2)) < c₀)
    (hA : ∀ z : M, ∀ y ∈ chartDomE ι z, aNormSq J (F ∘ chartInv ι z) y ≤ c₁ ^ 2)
    {p : M} (hp0 : F p = y₀) (hpv : verCL J U ∘L dFq ι F p = 0)
    {v : EuclideanSpace ℝ (Fin n)} (hv : ‖v‖ < 6) {b : ℝ} (hb : b ∈ Ioc (0 : ℝ) 1) {γ : ℝ → M}
    (hc : ContinuousOn γ (Ico 0 b)) (h0 : γ 0 = p)
    (hl : ∀ s ∈ Ico 0 b, γ s ∈ goodSet ι J U y₀ F 6 ∧ horC U y₀ (F (γ s)) = s • v) :
    ∀ s ∈ Ico 0 b, ‖slopeZ J U (dFq ι F (γ s))‖ < c₀ ∧
      ‖verC J U y₀ (F (γ s))‖ ≤ 6 * c₀ * s := by
  have hb0 : 0 < b := hb.1
  have key : ∀ s₁ ∈ Ico 0 b, ∃ φ : EuclideanSpace ℝ (Fin n) → ℝ, ContDiff ℝ ∞ φ ∧
      ∀ᶠ s in 𝓝[Ico 0 b] s₁,
      slopeZ J U (dFq ι F (γ s)) = jet 2 φ (s • v) ∧ verC J U y₀ (F (γ s)) = gradE φ (s • v) ∧
      ‖jet 3 φ (s • v)‖ ^ 2 ≤ (1 + ‖jet 2 φ (s • v)‖ ^ 2) ^ 3 * c₁ ^ 2 :=
    fun s₁ hs₁ => grf_local ι hJi hJJ hU hdim hF hA (hc s₁ hs₁) (hl s₁ hs₁).1
      (fun s hs => (hl s hs).2)
  have key' : ∀ s₁ ∈ Ioo 0 b, ∃ φ : EuclideanSpace ℝ (Fin n) → ℝ, ContDiff ℝ ∞ φ ∧
      ∀ᶠ s in 𝓝 s₁,
      slopeZ J U (dFq ι F (γ s)) = jet 2 φ (s • v) ∧ verC J U y₀ (F (γ s)) = gradE φ (s • v) ∧
      ‖jet 3 φ (s • v)‖ ^ 2 ≤ (1 + ‖jet 2 φ (s • v)‖ ^ 2) ^ 3 * c₁ ^ 2 := by
    intro s₁ hs₁
    obtain ⟨φ, hφ, ev⟩ := key s₁ ⟨hs₁.1.le, hs₁.2⟩
    exact ⟨φ, hφ, by rwa [nhdsWithin_eq_nhds.2 (Ico_mem_nhds hs₁.1 hs₁.2)] at ev⟩
  -- continuity of `u` and `w`
  have hucont : ContinuousOn (fun s => slopeZ J U (dFq ι F (γ s))) (Ico 0 b) := by
    intro s₁ hs₁
    obtain ⟨φ, hφ, ev⟩ := key s₁ hs₁
    have hcont : ContinuousAt (fun s : ℝ => jet 2 φ (s • v)) s₁ :=
      ((contDiff_jet hφ 2).continuous.comp (continuous_id.smul continuous_const)).continuousAt
    exact hcont.continuousWithinAt.congr_of_eventuallyEq (ev.mono fun s h => h.1)
      (ev.self_of_nhdsWithin hs₁).1
  have hwcont : ContinuousOn (fun s => verC J U y₀ (F (γ s))) (Ico 0 b) := by
    intro s₁ hs₁
    obtain ⟨φ, hφ, ev⟩ := key s₁ hs₁
    have hcont : ContinuousAt (fun s : ℝ => gradE φ (s • v)) s₁ :=
      ((contDiff_gradE hφ).continuous.comp (continuous_id.smul continuous_const)).continuousAt
    exact hcont.continuousWithinAt.congr_of_eventuallyEq (ev.mono fun s h => h.2.1)
      (ev.self_of_nhdsWithin hs₁).2.1
  -- the slope at `0`
  have hu0 : slopeZ J U (dFq ι F (γ 0)) = 0 := by
    rw [h0]; exact grf_slopeZ_zero hpv
  have hw0 : verC J U y₀ (F (γ 0)) = 0 := by
    rw [h0, hp0, verC_eq, sub_self, map_zero]
  -- derivative of `u`
  have hderu : ∀ s ∈ Ioo 0 b, ∃ u' : Z n 2,
      HasDerivAt (fun s => slopeZ J U (dFq ι F (γ s))) u' s ∧
      (‖slopeZ J U (dFq ι F (γ s))‖ ≤ c₀ →
        ‖u'‖ ≤ 6 * c₁ * (1 + ‖slopeZ J U (dFq ι F (γ s))‖ ^ 2) ^ ((3 : ℝ) / 2)) := by
    intro s hs
    obtain ⟨φ, hφ, ev⟩ := key' s hs
    refine ⟨fderiv ℝ (jet 2 φ) (s • v) v, ?_, ?_⟩
    · have h1 : HasFDerivAt (jet 2 φ) (fderiv ℝ (jet 2 φ) (s • v)) (s • v) :=
        ((contDiff_jet hφ 2).differentiable (by simp) _).hasFDerivAt
      have h2 : HasDerivAt (fun t : ℝ => t • v) v s := by
        simpa using (hasDerivAt_id s).smul_const v
      have h3 := h1.comp_hasDerivAt s h2
      exact h3.congr_of_eventuallyEq (ev.mono fun t ht => ht.1)
    · intro _
      obtain ⟨e1, -, e3⟩ := ev.self_of_nhds
      rw [e1]
      have h3 := grf_real_le (norm_nonneg (jet 3 φ (s • v))) hc₁ e3
      calc ‖fderiv ℝ (jet 2 φ) (s • v) v‖ ≤ ‖jet 3 φ (s • v)‖ * ‖v‖ := grf_fderiv_jet2 hφ _ _
        _ ≤ ‖jet 3 φ (s • v)‖ * 6 := mul_le_mul_of_nonneg_left hv.le (norm_nonneg _)
        _ ≤ (c₁ * (1 + ‖jet 2 φ (s • v)‖ ^ 2) ^ ((3 : ℝ) / 2)) * 6 :=
          mul_le_mul_of_nonneg_right h3 (by norm_num)
        _ = 6 * c₁ * (1 + ‖jet 2 φ (s • v)‖ ^ 2) ^ ((3 : ℝ) / 2) := by ring
  have hsmall' : b * (6 * c₁ * (1 + c₀ ^ 2) ^ ((3 : ℝ) / 2)) < c₀ := by
    calc b * (6 * c₁ * (1 + c₀ ^ 2) ^ ((3 : ℝ) / 2))
        ≤ 1 * (6 * c₁ * (1 + c₀ ^ 2) ^ ((3 : ℝ) / 2)) :=
          mul_le_mul_of_nonneg_right hb.2 (by positivity)
      _ = 6 * (c₁ * (1 + c₀ ^ 2) ^ ((3 : ℝ) / 2)) := by ring
      _ < c₀ := hsmall
  have hlt := norm_lt_of_deriv_le_Ico (u := fun s => slopeZ J U (dFq ι F (γ s))) (a := 6 * c₁)
    (c := c₀) (b := b) (by positivity) hc₀ hucont hderu hu0 hsmall'
  -- derivative of `w`
  have hderw : ∀ s ∈ Ioo 0 b, ∃ w' : EuclideanSpace ℝ (Fin n),
      HasDerivAt (fun s => verC J U y₀ (F (γ s))) w' s ∧ ‖w'‖ ≤ 6 * c₀ := by
    intro s hs
    obtain ⟨φ, hφ, ev⟩ := key' s hs
    refine ⟨fderiv ℝ (gradE φ) (s • v) v, ?_, ?_⟩
    · have h1 : HasFDerivAt (gradE φ) (fderiv ℝ (gradE φ) (s • v)) (s • v) :=
        ((contDiff_gradE hφ).differentiable (by simp) _).hasFDerivAt
      have h2 : HasDerivAt (fun t : ℝ => t • v) v s := by
        simpa using (hasDerivAt_id s).smul_const v
      have h3 := h1.comp_hasDerivAt s h2
      exact h3.congr_of_eventuallyEq (ev.mono fun t ht => ht.2.1)
    · obtain ⟨e1, -, -⟩ := ev.self_of_nhds
      have h4 : ‖jet 2 φ (s • v)‖ < c₀ := by
        rw [← e1]; exact hlt s ⟨hs.1.le, hs.2⟩
      calc ‖fderiv ℝ (gradE φ) (s • v) v‖ ≤ ‖jet 2 φ (s • v)‖ * ‖v‖ := grf_fderiv_gradE hφ _ _
        _ ≤ c₀ * 6 := mul_le_mul h4.le hv.le (norm_nonneg _) hc₀.le
        _ = 6 * c₀ := by ring
  have hmvt := grf_mvt hb0 hwcont hderw hw0
  intro s hs
  exact ⟨hlt s hs, hmvt s hs⟩

/-! ### Part 6: the radial section -/

theorem exists_radial_section [T2Space M] [CompactSpace M]
    (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v) (hU : IsLagFrame J U)
    [FiniteDimensional ℝ V] (hdim : Module.finrank ℝ V = 2 * n) (hF : IsLagImm ι J F)
    {c₀ c₁ : ℝ} (hc₀ : 0 < c₀) (hc₀1 : 6 * c₀ < 1) (hc₁ : 0 ≤ c₁)
    (hsmall : 6 * (c₁ * (1 + c₀ ^ 2) ^ ((3 : ℝ) / 2)) < c₀)
    (hA : ∀ z : M, ∀ y ∈ chartDomE ι z, aNormSq J (F ∘ chartInv ι z) y ≤ c₁ ^ 2)
    {p : M} (hp0 : F p = y₀) (hpv : verCL J U ∘L dFq ι F p = 0) :
    ∃ σ : EuclideanSpace ℝ (Fin n) → M, ContinuousOn σ (ball 0 6) ∧ σ 0 = p ∧
      (∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6,
        σ x ∈ goodSet ι J U y₀ F 6 ∧ horC U y₀ (F (σ x)) = x ∧
        ‖verC J U y₀ (F (σ x))‖ ≤ 6 * c₀ ∧ ‖slopeZ J U (dFq ι F (σ x))‖ < c₀) ∧
      σ '' ball 0 6 = connectedComponentIn (goodSet ι J U y₀ F 6) p := by
  classical
  have hFc := grl_continuous_F ι hF
  have hGo : IsOpen (goodSet ι J U y₀ F 6) := isOpen_goodSet ι hF 6
  have hf : IsLocalHomeomorphOn (fun q => horC U y₀ (F q)) (goodSet ι J U y₀ F 6) :=
    isLocalHomeomorphOn_horC_goodSet ι hJi hJJ hU hdim hF 6
  have hpG : p ∈ goodSet ι J U y₀ F 6 := grf_mem_good ι hJi hJJ hU hdim hF hp0 hpv
  have hfp : horC U y₀ (F p) = 0 := by rw [hp0, horC_eq, sub_self, map_zero]
  have hmaps : MapsTo (fun q => horC U y₀ (F q)) (goodSet ι J U y₀ F 6)
      (ball (0 : EuclideanSpace ℝ (Fin n)) 6) := fun q hq => mem_ball_zero_iff.2 hq.1
  have hlift : ∀ v ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6, ∃ γ : ℝ → M,
      ContinuousOn γ (Icc 0 1) ∧ γ 0 = p ∧
      ∀ s ∈ Icc (0 : ℝ) 1, γ s ∈ goodSet ι J U y₀ F 6 ∧ horC U y₀ (F (γ s)) = s • v := by
    intro v hv
    have hv' : ‖v‖ < 6 := mem_ball_zero_iff.1 hv
    set K : Set M := {q | ‖horC U y₀ (F q)‖ ≤ ‖v‖ ∧ ‖verC J U y₀ (F q)‖ ≤ 6 * c₀ ∧
      slopeLE J U c₀ (dFq ι F q)} with hK
    have hKc : IsClosed K :=
      (isClosed_le (continuous_norm.comp (grl_continuous_horC.comp hFc)) continuous_const).inter
        ((isClosed_le (continuous_norm.comp (grl_continuous_verC.comp hFc)) continuous_const).inter
          (isClosed_slopeLE_dFq ι hF c₀))
    have hKG : K ⊆ goodSet ι J U y₀ F 6 := by
      rintro q ⟨h1, h2, h3⟩
      exact ⟨lt_of_le_of_lt h1 hv', lt_of_le_of_lt h2 (by linarith),
        injective_horCL_of_slopeLE hJi hJJ hU hdim (grf_dFq_inj ι hF q) h3⟩
    refine exists_ray_lift hGo hf hpG hfp hKc.isCompact hKG ?_
    intro b hb γ hc h0 hl s hs
    obtain ⟨h1, h2⟩ := grf_apriori ι hJi hJJ hU hdim hF hc₀ hc₁ hsmall hA hp0 hpv hv' hb hc h0 hl s hs
    have hs1 : s ≤ 1 := hs.2.le.trans hb.2
    refine ⟨?_, ?_, slopeLE_of_norm_slopeZ_le hJi hJJ hU hdim (hl s hs).1.2.2 h1.le⟩
    · rw [(hl s hs).2, norm_smul, Real.norm_eq_abs, abs_of_nonneg hs.1]
      exact mul_le_of_le_one_left (norm_nonneg v) hs1
    · exact h2.trans (mul_le_of_le_one_right (by positivity) hs1)
  obtain ⟨σ, hσc, hσ0, hσs, hσim⟩ := radial_lift hGo hf hpG hfp hmaps hlift
  refine ⟨σ, hσc, hσ0, fun x hx => ?_, hσim⟩
  obtain ⟨hxG, hxf⟩ := hσs x hx
  have hxn : ‖x‖ < 6 := mem_ball_zero_iff.1 hx
  set t : ℝ := 1 + (6 - ‖x‖) / 12 with ht
  have ht1 : 1 < t := by rw [ht]; linarith
  have hvx : ‖t • x‖ < 6 := by
    rw [norm_smul, Real.norm_eq_abs, abs_of_pos (by linarith)]
    rw [ht]
    nlinarith [norm_nonneg x]
  have hsv : ∀ s ∈ Ico (0 : ℝ) 1, s • (t • x) ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6 := by
    intro s hs
    rw [mem_ball_zero_iff, norm_smul, Real.norm_eq_abs, abs_of_nonneg hs.1]
    calc s * ‖t • x‖ ≤ 1 * ‖t • x‖ := mul_le_mul_of_nonneg_right hs.2.le (norm_nonneg _)
      _ < 6 := by linarith
  have hγc : ContinuousOn (fun s : ℝ => σ (s • (t • x))) (Ico 0 1) :=
    hσc.comp (continuous_id.smul continuous_const).continuousOn hsv
  have h0 : (fun s : ℝ => σ (s • (t • x))) 0 = p := by simpa using hσ0
  have hl : ∀ s ∈ Ico (0 : ℝ) 1, (fun s : ℝ => σ (s • (t • x))) s ∈ goodSet ι J U y₀ F 6 ∧
      horC U y₀ (F ((fun s : ℝ => σ (s • (t • x))) s)) = s • (t • x) :=
    fun s hs => hσs _ (hsv s hs)
  have ht0 : 0 < t := by linarith
  have hs1 : 1 / t ∈ Ico (0 : ℝ) 1 := ⟨by positivity, by rw [div_lt_one ht0]; exact ht1⟩
  obtain ⟨h1, h2⟩ := grf_apriori ι hJi hJJ hU hdim hF hc₀ hc₁ hsmall hA hp0 hpv hvx
    (b := 1) ⟨one_pos, le_rfl⟩ hγc h0 hl (1 / t) hs1
  have hx' : (1 / t) • (t • x) = x := by
    rw [smul_smul, one_div, inv_mul_cancel₀ ht0.ne', one_smul]
  simp only [hx'] at h1 h2
  refine ⟨hxG, hxf, ?_, h1⟩
  exact h2.trans (mul_le_of_le_one_right (by positivity) ((div_le_one ht0).2 ht1.le))

/-- info: 'CWFlow.Extension.exists_radial_section' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_radial_section

end CWFlow.Extension
