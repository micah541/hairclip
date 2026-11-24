import CWFlow.Extension.GraphRadiusLift

/-!
# The graph radius theorem

A Lagrangian immersion `F : M → V` of a compact manifold with small second fundamental form is, near any
point `p` whose tangent plane is horizontal, the graph of the gradient of a function `φ` over the ball of radius
`6` in the Lagrangian plane `y₀ + U ℝⁿ`.  The parametrization `σ` of the corresponding open piece of `M` is a
chart-compatible parametrization (`IsParam`), and the images of the smaller balls are the connected components
of the sets `{‖horC‖ < ρ, ‖verC‖ < 1}`.

The continuous section `σ` with the slope bounds comes from `exists_radial_section`; here it is upgraded:

* `gr_local`: near each point of the ball, `σ = chartInv ∘ s_loc` for a local graph representation;
* `gr_isParam`: `σ` is a parametrization;
* `gr_component`: images of balls are connected components;
* `graph_radius`: the final statement.
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

local notation "𝔼" => EuclideanSpace ℝ (Fin n)

/-! ### Part 1: jets depend on the first derivative germ -/

theorem gr_fderiv_eq_of_gradE {φ ψ : 𝔼 → ℝ} {y : 𝔼} (h : gradE φ y = gradE ψ y) :
    fderiv ℝ φ y = fderiv ℝ ψ y := by
  refine ContinuousLinearMap.ext fun v => ?_
  rw [ir_lin_expand (fderiv ℝ φ y) v, ir_lin_expand (fderiv ℝ ψ y) v]
  refine Finset.sum_congr rfl fun i _ => ?_
  have : fderiv ℝ φ y (eb i) = fderiv ℝ ψ y (eb i) :=
    congrArg (fun w : 𝔼 => w i) h
  rw [this]

theorem gr_jet2_congr {φ ψ : 𝔼 → ℝ} {x : 𝔼} (h : fderiv ℝ φ =ᶠ[𝓝 x] fderiv ℝ ψ) :
    jet 2 φ x = jet 2 ψ x := by
  have h1 := (Filter.EventuallyEq.iteratedFDeriv ℝ h 1).eq_of_nhds
  have h2 : iteratedFDeriv ℝ 2 φ x = iteratedFDeriv ℝ 2 ψ x := by
    show iteratedFDeriv ℝ (1 + 1) φ x = iteratedFDeriv ℝ (1 + 1) ψ x
    rw [iteratedFDeriv_succ_eq_comp_right, iteratedFDeriv_succ_eq_comp_right]
    simp only [Function.comp_apply, h1]
  unfold jet
  rw [h2]

/-! ### Part 2: the local form of the section -/

theorem gr_local (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    (hU : IsLagFrame J U) [FiniteDimensional ℝ V] (hdim : Module.finrank ℝ V = 2 * n)
    (hF : IsLagImm ι J F) {σ : 𝔼 → M} (hσc : ContinuousOn σ (ball 0 6))
    (hσG : ∀ x ∈ ball (0 : 𝔼) 6, σ x ∈ goodSet ι J U y₀ F 6)
    (hσf : ∀ x ∈ ball (0 : 𝔼) 6, horC U y₀ (F (σ x)) = x) {x : 𝔼} (hx : x ∈ ball (0 : 𝔼) 6) :
    ∃ (r : ℝ) (s : 𝔼 → 𝔼) (φ : 𝔼 → ℝ), 0 < r ∧ ContDiff ℝ ∞ φ ∧ ContDiffOn ℝ ∞ s (ball x r) ∧
      MapsTo s (ball x r) (chartDomE ι (σ x)) ∧
      (∀ x' ∈ ball x r, F (chartInv ι (σ x) (s x')) = graphX J U y₀ φ x') ∧
      slopeZ J U (dFq ι F (σ x)) = jet 2 φ x ∧
      ∀ᶠ x' in 𝓝 x, x' ∈ ball x r ∧ x' ∈ ball (0 : 𝔼) 6 ∧
        σ x' ∈ (extChartAt 𝓘(ℝ, E) (σ x)).source ∧ σ x' = chartInv ι (σ x) (s x') := by
  classical
  have hc : horC U y₀ (F (σ x)) = x := hσf x hx
  obtain ⟨r, hr, s, φ, hφ, hs, hmaps, hx', hev⟩ :=
    exists_local_graph_M ι hJi hJJ hU hdim hF (hσG x hx)
  rw [hc] at hs hmaps hx' hev
  have hB : ball (0 : 𝔼) 6 ∈ 𝓝 x := isOpen_ball.mem_nhds hx
  have hev' := Filter.Tendsto.eventually (hσc.continuousAt hB) hev
  have hev2 : ∀ᶠ x' in 𝓝 x, x' ∈ ball x r ∧ x' ∈ ball (0 : 𝔼) 6 ∧
      σ x' ∈ (extChartAt 𝓘(ℝ, E) (σ x)).source ∧ σ x' = chartInv ι (σ x) (s x') := by
    filter_upwards [hev', hB] with x' h h'
    obtain ⟨h1, h2, h3⟩ := h
    rw [hσf x' h'] at h2 h3
    exact ⟨h2, h', h1, h3.symm⟩
  refine ⟨r, s, φ, hr, hφ, hs, hmaps, ?_, ?_, hev2⟩
  · intro x' hx'b
    exact (hx' x' hx'b).2.1
  · obtain ⟨-, -, -, h4⟩ := hev2.self_of_nhds
    have h5 := (hx' x (mem_ball_self hr)).2.2.2.1
    rw [← h4] at h5
    exact h5

theorem gr_contDiffAt (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    (hU : IsLagFrame J U) [FiniteDimensional ℝ V] (hdim : Module.finrank ℝ V = 2 * n)
    (hF : IsLagImm ι J F) {σ : 𝔼 → M} (hσc : ContinuousOn σ (ball 0 6))
    (hσG : ∀ x ∈ ball (0 : 𝔼) 6, σ x ∈ goodSet ι J U y₀ F 6)
    (hσf : ∀ x ∈ ball (0 : 𝔼) 6, horC U y₀ (F (σ x)) = x) {y : 𝔼} (hy : y ∈ ball (0 : 𝔼) 6)
    (z : M) (hz : σ y ∈ (extChartAt 𝓘(ℝ, E) z).source) :
    ContDiffAt ℝ ∞ (chartFwd ι z ∘ σ) y := by
  obtain ⟨r, s, φ, hr, -, hs, hmaps, -, -, hev⟩ := gr_local ι hJi hJJ hU hdim hF hσc hσG hσf hy
  have hsy : s y ∈ chartDomE ι (σ y) := hmaps (mem_ball_self hr)
  obtain ⟨-, -, -, h4⟩ := hev.self_of_nhds
  have hWo : IsOpen (chartDomE ι (σ y) ∩
      chartInv ι (σ y) ⁻¹' (extChartAt 𝓘(ℝ, E) z).source) :=
    (im_continuousOn_chartInv ι _).isOpen_inter_preimage (im_isOpen_chartDomE ι _)
      (isOpen_extChartAt_source z)
  have hmem : s y ∈ chartDomE ι (σ y) ∩
      chartInv ι (σ y) ⁻¹' (extChartAt 𝓘(ℝ, E) z).source :=
    ⟨hsy, by
      show chartInv ι (σ y) (s y) ∈ (extChartAt 𝓘(ℝ, E) z).source
      rw [← h4]; exact hz⟩
  have h1 : ContDiffAt ℝ ∞ (chartFwd ι z ∘ chartInv ι (σ y)) (s y) :=
    (im_contDiffOn_trans ι (σ y) z).contDiffAt (hWo.mem_nhds hmem)
  have h2 : ContDiffAt ℝ ∞ s y := hs.contDiffAt (isOpen_ball.mem_nhds (mem_ball_self hr))
  have h3 := h1.comp y h2
  refine h3.congr_of_eventuallyEq ?_
  filter_upwards [hev] with w hw
  simp only [Function.comp_apply]
  rw [hw.2.2.2]

/-! ### Part 3: the section is a parametrization -/

theorem gr_isParam (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    (hU : IsLagFrame J U) [FiniteDimensional ℝ V] (hdim : Module.finrank ℝ V = 2 * n)
    (hF : IsLagImm ι J F) {σ : 𝔼 → M} (hσc : ContinuousOn σ (ball 0 6))
    (hσG : ∀ x ∈ ball (0 : 𝔼) 6, σ x ∈ goodSet ι J U y₀ F 6)
    (hσf : ∀ x ∈ ball (0 : 𝔼) 6, horC U y₀ (F (σ x)) = x) {p : M}
    (hσim : σ '' ball (0 : 𝔼) 6 = connectedComponentIn (goodSet ι J U y₀ F 6) p) :
    IsParam ι (ball (0 : 𝔼) 6) σ where
  isOpen := isOpen_ball
  injOn := fun a ha b hb h =>
    calc a = horC U y₀ (F (σ a)) := (hσf a ha).symm
      _ = horC U y₀ (F (σ b)) := by rw [h]
      _ = b := hσf b hb
  isOpen_image := by
    rw [hσim]
    have : LocallyConnectedSpace E := inferInstance
    have : LocallyConnectedSpace M := ChartedSpace.locallyConnectedSpace E M
    exact (isOpen_goodSet ι hF 6).connectedComponentIn
  cont := hσc
  inv := ⟨fun q => horC U y₀ (F q),
    (grl_continuous_horC.comp (grl_continuous_F ι hF)).continuousOn, hσf⟩
  trans_smooth z := fun y hy =>
    (gr_contDiffAt ι hJi hJJ hU hdim hF hσc hσG hσf hy.1 z hy.2).contDiffWithinAt
  trans_det z y hy := by
    obtain ⟨hyB, hyS⟩ := hy
    have hτd : DifferentiableAt ℝ (chartFwd ι z ∘ σ) y :=
      (gr_contDiffAt ι hJi hJJ hU hdim hF hσc hσG hσf hyB z hyS).differentiableAt (by simp)
    have hyD : (chartFwd ι z ∘ σ) y ∈ chartDomE ι z := im_chartFwd_mem ι hyS
    have hFd : DifferentiableAt ℝ (F ∘ chartInv ι z) ((chartFwd ι z ∘ σ) y) :=
      ((hF.smooth z).contDiffAt ((im_isOpen_chartDomE ι z).mem_nhds hyD)).differentiableAt
        (by simp)
    have hΦ : DifferentiableAt ℝ (horC U y₀ ∘ (F ∘ chartInv ι z)) ((chartFwd ι z ∘ σ) y) :=
      (hasFDerivAt_horC ((F ∘ chartInv ι z) ((chartFwd ι z ∘ σ) y))).differentiableAt.comp _ hFd
    have hev : (horC U y₀ ∘ (F ∘ chartInv ι z)) ∘ (chartFwd ι z ∘ σ) =ᶠ[𝓝 y] id := by
      have h1 : ∀ᶠ w in 𝓝 y, w ∈ ball (0 : 𝔼) 6 := isOpen_ball.mem_nhds hyB
      have h2 : ∀ᶠ w in 𝓝 y, σ w ∈ (extChartAt 𝓘(ℝ, E) z).source :=
        (hσc.continuousAt (isOpen_ball.mem_nhds hyB)).preimage_mem_nhds
          ((isOpen_extChartAt_source z).mem_nhds hyS)
      filter_upwards [h1, h2] with w a b
      simp only [Function.comp_apply, id]
      rw [im_chartInv_chartFwd ι b, hσf w a]
    have h1 := fderiv_comp y hΦ hτd
    have h2 : fderiv ℝ ((horC U y₀ ∘ (F ∘ chartInv ι z)) ∘ (chartFwd ι z ∘ σ)) y =
        ContinuousLinearMap.id ℝ _ := by
      rw [hev.fderiv_eq]; exact fderiv_id
    refine sl_det_ne ?_
    intro a b hab
    have key : ∀ u, fderiv ℝ (horC U y₀ ∘ (F ∘ chartInv ι z)) ((chartFwd ι z ∘ σ) y)
        (fderiv ℝ (chartFwd ι z ∘ σ) y u) = u := fun u => by
      have := congrArg (fun L : 𝔼 →L[ℝ] 𝔼 => L u) (h1.symm.trans h2)
      simpa using this
    rw [← key a, ← key b, hab]

/-! ### Part 4: components -/

theorem gr_component {M' : Type*} [TopologicalSpace M'] [T2Space M'] {f : M' → 𝔼}
    (hf : Continuous f) {σ : 𝔼 → M'} (hσc : ContinuousOn σ (ball 0 6))
    (hσf : ∀ x ∈ ball (0 : 𝔼) 6, f (σ x) = x) (hopen : IsOpen (σ '' ball (0 : 𝔼) 6))
    {S : Set M'} {ρ : ℝ} (hρ : 0 < ρ) (hρ6 : ρ ≤ 6) (hSf : ∀ q ∈ S, ‖f q‖ < ρ)
    (hCS : σ '' ball (0 : 𝔼) ρ ⊆ S) :
    σ '' ball (0 : 𝔼) ρ = connectedComponentIn S (σ 0) := by
  set C := σ '' ball (0 : 𝔼) ρ with hC
  have hρ6' : ball (0 : 𝔼) ρ ⊆ ball 0 6 := ball_subset_ball hρ6
  have hCo : IsOpen C := by
    have : C = σ '' ball 0 6 ∩ f ⁻¹' ball 0 ρ := by
      ext q
      constructor
      · rintro ⟨x, hx, rfl⟩
        exact ⟨⟨x, hρ6' hx, rfl⟩, by rw [mem_preimage, hσf x (hρ6' hx)]; exact hx⟩
      · rintro ⟨⟨x, hx, rfl⟩, h2⟩
        refine ⟨x, ?_, rfl⟩
        rw [mem_preimage, hσf x hx] at h2
        exact h2
    rw [this]
    exact hopen.inter (isOpen_ball.preimage hf)
  have hCconn : IsPreconnected C :=
    ((convex_ball (0 : 𝔼) ρ).isPreconnected).image σ (hσc.mono hρ6')
  have hp0 : σ 0 ∈ C := ⟨0, mem_ball_self hρ, rfl⟩
  have hCK : C ⊆ connectedComponentIn S (σ 0) := hCconn.subset_connectedComponentIn hp0 hCS
  have hcl : ∀ q ∈ S, q ∈ closure C → q ∈ C := by
    intro q hqS hq'
    have hq6 : f q ∈ ball (0 : 𝔼) 6 := hρ6' (mem_ball_zero_iff.2 (hSf q hqS))
    have hg : ContinuousAt (fun q' => σ (f q')) q :=
      (hσc.continuousAt (isOpen_ball.mem_nhds hq6)).comp hf.continuousAt
    have : (𝓝[C] q).NeBot := mem_closure_iff_nhdsWithin_neBot.1 hq'
    have h1 : Tendsto (fun q' => σ (f q')) (𝓝[C] q) (𝓝 (σ (f q))) :=
      hg.tendsto.mono_left nhdsWithin_le_nhds
    have h2 : Tendsto id (𝓝[C] q) (𝓝 q) := tendsto_id.mono_left nhdsWithin_le_nhds
    have h3 : (fun q' => σ (f q')) =ᶠ[𝓝[C] q] id := by
      filter_upwards [self_mem_nhdsWithin] with q' hq''
      obtain ⟨x', hx', rfl⟩ := hq''
      simp only [id]
      rw [hσf x' (hρ6' hx')]
    have h4 := tendsto_nhds_unique (Filter.Tendsto.congr' h3 h1) h2
    exact ⟨f q, mem_ball_zero_iff.2 (hSf q hqS), h4⟩
  refine hCK.antisymm ?_
  have hv : IsOpen (closure C)ᶜ := isClosed_closure.isOpen_compl
  have hdisj : Disjoint C (closure C)ᶜ :=
    Set.disjoint_left.2 fun q hq hq' => hq' (subset_closure hq)
  have hsuv : connectedComponentIn S (σ 0) ⊆ C ∪ (closure C)ᶜ := by
    intro q hq
    by_cases h : q ∈ closure C
    · exact Or.inl (hcl q (connectedComponentIn_subset S (σ 0) hq) h)
    · exact Or.inr h
  exact (isPreconnected_connectedComponentIn).subset_left_of_subset_union hCo hv hdisj hsuv
    ⟨σ 0, mem_connectedComponentIn (hCS hp0), hp0⟩

/-! ### Part 5: the theorem -/

theorem graph_radius [T2Space M] [CompactSpace M]
    (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v) (hU : IsLagFrame J U)
    [FiniteDimensional ℝ V] (hdim : Module.finrank ℝ V = 2 * n) (hF : IsLagImm ι J F)
    {c₀ c₁ : ℝ} (hc₀ : 0 < c₀) (hc₀1 : 6 * c₀ < 1) (hc₁ : 0 ≤ c₁)
    (hsmall : 6 * (c₁ * (1 + c₀ ^ 2) ^ ((3 : ℝ) / 2)) < c₀)
    (hA : ∀ z : M, ∀ y ∈ chartDomE ι z, aNormSq J (F ∘ chartInv ι z) y ≤ c₁ ^ 2)
    {p : M} (hp0 : F p = y₀) (hpv : verCL J U ∘L dFq ι F p = 0) :
    ∃ (σ : EuclideanSpace ℝ (Fin n) → M) (φ : EuclideanSpace ℝ (Fin n) → ℝ),
      IsParam ι (ball 0 6) σ ∧ σ 0 = p ∧ ContDiffOn ℝ ∞ φ (ball 0 6) ∧ φ 0 = 0 ∧
      (∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6,
        F (σ x) = graphX J U y₀ φ x ∧ horC U y₀ (F (σ x)) = x ∧
        ‖jet 2 φ x‖ ≤ c₀ ∧ ‖gradE φ x‖ ≤ 6 * c₀) ∧
      ∀ ρ, 0 < ρ → ρ ≤ 6 →
        σ '' ball 0 ρ = connectedComponentIn
          {q | ‖horC U y₀ (F q)‖ < ρ ∧ ‖verC J U y₀ (F q)‖ < 1} p := by
  classical
  obtain ⟨σ, hσc, hσ0, hσb, hσim⟩ :=
    exists_radial_section ι hJi hJJ hU hdim hF hc₀ hc₀1 hc₁ hsmall hA hp0 hpv
  have hσG : ∀ x ∈ ball (0 : 𝔼) 6, σ x ∈ goodSet ι J U y₀ F 6 := fun x hx => (hσb x hx).1
  have hσf : ∀ x ∈ ball (0 : 𝔼) 6, horC U y₀ (F (σ x)) = x := fun x hx => (hσb x hx).2.1
  have hP : IsParam ι (ball (0 : 𝔼) 6) σ :=
    gr_isParam ι hJi hJJ hU hdim hF hσc hσG hσf hσim
  obtain ⟨h1, h2, h3⟩ := hF.comp_param ι hP
  have hfd : ∀ y ∈ ball (0 : 𝔼) 6, horCL U ∘L fderiv ℝ (F ∘ σ) y = ContinuousLinearMap.id ℝ 𝔼 := by
    intro y hy
    have hXd : DifferentiableAt ℝ (F ∘ σ) y :=
      (h1.contDiffAt (isOpen_ball.mem_nhds hy)).differentiableAt (by simp)
    have h4 : HasFDerivAt (fun w => horC U y₀ ((F ∘ σ) w)) (horCL U ∘L fderiv ℝ (F ∘ σ) y) y :=
      (hasFDerivAt_horC ((F ∘ σ) y)).comp y hXd.hasFDerivAt
    have h5 : HasFDerivAt (fun w => horC U y₀ ((F ∘ σ) w)) (ContinuousLinearMap.id ℝ 𝔼) y :=
      (hasFDerivAt_id y).congr_of_eventuallyEq (by
        filter_upwards [isOpen_ball.mem_nhds hy] with w hw
        exact hσf w hw)
    exact h4.unique h5
  have hpinj : ∀ y ∈ ball (0 : 𝔼) 6, Function.Injective (horCL U ∘L fderiv ℝ (F ∘ σ) y) := by
    intro y hy a b hab
    have e := hfd y hy
    rw [e] at hab
    simpa using hab
  have hinjOn : Set.InjOn (fun y => horC U y₀ ((F ∘ σ) y)) (ball (0 : 𝔼) 6) := by
    intro a ha b hb h
    have h' : horC U y₀ (F (σ a)) = horC U y₀ (F (σ b)) := h
    rw [hσf a ha, hσf b hb] at h'
    exact h'
  have hball : ball (0 : 𝔼) 6 ⊆ (fun y => horC U y₀ ((F ∘ σ) y)) '' ball (0 : 𝔼) 6 :=
    fun x hx => ⟨x, hx, hσf x hx⟩
  obtain ⟨s, φ, hsC, hsmaps, -, hsr, hφC, hφ0, hXs⟩ :=
    exists_graph_of_proj hJi hJJ hU hdim (X := F ∘ σ) (W := ball (0 : 𝔼) 6) isOpen_ball h1 h3
      hpinj hinjOn (c := 0) (R := 6) hball
  have hsx : ∀ x ∈ ball (0 : 𝔼) 6, s x = x := by
    intro x hx
    have e : horC U y₀ ((F ∘ σ) x) = x := hσf x hx
    have := hsr x hx (by rw [e]; exact hx)
    rwa [e] at this
  have hFσ : ∀ x ∈ ball (0 : 𝔼) 6, F (σ x) = graphX J U y₀ φ x := by
    intro x hx
    have := hXs x hx
    rw [hsx x hx] at this
    exact this
  have hjet : ∀ x ∈ ball (0 : 𝔼) 6, jet 2 φ x = slopeZ J U (dFq ι F (σ x)) := by
    intro x hx
    obtain ⟨r, s', ψ, hr, hψ, hs', hmaps', hgr, hsl', hev⟩ :=
      gr_local ι hJi hJJ hU hdim hF hσc hσG hσf hx
    rw [hsl']
    have hgrad : fderiv ℝ φ =ᶠ[𝓝 x] fderiv ℝ ψ := by
      filter_upwards [hev] with x' hx'
      obtain ⟨k1, k2, k3, k4⟩ := hx'
      have e1 : graphX J U y₀ φ x' = graphX J U y₀ ψ x' := by
        rw [← hFσ x' k2, ← hgr x' k1, ← k4]
      have e2 : gradE φ x' = gradE ψ x' :=
        calc gradE φ x' = verC J U y₀ (graphX J U y₀ φ x') := (verC_graphX hJi hJJ hU φ x').symm
          _ = verC J U y₀ (graphX J U y₀ ψ x') := by rw [e1]
          _ = gradE ψ x' := verC_graphX hJi hJJ hU ψ x'
      exact gr_fderiv_eq_of_gradE e2
    exact gr_jet2_congr hgrad
  have hbig : ∀ ρ : ℝ, ρ ≤ 6 → ball (0 : 𝔼) ρ ⊆ ball 0 6 := fun ρ hρ => ball_subset_ball hρ
  refine ⟨σ, φ, hP, hσ0, hφC, hφ0, ?_, ?_⟩
  · intro x hx
    obtain ⟨hG, hf, hv, hs⟩ := hσb x hx
    have hgr : gradE φ x = verC J U y₀ (F (σ x)) := by
      rw [hFσ x hx, verC_graphX hJi hJJ hU]
    refine ⟨hFσ x hx, hf, ?_, ?_⟩
    · rw [hjet x hx]; exact hs.le
    · rw [hgr]; exact hv
  · intro ρ hρ hρ6
    have hfc : Continuous (fun q => horC U y₀ (F q)) :=
      grl_continuous_horC.comp (grl_continuous_F ι hF)
    have key := gr_component (f := fun q => horC U y₀ (F q)) (M' := M) hfc hσc hσf hP.isOpen_image
      (S := {q | ‖horC U y₀ (F q)‖ < ρ ∧ ‖verC J U y₀ (F q)‖ < 1}) hρ hρ6
      (fun q hq => hq.1) (by
        rintro _ ⟨x, hx, rfl⟩
        refine ⟨?_, lt_of_le_of_lt (hσb x (hbig ρ hρ6 hx)).2.2.1 hc₀1⟩
        rw [hσf x (hbig ρ hρ6 hx)]
        exact mem_ball_zero_iff.1 hx)
    rw [hσ0] at key
    exact key

/-- info: 'CWFlow.Extension.graph_radius' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms graph_radius

end CWFlow.Extension
