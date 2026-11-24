import CWFlow.Flow.GlobalGauge
import CWFlow.Flow.GaugePotential

/-!
# A flow is a graph over its initial immersion for a short time

`exists_gauge_graph` (and `exists_gauge_graphS` for velocity `κ • flowVel`, with `κ * gaugePot` in place of
`gaugePot`): let `Φ` be a flow with `Φ 0 = F`. For a short time `δ` there is a jointly smooth family
of maps `ψ_t`, `ψ_0 = id`, such that `Φ_t ∘ ψ_t = graphImm (f_t)` for the jointly smooth potential
`f_t = gaugePot E J Φ ψ t`, with `f_0 = 0` and `∂ₜf_t = −div(JH)[Φ_t] ∘ ψ_t`.
-/

open Set Filter
open scoped InnerProductSpace RealInnerProductSpace Topology ContDiff Manifold

noncomputable section

set_option linter.unusedSectionVars false

namespace CWFlow

namespace Bridge

universe u

variable {E : Type u} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
  [NeZero (Module.finrank ℝ E)] [CompleteSpace E]
variable {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
  [CompactSpace M] [T2Space M] [SigmaCompactSpace M]
variable {V : Type u} [NormedAddCommGroup V] [InnerProductSpace ℝ V] [FiniteDimensional ℝ V]
variable {J : V →L[ℝ] V} {F : M → V}

/-- **A flow from `F` with velocity `κ • flowVel` is, for a short time and up to reparametrization, the graph of
an exact section** (the differential of `κ * gaugePot`). -/
theorem exists_gauge_graphS {κ : ℝ} (hJi : ∀ u v, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v, J (J v) = -v)
    (hF : IsLagImmersion E J F) (hdim : Module.finrank ℝ V = 2 * Module.finrank ℝ E)
    {T : ℝ} (hT : 0 < T) {Φ : ℝ → M → V} (hΦ : IsFlowS E J κ T Φ) (hΦ0 : Φ 0 = F) :
    ∃ δ, 0 < δ ∧ δ ≤ T ∧ ∃ ψ : ℝ → M → M,
      ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, E) ∞ (fun q : ℝ × M => ψ q.1 q.2) (Icc 0 δ ×ˢ univ) ∧
      (∀ x, ψ 0 x = x) ∧
      ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, ℝ) ∞
        (fun q : ℝ × M => κ * gaugePot E J Φ ψ q.1 q.2) (Icc 0 δ ×ˢ univ) ∧
      (∀ t ∈ Icc (0 : ℝ) δ, ∀ x : M, HasDerivWithinAt (fun s => κ * gaugePot E J Φ ψ s x)
        (-(κ * flowPotAt E J Φ t (ψ t x))) (Icc 0 δ) t) ∧
      ∀ t ∈ Icc (0 : ℝ) δ, ∀ x : M, (∃ α, Covers E J F (fun x => κ * gaugePot E J Φ ψ t x) α x) ∧
        Φ t (ψ t x) = graphImm E J F (fun x => κ * gaugePot E J Φ ψ t x) x := by
  obtain ⟨δ, hδ, hδT, ψ, hψ, hψ0, hloc⟩ := exists_global_gaugeS hJi hJJ hF hdim hT hΦ hΦ0
  have hΦδ := hΦ.mono hδT
  -- the chart statement at every point
  have hch : ∀ x : M, ∃ (α : M) (B : Set E) (ξ : ℝ × E → (E →L[ℝ] ℝ)), IsOpen B ∧
      B ⊆ (chart E α).target ∧ x ∈ (chart E α).source ∧ chart E α x ∈ B ∧
      ContDiffOn ℝ ∞ (fun q : ℝ × E => κ * gaugePot E J Φ ψ q.1 ((chart E α).symm q.2))
        (Icc 0 δ ×ˢ B) ∧
      ∀ t ∈ Icc (0 : ℝ) δ, ∀ y ∈ B,
        ((y, ξ (t, y)) ∈ chartDom E J F α ∧
          Φ t (ψ t ((chart E α).symm y)) = chartPsi E J F α (y, ξ (t, y))) ∧
        HasFDerivAt (fun z => κ * gaugePot E J Φ ψ t ((chart E α).symm z)) (ξ (t, y)) y ∧
        HasDerivWithinAt (fun s => κ * gaugePot E J Φ ψ s ((chart E α).symm y))
          (-(κ * flowPotAt E J Φ t (ψ t ((chart E α).symm y)))) (Icc 0 δ) t := by
    intro x
    obtain ⟨α, B, ξ, hBo, hBU, hxs, hxB, hξ, hξ0, hG⟩ := hloc x
    have hG' : ∀ t ∈ Icc (0 : ℝ) δ, ∀ y ∈ B, (y, ξ (t, y)) ∈ chartDom E J F α ∧
        Φ t (ψ t ((chart E α).symm y)) = chartPsi E J F α (y, ξ (t, y)) :=
      fun t ht y hy => (hG t ht y hy).2
    obtain ⟨hsm, hder⟩ := gaugePot_chartS hδ hΦδ hψ (isOpen_chartDom hJi hJJ hF α)
      (contDiffOn_chartPsi hJi hJJ hF α) (chartPsi_symplectic hJi hJJ hF α) hBo hBU hξ hG' hJJ hξ0
    exact ⟨α, B, ξ, hBo, hBU, hxs, hxB, hsm, fun t ht y hy => ⟨hG' t ht y hy, hder t ht y hy⟩⟩
  -- joint smoothness
  have hfsm : ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, ℝ) ∞
      (fun q : ℝ × M => κ * gaugePot E J Φ ψ q.1 q.2) (Icc 0 δ ×ˢ univ) := by
    rintro ⟨t, x⟩ hq
    obtain ⟨α, B, ξ, hBo, hBU, hxs, hxB, hsm, -⟩ := hch x
    have hDo : IsOpen ((chart E α).source ∩ chart E α ⁻¹' B) :=
      (continuousOn_extChartAt (I := 𝓘(ℝ, E)) α).isOpen_inter_preimage
        (isOpen_extChartAt_source (I := 𝓘(ℝ, E)) α) hBo
    have c1 : ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, ℝ × E) ∞
        (fun q : ℝ × M => ((q.1, chart E α q.2) : ℝ × E))
        (Icc 0 δ ×ˢ ((chart E α).source ∩ chart E α ⁻¹' B)) := by
      have ca : ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, ℝ) ∞ (fun q : ℝ × M => q.1)
          (Icc 0 δ ×ˢ ((chart E α).source ∩ chart E α ⁻¹' B)) := contMDiffOn_fst
      have cb : ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, E) ∞ (fun q : ℝ × M => chart E α q.2)
          (Icc 0 δ ×ˢ ((chart E α).source ∩ chart E α ⁻¹' B)) :=
        ContMDiffOn.comp (g := chart E α) (f := fun q : ℝ × M => q.2)
          (t := (chart E α).source)
          (by rw [extChartAt_source]; exact contMDiffOn_extChartAt) contMDiffOn_snd
          (fun q hq => hq.2.1)
      exact ca.prodMk_space cb
    have c3 := ContMDiffOn.comp
      (g := fun q : ℝ × E => κ * gaugePot E J Φ ψ q.1 ((chart E α).symm q.2))
      (f := fun q : ℝ × M => ((q.1, chart E α q.2) : ℝ × E))
      (contMDiffOn_iff_contDiffOn.mpr hsm) c1 (fun q hq => ⟨hq.1, hq.2.2⟩)
    have c5 : ContMDiffWithinAt (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, ℝ) ∞
        (fun q : ℝ × M => κ * gaugePot E J Φ ψ q.1 q.2)
        (Icc 0 δ ×ˢ ((chart E α).source ∩ chart E α ⁻¹' B)) (t, x) := by
      refine (c3 (t, x) ⟨hq.1, hxs, hxB⟩).congr (fun q hq' => ?_) ?_
      · show κ * gaugePot E J Φ ψ q.1 q.2 = κ * gaugePot E J Φ ψ q.1 ((chart E α).symm (chart E α q.2))
        rw [(chart E α).left_inv hq'.2.1]
      · show κ * gaugePot E J Φ ψ t x = κ * gaugePot E J Φ ψ t ((chart E α).symm (chart E α x))
        rw [(chart E α).left_inv hxs]
    exact c5.mono_of_mem_nhdsWithin (mem_nhdsWithin.2 ⟨univ ×ˢ ((chart E α).source ∩ chart E α ⁻¹' B),
      isOpen_univ.prod hDo, ⟨mem_univ _, hxs, hxB⟩, fun p hp => ⟨hp.2.1, hp.1.2⟩⟩)
  refine ⟨δ, hδ, hδT, ψ, hψ, hψ0, hfsm, fun t ht x => ?_, fun t ht x => ?_⟩
  · obtain ⟨α, B, ξ, hBo, hBU, hxs, hxB, hsm, hder⟩ := hch x
    have h := (hder t ht _ hxB).2.2
    rw [(chart E α).left_inv hxs] at h
    exact h
  · obtain ⟨α, B, ξ, hBo, hBU, hxs, hxB, hsm, hder⟩ := hch x
    -- the slice `f_t` is smooth on `M`
    have hslice : ContMDiff 𝓘(ℝ, E) 𝓘(ℝ, ℝ) ∞ (fun x => κ * gaugePot E J Φ ψ t x) := by
      intro x'
      have hc : ContMDiffAt 𝓘(ℝ, E) (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) ∞ (fun y : M => ((t, y) : ℝ × M)) x' :=
        contMDiffAt_const.prodMk contMDiffAt_id
      have h := ContMDiffWithinAt.comp (g := fun q : ℝ × M => κ * gaugePot E J Φ ψ q.1 q.2)
        (f := fun y : M => ((t, y) : ℝ × M)) (s := univ) x' (hfsm (t, x') ⟨ht, mem_univ _⟩)
        hc.contMDiffWithinAt (fun y _ => ⟨ht, mem_univ _⟩)
      exact contMDiffWithinAt_univ.mp h
    have hf : ∀ β : M, ContDiffOn ℝ ∞ (chartRep E (fun x => κ * gaugePot E J Φ ψ t x) β) (chart E β).target :=
      fun β => contMDiffOn_iff_contDiffOn.mp
        (hslice.comp_contMDiffOn (contMDiffOn_extChartAt_symm (I := 𝓘(ℝ, E)) β))
    have hd := hder t ht _ hxB
    have hfd : fderiv ℝ (chartRep E (fun x => κ * gaugePot E J Φ ψ t x) α) (chart E α x) =
        ξ (t, chart E α x) := hd.2.1.fderiv
    have hcov : Covers E J F (fun x => κ * gaugePot E J Φ ψ t x) α x := by
      refine ⟨hxs, hBU hxB, ?_⟩
      show (chart E α x, fderiv ℝ (chartRep E (fun x => κ * gaugePot E J Φ ψ t x) α) (chart E α x)) ∈
        chartDom E J F α
      rw [hfd]
      exact hd.1.1
    refine ⟨⟨α, hcov⟩, ?_⟩
    rw [graphImm_eq hf hJi hJJ hF hcov]
    show Φ t (ψ t x) = chartPsi E J F α
      (chart E α x, fderiv ℝ (chartRep E (fun x => κ * gaugePot E J Φ ψ t x) α) (chart E α x))
    rw [hfd, ← hd.1.2, (chart E α).left_inv hxs]

/-- **A flow from `F` is, for a short time and up to reparametrization, the graph of an exact section.** -/
theorem exists_gauge_graph (hJi : ∀ u v, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v, J (J v) = -v)
    (hF : IsLagImmersion E J F) (hdim : Module.finrank ℝ V = 2 * Module.finrank ℝ E)
    {T : ℝ} (hT : 0 < T) {Φ : ℝ → M → V} (hΦ : IsFlow E J T Φ) (hΦ0 : Φ 0 = F) :
    ∃ δ, 0 < δ ∧ δ ≤ T ∧ ∃ ψ : ℝ → M → M,
      ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, E) ∞ (fun q : ℝ × M => ψ q.1 q.2) (Icc 0 δ ×ˢ univ) ∧
      (∀ x, ψ 0 x = x) ∧
      ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, ℝ) ∞
        (fun q : ℝ × M => gaugePot E J Φ ψ q.1 q.2) (Icc 0 δ ×ˢ univ) ∧
      (∀ t ∈ Icc (0 : ℝ) δ, ∀ x : M, HasDerivWithinAt (fun s => gaugePot E J Φ ψ s x)
        (-flowPotAt E J Φ t (ψ t x)) (Icc 0 δ) t) ∧
      ∀ t ∈ Icc (0 : ℝ) δ, ∀ x : M, (∃ α, Covers E J F (gaugePot E J Φ ψ t) α x) ∧
        Φ t (ψ t x) = graphImm E J F (gaugePot E J Φ ψ t) x := by
  obtain ⟨δ, hδ, hδT, ψ, hψ, hψ0, h1, h2, h3⟩ :=
    exists_gauge_graphS (κ := 1) hJi hJJ hF hdim hT hΦ.toS hΦ0
  simp only [one_mul] at h1 h2 h3
  exact ⟨δ, hδ, hδT, ψ, hψ, hψ0, h1, h2, h3⟩

end Bridge

end CWFlow

/-- info: 'CWFlow.Bridge.exists_gauge_graphS' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms CWFlow.Bridge.exists_gauge_graphS

/-- info: 'CWFlow.Bridge.exists_gauge_graph' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms CWFlow.Bridge.exists_gauge_graph
