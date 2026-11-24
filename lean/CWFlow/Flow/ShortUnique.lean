import CWFlow.Flow.GaugeGraph
import CWFlow.Bridge.ScalarSolution

/-!
# Short-time uniqueness of the flow

`flow_unique_short`: two flows (`IsFlow`) on `[0, T]` with the same initial immersion `F` agree on `[0, δ]` for
some `δ > 0`.

Proof: by `exists_gauge_graph` each flow is, up to a reparametrization `ψ_t`, the graph of `d(gaugePot t)`.
For a short time the potential is a solution of the Mantegazza–Martinazzi equation
(`isSolution_gaugePot`: the jets stay in the domain by the tube lemma, and `∂ₜf = −div(JH)` of the graph by
invariance of `div(JH)` under reparametrization). The two potentials have the same initial datum `0`, so
`flow_unique_of_graph` applies.
-/

open Set Filter Bundle
open DifferentialGeometry DifferentialGeometry.Tensor0SBundle ShortTime ShortTime.HFlow CWFlow.JetRec
open DifferentialGeometry.Geometry.Curvature DifferentialGeometry.Geometry.Connection
open DifferentialGeometry.Integral.L2
open DifferentialGeometry.Analysis.Sobolev DifferentialGeometry.Analysis.Spectral
open DifferentialGeometry.Analysis.Parabolic.TensorSpectral DifferentialGeometry.Analysis.Elliptic
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

/-- A positive lower bound for finitely many positive numbers. -/
theorem exists_pos_le_of_fintype {ι : Type*} [Fintype ι] {δ : ℝ} (hδ : 0 < δ) {ε : ι → ℝ}
    (hε : ∀ i, 0 < ε i) : ∃ δ', 0 < δ' ∧ δ' ≤ δ ∧ ∀ i, δ' ≤ ε i := by
  classical
  have h : ∀ s : Finset ι, ∃ δ', 0 < δ' ∧ δ' ≤ δ ∧ ∀ i ∈ s, δ' ≤ ε i := by
    intro s
    induction s using Finset.induction_on with
    | empty => exact ⟨δ, hδ, le_rfl, by simp⟩
    | insert a s ha ih =>
      obtain ⟨d, h1, h2, h3⟩ := ih
      refine ⟨min d (ε a), lt_min h1 (hε a), (min_le_left _ _).trans h2, ?_⟩
      intro k hk
      rcases Finset.mem_insert.1 hk with hka | hk
      · rw [hka]
        exact min_le_right _ _
      · exact (min_le_left _ _).trans (h3 k hk)
  obtain ⟨d, h1, h2, h3⟩ := h Finset.univ
  exact ⟨d, h1, h2, fun i => h3 i (Finset.mem_univ i)⟩

variable (g : SmoothRiemannianMetric 𝓘(ℝ, E) M)

/-- `scalarFam g f t` only depends on `f t`. -/
theorem scalarFam_congr {f f' : ℝ → M → ℝ} {t : ℝ} (h : f t = f' t) :
    scalarFam g f t = scalarFam g f' t := by
  have key : ∀ (a b : M → ℝ) (_ : a = b) (ha : ContMDiff 𝓘(ℝ, E) 𝓘(ℝ, ℝ) ∞ a)
      (hb : ContMDiff 𝓘(ℝ, E) 𝓘(ℝ, ℝ) ∞ b), scalarTensor g a ha = scalarTensor g b hb := by
    rintro a _ rfl ha hb
    rfl
  by_cases hs : ContMDiff 𝓘(ℝ, E) 𝓘(ℝ, ℝ) ∞ (f t)
  · have hs' : ContMDiff 𝓘(ℝ, E) 𝓘(ℝ, ℝ) ∞ (f' t) := h ▸ hs
    rw [scalarFam_eq g hs, scalarFam_eq g hs']
    exact key _ _ h hs hs'
  · have hs' : ¬ContMDiff 𝓘(ℝ, E) 𝓘(ℝ, ℝ) ∞ (f' t) := fun h' => hs (h.symm ▸ h')
    unfold scalarFam
    rw [dif_neg hs, dif_neg hs']

variable {J : V →L[ℝ] V} {F : M → V}

/-- **Invariance of `div(JH)` under a reparametrization of `M`**, in charts. -/
theorem flowPot_chartRep_comp {Φt : M → V} {ψt : M → M} (hΦt : IsLagImmersion E J Φt)
    (hψt : ContMDiff 𝓘(ℝ, E) 𝓘(ℝ, E) ∞ ψt) (hG : IsLagImmersion E J (Φt ∘ ψt)) {α x : M}
    (hx : x ∈ (chart E α).source) :
    flowPot J (chartRep E (Φt ∘ ψt) α) (chart E α x) =
      flowPot J (chartRep E Φt (ψt x)) (chart E (ψt x) (ψt x)) := by
  obtain ⟨p, hp⟩ : ∃ p, p = ψt x := ⟨_, rfl⟩
  have hpsrc : p ∈ (chart E p).source := mem_extChartAt_source (I := 𝓘(ℝ, E)) p
  have hβ : ψt x ∈ (chart E p).source := by rw [← hp]; exact hpsrc
  have hNo : IsOpen (diffeoDom E ψt α p) := isOpen_diffeoDom' hψt.continuous α p
  have hτ : ContDiffOn ℝ ∞ (diffeoRep E ψt α p) (diffeoDom E ψt α p) :=
    contDiffOn_diffeoRep' hψt α p
  have hev : ∀ z ∈ diffeoDom E ψt α p,
      chartRep E (Φt ∘ ψt) α =ᶠ[𝓝 z] chartRep E Φt p ∘ diffeoRep E ψt α p := fun z hz => by
    filter_upwards [hNo.mem_nhds hz] with w hw using chartRep_comp_diffeo Φt hw
  have hXc : ∀ z ∈ diffeoDom E ψt α p,
      ContDiffAt ℝ ∞ (chartRep E Φt p) (diffeoRep E ψt α p z) := fun z hz =>
    (hΦt.smooth p).contDiffAt ((isOpen_chart_target (E := E) p).mem_nhds
      ((chart E p).map_source hz.2))
  have hτinv : ∀ z ∈ diffeoDom E ψt α p,
      ∃ e : E ≃L[ℝ] E, (e : E →L[ℝ] E) = fderiv ℝ (diffeoRep E ψt α p) z := by
    intro z hz
    have hρ : DifferentiableAt ℝ (diffeoRep E ψt α p) z :=
      (hτ.contDiffAt (hNo.mem_nhds hz)).differentiableAt (by simp)
    have hcomp := fderiv_comp z ((hXc z hz).differentiableAt (by simp)) hρ
    have hinj : Function.Injective (fderiv ℝ (diffeoRep E ψt α p) z) := by
      intro b b' hab
      apply hG.inj α z hz.1
      rw [(hev z hz).fderiv_eq, hcomp, ContinuousLinearMap.comp_apply,
        ContinuousLinearMap.comp_apply, hab]
    exact ⟨(LinearMap.linearEquivOfInjective (fderiv ℝ (diffeoRep E ψt α p) z : E →ₗ[ℝ] E) hinj
      rfl).toContinuousLinearEquiv, ContinuousLinearMap.ext fun _ => rfl⟩
  have hyN : chart E α x ∈ diffeoDom E ψt α p := mem_diffeoDom hx hβ
  have h := flowPot_comp hNo hτ hτinv hXc
    (fun z hz => hΦt.inj p _ ((chart E p).map_source hz.2)) J hyN
  rw [diffeoRep_chart hx] at h
  rw [flowPot_congr_of_eventuallyEq (hev _ hyN), h, ← hp]

variable (hJi : ∀ u v, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v, J (J v) = -v) (hF : IsLagImmersion E J F)
include hJi hJJ hF

/-- **The graphs stay in the Weinstein domains for a short time.** -/
theorem exists_time_graphDom {f : ℝ → M → ℝ} {δ : ℝ} (hδ : 0 < δ)
    (hf : ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, ℝ) ∞ (fun q : ℝ × M => f q.1 q.2)
      (Icc 0 δ ×ˢ univ)) (hf0 : f 0 = 0) :
    ∃ δ', 0 < δ' ∧ δ' ≤ δ ∧ ∀ t ∈ Icc (0 : ℝ) δ', ∀ i : PIdx E M,
      ∀ x ∈ tsupport (patchρ E M i), chart E i.1 x ∈ graphDom E J F (f t) i.1 := by
  have hW := weinstein_global (E := E) (M := M) hJi hJJ hF
  have hstep : ∀ i : PIdx E M, ∃ ε, 0 < ε ∧ ε ≤ δ ∧ ∀ t ∈ Icc (0 : ℝ) ε,
      ∀ x ∈ tsupport (patchρ E M i),
        ((chart E i.1 x, fderiv ℝ (fun z => f t ((chart E i.1).symm z)) (chart E i.1 x)) :
          E × (E →L[ℝ] ℝ)) ∈ chartDom E J F i.1 := by
    intro i
    have hsrc : ∀ x ∈ tsupport (patchρ E M i), x ∈ (chart E i.1).source := fun x hx =>
      (tsupport_patchρ_subset i hx).1
    have hw : ContDiffOn ℝ ∞ (fun q : ℝ × E => f q.1 ((chart E i.1).symm q.2))
        (Icc 0 δ ×ˢ (chart E i.1).target) :=
      contMDiffOn_iff_contDiffOn.mp (ContMDiffOn.comp (g := fun q : ℝ × M => f q.1 q.2) hf
        (contMDiffOn_chartPre (E := E) i.1 δ) (fun p hp => ⟨hp.1, mem_univ _⟩))
    have hD : ContinuousOn (fun p : ℝ × E =>
        fderiv ℝ (fun z => f p.1 ((chart E i.1).symm z)) p.2) (Icc 0 δ ×ˢ (chart E i.1).target) :=
      (contDiffOn_slice_fderiv hδ (isOpen_chart_target (E := E) i.1) hw).continuousOn
    have hκ : ContinuousOn (fun q : ℝ × M => ((q.1, chart E i.1 q.2) : ℝ × E))
        (Icc 0 δ ×ˢ tsupport (patchρ E M i)) :=
      continuous_fst.continuousOn.prodMk (ContinuousOn.comp (g := chart E i.1)
        (f := fun q : ℝ × M => q.2) (continuousOn_extChartAt (I := 𝓘(ℝ, E)) i.1)
        continuous_snd.continuousOn (fun q hq => hsrc q.2 hq.2))
    have hmap : MapsTo (fun q : ℝ × M => ((q.1, chart E i.1 q.2) : ℝ × E))
        (Icc 0 δ ×ˢ tsupport (patchρ E M i)) (Icc 0 δ ×ˢ (chart E i.1).target) :=
      fun q hq => ⟨hq.1, (chart E i.1).map_source (hsrc q.2 hq.2)⟩
    have hG : ContinuousOn (fun q : ℝ × M =>
        ((chart E i.1 q.2, fderiv ℝ (fun z => f q.1 ((chart E i.1).symm z)) (chart E i.1 q.2)) :
          E × (E →L[ℝ] ℝ))) (Icc 0 δ ×ˢ tsupport (patchρ E M i)) :=
      (continuous_snd.comp_continuousOn hκ).prodMk
        (ContinuousOn.comp (g := fun p : ℝ × E =>
          fderiv ℝ (fun z => f p.1 ((chart E i.1).symm z)) p.2)
          (f := fun q : ℝ × M => ((q.1, chart E i.1 q.2) : ℝ × E)) hD hκ hmap)
    refine exists_pos_forall_mem_of_isCompact (isClosed_tsupport _).isCompact hδ hG
      (isOpen_chartDom hJi hJJ hF i.1) fun x hx => ?_
    have hc : (fun z => f 0 ((chart E i.1).symm z)) = fun _ : E => (0 : ℝ) := by
      funext z
      rw [hf0]
      rfl
    have h0 : fderiv ℝ (fun z => f 0 ((chart E i.1).symm z)) (chart E i.1 x) = 0 := by
      rw [hc]
      exact fderiv_const_apply (0 : ℝ)
    show ((chart E i.1 x, fderiv ℝ (fun z => f 0 ((chart E i.1).symm z)) (chart E i.1 x)) :
      E × (E →L[ℝ] ℝ)) ∈ chartDom E J F i.1
    rw [h0]
    exact (hW.1 i.1).2 _ ((chart E i.1).map_source (hsrc x hx))
  choose ε hε hεδ hεmem using hstep
  obtain ⟨δ', hδ', hδ'δ, hδ'ε⟩ := exists_pos_le_of_fintype hδ hε
  exact ⟨δ', hδ', hδ'δ, fun t ht i x hx =>
    ⟨(chart E i.1).map_source (tsupport_patchρ_subset i hx).1,
      hεmem i t ⟨ht.1, ht.2.trans (hδ'ε i)⟩ x hx⟩⟩

/-- **The potential of a gauged flow is a Mantegazza–Martinazzi solution.** -/
theorem isSolution_gaugePot {δ : ℝ} (hδ : 0 < δ) (hδ1 : δ ≤ 1) {Φ : ℝ → M → V}
    (hΦ : IsFlow E J δ Φ) {ψ : ℝ → M → M}
    (hψ : ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, E) ∞ (fun q : ℝ × M => ψ q.1 q.2)
      (Icc 0 δ ×ˢ univ))
    (hf : ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, ℝ) ∞
      (fun q : ℝ × M => gaugePot E J Φ ψ q.1 q.2) (Icc 0 δ ×ˢ univ))
    (hd : ∀ t ∈ Icc (0 : ℝ) δ, ∀ x : M, HasDerivWithinAt (fun s => gaugePot E J Φ ψ s x)
      (-flowPotAt E J Φ t (ψ t x)) (Icc 0 δ) t)
    (hg : ∀ t ∈ Icc (0 : ℝ) δ, ∀ x : M, (∃ α, Covers E J F (gaugePot E J Φ ψ t) α x) ∧
      Φ t (ψ t x) = graphImm E J F (gaugePot E J Φ ψ t) x)
    (hdom : ∀ t ∈ Icc (0 : ℝ) δ, ∀ i : PIdx E M, ∀ x ∈ tsupport (patchρ E M i),
      chart E i.1 x ∈ graphDom E J F (gaugePot E J Φ ψ t) i.1) :
    (cwData g hJi hJJ hF).toEquationData.IsSolution g (scalarFam g (gaugePot E J Φ ψ) 0) δ
      (scalarFam g (gaugePot E J Φ ψ)) := by
  refine isSolution_scalarFam g hJi hJJ hF hδ hδ1 hf hdom fun t ht x α hα hgα => ?_
  have hsl := contMDiff_slice_fn hf ht
  have hfsl : ∀ β : M, ContDiffOn ℝ ∞ (chartRep E (gaugePot E J Φ ψ t) β) (chart E β).target :=
    fun β => contMDiffOn_iff_contDiffOn.mp
      (hsl.comp_contMDiffOn (contMDiffOn_extChartAt_symm (I := 𝓘(ℝ, E)) β))
  have hGeq : graphImm E J F (gaugePot E J Φ ψ t) = Φ t ∘ ψ t :=
    funext fun x' => ((hg t ht x').2).symm
  have hGimm : IsLagImmersion E J (Φ t ∘ ψ t) := by
    rw [← hGeq]
    exact isLagImmersion_graphImm hfsl hJi hJJ hF (fun x' => (hg t ht x').1)
  have hev : chartRep E (Φ t ∘ ψ t) α =ᶠ[𝓝 (chart E α x)]
      chartGraph E J F (gaugePot E J Φ ψ t) α := by
    filter_upwards [(isOpen_graphDom' hfsl hJi hJJ hF α).mem_nhds hgα] with w hw
    rw [← hGeq]
    show graphImm E J F (gaugePot E J Φ ψ t) ((chart E α).symm w) = _
    have hcov : Covers E J F (gaugePot E J Φ ψ t) α ((chart E α).symm w) := by
      refine ⟨(chart E α).map_target hw.1, ?_⟩
      rw [(chart E α).right_inv hw.1]
      exact hw
    rw [graphImm_eq hfsl hJi hJJ hF hcov, (chart E α).right_inv hw.1]
  have key : flowPot J (chartGraph E J F (gaugePot E J Φ ψ t) α) (chart E α x) =
      flowPotAt E J Φ t (ψ t x) := by
    rw [← flowPot_congr_of_eventuallyEq (J := J) hev]
    exact flowPot_chartRep_comp (hΦ.imm t ht) (contMDiff_slice hψ ht) hGimm hα.1
  rw [key]
  exact hd t ht x

include g in
/-- **Short-time uniqueness of the flow.** Two flows with the same initial immersion agree for a short
time. -/
theorem flow_unique_short (hdim : Module.finrank ℝ V = 2 * Module.finrank ℝ E) {T : ℝ} (hT : 0 < T)
    {Φ₁ Φ₂ : ℝ → M → V} (hΦ₁ : IsFlow E J T Φ₁) (hΦ₂ : IsFlow E J T Φ₂)
    (h₁ : Φ₁ 0 = F) (h₂ : Φ₂ 0 = F) :
    ∃ δ, 0 < δ ∧ δ ≤ T ∧ ∀ t ∈ Icc (0 : ℝ) δ, Φ₁ t = Φ₂ t := by
  obtain ⟨δ₁, hδ₁, hδ₁T, ψ₁, hψ₁, hψ₁0, hf₁, hd₁, hg₁⟩ :=
    exists_gauge_graph hJi hJJ hF hdim hT hΦ₁ h₁
  obtain ⟨δ₂, hδ₂, hδ₂T, ψ₂, hψ₂, hψ₂0, hf₂, hd₂, hg₂⟩ :=
    exists_gauge_graph hJi hJJ hF hdim hT hΦ₂ h₂
  -- a common time `δ₀ ≤ 1`
  obtain ⟨δ₀, hδ₀, hδ₀1, hδ₀2, hδ₀3⟩ : ∃ δ₀, 0 < δ₀ ∧ δ₀ ≤ δ₁ ∧ δ₀ ≤ δ₂ ∧ δ₀ ≤ 1 :=
    ⟨min (min δ₁ δ₂) 1, lt_min (lt_min hδ₁ hδ₂) one_pos, (min_le_left _ _).trans (min_le_left _ _),
      (min_le_left _ _).trans (min_le_right _ _), min_le_right _ _⟩
  have hsub₁ : Icc (0 : ℝ) δ₀ ×ˢ (univ : Set M) ⊆ Icc (0 : ℝ) δ₁ ×ˢ univ :=
    prod_mono (Icc_subset_Icc_right hδ₀1) subset_rfl
  have hsub₂ : Icc (0 : ℝ) δ₀ ×ˢ (univ : Set M) ⊆ Icc (0 : ℝ) δ₂ ×ˢ univ :=
    prod_mono (Icc_subset_Icc_right hδ₀2) subset_rfl
  obtain ⟨δa, hδa, hδa0, hdom₁⟩ := exists_time_graphDom hJi hJJ hF hδ₀ (hf₁.mono hsub₁)
    (funext fun x => gaugePot_zero J Φ₁ ψ₁ x)
  obtain ⟨δb, hδb, hδb0, hdom₂⟩ := exists_time_graphDom hJi hJJ hF hδ₀ (hf₂.mono hsub₂)
    (funext fun x => gaugePot_zero J Φ₂ ψ₂ x)
  obtain ⟨δ, hδ, hδa', hδb'⟩ : ∃ δ, 0 < δ ∧ δ ≤ δa ∧ δ ≤ δb :=
    ⟨min δa δb, lt_min hδa hδb, min_le_left _ _, min_le_right _ _⟩
  have hδδ₀ : δ ≤ δ₀ := hδa'.trans hδa0
  have hI₁ : Icc (0 : ℝ) δ ⊆ Icc (0 : ℝ) δ₁ := Icc_subset_Icc_right (hδδ₀.trans hδ₀1)
  have hI₂ : Icc (0 : ℝ) δ ⊆ Icc (0 : ℝ) δ₂ := Icc_subset_Icc_right (hδδ₀.trans hδ₀2)
  have hS₁ : Icc (0 : ℝ) δ ×ˢ (univ : Set M) ⊆ Icc (0 : ℝ) δ₁ ×ˢ univ := prod_mono hI₁ subset_rfl
  have hS₂ : Icc (0 : ℝ) δ ×ˢ (univ : Set M) ⊆ Icc (0 : ℝ) δ₂ ×ˢ univ := prod_mono hI₂ subset_rfl
  have hΦ₁δ : IsFlow E J δ Φ₁ := hΦ₁.mono ((hδδ₀.trans hδ₀1).trans hδ₁T)
  have hΦ₂δ : IsFlow E J δ Φ₂ := hΦ₂.mono ((hδδ₀.trans hδ₀2).trans hδ₂T)
  have hIcc : Icc (0 : ℝ) δ ⊆ Icc (0 : ℝ) δ := subset_rfl
  have hu₁ := isSolution_gaugePot g hJi hJJ hF hδ (hδδ₀.trans hδ₀3) hΦ₁δ (hψ₁.mono hS₁)
    (hf₁.mono hS₁) (fun t ht x => (hd₁ t (hI₁ ht) x).mono hI₁) (fun t ht x => hg₁ t (hI₁ ht) x)
    (fun t ht i x hx => hdom₁ t ⟨ht.1, ht.2.trans hδa'⟩ i x hx)
  have hu₂ := isSolution_gaugePot g hJi hJJ hF hδ (hδδ₀.trans hδ₀3) hΦ₂δ (hψ₂.mono hS₂)
    (hf₂.mono hS₂) (fun t ht x => (hd₂ t (hI₂ ht) x).mono hI₂) (fun t ht x => hg₂ t (hI₂ ht) x)
    (fun t ht i x hx => hdom₂ t ⟨ht.1, ht.2.trans hδb'⟩ i x hx)
  have hu0 : scalarFam g (gaugePot E J Φ₁ ψ₁) 0 = scalarFam g (gaugePot E J Φ₂ ψ₂) 0 :=
    scalarFam_congr g (funext fun x => by rw [gaugePot_zero, gaugePot_zero])
  rw [hu0] at hu₁
  have hG₁ : ∀ t ∈ Icc (0 : ℝ) δ, ∀ x : M, Φ₁ t (ψ₁ t x) =
      graphImm E J F (fnOf g (scalarFam g (gaugePot E J Φ₁ ψ₁) t)) x := fun t ht x => by
    rw [fnOf_scalarFam g (contMDiff_slice_fn (hf₁.mono hS₁) ht)]
    exact (hg₁ t (hI₁ ht) x).2
  have hG₂ : ∀ t ∈ Icc (0 : ℝ) δ, ∀ x : M, Φ₂ t (ψ₂ t x) =
      graphImm E J F (fnOf g (scalarFam g (gaugePot E J Φ₂ ψ₂) t)) x := fun t ht x => by
    rw [fnOf_scalarFam g (contMDiff_slice_fn (hf₂.mono hS₂) ht)]
    exact (hg₂ t (hI₂ ht) x).2
  exact ⟨δ, hδ, (hδδ₀.trans hδ₀1).trans hδ₁T,
    flow_unique_of_graph g hJi hJJ hF hdim hu₁ hu₂ hΦ₁δ hΦ₂δ (hψ₁.mono hS₁) (hψ₂.mono hS₂)
      hψ₁0 hψ₂0 hG₁ hG₂⟩

end Bridge

end CWFlow
