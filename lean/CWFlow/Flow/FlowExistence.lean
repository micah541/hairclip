import CWFlow.Flow.TanField
import CWFlow.Flow.Reparam
import CWFlow.Flow.DiffeoChart
import CWFlow.Scalar.Smooth
import DifferentialGeometry.Analysis.ODE.TimeDependentFlow.BoundaryExtension.SeeleyTimeExtension

/-!
# Short-time existence for the flow `∂ₜF = J∇div(JH)` of Lagrangian immersions

Let `u` solve the scalar equation, `G_t = graphImm (fnOf g (u t))`.

* `dtW_add_fderiv_eRep`: the normal part of the velocity of `G` is the flow velocity,
  `∂ₜG − DG(tan) = J∇div(JH)[G_t]`, in a good chart;
* `hasDerivWithinAt_reparam`: along a curve `χ s x` with velocity `tanField` (minus the tangential part),
  `d/ds G_s(χ_s x) = J∇div(JH)[G_t]` at `χ_t x`;
* `exists_flow`: **short-time existence**: there are `T > 0` and a family `Φ t` of Lagrangian immersions with
  `Φ 0 = F` and `∂ₜΦ = J∇div(JH)[Φ_t]` on `[0, T]` (Chen–Warren, Claim 3.3 and the existence part of the theorem).
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

namespace CWFlow.Bridge

universe u

variable {E : Type u} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
  [NeZero (Module.finrank ℝ E)] [CompleteSpace E]
variable {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
  [CompactSpace M] [T2Space M] [SigmaCompactSpace M]
variable {V : Type u} [NormedAddCommGroup V] [InnerProductSpace ℝ V] [FiniteDimensional ℝ V]

variable (g : SmoothRiemannianMetric 𝓘(ℝ, E) M) {J : V →L[ℝ] V} {F : M → V}
variable (hJi : ∀ u v, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v, J (J v) = -v) (hF : IsLagImmersion E J F)
variable (hdim : Module.finrank ℝ V = 2 * Module.finrank ℝ E)

omit [NeZero (Module.finrank ℝ E)] [CompleteSpace E] [CompactSpace M] [SigmaCompactSpace M]
  [FiniteDimensional ℝ V] in
theorem tanPart_neg' (L : E →L[ℝ] V) (w : V) : tanPart L (-w) = -tanPart L w := by
  have h0 : tanPart L 0 = 0 := by simpa using tanPart_sub L w w
  have h := tanPart_sub L 0 w
  rwa [zero_sub, h0, zero_sub] at h

variable {u₀ : SmoothCcTensor g 0 0} {T : ℝ} {u : ℝ → SmoothCcTensor g 0 0}
variable (hu : (cwData g hJi hJJ hF).toEquationData.IsSolution g u₀ T u)
include hu hdim

/-- **The normal part of the velocity of the family is the flow velocity** `J∇div(JH)`, at a good point. -/
theorem dtW_add_fderiv_eRep {t : ℝ} (ht : t ∈ Icc (0 : ℝ) T) {α x : M} (hα : x ∈ patchU E α)
    (hgα : good g J F α x (covJet g 0 0 (2 * 2) (u t) x)) :
    dtW (famRep g J F u α) T (t, chart E α x) +
        fderiv ℝ (chartRep E (graphImm E J F (fnOf g (u t))) α) (chart E α x)
          (eRep g J F u T α (t, chart E α x)) =
      flowVel J (chartRep E (graphImm E J F (fnOf g (u t))) α) (chart E α x) := by
  have hf := contDiffOn_chartRep_fnOf g (u t)
  have hy := mem_graphDom_of_good g J F (u t) hα hgα
  have hfun : (fun s => famRep g J F u α (s, chart E α x)) =
      fun s => graphImm E J F (fnOf g (u s)) x := by
    funext s
    simp only [famRep, chartRep, Function.comp_apply]
    rw [(chart E α).left_inv hα.1]
  have hw : dtW (famRep g J F u α) T (t, chart E α x) =
      graphVel J F (fnOf g (u t)) α (chart E α x) := by
    show derivWithin (fun s => famRep g J F u α (s, chart E α x)) (Icc 0 T) t = _
    rw [hfun]
    exact (hasDerivWithinAt_graphImm g hJi hJJ hF hu ht hα hgα).derivWithin
      (uniqueDiffOn_Icc hu.pos t ht)
  have heq : chartRep E (graphImm E J F (fnOf g (u t))) α =ᶠ[𝓝 (chart E α x)]
      chartGraph E J F (fnOf g (u t)) α := by
    filter_upwards [(isOpen_graphDom' hf hJi hJJ hF α).mem_nhds hy] with w hw
    have hc : Covers E J F (fnOf g (u t)) α ((chart E α).symm w) :=
      ⟨(chart E α).map_target hw.1, by rw [(chart E α).right_inv hw.1]; exact hw⟩
    have h := graphImm_eq hf hJi hJJ hF hc
    rw [(chart E α).right_inv hw.1] at h
    exact h
  have hgv := graphVel_eq hJi hJJ hF hdim hf hy
  have hinj := injective_fderiv_chartGraph hf hJi hJJ hF hy
  have hlag := lagrangian_fderiv_chartGraph hf hJi hJJ hF hy
  have h7 := sub_apply_tanPart_of_eq hinj hlag
    (b := ginvOf (fderiv ℝ (chartGraph E J F (fnOf g (u t)) α) (chart E α x))
      (fderiv ℝ (flowPot J (chartGraph E J F (fnOf g (u t)) α)) (chart E α x))) hgv
  show dtW (famRep g J F u α) T (t, chart E α x) +
      fderiv ℝ (chartRep E (graphImm E J F (fnOf g (u t))) α) (chart E α x)
        (tanPart (fderiv ℝ (chartRep E (graphImm E J F (fnOf g (u t))) α) (chart E α x))
          (-dtW (famRep g J F u α) T (t, chart E α x))) = _
  rw [flowVel_congr_of_eventuallyEq heq, heq.fderiv_eq, hw, tanPart_neg', map_neg,
    ← sub_eq_add_neg, h7]
  rfl

/-- **Along the reparametrization, the family moves by the flow velocity.** -/
theorem hasDerivWithinAt_reparam {χ : ℝ → M → M} {t : ℝ} (ht : t ∈ Icc (0 : ℝ) T) {x : M}
    (hχ : HasMFDerivAt 𝓘(ℝ, ℝ) 𝓘(ℝ, E) (fun s => χ s x) t
      ((1 : ℝ →L[ℝ] ℝ).smulRight (tanField g J F u T t (χ t x))))
    {α : M} (hα : χ t x ∈ patchU E α)
    (hgα : good g J F α (χ t x) (covJet g 0 0 (2 * 2) (u t) (χ t x))) :
    HasDerivWithinAt (fun s => graphImm E J F (fnOf g (u s)) (χ s x))
      (flowVel J (chartRep E (graphImm E J F (fnOf g (u t))) α) (chart E α (χ t x))) (Icc 0 T) t := by
  have hs : χ t x ∈ (chartAt E α).source := by
    have h := hα.1
    rwa [extChartAt_source] at h
  have hch := (mdifferentiableAt_extChartAt (I := 𝓘(ℝ, E)) hs).hasMFDerivAt
  have hcomp := hch.comp t hχ
  have hfd : HasFDerivAt (fun s => chart E α (χ s x)) (_ : ℝ →L[ℝ] E) t :=
    hasMFDerivAt_iff_hasFDerivAt.mp hcomp
  have hread := trivAt_tanField g hJi hJJ hF hu ht hα.1
  rw [TangentBundle.continuousLinearMapAt_trivializationAt hs] at hread
  have hd : HasDerivAt (fun s => chart E α (χ s x)) (eRep g J F u T α (t, chart E α (χ t x))) t := by
    refine hfd.hasDerivAt.congr_deriv ?_
    rw [← hread]
    exact congrArg _ (one_smul ℝ _)
  have hc : HasDerivAt (fun s => (s, chart E α (χ s x)))
      ((1 : ℝ), eRep g J F u T α (t, chart E α (χ t x))) t := (hasDerivAt_id t).prodMk hd
  have hN : {s | χ s x ∈ (chart E α).source} ∈ 𝓝 t :=
    hχ.continuousAt.preimage_mem_nhds ((isOpen_extChartAt_source (I := 𝓘(ℝ, E)) α).mem_nhds hα.1)
  have hG := hasFDerivWithinAt_prod hu.pos (isOpen_chart_target (E := E) α)
    (contDiffOn_famRep g hJi hJJ hF hu α) ht ((chart E α).map_source hα.1)
  have hm : MapsTo (fun s => (s, chart E α (χ s x))) (Icc 0 T ∩ {s | χ s x ∈ (chart E α).source})
      (Icc 0 T ×ˢ (chart E α).target) := fun s hs => ⟨hs.1, (chart E α).map_source hs.2⟩
  have h1 := HasFDerivWithinAt.comp_hasDerivWithinAt (f := fun s => (s, chart E α (χ s x))) t hG
    (hc.hasDerivWithinAt (s := Icc 0 T ∩ {s | χ s x ∈ (chart E α).source})) hm
  have h2 := h1.mono_of_mem_nhdsWithin (inter_mem self_mem_nhdsWithin (nhdsWithin_le_nhds hN))
  have hev : (fun s => graphImm E J F (fnOf g (u s)) (χ s x)) =ᶠ[𝓝[Icc 0 T] t]
      famRep g J F u α ∘ fun s => (s, chart E α (χ s x)) := by
    filter_upwards [nhdsWithin_le_nhds hN] with s hs
    simp only [famRep, chartRep, Function.comp_apply]
    rw [(chart E α).left_inv hs]
  refine (h2.congr_of_eventuallyEq hev (hev.self_of_nhdsWithin ht)).congr_deriv ?_
  rw [fderivWithin_prod_apply hu.pos (isOpen_chart_target (E := E) α)
    (contDiffOn_famRep g hJi hJJ hF hu α) ht ((chart E α).map_source hα.1), one_smul]
  exact dtW_add_fderiv_eRep g hJi hJJ hF hdim hu ht hα hgα

omit hdim in
/-- **The reparametrized family is jointly smooth** on `Icc 0 T × M`. -/
theorem contMDiffOn_reparam {χ : ℝ → M → M}
    (hχ : ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, E) ∞ (fun q : ℝ × M => χ q.1 q.2)
      (Ioo (-1) (T + 1) ×ˢ univ)) :
    ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, V) ∞
      (fun q : ℝ × M => graphImm E J F (fnOf g (u q.1)) (χ q.1 q.2)) (Icc 0 T ×ˢ univ) := by
  rintro ⟨t, x⟩ ⟨ht, -⟩
  obtain ⟨α, hα⟩ : ∃ α, α = χ t x := ⟨_, rfl⟩
  have hsrc : χ t x ∈ (chart E α).source := by
    rw [hα]
    exact mem_extChartAt_source (I := 𝓘(ℝ, E)) _
  have hχt : ContMDiffWithinAt (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, E) ∞ (fun q : ℝ × M => χ q.1 q.2)
      (Icc 0 T ×ˢ univ) (t, x) :=
    (hχ.mono (prod_mono Icc_subset_Ioo_reparam subset_rfl)) (t, x) ⟨ht, mem_univ _⟩
  have hA : {q : ℝ × M | χ q.1 q.2 ∈ (chart E α).source} ∈ 𝓝[Icc 0 T ×ˢ univ] (t, x) :=
    hχt.continuousWithinAt.preimage_mem_nhdsWithin
      ((isOpen_extChartAt_source (I := 𝓘(ℝ, E)) α).mem_nhds hsrc)
  have hsrc' : χ t x ∈ (chartAt E α).source := by rwa [extChartAt_source] at hsrc
  have hmc : MapsTo (fun q : ℝ × M => χ q.1 q.2)
      ((Icc 0 T ×ˢ univ) ∩ {q : ℝ × M | χ q.1 q.2 ∈ (chart E α).source}) (chartAt E α).source := by
    intro q hq
    have h := hq.2
    rwa [mem_ofPred_eq, extChartAt_source] at h
  have hc : ContMDiffWithinAt (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, E) ∞
      (fun q : ℝ × M => chart E α (χ q.1 q.2))
      ((Icc 0 T ×ˢ univ) ∩ {q : ℝ × M | χ q.1 q.2 ∈ (chart E α).source}) (t, x) :=
    ContMDiffWithinAt.comp (g := chart E α) (f := fun q : ℝ × M => χ q.1 q.2) (t, x)
      (contMDiffOn_extChartAt (I := 𝓘(ℝ, E)) (n := ∞) (x := α) _ hsrc')
      (hχt.mono inter_subset_left) hmc
  have h1 : ContMDiffWithinAt (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, ℝ × E) ∞
      (fun q : ℝ × M => (q.1, chart E α (χ q.1 q.2)))
      ((Icc 0 T ×ˢ univ) ∩ {q : ℝ × M | χ q.1 q.2 ∈ (chart E α).source}) (t, x) :=
    contMDiffWithinAt_fst.prodMk_space hc
  have h2 : ContMDiffWithinAt 𝓘(ℝ, ℝ × E) 𝓘(ℝ, V) ∞ (famRep g J F u α)
      (Icc 0 T ×ˢ (chart E α).target) (t, chart E α (χ t x)) :=
    (contDiffOn_famRep g hJi hJJ hF hu α (t, chart E α (χ t x))
      ⟨ht, (chart E α).map_source hsrc⟩).contMDiffWithinAt
  have hm : MapsTo (fun q : ℝ × M => (q.1, chart E α (χ q.1 q.2)))
      ((Icc 0 T ×ˢ univ) ∩ {q : ℝ × M | χ q.1 q.2 ∈ (chart E α).source})
      (Icc 0 T ×ˢ (chart E α).target) := fun q hq => ⟨hq.1.1, (chart E α).map_source hq.2⟩
  have h3 := ContMDiffWithinAt.comp (g := famRep g J F u α)
    (f := fun q : ℝ × M => (q.1, chart E α (χ q.1 q.2))) (t, x) h2 h1 hm
  have hev : (fun q : ℝ × M => graphImm E J F (fnOf g (u q.1)) (χ q.1 q.2))
      =ᶠ[𝓝[(Icc 0 T ×ˢ univ) ∩ {q : ℝ × M | χ q.1 q.2 ∈ (chart E α).source}] (t, x)]
      famRep g J F u α ∘ fun q : ℝ × M => (q.1, chart E α (χ q.1 q.2)) := by
    filter_upwards [self_mem_nhdsWithin] with q hq
    simp only [famRep, chartRep, Function.comp_apply]
    rw [(chart E α).left_inv hq.2]
  exact (h3.congr_of_eventuallyEq hev (hev.self_of_nhdsWithin ⟨⟨ht, mem_univ _⟩, hsrc⟩)).mono_of_mem_nhdsWithin
    (inter_mem self_mem_nhdsWithin hA)

omit hu
include g hJi hJJ hF

/-- **Short-time existence for the flow `∂ₜF = J∇div(JH)`** of a Lagrangian immersion `F` of a closed manifold:
a family `Φ t` of Lagrangian immersions on `[0, T]`, `T > 0`, jointly smooth in `(t, x)`, with `Φ 0 = F`, whose velocity at every point, read in
any chart, is `J∇div(JH)` of `Φ t`. -/
theorem exists_flow :
    ∃ T > 0, ∃ Φ : ℝ → M → V, Φ 0 = F ∧
      ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, V) ∞ (fun q : ℝ × M => Φ q.1 q.2) (Icc 0 T ×ˢ univ) ∧
      (∀ t ∈ Icc (0 : ℝ) T, IsLagImmersion E J (Φ t)) ∧
      ∀ t ∈ Icc (0 : ℝ) T, ∀ x z : M, x ∈ (chart E z).source →
        HasDerivWithinAt (fun s => Φ s x) (flowVel J (chartRep E (Φ t) z) (chart E z x)) (Icc 0 T) t := by
  obtain ⟨T, hT, u, hu⟩ := exists_scalar_solution g hJi hJJ hF
  obtain ⟨Y, hYsm, hYeq⟩ := DifferentialGeometry.Analysis.ODE.seeley_time_extend (I := 𝓘(ℝ, E))
    (tanField g J F u T) T hT (contMDiffOn_tanField g hJi hJJ hF hu)
  obtain ⟨χ, h0, hχsm, -, hode⟩ := exists_reparam_flow Y hYsm hT.le
  refine ⟨T, hT, fun t x => graphImm E J F (fnOf g (u t)) (χ t x), ?_, ?_, ?_, ?_⟩
  · funext x
    show graphImm E J F (fnOf g (u 0)) (χ 0 x) = F x
    rw [h0, hu.initial, fnOf_zero, graphImm_zero hJi hJJ hF]
  · exact contMDiffOn_reparam g hJi hJJ hF hu (χ := fun t x => χ t x) hχsm
  · intro t ht
    exact isLagImmersion_comp_diffeo (isLagImmersion_solution g hJi hJJ hF hu ht) (χ t)
  · intro t ht x z hz
    obtain ⟨i, hi, hgi⟩ := exists_good g hJi hJJ hF hu ht (χ t x)
    have hχ := hode t (Icc_subset_Ioo_reparam ht) x
    rw [hYeq t ht] at hχ
    have h := hasDerivWithinAt_reparam g hJi hJJ hF hdim hu ht hχ hi hgi
    have hG := isLagImmersion_solution g hJi hJJ hF hu ht
    have hmem : chart E z x ∈ diffeoDom E (χ t) z i.1 := mem_diffeoDom hz hi.1
    have heq : chartRep E (fun x => graphImm E J F (fnOf g (u t)) (χ t x)) z =ᶠ[𝓝 (chart E z x)]
        chartRep E (graphImm E J F (fnOf g (u t))) i.1 ∘ diffeoRep E (χ t) z i.1 := by
      filter_upwards [(isOpen_diffeoDom (χ t) z i.1).mem_nhds hmem] with w hw
      exact chartRep_comp_diffeo (graphImm E J F (fnOf g (u t))) hw
    have hfv := flowVel_comp' (isOpen_chart_target (E := E) i.1) (hG.smooth i.1) (hG.inj i.1) J
      (isOpen_diffeoDom (χ t) z i.1) (contDiffOn_diffeoRep (χ t) z i.1)
      (fun w hw => exists_equiv_fderiv_diffeoRep (χ t) hw)
      (fun w hw => (chart E i.1).map_source hw.2) hmem
    rw [diffeoRep_chart hz] at hfv
    show HasDerivWithinAt _ (flowVel J (chartRep E (fun x => graphImm E J F (fnOf g (u t)) (χ t x)) z)
      (chart E z x)) _ _
    rw [flowVel_congr_of_eventuallyEq heq, hfv]
    exact h

/-- info: 'CWFlow.Bridge.exists_flow' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_flow

end CWFlow.Bridge
