import CWFlow.Extension.Slope
import CWFlow.Extension.IntM

/-!
# Local structure of a Lagrangian immersion of a manifold over a fixed Lagrangian plane

For a Lagrangian immersion `F : M → V` read in charts: the derivative `dFq` of `F` at `q` in the chart at
`q`, its behavior under change of chart (`dFq_eq_comp`), the open set `goodSet` where the horizontal projection
of the tangent plane is injective and `F` is in a box (`isOpen_goodSet`), the horizontal projection as a local
homeomorphism there, closedness of the slope condition (`isClosed_slopeLE_dFq`) and local graph
representations on `M` (`exists_local_graph_M`).
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

/-- The derivative of `F` at `q`, read in the chart at `q`. -/
noncomputable def dFq (F : M → V) (q : M) : EuclideanSpace ℝ (Fin n) →L[ℝ] V :=
  fderiv ℝ (F ∘ chartInv ι q) (chartFwd ι q q)

/-- The points where `F` lies over the ball of radius `ρ`, below height `1`, with horizontal projection of the
tangent plane injective. -/
def goodSet (J : V →L[ℝ] V) (U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V) (y₀ : V) (F : M → V) (ρ : ℝ) :
    Set M :=
  {q | ‖horC U y₀ (F q)‖ < ρ ∧ ‖verC J U y₀ (F q)‖ < 1 ∧
    Function.Injective (horCL U ∘L dFq ι F q)}

/-! ### Part 1: change of chart -/

theorem dFq_eq_comp (hF : IsLagImm ι J F) {q z : M}
    (hz : q ∈ (extChartAt 𝓘(ℝ, E) z).source) :
    ∃ T : EuclideanSpace ℝ (Fin n) ≃L[ℝ] EuclideanSpace ℝ (Fin n),
      dFq ι F q = fderiv ℝ (F ∘ chartInv ι z) (chartFwd ι z q) ∘L (T : _ →L[ℝ] _) := by
  have hqq : q ∈ (extChartAt 𝓘(ℝ, E) q).source := mem_extChartAt_source q
  have hw0D : chartFwd ι q q ∈ chartDomE ι q := im_chartFwd_mem ι hqq
  have hback : chartInv ι q (chartFwd ι q q) = q := im_chartInv_chartFwd ι hqq
  have hWo : IsOpen (chartDomE ι q ∩ chartInv ι q ⁻¹' (extChartAt 𝓘(ℝ, E) z).source) :=
    (im_continuousOn_chartInv ι q).isOpen_inter_preimage (im_isOpen_chartDomE ι q)
      (isOpen_extChartAt_source z)
  have hw0W : chartFwd ι q q ∈
      chartDomE ι q ∩ chartInv ι q ⁻¹' (extChartAt 𝓘(ℝ, E) z).source :=
    ⟨hw0D, by show chartInv ι q (chartFwd ι q q) ∈ (extChartAt 𝓘(ℝ, E) z).source; rw [hback]; exact hz⟩
  have hτ : DifferentiableAt ℝ (chartFwd ι z ∘ chartInv ι q) (chartFwd ι q q) :=
    ((im_contDiffOn_trans ι q z).contDiffAt (hWo.mem_nhds hw0W)).differentiableAt (by simp)
  have hτw : (chartFwd ι z ∘ chartInv ι q) (chartFwd ι q q) = chartFwd ι z q := by
    simp only [Function.comp_apply, hback]
  have hτD : chartFwd ι z q ∈ chartDomE ι z := im_chartFwd_mem ι hz
  have hFd : DifferentiableAt ℝ (F ∘ chartInv ι z)
      ((chartFwd ι z ∘ chartInv ι q) (chartFwd ι q q)) := by
    rw [hτw]
    exact ((hF.smooth z).contDiffAt ((im_isOpen_chartDomE ι z).mem_nhds hτD)).differentiableAt
      (by simp)
  have hev : (F ∘ chartInv ι q) =ᶠ[𝓝 (chartFwd ι q q)] (F ∘ chartInv ι z) ∘
      (chartFwd ι z ∘ chartInv ι q) := by
    filter_upwards [hWo.mem_nhds hw0W] with w hw
    simp only [Function.comp_apply]
    rw [im_chartInv_chartFwd ι hw.2]
  have hcomp : dFq ι F q = (fderiv ℝ (F ∘ chartInv ι z) (chartFwd ι z q)).comp
      (fderiv ℝ (chartFwd ι z ∘ chartInv ι q) (chartFwd ι q q)) := by
    unfold dFq
    rw [hev.fderiv_eq, fderiv_comp _ hFd hτ, hτw]
  -- left inverse of the transition
  have hτ'D : chartFwd ι z q ∈ chartDomE ι z ∩ chartInv ι z ⁻¹' (extChartAt 𝓘(ℝ, E) q).source :=
    ⟨hτD, by show chartInv ι z (chartFwd ι z q) ∈ (extChartAt 𝓘(ℝ, E) q).source; rw [im_chartInv_chartFwd ι hz]; exact hqq⟩
  have hW'o : IsOpen (chartDomE ι z ∩ chartInv ι z ⁻¹' (extChartAt 𝓘(ℝ, E) q).source) :=
    (im_continuousOn_chartInv ι z).isOpen_inter_preimage (im_isOpen_chartDomE ι z)
      (isOpen_extChartAt_source q)
  have hτ' : DifferentiableAt ℝ (chartFwd ι q ∘ chartInv ι z)
      ((chartFwd ι z ∘ chartInv ι q) (chartFwd ι q q)) := by
    rw [hτw]
    exact ((im_contDiffOn_trans ι z q).contDiffAt (hW'o.mem_nhds hτ'D)).differentiableAt (by simp)
  have hlinv : (chartFwd ι q ∘ chartInv ι z) ∘ (chartFwd ι z ∘ chartInv ι q) =ᶠ[𝓝 (chartFwd ι q q)]
      id := by
    filter_upwards [hWo.mem_nhds hw0W] with w hw
    simp only [Function.comp_apply, id]
    rw [im_chartInv_chartFwd ι hw.2, im_chartFwd_chartInv ι hw.1]
  have h1 := fderiv_comp (chartFwd ι q q) hτ' hτ
  have h2 : fderiv ℝ ((chartFwd ι q ∘ chartInv ι z) ∘ (chartFwd ι z ∘ chartInv ι q))
      (chartFwd ι q q) = ContinuousLinearMap.id ℝ _ := by
    rw [hlinv.fderiv_eq]; exact fderiv_id
  have key : ∀ u, fderiv ℝ (chartFwd ι q ∘ chartInv ι z)
      ((chartFwd ι z ∘ chartInv ι q) (chartFwd ι q q))
      (fderiv ℝ (chartFwd ι z ∘ chartInv ι q) (chartFwd ι q q) u) = u := fun u => by
    have := congrArg (fun L : EuclideanSpace ℝ (Fin n) →L[ℝ] EuclideanSpace ℝ (Fin n) => L u)
      (h1.symm.trans h2)
    simpa using this
  have hinj : Function.Injective (fderiv ℝ (chartFwd ι z ∘ chartInv ι q) (chartFwd ι q q)) := by
    intro a b hab
    rw [← key a, ← key b, hab]
  obtain ⟨T, hT⟩ := gl_exists_equiv hinj
  exact ⟨T, by rw [hcomp, hT]⟩

theorem injective_horCL_dFq_iff (hF : IsLagImm ι J F) {q z : M}
    (hz : q ∈ (extChartAt 𝓘(ℝ, E) z).source) :
    Function.Injective (horCL U ∘L dFq ι F q) ↔
      Function.Injective (horCL U ∘L fderiv ℝ (F ∘ chartInv ι z) (chartFwd ι z q)) := by
  obtain ⟨T, hT⟩ := dFq_eq_comp ι hF hz
  have hfun : ⇑(horCL U ∘L dFq ι F q) =
      ⇑(horCL U ∘L fderiv ℝ (F ∘ chartInv ι z) (chartFwd ι z q)) ∘ T := by
    rw [hT]; rfl
  rw [hfun]
  constructor
  · intro h a b hab
    obtain ⟨a', rfl⟩ := T.surjective a
    obtain ⟨b', rfl⟩ := T.surjective b
    rw [h hab]
  · intro h
    exact h.comp T.injective

theorem slopeZ_dFq_eq (hF : IsLagImm ι J F) {q z : M}
    (hz : q ∈ (extChartAt 𝓘(ℝ, E) z).source)
    (hinj : Function.Injective (horCL U ∘L dFq ι F q)) :
    slopeZ J U (dFq ι F q) = slopeZ J U (fderiv ℝ (F ∘ chartInv ι z) (chartFwd ι z q)) := by
  obtain ⟨T, hT⟩ := dFq_eq_comp ι hF hz
  rw [hT]
  exact slopeZ_comp T ((injective_horCL_dFq_iff ι hF hz).1 hinj)

theorem slopeLE_dFq_iff (hF : IsLagImm ι J F) {q z : M}
    (hz : q ∈ (extChartAt 𝓘(ℝ, E) z).source) (c : ℝ) :
    slopeLE J U c (dFq ι F q) ↔
      slopeLE J U c (fderiv ℝ (F ∘ chartInv ι z) (chartFwd ι z q)) := by
  obtain ⟨T, hT⟩ := dFq_eq_comp ι hF hz
  rw [hT]
  constructor
  · intro h a
    obtain ⟨b, rfl⟩ := T.surjective a
    exact h b
  · intro h a
    exact h (T a)

/-! ### Part 2: continuity, openness, local homeomorphism -/

theorem grl_continuousOn_F (hF : IsLagImm ι J F) (z : M) :
    ContinuousOn F (extChartAt 𝓘(ℝ, E) z).source :=
  (((hF.smooth z).continuousOn.comp (im_continuousOn_chartFwd ι z)
    (fun p hp => im_chartFwd_mem ι hp))).congr fun p hp => by
      simp only [Function.comp_apply, im_chartInv_chartFwd ι hp]

theorem grl_continuous_F (hF : IsLagImm ι J F) : Continuous F :=
  continuous_iff_continuousAt.2 fun p => (grl_continuousOn_F ι hF p).continuousAt
    ((isOpen_extChartAt_source p).mem_nhds (mem_extChartAt_source p))

theorem grl_continuous_horC : Continuous (horC U y₀) :=
  (horCL U).continuous.comp (continuous_id.sub continuous_const)

theorem grl_continuous_verC : Continuous (verC J U y₀) :=
  (verCL J U).continuous.comp (continuous_id.sub continuous_const)

theorem continuousOn_fderiv_chart (hF : IsLagImm ι J F) (z : M) :
    ContinuousOn (fun q => fderiv ℝ (F ∘ chartInv ι z) (chartFwd ι z q))
      (extChartAt 𝓘(ℝ, E) z).source :=
  ((ims_contDiffOn_fderiv (im_isOpen_chartDomE ι z) (hF.smooth z)).continuousOn).comp
    (im_continuousOn_chartFwd ι z) (fun p hp => im_chartFwd_mem ι hp)

theorem grl_isOpen_inj :
    IsOpen {L : EuclideanSpace ℝ (Fin n) →L[ℝ] V | Function.Injective (horCL U ∘L L)} := by
  have h : {L : EuclideanSpace ℝ (Fin n) →L[ℝ] V | Function.Injective (horCL U ∘L L)} =
      (fun L : EuclideanSpace ℝ (Fin n) →L[ℝ] V => horCL U ∘L L) ⁻¹'
        Set.range ((↑) : (EuclideanSpace ℝ (Fin n) ≃L[ℝ] EuclideanSpace ℝ (Fin n)) →
          EuclideanSpace ℝ (Fin n) →L[ℝ] EuclideanSpace ℝ (Fin n)) := by
    ext L
    simp only [Set.mem_ofPred_eq, Set.mem_preimage, Set.mem_range]
    constructor
    · intro hL
      obtain ⟨E', hE'⟩ := gl_exists_equiv hL
      exact ⟨E', hE'⟩
    · rintro ⟨E', hE'⟩
      rw [← hE']
      exact E'.injective
  rw [h]
  exact ContinuousLinearEquiv.isOpen.preimage
    (ContinuousLinearMap.compL ℝ (EuclideanSpace ℝ (Fin n)) V (EuclideanSpace ℝ (Fin n))
      (horCL U)).continuous

theorem isOpen_goodSet (hF : IsLagImm ι J F) (ρ : ℝ) : IsOpen (goodSet ι J U y₀ F ρ) := by
  rw [isOpen_iff_mem_nhds]
  intro q hq
  obtain ⟨h1, h2, h3⟩ := hq
  have hFc := grl_continuous_F ι hF
  have e1 : ∀ᶠ q' in 𝓝 q, ‖horC U y₀ (F q')‖ < ρ :=
    (isOpen_lt (continuous_norm.comp (grl_continuous_horC.comp hFc)) continuous_const).mem_nhds h1
  have e2 : ∀ᶠ q' in 𝓝 q, ‖verC J U y₀ (F q')‖ < 1 :=
    (isOpen_lt (continuous_norm.comp (grl_continuous_verC.comp hFc)) continuous_const).mem_nhds h2
  have hsrc : (extChartAt 𝓘(ℝ, E) q).source ∈ 𝓝 q :=
    (isOpen_extChartAt_source q).mem_nhds (mem_extChartAt_source q)
  have hg : ContinuousAt (fun p => fderiv ℝ (F ∘ chartInv ι q) (chartFwd ι q p)) q :=
    (continuousOn_fderiv_chart ι hF q).continuousAt hsrc
  have h3' : Function.Injective
      (horCL U ∘L fderiv ℝ (F ∘ chartInv ι q) (chartFwd ι q q)) :=
    (injective_horCL_dFq_iff ι hF (mem_extChartAt_source q)).1 h3
  have e3 := hg.eventually_mem (grl_isOpen_inj.mem_nhds h3')
  filter_upwards [e1, e2, e3, hsrc] with q' a b c d
  exact ⟨a, b, (injective_horCL_dFq_iff ι hF d).2 c⟩

/-! ### Part 4: local graphs on `M` (placed before the local homeomorphism, which uses it) -/

theorem exists_local_graph_M (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    (hU : IsLagFrame J U) [FiniteDimensional ℝ V] (hdim : Module.finrank ℝ V = 2 * n)
    (hF : IsLagImm ι J F) {ρ : ℝ} {q : M} (hq : q ∈ goodSet ι J U y₀ F ρ) :
    ∃ r > 0, ∃ (s : EuclideanSpace ℝ (Fin n) → EuclideanSpace ℝ (Fin n))
      (φ : EuclideanSpace ℝ (Fin n) → ℝ),
      ContDiff ℝ ∞ φ ∧ ContDiffOn ℝ ∞ s (ball (horC U y₀ (F q)) r) ∧
      Set.MapsTo s (ball (horC U y₀ (F q)) r) (chartDomE ι q) ∧
      (∀ x ∈ ball (horC U y₀ (F q)) r,
        horC U y₀ (F (chartInv ι q (s x))) = x ∧ F (chartInv ι q (s x)) = graphX J U y₀ φ x ∧
        Function.Injective (horCL U ∘L dFq ι F (chartInv ι q (s x))) ∧
        slopeZ J U (dFq ι F (chartInv ι q (s x))) = jet 2 φ x ∧
        aNormSq J (F ∘ chartInv ι q) (s x) = gNormSq (ginvFun (jet 2 φ x)) (jet 3 φ x)) ∧
      ∀ᶠ q' in 𝓝 q, q' ∈ (extChartAt 𝓘(ℝ, E) q).source ∧
        horC U y₀ (F q') ∈ ball (horC U y₀ (F q)) r ∧
        chartInv ι q (s (horC U y₀ (F q'))) = q' := by
  have hqq : q ∈ (extChartAt 𝓘(ℝ, E) q).source := mem_extChartAt_source q
  have hw0D : chartFwd ι q q ∈ chartDomE ι q := im_chartFwd_mem ι hqq
  have hback : chartInv ι q (chartFwd ι q q) = q := im_chartInv_chartFwd ι hqq
  obtain ⟨r, hr, s, φ, hφ, hs, hmaps, hsy, hx, hev⟩ :=
    exists_local_graph hJi hJJ hU hdim y₀ (im_isOpen_chartDomE ι q) (hF.smooth q) (hF.inj q)
      (hF.lag q) hw0D hq.2.2
  have hq0 : horC U y₀ ((F ∘ chartInv ι q) (chartFwd ι q q)) = horC U y₀ (F q) := by
    rw [Function.comp_apply, hback]
  rw [hq0] at hs hmaps hx hev
  have hXd : ∀ x ∈ ball (horC U y₀ (F q)) r, DifferentiableAt ℝ (F ∘ chartInv ι q) (s x) :=
    fun x hxb => ((hF.smooth q).contDiffAt ((im_isOpen_chartDomE ι q).mem_nhds (hmaps hxb))).differentiableAt
      (by simp)
  refine ⟨r, hr, s, φ, hφ, hs, hmaps, ?_, ?_⟩
  · intro x hxb
    obtain ⟨h1, h2, h3, h4⟩ := hx x hxb
    have hsx := hmaps hxb
    have hp : chartInv ι q (s x) ∈ (extChartAt 𝓘(ℝ, E) q).source := im_chartInv_mem_source ι hsx
    have hfw : chartFwd ι q (chartInv ι q (s x)) = s x := im_chartFwd_chartInv ι hsx
    have hsd : DifferentiableAt ℝ s x :=
      (hs.contDiffAt (isOpen_ball.mem_nhds hxb)).differentiableAt (by simp)
    have hev' : ((F ∘ chartInv ι q) ∘ s) =ᶠ[𝓝 x] graphX J U y₀ φ :=
      Filter.eventuallyEq_of_mem (isOpen_ball.mem_nhds hxb) fun w hw => (hx w hw).2.1
    have hDG : fderiv ℝ (graphX J U y₀ φ) x =
        (fderiv ℝ (F ∘ chartInv ι q) (s x)).comp (fderiv ℝ s x) := by
      rw [← hev'.fderiv_eq]
      exact fderiv_comp x (hXd x hxb) hsd
    have hid := horCL_comp_fderiv_graphX hJi hJJ hU hφ y₀ x
    have hsurj : Function.Surjective (horCL U ∘L fderiv ℝ (F ∘ chartInv ι q) (s x)) := by
      intro u
      refine ⟨fderiv ℝ s x u, ?_⟩
      have := congrArg (fun L : EuclideanSpace ℝ (Fin n) →L[ℝ] EuclideanSpace ℝ (Fin n) => L u) hid
      simpa [hDG] using this
    have hinj' : Function.Injective (horCL U ∘L fderiv ℝ (F ∘ chartInv ι q) (s x)) :=
      (LinearMap.injective_iff_surjective
        (f := ((horCL U ∘L fderiv ℝ (F ∘ chartInv ι q) (s x) :
          EuclideanSpace ℝ (Fin n) →L[ℝ] EuclideanSpace ℝ (Fin n)) :
          EuclideanSpace ℝ (Fin n) →ₗ[ℝ] EuclideanSpace ℝ (Fin n))).2 hsurj)
    have hI : Function.Injective (horCL U ∘L dFq ι F (chartInv ι q (s x))) :=
      (injective_horCL_dFq_iff ι hF hp).2 (by rw [hfw]; exact hinj')
    refine ⟨h1, h2, hI, ?_, h4⟩
    rw [slopeZ_dFq_eq ι hF hp hI, hfw]
    exact h3
  · have hc : ContinuousAt (chartFwd ι q) q :=
      (im_continuousOn_chartFwd ι q).continuousAt
        ((isOpen_extChartAt_source q).mem_nhds hqq)
    have hsrc : (extChartAt 𝓘(ℝ, E) q).source ∈ 𝓝 q :=
      (isOpen_extChartAt_source q).mem_nhds hqq
    filter_upwards [hc.eventually hev, hsrc] with q' h h'
    obtain ⟨-, h1, h2⟩ := h
    have e := im_chartInv_chartFwd ι h'
    simp only [Function.comp_apply, e] at h1 h2
    refine ⟨h', h1, ?_⟩
    rw [h2, e]

theorem isLocalHomeomorphOn_horC_goodSet (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫)
    (hJJ : ∀ v : V, J (J v) = -v) (hU : IsLagFrame J U) [FiniteDimensional ℝ V]
    (hdim : Module.finrank ℝ V = 2 * n) (hF : IsLagImm ι J F) (ρ : ℝ) :
    IsLocalHomeomorphOn (fun q => horC U y₀ (F q)) (goodSet ι J U y₀ F ρ) := by
  apply IsLocalHomeomorphOn.mk
  intro q hq
  obtain ⟨r, hr, s, φ, hφ, hs, hmaps, hx, hev⟩ := exists_local_graph_M ι hJi hJJ hU hdim hF hq
  set c := horC U y₀ (F q) with hc
  set N : Set M := {q' | q' ∈ (extChartAt 𝓘(ℝ, E) q).source ∧ horC U y₀ (F q') ∈ ball c r ∧
    chartInv ι q (s (horC U y₀ (F q'))) = q'} with hNdef
  have hN : N ∈ 𝓝 q := hev
  have hqS : q ∈ interior N := mem_interior_iff_mem_nhds.2 hN
  have hcs : ContinuousOn (fun x => chartInv ι q (s x)) (ball c r) :=
    (im_continuousOn_chartInv ι q).comp hs.continuousOn hmaps
  have hT : IsOpen (ball c r ∩ (fun x => chartInv ι q (s x)) ⁻¹' interior N) :=
    hcs.isOpen_inter_preimage isOpen_ball isOpen_interior
  have hFc := grl_continuous_F ι hF
  refine ⟨{ toFun := fun q => horC U y₀ (F q)
            invFun := fun x => chartInv ι q (s x)
            source := interior N
            target := ball c r ∩ (fun x => chartInv ι q (s x)) ⁻¹' interior N
            map_source' := ?_
            map_target' := ?_
            left_inv' := ?_
            right_inv' := ?_
            open_source := isOpen_interior
            open_target := hT
            continuousOn_toFun := (grl_continuous_horC.comp hFc).continuousOn
            continuousOn_invFun := hcs.mono inter_subset_left }, hqS, fun _ _ => rfl⟩
  · intro q' hq'
    have h := interior_subset hq'
    exact ⟨h.2.1, by show chartInv ι q (s (horC U y₀ (F q'))) ∈ interior N; rw [h.2.2]; exact hq'⟩
  · intro x hxt
    exact hxt.2
  · intro q' hq'
    exact (interior_subset hq').2.2
  · intro x hxt
    exact (hx x hxt.1).1

/-! ### Part 3: closedness of the slope condition -/

theorem isClosed_slopeLE_dFq (hF : IsLagImm ι J F) (c : ℝ) :
    IsClosed {q : M | slopeLE J U c (dFq ι F q)} := by
  rw [← isOpen_compl_iff, isOpen_iff_mem_nhds]
  intro q hq
  have hqq : q ∈ (extChartAt 𝓘(ℝ, E) q).source := mem_extChartAt_source q
  have hsrc : (extChartAt 𝓘(ℝ, E) q).source ∈ 𝓝 q :=
    (isOpen_extChartAt_source q).mem_nhds hqq
  have hg : ContinuousAt (fun p => fderiv ℝ (F ∘ chartInv ι q) (chartFwd ι q p)) q :=
    (continuousOn_fderiv_chart ι hF q).continuousAt hsrc
  have hq' : fderiv ℝ (F ∘ chartInv ι q) (chartFwd ι q q) ∈
      {L : EuclideanSpace ℝ (Fin n) →L[ℝ] V | slopeLE J U c L}ᶜ := fun h =>
    hq ((slopeLE_dFq_iff ι hF hqq c).2 h)
  have e := hg.eventually_mem ((isClosed_slopeLE (J := J) (U := U) c).isOpen_compl.mem_nhds hq')
  filter_upwards [e, hsrc] with q' a b h
  exact a ((slopeLE_dFq_iff ι hF b c).1 h)

/-- info: 'CWFlow.Extension.dFq_eq_comp' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms dFq_eq_comp

/-- info: 'CWFlow.Extension.isOpen_goodSet' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms isOpen_goodSet

/-- info: 'CWFlow.Extension.isLocalHomeomorphOn_horC_goodSet' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms isLocalHomeomorphOn_horC_goodSet

/-- info: 'CWFlow.Extension.isClosed_slopeLE_dFq' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms isClosed_slopeLE_dFq

/-- info: 'CWFlow.Extension.exists_local_graph_M' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_local_graph_M

end CWFlow.Extension
