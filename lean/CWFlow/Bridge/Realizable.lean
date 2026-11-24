import CWFlow.Bridge.Data
import CWFlow.Scalar.Symm4

/-!
# The Mantegazza–Martinazzi equation on jets of functions

For a smooth function `u` on `M` (a `(0,0)`-tensor field, `fnOf g u` its scalar function), the right side of the
MM equation `cwData` at the covariant jet of `u` is the operator of the flow:

`lam x ((cwData …).F t x (covJet u x)) = −div(JH)[F_u]`, computed in any chart (`lam_cwData_F`).

Steps, for the covariant jet `ζ = covJet g 0 0 4 u x`:
* `cJet_covJet`: the coordinate jet of `ζ` in the chart at `α` is the jet of `u ∘ φ_α⁻¹`
  (`jetOf_scalRd`, the exact jet recursion);
* `permTop_jetOf`: the slot permutation does nothing on jets of functions (symmetry of `D⁴`);
* `Fi_covJet`: so the chart operator is `−div(JH)` of the chart graph (`neg_flowPot_chartGraph`), which does not
  depend on the chart (`flowPot_chartGraph_transition`);
* `lam_frameSum_cometricAt_eq`: for every jet with symmetric second-order part, MM's contraction of the top
  component with the cometric of one chart is `Fi ζ − Fi (zeroTop ζ)` (`Pjet_add_top`, `lam_frameSum_cometricAt`);
* `cometricAt_chart_indep`, `mmEf_covJet`: the cometrics of the charts agree on jets of functions, so the glued
  cometric is the cometric of any one chart.
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

universe u

section Calculus

variable {E : Type u} [NormedAddCommGroup E] [NormedSpace ℝ E]

/-- The jet of a function depends only on its germ. -/
theorem jetOf_congr {w w' : E → ℝ} {y : E} (h : w =ᶠ[𝓝 y] w') : jetOf w y = jetOf w' y := by
  unfold jetOf
  rw [h.fderiv_eq, (h.iteratedFDeriv ℝ 2).self_of_nhds, (h.iteratedFDeriv ℝ 3).self_of_nhds,
    (h.iteratedFDeriv ℝ 4).self_of_nhds]

theorem jetOf_two_symm {w : E → ℝ} {y : E} (hw : ContDiffAt ℝ ∞ w y) (a b : E) :
    (jetOf w y).2.1 ![a, b] = (jetOf w y).2.1 ![b, a] :=
  IsSymmSndFDerivAt.iteratedFDeriv_cons
    (hf := (hw.of_le (WithTop.coe_le_coe.2 le_top) : ContDiffAt ℝ 2 w y).isSymmSndFDerivAt (by simp))

end Calculus

namespace Bridge

variable {E : Type u} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
  [NeZero (Module.finrank ℝ E)] [CompleteSpace E]
variable {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
  [CompactSpace M] [T2Space M] [SigmaCompactSpace M]
variable {V : Type u} [NormedAddCommGroup V] [InnerProductSpace ℝ V]

/-- The slot permutation does nothing on jets of functions. -/
theorem permTop_jetOf {w : E → ℝ} {y : E} (hw : ContDiffAt ℝ ∞ w y) :
    permTop (jetOf w y) = jetOf w y := by
  refine Prod.ext rfl (Prod.ext rfl (Prod.ext rfl ?_))
  refine ContinuousMultilinearMap.ext fun m => ?_
  rw [permTop_apply_top]
  exact iteratedFDeriv_four_perm hw m

theorem permTop_addTop (z : SJet E) (T : E [×4]→L[ℝ] ℝ) :
    permTop (addTop z T) = addTop (permTop z) (T.domDomCongr σ4) := by
  unfold permTop addTop
  refine Prod.ext rfl (Prod.ext rfl (Prod.ext rfl ?_))
  exact ContinuousMultilinearMap.ext fun m => rfl

variable (g : SmoothRiemannianMetric 𝓘(ℝ, E) M)

/-- The scalar function of a `(0,0)`-tensor field. -/
def fnOf (u : SmoothCcTensor g 0 0) : M → ℝ := fun x => lam (I := 𝓘(ℝ, E)) x (u.toSection x)

theorem scalRd_eq_chartRep (α : M) (u : SmoothCcTensor g 0 0) {y : E} (hy : y ∈ (chart E α).target) :
    scalRd g α u y = chartRep E (fnOf g u) α y :=
  scal_chartRd g α u hy

theorem contDiffOn_chartRep_fnOf (u : SmoothCcTensor g 0 0) (α : M) :
    ContDiffOn ℝ ∞ (chartRep E (fnOf g u) α) (chart E α).target :=
  (contDiffOn_scalRd g α u).congr fun _ hy => (scalRd_eq_chartRep g α u hy).symm

theorem contDiffAt_chartRep_fnOf (u : SmoothCcTensor g 0 0) {α x : M} (hx : x ∈ (chart E α).source) :
    ContDiffAt ℝ ∞ (chartRep E (fnOf g u) α) (chart E α x) :=
  (contDiffOn_chartRep_fnOf g u α).contDiffAt
    ((isOpen_chart_target (E := E) α).mem_nhds ((chart E α).map_source hx))

/-- **The coordinate jet of the covariant jet of `u` is the jet of `u` in the chart.** -/
theorem cJet_covJet (u : SmoothCcTensor g 0 0) {α x : M} (hx : x ∈ patchU E α) :
    cJet g α x (covJet g 0 0 (2 * 2) u x) = jetOf (chartRep E (fnOf g u) α) (chart E α x) := by
  have h1 : cJet g α x (covJet g 0 0 (2 * 2) u x) = coordJet g α (chart E α x)
      (fun k : Fin (4 + 1) => chartRd g α (iteratedCovGrad g 0 0 k u) (chart E α x)) :=
    congrArg (coordJet g α (chart E α x)) (trivJet_covJet_eq g α u hx.1)
  have h2 : scalRd g α u =ᶠ[𝓝 (chart E α x)] chartRep E (fnOf g u) α := by
    filter_upwards [(isOpen_chart_target (E := E) α).mem_nhds ((chart E α).map_source hx.1)] with y hy
    exact scalRd_eq_chartRep g α u hy
  rw [h1, ← jetOf_scalRd g α u hx.2, jetOf_congr h2]

variable (J : V →L[ℝ] V) (F : M → V)

/-- The chart graph of a good jet of a function lies in the Weinstein domain. -/
theorem mem_graphDom_of_good (u : SmoothCcTensor g 0 0) {α x : M} (hx : x ∈ patchU E α)
    (hgood : good g J F α x (covJet g 0 0 (2 * 2) u x)) :
    chart E α x ∈ graphDom E J F (fnOf g u) α := by
  refine ⟨(chart E α).map_source hx.1, ?_⟩
  unfold good at hgood
  rw [cJet_covJet g u hx] at hgood
  exact hgood

variable {g J F}
variable (hJi : ∀ u v, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v, J (J v) = -v) (hF : IsLagImmersion E J F)
include hJi hJJ hF

/-- **The chart operator on the jet of a function is `−div(JH)` of the chart graph.** -/
theorem Fi_covJet (u : SmoothCcTensor g 0 0) {α x : M} (hx : x ∈ patchU E α)
    (hgood : good g J F α x (covJet g 0 0 (2 * 2) u x)) :
    Fi g J F α x (covJet g 0 0 (2 * 2) u x) =
      -flowPot J (chartGraph E J F (fnOf g u) α) (chart E α x) := by
  unfold Fi
  rw [cJet_covJet g u hx, permTop_jetOf (contDiffAt_chartRep_fnOf g u hx.1)]
  exact (neg_flowPot_chartGraph (contDiffOn_chartRep_fnOf g u) hJi hJJ hF
    (mem_graphDom_of_good g J F u hx hgood)).symm

/-- **Chart independence of the operator.** -/
theorem flowPot_chart_indep (u : SmoothCcTensor g 0 0) {α β x : M} (hα : x ∈ patchU E α)
    (hβ : x ∈ patchU E β) (hgα : good g J F α x (covJet g 0 0 (2 * 2) u x))
    (hgβ : good g J F β x (covJet g 0 0 (2 * 2) u x)) :
    flowPot J (chartGraph E J F (fnOf g u) α) (chart E α x) =
      flowPot J (chartGraph E J F (fnOf g u) β) (chart E β x) := by
  have hy : chart E α x ∈ overlap E α β :=
    ⟨(chart E α).map_source hα.1, by
      show (chart E α).symm (chart E α x) ∈ (chart E β).source
      rw [(chart E α).left_inv hα.1]
      exact hβ.1⟩
  have hτ : transition E α β (chart E α x) = chart E β x := by
    simp only [transition, Function.comp_apply]
    rw [(chart E α).left_inv hα.1]
  have h := flowPot_chartGraph_transition (contDiffOn_chartRep_fnOf g u) hJi hJJ hF hy
    (mem_graphDom_of_good g J F u hα hgα)
    (by rw [hτ]; exact mem_graphDom_of_good g J F u hβ hgβ)
  rw [hτ] at h
  exact h

/-- The differential of the graph at the jet of a function is the differential of the chart graph. -/
theorem Li_covJet (u : SmoothCcTensor g 0 0) {α x : M} (hx : x ∈ patchU E α)
    (hgood : good g J F α x (covJet g 0 0 (2 * 2) u x)) :
    Li g J F α x (covJet g 0 0 (2 * 2) u x) =
      fderiv ℝ (chartGraph E J F (fnOf g u) α) (chart E α x) := by
  have hf := contDiffOn_chartRep_fnOf g u
  have hy := mem_graphDom_of_good g J F u hx hgood
  have h1 := iteratedFDeriv_graphFu_eq_TkOf (isOpen_chartDom hJi hJJ hF α)
    (contDiffOn_chartPsi hJi hJJ hF α) (isOpen_graphDom' hf hJi hJJ hF α)
    ((hf α).mono fun w hw => hw.1) (fun w hw => hw.2) hy (k := 1) (by norm_num)
  unfold Li LOf
  rw [cJet_covJet g u hx, ← h1]
  ext a
  rw [continuousMultilinearCurryFin1_apply, iteratedFDeriv_one_apply]
  rfl

omit hJi hJJ hF in
/-- The symmetrized cometric transforms correctly under a change of chart. -/
theorem cometricAt_chart_indep {α β x : M} (hα : x ∈ (chart E α).source)
    (hβ : x ∈ (chart E β).source) {Lα Lβ : E →L[ℝ] V} (hinj : Function.Injective Lβ)
    (hL : Lα = Lβ.comp (fderiv ℝ (transition E α β) (chart E α x))) :
    cometricAt α x Lα = cometricAt β x Lβ := by
  have hy : chart E α x ∈ overlap E α β :=
    ⟨(chart E α).map_source hα, by
      show (chart E α).symm (chart E α x) ∈ (chart E β).source
      rw [(chart E α).left_inv hα]
      exact hβ⟩
  obtain ⟨e, he⟩ := exists_equiv_fderiv_transition (E := E) hy
  have hα' : x ∈ (chartAt E α).source := by rwa [extChartAt_source] at hα
  have hβ' : x ∈ (chartAt E β).source := by rwa [extChartAt_source] at hβ
  -- the tangent coordinate change is the derivative of the transition
  have hsymm : ∀ v : E, (trivializationAt E (TangentSpace 𝓘(ℝ, E)) α).symmL ℝ x v =
      (trivializationAt E (TangentSpace 𝓘(ℝ, E)) β).symmL ℝ x (e v) := by
    intro v
    have hcore : (tangentBundleCore 𝓘(ℝ, E) M).coordChange (achart E α) (achart E β) x v = e v := by
      rw [tangentBundleCore_coordChange_achart, ModelWithCorners.range_eq_univ, fderivWithin_univ]
      show fderiv ℝ (transition E α β) (chart E α x) v = e v
      rw [← he]
      rfl
    rw [TangentBundle.symmL_trivializationAt_eq_core hα',
      TangentBundle.symmL_trivializationAt_eq_core hβ', ← hcore]
    exact ((tangentBundleCore 𝓘(ℝ, E) M).coordChange_comp (achart E α) (achart E β) (achart E x) x
      ⟨⟨hα', hβ'⟩, mem_chart_source E x⟩ v).symm
  have hcov : ∀ ξ : TangentSpace 𝓘(ℝ, E) x →L[ℝ] ℝ,
      chartCovec α x ξ = (chartCovec β x ξ).comp (e : E →L[ℝ] E) := by
    intro ξ
    ext v
    simp only [chartCovec, ContinuousLinearMap.comp_apply, ContinuousLinearEquiv.coe_coe]
    rw [hsymm v]
  have hLe : Lα = Lβ.comp (e : E →L[ℝ] E) := by rw [hL, he]
  have hcancel : ∀ η : E →L[ℝ] ℝ, (η.comp (e : E →L[ℝ] E)).comp (e.symm : E →L[ℝ] E) = η := by
    intro η
    ext v
    simp
  have hkey : ∀ ξ η : TangentSpace 𝓘(ℝ, E) x →L[ℝ] ℝ,
      chartCovec α x ξ (ginvOf Lα (chartCovec α x η)) =
        chartCovec β x ξ (ginvOf Lβ (chartCovec β x η)) := by
    intro ξ η
    rw [hcov ξ, hcov η, hLe, ContinuousLinearMap.comp_apply, ContinuousLinearEquiv.coe_coe,
      apply_ginvOf_comp hinj e, hcancel]
  refine LinearMap.ext fun ξ => LinearMap.ext fun η => ?_
  rw [cometricAt_apply, cometricAt_apply, hkey ξ η, hkey η ξ]

/-- **The cometrics of two charts agree on jets of functions.** -/
theorem cometricAt_covJet_indep (u : SmoothCcTensor g 0 0) {α β x : M} (hα : x ∈ patchU E α)
    (hβ : x ∈ patchU E β) (hgα : good g J F α x (covJet g 0 0 (2 * 2) u x))
    (hgβ : good g J F β x (covJet g 0 0 (2 * 2) u x)) :
    cometricAt α x (Li g J F α x (covJet g 0 0 (2 * 2) u x)) =
      cometricAt β x (Li g J F β x (covJet g 0 0 (2 * 2) u x)) := by
  have hf := contDiffOn_chartRep_fnOf g u
  have hyα := mem_graphDom_of_good g J F u hα hgα
  have hyβ := mem_graphDom_of_good g J F u hβ hgβ
  have hy : chart E α x ∈ overlap E α β :=
    ⟨(chart E α).map_source hα.1, by
      show (chart E α).symm (chart E α x) ∈ (chart E β).source
      rw [(chart E α).left_inv hα.1]
      exact hβ.1⟩
  have hτ : transition E α β (chart E α x) = chart E β x := by
    simp only [transition, Function.comp_apply]
    rw [(chart E α).left_inv hα.1]
  set τ := transition E α β with hτdef
  have hτc : ContinuousOn τ (overlap E α β) := (contDiffOn_transition (E := E) α β).continuousOn
  have hN : IsOpen ((overlap E α β ∩ graphDom E J F (fnOf g u) α) ∩ τ ⁻¹' graphDom E J F (fnOf g u) β) :=
    (hτc.mono inter_subset_left).isOpen_inter_preimage
      ((isOpen_overlap (E := E) α β).inter (isOpen_graphDom' hf hJi hJJ hF α))
      (isOpen_graphDom' hf hJi hJJ hF β)
  have heq : chartGraph E J F (fnOf g u) α =ᶠ[𝓝 (chart E α x)]
      chartGraph E J F (fnOf g u) β ∘ τ := by
    filter_upwards [hN.mem_nhds ⟨⟨hy, hyα⟩, by rw [mem_preimage, hτ]; exact hyβ⟩] with w hw
    exact (chartGraph_transition hf hJi hJJ hF hw.1.1 hw.1.2 hw.2).symm
  have hτd : DifferentiableAt ℝ τ (chart E α x) :=
    ((contDiffOn_transition (E := E) α β).contDiffAt
      ((isOpen_overlap (E := E) α β).mem_nhds hy)).differentiableAt (by simp)
  have hGd : DifferentiableAt ℝ (chartGraph E J F (fnOf g u) β) (τ (chart E α x)) := by
    rw [hτ]
    exact (contDiffAt_chartGraph hf hJi hJJ hF hyβ).differentiableAt (by simp)
  have hL : fderiv ℝ (chartGraph E J F (fnOf g u) α) (chart E α x) =
      (fderiv ℝ (chartGraph E J F (fnOf g u) β) (chart E β x)).comp (fderiv ℝ τ (chart E α x)) := by
    rw [heq.fderiv_eq, fderiv_comp _ hGd hτd, hτ]
  rw [Li_covJet hJi hJJ hF u hα hgα, Li_covJet hJi hJJ hF u hβ hgβ]
  exact cometricAt_chart_indep hα.1 hβ.1 (injective_fderiv_chartGraph hf hJi hJJ hF hyβ) hL

/-- **The glued cometric on the jet of a function is the cometric of any chart.** -/
theorem mmEf_covJet (u : SmoothCcTensor g 0 0) (t : ℝ) {α x : M}
    (hζ : covJet g 0 0 (2 * 2) u x ∈ mmDom g J F t x) (hα : x ∈ patchU E α)
    (hgα : good g J F α x (covJet g 0 0 (2 * 2) u x)) :
    mmEf g J F t x (covJet g 0 0 (2 * 2) u x) =
      cometricAt α x (Li g J F α x (covJet g 0 0 (2 * 2) u x)) := by
  unfold mmEf
  have h : ∀ i : PIdx E M, patchρ E M i x •
      cometricAt i.1 x (Li g J F i.1 x (covJet g 0 0 (2 * 2) u x)) =
      patchρ E M i x • cometricAt α x (Li g J F α x (covJet g 0 0 (2 * 2) u x)) := by
    intro i
    by_cases h0 : patchρ E M i x = 0
    · rw [h0, zero_smul, zero_smul]
    · have hx : x ∈ tsupport (patchρ E M i) := subset_tsupport _ (Function.mem_support.2 h0)
      rw [cometricAt_covJet_indep hJi hJJ hF u (tsupport_patchρ_subset i hx) hα (hζ i hx) hgα]
  rw [Finset.sum_congr rfl fun i _ => h i, ← Finset.sum_smul]
  have h1 : ∑ i : PIdx E M, patchρ E M i x = 1 := by
    rw [← finsum_eq_sum_of_fintype]
    exact (patchρ E M).sum_eq_one (mem_univ x)
  rw [h1, one_smul]

/-- **MM's contraction of the top component is the top-order part of the chart operator**, for every jet with
symmetric second-order part. -/
theorem lam_frameSum_cometricAt_eq {α x : M} (hx : x ∈ (chartAt E α).source)
    {ζ : CovJetFiber 𝓘(ℝ, E) 0 0 (2 * 2) x} (hgood : good g J F α x ζ)
    (hsym : ∀ a b, (cJet g α x ζ).2.1 ![a, b] = (cJet g α x ζ).2.1 ![b, a]) :
    lam x (frameSum 2 x (fun _ => cometricAt α x (Li g J F α x ζ)) (stdFrame x) (stdCoframe x)
      (ζ (Fin.last (2 * 2)))) = Fi g J F α x ζ - Fi g J F α x (zeroTopC ζ) := by
  have hlow := cJet_lower g α x (ζ := zeroTopC ζ) (ζ' := ζ) fun _ hk => zeroTopC_apply_of_ne ζ hk
  set z0 := cJet g α x (zeroTopC ζ) with hz0
  set T := topRead (rdJet α x ζ (Fin.last 4)) with hT
  have hc : cJet g α x ζ = addTop z0 T := by
    rw [hz0, hT]
    unfold cJet
    rw [rdJet_zeroTopC]
    exact coordJet_eq_addTop g α (chart E α x) (rdJet α x ζ)
  have hq : (chart E α x, (permTop z0).1) ∈ chartDom E J F α :=
    (good_zeroTopC g J F α x ζ).2 hgood
  have hsym0 : ∀ a b, (permTop z0).2.1 ![a, b] = (permTop z0).2.1 ![b, a] := by
    intro a b
    show z0.2.1 ![a, b] = z0.2.1 ![b, a]
    rw [hlow.2.1]
    exact hsym a b
  have hL0 : LOf (chartPsi E J F α) (chart E α x) (permTop z0) = Li g J F α x ζ :=
    LOf_congr hlow.1 hlow.2.1
  have hinj : Function.Injective (LOf (chartPsi E J F α) (chart E α x) (permTop z0)) :=
    injective_LOf (chartPsi_symplectic hJi hJJ hF α) hq
  have hFi : Fi g J F α x ζ = Fi g J F α x (zeroTopC ζ) -
      trOf (Li g J F α x ζ) fun a c => trOf (Li g J F α x ζ) fun b d =>
        (T.domDomCongr σ4) ![a, b, d, c] := by
    have h := Pjet_add_top (J := J) (chartPsi_symplectic hJi hJJ hF α) hq hsym0 hinj (T.domDomCongr σ4)
    rw [hL0] at h
    unfold Fi
    rw [hc, permTop_addTop]
    exact h
  have hLinj : Function.Injective (Li g J F α x ζ) := injective_Li hJi hJJ hF hgood
  rw [lam_frameSum_cometricAt hx hLinj, hFi]
  have hrd : ∀ a b c d : E, (T.domDomCongr σ4) ![a, b, d, c] =
      rd4 α x (ζ (Fin.last (2 * 2))) ![a, c, b, d] := by
    intro a b c d
    rw [rd4_eq_trivAt hx, hT, ContinuousMultilinearMap.domDomCongr_apply, topRead_apply]
    exact congrArg _ (funext fun i => by fin_cases i <;> rfl)
  simp only [hrd]
  ring

/-- **The right side of the Mantegazza–Martinazzi equation at the jet of a function is the operator of the
flow**, `−div(JH)` of the graph, computed in any chart `α` around `x` in whose Weinstein domain the jet lies. -/
theorem lam_cwData_F (u : SmoothCcTensor g 0 0) (t : ℝ) {α x : M}
    (hζ : covJet g 0 0 (2 * 2) u x ∈ mmDom g J F t x) (hα : x ∈ patchU E α)
    (hgα : good g J F α x (covJet g 0 0 (2 * 2) u x)) :
    lam x ((cwData g hJi hJJ hF).F t x (covJet g 0 0 (2 * 2) u x)) =
      -flowPot J (chartGraph E J F (fnOf g u) α) (chart E α x) := by
  set ζ := covJet g 0 0 (2 * 2) u x with hζdef
  set Q := -flowPot J (chartGraph E J F (fnOf g u) α) (chart E α x) with hQ
  have hF' : (cwData g hJi hJJ hF).F t x ζ =
      frameSum 2 x (fun _ => mmEf g J F t x ζ) (stdFrame x) (stdCoframe x) (ζ (Fin.last (2 * 2))) +
        mmb g J F t x ζ := rfl
  rw [hF', map_add, lam_mmb, mmEf_covJet hJi hJJ hF u t hζ hα hgα]
  have hα' : x ∈ (chartAt E α).source := by
    have := hα.1
    rwa [extChartAt_source] at this
  have hsymα : ∀ a b, (cJet g α x ζ).2.1 ![a, b] = (cJet g α x ζ).2.1 ![b, a] := by
    intro a b
    rw [hζdef, cJet_covJet g u hα]
    exact jetOf_two_symm (contDiffAt_chartRep_fnOf g u hα.1) a b
  set A := lam x (frameSum 2 x (fun _ => cometricAt α x (Li g J F α x ζ)) (stdFrame x)
    (stdCoframe x) (ζ (Fin.last (2 * 2)))) with hA
  -- each active chart contributes `Q − A` to the lower-order term
  have hterm : ∀ i : PIdx E M, patchρ E M i x * Fi g J F i.1 x (zeroTopC ζ) =
      patchρ E M i x * (Q - A) := by
    intro i
    by_cases h0 : patchρ E M i x = 0
    · rw [h0, zero_mul, zero_mul]
    · have hx : x ∈ tsupport (patchρ E M i) := subset_tsupport _ (Function.mem_support.2 h0)
      have hi := tsupport_patchρ_subset i hx
      have hgi := hζ i hx
      have hsymi : ∀ a b, (cJet g i.1 x ζ).2.1 ![a, b] = (cJet g i.1 x ζ).2.1 ![b, a] := by
        intro a b
        rw [hζdef, cJet_covJet g u hi]
        exact jetOf_two_symm (contDiffAt_chartRep_fnOf g u hi.1) a b
      have h1 := lam_frameSum_cometricAt_eq hJi hJJ hF (mem_source_of_mem_tsupport hx) hgi hsymi
      rw [cometricAt_covJet_indep hJi hJJ hF u hi hα hgi hgα] at h1
      have h2 : Fi g J F i.1 x ζ = Q := by
        rw [hQ, Fi_covJet hJi hJJ hF u hi hgi, flowPot_chart_indep hJi hJJ hF u hi hα hgi hgα]
      rw [h2] at h1
      congr 1
      linarith
  unfold mmbS
  rw [Finset.sum_congr rfl fun i _ => hterm i, ← Finset.sum_mul]
  have h1 : ∑ i : PIdx E M, patchρ E M i x = 1 := by
    rw [← finsum_eq_sum_of_fintype]
    exact (patchρ E M).sum_eq_one (mem_univ x)
  rw [h1, one_mul]
  ring

/-- info: 'CWFlow.Bridge.lam_cwData_F' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms lam_cwData_F

end Bridge

end CWFlow
