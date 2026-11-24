import CWFlow.Flow.FlowS
import CWFlow.Flow.NormalFlow

/-!
# Transport of the reparametrization for a normal family

Let `Φ` be an `IsFlowS E J κ T` (its velocity `κ • flowVel` is normal to the image) and `G` a smooth family of
Lagrangian immersions with `Φ_t ∘ ψ_t = G_t` for a smooth family of maps `ψ_t` with `ψ_0 = id`. Then:

* `nt_transport`: `∂ₜψ + Dψ(X) = 0` in charts, `X = tanFieldN G` (minus the tangential part of the velocity of
  `G`); only the normality of `∂ₜΦ` is used;
* `psi_comp_chiN`: hence `ψ_t ∘ χ_t = id` for every flow `χ` of `tanFieldN G` starting at the identity;
* `IsFlow.comp_diffeo`: a flow composed with a fixed diffeomorphism is a flow.
-/

open Set Filter Bundle
open DifferentialGeometry DifferentialGeometry.Tensor0SBundle ShortTime ShortTime.HFlow CWFlow.JetRec
open DifferentialGeometry.Geometry.Curvature DifferentialGeometry.Geometry.Connection
open DifferentialGeometry.Integral.L2
open DifferentialGeometry.Analysis.Sobolev DifferentialGeometry.Analysis.Spectral
open DifferentialGeometry.Analysis.Parabolic.TensorSpectral DifferentialGeometry.Analysis.Elliptic
open scoped InnerProductSpace RealInnerProductSpace Topology ContDiff Manifold
open CWFlow.Bridge

noncomputable section

set_option linter.unusedSectionVars false

namespace CWFlow

universe u

variable {E : Type u} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
  [NeZero (Module.finrank ℝ E)] [CompleteSpace E]
variable {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
  [CompactSpace M] [T2Space M] [SigmaCompactSpace M]
variable {V : Type u} [NormedAddCommGroup V] [InnerProductSpace ℝ V] [FiniteDimensional ℝ V]

/-- **The transport equation** `∂ₜψ + Dψ(X) = 0`, `X = tanFieldN G`, in charts. -/
theorem nt_transport {J : V →L[ℝ] V} {κ T : ℝ} (hT : 0 < T) {G Φ : ℝ → M → V}
    (hsm : ∀ α : M, ContDiffOn ℝ ∞ (flowRep E G α) (Icc 0 T ×ˢ (chart E α).target))
    (himm : ∀ t ∈ Icc (0 : ℝ) T, IsLagImmersion E J (G t))
    (hΦ : IsFlowS E J κ T Φ) {ψ : ℝ → M → M}
    (hψ : ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, E) ∞ (fun q : ℝ × M => ψ q.1 q.2) (Icc 0 T ×ˢ univ))
    (hG : ∀ t ∈ Icc (0 : ℝ) T, ∀ x : M, Φ t (ψ t x) = G t x)
    {t : ℝ} (ht : t ∈ Icc (0 : ℝ) T) {α y β : M} (hα : y ∈ (chart E α).source)
    (hβ : ψ t y ∈ (chart E β).source) :
    fderivWithin ℝ (psiRep E ψ α β) (Icc 0 T ×ˢ (chart E α).target) (t, chart E α y)
      (1, eRepN E G T α (t, chart E α y)) = 0 := by
  have hzt : chart E α y ∈ (chart E α).target := (chart E α).map_source hα
  have hsymm : (chart E α).symm (chart E α y) = y := (chart E α).left_inv hα
  have hβz : ψ t ((chart E α).symm (chart E α y)) ∈ (chart E β).source := by rw [hsymm]; exact hβ
  have hceq : psiRep E ψ α β (t, chart E α y) = chart E β (ψ t y) := by
    simp only [psiRep]
    rw [hsymm]
  have hct : chart E β (ψ t y) ∈ (chart E β).target := (chart E β).map_source hβ
  have hP := contDiffWithinAt_psiRep hψ ht hzt hβz
  have hUD : UniqueDiffWithinAt ℝ (Icc (0 : ℝ) T ×ˢ (chart E α).target) (t, chart E α y) :=
    ((uniqueDiffOn_Icc hT).prod (isOpen_chart_target (E := E) α).uniqueDiffOn) _ ⟨ht, hzt⟩
  have hPd : HasFDerivWithinAt (psiRep E ψ α β)
      (fderivWithin ℝ (psiRep E ψ α β) (Icc 0 T ×ˢ (chart E α).target) (t, chart E α y))
      (Icc 0 T ×ˢ (chart E α).target) (t, chart E α y) :=
    (hP.differentiableWithinAt (by simp)).hasFDerivWithinAt
  have hN : {q : ℝ × E | ψ q.1 ((chart E α).symm q.2) ∈ (chart E β).source} ∈
      𝓝[Icc 0 T ×ˢ (chart E α).target] (t, chart E α y) :=
    (contMDiffOn_psiPre hψ α (t, chart E α y) ⟨ht, hzt⟩).continuousWithinAt.preimage_mem_nhdsWithin
      ((isOpen_extChartAt_source (I := 𝓘(ℝ, E)) β).mem_nhds hβz)
  have hΦrep := contDiffOn_flowRep hΦ.smooth β
  have hΦd := hasFDerivWithinAt_prod hT (isOpen_chart_target (E := E) β) hΦrep ht hct
  rw [← hceq] at hΦd
  have hκ : HasFDerivWithinAt (fun q : ℝ × E => (q.1, psiRep E ψ α β q))
      ((ContinuousLinearMap.fst ℝ ℝ E).prod
        (fderivWithin ℝ (psiRep E ψ α β) (Icc 0 T ×ˢ (chart E α).target) (t, chart E α y)))
      (Icc 0 T ×ˢ (chart E α).target) (t, chart E α y) := hasFDerivWithinAt_fst.prodMk hPd
  have hmaps : MapsTo (fun q : ℝ × E => (q.1, psiRep E ψ α β q))
      ((Icc 0 T ×ˢ (chart E α).target) ∩
        {q : ℝ × E | ψ q.1 ((chart E α).symm q.2) ∈ (chart E β).source})
      (Icc 0 T ×ˢ (chart E β).target) := fun q hq => ⟨hq.1.1, (chart E β).map_source hq.2⟩
  have h1 := HasFDerivWithinAt.comp (g := flowRep E Φ β)
    (f := fun q : ℝ × E => (q.1, psiRep E ψ α β q)) (t, chart E α y) hΦd
    (hκ.mono inter_subset_left) hmaps
  have h2 := h1.mono_of_mem_nhdsWithin (inter_mem self_mem_nhdsWithin hN)
  have hfam := hsm α
  have hfd := hasFDerivWithinAt_prod hT (isOpen_chart_target (E := E) α) hfam ht hzt
  have hev : flowRep E G α =ᶠ[𝓝[Icc 0 T ×ˢ (chart E α).target] (t, chart E α y)]
      flowRep E Φ β ∘ fun q : ℝ × E => (q.1, psiRep E ψ α β q) := by
    filter_upwards [inter_mem self_mem_nhdsWithin hN] with q hq
    simp only [flowRep, psiRep, Function.comp_apply]
    rw [(chart E β).left_inv hq.2]
    exact (hG q.1 hq.1.1 _).symm
  have h3 := h2.congr_of_eventuallyEq hev (hev.self_of_nhdsWithin ⟨ht, hzt⟩)
  have hEq := hUD.eq hfd h3
  have hEq1 := congrArg (fun L : ℝ × E →L[ℝ] V => L (1, eRepN E G T α (t, chart E α y))) hEq
  simp only [ContinuousLinearMap.comp_apply, ContinuousLinearMap.prod_apply,
    ContinuousLinearMap.coe_fst'] at hEq1
  rw [fderivWithin_prod_apply hT (isOpen_chart_target (E := E) α) hfam ht hzt, hceq,
    fderivWithin_prod_apply hT (isOpen_chart_target (E := E) β) hΦrep ht hct, one_smul,
    one_smul] at hEq1
  -- the time derivative of the flow
  have hdt : dtW (flowRep E Φ β) T (t, chart E β (ψ t y)) =
      κ • flowVel J (chartRep E (Φ t) β) (chart E β (ψ t y)) := by
    have hfun : (fun s => flowRep E Φ β (s, chart E β (ψ t y))) = fun s => Φ s (ψ t y) := by
      funext s
      simp only [flowRep]
      rw [(chart E β).left_inv hβ]
    show derivWithin (fun s => flowRep E Φ β (s, chart E β (ψ t y))) (Icc 0 T) t = _
    rw [hfun]
    exact (hΦ.vel t ht (ψ t y) β hβ).derivWithin (uniqueDiffOn_Icc hT t ht)
  rw [hdt] at hEq1
  generalize fderivWithin ℝ (psiRep E ψ α β) (Icc 0 T ×ˢ (chart E α).target) (t, chart E α y)
    (1, eRepN E G T α (t, chart E α y)) = P at hEq1 ⊢
  -- chart readings
  have hGt := himm t ht
  have hΦt := hΦ.imm t ht
  have hψt := contMDiff_slice hψ ht
  have hmem : chart E α y ∈ diffeoDom E (ψ t) α β := mem_diffeoDom hα hβ
  have hdr : diffeoRep E (ψ t) α β (chart E α y) = chart E β (ψ t y) := diffeoRep_chart hα
  have heq : chartRep E (G t) α =ᶠ[𝓝 (chart E α y)]
      chartRep E (Φ t) β ∘ diffeoRep E (ψ t) α β := by
    filter_upwards [(isOpen_diffeoDom' hψt.continuous α β).mem_nhds hmem] with w' hw'
    have h := chartRep_comp_diffeo (Φ t) hw'
    rw [Function.comp_apply, ← h]
    simp only [chartRep, Function.comp_apply]
    exact (hG t ht _).symm
  have hρ : DifferentiableAt ℝ (diffeoRep E (ψ t) α β) (chart E α y) :=
    ((contDiffOn_diffeoRep' hψt α β).contDiffAt
      ((isOpen_diffeoDom' hψt.continuous α β).mem_nhds hmem)).differentiableAt (by simp)
  have hX : DifferentiableAt ℝ (chartRep E (Φ t) β) (diffeoRep E (ψ t) α β (chart E α y)) := by
    rw [hdr]
    exact ((hΦt.smooth β).contDiffAt
      ((isOpen_chart_target (E := E) β).mem_nhds hct)).differentiableAt (by simp)
  have hLG : fderiv ℝ (chartRep E (G t) α) (chart E α y) =
      (fderiv ℝ (chartRep E (Φ t) β) (chart E β (ψ t y))).comp
        (fderiv ℝ (diffeoRep E (ψ t) α β) (chart E α y)) := by
    rw [heq.fderiv_eq, fderiv_comp _ hX hρ, hdr]
  have hinjG : Function.Injective (fderiv ℝ (chartRep E (G t) α) (chart E α y)) :=
    hGt.inj α _ hzt
  have hinjD : Function.Injective (fderiv ℝ (diffeoRep E (ψ t) α β) (chart E α y)) := by
    intro a b hab
    apply hinjG
    rw [hLG, ContinuousLinearMap.comp_apply, ContinuousLinearMap.comp_apply, hab]
  have hsurjD : Function.Surjective (fderiv ℝ (diffeoRep E (ψ t) α β) (chart E α y)) :=
    (LinearMap.injective_iff_surjective
      (f := (fderiv ℝ (diffeoRep E (ψ t) α β) (chart E α y) : E →ₗ[ℝ] E))).mp hinjD
  -- the normal part of the velocity of `G`
  have hA : ∀ a, ⟪dtW (flowRep E G α) T (t, chart E α y) +
      fderiv ℝ (chartRep E (G t) α) (chart E α y) (eRepN E G T α (t, chart E α y)),
      fderiv ℝ (chartRep E (G t) α) (chart E α y) a⟫ = 0 := by
    intro a
    have h := inner_sub_tanPart hinjG (dtW (flowRep E G α) T (t, chart E α y)) a
    show ⟪dtW (flowRep E G α) T (t, chart E α y) +
      fderiv ℝ (chartRep E (G t) α) (chart E α y) (tanPart (fderiv ℝ (chartRep E (G t) α) (chart E α y))
        (-dtW (flowRep E G α) T (t, chart E α y))), _⟫ = 0
    rw [tanPart_neg', map_neg, ← sub_eq_add_neg]
    exact h
  set LΦ := fderiv ℝ (chartRep E (Φ t) β) (chart E β (ψ t y)) with hLΦ
  have hAΦ : ∀ b, ⟪dtW (flowRep E G α) T (t, chart E α y) +
      fderiv ℝ (chartRep E (G t) α) (chart E α y) (eRepN E G T α (t, chart E α y)), LΦ b⟫ = 0 := by
    intro b
    obtain ⟨a, rfl⟩ := hsurjD b
    have h := congrArg (fun L : E →L[ℝ] V => L a) hLG
    simp only [ContinuousLinearMap.comp_apply] at h
    rw [← h]
    exact hA a
  have hfv : ∀ b, ⟪κ • flowVel J (chartRep E (Φ t) β) (chart E β (ψ t y)), LΦ b⟫ = 0 := by
    intro b
    rw [inner_smul_left]
    have : ⟪flowVel J (chartRep E (Φ t) β) (chart E β (ψ t y)), LΦ b⟫ = 0 :=
      hΦt.lag β _ hct _ b
    rw [this]
    simp
  have hLP : LΦ P = dtW (flowRep E G α) T (t, chart E α y) +
      fderiv ℝ (chartRep E (G t) α) (chart E α y) (eRepN E G T α (t, chart E α y)) -
      κ • flowVel J (chartRep E (Φ t) β) (chart E β (ψ t y)) := by
    exact eq_sub_of_add_eq' hEq1.symm
  have h0 : ⟪LΦ P, LΦ P⟫ = 0 := by
    have h : ⟪LΦ P, LΦ P⟫ = _ := congrArg (fun v => ⟪v, LΦ P⟫) hLP
    rw [h, inner_sub_left, hAΦ, hfv, sub_zero]
  exact hΦt.inj β _ hct ((inner_self_eq_zero.mp h0).trans (map_zero _).symm)

/-- Along a curve of the field `tanFieldN G`, the chart readings of `ψ` have zero derivative. -/
theorem nt_hasDerivWithinAt_chart_psi_chi {J : V →L[ℝ] V} {κ T : ℝ} (hT : 0 < T) {G Φ : ℝ → M → V}
    (hsm : ∀ α : M, ContDiffOn ℝ ∞ (flowRep E G α) (Icc 0 T ×ˢ (chart E α).target))
    (himm : ∀ t ∈ Icc (0 : ℝ) T, IsLagImmersion E J (G t))
    (hΦ : IsFlowS E J κ T Φ) {ψ : ℝ → M → M}
    (hψ : ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, E) ∞ (fun q : ℝ × M => ψ q.1 q.2) (Icc 0 T ×ˢ univ))
    (hG : ∀ t ∈ Icc (0 : ℝ) T, ∀ x : M, Φ t (ψ t x) = G t x)
    {χ : ℝ → M → M} {t : ℝ} (ht : t ∈ Icc (0 : ℝ) T) {x : M}
    (hχ : HasMFDerivAt 𝓘(ℝ, ℝ) 𝓘(ℝ, E) (fun s => χ s x) t
      ((1 : ℝ →L[ℝ] ℝ).smulRight (tanFieldN E G T t (χ t x))))
    {β : M} (hβ : ψ t (χ t x) ∈ (chart E β).source) :
    HasDerivWithinAt (fun s => chart E β (ψ s (χ s x))) 0 (Icc 0 T) t := by
  have hα : χ t x ∈ (chart E (χ t x)).source := mem_extChartAt_source (I := 𝓘(ℝ, E)) _
  have hs : χ t x ∈ (chartAt E (χ t x)).source := by
    have h := hα
    rwa [extChartAt_source] at h
  have hch := (mdifferentiableAt_extChartAt (I := 𝓘(ℝ, E)) hs).hasMFDerivAt
  have hcomp := hch.comp t hχ
  have hfd : HasFDerivAt (fun s => chart E (χ t x) (χ s x)) (_ : ℝ →L[ℝ] E) t :=
    hasMFDerivAt_iff_hasFDerivAt.mp hcomp
  have hread := trivAt_tanFieldN himm ht hα
  rw [TangentBundle.continuousLinearMapAt_trivializationAt hs] at hread
  have hd : HasDerivAt (fun s => chart E (χ t x) (χ s x))
      (eRepN E G T (χ t x) (t, chart E (χ t x) (χ t x))) t := by
    refine hfd.hasDerivAt.congr_deriv ?_
    rw [← hread]
    exact congrArg _ (one_smul ℝ _)
  have hc : HasDerivAt (fun s => (s, chart E (χ t x) (χ s x)))
      ((1 : ℝ), eRepN E G T (χ t x) (t, chart E (χ t x) (χ t x))) t := (hasDerivAt_id t).prodMk hd
  have hN : {s | χ s x ∈ (chart E (χ t x)).source} ∈ 𝓝 t :=
    hχ.continuousAt.preimage_mem_nhds ((isOpen_extChartAt_source (I := 𝓘(ℝ, E)) (χ t x)).mem_nhds hα)
  have hzt : chart E (χ t x) (χ t x) ∈ (chart E (χ t x)).target := (chart E (χ t x)).map_source hα
  have hβz : ψ t ((chart E (χ t x)).symm (chart E (χ t x) (χ t x))) ∈ (chart E β).source := by
    rw [(chart E (χ t x)).left_inv hα]
    exact hβ
  have hP := contDiffWithinAt_psiRep hψ (α := χ t x) (β := β) ht hzt hβz
  have hPd := (hP.differentiableWithinAt (by simp)).hasFDerivWithinAt
  have hm : MapsTo (fun s => (s, chart E (χ t x) (χ s x)))
      (Icc 0 T ∩ {s | χ s x ∈ (chart E (χ t x)).source}) (Icc 0 T ×ˢ (chart E (χ t x)).target) :=
    fun s hs => ⟨hs.1, (chart E (χ t x)).map_source hs.2⟩
  have h1 := HasFDerivWithinAt.comp_hasDerivWithinAt (f := fun s => (s, chart E (χ t x) (χ s x))) t hPd
    (hc.hasDerivWithinAt (s := Icc 0 T ∩ {s | χ s x ∈ (chart E (χ t x)).source})) hm
  have h2 := h1.mono_of_mem_nhdsWithin (inter_mem self_mem_nhdsWithin (nhdsWithin_le_nhds hN))
  have hev : (fun s => chart E β (ψ s (χ s x))) =ᶠ[𝓝[Icc 0 T] t]
      psiRep E ψ (χ t x) β ∘ fun s => (s, chart E (χ t x) (χ s x)) := by
    filter_upwards [nhdsWithin_le_nhds hN] with s hs
    simp only [psiRep, Function.comp_apply]
    rw [(chart E (χ t x)).left_inv hs]
  refine (h2.congr_of_eventuallyEq hev (hev.self_of_nhdsWithin ht)).congr_deriv ?_
  exact nt_transport hT hsm himm hΦ hψ hG ht hα hβ

/-- **`ψ_t` is the inverse of the reparametrization `χ_t`**: `ψ_t ∘ χ_t = id` for every flow `χ` of the field
`tanFieldN G` starting at the identity, when `Φ_t ∘ ψ_t = G_t` and `Φ` has normal velocity. -/
theorem psi_comp_chiN {J : V →L[ℝ] V} (hJi : ∀ u v, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v, J (J v) = -v)
    {κ T : ℝ} (hT : 0 < T) {G Φ : ℝ → M → V}
    (hsm : ∀ α : M, ContDiffOn ℝ ∞ (flowRep E G α) (Icc 0 T ×ˢ (chart E α).target))
    (himm : ∀ t ∈ Icc (0 : ℝ) T, IsLagImmersion E J (G t))
    (hΦ : IsFlowS E J κ T Φ) {ψ : ℝ → M → M}
    (hψ : ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, E) ∞ (fun q : ℝ × M => ψ q.1 q.2) (Icc 0 T ×ˢ univ))
    (hψ0 : ∀ x, ψ 0 x = x)
    (hG : ∀ t ∈ Icc (0 : ℝ) T, ∀ x : M, Φ t (ψ t x) = G t x)
    {χ : ℝ → M → M} (hχ0 : ∀ x, χ 0 x = x)
    (hχ : ∀ t ∈ Icc (0 : ℝ) T, ∀ x, HasMFDerivAt 𝓘(ℝ, ℝ) 𝓘(ℝ, E) (fun s => χ s x) t
      ((1 : ℝ →L[ℝ] ℝ).smulRight (tanFieldN E G T t (χ t x))))
    {t : ℝ} (ht : t ∈ Icc (0 : ℝ) T) (x : M) : ψ t (χ t x) = x := by
  have hc : ContinuousOn (fun s => ψ s (χ s x)) (Icc 0 T) := by
    intro s hs
    have h1 : ContinuousWithinAt (fun q : ℝ × M => ψ q.1 q.2) (Icc 0 T ×ˢ univ) (s, χ s x) :=
      (hψ (s, χ s x) ⟨hs, mem_univ _⟩).continuousWithinAt
    have h2 : ContinuousWithinAt (fun r => ((r, χ r x) : ℝ × M)) (Icc 0 T) s :=
      (continuousAt_id.prodMk (hχ s hs x).continuousAt).continuousWithinAt
    exact ContinuousWithinAt.comp (g := fun q : ℝ × M => ψ q.1 q.2)
      (f := fun r => ((r, χ r x) : ℝ × M)) h1 h2 (fun r hr => ⟨hr, mem_univ _⟩)
  have h := eq_of_chart_hasDerivWithinAt_zero (E := E) hc (fun s hs β hβ =>
    nt_hasDerivWithinAt_chart_psi_chi hT hsm himm hΦ hψ hG hs (hχ s hs x) hβ) t ht
  have h' : ψ t (χ t x) = ψ 0 (χ 0 x) := h
  rw [h', hχ0, hψ0]

/-- info: 'CWFlow.psi_comp_chiN' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms psi_comp_chiN

/-- **A flow composed with a fixed diffeomorphism is a flow.** -/
theorem IsFlow.comp_diffeo {J : V →L[ℝ] V} {T : ℝ} {Φ : ℝ → M → V} (hΦ : IsFlow E J T Φ)
    (φ : M ≃ₘ⟮𝓘(ℝ, E), 𝓘(ℝ, E)⟯ M) : IsFlow E J T (fun t x => Φ t (φ x)) where
  smooth := by
    have hf : ContMDiff ((𝓘(ℝ, ℝ).prod 𝓘(ℝ, E))) ((𝓘(ℝ, ℝ).prod 𝓘(ℝ, E))) ∞
        (fun q : ℝ × M => ((q.1, φ q.2) : ℝ × M)) :=
      contMDiff_fst.prodMk (φ.contMDiff.comp contMDiff_snd)
    exact ContMDiffOn.comp (g := fun q : ℝ × M => Φ q.1 q.2) hΦ.smooth hf.contMDiffOn
      (fun q hq => ⟨hq.1, mem_univ _⟩)
  imm := fun t ht => isLagImmersion_comp_diffeo (hΦ.imm t ht) φ
  vel := by
    intro t ht x z hz
    have hi : φ x ∈ (chart E (φ x)).source := mem_extChartAt_source (I := 𝓘(ℝ, E)) _
    have h := hΦ.vel t ht (φ x) (φ x) hi
    have hG := hΦ.imm t ht
    have hmem : chart E z x ∈ diffeoDom E φ z (φ x) := mem_diffeoDom hz hi
    have heq : chartRep E (fun x => Φ t (φ x)) z =ᶠ[𝓝 (chart E z x)]
        chartRep E (Φ t) (φ x) ∘ diffeoRep E φ z (φ x) := by
      filter_upwards [(isOpen_diffeoDom φ z (φ x)).mem_nhds hmem] with w hw
      exact chartRep_comp_diffeo (Φ t) hw
    have hfv := flowVel_comp' (isOpen_chart_target (E := E) (φ x)) (hG.smooth (φ x))
      (hG.inj (φ x)) J
      (isOpen_diffeoDom φ z (φ x)) (contDiffOn_diffeoRep φ z (φ x))
      (fun w hw => exists_equiv_fderiv_diffeoRep φ hw)
      (fun w hw => (chart E (φ x)).map_source hw.2) hmem
    rw [diffeoRep_chart hz] at hfv
    show HasDerivWithinAt _ (flowVel J (chartRep E (fun x => Φ t (φ x)) z) (chart E z x)) _ _
    rw [flowVel_congr_of_eventuallyEq heq, hfv]
    exact h

/-- info: 'CWFlow.IsFlow.comp_diffeo' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms IsFlow.comp_diffeo

end CWFlow
