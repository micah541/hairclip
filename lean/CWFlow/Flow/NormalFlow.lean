import CWFlow.Flow.FlowExistence
import CWFlow.Flow.GraphUnique

/-!
# Existence of a flow from an arbitrary family with the right normal velocity

Chen-Warren, Claim 3.3 for an arbitrary family: if `G_t` is a jointly smooth family of Lagrangian immersions
whose velocity, in every chart, is the flow velocity `J∇div(JH)[G_t]` plus a tangent vector, then composing with
the flow `χ_t` of minus the tangential part of the velocity gives a solution `G_t ∘ χ_t` of the flow.
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

variable (E) in
/-- Minus the tangential coordinates of the velocity of the family `G`, in the chart at `α`. -/
def eRepN (G : ℝ → M → V) (T : ℝ) (α : M) (p : ℝ × E) : E :=
  tanPart (fderiv ℝ (fun z => flowRep E G α (p.1, z)) p.2) (-dtW (flowRep E G α) T p)

variable (E) in
/-- Minus the tangential part of the velocity of `G`, as a tangent vector at `x`. -/
def tanFieldN (G : ℝ → M → V) (T : ℝ) (t : ℝ) (x : M) : TangentSpace 𝓘(ℝ, E) x :=
  (show E from eRepN E G T x (t, chart E x x))

theorem contDiffOn_eRepN {J : V →L[ℝ] V} {G : ℝ → M → V} {T : ℝ} (hT : 0 < T)
    (hsm : ∀ α : M, ContDiffOn ℝ ∞ (flowRep E G α) (Icc 0 T ×ˢ (chart E α).target))
    (himm : ∀ t ∈ Icc (0 : ℝ) T, IsLagImmersion E J (G t)) (α : M) :
    ContDiffOn ℝ ∞ (eRepN E G T α) (Icc 0 T ×ˢ (chart E α).target) := by
  have hG := hsm α
  have hsl := contDiffOn_slice_fderiv hT (isOpen_chart_target (E := E) α) hG
  have hdt := contDiffOn_dtW hT (isOpen_chart_target (E := E) α) hG
  intro p hp
  have hinj : Function.Injective (fderiv ℝ (fun z => flowRep E G α (p.1, z)) p.2) :=
    (himm p.1 hp.1).inj α p.2 hp.2
  exact (contDiffAt_tanPart hinj _).comp_contDiffWithinAt p ((hsl p hp).prodMk (hdt p hp).neg)

theorem eRepN_transition {J : V →L[ℝ] V} {G : ℝ → M → V} {T : ℝ}
    (himm : ∀ t ∈ Icc (0 : ℝ) T, IsLagImmersion E J (G t))
    {t : ℝ} (ht : t ∈ Icc (0 : ℝ) T) {α β x : M} (hα : x ∈ (chart E α).source)
    (hβ : x ∈ (chart E β).source) :
    fderiv ℝ (transition E α β) (chart E α x) (eRepN E G T α (t, chart E α x)) =
      eRepN E G T β (t, chart E β x) := by
  have hG := himm t ht
  have hy : chart E α x ∈ overlap E α β := chart_mem_overlap hα hβ
  have hτ : transition E α β (chart E α x) = chart E β x := transition_chart hα
  obtain ⟨e, he⟩ := exists_equiv_fderiv_transition (E := E) hy
  have heq : chartRep E (G t) α =ᶠ[𝓝 (chart E α x)]
      chartRep E (G t) β ∘ transition E α β := by
    filter_upwards [(isOpen_overlap (E := E) α β).mem_nhds hy] with w hw
    exact chartRep_transition _ hw
  have hdβ : DifferentiableAt ℝ (chartRep E (G t) β)
      (transition E α β (chart E α x)) :=
    ((hG.smooth β).contDiffAt ((isOpen_chart_target (E := E) β).mem_nhds
      (mapsTo_transition (E := E) α β hy))).differentiableAt (by simp)
  have hdτ : DifferentiableAt ℝ (transition E α β) (chart E α x) :=
    ((contDiffOn_transition (E := E) α β).contDiffAt
      ((isOpen_overlap (E := E) α β).mem_nhds hy)).differentiableAt (by simp)
  have hL : fderiv ℝ (chartRep E (G t) α) (chart E α x) =
      (fderiv ℝ (chartRep E (G t) β) (chart E β x)).comp
        (e : E →L[ℝ] E) := by
    rw [heq.fderiv_eq, fderiv_comp _ hdβ hdτ, hτ, he]
  have hv : dtW (flowRep E G α) T (t, chart E α x) = dtW (flowRep E G β) T (t, chart E β x) := by
    show derivWithin (fun s => flowRep E G α (s, chart E α x)) (Icc 0 T) t =
      derivWithin (fun s => flowRep E G β (s, chart E β x)) (Icc 0 T) t
    congr 1
    funext s
    simp only [flowRep]
    rw [(chart E α).left_inv hα, (chart E β).left_inv hβ]
  have key : eRepN E G T α (t, chart E α x) = e.symm (eRepN E G T β (t, chart E β x)) := by
    show tanPart (fderiv ℝ (chartRep E (G t) α) (chart E α x))
        (-dtW (flowRep E G α) T (t, chart E α x)) =
      e.symm (tanPart (fderiv ℝ (chartRep E (G t) β) (chart E β x))
        (-dtW (flowRep E G β) T (t, chart E β x)))
    rw [hL, hv]
    exact tanPart_comp_equiv (hG.inj β _ ((chart E β).map_source hβ)) e _
  rw [← he, key]
  exact e.apply_symm_apply _

theorem trivAt_tanFieldN {J : V →L[ℝ] V} {G : ℝ → M → V} {T : ℝ}
    (himm : ∀ t ∈ Icc (0 : ℝ) T, IsLagImmersion E J (G t))
    {t : ℝ} (ht : t ∈ Icc (0 : ℝ) T) {α x : M} (hα : x ∈ (chart E α).source) :
    (trivializationAt E (TangentSpace 𝓘(ℝ, E)) α).continuousLinearMapAt ℝ x (tanFieldN E G T t x) =
      eRepN E G T α (t, chart E α x) := by
  have hs : x ∈ (chartAt E α).source := by rwa [extChartAt_source] at hα
  rw [TangentBundle.continuousLinearMapAt_trivializationAt_eq_core hs,
    tangentBundleCore_coordChange_achart, ModelWithCorners.range_eq_univ, fderivWithin_univ]
  exact eRepN_transition himm ht (mem_extChartAt_source (I := 𝓘(ℝ, E)) x) hα

theorem contMDiffOn_tanFieldN {J : V →L[ℝ] V} {G : ℝ → M → V} {T : ℝ} (hT : 0 < T)
    (hsm : ∀ α : M, ContDiffOn ℝ ∞ (flowRep E G α) (Icc 0 T ×ˢ (chart E α).target))
    (himm : ∀ t ∈ Icc (0 : ℝ) T, IsLagImmersion E J (G t)) :
    ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) (𝓘(ℝ, E).prod 𝓘(ℝ, E)) ∞
      (fun q : ℝ × M => (TotalSpace.mk' E q.2 (tanFieldN E G T q.1 q.2) : TangentBundle 𝓘(ℝ, E) M))
      (Icc 0 T ×ˢ univ) := by
  rintro ⟨t, x⟩ ⟨ht, -⟩
  rw [contMDiffWithinAt_totalSpace]
  refine ⟨contMDiffWithinAt_snd, ?_⟩
  have hsrc : ∀ q : ℝ × M, q ∈ Icc 0 T ×ˢ (chart E x).source → q.2 ∈ (chartAt E x).source := by
    intro q hq
    have h := hq.2
    rwa [extChartAt_source] at h
  have h1 : ContMDiffWithinAt (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, ℝ × E) ∞
      (fun q : ℝ × M => (q.1, chart E x q.2)) (Icc 0 T ×ˢ (chart E x).source) (t, x) :=
    contMDiffWithinAt_fst.prodMk_space
      ((contMDiffOn_extChartAt (I := 𝓘(ℝ, E)) (n := ∞) (x := x)).comp contMDiffOn_snd hsrc (t, x)
        ⟨ht, mem_extChartAt_source (I := 𝓘(ℝ, E)) x⟩)
  have h2 : ContMDiffWithinAt 𝓘(ℝ, ℝ × E) 𝓘(ℝ, E) ∞ (eRepN E G T x)
      (Icc 0 T ×ˢ (chart E x).target) (t, chart E x x) :=
    (contDiffOn_eRepN hT hsm himm x (t, chart E x x)
      ⟨ht, mem_extChartAt_target (I := 𝓘(ℝ, E)) x⟩).contMDiffWithinAt
  have h3 : ContMDiffWithinAt (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, E) ∞
      (fun q : ℝ × M => eRepN E G T x (q.1, chart E x q.2)) (Icc 0 T ×ˢ (chart E x).source) (t, x) := by
    have hm : MapsTo (fun q : ℝ × M => (q.1, chart E x q.2)) (Icc 0 T ×ˢ (chart E x).source)
        (Icc 0 T ×ˢ (chart E x).target) := fun q hq => ⟨hq.1, (chart E x).map_source hq.2⟩
    exact ContMDiffWithinAt.comp (g := eRepN E G T x)
      (f := fun q : ℝ × M => (q.1, chart E x q.2)) (t, x) h2 h1 hm
  have hev : (fun q : ℝ × M => (trivializationAt E (TangentSpace 𝓘(ℝ, E)) x
        (TotalSpace.mk' E q.2 (tanFieldN E G T q.1 q.2) : TangentBundle 𝓘(ℝ, E) M)).2)
      =ᶠ[𝓝[Icc 0 T ×ˢ (chart E x).source] (t, x)]
      fun q : ℝ × M => eRepN E G T x (q.1, chart E x q.2) := by
    filter_upwards [self_mem_nhdsWithin] with q hq
    have hb : q.2 ∈ (trivializationAt E (TangentSpace 𝓘(ℝ, E)) x).baseSet := by
      rw [TangentBundle.trivializationAt_baseSet]
      exact hsrc q hq
    rw [← trivAt_tanFieldN himm hq.1 hq.2]
    show _ = (trivializationAt E (TangentSpace 𝓘(ℝ, E)) x).linearMapAt ℝ q.2 _
    rw [Trivialization.coe_linearMapAt_of_mem _ hb]
  refine (h3.congr_of_eventuallyEq hev (hev.self_of_nhdsWithin
    ⟨ht, mem_extChartAt_source (I := 𝓘(ℝ, E)) x⟩)).mono_of_mem_nhdsWithin ?_
  exact mem_nhdsWithin.2 ⟨univ ×ˢ (chart E x).source,
    isOpen_univ.prod (isOpen_extChartAt_source (I := 𝓘(ℝ, E)) x),
    ⟨mem_univ _, mem_extChartAt_source (I := 𝓘(ℝ, E)) x⟩, fun p hp => ⟨hp.2.1, hp.1.2⟩⟩

/-- The normal part of the velocity of the family is the flow velocity. -/
theorem dtW_add_fderiv_eRepN {J : V →L[ℝ] V} {G : ℝ → M → V} {T : ℝ}
    (himm : ∀ t ∈ Icc (0 : ℝ) T, IsLagImmersion E J (G t))
    (hvel : ∀ t ∈ Icc (0 : ℝ) T, ∀ α x : M, x ∈ (chart E α).source → ∃ c : E,
      dtW (flowRep E G α) T (t, chart E α x) =
        flowVel J (chartRep E (G t) α) (chart E α x) +
          fderiv ℝ (chartRep E (G t) α) (chart E α x) c)
    {t : ℝ} (ht : t ∈ Icc (0 : ℝ) T) {α x : M} (hα : x ∈ (chart E α).source) :
    dtW (flowRep E G α) T (t, chart E α x) +
        fderiv ℝ (chartRep E (G t) α) (chart E α x) (eRepN E G T α (t, chart E α x)) =
      flowVel J (chartRep E (G t) α) (chart E α x) := by
  have hG := himm t ht
  obtain ⟨c, hc⟩ := hvel t ht α x hα
  have hy := (chart E α).map_source hα
  have h7 := sub_apply_tanPart_of_eq (hG.inj α _ hy) (hG.lag α _ hy)
    (b := ginvOf (fderiv ℝ (chartRep E (G t) α) (chart E α x))
      (fderiv ℝ (flowPot J (chartRep E (G t) α)) (chart E α x))) hc
  show dtW (flowRep E G α) T (t, chart E α x) +
      fderiv ℝ (chartRep E (G t) α) (chart E α x)
        (tanPart (fderiv ℝ (chartRep E (G t) α) (chart E α x))
          (-dtW (flowRep E G α) T (t, chart E α x))) = _
  rw [tanPart_neg', map_neg, ← sub_eq_add_neg, h7]
  try rfl

/-- Along the reparametrization, the family moves by the flow velocity. -/
theorem hasDerivWithinAt_reparamN {J : V →L[ℝ] V} {G : ℝ → M → V} {T : ℝ} (hT : 0 < T)
    (hsm : ∀ α : M, ContDiffOn ℝ ∞ (flowRep E G α) (Icc 0 T ×ˢ (chart E α).target))
    (himm : ∀ t ∈ Icc (0 : ℝ) T, IsLagImmersion E J (G t))
    (hvel : ∀ t ∈ Icc (0 : ℝ) T, ∀ α x : M, x ∈ (chart E α).source → ∃ c : E,
      dtW (flowRep E G α) T (t, chart E α x) =
        flowVel J (chartRep E (G t) α) (chart E α x) +
          fderiv ℝ (chartRep E (G t) α) (chart E α x) c)
    {χ : ℝ → M → M} {t : ℝ} (ht : t ∈ Icc (0 : ℝ) T) {x : M}
    (hχ : HasMFDerivAt 𝓘(ℝ, ℝ) 𝓘(ℝ, E) (fun s => χ s x) t
      ((1 : ℝ →L[ℝ] ℝ).smulRight (tanFieldN E G T t (χ t x))))
    {α : M} (hα : χ t x ∈ (chart E α).source) :
    HasDerivWithinAt (fun s => G s (χ s x))
      (flowVel J (chartRep E (G t) α) (chart E α (χ t x))) (Icc 0 T) t := by
  have hs : χ t x ∈ (chartAt E α).source := by
    have h := hα
    rwa [extChartAt_source] at h
  have hch := (mdifferentiableAt_extChartAt (I := 𝓘(ℝ, E)) hs).hasMFDerivAt
  have hcomp := hch.comp t hχ
  have hfd : HasFDerivAt (fun s => chart E α (χ s x)) (_ : ℝ →L[ℝ] E) t :=
    hasMFDerivAt_iff_hasFDerivAt.mp hcomp
  have hread := trivAt_tanFieldN himm ht hα
  rw [TangentBundle.continuousLinearMapAt_trivializationAt hs] at hread
  have hd : HasDerivAt (fun s => chart E α (χ s x)) (eRepN E G T α (t, chart E α (χ t x))) t := by
    refine hfd.hasDerivAt.congr_deriv ?_
    rw [← hread]
    exact congrArg _ (one_smul ℝ _)
  have hc : HasDerivAt (fun s => (s, chart E α (χ s x)))
      ((1 : ℝ), eRepN E G T α (t, chart E α (χ t x))) t := (hasDerivAt_id t).prodMk hd
  have hN : {s | χ s x ∈ (chart E α).source} ∈ 𝓝 t :=
    hχ.continuousAt.preimage_mem_nhds ((isOpen_extChartAt_source (I := 𝓘(ℝ, E)) α).mem_nhds hα)
  have hG := hasFDerivWithinAt_prod hT (isOpen_chart_target (E := E) α)
    (hsm α) ht ((chart E α).map_source hα)
  have hm : MapsTo (fun s => (s, chart E α (χ s x))) (Icc 0 T ∩ {s | χ s x ∈ (chart E α).source})
      (Icc 0 T ×ˢ (chart E α).target) := fun s hs => ⟨hs.1, (chart E α).map_source hs.2⟩
  have h1 := HasFDerivWithinAt.comp_hasDerivWithinAt (f := fun s => (s, chart E α (χ s x))) t hG
    (hc.hasDerivWithinAt (s := Icc 0 T ∩ {s | χ s x ∈ (chart E α).source})) hm
  have h2 := h1.mono_of_mem_nhdsWithin (inter_mem self_mem_nhdsWithin (nhdsWithin_le_nhds hN))
  have hev : (fun s => G s (χ s x)) =ᶠ[𝓝[Icc 0 T] t]
      flowRep E G α ∘ fun s => (s, chart E α (χ s x)) := by
    filter_upwards [nhdsWithin_le_nhds hN] with s hs
    simp only [flowRep, Function.comp_apply]
    rw [(chart E α).left_inv hs]
  refine (h2.congr_of_eventuallyEq hev (hev.self_of_nhdsWithin ht)).congr_deriv ?_
  rw [fderivWithin_prod_apply hT (isOpen_chart_target (E := E) α)
    (hsm α) ht ((chart E α).map_source hα), one_smul]
  exact dtW_add_fderiv_eRepN himm hvel ht hα

/-- The reparametrized family is jointly smooth on `Icc 0 T × M`. -/
theorem contMDiffOn_reparamN {G : ℝ → M → V} {T : ℝ}
    (hsm : ∀ α : M, ContDiffOn ℝ ∞ (flowRep E G α) (Icc 0 T ×ˢ (chart E α).target))
    {χ : ℝ → M → M}
    (hχ : ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, E) ∞ (fun q : ℝ × M => χ q.1 q.2)
      (Ioo (-1) (T + 1) ×ˢ univ)) :
    ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, V) ∞
      (fun q : ℝ × M => G q.1 (χ q.1 q.2)) (Icc 0 T ×ˢ univ) := by
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
  have h2 : ContMDiffWithinAt 𝓘(ℝ, ℝ × E) 𝓘(ℝ, V) ∞ (flowRep E G α)
      (Icc 0 T ×ˢ (chart E α).target) (t, chart E α (χ t x)) :=
    (hsm α (t, chart E α (χ t x))
      ⟨ht, (chart E α).map_source hsrc⟩).contMDiffWithinAt
  have hm : MapsTo (fun q : ℝ × M => (q.1, chart E α (χ q.1 q.2)))
      ((Icc 0 T ×ˢ univ) ∩ {q : ℝ × M | χ q.1 q.2 ∈ (chart E α).source})
      (Icc 0 T ×ˢ (chart E α).target) := fun q hq => ⟨hq.1.1, (chart E α).map_source hq.2⟩
  have h3 := ContMDiffWithinAt.comp (g := flowRep E G α)
    (f := fun q : ℝ × M => (q.1, chart E α (χ q.1 q.2))) (t, x) h2 h1 hm
  have hev : (fun q : ℝ × M => G q.1 (χ q.1 q.2))
      =ᶠ[𝓝[(Icc 0 T ×ˢ univ) ∩ {q : ℝ × M | χ q.1 q.2 ∈ (chart E α).source}] (t, x)]
      flowRep E G α ∘ fun q : ℝ × M => (q.1, chart E α (χ q.1 q.2)) := by
    filter_upwards [self_mem_nhdsWithin] with q hq
    simp only [flowRep, Function.comp_apply]
    rw [(chart E α).left_inv hq.2]
  exact (h3.congr_of_eventuallyEq hev (hev.self_of_nhdsWithin ⟨⟨ht, mem_univ _⟩, hsrc⟩)).mono_of_mem_nhdsWithin
    (inter_mem self_mem_nhdsWithin hA)

/-- **Chen-Warren, Claim 3.3 for an arbitrary family**: a jointly smooth family `G_t` of Lagrangian immersions
whose velocity, in every chart, is the flow velocity `J∇div(JH)[G_t]` plus a tangent vector, becomes a solution
of the flow after composing with the flow `χ_t` of minus the tangential part of the velocity. -/
theorem exists_flow_of_normal_family {J : V →L[ℝ] V} {G : ℝ → M → V} {T : ℝ} (hT : 0 < T)
    (hsm : ∀ α : M, ContDiffOn ℝ ∞ (flowRep E G α) (Icc 0 T ×ˢ (chart E α).target))
    (himm : ∀ t ∈ Icc (0 : ℝ) T, IsLagImmersion E J (G t))
    (hvel : ∀ t ∈ Icc (0 : ℝ) T, ∀ α x : M, x ∈ (chart E α).source → ∃ c : E,
      dtW (flowRep E G α) T (t, chart E α x) =
        flowVel J (chartRep E (G t) α) (chart E α x) +
          fderiv ℝ (chartRep E (G t) α) (chart E α x) c) :
    ∃ χ : ℝ → (M ≃ₘ⟮𝓘(ℝ, E), 𝓘(ℝ, E)⟯ M), (∀ x, χ 0 x = x) ∧
      IsFlow E J T (fun t x => G t (χ t x)) := by
  obtain ⟨Y, hYsm, hYeq⟩ := DifferentialGeometry.Analysis.ODE.seeley_time_extend (I := 𝓘(ℝ, E))
    (tanFieldN E G T) T hT (contMDiffOn_tanFieldN hT hsm himm)
  obtain ⟨χ, h0, hχsm, -, hode⟩ := exists_reparam_flow Y hYsm hT.le
  refine ⟨χ, h0, ?_, ?_, ?_⟩
  · exact contMDiffOn_reparamN hsm (χ := fun t x => χ t x) hχsm
  · intro t ht
    exact isLagImmersion_comp_diffeo (himm t ht) (χ t)
  · intro t ht x z hz
    have hi : χ t x ∈ (chart E (χ t x)).source := mem_extChartAt_source (I := 𝓘(ℝ, E)) _
    have hχ := hode t (Icc_subset_Ioo_reparam ht) x
    rw [hYeq t ht] at hχ
    have h := hasDerivWithinAt_reparamN hT hsm himm hvel ht hχ hi
    have hG := himm t ht
    have hmem : chart E z x ∈ diffeoDom E (χ t) z (χ t x) := mem_diffeoDom hz hi
    have heq : chartRep E (fun x => G t (χ t x)) z =ᶠ[𝓝 (chart E z x)]
        chartRep E (G t) (χ t x) ∘ diffeoRep E (χ t) z (χ t x) := by
      filter_upwards [(isOpen_diffeoDom (χ t) z (χ t x)).mem_nhds hmem] with w hw
      exact chartRep_comp_diffeo (G t) hw
    have hfv := flowVel_comp' (isOpen_chart_target (E := E) (χ t x)) (hG.smooth (χ t x))
      (hG.inj (χ t x)) J
      (isOpen_diffeoDom (χ t) z (χ t x)) (contDiffOn_diffeoRep (χ t) z (χ t x))
      (fun w hw => exists_equiv_fderiv_diffeoRep (χ t) hw)
      (fun w hw => (chart E (χ t x)).map_source hw.2) hmem
    rw [diffeoRep_chart hz] at hfv
    show HasDerivWithinAt _ (flowVel J (chartRep E (fun x => G t (χ t x)) z)
      (chart E z x)) _ _
    rw [flowVel_congr_of_eventuallyEq heq, hfv]
    exact h

/-- info: 'CWFlow.exists_flow_of_normal_family' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_flow_of_normal_family

end CWFlow
