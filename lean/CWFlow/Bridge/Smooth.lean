import CWFlow.Bridge.Fields

/-!
# The Mantegazza–Martinazzi data of the scalar equation: openness and smoothness

In MM's local coordinates `q = (t, x, z)` at `x₀` (`z` the jet in the trivialisations at `x₀`):

* `rdQ x₀ α q`: the jet read in the chart at `α`; smooth on the overlap (`contMDiffAt_rdQ`, from MM's
  `contMDiffAt_trivJet_jetOfModel`);
* `cJQ g φ x₀ α q = coordJet g α (chart α x) (φ (rdQ x₀ α q))`: the coordinate jet, smooth over the patch of `α`;
* `mmDom_relOpen`: the domain is relatively open;
* `contMDiffOn_mmEf`, `contMDiffOn_mmb`: the coefficients are smooth on the domain.
-/

open Set Filter Bundle
open DifferentialGeometry DifferentialGeometry.Tensor0SBundle ShortTime ShortTime.HFlow CWFlow.JetRec
open DifferentialGeometry.Geometry.Curvature DifferentialGeometry.Geometry.Connection
open scoped InnerProductSpace RealInnerProductSpace Topology ContDiff Manifold

noncomputable section

set_option linter.unusedSectionVars false

namespace CWFlow

universe u

section LOf

variable {E V : Type u} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
  [NormedAddCommGroup V] [InnerProductSpace ℝ V]
  {Ψ : E × (E →L[ℝ] ℝ) → V} {O : Set (E × (E →L[ℝ] ℝ))}

theorem contDiffAt_LOf (hO : IsOpen O) (hΨ : ContDiffOn ℝ ∞ Ψ O) {p : E × SJet E}
    (hp : (p.1, p.2.1) ∈ O) : ContDiffAt ℝ ∞ (fun p : E × SJet E => LOf Ψ p.1 p.2) p := by
  have he : ContDiff ℝ ∞ (continuousMultilinearCurryFin1 ℝ E V) :=
    LinearIsometryEquiv.contDiff (𝕜 := ℝ) _
  exact he.contDiffAt.comp p (contDiffAt_TkOf hO hΨ hp 1)

end LOf

namespace Bridge


variable {E : Type u} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
  [NeZero (Module.finrank ℝ E)] [CompleteSpace E]
variable {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
  [CompactSpace M] [T2Space M] [SigmaCompactSpace M]
variable {V : Type u} [NormedAddCommGroup V] [InnerProductSpace ℝ V]

set_option quotPrecheck false in
/-- The model of MM's local coordinates `(t, x, z)`. -/
local notation "𝓠" => (𝓘(ℝ, ℝ).prod ((𝓘(ℝ, E)).prod 𝓘(ℝ, JetModel E 0 0 (2 * 2))))

/-! ### Linear pieces -/

theorem contDiff_zeroTop : ContDiff ℝ ∞ (zeroTop : JetModel E 0 0 4 → JetModel E 0 0 4) := by
  rw [contDiff_pi]
  intro k
  by_cases hk : k = Fin.last 4
  · subst hk
    simp only [zeroTop, zeroLast, Function.update_self]
    exact contDiff_const
  · simp only [zeroTop_apply_of_ne _ hk]
    exact (ContinuousLinearMap.proj (R := ℝ)
      (φ := fun k : Fin (4 + 1) => TensorRSModel 0 (0 + (k : ℕ)) ℝ E) k).contDiff

theorem contDiff_permTop : ContDiff ℝ ∞ (permTop : SJet E → SJet E) := by
  have h4 : ContDiff ℝ ∞ (fun T : E [×4]→L[ℝ] ℝ => T.domDomCongr σ4) :=
    LinearIsometryEquiv.contDiff (𝕜 := ℝ)
      (ContinuousMultilinearMap.domDomCongrₗᵢ ℝ E ℝ σ4)
  exact contDiff_fst.prodMk ((contDiff_fst.comp contDiff_snd).prodMk
    ((contDiff_fst.comp (contDiff_snd.comp contDiff_snd)).prodMk
      (h4.comp (contDiff_snd.comp (contDiff_snd.comp contDiff_snd)))))

theorem rdJet_zeroTopC (α x : M) (ζ : CovJetFiber 𝓘(ℝ, E) 0 0 (2 * 2) x) :
    rdJet α x (zeroTopC ζ) = zeroTop (rdJet α x ζ) := by
  funext k
  by_cases hk : k = Fin.last 4
  · subst hk
    rw [zeroTop, zeroLast, Function.update_self]
    show (trivAt 𝓘(ℝ, E) M 0 (0 + _) α).continuousLinearMapAt ℝ x (zeroTopC ζ (Fin.last (2 * 2))) = 0
    rw [zeroTopC_last, map_zero]
  · rw [zeroTop_apply_of_ne _ hk]
    show (trivAt 𝓘(ℝ, E) M 0 (0 + k) α).continuousLinearMapAt ℝ x (zeroTopC ζ k) =
      (trivAt 𝓘(ℝ, E) M 0 (0 + k) α).continuousLinearMapAt ℝ x (ζ k)
    rw [zeroTopC_apply_of_ne ζ hk]

/-- A model `(0,0)`-tensor is its scalar value times the identity. -/
theorem eq_scal_smul_id (T : TensorRSModel 0 0 ℝ E) :
    T = scal E T • (ContinuousLinearMap.id ℝ (Tensor0SModel 0 ℝ E) : TensorRSModel 0 0 ℝ E) := by
  refine ContinuousLinearMap.ext fun m => ?_
  have hm : m = m (fun i => i.elim0) • unitM E := by
    refine ContinuousMultilinearMap.ext fun w => ?_
    have hw : w = fun i => i.elim0 := Subsingleton.elim _ _
    subst hw
    show _ = m (fun i => i.elim0) * unitM E _
    rw [unitM_apply, mul_one]
  refine ContinuousMultilinearMap.ext fun v => ?_
  have hv : v = fun i => i.elim0 := Subsingleton.elim _ _
  subst hv
  show T m (fun i => i.elim0) = scal E T * m (fun i => i.elim0)
  conv_lhs => rw [hm, map_smul]
  show m (fun i => i.elim0) * T (unitM E) (fun i => i.elim0) = _
  rw [scal_apply, mul_comm]

/-! ### The jet read in a chart -/

/-- The jet with coordinates `z` in the trivialisations at `x₀`, read in the chart at `α`. -/
def rdQ (x₀ α : M) (q : ℝ × M × JetModel E 0 0 (2 * 2)) : JetModel E 0 0 4 :=
  rdJet α q.2.1 (jetOfModel 0 0 (2 * 2) x₀ q.2.1 q.2.2)

theorem contMDiffAt_snd_fst (q : ℝ × M × JetModel E 0 0 (2 * 2)) :
    ContMDiffAt 𝓠 𝓘(ℝ, E) ∞ (fun q : ℝ × M × JetModel E 0 0 (2 * 2) => q.2.1) q :=
  contMDiffAt_fst.comp q contMDiffAt_snd

theorem contMDiffAt_rdQ (x₀ α : M) {q : ℝ × M × JetModel E 0 0 (2 * 2)}
    (h₀ : q.2.1 ∈ jetBase (I := 𝓘(ℝ, E)) 0 0 (2 * 2) x₀) (hα : q.2.1 ∈ (chartAt E α).source) :
    ContMDiffAt 𝓠 𝓘(ℝ, JetModel E 0 0 4) ∞ (rdQ x₀ α) q :=
  contMDiffAt_trivJet_jetOfModel (I := 𝓘(ℝ, E)) (s := 0) (N := 2 * 2)
    (Nn := ℝ × M × JetModel E 0 0 (2 * 2)) x₀ α
    (fun q : ℝ × M × JetModel E 0 0 (2 * 2) => q.2.1)
    (fun q : ℝ × M × JetModel E 0 0 (2 * 2) => q.2.2) (contMDiffAt_snd_fst q) (contMDiffAt_snd.comp q contMDiffAt_snd) h₀
    (mem_jetBase_of_source (by rwa [extChartAt_source]))

theorem contMDiffAt_chart (α : M) {q : ℝ × M × JetModel E 0 0 (2 * 2)}
    (hα : q.2.1 ∈ (chartAt E α).source) :
    ContMDiffAt 𝓠 𝓘(ℝ, E) ∞ (fun q : ℝ × M × JetModel E 0 0 (2 * 2) => chart E α q.2.1) q :=
  (contMDiffAt_extChartAt' (I := 𝓘(ℝ, E)) (n := ∞) hα).comp q (contMDiffAt_snd_fst q)

variable (g : SmoothRiemannianMetric 𝓘(ℝ, E) M) (J : V →L[ℝ] V) (F : M → V)

/-- The chart point and the coordinate jet of `φ (rdQ x₀ α q)`. -/
def cJQ (φ : JetModel E 0 0 4 → JetModel E 0 0 4) (x₀ α : M)
    (q : ℝ × M × JetModel E 0 0 (2 * 2)) : E × SJet E :=
  (chart E α q.2.1, coordJet g α (chart E α q.2.1) (φ (rdQ x₀ α q)))

theorem contMDiffAt_cJQ {φ : JetModel E 0 0 4 → JetModel E 0 0 4} (hφ : ContDiff ℝ ∞ φ) (x₀ α : M)
    {q : ℝ × M × JetModel E 0 0 (2 * 2)}
    (h₀ : q.2.1 ∈ jetBase (I := 𝓘(ℝ, E)) 0 0 (2 * 2) x₀) (hα : q.2.1 ∈ patchU E α) :
    ContMDiffAt 𝓠 𝓘(ℝ, E × SJet E) ∞ (cJQ g φ x₀ α) q := by
  have hs : q.2.1 ∈ (chartAt E α).source := by
    have := patchU_subset_source α hα
    rwa [extChartAt_source] at this
  have hc : ContDiffAt ℝ ∞ (fun p : E × JetModel E 0 0 4 => coordJet g α p.1 p.2)
      (chart E α q.2.1, φ (rdQ x₀ α q)) :=
    (contDiffOn_coordJet g α).contDiffAt
      (((isOpen_patchΩ (I := 𝓘(ℝ, E)) α).prod isOpen_univ).mem_nhds ⟨hα.2, mem_univ _⟩)
  have h1 := contMDiffAt_chart α hs
  have h2 : ContMDiffAt 𝓠 𝓘(ℝ, JetModel E 0 0 4) ∞ (fun q => φ (rdQ x₀ α q)) q :=
    hφ.contDiffAt.comp_contMDiffAt (contMDiffAt_rdQ x₀ α h₀ hs)
  have h3 : ContMDiffAt 𝓠 𝓘(ℝ, SJet E) ∞
      (fun q => coordJet g α (chart E α q.2.1) (φ (rdQ x₀ α q))) q :=
    hc.comp_contMDiffAt (f := fun q => (chart E α q.2.1, φ (rdQ x₀ α q))) (h1.prodMk_space h2)
  exact h1.prodMk_space h3

theorem cJQ_id (x₀ α : M) (q : ℝ × M × JetModel E 0 0 (2 * 2)) :
    cJQ g id x₀ α q =
      (chart E α q.2.1, cJet g α q.2.1 (jetOfModel 0 0 (2 * 2) x₀ q.2.1 q.2.2)) := rfl

theorem cJQ_zeroTop (x₀ α : M) (q : ℝ × M × JetModel E 0 0 (2 * 2)) :
    cJQ g zeroTop x₀ α q =
      (chart E α q.2.1, cJet g α q.2.1 (zeroTopC (jetOfModel 0 0 (2 * 2) x₀ q.2.1 q.2.2))) := by
  rw [cJQ, cJet, rdJet_zeroTopC]
  rfl

/-! ### Openness of the domain -/

/-- The good points of the chart at `α`, in the local coordinates at `x₀`. -/
def goodQ (x₀ α : M) : Set (ℝ × M × JetModel E 0 0 (2 * 2)) :=
  {q | q.2.1 ∈ jetBase (I := 𝓘(ℝ, E)) 0 0 (2 * 2) x₀ ∧ q.2.1 ∈ patchU E α ∧
    good g J F α q.2.1 (jetOfModel 0 0 (2 * 2) x₀ q.2.1 q.2.2)}

variable {g J F}
variable (hJi : ∀ u v, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v, J (J v) = -v) (hF : IsLagImmersion E J F)
include hJi hJJ hF

theorem isOpen_goodQ (x₀ α : M) : IsOpen (goodQ g J F x₀ α) := by
  rw [isOpen_iff_mem_nhds]
  rintro q ⟨h₀, hα, hg⟩
  have hc : ContinuousAt (fun q : ℝ × M × JetModel E 0 0 (2 * 2) =>
      ((cJQ g id x₀ α q).1, (cJQ g id x₀ α q).2.1)) q := by
    have h := (contMDiffAt_cJQ g contDiff_id x₀ α h₀ hα).continuousAt
    exact (continuous_fst.prodMk (continuous_fst.comp continuous_snd)).continuousAt.comp h
  have hx : ContinuousAt (fun q : ℝ × M × JetModel E 0 0 (2 * 2) => q.2.1) q :=
    (continuous_fst.comp continuous_snd).continuousAt
  filter_upwards [hx.preimage_mem_nhds ((isOpen_jetBase (I := 𝓘(ℝ, E)) (s := 0) (2 * 2) x₀).mem_nhds h₀),
    hx.preimage_mem_nhds ((isOpen_patchU α).mem_nhds hα),
    hc.preimage_mem_nhds ((isOpen_chartDom hJi hJJ hF α).mem_nhds hg)] with q' h1 h2 h3
  exact ⟨h1, h2, h3⟩

/-- **The domain is relatively open** (`MMData.dom_relOpen`). -/
theorem mmDom_relOpen (T₀ : ℝ) (x₀ : M) : ∃ U : Set (ℝ × M × JetModel E 0 0 (2 * 2)), IsOpen U ∧
    localDom 0 0 (2 * 2) T₀ (mmDom g J F) x₀ = U ∩ localSlab (I := 𝓘(ℝ, E)) 0 0 (2 * 2) T₀ x₀ := by
  refine ⟨⋂ i : PIdx E M, ({q | q.2.1 ∉ tsupport (patchρ E M i)} ∪ goodQ g J F x₀ i.1),
    isOpen_iInter_of_finite fun i => IsOpen.union ?_ (isOpen_goodQ hJi hJJ hF x₀ i.1), ?_⟩
  · exact (isClosed_tsupport _).isOpen_compl.preimage (continuous_fst.comp continuous_snd)
  · ext q
    simp only [localDom, mem_ofPred_eq, mem_inter_iff, mem_iInter, mem_union]
    constructor
    · rintro ⟨hs, hd⟩
      refine ⟨fun i => ?_, hs⟩
      by_cases hx : q.2.1 ∈ tsupport (patchρ E M i)
      · exact Or.inr ⟨hs.2, tsupport_patchρ_subset i hx, hd i hx⟩
      · exact Or.inl hx
    · rintro ⟨hU, hs⟩
      refine ⟨hs, fun i hx => ?_⟩
      rcases hU i with h | h
      · exact absurd hx h
      · exact h.2.2

omit hJi hJJ hF in
/-- A function which is smooth near the support of a bump function, times the bump function. -/
theorem contMDiffAt_bump_mul (i : PIdx E M) {f : ℝ × M × JetModel E 0 0 (2 * 2) → ℝ}
    {q : ℝ × M × JetModel E 0 0 (2 * 2)}
    (hf : q.2.1 ∈ tsupport (patchρ E M i) → ContMDiffAt 𝓠 𝓘(ℝ, ℝ) ∞ f q) :
    ContMDiffAt 𝓠 𝓘(ℝ, ℝ) ∞ (fun q => patchρ E M i q.2.1 * f q) q := by
  by_cases hx : q.2.1 ∈ tsupport (patchρ E M i)
  · exact ((patchρ E M i).contMDiff.contMDiffAt.comp q (contMDiffAt_snd_fst q)).mul (hf hx)
  · have hopen : IsOpen {q : ℝ × M × JetModel E 0 0 (2 * 2) | q.2.1 ∉ tsupport (patchρ E M i)} :=
      (isClosed_tsupport _).isOpen_compl.preimage (continuous_fst.comp continuous_snd)
    have h0 : (fun q : ℝ × M × JetModel E 0 0 (2 * 2) => patchρ E M i q.2.1 * f q) =ᶠ[𝓝 q]
        fun _ => (0 : ℝ) := by
      filter_upwards [hopen.mem_nhds hx] with q' hq'
      rw [image_eq_zero_of_notMem_tsupport hq', zero_mul]
    exact contMDiffAt_const.congr_of_eventuallyEq h0

/-! ### Smoothness of the principal coefficient -/

omit hJi hJJ hF in
/-- A coordinate covector of `x₀`, read in the chart at `α`. -/
theorem chartCovec_coordCoframe_eq (x₀ α : M) {x : M}
    (h₀ : x ∈ (trivializationAt E (TangentSpace 𝓘(ℝ, E)) x₀).baseSet)
    (hα : x ∈ (trivializationAt E (TangentSpace 𝓘(ℝ, E)) α).baseSet)
    (a : Fin (Module.finrank ℝ E)) :
    chartCovec α x (coordCoframe (I := 𝓘(ℝ, E)) x₀ x a) =
      (LinearMap.toContinuousLinearMap ((Module.finBasis ℝ E).coord a)).comp
        ((trivializationAt E (TangentSpace 𝓘(ℝ, E)) α).coordChangeL ℝ
          (trivializationAt E (TangentSpace 𝓘(ℝ, E)) x₀) x : E →L[ℝ] E) := by
  ext v
  simp only [chartCovec, coordCoframe, ContinuousLinearMap.comp_apply, ContinuousLinearEquiv.coe_coe]
  rw [Trivialization.coordChangeL_apply _ _ ⟨hα, h₀⟩, Trivialization.continuousLinearMapAt_apply,
    Trivialization.coe_linearMapAt_of_mem _ h₀, Trivialization.symmL_apply _ hα]

omit hJi hJJ hF in
theorem contMDiffAt_chartCovec (x₀ α : M) {q : ℝ × M × JetModel E 0 0 (2 * 2)}
    (h₀ : q.2.1 ∈ (trivializationAt E (TangentSpace 𝓘(ℝ, E)) x₀).baseSet)
    (hα : q.2.1 ∈ (trivializationAt E (TangentSpace 𝓘(ℝ, E)) α).baseSet)
    (a : Fin (Module.finrank ℝ E)) :
    ContMDiffAt 𝓠 𝓘(ℝ, E →L[ℝ] ℝ) ∞ (fun q : ℝ × M × JetModel E 0 0 (2 * 2) =>
      chartCovec α q.2.1 (coordCoframe (I := 𝓘(ℝ, E)) x₀ q.2.1 a)) q := by
  have h1 : ContMDiffAt 𝓘(ℝ, E) 𝓘(ℝ, E →L[ℝ] E) ∞ (fun x : M =>
      ((trivializationAt E (TangentSpace 𝓘(ℝ, E)) α).coordChangeL ℝ
        (trivializationAt E (TangentSpace 𝓘(ℝ, E)) x₀) x : E →L[ℝ] E)) q.2.1 :=
    contMDiffAt_coordChangeL (IB := 𝓘(ℝ, E)) (n := ∞) hα h₀
  have h2 : ContDiff ℝ ∞ (fun A : E →L[ℝ] E =>
      (LinearMap.toContinuousLinearMap ((Module.finBasis ℝ E).coord a)).comp A) :=
    (ContinuousLinearMap.compL ℝ E E ℝ
      (LinearMap.toContinuousLinearMap ((Module.finBasis ℝ E).coord a))).contDiff
  have h3 := h2.contDiffAt.comp_contMDiffAt (h1.comp q (contMDiffAt_snd_fst q))
  have hopen : IsOpen {q : ℝ × M × JetModel E 0 0 (2 * 2) |
      q.2.1 ∈ (trivializationAt E (TangentSpace 𝓘(ℝ, E)) x₀).baseSet ∩
        (trivializationAt E (TangentSpace 𝓘(ℝ, E)) α).baseSet} :=
    ((trivializationAt E (TangentSpace 𝓘(ℝ, E)) x₀).open_baseSet.inter
      (trivializationAt E (TangentSpace 𝓘(ℝ, E)) α).open_baseSet).preimage
      (continuous_fst.comp continuous_snd)
  refine h3.congr_of_eventuallyEq ?_
  filter_upwards [hopen.mem_nhds ⟨h₀, hα⟩] with q' hq'
  exact chartCovec_coordCoframe_eq x₀ α hq'.1 hq'.2 a

/-- The differential of the graph is smooth in MM's local coordinates, at good points. -/
theorem contMDiffAt_Li (x₀ α : M) {q : ℝ × M × JetModel E 0 0 (2 * 2)} (hq : q ∈ goodQ g J F x₀ α) :
    ContMDiffAt 𝓠 𝓘(ℝ, E →L[ℝ] V) ∞ (fun q : ℝ × M × JetModel E 0 0 (2 * 2) =>
      Li g J F α q.2.1 (jetOfModel 0 0 (2 * 2) x₀ q.2.1 q.2.2)) q := by
  have hL : ContDiffAt ℝ ∞ (fun p : E × SJet E => LOf (chartPsi E J F α) p.1 p.2) (cJQ g id x₀ α q) :=
    contDiffAt_LOf (isOpen_chartDom hJi hJJ hF α) (contDiffOn_chartPsi hJi hJJ hF α) hq.2.2
  exact hL.comp_contMDiffAt (contMDiffAt_cJQ g contDiff_id x₀ α hq.1 hq.2.1)

/-- An entry of the cometric of one chart is smooth at good points. -/
theorem contMDiffAt_cometricAt (x₀ α : M) (a b : Fin (Module.finrank ℝ E))
    {q : ℝ × M × JetModel E 0 0 (2 * 2)} (hq : q ∈ goodQ g J F x₀ α) :
    ContMDiffAt 𝓠 𝓘(ℝ, ℝ) ∞ (fun q : ℝ × M × JetModel E 0 0 (2 * 2) =>
      cometricAt α q.2.1 (Li g J F α q.2.1 (jetOfModel 0 0 (2 * 2) x₀ q.2.1 q.2.2))
        (coordCoframe (I := 𝓘(ℝ, E)) x₀ q.2.1 a) (coordCoframe (I := 𝓘(ℝ, E)) x₀ q.2.1 b)) q := by
  have hs : q.2.1 ∈ (chartAt E α).source := by
    have := patchU_subset_source α hq.2.1
    rwa [extChartAt_source] at this
  have h₀ : q.2.1 ∈ (trivializationAt E (TangentSpace 𝓘(ℝ, E)) x₀).baseSet :=
    (mem_jetBase_iff (I := 𝓘(ℝ, E))).1 hq.1
  have hα := mem_tangent_baseSet hs
  have hG : ContMDiffAt 𝓠 𝓘(ℝ, (E →L[ℝ] ℝ) →L[ℝ] E) ∞ (fun q : ℝ × M × JetModel E 0 0 (2 * 2) =>
      ginvOf (Li g J F α q.2.1 (jetOfModel 0 0 (2 * 2) x₀ q.2.1 q.2.2))) q :=
    (contDiffAt_ginvOf (injective_Li hJi hJJ hF hq.2.2)).comp_contMDiffAt
      (f := fun q : ℝ × M × JetModel E 0 0 (2 * 2) =>
        Li g J F α q.2.1 (jetOfModel 0 0 (2 * 2) x₀ q.2.1 q.2.2))
      (contMDiffAt_Li hJi hJJ hF x₀ α hq)
  have ha := contMDiffAt_chartCovec x₀ α h₀ hα a
  have hb := contMDiffAt_chartCovec x₀ α h₀ hα b
  have h1 := ha.clm_apply (hG.clm_apply hb)
  have h2 := hb.clm_apply (hG.clm_apply ha)
  have h3 : ContMDiffAt 𝓠 𝓘(ℝ, ℝ) ∞ (fun _ : ℝ × M × JetModel E 0 0 (2 * 2) => (1 / 2 : ℝ)) q :=
    contMDiffAt_const
  exact h3.mul (h1.add h2)

/-- **The principal coefficient is smooth** (`MMData.Ef_smooth`). -/
theorem contMDiffOn_mmEf (T₀ : ℝ) (x₀ : M) (a b : Fin (Module.finrank ℝ E)) :
    ContMDiffOn 𝓠 𝓘(ℝ, ℝ) ∞
      (fun q : ℝ × M × JetModel E 0 0 (2 * 2) =>
        mmEf g J F q.1 q.2.1 (jetOfModel 0 0 (2 * 2) x₀ q.2.1 q.2.2)
          (coordCoframe (I := 𝓘(ℝ, E)) x₀ q.2.1 a) (coordCoframe (I := 𝓘(ℝ, E)) x₀ q.2.1 b))
      (localDom 0 0 (2 * 2) T₀ (mmDom g J F) x₀) := by
  intro q hq
  refine ContMDiffAt.contMDiffWithinAt ?_
  simp only [mmEf_apply]
  refine contMDiffAt_finsetSum fun i _ => contMDiffAt_bump_mul i fun hx => ?_
  exact contMDiffAt_cometricAt hJi hJJ hF x₀ i.1 a b ⟨hq.1.2, tsupport_patchρ_subset i hx, hq.2 i hx⟩

/-! ### Smoothness of the lower-order term -/

/-- The chart operator at the jet with top component `0` is smooth at good points. -/
theorem contMDiffAt_Fi (x₀ α : M) {q : ℝ × M × JetModel E 0 0 (2 * 2)} (hq : q ∈ goodQ g J F x₀ α) :
    ContMDiffAt 𝓠 𝓘(ℝ, ℝ) ∞ (fun q : ℝ × M × JetModel E 0 0 (2 * 2) =>
      Fi g J F α q.2.1 (zeroTopC (jetOfModel 0 0 (2 * 2) x₀ q.2.1 q.2.2))) q := by
  have hgood : good g J F α q.2.1 (zeroTopC (jetOfModel 0 0 (2 * 2) x₀ q.2.1 q.2.2)) :=
    (good_zeroTopC g J F α q.2.1 _).2 hq.2.2
  have hpt : ((cJQ g zeroTop x₀ α q).1, (permTop (cJQ g zeroTop x₀ α q).2).1) ∈ chartDom E J F α := by
    rw [cJQ_zeroTop]
    exact hgood
  have hP : ContDiffAt ℝ ∞ (Pjet J (chartPsi E J F α))
      ((cJQ g zeroTop x₀ α q).1, permTop (cJQ g zeroTop x₀ α q).2) :=
    contDiffAt_Pjet (isOpen_chartDom hJi hJJ hF α) (contDiffOn_chartPsi hJi hJJ hF α) hpt
      (injective_LOf (chartPsi_symplectic hJi hJJ hF α) hpt)
  have hperm : ContDiff ℝ ∞ (fun p : E × SJet E => (p.1, permTop p.2)) :=
    contDiff_fst.prodMk (contDiff_permTop.comp contDiff_snd)
  have h := hP.comp_contMDiffAt (f := fun q => ((cJQ g zeroTop x₀ α q).1, permTop (cJQ g zeroTop x₀ α q).2))
    (hperm.contDiffAt.comp_contMDiffAt (contMDiffAt_cJQ g contDiff_zeroTop x₀ α hq.1 hq.2.1))
  refine h.congr_of_eventuallyEq (Eventually.of_forall fun q' => ?_)
  show Fi g J F α q'.2.1 (zeroTopC (jetOfModel 0 0 (2 * 2) x₀ q'.2.1 q'.2.2)) =
    Pjet J (chartPsi E J F α) ((cJQ g zeroTop x₀ α q').1, permTop (cJQ g zeroTop x₀ α q').2)
  rw [cJQ_zeroTop]
  rfl

/-- The scalar lower-order term is smooth on the domain. -/
theorem contMDiffAt_mmbS (T₀ : ℝ) (x₀ : M) {q : ℝ × M × JetModel E 0 0 (2 * 2)}
    (hq : q ∈ localDom 0 0 (2 * 2) T₀ (mmDom g J F) x₀) :
    ContMDiffAt 𝓠 𝓘(ℝ, ℝ) ∞ (fun q : ℝ × M × JetModel E 0 0 (2 * 2) =>
      mmbS g J F q.2.1 (jetOfModel 0 0 (2 * 2) x₀ q.2.1 q.2.2)) q := by
  simp only [mmbS]
  refine contMDiffAt_finsetSum fun i _ => contMDiffAt_bump_mul i fun hx => ?_
  exact contMDiffAt_Fi hJi hJJ hF x₀ i.1 ⟨hq.1.2, tsupport_patchρ_subset i hx, hq.2 i hx⟩

omit hJi hJJ hF in
/-- The local representative of `b` is its scalar value times the identity. -/
theorem localRep_mmb (x₀ : M) {q : ℝ × M × JetModel E 0 0 (2 * 2)}
    (hq : q.2.1 ∈ jetBase (I := 𝓘(ℝ, E)) 0 0 (2 * 2) x₀) :
    localRep 0 0 (2 * 2) (mmb g J F) x₀ q =
      mmbS g J F q.2.1 (jetOfModel 0 0 (2 * 2) x₀ q.2.1 q.2.2) •
        (ContinuousLinearMap.id ℝ (Tensor0SModel 0 ℝ E) : TensorRSModel 0 0 ℝ E) := by
  have hs : q.2.1 ∈ (extChartAt 𝓘(ℝ, E) x₀).source := by
    rw [extChartAt_source, ← TangentBundle.trivializationAt_baseSet (I := 𝓘(ℝ, E))]
    exact (mem_jetBase_iff (I := 𝓘(ℝ, E))).1 hq
  rw [eq_scal_smul_id (localRep 0 0 (2 * 2) (mmb g J F) x₀ q)]
  congr 1
  show scal E ((trivAt 𝓘(ℝ, E) M 0 0 x₀).continuousLinearMapAt ℝ q.2.1 (mmb g J F q.1 q.2.1 _)) = _
  rw [scal_trivAt hs, lam_mmb]

/-- **The lower-order term is smooth** (`MMData.b_smooth`). -/
theorem contMDiffOn_mmb (T₀ : ℝ) (x₀ : M) :
    ContMDiffOn 𝓠 𝓘(ℝ, TensorRSModel 0 0 ℝ E) ∞ (localRep 0 0 (2 * 2) (mmb g J F) x₀)
      (localDom 0 0 (2 * 2) T₀ (mmDom g J F) x₀) := by
  have h : ContMDiffOn 𝓠 𝓘(ℝ, TensorRSModel 0 0 ℝ E) ∞
      (fun q : ℝ × M × JetModel E 0 0 (2 * 2) =>
        mmbS g J F q.2.1 (jetOfModel 0 0 (2 * 2) x₀ q.2.1 q.2.2) •
          (ContinuousLinearMap.id ℝ (Tensor0SModel 0 ℝ E) : TensorRSModel 0 0 ℝ E))
      (localDom 0 0 (2 * 2) T₀ (mmDom g J F) x₀) := fun q hq =>
    ((contMDiffAt_mmbS hJi hJJ hF T₀ x₀ hq).smul contMDiffAt_const).contMDiffWithinAt
  exact h.congr fun q hq => localRep_mmb x₀ hq.1.2

end Bridge

end CWFlow
