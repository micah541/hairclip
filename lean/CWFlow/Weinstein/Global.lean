import CWFlow.Moser.Transfer
import CWFlow.Weinstein.Overlap

/-!
# The immersed Weinstein theorem on a manifold, as a compatible family of chart maps

Let `M` be a Hausdorff smooth manifold modelled on `E`, and `F : M → V` a Lagrangian immersion into `ℂⁿ`, given
chartwise (`IsLagImmersion`). We use the charts `chart z = extChartAt 𝓘(ℝ, E) z`. For each `z`, the canonical
local Weinstein map `chartPsi z` lives on the cotangent coordinates of `chart z`. A covector at `x` is written in
the coordinates of a chart; the coordinates of `chartAt x` match Mathlib's `TangentSpace` and `mfderiv` at `x`.

`weinstein_global`:
* each `chartPsi z` is a local Weinstein map for `ω₀`, defined on an open set containing the zero section of the
  chart;
* the family is compatible: `chartPsi w (cotLift τ q) = chartPsi z q` for `τ = transition E z w`, whenever both
  sides are defined.

Together these are the data of a symplectic map from a neighbourhood of the zero section of `T*M` to `ℂⁿ` that
extends `F`. Compatibility needs no assumption that the Moser curves stay in the overlap of the two charts:
`moserCurve_transfer`, with the closure condition supplied by `T2Space M`, shows that they do.
-/

open Set Function Filter

open scoped Manifold Topology ContDiff InnerProductSpace RealInnerProductSpace

noncomputable section

namespace CWFlow

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
variable {M : Type*} [TopologicalSpace M] [ChartedSpace E M]
variable {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]

variable (E) in
/-- The chart at `z`. -/
abbrev chart (z : M) : PartialEquiv M E := extChartAt 𝓘(ℝ, E) z

variable (E) in
/-- `F` in the chart at `z`. -/
def chartRep (F : M → V) (z : M) : E → V := F ∘ (chart E z).symm

variable (E) in
/-- The change of chart `τ_{z→w} = φ_w ∘ φ_z⁻¹`. -/
def transition (z w : M) : E → E := chart E w ∘ (chart E z).symm

variable (E) in
/-- The overlap `φ_z(U_z ∩ U_w)`, in the coordinates of `z`. -/
def overlap (z w : M) : Set E := (chart E z).target ∩ (chart E z).symm ⁻¹' (chart E w).source

variable (E) in
/-- A Lagrangian immersion, chartwise: smooth, injective derivative, `F*ω₀ = 0`. -/
structure IsLagImmersion (J : V →L[ℝ] V) (F : M → V) : Prop where
  smooth : ∀ z : M, ContDiffOn ℝ ∞ (chartRep E F z) (chart E z).target
  inj : ∀ z : M, ∀ y ∈ (chart E z).target, Function.Injective (fderiv ℝ (chartRep E F z) y)
  lag : ∀ z : M, ∀ y ∈ (chart E z).target, ∀ u w,
    ⟪J (fderiv ℝ (chartRep E F z) y u), fderiv ℝ (chartRep E F z) y w⟫ = 0

variable (E) (J : V →L[ℝ] V) (F : M → V)

/-- The primitive `β` in the chart at `z`. -/
def chartBeta (z : M) : E × (E →L[ℝ] ℝ) → (E × (E →L[ℝ] ℝ)) →L[ℝ] ℝ :=
  beta J (chartRep E F z) (liftA (chartRep E F z))

/-- The domain of the canonical Weinstein map in the chart at `z`. -/
def chartDom (z : M) : Set (E × (E →L[ℝ] ℝ)) :=
  moserDom lamCan (chartBeta E J F z) (Prod.fst ⁻¹' (chart E z).target)

/-- The canonical Weinstein map in the chart at `z`. -/
def chartPsi (z : M) : E × (E →L[ℝ] ℝ) → V :=
  weinsteinMap (Phi0 J (chartRep E F z) (liftA (chartRep E F z))) (chartBeta E J F z)
    (Prod.fst ⁻¹' (chart E z).target)

variable {E J F}

/-! ### Charts -/

theorem isOpen_chart_target (z : M) : IsOpen (chart E z).target :=
  isOpen_extChartAt_target z

theorem isOpen_overlap (z w : M) : IsOpen (overlap E z w) :=
  (continuousOn_extChartAt_symm z).isOpen_inter_preimage (isOpen_extChartAt_target z)
    (isOpen_extChartAt_source w)

theorem overlap_subset (z w : M) : overlap E z w ⊆ (chart E z).target := inter_subset_left

theorem mapsTo_transition (z w : M) : MapsTo (transition E z w) (overlap E z w) (chart E w).target :=
  fun _ hy => (chart E w).map_source hy.2

theorem chartRep_transition (F : M → V) {z w : M} {y : E} (hy : y ∈ overlap E z w) :
    chartRep E F z y = chartRep E F w (transition E z w y) := by
  simp only [chartRep, transition, Function.comp_apply]
  rw [(chart E w).left_inv hy.2]

theorem transition_transition {z w : M} {y : E} (hy : y ∈ overlap E z w) :
    transition E w z (transition E z w y) = y := by
  simp only [transition, Function.comp_apply]
  rw [(chart E w).left_inv hy.2, (chart E z).right_inv hy.1]

theorem transition_mem_overlap {z w : M} {y : E} (hy : y ∈ overlap E z w) :
    transition E z w y ∈ overlap E w z := by
  refine ⟨mapsTo_transition (E := E) z w hy, ?_⟩
  simp only [transition, Function.comp_apply, mem_preimage]
  rw [(chart E w).left_inv hy.2]
  exact (chart E z).map_target hy.1

variable [IsManifold 𝓘(ℝ, E) ∞ M]

theorem contDiffOn_transition (z w : M) : ContDiffOn ℝ ∞ (transition E z w) (overlap E z w) := by
  have h := contDiffOn_ext_coord_change (I := 𝓘(ℝ, E)) (n := ∞) w z
  have e : ((extChartAt 𝓘(ℝ, E) z).symm ≫ extChartAt 𝓘(ℝ, E) w).source = overlap E z w := by
    simp [overlap, PartialEquiv.trans_source]
  rw [e] at h
  exact h

theorem exists_equiv_fderiv_transition {z w : M} {y : E} (hy : y ∈ overlap E z w) :
    ∃ e : E ≃L[ℝ] E, (e : E →L[ℝ] E) = fderiv ℝ (transition E z w) y := by
  have hτ := ((contDiffOn_transition (E := E) z w).contDiffAt ((isOpen_overlap (E := E) z w).mem_nhds hy)).differentiableAt
    (by simp)
  have hy' := transition_mem_overlap hy
  have hτ' := ((contDiffOn_transition (E := E) w z).contDiffAt
    ((isOpen_overlap (E := E) w z).mem_nhds hy')).differentiableAt (by simp)
  have hid : transition E w z ∘ transition E z w =ᶠ[𝓝 y] id := by
    filter_upwards [(isOpen_overlap (E := E) z w).mem_nhds hy] with y' hy'' using transition_transition hy''
  have hcomp := fderiv_comp y hτ' hτ
  rw [hid.fderiv_eq, fderiv_id] at hcomp
  have hinj : Function.Injective (fderiv ℝ (transition E z w) y) := by
    intro a b hab
    have := congrArg (fderiv ℝ (transition E w z) (transition E z w y)) hab
    rw [← ContinuousLinearMap.comp_apply, ← ContinuousLinearMap.comp_apply, ← hcomp] at this
    simpa using this
  exact ⟨(LinearMap.linearEquivOfInjective (fderiv ℝ (transition E z w) y : E →ₗ[ℝ] E) hinj
    rfl).toContinuousLinearEquiv, ContinuousLinearMap.ext fun _ => rfl⟩

/-! ### The compatible family -/

variable [T2Space M]

/-- **Compatibility of the chart Weinstein maps.** -/
theorem chartPsi_compat (hF : IsLagImmersion E J F) {z w : M} {q : E × (E →L[ℝ] ℝ)}
    (hq : q ∈ chartDom E J F z) (hqN : q.1 ∈ overlap E z w)
    (hq' : cotLift (transition E z w) q ∈ chartDom E J F w) :
    chartPsi E J F w (cotLift (transition E z w) q) = chartPsi E J F z q := by
  obtain ⟨γ₁, h₁⟩ := hq
  obtain ⟨γ₂, h₂⟩ := hq'
  set τ := transition E z w with hτdef
  have hS₁ := isOpen_chart_target (E := E) (M := M) z
  have hS₂ := isOpen_chart_target (E := E) (M := M) w
  have hS₁₂ := isOpen_overlap (E := E) z w
  have hτ : ContDiffOn ℝ ∞ τ (overlap E z w) := contDiffOn_transition (E := E) z w
  have hτA : ∀ y ∈ overlap E z w, ContDiffAt ℝ ∞ τ y := fun y hy => hτ.contDiffAt (hS₁₂.mem_nhds hy)
  have hτinv : ∀ y ∈ overlap E z w, ∃ e : E ≃L[ℝ] E, (e : E →L[ℝ] E) = fderiv ℝ τ y :=
    fun y hy => exists_equiv_fderiv_transition hy
  have hFτ : ∀ y ∈ overlap E z w, chartRep E F z y = chartRep E F w (τ y) := fun y hy =>
    chartRep_transition F hy
  have hF₂A : ∀ y ∈ overlap E z w, ContDiffAt ℝ ∞ (chartRep E F w) (τ y) := fun y hy =>
    (hF.smooth w).contDiffAt (hS₂.mem_nhds (mapsTo_transition (E := E) z w hy))
  have hinj₁' : ∀ y ∈ overlap E z w, Function.Injective (fderiv ℝ (chartRep E F z) y) := fun y hy =>
    hF.inj z y hy.1
  have hinj₂' : ∀ y ∈ overlap E z w, Function.Injective (fderiv ℝ (chartRep E F w) (τ y)) :=
    fun y hy => hF.inj w _ (mapsTo_transition (E := E) z w hy)
  have hO₂ : IsOpen (Prod.fst ⁻¹' (chart E w).target : Set (E × (E →L[ℝ] ℝ))) :=
    hS₂.preimage continuous_fst
  have hN : IsOpen (Prod.fst ⁻¹' overlap E z w : Set (E × (E →L[ℝ] ℝ))) :=
    hS₁₂.preimage continuous_fst
  have hTN : MapsTo (cotLift τ) (Prod.fst ⁻¹' overlap E z w) (Prod.fst ⁻¹' (chart E w).target) :=
    fun p hp => mapsTo_transition (E := E) z w hp
  have hT : ContDiffOn ℝ ∞ (cotLift τ) (Prod.fst ⁻¹' overlap E z w) := fun p hp => by
    obtain ⟨e, he⟩ := hτinv p.1 hp
    exact (contDiffAt_cotLift (hτA p.1 hp) he).contDiffWithinAt
  have hTinv : ∀ p ∈ (Prod.fst ⁻¹' overlap E z w : Set (E × (E →L[ℝ] ℝ))),
      ∃ e : (E × (E →L[ℝ] ℝ)) ≃L[ℝ] (E × (E →L[ℝ] ℝ)),
        (e : _ →L[ℝ] _) = fderiv ℝ (cotLift τ) p := fun p hp => by
    obtain ⟨e, he⟩ := hτinv p.1 hp
    exact exists_equiv_fderiv_cotLift (hτA p.1 hp) he
  have hlamT : ∀ p ∈ (Prod.fst ⁻¹' overlap E z w : Set (E × (E →L[ℝ] ℝ))),
      lamCan p = pb1form (cotLift τ) lamCan p := fun p hp => by
    obtain ⟨e, he⟩ := hτinv p.1 hp
    exact lamCan_natural (hτA p.1 hp) he
  have hβT : ∀ p ∈ (Prod.fst ⁻¹' overlap E z w : Set (E × (E →L[ℝ] ℝ))),
      chartBeta E J F z p = pb1form (cotLift τ) (chartBeta E J F w) p := fun p hp =>
    beta_natural hS₁₂ hτA hτinv hFτ hF₂A hinj₁' hinj₂' hp
  -- the closure condition, from the Hausdorff property
  have hcl : ∀ t ∈ Ioc (0 : ℝ) 1, (∀ r ∈ Ico (0 : ℝ) t,
      γ₁ r ∈ (Prod.fst ⁻¹' overlap E z w : Set (E × (E →L[ℝ] ℝ))) ∧ γ₂ r = cotLift τ (γ₁ r)) →
      γ₁ t ∈ (Prod.fst ⁻¹' overlap E z w : Set (E × (E →L[ℝ] ℝ))) := by
    intro t ht hbefore
    have htI : t ∈ Icc (0 : ℝ) 1 := Ioc_subset_Icc_self ht
    have h1t : (γ₁ t).1 ∈ (chart E z).target := h₁.mem htI
    have h2t : (γ₂ t).1 ∈ (chart E w).target := h₂.mem htI
    have : (𝓝[Ico 0 t] t).NeBot := mem_closure_iff_nhdsWithin_neBot.1
      (by rw [closure_Ico ht.1.ne]; exact right_mem_Icc.2 ht.1.le)
    have hsub : Ico (0 : ℝ) t ⊆ Icc 0 1 := fun r hr => ⟨hr.1, hr.2.le.trans ht.2⟩
    have hx₁ : ContinuousOn (fun r => (γ₁ r).1) (Icc 0 1) :=
      continuous_fst.comp_continuousOn h₁.continuousOn
    have hx₂ : ContinuousOn (fun r => (γ₂ r).1) (Icc 0 1) :=
      continuous_fst.comp_continuousOn h₂.continuousOn
    have hm₁ : MapsTo (fun r => (γ₁ r).1) (Icc 0 1) (chart E z).target := fun r hr => h₁.mem hr
    have hm₂ : MapsTo (fun r => (γ₂ r).1) (Icc 0 1) (chart E w).target := fun r hr => h₂.mem hr
    have hA : Tendsto (fun r => (chart E z).symm (γ₁ r).1) (𝓝[Ico 0 t] t)
        (𝓝 ((chart E z).symm (γ₁ t).1)) :=
      ((ContinuousWithinAt.comp (f := fun r => (γ₁ r).1) (x := t)
        (continuousOn_extChartAt_symm (I := 𝓘(ℝ, E)) z _ h1t) (hx₁ t htI) hm₁).mono
        hsub).tendsto
    have hB : Tendsto (fun r => (chart E w).symm (γ₂ r).1) (𝓝[Ico 0 t] t)
        (𝓝 ((chart E w).symm (γ₂ t).1)) :=
      ((ContinuousWithinAt.comp (f := fun r => (γ₂ r).1) (x := t)
        (continuousOn_extChartAt_symm (I := 𝓘(ℝ, E)) w _ h2t) (hx₂ t htI) hm₂).mono
        hsub).tendsto
    have heq : (chart E z).symm (γ₁ t).1 = (chart E w).symm (γ₂ t).1 := by
      refine tendsto_nhds_unique_of_eventuallyEq hA hB
        (eventually_nhdsWithin_of_forall fun r hr => ?_)
      obtain ⟨hrN, hr2⟩ := hbefore r hr
      simp only
      rw [hr2, cotLift_fst, hτdef, transition, Function.comp_apply, (chart E w).left_inv hrN.2]
    refine ⟨h1t, ?_⟩
    simp only [mem_preimage]
    rw [heq]
    exact (chart E w).map_target h2t
  have htr := moserCurve_transfer hO₂ contDiff_lamCan.contDiffOn
    (contDiffOn_beta_liftA hS₂ (hF.smooth w) (hF.inj w)) hN hTN hT hTinv hlamT hβT hqN h₁ h₂ hcl
  exact weinsteinMap_natural hS₁ hS₂ hS₁₂ (overlap_subset z w) (mapsTo_transition (E := E) z w)
    (hF.smooth z) (hF.smooth w) (hF.inj z) (hF.inj w) hτ hτinv hFτ h₁ fun t ht => (htr t ht).1

/-- **The immersed Weinstein theorem on a manifold.** For a (not necessarily injective) Lagrangian immersion
`F : M → ℂⁿ`, the canonical chart maps form a compatible family of local Weinstein maps, each defined on an open
neighbourhood of the zero section of its chart. -/
theorem weinstein_global (hJi : ∀ u v, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v, J (J v) = -v)
    (hF : IsLagImmersion E J F) :
    (∀ z : M, IsLocalWeinstein (omega0 J) (chartRep E F z) (chartDom E J F z) (chartPsi E J F z) ∧
      ∀ y ∈ (chart E z).target, (y, (0 : E →L[ℝ] ℝ)) ∈ chartDom E J F z) ∧
    ∀ (z w : M) (q : E × (E →L[ℝ] ℝ)), q ∈ chartDom E J F z → q.1 ∈ overlap E z w →
      cotLift (transition E z w) q ∈ chartDom E J F w →
      chartPsi E J F w (cotLift (transition E z w) q) = chartPsi E J F z q :=
  ⟨fun z => isLocalWeinstein_Cn hJi hJJ (isOpen_chart_target z) (hF.smooth z) (hF.inj z) (hF.lag z),
    fun _ _ _ hq hqN hq' => chartPsi_compat hF hq hqN hq'⟩

end CWFlow
