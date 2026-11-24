import CWFlow.Extension.ImmReparam
import Mathlib.Geometry.Manifold.PartitionOfUnity
import Mathlib.Geometry.Manifold.ContMDiff.Basic
import Mathlib.Geometry.Manifold.ContMDiff.NormedSpace
import Mathlib.Geometry.Manifold.ContMDiff.Atlas
import Mathlib.MeasureTheory.Function.Jacobian
import Mathlib.MeasureTheory.Integral.Bochner.Set

/-!
# Integration of natural densities over a compact manifold

A fixed finite atlas with a smooth partition of unity (`IntAtlas`) defines the integral `intM` of a
natural density (`IsNatDens`) of a Lagrangian immersion `F` of `M` and a function `h`. The transfer lemma
`intM_eq_param` computes this integral in any parametrization of an open subset of `M` whose image
contains the support of `h`.
-/

open scoped ContDiff Manifold Topology RealInnerProductSpace
open MeasureTheory

set_option linter.unusedSectionVars false

namespace CWFlow.Extension

variable {n : ℕ} {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
  {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
  {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
  (ι : EuclideanSpace ℝ (Fin n) ≃L[ℝ] E)

/-- The inverse chart at `z`, read on `ℝⁿ`. -/
noncomputable def chartInv (z : M) (y : EuclideanSpace ℝ (Fin n)) : M :=
  (extChartAt 𝓘(ℝ, E) z).symm (ι y)

/-- The chart at `z`, read in `ℝⁿ`. -/
noncomputable def chartFwd (z : M) (p : M) : EuclideanSpace ℝ (Fin n) :=
  ι.symm (extChartAt 𝓘(ℝ, E) z p)

/-- The domain of the chart at `z` in `ℝⁿ`. -/
def chartDomE (z : M) : Set (EuclideanSpace ℝ (Fin n)) := ι ⁻¹' (extChartAt 𝓘(ℝ, E) z).target

/-! ### Part 0: determinant of the derivative -/

theorem im_det_eq_jacM (τ : EuclideanSpace ℝ (Fin n) → EuclideanSpace ℝ (Fin n))
    (w : EuclideanSpace ℝ (Fin n)) : (fderiv ℝ τ w).det = (jacM τ w).det := by
  classical
  rw [ContinuousLinearMap.det, ← LinearMap.det_toMatrix (EuclideanSpace.basisFun (Fin n) ℝ).toBasis]
  congr 1

theorem im_det_ne_zero_of_left_inv {τ τ' : EuclideanSpace ℝ (Fin n) → EuclideanSpace ℝ (Fin n)}
    {w : EuclideanSpace ℝ (Fin n)} (hτ : DifferentiableAt ℝ τ w) (hτ' : DifferentiableAt ℝ τ' (τ w))
    (h : τ' ∘ τ =ᶠ[𝓝 w] id) : (jacM τ w).det ≠ 0 := by
  have h1 : fderiv ℝ (τ' ∘ τ) w = (fderiv ℝ τ' (τ w)).comp (fderiv ℝ τ w) := fderiv_comp w hτ' hτ
  have h2 : fderiv ℝ (τ' ∘ τ) w = ContinuousLinearMap.id ℝ _ := by
    rw [h.fderiv_eq]; exact fderiv_id
  have h3 : (fderiv ℝ τ' (τ w)).det * (fderiv ℝ τ w).det = 1 := by
    simp only [ContinuousLinearMap.det, ← LinearMap.det_comp]
    have := congrArg (fun f : EuclideanSpace ℝ (Fin n) →L[ℝ] EuclideanSpace ℝ (Fin n) =>
      LinearMap.det (f : EuclideanSpace ℝ (Fin n) →ₗ[ℝ] EuclideanSpace ℝ (Fin n))) (h1.symm.trans h2)
    simpa using this
  rw [← im_det_eq_jacM]
  exact right_ne_zero_of_mul_eq_one h3

/-! ### Part 1: charts -/

theorem im_isOpen_chartDomE (z : M) : IsOpen (chartDomE ι z) :=
  (isOpen_extChartAt_target z).preimage ι.continuous

theorem im_chartInv_mem_source {z : M} {y : EuclideanSpace ℝ (Fin n)} (hy : y ∈ chartDomE ι z) :
    chartInv ι z y ∈ (extChartAt 𝓘(ℝ, E) z).source :=
  (extChartAt 𝓘(ℝ, E) z).map_target hy

theorem im_chartFwd_chartInv {z : M} {y : EuclideanSpace ℝ (Fin n)} (hy : y ∈ chartDomE ι z) :
    chartFwd ι z (chartInv ι z y) = y := by
  simp only [chartFwd, chartInv]
  rw [(extChartAt 𝓘(ℝ, E) z).right_inv hy]; simp

theorem im_chartInv_chartFwd {z p : M} (hp : p ∈ (extChartAt 𝓘(ℝ, E) z).source) :
    chartInv ι z (chartFwd ι z p) = p := by
  simp only [chartFwd, chartInv]
  rw [ContinuousLinearEquiv.apply_symm_apply]; exact (extChartAt 𝓘(ℝ, E) z).left_inv hp

theorem im_chartFwd_mem {z p : M} (hp : p ∈ (extChartAt 𝓘(ℝ, E) z).source) :
    chartFwd ι z p ∈ chartDomE ι z := by
  simpa [chartFwd, chartDomE] using (extChartAt 𝓘(ℝ, E) z).map_source hp

theorem im_continuousOn_chartInv (z : M) : ContinuousOn (chartInv ι z) (chartDomE ι z) :=
  (continuousOn_extChartAt_symm z).comp ι.continuous.continuousOn fun _ hy => hy

theorem im_continuousOn_chartFwd (z : M) :
    ContinuousOn (chartFwd ι z) (extChartAt 𝓘(ℝ, E) z).source :=
  ι.symm.continuous.comp_continuousOn (continuousOn_extChartAt z)

/-- A parametrization of an open subset of `M` by an open subset `B` of `ℝⁿ`: injective, open image,
a homeomorphism onto its image, smooth transitions with nonzero Jacobian to all charts. -/
structure IsParam (B : Set (EuclideanSpace ℝ (Fin n))) (σ : EuclideanSpace ℝ (Fin n) → M) : Prop where
  isOpen : IsOpen B
  injOn : Set.InjOn σ B
  isOpen_image : IsOpen (σ '' B)
  cont : ContinuousOn σ B
  inv : ∃ κ : M → EuclideanSpace ℝ (Fin n), ContinuousOn κ (σ '' B) ∧ ∀ y ∈ B, κ (σ y) = y
  trans_smooth : ∀ z : M, ContDiffOn ℝ ∞ (chartFwd ι z ∘ σ)
    (B ∩ σ ⁻¹' (extChartAt 𝓘(ℝ, E) z).source)
  trans_det : ∀ z : M, ∀ y ∈ B ∩ σ ⁻¹' (extChartAt 𝓘(ℝ, E) z).source,
    (jacM (chartFwd ι z ∘ σ) y).det ≠ 0

theorem im_isOpen_W {B : Set (EuclideanSpace ℝ (Fin n))} {σ : EuclideanSpace ℝ (Fin n) → M}
    (hσ : IsParam ι B σ) (z : M) : IsOpen (B ∩ σ ⁻¹' (extChartAt 𝓘(ℝ, E) z).source) :=
  hσ.cont.isOpen_inter_preimage hσ.isOpen (isOpen_extChartAt_source z)

theorem im_contDiffOn_trans (z z' : M) :
    ContDiffOn ℝ ∞ (chartFwd ι z' ∘ chartInv ι z)
      (chartDomE ι z ∩ chartInv ι z ⁻¹' (extChartAt 𝓘(ℝ, E) z').source) := by
  have h1 := contDiffOn_ext_coord_change (I := 𝓘(ℝ, E)) (n := ∞) z' z
  have h2 : ContDiffOn ℝ ∞ ((extChartAt 𝓘(ℝ, E) z' ∘ (extChartAt 𝓘(ℝ, E) z).symm) ∘ ι)
      (chartDomE ι z ∩ chartInv ι z ⁻¹' (extChartAt 𝓘(ℝ, E) z').source) := by
    refine h1.comp ι.contDiff.contDiffOn ?_
    intro y hy
    simp only [PartialEquiv.trans_source, PartialEquiv.symm_source, Set.mem_inter_iff,
      Set.mem_preimage]
    exact ⟨hy.1, hy.2⟩
  exact ι.symm.contDiff.comp_contDiffOn h2

theorem isParam_chart (z : M) : IsParam ι (chartDomE ι z) (chartInv ι z) where
  isOpen := im_isOpen_chartDomE ι z
  injOn := by
    intro a ha b hb hab
    exact ι.injective ((extChartAt 𝓘(ℝ, E) z).symm.injOn ha hb hab)
  isOpen_image := by
    have : chartInv ι z '' chartDomE ι z = (extChartAt 𝓘(ℝ, E) z).source := by
      ext p
      constructor
      · rintro ⟨y, hy, rfl⟩; exact im_chartInv_mem_source ι hy
      · intro hp
        exact ⟨chartFwd ι z p, im_chartFwd_mem ι hp, im_chartInv_chartFwd ι hp⟩
    rw [this]; exact isOpen_extChartAt_source z
  cont := im_continuousOn_chartInv ι z
  inv := by
    refine ⟨chartFwd ι z, ?_, fun y hy => im_chartFwd_chartInv ι hy⟩
    refine (im_continuousOn_chartFwd ι z).mono ?_
    rintro _ ⟨y, hy, rfl⟩; exact im_chartInv_mem_source ι hy
  trans_smooth z' := im_contDiffOn_trans ι z z'
  trans_det z' := by
    rintro y ⟨hy, hy'⟩
    have hW : IsOpen (chartDomE ι z ∩ chartInv ι z ⁻¹' (extChartAt 𝓘(ℝ, E) z').source) :=
      (im_continuousOn_chartInv ι z).isOpen_inter_preimage (im_isOpen_chartDomE ι z)
        (isOpen_extChartAt_source z')
    have hyW : y ∈ chartDomE ι z ∩ chartInv ι z ⁻¹' (extChartAt 𝓘(ℝ, E) z').source := ⟨hy, hy'⟩
    have hτ : DifferentiableAt ℝ (chartFwd ι z' ∘ chartInv ι z) y :=
      ((im_contDiffOn_trans ι z z').contDiffAt (hW.mem_nhds hyW)).differentiableAt (by simp)
    -- the inverse transition
    set p := chartInv ι z y with hp
    have hτy : chartFwd ι z' p ∈ chartDomE ι z' := im_chartFwd_mem ι hy'
    have hinv : chartInv ι z' (chartFwd ι z' p) = p := im_chartInv_chartFwd ι hy'
    have hy2 : chartInv ι z' (chartFwd ι z' p) ∈ (extChartAt 𝓘(ℝ, E) z).source := by
      rw [hinv]; exact im_chartInv_mem_source ι hy
    have hW' : IsOpen (chartDomE ι z' ∩ chartInv ι z' ⁻¹' (extChartAt 𝓘(ℝ, E) z).source) :=
      (im_continuousOn_chartInv ι z').isOpen_inter_preimage (im_isOpen_chartDomE ι z')
        (isOpen_extChartAt_source z)
    have hτ' : DifferentiableAt ℝ (chartFwd ι z ∘ chartInv ι z') ((chartFwd ι z' ∘ chartInv ι z) y) :=
      ((im_contDiffOn_trans ι z' z).contDiffAt
        (hW'.mem_nhds ⟨hτy, hy2⟩)).differentiableAt (by simp)
    refine im_det_ne_zero_of_left_inv hτ hτ' ?_
    filter_upwards [hW.mem_nhds hyW] with w hw
    simp only [Function.comp_apply, id]
    rw [im_chartInv_chartFwd ι hw.2, im_chartFwd_chartInv ι hw.1]

/-! ### Part 2: Lagrangian immersions of `M` and parametrizations -/

/-- A Lagrangian immersion of `M`, read in the charts. -/
structure IsLagImm (J : V →L[ℝ] V) (F : M → V) : Prop where
  smooth : ∀ z : M, ContDiffOn ℝ ∞ (F ∘ chartInv ι z) (chartDomE ι z)
  inj : ∀ z : M, ∀ y ∈ chartDomE ι z, Function.Injective (fderiv ℝ (F ∘ chartInv ι z) y)
  lag : ∀ z : M, ∀ y ∈ chartDomE ι z, ∀ u w,
    ⟪J (fderiv ℝ (F ∘ chartInv ι z) y u), fderiv ℝ (F ∘ chartInv ι z) y w⟫ = 0

theorem im_mem_W {B : Set (EuclideanSpace ℝ (Fin n))} {σ : EuclideanSpace ℝ (Fin n) → M}
    {y : EuclideanSpace ℝ (Fin n)} (hy : y ∈ B) :
    y ∈ B ∩ σ ⁻¹' (extChartAt 𝓘(ℝ, E) (σ y)).source :=
  ⟨hy, mem_extChartAt_source (σ y)⟩

theorem im_contDiffOn_comp_param {β : Type*} [NormedAddCommGroup β] [NormedSpace ℝ β] {G : M → β}
    (hG : ∀ z : M, ContDiffOn ℝ ∞ (G ∘ chartInv ι z) (chartDomE ι z))
    {B : Set (EuclideanSpace ℝ (Fin n))} {σ : EuclideanSpace ℝ (Fin n) → M} (hσ : IsParam ι B σ) :
    ContDiffOn ℝ ∞ (G ∘ σ) B := by
  intro y hy
  set z := σ y
  have hW := im_isOpen_W ι hσ z
  have h1 : ContDiffOn ℝ ∞ ((G ∘ chartInv ι z) ∘ (chartFwd ι z ∘ σ))
      (B ∩ σ ⁻¹' (extChartAt 𝓘(ℝ, E) z).source) :=
    (hG z).comp (hσ.trans_smooth z) fun w hw => im_chartFwd_mem ι hw.2
  have h2 : ContDiffOn ℝ ∞ (G ∘ σ) (B ∩ σ ⁻¹' (extChartAt 𝓘(ℝ, E) z).source) :=
    h1.congr fun w hw => by simp [Function.comp_def, im_chartInv_chartFwd ι hw.2]
  exact (h2.contDiffAt (hW.mem_nhds (im_mem_W (σ := σ) hy))).contDiffWithinAt

theorem im_fderiv_eq_comp {J : V →L[ℝ] V} {F : M → V} (hF : IsLagImm ι J F)
    {B : Set (EuclideanSpace ℝ (Fin n))} {σ : EuclideanSpace ℝ (Fin n) → M} (hσ : IsParam ι B σ)
    {y : EuclideanSpace ℝ (Fin n)} (hy : y ∈ B) :
    (Function.Injective (fderiv ℝ (F ∘ σ) y)) ∧
      ∀ u w, ⟪J (fderiv ℝ (F ∘ σ) y u), fderiv ℝ (F ∘ σ) y w⟫ = 0 := by
  set z := σ y
  have hW := im_isOpen_W ι hσ z
  have hyW := im_mem_W (E := E) (σ := σ) hy
  have hmaps : ∀ w ∈ B ∩ σ ⁻¹' (extChartAt 𝓘(ℝ, E) z).source,
      (chartFwd ι z ∘ σ) w ∈ chartDomE ι z := fun w hw => im_chartFwd_mem ι hw.2
  have hev : (F ∘ σ) =ᶠ[𝓝 y] ((F ∘ chartInv ι z) ∘ (chartFwd ι z ∘ σ)) := by
    filter_upwards [hW.mem_nhds hyW] with w hw
    simp [Function.comp_def, im_chartInv_chartFwd ι hw.2]
  rw [hev.fderiv_eq]
  exact ⟨injective_fderiv_comp (im_isOpen_chartDomE ι z) (hF.smooth z) (hF.inj z) hW
      (hσ.trans_smooth z) hmaps (hσ.trans_det z) hyW,
    fun u w => lagrangian_comp (im_isOpen_chartDomE ι z) (hF.smooth z) (hF.lag z) hW
      (hσ.trans_smooth z) hmaps hyW u w⟩

theorem IsLagImm.comp_param {J : V →L[ℝ] V} {F : M → V} (hF : IsLagImm ι J F)
    {B : Set (EuclideanSpace ℝ (Fin n))} {σ : EuclideanSpace ℝ (Fin n) → M} (hσ : IsParam ι B σ) :
    ContDiffOn ℝ ∞ (F ∘ σ) B ∧ (∀ y ∈ B, Function.Injective (fderiv ℝ (F ∘ σ) y)) ∧
      ∀ y ∈ B, ∀ u w, ⟪J (fderiv ℝ (F ∘ σ) y u), fderiv ℝ (F ∘ σ) y w⟫ = 0 :=
  ⟨im_contDiffOn_comp_param ι hF.smooth hσ, fun y hy => (im_fderiv_eq_comp ι hF hσ hy).1,
    fun y hy => (im_fderiv_eq_comp ι hF hσ hy).2⟩

/-- A smooth function on `M`, read in a chart. -/
theorem im_contDiffOn_chart {h : M → ℝ} (hh : ContMDiff 𝓘(ℝ, E) 𝓘(ℝ) ∞ h) (z : M) :
    ContDiffOn ℝ ∞ (h ∘ chartInv ι z) (chartDomE ι z) := by
  have h1 : ContMDiffOn 𝓘(ℝ, E) 𝓘(ℝ) ∞ (h ∘ (extChartAt 𝓘(ℝ, E) z).symm)
      (extChartAt 𝓘(ℝ, E) z).target :=
    hh.comp_contMDiffOn (contMDiffOn_extChartAt_symm z)
  have h2 := h1.contDiffOn
  exact h2.comp ι.contDiff.contDiffOn fun y hy => hy

/-! ### Part 3: natural densities -/

/-- A natural density: local, vanishing where `h` vanishes, natural under changes of coordinates,
continuous. -/
structure IsNatDens (J : V →L[ℝ] V)
    (e : (EuclideanSpace ℝ (Fin n) → V) → (EuclideanSpace ℝ (Fin n) → ℝ) →
      EuclideanSpace ℝ (Fin n) → ℝ) : Prop where
  loc : ∀ X X' h h' z, X =ᶠ[𝓝 z] X' → h =ᶠ[𝓝 z] h' → e X h z = e X' h' z
  vanish : ∀ X h z, h =ᶠ[𝓝 z] 0 → e X h z = 0
  nat : ∀ (X : EuclideanSpace ℝ (Fin n) → V) (h : EuclideanSpace ℝ (Fin n) → ℝ)
    (Ω N : Set (EuclideanSpace ℝ (Fin n))) (τ : EuclideanSpace ℝ (Fin n) → EuclideanSpace ℝ (Fin n)) z,
    IsOpen Ω → ContDiffOn ℝ ∞ X Ω → (∀ x ∈ Ω, Function.Injective (fderiv ℝ X x)) →
    (∀ x ∈ Ω, ∀ u w, ⟪J (fderiv ℝ X x u), fderiv ℝ X x w⟫ = 0) → ContDiffOn ℝ ∞ h Ω →
    IsOpen N → ContDiffOn ℝ ∞ τ N → Set.MapsTo τ N Ω → (∀ z ∈ N, (jacM τ z).det ≠ 0) → z ∈ N →
    e (X ∘ τ) (h ∘ τ) z = e X h (τ z) * |(jacM τ z).det|
  cont : ∀ (X : EuclideanSpace ℝ (Fin n) → V) (h : EuclideanSpace ℝ (Fin n) → ℝ)
    (Ω : Set (EuclideanSpace ℝ (Fin n))),
    IsOpen Ω → ContDiffOn ℝ ∞ X Ω → (∀ x ∈ Ω, Function.Injective (fderiv ℝ X x)) →
    (∀ x ∈ Ω, ∀ u w, ⟪J (fderiv ℝ X x u), fderiv ℝ X x w⟫ = 0) → ContDiffOn ℝ ∞ h Ω →
    ContinuousOn (e X h) Ω

/-- The density `h · d`. -/
def densOf (d : (EuclideanSpace ℝ (Fin n) → V) → EuclideanSpace ℝ (Fin n) → ℝ) :
    (EuclideanSpace ℝ (Fin n) → V) → (EuclideanSpace ℝ (Fin n) → ℝ) →
      EuclideanSpace ℝ (Fin n) → ℝ :=
  fun X h z => h z * d X z

theorem isNatDens_densOf_eDens (J : V →L[ℝ] V) (k : ℕ) :
    IsNatDens (n := n) J (densOf fun X => eDens J X k) where
  loc X X' h h' z hX hh := by
    simp only [densOf]
    rw [eDens_congr hX k, hh.eq_of_nhds]
  vanish X h z hh := by
    simp only [densOf]
    rw [hh.eq_of_nhds]; simp
  nat X h Ω N τ z hΩ hX hinj hlag hh hN hτ hmaps hdet hz := by
    simp only [densOf, Function.comp_apply]
    rw [eDens_comp hΩ hX hinj hlag hN hτ hmaps hdet k hz]
    ring
  cont X h Ω hΩ hX hinj hlag hh :=
    hh.continuousOn.mul (contDiffOn_eDens hΩ hX hinj k).continuousOn

theorem isNatDens_densOf_volX (J : V →L[ℝ] V) : IsNatDens (n := n) J (densOf volX) where
  loc X X' h h' z hX hh := by
    simp only [densOf]
    rw [volX_congr hX, hh.eq_of_nhds]
  vanish X h z hh := by
    simp only [densOf]
    rw [hh.eq_of_nhds]; simp
  nat X h Ω N τ z hΩ hX hinj hlag hh hN hτ hmaps hdet hz := by
    simp only [densOf, Function.comp_apply]
    rw [volX_comp hΩ hX hN hτ hmaps hz]
    ring
  cont X h Ω hΩ hX hinj hlag hh :=
    hh.continuousOn.mul (contDiffOn_volX hΩ hX hinj).continuousOn

/-! ### Part 4: the atlas -/

/-- A finite atlas with a smooth partition of unity subordinate to it. -/
structure IntAtlas (M : Type*) [TopologicalSpace M] [ChartedSpace E M] where
  pts : Finset M
  θ : M → M → ℝ
  smooth : ∀ z ∈ pts, ContMDiff 𝓘(ℝ, E) 𝓘(ℝ) ∞ (θ z)
  nonneg : ∀ z p, 0 ≤ θ z p
  supp : ∀ z ∈ pts, tsupport (θ z) ⊆ (extChartAt 𝓘(ℝ, E) z).source
  sum_eq_one : ∀ p, ∑ z ∈ pts, θ z p = 1

theorem exists_intAtlas [T2Space M] [CompactSpace M] : Nonempty (IntAtlas (E := E) M) := by
  classical
  obtain ⟨t, ht⟩ := isCompact_univ.elim_nhds_subcover (fun z : M => (chartAt E z).source)
    (fun z _ => (chartAt E z).open_source.mem_nhds (mem_chart_source E z))
  obtain ⟨f, hf⟩ := SmoothPartitionOfUnity.exists_isSubordinate (ι := ↥t) 𝓘(ℝ, E) (M := M)
    isClosed_univ (fun i => (chartAt E (i : M)).source) (fun i => (chartAt E (i : M)).open_source)
    (by
      intro p _
      have := ht.2 (Set.mem_univ p)
      simp only [Set.mem_iUnion] at this ⊢
      obtain ⟨i, hi, hp⟩ := this
      exact ⟨⟨i, hi⟩, hp⟩)
  refine ⟨⟨t, fun z p => if hz : z ∈ t then f ⟨z, hz⟩ p else 0, ?_, ?_, ?_, ?_⟩⟩
  · intro z hz
    simp only [hz, dite_true]
    exact (f ⟨z, hz⟩).contMDiff
  · intro z p
    by_cases hz : z ∈ t
    · simp only [hz, dite_true]; exact f.nonneg _ _
    · simp [hz]
  · intro z hz
    simp only [hz, dite_true]
    rw [extChartAt_source]
    exact hf ⟨z, hz⟩
  · intro p
    have h1 := f.sum_eq_one (Set.mem_univ p)
    rw [finsum_eq_sum_of_fintype] at h1
    rw [← h1]
    exact (Finset.sum_coe_sort t fun z => if hz : z ∈ t then f ⟨z, hz⟩ p else 0).symm.trans
      (Finset.sum_congr rfl fun i _ => by simp)

/-! ### Part 5: integration -/

/-- The integral over `M` of the natural density `e` of `F` and `h`. -/
noncomputable def intM (A : IntAtlas (E := E) M) (F : M → V) (h : M → ℝ)
    (e : (EuclideanSpace ℝ (Fin n) → V) → (EuclideanSpace ℝ (Fin n) → ℝ) →
      EuclideanSpace ℝ (Fin n) → ℝ) : ℝ :=
  ∑ z ∈ A.pts, ∫ y in chartDomE ι z,
    A.θ z (chartInv ι z y) * e (F ∘ chartInv ι z) (h ∘ chartInv ι z) y

theorem im_integrableOn_of_vanish {f : EuclideanSpace ℝ (Fin n) → ℝ}
    {B K : Set (EuclideanSpace ℝ (Fin n))}
    (hB : IsOpen B) (hf : ContinuousOn f B) (hK : IsCompact K) (hKB : K ⊆ B)
    (hz : ∀ y ∈ B, y ∉ K → f y = 0) : IntegrableOn f B := by
  have h1 : IntegrableOn f K := (hf.mono hKB).integrableOn_compact hK
  exact h1.of_forall_sdiff_eq_zero hB.measurableSet fun y hy => hz y hy.1 hy.2

theorem im_theta_le_one (A : IntAtlas (E := E) M) {z : M} (hz : z ∈ A.pts) (p : M) :
    A.θ z p ≤ 1 := by
  rw [← A.sum_eq_one p]
  exact Finset.single_le_sum (f := fun z => A.θ z p) (fun w _ => A.nonneg w p) hz

/-- `h ∘ σ` vanishes near `y ∈ B` when `σ y` is outside the support of `h`. -/
theorem im_eventuallyEq_zero_param {h : M → ℝ} {B : Set (EuclideanSpace ℝ (Fin n))}
    {σ : EuclideanSpace ℝ (Fin n) → M} (hσc : ContinuousOn σ B) (hB : IsOpen B)
    {y : EuclideanSpace ℝ (Fin n)} (hy : y ∈ B) (hn : σ y ∉ tsupport h) :
    (h ∘ σ) =ᶠ[𝓝 y] 0 := by
  have h1 : ∀ᶠ w in 𝓝 y, σ w ∈ (tsupport h)ᶜ :=
    (hσc.continuousAt (hB.mem_nhds hy)).preimage_mem_nhds
      ((isClosed_tsupport h).isOpen_compl.mem_nhds hn)
  filter_upwards [h1] with w hw
  show h (σ w) = 0
  exact image_eq_zero_of_notMem_tsupport hw

theorem im_integrableOn_param [T2Space M] [CompactSpace M] {J : V →L[ℝ] V} {F : M → V}
    (hF : IsLagImm ι J F) {h : M → ℝ} (hh : ContMDiff 𝓘(ℝ, E) 𝓘(ℝ) ∞ h) {e} (he : IsNatDens J e)
    {B : Set (EuclideanSpace ℝ (Fin n))} {σ : EuclideanSpace ℝ (Fin n) → M} (hσ : IsParam ι B σ)
    (hsupp : tsupport h ⊆ σ '' B) : IntegrableOn (e (F ∘ σ) (h ∘ σ)) B := by
  obtain ⟨κ, hκc, hκ⟩ := hσ.inv
  have hKc : IsCompact (κ '' tsupport h) :=
    (isClosed_tsupport h).isCompact.image_of_continuousOn (hκc.mono hsupp)
  have hKB : κ '' tsupport h ⊆ B := by
    rintro _ ⟨p, hp, rfl⟩
    obtain ⟨y, hy, rfl⟩ := hsupp hp
    rw [hκ y hy]; exact hy
  obtain ⟨h1, h2, h3⟩ := hF.comp_param ι hσ
  have hcont : ContinuousOn (e (F ∘ σ) (h ∘ σ)) B :=
    he.cont _ _ B hσ.isOpen h1 h2 h3 (im_contDiffOn_comp_param ι (im_contDiffOn_chart ι hh) hσ)
  refine im_integrableOn_of_vanish hσ.isOpen hcont hKc hKB fun y hy hyK => ?_
  have hn : σ y ∉ tsupport h := fun hmem => hyK ⟨σ y, hmem, hκ y hy⟩
  exact he.vanish _ _ _ (im_eventuallyEq_zero_param hσ.cont hσ.isOpen hy hn)

theorem im_integrableOn_theta_mul {A : IntAtlas (E := E) M} {z : M} (hz : z ∈ A.pts)
    {B : Set (EuclideanSpace ℝ (Fin n))} {σ : EuclideanSpace ℝ (Fin n) → M} (hσ : IsParam ι B σ)
    {f : EuclideanSpace ℝ (Fin n) → ℝ} (hf : IntegrableOn f B) (hfc : ContinuousOn f B) :
    IntegrableOn (fun w => A.θ z (σ w) * f w) B := by
  refine Integrable.mono' hf.norm ?_ ?_
  · exact (((A.smooth z hz).continuous.comp_continuousOn hσ.cont).mul hfc).aestronglyMeasurable
      hσ.isOpen.measurableSet
  · refine ae_restrict_of_forall_mem hσ.isOpen.measurableSet fun w _ => ?_
    rw [norm_mul]
    have h1 : ‖A.θ z (σ w)‖ ≤ 1 := by
      rw [Real.norm_eq_abs, abs_of_nonneg (A.nonneg _ _)]; exact im_theta_le_one A hz _
    calc ‖A.θ z (σ w)‖ * ‖f w‖ ≤ 1 * ‖f w‖ := by gcongr
      _ = ‖f w‖ := one_mul _

/-- The term of a single chart, computed in `σ`. -/
theorem im_chart_term {J : V →L[ℝ] V} {F : M → V}
    (hF : IsLagImm ι J F) {h : M → ℝ} (hh : ContMDiff 𝓘(ℝ, E) 𝓘(ℝ) ∞ h) {e} (he : IsNatDens J e)
    {B : Set (EuclideanSpace ℝ (Fin n))} {σ : EuclideanSpace ℝ (Fin n) → M} (hσ : IsParam ι B σ)
    (hsupp : tsupport h ⊆ σ '' B) (z : M) (θ : M → ℝ)
    (hθ : tsupport θ ⊆ (extChartAt 𝓘(ℝ, E) z).source) :
    ∫ y in chartDomE ι z, θ (chartInv ι z y) * e (F ∘ chartInv ι z) (h ∘ chartInv ι z) y =
      ∫ w in B, θ (σ w) * e (F ∘ σ) (h ∘ σ) w := by
  set S := (extChartAt 𝓘(ℝ, E) z).source with hS
  set W := B ∩ σ ⁻¹' S with hWdef
  set τ := chartFwd ι z ∘ σ with hτdef
  have hWo : IsOpen W := im_isOpen_W ι hσ z
  have hτsm : ContDiffOn ℝ ∞ τ W := hσ.trans_smooth z
  have hmaps : ∀ w ∈ W, τ w ∈ chartDomE ι z := fun w hw => im_chartFwd_mem ι hw.2
  have hback : ∀ w ∈ W, chartInv ι z (τ w) = σ w := fun w hw => im_chartInv_chartFwd ι hw.2
  have hτinj : Set.InjOn τ W := by
    intro a ha b hb hab
    have h1 : chartInv ι z (τ a) = chartInv ι z (τ b) := by rw [hab]
    rw [hback a ha, hback b hb] at h1
    exact hσ.injOn ha.1 hb.1 h1
  have hdiff : ∀ w ∈ W, DifferentiableAt ℝ τ w := fun w hw => ir_dτ hWo hτsm hw
  -- Step 1
  have step1 : ∫ y in chartDomE ι z, θ (chartInv ι z y) * e (F ∘ chartInv ι z) (h ∘ chartInv ι z) y =
      ∫ y in τ '' W, θ (chartInv ι z y) * e (F ∘ chartInv ι z) (h ∘ chartInv ι z) y := by
    refine setIntegral_eq_of_subset_of_forall_sdiff_eq_zero (im_isOpen_chartDomE ι z).measurableSet
      ?_ ?_
    · rintro _ ⟨w, hw, rfl⟩; exact hmaps w hw
    · rintro y ⟨hy, hyn⟩
      have hnot : chartInv ι z y ∉ σ '' B := by
        rintro ⟨w, hw, hwy⟩
        have hwW : w ∈ W := ⟨hw, by
          show σ w ∈ S
          rw [hwy]; exact im_chartInv_mem_source ι hy⟩
        apply hyn
        refine ⟨w, hwW, ?_⟩
        show chartFwd ι z (σ w) = y
        rw [hwy]; exact im_chartFwd_chartInv ι hy
      have hn : chartInv ι z y ∉ tsupport h := fun hm => hnot (hsupp hm)
      have hev : (h ∘ chartInv ι z) =ᶠ[𝓝 y] 0 :=
        im_eventuallyEq_zero_param (im_continuousOn_chartInv ι z) (im_isOpen_chartDomE ι z) hy hn
      rw [he.vanish _ _ _ hev, mul_zero]
  -- Step 2
  have step2 := integral_image_eq_integral_abs_det_fderiv_smul (volume : Measure (EuclideanSpace ℝ (Fin n)))
    hWo.measurableSet (f' := fderiv ℝ τ)
    (fun x hx => (hdiff x hx).hasFDerivAt.hasFDerivWithinAt) hτinj
    (fun y => θ (chartInv ι z y) * e (F ∘ chartInv ι z) (h ∘ chartInv ι z) y)
  -- Step 3
  have step3 : ∫ w in W, |(fderiv ℝ τ w).det| •
      (θ (chartInv ι z (τ w)) * e (F ∘ chartInv ι z) (h ∘ chartInv ι z) (τ w)) =
      ∫ w in W, θ (σ w) * e (F ∘ σ) (h ∘ σ) w := by
    refine setIntegral_congr_fun hWo.measurableSet fun w hw => ?_
    simp only [smul_eq_mul]
    rw [im_det_eq_jacM, hback w hw]
    have hev1 : (F ∘ σ) =ᶠ[𝓝 w] ((F ∘ chartInv ι z) ∘ τ) := by
      filter_upwards [hWo.mem_nhds hw] with w' hw'
      simp [Function.comp_def, hback w' hw']
    have hev2 : (h ∘ σ) =ᶠ[𝓝 w] ((h ∘ chartInv ι z) ∘ τ) := by
      filter_upwards [hWo.mem_nhds hw] with w' hw'
      simp [Function.comp_def, hback w' hw']
    have hnat := he.nat (F ∘ chartInv ι z) (h ∘ chartInv ι z) (chartDomE ι z) W τ w
      (im_isOpen_chartDomE ι z) (hF.smooth z) (hF.inj z) (hF.lag z) (im_contDiffOn_chart ι hh z)
      hWo hτsm hmaps (hσ.trans_det z) hw
    rw [he.loc _ _ _ _ w hev1 hev2, hnat]
    ring
  -- Step 4
  have step4 : ∫ w in B, θ (σ w) * e (F ∘ σ) (h ∘ σ) w =
      ∫ w in W, θ (σ w) * e (F ∘ σ) (h ∘ σ) w := by
    refine setIntegral_eq_of_subset_of_forall_sdiff_eq_zero hσ.isOpen.measurableSet
      Set.inter_subset_left ?_
    rintro w ⟨hw, hwn⟩
    have : σ w ∉ tsupport θ := fun hm => hwn ⟨hw, hθ hm⟩
    rw [image_eq_zero_of_notMem_tsupport this, zero_mul]
  rw [step1, step2, step3, step4]

theorem intM_eq_param [T2Space M] [CompactSpace M] (A : IntAtlas (E := E) M) {J : V →L[ℝ] V} {F : M → V}
    (hF : IsLagImm ι J F) {h : M → ℝ} (hh : ContMDiff 𝓘(ℝ, E) 𝓘(ℝ) ∞ h) {e} (he : IsNatDens J e)
    {B : Set (EuclideanSpace ℝ (Fin n))} {σ : EuclideanSpace ℝ (Fin n) → M} (hσ : IsParam ι B σ)
    (hsupp : tsupport h ⊆ σ '' B) :
    intM ι A F h e = ∫ y in B, e (F ∘ σ) (h ∘ σ) y := by
  unfold intM
  have h0 : ∀ z ∈ A.pts, ∫ y in chartDomE ι z,
      A.θ z (chartInv ι z y) * e (F ∘ chartInv ι z) (h ∘ chartInv ι z) y =
      ∫ w in B, A.θ z (σ w) * e (F ∘ σ) (h ∘ σ) w :=
    fun z hz => im_chart_term ι hF hh he hσ hsupp z (A.θ z) (A.supp z hz)
  rw [Finset.sum_congr rfl h0]
  have hint := im_integrableOn_param ι hF hh he hσ hsupp
  obtain ⟨h1, h2, h3⟩ := hF.comp_param ι hσ
  have hcont : ContinuousOn (e (F ∘ σ) (h ∘ σ)) B :=
    he.cont _ _ B hσ.isOpen h1 h2 h3 (im_contDiffOn_comp_param ι (im_contDiffOn_chart ι hh) hσ)
  rw [← integral_finsetSum _ fun z hz => im_integrableOn_theta_mul ι hz hσ hint hcont]
  refine setIntegral_congr_fun hσ.isOpen.measurableSet fun w _ => ?_
  rw [← Finset.sum_mul, A.sum_eq_one, one_mul]

theorem intM_densOf_nonneg (A : IntAtlas (E := E) M) {F : M → V} {h : M → ℝ} (hh : ∀ p, 0 ≤ h p) {d}
    (hd : ∀ X y, 0 ≤ d X y) : 0 ≤ intM ι A F h (densOf d) := by
  refine Finset.sum_nonneg fun z _ => setIntegral_nonneg (im_isOpen_chartDomE ι z).measurableSet
    fun y _ => ?_
  exact mul_nonneg (A.nonneg _ _) (mul_nonneg (hh _) (hd _ _))

theorem intM_densOf_sum [T2Space M] [CompactSpace M] (A : IntAtlas (E := E) M) {J : V →L[ℝ] V}
    {F : M → V} (hF : IsLagImm ι J F) {d} (hd : IsNatDens J (densOf d)) {α : Type*} (s : Finset α)
    (h : α → M → ℝ) (hh : ∀ a ∈ s, ContMDiff 𝓘(ℝ, E) 𝓘(ℝ) ∞ (h a)) :
    intM ι A F (fun p => ∑ a ∈ s, h a p) (densOf d) = ∑ a ∈ s, intM ι A F (h a) (densOf d) := by
  unfold intM
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun z hz => ?_
  have hdc : ContinuousOn (d (F ∘ chartInv ι z)) (chartDomE ι z) := by
    have : ContinuousOn (fun y => 1 * d (F ∘ chartInv ι z) y) (chartDomE ι z) :=
      hd.cont (F ∘ chartInv ι z) (fun _ => 1) (chartDomE ι z) (im_isOpen_chartDomE ι z)
        (hF.smooth z) (hF.inj z) (hF.lag z) contDiffOn_const
    simpa using this
  have hKc : IsCompact ((chartFwd ι z) '' tsupport (A.θ z)) :=
    ((isClosed_tsupport _).isCompact).image_of_continuousOn
      ((im_continuousOn_chartFwd ι z).mono (A.supp z hz))
  have hKB : (chartFwd ι z) '' tsupport (A.θ z) ⊆ chartDomE ι z := by
    rintro _ ⟨p, hp, rfl⟩; exact im_chartFwd_mem ι (A.supp z hz hp)
  have hint : ∀ a ∈ s, Integrable (fun y => A.θ z (chartInv ι z y) *
      (h a (chartInv ι z y) * d (F ∘ chartInv ι z) y)) (volume.restrict (chartDomE ι z)) := by
    intro a ha
    refine im_integrableOn_of_vanish (im_isOpen_chartDomE ι z) ?_ hKc hKB ?_
    · exact (((A.smooth z hz).continuous.comp_continuousOn (im_continuousOn_chartInv ι z)).mul
        (((hh a ha).continuous.comp_continuousOn (im_continuousOn_chartInv ι z)).mul hdc))
    · intro y hy hyK
      have : A.θ z (chartInv ι z y) = 0 := by
        by_contra hne
        exact hyK ⟨chartInv ι z y, subset_tsupport _ (Function.mem_support.2 hne),
          im_chartFwd_chartInv ι hy⟩
      rw [this, zero_mul]
  simp only [densOf, Function.comp_apply]
  rw [← integral_finsetSum _ hint]
  refine setIntegral_congr_fun (im_isOpen_chartDomE ι z).measurableSet fun y _ => ?_
  rw [Finset.sum_mul, Finset.mul_sum]

/-! ### Axiom checks -/

/-- info: 'CWFlow.Extension.exists_intAtlas' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_intAtlas

/-- info: 'CWFlow.Extension.intM_eq_param' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms intM_eq_param

/-- info: 'CWFlow.Extension.isNatDens_densOf_eDens' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms isNatDens_densOf_eDens

/-- info: 'CWFlow.Extension.isNatDens_densOf_volX' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms isNatDens_densOf_volX

/-- info: 'CWFlow.Extension.intM_densOf_sum' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms intM_densOf_sum

end CWFlow.Extension
