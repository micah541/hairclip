import CWFlow.Extension.FlowM
import CWFlow.Extension.VolChart
import CWFlow.Extension.PotReparam

/-!
# The first variation of volume on a compact manifold

For a flow `Φ` of Lagrangian immersions of a compact manifold `M` (`IsFlowM`), the volume `volM` has
derivative `-∫_M div(JH)² dV` (`hasDerivAt_volM`), hence is nonincreasing (`volM_antitoneOn`).

In each chart, `d/dt ∫ w vol = -∫ w p² vol - ∫ divDens w` (`hasDerivAt_integral_volX`), where
`divDens J X h = p vol ⟨μ, dh⟩_g`. Summing over the atlas with the partition of unity `θ_z`, the second
terms add up to `∫ divDens (Σ θ_z) = ∫ divDens 1 = 0`: the terms involving derivatives of the partition
functions cancel. For this, `divDens` is a natural density (`isNatDens_divDens`), linear in `h`
(`intM_divDens_sum`) and vanishes on constants (`intM_divDens_const`).
-/

open MeasureTheory Set
open scoped ContDiff Manifold Topology RealInnerProductSpace

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false

namespace CWFlow.Extension

variable {n : ℕ} {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
  {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
  {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
  (ι : EuclideanSpace ℝ (Fin n) ≃L[ℝ] E)

/-! ### Part 1: two natural densities -/

/-- `div(JH)² dV`. -/
noncomputable def potSqVol (J : V →L[ℝ] V) :
    (EuclideanSpace ℝ (Fin n) → V) → EuclideanSpace ℝ (Fin n) → ℝ :=
  fun X y => potX J X y ^ 2 * volX X y

/-- `div(JH) ⟨μ, dh⟩_g dV`, linear in `h`. -/
noncomputable def divDens (J : V →L[ℝ] V) :
    (EuclideanSpace ℝ (Fin n) → V) → (EuclideanSpace ℝ (Fin n) → ℝ) →
      EuclideanSpace ℝ (Fin n) → ℝ :=
  fun X h y => potX J X y * volX X y *
    ∑ d, ∑ c, (gX X y)⁻¹ d c * muC J X y c * fderiv ℝ h y (eb d)

theorem vm_potX_congr {J : V →L[ℝ] V} {X X' : EuclideanSpace ℝ (Fin n) → V}
    {x : EuclideanSpace ℝ (Fin n)} (h : X =ᶠ[𝓝 x] X') : potX J X x = potX J X' x := by
  unfold potX
  rw [ginvZ_congr h, SX_congr h 1]

theorem vm_muC_congr {J : V →L[ℝ] V} {X X' : EuclideanSpace ℝ (Fin n) → V}
    {x : EuclideanSpace ℝ (Fin n)} (h : X =ᶠ[𝓝 x] X') (c : Fin n) :
    muC J X x c = muC J X' x c := by
  unfold muC
  rw [gX_congr h, cX_congr h]

/-- Finite-sum algebra for the naturality of `divDens`. -/
theorem vm_contract {A g g' : Matrix (Fin n) (Fin n) ℝ}
    (hk : ∀ a' b', ∑ a, ∑ b, A a' a * g'⁻¹ a b * A b' b = g⁻¹ a' b') (μ hv : Fin n → ℝ) :
    ∑ d, ∑ c, g'⁻¹ d c * (∑ c', A c' c * μ c') * (∑ d', A d' d * hv d') =
      ∑ d, ∑ c, g⁻¹ d c * μ c * hv d := by
  have h1 : ∀ d c, g'⁻¹ d c * (∑ c', A c' c * μ c') * (∑ d', A d' d * hv d') =
      ∑ c', ∑ d', (A d' d * g'⁻¹ d c * A c' c) * (μ c' * hv d') := by
    intro d c
    simp only [Finset.mul_sum, Finset.sum_mul]
    refine Finset.sum_comm.trans ?_
    refine Finset.sum_congr rfl fun c' _ => Finset.sum_congr rfl fun d' _ => ?_
    ring
  simp only [h1]
  refine Eq.trans (pr_rev4 (fun d c c' d' => (A d' d * g'⁻¹ d c * A c' c) * (μ c' * hv d'))) ?_
  refine Finset.sum_congr rfl fun d' _ => Finset.sum_congr rfl fun c' _ => ?_
  rw [← hk d' c']
  simp only [Finset.sum_mul]
  refine Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun b _ => ?_
  ring

theorem isNatDens_potSqVol (J : V →L[ℝ] V) (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫)
    (hJJ : ∀ v : V, J (J v) = -v) : IsNatDens (n := n) J (densOf (potSqVol J)) where
  loc X X' h h' z hX hh := by
    simp only [densOf, potSqVol]
    rw [vm_potX_congr hX, volX_congr hX, hh.eq_of_nhds]
  vanish X h z hh := by
    simp only [densOf]
    rw [hh.eq_of_nhds]; simp
  nat X h Ω N τ z hΩ hX hinj hlag hh hN hτ hmaps hdet hz := by
    simp only [densOf, potSqVol, Function.comp_apply]
    rw [potX_comp hΩ hX hinj hlag hN hτ hmaps hdet hz, volX_comp hΩ hX hN hτ hmaps hz]
    ring
  cont X h Ω hΩ hX hinj hlag hh := by
    exact hh.continuousOn.mul
      (((vv_contDiffOn_potX (J := J) hΩ hX hinj).continuousOn.pow 2).mul
        (contDiffOn_volX hΩ hX hinj).continuousOn)

theorem isNatDens_divDens (J : V →L[ℝ] V) (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫)
    (hJJ : ∀ v : V, J (J v) = -v) : IsNatDens (n := n) J (divDens J) where
  loc X X' h h' z hX hh := by
    simp only [divDens]
    rw [vm_potX_congr hX, volX_congr hX, gX_congr hX, hh.fderiv_eq]
    simp only [vm_muC_congr hX]
  vanish X h z hh := by
    have : fderiv ℝ h z = 0 := vc_fderiv_eq_zero hh
    simp [divDens, this]
  nat X h Ω N τ z hΩ hX hinj hlag hh hN hτ hmaps hdet hz := by
    have hA : IsUnit (jacM τ z).det := isUnit_iff_ne_zero.2 (hdet z hz)
    have hk := pr_key hA (gX_comp hΩ hX hN hτ hmaps hz)
    have hmu : ∀ c, muC J (X ∘ τ) z c = ∑ c', jacM τ z c' c * muC J X (τ z) c' := fun c =>
      muC_comp hΩ hX hlag hN hτ hmaps hdet hz c
    have hpd : DifferentiableAt ℝ h (τ z) := ir_dX hΩ hh (hmaps hz)
    have hfd : ∀ d, fderiv ℝ (h ∘ τ) z (eb d) =
        ∑ j, jacM τ z j d * fderiv ℝ h (τ z) (eb j) := by
      intro d
      rw [fderiv_comp z hpd (ir_dτ hN hτ hz)]
      show fderiv ℝ h (τ z) (fderiv ℝ τ z (eb d)) = _
      rw [ir_fderiv_tau]
      simp [map_sum, map_smul]
    simp only [divDens, hmu, hfd]
    rw [potX_comp hΩ hX hinj hlag hN hτ hmaps hdet hz, volX_comp hΩ hX hN hτ hmaps hz]
    have hc := vm_contract hk (fun c' => muC J X (τ z) c') (fun d' => fderiv ℝ h (τ z) (eb d'))
    rw [hc]
    ring
  cont X h Ω hΩ hX hinj hlag hh := by
    have hp := (vv_contDiffOn_potX (J := J) hΩ hX hinj).continuousOn
    have hv := (contDiffOn_volX hΩ hX hinj).continuousOn
    have hmu := vc_contDiffOn_muC (J := J) hΩ hX hinj
    have hgi := contDiffOn_gX_inv hΩ hX hinj
    have hfh : ContDiffOn ℝ ∞ (fderiv ℝ h) Ω := hh.fderiv_of_isOpen hΩ (by simp)
    show ContinuousOn (fun y => potX J X y * volX X y *
      ∑ d, ∑ c, (gX X y)⁻¹ d c * muC J X y c * fderiv ℝ h y (eb d)) Ω
    refine (hp.mul hv).mul ?_
    refine continuousOn_finsetSum _ fun d _ => continuousOn_finsetSum _ fun c _ => ?_
    exact (((hgi d c).continuousOn).mul (hmu c).continuousOn).mul
      (hfh.continuousOn.clm_apply continuousOn_const)

/-! ### Part 2: linearity and vanishing for `divDens` -/

theorem vm_divDens_sum {J : V →L[ℝ] V} {X : EuclideanSpace ℝ (Fin n) → V} {α : Type*}
    (s : Finset α) (φ : α → EuclideanSpace ℝ (Fin n) → ℝ) {y : EuclideanSpace ℝ (Fin n)}
    (hd : ∀ a ∈ s, DifferentiableAt ℝ (φ a) y) :
    divDens J X (fun y => ∑ a ∈ s, φ a y) y = ∑ a ∈ s, divDens J X (φ a) y := by
  simp only [divDens]
  rw [fderiv_fun_sum hd]
  simp only [sum_apply, Finset.mul_sum]
  have : ∀ F : Fin n → Fin n → α → ℝ,
      ∑ d, ∑ c, ∑ a ∈ s, F d c a = ∑ a ∈ s, ∑ d, ∑ c, F d c a := by
    intro F
    calc ∑ d, ∑ c, ∑ a ∈ s, F d c a = ∑ d, ∑ a ∈ s, ∑ c, F d c a :=
          Finset.sum_congr rfl fun d _ => Finset.sum_comm
      _ = _ := Finset.sum_comm
  exact this _

theorem intM_divDens_sum [T2Space M] [CompactSpace M] (A : IntAtlas (E := E) M) {J : V →L[ℝ] V}
    (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v) {F : M → V}
    (hF : IsLagImm ι J F) {α : Type*} (s : Finset α) (h : α → M → ℝ)
    (hh : ∀ a ∈ s, ContMDiff 𝓘(ℝ, E) 𝓘(ℝ) ∞ (h a)) :
    intM ι A F (fun p => ∑ a ∈ s, h a p) (divDens J) =
      ∑ a ∈ s, intM ι A F (h a) (divDens J) := by
  have hd := isNatDens_divDens (n := n) J hJi hJJ
  unfold intM
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun z hz => ?_
  have hKc : IsCompact ((chartFwd ι z) '' tsupport (A.θ z)) :=
    ((isClosed_tsupport _).isCompact).image_of_continuousOn
      ((im_continuousOn_chartFwd ι z).mono (A.supp z hz))
  have hKB : (chartFwd ι z) '' tsupport (A.θ z) ⊆ chartDomE ι z := by
    rintro _ ⟨p, hp, rfl⟩; exact im_chartFwd_mem ι (A.supp z hz hp)
  have hint : ∀ a ∈ s, Integrable (fun y => A.θ z (chartInv ι z y) *
      divDens J (F ∘ chartInv ι z) (h a ∘ chartInv ι z) y) (volume.restrict (chartDomE ι z)) := by
    intro a ha
    have hc := hd.cont (F ∘ chartInv ι z) (h a ∘ chartInv ι z) (chartDomE ι z)
      (im_isOpen_chartDomE ι z) (hF.smooth z) (hF.inj z) (hF.lag z)
      (im_contDiffOn_chart ι (hh a ha) z)
    refine im_integrableOn_of_vanish (im_isOpen_chartDomE ι z) ?_ hKc hKB ?_
    · exact ((A.smooth z hz).continuous.comp_continuousOn (im_continuousOn_chartInv ι z)).mul hc
    · intro y hy hyK
      have : A.θ z (chartInv ι z y) = 0 := by
        by_contra hne
        exact hyK ⟨chartInv ι z y, subset_tsupport _ (Function.mem_support.2 hne),
          im_chartFwd_chartInv ι hy⟩
      rw [this, zero_mul]
  rw [← integral_finsetSum _ hint]
  refine setIntegral_congr_fun (im_isOpen_chartDomE ι z).measurableSet fun y hy => ?_
  have hdiff : ∀ a ∈ s, DifferentiableAt ℝ (h a ∘ chartInv ι z) y := fun a ha =>
    (((im_contDiffOn_chart ι (hh a ha) z).contDiffAt
      ((im_isOpen_chartDomE ι z).mem_nhds hy)).differentiableAt (by simp))
  have := vm_divDens_sum (J := J) (X := F ∘ chartInv ι z) s (fun a => h a ∘ chartInv ι z) hdiff
  show A.θ z (chartInv ι z y) *
      divDens J (F ∘ chartInv ι z) (fun y => ∑ a ∈ s, (h a ∘ chartInv ι z) y) y = _
  rw [this, Finset.mul_sum]

theorem intM_divDens_const (A : IntAtlas (E := E) M) {J : V →L[ℝ] V} {F : M → V} (c : ℝ) :
    intM ι A F (fun _ => c) (divDens J) = 0 := by
  unfold intM
  refine Finset.sum_eq_zero fun z _ => ?_
  have : ∀ y, A.θ z (chartInv ι z y) *
      divDens J (F ∘ chartInv ι z) ((fun _ => c) ∘ chartInv ι z) y = 0 := by
    intro y
    simp [divDens, Function.comp_def]
  simp [this]

/-! ### Part 3: the derivative of the volume -/

theorem vm_setIntegral_eq (θ : M → ℝ) (σ : EuclideanSpace ℝ (Fin n) → M)
    {B : Set (EuclideanSpace ℝ (Fin n))} (hB : MeasurableSet B) (f : EuclideanSpace ℝ (Fin n) → ℝ) :
    ∫ y in B, θ (σ y) * (1 * f y) = ∫ y, B.indicator (fun y => θ (σ y)) y * f y := by
  refine (integral_indicator hB).symm.trans ?_
  refine integral_congr_ae (Filter.Eventually.of_forall fun y => ?_)
  by_cases hy : y ∈ B <;> simp [Set.indicator, hy]

/-- The cut-off `θ_z ∘ chartInv ι z`, extended by zero outside the chart domain. -/
noncomputable def vm_w (A : IntAtlas (E := E) M) (z : M) : EuclideanSpace ℝ (Fin n) → ℝ :=
  (chartDomE ι z).indicator (fun y => A.θ z (chartInv ι z y))

theorem vm_w_eq_zero (A : IntAtlas (E := E) M) {z : M} {y : EuclideanSpace ℝ (Fin n)}
    (hy : y ∉ chartFwd ι z '' tsupport (A.θ z)) : vm_w ι A z y = 0 := by
  unfold vm_w
  by_cases hB : y ∈ chartDomE ι z
  · rw [Set.indicator_of_mem hB]
    by_contra hne
    exact hy ⟨chartInv ι z y, subset_tsupport _ (Function.mem_support.2 hne),
      im_chartFwd_chartInv ι hB⟩
  · exact Set.indicator_of_notMem hB _

theorem vm_K_compact [CompactSpace M] (A : IntAtlas (E := E) M) {z : M} (hz : z ∈ A.pts) :
    IsCompact (chartFwd ι z '' tsupport (A.θ z)) :=
  ((isClosed_tsupport _).isCompact).image_of_continuousOn
    ((im_continuousOn_chartFwd ι z).mono (A.supp z hz))

theorem vm_K_sub (A : IntAtlas (E := E) M) {z : M} (hz : z ∈ A.pts) :
    chartFwd ι z '' tsupport (A.θ z) ⊆ chartDomE ι z := by
  rintro _ ⟨p, hp, rfl⟩; exact im_chartFwd_mem ι (A.supp z hz hp)

theorem vm_w_tsupport [CompactSpace M] (A : IntAtlas (E := E) M) {z : M} (hz : z ∈ A.pts) :
    tsupport (vm_w ι A z) ⊆ chartFwd ι z '' tsupport (A.θ z) :=
  closure_minimal (fun y hy => by
    by_contra h
    exact hy (vm_w_eq_zero ι A h)) (vm_K_compact ι A hz).isClosed

theorem vm_w_compact [CompactSpace M] (A : IntAtlas (E := E) M) {z : M} (hz : z ∈ A.pts) :
    HasCompactSupport (vm_w ι A z) :=
  HasCompactSupport.intro (vm_K_compact ι A hz) fun _ hy => vm_w_eq_zero ι A hy

theorem vm_w_contDiff [CompactSpace M] (A : IntAtlas (E := E) M) {z : M} (hz : z ∈ A.pts) :
    ContDiff ℝ ∞ (vm_w ι A z) := by
  rw [contDiff_iff_contDiffAt]
  intro y
  by_cases hy : y ∈ chartFwd ι z '' tsupport (A.θ z)
  · have hB := vm_K_sub ι A hz hy
    have hev : vm_w ι A z =ᶠ[𝓝 y] fun y => A.θ z (chartInv ι z y) := by
      filter_upwards [(im_isOpen_chartDomE ι z).mem_nhds hB] with w hw
      exact Set.indicator_of_mem hw _
    exact ((im_contDiffOn_chart ι (A.smooth z hz) z).contDiffAt
      ((im_isOpen_chartDomE ι z).mem_nhds hB)).congr_of_eventuallyEq hev
  · have hev : vm_w ι A z =ᶠ[𝓝 y] fun _ => 0 := by
      filter_upwards [(vm_K_compact ι A hz).isClosed.isOpen_compl.mem_nhds hy] with w hw
      exact vm_w_eq_zero ι A hw
    exact contDiffAt_const.congr_of_eventuallyEq hev

/-- The derivative of the term of one chart. -/
theorem vm_hasDerivAt_chart [T2Space M] [CompactSpace M] (A : IntAtlas (E := E) M) {J : V →L[ℝ] V}
    (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    {I : Set ℝ} {Φ : ℝ → M → V} (hΦ : IsFlowM ι J I Φ) {t : ℝ} (ht : t ∈ I) {z : M}
    (hz : z ∈ A.pts) :
    HasDerivAt (fun s => ∫ y in chartDomE ι z, A.θ z (chartInv ι z y) *
        densOf volX (Φ s ∘ chartInv ι z) ((fun _ => (1 : ℝ)) ∘ chartInv ι z) y)
      (-(∫ y in chartDomE ι z, A.θ z (chartInv ι z y) *
          densOf (potSqVol J) (Φ t ∘ chartInv ι z) ((fun _ => (1 : ℝ)) ∘ chartInv ι z) y)
        - intM ι A (Φ t) (A.θ z) (divDens J)) t := by
  have hB := im_isOpen_chartDomE ι z
  have h1 := hasDerivAt_integral_volX (J := J) hJi hJJ
    (Y := fun q : ℝ × EuclideanSpace ℝ (Fin n) => Φ q.1 (chartInv ι z q.2)) (I := I)
    (B := chartDomE ι z) hΦ.isOpen hB (hΦ.smooth z)
    (fun t ht y hy => (hΦ.imm t ht).inj z y hy) (fun t ht y hy => (hΦ.imm t ht).lag z y hy)
    (fun t ht y hy => hΦ.vel t ht z y hy) (vm_w_contDiff ι A hz) (vm_w_compact ι A hz)
    ((vm_w_tsupport ι A hz).trans (vm_K_sub ι A hz)) ht
  have hnd := isNatDens_divDens (n := n) J hJi hJJ
  have hsupp : tsupport (A.θ z) ⊆ chartInv ι z '' chartDomE ι z := fun p hp =>
    ⟨chartFwd ι z p, im_chartFwd_mem ι (A.supp z hz hp), im_chartInv_chartFwd ι (A.supp z hz hp)⟩
  have hvan : ∀ y, y ∉ chartDomE ι z →
      divDens J (Φ t ∘ chartInv ι z) (vm_w ι A z) y = 0 := by
    intro y hy
    refine hnd.vanish _ _ _ ?_
    filter_upwards [(vm_K_compact ι A hz).isClosed.isOpen_compl.mem_nhds
      (fun hyK => hy (vm_K_sub ι A hz hyK))] with w hw
    exact vm_w_eq_zero ι A hw
  have hb : intM ι A (Φ t) (A.θ z) (divDens J) =
      ∫ y, divDens J (Φ t ∘ chartInv ι z) (vm_w ι A z) y := by
    rw [intM_eq_param ι A (hΦ.imm t ht) (A.smooth z hz) hnd (isParam_chart ι z) hsupp]
    refine Eq.trans ?_ (setIntegral_eq_integral_of_forall_compl_eq_zero hvan)
    refine setIntegral_congr_fun hB.measurableSet fun y hy => ?_
    refine hnd.loc _ _ _ _ y Filter.EventuallyEq.rfl ?_
    filter_upwards [hB.mem_nhds hy] with w hw
    exact (Set.indicator_of_mem hw _).symm
  have hfun : (fun s => ∫ y in chartDomE ι z, A.θ z (chartInv ι z y) *
        densOf volX (Φ s ∘ chartInv ι z) ((fun _ => (1 : ℝ)) ∘ chartInv ι z) y) =
      fun s => ∫ y, vm_w ι A z y * volX (fun y' => Φ s (chartInv ι z y')) y := by
    funext s
    exact vm_setIntegral_eq (A.θ z) (chartInv ι z) hB.measurableSet (volX (Φ s ∘ chartInv ι z))
  have ha : ∫ y in chartDomE ι z, A.θ z (chartInv ι z y) *
        densOf (potSqVol J) (Φ t ∘ chartInv ι z) ((fun _ => (1 : ℝ)) ∘ chartInv ι z) y =
      ∫ y, vm_w ι A z y * (potX J (fun y' => Φ t (chartInv ι z y')) y ^ 2 *
        volX (fun y' => Φ t (chartInv ι z y')) y) :=
    vm_setIntegral_eq (A.θ z) (chartInv ι z) hB.measurableSet
      (potSqVol J (Φ t ∘ chartInv ι z))
  rw [hfun, ha, hb]
  exact h1

theorem hasDerivAt_volM [T2Space M] [CompactSpace M] (A : IntAtlas (E := E) M) {J : V →L[ℝ] V}
    (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    {I : Set ℝ} {Φ : ℝ → M → V} (hΦ : IsFlowM ι J I Φ) {t : ℝ} (ht : t ∈ I) :
    HasDerivAt (fun s => volM ι A (Φ s))
      (-intM ι A (Φ t) (fun _ => 1) (densOf (potSqVol J))) t := by
  have hsum := HasDerivAt.fun_sum (u := A.pts)
    (A := fun z s => ∫ y in chartDomE ι z, A.θ z (chartInv ι z y) *
        densOf volX (Φ s ∘ chartInv ι z) ((fun _ => (1 : ℝ)) ∘ chartInv ι z) y)
    (A' := fun z => -(∫ y in chartDomE ι z, A.θ z (chartInv ι z y) *
          densOf (potSqVol J) (Φ t ∘ chartInv ι z) ((fun _ => (1 : ℝ)) ∘ chartInv ι z) y)
        - intM ι A (Φ t) (A.θ z) (divDens J))
    (fun z hz => vm_hasDerivAt_chart ι A hJi hJJ hΦ ht hz)
  have hdiv : ∑ z ∈ A.pts, intM ι A (Φ t) (A.θ z) (divDens J) = 0 := by
    rw [← intM_divDens_sum ι A hJi hJJ (hΦ.imm t ht) A.pts A.θ A.smooth]
    have : (fun p => ∑ z ∈ A.pts, A.θ z p) = fun _ => (1 : ℝ) := funext A.sum_eq_one
    rw [this]
    exact intM_divDens_const ι A 1
  have hE : ∑ z ∈ A.pts, (-(∫ y in chartDomE ι z, A.θ z (chartInv ι z y) *
          densOf (potSqVol J) (Φ t ∘ chartInv ι z) ((fun _ => (1 : ℝ)) ∘ chartInv ι z) y)
        - intM ι A (Φ t) (A.θ z) (divDens J)) =
      -intM ι A (Φ t) (fun _ => 1) (densOf (potSqVol J)) := by
    rw [Finset.sum_sub_distrib, hdiv, sub_zero, Finset.sum_neg_distrib]
    rfl
  rw [hE] at hsum
  exact hsum

/-! ### Part 4: monotonicity -/

theorem volM_antitoneOn [T2Space M] [CompactSpace M] (A : IntAtlas (E := E) M) {J : V →L[ℝ] V}
    (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    {a b : ℝ} {Φ : ℝ → M → V} (hΦ : IsFlowM ι J (Ioo a b) Φ) :
    AntitoneOn (fun s => volM ι A (Φ s)) (Ioo a b) := by
  have hder : ∀ t ∈ Ioo a b, HasDerivAt (fun s => volM ι A (Φ s))
      (-intM ι A (Φ t) (fun _ => 1) (densOf (potSqVol J))) t :=
    fun t ht => hasDerivAt_volM ι A hJi hJJ hΦ ht
  have hnn : ∀ t, 0 ≤ intM ι A (Φ t) (fun _ => (1 : ℝ)) (densOf (potSqVol J)) := fun t =>
    intM_densOf_nonneg ι A (fun _ => zero_le_one) fun X y =>
      mul_nonneg (sq_nonneg _) (Real.sqrt_nonneg _)
  refine antitoneOn_of_deriv_nonpos (convex_Ioo a b)
    (fun t ht => (hder t ht).continuousAt.continuousWithinAt) ?_ ?_
  · rw [interior_Ioo]
    exact fun t ht => (hder t ht).differentiableAt.differentiableWithinAt
  · rw [interior_Ioo]
    intro t ht
    rw [(hder t ht).deriv]
    exact neg_nonpos.2 (hnn t)

/-! ### Axiom checks -/

/-- info: 'CWFlow.Extension.isNatDens_potSqVol' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms isNatDens_potSqVol

/-- info: 'CWFlow.Extension.isNatDens_divDens' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms isNatDens_divDens

/-- info: 'CWFlow.Extension.hasDerivAt_volM' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms hasDerivAt_volM

/-- info: 'CWFlow.Extension.volM_antitoneOn' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms volM_antitoneOn

end CWFlow.Extension
