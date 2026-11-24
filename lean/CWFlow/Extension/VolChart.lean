import CWFlow.Extension.VolEvol
import CWFlow.Extension.SliceSmooth
import Mathlib.Analysis.Calculus.LineDeriv.IntegrationByParts

/-!
# The time derivative of a cut-off volume integral, in one chart

For a smooth family `Y_t : ℝⁿ → V` of Lagrangian immersions (on `I ×ˢ B`) moving by the velocity
`velX J Y_t = J ∇ div(JH)`, and a smooth compactly supported cut-off `w` with `tsupport w ⊆ B`,
`d/dt ∫ w √det g_t = -∫ w p² √det g - ∫ p √det g g^{dc} μ_c ∂_d w`, where `p = potX J Y_t` and
`μ_c = muC J Y_t · c`. The proof differentiates under the integral sign (`hasDerivAt_integral_mul_slice`)
and integrates by parts once in each coordinate direction.
-/

open MeasureTheory
open scoped ContDiff Matrix Topology RealInnerProductSpace

namespace CWFlow.Extension

variable {n : ℕ} {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]

/-! ### Helpers -/

/-- Joint smoothness of `(s, y) ↦ volX Y_s y`. -/
theorem vc_contDiffOn_volX_joint {Y : ℝ × EuclideanSpace ℝ (Fin n) → V}
    {O : Set (ℝ × EuclideanSpace ℝ (Fin n))} (hO : IsOpen O) (hY : ContDiffOn ℝ ∞ Y O)
    (hinj : ∀ q ∈ O, Function.Injective (fderiv ℝ (fun y => Y (q.1, y)) q.2)) :
    ContDiffOn ℝ ∞ (fun q : ℝ × EuclideanSpace ℝ (Fin n) => volX (fun y' => Y (q.1, y')) q.2) O := by
  have hd : ∀ (q : ℝ × EuclideanSpace ℝ (Fin n)) (a : Fin n), q ∈ O →
      dX (fun y' => Y (q.1, y')) q.2 a = fderiv ℝ Y q (0, eb a) := by
    intro q a hq
    have h := vv_hasFDerivAt_slice_y (t := q.1) (y := q.2)
      ((hY.contDiffAt (hO.mem_nhds hq)).differentiableAt (by simp))
    unfold dX
    rw [h.fderiv]
    rfl
  have hDY : ContDiffOn ℝ ∞ (fderiv ℝ Y) O := hY.fderiv_of_isOpen hO (by simp)
  have hg : ∀ a b, ContDiffOn ℝ ∞ (fun q : ℝ × EuclideanSpace ℝ (Fin n) =>
      gX (fun y' => Y (q.1, y')) q.2 a b) O := fun a b =>
    (ContDiffOn.inner ℝ (hDY.clm_apply (contDiffOn_const (c := ((0 : ℝ), eb a))))
      (hDY.clm_apply (contDiffOn_const (c := ((0 : ℝ), eb b))))).congr (fun q hq => by
        show ⟪dX (fun y' => Y (q.1, y')) q.2 a, dX (fun y' => Y (q.1, y')) q.2 b⟫ = _
        rw [hd q a hq, hd q b hq])
  exact (ims_contDiffOn_det hg).sqrt fun q hq => (gX_det_pos (hinj q hq)).ne'

theorem vc_contDiffOn_muC {J : V →L[ℝ] V} {X : EuclideanSpace ℝ (Fin n) → V}
    {Ω : Set (EuclideanSpace ℝ (Fin n))} (hΩ : IsOpen Ω) (hX : ContDiffOn ℝ ∞ X Ω)
    (hinj : ∀ x ∈ Ω, Function.Injective (fderiv ℝ X x)) (c : Fin n) :
    ContDiffOn ℝ ∞ (fun y => muC J X y c) Ω := by
  unfold muC
  refine ContDiffOn.neg ?_
  refine ContDiffOn.sum fun a _ => ContDiffOn.sum fun b _ => ?_
  exact (contDiffOn_gX_inv hΩ hX hinj a b).mul
    ((contDiff_coordG (![a, b, c] : Fin 3 → Fin n)).comp_contDiffOn (contDiffOn_cX (J := J) hΩ hX))

/-- The vector field `Q_d = √det g · Σ_c g^{dc} μ_c`. -/
noncomputable def vc_Q (J : V →L[ℝ] V) (X : EuclideanSpace ℝ (Fin n) → V) (d : Fin n)
    (y : EuclideanSpace ℝ (Fin n)) : ℝ :=
  volX X y * ∑ c, (gX X y)⁻¹ d c * muC J X y c

/-- A cut-off times a function smooth on a neighbourhood of its support is globally smooth. -/
theorem vc_contDiff_mul {w : EuclideanSpace ℝ (Fin n) → ℝ} (hw : ContDiff ℝ ∞ w)
    {B : Set (EuclideanSpace ℝ (Fin n))} (hB : IsOpen B) (hwB : tsupport w ⊆ B)
    {h : EuclideanSpace ℝ (Fin n) → ℝ} (hh : ContDiffOn ℝ ∞ h B) :
    ContDiff ℝ ∞ (fun y => w y * h y) := by
  rw [contDiff_iff_contDiffAt]
  intro y
  by_cases hy : y ∈ tsupport w
  · exact hw.contDiffAt.mul (hh.contDiffAt (hB.mem_nhds (hwB hy)))
  · have : (fun y => w y * h y) =ᶠ[𝓝 y] fun _ => 0 := by
      filter_upwards [(isClosed_tsupport w).isOpen_compl.mem_nhds hy] with z hz
      rw [image_eq_zero_of_notMem_tsupport hz, zero_mul]
    exact contDiffAt_const.congr_of_eventuallyEq this

theorem vc_fderiv_eq_zero {h : EuclideanSpace ℝ (Fin n) → ℝ} {y : EuclideanSpace ℝ (Fin n)}
    (hh : h =ᶠ[𝓝 y] fun _ => 0) : fderiv ℝ h y = 0 := by
  rw [hh.fderiv_eq]
  simp

theorem vc_integrable {w : EuclideanSpace ℝ (Fin n) → ℝ} (hwc : HasCompactSupport w)
    {B : Set (EuclideanSpace ℝ (Fin n))} (hwB : tsupport w ⊆ B)
    {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContinuousOn φ B)
    (h0 : ∀ y, y ∉ tsupport w → φ y = 0) : Integrable φ := by
  have hK : IsCompact (tsupport w) := hwc
  exact ((hφ.mono hwB).integrableOn_compact hK).integrable_of_forall_notMem_eq_zero h0

/-! ### The main theorem -/

theorem hasDerivAt_integral_volX {J : V →L[ℝ] V} (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫)
    (hJJ : ∀ v : V, J (J v) = -v)
    {Y : ℝ × EuclideanSpace ℝ (Fin n) → V} {I : Set ℝ} {B : Set (EuclideanSpace ℝ (Fin n))}
    (hI : IsOpen I) (hB : IsOpen B) (hY : ContDiffOn ℝ ∞ Y (I ×ˢ B))
    (hinj : ∀ t ∈ I, ∀ y ∈ B, Function.Injective (fderiv ℝ (fun y' => Y (t, y')) y))
    (hlag : ∀ t ∈ I, ∀ y ∈ B, ∀ u w,
      ⟪J (fderiv ℝ (fun y' => Y (t, y')) y u), fderiv ℝ (fun y' => Y (t, y')) y w⟫ = 0)
    (hvel : ∀ t ∈ I, ∀ y ∈ B,
      HasDerivAt (fun s => Y (s, y)) (velX J (fun y' => Y (t, y')) y) t)
    {w : EuclideanSpace ℝ (Fin n) → ℝ} (hw : ContDiff ℝ ∞ w) (hwc : HasCompactSupport w)
    (hwB : tsupport w ⊆ B) {t : ℝ} (ht : t ∈ I) :
    HasDerivAt (fun s => ∫ y, w y * volX (fun y' => Y (s, y')) y)
      (-(∫ y, w y * (potX J (fun y' => Y (t, y')) y ^ 2 * volX (fun y' => Y (t, y')) y))
        - ∫ y, potX J (fun y' => Y (t, y')) y * volX (fun y' => Y (t, y')) y *
            ∑ d, ∑ c, (gX (fun y' => Y (t, y')) y)⁻¹ d c * muC J (fun y' => Y (t, y')) y c *
              fderiv ℝ w y (eb d)) t := by
  have hO : IsOpen (I ×ˢ B) := hI.prod hB
  have hinjO : ∀ q ∈ I ×ˢ B, Function.Injective (fderiv ℝ (fun y => Y (q.1, y)) q.2) :=
    fun q hq => hinj q.1 hq.1 q.2 hq.2
  have hYt : ContDiffOn ℝ ∞ (fun y' => Y (t, y')) B :=
    hY.comp (contDiff_const.prodMk contDiff_id).contDiffOn (fun y hy => ⟨ht, hy⟩)
  have hinjt := hinj t ht
  have hlagt := hlag t ht
  have hFs := vc_contDiffOn_volX_joint hO hY hinjO
  have h1 := hasDerivAt_integral_mul_slice hO hFs hw.continuous hwc (t := t)
    (fun x hx => ⟨ht, hwB hx⟩)
  -- pointwise time derivative
  have hP1 : ∀ y ∈ B,
      fderiv ℝ (fun q : ℝ × EuclideanSpace ℝ (Fin n) => volX (fun y' => Y (q.1, y')) q.2) (t, y)
        (1, 0) =
      volX (fun y' => Y (t, y')) y * ∑ c, ∑ d, (gX (fun y' => Y (t, y')) y)⁻¹ c d *
        muC J (fun y' => Y (t, y')) y c *
        fderiv ℝ (potX J (fun y' => Y (t, y'))) y (eb d) := by
    intro y hy
    have hq : (t, y) ∈ I ×ˢ B := ⟨ht, hy⟩
    have hd1 := hasDerivAt_volX_slice hO hY hinjO hq
    have hd2 : HasDerivAt (fun s => volX (fun y' => Y (s, y')) y)
        (fderiv ℝ (fun q : ℝ × EuclideanSpace ℝ (Fin n) => volX (fun y' => Y (q.1, y')) q.2)
          (t, y) (1, 0)) t :=
      vv_hasDerivAt_slice_t (G := EuclideanSpace ℝ (Fin n))
        ((hFs.contDiffAt (hO.mem_nhds hq)).differentiableAt (by simp))
    rw [hd2.unique hd1]
    have hev : (fun y' => deriv (fun s => Y (s, y')) t) =ᶠ[𝓝 y]
        velX J (fun y' => Y (t, y')) := by
      filter_upwards [hB.mem_nhds hy] with y' hy'
      exact (hvel t ht y' hy').deriv
    rw [hev.fderiv_eq, sum_gX_inv_inner_fderiv_velX hJi hJJ hB hYt hinjt hlagt hy]
  -- smoothness on `B`
  have hp : ContDiffOn ℝ ∞ (potX J (fun y' => Y (t, y'))) B := vv_contDiffOn_potX hB hYt hinjt
  have hpc : ContinuousOn (potX J (fun y' => Y (t, y'))) B := hp.continuousOn
  have hfp : ContDiffOn ℝ ∞ (fderiv ℝ (potX J (fun y' => Y (t, y')))) B :=
    hp.fderiv_of_isOpen hB (by simp)
  have hdp : ∀ d, ContinuousOn
      (fun y => fderiv ℝ (potX J (fun y' => Y (t, y'))) y (eb d)) B := fun d =>
    hfp.continuousOn.clm_apply continuousOn_const
  have hvol : ContDiffOn ℝ ∞ (volX (fun y' => Y (t, y'))) B := contDiffOn_volX hB hYt hinjt
  have hmu : ∀ c, ContDiffOn ℝ ∞ (fun y => muC J (fun y' => Y (t, y')) y c) B :=
    vc_contDiffOn_muC hB hYt hinjt
  have hgi := contDiffOn_gX_inv hB hYt hinjt
  have hQs : ∀ d, ContDiffOn ℝ ∞ (vc_Q J (fun y' => Y (t, y')) d) B := fun d =>
    hvol.mul (ContDiffOn.sum fun c _ => (hgi d c).mul (hmu c))
  have hG : ∀ d, ContDiff ℝ ∞ (fun y => w y * vc_Q J (fun y' => Y (t, y')) d y) := fun d =>
    vc_contDiff_mul hw hB hwB (hQs d)
  have hw0 : ∀ y, y ∉ tsupport w → w y = 0 := fun y hy => image_eq_zero_of_notMem_tsupport hy
  have hGz : ∀ d y, y ∉ tsupport w →
      fderiv ℝ (fun y => w y * vc_Q J (fun y' => Y (t, y')) d y) y = 0 := by
    intro d y hy
    refine vc_fderiv_eq_zero ?_
    filter_upwards [(isClosed_tsupport w).isOpen_compl.mem_nhds hy] with z hz
    rw [hw0 z hz, zero_mul]
  have hwz : ∀ y, y ∉ tsupport w → fderiv ℝ w y = 0 := by
    intro y hy
    refine vc_fderiv_eq_zero ?_
    filter_upwards [(isClosed_tsupport w).isOpen_compl.mem_nhds hy] with z hz
    exact hw0 z hz
  -- integration by parts in each direction
  have hibp : ∀ d, ∫ y, fderiv ℝ (potX J (fun y' => Y (t, y'))) y (eb d) *
        (w y * vc_Q J (fun y' => Y (t, y')) d y) =
      -∫ y, potX J (fun y' => Y (t, y')) y *
        fderiv ℝ (fun y => w y * vc_Q J (fun y' => Y (t, y')) d y) y (eb d) := by
    intro d
    have hi1 : Integrable (fun y => fderiv ℝ (potX J (fun y' => Y (t, y'))) y (eb d) *
        (w y * vc_Q J (fun y' => Y (t, y')) d y)) :=
      vc_integrable hwc hwB ((hdp d).mul (hw.continuous.continuousOn.mul (hQs d).continuousOn))
        (fun y hy => by rw [hw0 y hy]; simp)
    have hi2 : Integrable (fun y => potX J (fun y' => Y (t, y')) y *
        fderiv ℝ (fun y => w y * vc_Q J (fun y' => Y (t, y')) d y) y (eb d)) :=
      vc_integrable hwc hwB (hpc.mul
        ((((hG d).continuous_fderiv (by simp)).clm_apply continuous_const).continuousOn))
        (fun y hy => by rw [hGz d y hy]; simp)
    have hi3 : Integrable (fun y => potX J (fun y' => Y (t, y')) y *
        (w y * vc_Q J (fun y' => Y (t, y')) d y)) :=
      vc_integrable hwc hwB (hpc.mul (hw.continuous.continuousOn.mul (hQs d).continuousOn))
        (fun y hy => by rw [hw0 y hy]; simp)
    have := integral_mul_fderiv_eq_neg_fderiv_mul_of_integrable
      (f := potX J (fun y' => Y (t, y')))
      (g := fun y => w y * vc_Q J (fun y' => Y (t, y')) d y) (v := eb d) hi1 hi2 hi3
      (fun x hx => (hp.contDiffAt (hB.mem_nhds (hwB (tsupport_mul_subset_left hx)))).differentiableAt
        (by simp))
      (fun x _ => (hG d).differentiable (by simp) x)
    rw [this, neg_neg]
  have hi2' : ∀ d, Integrable (fun y => potX J (fun y' => Y (t, y')) y *
        fderiv ℝ (fun y => w y * vc_Q J (fun y' => Y (t, y')) d y) y (eb d)) := fun d =>
    vc_integrable hwc hwB (hpc.mul
      ((((hG d).continuous_fderiv (by simp)).clm_apply continuous_const).continuousOn))
      (fun y hy => by rw [hGz d y hy]; simp)
  have hi1' : ∀ d, Integrable (fun y => fderiv ℝ (potX J (fun y' => Y (t, y'))) y (eb d) *
        (w y * vc_Q J (fun y' => Y (t, y')) d y)) := fun d =>
    vc_integrable hwc hwB ((hdp d).mul (hw.continuous.continuousOn.mul (hQs d).continuousOn))
      (fun y hy => by rw [hw0 y hy]; simp)
  -- pointwise rewriting of the integrand
  have e1 : ∀ y, w y * fderiv ℝ (fun q : ℝ × EuclideanSpace ℝ (Fin n) =>
        volX (fun y' => Y (q.1, y')) q.2) (t, y) (1, 0) =
      ∑ d, fderiv ℝ (potX J (fun y' => Y (t, y'))) y (eb d) *
        (w y * vc_Q J (fun y' => Y (t, y')) d y) := by
    intro y
    by_cases hy : y ∈ tsupport w
    · rw [hP1 y (hwB hy)]
      simp only [vc_Q, Finset.mul_sum]
      rw [Finset.sum_comm]
      refine Finset.sum_congr rfl fun d _ => Finset.sum_congr rfl fun c _ => ?_
      rw [gX_inv_symm (hinjt y (hwB hy)) c d]
      ring
    · rw [hw0 y hy]
      simp
  have e2 : ∀ y, ∑ d, potX J (fun y' => Y (t, y')) y *
        fderiv ℝ (fun y => w y * vc_Q J (fun y' => Y (t, y')) d y) y (eb d) =
      w y * (potX J (fun y' => Y (t, y')) y ^ 2 * volX (fun y' => Y (t, y')) y) +
        potX J (fun y' => Y (t, y')) y * volX (fun y' => Y (t, y')) y *
          ∑ d, ∑ c, (gX (fun y' => Y (t, y')) y)⁻¹ d c * muC J (fun y' => Y (t, y')) y c *
            fderiv ℝ w y (eb d) := by
    intro y
    by_cases hy : y ∈ tsupport w
    · have hyB := hwB hy
      have hprod : ∀ d, fderiv ℝ (fun y => w y * vc_Q J (fun y' => Y (t, y')) d y) y (eb d) =
          w y * fderiv ℝ (vc_Q J (fun y' => Y (t, y')) d) y (eb d) +
            vc_Q J (fun y' => Y (t, y')) d y * fderiv ℝ w y (eb d) := by
        intro d
        have hwd : DifferentiableAt ℝ w y := hw.differentiable (by simp) y
        have hQd : DifferentiableAt ℝ (vc_Q J (fun y' => Y (t, y')) d) y :=
          ((hQs d).contDiffAt (hB.mem_nhds hyB)).differentiableAt (by simp)
        rw [fderiv_fun_mul hwd hQd]
        simp only [ContinuousLinearMap.add_apply, ContinuousLinearMap.smul_apply, smul_eq_mul]
      have hsq : ∑ d, fderiv ℝ (vc_Q J (fun y' => Y (t, y')) d) y (eb d) =
          volX (fun y' => Y (t, y')) y * potX J (fun y' => Y (t, y')) y :=
        sum_fderiv_volX_muC hB hYt hinjt hyB
      have hQw : ∀ d, potX J (fun y' => Y (t, y')) y *
          (vc_Q J (fun y' => Y (t, y')) d y * fderiv ℝ w y (eb d)) =
          potX J (fun y' => Y (t, y')) y * volX (fun y' => Y (t, y')) y *
            ∑ c, (gX (fun y' => Y (t, y')) y)⁻¹ d c * muC J (fun y' => Y (t, y')) y c *
              fderiv ℝ w y (eb d) := by
        intro d
        rw [← Finset.sum_mul]
        unfold vc_Q
        ring
      calc ∑ d, potX J (fun y' => Y (t, y')) y *
            fderiv ℝ (fun y => w y * vc_Q J (fun y' => Y (t, y')) d y) y (eb d)
          = ∑ d, (potX J (fun y' => Y (t, y')) y * w y *
              fderiv ℝ (vc_Q J (fun y' => Y (t, y')) d) y (eb d) +
            potX J (fun y' => Y (t, y')) y *
              (vc_Q J (fun y' => Y (t, y')) d y * fderiv ℝ w y (eb d))) :=
            Finset.sum_congr rfl fun d _ => by rw [hprod d]; ring
        _ = potX J (fun y' => Y (t, y')) y * w y *
              ∑ d, fderiv ℝ (vc_Q J (fun y' => Y (t, y')) d) y (eb d) +
            ∑ d, potX J (fun y' => Y (t, y')) y *
              (vc_Q J (fun y' => Y (t, y')) d y * fderiv ℝ w y (eb d)) := by
            rw [Finset.sum_add_distrib, Finset.mul_sum]
        _ = _ := by
            rw [hsq, Finset.sum_congr rfl fun d _ => hQw d, ← Finset.mul_sum]
            ring
    · rw [hw0 y hy, hwz y hy]
      simp [hGz _ y hy]
  have hA : Integrable (fun y => w y * (potX J (fun y' => Y (t, y')) y ^ 2 *
      volX (fun y' => Y (t, y')) y)) :=
    vc_integrable hwc hwB (hw.continuous.continuousOn.mul ((hpc.pow 2).mul hvol.continuousOn))
      (fun y hy => by rw [hw0 y hy]; simp)
  have hBi : Integrable (fun y => potX J (fun y' => Y (t, y')) y *
      volX (fun y' => Y (t, y')) y * ∑ d, ∑ c, (gX (fun y' => Y (t, y')) y)⁻¹ d c *
        muC J (fun y' => Y (t, y')) y c * fderiv ℝ w y (eb d)) := by
    refine vc_integrable hwc hwB ?_ (fun y hy => by rw [hwz y hy]; simp)
    refine (hpc.mul hvol.continuousOn).mul ?_
    refine continuousOn_finsetSum _ fun d _ => continuousOn_finsetSum _ fun c _ => ?_
    exact (((hgi d c).continuousOn).mul (hmu c).continuousOn).mul
      (((hw.continuous_fderiv (by simp)).clm_apply continuous_const).continuousOn)
  have hval : ∫ y, w y * fderiv ℝ (fun q : ℝ × EuclideanSpace ℝ (Fin n) =>
        volX (fun y' => Y (q.1, y')) q.2) (t, y) (1, 0) =
      -(∫ y, w y * (potX J (fun y' => Y (t, y')) y ^ 2 * volX (fun y' => Y (t, y')) y))
        - ∫ y, potX J (fun y' => Y (t, y')) y * volX (fun y' => Y (t, y')) y *
            ∑ d, ∑ c, (gX (fun y' => Y (t, y')) y)⁻¹ d c * muC J (fun y' => Y (t, y')) y c *
              fderiv ℝ w y (eb d) := by
    calc ∫ y, w y * fderiv ℝ (fun q : ℝ × EuclideanSpace ℝ (Fin n) =>
            volX (fun y' => Y (q.1, y')) q.2) (t, y) (1, 0)
        = ∫ y, ∑ d, fderiv ℝ (potX J (fun y' => Y (t, y'))) y (eb d) *
            (w y * vc_Q J (fun y' => Y (t, y')) d y) :=
          integral_congr_ae (Filter.Eventually.of_forall e1)
      _ = ∑ d, ∫ y, fderiv ℝ (potX J (fun y' => Y (t, y'))) y (eb d) *
            (w y * vc_Q J (fun y' => Y (t, y')) d y) :=
          integral_finsetSum _ fun d _ => hi1' d
      _ = ∑ d, -∫ y, potX J (fun y' => Y (t, y')) y *
            fderiv ℝ (fun y => w y * vc_Q J (fun y' => Y (t, y')) d y) y (eb d) :=
          Finset.sum_congr rfl fun d _ => hibp d
      _ = -∫ y, ∑ d, potX J (fun y' => Y (t, y')) y *
            fderiv ℝ (fun y => w y * vc_Q J (fun y' => Y (t, y')) d y) y (eb d) := by
          rw [Finset.sum_neg_distrib, integral_finsetSum _ fun d _ => hi2' d]
      _ = -∫ y, (w y * (potX J (fun y' => Y (t, y')) y ^ 2 * volX (fun y' => Y (t, y')) y) +
            potX J (fun y' => Y (t, y')) y * volX (fun y' => Y (t, y')) y *
              ∑ d, ∑ c, (gX (fun y' => Y (t, y')) y)⁻¹ d c * muC J (fun y' => Y (t, y')) y c *
                fderiv ℝ w y (eb d)) := by
          rw [integral_congr_ae (Filter.Eventually.of_forall e2)]
      _ = _ := by rw [integral_add hA hBi]; ring
  exact h1.congr_deriv hval

/-- info: 'CWFlow.Extension.hasDerivAt_integral_volX' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms hasDerivAt_integral_volX

end CWFlow.Extension
