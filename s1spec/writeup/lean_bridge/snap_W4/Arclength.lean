import CWFlow.Curve.Periodic
import Mathlib.Analysis.Calculus.FDeriv.Symmetric
import Mathlib.Analysis.Calculus.ParametricIntervalIntegral
import Mathlib.Analysis.Calculus.MeanValue

/-!
# Arclength calculus for a time-dependent closed curve (SCOPE W3)

`γ : ℝ × ℝ → ℂ`, `(t, θ) ↦ γ(t, θ)`, smooth on `I × ℝ` (`I` open), `2π`-periodic in `θ`, regular, moving by
curve diffusion `∂ₜγ = −κ_ss N` (`IsCDF`).

* partial derivatives `pθ f`, `pt f` (derivatives of the slices), their smoothness, the symmetry of mixed partials;
* the curve quantities as functions of `(t, θ)`: `G1 = ∂θγ`, `v = |∂θγ|`, `T`, `N = iT`, `κ`, `Ds f = ∂θf / v`;
  they are *definitionally* the quantities `deriv`, `curv`, `dS` of the slice `θ ↦ γ(t, θ)`;
* Frenet formulas and the evolution equations (E1)–(E4) for a general normal speed;
* differentiation under the integral sign.
-/

open Set Function Filter intervalIntegral Real
open scoped RealInnerProductSpace ContDiff Topology

noncomputable section

namespace CWFlow.Curve

variable {F : Type*} [NormedAddCommGroup F] [NormedSpace ℝ F]

/-- The `θ`-partial derivative: the derivative of the slice `θ ↦ f (t, θ)`. -/
def pθ (f : ℝ × ℝ → F) (q : ℝ × ℝ) : F := deriv (fun θ => f (q.1, θ)) q.2

/-- The `t`-partial derivative: the derivative of the slice `s ↦ f (s, θ)`. -/
def pt (f : ℝ × ℝ → F) (q : ℝ × ℝ) : F := deriv (fun s => f (s, q.2)) q.1

section Partials

variable {U : Set (ℝ × ℝ)} {f g : ℝ × ℝ → F}

theorem hasDerivAt_slice_θ {t θ : ℝ} (hf : DifferentiableAt ℝ f (t, θ)) :
    HasDerivAt (fun x => f (t, x)) (fderiv ℝ f (t, θ) (0, 1)) θ := by
  have h : HasDerivAt (fun x : ℝ => ((t, x) : ℝ × ℝ)) ((0, 1) : ℝ × ℝ) θ :=
    (hasDerivAt_const θ t).prodMk (hasDerivAt_id θ)
  exact hf.hasFDerivAt.comp_hasDerivAt θ h

theorem hasDerivAt_slice_t {t θ : ℝ} (hf : DifferentiableAt ℝ f (t, θ)) :
    HasDerivAt (fun s => f (s, θ)) (fderiv ℝ f (t, θ) (1, 0)) t := by
  have h : HasDerivAt (fun s : ℝ => ((s, θ) : ℝ × ℝ)) ((1, 0) : ℝ × ℝ) t :=
    (hasDerivAt_id t).prodMk (hasDerivAt_const t θ)
  exact hf.hasFDerivAt.comp_hasDerivAt t h

theorem pθ_eq_fderiv {q : ℝ × ℝ} (hf : DifferentiableAt ℝ f q) : pθ f q = fderiv ℝ f q (0, 1) :=
  (hasDerivAt_slice_θ (t := q.1) (θ := q.2) hf).deriv

theorem pt_eq_fderiv {q : ℝ × ℝ} (hf : DifferentiableAt ℝ f q) : pt f q = fderiv ℝ f q (1, 0) :=
  (hasDerivAt_slice_t (t := q.1) (θ := q.2) hf).deriv

theorem hasDerivAt_pθ {t θ : ℝ} (hf : DifferentiableAt ℝ f (t, θ)) :
    HasDerivAt (fun x => f (t, x)) (pθ f (t, θ)) θ := by
  rw [pθ_eq_fderiv hf]; exact hasDerivAt_slice_θ hf

theorem hasDerivAt_pt {t θ : ℝ} (hf : DifferentiableAt ℝ f (t, θ)) :
    HasDerivAt (fun s => f (s, θ)) (pt f (t, θ)) t := by
  rw [pt_eq_fderiv hf]; exact hasDerivAt_slice_t hf

theorem _root_.ContDiffOn.differentiableAt' (hU : IsOpen U) (hf : ContDiffOn ℝ ∞ f U) {q : ℝ × ℝ}
    (hq : q ∈ U) : DifferentiableAt ℝ f q :=
  (hf.contDiffAt (hU.mem_nhds hq)).differentiableAt (by simp)

theorem contDiffOn_pθ (hU : IsOpen U) (hf : ContDiffOn ℝ ∞ f U) : ContDiffOn ℝ ∞ (pθ f) U := by
  have h : ContDiffOn ℝ ∞ (fun q => fderiv ℝ f q ((0, 1) : ℝ × ℝ)) U :=
    (hf.fderiv_of_isOpen hU (m := ∞) (by simp)).clm_apply contDiffOn_const
  exact h.congr fun q hq => pθ_eq_fderiv (hf.differentiableAt' hU hq)

theorem contDiffOn_pt (hU : IsOpen U) (hf : ContDiffOn ℝ ∞ f U) : ContDiffOn ℝ ∞ (pt f) U := by
  have h : ContDiffOn ℝ ∞ (fun q => fderiv ℝ f q ((1, 0) : ℝ × ℝ)) U :=
    (hf.fderiv_of_isOpen hU (m := ∞) (by simp)).clm_apply contDiffOn_const
  exact h.congr fun q hq => pt_eq_fderiv (hf.differentiableAt' hU hq)

/-- Locality of `pθ`. -/
theorem pθ_congr (hU : IsOpen U) (h : ∀ q ∈ U, f q = g q) {q : ℝ × ℝ} (hq : q ∈ U) :
    pθ f q = pθ g q := by
  unfold pθ
  apply Filter.EventuallyEq.deriv_eq
  have hc : Continuous fun x : ℝ => ((q.1, x) : ℝ × ℝ) := by fun_prop
  have := hc.continuousAt.preimage_mem_nhds (x := q.2) (hU.mem_nhds (by simpa using hq))
  filter_upwards [this] with x hx using h _ hx

/-- Locality of `pt`. -/
theorem pt_congr (hU : IsOpen U) (h : ∀ q ∈ U, f q = g q) {q : ℝ × ℝ} (hq : q ∈ U) :
    pt f q = pt g q := by
  unfold pt
  apply Filter.EventuallyEq.deriv_eq
  have hc : Continuous fun s : ℝ => ((s, q.2) : ℝ × ℝ) := by fun_prop
  have := hc.continuousAt.preimage_mem_nhds (x := q.1) (hU.mem_nhds (by simpa using hq))
  filter_upwards [this] with x hx using h _ hx

/-- Symmetry of mixed partials. -/
theorem pt_pθ (hU : IsOpen U) (hf : ContDiffOn ℝ ∞ f U) {q : ℝ × ℝ} (hq : q ∈ U) :
    pt (pθ f) q = pθ (pt f) q := by
  have hfd : ContDiffOn ℝ ∞ (fderiv ℝ f) U := hf.fderiv_of_isOpen hU (m := ∞) (by simp)
  have hdd : DifferentiableAt ℝ (fderiv ℝ f) q := hfd.differentiableAt' hU hq
  have e1 : pt (pθ f) q = pt (fun q => fderiv ℝ f q ((0, 1) : ℝ × ℝ)) q :=
    pt_congr hU (fun q hq => pθ_eq_fderiv (hf.differentiableAt' hU hq)) hq
  have e2 : pθ (pt f) q = pθ (fun q => fderiv ℝ f q ((1, 0) : ℝ × ℝ)) q :=
    pθ_congr hU (fun q hq => pt_eq_fderiv (hf.differentiableAt' hU hq)) hq
  have hA : HasFDerivAt (fun q => fderiv ℝ f q ((0, 1) : ℝ × ℝ))
      ((fderiv ℝ (fderiv ℝ f) q).flip ((0, 1) : ℝ × ℝ)) q := by
    have := hdd.hasFDerivAt.clm_apply (hasFDerivAt_const ((0, 1) : ℝ × ℝ) q)
    simpa using this
  have hB : HasFDerivAt (fun q => fderiv ℝ f q ((1, 0) : ℝ × ℝ))
      ((fderiv ℝ (fderiv ℝ f) q).flip ((1, 0) : ℝ × ℝ)) q := by
    have := hdd.hasFDerivAt.clm_apply (hasFDerivAt_const ((1, 0) : ℝ × ℝ) q)
    simpa using this
  rw [e1, e2, pt_eq_fderiv hA.differentiableAt, pθ_eq_fderiv hB.differentiableAt, hA.fderiv,
    hB.fderiv]
  simp only [ContinuousLinearMap.flip_apply]
  have hs := (hf.contDiffAt (hU.mem_nhds hq)).isSymmSndFDerivAt (n := ∞)
    (by simp; exact WithTop.coe_le_coe.2 le_top)
  exact hs _ _

theorem isOpen_U {I : Set ℝ} (hI : IsOpen I) : IsOpen (I ×ˢ (univ : Set ℝ)) := hI.prod isOpen_univ

/-- A slice of a function smooth on `I × ℝ` is smooth. -/
theorem contDiff_slice {I : Set ℝ} {f : ℝ × ℝ → F} (hf : ContDiffOn ℝ ∞ f (I ×ˢ univ)) {t : ℝ}
    (ht : t ∈ I) : ContDiff ℝ ∞ (fun x => f (t, x)) :=
  hf.comp_contDiff (by fun_prop : ContDiff ℝ ∞ fun x : ℝ => ((t, x) : ℝ × ℝ))
    fun x => ⟨ht, trivial⟩

end Partials

/-! ### Differentiation under the integral sign -/

theorem hasDerivAt_integral_of_contDiffOn {Φ : ℝ × ℝ → ℝ} {I : Set ℝ} (hI : IsOpen I)
    (hΦ : ContDiffOn ℝ ∞ Φ (I ×ˢ univ)) {t : ℝ} (ht : t ∈ I) (a b : ℝ) :
    HasDerivAt (fun s => ∫ θ in a..b, Φ (s, θ)) (∫ θ in a..b, pt Φ (t, θ)) t := by
  have hU : IsOpen (I ×ˢ (univ : Set ℝ)) := hI.prod isOpen_univ
  obtain ⟨ε, hε, hball⟩ := Metric.isOpen_iff.1 hI t ht
  set δ := ε / 2 with hδ
  have hδ0 : 0 < δ := by positivity
  have hIcc : Icc (t - δ) (t + δ) ⊆ I := by
    intro x hx
    apply hball
    rw [Metric.mem_ball, Real.dist_eq, abs_lt]
    constructor <;> linarith [hx.1, hx.2]
  have hptc : ContinuousOn (pt Φ) (I ×ˢ univ) := (contDiffOn_pt hU hΦ).continuousOn
  have hK : IsCompact (Icc (t - δ) (t + δ) ×ˢ uIcc a b) := isCompact_Icc.prod isCompact_uIcc
  obtain ⟨C, hC⟩ := hK.exists_bound_of_continuousOn
    (hptc.mono (prod_mono hIcc (subset_univ _)))
  have hslice : ∀ x ∈ I, Continuous fun θ => Φ (x, θ) := fun x hx =>
    continuous_iff_continuousAt.2 fun θ =>
      ((hΦ.differentiableAt' hU (show (x, θ) ∈ I ×ˢ univ from ⟨hx, trivial⟩)).continuousAt).comp
        (by fun_prop : Continuous fun θ : ℝ => ((x, θ) : ℝ × ℝ)).continuousAt
  have hslice' : Continuous fun θ => pt Φ (t, θ) :=
    continuous_iff_continuousAt.2 fun θ =>
      ((hptc.continuousAt (hU.mem_nhds (show (t, θ) ∈ I ×ˢ univ from ⟨ht, trivial⟩)))).comp
        (by fun_prop : Continuous fun θ : ℝ => ((t, θ) : ℝ × ℝ)).continuousAt
  have hballδ : Metric.ball t δ ⊆ Icc (t - δ) (t + δ) := by
    intro x hx
    rw [Metric.mem_ball, Real.dist_eq, abs_lt] at hx
    constructor <;> linarith [hx.1, hx.2]
  have key := hasDerivAt_integral_of_dominated_loc_of_deriv_le (μ := MeasureTheory.volume) (a := a) (b := b)
    (F := fun s θ => Φ (s, θ)) (F' := fun s θ => pt Φ (s, θ)) (x₀ := t) (bound := fun _ => C)
    (Metric.ball_mem_nhds t hδ0)
    (by
      filter_upwards [hI.mem_nhds ht] with x hx
      exact (hslice x hx).aestronglyMeasurable)
    ((hslice t ht).intervalIntegrable a b)
    hslice'.aestronglyMeasurable
    (Filter.Eventually.of_forall fun θ hθ x hx => hC (x, θ) ⟨hballδ hx, uIoc_subset_uIcc hθ⟩)
    intervalIntegrable_const
    (Filter.Eventually.of_forall fun θ _ x hx =>
      hasDerivAt_pt (hΦ.differentiableAt' hU ⟨hIcc (hballδ hx), trivial⟩))
  exact key.2

/-! ### Plane algebra in the frame `(T, iT)` -/

theorem inner_frame {T : ℂ} (hT : ‖T‖ = 1) (a b c d : ℝ) :
    ⟪(a : ℂ) * (Complex.I * T) + (b : ℂ) * T, (c : ℂ) * (Complex.I * T) + (d : ℂ) * T⟫ = a * c + b * d := by
  have h : T.re * T.re + T.im * T.im = 1 := by
    rw [← Complex.normSq_apply, ← Complex.sq_norm, hT, one_pow]
  simp only [Complex.inner, Complex.mul_re, Complex.mul_im, Complex.add_re, Complex.add_im,
    Complex.I_re, Complex.I_im, Complex.ofReal_re, Complex.ofReal_im, Complex.conj_re, Complex.conj_im]
  linear_combination (a * c + b * d) * h

/-! ### The curve quantities as functions of `(t, θ)` -/

section CurveQ

variable (γ : ℝ × ℝ → ℂ)

/-- Speed `v = |∂θγ|`. -/
def cSpd (q : ℝ × ℝ) : ℝ := ‖pθ γ q‖
/-- Unit tangent `T = ∂θγ / v`. -/
def cTan (q : ℝ × ℝ) : ℂ := pθ γ q / ((cSpd γ q : ℝ) : ℂ)
/-- Unit normal `N = iT`. -/
def cNor (q : ℝ × ℝ) : ℂ := Complex.I * cTan γ q
/-- Curvature `κ = ⟪∂θ²γ, i∂θγ⟫ / v³`. -/
def cKap (q : ℝ × ℝ) : ℝ := ⟪pθ (pθ γ) q, Complex.I * pθ γ q⟫ / cSpd γ q ^ 3
/-- Arclength derivative `∂_s f = ∂θ f / v`. -/
def cDs (f : ℝ × ℝ → ℝ) (q : ℝ × ℝ) : ℝ := pθ f q / cSpd γ q

/-- The joint quantities are those of the slice, definitionally. -/
theorem cKap_eq_curv (t θ : ℝ) : cKap γ (t, θ) = curv (fun x => γ (t, x)) θ := rfl
theorem cDs_eq_dS (f : ℝ × ℝ → ℝ) (t θ : ℝ) :
    cDs γ f (t, θ) = dS (fun x => γ (t, x)) (fun x => f (t, x)) θ := rfl
theorem cSpd_eq (t θ : ℝ) : cSpd γ (t, θ) = ‖deriv (fun x => γ (t, x)) θ‖ := rfl

variable {γ} {I : Set ℝ} (hI : IsOpen I) (hγ : ContDiffOn ℝ ∞ γ (I ×ˢ univ))
  (hreg : ∀ q ∈ I ×ˢ (univ : Set ℝ), pθ γ q ≠ 0)
include hI hγ hreg

theorem cSpd_pos {q : ℝ × ℝ} (hq : q ∈ I ×ˢ (univ : Set ℝ)) : 0 < cSpd γ q :=
  norm_pos_iff.2 (hreg q hq)

omit hreg in
theorem contDiffOn_G1 : ContDiffOn ℝ ∞ (pθ γ) (I ×ˢ univ) := contDiffOn_pθ (isOpen_U hI) hγ

omit hreg in
theorem contDiffOn_G2 : ContDiffOn ℝ ∞ (pθ (pθ γ)) (I ×ˢ univ) :=
  contDiffOn_pθ (isOpen_U hI) (contDiffOn_G1 hI hγ)

theorem contDiffOn_cSpd : ContDiffOn ℝ ∞ (cSpd γ) (I ×ˢ univ) :=
  (contDiffOn_G1 hI hγ).norm ℝ hreg

theorem contDiffOn_cTan : ContDiffOn ℝ ∞ (cTan γ) (I ×ˢ univ) := by
  have h : ContDiffOn ℝ ∞ (fun q => (cSpd γ q)⁻¹ • pθ γ q) (I ×ˢ univ) :=
    ((contDiffOn_cSpd hI hγ hreg).inv fun q hq => (cSpd_pos hI hγ hreg hq).ne').smul
      (contDiffOn_G1 hI hγ)
  refine h.congr fun q _ => ?_
  simp only [cTan, Complex.real_smul, Complex.ofReal_inv, div_eq_inv_mul]

theorem contDiffOn_cNor : ContDiffOn ℝ ∞ (cNor γ) (I ×ˢ univ) :=
  contDiffOn_const.mul (contDiffOn_cTan hI hγ hreg)

theorem contDiffOn_cKap : ContDiffOn ℝ ∞ (cKap γ) (I ×ˢ univ) :=
  ((contDiffOn_G2 hI hγ).inner ℝ (contDiffOn_const.mul (contDiffOn_G1 hI hγ))).div
    ((contDiffOn_cSpd hI hγ hreg).pow 3) fun q hq => pow_ne_zero 3 (cSpd_pos hI hγ hreg hq).ne'

theorem contDiffOn_cDs {f : ℝ × ℝ → ℝ} (hf : ContDiffOn ℝ ∞ f (I ×ˢ univ)) :
    ContDiffOn ℝ ∞ (cDs γ f) (I ×ˢ univ) :=
  (contDiffOn_pθ (isOpen_U hI) hf).div (contDiffOn_cSpd hI hγ hreg)
    fun q hq => (cSpd_pos hI hγ hreg hq).ne'

theorem norm_cTan {q : ℝ × ℝ} (hq : q ∈ I ×ˢ (univ : Set ℝ)) : ‖cTan γ q‖ = 1 :=
  norm_unitTangent (hreg q hq)

theorem G1_eq {q : ℝ × ℝ} (hq : q ∈ I ×ˢ (univ : Set ℝ)) :
    pθ γ q = ((cSpd γ q : ℝ) : ℂ) * cTan γ q := by
  have h : ((cSpd γ q : ℝ) : ℂ) ≠ 0 := by exact_mod_cast (cSpd_pos hI hγ hreg hq).ne'
  simp only [cTan]; field_simp

/-! #### Frenet formulas in `θ` -/

theorem pθ_cTan {t θ : ℝ} (ht : t ∈ I) :
    pθ (cTan γ) (t, θ) = ((cKap γ (t, θ) * cSpd γ (t, θ) : ℝ) : ℂ) * cNor γ (t, θ) := by
  have h := hasDerivAt_unitTangent (contDiff_slice hγ ht) (fun x => hreg (t, x) ⟨ht, trivial⟩) θ
  have e : pθ (cTan γ) (t, θ) = Complex.I * cTan γ (t, θ) *
      ((cKap γ (t, θ) * cSpd γ (t, θ) : ℝ) : ℂ) := h.deriv
  rw [e, cNor]; ring

theorem pθ_cNor {t θ : ℝ} (ht : t ∈ I) :
    pθ (cNor γ) (t, θ) = -(((cKap γ (t, θ) * cSpd γ (t, θ) : ℝ) : ℂ) * cTan γ (t, θ)) := by
  have hd : DifferentiableAt ℝ (cTan γ) (t, θ) :=
    (contDiffOn_cTan hI hγ hreg).differentiableAt' (isOpen_U hI) ⟨ht, trivial⟩
  have h := (hasDerivAt_pθ hd).const_mul Complex.I
  have e : pθ (cNor γ) (t, θ) = Complex.I * pθ (cTan γ) (t, θ) := h.deriv
  rw [e, pθ_cTan hI hγ hreg ht, cNor]
  ring_nf; rw [Complex.I_sq]; ring

/-! #### Evolution under a normal velocity `∂ₜγ = V N` -/

variable {V : ℝ × ℝ → ℝ} (hV : ContDiffOn ℝ ∞ V (I ×ˢ univ))
  (hvel : ∀ q ∈ I ×ˢ (univ : Set ℝ), pt γ q = ((V q : ℝ) : ℂ) * cNor γ q)
include hV hvel

theorem pt_G1 {t θ : ℝ} (ht : t ∈ I) :
    pt (pθ γ) (t, θ) = ((pθ V (t, θ) : ℝ) : ℂ) * (Complex.I * cTan γ (t, θ)) +
      ((-(V (t, θ) * (cKap γ (t, θ) * cSpd γ (t, θ))) : ℝ) : ℂ) * cTan γ (t, θ) := by
  have hU := isOpen_U hI
  have hq : ((t, θ) : ℝ × ℝ) ∈ I ×ˢ (univ : Set ℝ) := ⟨ht, trivial⟩
  rw [pt_pθ hU hγ hq, pθ_congr hU hvel hq]
  have hVd := (hasDerivAt_pθ (hV.differentiableAt' hU hq)).ofReal_comp
  have hNd := hasDerivAt_pθ ((contDiffOn_cNor hI hγ hreg).differentiableAt' hU hq)
  have e : pθ (fun q => ((V q : ℝ) : ℂ) * cNor γ q) (t, θ) =
      ((pθ V (t, θ) : ℝ) : ℂ) * cNor γ (t, θ) + ((V (t, θ) : ℝ) : ℂ) * pθ (cNor γ) (t, θ) :=
    (hVd.mul hNd).deriv
  rw [e, pθ_cNor hI hγ hreg ht, cNor]
  push_cast; ring

/-- (E1) `∂ₜ v = −κ V v`. -/
theorem pt_cSpd {t θ : ℝ} (ht : t ∈ I) :
    pt (cSpd γ) (t, θ) = -(cKap γ (t, θ) * V (t, θ) * cSpd γ (t, θ)) := by
  have hU := isOpen_U hI
  have hq : ((t, θ) : ℝ × ℝ) ∈ I ×ˢ (univ : Set ℝ) := ⟨ht, trivial⟩
  have h := hasDerivAt_norm_of_ne
    (hasDerivAt_pt ((contDiffOn_G1 hI hγ).differentiableAt' hU hq)) (hreg _ hq)
  have e : pt (cSpd γ) (t, θ) = ⟪pθ γ (t, θ), pt (pθ γ) (t, θ)⟫ / cSpd γ (t, θ) := h.deriv
  have hv := cSpd_pos hI hγ hreg hq
  rw [e, pt_G1 hI hγ hreg hV hvel ht]
  have hG : pθ γ (t, θ) = ((0 : ℝ) : ℂ) * (Complex.I * cTan γ (t, θ)) +
      ((cSpd γ (t, θ) : ℝ) : ℂ) * cTan γ (t, θ) := by
    rw [G1_eq hI hγ hreg hq]; push_cast; ring
  rw [hG, inner_frame (norm_cTan hI hγ hreg hq)]
  field_simp
  ring

/-- (E3) `∂ₜ T = (∂_s V) N`. -/
theorem pt_cTan {t θ : ℝ} (ht : t ∈ I) :
    pt (cTan γ) (t, θ) = ((cDs γ V (t, θ) : ℝ) : ℂ) * cNor γ (t, θ) := by
  have hU := isOpen_U hI
  have hq : ((t, θ) : ℝ × ℝ) ∈ I ×ˢ (univ : Set ℝ) := ⟨ht, trivial⟩
  have hv := cSpd_pos hI hγ hreg hq
  have hvC : ((cSpd γ (t, θ) : ℝ) : ℂ) ≠ 0 := by exact_mod_cast hv.ne'
  have hG := hasDerivAt_pt ((contDiffOn_G1 hI hγ).differentiableAt' hU hq)
  have hS := (hasDerivAt_pt ((contDiffOn_cSpd hI hγ hreg).differentiableAt' hU hq)).ofReal_comp
  have e : pt (cTan γ) (t, θ) = (pt (pθ γ) (t, θ) * ((cSpd γ (t, θ) : ℝ) : ℂ) -
      pθ γ (t, θ) * ((pt (cSpd γ) (t, θ) : ℝ) : ℂ)) / ((cSpd γ (t, θ) : ℝ) : ℂ) ^ 2 :=
    (hG.div hS hvC).deriv
  rw [e, pt_G1 hI hγ hreg hV hvel ht, pt_cSpd hI hγ hreg hV hvel ht, G1_eq hI hγ hreg hq, cNor, cDs]
  push_cast
  field_simp
  ring

/-- `∂ₜ N = −(∂_s V) T`. -/
theorem pt_cNor {t θ : ℝ} (ht : t ∈ I) :
    pt (cNor γ) (t, θ) = -(((cDs γ V (t, θ) : ℝ) : ℂ) * cTan γ (t, θ)) := by
  have hU := isOpen_U hI
  have hq : ((t, θ) : ℝ × ℝ) ∈ I ×ˢ (univ : Set ℝ) := ⟨ht, trivial⟩
  have h := (hasDerivAt_pt ((contDiffOn_cTan hI hγ hreg).differentiableAt' hU hq)).const_mul Complex.I
  have e : pt (cNor γ) (t, θ) = Complex.I * pt (cTan γ) (t, θ) := h.deriv
  rw [e, pt_cTan hI hγ hreg hV hvel ht, cNor]
  ring_nf; rw [Complex.I_sq]; ring

/-- (E2) the commutator `∂ₜ ∂_s f = ∂_s ∂ₜ f + κ V ∂_s f`. -/
theorem pt_cDs {f : ℝ × ℝ → ℝ} (hf : ContDiffOn ℝ ∞ f (I ×ˢ univ)) {t θ : ℝ} (ht : t ∈ I) :
    pt (cDs γ f) (t, θ) = cDs γ (pt f) (t, θ) + cKap γ (t, θ) * V (t, θ) * cDs γ f (t, θ) := by
  have hU := isOpen_U hI
  have hq : ((t, θ) : ℝ × ℝ) ∈ I ×ˢ (univ : Set ℝ) := ⟨ht, trivial⟩
  have hv := cSpd_pos hI hγ hreg hq
  have hP := hasDerivAt_pt ((contDiffOn_pθ hU hf).differentiableAt' hU hq)
  have hS := hasDerivAt_pt ((contDiffOn_cSpd hI hγ hreg).differentiableAt' hU hq)
  have e : pt (cDs γ f) (t, θ) = (pt (pθ f) (t, θ) * cSpd γ (t, θ) -
      pθ f (t, θ) * pt (cSpd γ) (t, θ)) / cSpd γ (t, θ) ^ 2 := (hP.div hS hv.ne').deriv
  rw [e, pt_pθ hU hf hq, pt_cSpd hI hγ hreg hV hvel ht]
  simp only [cDs]
  field_simp
  ring

/-- (E4) `∂ₜ κ = ∂_s² V + κ² V`. -/
theorem pt_cKap {t θ : ℝ} (ht : t ∈ I) :
    pt (cKap γ) (t, θ) = cDs γ (cDs γ V) (t, θ) + cKap γ (t, θ) ^ 2 * V (t, θ) := by
  have hU := isOpen_U hI
  have hq : ((t, θ) : ℝ × ℝ) ∈ I ×ˢ (univ : Set ℝ) := ⟨ht, trivial⟩
  have hv := cSpd_pos hI hγ hreg hq
  -- `κ = ⟪∂θ T, N⟫ / v` on `I × ℝ`
  have hk : ∀ q ∈ I ×ˢ (univ : Set ℝ), cKap γ q = ⟪pθ (cTan γ) q, cNor γ q⟫ / cSpd γ q := by
    rintro ⟨t', θ'⟩ hq'
    have hv' := cSpd_pos hI hγ hreg hq'
    rw [pθ_cTan hI hγ hreg hq'.1]
    have e1 : ((cKap γ (t', θ') * cSpd γ (t', θ') : ℝ) : ℂ) * cNor γ (t', θ') =
        ((cKap γ (t', θ') * cSpd γ (t', θ') : ℝ) : ℂ) * (Complex.I * cTan γ (t', θ')) +
          ((0 : ℝ) : ℂ) * cTan γ (t', θ') := by rw [cNor]; push_cast; ring
    have e2 : cNor γ (t', θ') = ((1 : ℝ) : ℂ) * (Complex.I * cTan γ (t', θ')) +
          ((0 : ℝ) : ℂ) * cTan γ (t', θ') := by rw [cNor]; push_cast; ring
    rw [e1, e2, inner_frame (norm_cTan hI hγ hreg hq')]
    field_simp
    ring
  have hTsm : ContDiffOn ℝ ∞ (pθ (cTan γ)) (I ×ˢ univ) := contDiffOn_pθ hU (contDiffOn_cTan hI hγ hreg)
  have hA := hasDerivAt_pt (hTsm.differentiableAt' hU hq)
  have hB := hasDerivAt_pt ((contDiffOn_cNor hI hγ hreg).differentiableAt' hU hq)
  have hS := hasDerivAt_pt ((contDiffOn_cSpd hI hγ hreg).differentiableAt' hU hq)
  have e0 : pt (fun q => ⟪pθ (cTan γ) q, cNor γ q⟫ / cSpd γ q) (t, θ) = _ :=
    ((hA.inner ℝ hB).div hS hv.ne').deriv
  have e : pt (cKap γ) (t, θ) = ((⟪pθ (cTan γ) (t, θ), pt (cNor γ) (t, θ)⟫ +
      ⟪pt (pθ (cTan γ)) (t, θ), cNor γ (t, θ)⟫) * cSpd γ (t, θ) -
      ⟪pθ (cTan γ) (t, θ), cNor γ (t, θ)⟫ * pt (cSpd γ) (t, θ)) / cSpd γ (t, θ) ^ 2 := by
    rw [pt_congr hU hk hq, e0]
  -- `∂ₜ ∂θ T = ∂θ ∂ₜ T = ∂θ ((∂_s V) N)`
  have hDsV : ContDiffOn ℝ ∞ (cDs γ V) (I ×ˢ univ) := contDiffOn_cDs hI hγ hreg hV
  have hptT : pt (pθ (cTan γ)) (t, θ) =
      ((cSpd γ (t, θ) * cDs γ (cDs γ V) (t, θ) : ℝ) : ℂ) * (Complex.I * cTan γ (t, θ)) +
        ((-(cDs γ V (t, θ) * (cKap γ (t, θ) * cSpd γ (t, θ))) : ℝ) : ℂ) * cTan γ (t, θ) := by
    rw [pt_pθ hU (contDiffOn_cTan hI hγ hreg) hq,
      pθ_congr (g := fun q => ((cDs γ V q : ℝ) : ℂ) * cNor γ q) hU (fun q hq' => by
        obtain ⟨t', θ'⟩ := q
        exact pt_cTan hI hγ hreg hV hvel (t := t') (θ := θ') hq'.1) hq]
    have h1 := (hasDerivAt_pθ (hDsV.differentiableAt' hU hq)).ofReal_comp
    have h2 := hasDerivAt_pθ ((contDiffOn_cNor hI hγ hreg).differentiableAt' hU hq)
    have e3 : pθ (fun q => ((cDs γ V q : ℝ) : ℂ) * cNor γ q) (t, θ) =
        ((pθ (cDs γ V) (t, θ) : ℝ) : ℂ) * cNor γ (t, θ) +
          ((cDs γ V (t, θ) : ℝ) : ℂ) * pθ (cNor γ) (t, θ) := (h1.mul h2).deriv
    rw [e3, pθ_cNor hI hγ hreg ht, cNor]
    have e4 : pθ (cDs γ V) (t, θ) = cSpd γ (t, θ) * cDs γ (cDs γ V) (t, θ) := by
      simp only [cDs]; field_simp
    rw [e4]; push_cast; ring
  have hN : cNor γ (t, θ) = ((1 : ℝ) : ℂ) * (Complex.I * cTan γ (t, θ)) +
      ((0 : ℝ) : ℂ) * cTan γ (t, θ) := by rw [cNor]; push_cast; ring
  have hT' : pθ (cTan γ) (t, θ) = ((cKap γ (t, θ) * cSpd γ (t, θ) : ℝ) : ℂ) *
      (Complex.I * cTan γ (t, θ)) + ((0 : ℝ) : ℂ) * cTan γ (t, θ) := by
    rw [pθ_cTan hI hγ hreg ht, cNor]; push_cast; ring
  have hN' : pt (cNor γ) (t, θ) = ((0 : ℝ) : ℂ) * (Complex.I * cTan γ (t, θ)) +
      (((-cDs γ V (t, θ)) : ℝ) : ℂ) * cTan γ (t, θ) := by
    rw [pt_cNor hI hγ hreg hV hvel ht]; push_cast; ring
  have hTn := norm_cTan hI hγ hreg hq
  rw [e, hptT, hT', hN', pt_cSpd hI hγ hreg hV hvel ht]
  rw [hN, inner_frame hTn, inner_frame hTn, inner_frame hTn]
  field_simp
  ring

end CurveQ

/-! ### Curve diffusion flow -/

/-- The first energy `E₁ = ∫ κ_s² ds` of a curve, over one period. -/
def energy1 (c : ℝ → ℂ) : ℝ := ∫ θ in (0 : ℝ)..(2 * π), dS c (curv c) θ ^ 2 * ‖deriv c θ‖

/-- The integrand of `d/dt E₁` under curve diffusion, written with the slice at time `t`: with
`k₀ = κ`, `k₁ = κ_s`, `k₂ = κ_ss`, `V = −k₂`,
`2 k₁ (∂_s(∂_s² V + κ² V) + κ V k₁) v + k₁² (−κ V v)`. -/
def rate (c : ℝ → ℂ) (θ : ℝ) : ℝ :=
  2 * dS c (curv c) θ *
      (dS c (fun x => dS c (dS c (fun y => -dS c (dS c (curv c)) y)) x +
          curv c x ^ 2 * -dS c (dS c (curv c)) x) θ +
        curv c θ * -dS c (dS c (curv c)) θ * dS c (curv c) θ) * ‖deriv c θ‖ +
    dS c (curv c) θ ^ 2 * -(curv c θ * -dS c (dS c (curv c)) θ * ‖deriv c θ‖)

/-- **A curve diffusion flow** on the open time set `I`: `γ(t, θ)` smooth on `I × ℝ`, `2π`-periodic in `θ`,
regular, with `∂ₜγ = −κ_ss N`, `N = i γ_θ/|γ_θ|` (the form of `flowVel_eq_neg_kss`). -/
structure IsCDF (I : Set ℝ) (γ : ℝ × ℝ → ℂ) : Prop where
  isOpen : IsOpen I
  smooth : ContDiffOn ℝ ∞ γ (I ×ˢ univ)
  periodic : ∀ t θ, γ (t, θ + 2 * π) = γ (t, θ)
  reg : ∀ t ∈ I, ∀ θ, deriv (fun θ' => γ (t, θ')) θ ≠ 0
  vel : ∀ t ∈ I, ∀ θ, HasDerivAt (fun s => γ (s, θ))
    (-((dS (fun θ' => γ (t, θ')) (dS (fun θ' => γ (t, θ')) (curv (fun θ' => γ (t, θ')))) θ : ℝ) : ℂ) *
      (Complex.I * (deriv (fun θ' => γ (t, θ')) θ / (‖deriv (fun θ' => γ (t, θ')) θ‖ : ℂ)))) t

namespace IsCDF

variable {I : Set ℝ} {γ : ℝ × ℝ → ℂ} (h : IsCDF I γ)
include h

theorem reg' : ∀ q ∈ I ×ˢ (univ : Set ℝ), pθ γ q ≠ 0 := fun q hq => h.reg q.1 hq.1 q.2

/-- The normal speed `V = −κ_ss`. -/
abbrev nspd (γ : ℝ × ℝ → ℂ) (q : ℝ × ℝ) : ℝ := -cDs γ (cDs γ (cKap γ)) q

theorem contDiffOn_nspd : ContDiffOn ℝ ∞ (nspd γ) (I ×ˢ univ) :=
  (contDiffOn_cDs h.isOpen h.smooth h.reg' (contDiffOn_cDs h.isOpen h.smooth h.reg'
    (contDiffOn_cKap h.isOpen h.smooth h.reg'))).neg

theorem pt_eq : ∀ q ∈ I ×ˢ (univ : Set ℝ), pt γ q = ((nspd γ q : ℝ) : ℂ) * cNor γ q := by
  rintro ⟨t, θ⟩ hq
  show deriv (fun s => γ (s, θ)) t = _
  rw [(h.vel t hq.1 θ).deriv]
  simp only [Complex.ofReal_neg, neg_mul]
  rfl

theorem slice_contDiff {t : ℝ} (ht : t ∈ I) : ContDiff ℝ ∞ (fun x => γ (t, x)) :=
  contDiff_slice h.smooth ht

theorem slice_periodic (t : ℝ) : Periodic (fun x => γ (t, x)) (2 * π) := fun θ => h.periodic t θ

/-- **`d/dt ∫ κ_s² ds = ∫ rate`** along a curve diffusion flow. -/
theorem hasDerivAt_energy1 {t : ℝ} (ht : t ∈ I) :
    HasDerivAt (fun s => energy1 (fun x => γ (s, x)))
      (∫ θ in (0 : ℝ)..(2 * π), rate (fun x => γ (t, x)) θ) t := by
  have hI := h.isOpen
  have hγ := h.smooth
  have hreg := h.reg'
  have hU := isOpen_U hI
  have hk := contDiffOn_cKap hI hγ hreg
  have hK1 := contDiffOn_cDs hI hγ hreg hk
  have hS := contDiffOn_cSpd hI hγ hreg
  set Φ : ℝ × ℝ → ℝ := fun q => cDs γ (cKap γ) q ^ 2 * cSpd γ q with hΦ
  have hΦs : ContDiffOn ℝ ∞ Φ (I ×ˢ univ) := (hK1.pow 2).mul hS
  have key := hasDerivAt_integral_of_contDiffOn hI hΦs ht 0 (2 * π)
  have hV := h.contDiffOn_nspd
  have hvel := h.pt_eq
  refine key.congr_deriv ?_
  refine integral_congr fun θ _ => ?_
  have hq : ((t, θ) : ℝ × ℝ) ∈ I ×ˢ (univ : Set ℝ) := ⟨ht, trivial⟩
  have h1 := hasDerivAt_pt (hK1.differentiableAt' hU hq)
  have h2 := hasDerivAt_pt (hS.differentiableAt' hU hq)
  have e : pt Φ (t, θ) = 2 * cDs γ (cKap γ) (t, θ) * pt (cDs γ (cKap γ)) (t, θ) * cSpd γ (t, θ) +
      cDs γ (cKap γ) (t, θ) ^ 2 * pt (cSpd γ) (t, θ) := by
    have e0 : pt Φ (t, θ) = _ := ((h1.pow 2).mul h2).deriv
    rw [e0]
    norm_num <;> ring
  rw [e, pt_cDs hI hγ hreg hV hvel hk ht, pt_cSpd hI hγ hreg hV hvel ht]
  have e2 : cDs γ (pt (cKap γ)) (t, θ) =
      cDs γ (fun q => cDs γ (cDs γ (nspd γ)) q + cKap γ q ^ 2 * nspd γ q) (t, θ) := by
    simp only [cDs]
    rw [pθ_congr (g := fun q => cDs γ (cDs γ (nspd γ)) q + cKap γ q ^ 2 * nspd γ q) hU (fun q hq' => by
      obtain ⟨t', θ'⟩ := q
      exact pt_cKap hI hγ hreg hV hvel (t := t') (θ := θ') hq'.1) hq]
    rfl
  rw [e2]
  rfl

end IsCDF

end CWFlow.Curve
