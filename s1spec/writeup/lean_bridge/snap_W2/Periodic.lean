import CWFlow.Curve.Defs
import Mathlib.MeasureTheory.Integral.IntervalIntegral.IntegrationByParts
import Mathlib.MeasureTheory.Integral.IntervalIntegral.Periodic
import Mathlib.Algebra.QuadraticDiscriminant
import Mathlib.Analysis.SpecialFunctions.Sqrt
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic
import Mathlib.Analysis.InnerProductSpace.Calculus
import Mathlib.Analysis.Calculus.ContDiff.Operations
import Mathlib.Analysis.Calculus.Deriv.Shift

/-!
# Periodic inequalities for the curvature bound (SCOPE W2)

Elementary one-dimensional analysis over one period `[0, 2π]`, with a positive continuous weight `w`
(later `w = |c'|`, so that `w dθ = ds`):

* `exists_le_zero_of_integral_nonpos`: a continuous function with nonpositive integral is `≤ 0` somewhere;
* `sq_integral_mul_le`: Cauchy–Schwarz in `L²(w dθ)` (P2);
* `integral_deriv_mul_periodic`: integration by parts over a period (P1);
* `exists_deriv_eq_zero_of_periodic`: the derivative of a periodic function vanishes somewhere;
* `sq_le_sq_add_of_periodic`: `f(θ)² ≤ f(θ₀)² + 2 ‖f‖ ‖∂_s f‖` (P3, P4);
* `two_le_length_mul_bendingEnergy`: `L · ∫ κ² ds ≥ 2` for a closed regular curve (P5);
* `rate_le`: the real-number lemma (C) closing the energy inequality.
-/

open Set Function intervalIntegral Real
open scoped RealInnerProductSpace ContDiff

noncomputable section

namespace CWFlow.Curve

/-! ### Generalities -/

theorem exists_le_zero_of_integral_nonpos {g : ℝ → ℝ} {a b : ℝ} (hab : a < b) (hg : Continuous g)
    (h : ∫ x in a..b, g x ≤ 0) : ∃ x ∈ Icc a b, g x ≤ 0 := by
  by_contra hcon
  push_neg at hcon
  have hpos : 0 < ∫ x in a..b, g x :=
    intervalIntegral_pos_of_pos_on (hg.intervalIntegrable a b)
      (fun x hx => hcon x (Ioo_subset_Icc_self hx)) hab
  linarith

/-- Cauchy–Schwarz with a nonnegative weight (P2). -/
theorem sq_integral_mul_le {f g w : ℝ → ℝ} {a b : ℝ} (hab : a ≤ b) (hf : Continuous f)
    (hg : Continuous g) (hw : Continuous w) (hw0 : ∀ x, 0 ≤ w x) :
    (∫ x in a..b, f x * g x * w x) ^ 2 ≤
      (∫ x in a..b, f x ^ 2 * w x) * ∫ x in a..b, g x ^ 2 * w x := by
  set A := ∫ x in a..b, f x ^ 2 * w x
  set B := ∫ x in a..b, f x * g x * w x
  set C := ∫ x in a..b, g x ^ 2 * w x
  have key : ∀ l : ℝ, 0 ≤ A * (l * l) + (-2 * B) * l + C := by
    intro l
    have h1 : 0 ≤ ∫ x in a..b, (l * f x - g x) ^ 2 * w x :=
      integral_nonneg hab (fun x _ => mul_nonneg (sq_nonneg _) (hw0 x))
    have h2 : ∫ x in a..b, (l * f x - g x) ^ 2 * w x = A * (l * l) + (-2 * B) * l + C := by
      have e : ∀ x, (l * f x - g x) ^ 2 * w x =
          (l * l) * (f x ^ 2 * w x) + (-2 * l) * (f x * g x * w x) + g x ^ 2 * w x :=
        fun x => by ring
      simp_rw [e]
      rw [intervalIntegral.integral_add, intervalIntegral.integral_add,
        intervalIntegral.integral_const_mul, intervalIntegral.integral_const_mul]
      · ring
      all_goals (apply Continuous.intervalIntegrable; fun_prop)
    linarith
  have := discrim_le_zero key
  unfold discrim at this
  nlinarith

/-- `|∫ f g w| ≤ √(∫ f² w) √(∫ g² w)`. -/
theorem abs_integral_mul_le {f g w : ℝ → ℝ} {a b : ℝ} (hab : a ≤ b) (hf : Continuous f)
    (hg : Continuous g) (hw : Continuous w) (hw0 : ∀ x, 0 ≤ w x) :
    |∫ x in a..b, f x * g x * w x| ≤
      √(∫ x in a..b, f x ^ 2 * w x) * √(∫ x in a..b, g x ^ 2 * w x) := by
  rw [← Real.sqrt_mul (integral_nonneg hab (fun x _ => mul_nonneg (sq_nonneg _) (hw0 x))),
    ← Real.sqrt_sq_eq_abs]
  exact Real.sqrt_le_sqrt (sq_integral_mul_le hab hf hg hw hw0)

/-! ### Periodic functions -/

theorem periodic_deriv {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] {f : ℝ → E} {T : ℝ}
    (h : Periodic f T) : Periodic (deriv f) T := fun x => by
  have e : (fun y => f (y + T)) = f := funext h
  rw [← deriv_comp_add_const, e]

/-- Integration by parts over one period (P1). -/
theorem integral_deriv_mul_periodic {f g : ℝ → ℝ} (hf : ContDiff ℝ 1 f) (hg : ContDiff ℝ 1 g)
    (hfp : Periodic f (2 * π)) (hgp : Periodic g (2 * π)) :
    ∫ θ in (0 : ℝ)..(2 * π), deriv f θ * g θ = -∫ θ in (0 : ℝ)..(2 * π), f θ * deriv g θ := by
  have hf' : Continuous (deriv f) := hf.continuous_deriv le_rfl
  have hg' : Continuous (deriv g) := hg.continuous_deriv le_rfl
  have hfd : Differentiable ℝ f := hf.differentiable (by simp)
  have hgd : Differentiable ℝ g := hg.differentiable (by simp)
  have h := intervalIntegral.integral_mul_deriv_eq_deriv_mul (a := 0) (b := 2 * π) (u := f) (v := g)
    (u' := deriv f) (v' := deriv g) (fun x _ => (hfd x).hasDerivAt) (fun x _ => (hgd x).hasDerivAt)
    (hf'.intervalIntegrable _ _) (hg'.intervalIntegrable _ _)
  have e1 : f (2 * π) = f 0 := by simpa using hfp 0
  have e2 : g (2 * π) = g 0 := by simpa using hgp 0
  rw [e1, e2, sub_self, zero_sub] at h
  rw [h, neg_neg]

/-- A periodic differentiable function has a critical point. -/
theorem exists_deriv_eq_zero_of_periodic {f : ℝ → ℝ} (hf : Differentiable ℝ f)
    (hfp : Periodic f (2 * π)) : ∃ θ, deriv f θ = 0 := by
  obtain ⟨θ, hθ, hmax⟩ := isCompact_Icc.exists_isMaxOn (nonempty_Icc.2 (by positivity : (0 : ℝ) ≤ 2 * π))
    hf.continuous.continuousOn
  refine ⟨θ, IsLocalMax.deriv_eq_zero (Filter.Eventually.of_forall fun x => ?_)⟩
  obtain ⟨y, hy, hfy⟩ : ∃ y ∈ Icc (0 : ℝ) (2 * π), f y = f x := by
    obtain ⟨y, hy, e⟩ := hfp.exists_mem_Ico₀ (by positivity) x
    exact ⟨y, Ico_subset_Icc_self hy, e.symm⟩
  rw [← hfy]
  exact hmax hy

/-- `f(θ)² ≤ f(θ₀)² + 2 ‖f‖_{L²(w)} ‖f'/w‖_{L²(w)}` for a periodic `C¹` function (P3, P4). -/
theorem sq_le_sq_add_of_periodic {f w : ℝ → ℝ} (hf : ContDiff ℝ 1 f) (hw : Continuous w)
    (hw0 : ∀ θ, 0 < w θ) (hfp : Periodic f (2 * π)) (hwp : Periodic w (2 * π)) (θ₀ θ : ℝ) :
    f θ ^ 2 ≤ f θ₀ ^ 2 + 2 * (√(∫ x in (0 : ℝ)..(2 * π), f x ^ 2 * w x) *
      √(∫ x in (0 : ℝ)..(2 * π), (deriv f x / w x) ^ 2 * w x)) := by
  have hf' : Continuous (deriv f) := hf.continuous_deriv le_rfl
  have hfd : Differentiable ℝ f := hf.differentiable (by simp)
  have hwne : ∀ x, w x ≠ 0 := fun x => (hw0 x).ne'
  -- reduce `θ` to `[θ₀, θ₀ + 2π]`
  obtain ⟨θ', hθ', e⟩ := hfp.exists_mem_Ico (by positivity) θ θ₀
  rw [e]
  have hdiv : Continuous fun x => deriv f x / w x := hf'.div hw hwne
  -- the integrand `2 f f'` written with the weight
  have hrw : ∀ x, 2 * (f x * deriv f x) = 2 * (f x * (deriv f x / w x) * w x) := fun x => by
    field_simp [hwne x]
  have hftc : f θ' ^ 2 - f θ₀ ^ 2 = ∫ x in θ₀..θ', 2 * (f x * deriv f x) := by
    rw [integral_eq_sub_of_hasDerivAt (f := fun x => f x * f x)
      (fun x _ => ((hfd x).hasDerivAt.mul (hfd x).hasDerivAt).congr_deriv (by ring))
      (by apply Continuous.intervalIntegrable; fun_prop)]
    ring
  -- absolute values
  set g : ℝ → ℝ := fun x => |f x| * |deriv f x / w x| * w x with hg
  have hgc : Continuous g := by fun_prop
  have hg0 : ∀ x, 0 ≤ g x := fun x => mul_nonneg (mul_nonneg (abs_nonneg _) (abs_nonneg _)) (hw0 x).le
  have hgp : Periodic g (2 * π) := by
    intro x
    have hd : deriv f (x + 2 * π) = deriv f x := periodic_deriv hfp x
    simp only [hg, hfp x, hwp x, hd]
  have h1 : |∫ x in θ₀..θ', 2 * (f x * deriv f x)| ≤ 2 * ∫ x in θ₀..θ', g x := by
    simp_rw [hrw]
    rw [integral_const_mul, abs_mul, abs_two]
    gcongr
    refine (abs_integral_le_integral_abs hθ'.1).trans (le_of_eq ?_)
    congr 1; funext x
    simp only [hg, abs_mul, abs_of_pos (hw0 x)]
  have h2 : ∫ x in θ₀..θ', g x ≤ ∫ x in θ₀..(θ₀ + 2 * π), g x :=
    integral_mono_interval le_rfl hθ'.1 hθ'.2.le (Filter.Eventually.of_forall hg0)
      (hgc.intervalIntegrable _ _)
  have h3 : ∫ x in θ₀..(θ₀ + 2 * π), g x = ∫ x in (0 : ℝ)..(2 * π), g x := by
    have := hgp.intervalIntegral_add_eq θ₀ 0
    simpa [add_comm] using this
  have h4 : ∫ x in (0 : ℝ)..(2 * π), g x ≤
      √(∫ x in (0 : ℝ)..(2 * π), f x ^ 2 * w x) *
        √(∫ x in (0 : ℝ)..(2 * π), (deriv f x / w x) ^ 2 * w x) := by
    have := abs_integral_mul_le (a := 0) (b := 2 * π) (by positivity)
      (continuous_abs.comp hf.continuous) (continuous_abs.comp hdiv) hw (fun x => (hw0 x).le)
    simp only [Function.comp_def, sq_abs] at this
    exact (le_abs_self _).trans this
  have : f θ' ^ 2 - f θ₀ ^ 2 ≤ 2 * (√(∫ x in (0 : ℝ)..(2 * π), f x ^ 2 * w x) *
      √(∫ x in (0 : ℝ)..(2 * π), (deriv f x / w x) ^ 2 * w x)) := by
    rw [hftc]
    refine (le_abs_self _).trans (h1.trans ?_)
    linarith
  linarith

/-! ### (P5) `L · ∫ κ² ds ≥ 2` -/

/-- The plane identity `b |a|² − ⟪a, b⟫ a = ⟪b, i a⟫ i a`. -/
theorem plane_identity (a b : ℂ) :
    b * ((‖a‖ ^ 2 : ℝ) : ℂ) - a * ((⟪a, b⟫ : ℝ) : ℂ) = Complex.I * a * ((⟪b, Complex.I * a⟫ : ℝ) : ℂ) := by
  rw [Complex.sq_norm, Complex.normSq_apply]
  apply Complex.ext <;>
  · simp only [Complex.inner, Complex.mul_re, Complex.mul_im, Complex.I_re, Complex.I_im,
      Complex.sub_re, Complex.sub_im, Complex.ofReal_re, Complex.ofReal_im, Complex.conj_re,
      Complex.conj_im]
    ring

/-- Derivative of the norm of a nonvanishing `ℂ`-valued curve. -/
theorem hasDerivAt_norm_of_ne {u : ℝ → ℂ} {u' : ℂ} {x : ℝ} (hu : HasDerivAt u u' x) (h0 : u x ≠ 0) :
    HasDerivAt (fun y => ‖u y‖) (⟪u x, u'⟫ / ‖u x‖) x := by
  have hn : DifferentiableAt ℝ (fun y => ‖u y‖) x := hu.differentiableAt.norm ℝ h0
  have h1 := hu.norm_sq
  have h2 := hn.hasDerivAt.pow 2
  have := h1.unique h2
  have hpos : 0 < ‖u x‖ := norm_pos_iff.2 h0
  simp only [Nat.cast_ofNat, Nat.add_one_sub_one, pow_one] at this
  have hD : deriv (fun y => ‖u y‖) x = ⟪u x, u'⟫ / ‖u x‖ := by
    rw [eq_div_iff hpos.ne']
    linear_combination (-1 / 2 : ℝ) * this
  rw [← hD]
  exact hn.hasDerivAt

/-- The unit tangent `T = c'/|c'|` and its derivative `T' = κ |c'| · i T`. -/
theorem hasDerivAt_unitTangent {c : ℝ → ℂ} (hc : ContDiff ℝ ∞ c) (hreg : ∀ θ, deriv c θ ≠ 0) (θ : ℝ) :
    HasDerivAt (fun x => deriv c x / ((‖deriv c x‖ : ℝ) : ℂ))
      (Complex.I * (deriv c θ / ((‖deriv c θ‖ : ℝ) : ℂ)) * ((curv c θ * ‖deriv c θ‖ : ℝ) : ℂ)) θ := by
  have hd1 : ContDiff ℝ ∞ (deriv c) := (contDiff_infty_iff_deriv.mp hc).2
  have hu : HasDerivAt (deriv c) (deriv (deriv c) θ) θ :=
    ((hd1.differentiable (by simp)) θ).hasDerivAt
  have hn := hasDerivAt_norm_of_ne hu (hreg θ)
  have hnC : HasDerivAt (fun x => ((‖deriv c x‖ : ℝ) : ℂ))
      (((⟪deriv c θ, deriv (deriv c) θ⟫ / ‖deriv c θ‖ : ℝ) : ℂ)) θ := hn.ofReal_comp
  have hne : ((‖deriv c θ‖ : ℝ) : ℂ) ≠ 0 := by exact_mod_cast norm_ne_zero_iff.2 (hreg θ)
  refine (hu.div hnC hne).congr_deriv ?_
  have key := plane_identity (deriv c θ) (deriv (deriv c) θ)
  have hκ : curv c θ = ⟪deriv (deriv c) θ, Complex.I * deriv c θ⟫ / ‖deriv c θ‖ ^ 3 := rfl
  rw [hκ]
  generalize deriv c θ = a at *
  generalize deriv (deriv c) θ = b at *
  generalize hn : ‖a‖ = n at *
  have hn0 : (n : ℂ) ≠ 0 := hne
  calc (b * (n : ℂ) - a * ((⟪a, b⟫ / n : ℝ) : ℂ)) / (n : ℂ) ^ 2
      = (b * ((n ^ 2 : ℝ) : ℂ) - a * ((⟪a, b⟫ : ℝ) : ℂ)) / (n : ℂ) ^ 3 := by
        push_cast; field_simp
    _ = (Complex.I * a * ((⟪b, Complex.I * a⟫ : ℝ) : ℂ)) / (n : ℂ) ^ 3 := by rw [key]
    _ = _ := by push_cast; field_simp

/-- `‖T‖ = 1`. -/
theorem norm_unitTangent {a : ℂ} (ha : a ≠ 0) : ‖a / ((‖a‖ : ℝ) : ℂ)‖ = 1 := by
  rw [norm_div, Complex.norm_real, Real.norm_eq_abs, abs_norm]
  exact div_self (norm_ne_zero_iff.2 ha)

/-- **(P5)** `L · ∫ κ² ds ≥ 2` for a closed regular smooth curve. -/
theorem two_le_length_mul_bendingEnergy {c : ℝ → ℂ} (hc : ContDiff ℝ ∞ c)
    (hreg : ∀ θ, deriv c θ ≠ 0) (hcp : Periodic c (2 * π)) :
    2 ≤ (∫ θ in (0 : ℝ)..(2 * π), ‖deriv c θ‖) * bendingEnergy c := by
  have hd1 : ContDiff ℝ ∞ (deriv c) := (contDiff_infty_iff_deriv.mp hc).2
  have hd2 : ContDiff ℝ ∞ (deriv (deriv c)) := (contDiff_infty_iff_deriv.mp hd1).2
  have hw : Continuous fun x => ‖deriv c x‖ := hd1.continuous.norm
  have hw0 : ∀ x, 0 < ‖deriv c x‖ := fun x => norm_pos_iff.2 (hreg x)
  have hκ : Continuous (curv c) := by
    unfold curv
    exact (hd2.continuous.inner (continuous_const.mul hd1.continuous)).div (hw.pow 3)
      (fun x => pow_ne_zero 3 (hw0 x).ne')
  set T : ℝ → ℂ := fun x => deriv c x / ((‖deriv c x‖ : ℝ) : ℂ) with hTdef
  set T' : ℝ → ℂ := fun x => Complex.I * T x * ((curv c x * ‖deriv c x‖ : ℝ) : ℂ) with hT'def
  have hT : ∀ x, HasDerivAt T (T' x) x := fun x => hasDerivAt_unitTangent hc hreg x
  have hTn : ∀ x, ‖T x‖ = 1 := fun x => norm_unitTangent (hreg x)
  have hTc : Continuous T := hd1.continuous.div (Complex.continuous_ofReal.comp hw)
    (fun x => by exact_mod_cast (hw0 x).ne')
  have hT'c : Continuous T' := by
    simp only [hT'def]
    exact (continuous_const.mul hTc).mul (Complex.continuous_ofReal.comp (hκ.mul hw))
  have hT'n : ∀ x, ‖T' x‖ = |curv c x| * ‖deriv c x‖ := fun x => by
    simp only [hT'def, norm_mul, Complex.norm_I, hTn, one_mul, Complex.norm_real, Real.norm_eq_abs,
      abs_norm]
  -- `∫ ⟪c', T 0⟫ = 0`
  have hcd : Differentiable ℝ c := hc.differentiable (by simp)
  have hint0 : ∫ x in (0 : ℝ)..(2 * π), ⟪deriv c x, T 0⟫ = 0 := by
    rw [integral_eq_sub_of_hasDerivAt (f := fun x => ⟪c x, T 0⟫)]
    · simp only [hcp.eq, sub_self]
    · intro x _
      exact ((hcd x).hasDerivAt).inner ℝ (hasDerivAt_const x (T 0)) |>.congr_deriv (by simp)
    · apply Continuous.intervalIntegrable
      exact hd1.continuous.inner continuous_const
  obtain ⟨θ₁, hθ₁, hle⟩ := exists_le_zero_of_integral_nonpos (by positivity)
    (hd1.continuous.inner continuous_const) hint0.le
  have hcT : deriv c θ₁ = ((‖deriv c θ₁‖ : ℝ) : ℂ) * T θ₁ := by
    simp only [hTdef]
    field_simp [show ((‖deriv c θ₁‖ : ℝ) : ℂ) ≠ 0 by exact_mod_cast (hw0 θ₁).ne']
  have hinner : ⟪T θ₁, T 0⟫ ≤ 0 := by
    rw [hcT, ← Complex.real_smul, real_inner_smul_left] at hle
    exact nonpos_of_mul_nonpos_right hle (hw0 θ₁) |> fun h => by linarith [h]
  have h2 : 2 ≤ ‖T θ₁ - T 0‖ ^ 2 := by
    rw [norm_sub_sq_real, hTn, hTn]; linarith
  have hftc : T θ₁ - T 0 = ∫ x in (0 : ℝ)..θ₁, T' x :=
    (integral_eq_sub_of_hasDerivAt (fun x _ => hT x) (hT'c.intervalIntegrable _ _)).symm
  set X := ∫ x in (0 : ℝ)..(2 * π), |curv c x| * ‖deriv c x‖ with hX
  have h3 : ‖T θ₁ - T 0‖ ≤ X := by
    rw [hftc]
    refine (norm_integral_le_integral_norm hθ₁.1).trans ?_
    simp_rw [hT'n]
    exact integral_mono_interval (f := fun x => |curv c x| * ‖deriv c x‖) le_rfl hθ₁.1 hθ₁.2
      (Filter.Eventually.of_forall fun x => mul_nonneg (abs_nonneg _) (hw0 x).le)
      (by apply Continuous.intervalIntegrable; exact (continuous_abs.comp hκ).mul hw)
  have h4 : X ^ 2 ≤ (∫ θ in (0 : ℝ)..(2 * π), ‖deriv c θ‖) * bendingEnergy c := by
    have := sq_integral_mul_le (f := fun _ => (1 : ℝ)) (g := fun x => |curv c x|) (a := 0)
      (b := 2 * π) (by positivity) continuous_const
      (continuous_abs.comp hκ) hw (fun x => (hw0 x).le)
    simp only [one_pow, one_mul, sq_abs] at this
    simpa [bendingEnergy, hX] using this
  have h5 : ‖T θ₁ - T 0‖ ^ 2 ≤ X ^ 2 := pow_le_pow_left₀ (norm_nonneg _) h3 2
  linarith

/-! ### The real-number lemma (C) -/

/-- **Lemma (C)**: if `a₁² ≤ A a₂`, `a₂² ≤ a₁ a₃` and `X ≤ −2a₃² + 4A²a₂a₃ + (2/3)a₁³a₂`, then
`X ≤ (14/3)⁶ A¹⁴`. -/
theorem rate_le {A a₁ a₂ a₃ X : ℝ} (hA : 0 ≤ A) (h1 : 0 ≤ a₁) (h2 : 0 ≤ a₂) (h3 : 0 ≤ a₃)
    (h12 : a₁ ^ 2 ≤ A * a₂) (h23 : a₂ ^ 2 ≤ a₁ * a₃)
    (hX : X ≤ -2 * a₃ ^ 2 + 4 * A ^ 2 * a₂ * a₃ + 2 / 3 * a₁ ^ 3 * a₂) :
    X ≤ (14 / 3) ^ 6 * A ^ 14 := by
  -- `a₁³ ≤ A² a₃`
  have e1 : a₁ ^ 3 ≤ A ^ 2 * a₃ := by
    rcases h1.eq_or_lt with h | h
    · rw [← h]; simp; positivity
    · have : a₁ ^ 4 ≤ A ^ 2 * a₂ ^ 2 := by
        have := pow_le_pow_left₀ (by positivity) h12 2
        nlinarith
      have h' : a₁ * a₁ ^ 3 ≤ a₁ * (A ^ 2 * a₃) := by nlinarith
      exact le_of_mul_le_mul_left h' h
  -- `a₂³ ≤ A a₃²`
  have e2 : a₂ ^ 3 ≤ A * a₃ ^ 2 := by
    rcases h2.eq_or_lt with h | h
    · rw [← h]; simp; positivity
    · have : a₂ ^ 4 ≤ a₁ ^ 2 * a₃ ^ 2 := by
        have := pow_le_pow_left₀ (by positivity) h23 2
        nlinarith
      have h' : a₂ * a₂ ^ 3 ≤ a₂ * (A * a₃ ^ 2) := by nlinarith
      exact le_of_mul_le_mul_left h' h
  set c := 14 / 3 * A ^ 2 with hc
  have hc0 : 0 ≤ c := by positivity
  have hX' : X ≤ -2 * a₃ ^ 2 + c * a₂ * a₃ := by
    have : 2 / 3 * a₁ ^ 3 * a₂ ≤ 2 / 3 * (A ^ 2 * a₃) * a₂ := by gcongr
    have hce : c * a₂ * a₃ = 4 * A ^ 2 * a₂ * a₃ + 2 / 3 * (A ^ 2 * a₃) * a₂ := by rw [hc]; ring
    linarith
  have hB : 0 ≤ A ^ 2 * c ^ 6 := by positivity
  have hBeq : A ^ 2 * c ^ 6 = (14 / 3) ^ 6 * A ^ 14 := by rw [hc]; ring
  rw [← hBeq]
  rcases le_or_gt (c * a₂) a₃ with hcase | hcase
  · nlinarith
  · have ha2 : 0 < a₂ := by
      by_contra h; push_neg at h
      have : a₂ = 0 := le_antisymm h h2
      rw [this] at hcase; simp at hcase; linarith
    have hc' : 0 < c := by
      by_contra h; push_neg at h
      have : c = 0 := le_antisymm h hc0
      rw [this] at hcase; simp at hcase; linarith
    -- `a₂ ≤ A c²`
    have hle : a₂ ≤ A * c ^ 2 := by
      have h' : a₂ ^ 2 * a₂ ≤ a₂ ^ 2 * (A * c ^ 2) := by
        have : a₃ ^ 2 ≤ (c * a₂) ^ 2 := pow_le_pow_left₀ h3 hcase.le 2
        nlinarith
      exact le_of_mul_le_mul_left h' (by positivity)
    have h5 : c * a₂ * a₃ ≤ A ^ 2 * c ^ 6 := by
      have hx : c * a₂ ≤ A * c ^ 3 := by nlinarith
      have hy : a₃ ≤ A * c ^ 3 := by linarith
      have : c * a₂ * a₃ ≤ (A * c ^ 3) * (A * c ^ 3) :=
        mul_le_mul hx hy h3 (by positivity)
      nlinarith
    nlinarith

end CWFlow.Curve
