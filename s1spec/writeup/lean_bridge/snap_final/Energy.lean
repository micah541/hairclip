import CWFlow.Curve.Arclength

/-!
# The energy inequality and the curvature bound (SCOPE W4)

For a smooth closed regular curve `c` (one time slice), with `w = |c'|`, `k₀ = κ`, `k_{j+1} = ∂_s k_j`:

* `integral_rate`: the energy identity (I),
  `∫ rate = −2 ∫ k₃² ds + 2 ∫ κ² k₂² ds + (1/3) ∫ k₁⁴ ds`;
* `integral_rate_le`: `∫ rate ≤ (14/3)⁶ K⁷` when `∫ κ² ds ≤ K`;
* `curv_sq_le`: `κ² ≤ K²/2 + K + E₁` when `∫ κ² ds ≤ K`;
* `IsCDF.curv_sq_le` (T5): along a curve diffusion flow with `∫ κ² ds ≤ K` on `[t₀, b)`, `κ²` is bounded on
  `[t₀, b)`.
-/

open Set Function Filter MeasureTheory intervalIntegral Real
open scoped RealInnerProductSpace ContDiff Topology

noncomputable section

namespace CWFlow.Curve

section Slice

variable {c : ℝ → ℂ} (hc : ContDiff ℝ ∞ c) (hreg : ∀ θ, deriv c θ ≠ 0) (hcp : Periodic c (2 * π))
include hc hreg hcp

omit hcp in
theorem contDiff_deriv_c : ContDiff ℝ ∞ (deriv c) := (contDiff_infty_iff_deriv.mp hc).2

omit hcp in
theorem contDiff_w : ContDiff ℝ ∞ (fun x => ‖deriv c x‖) :=
  (contDiff_deriv_c hc hreg).norm ℝ hreg

omit hcp in
theorem w_pos (x : ℝ) : 0 < ‖deriv c x‖ := norm_pos_iff.2 (hreg x)

omit hcp in
theorem contDiff_curv : ContDiff ℝ ∞ (curv c) := by
  have h1 := contDiff_deriv_c hc hreg
  have h2 : ContDiff ℝ ∞ (deriv (deriv c)) := (contDiff_infty_iff_deriv.mp h1).2
  unfold curv
  exact (h2.inner ℝ (contDiff_const.mul h1)).div ((contDiff_w hc hreg).pow 3)
    fun x => pow_ne_zero 3 (w_pos hc hreg x).ne'

omit hcp in
theorem contDiff_dS {f : ℝ → ℝ} (hf : ContDiff ℝ ∞ f) : ContDiff ℝ ∞ (dS c f) := by
  unfold dS
  exact (contDiff_infty_iff_deriv.mp hf).2.div (contDiff_w hc hreg) fun x => (w_pos hc hreg x).ne'

theorem periodic_w : Periodic (fun x => ‖deriv c x‖) (2 * π) := fun x => by
  simp only [periodic_deriv hcp x]

theorem periodic_curv : Periodic (curv c) (2 * π) := fun x => by
  have h1 := periodic_deriv hcp
  have h2 := periodic_deriv h1
  simp only [curv, h1 x, h2 x]

theorem periodic_dS {f : ℝ → ℝ} (hfp : Periodic f (2 * π)) : Periodic (dS c f) (2 * π) := fun x => by
  simp only [dS, periodic_deriv hfp x, periodic_deriv hcp x]

omit hcp in
theorem dS_mul_w (f : ℝ → ℝ) (x : ℝ) : dS c f x * ‖deriv c x‖ = deriv f x := by
  unfold dS; field_simp [(w_pos hc hreg x).ne']

/-- Weighted integration by parts over a period (P1): `∫ (∂_s f) g ds = −∫ f ∂_s g ds`. -/
theorem integral_dS_mul {f g : ℝ → ℝ} (hf : ContDiff ℝ ∞ f) (hg : ContDiff ℝ ∞ g)
    (hfp : Periodic f (2 * π)) (hgp : Periodic g (2 * π)) :
    ∫ x in (0 : ℝ)..(2 * π), dS c f x * g x * ‖deriv c x‖ =
      -∫ x in (0 : ℝ)..(2 * π), f x * dS c g x * ‖deriv c x‖ := by
  have e1 : ∀ x, dS c f x * g x * ‖deriv c x‖ = deriv f x * g x := fun x => by
    rw [mul_right_comm, dS_mul_w hc hreg]
  have e2 : ∀ x, f x * dS c g x * ‖deriv c x‖ = f x * deriv g x := fun x => by
    rw [mul_assoc, dS_mul_w hc hreg]
  simp_rw [e1, e2]
  exact integral_deriv_mul_periodic (hf.of_le (by simp)) (hg.of_le (by simp)) hfp hgp

/-- `f(θ)² ≤ f(θ₀)² + 2 ‖f‖ ‖∂_s f‖`, in `L²(ds)` (P3, P4 for the slice). -/
theorem sq_le_slice {f : ℝ → ℝ} (hf : ContDiff ℝ ∞ f) (hfp : Periodic f (2 * π)) (θ₀ θ : ℝ) :
    f θ ^ 2 ≤ f θ₀ ^ 2 + 2 * (√(∫ x in (0 : ℝ)..(2 * π), f x ^ 2 * ‖deriv c x‖) *
      √(∫ x in (0 : ℝ)..(2 * π), dS c f x ^ 2 * ‖deriv c x‖)) :=
  sq_le_sq_add_of_periodic (w := fun x => ‖deriv c x‖) (hf.of_le (by simp))
    (contDiff_w hc hreg).continuous (w_pos hc hreg) hfp (periodic_w hc hreg hcp) θ₀ θ

/-- `∂_s f` vanishes somewhere. -/
theorem exists_dS_eq_zero {f : ℝ → ℝ} (hf : ContDiff ℝ ∞ f) (hfp : Periodic f (2 * π)) :
    ∃ θ, dS c f θ = 0 := by
  obtain ⟨θ, h⟩ := exists_deriv_eq_zero_of_periodic (hf.differentiable (by simp)) hfp
  exact ⟨θ, by simp [dS, h]⟩

/-- **The energy identity (I)**:
`∫ rate = −2 ∫ κ_sss² ds + 2 ∫ κ² κ_ss² ds + (1/3) ∫ κ_s⁴ ds`. -/
theorem integral_rate :
    ∫ θ in (0 : ℝ)..(2 * π), rate c θ =
      -2 * (∫ θ in (0 : ℝ)..(2 * π), dS c (dS c (dS c (curv c))) θ ^ 2 * ‖deriv c θ‖) +
      2 * (∫ θ in (0 : ℝ)..(2 * π), curv c θ ^ 2 * dS c (dS c (curv c)) θ ^ 2 * ‖deriv c θ‖) +
      1 / 3 * ∫ θ in (0 : ℝ)..(2 * π), dS c (curv c) θ ^ 4 * ‖deriv c θ‖ := by
  set w : ℝ → ℝ := fun x => ‖deriv c x‖ with hw
  set k0 := curv c with hk0
  set k1 := dS c k0 with hk1
  set k2 := dS c k1 with hk2
  set k3 := dS c k2 with hk3
  have s0 : ContDiff ℝ ∞ k0 := contDiff_curv hc hreg
  have s1 : ContDiff ℝ ∞ k1 := contDiff_dS hc hreg s0
  have s2 : ContDiff ℝ ∞ k2 := contDiff_dS hc hreg s1
  have s3 : ContDiff ℝ ∞ k3 := contDiff_dS hc hreg s2
  have p0 : Periodic k0 (2 * π) := periodic_curv hc hreg hcp
  have p1 : Periodic k1 (2 * π) := periodic_dS hc hreg hcp p0
  have p2 : Periodic k2 (2 * π) := periodic_dS hc hreg hcp p1
  have p3 : Periodic k3 (2 * π) := periodic_dS hc hreg hcp p2
  have cw : Continuous w := (contDiff_w hc hreg).continuous
  have c0 := s0.continuous
  have c1 := s1.continuous
  have c2 := s2.continuous
  have c3 := s3.continuous
  -- `V = −k₂`, `H = ∂_s V = −k₃`, `G = ∂_s² V + κ² V`, `P = k₁³/3`
  set Vc : ℝ → ℝ := fun y => -k2 y with hVc
  set H := dS c Vc with hH
  have hHe : ∀ x, H x = -k3 x := fun x => by
    simp only [hH, hVc, hk3, dS, deriv.fun_neg, neg_div]
  have sH : ContDiff ℝ ∞ H := contDiff_dS hc hreg s2.neg
  have pH : Periodic H (2 * π) := fun x => by rw [hHe, hHe, p3 x]
  set G : ℝ → ℝ := fun x => dS c H x + k0 x ^ 2 * Vc x with hG
  have sG : ContDiff ℝ ∞ G := (contDiff_dS hc hreg sH).add ((s0.pow 2).mul s2.neg)
  have pG : Periodic G (2 * π) := fun x => by
    simp only [hG, hVc, periodic_dS hc hreg hcp pH x, p0 x, p2 x]
  set P : ℝ → ℝ := fun x => k1 x ^ 3 / 3 with hP
  have sP : ContDiff ℝ ∞ P := (s1.pow 3).div_const 3
  have pP : Periodic P (2 * π) := fun x => by simp only [hP, p1 x]
  have hPe : ∀ x, dS c P x = k1 x ^ 2 * k2 x := fun x => by
    have hd : deriv P x = k1 x ^ 2 * deriv k1 x := by
      have h := ((s1.differentiable (by simp) x).hasDerivAt.pow 3).div_const 3
      have e : deriv P x = ((3 : ℕ) : ℝ) * k1 x ^ (3 - 1) * deriv k1 x / 3 := h.deriv
      rw [e, show (3 : ℕ) - 1 = 2 from rfl]
      push_cast; ring
    simp only [dS, hd, hk2]
    ring
  -- pointwise form of the integrand
  have hrate : ∀ θ, rate c θ = 2 * (k1 θ * dS c G θ * w θ) - k0 θ * dS c P θ * w θ := fun θ => by
    rw [hPe]; simp only [rate, hG, hH, hVc, hk2, hk1, hk0, hw]; ring
  have ii : ∀ {f : ℝ → ℝ}, Continuous f → IntervalIntegrable f MeasureTheory.volume 0 (2 * π) :=
    fun hf => hf.intervalIntegrable _ _
  have cdG : Continuous (dS c G) := (contDiff_dS hc hreg sG).continuous
  have cdP : Continuous (dS c P) := (contDiff_dS hc hreg sP).continuous
  have cG := sG.continuous
  have cP := sP.continuous
  have cH := sH.continuous
  have cdH : Continuous (dS c H) := (contDiff_dS hc hreg sH).continuous
  -- the three integrations by parts
  have ibp1 : ∫ x in (0 : ℝ)..(2 * π), k1 x * dS c G x * w x =
      -∫ x in (0 : ℝ)..(2 * π), k2 x * G x * w x := by
    have := integral_dS_mul hc hreg hcp s1 sG p1 pG
    rw [this, neg_neg]
  have ibp2 : ∫ x in (0 : ℝ)..(2 * π), k2 x * dS c H x * w x =
      -∫ x in (0 : ℝ)..(2 * π), k3 x * H x * w x := by
    have := integral_dS_mul hc hreg hcp s2 sH p2 pH
    rw [this, neg_neg]
  have ibp3 : ∫ x in (0 : ℝ)..(2 * π), k0 x * dS c P x * w x =
      -∫ x in (0 : ℝ)..(2 * π), k1 x * P x * w x := by
    have := integral_dS_mul hc hreg hcp s0 sP p0 pP
    rw [this, neg_neg]
  have eG : ∫ x in (0 : ℝ)..(2 * π), k2 x * G x * w x =
      (∫ x in (0 : ℝ)..(2 * π), k2 x * dS c H x * w x) -
        ∫ x in (0 : ℝ)..(2 * π), k0 x ^ 2 * k2 x ^ 2 * w x := by
    have hs := intervalIntegral.integral_sub (μ := MeasureTheory.volume)
      (f := fun x => k2 x * dS c H x * w x) (g := fun x => k0 x ^ 2 * k2 x ^ 2 * w x)
      (ii (by fun_prop)) (ii (by fun_prop))
    rw [← hs]
    congr 1; funext x; simp only [hG, hVc]; ring
  have eH : ∫ x in (0 : ℝ)..(2 * π), k3 x * H x * w x =
      -∫ x in (0 : ℝ)..(2 * π), k3 x ^ 2 * w x := by
    rw [← intervalIntegral.integral_neg]
    congr 1; funext x; rw [hHe]; ring
  have eP : ∫ x in (0 : ℝ)..(2 * π), k1 x * P x * w x =
      1 / 3 * ∫ x in (0 : ℝ)..(2 * π), k1 x ^ 4 * w x := by
    rw [← intervalIntegral.integral_const_mul]
    congr 1; funext x; simp only [hP]; ring
  have eR : ∫ θ in (0 : ℝ)..(2 * π), rate c θ =
      2 * (∫ x in (0 : ℝ)..(2 * π), k1 x * dS c G x * w x) -
        ∫ x in (0 : ℝ)..(2 * π), k0 x * dS c P x * w x := by
    have hs := intervalIntegral.integral_sub (μ := MeasureTheory.volume)
      (f := fun x => 2 * (k1 x * dS c G x * w x)) (g := fun x => k0 x * dS c P x * w x)
      (ii (by fun_prop)) (ii (by fun_prop))
    rw [← intervalIntegral.integral_const_mul, ← hs]
    congr 1; funext θ; rw [hrate]
  have e4 : ∫ x in (0 : ℝ)..(2 * π), k0 x ^ 2 * k2 x ^ 2 * w x =
      ∫ θ in (0 : ℝ)..(2 * π), curv c θ ^ 2 * dS c (dS c (curv c)) θ ^ 2 * ‖deriv c θ‖ := rfl
  rw [eR, ibp1, eG, ibp2, eH, ibp3, eP, ← e4]
  simp only [hk3, hk2, hk1, hk0, hw]
  ring

/-- **The energy inequality**: if `∫ κ² ds ≤ K` then `∫ rate ≤ (14/3)⁶ K⁷`. -/
theorem integral_rate_le {K : ℝ} (hK : bendingEnergy c ≤ K) :
    ∫ θ in (0 : ℝ)..(2 * π), rate c θ ≤ (14 / 3) ^ 6 * K ^ 7 := by
  set w : ℝ → ℝ := fun x => ‖deriv c x‖ with hw
  set k0 := curv c with hk0
  set k1 := dS c k0 with hk1
  set k2 := dS c k1 with hk2
  set k3 := dS c k2 with hk3
  have s0 : ContDiff ℝ ∞ k0 := contDiff_curv hc hreg
  have s1 : ContDiff ℝ ∞ k1 := contDiff_dS hc hreg s0
  have s2 : ContDiff ℝ ∞ k2 := contDiff_dS hc hreg s1
  have s3 : ContDiff ℝ ∞ k3 := contDiff_dS hc hreg s2
  have p0 : Periodic k0 (2 * π) := periodic_curv hc hreg hcp
  have p1 : Periodic k1 (2 * π) := periodic_dS hc hreg hcp p0
  have p2 : Periodic k2 (2 * π) := periodic_dS hc hreg hcp p1
  have cw : Continuous w := (contDiff_w hc hreg).continuous
  have hw0 : ∀ x, 0 ≤ w x := fun x => (w_pos hc hreg x).le
  have c0 := s0.continuous
  have c1 := s1.continuous
  have c2 := s2.continuous
  have c3 := s3.continuous
  have ii : ∀ {f : ℝ → ℝ}, Continuous f → IntervalIntegrable f MeasureTheory.volume 0 (2 * π) :=
    fun hf => hf.intervalIntegrable _ _
  have hpi : (0 : ℝ) ≤ 2 * π := by positivity
  have nn : ∀ f : ℝ → ℝ, 0 ≤ ∫ x in (0 : ℝ)..(2 * π), f x ^ 2 * w x := fun f =>
    intervalIntegral.integral_nonneg hpi fun x _ => mul_nonneg (sq_nonneg _) (hw0 x)
  set P0 := ∫ x in (0 : ℝ)..(2 * π), k0 x ^ 2 * w x with hP0
  set P1 := ∫ x in (0 : ℝ)..(2 * π), k1 x ^ 2 * w x with hP1
  set P2 := ∫ x in (0 : ℝ)..(2 * π), k2 x ^ 2 * w x with hP2
  set P3 := ∫ x in (0 : ℝ)..(2 * π), k3 x ^ 2 * w x with hP3
  set A := √P0
  set a₁ := √P1
  set a₂ := √P2
  set a₃ := √P3
  have hA2 : A ^ 2 = P0 := Real.sq_sqrt (nn k0)
  have h12 : a₁ ^ 2 = P1 := Real.sq_sqrt (nn k1)
  have h22 : a₂ ^ 2 = P2 := Real.sq_sqrt (nn k2)
  have h32 : a₃ ^ 2 = P3 := Real.sq_sqrt (nn k3)
  -- `a₁² ≤ A a₂` and `a₂² ≤ a₁ a₃`
  have i12 : a₁ ^ 2 ≤ A * a₂ := by
    have e : P1 = -∫ x in (0 : ℝ)..(2 * π), k0 x * k2 x * w x := by
      have := integral_dS_mul hc hreg hcp s0 s1 p0 p1
      rw [← this]; simp only [hP1, hk1, sq]; rfl
    have hb := abs_integral_mul_le (a := 0) (b := 2 * π) hpi c0 c2 cw hw0
    rw [h12, e]
    exact (neg_le_abs _).trans hb
  have i23 : a₂ ^ 2 ≤ a₁ * a₃ := by
    have e : P2 = -∫ x in (0 : ℝ)..(2 * π), k1 x * k3 x * w x := by
      have := integral_dS_mul hc hreg hcp s1 s2 p1 p2
      rw [← this]; simp only [hP2, hk2, sq]; rfl
    have hb := abs_integral_mul_le (a := 0) (b := 2 * π) hpi c1 c3 cw hw0
    rw [h22, e]
    exact (neg_le_abs _).trans hb
  -- sup bounds `k₂² ≤ 2 a₂ a₃`, `k₁² ≤ 2 a₁ a₂`
  obtain ⟨z2, hz2⟩ := exists_dS_eq_zero hc hreg hcp s1 p1
  obtain ⟨z1, hz1⟩ := exists_dS_eq_zero hc hreg hcp s0 p0
  have sup2 : ∀ x, k2 x ^ 2 ≤ 2 * (a₂ * a₃) := fun x => by
    have := sq_le_slice hc hreg hcp s2 p2 z2 x
    have h0 : k2 z2 = 0 := hz2
    calc k2 x ^ 2 ≤ k2 z2 ^ 2 + 2 * (a₂ * a₃) := this
      _ = 2 * (a₂ * a₃) := by rw [h0]; ring
  have sup1 : ∀ x, k1 x ^ 2 ≤ 2 * (a₁ * a₂) := fun x => by
    have := sq_le_slice hc hreg hcp s1 p1 z1 x
    have h0 : k1 z1 = 0 := hz1
    calc k1 x ^ 2 ≤ k1 z1 ^ 2 + 2 * (a₁ * a₂) := this
      _ = 2 * (a₁ * a₂) := by rw [h0]; ring
  have hQ : ∫ x in (0 : ℝ)..(2 * π), k0 x ^ 2 * k2 x ^ 2 * w x ≤ 2 * (a₂ * a₃) * P0 := by
    rw [hP0, ← intervalIntegral.integral_const_mul]
    refine intervalIntegral.integral_mono_on hpi (ii (by fun_prop)) (ii (by fun_prop)) fun x _ => ?_
    have := mul_le_mul_of_nonneg_left (sup2 x) (mul_nonneg (sq_nonneg (k0 x)) (hw0 x))
    nlinarith [this]
  have hR : ∫ x in (0 : ℝ)..(2 * π), k1 x ^ 4 * w x ≤ 2 * (a₁ * a₂) * P1 := by
    rw [hP1, ← intervalIntegral.integral_const_mul]
    refine intervalIntegral.integral_mono_on hpi (ii (by fun_prop)) (ii (by fun_prop)) fun x _ => ?_
    have := mul_le_mul_of_nonneg_left (sup1 x) (mul_nonneg (sq_nonneg (k1 x)) (hw0 x))
    nlinarith [this]
  have hI := integral_rate hc hreg hcp
  have hX : ∫ θ in (0 : ℝ)..(2 * π), rate c θ ≤
      -2 * a₃ ^ 2 + 4 * A ^ 2 * a₂ * a₃ + 2 / 3 * a₁ ^ 3 * a₂ := by
    rw [hI]
    have e3 : ∫ θ in (0 : ℝ)..(2 * π), dS c (dS c (dS c (curv c))) θ ^ 2 * ‖deriv c θ‖ = a₃ ^ 2 := by
      rw [h32]
    rw [e3]
    have hQ' : ∫ θ in (0 : ℝ)..(2 * π), curv c θ ^ 2 * dS c (dS c (curv c)) θ ^ 2 * ‖deriv c θ‖
        ≤ 2 * (a₂ * a₃) * A ^ 2 := by rw [hA2]; exact hQ
    have hR' : ∫ θ in (0 : ℝ)..(2 * π), dS c (curv c) θ ^ 4 * ‖deriv c θ‖ ≤ 2 * (a₁ * a₂) * a₁ ^ 2 := by
      rw [h12]; exact hR
    nlinarith [hQ', hR']
  have hfin := rate_le (Real.sqrt_nonneg _) (Real.sqrt_nonneg _) (Real.sqrt_nonneg _)
    (Real.sqrt_nonneg _) i12 i23 hX
  have hA14 : A ^ 14 ≤ K ^ 7 := by
    have : A ^ 14 = P0 ^ 7 := by rw [← hA2]; ring
    rw [this]
    exact pow_le_pow_left₀ (nn k0) hK 7
  calc _ ≤ (14 / 3) ^ 6 * A ^ 14 := hfin
    _ ≤ (14 / 3) ^ 6 * K ^ 7 := by gcongr

/-- **The curvature bound**: if `∫ κ² ds ≤ K` then `κ² ≤ K²/2 + K + E₁`. -/
theorem curv_sq_le {K : ℝ} (hK : bendingEnergy c ≤ K) (θ : ℝ) :
    curv c θ ^ 2 ≤ K ^ 2 / 2 + K + energy1 c := by
  set w : ℝ → ℝ := fun x => ‖deriv c x‖ with hw
  have s0 : ContDiff ℝ ∞ (curv c) := contDiff_curv hc hreg
  have p0 : Periodic (curv c) (2 * π) := periodic_curv hc hreg hcp
  have cw : Continuous w := (contDiff_w hc hreg).continuous
  have hw0 : ∀ x, 0 < w x := w_pos hc hreg
  have hpi : (0 : ℝ) < 2 * π := by positivity
  set P0 := bendingEnergy c with hP0
  set E1 := energy1 c with hE1
  set L := ∫ x in (0 : ℝ)..(2 * π), w x with hL
  have hP0nn : 0 ≤ P0 := intervalIntegral.integral_nonneg hpi.le fun x _ =>
    mul_nonneg (sq_nonneg _) (hw0 x).le
  have hE1nn : 0 ≤ E1 := intervalIntegral.integral_nonneg hpi.le fun x _ =>
    mul_nonneg (sq_nonneg _) (hw0 x).le
  have hLpos : 0 < L := intervalIntegral.intervalIntegral_pos_of_pos_on (cw.intervalIntegrable _ _)
    (fun x _ => hw0 x) hpi
  -- (P5)
  have h5 : 2 ≤ L * P0 := two_le_length_mul_bendingEnergy hc hreg hcp
  -- (P4): a point where `κ² ≤ P0 / L`
  have hint : ∫ x in (0 : ℝ)..(2 * π), (curv c x ^ 2 - P0 / L) * w x ≤ 0 := by
    have e : ∫ x in (0 : ℝ)..(2 * π), (curv c x ^ 2 - P0 / L) * w x =
        P0 - P0 / L * L := by
      have h1 : ∀ x, (curv c x ^ 2 - P0 / L) * w x = curv c x ^ 2 * w x - P0 / L * w x :=
        fun x => by ring
      simp_rw [h1]
      rw [intervalIntegral.integral_sub, intervalIntegral.integral_const_mul]
      · rfl
      · exact ((s0.continuous.pow 2).mul cw).intervalIntegrable _ _
      · exact (continuous_const.mul cw).intervalIntegrable _ _
    rw [e, div_mul_cancel₀ _ hLpos.ne', sub_self]
  obtain ⟨θ₀, -, hθ₀⟩ := exists_le_zero_of_integral_nonpos hpi
    (((s0.continuous.pow 2).sub continuous_const).mul cw) hint
  have hk0 : curv c θ₀ ^ 2 ≤ P0 / L := by
    have h' : (curv c θ₀ ^ 2 - P0 / L) * w θ₀ ≤ 0 := hθ₀
    by_contra hcon
    push_neg at hcon
    have := mul_pos (sub_pos.2 hcon) (hw0 θ₀)
    linarith
  have hsl := sq_le_slice hc hreg hcp s0 p0 θ₀ θ
  have hA : √P0 ^ 2 = P0 := Real.sq_sqrt hP0nn
  have ha : √E1 ^ 2 = E1 := Real.sq_sqrt hE1nn
  have hAa : 2 * (√P0 * √E1) ≤ P0 + E1 := by
    nlinarith [sq_nonneg (√P0 - √E1)]
  have hPL : P0 / L ≤ P0 ^ 2 / 2 := by
    rw [div_le_iff₀ hLpos]
    nlinarith
  have hPK : P0 ^ 2 ≤ K ^ 2 := pow_le_pow_left₀ hP0nn hK 2
  have : curv c θ ^ 2 ≤ P0 / L + 2 * (√P0 * √E1) := by
    calc curv c θ ^ 2 ≤ curv c θ₀ ^ 2 + 2 * (√P0 * √E1) := hsl
      _ ≤ P0 / L + 2 * (√P0 * √E1) := by linarith
  linarith

end Slice

/-! ### (T5) the curvature bound along curve diffusion -/

/-- **(T5)** Along a curve diffusion flow on `(a, b)`, a bound on `∫ κ² ds` on `[t₀, b)` gives a bound on
`κ²` on `[t₀, b)`. -/
theorem IsCDF.curv_sq_le {a b : ℝ} {γ : ℝ × ℝ → ℂ} (h : IsCDF (Ioo a b) γ) {t₀ K : ℝ}
    (ht₀ : t₀ ∈ Ioo a b) (hK : ∀ t ∈ Ico t₀ b, bendingEnergy (fun θ => γ (t, θ)) ≤ K) :
    ∃ C : ℝ, ∀ t ∈ Ico t₀ b, ∀ θ, curv (fun θ' => γ (t, θ')) θ ^ 2 ≤ C := by
  set E : ℝ → ℝ := fun s => energy1 (fun x => γ (s, x)) with hE
  set Cr := (14 / 3 : ℝ) ^ 6 * K ^ 7 with hCr
  have hsub : Ico t₀ b ⊆ Ioo a b := fun t ht => ⟨ht₀.1.trans_le ht.1, ht.2⟩
  have hder : ∀ t ∈ Ico t₀ b, HasDerivAt E (∫ θ in (0 : ℝ)..(2 * π), rate (fun x => γ (t, x)) θ) t :=
    fun t ht => h.hasDerivAt_energy1 (hsub ht)
  have hrate : ∀ t ∈ Ico t₀ b, ∫ θ in (0 : ℝ)..(2 * π), rate (fun x => γ (t, x)) θ ≤ Cr :=
    fun t ht => integral_rate_le (h.slice_contDiff (hsub ht)) (h.reg t (hsub ht))
      (h.slice_periodic t) (hK t ht)
  have hmvt : ∀ t ∈ Ico t₀ b, E t - E t₀ ≤ Cr * (t - t₀) := by
    intro t ht
    have := Convex.image_sub_le_mul_sub_of_deriv_le (convex_Ico t₀ b) (f := E) (C := Cr)
      (fun x hx => (hder x hx).continuousAt.continuousWithinAt)
      (fun x hx => (hder x (interior_subset hx)).differentiableAt.differentiableWithinAt)
      (fun x hx => by
        have hx' := interior_subset hx
        rw [(hder x hx').deriv]; exact hrate x hx')
      t₀ ⟨le_rfl, ht₀.2⟩ t ht ht.1
    exact this
  refine ⟨K ^ 2 / 2 + K + (E t₀ + max Cr 0 * (b - t₀)), fun t ht θ => ?_⟩
  have h1 := CWFlow.Curve.curv_sq_le (h.slice_contDiff (hsub ht)) (h.reg t (hsub ht)) (h.slice_periodic t)
    (hK t ht) θ
  have h2 : E t ≤ E t₀ + max Cr 0 * (b - t₀) := by
    have := hmvt t ht
    have hle : Cr * (t - t₀) ≤ max Cr 0 * (b - t₀) := by
      calc Cr * (t - t₀) ≤ max Cr 0 * (t - t₀) :=
            mul_le_mul_of_nonneg_right (le_max_left _ _) (by linarith [ht.1])
        _ ≤ max Cr 0 * (b - t₀) :=
            mul_le_mul_of_nonneg_left (by linarith [ht.2]) (le_max_right _ _)
    linarith
  have h3 : energy1 (fun x => γ (t, x)) = E t := rfl
  linarith

end CWFlow.Curve
