import Mathlib.Analysis.InnerProductSpace.PiL2
import Mathlib.Analysis.InnerProductSpace.Calculus
import Mathlib.Analysis.Calculus.LineDeriv.IntegrationByParts
import Mathlib.Analysis.SpecialFunctions.Pow.Deriv
import Mathlib.Analysis.SpecialFunctions.Pow.Continuity
import Mathlib.MeasureTheory.Integral.MeanInequalities
import Mathlib.MeasureTheory.Measure.Haar.InnerProductSpace
import Mathlib.Tactic

/-!
# Hamilton's interpolation step on `ℝⁿ`

Hamilton, *Three-manifolds with positive Ricci curvature*, Theorem 12.1, in `ℝⁿ`: for a `C²`
compactly supported `η : ℝⁿ → W` and `s ≥ 1`,
`∫ |Dη|^{2s} ≤ (n + 2s - 2) ∫ |η| |D²η| |Dη|^{2s-2}`
(`integral_gradSq_rpow_le`), and the Hölder forms (`hamilton_step`, `hamilton_step_sup`).

Here `gradSq η = |Dη|²` and `hessSq η = |D²η|²` are computed in the standard coordinates.
-/

open MeasureTheory Filter Topology

namespace CWFlow.Interp

variable {n : ℕ} {W : Type*} [NormedAddCommGroup W] [InnerProductSpace ℝ W] [FiniteDimensional ℝ W]


/-- `|Dη|² = Σᵢ |∂ᵢη|²`. -/
noncomputable def gradSq (η : EuclideanSpace ℝ (Fin n) → W) (x : EuclideanSpace ℝ (Fin n)) : ℝ :=
  ∑ i, ‖fderiv ℝ η x (EuclideanSpace.single i 1)‖ ^ 2

/-- `|D²η|² = Σᵢₖ |∂ᵢ∂ₖη|²`. -/
noncomputable def hessSq (η : EuclideanSpace ℝ (Fin n) → W) (x : EuclideanSpace ℝ (Fin n)) : ℝ :=
  ∑ i, ∑ k, ‖fderiv ℝ (fun y => fderiv ℝ η y (EuclideanSpace.single k 1)) x
    (EuclideanSpace.single i 1)‖ ^ 2

/-- The partial derivative `∂ᵢη` as a function. -/
noncomputable def pd (η : EuclideanSpace ℝ (Fin n) → W) (i : Fin n) : EuclideanSpace ℝ (Fin n) → W :=
  fun y => fderiv ℝ η y (EuclideanSpace.single i 1)

/-- The second partial derivative `∂ᵢ∂ₖη`. -/
noncomputable def pd2 (η : EuclideanSpace ℝ (Fin n) → W) (i k : Fin n) : EuclideanSpace ℝ (Fin n) → W :=
  fun x => fderiv ℝ (pd η k) x (EuclideanSpace.single i 1)

lemma gradSq_eq (η : EuclideanSpace ℝ (Fin n) → W) (x : EuclideanSpace ℝ (Fin n)) : gradSq η x = ∑ i, ‖pd η i x‖ ^ 2 := rfl

lemma hessSq_eq (η : EuclideanSpace ℝ (Fin n) → W) (x : EuclideanSpace ℝ (Fin n)) : hessSq η x = ∑ i, ∑ k, ‖pd2 η i k x‖ ^ 2 := rfl

/-! ### Algebraic lemmas -/

lemma norm_sum_le_sqrt_mul {ι : Type*} [Fintype ι] (H : ι → W) :
    ‖∑ i, H i‖ ≤ Real.sqrt (Fintype.card ι) * Real.sqrt (∑ i, ‖H i‖ ^ 2) := by
  calc ‖∑ i, H i‖ ≤ ∑ i, ‖H i‖ := norm_sum_le _ _
    _ = ∑ i, (1 : ℝ) * ‖H i‖ := by simp
    _ ≤ Real.sqrt (∑ _i : ι, (1 : ℝ) ^ 2) * Real.sqrt (∑ i, ‖H i‖ ^ 2) :=
        Real.sum_mul_le_sqrt_mul_sqrt _ _ _
    _ = _ := by simp

lemma abs_double_sum_le {ι : Type*} [Fintype ι] (η : W) (D : ι → W) (H : ι → ι → W) :
    |∑ i, inner ℝ η (D i) * ∑ k, inner ℝ (D k) (H i k)| ≤
      ‖η‖ * (∑ i, ‖D i‖ ^ 2) * Real.sqrt (∑ i, ∑ k, ‖H i k‖ ^ 2) := by
  have h0 : ∑ i, inner ℝ η (D i) * ∑ k, inner ℝ (D k) (H i k)
      = ∑ i, ∑ k, inner ℝ η (D i) * inner ℝ (D k) (H i k) := by
    simp [Finset.mul_sum]
  have h1 : |∑ i, ∑ k, inner ℝ η (D i) * inner ℝ (D k) (H i k)| ≤
      ∑ i, ∑ k, ‖η‖ * ((‖D i‖ * ‖D k‖) * ‖H i k‖) := by
    refine (Finset.abs_sum_le_sum_abs _ _).trans (Finset.sum_le_sum fun i _ => ?_)
    refine (Finset.abs_sum_le_sum_abs _ _).trans (Finset.sum_le_sum fun k _ => ?_)
    rw [abs_mul]
    have a1 := abs_real_inner_le_norm η (D i)
    have a2 := abs_real_inner_le_norm (D k) (H i k)
    calc |inner ℝ η (D i)| * |inner ℝ (D k) (H i k)|
        ≤ (‖η‖ * ‖D i‖) * (‖D k‖ * ‖H i k‖) :=
          mul_le_mul a1 a2 (abs_nonneg _) (by positivity)
      _ = _ := by ring
  have h2 : ∑ i, ∑ k, ‖η‖ * ((‖D i‖ * ‖D k‖) * ‖H i k‖)
      = ‖η‖ * ∑ p : ι × ι, (‖D p.1‖ * ‖D p.2‖) * ‖H p.1 p.2‖ := by
    rw [Fintype.sum_prod_type, Finset.mul_sum]
    simp [Finset.mul_sum]
  have h3 : ∑ p : ι × ι, (‖D p.1‖ * ‖D p.2‖) * ‖H p.1 p.2‖ ≤
      Real.sqrt (∑ p : ι × ι, (‖D p.1‖ * ‖D p.2‖) ^ 2) *
        Real.sqrt (∑ p : ι × ι, ‖H p.1 p.2‖ ^ 2) :=
    Real.sum_mul_le_sqrt_mul_sqrt _ _ _
  have h4 : ∑ p : ι × ι, (‖D p.1‖ * ‖D p.2‖) ^ 2 = (∑ i, ‖D i‖ ^ 2) ^ 2 := by
    rw [Fintype.sum_prod_type, sq, Finset.sum_mul_sum]
    refine Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun k _ => ?_
    ring
  have h5 : ∑ p : ι × ι, ‖H p.1 p.2‖ ^ 2 = ∑ i, ∑ k, ‖H i k‖ ^ 2 := by
    rw [Fintype.sum_prod_type]
  rw [h4, h5, Real.sqrt_sq (Finset.sum_nonneg fun _ _ => by positivity)] at h3
  rw [h0]
  calc _ ≤ _ := h1
    _ = _ := h2
    _ ≤ ‖η‖ * ((∑ i, ‖D i‖ ^ 2) * Real.sqrt (∑ i, ∑ k, ‖H i k‖ ^ 2)) :=
        mul_le_mul_of_nonneg_left h3 (norm_nonneg _)
    _ = _ := by ring

/-! ### Regularity -/

section Reg

variable {η : EuclideanSpace ℝ (Fin n) → W}

lemma pd_contDiff (hη : ContDiff ℝ 2 η) (i : Fin n) : ContDiff ℝ 1 (pd η i) :=
  (hη.fderiv_right (m := 1) (by norm_num)).clm_apply contDiff_const

lemma pd2_continuous (hη : ContDiff ℝ 2 η) (i k : Fin n) : Continuous (pd2 η i k) :=
  (((pd_contDiff hη k).continuous_fderiv (by norm_num))).clm_apply continuous_const

lemma gradSq_contDiff (hη : ContDiff ℝ 2 η) : ContDiff ℝ 1 (gradSq η) :=
  ContDiff.sum fun i _ => ContDiff.norm_sq ℝ (pd_contDiff hη i)

lemma hessSq_continuous (hη : ContDiff ℝ 2 η) : Continuous (hessSq η) :=
  continuous_finsetSum _ fun i _ => continuous_finsetSum _ fun k _ =>
    ((pd2_continuous hη i k).norm).pow 2

lemma pd_eq_zero (hc : HasCompactSupport η) {i : Fin n} {x : EuclideanSpace ℝ (Fin n)} (hx : x ∉ tsupport η) :
    pd η i x = 0 := by
  simp [pd, fderiv_of_notMem_tsupport ℝ hx]

lemma tsupport_pd_subset (i : Fin n) : tsupport (pd η i) ⊆ tsupport η :=
  tsupport_fderiv_apply_subset ℝ _

lemma pd2_eq_zero {i k : Fin n} {x : EuclideanSpace ℝ (Fin n)} (hx : x ∉ tsupport η) : pd2 η i k x = 0 := by
  have : x ∉ tsupport (pd η k) := fun h => hx (tsupport_pd_subset k h)
  simp [pd2, fderiv_of_notMem_tsupport ℝ this]

lemma gradSq_eq_zero {x : EuclideanSpace ℝ (Fin n)} (hc : HasCompactSupport η) (hx : x ∉ tsupport η) :
    gradSq η x = 0 := by
  simp [gradSq_eq, pd_eq_zero hc hx]

lemma hessSq_eq_zero {x : EuclideanSpace ℝ (Fin n)} (hx : x ∉ tsupport η) : hessSq η x = 0 := by
  simp [hessSq_eq, pd2_eq_zero hx]

lemma gradSq_nonneg (η : EuclideanSpace ℝ (Fin n) → W) (x : EuclideanSpace ℝ (Fin n)) : 0 ≤ gradSq η x :=
  Finset.sum_nonneg fun _ _ => by positivity

lemma hessSq_nonneg (η : EuclideanSpace ℝ (Fin n) → W) (x : EuclideanSpace ℝ (Fin n)) : 0 ≤ hessSq η x :=
  Finset.sum_nonneg fun _ _ => Finset.sum_nonneg fun _ _ => by positivity

lemma hasFDerivAt_gradSq (hη : ContDiff ℝ 2 η) (x : EuclideanSpace ℝ (Fin n)) (j : Fin n) :
    fderiv ℝ (gradSq η) x (EuclideanSpace.single j 1) =
      ∑ k, 2 * inner ℝ (pd η k x) (pd2 η j k x) := by
  have hd : ∀ k, DifferentiableAt ℝ (pd η k) x := fun k =>
    ((pd_contDiff hη k).differentiable (by norm_num)) x
  have h := HasFDerivAt.fun_sum (u := Finset.univ) (x := x)
    (A := fun k y => ‖pd η k y‖ ^ 2)
    (fun k _ => ((hd k).hasFDerivAt).norm_sq)
  have h2 : gradSq η = fun y => ∑ k, ‖pd η k y‖ ^ 2 := rfl
  rw [h2, h.fderiv]
  simp [pd2]

end Reg

/-! ### Integration by parts and the regularised inequality -/

section Reg2

variable {η : EuclideanSpace ℝ (Fin n) → W}

lemma integrable_fderiv_apply {G : EuclideanSpace ℝ (Fin n) → ℝ} (hG : ContDiff ℝ 1 G)
    (hc : HasCompactSupport G) (v : EuclideanSpace ℝ (Fin n)) :
    Integrable (fun x => fderiv ℝ G x v) :=
  ((hG.continuous_fderiv (by norm_num)).clm_apply continuous_const).integrable_of_hasCompactSupport
    (hc.fderiv_apply ℝ v)

lemma integral_fderiv_apply_eq_zero {G : EuclideanSpace ℝ (Fin n) → ℝ} (hG : ContDiff ℝ 1 G)
    (hc : HasCompactSupport G) (v : EuclideanSpace ℝ (Fin n)) :
    ∫ x, fderiv ℝ G x v = 0 := by
  have hint := integrable_fderiv_apply hG hc v
  have hGint : Integrable G := hG.continuous.integrable_of_hasCompactSupport hc
  have := integral_mul_fderiv_eq_neg_fderiv_mul_of_integrable (μ := volume)
    (f := fun _ => (1 : ℝ)) (g := G) (v := v) (by simp) (by simpa using hint)
    (by simpa using hGint) (fun x _ => differentiableAt_const _)
    (fun x _ => (hG.differentiable (by norm_num)) x)
  simpa using this

/-- Regularised weight `(|Dη|² + ε)^(s-1)`. -/
noncomputable def Phi (η : EuclideanSpace ℝ (Fin n) → W) (s ε : ℝ)
    (x : EuclideanSpace ℝ (Fin n)) : ℝ :=
  (gradSq η x + ε) ^ (s - 1)

lemma gradSq_add_pos (η : EuclideanSpace ℝ (Fin n) → W) {ε : ℝ} (hε : 0 < ε)
    (x : EuclideanSpace ℝ (Fin n)) : 0 < gradSq η x + ε := by
  have := gradSq_nonneg η x
  linarith

lemma Phi_contDiff (hη : ContDiff ℝ 2 η) {ε : ℝ} (hε : 0 < ε) (s : ℝ) :
    ContDiff ℝ 1 (Phi η s ε) :=
  ((gradSq_contDiff hη).add contDiff_const).rpow_const_of_ne
    (fun x => (gradSq_add_pos η hε x).ne')

lemma Phi_nonneg (η : EuclideanSpace ℝ (Fin n) → W) {ε : ℝ} (hε : 0 < ε) (s : ℝ)
    (x : EuclideanSpace ℝ (Fin n)) : 0 ≤ Phi η s ε x :=
  Real.rpow_nonneg (gradSq_add_pos η hε x).le _

lemma fderiv_Phi (hη : ContDiff ℝ 2 η) {ε : ℝ} (hε : 0 < ε) (s : ℝ)
    (x : EuclideanSpace ℝ (Fin n)) (j : Fin n) :
    fderiv ℝ (Phi η s ε) x (EuclideanSpace.single j 1) =
      (s - 1) * (gradSq η x + ε) ^ (s - 1 - 1) *
        ∑ k, 2 * inner ℝ (pd η k x) (pd2 η j k x) := by
  have hg : HasFDerivAt (gradSq η) (fderiv ℝ (gradSq η) x) x :=
    (((gradSq_contDiff hη).differentiable (by norm_num)) x).hasFDerivAt
  have h := (hg.add_const ε).rpow_const (p := s - 1) (Or.inl (gradSq_add_pos η hε x).ne')
  have h2 : Phi η s ε = fun y => (gradSq η y + ε) ^ (s - 1) := rfl
  rw [h2, h.fderiv]
  simp [hasFDerivAt_gradSq hη x j]

/-- The vector field component whose divergence-type integral vanishes. -/
noncomputable def Gi (η : EuclideanSpace ℝ (Fin n) → W) (s ε : ℝ) (i : Fin n)
    (x : EuclideanSpace ℝ (Fin n)) : ℝ :=
  inner ℝ (η x) (pd η i x) * Phi η s ε x

lemma Gi_contDiff (hη : ContDiff ℝ 2 η) {ε : ℝ} (hε : 0 < ε) (s : ℝ) (i : Fin n) :
    ContDiff ℝ 1 (Gi η s ε i) :=
  (ContDiff.inner ℝ (hη.of_le (by norm_num)) (pd_contDiff hη i)).mul (Phi_contDiff hη hε s)

lemma Gi_compact (hc : HasCompactSupport η) (s ε : ℝ) (i : Fin n) :
    HasCompactSupport (Gi η s ε i) := by
  refine HasCompactSupport.intro hc (fun x hx => ?_)
  simp [Gi, image_eq_zero_of_notMem_tsupport hx]

lemma fderiv_Gi (hη : ContDiff ℝ 2 η) {ε : ℝ} (hε : 0 < ε) (s : ℝ) (i : Fin n)
    (x : EuclideanSpace ℝ (Fin n)) :
    fderiv ℝ (Gi η s ε i) x (EuclideanSpace.single i 1) =
      (‖pd η i x‖ ^ 2 + inner ℝ (η x) (pd2 η i i x)) * Phi η s ε x +
        inner ℝ (η x) (pd η i x) * ((s - 1) * (gradSq η x + ε) ^ (s - 1 - 1) *
          ∑ k, 2 * inner ℝ (pd η k x) (pd2 η i k x)) := by
  have hη' : DifferentiableAt ℝ η x := (hη.differentiable (by norm_num)) x
  have hp : DifferentiableAt ℝ (pd η i) x := ((pd_contDiff hη i).differentiable (by norm_num)) x
  have hI : DifferentiableAt ℝ (fun y => inner ℝ (η y) (pd η i y)) x := hη'.inner ℝ hp
  have hP : DifferentiableAt ℝ (Phi η s ε) x :=
    ((Phi_contDiff hη hε s).differentiable (by norm_num)) x
  have h2 : Gi η s ε i = fun y => inner ℝ (η y) (pd η i y) * Phi η s ε y := rfl
  rw [h2, fderiv_fun_mul hI hP]
  simp only [ContinuousLinearMap.add_apply, ContinuousLinearMap.smul_apply, smul_eq_mul]
  rw [fderiv_inner_apply ℝ hη' hp, fderiv_Phi hη hε s x i]
  have h3 : fderiv ℝ η x (EuclideanSpace.single i 1) = pd η i x := rfl
  have h4 : fderiv ℝ (pd η i) x (EuclideanSpace.single i 1) = pd2 η i i x := rfl
  rw [h3, h4, real_inner_self_eq_norm_sq]
  ring

lemma pointwise_ineq (hη : ContDiff ℝ 2 η) {s : ℝ} (hs : 1 ≤ s) {ε : ℝ} (hε : 0 < ε)
    (x : EuclideanSpace ℝ (Fin n)) :
    gradSq η x * Phi η s ε x ≤
      (Real.sqrt n + 2 * (s - 1)) * (‖η x‖ * Real.sqrt (hessSq η x) * Phi η s ε x) +
        ∑ i, fderiv ℝ (Gi η s ε i) x (EuclideanSpace.single i 1) := by
  have hg0 := gradSq_nonneg η x
  have hpos := gradSq_add_pos η hε x
  have hB : 0 < (gradSq η x + ε) ^ (s - 1 - 1) := Real.rpow_pos_of_pos hpos _
  have hΦ : Phi η s ε x = (gradSq η x + ε) ^ (s - 1 - 1) * (gradSq η x + ε) := by
    rw [← Real.rpow_add_one hpos.ne']
    unfold Phi
    congr 1
    ring
  have hΦ0 := Phi_nonneg η hε s x
  -- abbreviations
  obtain ⟨B, hBdef⟩ : ∃ B, B = (gradSq η x + ε) ^ (s - 1 - 1) := ⟨_, rfl⟩
  rw [← hBdef] at hB hΦ
  obtain ⟨T, hT⟩ : ∃ T, T = ∑ i, inner ℝ (η x) (pd η i x) *
      ∑ k, inner ℝ (pd η k x) (pd2 η i k x) := ⟨_, rfl⟩
  obtain ⟨Sx, hSx⟩ : ∃ Sx, Sx = ∑ i, pd2 η i i x := ⟨_, rfl⟩
  have hsum : ∑ i, fderiv ℝ (Gi η s ε i) x (EuclideanSpace.single i 1) =
      (gradSq η x + inner ℝ (η x) Sx) * Phi η s ε x + (s - 1) * B * 2 * T := by
    simp only [fderiv_Gi hη hε s, ← hBdef]
    rw [Finset.sum_add_distrib, ← Finset.sum_mul, Finset.sum_add_distrib, hSx, inner_sum, hT,
      Finset.mul_sum, gradSq_eq]
    congr 1
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [← Finset.mul_sum]
    ring
  have hSb : |inner ℝ (η x) Sx| ≤ ‖η x‖ * (Real.sqrt n * Real.sqrt (hessSq η x)) := by
    refine (abs_real_inner_le_norm _ _).trans (mul_le_mul_of_nonneg_left ?_ (norm_nonneg _))
    rw [hSx]
    refine (norm_sum_le_sqrt_mul _).trans ?_
    rw [Fintype.card_fin]
    refine mul_le_mul_of_nonneg_left (Real.sqrt_le_sqrt ?_) (Real.sqrt_nonneg _)
    rw [hessSq_eq]
    exact Finset.sum_le_sum fun i _ =>
      Finset.single_le_sum (f := fun k => ‖pd2 η i k x‖ ^ 2) (fun _ _ => by positivity)
        (Finset.mem_univ i)
  have hTb : |T| ≤ ‖η x‖ * gradSq η x * Real.sqrt (hessSq η x) := by
    rw [hT]
    exact abs_double_sum_le (η x) (fun i => pd η i x) (fun i k => pd2 η i k x)
  have hR : 0 ≤ ‖η x‖ * Real.sqrt (hessSq η x) := by positivity
  have hs1 : 0 ≤ s - 1 := by linarith
  have k1 : -(‖η x‖ * (Real.sqrt n * Real.sqrt (hessSq η x))) * Phi η s ε x ≤
      inner ℝ (η x) Sx * Phi η s ε x :=
    mul_le_mul_of_nonneg_right (neg_le_of_abs_le hSb) hΦ0
  have k2 : (2 * (s - 1) * B) * (-(‖η x‖ * gradSq η x * Real.sqrt (hessSq η x))) ≤
      (2 * (s - 1) * B) * T :=
    mul_le_mul_of_nonneg_left (neg_le_of_abs_le hTb) (by positivity)
  have hBg : B * gradSq η x ≤ Phi η s ε x := by
    rw [hΦ]; nlinarith
  have k3 : (2 * (s - 1) * (‖η x‖ * Real.sqrt (hessSq η x))) * (B * gradSq η x) ≤
      (2 * (s - 1) * (‖η x‖ * Real.sqrt (hessSq η x))) * Phi η s ε x :=
    mul_le_mul_of_nonneg_left hBg (by positivity)
  rw [hsum]
  nlinarith [k1, k2, k3]

lemma eps_ineq (hη : ContDiff ℝ 2 η) (hc : HasCompactSupport η) {s : ℝ} (hs : 1 ≤ s) {ε : ℝ}
    (hε : 0 < ε) :
    ∫ x, gradSq η x * Phi η s ε x ≤
      (Real.sqrt n + 2 * (s - 1)) *
        ∫ x, ‖η x‖ * Real.sqrt (hessSq η x) * Phi η s ε x := by
  have hint_i : ∀ i : Fin n, Integrable
      (fun x => fderiv ℝ (Gi η s ε i) x (EuclideanSpace.single i 1)) := fun i =>
    integrable_fderiv_apply (Gi_contDiff hη hε s i) (Gi_compact hc s ε i) _
  have h0 : ∫ x, ∑ i, fderiv ℝ (Gi η s ε i) x (EuclideanSpace.single i 1) = 0 := by
    rw [integral_finset_sum _ (fun i _ => hint_i i)]
    exact Finset.sum_eq_zero fun i _ =>
      integral_fderiv_apply_eq_zero (Gi_contDiff hη hε s i) (Gi_compact hc s ε i) _
  have hΨ : Integrable (fun x => ∑ i, fderiv ℝ (Gi η s ε i) x (EuclideanSpace.single i 1)) :=
    integrable_finset_sum _ (fun i _ => hint_i i)
  have hPhi := (Phi_contDiff hη hε s).continuous
  have hL : Integrable (fun x => gradSq η x * Phi η s ε x) :=
    (((gradSq_contDiff hη).continuous).mul hPhi).integrable_of_hasCompactSupport
      (HasCompactSupport.intro hc fun x hx => by simp [gradSq_eq_zero hc hx])
  have hR : Integrable (fun x => ‖η x‖ * Real.sqrt (hessSq η x) * Phi η s ε x) :=
    (((hη.continuous.norm).mul (hessSq_continuous hη).sqrt).mul hPhi).integrable_of_hasCompactSupport
      (HasCompactSupport.intro hc fun x hx => by simp [image_eq_zero_of_notMem_tsupport hx])
  calc ∫ x, gradSq η x * Phi η s ε x
      ≤ ∫ x, ((Real.sqrt n + 2 * (s - 1)) *
          (‖η x‖ * Real.sqrt (hessSq η x) * Phi η s ε x) +
          ∑ i, fderiv ℝ (Gi η s ε i) x (EuclideanSpace.single i 1)) :=
        integral_mono hL ((hR.const_mul _).add hΨ) (fun x => pointwise_ineq hη hs hε x)
    _ = _ := by
        rw [integral_add (hR.const_mul _) hΨ, integral_const_mul, h0, add_zero]

lemma tendsto_Phi (η : EuclideanSpace ℝ (Fin n) → W) {s : ℝ} (hs : 1 ≤ s)
    (x : EuclideanSpace ℝ (Fin n)) :
    Tendsto (fun m : ℕ => Phi η s (1 / ((m : ℝ) + 1)) x) atTop (𝓝 (gradSq η x ^ (s - 1))) := by
  have h1 : Tendsto (fun m : ℕ => gradSq η x + 1 / ((m : ℝ) + 1)) atTop (𝓝 (gradSq η x)) := by
    simpa using (tendsto_one_div_add_atTop_nhds_zero_nat (𝕜 := ℝ)).const_add (gradSq η x)
  exact ((Real.continuousAt_rpow_const _ _ (Or.inr (by linarith))).tendsto).comp h1

lemma mul_rpow_sub_one {g s : ℝ} (hg : 0 ≤ g) (hs : 1 ≤ s) : g * g ^ (s - 1) = g ^ s := by
  rcases hg.eq_or_lt with h | h
  · subst h
    rw [Real.zero_rpow (by linarith : s ≠ 0)]
    simp
  · conv_rhs => rw [show s = 1 + (s - 1) by ring, Real.rpow_add h, Real.rpow_one]

theorem integral_gradSq_rpow_le {η : EuclideanSpace ℝ (Fin n) → W} (hη : ContDiff ℝ 2 η)
    (hc : HasCompactSupport η) {s : ℝ} (hs : 1 ≤ s) :
    ∫ x, gradSq η x ^ s ≤
      (n + 2 * s - 2) * ∫ x, ‖η x‖ * Real.sqrt (hessSq η x) * gradSq η x ^ (s - 1) := by
  have hε : ∀ m : ℕ, (0 : ℝ) < 1 / ((m : ℝ) + 1) := fun m => by positivity
  have hε1 : ∀ m : ℕ, 1 / ((m : ℝ) + 1) ≤ 1 := fun m => by
    rw [div_le_one (by positivity)]
    have := m.cast_nonneg (α := ℝ)
    linarith
  have hcont_g := (gradSq_contDiff hη).continuous
  have hbw : Continuous (fun x => (gradSq η x + 1) ^ (s - 1)) :=
    (hcont_g.add continuous_const).rpow_const (fun x => Or.inr (by linarith))
  have hbL : Integrable (fun x => gradSq η x * (gradSq η x + 1) ^ (s - 1)) :=
    (hcont_g.mul hbw).integrable_of_hasCompactSupport
      (HasCompactSupport.intro hc fun x hx => by simp [gradSq_eq_zero hc hx])
  have hbR : Integrable (fun x => ‖η x‖ * Real.sqrt (hessSq η x) * (gradSq η x + 1) ^ (s - 1)) :=
    (((hη.continuous.norm).mul (hessSq_continuous hη).sqrt).mul hbw).integrable_of_hasCompactSupport
      (HasCompactSupport.intro hc fun x hx => by simp [image_eq_zero_of_notMem_tsupport hx])
  have hL : Tendsto (fun m : ℕ => ∫ x, gradSq η x * Phi η s (1 / ((m : ℝ) + 1)) x) atTop
      (𝓝 (∫ x, gradSq η x ^ s)) := by
    refine tendsto_integral_of_dominated_convergence
      (fun x => gradSq η x * (gradSq η x + 1) ^ (s - 1)) (fun m => ?_) hbL (fun m => ?_) ?_
    · exact (hcont_g.mul (Phi_contDiff hη (hε m) s).continuous).aestronglyMeasurable
    · refine Filter.Eventually.of_forall fun x => ?_
      have h0 := gradSq_nonneg η x
      rw [Real.norm_of_nonneg (mul_nonneg h0 (Phi_nonneg η (hε m) s x))]
      exact mul_le_mul_of_nonneg_left (Real.rpow_le_rpow (gradSq_add_pos η (hε m) x).le
        (by linarith [hε1 m]) (by linarith)) h0
    · refine Filter.Eventually.of_forall fun x => ?_
      have h := (tendsto_Phi η hs x).const_mul (gradSq η x)
      rw [mul_rpow_sub_one (gradSq_nonneg η x) hs] at h
      exact h
  have hR : Tendsto (fun m : ℕ => ∫ x, ‖η x‖ * Real.sqrt (hessSq η x) *
        Phi η s (1 / ((m : ℝ) + 1)) x) atTop
      (𝓝 (∫ x, ‖η x‖ * Real.sqrt (hessSq η x) * gradSq η x ^ (s - 1))) := by
    refine tendsto_integral_of_dominated_convergence
      (fun x => ‖η x‖ * Real.sqrt (hessSq η x) * (gradSq η x + 1) ^ (s - 1))
      (fun m => ?_) hbR (fun m => ?_) ?_
    · exact (((hη.continuous.norm).mul (hessSq_continuous hη).sqrt).mul
        (Phi_contDiff hη (hε m) s).continuous).aestronglyMeasurable
    · refine Filter.Eventually.of_forall fun x => ?_
      have h0 : 0 ≤ ‖η x‖ * Real.sqrt (hessSq η x) := by positivity
      rw [Real.norm_of_nonneg (mul_nonneg h0 (Phi_nonneg η (hε m) s x))]
      exact mul_le_mul_of_nonneg_left (Real.rpow_le_rpow (gradSq_add_pos η (hε m) x).le
        (by linarith [hε1 m]) (by linarith)) h0
    · exact Filter.Eventually.of_forall fun x => (tendsto_Phi η hs x).const_mul _
  have hle := le_of_tendsto_of_tendsto' hL (hR.const_mul (Real.sqrt n + 2 * (s - 1)))
    (fun m => eps_ineq hη hc hs (hε m))
  have hRnn : 0 ≤ ∫ x, ‖η x‖ * Real.sqrt (hessSq η x) * gradSq η x ^ (s - 1) :=
    integral_nonneg fun x => mul_nonneg (mul_nonneg (norm_nonneg _) (Real.sqrt_nonneg _))
      (Real.rpow_nonneg (gradSq_nonneg η x) _)
  have hsq : Real.sqrt n ≤ n := by
    rw [Real.sqrt_le_iff]
    refine ⟨Nat.cast_nonneg n, ?_⟩
    rcases Nat.eq_zero_or_pos n with h | h
    · subst h; simp
    · have : (1 : ℝ) ≤ n := Nat.one_le_cast.mpr h
      nlinarith
  refine hle.trans (mul_le_mul_of_nonneg_right ?_ hRnn)
  linarith

end Reg2

/-! ### Hölder forms -/

section Holder

lemma integral_le_of_lintegral_le {α : Type*} [MeasurableSpace α] {μ : Measure α} {f : α → ℝ}
    (hf0 : ∀ x, 0 ≤ f x) {r : ℝ} (hr : 0 ≤ r)
    (h : ∫⁻ x, ENNReal.ofReal (f x) ∂μ ≤ ENNReal.ofReal r) : ∫ x, f x ∂μ ≤ r := by
  by_cases hf : Integrable f μ
  · rw [← ENNReal.ofReal_le_ofReal_iff hr,
      ofReal_integral_eq_lintegral_ofReal hf (Filter.Eventually.of_forall hf0)]
    exact h
  · rw [integral_undef hf]
    exact hr

/-- Hölder's inequality for three nonnegative functions with weights summing to one. -/
lemma integral_holder3 {α : Type*} [MeasurableSpace α] {μ : Measure α} {A B C : α → ℝ}
    (hA0 : ∀ x, 0 ≤ A x) (hB0 : ∀ x, 0 ≤ B x) (hC0 : ∀ x, 0 ≤ C x)
    (hA : Integrable A μ) (hB : Integrable B μ) (hC : Integrable C μ)
    {a b c : ℝ} (ha : 0 ≤ a) (hb : 0 ≤ b) (hc : 0 ≤ c) (habc : a + b + c = 1) :
    ∫ x, A x ^ a * B x ^ b * C x ^ c ∂μ ≤
      (∫ x, A x ∂μ) ^ a * (∫ x, B x ∂μ) ^ b * (∫ x, C x ∂μ) ^ c := by
  have hIA : 0 ≤ ∫ x, A x ∂μ := integral_nonneg hA0
  have hIB : 0 ≤ ∫ x, B x ∂μ := integral_nonneg hB0
  have hIC : 0 ≤ ∫ x, C x ∂μ := integral_nonneg hC0
  have key := ENNReal.lintegral_prod_norm_pow_le (μ := μ) (Finset.univ : Finset (Fin 3))
    (f := ![fun x => ENNReal.ofReal (A x), fun x => ENNReal.ofReal (B x),
      fun x => ENNReal.ofReal (C x)])
    (fun i _ => by
      fin_cases i
      · exact hA.aemeasurable.ennreal_ofReal
      · exact hB.aemeasurable.ennreal_ofReal
      · exact hC.aemeasurable.ennreal_ofReal)
    (p := ![a, b, c]) (by simp [Fin.sum_univ_three, habc])
    (fun i _ => by fin_cases i <;> simpa)
  simp only [Fin.prod_univ_three, Matrix.cons_val_zero, Matrix.cons_val_one,
    Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons] at key
  have hX := Real.rpow_nonneg hIA a
  have hY := Real.rpow_nonneg hIB b
  have hZ := Real.rpow_nonneg hIC c
  refine integral_le_of_lintegral_le
    (fun x => mul_nonneg (mul_nonneg (Real.rpow_nonneg (hA0 x) a) (Real.rpow_nonneg (hB0 x) b))
      (Real.rpow_nonneg (hC0 x) c))
    (mul_nonneg (mul_nonneg hX hY) hZ) ?_
  have e1 : ∀ x, ENNReal.ofReal (A x ^ a * B x ^ b * C x ^ c) =
      ENNReal.ofReal (A x) ^ a * ENNReal.ofReal (B x) ^ b * ENNReal.ofReal (C x) ^ c := by
    intro x
    rw [ENNReal.ofReal_mul (by have := hA0 x; have := hB0 x; positivity),
      ENNReal.ofReal_mul (by have := hA0 x; positivity),
      ENNReal.ofReal_rpow_of_nonneg (hA0 x) ha, ENNReal.ofReal_rpow_of_nonneg (hB0 x) hb,
      ENNReal.ofReal_rpow_of_nonneg (hC0 x) hc]
  simp_rw [e1]
  refine key.trans (le_of_eq ?_)
  rw [← ofReal_integral_eq_lintegral_ofReal hA (Filter.Eventually.of_forall hA0),
    ← ofReal_integral_eq_lintegral_ofReal hB (Filter.Eventually.of_forall hB0),
    ← ofReal_integral_eq_lintegral_ofReal hC (Filter.Eventually.of_forall hC0),
    ENNReal.ofReal_rpow_of_nonneg hIA ha, ENNReal.ofReal_rpow_of_nonneg hIB hb,
    ENNReal.ofReal_rpow_of_nonneg hIC hc,
    ← ENNReal.ofReal_mul hX, ← ENNReal.ofReal_mul (mul_nonneg hX hY)]

variable {η : EuclideanSpace ℝ (Fin n) → W}

lemma integrable_rpow_of_compact {f : EuclideanSpace ℝ (Fin n) → ℝ} (hf : Continuous f)
    (hc : HasCompactSupport f) {p : ℝ} (hp : 0 < p) :
    Integrable (fun x => f x ^ p) :=
  (hf.rpow_const (fun x => Or.inr hp.le)).integrable_of_hasCompactSupport
    (hc.comp_left (g := fun y : ℝ => y ^ p) (Real.zero_rpow hp.ne'))

lemma hasCompactSupport_norm (hc : HasCompactSupport η) :
    HasCompactSupport (fun x => ‖η x‖) :=
  HasCompactSupport.intro hc fun x hx => by simp [image_eq_zero_of_notMem_tsupport hx]

lemma hasCompactSupport_sqrt_hessSq (hc : HasCompactSupport η) :
    HasCompactSupport (fun x => Real.sqrt (hessSq η x)) :=
  HasCompactSupport.intro hc fun x hx => by simp [hessSq_eq_zero hx]

lemma hasCompactSupport_gradSq (hc : HasCompactSupport η) :
    HasCompactSupport (gradSq η) :=
  HasCompactSupport.intro hc fun x hx => gradSq_eq_zero hc hx

lemma integrable_main (hη : ContDiff ℝ 2 η) (hc : HasCompactSupport η) {s : ℝ} (hs : 1 ≤ s) :
    Integrable (fun x => ‖η x‖ * Real.sqrt (hessSq η x) * gradSq η x ^ (s - 1)) :=
  ((hη.continuous.norm.mul (hessSq_continuous hη).sqrt).mul
    ((gradSq_contDiff hη).continuous.rpow_const (fun x => Or.inr (by linarith)))
    ).integrable_of_hasCompactSupport
    (HasCompactSupport.intro hc fun x hx => by simp [image_eq_zero_of_notMem_tsupport hx])

lemma integrable_sqrt_mul (hη : ContDiff ℝ 2 η) (hc : HasCompactSupport η) {s : ℝ} (hs : 1 ≤ s) :
    Integrable (fun x => Real.sqrt (hessSq η x) * gradSq η x ^ (s - 1)) :=
  (((hessSq_continuous hη).sqrt).mul
    ((gradSq_contDiff hη).continuous.rpow_const (fun x => Or.inr (by linarith)))
    ).integrable_of_hasCompactSupport
    (HasCompactSupport.intro hc fun x hx => by simp [hessSq_eq_zero hx])

lemma rpow_cancel {y q : ℝ} (hy : 0 ≤ y) (hq : 0 < q) : (y ^ q) ^ (1 / q) = y := by
  rw [one_div, Real.rpow_rpow_inv hy hq.ne']

lemma const_nonneg {s : ℝ} (hs : 1 ≤ s) : (0 : ℝ) ≤ (n : ℝ) + 2 * s - 2 := by
  have := (Nat.cast_nonneg n : (0 : ℝ) ≤ n)
  linarith

/-- Final absorption step: `X ≤ C Y`, `Y ≤ Q X^((s-1)/s)` give `X^(1/s) ≤ C Q`. -/
lemma absorb {X C Q s : ℝ} (hs : 1 ≤ s) (hX : 0 ≤ X) (hC : 0 ≤ C) (hQ : 0 ≤ Q)
    (h : X ≤ C * (Q * X ^ ((s - 1) / s))) : X ^ (1 / s) ≤ C * Q := by
  have hs0 : 0 < s := by linarith
  rcases hX.eq_or_lt with h0 | hpos
  · rw [← h0, Real.zero_rpow (by positivity)]
    positivity
  · have hpow : 0 < X ^ ((s - 1) / s) := Real.rpow_pos_of_pos hpos _
    have hsplit : X = X ^ (1 / s) * X ^ ((s - 1) / s) := by
      rw [← Real.rpow_add hpos]
      have : 1 / s + (s - 1) / s = 1 := by field_simp; ring
      rw [this, Real.rpow_one]
    have : X ^ (1 / s) * X ^ ((s - 1) / s) ≤ (C * Q) * X ^ ((s - 1) / s) := by
      calc X ^ (1 / s) * X ^ ((s - 1) / s) = X := hsplit.symm
        _ ≤ _ := h
        _ = _ := by ring
    exact le_of_mul_le_mul_right this hpow

theorem hamilton_step {η : EuclideanSpace ℝ (Fin n) → W} (hη : ContDiff ℝ 2 η)
    (hc : HasCompactSupport η) {p q s : ℝ} (hs : 1 ≤ s) (hp : 0 < p) (hq : 0 < q)
    (hpq : 1 / p + 1 / q = 1 / s) :
    (∫ x, gradSq η x ^ s) ^ (1 / s) ≤
      (n + 2 * s - 2) * (∫ x, Real.sqrt (hessSq η x) ^ p) ^ (1 / p) *
        (∫ x, ‖η x‖ ^ q) ^ (1 / q) := by
  have hs0 : 0 < s := by linarith
  have h1 := integral_gradSq_rpow_le hη hc hs
  have hH := integral_holder3 (μ := volume)
    (A := fun x => ‖η x‖ ^ q) (B := fun x => Real.sqrt (hessSq η x) ^ p)
    (C := fun x => gradSq η x ^ s)
    (fun x => Real.rpow_nonneg (norm_nonneg _) _) (fun x => Real.rpow_nonneg (Real.sqrt_nonneg _) _)
    (fun x => Real.rpow_nonneg (gradSq_nonneg η x) _)
    (integrable_rpow_of_compact hη.continuous.norm (hasCompactSupport_norm hc) hq)
    (integrable_rpow_of_compact (hessSq_continuous hη).sqrt (hasCompactSupport_sqrt_hessSq hc) hp)
    (integrable_rpow_of_compact (gradSq_contDiff hη).continuous (hasCompactSupport_gradSq hc) hs0)
    (a := 1 / q) (b := 1 / p) (c := (s - 1) / s) (by positivity) (by positivity)
    (div_nonneg (by linarith) hs0.le) (by
      have : (s - 1) / s = 1 - 1 / s := by field_simp
      linarith)
  have hfun : (fun x => (‖η x‖ ^ q) ^ (1 / q) * (Real.sqrt (hessSq η x) ^ p) ^ (1 / p) *
      (gradSq η x ^ s) ^ ((s - 1) / s)) =
      fun x => ‖η x‖ * Real.sqrt (hessSq η x) * gradSq η x ^ (s - 1) := by
    funext x
    rw [rpow_cancel (norm_nonneg _) hq, rpow_cancel (Real.sqrt_nonneg _) hp,
      ← Real.rpow_mul (gradSq_nonneg η x), mul_div_cancel₀ _ hs0.ne']
  rw [hfun] at hH
  rw [mul_assoc]
  refine absorb hs (integral_nonneg fun x => Real.rpow_nonneg (gradSq_nonneg η x) _)
    (const_nonneg hs) (by positivity) ?_
  calc ∫ x, gradSq η x ^ s
      ≤ (n + 2 * s - 2) * ∫ x, ‖η x‖ * Real.sqrt (hessSq η x) * gradSq η x ^ (s - 1) := h1
    _ ≤ (n + 2 * s - 2) * ((∫ x, ‖η x‖ ^ q) ^ (1 / q) * (∫ x, Real.sqrt (hessSq η x) ^ p) ^ (1 / p) *
        (∫ x, gradSq η x ^ s) ^ ((s - 1) / s)) :=
        mul_le_mul_of_nonneg_left hH (const_nonneg hs)
    _ = (n + 2 * s - 2) * (((∫ x, Real.sqrt (hessSq η x) ^ p) ^ (1 / p) *
        (∫ x, ‖η x‖ ^ q) ^ (1 / q)) * (∫ x, gradSq η x ^ s) ^ ((s - 1) / s)) := by ring

theorem hamilton_step_sup {η : EuclideanSpace ℝ (Fin n) → W} (hη : ContDiff ℝ 2 η)
    (hc : HasCompactSupport η) {s M : ℝ} (hs : 1 ≤ s) (hM : ∀ x, ‖η x‖ ≤ M) :
    (∫ x, gradSq η x ^ s) ^ (1 / s) ≤
      (n + 2 * s - 2) * M * (∫ x, Real.sqrt (hessSq η x) ^ s) ^ (1 / s) := by
  have hs0 : 0 < s := by linarith
  have hM0 : 0 ≤ M := (norm_nonneg _).trans (hM 0)
  have h1 := integral_gradSq_rpow_le hη hc hs
  have h2 : ∫ x, ‖η x‖ * Real.sqrt (hessSq η x) * gradSq η x ^ (s - 1) ≤
      M * ∫ x, Real.sqrt (hessSq η x) * gradSq η x ^ (s - 1) := by
    rw [← integral_const_mul]
    refine integral_mono (integrable_main hη hc hs) ((integrable_sqrt_mul hη hc hs).const_mul M)
      fun x => ?_
    have h0 : 0 ≤ Real.sqrt (hessSq η x) * gradSq η x ^ (s - 1) :=
      mul_nonneg (Real.sqrt_nonneg _) (Real.rpow_nonneg (gradSq_nonneg η x) _)
    calc ‖η x‖ * Real.sqrt (hessSq η x) * gradSq η x ^ (s - 1)
        = ‖η x‖ * (Real.sqrt (hessSq η x) * gradSq η x ^ (s - 1)) := by ring
      _ ≤ M * (Real.sqrt (hessSq η x) * gradSq η x ^ (s - 1)) :=
        mul_le_mul_of_nonneg_right (hM x) h0
  have hH := integral_holder3 (μ := volume)
    (A := fun x => Real.sqrt (hessSq η x) ^ s) (B := fun x => gradSq η x ^ s)
    (C := fun x => Real.sqrt (hessSq η x) ^ s)
    (fun x => Real.rpow_nonneg (Real.sqrt_nonneg _) _)
    (fun x => Real.rpow_nonneg (gradSq_nonneg η x) _)
    (fun x => Real.rpow_nonneg (Real.sqrt_nonneg _) _)
    (integrable_rpow_of_compact (hessSq_continuous hη).sqrt (hasCompactSupport_sqrt_hessSq hc) hs0)
    (integrable_rpow_of_compact (gradSq_contDiff hη).continuous (hasCompactSupport_gradSq hc) hs0)
    (integrable_rpow_of_compact (hessSq_continuous hη).sqrt (hasCompactSupport_sqrt_hessSq hc) hs0)
    (a := 1 / s) (b := (s - 1) / s) (c := 0) (by positivity)
    (div_nonneg (by linarith) hs0.le) le_rfl (by
      have : (s - 1) / s = 1 - 1 / s := by field_simp
      linarith)
  have hfun : (fun x => (Real.sqrt (hessSq η x) ^ s) ^ (1 / s) * (gradSq η x ^ s) ^ ((s - 1) / s) *
      (Real.sqrt (hessSq η x) ^ s) ^ (0 : ℝ)) =
      fun x => Real.sqrt (hessSq η x) * gradSq η x ^ (s - 1) := by
    funext x
    rw [rpow_cancel (Real.sqrt_nonneg _) hs0, Real.rpow_zero, mul_one,
      ← Real.rpow_mul (gradSq_nonneg η x), mul_div_cancel₀ _ hs0.ne']
  rw [hfun, Real.rpow_zero, mul_one] at hH
  refine absorb hs (integral_nonneg fun x => Real.rpow_nonneg (gradSq_nonneg η x) _)
    (mul_nonneg (const_nonneg hs) hM0) (Real.rpow_nonneg (integral_nonneg fun x =>
      Real.rpow_nonneg (Real.sqrt_nonneg _) _) _) ?_
  calc ∫ x, gradSq η x ^ s
      ≤ (n + 2 * s - 2) * ∫ x, ‖η x‖ * Real.sqrt (hessSq η x) * gradSq η x ^ (s - 1) := h1
    _ ≤ (n + 2 * s - 2) * (M * ∫ x, Real.sqrt (hessSq η x) * gradSq η x ^ (s - 1)) :=
        mul_le_mul_of_nonneg_left h2 (const_nonneg hs)
    _ ≤ (n + 2 * s - 2) * (M * ((∫ x, Real.sqrt (hessSq η x) ^ s) ^ (1 / s) *
        (∫ x, gradSq η x ^ s) ^ ((s - 1) / s))) :=
        mul_le_mul_of_nonneg_left (mul_le_mul_of_nonneg_left hH hM0) (const_nonneg hs)
    _ = (n + 2 * s - 2) * M * ((∫ x, Real.sqrt (hessSq η x) ^ s) ^ (1 / s) *
        (∫ x, gradSq η x ^ s) ^ ((s - 1) / s)) := by ring

end Holder

/-- info: 'CWFlow.Interp.integral_gradSq_rpow_le' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms integral_gradSq_rpow_le

/-- info: 'CWFlow.Interp.hamilton_step' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms hamilton_step

/-- info: 'CWFlow.Interp.hamilton_step_sup' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms hamilton_step_sup

end CWFlow.Interp
