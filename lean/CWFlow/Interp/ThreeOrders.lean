import CWFlow.Interp.HamiltonStep
import CWFlow.Interp.Jets
import Mathlib.Analysis.Calculus.BumpFunction.Basic
import Mathlib.Analysis.Calculus.BumpFunction.FiniteDimension

/-!
# Interpolation between three consecutive orders on balls of `ℝⁿ`

For `ρ < R` and `ε > 0` there is a constant `C` such that for every `C²` map `η` from `ℝⁿ` to a
finite-dimensional inner product space,
`∫_{B_ρ} |Dη|² ≤ ε ∫_{B_R} |D²η|² + C ∫_{B_R} |η|²`
(`exists_const_gradSq_ball_le`), and, with `C` independent of the order `j`, the same for the jets of a
smooth function (`exists_const_jet_three_orders`).

The proof is one integration by parts against a fixed cutoff `χ`: with `Gᵢ = χ² ⟪η, ∂ᵢη⟫`,
`∫ ∑ᵢ ∂ᵢ Gᵢ = 0`, followed by Cauchy–Schwarz and `ab ≤ δ a² + b² / (4δ)`.
-/

open MeasureTheory Metric
open scoped ContDiff

namespace CWFlow.Interp

universe u v

section Real

lemma real_ineq {c N q r T P D L δ sn : ℝ} (hc0 : 0 ≤ c) (hc1 : c ≤ 1) (hN : 0 ≤ N) (hq : 0 ≤ q)
    (hr : 0 ≤ r) (hL : 0 ≤ L) (hδ : 0 < δ) (hsn : 0 ≤ sn)
    (hT : |T| ≤ L * N * q) (hP : |P| ≤ N * sn * r) (hD : D = 2 * c * T + c ^ 2 * q ^ 2 + c ^ 2 * P) :
    (1 / 2) * (c ^ 2 * q ^ 2) ≤ D + δ * r ^ 2 + (sn ^ 2 / (4 * δ) + 2 * L ^ 2) * N ^ 2 := by
  have h1 : 2 * c * (-T) ≤ 2 * c * (L * N * q) := by
    have : -T ≤ L * N * q := by linarith [neg_le_of_abs_le hT]
    exact mul_le_mul_of_nonneg_left this (by positivity)
  have h2 : 2 * c * (L * N * q) ≤ (1 / 2) * (c ^ 2 * q ^ 2) + 2 * (L * N) ^ 2 := by
    nlinarith [sq_nonneg (c * q - 2 * (L * N))]
  have h3 : c ^ 2 * (-P) ≤ c ^ 2 * (N * sn * r) := by
    have : -P ≤ N * sn * r := by linarith [neg_le_of_abs_le hP]
    exact mul_le_mul_of_nonneg_left this (sq_nonneg c)
  have h4 : N * sn * r ≤ δ * r ^ 2 + sn ^ 2 * N ^ 2 / (4 * δ) := by
    rw [← sub_nonneg]
    have : δ * r ^ 2 + sn ^ 2 * N ^ 2 / (4 * δ) - N * sn * r = (2 * δ * r - sn * N) ^ 2 / (4 * δ) := by
      field_simp
      ring
    rw [this]
    positivity
  have h5 : c ^ 2 * (δ * r ^ 2 + sn ^ 2 * N ^ 2 / (4 * δ)) ≤ δ * r ^ 2 + sn ^ 2 * N ^ 2 / (4 * δ) := by
    have h0 : 0 ≤ δ * r ^ 2 + sn ^ 2 * N ^ 2 / (4 * δ) := by positivity
    have hc2 : c ^ 2 ≤ 1 := by nlinarith
    exact mul_le_of_le_one_left h0 hc2
  have h6 : c ^ 2 * (N * sn * r) ≤ c ^ 2 * (δ * r ^ 2 + sn ^ 2 * N ^ 2 / (4 * δ)) :=
    mul_le_mul_of_nonneg_left h4 (sq_nonneg c)
  have e : (sn ^ 2 / (4 * δ) + 2 * L ^ 2) * N ^ 2 = sn ^ 2 * N ^ 2 / (4 * δ) + 2 * (L * N) ^ 2 := by
    ring
  rw [e]
  nlinarith [h1, h2, h3, h5, h6]

lemma abs_sum_mul_inner_le {W : Type*} [NormedAddCommGroup W] [InnerProductSpace ℝ W] {ι : Type*}
    [Fintype ι] (a : ι → ℝ) (η : W) (D : ι → W) {L : ℝ} (hL : ∑ i, a i ^ 2 ≤ L ^ 2) (hL0 : 0 ≤ L) :
    |∑ i, a i * inner ℝ η (D i)| ≤ L * ‖η‖ * Real.sqrt (∑ i, ‖D i‖ ^ 2) := by
  have hs : Real.sqrt (∑ i, |a i| ^ 2) ≤ L :=
    Real.sqrt_le_iff.mpr ⟨hL0, by simpa only [sq_abs] using hL⟩
  calc |∑ i, a i * inner ℝ η (D i)| ≤ ∑ i, |a i * inner ℝ η (D i)| := Finset.abs_sum_le_sum_abs _ _
    _ ≤ ∑ i, |a i| * (‖η‖ * ‖D i‖) := Finset.sum_le_sum fun i _ => by
        rw [abs_mul]
        exact mul_le_mul_of_nonneg_left (abs_real_inner_le_norm _ _) (abs_nonneg _)
    _ = ‖η‖ * ∑ i, |a i| * ‖D i‖ := by
        rw [Finset.mul_sum]
        exact Finset.sum_congr rfl fun i _ => by ring
    _ ≤ ‖η‖ * (Real.sqrt (∑ i, |a i| ^ 2) * Real.sqrt (∑ i, ‖D i‖ ^ 2)) :=
        mul_le_mul_of_nonneg_left (Real.sum_mul_le_sqrt_mul_sqrt _ _ _) (norm_nonneg _)
    _ ≤ ‖η‖ * (L * Real.sqrt (∑ i, ‖D i‖ ^ 2)) :=
        mul_le_mul_of_nonneg_left (mul_le_mul_of_nonneg_right hs (Real.sqrt_nonneg _))
          (norm_nonneg _)
    _ = _ := by ring

end Real

section Key

variable {n : ℕ} {W : Type*} [NormedAddCommGroup W] [InnerProductSpace ℝ W] [FiniteDimensional ℝ W]

lemma integrableOn_ball_of_continuous {f : EuclideanSpace ℝ (Fin n) → ℝ} (hf : Continuous f)
    (R : ℝ) : IntegrableOn f (ball (0 : EuclideanSpace ℝ (Fin n)) R) :=
  (hf.continuousOn.integrableOn_compact (isCompact_closedBall 0 R)).mono_set ball_subset_closedBall

lemma abs_inner_lap_le (η : EuclideanSpace ℝ (Fin n) → W) (x : EuclideanSpace ℝ (Fin n)) (v : W) :
    |inner ℝ v (∑ i, pd2 η i i x)| ≤ ‖v‖ * Real.sqrt n * Real.sqrt (hessSq η x) := by
  refine (abs_real_inner_le_norm _ _).trans ?_
  rw [mul_assoc]
  refine mul_le_mul_of_nonneg_left ?_ (norm_nonneg _)
  refine (norm_sum_le_sqrt_mul _).trans ?_
  rw [Fintype.card_fin]
  refine mul_le_mul_of_nonneg_left (Real.sqrt_le_sqrt ?_) (Real.sqrt_nonneg _)
  rw [hessSq_eq]
  exact Finset.sum_le_sum fun i _ =>
    Finset.single_le_sum (f := fun k => ‖pd2 η i k x‖ ^ 2) (fun _ _ => by positivity)
      (Finset.mem_univ i)

/-- The vector field component `χ² ⟪η, ∂ᵢη⟫` (written `χ (χ ⟪η, ∂ᵢη⟫)`). -/
noncomputable def Gc (χ : EuclideanSpace ℝ (Fin n) → ℝ) (η : EuclideanSpace ℝ (Fin n) → W) (i : Fin n)
    (x : EuclideanSpace ℝ (Fin n)) : ℝ :=
  χ x * (χ x * inner ℝ (η x) (pd η i x))

lemma Gc_contDiff {χ : EuclideanSpace ℝ (Fin n) → ℝ} (hχ : ContDiff ℝ 1 χ)
    {η : EuclideanSpace ℝ (Fin n) → W} (hη : ContDiff ℝ 2 η) (i : Fin n) :
    ContDiff ℝ 1 (Gc χ η i) :=
  hχ.mul (hχ.mul (ContDiff.inner ℝ (hη.of_le (by norm_num)) (pd_contDiff hη i)))

lemma Gc_compact {χ : EuclideanSpace ℝ (Fin n) → ℝ} (hχc : HasCompactSupport χ)
    (η : EuclideanSpace ℝ (Fin n) → W) (i : Fin n) : HasCompactSupport (Gc χ η i) :=
  HasCompactSupport.intro hχc fun x hx => by simp [Gc, image_eq_zero_of_notMem_tsupport hx]

lemma tsupport_Gc_subset (χ : EuclideanSpace ℝ (Fin n) → ℝ) (η : EuclideanSpace ℝ (Fin n) → W)
    (i : Fin n) : tsupport (Gc χ η i) ⊆ tsupport χ := by
  refine closure_mono fun x hx => ?_
  simp only [Function.mem_support, Gc] at hx ⊢
  intro h
  apply hx
  simp [h]

lemma fderiv_Gc {χ : EuclideanSpace ℝ (Fin n) → ℝ} (hχ : ContDiff ℝ 1 χ)
    {η : EuclideanSpace ℝ (Fin n) → W} (hη : ContDiff ℝ 2 η) (i : Fin n)
    (x : EuclideanSpace ℝ (Fin n)) :
    fderiv ℝ (Gc χ η i) x (EuclideanSpace.single i 1) =
      χ x * (χ x * (inner ℝ (η x) (pd2 η i i x) + ‖pd η i x‖ ^ 2) +
        inner ℝ (η x) (pd η i x) * fderiv ℝ χ x (EuclideanSpace.single i 1)) +
      (χ x * inner ℝ (η x) (pd η i x)) * fderiv ℝ χ x (EuclideanSpace.single i 1) := by
  have hη' : DifferentiableAt ℝ η x := (hη.differentiable (by norm_num)) x
  have hp : DifferentiableAt ℝ (pd η i) x := ((pd_contDiff hη i).differentiable (by norm_num)) x
  have hc : DifferentiableAt ℝ χ x := (hχ.differentiable (by norm_num)) x
  have hI : DifferentiableAt ℝ (fun y => inner ℝ (η y) (pd η i y)) x := hη'.inner ℝ hp
  have hJ : DifferentiableAt ℝ (fun y => χ y * inner ℝ (η y) (pd η i y)) x := hc.mul hI
  have h2 : Gc χ η i = fun y => χ y * (χ y * inner ℝ (η y) (pd η i y)) := rfl
  rw [h2, fderiv_fun_mul hc hJ, fderiv_fun_mul hc hI]
  simp only [ContinuousLinearMap.add_apply, ContinuousLinearMap.smul_apply, smul_eq_mul]
  rw [fderiv_inner_apply ℝ hη' hp]
  have h3 : fderiv ℝ η x (EuclideanSpace.single i 1) = pd η i x := rfl
  have h4 : fderiv ℝ (pd η i) x (EuclideanSpace.single i 1) = pd2 η i i x := rfl
  rw [h3, h4, real_inner_self_eq_norm_sq]

lemma sum_fderiv_Gc {χ : EuclideanSpace ℝ (Fin n) → ℝ} (hχ : ContDiff ℝ 1 χ)
    {η : EuclideanSpace ℝ (Fin n) → W} (hη : ContDiff ℝ 2 η) (x : EuclideanSpace ℝ (Fin n)) :
    ∑ i, fderiv ℝ (Gc χ η i) x (EuclideanSpace.single i 1) =
      2 * χ x * (∑ i, fderiv ℝ χ x (EuclideanSpace.single i 1) * inner ℝ (η x) (pd η i x)) +
      χ x ^ 2 * gradSq η x + χ x ^ 2 * inner ℝ (η x) (∑ i, pd2 η i i x) := by
  have e1 : 2 * χ x * (∑ i, fderiv ℝ χ x (EuclideanSpace.single i 1) * inner ℝ (η x) (pd η i x)) =
      ∑ i, 2 * χ x * (fderiv ℝ χ x (EuclideanSpace.single i 1) * inner ℝ (η x) (pd η i x)) :=
    Finset.mul_sum _ _ _
  have e2 : χ x ^ 2 * gradSq η x = ∑ i, χ x ^ 2 * ‖pd η i x‖ ^ 2 := by
    rw [gradSq_eq, Finset.mul_sum]
  have e3 : χ x ^ 2 * inner ℝ (η x) (∑ i, pd2 η i i x) =
      ∑ i, χ x ^ 2 * inner ℝ (η x) (pd2 η i i x) := by
    rw [inner_sum, Finset.mul_sum]
  rw [e1, e2, e3, ← Finset.sum_add_distrib, ← Finset.sum_add_distrib]
  exact Finset.sum_congr rfl fun i _ => by rw [fderiv_Gc hχ hη]; ring

/-- The cutoff-weighted inequality. -/
theorem key_ineq {χ : EuclideanSpace ℝ (Fin n) → ℝ} (hχ : ContDiff ℝ ∞ χ) (hχc : HasCompactSupport χ)
    (h0 : ∀ x, 0 ≤ χ x) (h1 : ∀ x, χ x ≤ 1) {L : ℝ} (hL0 : 0 ≤ L)
    (hL : ∀ x, ∑ i, (fderiv ℝ χ x (EuclideanSpace.single i 1)) ^ 2 ≤ L ^ 2) {ρ R : ℝ}
    (hρR : ρ ≤ R) (hsupp : tsupport χ ⊆ ball (0 : EuclideanSpace ℝ (Fin n)) R)
    (hχρ : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) ρ, χ x = 1) {δ : ℝ} (hδ : 0 < δ)
    {η : EuclideanSpace ℝ (Fin n) → W} (hη : ContDiff ℝ 2 η) :
    ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) ρ, gradSq η x ≤
      2 * δ * (∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) R, hessSq η x) +
      ((n : ℝ) / (2 * δ) + 4 * L ^ 2) *
        (∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) R, ‖η x‖ ^ 2) := by
  have hχ1 : ContDiff ℝ 1 χ := hχ.of_le (by simp)
  have hpt : ∀ x, (1 / 2) * (χ x ^ 2 * gradSq η x) ≤
      (∑ i, fderiv ℝ (Gc χ η i) x (EuclideanSpace.single i 1)) + δ * hessSq η x +
        ((n : ℝ) / (4 * δ) + 2 * L ^ 2) * ‖η x‖ ^ 2 := by
    intro x
    have hq2 : Real.sqrt (gradSq η x) ^ 2 = gradSq η x := Real.sq_sqrt (gradSq_nonneg _ _)
    have hr2 : Real.sqrt (hessSq η x) ^ 2 = hessSq η x := Real.sq_sqrt (hessSq_nonneg _ _)
    have hsn2 : Real.sqrt (n : ℝ) ^ 2 = n := Real.sq_sqrt (Nat.cast_nonneg n)
    have hT : |∑ i, fderiv ℝ χ x (EuclideanSpace.single i 1) * inner ℝ (η x) (pd η i x)| ≤
        L * ‖η x‖ * Real.sqrt (gradSq η x) := by
      rw [gradSq_eq]
      exact abs_sum_mul_inner_le _ (η x) (fun i => pd η i x) (hL x) hL0
    have key := real_ineq (c := χ x) (N := ‖η x‖) (q := Real.sqrt (gradSq η x))
      (r := Real.sqrt (hessSq η x)) (L := L) (δ := δ) (sn := Real.sqrt n) (h0 x) (h1 x)
      (norm_nonneg _) (Real.sqrt_nonneg _) (Real.sqrt_nonneg _) hL0 hδ (Real.sqrt_nonneg _) hT
      (abs_inner_lap_le η x (η x))
      ((sum_fderiv_Gc hχ1 hη x).trans (by rw [hq2]))
    rw [hq2, hr2, hsn2] at key
    exact key
  have hvan : ∀ x, x ∉ ball (0 : EuclideanSpace ℝ (Fin n)) R →
      ∑ i, fderiv ℝ (Gc χ η i) x (EuclideanSpace.single i 1) = 0 := by
    intro x hx
    refine Finset.sum_eq_zero fun i _ => ?_
    have hx' : x ∉ tsupport (Gc χ η i) := fun h => hx (hsupp (tsupport_Gc_subset χ η i h))
    simp [fderiv_of_notMem_tsupport ℝ hx']
  have hdiv0 : ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) R,
      ∑ i, fderiv ℝ (Gc χ η i) x (EuclideanSpace.single i 1) = 0 := by
    rw [setIntegral_eq_integral_of_forall_compl_eq_zero hvan,
      integral_finset_sum _ (fun i _ => integrable_fderiv_apply (Gc_contDiff hχ1 hη i)
        (Gc_compact hχc η i) _)]
    exact Finset.sum_eq_zero fun i _ =>
      integral_fderiv_apply_eq_zero (Gc_contDiff hχ1 hη i) (Gc_compact hχc η i) _
  have hdivC : Continuous (fun x => ∑ i, fderiv ℝ (Gc χ η i) x (EuclideanSpace.single i 1)) :=
    continuous_finsetSum _ fun i _ =>
      ((Gc_contDiff hχ1 hη i).continuous_fderiv (by norm_num)).clm_apply continuous_const
  have cg : Continuous (fun x => χ x ^ 2 * gradSq η x) :=
    (hχ.continuous.pow 2).mul (gradSq_contDiff hη).continuous
  have cH : Continuous (fun x => hessSq η x) := hessSq_continuous hη
  have cN : Continuous (fun x => ‖η x‖ ^ 2) := hη.continuous.norm.pow 2
  set K : ℝ := (n : ℝ) / (4 * δ) + 2 * L ^ 2 with hK
  have iD := integrableOn_ball_of_continuous hdivC R
  have iH := integrableOn_ball_of_continuous (cH.const_mul δ) R
  have iN := integrableOn_ball_of_continuous (cN.const_mul K) R
  have ig := integrableOn_ball_of_continuous cg R
  have iDH : IntegrableOn (fun x => ∑ i, fderiv ℝ (Gc χ η i) x (EuclideanSpace.single i 1) +
      δ * hessSq η x) (ball (0 : EuclideanSpace ℝ (Fin n)) R) := iD.add iH
  have hmono : ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) R, (1 / 2) * (χ x ^ 2 * gradSq η x) ≤
      ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) R,
        ((∑ i, fderiv ℝ (Gc χ η i) x (EuclideanSpace.single i 1) + δ * hessSq η x) +
          K * ‖η x‖ ^ 2) :=
    integral_mono (μ := volume.restrict (ball (0 : EuclideanSpace ℝ (Fin n)) R))
      (ig.const_mul (1 / 2)) (iDH.add iN) hpt
  rw [integral_add iDH iN, integral_add iD iH, integral_const_mul, integral_const_mul,
    integral_const_mul, hdiv0] at hmono
  have hA : ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) ρ, gradSq η x ≤
      ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) R, χ x ^ 2 * gradSq η x := by
    calc ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) ρ, gradSq η x
        = ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) ρ, χ x ^ 2 * gradSq η x :=
          setIntegral_congr_fun measurableSet_ball fun x hx => by
            simp only [hχρ x hx]; ring
      _ ≤ _ := setIntegral_mono_set ig
          (Filter.Eventually.of_forall fun x => mul_nonneg (sq_nonneg _) (gradSq_nonneg _ _))
          (ball_subset_ball hρR).eventuallyLE
  have hK2 : 2 * K = (n : ℝ) / (2 * δ) + 4 * L ^ 2 := by
    rw [hK]
    field_simp
    ring
  rw [← hK2]
  nlinarith [hmono, hA]

end Key

section Cutoff

variable {n : ℕ}

lemma exists_cutoff {ρ R : ℝ} (hρ : 0 < ρ) (hρR : ρ < R) :
    ∃ (χ : EuclideanSpace ℝ (Fin n) → ℝ) (L : ℝ), ContDiff ℝ ∞ χ ∧ HasCompactSupport χ ∧
      (∀ x, 0 ≤ χ x) ∧ (∀ x, χ x ≤ 1) ∧ 0 ≤ L ∧
      (∀ x, ∑ i, (fderiv ℝ χ x (EuclideanSpace.single i 1)) ^ 2 ≤ L ^ 2) ∧
      tsupport χ ⊆ ball (0 : EuclideanSpace ℝ (Fin n)) R ∧
      ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) ρ, χ x = 1 := by
  let f : ContDiffBump (0 : EuclideanSpace ℝ (Fin n)) :=
    ⟨ρ, (ρ + R) / 2, hρ, by linarith⟩
  have hc : HasCompactSupport (⇑f) := f.hasCompactSupport
  have hχ : ContDiff ℝ ∞ (⇑f) := f.contDiff
  have hg : Continuous (fun x => ∑ i, (fderiv ℝ (⇑f) x (EuclideanSpace.single i 1)) ^ 2) :=
    continuous_finsetSum _ fun i _ =>
      ((hχ.continuous_fderiv (by simp)).clm_apply continuous_const).pow 2
  have hgc : HasCompactSupport
      (fun x => ∑ i, (fderiv ℝ (⇑f) x (EuclideanSpace.single i 1)) ^ 2) :=
    HasCompactSupport.intro hc fun x hx => by simp [fderiv_of_notMem_tsupport ℝ hx]
  obtain ⟨C, hC⟩ := hg.bounded_above_of_compact_support hgc
  refine ⟨⇑f, Real.sqrt C, hχ, hc, fun x => f.nonneg, fun x => f.le_one, Real.sqrt_nonneg _,
    fun x => ?_, ?_, fun x hx => f.one_of_mem_closedBall (ball_subset_closedBall hx)⟩
  · have h := hC x
    rw [Real.norm_of_nonneg (Finset.sum_nonneg fun _ _ => sq_nonneg _)] at h
    rw [Real.sq_sqrt (le_trans (norm_nonneg _) (hC x))]
    exact h
  · rw [f.tsupport_eq]
    exact closedBall_subset_ball (by show (ρ + R) / 2 < R; linarith)

end Cutoff

section Main

variable {n : ℕ}

/-- Universe-uniform form: the constant is chosen before the target space. -/
theorem exists_const_gradSq_ball_le_univ {ρ R ε : ℝ} (hρ : 0 < ρ) (hρR : ρ < R) (hε : 0 < ε) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ (W : Type u) [NormedAddCommGroup W] [InnerProductSpace ℝ W]
      [FiniteDimensional ℝ W] (η : EuclideanSpace ℝ (Fin n) → W), ContDiff ℝ 2 η →
      ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) ρ, gradSq η x ≤
        ε * (∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) R, hessSq η x) +
        C * (∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) R, ‖η x‖ ^ 2) := by
  obtain ⟨χ, L, hχ, hχc, h0, h1, hL0, hL, hsupp, hχρ⟩ := exists_cutoff (n := n) hρ hρR
  refine ⟨(n : ℝ) / ε + 4 * L ^ 2, by positivity, fun W _ _ _ η hη => ?_⟩
  have h := key_ineq hχ hχc h0 h1 hL0 hL hρR.le hsupp hχρ (half_pos hε) hη
  have e1 : 2 * (ε / 2) = ε := by ring
  rw [e1] at h
  exact h

theorem exists_const_gradSq_ball_le {W : Type v} [NormedAddCommGroup W] [InnerProductSpace ℝ W]
    [FiniteDimensional ℝ W] {ρ R ε : ℝ} (hρ : 0 < ρ) (hρR : ρ < R) (hε : 0 < ε) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ η : EuclideanSpace ℝ (Fin n) → W, ContDiff ℝ 2 η →
      ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) ρ, gradSq η x ≤
        ε * (∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) R, hessSq η x) +
        C * (∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) R, ‖η x‖ ^ 2) := by
  obtain ⟨C, hC0, hC⟩ := exists_const_gradSq_ball_le_univ.{v} (n := n) hρ hρR hε
  exact ⟨C, hC0, fun η hη => hC W η hη⟩

theorem exists_const_jet_three_orders {V : Type v} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
    [FiniteDimensional ℝ V] {ρ R ε : ℝ} (hρ : 0 < ρ) (hρR : ρ < R) (hε : 0 < ε) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ (j : ℕ) (f : EuclideanSpace ℝ (Fin n) → V), ContDiff ℝ ∞ f →
      ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) ρ, ‖jet (j + 1) f x‖ ^ 2 ≤
        ε * (∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) R, ‖jet (j + 2) f x‖ ^ 2) +
        C * (∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) R, ‖jet j f x‖ ^ 2) := by
  obtain ⟨C, hC0, hC⟩ := exists_const_gradSq_ball_le_univ.{v} (n := n) hρ hρR hε
  refine ⟨C, hC0, fun j f hf => ?_⟩
  have hj : ContDiff ℝ 2 (jet j f) := by
    exact_mod_cast contDiff_infty.mp (contDiff_jet hf j) 2
  have h := hC (PiLp 2 (fun _ : Fin j → Fin n => V)) (jet j f) hj
  have e1 : ∀ x, gradSq (jet j f) x = ‖jet (j + 1) f x‖ ^ 2 := fun x =>
    sum_norm_fderiv_jet_sq hf j x
  have e2 : ∀ x, hessSq (jet j f) x = ‖jet (j + 2) f x‖ ^ 2 := fun x =>
    sum_norm_fderiv_fderiv_jet_sq hf j x
  simp only [e1, e2] at h
  exact h

end Main

/-- info: 'CWFlow.Interp.exists_const_gradSq_ball_le' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_const_gradSq_ball_le

/-- info: 'CWFlow.Interp.exists_const_jet_three_orders' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_const_jet_three_orders

end CWFlow.Interp
