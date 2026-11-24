import CWFlow.Interp.HamiltonStep
import CWFlow.Interp.Jets
import CWFlow.Interp.LogConvex

/-!
# Hamilton's interpolation inequality

For `ξ : ℝⁿ → V` smooth with compact support, `|ξ| ≤ M`, and `1 ≤ j ≤ m`:
`∫ |D^j ξ|^(2m/j) ≤ (n + 2m)^(m (m - j)) M^(2m/j - 2) ∫ |D^m ξ|²`
(`hamilton_interp`).  It is obtained from the one-step inequality `hamilton_step` applied to the
jets `jet k ξ`, whose norms satisfy `N i ^ 2 ≤ c N (i-1) N (i+1)` for
`N i = (∫ |D^i ξ|^(2m/i))^(i/(2m))`, and from `discrete_logconvex_le`.
-/

open MeasureTheory
open scoped ContDiff

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false

namespace CWFlow.Interp

variable {n : ℕ} {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V] [FiniteDimensional ℝ V]

lemma sq_rpow_eq {y : ℝ} (hy : 0 ≤ y) (s : ℝ) : (y ^ 2) ^ s = y ^ (2 * s) := by
  rw [← Real.rpow_natCast y 2, ← Real.rpow_mul hy]; norm_num

/-- The normalised `L^{2m/i}` norms of the `i`-jet, with `N 0 = M`. -/
noncomputable def Nseq (ξ : EuclideanSpace ℝ (Fin n) → V) (m : ℕ) (M : ℝ) (i : ℕ) : ℝ :=
  if i = 0 then M else (∫ x, ‖jet i ξ x‖ ^ (2 * (m : ℝ) / i)) ^ ((i : ℝ) / (2 * m))

lemma Nseq_nonneg (ξ : EuclideanSpace ℝ (Fin n) → V) (m : ℕ) {M : ℝ} (hM : 0 ≤ M) (i : ℕ) :
    0 ≤ Nseq ξ m M i := by
  unfold Nseq
  split_ifs
  · exact hM
  · exact Real.rpow_nonneg (integral_nonneg fun x => Real.rpow_nonneg (norm_nonneg _) _) _

lemma integral_jet_nonneg (ξ : EuclideanSpace ℝ (Fin n) → V) (m : ℕ) (i : ℕ) :
    0 ≤ ∫ x, ‖jet i ξ x‖ ^ (2 * (m : ℝ) / i) :=
  integral_nonneg fun x => Real.rpow_nonneg (norm_nonneg _) _

lemma Nseq_step_mid {ξ : EuclideanSpace ℝ (Fin n) → V} (hξ : ContDiff ℝ ∞ ξ)
    (hc : HasCompactSupport ξ) {m k : ℕ} (hk : 1 ≤ k) (hkm : k + 1 ≤ m) (M : ℝ) :
    Nseq ξ m M (k + 1) ^ 2 ≤ ((n : ℝ) + 2 * m) * Nseq ξ m M k * Nseq ξ m M (k + 2) := by
  have hk0 : (0 : ℝ) < k := by exact_mod_cast hk
  have hm0 : (0 : ℝ) < m := by exact_mod_cast (by omega : 0 < m)
  have hkm' : (k : ℝ) + 1 ≤ m := by exact_mod_cast hkm
  have hs : 1 ≤ (m : ℝ) / ((k + 1 : ℕ) : ℝ) := by
    rw [le_div_iff₀ (by positivity)]; push_cast; linarith
  have hsm : (m : ℝ) / ((k + 1 : ℕ) : ℝ) ≤ m :=
    div_le_self hm0.le (by push_cast; linarith)
  have hη : ContDiff ℝ 2 (jet k ξ) := (contDiff_jet hξ k).of_le (by exact WithTop.coe_le_coe.mpr le_top)
  have hp : 0 < 2 * (m : ℝ) / ((k + 2 : ℕ) : ℝ) := by positivity
  have hq : 0 < 2 * (m : ℝ) / (k : ℝ) := by positivity
  have hpq : 1 / (2 * (m : ℝ) / ((k + 2 : ℕ) : ℝ)) + 1 / (2 * (m : ℝ) / (k : ℝ)) =
      1 / ((m : ℝ) / ((k + 1 : ℕ) : ℝ)) := by
    push_cast; field_simp; ring
  have h := hamilton_step hη (hasCompactSupport_jet hc k) hs hp hq hpq
  have hg : ∀ x, gradSq (jet k ξ) x = ‖jet (k + 1) ξ x‖ ^ 2 :=
    fun x => sum_norm_fderiv_jet_sq hξ k x
  have hh : ∀ x, hessSq (jet k ξ) x = ‖jet (k + 2) ξ x‖ ^ 2 :=
    fun x => sum_norm_fderiv_fderiv_jet_sq hξ k x
  have e1 : ∀ x, (‖jet (k + 1) ξ x‖ ^ 2) ^ ((m : ℝ) / ((k + 1 : ℕ) : ℝ)) =
      ‖jet (k + 1) ξ x‖ ^ (2 * (m : ℝ) / ((k + 1 : ℕ) : ℝ)) := by
    intro x; rw [sq_rpow_eq (norm_nonneg _)]; congr 1; ring
  have e2 : ∀ x, Real.sqrt (‖jet (k + 2) ξ x‖ ^ 2) ^ (2 * (m : ℝ) / ((k + 2 : ℕ) : ℝ)) =
      ‖jet (k + 2) ξ x‖ ^ (2 * (m : ℝ) / ((k + 2 : ℕ) : ℝ)) := by
    intro x; rw [Real.sqrt_sq (norm_nonneg _)]
  simp only [hg, hh] at h
  simp only [e1, e2] at h
  rw [one_div_div, one_div_div, one_div_div] at h
  have hB := integral_jet_nonneg ξ m k
  have hC := integral_jet_nonneg ξ m (k + 2)
  have hA := integral_jet_nonneg ξ m (k + 1)
  have hB' := Real.rpow_nonneg hB ((k : ℝ) / (2 * m))
  have hC' := Real.rpow_nonneg hC (((k + 2 : ℕ) : ℝ) / (2 * m))
  simp only [Nseq, if_neg (show k + 1 ≠ 0 by omega), if_neg (show k ≠ 0 by omega),
    if_neg (show k + 2 ≠ 0 by omega)]
  have hAe : ((∫ x, ‖jet (k + 1) ξ x‖ ^ (2 * (m : ℝ) / ((k + 1 : ℕ) : ℝ))) ^
      (((k + 1 : ℕ) : ℝ) / (2 * m))) ^ 2 =
      (∫ x, ‖jet (k + 1) ξ x‖ ^ (2 * (m : ℝ) / ((k + 1 : ℕ) : ℝ))) ^
        (((k + 1 : ℕ) : ℝ) / (m : ℝ)) := by
    rw [← Real.rpow_natCast, ← Real.rpow_mul hA]
    congr 1; push_cast; field_simp
  rw [hAe]
  have hconst : (n : ℝ) + 2 * ((m : ℝ) / ((k + 1 : ℕ) : ℝ)) - 2 ≤ (n : ℝ) + 2 * m := by
    linarith
  calc _ ≤ _ := h
    _ = ((n : ℝ) + 2 * ((m : ℝ) / ((k + 1 : ℕ) : ℝ)) - 2) *
        ((∫ x, ‖jet k ξ x‖ ^ (2 * (m : ℝ) / (k : ℝ))) ^ ((k : ℝ) / (2 * m)) *
          (∫ x, ‖jet (k + 2) ξ x‖ ^ (2 * (m : ℝ) / ((k + 2 : ℕ) : ℝ))) ^
            (((k + 2 : ℕ) : ℝ) / (2 * m))) := by ring
    _ ≤ ((n : ℝ) + 2 * m) * _ := mul_le_mul_of_nonneg_right hconst (mul_nonneg hB' hC')
    _ = _ := by ring

lemma Nseq_step_one {ξ : EuclideanSpace ℝ (Fin n) → V} (hξ : ContDiff ℝ ∞ ξ)
    (hc : HasCompactSupport ξ) {m : ℕ} (hm : 2 ≤ m) {M : ℝ} (hM : ∀ x, ‖ξ x‖ ≤ M) :
    Nseq ξ m M 1 ^ 2 ≤ ((n : ℝ) + 2 * m) * Nseq ξ m M 0 * Nseq ξ m M 2 := by
  have hM0 : 0 ≤ M := (norm_nonneg _).trans (hM 0)
  have hm0 : (0 : ℝ) < m := by exact_mod_cast (by omega : 0 < m)
  have hm1 : (1 : ℝ) ≤ m := by exact_mod_cast (by omega : 1 ≤ m)
  have hη : ContDiff ℝ 2 (jet 0 ξ) := (contDiff_jet hξ 0).of_le (by exact WithTop.coe_le_coe.mpr le_top)
  have hM' : ∀ x, ‖jet 0 ξ x‖ ≤ M := fun x => (norm_jet_zero ξ x).le.trans (hM x)
  have h := hamilton_step_sup hη (hasCompactSupport_jet hc 0) hm1 hM'
  have hg : ∀ x, gradSq (jet 0 ξ) x = ‖jet 1 ξ x‖ ^ 2 :=
    fun x => sum_norm_fderiv_jet_sq hξ 0 x
  have hh : ∀ x, hessSq (jet 0 ξ) x = ‖jet 2 ξ x‖ ^ 2 :=
    fun x => sum_norm_fderiv_fderiv_jet_sq hξ 0 x
  have e1 : ∀ x, (‖jet 1 ξ x‖ ^ 2) ^ (m : ℝ) =
      ‖jet 1 ξ x‖ ^ (2 * (m : ℝ) / ((1 : ℕ) : ℝ)) := by
    intro x; rw [sq_rpow_eq (norm_nonneg _)]; congr 1; simp
  have e2 : ∀ x, Real.sqrt (‖jet 2 ξ x‖ ^ 2) ^ (m : ℝ) =
      ‖jet 2 ξ x‖ ^ (2 * (m : ℝ) / ((2 : ℕ) : ℝ)) := by
    intro x; rw [Real.sqrt_sq (norm_nonneg _)]; congr 1; push_cast; ring
  simp only [hg, hh] at h
  simp only [e1, e2] at h
  have hA := integral_jet_nonneg ξ m 1
  have hC := integral_jet_nonneg ξ m 2
  have hC' := Real.rpow_nonneg hC (1 / (m : ℝ))
  have r2 : ((2 : ℕ) : ℝ) / (2 * m) = 1 / (m : ℝ) := by push_cast; field_simp
  have r1 : ((1 : ℕ) : ℝ) / (2 * m) * (2 : ℕ) = 1 / (m : ℝ) := by push_cast; field_simp
  have hN0 : Nseq ξ m M 0 = M := by simp [Nseq]
  rw [hN0]
  simp only [Nseq, if_neg (show (1 : ℕ) ≠ 0 by norm_num), if_neg (show (2 : ℕ) ≠ 0 by norm_num)]
  rw [r2, ← Real.rpow_natCast, ← Real.rpow_mul hA, r1]
  calc _ ≤ _ := h
    _ = ((n : ℝ) + 2 * m - 2) * (M * (∫ x, ‖jet 2 ξ x‖ ^ (2 * (m : ℝ) / ((2 : ℕ) : ℝ))) ^
        (1 / (m : ℝ))) := by ring
    _ ≤ ((n : ℝ) + 2 * m) * (M * (∫ x, ‖jet 2 ξ x‖ ^ (2 * (m : ℝ) / ((2 : ℕ) : ℝ))) ^
        (1 / (m : ℝ))) :=
      mul_le_mul_of_nonneg_right (by linarith) (mul_nonneg hM0 hC')
    _ = _ := by ring

theorem hamilton_interp {ξ : EuclideanSpace ℝ (Fin n) → V} (hξ : ContDiff ℝ ∞ ξ)
    (hc : HasCompactSupport ξ) {m j : ℕ} (hj : 1 ≤ j) (hjm : j ≤ m) {M : ℝ}
    (hM : ∀ x, ‖ξ x‖ ≤ M) :
    ∫ x, ‖jet j ξ x‖ ^ (2 * (m : ℝ) / j) ≤
      ((n : ℝ) + 2 * m) ^ ((m : ℝ) * (m - j)) * M ^ (2 * (m : ℝ) / j - 2) *
        ∫ x, ‖jet m ξ x‖ ^ 2 := by
  have hM0 : 0 ≤ M := (norm_nonneg _).trans (hM 0)
  have hm : 0 < m := by omega
  have hm0 : (0 : ℝ) < m := by exact_mod_cast hm
  have hj0 : (0 : ℝ) < j := by exact_mod_cast hj
  have hc0 : (0 : ℝ) < (n : ℝ) + 2 * m := by positivity
  have hN := Nseq_nonneg ξ m hM0
  have hstep : ∀ i, 0 < i → i < m →
      Nseq ξ m M i ^ 2 ≤ ((n : ℝ) + 2 * m) * Nseq ξ m M (i - 1) * Nseq ξ m M (i + 1) := by
    intro i hi him
    obtain ⟨k, rfl⟩ : ∃ k, i = k + 1 := ⟨i - 1, by omega⟩
    rw [Nat.add_sub_cancel]
    rcases Nat.eq_zero_or_pos k with rfl | hk
    · exact Nseq_step_one hξ hc (by omega) hM
    · exact Nseq_step_mid hξ hc hk (by omega) M
  have key := discrete_logconvex_le hm hc0 hN hstep j hjm
  have hNj : Nseq ξ m M j = (∫ x, ‖jet j ξ x‖ ^ (2 * (m : ℝ) / j)) ^ ((j : ℝ) / (2 * m)) := by
    simp only [Nseq, if_neg (show j ≠ 0 by omega)]
  have hN0 : Nseq ξ m M 0 = M := by simp [Nseq]
  have hIj := integral_jet_nonneg ξ m j
  have hNm : Nseq ξ m M m ^ 2 = ∫ x, ‖jet m ξ x‖ ^ 2 := by
    have hI := integral_jet_nonneg ξ m m
    simp only [Nseq, if_neg (show m ≠ 0 by omega)]
    have e : ∀ x, ‖jet m ξ x‖ ^ (2 * (m : ℝ) / (m : ℝ)) = ‖jet m ξ x‖ ^ 2 := by
      intro x
      rw [show 2 * (m : ℝ) / (m : ℝ) = (2 : ℕ) by field_simp; norm_num, Real.rpow_natCast]
    simp only [e] at hI ⊢
    rw [← Real.rpow_natCast, ← Real.rpow_mul hI]
    have : ((m : ℝ) / (2 * m)) * ((2 : ℕ) : ℝ) = 1 := by push_cast; field_simp
    rw [this, Real.rpow_one]
  have hNmnn := hN m
  have hmain : (Nseq ξ m M j) ^ (2 * (m : ℝ) / j) =
      ∫ x, ‖jet j ξ x‖ ^ (2 * (m : ℝ) / j) := by
    rw [hNj, ← Real.rpow_mul hIj]
    have : (j : ℝ) / (2 * m) * (2 * (m : ℝ) / j) = 1 := by field_simp
    rw [this, Real.rpow_one]
  have he : 0 ≤ 2 * (m : ℝ) / j := by positivity
  rw [← hmain]
  calc Nseq ξ m M j ^ (2 * (m : ℝ) / j)
      ≤ (((n : ℝ) + 2 * m) ^ ((j : ℝ) * (m - j) / 2) * Nseq ξ m M 0 ^ (1 - (j : ℝ) / m) *
        Nseq ξ m M m ^ ((j : ℝ) / m)) ^ (2 * (m : ℝ) / j) :=
        Real.rpow_le_rpow (hN j) key he
    _ = ((n : ℝ) + 2 * m) ^ ((m : ℝ) * (m - j)) * M ^ (2 * (m : ℝ) / j - 2) *
        (Nseq ξ m M m ^ 2) := by
      have hjm' : (j : ℝ) ≤ m := by exact_mod_cast hjm
      have hA : 0 ≤ ((n : ℝ) + 2 * m) ^ ((j : ℝ) * (m - j) / 2) := Real.rpow_nonneg hc0.le _
      have hB : 0 ≤ Nseq ξ m M 0 ^ (1 - (j : ℝ) / m) := Real.rpow_nonneg (hN 0) _
      have hC : 0 ≤ Nseq ξ m M m ^ ((j : ℝ) / m) := Real.rpow_nonneg hNmnn _
      rw [Real.mul_rpow (mul_nonneg hA hB) hC, Real.mul_rpow hA hB, ← Real.rpow_mul hc0.le,
        hN0, ← Real.rpow_mul hM0, ← Real.rpow_mul hNmnn]
      have x1 : (j : ℝ) * (m - j) / 2 * (2 * (m : ℝ) / j) = (m : ℝ) * (m - j) := by
        field_simp
      have x2 : (1 - (j : ℝ) / m) * (2 * (m : ℝ) / j) = 2 * (m : ℝ) / j - 2 := by
        field_simp
      have x3 : (j : ℝ) / m * (2 * (m : ℝ) / j) = ((2 : ℕ) : ℝ) := by
        push_cast; field_simp
      rw [x1, x2, x3, Real.rpow_natCast]
    _ = _ := by rw [hNm]

/-- info: 'CWFlow.Interp.hamilton_interp' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms hamilton_interp

/-! ### Products of jets -/

lemma jet_eq_zero_of_notMem {ξ : EuclideanSpace ℝ (Fin n) → V} (j : ℕ)
    {x : EuclideanSpace ℝ (Fin n)} (hx : x ∉ tsupport ξ) : jet j ξ x = 0 := by
  have h0 : iteratedFDeriv ℝ j ξ x = 0 := by
    by_contra h
    exact hx (support_iteratedFDeriv_subset (𝕜 := ℝ) j h)
  ext I
  simp [jet_apply, h0]

/-- Hölder's inequality for finitely many nonnegative functions and the indicator of a set `K`
of finite measure, with weights summing to one. -/
lemma integral_prod_holder {α ι : Type*} [MeasurableSpace α] {μ : Measure α} [Fintype ι]
    {K : Set α} (hK : MeasurableSet K) (hKμ : μ K ≠ ⊤) {f : ι → α → ℝ}
    (hf0 : ∀ i x, 0 ≤ f i x) (hfi : ∀ i, Integrable (f i) μ) {p : ι → ℝ} (hp : ∀ i, 0 ≤ p i)
    {w : ℝ} (hw : 0 ≤ w) (hsum : w + ∑ i, p i = 1)
    (hK0 : ∀ x, x ∉ K → ∏ i, f i x ^ p i = 0) :
    ∫ x, ∏ i, f i x ^ p i ∂μ ≤ (μ K).toReal ^ w * ∏ i, (∫ x, f i x ∂μ) ^ p i := by
  have hI : ∀ i, 0 ≤ ∫ x, f i x ∂μ := fun i => integral_nonneg (hf0 i)
  have hprod0 : ∀ x, 0 ≤ ∏ i, f i x ^ p i := fun x =>
    Finset.prod_nonneg fun i _ => Real.rpow_nonneg (hf0 i x) _
  have hr : 0 ≤ (μ K).toReal ^ w * ∏ i, (∫ x, f i x ∂μ) ^ p i :=
    mul_nonneg (Real.rpow_nonneg ENNReal.toReal_nonneg _)
      (Finset.prod_nonneg fun i _ => Real.rpow_nonneg (hI i) _)
  refine integral_le_of_lintegral_le hprod0 hr ?_
  have hP : ∀ x, ∏ i, ENNReal.ofReal (f i x) ^ p i = ENNReal.ofReal (∏ i, f i x ^ p i) := by
    intro x
    rw [ENNReal.ofReal_prod_of_nonneg (fun i _ => Real.rpow_nonneg (hf0 i x) _)]
    exact Finset.prod_congr rfl fun i _ => (ENNReal.ofReal_rpow_of_nonneg (hf0 i x) (hp i))
  have e : ∀ x, ENNReal.ofReal (∏ i, f i x ^ p i) =
      (K.indicator (1 : α → ENNReal) x) ^ w * ∏ i, ENNReal.ofReal (f i x) ^ p i := by
    intro x
    rw [hP]
    by_cases hx : x ∈ K
    · simp [hx]
    · rw [hK0 x hx]; simp
  have hL : ∀ i, ∫⁻ x, ENNReal.ofReal (f i x) ∂μ = ENNReal.ofReal (∫ x, f i x ∂μ) := fun i =>
    (ofReal_integral_eq_lintegral_ofReal (hfi i) (Filter.Eventually.of_forall (hf0 i))).symm
  simp_rw [e]
  refine (ENNReal.lintegral_mul_prod_norm_pow_le Finset.univ (g := K.indicator 1)
    (f := fun i x => ENNReal.ofReal (f i x)) (measurable_one.indicator hK).aemeasurable
    (fun i _ => (hfi i).aemeasurable.ennreal_ofReal) w (p := p) hsum hw (fun i _ => hp i)).trans
    (le_of_eq ?_)
  rw [lintegral_indicator_one hK]
  simp only [hL]
  rw [ENNReal.ofReal_mul (Real.rpow_nonneg ENNReal.toReal_nonneg _),
    ← ENNReal.ofReal_rpow_of_nonneg ENNReal.toReal_nonneg hw, ENNReal.ofReal_toReal hKμ,
    ENNReal.ofReal_prod_of_nonneg (fun i _ => Real.rpow_nonneg (hI i) _)]
  congr 1
  exact Finset.prod_congr rfl fun i _ => (ENNReal.ofReal_rpow_of_nonneg (hI i) (hp i))

theorem integral_prod_norm_jet_le {ξ : EuclideanSpace ℝ (Fin n) → V} (hξ : ContDiff ℝ ∞ ξ)
    (hc : HasCompactSupport ξ) {m : ℕ} (hm : 1 ≤ m) {q : ℕ} (hq : 1 ≤ q) (j : Fin q → ℕ)
    (hj1 : ∀ i, 1 ≤ j i) (hjm : ∀ i, j i ≤ m) (hsum : ∑ i, j i ≤ 2 * m) {M : ℝ}
    (hM : ∀ x, ‖ξ x‖ ≤ M) :
    ∫ x, ∏ i, ‖jet (j i) ξ x‖ ≤
      ((n : ℝ) + 2 * m) ^ ((q : ℝ) * m ^ 2) * M ^ ((q : ℝ) - (∑ i, (j i : ℝ)) / m) *
        (MeasureTheory.volume (tsupport ξ)).toReal ^ (1 - (∑ i, (j i : ℝ)) / (2 * m)) *
        (∫ x, ‖jet m ξ x‖ ^ 2) ^ ((∑ i, (j i : ℝ)) / (2 * m)) := by
  classical
  have hM0 : 0 ≤ M := (norm_nonneg _).trans (hM 0)
  have hm0 : (0 : ℝ) < m := by exact_mod_cast hm
  have hm1 : (1 : ℝ) ≤ m := by exact_mod_cast hm
  have hc0 : (1 : ℝ) ≤ (n : ℝ) + 2 * m := by
    have := (Nat.cast_nonneg n : (0 : ℝ) ≤ n); linarith
  have hcpos : (0 : ℝ) < (n : ℝ) + 2 * m := by linarith
  have hj0 : ∀ i, (0 : ℝ) < j i := fun i => by exact_mod_cast hj1 i
  have hjm' : ∀ i, (j i : ℝ) ≤ m := fun i => by exact_mod_cast hjm i
  have hS2 : (∑ i, (j i : ℝ)) ≤ 2 * m := by exact_mod_cast hsum
  set S : ℝ := ∑ i, (j i : ℝ) with hS
  set Q : ℝ := ∫ x, ‖jet m ξ x‖ ^ 2 with hQ
  have hQ0 : 0 ≤ Q := integral_nonneg fun x => by positivity
  have hfp : ∀ i x, (‖jet (j i) ξ x‖ ^ (2 * (m : ℝ) / (j i))) ^ ((j i : ℝ) / (2 * m)) =
      ‖jet (j i) ξ x‖ := by
    intro i x
    rw [← Real.rpow_mul (norm_nonneg _)]
    have : 2 * (m : ℝ) / (j i) * ((j i : ℝ) / (2 * m)) = 1 := by
      have := hj0 i; field_simp
    rw [this, Real.rpow_one]
  have hsump : ∑ i : Fin q, (j i : ℝ) / (2 * m) = S / (2 * m) := by rw [hS, Finset.sum_div]
  have hK0 : ∀ x, x ∉ tsupport ξ →
      ∏ i : Fin q, (‖jet (j i) ξ x‖ ^ (2 * (m : ℝ) / (j i))) ^ ((j i : ℝ) / (2 * m)) = 0 := by
    intro x hx
    simp only [hfp]
    exact Finset.prod_eq_zero (Finset.mem_univ (⟨0, hq⟩ : Fin q))
      (by rw [jet_eq_zero_of_notMem _ hx, norm_zero])
  have hH := integral_prod_holder (μ := volume) (K := tsupport ξ)
    (isClosed_tsupport ξ).measurableSet (hc : IsCompact (tsupport ξ)).measure_lt_top.ne
    (f := fun i x => ‖jet (j i) ξ x‖ ^ (2 * (m : ℝ) / (j i)))
    (fun i x => Real.rpow_nonneg (norm_nonneg _) _)
    (fun i => integrable_rpow_of_compact (contDiff_jet hξ (j i)).continuous.norm
      (hasCompactSupport_norm (hasCompactSupport_jet hc (j i)))
      (by have := hj0 i; positivity))
    (p := fun i => (j i : ℝ) / (2 * m)) (fun i => by have := hj0 i; positivity)
    (w := 1 - S / (2 * m)) (by rw [sub_nonneg, div_le_one (by positivity)]; exact hS2)
    (by rw [hsump]; ring) hK0
  simp only [hfp] at hH
  have hper : ∀ i, (∫ x, ‖jet (j i) ξ x‖ ^ (2 * (m : ℝ) / (j i))) ^ ((j i : ℝ) / (2 * m)) ≤
      ((n : ℝ) + 2 * m) ^ ((m : ℝ) ^ 2) * M ^ (1 - (j i : ℝ) / m) * Q ^ ((j i : ℝ) / (2 * m)) := by
    intro i
    have hji := hj0 i
    have hjmi := hjm' i
    have hst := hamilton_interp hξ hc (hj1 i) (hjm i) hM
    have hI := integral_jet_nonneg ξ m (j i)
    have hp : 0 ≤ (j i : ℝ) / (2 * m) := by positivity
    have hA : 0 ≤ ((n : ℝ) + 2 * m) ^ ((m : ℝ) * (m - j i)) := Real.rpow_nonneg hcpos.le _
    have hB : 0 ≤ M ^ (2 * (m : ℝ) / (j i) - 2) := Real.rpow_nonneg hM0 _
    have e : (2 * (m : ℝ) / (j i) - 2) * ((j i : ℝ) / (2 * m)) = 1 - (j i : ℝ) / m := by
      field_simp
    have hexp : (m : ℝ) * (m - j i) * ((j i : ℝ) / (2 * m)) ≤ (m : ℝ) ^ 2 := by
      have e2 : (m : ℝ) * (m - j i) * ((j i : ℝ) / (2 * m)) = (m - j i) * j i / 2 := by
        field_simp
      have : ((m : ℝ) - j i) * j i ≤ m * m :=
        mul_le_mul (by linarith) hjmi hji.le hm0.le
      rw [e2]; nlinarith
    calc _ ≤ (((n : ℝ) + 2 * m) ^ ((m : ℝ) * (m - j i)) * M ^ (2 * (m : ℝ) / (j i) - 2) * Q) ^
          ((j i : ℝ) / (2 * m)) := Real.rpow_le_rpow hI hst hp
      _ = ((n : ℝ) + 2 * m) ^ ((m : ℝ) * (m - j i) * ((j i : ℝ) / (2 * m))) *
          M ^ ((2 * (m : ℝ) / (j i) - 2) * ((j i : ℝ) / (2 * m))) * Q ^ ((j i : ℝ) / (2 * m)) := by
        rw [Real.mul_rpow (mul_nonneg hA hB) hQ0, Real.mul_rpow hA hB,
          ← Real.rpow_mul hcpos.le, ← Real.rpow_mul hM0]
      _ ≤ _ := by
        rw [e]
        exact mul_le_mul_of_nonneg_right
          (mul_le_mul_of_nonneg_right (Real.rpow_le_rpow_of_exponent_le hc0 hexp)
            (Real.rpow_nonneg hM0 _)) (Real.rpow_nonneg hQ0 _)
  have hcp : ∏ _i : Fin q, ((n : ℝ) + 2 * m) ^ ((m : ℝ) ^ 2) =
      ((n : ℝ) + 2 * m) ^ ((q : ℝ) * m ^ 2) := by
    rw [Finset.prod_const, Finset.card_univ, Fintype.card_fin, ← Real.rpow_natCast,
      ← Real.rpow_mul hcpos.le]
    congr 1; ring
  have hMp : ∏ i : Fin q, M ^ (1 - (j i : ℝ) / m) = M ^ ((q : ℝ) - S / m) := by
    rw [← Real.rpow_sum_of_nonneg hM0 (fun i _ => by
      have := hjm' i; rw [sub_nonneg, div_le_one hm0]; exact this)]
    congr 1
    rw [Finset.sum_sub_distrib, ← Finset.sum_div]
    simp [hS]
  have hQp : ∏ i : Fin q, Q ^ ((j i : ℝ) / (2 * m)) = Q ^ (S / (2 * m)) := by
    rw [← Real.rpow_sum_of_nonneg hQ0 (fun i _ => by have := hj0 i; positivity), hsump]
  calc ∫ x, ∏ i, ‖jet (j i) ξ x‖
      ≤ (volume (tsupport ξ)).toReal ^ (1 - S / (2 * m)) *
        ∏ i, (∫ x, ‖jet (j i) ξ x‖ ^ (2 * (m : ℝ) / (j i))) ^ ((j i : ℝ) / (2 * m)) := hH
    _ ≤ (volume (tsupport ξ)).toReal ^ (1 - S / (2 * m)) *
        ∏ i, (((n : ℝ) + 2 * m) ^ ((m : ℝ) ^ 2) * M ^ (1 - (j i : ℝ) / m) *
          Q ^ ((j i : ℝ) / (2 * m))) :=
      mul_le_mul_of_nonneg_left (Finset.prod_le_prod
        (fun i _ => Real.rpow_nonneg (integral_jet_nonneg ξ m (j i)) _) fun i _ => hper i)
        (Real.rpow_nonneg ENNReal.toReal_nonneg _)
    _ = _ := by
      rw [Finset.prod_mul_distrib, Finset.prod_mul_distrib, hcp, hMp, hQp]
      ring

/-- info: 'CWFlow.Interp.integral_prod_norm_jet_le' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms integral_prod_norm_jet_le

end CWFlow.Interp
