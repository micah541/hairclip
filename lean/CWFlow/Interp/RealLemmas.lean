import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.Analysis.SpecificLimits.Basic
import Mathlib.Tactic

/-!
# Two elementary real lemmas for interpolation inequalities

* `discrete_convex_le`: a real sequence satisfying `2 a j ≤ a (j-1) + a (j+1) + c` on the interior
  of `{0, …, m}` lies below the linear interpolation of its end values plus the parabola
  `c j (m - j) / 2` (discrete maximum principle).
* `exists_iteration_const`: the standard absorbing iteration on nested intervals:
  `Z s ≤ ½ Z t + A / (t - s)^α + B` for all `ρ ≤ s < t ≤ R` and `Z` bounded implies
  `Z ρ ≤ C (A / (R - ρ)^α + B)` with `C` depending only on `α`.
-/

namespace CWFlow.Interp

theorem discrete_convex_le {m : ℕ} (hm : 0 < m) {a : ℕ → ℝ} {c : ℝ}
    (h : ∀ j, 0 < j → j < m → 2 * a j ≤ a (j - 1) + a (j + 1) + c) :
    ∀ j, j ≤ m → a j ≤ (1 - (j : ℝ) / m) * a 0 + ((j : ℝ) / m) * a m + c * j * (m - j) / 2 := by
  classical
  have hmr : (m : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr hm.ne'
  obtain ⟨b, hb⟩ : ∃ b : ℕ → ℝ, ∀ j, b j = a j - (1 - (j : ℝ) / m) * a 0 - ((j : ℝ) / m) * a m
      - c * j * (m - j) / 2 := ⟨_, fun _ => rfl⟩
  have hb0 : b 0 = 0 := by rw [hb]; simp
  have hbm : b m = 0 := by rw [hb, div_self hmr]; ring
  have hconv : ∀ k, k + 1 < m → 2 * b (k + 1) ≤ b k + b (k + 1 + 1) := by
    intro k hk
    have h1 := h (k + 1) (by omega) hk
    simp only [Nat.add_sub_cancel] at h1
    have key : b k + b (k + 1 + 1) - 2 * b (k + 1)
        = a k + a (k + 1 + 1) - 2 * a (k + 1) + c := by
      rw [hb, hb, hb]; push_cast; ring
    linarith
  intro j hj
  by_contra hlt
  push_neg at hlt
  have hpos : 0 < b j := by rw [hb]; linarith
  obtain ⟨j0, hj0, hmax⟩ := Finset.exists_max_image (Finset.range (m + 1)) b ⟨0, by simp⟩
  have hP : ∃ i, i ≤ m ∧ ∀ x, x ≤ m → b x ≤ b i :=
    ⟨j0, by simpa [Nat.lt_succ_iff] using hj0,
      fun x hx => hmax x (by simp; omega)⟩
  obtain ⟨i, ⟨hi, himax⟩, hmin⟩ : ∃ i, (i ≤ m ∧ ∀ x, x ≤ m → b x ≤ b i) ∧
      ∀ n, n < i → ¬ (n ≤ m ∧ ∀ x, x ≤ m → b x ≤ b n) :=
    ⟨Nat.find hP, Nat.find_spec hP, fun n hn => Nat.find_min hP hn⟩
  have hbi : 0 < b i := lt_of_lt_of_le hpos (himax j hj)
  have hi0 : i ≠ 0 := by rintro rfl; linarith
  have him : i ≠ m := by rintro rfl; linarith
  obtain ⟨k, rfl⟩ : ∃ k, i = k + 1 := ⟨i - 1, by omega⟩
  have hlt1 : b k < b (k + 1) := by
    by_contra hc
    push_neg at hc
    exact hmin k (by omega) ⟨by omega, fun x hx => le_trans (himax x hx) hc⟩
  have hle2 : b (k + 1 + 1) ≤ b (k + 1) := himax _ (by omega)
  have := hconv k (by omega)
  linarith

/-- Explicit-constant version of the iteration lemma. -/
theorem iteration_explicit {α : ℝ} (hα : 0 < α) {τ : ℝ} (hτ0 : 0 < τ) (hτ1 : τ < 1)
    (hτα : τ ^ α = 3 / 4) (ρ R A B : ℝ) (Z : ℝ → ℝ) (hρR : ρ < R) (hA : 0 ≤ A) (hB : 0 ≤ B)
    (hM : ∃ M, ∀ s ∈ Set.Icc ρ R, Z s ≤ M)
    (hZ : ∀ s t, ρ ≤ s → s < t → t ≤ R → Z s ≤ (1 / 2) * Z t + A / (t - s) ^ α + B) :
    Z ρ ≤ 3 * (A / ((1 - τ) ^ α * (R - ρ) ^ α) + B) := by
  obtain ⟨M, hM⟩ := hM
  set L : ℝ := R - ρ with hL
  have hLpos : 0 < L := by rw [hL]; linarith
  have h1τ : 0 < 1 - τ := by linarith
  have hP : 0 < (1 - τ) ^ α := Real.rpow_pos_of_pos h1τ α
  have hQ : 0 < L ^ α := Real.rpow_pos_of_pos hLpos α
  obtain ⟨t, ht⟩ : ∃ t : ℕ → ℝ, ∀ i, t i = ρ + (1 - τ ^ i) * L := ⟨_, fun _ => rfl⟩
  set D : ℝ := A / ((1 - τ) ^ α * L ^ α) with hD
  have hD0 : 0 ≤ D := by positivity
  set K : ℝ := D + B with hK
  have hK0 : 0 ≤ K := by positivity
  have hτi : ∀ i : ℕ, 0 < τ ^ i := fun i => pow_pos hτ0 i
  have hτi1 : ∀ i : ℕ, τ ^ i ≤ 1 := fun i => pow_le_one₀ hτ0.le hτ1.le
  have ht0 : t 0 = ρ := by rw [ht]; simp
  have htρ : ∀ i, ρ ≤ t i := by
    intro i; rw [ht]
    have : 0 ≤ (1 - τ ^ i) * L := mul_nonneg (by linarith [hτi1 i]) hLpos.le
    linarith
  have htR : ∀ i, t i ≤ R := by
    intro i; rw [ht]
    have : (1 - τ ^ i) * L ≤ 1 * L :=
      mul_le_mul_of_nonneg_right (by linarith [hτi i]) hLpos.le
    rw [hL] at this ⊢; linarith
  have htdiff : ∀ i, t (i + 1) - t i = (1 - τ) * τ ^ i * L := by
    intro i; rw [ht, ht]; ring
  have htlt : ∀ i, t i < t (i + 1) := by
    intro i
    have : 0 < (1 - τ) * τ ^ i * L := by have := hτi i; positivity
    linarith [htdiff i]
  have hτpow : ∀ i : ℕ, (τ ^ i) ^ α = (3 / 4 : ℝ) ^ i := by
    intro i
    rw [← Real.rpow_natCast τ i, ← Real.rpow_mul hτ0.le, mul_comm, Real.rpow_mul hτ0.le,
      hτα, Real.rpow_natCast]
  -- one step
  have hstep : ∀ i : ℕ, Z (t i) ≤ (1 / 2) * Z (t (i + 1)) + (4 / 3 : ℝ) ^ i * K := by
    intro i
    have h1 := hZ (t i) (t (i + 1)) (htρ i) (htlt i) (htR (i + 1))
    rw [htdiff i, Real.mul_rpow (mul_nonneg h1τ.le (hτi i).le) hLpos.le,
      Real.mul_rpow h1τ.le (hτi i).le, hτpow i] at h1
    have hq : (0 : ℝ) < (3 / 4) ^ i := by positivity
    have hinv : (4 / 3 : ℝ) ^ i = ((3 / 4 : ℝ) ^ i)⁻¹ := by
      apply eq_inv_of_mul_eq_one_right
      rw [← mul_pow]; norm_num
    have heq : A / ((1 - τ) ^ α * (3 / 4 : ℝ) ^ i * L ^ α) = (4 / 3 : ℝ) ^ i * D := by
      rw [hinv, hD]
      field_simp
    have hB' : B ≤ (4 / 3 : ℝ) ^ i * B := by
      have : (1 : ℝ) ≤ (4 / 3) ^ i := one_le_pow₀ (by norm_num)
      nlinarith
    rw [heq] at h1
    have : (4 / 3 : ℝ) ^ i * K = (4 / 3 : ℝ) ^ i * D + (4 / 3 : ℝ) ^ i * B := by rw [hK]; ring
    linarith
  -- induction
  have hind : ∀ k : ℕ, Z ρ ≤ (1 / 2 : ℝ) ^ k * Z (t k)
      + K * ∑ i ∈ Finset.range k, (2 / 3 : ℝ) ^ i := by
    intro k
    induction k with
    | zero => simp [ht0]
    | succ k ih =>
      have hs := hstep k
      have hp : (0 : ℝ) < (1 / 2) ^ k := by positivity
      have h2 : (1 / 2 : ℝ) ^ k * Z (t k) ≤ (1 / 2) ^ k * ((1 / 2) * Z (t (k + 1))
          + (4 / 3 : ℝ) ^ k * K) := mul_le_mul_of_nonneg_left hs hp.le
      have h3 : (1 / 2 : ℝ) ^ k * (4 / 3 : ℝ) ^ k = (2 / 3 : ℝ) ^ k := by
        rw [← mul_pow]; norm_num
      rw [Finset.sum_range_succ, pow_succ]
      have h4 : (1 / 2 : ℝ) ^ k * ((1 / 2) * Z (t (k + 1)) + (4 / 3 : ℝ) ^ k * K)
          = (1 / 2 : ℝ) ^ k * (1 / 2) * Z (t (k + 1)) + (2 / 3 : ℝ) ^ k * K := by
        rw [← h3]; ring
      nlinarith [h2, h4, ih]
  have hsum : ∀ k : ℕ, ∑ i ∈ Finset.range k, (2 / 3 : ℝ) ^ i = 3 - 3 * (2 / 3 : ℝ) ^ k := by
    intro k
    induction k with
    | zero => norm_num
    | succ k ih => rw [Finset.sum_range_succ, ih, pow_succ]; ring
  have hsum3 : ∀ k : ℕ, K * ∑ i ∈ Finset.range k, (2 / 3 : ℝ) ^ i ≤ 3 * K := by
    intro k
    rw [hsum]
    have : 0 ≤ (2 / 3 : ℝ) ^ k := by positivity
    nlinarith
  set M' : ℝ := max M 0 with hM'
  have hlim : Filter.Tendsto (fun k : ℕ => (1 / 2 : ℝ) ^ k * M' + 3 * K) Filter.atTop
      (nhds (0 * M' + 3 * K)) :=
    ((tendsto_pow_atTop_nhds_zero_of_lt_one (by norm_num) (by norm_num)).mul_const M').add_const _
  have hbound : ∀ k : ℕ, Z ρ ≤ (1 / 2 : ℝ) ^ k * M' + 3 * K := by
    intro k
    have hZk : Z (t k) ≤ M' := le_trans (hM _ ⟨htρ k, htR k⟩) (le_max_left _ _)
    have hp : (0 : ℝ) < (1 / 2) ^ k := by positivity
    have := mul_le_mul_of_nonneg_left hZk hp.le
    linarith [hind k, hsum3 k]
  have := ge_of_tendsto' hlim hbound
  simpa [hK, hD] using this

theorem exists_iteration_const {α : ℝ} (hα : 0 < α) :
    ∃ C : ℝ, 0 < C ∧ ∀ (ρ R A B : ℝ) (Z : ℝ → ℝ), ρ < R → 0 ≤ A → 0 ≤ B →
      (∃ M, ∀ s ∈ Set.Icc ρ R, Z s ≤ M) →
      (∀ s t, ρ ≤ s → s < t → t ≤ R → Z s ≤ (1 / 2) * Z t + A / (t - s) ^ α + B) →
      Z ρ ≤ C * (A / (R - ρ) ^ α + B) := by
  set τ : ℝ := (3 / 4 : ℝ) ^ (1 / α) with hτ
  have hτ0 : 0 < τ := Real.rpow_pos_of_pos (by norm_num) _
  have hτ1 : τ < 1 := Real.rpow_lt_one (by norm_num) (by norm_num) (by positivity)
  have hτα : τ ^ α = 3 / 4 := by
    rw [hτ, ← Real.rpow_mul (by norm_num), one_div_mul_cancel hα.ne', Real.rpow_one]
  have h1τ : 0 < 1 - τ := by linarith
  have hP : 0 < (1 - τ) ^ α := Real.rpow_pos_of_pos h1τ α
  set c1 : ℝ := ((1 - τ) ^ α)⁻¹ with hc1
  have hc10 : 0 < c1 := inv_pos.mpr hP
  refine ⟨3 * max 1 c1, by positivity, ?_⟩
  intro ρ R A B Z hρR hA hB hM hZ
  have h := iteration_explicit hα hτ0 hτ1 hτα ρ R A B Z hρR hA hB hM hZ
  have hQ : 0 < (R - ρ) ^ α := Real.rpow_pos_of_pos (by linarith) α
  have hE : 0 ≤ A / (R - ρ) ^ α := by positivity
  have hrw : A / ((1 - τ) ^ α * (R - ρ) ^ α) = c1 * (A / (R - ρ) ^ α) := by
    rw [hc1]; field_simp
  rw [hrw] at h
  have h1 : c1 * (A / (R - ρ) ^ α) ≤ max 1 c1 * (A / (R - ρ) ^ α) :=
    mul_le_mul_of_nonneg_right (le_max_right _ _) hE
  have h2 : B ≤ max 1 c1 * B := by
    have := mul_le_mul_of_nonneg_right (le_max_left 1 c1) hB
    linarith
  nlinarith

/-- info: 'CWFlow.Interp.discrete_convex_le' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms discrete_convex_le

/-- info: 'CWFlow.Interp.exists_iteration_const' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_iteration_const

end CWFlow.Interp
