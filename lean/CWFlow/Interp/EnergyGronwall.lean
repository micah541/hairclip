import Mathlib.Analysis.ODE.Gronwall
import Mathlib.Analysis.Calculus.Deriv.Basic
import Mathlib.Algebra.Order.BigOperators.Group.Finset

/-!
# Two elementary real-analysis lemmas for energy estimates

* `absorb_le`: the absorption step.  From (G1) `E' ≤ −λ𝔇 + ε(X + Y) + c₁(ε)(1 + S + Y)` for all
  `ε > 0`, (G2) `X ≤ A𝔇 + A(1 + Y + S)` and (G3) `Y ≤ ε X + c₃(ε)(1 + S)` for all `ε > 0`,
  conclude `E' ≤ C (1 + S)`.
* `bounded_of_deriv_le_sum`: Gronwall induction.  If `(E k)' ≤ C k (1 + Σ_{l ≤ k} E l)` for
  `k ≥ 1` and `E 0` is bounded, every `E k` is bounded on `[a, T)`.
-/

namespace CWFlow.Interp

open Set Finset

/-- From (G1) `E' ≤ −λ𝔇 + ε(X + Y) + c₁(ε)(1 + S + Y)` for all `ε > 0`, (G2) `X ≤ A𝔇 + A(1 + Y + S)` and
(G3) `Y ≤ ε X + c₃(ε)(1 + S)` for all `ε > 0`, with `𝔇, S ≥ 0`, conclude `E' ≤ C (1 + S)` with `C` depending only on
`λ, A, c₁, c₃`. -/
theorem absorb_le {lam A : ℝ} (hlam : 0 < lam) (hA : 0 < A) (c₁ c₃ : ℝ → ℝ) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ (E' X Y S D : ℝ), 0 ≤ D → 0 ≤ S → 0 ≤ X → 0 ≤ Y →
      (∀ ε > 0, 0 ≤ c₁ ε ∧ E' ≤ -lam * D + ε * (X + Y) + c₁ ε * (1 + S + Y)) →
      X ≤ A * D + A * (1 + Y + S) →
      (∀ ε > 0, 0 ≤ c₃ ε ∧ Y ≤ ε * X + c₃ ε * (1 + S)) →
      E' ≤ C * (1 + S) := by
  have hε₁ : 0 < lam / (8 * A) := by positivity
  have hε₂ : 0 < 1 / (2 * A) := by positivity
  set ε₁ := lam / (8 * A) with hε₁def
  set ε₂ := 1 / (2 * A) with hε₂def
  set m := ε₁ + c₁ ε₁ with hm
  set ε₃ := lam / (4 * A * m + 1) with hε₃def
  -- the constant, made nonnegative with absolute values
  refine ⟨|3 * lam / (8 * A) * (2 * A * (1 + c₃ ε₂)) + m * c₃ ε₃ + c₁ ε₁|, abs_nonneg _, ?_⟩
  intro E' X Y S D hD hS hX hY hG1 hG2 hG3
  obtain ⟨hc1, h1⟩ := hG1 ε₁ hε₁
  obtain ⟨hc3a, h3a⟩ := hG3 ε₂ hε₂
  have hm0 : 0 ≤ m := by rw [hm]; positivity
  have hε₃pos : 0 < ε₃ := by rw [hε₃def]; positivity
  obtain ⟨hc3b, h3b⟩ := hG3 ε₃ hε₃pos
  have hAε₂ : A * ε₂ = 1 / 2 := by rw [hε₂def]; field_simp
  have hAε₁ : 8 * A * ε₁ = lam := by rw [hε₁def]; field_simp
  -- X ≤ 2 A D + C₂ (1 + S)
  have hAY := mul_le_mul_of_nonneg_left h3a hA.le
  have hexp : A * (ε₂ * X + c₃ ε₂ * (1 + S)) = (A * ε₂) * X + A * c₃ ε₂ * (1 + S) := by ring
  rw [hexp, hAε₂] at hAY
  have hX2 : X ≤ 2 * A * D + 2 * A * (1 + c₃ ε₂) * (1 + S) := by linarith
  -- m * ε₃ ≤ lam / (4A)
  have hmε₃ : 4 * A * (m * ε₃) ≤ lam := by
    rw [hε₃def]
    have hpos : 0 < 4 * A * m + 1 := by positivity
    have : 4 * A * (m * (lam / (4 * A * m + 1))) = lam * ((4 * A * m) / (4 * A * m + 1)) := by
      field_simp
    rw [this]
    have : (4 * A * m) / (4 * A * m + 1) ≤ 1 := by
      rw [div_le_one hpos]; linarith
    exact mul_le_of_le_one_right hlam.le this
  have hmY := mul_le_mul_of_nonneg_left h3b hm0
  have hκ0 : 0 ≤ ε₁ + m * ε₃ := by positivity
  clear_value ε₁ ε₂ m ε₃
  have hκ : 8 * A * (ε₁ + m * ε₃) ≤ 3 * lam := by linarith
  have hmm : m * Y = ε₁ * Y + c₁ ε₁ * Y := by rw [hm]; ring
  have hE : E' ≤ -lam * D + (ε₁ + m * ε₃) * X + (m * c₃ ε₃ + c₁ ε₁) * (1 + S) := by
    linarith
  have hκX := mul_le_mul_of_nonneg_left hX2 hκ0
  have hC₂0 : 0 ≤ 2 * A * (1 + c₃ ε₂) := by positivity
  have hS1 : 0 ≤ 1 + S := by linarith
  have hκD := mul_le_mul_of_nonneg_right hκ hD
  have hk : ε₁ + m * ε₃ ≤ 3 * lam / (8 * A) := by
    rw [le_div_iff₀ (by positivity)]; linarith
  have hκC := mul_le_mul_of_nonneg_right hk (mul_nonneg hC₂0 hS1)
  have hlD := mul_nonneg hlam.le hD
  have hfin : E' ≤ (3 * lam / (8 * A) * (2 * A * (1 + c₃ ε₂)) + m * c₃ ε₃ + c₁ ε₁) * (1 + S) := by
    linarith
  exact hfin.trans (mul_le_mul_of_nonneg_right (le_abs_self _) hS1)

/-- One step of the Gronwall induction: a function continuous on `[a, T)`, differentiable on
`(a, T)` with `f' ≤ K f + ε` there (`K, ε ≥ 0`), is bounded on `[a, T)`. -/
theorem eg_bounded_step {a T : ℝ} (haT : a < T) {f D : ℝ → ℝ} {K ε : ℝ} (hK : 0 ≤ K) (hε : 0 ≤ ε)
    (hnn : ∀ t ∈ Ico a T, 0 ≤ f t) (hcont : ContinuousOn f (Ico a T)) (hder : ∀ t ∈ Ioo a T, HasDerivAt f (D t) t)
    (hle : ∀ t ∈ Ioo a T, D t ≤ K * f t + ε) : ∃ B, ∀ t ∈ Ico a T, f t ≤ B := by
  have hca : ContinuousWithinAt f (Ico a T) a := hcont a ⟨le_rfl, haT⟩
  obtain ⟨η, hη, hηf⟩ := Metric.continuousWithinAt_iff.1 hca 1 one_pos
  set a' := min (a + η / 2) ((a + T) / 2) with ha'
  have ha'a : a < a' := lt_min (by linarith) (by linarith)
  have ha'T : a' < T := lt_of_le_of_lt (min_le_right _ _) (by linarith)
  have hnear : ∀ x ∈ Icc a a', f x ≤ f a + 1 := by
    intro x hx
    have hxT : x ∈ Ico a T := ⟨hx.1, lt_of_le_of_lt hx.2 ha'T⟩
    have hd : dist x a < η := by
      rw [Real.dist_eq, abs_of_nonneg (by linarith [hx.1])]
      have := hx.2.trans (min_le_left _ _)
      linarith
    have := hηf hxT hd
    rw [Real.dist_eq, abs_lt] at this
    linarith [this.2]
  set δ := f a + 1 with hδ
  set L := gronwallBound δ K ε (T - a) with hL
  refine ⟨max δ L, fun t ht => ?_⟩
  by_cases hta : t ≤ a'
  · exact (hnear t ⟨ht.1, hta⟩).trans (le_max_left _ _)
  · replace hta := not_le.1 hta
    refine le_trans ?_ (le_max_right _ _)
    have hIcc : Icc a' t ⊆ Ico a T := fun x hx => ⟨by linarith [hx.1], lt_of_le_of_lt hx.2 ht.2⟩
    have hmain := le_gronwallBound_of_liminf_deriv_right_le (f := f) (f' := D) (δ := δ) (K := K)
      (ε := ε) (a := a') (b := t) (hcont.mono hIcc) ?_ ?_ ?_ t ⟨hta.le, le_rfl⟩
    · refine hmain.trans ?_
      have hδ0 : 0 ≤ δ := by
        have := hnn a ⟨le_rfl, haT⟩
        rw [hδ]; linarith
      exact gronwallBound_mono hδ0 hε hK (by linarith [ht.2])
    · intro x hx r hr
      have hxo : x ∈ Ioo a T := ⟨by linarith [hx.1], lt_of_lt_of_le hx.2 ht.2.le⟩
      have := (hder x hxo).hasDerivWithinAt.liminf_right_slope_le hr
      refine this.mono fun z hz => ?_
      rwa [slope_def_field, div_eq_inv_mul] at hz
    · exact hnear a' ⟨ha'a.le, le_rfl⟩
    · intro x hx
      exact hle x ⟨by linarith [hx.1], lt_of_lt_of_le hx.2 ht.2.le⟩

theorem bounded_of_deriv_le_sum {a T : ℝ} (haT : a < T) (E : ℕ → ℝ → ℝ) (D : ℕ → ℝ → ℝ) (C : ℕ → ℝ)
    (hC : ∀ k, 0 ≤ C k) (hE0 : ∃ B, ∀ t ∈ Ico a T, E 0 t ≤ B)
    (hnn : ∀ k, ∀ t ∈ Ico a T, 0 ≤ E k t)
    (hcont : ∀ k, ContinuousOn (E k) (Ico a T))
    (hder : ∀ k, 1 ≤ k → ∀ t ∈ Ioo a T, HasDerivAt (E k) (D k t) t)
    (hle : ∀ k, 1 ≤ k → ∀ t ∈ Ioo a T, D k t ≤ C k * (1 + ∑ l ∈ range (k + 1), E l t)) :
    ∀ k, ∃ B, ∀ t ∈ Ico a T, E k t ≤ B := by
  have key : ∀ k, ∃ B : ℝ, 0 ≤ B ∧ ∀ l ≤ k, ∀ t ∈ Ico a T, E l t ≤ B := by
    intro k
    induction k with
    | zero =>
      obtain ⟨B, hB⟩ := hE0
      refine ⟨max B 0, le_max_right _ _, fun l hl t ht => ?_⟩
      obtain rfl : l = 0 := by omega
      exact (hB t ht).trans (le_max_left _ _)
    | succ k ih =>
      obtain ⟨B, hB0, hB⟩ := ih
      have hbd : ∀ t ∈ Ioo a T, ∑ l ∈ range (k + 1), E l t ≤ (k + 1 : ℝ) * B := by
        intro t ht
        have := Finset.sum_le_card_nsmul (range (k + 1)) (fun l => E l t) B
          (fun l hl => hB l (by simpa [Nat.lt_succ_iff] using hl) t ⟨ht.1.le, ht.2⟩)
        simpa [nsmul_eq_mul] using this
      obtain ⟨B', hB'⟩ := eg_bounded_step haT (K := C (k + 1)) (ε := C (k + 1) * (1 + (k + 1 : ℝ) * B))
        (f := E (k + 1)) (D := D (k + 1)) (hC _) (mul_nonneg (hC _) (by
          have : (0 : ℝ) ≤ k := k.cast_nonneg
          nlinarith))
        (hnn _) (hcont _) (hder _ (by omega)) (by
          intro t ht
          refine (hle (k + 1) (by omega) t ht).trans ?_
          rw [Finset.sum_range_succ]
          have h2 : 1 + (∑ l ∈ range (k + 1), E l t + E (k + 1) t) ≤
              1 + (k + 1 : ℝ) * B + E (k + 1) t := by linarith [hbd t ht]
          have := mul_le_mul_of_nonneg_left h2 (hC (k + 1))
          linarith [this])
      refine ⟨max B B', hB0.trans (le_max_left _ _), fun l hl t ht => ?_⟩
      rcases Nat.lt_or_ge l (k + 1) with h | h
      · exact (hB l (by omega) t ht).trans (le_max_left _ _)
      · obtain rfl : l = k + 1 := by omega
        exact (hB' t ht).trans (le_max_right _ _)
  intro k
  obtain ⟨B, -, hB⟩ := key k
  exact ⟨B, hB k le_rfl⟩

/-- info: 'CWFlow.Interp.absorb_le' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms absorb_le

/-- info: 'CWFlow.Interp.bounded_of_deriv_le_sum' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms bounded_of_deriv_le_sum

end CWFlow.Interp
