import CWFlow.Interp.BallHamilton
import CWFlow.Extension.Weights
import Mathlib.Algebra.Order.Chebyshev

/-!
# Integral bounds for weight functions on jets

For `φ : ℝⁿ → ℝ` let `a_j(x) = |D^{3+j} φ(x)|` (`aJet`).  We combine Hamilton's interpolation
inequality on balls (`CWFlow.Interp.exists_const_jet_rpow_ball`) with the pointwise comparison
inequalities for the weight functions `Om` to bound
`∫_{B_ρ} Ω_{2m}^{≤m}`, `∫_{B_ρ} (Ω_v^{≤ℓ})²` and `∫_{B_ρ} Ω_{2m-1}^{≤m}` in terms of the energy
`∑_{i ≤ m} ∫_{B_R} |D^{3+i} φ|²`, for `φ` with `|D³ φ|` bounded on `B_R`.
-/

open Metric MeasureTheory CWFlow.Interp
open scoped ContDiff

namespace CWFlow.Extension

variable {n : ℕ}

/-- `a_j(x) = |D^{3+j} φ(x)|`: the size of the jet of weight `j`. -/
noncomputable def aJet (φ : EuclideanSpace ℝ (Fin n) → ℝ) (x : EuclideanSpace ℝ (Fin n)) (j : ℕ) :
    ℝ :=
  ‖jet (3 + j) φ x‖

/-- `Σ_{i ≤ m} ∫_{B_R} |D^{3+i} φ|²`. -/
noncomputable def energy3 (m : ℕ) (R : ℝ) (φ : EuclideanSpace ℝ (Fin n) → ℝ) : ℝ :=
  ∑ i ∈ Finset.range (m + 1), ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) R, ‖jet (3 + i) φ x‖ ^ 2

theorem aJet_nonneg (φ : EuclideanSpace ℝ (Fin n) → ℝ) (x : EuclideanSpace ℝ (Fin n)) (j : ℕ) :
    0 ≤ aJet φ x j := norm_nonneg _

theorem norm_jet_three {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ) (i : ℕ)
    (x : EuclideanSpace ℝ (Fin n)) : ‖jet i (jet 3 φ) x‖ = ‖jet (3 + i) φ x‖ := by
  rw [norm_jet_jet hφ, Nat.add_comm]

theorem energy3_eq_ballEnergy {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ) (m : ℕ)
    (R : ℝ) : energy3 m R φ = ballEnergy m R (jet 3 φ) := by
  unfold energy3 ballEnergy
  refine Finset.sum_congr rfl fun i _ => ?_
  refine setIntegral_congr_fun measurableSet_ball fun x _ => ?_
  simp only [norm_jet_three hφ]

theorem energy3_nonneg (m : ℕ) (R : ℝ) (φ : EuclideanSpace ℝ (Fin n) → ℝ) :
    0 ≤ energy3 m R φ :=
  Finset.sum_nonneg fun i _ =>
    setIntegral_nonneg measurableSet_ball fun x _ => by positivity

theorem continuous_aJet {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ) (j : ℕ) :
    Continuous fun x => aJet φ x j :=
  (contDiff_jet hφ (3 + j)).continuous.norm

theorem integrableOn_ball_of_continuous {f : EuclideanSpace ℝ (Fin n) → ℝ} (hf : Continuous f)
    (ρ : ℝ) : IntegrableOn f (ball (0 : EuclideanSpace ℝ (Fin n)) ρ) :=
  (hf.continuousOn.integrableOn_compact (isCompact_closedBall 0 ρ)).mono_set
    ball_subset_closedBall

theorem integrableOn_const_ball (c : ℝ) (ρ : ℝ) :
    IntegrableOn (fun _ : EuclideanSpace ℝ (Fin n) => c) (ball (0 : EuclideanSpace ℝ (Fin n)) ρ) :=
  integrableOn_const measure_ball_lt_top.ne

theorem continuous_aJet_rpow {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ) (j : ℕ)
    {p : ℝ} (hp : 0 ≤ p) : Continuous fun x => aJet φ x j ^ p :=
  (continuous_aJet hφ j).rpow_const fun _ => Or.inr hp

theorem continuous_Om {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ) (v ℓ : ℕ) :
    Continuous fun x => Om v ℓ (aJet φ x) := by
  unfold Om
  refine continuous_const.add (continuous_finsetSum _ fun j _ => ?_)
  exact continuous_aJet_rpow hφ j (by positivity)

/-- Single powers: `∫_{B_ρ} |D^{3+j}φ|^{2m/j} ≤ C · energy`. -/
theorem exists_const_aJet_rpow {m j : ℕ} (hj : 1 ≤ j) (hjm : j ≤ m) {ρ R : ℝ} (hρ : 0 < ρ)
    (hρR : ρ < R) (K : ℝ) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ φ : EuclideanSpace ℝ (Fin n) → ℝ, ContDiff ℝ ∞ φ →
      (∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) R, ‖jet 3 φ x‖ ≤ K) →
      ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) ρ, aJet φ x j ^ (2 * (m : ℝ) / j) ≤
        C * energy3 m R φ := by
  obtain ⟨C, hC, h⟩ := exists_const_jet_rpow_ball (n := n)
    (V := PiLp 2 (fun _ : Fin 3 → Fin n => ℝ)) hj hjm hρ hρR K
  refine ⟨C, hC, fun φ hφ hK => ?_⟩
  have h1 := h (jet 3 φ) (contDiff_jet hφ 3) hK
  rw [← energy3_eq_ballEnergy hφ] at h1
  have e : ∀ x, aJet φ x j = ‖jet j (jet 3 φ) x‖ := fun x => (norm_jet_three hφ j x).symm
  simp only [e]
  exact h1

theorem integral_Om_two_mul_eq {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ) (m : ℕ)
    (ρ : ℝ) :
    ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) ρ, Om (2 * m) m (aJet φ x) =
      (volume.real (ball (0 : EuclideanSpace ℝ (Fin n)) ρ)) +
        ∑ j ∈ Finset.Icc 1 m,
          ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) ρ, aJet φ x j ^ (2 * (m : ℝ) / j) := by
  simp only [Om_two_mul]
  have hint : ∀ j ∈ Finset.Icc 1 m, IntegrableOn (fun x => aJet φ x j ^ (2 * (m : ℝ) / j))
      (ball (0 : EuclideanSpace ℝ (Fin n)) ρ) volume := fun j _ =>
    integrableOn_ball_of_continuous (continuous_aJet_rpow hφ j (by positivity)) ρ
  rw [integral_add (integrableOn_const_ball _ ρ) (integrable_finsetSum _ hint),
    integral_finsetSum _ hint]
  simp

/-- `∫_{B_ρ} Ω_{2m}^{≤m} ≤ C (1 + energy)`. -/
theorem exists_const_integral_Om_two_mul (m : ℕ) {ρ R : ℝ} (hρ : 0 < ρ) (hρR : ρ < R) (K : ℝ) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ φ : EuclideanSpace ℝ (Fin n) → ℝ, ContDiff ℝ ∞ φ →
      (∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) R, ‖jet 3 φ x‖ ≤ K) →
      ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) ρ, Om (2 * m) m (aJet φ x) ≤
        C * (1 + energy3 m R φ) := by
  have hex : ∀ j : ℕ, ∃ C : ℝ, 0 ≤ C ∧ (j ∈ Finset.Icc 1 m →
      ∀ φ : EuclideanSpace ℝ (Fin n) → ℝ, ContDiff ℝ ∞ φ →
      (∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) R, ‖jet 3 φ x‖ ≤ K) →
      ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) ρ, aJet φ x j ^ (2 * (m : ℝ) / j) ≤
        C * energy3 m R φ) := by
    intro j
    by_cases hj : j ∈ Finset.Icc 1 m
    · obtain ⟨C, hC, h⟩ := exists_const_aJet_rpow (n := n) (m := m) (Finset.mem_Icc.1 hj).1
        (Finset.mem_Icc.1 hj).2 hρ hρR K
      exact ⟨C, hC, fun _ => h⟩
    · exact ⟨0, le_rfl, fun h => absurd h hj⟩
  choose Cj hCj hbd using hex
  set S := ∑ j ∈ Finset.Icc 1 m, Cj j with hS
  have hS0 : 0 ≤ S := Finset.sum_nonneg fun j _ => hCj j
  set Vol := volume.real (ball (0 : EuclideanSpace ℝ (Fin n)) ρ) with hVol
  have hV0 : 0 ≤ Vol := measureReal_nonneg
  refine ⟨Vol + S, add_nonneg hV0 hS0, fun φ hφ hK => ?_⟩
  rw [integral_Om_two_mul_eq hφ]
  have hE := energy3_nonneg m R φ
  have h1 : ∑ j ∈ Finset.Icc 1 m,
      ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) ρ, aJet φ x j ^ (2 * (m : ℝ) / j) ≤
      S * energy3 m R φ := by
    rw [hS, Finset.sum_mul]
    exact Finset.sum_le_sum fun j hj => hbd j hj φ hφ hK
  nlinarith [mul_nonneg hV0 hE]

/-- Pointwise: `(Ω_v^{≤ℓ})² ≤ 2(v+1) Ω_{2v}^{≤v}`. -/
theorem Om_sq_le (v ℓ : ℕ) {a : ℕ → ℝ} (ha : ∀ j, 0 ≤ a j) :
    (Om v ℓ a) ^ 2 ≤ (2 * ((v : ℝ) + 1)) * Om (2 * v) v a := by
  set s := Finset.Icc 1 (min v ℓ) with hs
  set S := ∑ j ∈ s, a j ^ ((v : ℝ) / j) with hSdef
  have hS0 : 0 ≤ S := Finset.sum_nonneg fun j _ => Real.rpow_nonneg (ha j) _
  have hcs : S ^ 2 ≤ (s.card : ℝ) * ∑ j ∈ s, (a j ^ ((v : ℝ) / j)) ^ 2 :=
    sq_sum_le_card_mul_sum_sq
  have hcard : (s.card : ℝ) ≤ v := by
    rw [hs, Nat.card_Icc]
    have : min v ℓ + 1 - 1 ≤ v := by omega
    exact_mod_cast this
  have hsq : ∀ j ∈ s, (a j ^ ((v : ℝ) / j)) ^ 2 = a j ^ ((2 * v : ℝ) / j) := by
    intro j _
    rw [← Real.rpow_natCast, ← Real.rpow_mul (ha j)]
    congr 1
    push_cast
    ring
  have hsub : ∑ j ∈ s, a j ^ ((2 * v : ℝ) / j) ≤ ∑ j ∈ Finset.Icc 1 v, a j ^ ((2 * v : ℝ) / j) :=
    Finset.sum_le_sum_of_subset_of_nonneg
      (Finset.Icc_subset_Icc le_rfl (min_le_left _ _))
      fun j _ _ => Real.rpow_nonneg (ha j) _
  rw [Finset.sum_congr rfl hsq] at hcs
  have hT0 : 0 ≤ ∑ j ∈ s, a j ^ ((2 * v : ℝ) / j) :=
    Finset.sum_nonneg fun j _ => Real.rpow_nonneg (ha j) _
  rw [Om_two_mul]
  unfold Om
  rw [← hs, ← hSdef]
  have hv : (0 : ℝ) ≤ v := by positivity
  nlinarith [sq_nonneg (1 - S), mul_le_mul_of_nonneg_right hcard hT0]

/-- `∫_{B_ρ} (Ω_v^{≤ℓ})² ≤ C (1 + energy of order v)`. -/
theorem exists_const_integral_Om_sq (v ℓ : ℕ) {ρ R : ℝ} (hρ : 0 < ρ) (hρR : ρ < R) (K : ℝ) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ φ : EuclideanSpace ℝ (Fin n) → ℝ, ContDiff ℝ ∞ φ →
      (∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) R, ‖jet 3 φ x‖ ≤ K) →
      ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) ρ, (Om v ℓ (aJet φ x)) ^ 2 ≤
        C * (1 + energy3 v R φ) := by
  obtain ⟨C, hC, h⟩ := exists_const_integral_Om_two_mul (n := n) v hρ hρR K
  have hc : 0 ≤ 2 * ((v : ℝ) + 1) := by positivity
  refine ⟨2 * ((v : ℝ) + 1) * C, mul_nonneg hc hC, fun φ hφ hK => ?_⟩
  calc ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) ρ, (Om v ℓ (aJet φ x)) ^ 2
      ≤ ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) ρ, 2 * ((v : ℝ) + 1) * Om (2 * v) v (aJet φ x) := by
        refine setIntegral_mono_on ?_ ?_ measurableSet_ball fun x _ => Om_sq_le v ℓ (aJet_nonneg φ x)
        · exact integrableOn_ball_of_continuous ((continuous_Om hφ v ℓ).pow 2) ρ
        · exact (integrableOn_ball_of_continuous (continuous_Om hφ (2 * v) v) ρ).const_mul _
    _ = 2 * ((v : ℝ) + 1) * ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) ρ, Om (2 * v) v (aJet φ x) :=
        integral_const_mul _ _
    _ ≤ 2 * ((v : ℝ) + 1) * (C * (1 + energy3 v R φ)) :=
        mul_le_mul_of_nonneg_left (h φ hφ hK) hc
    _ = _ := by ring

/-- The `ε`-bound for odd total weight: `∫_{B_ρ} Ω_{2m-1}^{≤m} ≤ ε · energy + C`. -/
theorem exists_const_integral_Om_odd {m : ℕ} (hm : 1 ≤ m) {ρ R : ℝ} (hρ : 0 < ρ) (hρR : ρ < R)
    (K : ℝ) {ε : ℝ} (hε : 0 < ε) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ φ : EuclideanSpace ℝ (Fin n) → ℝ, ContDiff ℝ ∞ φ →
      (∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) R, ‖jet 3 φ x‖ ≤ K) →
      ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) ρ, Om (2 * m - 1) m (aJet φ x) ≤
        ε * energy3 m R φ + C := by
  have hex : ∀ j : ℕ, ∃ C : ℝ, 0 ≤ C ∧ (j ∈ Finset.Icc 1 m →
      ∀ φ : EuclideanSpace ℝ (Fin n) → ℝ, ContDiff ℝ ∞ φ →
      (∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) R, ‖jet 3 φ x‖ ≤ K) →
      ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) ρ, aJet φ x j ^ (2 * (m : ℝ) / j) ≤
        C * energy3 m R φ) := by
    intro j
    by_cases hj : j ∈ Finset.Icc 1 m
    · obtain ⟨C, hC, h⟩ := exists_const_aJet_rpow (n := n) (m := m) (Finset.mem_Icc.1 hj).1
        (Finset.mem_Icc.1 hj).2 hρ hρR K
      exact ⟨C, hC, fun _ => h⟩
    · exact ⟨0, le_rfl, fun h => absurd h hj⟩
  choose Cj hCj hbd using hex
  set S := ∑ j ∈ Finset.Icc 1 m, Cj j with hS
  have hS0 : 0 ≤ S := Finset.sum_nonneg fun j _ => hCj j
  have hpos : 0 < 1 + S := by linarith
  set ε' := ε / (1 + S) with hε'
  have hε'0 : 0 < ε' := div_pos hε hpos
  have hε'S : ε' * S ≤ ε := by
    rw [hε', div_mul_eq_mul_div, div_le_iff₀ hpos]
    nlinarith
  obtain ⟨C', hC', hpt⟩ := exists_Om_le_eps hm hε'0
  set Vol := volume.real (ball (0 : EuclideanSpace ℝ (Fin n)) ρ) with hVol
  have hV0 : 0 ≤ Vol := measureReal_nonneg
  refine ⟨C' * Vol, mul_nonneg hC' hV0, fun φ hφ hK => ?_⟩
  have hint : ∀ j ∈ Finset.Icc 1 m, IntegrableOn (fun x => aJet φ x j ^ (2 * (m : ℝ) / j))
      (ball (0 : EuclideanSpace ℝ (Fin n)) ρ) volume := fun j _ =>
    integrableOn_ball_of_continuous (continuous_aJet_rpow hφ j (by positivity)) ρ
  have hE := energy3_nonneg m R φ
  have h1 : ∑ j ∈ Finset.Icc 1 m,
      ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) ρ, aJet φ x j ^ (2 * (m : ℝ) / j) ≤
      S * energy3 m R φ := by
    rw [hS, Finset.sum_mul]
    exact Finset.sum_le_sum fun j hj => hbd j hj φ hφ hK
  have hsumint : IntegrableOn (fun x => ∑ j ∈ Finset.Icc 1 m, aJet φ x j ^ (2 * (m : ℝ) / j))
      (ball (0 : EuclideanSpace ℝ (Fin n)) ρ) volume := integrable_finsetSum _ hint
  calc ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) ρ, Om (2 * m - 1) m (aJet φ x)
      ≤ ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) ρ,
          (ε' * ∑ j ∈ Finset.Icc 1 m, aJet φ x j ^ (2 * (m : ℝ) / j) + C') := by
        refine setIntegral_mono_on ?_ ?_ measurableSet_ball fun x _ => hpt _ (aJet_nonneg φ x)
        · exact integrableOn_ball_of_continuous (continuous_Om hφ _ _) ρ
        · exact (hsumint.const_mul _).add (integrableOn_const_ball _ ρ)
    _ = ε' * (∑ j ∈ Finset.Icc 1 m,
          ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) ρ, aJet φ x j ^ (2 * (m : ℝ) / j)) +
          C' * Vol := by
        rw [integral_add (hsumint.const_mul _) (integrableOn_const_ball _ ρ),
          integral_const_mul, integral_finsetSum _ hint]
        simp [hVol, mul_comm]
    _ ≤ ε * energy3 m R φ + C' * Vol := by
        have : ε' * (∑ j ∈ Finset.Icc 1 m,
          ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) ρ, aJet φ x j ^ (2 * (m : ℝ) / j)) ≤
            ε * energy3 m R φ := by
          calc _ ≤ ε' * (S * energy3 m R φ) := mul_le_mul_of_nonneg_left h1 hε'0.le
            _ = (ε' * S) * energy3 m R φ := by ring
            _ ≤ ε * energy3 m R φ := mul_le_mul_of_nonneg_right hε'S hE
        linarith

/-- info: 'CWFlow.Extension.exists_const_aJet_rpow' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_const_aJet_rpow

/-- info: 'CWFlow.Extension.exists_const_integral_Om_sq' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_const_integral_Om_sq

/-- info: 'CWFlow.Extension.exists_const_integral_Om_odd' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_const_integral_Om_odd

end CWFlow.Extension
