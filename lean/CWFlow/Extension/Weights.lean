import Mathlib.Analysis.MeanInequalities
import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.Algebra.BigOperators.Ring.Finset
import Mathlib.Order.Interval.Finset.Nat
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.GCongr

/-!
# Weight functions

For a sequence `a : ℕ → ℝ` of nonnegative numbers (think `a j = |D^{3+j} φ (x)|`, a quantity of
"weight `j`"), the weight function of total weight `v` with factors of weight at most `ℓ` is
`Ω_v^{≤ℓ}(a) = 1 + ∑_{1 ≤ j ≤ min v ℓ} a_j^{v/j}`.  We prove the elementary comparison inequalities
(monotonicity, products, absorption of a lower weight term) from Young's inequality.

In `exists_Om_le_eps` and `Om_two_mul` the real number `(2 * m : ℝ)` is `2 * (m : ℝ)`.
-/

namespace CWFlow.Extension

open Finset

/-- `Ω_v^{≤ℓ}(a) = 1 + Σ_{1 ≤ j ≤ min v ℓ} a_j^{v/j}`. -/
noncomputable def Om (v ℓ : ℕ) (a : ℕ → ℝ) : ℝ :=
  1 + ∑ j ∈ Finset.Icc 1 (min v ℓ), a j ^ ((v : ℝ) / j)

theorem Om_zero (ℓ : ℕ) (a : ℕ → ℝ) : Om 0 ℓ a = 1 := by
  simp [Om]

theorem one_le_Om (v ℓ : ℕ) {a : ℕ → ℝ} (ha : ∀ j, 0 ≤ a j) : 1 ≤ Om v ℓ a := by
  unfold Om
  have : 0 ≤ ∑ j ∈ Finset.Icc 1 (min v ℓ), a j ^ ((v : ℝ) / j) :=
    Finset.sum_nonneg fun j _ => Real.rpow_nonneg (ha j) _
  linarith

/-- `x ^ p ≤ 1 + x ^ p'` for `0 ≤ p ≤ p'`, `x ≥ 0`. -/
lemma rpow_le_one_add {x p p' : ℝ} (hx : 0 ≤ x) (hp : 0 ≤ p) (hpp : p ≤ p') :
    x ^ p ≤ 1 + x ^ p' := by
  have h' : 0 ≤ x ^ p' := Real.rpow_nonneg hx _
  rcases le_total x 1 with h | h
  · have := Real.rpow_le_one hx h hp
    linarith
  · have := Real.rpow_le_rpow_of_exponent_le h hpp
    linarith

/-- Young's inequality with exponents `(P+Q)/P`, `(P+Q)/Q`. -/
lemma young_aux {X Y P Q : ℝ} (hX : 0 ≤ X) (hY : 0 ≤ Y) (hP : 0 < P) (hQ : 0 < Q) :
    X * Y ≤ X ^ ((P + Q) / P) + Y ^ ((P + Q) / Q) := by
  have hs : 1 < (P + Q) / P := by rw [lt_div_iff₀ hP]; linarith
  have ht : 1 < (P + Q) / Q := by rw [lt_div_iff₀ hQ]; linarith
  have hne : P + Q ≠ 0 := (add_pos hP hQ).ne'
  have hc : ((P + Q) / P).HolderConjugate ((P + Q) / Q) := by
    rw [Real.holderConjugate_iff]
    refine ⟨hs, ?_⟩
    rw [inv_div, inv_div, ← add_div]
    exact div_self hne
  have h := Real.young_inequality_of_nonneg hX hY hc
  have h1 : X ^ ((P + Q) / P) / ((P + Q) / P) ≤ X ^ ((P + Q) / P) :=
    div_le_self (Real.rpow_nonneg hX _) hs.le
  have h2 : Y ^ ((P + Q) / Q) / ((P + Q) / Q) ≤ Y ^ ((P + Q) / Q) :=
    div_le_self (Real.rpow_nonneg hY _) ht.le
  linarith

lemma rpow_rpow_div {x : ℝ} (hx : 0 ≤ x) {P j R : ℝ} (hP : P ≠ 0) (hj : j ≠ 0) :
    (x ^ (P / j)) ^ (R / P) = x ^ (R / j) := by
  rw [← Real.rpow_mul hx]
  congr 1
  field_simp

lemma card_Icc_le (n ℓ : ℕ) : (((Finset.Icc 1 (min n ℓ)).card : ℕ) : ℝ) ≤ n := by
  rw [Nat.card_Icc]
  have : min n ℓ + 1 - 1 ≤ n := by omega
  exact_mod_cast this

lemma sum_rpow_le {I S : Finset ℕ} (hIS : I ⊆ S) {a : ℕ → ℝ} (ha : ∀ j, 0 ≤ a j)
    {p q : ℕ → ℝ} (hp : ∀ j ∈ I, 0 ≤ p j) (hpq : ∀ j ∈ I, p j ≤ q j) :
    ∑ j ∈ I, a j ^ p j ≤ I.card + ∑ j ∈ S, a j ^ q j := by
  calc ∑ j ∈ I, a j ^ p j ≤ ∑ j ∈ I, (1 + a j ^ q j) :=
        Finset.sum_le_sum fun j hj => rpow_le_one_add (ha j) (hp j hj) (hpq j hj)
    _ = I.card + ∑ j ∈ I, a j ^ q j := by rw [Finset.sum_add_distrib]; simp
    _ ≤ _ := add_le_add_right
        (Finset.sum_le_sum_of_subset_of_nonneg hIS fun j _ _ => Real.rpow_nonneg (ha j) _) _

theorem Om_mono {v v' ℓ ℓ' : ℕ} (hv : v ≤ v') (hℓ : ℓ ≤ ℓ') {a : ℕ → ℝ} (ha : ∀ j, 0 ≤ a j) :
    Om v ℓ a ≤ ((v : ℝ) + 1) * Om v' ℓ' a := by
  unfold Om
  have hsub : Finset.Icc 1 (min v ℓ) ⊆ Finset.Icc 1 (min v' ℓ') :=
    Finset.Icc_subset_Icc le_rfl (min_le_min hv hℓ)
  have h := sum_rpow_le hsub ha (p := fun j => (v : ℝ) / j) (q := fun j => (v' : ℝ) / j)
    (fun j _ => by positivity) (fun j _ => by
      have : (v : ℝ) ≤ v' := by exact_mod_cast hv
      gcongr)
  have hc := card_Icc_le v ℓ
  have hT : 0 ≤ ∑ j ∈ Finset.Icc 1 (min v' ℓ'), a j ^ ((v' : ℝ) / j) :=
    Finset.sum_nonneg fun j _ => Real.rpow_nonneg (ha j) _
  have := mul_nonneg (Nat.cast_nonneg v : (0 : ℝ) ≤ v) hT
  nlinarith

theorem exists_Om_mul (v v' : ℕ) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ (ℓ ℓ' : ℕ) (a : ℕ → ℝ), (∀ j, 0 ≤ a j) →
      Om v ℓ a * Om v' ℓ' a ≤ C * Om (v + v') (max ℓ ℓ') a := by
  refine ⟨2 + v + v', by positivity, fun ℓ ℓ' a ha => ?_⟩
  obtain ⟨W, hWdef⟩ : ∃ W : ℕ → ℝ, ∀ j, W j = a j ^ (((v + v' : ℕ) : ℝ) / j) :=
    ⟨_, fun _ => rfl⟩
  have hW : ∀ j, 0 ≤ W j := fun j => by rw [hWdef]; exact Real.rpow_nonneg (ha j) _
  have hI : Finset.Icc 1 (min v ℓ) ⊆ Finset.Icc 1 (min (v + v') (max ℓ ℓ')) :=
    Finset.Icc_subset_Icc le_rfl (by omega)
  have hI' : Finset.Icc 1 (min v' ℓ') ⊆ Finset.Icc 1 (min (v + v') (max ℓ ℓ')) :=
    Finset.Icc_subset_Icc le_rfl (by omega)
  have hT : ∑ j ∈ Finset.Icc 1 (min (v + v') (max ℓ ℓ')),
      a j ^ (((v + v' : ℕ) : ℝ) / j) =
      ∑ j ∈ Finset.Icc 1 (min (v + v') (max ℓ ℓ')), W j :=
    Finset.sum_congr rfl fun j _ => (hWdef j).symm
  have hv0 : (0 : ℝ) ≤ v := Nat.cast_nonneg v
  have hv0' : (0 : ℝ) ≤ v' := Nat.cast_nonneg v'
  have hTnn : 0 ≤ ∑ j ∈ Finset.Icc 1 (min (v + v') (max ℓ ℓ')), W j :=
    Finset.sum_nonneg fun j _ => hW j
  generalize hTdef : ∑ j ∈ Finset.Icc 1 (min (v + v') (max ℓ ℓ')), W j = T at hT hTnn
  -- single sums
  have hA : ∑ j ∈ Finset.Icc 1 (min v ℓ), a j ^ ((v : ℝ) / j) ≤ v + T := by
    have h := sum_rpow_le hI ha (p := fun j => (v : ℝ) / j)
      (q := fun j => (((v + v' : ℕ) : ℝ)) / j)
      (fun j _ => by positivity) (fun j _ => by
        have : (v : ℝ) ≤ ((v + v' : ℕ) : ℝ) := by push_cast; linarith
        gcongr)
    have hc := card_Icc_le v ℓ
    simp only [hT] at h
    linarith
  have hB : ∑ i ∈ Finset.Icc 1 (min v' ℓ'), a i ^ ((v' : ℝ) / i) ≤ v' + T := by
    have h := sum_rpow_le hI' ha (p := fun j => (v' : ℝ) / j)
      (q := fun j => (((v + v' : ℕ) : ℝ)) / j)
      (fun j _ => by positivity) (fun j _ => by
        have : (v' : ℝ) ≤ ((v + v' : ℕ) : ℝ) := by push_cast; linarith
        gcongr)
    have hc := card_Icc_le v' ℓ'
    simp only [hT] at h
    linarith
  -- cross terms
  have hX : ∀ j ∈ Finset.Icc 1 (min v ℓ), ∀ i ∈ Finset.Icc 1 (min v' ℓ'),
      a j ^ ((v : ℝ) / j) * a i ^ ((v' : ℝ) / i) ≤ W j + W i := by
    intro j hj i hi
    have hj' := Finset.mem_Icc.1 hj
    have hi' := Finset.mem_Icc.1 hi
    have hvpos : (0 : ℝ) < v := Nat.cast_pos.2 (by omega)
    have hvpos' : (0 : ℝ) < v' := Nat.cast_pos.2 (by omega)
    have hjne : (j : ℝ) ≠ 0 := Nat.cast_ne_zero.2 (by omega)
    have hine : (i : ℝ) ≠ 0 := Nat.cast_ne_zero.2 (by omega)
    have h := young_aux (Real.rpow_nonneg (ha j) ((v : ℝ) / j))
      (Real.rpow_nonneg (ha i) ((v' : ℝ) / i)) hvpos hvpos'
    rw [rpow_rpow_div (ha j) hvpos.ne' hjne, rpow_rpow_div (ha i) hvpos'.ne' hine] at h
    rw [hWdef, hWdef]
    push_cast
    exact h
  have hcr : ∑ j ∈ Finset.Icc 1 (min v ℓ), ∑ i ∈ Finset.Icc 1 (min v' ℓ'),
      a j ^ ((v : ℝ) / j) * a i ^ ((v' : ℝ) / i) ≤
      ∑ j ∈ Finset.Icc 1 (min v ℓ), ∑ i ∈ Finset.Icc 1 (min v' ℓ'), (W j + W i) :=
    Finset.sum_le_sum fun j hj => Finset.sum_le_sum fun i hi => hX j hj i hi
  have hexp : ∑ j ∈ Finset.Icc 1 (min v ℓ), ∑ i ∈ Finset.Icc 1 (min v' ℓ'), (W j + W i) =
      ((Finset.Icc 1 (min v' ℓ')).card : ℝ) * ∑ j ∈ Finset.Icc 1 (min v ℓ), W j +
      ((Finset.Icc 1 (min v ℓ)).card : ℝ) * ∑ i ∈ Finset.Icc 1 (min v' ℓ'), W i := by
    simp [Finset.sum_add_distrib, Finset.mul_sum]
  have hs1 : ∑ j ∈ Finset.Icc 1 (min v ℓ), W j ≤ T := by
    rw [← hTdef]; exact Finset.sum_le_sum_of_subset_of_nonneg hI fun j _ _ => hW j
  have hs2 : ∑ i ∈ Finset.Icc 1 (min v' ℓ'), W i ≤ T := by
    rw [← hTdef]; exact Finset.sum_le_sum_of_subset_of_nonneg hI' fun j _ _ => hW j
  have hs1' : 0 ≤ ∑ j ∈ Finset.Icc 1 (min v ℓ), W j := Finset.sum_nonneg fun j _ => hW j
  have hs2' : 0 ≤ ∑ j ∈ Finset.Icc 1 (min v' ℓ'), W j := Finset.sum_nonneg fun j _ => hW j
  have hm1 : ((Finset.Icc 1 (min v' ℓ')).card : ℝ) * ∑ j ∈ Finset.Icc 1 (min v ℓ), W j ≤ v' * T :=
    mul_le_mul (card_Icc_le v' ℓ') hs1 hs1' hv0'
  have hm2 : ((Finset.Icc 1 (min v ℓ)).card : ℝ) * ∑ i ∈ Finset.Icc 1 (min v' ℓ'), W i ≤ v * T :=
    mul_le_mul (card_Icc_le v ℓ) hs2 hs2' hv0
  unfold Om
  rw [hT]
  rw [add_mul, mul_add, mul_add, Finset.sum_mul_sum]
  nlinarith [mul_nonneg hv0 hTnn, mul_nonneg hv0' hTnn]

theorem exists_mul_Om (v : ℕ) {j : ℕ} (hj : 1 ≤ j) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ (ℓ : ℕ) (a : ℕ → ℝ), (∀ i, 0 ≤ a i) →
      a j * Om v ℓ a ≤ C * Om (v + j) (max ℓ j) a := by
  refine ⟨(v : ℝ) + 2, by positivity, fun ℓ a ha => ?_⟩
  obtain ⟨W, hWdef⟩ : ∃ W : ℕ → ℝ, ∀ i, W i = a i ^ (((v + j : ℕ) : ℝ) / i) :=
    ⟨_, fun _ => rfl⟩
  have hW : ∀ i, 0 ≤ W i := fun i => by rw [hWdef]; exact Real.rpow_nonneg (ha i) _
  have hI : Finset.Icc 1 (min v ℓ) ⊆ Finset.Icc 1 (min (v + j) (max ℓ j)) :=
    Finset.Icc_subset_Icc le_rfl (by omega)
  have hjS : j ∈ Finset.Icc 1 (min (v + j) (max ℓ j)) := by
    rw [Finset.mem_Icc]; omega
  have hT : ∑ i ∈ Finset.Icc 1 (min (v + j) (max ℓ j)), a i ^ (((v + j : ℕ) : ℝ) / i) =
      ∑ i ∈ Finset.Icc 1 (min (v + j) (max ℓ j)), W i :=
    Finset.sum_congr rfl fun i _ => (hWdef i).symm
  have hTnn : 0 ≤ ∑ i ∈ Finset.Icc 1 (min (v + j) (max ℓ j)), W i :=
    Finset.sum_nonneg fun i _ => hW i
  have hWj : W j ≤ ∑ i ∈ Finset.Icc 1 (min (v + j) (max ℓ j)), W i :=
    Finset.single_le_sum (f := W) (fun i _ => hW i) hjS
  have hv0 : (0 : ℝ) ≤ v := Nat.cast_nonneg v
  have hj0 : (0 : ℝ) < j := Nat.cast_pos.2 (by omega)
  have haj : a j ≤ 1 + W j := by
    have h := rpow_le_one_add (ha j) (p := 1) (p' := ((v + j : ℕ) : ℝ) / j) zero_le_one (by
      rw [le_div_iff₀ hj0]; push_cast; linarith)
    rw [Real.rpow_one] at h
    rw [hWdef]; exact h
  have hterm : ∀ i ∈ Finset.Icc 1 (min v ℓ),
      a j * a i ^ ((v : ℝ) / i) ≤ W j + W i := by
    intro i hi
    have hi' := Finset.mem_Icc.1 hi
    have hvpos : (0 : ℝ) < v := Nat.cast_pos.2 (by omega)
    have hine : (i : ℝ) ≠ 0 := Nat.cast_ne_zero.2 (by omega)
    have h := young_aux (Real.rpow_nonneg (ha i) ((v : ℝ) / i)) (ha j) hvpos hj0
    rw [rpow_rpow_div (ha i) hvpos.ne' hine] at h
    rw [hWdef, hWdef]
    push_cast
    linarith
  have hsum : ∑ i ∈ Finset.Icc 1 (min v ℓ), a j * a i ^ ((v : ℝ) / i) ≤
      ∑ i ∈ Finset.Icc 1 (min v ℓ), (W j + W i) := Finset.sum_le_sum hterm
  have hexp : ∑ i ∈ Finset.Icc 1 (min v ℓ), (W j + W i) =
      ((Finset.Icc 1 (min v ℓ)).card : ℝ) * W j + ∑ i ∈ Finset.Icc 1 (min v ℓ), W i := by
    simp [Finset.sum_add_distrib]
  have hs : ∑ i ∈ Finset.Icc 1 (min v ℓ), W i ≤
      ∑ i ∈ Finset.Icc 1 (min (v + j) (max ℓ j)), W i :=
    Finset.sum_le_sum_of_subset_of_nonneg hI fun i _ _ => hW i
  have hm : ((Finset.Icc 1 (min v ℓ)).card : ℝ) * W j ≤ v * W j :=
    mul_le_mul_of_nonneg_right (card_Icc_le v ℓ) (hW j)
  unfold Om
  rw [← Finset.mul_sum] at hsum
  rw [hT, mul_add, mul_one]
  nlinarith [mul_nonneg hv0 hTnn, mul_le_mul_of_nonneg_left hWj hv0]

/-- `x ^ p ≤ ε x ^ p' + C` for `0 < p < p'`. -/
lemma eps_aux {p p' ε : ℝ} (hp : 0 < p) (hpp : p < p') (hε : 0 < ε) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ x : ℝ, 0 ≤ x → x ^ p ≤ ε * x ^ p' + C := by
  have hd : p' - p ≠ 0 := by linarith
  set T : ℝ := ε ^ (-1 / (p' - p)) with hTdef
  have hT : 0 < T := Real.rpow_pos_of_pos hε _
  have hTd : T ^ (p - p') = ε := by
    rw [hTdef, ← Real.rpow_mul hε.le]
    have : (-1 / (p' - p)) * (p - p') = 1 := by
      calc (-1 / (p' - p)) * (p - p') = (p' - p) / (p' - p) := by ring
        _ = 1 := div_self hd
    rw [this, Real.rpow_one]
  refine ⟨T ^ p, Real.rpow_nonneg hT.le p, fun x hx => ?_⟩
  have hxp' : 0 ≤ x ^ p' := Real.rpow_nonneg hx _
  by_cases hxT : x ≤ T
  · have := Real.rpow_le_rpow hx hxT hp.le
    nlinarith [mul_nonneg hε.le hxp']
  · replace hxT := not_le.1 hxT
    have hxpos : 0 < x := hT.trans hxT
    have e : x ^ p = x ^ p' * x ^ (p - p') := by
      rw [← Real.rpow_add hxpos]; congr 1; ring
    have h1 : x ^ (p - p') ≤ T ^ (p - p') :=
      Real.rpow_le_rpow_of_nonpos hT hxT.le (by linarith)
    rw [hTd] at h1
    have := mul_le_mul_of_nonneg_left h1 hxp'
    have h3 : 0 ≤ T ^ p := Real.rpow_nonneg hT.le p
    nlinarith

theorem exists_Om_le_eps {m : ℕ} (hm : 1 ≤ m) {ε : ℝ} (hε : 0 < ε) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ a : ℕ → ℝ, (∀ j, 0 ≤ a j) →
      Om (2 * m - 1) m a ≤ ε * ∑ j ∈ Finset.Icc 1 m, a j ^ ((2 * m : ℝ) / j) + C := by
  have hm' : (1 : ℝ) ≤ m := by exact_mod_cast hm
  have h : ∀ j : ℕ, ∃ C : ℝ, 0 ≤ C ∧ (1 ≤ j → ∀ x : ℝ, 0 ≤ x →
      x ^ ((2 * (m : ℝ) - 1) / j) ≤ ε * x ^ ((2 * (m : ℝ)) / j) + C) := by
    intro j
    by_cases hj : 1 ≤ j
    · have hj0 : (0 : ℝ) < j := Nat.cast_pos.2 (by omega)
      obtain ⟨C, hC, hC'⟩ := eps_aux (p := (2 * (m : ℝ) - 1) / j) (p' := (2 * (m : ℝ)) / j)
        (ε := ε) (div_pos (by linarith) hj0)
        (by rw [div_lt_div_iff_of_pos_right hj0]; linarith) hε
      exact ⟨C, hC, fun _ => hC'⟩
    · exact ⟨0, le_rfl, fun h => absurd h hj⟩
  choose Cf hC0 hCf using h
  refine ⟨1 + ∑ j ∈ Finset.Icc 1 m, Cf j,
    add_nonneg zero_le_one (Finset.sum_nonneg fun j _ => hC0 j), fun a ha => ?_⟩
  have hmin : min (2 * m - 1) m = m := by omega
  have hcast : ((2 * m - 1 : ℕ) : ℝ) = 2 * (m : ℝ) - 1 := by
    rw [Nat.cast_sub (by omega)]; push_cast; ring
  unfold Om
  rw [hmin, hcast]
  have hs : ∑ j ∈ Finset.Icc 1 m, a j ^ ((2 * (m : ℝ) - 1) / j) ≤
      ∑ j ∈ Finset.Icc 1 m, (ε * a j ^ ((2 * (m : ℝ)) / j) + Cf j) :=
    Finset.sum_le_sum fun j hj => hCf j (Finset.mem_Icc.1 hj).1 (a j) (ha j)
  rw [Finset.sum_add_distrib, ← Finset.mul_sum] at hs
  linarith

theorem Om_two_mul (m : ℕ) (a : ℕ → ℝ) :
    Om (2 * m) m a = 1 + ∑ j ∈ Finset.Icc 1 m, a j ^ ((2 * m : ℝ) / j) := by
  have hmin : min (2 * m) m = m := by omega
  unfold Om
  rw [hmin]
  push_cast
  rfl

/-- info: 'CWFlow.Extension.exists_Om_mul' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_Om_mul

/-- info: 'CWFlow.Extension.exists_mul_Om' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_mul_Om

/-- info: 'CWFlow.Extension.exists_Om_le_eps' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_Om_le_eps

end CWFlow.Extension
