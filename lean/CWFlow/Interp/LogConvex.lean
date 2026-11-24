import CWFlow.Interp.RealLemmas

/-!
# Multiplicative discrete convexity

`discrete_logconvex_le`: a nonnegative sequence with `N j ^ 2 ≤ c * N (j-1) * N (j+1)` on the
interior of `{0, …, m}` is bounded by the geometric interpolation of its end values times
`c ^ (j (m - j) / 2)`.  This is the multiplicative form of `discrete_convex_le`.
-/

namespace CWFlow.Interp

theorem discrete_logconvex_le {m : ℕ} (hm : 0 < m) {N : ℕ → ℝ} {c : ℝ} (hc : 0 < c)
    (hN : ∀ j, 0 ≤ N j)
    (h : ∀ j, 0 < j → j < m → N j ^ 2 ≤ c * N (j - 1) * N (j + 1)) :
    ∀ j, j ≤ m →
      N j ≤ c ^ ((j : ℝ) * (m - j) / 2) * N 0 ^ (1 - (j : ℝ) / m) * N m ^ ((j : ℝ) / m) := by
  by_cases hpos : ∀ j, j ≤ m → 0 < N j
  · intro j hj
    have h' : ∀ j, 0 < j → j < m →
        2 * Real.log (N j) ≤ Real.log (N (j - 1)) + Real.log (N (j + 1)) + Real.log c := by
      intro j h0 h1
      have p0 := hpos j h1.le
      have p1 := hpos (j - 1) (by omega)
      have p2 := hpos (j + 1) (by omega)
      have e : 2 * Real.log (N j) = Real.log (N j ^ 2) := by
        rw [Real.log_pow]; norm_num
      rw [e]
      calc Real.log (N j ^ 2) ≤ Real.log (c * N (j - 1) * N (j + 1)) :=
            Real.log_le_log (by positivity) (h j h0 h1)
        _ = _ := by
            rw [Real.log_mul (by positivity) p2.ne', Real.log_mul hc.ne' p1.ne']; ring
    have key := discrete_convex_le hm h' j hj
    have hm' : (0 : ℝ) < m := by exact_mod_cast hm
    have p0 := hpos 0 (Nat.zero_le _)
    have pm := hpos m le_rfl
    calc N j = Real.exp (Real.log (N j)) := (Real.exp_log (hpos j hj)).symm
      _ ≤ Real.exp ((1 - (j : ℝ) / m) * Real.log (N 0) + ((j : ℝ) / m) * Real.log (N m)
            + Real.log c * j * (m - j) / 2) := Real.exp_le_exp.mpr key
      _ = _ := by
          rw [Real.rpow_def_of_pos hc, Real.rpow_def_of_pos p0, Real.rpow_def_of_pos pm,
            ← Real.exp_add, ← Real.exp_add]
          congr 1; ring
  · push_neg at hpos
    obtain ⟨j0, hj0, hz⟩ := hpos
    have hz0 : N j0 = 0 := le_antisymm hz (hN _)
    have step : ∀ j, 0 < j → j < m → (N (j - 1) = 0 ∨ N (j + 1) = 0) → N j = 0 := by
      intro j h0 h1 hh
      have := h j h0 h1
      have : N j ^ 2 ≤ 0 := by
        rcases hh with hh | hh <;> simpa [hh] using this
      exact pow_eq_zero_iff (two_ne_zero) |>.mp (le_antisymm this (sq_nonneg _))
    have up : ∀ d, j0 + d < m → N (j0 + d) = 0 := by
      intro d
      induction d with
      | zero => intro _; simpa using hz0
      | succ d ih =>
        intro hd
        apply step (j0 + (d + 1)) (by omega) hd
        left
        have := ih (by omega)
        have e : j0 + (d + 1) - 1 = j0 + d := by omega
        rw [e]; exact this
    have down : ∀ d, 0 < j0 - d → N (j0 - d) = 0 := by
      intro d
      induction d with
      | zero => intro _; simpa using hz0
      | succ d ih =>
        intro hd
        apply step (j0 - (d + 1)) hd (by omega)
        right
        have := ih (by omega)
        have e : j0 - (d + 1) + 1 = j0 - d := by omega
        rw [e]; exact this
    have interior : ∀ j, 0 < j → j < m → N j = 0 := by
      intro j h0 h1
      by_cases hj : j0 ≤ j
      · have := up (j - j0) (by omega)
        have e : j0 + (j - j0) = j := by omega
        rwa [e] at this
      · have := down (j0 - j) (by omega)
        have e : j0 - (j0 - j) = j := by omega
        rwa [e] at this
    have hm' : (0 : ℝ) < m := by exact_mod_cast hm
    intro j hj
    rcases Nat.eq_zero_or_pos j with rfl | h0
    · simp
    rcases hj.lt_or_eq with h1 | rfl
    · rw [interior j h0 h1]
      have := hN 0
      have := hN m
      positivity
    · have e : (j : ℝ) / j = 1 := div_self hm'.ne'
      simp [e]

/-- info: 'CWFlow.Interp.discrete_logconvex_le' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms discrete_logconvex_le

end CWFlow.Interp
