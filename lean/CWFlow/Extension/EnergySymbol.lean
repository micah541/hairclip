import CWFlow.Extension.CovSymbols
import CWFlow.Extension.GNorm

/-!
# The energy density as a symbol

For a Lagrangian graph with `u = D²φ`, the energy density `|∇^k C|²_g √det g` is, on jet space,
`eSym k = gNormSq (ginvFun ζ₂) (tSym k) * volFun ζ₂`.  We prove the lower bound
`(1 + ‖u‖²)⁻¹ |v|² ≤ vᵀ g⁻¹ v` for the inverse metric, that `qSym m = |·|²_{g⁻¹}` and `eSym k`
are symbols (of weights `(0,0)` and `(k+k, k)`), and the two-sided comparison of `eSym k` with
`‖tSym k‖² √det g`.
-/

open scoped ContDiff Matrix

namespace CWFlow.Extension

variable {n N : ℕ}

/-! ### Part 1: the inverse metric as a form -/

theorem ginvFun_apply (u : Z n 2) (a b : Fin n) : ginvFun u ![a, b] = (gMat u)⁻¹ a b := by
  simp [ginvFun]

theorem gMat_symm (u : Z n 2) (a b : Fin n) : gMat u a b = gMat u b a := by
  unfold gMat
  have : (if a = b then (1 : ℝ) else 0) = if b = a then 1 else 0 := by
    by_cases h : a = b <;> simp [h, eq_comm]
  rw [this]
  congr 1
  exact Finset.sum_congr rfl fun c _ => mul_comm _ _

theorem ginvFun_symm (u : Z n 2) (a b : Fin n) : ginvFun u ![a, b] = ginvFun u ![b, a] := by
  rw [ginvFun_apply, ginvFun_apply]
  have hT : (gMat u)ᵀ = gMat u := by
    ext i j; exact gMat_symm u j i
  have := Matrix.transpose_nonsing_inv (gMat u)
  rw [hT] at this
  have h2 := congrFun (congrFun this a) b
  simpa [Matrix.transpose_apply] using h2.symm

theorem sum_sq_eq (u : Z n 2) : ‖u‖ ^ 2 = ∑ a, ∑ c, u ![a, c] ^ 2 := by
  rw [normSq_eq_sum]
  rw [← (finTwoArrowEquiv (Fin n)).symm.sum_comp, Fintype.sum_prod_type]
  rfl

theorem ginvFun_quadratic_ge (u : Z n 2) (v : Fin n → ℝ) :
    (1 + ‖u‖ ^ 2)⁻¹ * ∑ a, v a ^ 2 ≤ ∑ a, ∑ b, (gMat u)⁻¹ a b * v a * v b := by
  set s : ℝ := ‖u‖ ^ 2 with hs
  have hs0 : 0 ≤ s := by positivity
  have hsum : s = ∑ c, ∑ a, uMat u a c ^ 2 := by
    rw [hs, sum_sq_eq, Finset.sum_comm]; rfl
  have hsum' : s = ∑ a, ∑ c, uMat u a c ^ 2 := by
    rw [hs, sum_sq_eq]; rfl
  set w := (gMat u)⁻¹ *ᵥ v with hw
  set y := w ᵥ* uMat u with hy
  have hv : gMat u *ᵥ w = v := by
    rw [hw, Matrix.mulVec_mulVec, gMat_mul_inv, Matrix.one_mulVec]
  have hv2 : v = w + uMat u *ᵥ y := by
    rw [← hv, gMat_eq, Matrix.add_mulVec, Matrix.one_mulVec, ← Matrix.mulVec_mulVec,
      Matrix.mulVec_transpose]
  have hcross : w ⬝ᵥ (uMat u *ᵥ y) = y ⬝ᵥ y := by
    rw [Matrix.dotProduct_mulVec]
  have hvv : v ⬝ᵥ v = w ⬝ᵥ w + 2 * (y ⬝ᵥ y) + (uMat u *ᵥ y) ⬝ᵥ (uMat u *ᵥ y) := by
    have : v ⬝ᵥ v = (w + uMat u *ᵥ y) ⬝ᵥ (w + uMat u *ᵥ y) := by rw [← hv2]
    rw [this, add_dotProduct, dotProduct_add, dotProduct_add, hcross,
      dotProduct_comm (uMat u *ᵥ y) w, hcross]
    ring
  -- Cauchy–Schwarz 1: |y|² ≤ s |w|²
  have cs1 : y ⬝ᵥ y ≤ s * (w ⬝ᵥ w) := by
    have : y ⬝ᵥ y = ∑ c, (∑ a, w a * uMat u a c) ^ 2 := by
      simp [dotProduct, hy, Matrix.vecMul, sq]
    rw [this, hsum, Finset.sum_mul]
    refine Finset.sum_le_sum fun c _ => ?_
    have := Finset.sum_mul_sq_le_sq_mul_sq Finset.univ w (fun a => uMat u a c)
    have hww : w ⬝ᵥ w = ∑ a, w a ^ 2 := by simp [dotProduct, sq]
    rw [hww]; linarith
  have cs2 : (uMat u *ᵥ y) ⬝ᵥ (uMat u *ᵥ y) ≤ s * (y ⬝ᵥ y) := by
    have : (uMat u *ᵥ y) ⬝ᵥ (uMat u *ᵥ y) = ∑ a, (∑ c, uMat u a c * y c) ^ 2 := by
      simp [dotProduct, Matrix.mulVec, sq]
    rw [this, hsum', Finset.sum_mul]
    refine Finset.sum_le_sum fun a _ => ?_
    have := Finset.sum_mul_sq_le_sq_mul_sq Finset.univ (fun c => uMat u a c) y
    have hyy : y ⬝ᵥ y = ∑ c, y c ^ 2 := by simp [dotProduct, sq]
    rw [hyy]; linarith
  have hq : v ⬝ᵥ w = w ⬝ᵥ w + y ⬝ᵥ y := by
    have := quad_g u w
    rw [hv] at this
    rw [dotProduct_comm]
    exact this
  have hfin : v ⬝ᵥ v ≤ (1 + s) * (v ⬝ᵥ w) := by
    rw [hq, hvv]
    nlinarith [cs1, cs2]
  rw [quad_inv_eq, inv_mul_le_iff₀ (by positivity)]
  have hV : ∑ a, v a ^ 2 = v ⬝ᵥ v := by simp [dotProduct, sq]
  rw [hV]
  exact hfin

theorem one_le_volFun (u : Z n 2) : 1 ≤ volFun u := by
  unfold volFun
  exact Real.one_le_sqrt.2 (one_le_det_gMat u)

/-! ### Part 2: the density and its weight -/

/-- The form `|·|²_{g⁻¹}` on `Z n m`, on jet space. -/
noncomputable def qSym (m : ℕ) (p : P n N) : Z n m →L[ℝ] Z n m →L[ℝ] ℝ :=
  gForm (ginvFun (low23 p).1)

/-- The energy density `|∇^k C|²_g √det g` on jet space. -/
noncomputable def eSym (k : ℕ) (p : P n N) : ℝ :=
  gNormSq (ginvFun (low23 p).1) (tSym k p) * volFun (low23 p).1

theorem isSymbol_qSym {A : Set (P n N)} (hA : IsOpen A)
    (hlow : ∀ m : Fin (N + 1), (m : ℕ) ≤ 3 → ∃ K : ℝ, ∀ p ∈ A, ‖p.2 m‖ ≤ K) (m : ℕ) :
    IsSymbol A 0 0 (qSym (N := N) m) :=
  isSymbol_comp_low hA low23 low23_high (bdd_low23 hlow)
    (c := fun q : Z n 2 × Z n 3 => (gForm (ginvFun q.1) : Z n m →L[ℝ] Z n m →L[ℝ] ℝ))
    (contDiff_gForm.comp (contDiff_ginvFun.comp contDiff_fst))

theorem eSym_eq (k : ℕ) (p : P n N) :
    eSym k p = qSym (k + 3) p (tSym k p) (tSym k p) * volFun (low23 p).1 := by
  unfold eSym qSym
  rw [gForm_apply, gNormSq_eq_gInner]

theorem isSymbol_eSym {A : Set (P n N)} (hA : IsOpen A)
    (hlow : ∀ m : Fin (N + 1), (m : ℕ) ≤ 3 → ∃ K : ℝ, ∀ p ∈ A, ‖p.2 m‖ ≤ K)
    (k : ℕ) (hN : k + 3 ≤ N) : IsSymbol A (k + k) k (eSym k) := by
  have hq := isSymbol_qSym hA hlow (n := n) (k + 3)
  have ht := isSymbol_tSym hA hlow k hN
  have h1 := IsSymbol.bilin
    (ContinuousLinearMap.apply ℝ (Z n (k + 3) →L[ℝ] ℝ)).flip hq ht
  have h2 := IsSymbol.bilin
    (ContinuousLinearMap.apply ℝ ℝ : Z n (k + 3) →L[ℝ] (Z n (k + 3) →L[ℝ] ℝ) →L[ℝ] ℝ).flip h1 ht
  have hv := isSymbol_vol hA (bdd_low23 hlow)
  have h3 := (h2.mul hv).mono (w' := k + k) (ℓ' := k) (by omega) (by simp)
  exact h3.congr fun p _ => by rw [eSym_eq]; rfl

/-! ### Part 3: comparison with the flat norm -/

theorem ginv_hsymm (u : Z n 2) : ∀ a b, ginvFun u ![a, b] = ginvFun u ![b, a] := ginvFun_symm u

theorem ginv_hnn (u : Z n 2) :
    ∀ v : Fin n → ℝ, 0 ≤ ∑ a, ∑ b, ginvFun u ![a, b] * v a * v b := by
  intro v
  simp only [ginvFun_apply]
  exact ginvFun_quadratic_nonneg u v

theorem eSym_nonneg (k : ℕ) (p : P n N) : 0 ≤ eSym k p :=
  mul_nonneg (gNormSq_nonneg (ginv_hsymm _) (ginv_hnn _) _)
    (zero_le_one.trans (one_le_volFun _))

theorem eSym_le (k : ℕ) (p : P n N) : eSym k p ≤ ‖tSym k p‖ ^ 2 * volFun (low23 p).1 := by
  refine mul_le_mul_of_nonneg_right (gNormSq_le (ginv_hsymm _) (ginv_hnn _) ?_ _)
    (zero_le_one.trans (one_le_volFun _))
  intro v
  simp only [ginvFun_apply]
  exact ginvFun_quadratic_le _ v

theorem le_eSym (k : ℕ) (p : P n N) :
    ((1 + ‖(low23 p).1‖ ^ 2)⁻¹) ^ (k + 3) * ‖tSym k p‖ ^ 2 ≤ eSym k p := by
  have h1 : ((1 + ‖(low23 p).1‖ ^ 2)⁻¹) ^ (k + 3) * ‖tSym k p‖ ^ 2 ≤
      gNormSq (ginvFun (low23 p).1) (tSym k p) := by
    refine le_gNormSq (by positivity) (ginv_hsymm _) (ginv_hnn _) ?_ _
    intro v
    simp only [ginvFun_apply]
    exact ginvFun_quadratic_ge _ v
  exact h1.trans (le_mul_of_one_le_right (gNormSq_nonneg (ginv_hsymm _) (ginv_hnn _) _)
    (one_le_volFun _))

/-- info: 'CWFlow.Extension.ginvFun_quadratic_ge' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms ginvFun_quadratic_ge

/-- info: 'CWFlow.Extension.isSymbol_eSym' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms isSymbol_eSym

/-- info: 'CWFlow.Extension.eSym_le' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms eSym_le

/-- info: 'CWFlow.Extension.le_eSym' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms le_eSym

end CWFlow.Extension
