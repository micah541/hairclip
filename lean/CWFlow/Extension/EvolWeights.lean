import CWFlow.Extension.EnergySymbol
import CWFlow.Extension.SymbolLow

/-!
# The weight tensor `√det g · ♯(∇^k C)` as a symbol

`sharp h` raises all indices of an `m`-tensor with the form `h` (the inverse metric).  The weight
tensor `wSym k = √det g · ♯(∇^k C)` on jet space satisfies
`|∇^k C|²-pairing · √det g = Σ_I S_I (wSym k)_I`, and `wSym k`, `Dz (wSym k)`, `Dz (Dz (wSym k))`
are symbols of weights `(k,k)`, `(k+1,k+1)`, `(k+2,k+2)`.
-/

open scoped ContDiff

namespace CWFlow.Extension

variable {n N : ℕ}

/-! ### Part 1: raising all indices -/

/-- The linear map underlying `sharp`. -/
noncomputable def sharpLin {m : ℕ} (h : Z n 2) : Z n m →ₗ[ℝ] Z n m where
  toFun T := WithLp.toLp 2 (fun I : Fin m → Fin n => ∑ J : Fin m → Fin n, (∏ i, h ![I i, J i]) * T J)
  map_add' S T := by
    ext I
    simp [mul_add, Finset.sum_add_distrib]
  map_smul' c T := by
    ext I
    simp [Finset.mul_sum, mul_left_comm]

/-- `(♯T)_I = Σ_J (Π_i h(I_i, J_i)) T_J`. -/
noncomputable def sharp {m : ℕ} (h : Z n 2) : Z n m →L[ℝ] Z n m :=
  LinearMap.toContinuousLinearMap (sharpLin h)

theorem sharp_apply {m : ℕ} (h : Z n 2) (T : Z n m) (I : Fin m → Fin n) :
    sharp h T I = ∑ J : Fin m → Fin n, (∏ i, h ![I i, J i]) * T J := rfl

theorem gInner_eq_sum_sharp {m : ℕ} (h : Z n 2) (S T : Z n m) :
    gInner h S T = ∑ I : Fin m → Fin n, S I * sharp h T I := by
  unfold gInner
  refine Finset.sum_congr rfl fun I _ => ?_
  rw [sharp_apply, Finset.mul_sum]
  refine Finset.sum_congr rfl fun J _ => ?_
  ring

theorem ew_contDiff_sharp_apply {m : ℕ} (T : Z n m) :
    ContDiff ℝ ∞ (fun h : Z n 2 => sharp h T) := by
  refine contDiff_piLp' (𝕜 := ℝ) 2 fun I => ?_
  simp only [sharp_apply]
  refine ContDiff.sum fun J _ => ?_
  exact (contDiff_prod fun i _ => contDiff_coordG ![I i, J i]).mul contDiff_const

theorem contDiff_sharp {m : ℕ} : ContDiff ℝ ∞ (fun h : Z n 2 => (sharp h : Z n m →L[ℝ] Z n m)) := by
  rw [contDiff_clm_apply_iff]
  exact fun T => ew_contDiff_sharp_apply T

/-! ### Part 2: the weight tensor -/

/-- The weight tensor `√det g · ♯(∇^k C)` on jet space. -/
noncomputable def wSym (k : ℕ) (p : P n N) : Z n (k + 3) :=
  volFun (low23 p).1 • sharp (ginvFun (low23 p).1) (tSym k p)

/-- The density pairing in terms of the weight tensor. -/
theorem qSym_mul_vol_eq (k : ℕ) (p : P n N) (S : Z n (k + 3)) :
    qSym (k + 3) p S (tSym k p) * volFun (low23 p).1 = ∑ I : Fin (k + 3) → Fin n, S I * wSym k p I := by
  unfold qSym
  rw [gForm_apply, gInner_eq_sum_sharp, Finset.sum_mul]
  refine Finset.sum_congr rfl fun I _ => ?_
  simp only [wSym, PiLp.smul_apply, smul_eq_mul]
  ring

/-- The coefficient map `u ↦ √det g(u) · ♯_{g(u)⁻¹}`, a globally smooth function of the second
jet. -/
noncomputable def wCoef (m : ℕ) (u : Z n 2) : Z n m →L[ℝ] Z n m := volFun u • sharp (ginvFun u)

theorem contDiff_wCoef (m : ℕ) : ContDiff ℝ ∞ (wCoef (n := n) m) :=
  contDiff_volFun.smul (contDiff_sharp.comp contDiff_ginvFun)

theorem wSym_eq (k : ℕ) (p : P n N) : wSym k p = wCoef (k + 3) (low23 p).1 (tSym k p) := rfl

/-! ### Part 3: weights -/

theorem isSymbol_wCoef {A : Set (P n N)} (hA : IsOpen A)
    (hlow : ∀ m : Fin (N + 1), (m : ℕ) ≤ 3 → ∃ K : ℝ, ∀ p ∈ A, ‖p.2 m‖ ≤ K) (m : ℕ) :
    IsSymbol A 0 0 (fun p : P n N => wCoef m (low23 p).1) :=
  isSymbol_comp_low hA low23 low23_high (bdd_low23 hlow)
    (c := fun q : Z n 2 × Z n 3 => wCoef m q.1) ((contDiff_wCoef m).comp contDiff_fst)

theorem isSymbol_wSym {A : Set (P n N)} (hA : IsOpen A)
    (hlow : ∀ m : Fin (N + 1), (m : ℕ) ≤ 3 → ∃ K : ℝ, ∀ p ∈ A, ‖p.2 m‖ ≤ K)
    (k : ℕ) (hN : k + 3 ≤ N) : IsSymbol A k k (wSym k) := by
  have hc := isSymbol_wCoef hA hlow (n := n) (k + 3)
  have ht := isSymbol_tSym hA hlow k hN
  have h1 := IsSymbol.bilin
    (ContinuousLinearMap.apply ℝ (Z n (k + 3)) : Z n (k + 3) →L[ℝ] (Z n (k + 3) →L[ℝ] Z n (k + 3)) →L[ℝ] Z n (k + 3)).flip
    hc ht
  have h2 := h1.mono (w' := k) (ℓ' := k) (by omega) (by simp)
  exact h2.congr fun p _ => rfl

theorem isSymbol_Dz_wSym {A : Set (P n N)} (hA : IsOpen A)
    (hlow : ∀ m : Fin (N + 1), (m : ℕ) ≤ 3 → ∃ K : ℝ, ∀ p ∈ A, ‖p.2 m‖ ≤ K)
    (k : ℕ) (hN : k + 4 ≤ N) : IsSymbol A (k + 1) (k + 1) (Dz (wSym (N := N) k)) := by
  have h := (isSymbol_wSym hA hlow k (by omega)).Dz (by omega) hlow
  exact h.mono (le_refl _) (by simp)

theorem isSymbol_Dz_Dz_wSym {A : Set (P n N)} (hA : IsOpen A)
    (hlow : ∀ m : Fin (N + 1), (m : ℕ) ≤ 3 → ∃ K : ℝ, ∀ p ∈ A, ‖p.2 m‖ ≤ K)
    (k : ℕ) (hN : k + 5 ≤ N) : IsSymbol A (k + 2) (k + 2) (Dz (Dz (wSym (N := N) k))) := by
  have h := (isSymbol_Dz_wSym hA hlow k (by omega)).Dz (by omega) hlow
  exact h.mono (le_refl _) (by simp)

/-- Coordinates of the weight tensor are real-valued symbols. -/
theorem isSymbol_wSym_apply {A : Set (P n N)} (hA : IsOpen A)
    (hlow : ∀ m : Fin (N + 1), (m : ℕ) ≤ 3 → ∃ K : ℝ, ∀ p ∈ A, ‖p.2 m‖ ≤ K)
    (k : ℕ) (hN : k + 3 ≤ N) (I : Fin (k + 3) → Fin n) :
    IsSymbol A k k (fun p : P n N => wSym k p I) :=
  (isSymbol_wSym hA hlow k hN).clm_comp (PiLp.proj 2 (fun _ : Fin (k + 3) → Fin n => ℝ) I)

/-- info: 'CWFlow.Extension.gInner_eq_sum_sharp' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms gInner_eq_sum_sharp

/-- info: 'CWFlow.Extension.qSym_mul_vol_eq' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms qSym_mul_vol_eq

/-- info: 'CWFlow.Extension.isSymbol_wSym' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms isSymbol_wSym

/-- info: 'CWFlow.Extension.isSymbol_Dz_wSym' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms isSymbol_Dz_wSym

/-- info: 'CWFlow.Extension.isSymbol_Dz_Dz_wSym' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms isSymbol_Dz_Dz_wSym

end CWFlow.Extension
