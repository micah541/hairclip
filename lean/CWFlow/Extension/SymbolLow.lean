import CWFlow.Extension.CovSymbols

/-!
# Algebra of the total derivative, and functions of the second jet

Linearity and the product rule for `totalD`, and the sharp weight statement: the total derivative of
a globally smooth function of the 2-jet slot is a symbol of weight `0` (it involves the 3-jet slot,
which has weight `0`).
-/

open scoped ContDiff

namespace CWFlow.Extension

variable {n N : ℕ}

section Algebra

variable {W : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W]
  {𝔭 𝔮 : P n N → W} {p : P n N}

theorem totalD_add (h1 : DifferentiableAt ℝ 𝔭 p) (h2 : DifferentiableAt ℝ 𝔮 p) :
    totalD (fun q => 𝔭 q + 𝔮 q) p = totalD 𝔭 p + totalD 𝔮 p := by
  ext e
  simp only [totalD_apply, ContinuousLinearMap.add_apply]
  rw [fderiv_fun_add h1 h2]
  rfl

theorem totalD_sub (h1 : DifferentiableAt ℝ 𝔭 p) (h2 : DifferentiableAt ℝ 𝔮 p) :
    totalD (fun q => 𝔭 q - 𝔮 q) p = totalD 𝔭 p - totalD 𝔮 p := by
  ext e
  simp only [totalD_apply, ContinuousLinearMap.sub_apply]
  rw [fderiv_fun_sub h1 h2]
  rfl

theorem totalD_neg (p : P n N) : totalD (fun q => -𝔭 q) p = -totalD 𝔭 p := by
  ext e
  simp only [totalD_apply, ContinuousLinearMap.neg_apply]
  rw [fderiv_fun_neg]
  rfl

theorem totalD_const (c : W) (p : P n N) : totalD (fun _ : P n N => c) p = 0 := by
  ext e
  simp [totalD_apply]

theorem totalD_clm_comp {W' : Type*} [NormedAddCommGroup W'] [NormedSpace ℝ W']
    (L : W →L[ℝ] W') (h1 : DifferentiableAt ℝ 𝔭 p) :
    totalD (fun q => L (𝔭 q)) p = L.comp (totalD 𝔭 p) := by
  ext e
  simp only [totalD_apply, ContinuousLinearMap.comp_apply]
  have := (L.hasFDerivAt.comp p h1.hasFDerivAt).fderiv
  rw [show (fun q => L (𝔭 q)) = L ∘ 𝔭 from rfl, this]
  rfl

theorem totalD_sum {ι : Type*} (S : Finset ι) {𝔭 : ι → P n N → W}
    (h : ∀ i ∈ S, DifferentiableAt ℝ (𝔭 i) p) :
    totalD (fun q => ∑ i ∈ S, 𝔭 i q) p = ∑ i ∈ S, totalD (𝔭 i) p := by
  ext e
  simp only [totalD_apply, ContinuousLinearMap.sum_apply]
  rw [fderiv_fun_sum h]
  simp

end Algebra

section Product

theorem totalD_bilin {W₁ W₂ W₃ : Type*} [NormedAddCommGroup W₁] [NormedSpace ℝ W₁]
    [NormedAddCommGroup W₂] [NormedSpace ℝ W₂] [NormedAddCommGroup W₃] [NormedSpace ℝ W₃]
    (B : W₁ →L[ℝ] W₂ →L[ℝ] W₃) {𝔭 : P n N → W₁} {𝔮 : P n N → W₂} {p : P n N}
    (h1 : DifferentiableAt ℝ 𝔭 p) (h2 : DifferentiableAt ℝ 𝔮 p)
    (e : EuclideanSpace ℝ (Fin n)) :
    totalD (fun q => B (𝔭 q) (𝔮 q)) p e = B (totalD 𝔭 p e) (𝔮 p) + B (𝔭 p) (totalD 𝔮 p e) := by
  simp only [totalD_apply]
  rw [ContinuousLinearMap.fderiv_of_bilinear B h1 h2]
  simp [add_comm]

theorem totalD_mul {𝔭 𝔮 : P n N → ℝ} {p : P n N} (h1 : DifferentiableAt ℝ 𝔭 p)
    (h2 : DifferentiableAt ℝ 𝔮 p) (e : EuclideanSpace ℝ (Fin n)) :
    totalD (fun q => 𝔭 q * 𝔮 q) p e = totalD 𝔭 p e * 𝔮 p + 𝔭 p * totalD 𝔮 p e :=
  totalD_bilin (ContinuousLinearMap.mul ℝ ℝ) h1 h2 e

theorem totalD_smul {W : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W] {𝔭 : P n N → ℝ}
    {𝔮 : P n N → W} {p : P n N} (h1 : DifferentiableAt ℝ 𝔭 p)
    (h2 : DifferentiableAt ℝ 𝔮 p) (e : EuclideanSpace ℝ (Fin n)) :
    totalD (fun q => 𝔭 q • 𝔮 q) p e = totalD 𝔭 p e • 𝔮 p + 𝔭 p • totalD 𝔮 p e :=
  totalD_bilin (ContinuousLinearMap.lsmul ℝ ℝ : ℝ →L[ℝ] W →L[ℝ] W) h1 h2 e

end Product

section Low2

variable {W : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W]

theorem low23_fst_shift (p : P n N) (e : EuclideanSpace ℝ (Fin n)) :
    (low23 ((e, shift e p.2) : P n N)).1 = contract e (low23 p).2 := by
  by_cases h3 : 3 < N + 1
  · have h2 : 2 < N + 1 := by omega
    rw [low23_fst _ h2, low23_snd _ h3]
    simp only [shift]
    rw [dif_pos (by simp; omega)]
  · have hz : (low23 p).2 = 0 := by
      simp [low23, slotOrZero, h3]
    rw [hz, map_zero]
    by_cases h2 : 2 < N + 1
    · rw [low23_fst _ h2]
      simp only [shift]
      rw [dif_neg (by simp; omega)]
    · simp [low23, slotOrZero, h2]

theorem totalD_comp_low2 {c : Z n 2 → W} (hc : ContDiff ℝ ∞ c) (p : P n N)
    (e : EuclideanSpace ℝ (Fin n)) :
    totalD (fun q : P n N => c (low23 q).1) p e =
      fderiv ℝ c (low23 p).1 (contract e (low23 p).2) := by
  have hL : (fun q : P n N => c (low23 q).1) =
      c ∘ ⇑((ContinuousLinearMap.fst ℝ (Z n 2) (Z n 3)).comp low23) := rfl
  rw [totalD_apply, hL, fderiv_comp p
    (hc.differentiable (by simp) _) (ContinuousLinearMap.differentiableAt _),
    ContinuousLinearMap.fderiv]
  simp only [ContinuousLinearMap.comp_apply, ContinuousLinearMap.coe_fst']
  rw [low23_fst_shift]

theorem isSymbol_totalD_comp_low2 {A : Set (P n N)} (hA : IsOpen A)
    (hb : Bornology.IsBounded (low23 '' A)) {c : Z n 2 → W} (hc : ContDiff ℝ ∞ c) :
    IsSymbol A 0 0 (totalD (fun q : P n N => c (low23 q).1)) := by
  let ct : Z n 2 × Z n 3 → (EuclideanSpace ℝ (Fin n) →L[ℝ] W) := fun q =>
    (fderiv ℝ c q.1).comp (contractFlip n 2 q.2)
  have hct : ContDiff ℝ ∞ ct := by
    have h1 : ContDiff ℝ ∞ (fun q : Z n 2 × Z n 3 => fderiv ℝ c q.1) :=
      (hc.fderiv_right (m := ∞) le_rfl).comp contDiff_fst
    have h2 : ContDiff ℝ ∞ (fun q : Z n 2 × Z n 3 => contractFlip n 2 q.2) :=
      (contractFlip n 2).contDiff.comp contDiff_snd
    exact (ContinuousLinearMap.compL ℝ (EuclideanSpace ℝ (Fin n)) (Z n 2) W).isBoundedBilinearMap.contDiff.comp
      (h1.prodMk h2)
  have := isSymbol_comp_low hA low23 low23_high hb hct
  refine this.congr fun q _ => ?_
  ext e
  rw [totalD_comp_low2 hc]
  rfl

end Low2

section Coeff

variable {A : Set (P n N)}

theorem isSymbol_totalD_ginv (hA : IsOpen A)
    (hlow : ∀ m : Fin (N + 1), (m : ℕ) ≤ 3 → ∃ K : ℝ, ∀ p ∈ A, ‖p.2 m‖ ≤ K) :
    IsSymbol A 0 0 (totalD (fun p : P n N => ginvFun (low23 p).1)) :=
  isSymbol_totalD_comp_low2 hA (bdd_low23 hlow) contDiff_ginvFun

theorem isSymbol_totalD_g (hA : IsOpen A)
    (hlow : ∀ m : Fin (N + 1), (m : ℕ) ≤ 3 → ∃ K : ℝ, ∀ p ∈ A, ‖p.2 m‖ ≤ K) :
    IsSymbol A 0 0 (totalD (fun p : P n N => gFun (low23 p).1)) :=
  isSymbol_totalD_comp_low2 hA (bdd_low23 hlow) contDiff_gFun

theorem isSymbol_totalD_vol (hA : IsOpen A)
    (hlow : ∀ m : Fin (N + 1), (m : ℕ) ≤ 3 → ∃ K : ℝ, ∀ p ∈ A, ‖p.2 m‖ ≤ K) :
    IsSymbol A 0 0 (totalD (fun p : P n N => volFun (low23 p).1)) :=
  isSymbol_totalD_comp_low2 hA (bdd_low23 hlow) contDiff_volFun

end Coeff

/-- info: 'CWFlow.Extension.totalD_bilin' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms totalD_bilin

/-- info: 'CWFlow.Extension.totalD_comp_low2' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms totalD_comp_low2

/-- info: 'CWFlow.Extension.isSymbol_totalD_comp_low2' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms isSymbol_totalD_comp_low2

end CWFlow.Extension
