import CWFlow.Extension.SymbolTotalD

/-!
# Index-level tensor operations on jet space

* `slotN m`: the `m`-jet coordinate with a natural-number index (`0` for `m > N`);
* `curryZ`: turn a `Z n k`-valued linear map on `ℝⁿ` into a `(k+1)`-tensor (new index first);
* `Dz 𝔭`: the total derivative of a tensor-valued function as a tensor with one more index, with
  its symbol closure and its realization along jets;
* `actG Γ T`: the Christoffel correction term of the covariant derivative of a `k`-tensor.
-/

open scoped ContDiff

namespace CWFlow.Extension

variable {n N : ℕ}

/-! ### Part 1: coordinates with a natural-number index -/

/-- The `m`-jet coordinate (`0` if `m > N`). -/
noncomputable def slotN (m : ℕ) : P n N →L[ℝ] Z n m :=
  if h : m < N + 1 then slotProjL (⟨m, h⟩ : Fin (N + 1)) else 0

theorem slotN_apply (m : ℕ) (h : m < N + 1) (p : P n N) : slotN m p = p.2 ⟨m, h⟩ := by
  unfold slotN
  rw [dif_pos h]
  rfl

theorem slotN_of_lt (m : ℕ) (h : N < m) (p : P n N) : slotN m p = 0 := by
  unfold slotN
  rw [dif_neg (by omega)]
  rfl

theorem slotN_jp (m : ℕ) (h : m ≤ N) (φ : EuclideanSpace ℝ (Fin n) → ℝ)
    (x : EuclideanSpace ℝ (Fin n)) :
    slotN m (jp N φ x) = CWFlow.Interp.jet m φ x := by
  rw [slotN_apply m (by omega)]
  rfl

theorem isSymbol_slotN {A : Set (P n N)} (hA : IsOpen A) (m : ℕ) (hm : 4 ≤ m) :
    IsSymbol A (m - 3) (m - 3) (fun p => slotN m p) := by
  by_cases h : m < N + 1
  · have := isSymbol_slotProj (n := n) hA (⟨m, h⟩ : Fin (N + 1)) hm
    refine this.congr fun p _ => ?_
    rw [slotN_apply m h]
    rfl
  · exact ((isSymbol_const hA (0 : Z n m)).mono (Nat.zero_le _) (Nat.zero_le _)).congr
      fun p _ => (slotN_of_lt m (by omega) p).symm ▸ rfl

theorem isSymbol_slotN_low {A : Set (P n N)} (hA : IsOpen A) (m : ℕ) (hm : m ≤ 3)
    (hbd : ∃ K : ℝ, ∀ p ∈ A, ‖slotN m p‖ ≤ K) : IsSymbol A 0 0 (fun p => slotN m p) := by
  by_cases h : m < N + 1
  · have := isSymbol_slotProj_low (n := n) hA (⟨m, h⟩ : Fin (N + 1)) hm (by
      obtain ⟨K, hK⟩ := hbd
      exact ⟨K, fun p hp => by simpa [slotN_apply m h] using hK p hp⟩)
    refine this.congr fun p _ => ?_
    rw [slotN_apply m h]
    rfl
  · exact ((isSymbol_const hA (0 : Z n m)).congr
      fun p _ => (slotN_of_lt m (by omega) p).symm ▸ rfl)

/-! ### Part 2: currying a derivative into the first index -/

/-- The linear map underlying `curryZ`. -/
noncomputable def curryZLin {k : ℕ} :
    (EuclideanSpace ℝ (Fin n) →L[ℝ] Z n k) →ₗ[ℝ] Z n (k + 1) where
  toFun L := WithLp.toLp 2 (fun I => L (EuclideanSpace.single (I 0) 1) (Fin.tail I))
  map_add' L L' := by
    ext I
    simp
  map_smul' c L := by
    ext I
    simp

/-- Turn a `Z n k`-valued linear map on `ℝⁿ` into a `(k+1)`-tensor: the new index comes first. -/
noncomputable def curryZ {k : ℕ} : (EuclideanSpace ℝ (Fin n) →L[ℝ] Z n k) →L[ℝ] Z n (k + 1) :=
  LinearMap.toContinuousLinearMap curryZLin

theorem curryZ_apply {k : ℕ} (L : EuclideanSpace ℝ (Fin n) →L[ℝ] Z n k) (I : Fin (k + 1) → Fin n) :
    curryZ L I = L (EuclideanSpace.single (I 0) 1) (Fin.tail I) := rfl

theorem curryZ_contractFlip {k : ℕ} (z : Z n (k + 1)) : curryZ (contractFlip n k z) = z := by
  ext I
  rw [curryZ_apply, contractFlip_apply, contract_apply]
  simp

theorem norm_curryZ_sq {k : ℕ} (L : EuclideanSpace ℝ (Fin n) →L[ℝ] Z n k) :
    ‖curryZ L‖ ^ 2 = ∑ a, ‖L (EuclideanSpace.single a 1)‖ ^ 2 := by
  rw [PiLp.norm_sq_eq_of_L2, CWFlow.Interp.sum_fin_cons]
  refine Finset.sum_congr rfl fun a _ => ?_
  rw [PiLp.norm_sq_eq_of_L2]
  refine Finset.sum_congr rfl fun I _ => ?_
  have : curryZ L (Fin.cons a I) = L (EuclideanSpace.single a 1) I := by
    rw [curryZ_apply]
    simp
  rw [this]

/-! ### Part 3: the total derivative as a tensor -/

/-- The total derivative of a tensor-valued function on jet space, as a tensor with one more
index. -/
noncomputable def Dz {k : ℕ} (𝔭 : P n N → Z n k) (p : P n N) : Z n (k + 1) :=
  curryZ (totalD 𝔭 p)

theorem IsSymbol.Dz {k : ℕ} {A : Set (P n N)} {w ℓ : ℕ} {𝔭 : P n N → Z n k}
    (h : IsSymbol A w ℓ 𝔭) (hN : w + 4 ≤ N)
    (hlow : ∀ m : Fin (N + 1), (m : ℕ) ≤ 3 → ∃ K : ℝ, ∀ p ∈ A, ‖p.2 m‖ ≤ K) :
    IsSymbol A (w + 1) (max ℓ w + 1) (CWFlow.Extension.Dz 𝔭) :=
  (h.totalD hN hlow).clm_comp curryZ

/-- The total derivative of the `m`-jet coordinate is the `(m+1)`-jet coordinate. -/
theorem Dz_slotN (m : ℕ) (hm : m + 1 ≤ N) (p : P n N) :
    Dz (fun q : P n N => slotN m q) p = slotN (m + 1) p := by
  have h1 : totalD (fun q : P n N => slotN m q) p = contractFlip n m (slotN (m + 1) p) := by
    ext e : 1
    rw [totalD_apply, contractFlip_apply, ContinuousLinearMap.fderiv, slotN_apply m (by omega),
      slotN_apply (m + 1) (by omega)]
    show shift e p.2 ⟨m, _⟩ = _
    unfold shift
    rw [dif_pos (by show m + 1 < N + 1; omega)]
  unfold Dz
  rw [h1, curryZ_contractFlip]

/-- Along the jets of a function, `Dz` is the partial derivative in the new (first) index. -/
theorem Dz_jp_apply {k : ℕ} {A : Set (P n N)} {w ℓ : ℕ} {𝔭 : P n N → Z n k}
    (h : IsSymbol A w ℓ 𝔭) (hN : w + 4 ≤ N) {φ : EuclideanSpace ℝ (Fin n) → ℝ}
    (hφ : ContDiff ℝ ∞ φ) {x : EuclideanSpace ℝ (Fin n)} (hx : jp N φ x ∈ A)
    (I : Fin (k + 1) → Fin n) :
    Dz 𝔭 (jp N φ x) I =
      fderiv ℝ (fun y => 𝔭 (jp N φ y) (Fin.tail I)) x (EuclideanSpace.single (I 0) 1) := by
  have hd := hasFDerivAt_symbol_jp h hN hφ hx
  rw [← CWFlow.Interp.fderiv_piLp_apply hd.differentiableAt, hd.fderiv]
  rfl

/-! ### Part 4: the action of Christoffel symbols -/

/-- The linear map `T ↦ Γ · T`. -/
noncomputable def actGLin {k : ℕ} (G : Z n 3) : Z n k →L[ℝ] Z n (k + 1) :=
  LinearMap.toContinuousLinearMap
    { toFun := fun T => WithLp.toLp 2 (fun I : Fin (k + 1) → Fin n =>
        ∑ i : Fin k, ∑ c : Fin n,
          G ![I 0, Fin.tail I i, c] * T (Function.update (Fin.tail I) i c))
      map_add' := fun T T' => by
        ext I
        simp [mul_add, Finset.sum_add_distrib]
      map_smul' := fun a T => by
        ext I
        simp [Finset.mul_sum, mul_left_comm] }

/-- `(Γ · T)(a, I) = Σ_i Σ_c Γ(a, I_i, c) T(I with the i-th index replaced by c)`:
the correction term of the covariant derivative of a `k`-tensor (`Γ` has index order `(a, b, c)`
for `Γ^c_{ab}`). -/
noncomputable def actG {k : ℕ} : Z n 3 →L[ℝ] Z n k →L[ℝ] Z n (k + 1) :=
  LinearMap.toContinuousLinearMap
    { toFun := fun G => actGLin G
      map_add' := fun G G' => by
        ext T I
        simp [actGLin, add_mul, Finset.sum_add_distrib]
      map_smul' := fun a G => by
        ext T I
        simp [actGLin, Finset.mul_sum, mul_assoc] }

theorem actG_apply {k : ℕ} (G : Z n 3) (T : Z n k) (I : Fin (k + 1) → Fin n) :
    actG G T I = ∑ i : Fin k, ∑ c : Fin n,
      G ![I 0, Fin.tail I i, c] * T (Function.update (Fin.tail I) i c) := rfl

/-- info: 'CWFlow.Extension.isSymbol_slotN' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms isSymbol_slotN

/-- info: 'CWFlow.Extension.IsSymbol.Dz' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms IsSymbol.Dz

/-- info: 'CWFlow.Extension.Dz_slotN' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms Dz_slotN

/-- info: 'CWFlow.Extension.Dz_jp_apply' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms Dz_jp_apply

/-- info: 'CWFlow.Extension.actG_apply' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms actG_apply

end CWFlow.Extension
