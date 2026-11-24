import CWFlow.Extension.Weights
import CWFlow.Interp.Jets
import Mathlib.Analysis.Calculus.ContDiff.Defs
import Mathlib.Analysis.Calculus.ContDiff.Operations
import Mathlib.Analysis.Calculus.IteratedDeriv.Defs

/-!
# Symbol classes

Smooth functions of the base point and the jets up to order `N`, whose derivatives obey weighted
growth bounds in terms of the weight function `Om`.
-/

open scoped ContDiff

namespace CWFlow.Extension

/-- The space of `m`-jets of a scalar function on `ℝⁿ`. -/
abbrev Z (n m : ℕ) : Type := PiLp 2 (fun _ : Fin m → Fin n => ℝ)

/-- Jets of orders `0, …, N`. -/
abbrev J (n N : ℕ) : Type := (m : Fin (N + 1)) → Z n m

/-- Base point and jets. -/
abbrev P (n N : ℕ) : Type := EuclideanSpace ℝ (Fin n) × J n N

variable {n N : ℕ}

/-- The jets of `φ` at `x`, as a point of `P n N`. -/
noncomputable def jp (N : ℕ) (φ : EuclideanSpace ℝ (Fin n) → ℝ) (x : EuclideanSpace ℝ (Fin n)) :
    P n N :=
  (x, fun m => CWFlow.Interp.jet (m : ℕ) φ x)

/-- Slots: `none` is the base point, `some m` is the `m`-jet. -/
abbrev Slot (N : ℕ) : Type := Option (Fin (N + 1))

/-- Weight of a slot. -/
def wt : Slot N → ℕ
  | none => 0
  | some m => (m : ℕ) - 3

/-- Size of the component of a variation in a slot. -/
noncomputable def comp : Slot N → P n N → ℝ
  | none, δ => ‖δ.1‖
  | some m, δ => ‖δ.2 m‖

/-- `a_j = ‖ζ_{3+j}‖` (and `0` if `3 + j > N`). -/
noncomputable def aOf (p : P n N) (j : ℕ) : ℝ :=
  if h : 3 + j < N + 1 then ‖p.2 ⟨3 + j, h⟩‖ else 0

/-- The weight function with an integer index; zero for negative index. -/
noncomputable def OmZ (v : ℤ) (ℓ : ℕ) (a : ℕ → ℝ) : ℝ := if v < 0 then 0 else Om v.toNat ℓ a

/-- `𝔭` is a symbol of weight `w` with factors of weight at most `ℓ` on the open set `A`. -/
structure IsSymbol {W : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W]
    (A : Set (P n N)) (w ℓ : ℕ) (𝔭 : P n N → W) : Prop where
  isOpen : IsOpen A
  smooth : ContDiffOn ℝ ∞ 𝔭 A
  bound : ∀ r : ℕ, ∃ C : ℝ, ∀ p ∈ A, ∀ δ : Fin r → P n N,
    ‖iteratedFDeriv ℝ r 𝔭 p δ‖ ≤
      C * ∑ s : Fin r → Slot N,
        OmZ ((w : ℤ) - ∑ i, (wt (s i) : ℤ)) ℓ (aOf p) * ∏ i, comp (s i) (δ i)

/-! ### Lemmas about `OmZ` -/

theorem OmZ_nonneg {v : ℤ} {ℓ : ℕ} {a : ℕ → ℝ} (ha : ∀ j, 0 ≤ a j) : 0 ≤ OmZ v ℓ a := by
  unfold OmZ
  split_ifs
  · exact le_rfl
  · exact zero_le_one.trans (one_le_Om _ _ ha)

theorem OmZ_natCast (v : ℕ) (ℓ : ℕ) (a : ℕ → ℝ) : OmZ (v : ℤ) ℓ a = Om v ℓ a := by
  unfold OmZ
  rw [if_neg (by omega)]
  simp

theorem OmZ_of_neg {v : ℤ} (h : v < 0) (ℓ : ℕ) (a : ℕ → ℝ) : OmZ v ℓ a = 0 := by
  unfold OmZ
  rw [if_pos h]

theorem OmZ_of_nonneg {v : ℤ} (h : 0 ≤ v) (ℓ : ℕ) (a : ℕ → ℝ) : OmZ v ℓ a = Om v.toNat ℓ a := by
  unfold OmZ
  rw [if_neg (by omega)]

theorem aOf_nonneg (p : P n N) (j : ℕ) : 0 ≤ aOf p j := by
  unfold aOf
  split_ifs
  · exact norm_nonneg _
  · exact le_rfl

theorem comp_nonneg (s : Slot N) (δ : P n N) : 0 ≤ comp s δ := by
  cases s <;> simp [comp]

theorem exists_OmZ_mono (v v' : ℤ) (hv : v ≤ v') :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ (ℓ ℓ' : ℕ) (a : ℕ → ℝ), ℓ ≤ ℓ' → (∀ j, 0 ≤ a j) →
      OmZ v ℓ a ≤ C * OmZ v' ℓ' a := by
  refine ⟨(v.toNat : ℝ) + 1, by positivity, fun ℓ ℓ' a hℓ ha => ?_⟩
  by_cases h : v < 0
  · rw [OmZ_of_neg h]
    exact mul_nonneg (by positivity) (OmZ_nonneg ha)
  · rw [OmZ_of_nonneg (by omega), OmZ_of_nonneg (by omega)]
    exact Om_mono (by omega) hℓ ha

theorem exists_OmZ_mul (v v' : ℤ) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ (ℓ ℓ' : ℕ) (a : ℕ → ℝ), (∀ j, 0 ≤ a j) →
      OmZ v ℓ a * OmZ v' ℓ' a ≤ C * OmZ (v + v') (max ℓ ℓ') a := by
  obtain ⟨C, hC, h⟩ := exists_Om_mul v.toNat v'.toNat
  refine ⟨C, hC, fun ℓ ℓ' a ha => ?_⟩
  by_cases h1 : v < 0
  · rw [OmZ_of_neg h1]
    simpa using mul_nonneg hC (OmZ_nonneg (v := v + v') (ℓ := max ℓ ℓ') ha)
  by_cases h2 : v' < 0
  · rw [OmZ_of_neg h2]
    simpa using mul_nonneg hC (OmZ_nonneg (v := v + v') (ℓ := max ℓ ℓ') ha)
  rw [OmZ_of_nonneg (by omega), OmZ_of_nonneg (by omega), OmZ_of_nonneg (by omega)]
  have : (v + v').toNat = v.toNat + v'.toNat := by omega
  rw [this]
  exact h ℓ ℓ' a ha

theorem exists_mul_OmZ (v : ℤ) {j : ℕ} (hj : 1 ≤ j) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ (ℓ : ℕ) (a : ℕ → ℝ), (∀ i, 0 ≤ a i) →
      a j * OmZ v ℓ a ≤ C * OmZ (v + j) (max ℓ j) a := by
  obtain ⟨C, hC, h⟩ := exists_mul_Om v.toNat hj
  refine ⟨C, hC, fun ℓ a ha => ?_⟩
  by_cases h1 : v < 0
  · rw [OmZ_of_neg h1]
    simpa using mul_nonneg hC (OmZ_nonneg (v := v + j) (ℓ := max ℓ j) ha)
  rw [OmZ_of_nonneg (by omega), OmZ_of_nonneg (by omega)]
  have : (v + j).toNat = v.toNat + j := by omega
  rw [this]
  exact h ℓ a ha

theorem bdSum_nonneg (w ℓ r : ℕ) (p : P n N) (δ : Fin r → P n N) :
    0 ≤ ∑ s : Fin r → Slot N,
        OmZ ((w : ℤ) - ∑ i, (wt (s i) : ℤ)) ℓ (aOf p) * ∏ i, comp (s i) (δ i) :=
  Finset.sum_nonneg fun s _ =>
    mul_nonneg (OmZ_nonneg (aOf_nonneg p)) (Finset.prod_nonneg fun i _ => comp_nonneg _ _)

/-! ### Basic consequences -/

section Basic

variable {W : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W]
  {A : Set (P n N)} {w ℓ : ℕ} {𝔭 : P n N → W}

theorem IsSymbol.norm_le (h : IsSymbol A w ℓ 𝔭) :
    ∃ C : ℝ, ∀ p ∈ A, ‖𝔭 p‖ ≤ C * Om w ℓ (aOf p) := by
  obtain ⟨C, hC⟩ := h.bound 0
  refine ⟨C, fun p hp => ?_⟩
  have := hC p hp (fun i => i.elim0)
  simpa [OmZ_natCast, iteratedFDeriv_zero_apply] using this

theorem IsSymbol.fderiv_le (h : IsSymbol A w ℓ 𝔭) :
    ∃ C : ℝ, ∀ p ∈ A, ∀ δ : P n N, ‖fderiv ℝ 𝔭 p δ‖ ≤
      C * ∑ s : Slot N, OmZ ((w : ℤ) - wt s) ℓ (aOf p) * comp s δ := by
  obtain ⟨C, hC⟩ := h.bound 1
  refine ⟨C, fun p hp δ => ?_⟩
  have h1 := hC p hp (fun _ => δ)
  rw [iteratedFDeriv_one_apply] at h1
  have e : ∑ s : Fin 1 → Slot N,
        OmZ ((w : ℤ) - ∑ i, (wt (s i) : ℤ)) ℓ (aOf p) * ∏ i, comp (s i) δ =
      ∑ s : Slot N, OmZ ((w : ℤ) - wt s) ℓ (aOf p) * comp s δ := by
    refine Fintype.sum_equiv (Equiv.funUnique (Fin 1) (Slot N)) _ _ fun s => ?_
    simp
  rw [e] at h1
  exact h1

theorem IsSymbol.fderiv_eq_zero (h : IsSymbol A w ℓ 𝔭) :
    ∀ p ∈ A, ∀ δ : P n N, (∀ s : Slot N, wt s ≤ w → comp s δ = 0) → fderiv ℝ 𝔭 p δ = 0 := by
  obtain ⟨C, hC⟩ := h.fderiv_le
  intro p hp δ hδ
  have hz : ∑ s : Slot N, OmZ ((w : ℤ) - wt s) ℓ (aOf p) * comp s δ = 0 := by
    refine Finset.sum_eq_zero fun s _ => ?_
    by_cases hs : wt s ≤ w
    · rw [hδ s hs, mul_zero]
    · rw [OmZ_of_neg (by omega), zero_mul]
  have := hC p hp δ
  rw [hz, mul_zero] at this
  exact norm_le_zero_iff.1 this

theorem IsSymbol.mono (h : IsSymbol A w ℓ 𝔭) {w' ℓ' : ℕ} (hw : w ≤ w') (hℓ : ℓ ≤ ℓ') :
    IsSymbol A w' ℓ' 𝔭 := by
  refine ⟨h.isOpen, h.smooth, fun r => ?_⟩
  obtain ⟨C, hC⟩ := h.bound r
  have hK : ∀ s : Fin r → Slot N, ∃ K : ℝ, 0 ≤ K ∧ ∀ (ℓ ℓ' : ℕ) (a : ℕ → ℝ), ℓ ≤ ℓ' →
      (∀ j, 0 ≤ a j) → OmZ ((w : ℤ) - ∑ i, (wt (s i) : ℤ)) ℓ a ≤
        K * OmZ ((w' : ℤ) - ∑ i, (wt (s i) : ℤ)) ℓ' a := fun s =>
    exists_OmZ_mono _ _ (by omega)
  choose K hK0 hK using hK
  set Kt : ℝ := ∑ s, K s with hKt
  have hKt0 : 0 ≤ Kt := Finset.sum_nonneg fun s _ => hK0 s
  refine ⟨|C| * Kt, fun p hp δ => ?_⟩
  refine (hC p hp δ).trans ?_
  have hS := bdSum_nonneg w ℓ r p δ
  have h1 : C * _ ≤ |C| * _ := mul_le_mul_of_nonneg_right (le_abs_self C) hS
  refine h1.trans ?_
  rw [mul_assoc]
  refine mul_le_mul_of_nonneg_left ?_ (abs_nonneg C)
  rw [Finset.mul_sum]
  refine Finset.sum_le_sum fun s _ => ?_
  have hP : 0 ≤ ∏ i, comp (s i) (δ i) := Finset.prod_nonneg fun i _ => comp_nonneg _ _
  have hle : K s ≤ Kt := Finset.single_le_sum (f := K) (fun s _ => hK0 s) (Finset.mem_univ s)
  calc OmZ ((w : ℤ) - ∑ i, (wt (s i) : ℤ)) ℓ (aOf p) * ∏ i, comp (s i) (δ i)
      ≤ (K s * OmZ ((w' : ℤ) - ∑ i, (wt (s i) : ℤ)) ℓ' (aOf p)) * ∏ i, comp (s i) (δ i) :=
        mul_le_mul_of_nonneg_right (hK s ℓ ℓ' _ hℓ (aOf_nonneg p)) hP
    _ ≤ (Kt * OmZ ((w' : ℤ) - ∑ i, (wt (s i) : ℤ)) ℓ' (aOf p)) * ∏ i, comp (s i) (δ i) :=
        mul_le_mul_of_nonneg_right
          (mul_le_mul_of_nonneg_right hle (OmZ_nonneg (aOf_nonneg p))) hP
    _ = Kt * (OmZ ((w' : ℤ) - ∑ i, (wt (s i) : ℤ)) ℓ' (aOf p) * ∏ i, comp (s i) (δ i)) := by
        ring

theorem IsSymbol.add {𝔮 : P n N → W} (h : IsSymbol A w ℓ 𝔭) (h' : IsSymbol A w ℓ 𝔮) :
    IsSymbol A w ℓ (fun p => 𝔭 p + 𝔮 p) := by
  refine ⟨h.isOpen, h.smooth.add h'.smooth, fun r => ?_⟩
  obtain ⟨C, hC⟩ := h.bound r
  obtain ⟨C', hC'⟩ := h'.bound r
  refine ⟨C + C', fun p hp δ => ?_⟩
  have e : iteratedFDeriv ℝ r (fun p => 𝔭 p + 𝔮 p) p =
      iteratedFDeriv ℝ r 𝔭 p + iteratedFDeriv ℝ r 𝔮 p :=
    iteratedFDeriv_add_apply
      ((h.smooth.contDiffAt (h.isOpen.mem_nhds hp)).of_le (by exact_mod_cast le_top))
      ((h'.smooth.contDiffAt (h.isOpen.mem_nhds hp)).of_le (by exact_mod_cast le_top))
  rw [e, ContinuousMultilinearMap.add_apply, add_mul]
  exact (norm_add_le _ _).trans (add_le_add (hC p hp δ) (hC' p hp δ))

theorem IsSymbol.clm_comp {W' : Type*} [NormedAddCommGroup W'] [NormedSpace ℝ W']
    (h : IsSymbol A w ℓ 𝔭) (L : W →L[ℝ] W') : IsSymbol A w ℓ (fun p => L (𝔭 p)) := by
  refine ⟨h.isOpen, L.contDiff.comp_contDiffOn h.smooth, fun r => ?_⟩
  obtain ⟨C, hC⟩ := h.bound r
  refine ⟨‖L‖ * C, fun p hp δ => ?_⟩
  have e : iteratedFDeriv ℝ r (fun p => L (𝔭 p)) p =
      L.compContinuousMultilinearMap (iteratedFDeriv ℝ r 𝔭 p) :=
    L.iteratedFDeriv_comp_left
      ((h.smooth.contDiffAt (h.isOpen.mem_nhds hp)).of_le (by exact_mod_cast le_top)) le_rfl
  rw [e, ContinuousLinearMap.compContinuousMultilinearMap_coe, Function.comp_apply, mul_assoc]
  exact (L.le_opNorm _).trans (mul_le_mul_of_nonneg_left (hC p hp δ) (norm_nonneg _))

theorem IsSymbol.const_smul (h : IsSymbol A w ℓ 𝔭) (c : ℝ) :
    IsSymbol A w ℓ (fun p => c • 𝔭 p) :=
  h.clm_comp (c • ContinuousLinearMap.id ℝ W)

theorem IsSymbol.neg (h : IsSymbol A w ℓ 𝔭) : IsSymbol A w ℓ (fun p => -𝔭 p) := by
  simpa using h.const_smul (-1)

theorem isSymbol_const (hA : IsOpen A) (c : W) : IsSymbol A 0 0 (fun _ : P n N => c) := by
  refine ⟨hA, contDiffOn_const, fun r => ?_⟩
  by_cases hr : r = 0
  · subst hr
    refine ⟨‖c‖, fun p hp δ => ?_⟩
    have h0 : OmZ 0 0 (aOf p) = 1 := by
      simpa [Om_zero] using OmZ_natCast 0 0 (aOf p)
    simp [iteratedFDeriv_zero_apply, h0]
  · refine ⟨0, fun p hp δ => ?_⟩
    rw [iteratedFDeriv_const_of_ne hr]
    simp

theorem IsSymbol.mono_set {B : Set (P n N)} (h : IsSymbol A w ℓ 𝔭) (hB : IsOpen B) (hBA : B ⊆ A) :
    IsSymbol B w ℓ 𝔭 := by
  refine ⟨hB, h.smooth.mono hBA, fun r => ?_⟩
  obtain ⟨C, hC⟩ := h.bound r
  exact ⟨C, fun p hp δ => hC p (hBA hp) δ⟩

end Basic

/-! ### Coordinate symbols -/

theorem OmZ_zero (ℓ : ℕ) (a : ℕ → ℝ) : OmZ 0 ℓ a = 1 := by
  simpa [Om_zero] using OmZ_natCast 0 ℓ a

theorem le_Om_self {j : ℕ} (hj : 1 ≤ j) {a : ℕ → ℝ} (ha : ∀ i, 0 ≤ a i) : a j ≤ Om j j a := by
  unfold Om
  have hj0 : (j : ℝ) ≠ 0 := Nat.cast_ne_zero.2 (by omega)
  have hmem : j ∈ Finset.Icc 1 (min j j) := by rw [Finset.mem_Icc]; omega
  have h := Finset.single_le_sum (f := fun i : ℕ => a i ^ ((j : ℝ) / i))
    (fun i _ => Real.rpow_nonneg (ha i) _) hmem
  simp only [div_self hj0, Real.rpow_one] at h
  linarith

theorem iteratedFDeriv_clm_succ_succ {E F : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    [NormedAddCommGroup F] [NormedSpace ℝ F] (L : E →L[ℝ] F) (r : ℕ) :
    iteratedFDeriv ℝ (r + 2) L = 0 := by
  have h1 : iteratedFDeriv ℝ 1 L = fun _ => iteratedFDeriv ℝ 1 L 0 := by
    funext x
    ext m
    rw [iteratedFDeriv_one_apply, iteratedFDeriv_one_apply, L.fderiv, L.fderiv]
  induction r with
  | zero =>
    funext x
    ext m
    show iteratedFDeriv ℝ (1 + 1) L x m = _
    rw [iteratedFDeriv_succ_apply_left, h1]
    simp
  | succ r ih =>
    funext x
    ext m
    show iteratedFDeriv ℝ ((r + 2) + 1) L x m = _
    rw [iteratedFDeriv_succ_apply_left, ih]
    simp

/-- The `m`-jet coordinate. -/
def slotProj (m : Fin (N + 1)) : P n N → Z n m := fun p => p.2 m

/-- `slotProj` as a continuous linear map. -/
noncomputable def slotProjL (m : Fin (N + 1)) : P n N →L[ℝ] Z n m :=
  (ContinuousLinearMap.proj m).comp (ContinuousLinearMap.snd ℝ _ _)

theorem slotProjL_apply (m : Fin (N + 1)) (p : P n N) : slotProjL m p = slotProj m p := rfl

theorem slotProj_eq (m : Fin (N + 1)) : (slotProj m : P n N → Z n m) = ⇑(slotProjL m) := rfl

/-- The derivative bounds of order `≥ 1` for the coordinate symbol. -/
theorem slotProj_bound_succ (m : Fin (N + 1)) {w ℓ : ℕ} (hw : w = wt (some m : Slot N)) (r : ℕ) :
    ∃ C : ℝ, ∀ (p : P n N) (δ : Fin (r + 1) → P n N),
      ‖iteratedFDeriv ℝ (r + 1) (slotProj m) p δ‖ ≤
        C * ∑ s : Fin (r + 1) → Slot N,
          OmZ ((w : ℤ) - ∑ i, (wt (s i) : ℤ)) ℓ (aOf p) * ∏ i, comp (s i) (δ i) := by
  rw [slotProj_eq]
  rcases r with _ | r
  · refine ⟨1, fun p δ => ?_⟩
    rw [iteratedFDeriv_one_apply, ContinuousLinearMap.fderiv]
    set s0 : Fin 1 → Slot N := fun _ => some m with hs0
    have hle := Finset.single_le_sum
      (f := fun s : Fin 1 → Slot N =>
        OmZ ((w : ℤ) - ∑ i, (wt (s i) : ℤ)) ℓ (aOf p) * ∏ i, comp (s i) (δ i))
      (fun s _ => mul_nonneg (OmZ_nonneg (aOf_nonneg p))
        (Finset.prod_nonneg fun i _ => comp_nonneg _ _)) (Finset.mem_univ s0)
    rw [one_mul]
    have e : OmZ ((w : ℤ) - ∑ i : Fin 1, (wt (s0 i) : ℤ)) ℓ (aOf p) *
        ∏ i : Fin 1, comp (s0 i) (δ i) = ‖(δ 0).2 m‖ := by
      simp [hs0, comp, hw, wt, OmZ_zero]
    have e2 : ‖slotProjL m (δ 0)‖ = ‖(δ 0).2 m‖ := rfl
    rw [e2, ← e]
    exact hle
  · refine ⟨0, fun p δ => ?_⟩
    rw [iteratedFDeriv_clm_succ_succ]
    simp

theorem isSymbol_slotProj {A : Set (P n N)} (hA : IsOpen A) (m : Fin (N + 1)) (hm : 4 ≤ (m : ℕ)) :
    IsSymbol A ((m : ℕ) - 3) ((m : ℕ) - 3) (slotProj m) := by
  refine ⟨hA, ?_, fun r => ?_⟩
  · rw [slotProj_eq]; exact (slotProjL m).contDiff.contDiffOn
  rcases r with _ | r
  · refine ⟨1, fun p hp δ => ?_⟩
    have hj : 3 + ((m : ℕ) - 3) < N + 1 := by have := m.2; omega
    have e : ‖p.2 m‖ = aOf p ((m : ℕ) - 3) := by
      unfold aOf
      have hmeq : (⟨3 + ((m : ℕ) - 3), hj⟩ : Fin (N + 1)) = m := Fin.ext (by simp; omega)
      rw [dif_pos hj, hmeq]
    have := le_Om_self (j := (m : ℕ) - 3) (by omega) (a := aOf p) (aOf_nonneg p)
    simpa [OmZ_natCast, iteratedFDeriv_zero_apply, slotProj, e] using this
  · obtain ⟨C, hC⟩ := slotProj_bound_succ (n := n) m (w := (m : ℕ) - 3) (ℓ := (m : ℕ) - 3)
      (by simp [wt]) r
    exact ⟨C, fun p _ δ => hC p δ⟩

theorem isSymbol_slotProj_low {A : Set (P n N)} (hA : IsOpen A) (m : Fin (N + 1))
    (hm : (m : ℕ) ≤ 3) (hbd : ∃ K : ℝ, ∀ p ∈ A, ‖p.2 m‖ ≤ K) : IsSymbol A 0 0 (slotProj m) := by
  refine ⟨hA, ?_, fun r => ?_⟩
  · rw [slotProj_eq]; exact (slotProjL m).contDiff.contDiffOn
  rcases r with _ | r
  · obtain ⟨K, hK⟩ := hbd
    refine ⟨K, fun p hp δ => ?_⟩
    have := hK p hp
    simpa [OmZ_natCast, iteratedFDeriv_zero_apply, slotProj, Om_zero, OmZ_zero] using this
  · obtain ⟨C, hC⟩ := slotProj_bound_succ (n := n) m (w := 0) (ℓ := 0) (by simp [wt]; omega) r
    exact ⟨C, fun p _ δ => hC p δ⟩

/-! ### Weight-zero symbols from functions with bounded derivatives -/

/-- The low part of a variation: base component and slots `m ≤ 3`. -/
def lo (δ : P n N) : P n N := (δ.1, fun m => if (m : ℕ) ≤ 3 then δ.2 m else 0)

/-- The high part of a variation: slots `m ≥ 4`. -/
def hi (δ : P n N) : P n N := (0, fun m => if (m : ℕ) ≤ 3 then 0 else δ.2 m)

theorem hi_add_lo (δ : P n N) : hi δ + lo δ = δ := by
  refine Prod.ext (by simp [hi, lo]) ?_
  funext m
  by_cases h : (m : ℕ) ≤ 3 <;> simp [hi, lo, h]

theorem hi_isHigh (δ : P n N) : (hi δ).1 = 0 ∧ ∀ m : Fin (N + 1), (m : ℕ) ≤ 3 → (hi δ).2 m = 0 := by
  refine ⟨rfl, fun m hm => ?_⟩
  simp [hi, hm]

theorem map_eq_lo {r : ℕ} {W : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W]
    (M : ContinuousMultilinearMap ℝ (fun _ : Fin r => P n N) W)
    (hM : ∀ δ : Fin r → P n N, (∃ i, (δ i).1 = 0 ∧ ∀ m : Fin (N + 1), (m : ℕ) ≤ 3 → (δ i).2 m = 0) →
      M δ = 0) (δ : Fin r → P n N) : M δ = M (fun i => lo (δ i)) := by
  have h : δ = (fun i => hi (δ i)) + (fun i => lo (δ i)) :=
    funext fun i => (hi_add_lo (δ i)).symm
  have h2 : M δ = M ((fun i => hi (δ i)) + (fun i => lo (δ i))) := congrArg _ h
  rw [h2, M.map_add_univ, Finset.sum_eq_single (∅ : Finset (Fin r))]
  · simp
  · intro s _ hs
    obtain ⟨i, hi'⟩ := Finset.nonempty_iff_ne_empty.2 hs
    apply hM
    refine ⟨i, ?_⟩
    rw [Finset.piecewise_eq_of_mem _ _ _ hi']
    exact hi_isHigh (δ i)
  · intro h; exact absurd (Finset.mem_univ _) h

/-- Indicator of the slots of weight zero. -/
noncomputable def g0 (t : Slot N) (δ : P n N) : ℝ := if wt t = 0 then comp t δ else 0

theorem g0_nonneg (t : Slot N) (δ : P n N) : 0 ≤ g0 t δ := by
  unfold g0; split_ifs
  · exact comp_nonneg _ _
  · exact le_rfl

theorem norm_lo_le (δ : P n N) : ‖lo δ‖ ≤ ∑ t : Slot N, g0 t δ := by
  have hnn : 0 ≤ ∑ t : Slot N, g0 t δ := Finset.sum_nonneg fun t _ => g0_nonneg t δ
  rw [lo, Prod.norm_def, max_le_iff]
  refine ⟨?_, ?_⟩
  · have := Finset.single_le_sum (f := fun t : Slot N => g0 t δ) (fun t _ => g0_nonneg t δ)
      (Finset.mem_univ (none : Slot N))
    simpa [g0, comp, wt] using this
  · refine (pi_norm_le_iff_of_nonneg hnn).2 fun m => ?_
    by_cases hm : (m : ℕ) ≤ 3
    · have := Finset.single_le_sum (f := fun t : Slot N => g0 t δ) (fun t _ => g0_nonneg t δ)
        (Finset.mem_univ (some m : Slot N))
      have h0 : wt (some m : Slot N) = 0 := by simp [wt]; omega
      simpa [g0, comp, h0, hm] using this
    · simpa [hm] using hnn

theorem isSymbol_of_low {W : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W]
    {A : Set (P n N)} (hA : IsOpen A) {c : P n N → W} (hc : ContDiffOn ℝ ∞ c A)
    (hbd : ∀ r : ℕ, ∃ C : ℝ, ∀ p ∈ A, ‖iteratedFDeriv ℝ r c p‖ ≤ C)
    (hlow : ∀ r : ℕ, ∀ p ∈ A, ∀ δ : Fin r → P n N,
      (∃ i, (δ i).1 = 0 ∧ ∀ m : Fin (N + 1), (m : ℕ) ≤ 3 → (δ i).2 m = 0) →
      iteratedFDeriv ℝ r c p δ = 0) :
    IsSymbol A 0 0 c := by
  refine ⟨hA, hc, fun r => ?_⟩
  obtain ⟨C, hC⟩ := hbd r
  refine ⟨C, fun p hp δ => ?_⟩
  set M := iteratedFDeriv ℝ r c p with hM
  have hCn : 0 ≤ C := (norm_nonneg _).trans (hC p hp)
  rw [map_eq_lo M (hlow r p hp) δ]
  have h1 : ‖M (fun i => lo (δ i))‖ ≤ ‖M‖ * ∏ i, ‖lo (δ i)‖ := M.le_opNorm _
  have h2 : ∏ i, ‖lo (δ i)‖ ≤ ∏ i : Fin r, ∑ t : Slot N, g0 t (δ i) :=
    Finset.prod_le_prod (fun i _ => norm_nonneg _) fun i _ => norm_lo_le (δ i)
  have h3 : ∏ i : Fin r, ∑ t : Slot N, g0 t (δ i) = ∑ s : Fin r → Slot N, ∏ i, g0 (s i) (δ i) := by
    rw [Finset.prod_univ_sum]
    simp
  have h4 : ∑ s : Fin r → Slot N, ∏ i, g0 (s i) (δ i) ≤ ∑ s : Fin r → Slot N,
      OmZ (((0 : ℕ) : ℤ) - ∑ i, (wt (s i) : ℤ)) 0 (aOf p) * ∏ i, comp (s i) (δ i) := by
    refine Finset.sum_le_sum fun s _ => ?_
    by_cases hs : ∀ i, wt (s i) = 0
    · have e1 : ∑ i, (wt (s i) : ℤ) = 0 := Finset.sum_eq_zero fun i _ => by simp [hs i]
      have e2 : ∏ i, g0 (s i) (δ i) = ∏ i, comp (s i) (δ i) :=
        Finset.prod_congr rfl fun i _ => by simp [g0, hs i]
      rw [e1, e2]
      simp [OmZ_zero]
    · push_neg at hs
      obtain ⟨i, hi'⟩ := hs
      have : ∏ i, g0 (s i) (δ i) = 0 :=
        Finset.prod_eq_zero (Finset.mem_univ i) (by simp [g0, hi'])
      rw [this]
      exact mul_nonneg (OmZ_nonneg (aOf_nonneg p))
        (Finset.prod_nonneg fun i _ => comp_nonneg _ _)
  have hP : 0 ≤ ∏ i, ‖lo (δ i)‖ := Finset.prod_nonneg fun i _ => norm_nonneg _
  calc ‖M (fun i => lo (δ i))‖ ≤ ‖M‖ * ∏ i, ‖lo (δ i)‖ := h1
    _ ≤ C * ∏ i, ‖lo (δ i)‖ := mul_le_mul_of_nonneg_right (hC p hp) hP
    _ ≤ C * ∑ s : Fin r → Slot N, ∏ i, g0 (s i) (δ i) :=
        mul_le_mul_of_nonneg_left (h2.trans h3.le) hCn
    _ ≤ _ := mul_le_mul_of_nonneg_left h4 hCn

/-- info: 'CWFlow.Extension.IsSymbol.norm_le' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms IsSymbol.norm_le

/-- info: 'CWFlow.Extension.IsSymbol.fderiv_le' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms IsSymbol.fderiv_le

/-- info: 'CWFlow.Extension.IsSymbol.fderiv_eq_zero' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms IsSymbol.fderiv_eq_zero

/-- info: 'CWFlow.Extension.IsSymbol.mono' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms IsSymbol.mono

/-- info: 'CWFlow.Extension.IsSymbol.add' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms IsSymbol.add

/-- info: 'CWFlow.Extension.IsSymbol.clm_comp' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms IsSymbol.clm_comp

/-- info: 'CWFlow.Extension.isSymbol_slotProj' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms isSymbol_slotProj

/-- info: 'CWFlow.Extension.isSymbol_slotProj_low' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms isSymbol_slotProj_low

/-- info: 'CWFlow.Extension.isSymbol_of_low' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms isSymbol_of_low

end CWFlow.Extension
