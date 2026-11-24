import CWFlow.Extension.EvolWeights
import CWFlow.Extension.JetOfSymbol

/-!
# The top order part of the second total derivative of the weight tensor

`liftTail c` applies a linear map `c` in the last `m` indices of an `(m+1)`-tensor.  The product rule
`Dz (c · T) = liftTail c (Dz T) + (D c) T` (`Dz_clm_apply`) lets the derivative index pass through the
coefficient `c = wCoef`, and `Dz (Dz (wSym k))` equals `wTop k` (the coefficient applied to the
`(k+5)`-jet coordinate) up to a symbol of weight `k + 1` (`isSymbol_Dz_Dz_wSym_sub`).
-/

open scoped ContDiff

namespace CWFlow.Extension

variable {n N : ℕ}

/-! ### Part 1: slices and lifting a linear map past derivative indices -/

/-- Fix the first index of a tensor. -/
noncomputable def sliceZ {m : ℕ} (a : Fin n) : Z n (m + 1) →L[ℝ] Z n m :=
  LinearMap.toContinuousLinearMap
    { toFun := fun X => WithLp.toLp 2 (fun I => X (Fin.cons a I))
      map_add' := fun X Y => by
        ext I
        simp
      map_smul' := fun c X => by
        ext I
        simp }

theorem sliceZ_apply {m : ℕ} (a : Fin n) (X : Z n (m + 1)) (I : Fin m → Fin n) :
    sliceZ a X I = X (Fin.cons a I) := rfl

/-- Apply a linear map in the last `m` indices of an `(m+1)`-tensor. -/
noncomputable def liftTail {m : ℕ} (c : Z n m →L[ℝ] Z n m) : Z n (m + 1) →L[ℝ] Z n (m + 1) :=
  LinearMap.toContinuousLinearMap
    { toFun := fun X => WithLp.toLp 2 (fun J => c (sliceZ (J 0) X) (Fin.tail J))
      map_add' := fun X Y => by
        ext J
        simp
      map_smul' := fun a X => by
        ext J
        simp }

theorem liftTail_apply {m : ℕ} (c : Z n m →L[ℝ] Z n m) (X : Z n (m + 1))
    (J : Fin (m + 1) → Fin n) : liftTail c X J = c (sliceZ (J 0) X) (Fin.tail J) := rfl

/-- `liftTail` is linear in `c`. -/
noncomputable def liftTailLin {m : ℕ} :
    (Z n m →L[ℝ] Z n m) →ₗ[ℝ] (Z n (m + 1) →L[ℝ] Z n (m + 1)) where
  toFun c := liftTail c
  map_add' c d := by
    ext X J
    simp [liftTail_apply]
  map_smul' a c := by
    ext X J
    simp [liftTail_apply]

/-- `liftTail` as a continuous linear map in `c`. -/
noncomputable def liftTailL {m : ℕ} :
    (Z n m →L[ℝ] Z n m) →L[ℝ] Z n (m + 1) →L[ℝ] Z n (m + 1) :=
  LinearMap.toContinuousLinearMap liftTailLin

theorem liftTailL_apply {m : ℕ} (c : Z n m →L[ℝ] Z n m) : liftTailL c = liftTail c := rfl

/-! ### Part 2: product rule for `Dz` -/

theorem wt_sliceZ_Dz {m : ℕ} (𝔗 : P n N → Z n m) (p : P n N) (a : Fin n) :
    sliceZ a (Dz 𝔗 p) = totalD 𝔗 p (EuclideanSpace.single a 1) := by
  ext I
  rw [sliceZ_apply]
  unfold Dz
  rw [curryZ_apply]
  simp

theorem Dz_apply_eq {m : ℕ} (𝔭 : P n N → Z n m) (p : P n N) (J : Fin (m + 1) → Fin n) :
    Dz 𝔭 p J = totalD 𝔭 p (EuclideanSpace.single (J 0) 1) (Fin.tail J) := by
  unfold Dz
  rw [curryZ_apply]

theorem Dz_clm_apply {m : ℕ} {𝔠 : P n N → (Z n m →L[ℝ] Z n m)} {𝔗 : P n N → Z n m} {p : P n N}
    (hc : DifferentiableAt ℝ 𝔠 p) (hT : DifferentiableAt ℝ 𝔗 p) (J : Fin (m + 1) → Fin n) :
    Dz (fun q => 𝔠 q (𝔗 q)) p J =
      liftTail (𝔠 p) (Dz 𝔗 p) J + (totalD 𝔠 p (EuclideanSpace.single (J 0) 1)) (𝔗 p) (Fin.tail J) := by
  rw [Dz_apply_eq]
  have h := totalD_bilin
    ((ContinuousLinearMap.apply ℝ (Z n m)).flip : (Z n m →L[ℝ] Z n m) →L[ℝ] Z n m →L[ℝ] Z n m)
    hc hT (EuclideanSpace.single (J 0) 1)
  have h' : totalD (fun q => 𝔠 q (𝔗 q)) p (EuclideanSpace.single (J 0) 1) =
      totalD 𝔠 p (EuclideanSpace.single (J 0) 1) (𝔗 p) +
        𝔠 p (totalD 𝔗 p (EuclideanSpace.single (J 0) 1)) := h
  rw [h', PiLp.add_apply, liftTail_apply, wt_sliceZ_Dz]
  ring

theorem wt_differentiableAt_wCoef (m : ℕ) :
    Differentiable ℝ (fun q : P n N => wCoef m (low23 q).1) :=
  ((contDiff_wCoef (n := n) m).comp
    (ContinuousLinearMap.contDiff ((ContinuousLinearMap.fst ℝ (Z n 2) (Z n 3)).comp low23))).differentiable
    (by simp)

/-- The special case of the weight coefficient. -/
theorem Dz_wCoef_apply {m : ℕ} {𝔗 : P n N → Z n m} {p : P n N} (hT : DifferentiableAt ℝ 𝔗 p)
    (J : Fin (m + 1) → Fin n) :
    Dz (fun q => wCoef m (low23 q).1 (𝔗 q)) p J =
      liftTail (wCoef m (low23 p).1) (Dz 𝔗 p) J +
        (totalD (fun q : P n N => wCoef m (low23 q).1) p (EuclideanSpace.single (J 0) 1)) (𝔗 p)
          (Fin.tail J) :=
  Dz_clm_apply (𝔠 := fun q : P n N => wCoef m (low23 q).1) (wt_differentiableAt_wCoef m p) hT J

/-! ### Part 3: the top part of `Dz (Dz (wSym k))` -/

/-- The top order part of the second total derivative of the weight tensor: the weight coefficient applied,
in the last `k+3` indices, to the `(k+5)`-jet coordinate. -/
noncomputable def wTop (k : ℕ) (p : P n N) : Z n (k + 3 + 1 + 1) :=
  liftTail (liftTail (wCoef (k + 3) (low23 p).1)) (slotN (k + 3 + 1 + 1) p)

/-- The weight coefficient as a function on jet space. -/
noncomputable def wt_c (m : ℕ) (q : P n N) : Z n m →L[ℝ] Z n m := wCoef m (low23 q).1

/-- Pairing of a derivative direction (first index) of a coefficient derivative with a tensor. -/
noncomputable def wt_pair {m m' : ℕ} (L : (Z n m →L[ℝ] Z n m) →L[ℝ] (Z n m' →L[ℝ] Z n m'))
    (D : EuclideanSpace ℝ (Fin n) →L[ℝ] (Z n m →L[ℝ] Z n m)) (T : Z n m') : Z n (m' + 1) :=
  WithLp.toLp 2 (fun J => L (D (EuclideanSpace.single (J 0) 1)) T (Fin.tail J))

theorem Dz_congr_of_eventuallyEq {k : ℕ} {f g : P n N → Z n k} {p : P n N}
    (h : f =ᶠ[nhds p] g) : Dz f p = Dz g p := by
  unfold Dz totalD
  rw [h.fderiv_eq]

theorem wt_diff {W : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W] {A : Set (P n N)}
    {w ℓ : ℕ} {f : P n N → W} (h : IsSymbol A w ℓ f) {p : P n N} (hp : p ∈ A) :
    DifferentiableAt ℝ f p :=
  (h.smooth.contDiffAt (h.isOpen.mem_nhds hp)).differentiableAt (by simp)

theorem wt_isSymbol_pair {A : Set (P n N)} {m m' w ℓ : ℕ} (hA : IsOpen A)
    (L : (Z n m →L[ℝ] Z n m) →L[ℝ] (Z n m' →L[ℝ] Z n m'))
    {D : P n N → (EuclideanSpace ℝ (Fin n) →L[ℝ] (Z n m →L[ℝ] Z n m))} {T : P n N → Z n m'}
    (hD : IsSymbol A 0 0 D) (hT : IsSymbol A w ℓ T) :
    IsSymbol A w (max 0 ℓ) (fun p => wt_pair L (D p) (T p)) := by
  refine isSymbol_of_coord hA fun J => ?_
  have h1 : IsSymbol A 0 0 (fun p => L (D p (EuclideanSpace.single (J 0) 1))) :=
    (hD.clm_comp (ContinuousLinearMap.apply ℝ (Z n m →L[ℝ] Z n m)
      (EuclideanSpace.single (J 0) 1))).clm_comp L
  have h2 := IsSymbol.bilin ((ContinuousLinearMap.apply ℝ (Z n m')).flip) h1 hT
  have h3 := h2.clm_comp (PiLp.proj 2 (fun _ : Fin m' → Fin n => ℝ) (Fin.tail J))
  rw [zero_add] at h3
  exact h3.congr fun p _ => rfl

theorem wt_Dz_zero {m : ℕ} : Dz (fun _ : P n N => (0 : Z n m)) = fun _ => 0 := by
  funext p
  unfold Dz
  rw [totalD_const]
  simp

theorem wt_sSym_zero : (sSym (N := N) (n := n) 0) = fun _ => 0 := by
  funext p
  simp [sSym, tSym_zero]

section Sym

variable {A : Set (P n N)}

theorem wt_isSymbol_sSym (hA : IsOpen A)
    (hlow : ∀ m : Fin (N + 1), (m : ℕ) ≤ 3 → ∃ K : ℝ, ∀ p ∈ A, ‖p.2 m‖ ≤ K)
    (k : ℕ) (hN : k + 3 ≤ N) : IsSymbol A k k (sSym (N := N) (n := n) k) := by
  rcases Nat.eq_zero_or_pos k with rfl | hk
  · rw [wt_sSym_zero]
    exact isSymbol_const hA 0
  · exact (isSymbol_sSym hA hlow k hk hN).mono (by omega) (by omega)

theorem wt_isSymbol_Dz_sSym (hA : IsOpen A)
    (hlow : ∀ m : Fin (N + 1), (m : ℕ) ≤ 3 → ∃ K : ℝ, ∀ p ∈ A, ‖p.2 m‖ ≤ K)
    (k : ℕ) (hN : k + 4 ≤ N) : IsSymbol A (k + 1) (k + 1) (Dz (sSym (N := N) (n := n) k)) := by
  rcases Nat.eq_zero_or_pos k with rfl | hk
  · rw [wt_sSym_zero, wt_Dz_zero]
    exact (isSymbol_const hA 0).mono (Nat.zero_le _) (Nat.zero_le _)
  · have h := isSymbol_sSym hA hlow k hk (by omega)
    exact (h.Dz (by omega) hlow).mono (by omega) (by omega)

theorem wt_isSymbol_Dz_Dz_sSym (hA : IsOpen A)
    (hlow : ∀ m : Fin (N + 1), (m : ℕ) ≤ 3 → ∃ K : ℝ, ∀ p ∈ A, ‖p.2 m‖ ≤ K)
    (k : ℕ) (hN : k + 5 ≤ N) :
    IsSymbol A (k + 1) (k + 1) (Dz (Dz (sSym (N := N) (n := n) k))) := by
  rcases Nat.eq_zero_or_pos k with rfl | hk
  · rw [wt_sSym_zero, wt_Dz_zero, wt_Dz_zero]
    exact (isSymbol_const hA 0).mono (Nat.zero_le _) (Nat.zero_le _)
  · have h := isSymbol_sSym hA hlow k hk (by omega)
    have h1 := h.Dz (by omega) hlow
    have h2 := h1.Dz (by omega) hlow
    refine h2.mono ?_ ?_ <;> omega

theorem wt_Dz_Dz_tSym (hA : IsOpen A)
    (hlow : ∀ m : Fin (N + 1), (m : ℕ) ≤ 3 → ∃ K : ℝ, ∀ p ∈ A, ‖p.2 m‖ ≤ K)
    (k : ℕ) (hN : k + 5 ≤ N) {p : P n N} (hp : p ∈ A) :
    Dz (Dz (tSym (N := N) (n := n) k)) p =
      Dz (Dz (sSym (N := N) (n := n) k)) p + slotN (k + 3 + 1 + 1) p := by
  have hS := wt_isSymbol_sSym hA hlow k (n := n) (by omega)
  have eT : tSym (N := N) (n := n) k = fun q => sSym k q + slotN (k + 3) q := by
    funext q
    simp [sSym]
  have eDT : Dz (tSym (N := N) (n := n) k) =ᶠ[nhds p]
      fun q => Dz (sSym (N := N) (n := n) k) q + slotN (k + 3 + 1) q := by
    refine Filter.eventuallyEq_of_mem (hA.mem_nhds hp) fun q hq => ?_
    have hadd := Dz_add (wt_diff hS hq) (slotN (n := n) (N := N) (k + 3)).differentiableAt
    rw [Dz_slotN (n := n) (k + 3) (by omega)] at hadd
    rw [eT]
    exact hadd
  rw [Dz_congr_of_eventuallyEq eDT]
  have hadd := Dz_add (wt_diff (wt_isSymbol_Dz_sSym hA hlow k (by omega)) hp)
    (slotN (n := n) (N := N) (k + 3 + 1)).differentiableAt
  rw [Dz_slotN (n := n) (k + 3 + 1) (by omega)] at hadd
  exact hadd

end Sym

theorem isSymbol_Dz_Dz_wSym_sub {A : Set (P n N)} (hA : IsOpen A)
    (hlow : ∀ m : Fin (N + 1), (m : ℕ) ≤ 3 → ∃ K : ℝ, ∀ p ∈ A, ‖p.2 m‖ ≤ K)
    (k : ℕ) (hN : k + 5 ≤ N) :
    IsSymbol A (k + 1) (k + 1) (fun p => Dz (Dz (wSym (N := N) k)) p - wTop k p) := by
  have hc0 : IsSymbol A 0 0 (wt_c (n := n) (N := N) (k + 3)) := isSymbol_wCoef hA hlow (k + 3)
  have hdc : IsSymbol A 0 0 (totalD (wt_c (n := n) (N := N) (k + 3))) :=
    isSymbol_totalD_comp_low2 hA (bdd_low23 hlow) (contDiff_wCoef (k + 3))
  have hcd : Differentiable ℝ (wt_c (n := n) (N := N) (k + 3)) := wt_differentiableAt_wCoef (k + 3)
  have hT : IsSymbol A k k (tSym (N := N) (n := n) k) := isSymbol_tSym hA hlow k (by omega)
  have hDT : IsSymbol A (k + 1) (k + 1) (Dz (tSym (N := N) (n := n) k)) :=
    (hT.Dz (by omega) hlow).mono le_rfl (by omega)
  -- the two error terms
  let E1 : P n N → Z n (k + 3 + 1) := fun p =>
    wt_pair (ContinuousLinearMap.id ℝ _) (totalD (wt_c (k + 3)) p) (tSym k p)
  let E2 : P n N → Z n (k + 3 + 1 + 1) := fun p =>
    wt_pair liftTailL (totalD (wt_c (k + 3)) p) (Dz (tSym k) p)
  have hE1 : IsSymbol A k k E1 :=
    (wt_isSymbol_pair hA (ContinuousLinearMap.id ℝ _) hdc hT).mono (by omega) (by omega)
  have hE2 : IsSymbol A (k + 1) (k + 1) E2 :=
    (wt_isSymbol_pair hA liftTailL hdc hDT).mono (by omega) (by omega)
  have hDE1 : IsSymbol A (k + 1) (k + 1) (Dz E1) :=
    (hE1.Dz (by omega) hlow).mono le_rfl (by omega)
  have hG : IsSymbol A (k + 1) (k + 1)
      (fun q => liftTail (wt_c (k + 3) q) (Dz (tSym k) q)) := by
    have h1 := IsSymbol.bilin
      ((ContinuousLinearMap.apply ℝ (Z n (k + 3 + 1))).flip) (hc0.clm_comp liftTailL) hDT
    rw [zero_add] at h1
    exact (h1.mono le_rfl (by omega)).congr fun p _ => rfl
  have hLL : IsSymbol A 0 0 (fun q => liftTail (liftTail (wt_c (n := n) (N := N) (k + 3) q))) :=
    (hc0.clm_comp liftTailL).clm_comp liftTailL
  have hR1 : IsSymbol A (k + 1) (k + 1)
      (fun p => liftTail (liftTail (wt_c (k + 3) p)) (Dz (Dz (sSym (N := N) (n := n) k)) p)) := by
    have h1 := IsSymbol.bilin
      ((ContinuousLinearMap.apply ℝ (Z n (k + 3 + 1 + 1))).flip) hLL
      (wt_isSymbol_Dz_Dz_sSym hA hlow k hN)
    rw [zero_add] at h1
    exact (h1.mono le_rfl (by omega)).congr fun p _ => rfl
  have hR := (hR1.add hE2).add hDE1
  have hW : (wSym (N := N) (n := n) k) = fun q => wCoef (k + 3) (low23 q).1 (tSym k q) :=
    funext (wSym_eq k)
  refine hR.congr fun p hp => ?_
  -- step 1
  have hW1 : ∀ q ∈ A, Dz (wSym (N := N) (n := n) k) q =
      liftTail (wt_c (k + 3) q) (Dz (tSym k) q) + E1 q := by
    intro q hq
    ext J
    rw [hW]
    exact Dz_wCoef_apply (wt_diff hT hq) J
  have hev : Dz (wSym (N := N) (n := n) k) =ᶠ[nhds p]
      fun q => liftTail (wt_c (k + 3) q) (Dz (tSym k) q) + E1 q :=
    Filter.eventuallyEq_of_mem (hA.mem_nhds hp) hW1
  have hA1 := Dz_congr_of_eventuallyEq hev
  have hA2 := Dz_add (wt_diff hG hp) (wt_diff hE1 hp)
  -- step 2
  have hc' : DifferentiableAt ℝ (fun q : P n N => liftTail (wt_c (k + 3) q)) p :=
    (liftTailL.differentiable.comp hcd) p
  have hTotal : ∀ e : EuclideanSpace ℝ (Fin n),
      totalD (fun q : P n N => liftTail (wt_c (k + 3) q)) p e =
        liftTailL (totalD (wt_c (k + 3)) p e) := by
    intro e
    have := totalD_clm_comp (𝔭 := wt_c (n := n) (N := N) (k + 3)) liftTailL (hcd p)
    exact congrArg (fun L => L e) this
  have hA3 : Dz (fun q => liftTail (wt_c (k + 3) q) (Dz (tSym k) q)) p =
      liftTail (liftTail (wt_c (k + 3) p)) (Dz (Dz (tSym k)) p) + E2 p := by
    ext J
    have := Dz_clm_apply (𝔠 := fun q : P n N => liftTail (wt_c (k + 3) q))
      (𝔗 := Dz (tSym (N := N) (n := n) k)) hc' (wt_diff hDT hp) J
    rw [this, hTotal]
    rfl
  have hA4 := wt_Dz_Dz_tSym hA hlow k hN hp (n := n)
  have hw : wTop k p = liftTail (liftTail (wt_c (k + 3) p)) (slotN (k + 3 + 1 + 1) p) := rfl
  rw [hA1, hA2, hA3, hA4, map_add, hw]
  abel

/-- info: 'CWFlow.Extension.Dz_wCoef_apply' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms Dz_wCoef_apply

/-- info: 'CWFlow.Extension.isSymbol_Dz_Dz_wSym_sub' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms isSymbol_Dz_Dz_wSym_sub

end CWFlow.Extension
