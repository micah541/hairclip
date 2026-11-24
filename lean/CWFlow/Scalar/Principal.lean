import CWFlow.Scalar.JetFunction

/-!
# The principal part of the scalar operator

On a graph of a closed form (`z₂` symmetric), the jet function is affine in the top jet `z₄`:
`Pjet(y, z₁, z₂, z₃, z₄ + ζ) = Pjet(y, z) − tr_g^{(a,c)} tr_g^{(b,d)} ζ(a, b, d, c)` (`Pjet_add_top`), where `g` is
the induced metric of `F_u`. The ingredients:

* Faà di Bruno in degree one: `T₁ = q₁ ∘ p₁` (`taylorComp_one_apply`), so `DF_u(a) = DΨ(a, z₂(a, ·))` (`LOf_apply`);
* in degree three a top-degree perturbation of `p` enters only through the one-block partition
  (`taylorComp_add_top`, `taylorComp_top_apply`);
* `PotM` is affine in `T₃` (`PotM_add`), and the linear part reduces, by the Lagrangian property and the vertical
  pairing `⟨J DΨ(0, η), DF_u c⟩ = η(c)`, to `tr tr ζ`.
-/

open Filter Set Function

open scoped InnerProductSpace RealInnerProductSpace Topology ContDiff

noncomputable section

namespace CWFlow

/-! ### Ordered finpartitions with a block of full size -/

theorem OrderedFinpartition.length_eq_one_of_partSize_eq {n : ℕ} (c : OrderedFinpartition n)
    {m : Fin c.length} (h : c.partSize m = n) : c.length = 1 := by
  have hbij : Function.Bijective (c.emb m) :=
    (Fintype.bijective_iff_injective_and_card (c.emb m)).2
      ⟨(c.emb_strictMono m).injective, by simp [h]⟩
  have hall : ∀ m' : Fin c.length, m' = m := by
    intro m'
    by_contra hne
    have hpos := c.partSize_pos m'
    obtain ⟨j, hj⟩ := hbij.2 (c.emb m' ⟨0, hpos⟩)
    have hdisj := c.disjoint (mem_univ m') (mem_univ m) hne
    exact Set.disjoint_left.1 hdisj (mem_range_self _) ⟨j, hj⟩
  have h1 : 0 < c.length := Fin.pos m
  by_contra hne
  have h2 : 1 < c.length := by omega
  exact absurd (hall ⟨0, h1⟩) (by
    intro h0
    have h1' := hall ⟨1, h2⟩
    rw [← h0] at h1'
    simp [Fin.ext_iff] at h1')

/-- The partition of `Fin 3` into one block. -/
def onePart3 : OrderedFinpartition 3 where
  length := 1
  partSize _ := 3
  partSize_pos _ := by norm_num
  emb _ := id
  emb_strictMono _ := strictMono_id
  parts_strictMono := Subsingleton.strictMono _
  disjoint i _ j _ h := (h (Subsingleton.elim i j)).elim
  cover x := ⟨0, x, rfl⟩

theorem eq_onePart3 (c : OrderedFinpartition 3) (h : c.length = 1) : c = onePart3 := by
  have hsize : ∀ i, c.partSize i = 3 := by
    intro i
    have hsurj : Function.Surjective (c.emb i) := by
      intro x
      obtain ⟨m, hm⟩ := c.cover x
      have : m = i := by
        have := Fin.subsingleton_iff_le_one.2 h.le
        exact Subsingleton.elim m i
      subst this
      exact hm
    have := Fintype.card_le_of_surjective _ hsurj
    simp only [Fintype.card_fin] at this
    exact le_antisymm (c.partSize_le i) this
  rcases c with ⟨length, partSize, hpos, emb, hmono, hpm, hdisj, hcover⟩
  simp only at h hsize
  subst h
  obtain rfl : partSize = fun _ => 3 := funext hsize
  have hid : ∀ i, emb i = id := by
    intro i
    have hm := hmono i
    funext x
    have h0 : (emb i 0 : ℕ) < emb i 1 := hm (show (0 : Fin 3) < 1 by decide)
    have h1 : (emb i 1 : ℕ) < emb i 2 := hm (show (1 : Fin 3) < 2 by decide)
    have h2 : (emb i 2 : ℕ) < 3 := (emb i 2).2
    fin_cases x <;> ext <;> simp <;> omega
  have : emb = fun _ => id := funext hid
  subst this
  rfl

/-! ### Faà di Bruno in degrees one and three -/

universe u

variable {E F G : Type u} [NormedAddCommGroup E] [NormedSpace ℝ E] [NormedAddCommGroup F]
  [NormedSpace ℝ F] [NormedAddCommGroup G] [NormedSpace ℝ G]

theorem taylorComp_one_apply (q : FormalMultilinearSeries ℝ F G) (p : FormalMultilinearSeries ℝ E F)
    (v : Fin 1 → E) : q.taylorComp p 1 v = q 1 (fun m => p 1 (fun _ => v m)) := by
  simp only [FormalMultilinearSeries.taylorComp, Fintype.sum_unique,
    ContinuousMultilinearMap.sum_apply]
  rfl

/-- A formal series concentrated in degree three. -/
def topSeries (ζ : E [×3]→L[ℝ] F) : FormalMultilinearSeries ℝ E F
  | 3 => ζ
  | _ => 0

theorem topSeries_ne {ζ : E [×3]→L[ℝ] F} {j : ℕ} (hj : j ≠ 3) : topSeries ζ j = 0 := by
  match j, hj with
  | 0, _ => rfl
  | 1, _ => rfl
  | 2, _ => rfl
  | n + 4, _ => rfl

theorem compAlong_zero_of_length_ne {n : ℕ} (hn : 0 < n) (c : OrderedFinpartition n)
    (f : F [×c.length]→L[ℝ] G) : c.compAlongOrderedFinpartition f (fun m => (0 : E [×c.partSize m]→L[ℝ] F)) = 0 := by
  have m₀ : Fin c.length := ⟨0, c.length_pos hn⟩
  rw [← OrderedFinpartition.compAlongOrderedFinpartitionL_apply]
  exact (c.compAlongOrderedFinpartitionL ℝ E F G f).map_coord_zero m₀ rfl

/-- In degree three, a top-degree perturbation enters additively. -/
theorem taylorComp_add_top (q : FormalMultilinearSeries ℝ F G) (p : FormalMultilinearSeries ℝ E F)
    (ζ : E [×3]→L[ℝ] F) :
    q.taylorComp (p + topSeries ζ) 3 = q.taylorComp p 3 + q.taylorComp (topSeries ζ) 3 := by
  simp only [FormalMultilinearSeries.taylorComp, ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun c _ => ?_
  simp only [FormalMultilinearSeries.compAlongOrderedFinpartition]
  by_cases hc : c.length = 1
  · have : Subsingleton (Fin c.length) := by rw [hc]; infer_instance
    let i0 : Fin c.length := ⟨0, by omega⟩
    set A : ∀ m : Fin c.length, E [×c.partSize m]→L[ℝ] F := fun m => p (c.partSize m)
    set B : ∀ m : Fin c.length, E [×c.partSize m]→L[ℝ] F := fun m => topSeries ζ (c.partSize m)
    have hAB : (fun m => (p + topSeries ζ) (c.partSize m)) = Function.update A i0 (A i0 + B i0) := by
      funext m
      rw [Subsingleton.elim m i0, Function.update_self]
      rfl
    have hA : A = Function.update A i0 (A i0) := by simp
    have hB : B = Function.update A i0 (B i0) := by
      funext m
      rw [Subsingleton.elim m i0, Function.update_self]
    rw [hAB, ← OrderedFinpartition.compAlongOrderedFinpartitionL_apply,
      ← OrderedFinpartition.compAlongOrderedFinpartitionL_apply,
      ← OrderedFinpartition.compAlongOrderedFinpartitionL_apply,
      ContinuousMultilinearMap.map_update_add, ← hA, ← hB]
  · have hs : ∀ m, topSeries ζ (c.partSize m) = 0 := fun m =>
      topSeries_ne fun h => hc (OrderedFinpartition.length_eq_one_of_partSize_eq c h)
    have e1 : (fun m => (p + topSeries ζ) (c.partSize m)) = fun m => p (c.partSize m) := by
      funext m
      change p (c.partSize m) + topSeries ζ (c.partSize m) = _
      rw [hs, add_zero]
    have e2 : (fun m => topSeries ζ (c.partSize m)) = fun m => (0 : E [×c.partSize m]→L[ℝ] F) :=
      funext hs
    rw [e1, e2, compAlong_zero_of_length_ne (by norm_num), add_zero]

/-- The linear part: only the one-block partition contributes. -/
theorem taylorComp_top_apply (q : FormalMultilinearSeries ℝ F G) (ζ : E [×3]→L[ℝ] F)
    (v : Fin 3 → E) : q.taylorComp (topSeries ζ) 3 v = q 1 (fun _ => ζ v) := by
  simp only [FormalMultilinearSeries.taylorComp]
  rw [Finset.sum_eq_single onePart3]
  · rfl
  · intro c _ hne
    have hc : c.length ≠ 1 := fun h => hne (eq_onePart3 c h)
    have e2 : (fun m => topSeries ζ (c.partSize m)) = fun m => (0 : E [×c.partSize m]→L[ℝ] F) :=
      funext fun m => topSeries_ne fun h => hc (OrderedFinpartition.length_eq_one_of_partSize_eq c h)
    simp only [FormalMultilinearSeries.compAlongOrderedFinpartition]
    rw [e2, compAlong_zero_of_length_ne (by norm_num)]
  · intro h
    exact absurd (Finset.mem_univ _) h

/-! ### The top jet of the scalar operator -/

section Geometry

variable {E V : Type u} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
  [NormedAddCommGroup V] [InnerProductSpace ℝ V] {J : V →L[ℝ] V} {Ψ : E × (E →L[ℝ] ℝ) → V}

/-- Change the top jet: `z₄ ↦ z₄ + ζ`. -/
def addTop (z : SJet E) (ζ : E [×4]→L[ℝ] ℝ) : SJet E := (z.1, z.2.1, z.2.2.1, z.2.2.2 + ζ)

/-- The corresponding change of the degree-three Taylor coefficient of the graph map. -/
def topCorr (ζ : E [×4]→L[ℝ] ℝ) : E [×3]→L[ℝ] (E × (E →L[ℝ] ℝ)) :=
  (0 : E [×3]→L[ℝ] E).prod (continuousMultilinearCurryRightEquiv' ℝ 3 E ℝ ζ)

omit [FiniteDimensional ℝ E] in
theorem gSeries_addTop_ne (y : E) (z : SJet E) (ζ : E [×4]→L[ℝ] ℝ) {j : ℕ} (hj : j ≠ 3) :
    gSeries y (addTop z ζ) j = gSeries y z j := by
  match j, hj with
  | 0, _ => rfl
  | 1, _ => rfl
  | 2, _ => rfl
  | n + 4, _ => rfl

omit [FiniteDimensional ℝ E] in
theorem gSeries_addTop_three (y : E) (z : SJet E) (ζ : E [×4]→L[ℝ] ℝ) :
    gSeries y (addTop z ζ) 3 = gSeries y z 3 + topCorr ζ := by
  change (0 : E [×3]→L[ℝ] E).prod (continuousMultilinearCurryRightEquiv' ℝ 3 E ℝ (z.2.2.2 + ζ)) =
    (0 : E [×3]→L[ℝ] E).prod (continuousMultilinearCurryRightEquiv' ℝ 3 E ℝ z.2.2.2) +
      (0 : E [×3]→L[ℝ] E).prod (continuousMultilinearCurryRightEquiv' ℝ 3 E ℝ ζ)
  rw [map_add]
  ext m
  · simp
  · simp

/-- `taylorComp_add_top`, with the perturbation given coefficientwise. -/
theorem taylorComp_add_top' (q : FormalMultilinearSeries ℝ F G)
    {p p' : FormalMultilinearSeries ℝ E F} {δ : E [×3]→L[ℝ] F} (h3 : p' 3 = p 3 + δ)
    (hne : ∀ j, j ≠ 3 → p' j = p j) :
    q.taylorComp p' 3 = q.taylorComp p 3 + q.taylorComp (topSeries δ) 3 := by
  have e : p' = p + topSeries δ := by
    funext j
    by_cases hj : j = 3
    · subst hj
      exact h3
    · rw [hne j hj]
      change p j = p j + topSeries δ j
      rw [topSeries_ne hj, add_zero]
  rw [e]
  exact taylorComp_add_top q p δ

theorem TkOf_addTop_of_le {k : ℕ} (hk : k ≤ 2) (y : E) (z : SJet E) (ζ : E [×4]→L[ℝ] ℝ) :
    TkOf Ψ k y (addTop z ζ) = TkOf Ψ k y z :=
  taylorComp_congr (fun _ _ => rfl) fun j _ hjk => gSeries_addTop_ne y z ζ (by omega)

theorem LOf_addTop (y : E) (z : SJet E) (ζ : E [×4]→L[ℝ] ℝ) :
    LOf Ψ y (addTop z ζ) = LOf Ψ y z := by
  rw [LOf, LOf, TkOf_addTop_of_le (by norm_num)]

theorem TkOf_three_addTop (y : E) (z : SJet E) (ζ : E [×4]→L[ℝ] ℝ) :
    TkOf Ψ 3 y (addTop z ζ) =
      TkOf Ψ 3 y z + (ftaylorSeries ℝ Ψ (y, z.1)).taylorComp (topSeries (topCorr ζ)) 3 :=
  taylorComp_add_top' _ (gSeries_addTop_three y z ζ) fun j hj => gSeries_addTop_ne y z ζ hj

theorem topTerm_apply (y : E) (z : SJet E) (ζ : E [×4]→L[ℝ] ℝ) (v : Fin 3 → E) :
    (ftaylorSeries ℝ Ψ (y, z.1)).taylorComp (topSeries (topCorr ζ)) 3 v =
      fderiv ℝ Ψ (y, z.1) (0, continuousMultilinearCurryRightEquiv' ℝ 3 E ℝ ζ v) := by
  rw [taylorComp_top_apply]
  simp [ftaylorSeries, iteratedFDeriv_one_apply, topCorr]

/-- `DF_u(a) = DΨ(a, z₂(a, ·))`. -/
theorem LOf_apply (y : E) (z : SJet E) (a : E) :
    LOf Ψ y z a = fderiv ℝ Ψ (y, z.1) (a, slot1L z.2.1 a) := by
  simp only [LOf, continuousMultilinearCurryFin1_apply, TkOf, taylorComp_one_apply]
  simp [ftaylorSeries, iteratedFDeriv_one_apply, gSeries, Fin.snoc]

theorem HOfM_add_smul (L : E →L[ℝ] V) (T Q : E [×2]→L[ℝ] V) (s : ℝ) :
    HOfM L (T + s • Q) = HOfM L T + s • HOfM L Q := by
  simp only [HOfM, trOf, ContinuousMultilinearMap.add_apply, ContinuousMultilinearMap.smul_apply,
    perpOf_add, perpOf_smul, Finset.sum_add_distrib, Finset.smul_sum]

/-- `HOfM` is linear in the second-derivative slot, so its vertical derivative is itself. -/
theorem fderiv_HOfMP_vertical {L : E →L[ℝ] V} (hL : Function.Injective L) (T Q : E [×2]→L[ℝ] V) :
    fderiv ℝ HOfMP (L, T) (0, Q) = HOfM L Q := by
  have hd := (contDiffAt_HOfMP (p := (L, T)) hL).differentiableAt (by simp)
  have h1 : HasDerivAt (fun s : ℝ => T + s • Q) Q 0 := by
    simpa using ((hasDerivAt_id (0 : ℝ)).smul_const Q).const_add T
  have h2 := (hasDerivAt_const (0 : ℝ) L).prodMk h1
  have h3 := hd.hasFDerivAt.comp_hasDerivAt_of_eq (0 : ℝ) h2 (by simp)
  have h4 : HasDerivAt (fun s : ℝ => HOfM L T + s • HOfM L Q) (HOfM L Q) 0 := by
    simpa using ((hasDerivAt_id (0 : ℝ)).smul_const (HOfM L Q)).const_add (HOfM L T)
  have e : (HOfMP ∘ fun s : ℝ => (L, T + s • Q)) = fun s => HOfM L T + s • HOfM L Q :=
    funext fun s => HOfM_add_smul L T Q s
  rw [e] at h3
  exact h3.unique h4

/-- `PotM` is affine in the third-derivative slot. -/
theorem PotM_add (J : V →L[ℝ] V) {L : E →L[ℝ] V} (hL : Function.Injective L) (T₂ : E [×2]→L[ℝ] V)
    (T₃ Λ : E [×3]→L[ℝ] V) :
    PotM J L T₂ (T₃ + Λ) = PotM J L T₂ T₃ + trOf L fun a c => ⟪J (HOfM L (Λ.curryLeft a)), L c⟫ := by
  simp only [PotM, trOf, ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun i _ => ?_
  have hc : (T₃ + Λ).curryLeft (bas E i) = T₃.curryLeft (bas E i) + Λ.curryLeft (bas E i) := by
    ext m
    simp
  have hp : ((slot1 T₂ (bas E i), T₃.curryLeft (bas E i) + Λ.curryLeft (bas E i)) :
      (E →L[ℝ] V) × (E [×2]→L[ℝ] V)) =
      (slot1 T₂ (bas E i), T₃.curryLeft (bas E i)) + (0, Λ.curryLeft (bas E i)) := by simp
  rw [hc, hp, map_add, fderiv_HOfMP_vertical hL, map_add, inner_add_left]

variable {O : Set (E × (E →L[ℝ] ℝ))}
  (hs : ∀ q ∈ O, ∀ a b, ⟪J (fderiv ℝ Ψ q a), fderiv ℝ Ψ q b⟫ = omegaCanVal a b)
include hs

/-- On a graph of a closed form, `DF_u` is Lagrangian, so `J perp` pairs with `DF_u` like `J`. -/
theorem inner_J_perpOf_LOf {y : E} {z : SJet E} (hq : (y, z.1) ∈ O)
    (hsym : ∀ a b, z.2.1 ![a, b] = z.2.1 ![b, a]) (v : V) (c : E) :
    ⟪J (perpOf (LOf Ψ y z) v), LOf Ψ y z c⟫ = ⟪J v, LOf Ψ y z c⟫ := by
  have hlag : ∀ w, ⟪J (LOf Ψ y z w), LOf Ψ y z c⟫ = 0 := fun w => by
    rw [LOf_apply, LOf_apply, hs _ hq]
    simp only [omegaCanVal]
    rw [slot1L_apply, slot1L_apply, hsym, sub_self]
  rw [perpOf, map_sub, inner_sub_left, hlag, sub_zero]

theorem inner_J_topTerm {y : E} {z : SJet E} (hq : (y, z.1) ∈ O) (ζ : E [×4]→L[ℝ] ℝ)
    (v : Fin 3 → E) (c : E) :
    ⟪J (fderiv ℝ Ψ (y, z.1) (0, continuousMultilinearCurryRightEquiv' ℝ 3 E ℝ ζ v)),
      LOf Ψ y z c⟫ = ζ (Fin.snoc v c) := by
  rw [LOf_apply, inner_J_fderiv_vertical (hs _ hq)]
  rfl

/-- **The principal part of the scalar operator.** On the graph of a closed form, the jet function is affine in
the top jet: `Pjet(y, z₁, z₂, z₃, z₄ + ζ) = Pjet(y, z) − tr_g^{(a,c)} tr_g^{(b,d)} ζ(a, b, d, c)`. -/
theorem Pjet_add_top {y : E} {z : SJet E} (hq : (y, z.1) ∈ O)
    (hsym : ∀ a b, z.2.1 ![a, b] = z.2.1 ![b, a]) (hL : Function.Injective (LOf Ψ y z))
    (ζ : E [×4]→L[ℝ] ℝ) :
    Pjet J Ψ (y, addTop z ζ) = Pjet J Ψ (y, z) -
      trOf (LOf Ψ y z) fun a c => trOf (LOf Ψ y z) fun b d => ζ ![a, b, d, c] := by
  simp only [Pjet, PotMP]
  rw [LOf_addTop, TkOf_addTop_of_le (k := 2) (by norm_num), TkOf_three_addTop, PotM_add J hL,
    neg_add, sub_eq_add_neg]
  congr 2
  simp only [trOf]
  refine Finset.sum_congr rfl fun i _ => ?_
  simp only [HOfM, trOf, map_sum, sum_inner]
  refine Finset.sum_congr rfl fun k _ => ?_
  rw [inner_J_perpOf_LOf hs hq hsym, ContinuousMultilinearMap.curryLeft_apply, topTerm_apply,
    inner_J_topTerm hs hq]
  congr 1
  ext l
  fin_cases l <;> rfl

end Geometry

/-! ### The cometric `g⁻¹`: symmetric and positive definite (MM's `Ef_symm`, `Ef_pos`) -/

section Cometric

variable {E V : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
  [NormedAddCommGroup V] [InnerProductSpace ℝ V] {L : E →L[ℝ] V}

theorem cometric_eq_inner (hL : Function.Injective L) (α β : E →L[ℝ] ℝ) :
    α (ginvOf L β) = ⟪L (ginvOf L α), L (ginvOf L β)⟫ := by
  rw [inner_ginvOf hL]

theorem cometric_symm (hL : Function.Injective L) (α β : E →L[ℝ] ℝ) :
    α (ginvOf L β) = β (ginvOf L α) := by
  rw [cometric_eq_inner hL, cometric_eq_inner hL β α, real_inner_comm]

theorem cometric_pos (hL : Function.Injective L) {α : E →L[ℝ] ℝ} (hα : α ≠ 0) :
    0 < α (ginvOf L α) := by
  rw [cometric_eq_inner hL, real_inner_self_eq_norm_sq]
  have hne : L (ginvOf L α) ≠ 0 := by
    intro h0
    have h1 : ginvOf L α = 0 := hL (by rw [h0, map_zero])
    apply hα
    rw [← gOf_ginvOf hL α, h1, map_zero]
  positivity

end Cometric

end CWFlow
