import CWFlow.Extension.GraphCoeff
import CWFlow.Extension.TensorOps

/-!
# Iterated covariant derivatives of the second fundamental form, as symbols

`tSym k` is `∇^k C` of a Lagrangian graph as a function on jet space: `tSym 0 = ζ₃` and
`tSym (k+1) = D (tSym k) − Γ · tSym k`.  It is a symbol of weight `k` and equals the top jet
`ζ_{k+3}` plus a symbol of weight `k - 1` (CW Lemma 4.7).  Along the jets of a function the
recursion is the usual covariant derivative.  Finally the potential `div (J H)` is a symbol of
weight `1` whose top part is `-g^{cd} g^{ab} ζ₄(d,a,b,c)`.
-/

open scoped ContDiff

namespace CWFlow.Extension

variable {n N : ℕ}

/-! ### Part 1: definitions -/

/-- The Christoffel symbols of the graph, on jet space. -/
noncomputable def gamSym (p : P n N) : Z n 3 := gamFun (low23 p)

/-- `∇^k C` on jet space: `tSym 0 = ζ₃` and `tSym (k+1) = D(tSym k) − Γ · tSym k`. It has `k + 3`
indices, the derivative indices first. -/
noncomputable def tSym : (k : ℕ) → P n N → Z n (k + 3)
  | 0 => fun p => slotN 3 p
  | k + 1 => fun p => Dz (tSym k) p - actG (gamSym p) (tSym k p)

/-- The lower order part: `sSym k = tSym k − ζ_{k+3}`. -/
noncomputable def sSym (k : ℕ) (p : P n N) : Z n (k + 3) := tSym k p - slotN (k + 3) p

theorem tSym_zero (p : P n N) : tSym 0 p = slotN 3 p := rfl

theorem tSym_succ (k : ℕ) (p : P n N) :
    tSym (k + 1) p = Dz (tSym k) p - actG (gamSym p) (tSym k p) := rfl

/-! ### Part 0: boundedness of low jets -/

section Bdd

variable {A : Set (P n N)}

theorem bdd_slotN_low
    (hlow : ∀ m : Fin (N + 1), (m : ℕ) ≤ 3 → ∃ K : ℝ, ∀ p ∈ A, ‖p.2 m‖ ≤ K)
    (m : ℕ) (hm : m ≤ 3) : ∃ K : ℝ, ∀ p ∈ A, ‖slotN m p‖ ≤ K := by
  by_cases h : m < N + 1
  · obtain ⟨K, hK⟩ := hlow ⟨m, h⟩ hm
    exact ⟨K, fun p hp => by rw [slotN_apply m h]; exact hK p hp⟩
  · exact ⟨0, fun p _ => by rw [slotN_of_lt m (by omega)]; simp⟩

theorem bdd_low23
    (hlow : ∀ m : Fin (N + 1), (m : ℕ) ≤ 3 → ∃ K : ℝ, ∀ p ∈ A, ‖p.2 m‖ ≤ K) :
    Bornology.IsBounded (low23 '' A : Set (Z n 2 × Z n 3)) := by
  obtain ⟨K2, h2⟩ := bdd_slotN_low hlow 2 (by norm_num)
  obtain ⟨K3, h3⟩ := bdd_slotN_low hlow 3 (by norm_num)
  rw [isBounded_iff_forall_norm_le]
  refine ⟨max K2 K3, ?_⟩
  rintro _ ⟨p, hp, rfl⟩
  have e2 : (low23 p).1 = slotN 2 p := rfl
  have e3 : (low23 p).2 = slotN 3 p := rfl
  rw [Prod.norm_def, e2, e3]
  exact max_le_max (h2 p hp) (h3 p hp)

end Bdd

/-! ### Part 2: weights -/

section Weights

variable {A : Set (P n N)}

theorem IsSymbol.sub {W : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W] {w ℓ : ℕ}
    {𝔭 𝔮 : P n N → W} (h : IsSymbol A w ℓ 𝔭) (h' : IsSymbol A w ℓ 𝔮) :
    IsSymbol A w ℓ (fun p => 𝔭 p - 𝔮 p) :=
  (h.add h'.neg).congr fun p _ => sub_eq_add_neg _ _

theorem isSymbol_gamSym (hA : IsOpen A)
    (hlow : ∀ m : Fin (N + 1), (m : ℕ) ≤ 3 → ∃ K : ℝ, ∀ p ∈ A, ‖p.2 m‖ ≤ K) :
    IsSymbol A 0 0 (gamSym : P n N → Z n 3) :=
  isSymbol_gam hA (bdd_low23 hlow)

theorem isSymbol_act (hA : IsOpen A)
    (hlow : ∀ m : Fin (N + 1), (m : ℕ) ≤ 3 → ∃ K : ℝ, ∀ p ∈ A, ‖p.2 m‖ ≤ K)
    {k w ℓ : ℕ} {𝔭 : P n N → Z n k} (h : IsSymbol A w ℓ 𝔭) :
    IsSymbol A w (max 0 ℓ) (fun p => actG (gamSym p) (𝔭 p)) := by
  have := IsSymbol.bilin (actG (n := n) (k := k)) (isSymbol_gamSym hA hlow) h
  simpa using this

/-- `∇^k C` has weight `k`. -/
theorem isSymbol_tSym (hA : IsOpen A)
    (hlow : ∀ m : Fin (N + 1), (m : ℕ) ≤ 3 → ∃ K : ℝ, ∀ p ∈ A, ‖p.2 m‖ ≤ K)
    (k : ℕ) (hN : k + 3 ≤ N) : IsSymbol A k k (tSym k) := by
  revert hN
  induction k with
  | zero =>
    intro _
    exact isSymbol_slotN_low hA 3 le_rfl (bdd_slotN_low hlow 3 le_rfl)
  | succ k ih =>
    intro hN
    have ih := ih (by omega)
    have h1 : IsSymbol A (k + 1) (k + 1) (Dz (tSym k)) :=
      (ih.Dz (by omega) hlow).mono le_rfl (by omega)
    have h2 : IsSymbol A (k + 1) (k + 1) (fun p => actG (gamSym p) (tSym k p)) :=
      (isSymbol_act hA hlow ih).mono (by omega) (by omega)
    exact (h1.sub h2).congr fun p _ => rfl

theorem sSym_one (hN : 4 ≤ N) (p : P n N) : sSym 1 p = -actG (gamSym p) (slotN 3 p) := by
  have h := Dz_slotN (n := n) 3 hN p
  have e : sSym 1 p = Dz (fun q : P n N => slotN 3 q) p - actG (gamSym p) (slotN 3 p) -
      slotN 4 p := rfl
  rw [e, h]
  abel

theorem Dz_add {k : ℕ} {𝔭 𝔮 : P n N → Z n k} {p : P n N} (h1 : DifferentiableAt ℝ 𝔭 p)
    (h2 : DifferentiableAt ℝ 𝔮 p) :
    Dz (fun q => 𝔭 q + 𝔮 q) p = Dz 𝔭 p + Dz 𝔮 p := by
  unfold Dz
  have : totalD (fun q => 𝔭 q + 𝔮 q) p = totalD 𝔭 p + totalD 𝔮 p := by
    ext e : 1
    simp only [totalD_apply, fderiv_fun_add h1 h2, add_apply]
  rw [this, map_add]

/-- CW Lemma 4.7: `∇^k C = D^{k+3}φ +` terms of lower weight. -/
theorem isSymbol_sSym (hA : IsOpen A)
    (hlow : ∀ m : Fin (N + 1), (m : ℕ) ≤ 3 → ∃ K : ℝ, ∀ p ∈ A, ‖p.2 m‖ ≤ K)
    (k : ℕ) (hk : 1 ≤ k) (hN : k + 3 ≤ N) : IsSymbol A (k - 1) (k - 1) (sSym k) := by
  revert hN
  induction k, hk using Nat.le_induction with
  | base =>
    intro hN
    have h0 : IsSymbol A 0 0 (fun p : P n N => slotN 3 p) :=
      isSymbol_slotN_low hA 3 le_rfl (bdd_slotN_low hlow 3 le_rfl)
    have := (isSymbol_act hA hlow h0).neg
    exact (this.mono le_rfl (by omega)).congr fun p _ => sSym_one (by omega) p
  | succ k hk ih =>
    intro hN
    have ih := ih (by omega)
    have hDz : IsSymbol A k k (Dz (sSym k)) :=
      (ih.Dz (by omega) hlow).mono (by omega) (by omega)
    have hT := isSymbol_tSym hA hlow k (by omega)
    have h2 : IsSymbol A k k (fun p => actG (gamSym p) (tSym k p)) :=
      (isSymbol_act hA hlow hT).mono (by omega) (by omega)
    have h3 := hDz.sub h2
    show IsSymbol A k k (sSym (k + 1))
    refine h3.congr fun p hp => ?_
    have hd1 : DifferentiableAt ℝ (sSym k) p :=
      (ih.smooth.contDiffAt (ih.isOpen.mem_nhds hp)).differentiableAt (by simp)
    have hd2 : DifferentiableAt ℝ (fun q : P n N => slotN (k + 3) q) p :=
      (slotN (k + 3)).differentiableAt
    have e : tSym k = fun q : P n N => sSym k q + slotN (k + 3) q := by
      funext q
      simp [sSym]
    have hz := Dz_slotN (n := n) (k + 3) (by omega) p
    have hadd := Dz_add hd1 hd2
    rw [hz] at hadd
    have e2 : sSym (k + 1) p = Dz (tSym k) p - actG (gamSym p) (tSym k p) -
        slotN (k + 3 + 1) p := rfl
    rw [e2]
    have e3 : Dz (tSym k) p = Dz (sSym k) p + slotN (k + 3 + 1) p := by
      rw [← hadd, e]
    rw [e3]
    abel

/-! ### Part 3: realization along the jets of a function -/

theorem tSym_zero_jp (hN : 3 ≤ N) (φ : EuclideanSpace ℝ (Fin n) → ℝ)
    (x : EuclideanSpace ℝ (Fin n)) : tSym 0 (jp N φ x) = CWFlow.Interp.jet 3 φ x :=
  slotN_jp 3 hN φ x

theorem tSym_succ_jp_apply (hA : IsOpen A)
    (hlow : ∀ m : Fin (N + 1), (m : ℕ) ≤ 3 → ∃ K : ℝ, ∀ p ∈ A, ‖p.2 m‖ ≤ K)
    (k : ℕ) (hN : k + 4 ≤ N) {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ)
    {x : EuclideanSpace ℝ (Fin n)} (hx : jp N φ x ∈ A) (I : Fin (k + 3 + 1) → Fin n) :
    tSym (k + 1) (jp N φ x) I =
      fderiv ℝ (fun y => tSym k (jp N φ y) (Fin.tail I)) x (EuclideanSpace.single (I 0) 1)
      - ∑ i : Fin (k + 3), ∑ c : Fin n,
          gamSym (jp N φ x) ![I 0, Fin.tail I i, c] *
            tSym k (jp N φ x) (Function.update (Fin.tail I) i c) := by
  have hT := isSymbol_tSym hA hlow k (by omega)
  have h1 := Dz_jp_apply hT (by omega) hφ hx I
  have h2 := actG_apply (gamSym (jp N φ x)) (tSym k (jp N φ x)) I
  rw [tSym_succ]
  show (Dz (tSym k) (jp N φ x) I - actG (gamSym (jp N φ x)) (tSym k (jp N φ x)) I) = _
  rw [h1, h2]

/-! ### Part 4: the potential -/

/-- `div(JH)` of the graph, on jet space: `−g^{cd} g^{ab} (∇C)(d, a, b, c)`. -/
noncomputable def potSym (p : P n N) : ℝ :=
  -∑ a, ∑ b, ∑ c, ∑ d, ginvFun (low23 p).1 ![c, d] * ginvFun (low23 p).1 ![a, b] *
    tSym 1 p ![d, a, b, c]

theorem isSymbol_potSym (hA : IsOpen A)
    (hlow : ∀ m : Fin (N + 1), (m : ℕ) ≤ 3 → ∃ K : ℝ, ∀ p ∈ A, ‖p.2 m‖ ≤ K)
    (hN : 4 ≤ N) : IsSymbol A 1 1 (potSym : P n N → ℝ) := by
  have hG := isSymbol_ginv hA (bdd_low23 hlow)
  have hg : ∀ c d : Fin n, IsSymbol A 0 0 (fun p : P n N => ginvFun (low23 p).1 ![c, d]) :=
    fun c d => hG.clm_comp (PiLp.proj (𝕜 := ℝ) 2 (fun _ : Fin 2 → Fin n => ℝ) ![c, d])
  have hT0 := isSymbol_tSym hA hlow 1 (by omega)
  have hT : ∀ I : Fin 4 → Fin n, IsSymbol A 1 1 (fun p : P n N => tSym 1 p I) :=
    fun I => hT0.clm_comp (PiLp.proj (𝕜 := ℝ) 2 (fun _ : Fin 4 → Fin n => ℝ) I)
  have hterm : ∀ a b c d : Fin n, IsSymbol A 1 1 (fun p : P n N =>
      ginvFun (low23 p).1 ![c, d] * ginvFun (low23 p).1 ![a, b] * tSym 1 p ![d, a, b, c]) :=
    fun a b c d => (((hg c d).mul (hg a b)).mul (hT ![d, a, b, c])).mono (by omega) (by omega)
  have hs : IsSymbol A 1 1 (fun p : P n N => ∑ a, ∑ b, ∑ c, ∑ d,
      ginvFun (low23 p).1 ![c, d] * ginvFun (low23 p).1 ![a, b] * tSym 1 p ![d, a, b, c]) :=
    IsSymbol.sum hA _ fun a _ => IsSymbol.sum hA _ fun b _ => IsSymbol.sum hA _ fun c _ =>
      IsSymbol.sum hA _ fun d _ => hterm a b c d
  exact hs.neg

/-- The top order part of the potential. -/
theorem isSymbol_potSym_sub (hA : IsOpen A)
    (hlow : ∀ m : Fin (N + 1), (m : ℕ) ≤ 3 → ∃ K : ℝ, ∀ p ∈ A, ‖p.2 m‖ ≤ K)
    (hN : 4 ≤ N) :
    IsSymbol A 0 0 (fun p : P n N => potSym p +
      ∑ a, ∑ b, ∑ c, ∑ d, ginvFun (low23 p).1 ![c, d] * ginvFun (low23 p).1 ![a, b] *
        slotN 4 p ![d, a, b, c]) := by
  have hG := isSymbol_ginv hA (bdd_low23 hlow)
  have hg : ∀ c d : Fin n, IsSymbol A 0 0 (fun p : P n N => ginvFun (low23 p).1 ![c, d]) :=
    fun c d => hG.clm_comp (PiLp.proj (𝕜 := ℝ) 2 (fun _ : Fin 2 → Fin n => ℝ) ![c, d])
  have hS0 := isSymbol_sSym hA hlow 1 le_rfl (by omega)
  have hS : ∀ I : Fin 4 → Fin n, IsSymbol A 0 0 (fun p : P n N => sSym 1 p I) :=
    fun I => hS0.clm_comp (PiLp.proj (𝕜 := ℝ) 2 (fun _ : Fin 4 → Fin n => ℝ) I)
  have hterm : ∀ a b c d : Fin n, IsSymbol A 0 0 (fun p : P n N =>
      ginvFun (low23 p).1 ![c, d] * ginvFun (low23 p).1 ![a, b] * sSym 1 p ![d, a, b, c]) :=
    fun a b c d => (((hg c d).mul (hg a b)).mul (hS ![d, a, b, c])).mono (by omega) (by omega)
  have hs : IsSymbol A 0 0 (fun p : P n N => ∑ a, ∑ b, ∑ c, ∑ d,
      ginvFun (low23 p).1 ![c, d] * ginvFun (low23 p).1 ![a, b] * sSym 1 p ![d, a, b, c]) :=
    IsSymbol.sum hA _ fun a _ => IsSymbol.sum hA _ fun b _ => IsSymbol.sum hA _ fun c _ =>
      IsSymbol.sum hA _ fun d _ => hterm a b c d
  refine hs.neg.congr fun p _ => ?_
  have e : ∀ I : Fin 4 → Fin n, sSym 1 p I = tSym 1 p I - slotN 4 p I := fun I => rfl
  simp only [potSym, e, mul_sub, Finset.sum_sub_distrib]
  ring

/-- info: 'CWFlow.Extension.isSymbol_tSym' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms isSymbol_tSym

/-- info: 'CWFlow.Extension.isSymbol_sSym' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms isSymbol_sSym

/-- info: 'CWFlow.Extension.tSym_succ_jp_apply' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms tSym_succ_jp_apply

/-- info: 'CWFlow.Extension.tSym_zero_jp' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms tSym_zero_jp

/-- info: 'CWFlow.Extension.isSymbol_potSym' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms isSymbol_potSym

/-- info: 'CWFlow.Extension.isSymbol_potSym_sub' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms isSymbol_potSym_sub

end Weights

end CWFlow.Extension
