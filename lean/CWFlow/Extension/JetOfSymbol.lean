import CWFlow.Extension.SymbolLow
import CWFlow.Extension.JetSymm

/-!
# Coordinate jets of an operator given by a symbol

If `f = fun y => 𝔭 (jp N φ y)` for a scalar symbol `𝔭`, then the `m`-jet of `f` is again given by a function
on jet space, `jetSym 𝔭 m` (the `m`-fold iterated total derivative `Dz`), of weight `w + m`
(`isSymbol_jetSym`, `jet_symbol_jp`).  For the potential `potSym` the top order part of
`jetSym potSym m` is `topPot m` (`isSymbol_jetSym_potSym_sub`).
-/

set_option linter.unusedSectionVars false

open scoped ContDiff
open CWFlow.Interp

namespace CWFlow.Extension

variable {n N : ℕ}

/-! ### Part 1: the iterated total derivative -/

/-- The `m`-jet, on jet space, of the operator defined by a scalar symbol. -/
noncomputable def jetSym (𝔭 : P n N → ℝ) : (m : ℕ) → P n N → Z n m
  | 0 => fun p => WithLp.toLp 2 (fun _ => 𝔭 p)
  | m + 1 => fun p => Dz (jetSym 𝔭 m) p

/-- The constant tensor with a given value, as a continuous linear map. -/
noncomputable def constZ (n : ℕ) : ℝ →L[ℝ] Z n 0 :=
  LinearMap.toContinuousLinearMap
    { toFun := fun c => WithLp.toLp 2 (fun _ => c)
      map_add' := fun c c' => by ext I; simp
      map_smul' := fun a c => by ext I; simp }

theorem isSymbol_jetSym {A : Set (P n N)} {w ℓ : ℕ} {𝔭 : P n N → ℝ} (h : IsSymbol A w ℓ 𝔭)
    (hlow : ∀ m : Fin (N + 1), (m : ℕ) ≤ 3 → ∃ K : ℝ, ∀ p ∈ A, ‖p.2 m‖ ≤ K)
    (m : ℕ) (hN : w + m + 3 ≤ N) : IsSymbol A (w + m) (max ℓ w + m) (jetSym 𝔭 m) := by
  induction m with
  | zero =>
    have := (h.clm_comp (constZ n)).mono (le_refl w) (le_max_left ℓ w)
    exact this
  | succ m ih =>
    have ih := ih (by omega)
    have h2 := ih.Dz (by omega) hlow
    have e : max (max ℓ w + m) (w + m) = max ℓ w + m := by omega
    rw [e] at h2
    have e2 : w + (m + 1) = w + m + 1 := by omega
    have e3 : max ℓ w + (m + 1) = max ℓ w + m + 1 := by omega
    rw [e2, e3]
    exact h2

/-! ### Part 2: realization -/

/-- One partial derivative of a coordinate of the jet, using only smoothness near `x`. -/
theorem jet_cons_eq_fderiv_of_contDiffAt {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
    [FiniteDimensional ℝ V] {f : EuclideanSpace ℝ (Fin n) → V} {x : EuclideanSpace ℝ (Fin n)}
    (hf : ContDiffAt ℝ ∞ f x) (m : ℕ) (a : Fin n) (J : Fin m → Fin n) :
    jet (m + 1) f x (Fin.cons a J) =
      fderiv ℝ (fun y => jet m f y J) x (EuclideanSpace.single a 1) := by
  have hd : DifferentiableAt ℝ (iteratedFDeriv ℝ m f) x :=
    hf.differentiableAt_iteratedFDeriv (by exact_mod_cast ENat.natCast_lt_top m)
  have e : (fun y => jet m f y J) =
      fun y => iteratedFDeriv ℝ m f y (fun k => EuclideanSpace.single (J k) (1 : ℝ)) := rfl
  rw [e, fderiv_continuousMultilinear_apply_const_apply hd, jet_apply,
    iteratedFDeriv_succ_apply_left]
  rfl

theorem jet_symbol_jp {A : Set (P n N)} {w ℓ : ℕ} {𝔭 : P n N → ℝ} (h : IsSymbol A w ℓ 𝔭)
    (hlow : ∀ m : Fin (N + 1), (m : ℕ) ≤ 3 → ∃ K : ℝ, ∀ p ∈ A, ‖p.2 m‖ ≤ K)
    (m : ℕ) (hN : w + m + 3 ≤ N) {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ)
    {x : EuclideanSpace ℝ (Fin n)} (hx : jp N φ x ∈ A) :
    CWFlow.Interp.jet m (fun y => 𝔭 (jp N φ y)) x = jetSym 𝔭 m (jp N φ x) := by
  induction m generalizing x with
  | zero =>
    ext I
    simp [jet_apply, iteratedFDeriv_zero_apply, jetSym]
  | succ m ih =>
    have hU : IsOpen {y | jp N φ y ∈ A} := isOpen_jp_preimage h hφ
    have hf : ContDiffAt ℝ ∞ (fun y => 𝔭 (jp N φ y)) x :=
      (contDiffOn_symbol_jp h hφ).contDiffAt (hU.mem_nhds hx)
    have hs := isSymbol_jetSym h hlow m (by omega)
    ext I
    obtain ⟨a, J, rfl⟩ : ∃ (a : Fin n) (J : Fin m → Fin n), I = Fin.cons a J :=
      ⟨I 0, Fin.tail I, (Fin.cons_self_tail I).symm⟩
    rw [jet_cons_eq_fderiv_of_contDiffAt hf]
    have hev : (fun y => jet m (fun y => 𝔭 (jp N φ y)) y J) =ᶠ[nhds x]
        (fun y => jetSym 𝔭 m (jp N φ y) J) := by
      filter_upwards [hU.mem_nhds hx] with y hy
      rw [ih (by omega) hy]
    rw [hev.fderiv_eq]
    have := Dz_jp_apply hs (by omega) hφ hx (Fin.cons a J)
    simp only [Fin.tail_cons, Fin.cons_zero] at this
    exact this.symm

/-! ### Part 3: the top order part of the jets of the potential -/

section Top

variable {A : Set (P n N)}

/-- A tensor-valued function is a symbol if all its coordinates are. -/
theorem isSymbol_of_coord {k : ℕ} {w ℓ : ℕ} {𝔭 : P n N → Z n k} (hA : IsOpen A)
    (h : ∀ I, IsSymbol A w ℓ (fun p => 𝔭 p I)) : IsSymbol A w ℓ 𝔭 := by
  have hs : IsSymbol A w ℓ (fun p => ∑ I : Fin k → Fin n,
      ((ContinuousLinearMap.id ℝ ℝ).smulRight
        (WithLp.toLp 2 (Pi.single I (1 : ℝ)) : Z n k)) (𝔭 p I)) :=
    IsSymbol.sum hA _ fun I _ => (h I).clm_comp _
  refine hs.congr fun p _ => ?_
  ext J
  simp [Pi.single_apply]

/-- `−Σ g^{cd} g^{ab} ζ_{m+4}(I, d, a, b, c)`: the top order part of the `m`-jet of the potential. -/
noncomputable def topPot (m : ℕ) (p : P n N) : Z n m :=
  WithLp.toLp 2 (fun I => -∑ a, ∑ b, ∑ c, ∑ d,
    ginvFun (low23 p).1 ![c, d] * ginvFun (low23 p).1 ![a, b] *
      slotN (m + 4) p (Fin.append I ![d, a, b, c]))

/-- The product of two coordinates of the inverse metric. -/
@[reducible] noncomputable def gg (a b c d : Fin n) (q : P n N) : ℝ :=
  ginvFun (low23 q).1 ![c, d] * ginvFun (low23 q).1 ![a, b]

theorem contDiff_gg (a b c d : Fin n) : ContDiff ℝ ∞ (gg a b c d : P n N → ℝ) := by
  have h : ∀ c d : Fin n, ContDiff ℝ ∞ (fun q : P n N => ginvFun (low23 q).1 ![c, d]) := fun c d =>
    (PiLp.proj (𝕜 := ℝ) 2 (fun _ : Fin 2 → Fin n => ℝ) ![c, d]).contDiff.comp
      (contDiff_ginvFun.comp (contDiff_fst.comp low23.contDiff))
  exact (h c d).mul (h a b)

theorem contDiff_slotN_apply (j : ℕ) (K : Fin j → Fin n) :
    ContDiff ℝ ∞ (fun q : P n N => slotN j q K) :=
  (PiLp.proj (𝕜 := ℝ) 2 (fun _ : Fin j → Fin n => ℝ) K).contDiff.comp (slotN j).contDiff

theorem totalD_apply_coord {k : ℕ} {𝔭 : P n N → Z n k} {p : P n N} (hd : DifferentiableAt ℝ 𝔭 p)
    (v : EuclideanSpace ℝ (Fin n)) (I : Fin k → Fin n) :
    totalD 𝔭 p v I = totalD (fun q => 𝔭 q I) p v := by
  have h := (PiLp.proj (𝕜 := ℝ) 2 (fun _ : Fin k → Fin n => ℝ) I).hasFDerivAt.comp p hd.hasFDerivAt
  have h2 : fderiv ℝ (fun y => 𝔭 y I) p =
      (PiLp.proj (𝕜 := ℝ) 2 (fun _ : Fin k → Fin n => ℝ) I).comp (fderiv ℝ 𝔭 p) := h.fderiv
  rw [totalD_apply, totalD_apply, h2]
  rfl

/-- The total derivative of a coordinate of the jet slot is a coordinate of the next jet slot. -/
theorem totalD_slotN_apply (j : ℕ) (hj : j + 1 ≤ N) (p : P n N) (e : Fin n) (K : Fin j → Fin n) :
    totalD (fun q : P n N => slotN j q K) p (EuclideanSpace.single e 1) =
      slotN (j + 1) p (Fin.cons e K) := by
  have h := Dz_slotN (n := n) j hj p
  have h2 := congrArg (fun z : Z n (j + 1) => z (Fin.cons e K)) h
  simp only [Dz, curryZ_apply, Fin.tail_cons, Fin.cons_zero] at h2
  rw [← h2]
  exact (totalD_apply_coord (slotN j).differentiableAt _ K).symm

/-- The cast identity for the slot index. -/
theorem slotN_append_cons (m : ℕ) (p : P n N) (e : Fin n) (I : Fin m → Fin n) (L : Fin 4 → Fin n) :
    slotN (m + 1 + 4) p (Fin.append (Fin.cons e I) L) =
      slotN (m + 4 + 1) p (Fin.cons e (Fin.append I L)) := by
  rw [Fin.append_cons]
  rfl

/-- The derivative of the top order part, minus the top order part of the next jet. -/
noncomputable def rPot (m : ℕ) (p : P n N) : Z n (m + 1) :=
  WithLp.toLp 2 (fun J => -∑ a, ∑ b, ∑ c, ∑ d,
    totalD (gg a b c d) p (EuclideanSpace.single (J 0) 1) *
      slotN (m + 4) p (Fin.append (Fin.tail J) ![d, a, b, c]))

theorem totalD_sum_apply {W : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W] {ι : Type*}
    (S : Finset ι) {𝔭 : ι → P n N → W} {p : P n N}
    (h : ∀ i ∈ S, DifferentiableAt ℝ (𝔭 i) p) (v : EuclideanSpace ℝ (Fin n)) :
    totalD (fun q => ∑ i ∈ S, 𝔭 i q) p v = ∑ i ∈ S, totalD (𝔭 i) p v := by
  rw [totalD_sum S h]
  simp

theorem topPot_apply (m : ℕ) (p : P n N) (I : Fin m → Fin n) :
    topPot m p I = -∑ a, ∑ b, ∑ c, ∑ d,
      gg a b c d p * slotN (m + 4) p (Fin.append I ![d, a, b, c]) := rfl

theorem rPot_apply (m : ℕ) (p : P n N) (J : Fin (m + 1) → Fin n) :
    rPot m p J = -∑ a, ∑ b, ∑ c, ∑ d,
      totalD (gg a b c d) p (EuclideanSpace.single (J 0) 1) *
        slotN (m + 4) p (Fin.append (Fin.tail J) ![d, a, b, c]) := rfl

theorem Dz_topPot_sub (m : ℕ) (hN : m + 5 ≤ N) (p : P n N)
    (hd : DifferentiableAt ℝ (topPot m) p) :
    Dz (topPot m) p - topPot (m + 1) p = rPot m p := by
  ext J
  obtain ⟨e, I, rfl⟩ : ∃ (e : Fin n) (I : Fin m → Fin n), J = Fin.cons e I :=
    ⟨J 0, Fin.tail J, (Fin.cons_self_tail J).symm⟩
  have hD : Dz (topPot m) p (Fin.cons e I) =
      totalD (fun q => topPot m q I) p (EuclideanSpace.single e 1) := by
    rw [Dz, curryZ_apply]
    simp only [Fin.tail_cons, Fin.cons_zero]
    exact totalD_apply_coord hd _ I
  have hf : (fun q => topPot m q I) = fun q : P n N => -∑ a, ∑ b, ∑ c, ∑ d,
      gg a b c d q * slotN (m + 4) q (Fin.append I ![d, a, b, c]) := funext fun q => topPot_apply m q I
  have hdg : ∀ a b c d : Fin n, DifferentiableAt ℝ (gg a b c d : P n N → ℝ) p :=
    fun a b c d => ((contDiff_gg a b c d).differentiable (by simp)) p
  have hds : ∀ K : Fin (m + 4) → Fin n, DifferentiableAt ℝ (fun q : P n N => slotN (m + 4) q K) p :=
    fun K => ((contDiff_slotN_apply (m + 4) K).differentiable (by simp)) p
  have hT : totalD (fun q => topPot m q I) p (EuclideanSpace.single e 1) =
      -∑ a, ∑ b, ∑ c, ∑ d, (totalD (gg a b c d) p (EuclideanSpace.single e 1) *
          slotN (m + 4) p (Fin.append I ![d, a, b, c]) +
        gg a b c d p * slotN (m + 4 + 1) p (Fin.cons e (Fin.append I ![d, a, b, c]))) := by
    rw [hf, totalD_neg, neg_apply]
    congr 1
    rw [totalD_sum_apply _ (fun a _ => by fun_prop)]
    refine Finset.sum_congr rfl fun a _ => ?_
    rw [totalD_sum_apply _ (fun b _ => by fun_prop)]
    refine Finset.sum_congr rfl fun b _ => ?_
    rw [totalD_sum_apply _ (fun c _ => by fun_prop)]
    refine Finset.sum_congr rfl fun c _ => ?_
    rw [totalD_sum_apply _ (fun d _ => by fun_prop)]
    refine Finset.sum_congr rfl fun d _ => ?_
    rw [totalD_mul (hdg a b c d) (hds _), totalD_slotN_apply (m + 4) (by omega)]
  have hR : topPot (m + 1) p (Fin.cons e I) =
      -∑ a, ∑ b, ∑ c, ∑ d, gg a b c d p *
        slotN (m + 4 + 1) p (Fin.cons e (Fin.append I ![d, a, b, c])) := by
    rw [topPot_apply]
    simp only [slotN_append_cons]
  rw [PiLp.sub_apply, hD, hT, hR, rPot_apply]
  simp only [Fin.tail_cons, Fin.cons_zero, Finset.sum_add_distrib]
  ring

theorem isSymbol_slotN_apply (hA : IsOpen A) (m : ℕ) (K : Fin (m + 4) → Fin n) :
    IsSymbol A (m + 1) (m + 1) (fun p : P n N => slotN (m + 4) p K) :=
  ((isSymbol_slotN hA (m + 4) (by omega)).clm_comp
    (PiLp.proj (𝕜 := ℝ) 2 (fun _ : Fin (m + 4) → Fin n => ℝ) K)).mono (by omega) (by omega)

theorem isSymbol_gg (hA : IsOpen A)
    (hlow : ∀ m : Fin (N + 1), (m : ℕ) ≤ 3 → ∃ K : ℝ, ∀ p ∈ A, ‖p.2 m‖ ≤ K) (a b c d : Fin n) :
    IsSymbol A 0 0 (gg a b c d : P n N → ℝ) := by
  have hG := isSymbol_ginv hA (bdd_low23 hlow)
  have hg : ∀ c d : Fin n, IsSymbol A 0 0 (fun p : P n N => ginvFun (low23 p).1 ![c, d]) :=
    fun c d => hG.clm_comp (PiLp.proj (𝕜 := ℝ) 2 (fun _ : Fin 2 → Fin n => ℝ) ![c, d])
  exact ((hg c d).mul (hg a b)).mono (by omega) (by omega)

theorem isSymbol_topPot (hA : IsOpen A)
    (hlow : ∀ m : Fin (N + 1), (m : ℕ) ≤ 3 → ∃ K : ℝ, ∀ p ∈ A, ‖p.2 m‖ ≤ K) (m : ℕ) :
    IsSymbol A (m + 1) (m + 1) (topPot m : P n N → Z n m) := by
  refine isSymbol_of_coord hA fun I => ?_
  have hs : IsSymbol A (m + 1) (m + 1) (fun p : P n N => ∑ a, ∑ b, ∑ c, ∑ d,
      gg a b c d p * slotN (m + 4) p (Fin.append I ![d, a, b, c])) :=
    IsSymbol.sum hA _ fun a _ => IsSymbol.sum hA _ fun b _ => IsSymbol.sum hA _ fun c _ =>
      IsSymbol.sum hA _ fun d _ =>
        ((isSymbol_gg hA hlow a b c d).mul (isSymbol_slotN_apply hA m _)).mono (by omega)
          (by omega)
  exact hs.neg.congr fun p _ => topPot_apply m p I

theorem isSymbol_totalD_gg (hA : IsOpen A)
    (hlow : ∀ m : Fin (N + 1), (m : ℕ) ≤ 3 → ∃ K : ℝ, ∀ p ∈ A, ‖p.2 m‖ ≤ K) (a b c d e : Fin n) :
    IsSymbol A 0 0 (fun p : P n N => totalD (gg a b c d) p (EuclideanSpace.single e 1)) := by
  have hc : ContDiff ℝ ∞ (fun u : Z n 2 => ginvFun u ![c, d] * ginvFun u ![a, b]) := by
    have h : ∀ c d : Fin n, ContDiff ℝ ∞ (fun u : Z n 2 => ginvFun u ![c, d]) := fun c d =>
      (PiLp.proj (𝕜 := ℝ) 2 (fun _ : Fin 2 → Fin n => ℝ) ![c, d]).contDiff.comp contDiff_ginvFun
    exact (h c d).mul (h a b)
  exact (isSymbol_totalD_comp_low2 hA (bdd_low23 hlow) hc).clm_comp
    (ContinuousLinearMap.apply ℝ ℝ (EuclideanSpace.single e (1 : ℝ)))

theorem isSymbol_rPot (hA : IsOpen A)
    (hlow : ∀ m : Fin (N + 1), (m : ℕ) ≤ 3 → ∃ K : ℝ, ∀ p ∈ A, ‖p.2 m‖ ≤ K) (m : ℕ) :
    IsSymbol A (m + 1) (m + 1) (rPot m : P n N → Z n (m + 1)) := by
  refine isSymbol_of_coord hA fun J => ?_
  have hs : IsSymbol A (m + 1) (m + 1) (fun p : P n N => ∑ a, ∑ b, ∑ c, ∑ d,
      totalD (gg a b c d) p (EuclideanSpace.single (J 0) 1) *
        slotN (m + 4) p (Fin.append (Fin.tail J) ![d, a, b, c])) :=
    IsSymbol.sum hA _ fun a _ => IsSymbol.sum hA _ fun b _ => IsSymbol.sum hA _ fun c _ =>
      IsSymbol.sum hA _ fun d _ =>
        ((isSymbol_totalD_gg hA hlow a b c d (J 0)).mul (isSymbol_slotN_apply hA m _)).mono
          (by omega) (by omega)
  exact hs.neg.congr fun p _ => rPot_apply m p J

theorem isSymbol_jetSym_potSym_sub {A : Set (P n N)} (hA : IsOpen A)
    (hlow : ∀ m : Fin (N + 1), (m : ℕ) ≤ 3 → ∃ K : ℝ, ∀ p ∈ A, ‖p.2 m‖ ≤ K)
    (m : ℕ) (hN : m + 4 ≤ N) :
    IsSymbol A m m (fun p => jetSym potSym m p - topPot m p) := by
  induction m with
  | zero =>
    have h0 := (isSymbol_potSym_sub hA hlow (by omega)).clm_comp (constZ n)
    refine h0.congr fun p _ => ?_
    ext I
    have hI : ∀ (I : Fin 0 → Fin n) (v : Fin 4 → Fin n),
        (Fin.append I v : Fin (0 + 4) → Fin n) = v := by
      intro I v
      rw [Fin.append_left_nil I v rfl]
      rfl
    rw [PiLp.sub_apply, topPot_apply]
    simp only [hI]
    show potSym p - -∑ a, ∑ b, ∑ c, ∑ d, ginvFun (low23 p).1 ![c, d] * ginvFun (low23 p).1 ![a, b] *
      slotN 4 p ![d, a, b, c] = _
    rw [sub_neg_eq_add]
    rfl
  | succ m ih =>
    have ih := ih (by omega)
    have h1 : IsSymbol A (m + 1) (m + 1) (Dz (fun p => jetSym potSym m p - topPot m p)) :=
      (ih.Dz (by omega) hlow).mono le_rfl (by omega)
    have h2 := isSymbol_rPot hA hlow m
    refine (h1.add h2).congr fun p hp => ?_
    have hj : DifferentiableAt ℝ (fun q => jetSym potSym m q - topPot m q) p :=
      (ih.smooth.contDiffAt (hA.mem_nhds hp)).differentiableAt (by simp)
    have ht : DifferentiableAt ℝ (topPot m) p :=
      ((isSymbol_topPot hA hlow m).smooth.contDiffAt (hA.mem_nhds hp)).differentiableAt (by simp)
    have hadd := Dz_add hj ht
    have e : (fun q : P n N => (jetSym (potSym : P n N → ℝ) m q - topPot m q) + topPot m q) =
        jetSym (potSym : P n N → ℝ) m :=
      funext fun q => sub_add_cancel _ _
    rw [e] at hadd
    rw [← Dz_topPot_sub m (by omega) p ht]
    show Dz (jetSym potSym m) p - topPot (m + 1) p = _
    rw [hadd]
    abel

/-- info: 'CWFlow.Extension.isSymbol_jetSym' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms isSymbol_jetSym

/-- info: 'CWFlow.Extension.jet_symbol_jp' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms jet_symbol_jp

/-- info: 'CWFlow.Extension.isSymbol_jetSym_potSym_sub' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms isSymbol_jetSym_potSym_sub

end Top

end CWFlow.Extension
