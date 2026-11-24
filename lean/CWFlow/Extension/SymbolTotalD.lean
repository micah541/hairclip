import CWFlow.Extension.SymbolProd
import CWFlow.Extension.SymbolDeriv

/-!
# Further closure properties of symbol classes

Congruence on the domain, finite sums, composition of globally smooth functions with low-order
projections, the slot derivative `Dsl`, and closure under the total derivative `totalD`.
-/

open scoped ContDiff

namespace CWFlow.Extension

variable {n N : ℕ}

/-! ### Part 1: congruence and finite sums -/

theorem IsSymbol.congr {W : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W] {A : Set (P n N)}
    {w ℓ : ℕ} {𝔭 𝔮 : P n N → W} (h : IsSymbol A w ℓ 𝔭) (heq : ∀ p ∈ A, 𝔮 p = 𝔭 p) :
    IsSymbol A w ℓ 𝔮 := by
  refine ⟨h.isOpen, h.smooth.congr heq, fun r => ?_⟩
  obtain ⟨C, hC⟩ := h.bound r
  refine ⟨C, fun p hp δ => ?_⟩
  have hev : 𝔭 =ᶠ[nhds p] 𝔮 :=
    Filter.eventuallyEq_of_mem (h.isOpen.mem_nhds hp) fun q hq => (heq q hq).symm
  have : iteratedFDeriv ℝ r 𝔭 p = iteratedFDeriv ℝ r 𝔮 p :=
    (Filter.EventuallyEq.iteratedFDeriv ℝ hev r).eq_of_nhds
  rw [← this]
  exact hC p hp δ

theorem IsSymbol.sum {W : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W] {A : Set (P n N)}
    (hA : IsOpen A) {w ℓ : ℕ} {ι : Type*} (S : Finset ι) {𝔭 : ι → P n N → W}
    (h : ∀ i ∈ S, IsSymbol A w ℓ (𝔭 i)) :
    IsSymbol A w ℓ (fun p => ∑ i ∈ S, 𝔭 i p) := by
  classical
  induction S using Finset.induction_on with
  | empty =>
    simp only [Finset.sum_empty]
    exact (isSymbol_const hA (0 : W)).mono (Nat.zero_le _) (Nat.zero_le _)
  | insert a s ha ih =>
    have h1 := (h a (Finset.mem_insert_self a s)).add (ih fun i hi => h i (Finset.mem_insert_of_mem hi))
    refine h1.congr fun p _ => ?_
    simp [Finset.sum_insert ha]

/-! ### Part 2: globally smooth functions of low slots -/

theorem isSymbol_comp_low {W Q : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W]
    [NormedAddCommGroup Q] [NormedSpace ℝ Q] [FiniteDimensional ℝ Q]
    {A : Set (P n N)} (hA : IsOpen A) (π : P n N →L[ℝ] Q)
    (hπ : ∀ δ : P n N, (δ.1 = 0 ∧ ∀ m : Fin (N + 1), (m : ℕ) ≤ 3 → δ.2 m = 0) → π δ = 0)
    (hbdd : Bornology.IsBounded (π '' A)) {c : Q → W} (hc : ContDiff ℝ ∞ c) :
    IsSymbol A 0 0 (fun p => c (π p)) := by
  have hcπ : ContDiff ℝ ∞ (fun p : P n N => c (π p)) := hc.comp π.contDiff
  have key : ∀ (r : ℕ) (p : P n N), iteratedFDeriv ℝ r (fun p => c (π p)) p =
      (iteratedFDeriv ℝ r c (π p)).compContinuousLinearMap fun _ => π := fun r p =>
    π.iteratedFDeriv_comp_right hc p (by exact_mod_cast le_top)
  refine isSymbol_of_low hA hcπ.contDiffOn (fun r => ?_) (fun r p _ δ hδ => ?_)
  · have hcont : Continuous (iteratedFDeriv ℝ r c) :=
      hc.continuous_iteratedFDeriv (by exact_mod_cast le_top)
    obtain ⟨M, hM⟩ := hbdd.isCompact_closure.exists_bound_of_continuousOn hcont.continuousOn
    refine ⟨M * ‖π‖ ^ r, fun p hp => ?_⟩
    have hmem : π p ∈ closure (π '' A) := subset_closure ⟨p, hp, rfl⟩
    rw [key]
    refine (ContinuousMultilinearMap.norm_compContinuousLinearMap_le _ _).trans ?_
    simp only [Finset.prod_const, Finset.card_univ, Fintype.card_fin]
    exact mul_le_mul_of_nonneg_right (hM _ hmem) (by positivity)
  · obtain ⟨i, hi1, hi2⟩ := hδ
    rw [key]
    simp only [ContinuousMultilinearMap.compContinuousLinearMap_apply]
    exact ContinuousMultilinearMap.map_coord_zero _ i (hπ _ ⟨hi1, hi2⟩)

/-! ### Part 3: the slot derivative of a symbol -/

theorem IsSymbol.Dsl {W : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W] {A : Set (P n N)}
    {w ℓ : ℕ} {𝔭 : P n N → W} (h : IsSymbol A w ℓ 𝔭) (s : Slot N) (hs : wt s ≤ w) :
    IsSymbol A (w - wt s) ℓ (CWFlow.Extension.Dsl s 𝔭) := by
  refine ⟨h.isOpen, contDiffOn_Dsl h.isOpen h.smooth s, fun r => ?_⟩
  obtain ⟨C, hC⟩ := bu_deriv h.isOpen h.smooth (h.bu (r + 1)) s r le_rfl
  have e : ((w - wt s : ℕ) : ℤ) = (w : ℤ) - wt s := by omega
  refine ⟨C, fun p hp δ => ?_⟩
  rw [e]
  exact hC p hp δ

theorem IsSymbol.Dsl_eq_zero {W : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W]
    {A : Set (P n N)} {w ℓ : ℕ} {𝔭 : P n N → W} (h : IsSymbol A w ℓ 𝔭) (s : Slot N)
    (hs : w < wt s) : ∀ p ∈ A, CWFlow.Extension.Dsl s 𝔭 p = 0 := by
  obtain ⟨C, hC⟩ := bu_deriv h.isOpen h.smooth (h.bu 1) s 0 (Nat.zero_le _)
  intro p hp
  have h1 := hC p hp (fun i => i.elim0)
  rw [iteratedFDeriv_zero_apply, bdS_zero, OmZ_of_neg (by omega), mul_zero] at h1
  exact norm_le_zero_iff.1 h1

/-! ### Part 4: closure under the total derivative -/

/-- Contraction as a continuous bilinear map. -/
noncomputable def contractCLM (n m : ℕ) :
    EuclideanSpace ℝ (Fin n) →L[ℝ] (Z n (m + 1) →L[ℝ] Z n m) :=
  LinearMap.toContinuousLinearMap
    { toFun := fun e => contract e
      map_add' := fun e e' => ContinuousLinearMap.ext fun z => contract_add e e' z
      map_smul' := fun c e => ContinuousLinearMap.ext fun z => contract_smul c e z }

/-- Contraction with the jet argument first. -/
noncomputable def contractFlip (n m : ℕ) :
    Z n (m + 1) →L[ℝ] (EuclideanSpace ℝ (Fin n) →L[ℝ] Z n m) :=
  (contractCLM n m).flip

theorem contractFlip_apply (n m : ℕ) (z : Z n (m + 1)) (e : EuclideanSpace ℝ (Fin n)) :
    contractFlip n m z e = contract e z := rfl

/-- Inclusion of the `m`-jet slot into `P n N`. -/
noncomputable def injSlot (n N : ℕ) (m : Fin (N + 1)) : Z n m →L[ℝ] P n N :=
  (ContinuousLinearMap.inr ℝ (EuclideanSpace ℝ (Fin n)) (J n N)).comp
    (ContinuousLinearMap.single ℝ (fun m : Fin (N + 1) => Z n m) m)

theorem injSlot_apply (n N : ℕ) (m : Fin (N + 1)) (x : Z n m) :
    injSlot n N m x = (0, Pi.single m x) := rfl

/-- The bilinear map `(L, z) ↦ (e ↦ L (0, Pi.single m (contract e z)))`. -/
noncomputable def Bm (W : Type*) [NormedAddCommGroup W] [NormedSpace ℝ W] (n N : ℕ)
    (m : Fin (N + 1)) :
    (P n N →L[ℝ] W) →L[ℝ] Z n ((m : ℕ) + 1) →L[ℝ] (EuclideanSpace ℝ (Fin n) →L[ℝ] W) :=
  ((ContinuousLinearMap.compL ℝ (EuclideanSpace ℝ (Fin n)) (P n N) W).flip.comp
    (((ContinuousLinearMap.compL ℝ (EuclideanSpace ℝ (Fin n)) (Z n m) (P n N)) (injSlot n N m)).comp
      (contractFlip n (m : ℕ)))).flip

theorem Bm_apply (W : Type*) [NormedAddCommGroup W] [NormedSpace ℝ W] (n N : ℕ) (m : Fin (N + 1))
    (L : P n N →L[ℝ] W) (z : Z n ((m : ℕ) + 1)) (e : EuclideanSpace ℝ (Fin n)) :
    Bm W n N m L z e = L (injSlot n N m (contract e z)) := rfl

section Total

variable {W : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W]

/-- The `m`-th term of the total derivative. -/
noncomputable def termM (𝔭 : P n N → W) (w : ℕ) (m : Fin (N + 1)) (p : P n N) :
    EuclideanSpace ℝ (Fin n) →L[ℝ] W :=
  if h : (m : ℕ) + 1 < N + 1 ∧ wt (some m : Slot N) ≤ w then
    Bm W n N m (Dsl (some m) 𝔭 p) (p.2 ⟨(m : ℕ) + 1, h.1⟩)
  else 0

theorem isSymbol_termM {A : Set (P n N)} {w ℓ : ℕ} {𝔭 : P n N → W} (h : IsSymbol A w ℓ 𝔭)
    (hlow : ∀ m : Fin (N + 1), (m : ℕ) ≤ 3 → ∃ K : ℝ, ∀ p ∈ A, ‖p.2 m‖ ≤ K)
    (m : Fin (N + 1)) :
    IsSymbol A (w + 1) (max ℓ w + 1) (termM 𝔭 w m) := by
  by_cases hh : (m : ℕ) + 1 < N + 1 ∧ wt (some m : Slot N) ≤ w
  · obtain ⟨hlt, hs⟩ := hh
    have hD := h.Dsl (some m) hs
    have hwt : wt (some m : Slot N) = (m : ℕ) - 3 := rfl
    have hs' : (m : ℕ) - 3 ≤ w := hwt ▸ hs
    have hfin : IsSymbol A (w + 1) (max ℓ w + 1) (fun p =>
        Bm W n N m (Dsl (some m) 𝔭 p) (slotProj (⟨(m : ℕ) + 1, hlt⟩ : Fin (N + 1)) p)) := by
      by_cases h4 : 4 ≤ (m : ℕ) + 1
      · have h2 : IsSymbol A ((m : ℕ) + 1 - 3) ((m : ℕ) + 1 - 3)
            (slotProj (⟨(m : ℕ) + 1, hlt⟩ : Fin (N + 1))) :=
          isSymbol_slotProj h.isOpen (⟨(m : ℕ) + 1, hlt⟩ : Fin (N + 1)) h4
        exact (IsSymbol.bilin (Bm W n N m) hD h2).mono (by rw [hwt]; omega) (by omega)
      · have h2 : IsSymbol A 0 0 (slotProj (⟨(m : ℕ) + 1, hlt⟩ : Fin (N + 1))) :=
          isSymbol_slotProj_low h.isOpen _ (by show (m : ℕ) + 1 ≤ 3; omega)
            (hlow ⟨(m : ℕ) + 1, hlt⟩ (by show (m : ℕ) + 1 ≤ 3; omega))
        exact (IsSymbol.bilin (Bm W n N m) hD h2).mono (by rw [hwt]; omega) (by omega)
    refine hfin.congr fun p _ => ?_
    unfold termM
    rw [dif_pos ⟨hlt, hs⟩]
    rfl
  · exact ((isSymbol_const h.isOpen (0 : EuclideanSpace ℝ (Fin n) →L[ℝ] W)).mono
      (Nat.zero_le _) (Nat.zero_le _)).congr fun p _ => by
        unfold termM
        rw [dif_neg hh]

theorem totalD_eq_sum {A : Set (P n N)} {w ℓ : ℕ} {𝔭 : P n N → W} (h : IsSymbol A w ℓ 𝔭)
    {p : P n N} (hp : p ∈ A) :
    totalD 𝔭 p = (Dsl none 𝔭 p).comp (ContinuousLinearMap.inl ℝ (EuclideanSpace ℝ (Fin n)) (J n N))
      + ∑ m : Fin (N + 1), termM 𝔭 w m p := by
  ext e
  rw [totalD_apply]
  simp only [_root_.add_apply, ContinuousLinearMap.comp_apply, _root_.sum_apply]
  have h1 : fderiv ℝ 𝔭 p (e, shift e p.2) = ∑ s : Slot N, Dsl s 𝔭 p (e, shift e p.2) := by
    simp only [Dsl, ContinuousLinearMap.comp_apply]
    rw [← map_sum, iotaP_sum]
  rw [h1, Fintype.sum_option]
  congr 1
  refine Finset.sum_congr rfl fun m _ => ?_
  unfold termM
  by_cases hh : (m : ℕ) + 1 < N + 1 ∧ wt (some m : Slot N) ≤ w
  · rw [dif_pos hh, Bm_apply]
    have hsh : shift e p.2 m = contract e (p.2 ⟨(m : ℕ) + 1, hh.1⟩) := by
      unfold shift
      rw [dif_pos hh.1]
    simp [Dsl, iotaP_some_apply, injSlot_apply, hsh]
  · rw [dif_neg hh, _root_.zero_apply]
    by_cases hlt : (m : ℕ) + 1 < N + 1
    · have hw : w < wt (some m : Slot N) := not_le.1 fun hs => hh ⟨hlt, hs⟩
      rw [h.Dsl_eq_zero (some m) hw p hp]
      rfl
    · have hsh : shift e p.2 m = 0 := by
        unfold shift
        rw [dif_neg hlt]
      simp [Dsl, iotaP_some_apply, hsh]
      exact map_zero (fderiv ℝ 𝔭 p)

theorem IsSymbol.totalD {A : Set (P n N)} {w ℓ : ℕ} {𝔭 : P n N → W} (h : IsSymbol A w ℓ 𝔭)
    (hN : w + 4 ≤ N)
    (hlow : ∀ m : Fin (N + 1), (m : ℕ) ≤ 3 → ∃ K : ℝ, ∀ p ∈ A, ‖p.2 m‖ ≤ K) :
    IsSymbol A (w + 1) (max ℓ w + 1) (CWFlow.Extension.totalD 𝔭) := by
  have h0 : IsSymbol A (w + 1) (max ℓ w + 1) (fun p => (CWFlow.Extension.Dsl none 𝔭 p).comp
      (ContinuousLinearMap.inl ℝ (EuclideanSpace ℝ (Fin n)) (J n N))) := by
    have hD := h.Dsl none (by simp [wt])
    have hD' := hD.clm_comp ((ContinuousLinearMap.compL ℝ (EuclideanSpace ℝ (Fin n)) (P n N) W).flip
      (ContinuousLinearMap.inl ℝ (EuclideanSpace ℝ (Fin n)) (J n N)))
    exact hD'.mono (by simp [wt]) (by omega)
  have hsum := IsSymbol.sum h.isOpen (Finset.univ : Finset (Fin (N + 1)))
    (fun m _ => isSymbol_termM h hlow m)
  exact (h0.add hsum).congr fun p hp => totalD_eq_sum h hp

end Total

/-- info: 'CWFlow.Extension.IsSymbol.congr' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms IsSymbol.congr

/-- info: 'CWFlow.Extension.IsSymbol.sum' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms IsSymbol.sum

/-- info: 'CWFlow.Extension.isSymbol_comp_low' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms isSymbol_comp_low

/-- info: 'CWFlow.Extension.IsSymbol.Dsl' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms IsSymbol.Dsl

/-- info: 'CWFlow.Extension.IsSymbol.totalD' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms IsSymbol.totalD

end CWFlow.Extension
