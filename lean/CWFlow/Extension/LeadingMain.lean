import CWFlow.Extension.IBPBounds
import CWFlow.Extension.WeightTop

/-!
# The main (coercive) estimate

The top order part `2 Σ_I topPot_I · wSym_I` of the energy inequality is integrated by parts twice and
identified, up to errors controlled by `ε ∫_{B₂} |D^{k+5} φ|² + C (1 + energy)`, with
`-2 ∫ ρ² √det g · |D^{k+5} φ|²_{g⁻¹}`, which is bounded by `-2 ((1+K²)⁻¹)^{k+5} ∫ ρ² |D^{k+5} φ|²`
(`exists_const_main_le`).
-/

set_option linter.unusedSectionVars false

open CWFlow.Interp Metric MeasureTheory
open scoped ContDiff

namespace CWFlow.Extension

variable {n N : ℕ}

/-! ### Part 1: symmetry of jets via counting indices -/

/-- Two index tuples with the same number of occurrences of each index give the same jet coordinate. -/
theorem lm_jet_eq_of_count {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ) (m : ℕ)
    (x : EuclideanSpace ℝ (Fin n)) (L L' : Fin m → Fin n)
    (h : ∀ i : Fin n, (∑ j, if L j = i then 1 else 0 : ℕ) = ∑ j, if L' j = i then 1 else 0) :
    jet m φ x L = jet m φ x L' := by
  classical
  have hc : ∀ i : Fin n, Fintype.card {j // L j = i} = Fintype.card {j // L' j = i} := by
    intro i
    rw [Fintype.card_subtype, Fintype.card_subtype, Finset.card_filter, Finset.card_filter]
    exact h i
  let e : ∀ i, {j // L j = i} ≃ {j // L' j = i} := fun i => Fintype.equivOfCardEq (hc i)
  let σ : Equiv.Perm (Fin m) := (Equiv.sigmaFiberEquiv L).symm.trans
    ((Equiv.sigmaCongrRight e).trans (Equiv.sigmaFiberEquiv L'))
  have hσ : L' ∘ σ = L := by
    funext j
    have h1 : σ j = (e (L j) ⟨j, rfl⟩).1 := rfl
    simp only [Function.comp_apply, h1]
    exact (e (L j) ⟨j, rfl⟩).2
  calc jet m φ x L = jet m φ x (L' ∘ σ) := by rw [hσ]
    _ = jet m φ x L' := jet_comp_perm hφ m x σ L'

theorem lm_count_append {m r : ℕ} (I : Fin m → Fin n) (v : Fin r → Fin n) (i : Fin n) :
    (∑ j : Fin (m + r), if Fin.append I v j = i then 1 else 0 : ℕ) =
      (∑ j, if I j = i then 1 else 0) + ∑ j, if v j = i then 1 else 0 := by
  rw [Fin.sum_univ_add]
  simp [Fin.append_left, Fin.append_right]

theorem lm_count_cons {m : ℕ} (a : Fin n) (I : Fin m → Fin n) (i : Fin n) :
    (∑ j : Fin (m + 1), if (Fin.cons a I : Fin (m + 1) → Fin n) j = i then 1 else 0 : ℕ) =
      (if a = i then 1 else 0) + ∑ j, if I j = i then 1 else 0 := by
  rw [Fin.sum_univ_succ]
  simp

/-- Step 1: `φ_{I d a b c} = φ_{a c I d b}`. -/
theorem lm_jet_perm1 {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ) (k : ℕ)
    (x : EuclideanSpace ℝ (Fin n)) (I : Fin (k + 3) → Fin n) (a b c d : Fin n) :
    jet (k + 5 + 2) φ x (Fin.cons a (Fin.cons c (Fin.append I ![d, b] : Fin (k + 3 + 2) → Fin n))) =
      jet (k + 3 + 4) φ x (Fin.append I ![d, a, b, c] : Fin (k + 3 + 4) → Fin n) := by
  refine lm_jet_eq_of_count hφ (k + 7) x _ _ fun i => ?_
  have e1 : (∑ j : Fin (k + 3 + 2 + 1 + 1), if (Fin.cons a (Fin.cons c
      (Fin.append I ![d, b] : Fin (k + 3 + 2) → Fin n) : Fin (k + 3 + 2 + 1) → Fin n) :
        Fin (k + 3 + 2 + 1 + 1) → Fin n) j = i then 1 else 0 : ℕ) =
      (if a = i then 1 else 0) + ((if c = i then 1 else 0) +
        ((∑ j, if I j = i then 1 else 0) + ∑ j : Fin 2, if ![d, b] j = i then 1 else 0)) := by
    rw [lm_count_cons, lm_count_cons, add_right_inj, add_right_inj]
    exact lm_count_append I ![d, b] i
  have e2 : (∑ j : Fin (k + 3 + 4), if (Fin.append I ![d, a, b, c] : Fin (k + 3 + 4) → Fin n) j = i
      then 1 else 0 : ℕ) =
      (∑ j, if I j = i then 1 else 0) + ∑ j : Fin 4, if ![d, a, b, c] j = i then 1 else 0 :=
    lm_count_append _ _ _
  refine e1.trans (Eq.trans ?_ e2.symm)
  have h0 : (![d, a, b, c] : Fin 4 → Fin n) 0 = d := rfl
  have h1 : (![d, a, b, c] : Fin 4 → Fin n) 1 = a := rfl
  have h2 : (![d, a, b, c] : Fin 4 → Fin n) 2 = b := rfl
  have h3 : (![d, a, b, c] : Fin 4 → Fin n) 3 = c := rfl
  have g0 : (![d, b] : Fin 2 → Fin n) 0 = d := rfl
  have g1 : (![d, b] : Fin 2 → Fin n) 1 = b := rfl
  simp only [Fin.sum_univ_four, Fin.sum_univ_two, h0, h1, h2, h3, g0, g1]
  ring

/-- `φ_{a c J} = φ_{J c a}`. -/
theorem lm_jet_perm2 {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ) (k : ℕ)
    (x : EuclideanSpace ℝ (Fin n)) (J : Fin (k + 3) → Fin n) (a c : Fin n) :
    jet (k + 5) φ x (Fin.cons a (Fin.cons c J)) =
      jet (k + 5) φ x (Fin.append J ![c, a] : Fin (k + 3 + 2) → Fin n) := by
  refine lm_jet_eq_of_count hφ (k + 5) x _ _ fun i => ?_
  have e1 : (∑ j : Fin (k + 3 + 1 + 1), if (Fin.cons a (Fin.cons c J : Fin (k + 3 + 1) → Fin n) :
        Fin (k + 3 + 1 + 1) → Fin n) j = i then 1 else 0 : ℕ) =
      (if a = i then 1 else 0) + ((if c = i then 1 else 0) + ∑ j, if J j = i then 1 else 0) := by
    rw [lm_count_cons, lm_count_cons]
  have e2 : (∑ j : Fin (k + 3 + 2), if (Fin.append J ![c, a] : Fin (k + 3 + 2) → Fin n) j = i
      then 1 else 0 : ℕ) =
      (∑ j, if J j = i then 1 else 0) + ∑ j : Fin 2, if ![c, a] j = i then 1 else 0 :=
    lm_count_append _ _ _
  refine e1.trans (Eq.trans ?_ e2.symm)
  have g0 : (![c, a] : Fin 2 → Fin n) 0 = c := rfl
  have g1 : (![c, a] : Fin 2 → Fin n) 1 = a := rfl
  simp only [Fin.sum_univ_two, g0, g1]
  ring

/-! ### Part 2: the algebraic identity behind the main term -/

theorem lm_prod_split {k : ℕ} (g : Fin n → Fin n → ℝ) (I J : Fin (k + 3) → Fin n)
    (v w : Fin 2 → Fin n) :
    ∏ i : Fin (k + 3 + 2), g ((Fin.append I v : Fin (k + 3 + 2) → Fin n) i)
        ((Fin.append J w : Fin (k + 3 + 2) → Fin n) i) =
      (∏ i, g (I i) (J i)) * (g (v 0) (w 0) * g (v 1) (w 1)) := by
  rw [Fin.prod_univ_add]
  simp [Fin.append_left, Fin.append_right, Fin.prod_univ_two]

theorem lm_algebra {k : ℕ} (g : Fin n → Fin n → ℝ) (hg : ∀ a b, g a b = g b a)
    (T T' : (Fin (k + 3 + 2) → Fin n) → ℝ) :
    ∑ I : Fin (k + 3) → Fin n, ∑ a, ∑ b, ∑ c, ∑ d, ∑ J : Fin (k + 3) → Fin n,
      (g c d * g a b * ∏ i, g (I i) (J i)) * T (Fin.append I ![d, b]) * T' (Fin.append J ![c, a]) =
    ∑ L, ∑ L', (∏ i, g (L i) (L' i)) * T L * T' L' := by
  classical
  let e1 : ((Fin (k + 3) → Fin n) × Fin n × Fin n × Fin n × Fin n × (Fin (k + 3) → Fin n)) ≃
      ((Fin (k + 3) → Fin n) × (Fin n × Fin n)) × ((Fin (k + 3) → Fin n) × (Fin n × Fin n)) :=
    { toFun := fun x => ((x.1, (x.2.2.2.2.1, x.2.2.1)), (x.2.2.2.2.2, (x.2.2.2.1, x.2.1)))
      invFun := fun y => (y.1.1, y.2.2.2, y.1.2.2, y.2.2.1, y.1.2.1, y.2.1)
      left_inv := fun x => rfl
      right_inv := fun y => rfl }
  let e3 : ((Fin (k + 3) → Fin n) × (Fin n × Fin n)) ≃ (Fin (k + 3 + 2) → Fin n) :=
    (Equiv.prodCongr (Equiv.refl _) (finTwoArrowEquiv (Fin n)).symm).trans
      (Fin.appendEquiv (k + 3) 2)
  let e : ((Fin (k + 3) → Fin n) × Fin n × Fin n × Fin n × Fin n × (Fin (k + 3) → Fin n)) ≃
      (Fin (k + 3 + 2) → Fin n) × (Fin (k + 3 + 2) → Fin n) :=
    e1.trans (Equiv.prodCongr e3 e3)
  have hL : (∑ I : Fin (k + 3) → Fin n, ∑ a, ∑ b, ∑ c, ∑ d, ∑ J : Fin (k + 3) → Fin n,
      (g c d * g a b * ∏ i, g (I i) (J i)) * T (Fin.append I ![d, b]) * T' (Fin.append J ![c, a])) =
      ∑ x : (Fin (k + 3) → Fin n) × Fin n × Fin n × Fin n × Fin n × (Fin (k + 3) → Fin n),
        (g x.2.2.2.1 x.2.2.2.2.1 * g x.2.1 x.2.2.1 * ∏ i, g (x.1 i) (x.2.2.2.2.2 i)) *
          T (Fin.append x.1 ![x.2.2.2.2.1, x.2.2.1]) *
          T' (Fin.append x.2.2.2.2.2 ![x.2.2.2.1, x.2.1]) := by
    simp only [Fintype.sum_prod_type]
  have hR : (∑ L, ∑ L', (∏ i, g (L i) (L' i)) * T L * T' L') =
      ∑ y : (Fin (k + 3 + 2) → Fin n) × (Fin (k + 3 + 2) → Fin n),
        (∏ i, g (y.1 i) (y.2 i)) * T y.1 * T' y.2 := by
    simp only [Fintype.sum_prod_type]
  rw [hL, hR]
  refine Fintype.sum_equiv e _ _ fun x => ?_
  obtain ⟨I, a, b, c, d, J⟩ := x
  have h1 : e (I, a, b, c, d, J) = (Fin.append I ![d, b], Fin.append J ![c, a]) := rfl
  rw [h1]
  dsimp only
  rw [lm_prod_split g I J ![d, b] ![c, a]]
  have h2 : (![d, b] : Fin 2 → Fin n) 0 = d := rfl
  have h3 : (![d, b] : Fin 2 → Fin n) 1 = b := rfl
  have h4 : (![c, a] : Fin 2 → Fin n) 0 = c := rfl
  have h5 : (![c, a] : Fin 2 → Fin n) 1 = a := rfl
  rw [h2, h3, h4, h5, hg d c, hg b a]
  ring

/-! ### Part 3: derivatives of scalar symbols -/

/-- The partial derivative `D_c 𝔣` of a scalar function on jet space (total derivative in the direction
`e_c`). -/
noncomputable def lm_Dc (𝔣 : P n N → ℝ) (c : Fin n) : P n N → ℝ :=
  fun q => totalD 𝔣 q (EuclideanSpace.single c 1)

theorem lm_totalD_congr {f g : P n N → ℝ} {p : P n N} (h : f =ᶠ[nhds p] g)
    (e : EuclideanSpace ℝ (Fin n)) : totalD f p e = totalD g p e := by
  rw [totalD_apply, totalD_apply, h.fderiv_eq]

theorem lm_isSymbol_Dc {A : Set (P n N)} {w : ℕ} {f : P n N → ℝ} (h : IsSymbol A w w f)
    (hN : w + 4 ≤ N)
    (hlow : ∀ m : Fin (N + 1), (m : ℕ) ≤ 3 → ∃ K : ℝ, ∀ p ∈ A, ‖p.2 m‖ ≤ K) (c : Fin n) :
    IsSymbol A (w + 1) (w + 1) (lm_Dc f c) :=
  ((h.totalD hN hlow).clm_comp
    (ContinuousLinearMap.apply ℝ ℝ (EuclideanSpace.single c (1 : ℝ) :
      EuclideanSpace ℝ (Fin n)))).mono (le_refl _) (by omega)

theorem lm_DD_mul {A : Set (P n N)} (hA : IsOpen A) {f g : P n N → ℝ} (hf : ContDiffOn ℝ ∞ f A)
    (hg : ContDiffOn ℝ ∞ g A) (a c : Fin n) (hfc : ContDiffOn ℝ ∞ (lm_Dc f c) A)
    (hgc : ContDiffOn ℝ ∞ (lm_Dc g c) A) {p : P n N} (hp : p ∈ A) :
    lm_Dc (lm_Dc (fun q => f q * g q) c) a p =
      f p * lm_Dc (lm_Dc g c) a p + lm_Dc f c p * lm_Dc g a p + lm_Dc f a p * lm_Dc g c p +
        lm_Dc (lm_Dc f c) a p * g p := by
  have hdiff : ∀ {h : P n N → ℝ}, ContDiffOn ℝ ∞ h A → ∀ q ∈ A, DifferentiableAt ℝ h q :=
    fun hh q hq => (hh.contDiffAt (hA.mem_nhds hq)).differentiableAt (by simp)
  have hev : lm_Dc (fun q => f q * g q) c =ᶠ[nhds p]
      fun q => lm_Dc f c q * g q + f q * lm_Dc g c q := by
    filter_upwards [hA.mem_nhds hp] with q hq
    exact totalD_mul (hdiff hf q hq) (hdiff hg q hq) _
  have d1 := hdiff hf p hp
  have d2 := hdiff hg p hp
  have d3 := hdiff hfc p hp
  have d4 := hdiff hgc p hp
  show totalD (lm_Dc (fun q => f q * g q) c) p (EuclideanSpace.single a 1) = _
  have hadd : totalD (fun q => lm_Dc f c q * g q + f q * lm_Dc g c q) p =
      totalD (fun q => lm_Dc f c q * g q) p + totalD (fun q => f q * lm_Dc g c q) p :=
    totalD_add (d3.mul d2) (d1.mul d4)
  rw [lm_totalD_congr hev, hadd, add_apply, totalD_mul d3 d2, totalD_mul d1 d4]
  simp only [lm_Dc]
  ring

/-! ### Part 4: the coefficient of the top term -/

/-- The top order coefficient `√det g · Σ_J (Π_i g^{I_i J_i}) ζ_{k+5}(a, c, J)`. -/
noncomputable def lm_M (k : ℕ) (a c : Fin n) (I : Fin (k + 3) → Fin n) (p : P n N) : ℝ :=
  volFun (low23 p).1 * ∑ J : Fin (k + 3) → Fin n,
    (∏ i, ginvFun (low23 p).1 ![I i, J i]) * slotN (k + 3 + 1 + 1) p (Fin.cons a (Fin.cons c J))

theorem lm_wTop_eq (k : ℕ) (a c : Fin n) (I : Fin (k + 3) → Fin n) (p : P n N) :
    wTop k p (Fin.cons a (Fin.cons c I)) = lm_M k a c I p := by
  unfold wTop lm_M
  rw [liftTail_apply]
  simp only [Fin.cons_zero, Fin.tail_cons]
  rw [liftTail_apply]
  simp only [Fin.cons_zero, Fin.tail_cons]
  simp only [wCoef, smul_apply, PiLp.smul_apply, smul_eq_mul, sharp_apply,
    sliceZ_apply]

theorem lm_DD_wSym {K : ℝ} (k : ℕ) (hN : k + 5 ≤ N) (a c : Fin n) (I : Fin (k + 3) → Fin n)
    {p : P n N} (hp : p ∈ admSet n N K) :
    lm_Dc (lm_Dc (fun q => wSym k q I) c) a p =
      Dz (Dz (wSym (N := N) k)) p (Fin.cons a (Fin.cons c I)) := by
  have hA := isOpen_admSet (n := n) (N := N) K
  have hlow := admSet_low (n := n) (N := N) K
  have hdz : DifferentiableAt ℝ (Dz (wSym (N := N) k)) p :=
    wt_diff (isSymbol_Dz_wSym hA hlow k (by omega)) hp
  rw [Dz_apply_eq]
  simp only [Fin.cons_zero, Fin.tail_cons]
  rw [totalD_apply_coord hdz]
  show totalD (fun q => totalD (fun r => wSym k r I) q (EuclideanSpace.single c 1)) p
    (EuclideanSpace.single a 1) = totalD (fun q => Dz (wSym (N := N) k) q (Fin.cons c I)) p
    (EuclideanSpace.single a 1)
  refine lm_totalD_congr ?_ _
  filter_upwards [hA.mem_nhds hp] with q hq
  have hw : DifferentiableAt ℝ (wSym (N := N) k) q :=
    wt_diff (isSymbol_wSym hA hlow k (by omega)) hq
  rw [Dz_apply_eq]
  simp only [Fin.cons_zero, Fin.tail_cons]
  rw [totalD_apply_coord hw]

/-! ### Part 5: the error symbol -/

/-- The scalar symbol `γ ω` with `γ = g^{cd} g^{ab}` and `ω = w_I`. -/
noncomputable def lm_F (k : ℕ) (I : Fin (k + 3) → Fin n) (a b c d : Fin n) : P n N → ℝ :=
  fun q => gg a b c d q * wSym k q I

/-- `D_a D_c (γ ω)`. -/
noncomputable def lm_G (k : ℕ) (I : Fin (k + 3) → Fin n) (a b c d : Fin n) : P n N → ℝ :=
  lm_Dc (lm_Dc (lm_F k I a b c d) c) a

/-- `D_a D_c (γ ω) - γ · (top coefficient)`: the lower order error of the second derivative. -/
noncomputable def lm_err (k : ℕ) (I : Fin (k + 3) → Fin n) (a b c d : Fin n) : P n N → ℝ :=
  fun p => lm_G k I a b c d p - gg a b c d p * lm_M k a c I p

theorem lm_isSymbol_F {K : ℝ} (k : ℕ) (hN : k + 3 ≤ N) (I : Fin (k + 3) → Fin n) (a b c d : Fin n) :
    IsSymbol (admSet n N K) k k (lm_F k I a b c d) :=
  ((isSymbol_gg (isOpen_admSet K) (admSet_low K) a b c d).mul
    (isSymbol_wSym_apply (isOpen_admSet K) (admSet_low K) k hN I)).mono (by omega) (by omega)

theorem lm_isSymbol_err {K : ℝ} (k : ℕ) (hN : k + 5 ≤ N) (I : Fin (k + 3) → Fin n)
    (a b c d : Fin n) : IsSymbol (admSet n N K) (k + 1) (k + 1) (lm_err k I a b c d) := by
  have hA := isOpen_admSet (n := n) (N := N) K
  have hlow := admSet_low (n := n) (N := N) K
  have hγ : IsSymbol (admSet n N K) 0 0 (gg a b c d : P n N → ℝ) := isSymbol_gg hA hlow a b c d
  have hω : IsSymbol (admSet n N K) k k (fun q : P n N => wSym k q I) :=
    isSymbol_wSym_apply hA hlow k (by omega) I
  have hDγ : ∀ e : Fin n, IsSymbol (admSet n N K) 0 0 (lm_Dc (gg a b c d : P n N → ℝ) e) :=
    fun e => isSymbol_totalD_gg hA hlow a b c d e
  have hDDγ : IsSymbol (admSet n N K) 1 1 (lm_Dc (lm_Dc (gg a b c d : P n N → ℝ) c) a) :=
    lm_isSymbol_Dc (w := 0) (hDγ c) (by omega) hlow a
  have hDω : ∀ e : Fin n, IsSymbol (admSet n N K) (k + 1) (k + 1)
      (lm_Dc (fun q : P n N => wSym k q I) e) := fun e => lm_isSymbol_Dc hω (by omega) hlow e
  have hR : IsSymbol (admSet n N K) (k + 1) (k + 1) (fun p : P n N =>
      (Dz (Dz (wSym (N := N) k)) p - wTop k p) (Fin.cons a (Fin.cons c I))) :=
    (isSymbol_Dz_Dz_wSym_sub hA hlow k hN).clm_comp
      (PiLp.proj 2 (fun _ : Fin (k + 3 + 1 + 1) → Fin n => ℝ) (Fin.cons a (Fin.cons c I)))
  have S1 := (hγ.mul hR).mono (w' := k + 1) (ℓ' := k + 1) (by omega) (by omega)
  have S2 := ((hDγ c).mul (hDω a)).mono (w' := k + 1) (ℓ' := k + 1) (by omega) (by omega)
  have S3 := ((hDγ a).mul (hDω c)).mono (w' := k + 1) (ℓ' := k + 1) (by omega) (by omega)
  have S4 := (hDDγ.mul hω).mono (w' := k + 1) (ℓ' := k + 1) (by omega) (by omega)
  refine (((S1.add S2).add S3).add S4).congr fun p hp => ?_
  have key := lm_DD_mul hA hγ.smooth hω.smooth a c
    (hDγ c).smooth (hDω c).smooth hp
  have hG : lm_G k I a b c d p = _ := key
  unfold lm_err
  rw [hG, lm_DD_wSym k hN a c I hp, ← lm_wTop_eq k a c I p, PiLp.sub_apply]
  ring

/-! ### Part 6: error terms -/

/-- Generic bound for an error term: a symbol of weight `(k+1, k+1)` against a jet coordinate of order
`k + 5`. -/
theorem lm_error_bound {k : ℕ} (hN : k + 7 ≤ N) {K : ℝ} {𝔢 : P n N → ℝ}
    (h𝔢 : IsSymbol (admSet n N K) (k + 1) (k + 1) 𝔢) {ε : ℝ} (hε : 0 < ε) (Λ : ℝ) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ ρ : EuclideanSpace ℝ (Fin n) → ℝ, IsCutoff ρ Λ →
      ∀ φ : EuclideanSpace ℝ (Fin n) → ℝ, ContDiff ℝ ∞ φ →
      (∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 3, jp N φ x ∈ admSet n N K) →
      ∀ L : Fin (k + 5) → Fin n,
        |∫ x, jet (k + 5) φ x L * locWeight ρ 𝔢 φ x| ≤
          ε * (∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 2, ‖jet (k + 5) φ x‖ ^ 2) +
          C * (1 + energy3 (k + 1) 3 φ) := by
  obtain ⟨C₀, hC₀, hpt⟩ := ibb_pt h𝔢 (by omega : 3 + min (k + 1) (k + 1) ≤ N)
  obtain ⟨C₂, hC₂, hom⟩ := exists_const_integral_Om_sq (n := n) (k + 1) (k + 1) (ρ := 2) (R := 3)
    (by norm_num) (by norm_num) K
  refine ⟨C₀ ^ 2 / (4 * ε) * C₂, by positivity, fun ρ hρ φ hφ hadm L => ?_⟩
  have hK : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 3, ‖jet 3 φ x‖ ≤ K :=
    fun x hx => ibb_norm3 (by omega) (hadm x hx)
  have hsupp : ∀ x ∈ tsupport ρ, jp N φ x ∈ admSet n N K :=
    fun x hx => hadm x (ibb_ball_sub x (hρ.supp hx))
  have hw := contDiff_locWeight h𝔢 hφ hρ.smooth hsupp
  have := ibb_young (k := k) (ε := ε) (C₁ := C₀) (C₂ := C₂) hε hφ (hom φ hφ hK)
    (h := fun x => jet (k + 5) φ x L * locWeight ρ 𝔢 φ x)
    ((ibb_coord_contDiff hφ _ L).continuous.mul hw.continuous)
    (fun x hx => by
      have hxs : x ∉ tsupport ρ := fun h => hx (hρ.supp h)
      simp [locWeight_eq_zero 𝔢 φ hxs])
    (fun x hx => by
      by_cases hxs : x ∈ tsupport ρ
      · rw [abs_mul]
        have h1 := hpt φ x (hsupp x hxs)
        rw [Real.norm_eq_abs] at h1
        have h2 : |ρ x ^ 2| ≤ 1 := by
          rw [abs_pow]; exact pow_le_one₀ (abs_nonneg _) (hρ.le_one x)
        refine mul_le_mul (ibb_coord_le _ x L) ?_ (abs_nonneg _) (norm_nonneg _)
        unfold locWeight
        rw [abs_mul]
        calc |ρ x ^ 2| * |𝔢 (jp N φ x)| ≤ 1 * (C₀ * Om (k + 1) (k + 1) (aJet φ x)) :=
              mul_le_mul h2 h1 (abs_nonneg _) zero_le_one
          _ = _ := one_mul _
      · simp only [locWeight_eq_zero 𝔢 φ hxs, mul_zero, abs_zero]
        exact mul_nonneg (norm_nonneg _) (mul_nonneg hC₀
          (zero_le_one.trans (one_le_Om _ _ (aJet_nonneg φ x)))))
  refine this.trans (le_of_eq ?_)
  ring

theorem lm_integrable_smooth {F : P n N → ℝ} (hF : ContDiff ℝ ∞ F) {Λ : ℝ}
    {ρ : EuclideanSpace ℝ (Fin n) → ℝ} (hρ : IsCutoff ρ Λ) {φ : EuclideanSpace ℝ (Fin n) → ℝ}
    (hφ : ContDiff ℝ ∞ φ) (m : ℕ) (L : Fin m → Fin n) :
    Integrable (fun x => jet m φ x L * (ρ x ^ 2 * F (jp N φ x))) :=
  ((ibb_coord_contDiff hφ m L).continuous.mul
    ((hρ.smooth.continuous.pow 2).mul (hF.continuous.comp (contDiff_jp hφ).continuous))
    ).integrable_of_hasCompactSupport
    (hasCompactSupport_locWeight (ibb_hasCompactSupport hρ) F φ).mul_left

theorem lm_integrable_sym {A : Set (P n N)} {w ℓ : ℕ} {𝔣 : P n N → ℝ} (h : IsSymbol A w ℓ 𝔣)
    {Λ : ℝ} {ρ : EuclideanSpace ℝ (Fin n) → ℝ} (hρ : IsCutoff ρ Λ)
    {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ)
    (hsupp : ∀ x ∈ tsupport ρ, jp N φ x ∈ A) (m : ℕ) (L : Fin m → Fin n) :
    Integrable (fun x => jet m φ x L * locWeight ρ 𝔣 φ x) :=
  ((ibb_coord_contDiff hφ m L).continuous.mul
    (contDiff_locWeight h hφ hρ.smooth hsupp).continuous).integrable_of_hasCompactSupport
    (hasCompactSupport_locWeight (ibb_hasCompactSupport hρ) 𝔣 φ).mul_left

theorem lm_contDiff_M (k : ℕ) (a c : Fin n) (I : Fin (k + 3) → Fin n) :
    ContDiff ℝ ∞ (lm_M k a c I : P n N → ℝ) := by
  have hu : ContDiff ℝ ∞ (fun q : P n N => (low23 q).1) := contDiff_fst.comp low23.contDiff
  unfold lm_M
  refine (contDiff_volFun.comp hu).mul (ContDiff.sum fun J _ => ?_)
  refine (contDiff_prod fun i _ => ?_).mul (contDiff_slotN_apply _ _)
  exact (PiLp.proj (𝕜 := ℝ) 2 (fun _ : Fin 2 → Fin n => ℝ) ![I i, J i]).contDiff.comp
    (contDiff_ginvFun.comp hu)

theorem lm_contDiff_gM (k : ℕ) (a b c d : Fin n) (I : Fin (k + 3) → Fin n) :
    ContDiff ℝ ∞ (fun q : P n N => gg a b c d q * lm_M k a c I q) :=
  (contDiff_gg a b c d).mul (lm_contDiff_M k a c I)

/-- The bound for one term of the sum: two integrations by parts, up to controlled errors, give the main
term. -/
theorem lm_term_bound {k : ℕ} (hN : k + 7 ≤ N) {K : ℝ} (Λ : ℝ) {ε : ℝ} (hε : 0 < ε)
    (I : Fin (k + 3) → Fin n) (a b c d : Fin n) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ ρ : EuclideanSpace ℝ (Fin n) → ℝ, IsCutoff ρ Λ →
      ∀ φ : EuclideanSpace ℝ (Fin n) → ℝ, ContDiff ℝ ∞ φ →
      (∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 3, jp N φ x ∈ admSet n N K) →
      |(∫ x, jet (k + 5 + 2) φ x (Fin.cons a (Fin.cons c
            (Fin.append I ![d, b] : Fin (k + 3 + 2) → Fin n))) *
              locWeight ρ (lm_F (N := N) k I a b c d) φ x) -
        ∫ x, jet (k + 5) φ x (Fin.append I ![d, b] : Fin (k + 3 + 2) → Fin n) *
          (ρ x ^ 2 * (gg a b c d (jp N φ x) * lm_M k a c I (jp N φ x)))| ≤
        ε * (∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 2, ‖jet (k + 5) φ x‖ ^ 2) +
        C * (1 + energy3 (k + 1) 3 φ) := by
  obtain ⟨C₁, hC₁, h1⟩ := exists_const_integral_jet2_locWeight (by omega : k + 5 ≤ N) Λ
    (lm_isSymbol_F (K := K) k (by omega) I a b c d) (ε := ε / 2) (by positivity)
  obtain ⟨C₂, hC₂, h2⟩ := lm_error_bound hN (lm_isSymbol_err (K := K) k (by omega) I a b c d)
    (ε := ε / 2) (by positivity) Λ
  refine ⟨C₁ + C₂, by positivity, fun ρ hρ φ hφ hadm => ?_⟩
  set L : Fin (k + 3 + 2) → Fin n := Fin.append I ![d, b] with hL
  have t1 := h1 ρ hρ φ hφ hadm a c L
  have t2 := h2 ρ hρ φ hφ hadm L
  have hsupp : ∀ x ∈ tsupport ρ, jp N φ x ∈ admSet n N K :=
    fun x hx => hadm x (ibb_ball_sub x (hρ.supp hx))
  have i1 := lm_integrable_sym (lm_isSymbol_err (K := K) k (by omega) I a b c d) hρ hφ hsupp
    (k + 5) L
  have i2 := lm_integrable_smooth (lm_contDiff_gM (N := N) k a b c d I) hρ hφ (k + 5) L
  have hsplit : (∫ x, jet (k + 5) φ x L *
      (ρ x ^ 2 * totalD (fun q => totalD (lm_F (N := N) k I a b c d) q (EuclideanSpace.single c 1)) (jp N φ x)
        (EuclideanSpace.single a 1))) =
      (∫ x, jet (k + 5) φ x L * (ρ x ^ 2 * (gg a b c d (jp N φ x) * lm_M k a c I (jp N φ x)))) +
        ∫ x, jet (k + 5) φ x L * locWeight ρ (lm_err (N := N) k I a b c d) φ x := by
    rw [← integral_add i2 i1]
    refine integral_congr_ae (Filter.Eventually.of_forall fun x => ?_)
    show jet (k + 5) φ x L * (ρ x ^ 2 * lm_G k I a b c d (jp N φ x)) = _
    unfold locWeight lm_err
    ring
  have hdiff : (∫ x, jet (k + 5) φ x L * (ρ x ^ 2 * totalD (fun q => totalD (lm_F (N := N) k I a b c d) q
      (EuclideanSpace.single c 1)) (jp N φ x) (EuclideanSpace.single a 1))) -
      (∫ x, jet (k + 5) φ x L * (ρ x ^ 2 * (gg a b c d (jp N φ x) * lm_M k a c I (jp N φ x)))) =
      ∫ x, jet (k + 5) φ x L * locWeight ρ (lm_err (N := N) k I a b c d) φ x := by
    rw [hsplit]
    ring
  have key := abs_sub_le
    (∫ x, jet (k + 5 + 2) φ x (Fin.cons a (Fin.cons c L)) * locWeight ρ (lm_F (N := N) k I a b c d) φ x)
    (∫ x, jet (k + 5) φ x L * (ρ x ^ 2 * totalD (fun q => totalD (lm_F (N := N) k I a b c d) q
      (EuclideanSpace.single c 1)) (jp N φ x) (EuclideanSpace.single a 1)))
    (∫ x, jet (k + 5) φ x L * (ρ x ^ 2 * (gg a b c d (jp N φ x) * lm_M k a c I (jp N φ x))))
  rw [hdiff] at key
  have := key.trans (add_le_add t1 t2)
  refine this.trans (le_of_eq ?_)
  ring

/-! ### Part 7: pointwise identities for the sums -/

/-- The integrand of the left hand side as a sum of localized weights against `(k+7)`-jets. -/
theorem lm_lhs_pointwise {k : ℕ} (hN : k + 7 ≤ N) {φ : EuclideanSpace ℝ (Fin n) → ℝ}
    (hφ : ContDiff ℝ ∞ φ) (ρ : EuclideanSpace ℝ (Fin n) → ℝ) (x : EuclideanSpace ℝ (Fin n)) :
    ρ x ^ 2 * (2 * ∑ I : Fin (k + 3) → Fin n, topPot (k + 3) (jp N φ x) I * wSym k (jp N φ x) I) =
      ∑ I : Fin (k + 3) → Fin n, ∑ a, ∑ b, ∑ c, ∑ d, (-2) *
        (jet (k + 5 + 2) φ x (Fin.cons a (Fin.cons c
          (Fin.append I ![d, b] : Fin (k + 3 + 2) → Fin n))) *
          locWeight ρ (lm_F (N := N) k I a b c d) φ x) := by
  rw [Finset.mul_sum, Finset.mul_sum]
  refine Finset.sum_congr rfl fun I _ => ?_
  have hR : (∑ a, ∑ b, ∑ c, ∑ d, (-2 : ℝ) *
        (jet (k + 5 + 2) φ x (Fin.cons a (Fin.cons c
          (Fin.append I ![d, b] : Fin (k + 3 + 2) → Fin n))) *
          locWeight ρ (lm_F (N := N) k I a b c d) φ x)) =
      (-2 * ρ x ^ 2 * wSym k (jp N φ x) I) * ∑ a, ∑ b, ∑ c, ∑ d,
        gg a b c d (jp N φ x) * slotN (k + 3 + 4) (jp N φ x) (Fin.append I ![d, a, b, c]) := by
    simp only [Finset.mul_sum]
    refine Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun b _ =>
      Finset.sum_congr rfl fun c _ => Finset.sum_congr rfl fun d _ => ?_
    rw [lm_jet_perm1 hφ k x I a b c d, slotN_jp (k + 3 + 4) (by omega)]
    unfold locWeight lm_F
    ring
  rw [hR, topPot_apply]
  ring

/-- The main term: `Σ_{I,a,b,c,d} φ_L ρ² γ M = ρ² √det g |D^{k+5} φ|²_{g⁻¹}`. -/
theorem lm_main_pointwise {k : ℕ} (hN : k + 5 ≤ N) {φ : EuclideanSpace ℝ (Fin n) → ℝ}
    (hφ : ContDiff ℝ ∞ φ) (r : ℝ) (x : EuclideanSpace ℝ (Fin n)) :
    ∑ I : Fin (k + 3) → Fin n, ∑ a, ∑ b, ∑ c, ∑ d,
      jet (k + 5) φ x (Fin.append I ![d, b] : Fin (k + 3 + 2) → Fin n) *
        (r * (gg a b c d (jp N φ x) * lm_M k a c I (jp N φ x))) =
    r * (volFun (low23 (jp N φ x)).1 *
      gNormSq (ginvFun (low23 (jp N φ x)).1) (jet (k + 5) φ x)) := by
  set u := (low23 (jp N φ x)).1 with hu
  have h1 : ∀ (I : Fin (k + 3) → Fin n) (a b c d : Fin n),
      jet (k + 5) φ x (Fin.append I ![d, b] : Fin (k + 3 + 2) → Fin n) *
        (r * (gg a b c d (jp N φ x) * lm_M k a c I (jp N φ x))) =
      ∑ J : Fin (k + 3) → Fin n, (r * volFun u) *
        (((ginvFun u ![c, d] * ginvFun u ![a, b] * ∏ i, ginvFun u ![I i, J i]) *
          jet (k + 5) φ x (Fin.append I ![d, b] : Fin (k + 3 + 2) → Fin n)) *
          jet (k + 5) φ x (Fin.append J ![c, a] : Fin (k + 3 + 2) → Fin n)) := by
    intro I a b c d
    simp only [lm_M, Finset.mul_sum]
    refine Finset.sum_congr rfl fun J _ => ?_
    rw [← lm_jet_perm2 hφ k x J a c, slotN_jp (k + 3 + 1 + 1) (by omega)]
    simp only [gg, hu]
    ring
  simp only [h1, ← Finset.mul_sum]
  rw [lm_algebra (fun a b => ginvFun u ![a, b]) (fun a b => ginvFun_symm u a b)
    (fun L => jet (k + 5) φ x L) (fun L => jet (k + 5) φ x L)]
  rw [mul_assoc]
  rfl

/-! ### Part 8: coercivity and finite sums -/

theorem lm_coercive {K : ℝ} (hN : 2 ≤ N) (m : ℕ) {p : P n N} (hp : p ∈ admSet n N K) (T : Z n m) :
    ((1 + K ^ 2)⁻¹) ^ m * ‖T‖ ^ 2 ≤ volFun (low23 p).1 * gNormSq (ginvFun (low23 p).1) T := by
  have hu : ‖(low23 p).1‖ ≤ K := by
    rw [low23_fst p (by omega)]
    exact (hp ⟨2, by omega⟩ (by simp)).le
  have hlam : (1 + K ^ 2)⁻¹ ≤ (1 + ‖(low23 p).1‖ ^ 2)⁻¹ :=
    inv_anti₀ (by positivity) (by nlinarith [norm_nonneg (low23 p).1])
  have h1 : ((1 + K ^ 2)⁻¹) ^ m * ‖T‖ ^ 2 ≤ gNormSq (ginvFun (low23 p).1) T := by
    refine le_gNormSq (by positivity) (ginv_hsymm _) (ginv_hnn _) ?_ T
    intro v
    simp only [ginvFun_apply]
    exact le_trans (mul_le_mul_of_nonneg_right hlam (by positivity))
      (ginvFun_quadratic_ge _ v)
  exact h1.trans (le_mul_of_one_le_left (gNormSq_nonneg (ginv_hsymm _) (ginv_hnn _) _)
    (one_le_volFun _))

theorem lm_integrable_rho {Λ : ℝ} {ρ : EuclideanSpace ℝ (Fin n) → ℝ} (hρ : IsCutoff ρ Λ)
    {G : EuclideanSpace ℝ (Fin n) → ℝ} (hG : Continuous G) :
    Integrable (fun x => ρ x ^ 2 * G x) := by
  have hc : Continuous (fun x => ρ x ^ 2 * G x) := (hρ.smooth.continuous.pow 2).mul hG
  refine hc.integrable_of_hasCompactSupport ?_
  refine (ibb_hasCompactSupport hρ).mono' ?_
  intro x hx
  refine subset_tsupport ρ ?_
  rw [Function.mem_support] at hx ⊢
  intro h0
  apply hx
  simp [h0]

/-- The five-fold sum over `(I, a, b, c, d)`. -/
noncomputable def lm_S5 (k : ℕ)
    (f : (Fin (k + 3) → Fin n) → Fin n → Fin n → Fin n → Fin n → ℝ) : ℝ :=
  ∑ I, ∑ a, ∑ b, ∑ c, ∑ d, f I a b c d

variable {k : ℕ} {f g : (Fin (k + 3) → Fin n) → Fin n → Fin n → Fin n → Fin n → ℝ}

theorem lm_S5_mono (h : ∀ I a b c d, f I a b c d ≤ g I a b c d) : lm_S5 k f ≤ lm_S5 k g :=
  Finset.sum_le_sum fun I _ => Finset.sum_le_sum fun a _ => Finset.sum_le_sum fun b _ =>
    Finset.sum_le_sum fun c _ => Finset.sum_le_sum fun d _ => h I a b c d

theorem lm_S5_nonneg (h : ∀ I a b c d, 0 ≤ f I a b c d) : 0 ≤ lm_S5 k f :=
  Finset.sum_nonneg fun I _ => Finset.sum_nonneg fun a _ => Finset.sum_nonneg fun b _ =>
    Finset.sum_nonneg fun c _ => Finset.sum_nonneg fun d _ => h I a b c d

theorem lm_S5_add :
    lm_S5 k (fun I a b c d => f I a b c d + g I a b c d) = lm_S5 k f + lm_S5 k g := by
  simp only [lm_S5, Finset.sum_add_distrib]

theorem lm_S5_sub :
    lm_S5 k (fun I a b c d => f I a b c d - g I a b c d) = lm_S5 k f - lm_S5 k g := by
  simp only [lm_S5, Finset.sum_sub_distrib]

theorem lm_S5_const_mul (r : ℝ) :
    lm_S5 k (fun I a b c d => r * f I a b c d) = r * lm_S5 k f := by
  simp only [lm_S5, ← Finset.mul_sum]

theorem lm_S5_mul_const (r : ℝ) :
    lm_S5 k (fun I a b c d => f I a b c d * r) = lm_S5 k f * r := by
  simp only [lm_S5, ← Finset.sum_mul]

theorem lm_S5_const (r : ℝ) :
    lm_S5 k (fun (_ : Fin (k + 3) → Fin n) _ _ _ _ => r) =
      lm_S5 k (fun (_ : Fin (k + 3) → Fin n) _ _ _ _ => (1 : ℝ)) * r := by
  rw [← lm_S5_mul_const]
  simp only [one_mul]

theorem lm_integral_S5 {F : (Fin (k + 3) → Fin n) → Fin n → Fin n → Fin n → Fin n →
      EuclideanSpace ℝ (Fin n) → ℝ}
    (hF : ∀ I a b c d, Integrable (F I a b c d)) :
    ∫ x, lm_S5 k (fun I a b c d => F I a b c d x) = lm_S5 k (fun I a b c d => ∫ x, F I a b c d x) := by
  unfold lm_S5
  rw [integral_finsetSum _ fun I _ => integrable_finsetSum _ fun a _ =>
    integrable_finsetSum _ fun b _ => integrable_finsetSum _ fun c _ =>
      integrable_finsetSum _ fun d _ => hF I a b c d]
  refine Finset.sum_congr rfl fun I _ => ?_
  rw [integral_finsetSum _ fun a _ => integrable_finsetSum _ fun b _ =>
    integrable_finsetSum _ fun c _ => integrable_finsetSum _ fun d _ => hF I a b c d]
  refine Finset.sum_congr rfl fun a _ => ?_
  rw [integral_finsetSum _ fun b _ => integrable_finsetSum _ fun c _ =>
    integrable_finsetSum _ fun d _ => hF I a b c d]
  refine Finset.sum_congr rfl fun b _ => ?_
  rw [integral_finsetSum _ fun c _ => integrable_finsetSum _ fun d _ => hF I a b c d]
  refine Finset.sum_congr rfl fun c _ => ?_
  rw [integral_finsetSum _ fun d _ => hF I a b c d]

/-! ### Part 9: the main estimate -/

/-- The term `∫ φ_{a c L} · ρ² γ w_I`. -/
noncomputable def lm_Tf (k : ℕ) (ρ φ : EuclideanSpace ℝ (Fin n) → ℝ) :
    (Fin (k + 3) → Fin n) → Fin n → Fin n → Fin n → Fin n → ℝ := fun I a b c d =>
  ∫ x, jet (k + 5 + 2) φ x (Fin.cons a (Fin.cons c
      (Fin.append I ![d, b] : Fin (k + 3 + 2) → Fin n))) *
    locWeight ρ (lm_F (N := N) k I a b c d) φ x

/-- The main part `∫ φ_L ρ² γ M` of the term. -/
noncomputable def lm_Mf (k : ℕ) (ρ φ : EuclideanSpace ℝ (Fin n) → ℝ) :
    (Fin (k + 3) → Fin n) → Fin n → Fin n → Fin n → Fin n → ℝ := fun I a b c d =>
  ∫ x, jet (k + 5) φ x (Fin.append I ![d, b] : Fin (k + 3 + 2) → Fin n) *
    (ρ x ^ 2 * (gg a b c d (jp N φ x) * lm_M k a c I (jp N φ x)))

theorem lm_hLHS {k : ℕ} (hN : k + 7 ≤ N) {K : ℝ} {Λ : ℝ} {ρ : EuclideanSpace ℝ (Fin n) → ℝ}
    (hρ : IsCutoff ρ Λ) {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ)
    (hsupp : ∀ x ∈ tsupport ρ, jp N φ x ∈ admSet n N K) :
    (∫ x, ρ x ^ 2 * (2 * ∑ I : Fin (k + 3) → Fin n,
      topPot (k + 3) (jp N φ x) I * wSym k (jp N φ x) I)) =
      -2 * lm_S5 k (lm_Tf (N := N) k ρ φ) := by
  rw [integral_congr_ae (Filter.Eventually.of_forall (lm_lhs_pointwise hN hφ ρ))]
  have := lm_integral_S5 (k := k) (F := fun I a b c d x => (-2 : ℝ) *
    (jet (k + 5 + 2) φ x (Fin.cons a (Fin.cons c
      (Fin.append I ![d, b] : Fin (k + 3 + 2) → Fin n))) *
      locWeight ρ (lm_F (N := N) k I a b c d) φ x))
    (fun I a b c d => (lm_integrable_sym (lm_isSymbol_F (K := K) k (by omega) I a b c d) hρ hφ
      hsupp (k + 5 + 2) _).const_mul (-2))
  refine this.trans ?_
  simp only [integral_const_mul]
  rw [lm_S5_const_mul]
  rfl

theorem lm_hM {k : ℕ} (hN : k + 5 ≤ N) {Λ : ℝ} {ρ : EuclideanSpace ℝ (Fin n) → ℝ}
    (hρ : IsCutoff ρ Λ) {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ) :
    lm_S5 k (lm_Mf (N := N) k ρ φ) = ∫ x, ρ x ^ 2 * (volFun (low23 (jp N φ x)).1 *
      gNormSq (ginvFun (low23 (jp N φ x)).1) (jet (k + 5) φ x)) := by
  have := lm_integral_S5 (k := k) (F := fun I a b c d x =>
    jet (k + 5) φ x (Fin.append I ![d, b] : Fin (k + 3 + 2) → Fin n) *
      (ρ x ^ 2 * (gg a b c d (jp N φ x) * lm_M k a c I (jp N φ x))))
    (fun I a b c d => lm_integrable_smooth (lm_contDiff_gM (N := N) k a b c d I) hρ hφ (k + 5) _)
  refine this.symm.trans ?_
  exact integral_congr_ae (Filter.Eventually.of_forall fun x =>
    lm_main_pointwise hN hφ (ρ x ^ 2) x)

theorem lm_continuous_Q {k : ℕ} {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ) :
    Continuous (fun x => volFun (low23 (jp N φ x)).1 *
      gNormSq (ginvFun (low23 (jp N φ x)).1) (jet (k + 5) φ x)) := by
  have hu0 : ContDiff ℝ ∞ (fun q : P n N => (low23 q).1) := contDiff_fst.comp low23.contDiff
  have hu : Continuous fun x => (low23 (jp N φ x)).1 :=
    hu0.continuous.comp (contDiff_jp hφ).continuous
  have hg : Continuous (fun q : Z n 2 × Z n (k + 5) => gNormSq q.1 q.2) :=
    (contDiff_gNormSq (n := n) (m := k + 5)).continuous
  have hp : Continuous fun x => (ginvFun (low23 (jp N φ x)).1, jet (k + 5) φ x) :=
    ((contDiff_ginvFun (n := n)).continuous.comp hu).prodMk (contDiff_jet hφ (k + 5)).continuous
  have h5 := hg.comp hp
  have h6 : Continuous (fun x => gNormSq (ginvFun (low23 (jp N φ x)).1) (jet (k + 5) φ x)) :=
    h5.congr fun x => rfl
  have h7 : Continuous (fun x => volFun (low23 (jp N φ x)).1) :=
    (contDiff_volFun (n := n)).continuous.comp hu
  exact h7.mul h6

theorem lm_hlow {k : ℕ} (hN : k + 7 ≤ N) {K : ℝ} {Λ : ℝ} {ρ : EuclideanSpace ℝ (Fin n) → ℝ}
    (hρ : IsCutoff ρ Λ) {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ)
    (hsupp : ∀ x ∈ tsupport ρ, jp N φ x ∈ admSet n N K) :
    ((1 + K ^ 2)⁻¹) ^ (k + 5) * (∫ x, ρ x ^ 2 * ‖jet (k + 5) φ x‖ ^ 2) ≤
      lm_S5 k (lm_Mf (N := N) k ρ φ) := by
  have hQint : Integrable (fun x => ρ x ^ 2 * (volFun (low23 (jp N φ x)).1 *
      gNormSq (ginvFun (low23 (jp N φ x)).1) (jet (k + 5) φ x))) :=
    lm_integrable_rho hρ (lm_continuous_Q hφ)
  have hI1 : Integrable (fun x => ρ x ^ 2 * ‖jet (k + 5) φ x‖ ^ 2) :=
    lm_integrable_rho hρ ((contDiff_jet hφ (k + 5)).continuous.norm.pow 2)
  rw [lm_hM (by omega) hρ hφ, ← integral_const_mul]
  refine integral_mono (hI1.const_mul _) hQint fun x => ?_
  by_cases h : ρ x = 0
  · simp [h]
  · have hx : jp N φ x ∈ admSet n N K :=
      hsupp x (subset_tsupport ρ (Function.mem_support.2 h))
    have := lm_coercive (K := K) (by omega) (k + 5) hx (jet (k + 5) φ x)
    calc ((1 + K ^ 2)⁻¹) ^ (k + 5) * (ρ x ^ 2 * ‖jet (k + 5) φ x‖ ^ 2)
        = ρ x ^ 2 * (((1 + K ^ 2)⁻¹) ^ (k + 5) * ‖jet (k + 5) φ x‖ ^ 2) := by ring
      _ ≤ _ := mul_le_mul_of_nonneg_left this (sq_nonneg _)

theorem exists_const_main_le {N k : ℕ} (hN : k + 7 ≤ N) {K : ℝ} (hK : 0 ≤ K) (Λ : ℝ) {ε : ℝ}
    (hε : 0 < ε) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ ρ : EuclideanSpace ℝ (Fin n) → ℝ, IsCutoff ρ Λ →
      ∀ φ : EuclideanSpace ℝ (Fin n) → ℝ, ContDiff ℝ ∞ φ →
      (∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 3, jp N φ x ∈ admSet n N K) →
      ∫ x, ρ x ^ 2 * (2 * ∑ I : Fin (k + 3) → Fin n,
          topPot (k + 3) (jp N φ x) I * wSym k (jp N φ x) I) ≤
        -(2 * ((1 + K ^ 2)⁻¹) ^ (k + 5)) * (∫ x, ρ x ^ 2 * ‖jet (k + 5) φ x‖ ^ 2)
        + ε * (∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 2, ‖jet (k + 5) φ x‖ ^ 2)
        + C * (1 + energy3 (k + 1) 3 φ) := by
  classical
  obtain ⟨cd, hcd⟩ : ∃ cd : ℝ, cd = lm_S5 (n := n) k (fun _ _ _ _ _ => (1 : ℝ)) := ⟨_, rfl⟩
  have hcd0 : 0 ≤ cd := by
    rw [hcd]
    exact lm_S5_nonneg fun _ _ _ _ _ => zero_le_one
  obtain ⟨ε₀, hε₀⟩ : ∃ ε₀ : ℝ, ε₀ = ε / (2 * (cd + 1)) := ⟨_, rfl⟩
  have hε₀pos : 0 < ε₀ := by rw [hε₀]; positivity
  have h3 : 2 * cd * ε₀ ≤ ε := by
    have hpos : 0 < cd + 1 := by linarith
    rw [hε₀, show 2 * cd * (ε / (2 * (cd + 1))) = ε * (cd / (cd + 1)) by field_simp]
    calc ε * (cd / (cd + 1)) ≤ ε * 1 :=
          mul_le_mul_of_nonneg_left (div_le_one_of_le₀ (by linarith) hpos.le) hε.le
      _ = ε := mul_one _
  choose Cf hCf0 hCf using fun (I : Fin (k + 3) → Fin n) (a b c d : Fin n) =>
    lm_term_bound (N := N) hN (K := K) Λ hε₀pos I a b c d
  refine ⟨2 * lm_S5 k Cf, mul_nonneg zero_le_two (lm_S5_nonneg hCf0), fun ρ hρ φ hφ hadm => ?_⟩
  have hsupp : ∀ x ∈ tsupport ρ, jp N φ x ∈ admSet n N K :=
    fun x hx => hadm x (ibb_ball_sub x (hρ.supp hx))
  have hB0 : 0 ≤ ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 2, ‖jet (k + 5) φ x‖ ^ 2 :=
    setIntegral_nonneg measurableSet_ball fun x _ => by positivity
  obtain ⟨B, hB⟩ : ∃ B : ℝ, B = ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 2,
      ‖jet (k + 5) φ x‖ ^ 2 := ⟨_, rfl⟩
  rw [← hB] at hB0 ⊢
  have hE0 : 0 ≤ 1 + energy3 (k + 1) 3 φ := by
    have := energy3_nonneg (n := n) (k + 1) 3 φ
    linarith
  obtain ⟨E, hE⟩ : ∃ E : ℝ, E = 1 + energy3 (k + 1) 3 φ := ⟨_, rfl⟩
  rw [← hE] at hE0 ⊢
  have hdiff : lm_S5 k (lm_Mf (N := N) k ρ φ) - lm_S5 k (lm_Tf (N := N) k ρ φ) ≤
      cd * (ε₀ * B) + lm_S5 k Cf * E := by
    rw [← lm_S5_sub]
    have h1 : lm_S5 k (fun I a b c d => lm_Mf (N := N) k ρ φ I a b c d -
        lm_Tf (N := N) k ρ φ I a b c d) ≤
        lm_S5 k (fun I a b c d => ε₀ * B + Cf I a b c d * E) := by
      refine lm_S5_mono fun I a b c d => ?_
      have := hCf I a b c d ρ hρ φ hφ hadm
      rw [← hB, ← hE] at this
      have h2 : |lm_Tf (N := N) k ρ φ I a b c d - lm_Mf (N := N) k ρ φ I a b c d| ≤ _ := this
      rw [abs_sub_comm] at h2
      exact (le_abs_self _).trans h2
    refine h1.trans (le_of_eq ?_)
    rw [lm_S5_add, lm_S5_const, lm_S5_mul_const, ← hcd]
  have hfin := mul_le_mul_of_nonneg_right h3 hB0
  have hlow := lm_hlow hN hρ hφ hsupp
  rw [lm_hLHS hN hρ hφ hsupp]
  nlinarith [hlow, hdiff, hfin]

/-- info: 'CWFlow.Extension.exists_const_main_le' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_const_main_le

end CWFlow.Extension
