import CWFlow.Extension.IBPBounds
import CWFlow.Extension.WeightTop

/-!
# The almost top order terms of the energy inequality

Four groups of terms in the derivative of the energy, each of the form
`∫ ρ² · 2 Σ_I (…)_I w_I` with `w = wSym k` the weight tensor:

* (G2) `rPot (k+2)`,
* (G3) `Dz (jetSym potSym (k+2) − topPot (k+2))`,
* (G4) `cTop k (topPot (k+2))`,
* (G5) `cTop k (jetSym potSym (k+2) − topPot (k+2))`.

Each is bounded by `ε ∫_{B₂} |D^{k+5}φ|² + ε · energy + C (1 + lower energy)`.  (G5) needs no integration by
parts, (G2) and (G4) need one (the top jet coordinate gets a derivative moved onto the weight),
(G3) is one integration by parts moving the derivative from the remainder onto the weight.
`exists_const_lower_le` combines the four.
-/

open CWFlow.Interp Metric MeasureTheory
open scoped ContDiff

namespace CWFlow.Extension

variable {n N : ℕ}

/-! ### Part 0: generic bookkeeping -/

theorem ll_J_nonneg (m : ℕ) (φ : EuclideanSpace ℝ (Fin n) → ℝ) :
    0 ≤ ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 2, ‖jet m φ x‖ ^ 2 :=
  setIntegral_nonneg measurableSet_ball fun x _ => by positivity

/-- Finite sums of bounds. -/
theorem ll_sum_bound {ι Ω : Type*} [Fintype ι] {S : Set Ω} (F : ι → Ω → ℝ) (X Y : Ω → ℝ)
    (hX : ∀ ξ ∈ S, 0 ≤ X ξ) {ε ε' : ℝ} (hcard : ε' * (Fintype.card ι : ℝ) ≤ ε)
    (h : ∀ i, ∃ C : ℝ, 0 ≤ C ∧ ∀ ξ ∈ S, |F i ξ| ≤ ε' * X ξ + C * Y ξ) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ ξ ∈ S, |∑ i, F i ξ| ≤ ε * X ξ + C * Y ξ := by
  choose C hC0 hC using h
  refine ⟨∑ i, C i, Finset.sum_nonneg fun i _ => hC0 i, fun ξ hξ => ?_⟩
  calc |∑ i, F i ξ| ≤ ∑ i, |F i ξ| := Finset.abs_sum_le_sum_abs _ _
    _ ≤ ∑ i, (ε' * X ξ + C i * Y ξ) := Finset.sum_le_sum fun i _ => hC i ξ hξ
    _ = (Fintype.card ι : ℝ) * (ε' * X ξ) + (∑ i, C i) * Y ξ := by
        rw [Finset.sum_add_distrib, Finset.sum_const, Finset.card_univ, nsmul_eq_mul,
          Finset.sum_mul]
    _ ≤ ε * X ξ + (∑ i, C i) * Y ξ := by
        have := mul_le_mul_of_nonneg_right hcard (hX ξ hξ)
        nlinarith

theorem ll_abs4 (a b c d : ℝ) : |a + b + c + d| ≤ |a| + |b| + |c| + |d| := by
  calc |a + b + c + d| ≤ |a + b + c| + |d| := abs_add_le _ _
    _ ≤ |a| + |b| + |c| + |d| := by linarith [abs_add_three a b c]

theorem ll_fin_odd {I J E1 E2 C ε : ℝ} (hε : 0 < ε) (hJ : 0 ≤ J) (hE1 : 0 ≤ E1) (hC : 0 ≤ C)
    (h : |I| ≤ ε * E2 + C * 1) : |I| ≤ ε * J + ε * E2 + C * (1 + E1) := by
  nlinarith [mul_nonneg hε.le hJ, mul_nonneg hC hE1]

theorem ll_fin_jet {I J E1 E2 C ε : ℝ} (hε : 0 < ε) (hE2 : 0 ≤ E2) (h : |I| ≤ ε * J + C * (1 + E1)) :
    |I| ≤ ε * J + ε * E2 + C * (1 + E1) := by
  nlinarith [mul_nonneg hε.le hE2]

/-- The admissible data: a cutoff and a smooth function with admissible jets on `B₃`. -/
def LLAdm (N : ℕ) (K Λ : ℝ) (ξ : (EuclideanSpace ℝ (Fin n) → ℝ) × (EuclideanSpace ℝ (Fin n) → ℝ)) :
    Prop :=
  IsCutoff ξ.1 Λ ∧ ContDiff ℝ ∞ ξ.2 ∧
    ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 3, jp N ξ.2 x ∈ admSet n N K

/-- Sum of the integrals of `jet(k+6) · locWeight`, one integration by parts each. -/
theorem ll_sum_jet {N k : ℕ} (hN : k + 5 ≤ N) {K : ℝ} (Λ : ℝ) {ι : Type*} [Fintype ι]
    (𝔣 : ι → P n N → ℝ) (h𝔣 : ∀ i, IsSymbol (admSet n N K) k k (𝔣 i))
    (L : ι → Fin (k + 5 + 1) → Fin n) {ε : ℝ} (hε : 0 < ε) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ ρ : EuclideanSpace ℝ (Fin n) → ℝ, IsCutoff ρ Λ →
      ∀ φ : EuclideanSpace ℝ (Fin n) → ℝ, ContDiff ℝ ∞ φ →
      (∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 3, jp N φ x ∈ admSet n N K) →
      |∫ x, ρ x ^ 2 * ∑ i, (jet (k + 5 + 1) φ x (L i) * 𝔣 i (jp N φ x))| ≤
        ε * (∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 2, ‖jet (k + 5) φ x‖ ^ 2) +
          C * (1 + energy3 (k + 1) 3 φ) := by
  classical
  set ε' : ℝ := ε / ((Fintype.card ι : ℝ) + 1) with hε'
  have hε'0 : 0 < ε' := by positivity
  have hcard : ε' * (Fintype.card ι : ℝ) ≤ ε := by
    rw [hε', div_mul_eq_mul_div, div_le_iff₀ (by positivity)]
    nlinarith [Nat.cast_nonneg (α := ℝ) (Fintype.card ι)]
  have hex : ∀ i, ∃ C : ℝ, 0 ≤ C ∧ ∀ ξ ∈ {ξ | LLAdm N K Λ ξ},
      |∫ x, jet (k + 5 + 1) ξ.2 x (L i) * locWeight ξ.1 (𝔣 i) ξ.2 x| ≤
        ε' * (∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 2, ‖jet (k + 5) ξ.2 x‖ ^ 2) +
          C * (1 + energy3 (k + 1) 3 ξ.2) := by
    intro i
    obtain ⟨C, hC0, hC⟩ := exists_const_integral_jet_locWeight hN Λ (h𝔣 i) hε'0
    exact ⟨C, hC0, fun ξ hξ => hC ξ.1 hξ.1 ξ.2 hξ.2.1 hξ.2.2 (L i)⟩
  obtain ⟨C, hC0, hC⟩ := ll_sum_bound (S := {ξ | LLAdm N K Λ ξ})
    (fun i ξ => ∫ x, jet (k + 5 + 1) ξ.2 x (L i) * locWeight ξ.1 (𝔣 i) ξ.2 x)
    (fun ξ => ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 2, ‖jet (k + 5) ξ.2 x‖ ^ 2)
    (fun ξ => 1 + energy3 (k + 1) 3 ξ.2) (fun ξ _ => ll_J_nonneg _ _) hcard hex
  refine ⟨C, hC0, fun ρ hρ φ hφ hadm => ?_⟩
  have hsupp : ∀ x ∈ tsupport ρ, jp N φ x ∈ admSet n N K :=
    fun x hx => hadm x (ibb_ball_sub x (hρ.supp hx))
  have hint : ∀ i ∈ (Finset.univ : Finset ι),
      Integrable (fun x => jet (k + 5 + 1) φ x (L i) * locWeight ρ (𝔣 i) φ x) := by
    intro i _
    exact ((ibb_coord_contDiff hφ _ (L i)).continuous.mul
      (contDiff_locWeight (h𝔣 i) hφ hρ.smooth hsupp).continuous).integrable_of_hasCompactSupport
      (hasCompactSupport_locWeight (ibb_hasCompactSupport hρ) (𝔣 i) φ).mul_left
  have e : ∫ x, ρ x ^ 2 * ∑ i, (jet (k + 5 + 1) φ x (L i) * 𝔣 i (jp N φ x)) =
      ∑ i, ∫ x, jet (k + 5 + 1) φ x (L i) * locWeight ρ (𝔣 i) φ x := by
    rw [← integral_finsetSum _ hint]
    congr 1
    funext x
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl fun i _ => ?_
    simp only [locWeight]
    ring
  rw [e]
  exact hC (ρ, φ) ⟨hρ, hφ, hadm⟩

/-- From a pointwise bound by `Om (2k+3) (k+2)` to an `ε`-bound by the energy. -/
theorem ll_odd_of_pt {k : ℕ} (K : ℝ) {ε : ℝ} (hε : 0 < ε) (C₀ : ℝ) (hC₀ : 0 ≤ C₀) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ φ : EuclideanSpace ℝ (Fin n) → ℝ, ContDiff ℝ ∞ φ →
      (∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 3, ‖jet 3 φ x‖ ≤ K) →
      ∀ h : EuclideanSpace ℝ (Fin n) → ℝ, Continuous h →
      (∀ x, x ∉ ball (0 : EuclideanSpace ℝ (Fin n)) 2 → h x = 0) →
      (∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 2,
        |h x| ≤ C₀ * Om (2 * k + 3) (k + 2) (aJet φ x)) →
      |∫ x, h x| ≤ ε * energy3 (k + 2) 3 φ + C := by
  obtain ⟨C', hC', hint⟩ := exists_const_integral_Om_odd (n := n) (m := k + 2) (by omega)
    (ρ := 2) (R := 3) (by norm_num) (by norm_num) K (ε := ε / (C₀ + 1)) (by positivity)
  refine ⟨C₀ * C', mul_nonneg hC₀ hC', fun φ hφ hK h hcont hz hb => ?_⟩
  have hom : ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 2, Om (2 * k + 3) (k + 2) (aJet φ x) ≤
      ε / (C₀ + 1) * energy3 (k + 2) 3 φ + C' := by
    have := hint φ hφ hK
    have e : 2 * (k + 2) - 1 = 2 * k + 3 := by omega
    rwa [e] at this
  have hq : C₀ * (ε / (C₀ + 1)) ≤ ε := by
    calc C₀ * (ε / (C₀ + 1)) = ε * (C₀ / (C₀ + 1)) := by ring
      _ ≤ ε * 1 := mul_le_mul_of_nonneg_left (div_le_one_of_le₀ (by linarith) (by positivity)) hε.le
      _ = ε := mul_one _
  rw [← setIntegral_eq_integral_of_forall_compl_eq_zero hz]
  calc |∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 2, h x|
      ≤ ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 2, |h x| := by
        simpa [Real.norm_eq_abs] using norm_integral_le_integral_norm (μ := volume.restrict
          (ball (0 : EuclideanSpace ℝ (Fin n)) 2)) h
    _ ≤ ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 2, C₀ * Om (2 * k + 3) (k + 2) (aJet φ x) :=
        setIntegral_mono_on (integrableOn_ball_of_continuous hcont 2).abs
          ((integrableOn_ball_of_continuous (continuous_Om hφ _ _) 2).const_mul _)
          measurableSet_ball hb
    _ = C₀ * ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 2, Om (2 * k + 3) (k + 2) (aJet φ x) :=
        integral_const_mul _ _
    _ ≤ C₀ * (ε / (C₀ + 1) * energy3 (k + 2) 3 φ + C') := mul_le_mul_of_nonneg_left hom hC₀
    _ ≤ ε * energy3 (k + 2) 3 φ + C₀ * C' := by
        nlinarith [mul_le_mul_of_nonneg_right hq (energy3_nonneg (n := n) (k + 2) 3 φ)]

/-- The remainder `jet(k+2)` of the potential minus its top part. -/
@[reducible] noncomputable def ll_Rm (n N k : ℕ) (q : P n N) : Z n (k + 2) :=
  jetSym potSym (k + 2) q - topPot (k + 2) q

/-! ### Part 1: the coefficient of the `(k+2)`-slot in `cTop` -/

/-- The inclusion of the `(k+2)`-jet slot into `P n N`. -/
noncomputable def ll_injZ (n N k : ℕ) (h2 : k + 2 < N + 1) : Z n (k + 2) →L[ℝ] P n N :=
  injSlot n N ⟨k + 2, h2⟩

/-- `cTop k q` as a continuous linear map. -/
noncomputable def ll_cL (k : ℕ) (h2 : k + 2 < N + 1) (q : P n N) : Z n (k + 2) →L[ℝ] Z n (k + 3) :=
  (Dsl (some (⟨k + 2, h2⟩ : Fin (N + 1))) (sSym (N := N) (n := n) k) q).comp (ll_injZ n N k h2)

theorem ll_cTop_eq (k : ℕ) (h2 : k + 2 < N + 1) (q : P n N) (z : Z n (k + 2)) :
    cTop k q z = ll_cL k h2 q z := by
  have e : (fun m : Fin (N + 1) => if h : (m : ℕ) = k + 2 then h ▸ z else (0 : Z n m)) =
      (Pi.single (⟨k + 2, h2⟩ : Fin (N + 1)) z : J n N) := by
    funext m
    by_cases h : (m : ℕ) = k + 2
    · obtain rfl : m = ⟨k + 2, h2⟩ := Fin.ext h
      simp
    · have hm : m ≠ ⟨k + 2, h2⟩ := fun e => h (by rw [e])
      rw [dif_neg h, Pi.single_eq_of_ne hm]
  have e2 : iotaP (some (⟨k + 2, h2⟩ : Fin (N + 1))) (ll_injZ n N k h2 z) =
      ((0 : EuclideanSpace ℝ (Fin n)),
        (Pi.single (⟨k + 2, h2⟩ : Fin (N + 1)) z : J n N)) := by
    have h1 : ((ll_injZ n N k h2) z).2 ⟨k + 2, h2⟩ = z := by
      unfold ll_injZ
      rw [injSlot_apply]
      exact Pi.single_eq_same _ _
    rw [iotaP_some_apply, h1]
  unfold cTop ll_cL
  rw [ContinuousLinearMap.comp_apply]
  show _ = (fderiv ℝ (sSym (N := N) (n := n) k) q).comp (iotaP (some (⟨k + 2, h2⟩ : Fin (N + 1))))
    (ll_injZ n N k h2 z)
  rw [ContinuousLinearMap.comp_apply, e2, e]

section Sym

variable {A : Set (P n N)}

theorem ll_isSymbol_cL (hA : IsOpen A)
    (hlow : ∀ m : Fin (N + 1), (m : ℕ) ≤ 3 → ∃ K : ℝ, ∀ p ∈ A, ‖p.2 m‖ ≤ K)
    (k : ℕ) (hk : 1 ≤ k) (hN : k + 3 ≤ N) (h2 : k + 2 < N + 1) :
    IsSymbol A 0 (k - 1) (fun q : P n N => ll_cL k h2 q) := by
  have hS := isSymbol_sSym hA hlow k hk hN
  have hD := hS.Dsl (some (⟨k + 2, h2⟩ : Fin (N + 1))) (by simp [wt]; try omega)
  have hD0 := hD.mono (w' := 0) (by simp [wt]; try omega) le_rfl
  exact hD0.clm_comp ((ContinuousLinearMap.compL ℝ (Z n (k + 2)) (P n N) (Z n (k + 3))).flip
    (ll_injZ n N k h2))

theorem ll_isSymbol_cL_apply (hA : IsOpen A)
    (hlow : ∀ m : Fin (N + 1), (m : ℕ) ≤ 3 → ∃ K : ℝ, ∀ p ∈ A, ‖p.2 m‖ ≤ K)
    (k : ℕ) (hk : 1 ≤ k) (hN : k + 3 ≤ N) (h2 : k + 2 < N + 1) (z : Z n (k + 2))
    (I : Fin (k + 3) → Fin n) :
    IsSymbol A 0 (k - 1) (fun q : P n N => ll_cL k h2 q z I) :=
  ((ll_isSymbol_cL hA hlow k hk hN h2).clm_comp (ContinuousLinearMap.apply ℝ (Z n (k + 3)) z)).clm_comp
    (PiLp.proj (𝕜 := ℝ) 2 (fun _ : Fin (k + 3) → Fin n => ℝ) I)

theorem ll_isSymbol_Rm (hA : IsOpen A)
    (hlow : ∀ m : Fin (N + 1), (m : ℕ) ≤ 3 → ∃ K : ℝ, ∀ p ∈ A, ‖p.2 m‖ ≤ K)
    (k : ℕ) (hN : k + 6 ≤ N) : IsSymbol A (k + 2) (k + 2) (ll_Rm n N k) :=
  isSymbol_jetSym_potSym_sub hA hlow (k + 2) (by omega)

theorem ll_isSymbol_cL_Rm (hA : IsOpen A)
    (hlow : ∀ m : Fin (N + 1), (m : ℕ) ≤ 3 → ∃ K : ℝ, ∀ p ∈ A, ‖p.2 m‖ ≤ K)
    (k : ℕ) (hk : 1 ≤ k) (hN : k + 6 ≤ N) (h2 : k + 2 < N + 1) (I : Fin (k + 3) → Fin n) :
    IsSymbol A (0 + (k + 2)) (max (k - 1) (k + 2))
      (fun q : P n N => ll_cL k h2 q (ll_Rm n N k q) I) := by
  have hb := IsSymbol.bilin
    (ContinuousLinearMap.apply ℝ (Z n (k + 3)) :
      Z n (k + 2) →L[ℝ] (Z n (k + 2) →L[ℝ] Z n (k + 3)) →L[ℝ] Z n (k + 3)).flip
    (ll_isSymbol_cL hA hlow k hk (by omega) h2) (ll_isSymbol_Rm hA hlow k hN)
  exact hb.clm_comp (PiLp.proj (𝕜 := ℝ) 2 (fun _ : Fin (k + 3) → Fin n => ℝ) I)

/-- The scalar symbol of group (G5). -/
theorem ll_isSymbol_g5 (hA : IsOpen A)
    (hlow : ∀ m : Fin (N + 1), (m : ℕ) ≤ 3 → ∃ K : ℝ, ∀ p ∈ A, ‖p.2 m‖ ≤ K)
    (k : ℕ) (hk : 1 ≤ k) (hN : k + 6 ≤ N) :
    IsSymbol A (2 * k + 3) (k + 2) (fun q : P n N => 2 * ∑ I : Fin (k + 3) → Fin n,
      cTop k q (jetSym potSym (k + 2) q - topPot (k + 2) q) I * wSym k q I) := by
  have h2 : k + 2 < N + 1 := by omega
  have hterm : ∀ I : Fin (k + 3) → Fin n, IsSymbol A (2 * k + 3) (k + 2)
      (fun q : P n N => ll_cL k h2 q (ll_Rm n N k q) I * wSym k q I) := fun I =>
    ((ll_isSymbol_cL_Rm hA hlow k hk hN h2 I).mul
      (isSymbol_wSym_apply hA hlow k (by omega) I)).mono (by omega) (by omega)
  have hs := IsSymbol.sum hA (Finset.univ : Finset (Fin (k + 3) → Fin n)) (fun I _ => hterm I)
  have h3 := ((isSymbol_const hA (2 : ℝ)).mul hs).mono (w' := 2 * k + 3) (ℓ' := k + 2)
    (by omega) (by omega)
  refine h3.congr fun p _ => ?_
  simp only [ll_cTop_eq k h2 p]

end Sym

/-! ### Part 2: group (G5) -/

theorem exists_const_group5_le {N k : ℕ} (hk : 1 ≤ k) (hN : k + 7 ≤ N) {K : ℝ} (Λ : ℝ) {ε : ℝ}
    (hε : 0 < ε) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ ρ : EuclideanSpace ℝ (Fin n) → ℝ, IsCutoff ρ Λ →
      ∀ φ : EuclideanSpace ℝ (Fin n) → ℝ, ContDiff ℝ ∞ φ →
      (∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 3, jp N φ x ∈ admSet n N K) →
      |∫ x, ρ x ^ 2 * (2 * ∑ I : Fin (k + 3) → Fin n,
          cTop k (jp N φ x) (jetSym potSym (k + 2) (jp N φ x) - topPot (k + 2) (jp N φ x)) I *
            wSym k (jp N φ x) I)| ≤
        ε * (∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 2, ‖jet (k + 5) φ x‖ ^ 2)
        + ε * energy3 (k + 2) 3 φ + C * (1 + energy3 (k + 1) 3 φ) := by
  have hsym := ll_isSymbol_g5 (isOpen_admSet (n := n) (N := N) K) (admSet_low K) k hk (by omega)
  obtain ⟨C, hC0, hC⟩ := exists_const_integral_locWeight_odd (by omega : k + 5 ≤ N) hsym hε
  refine ⟨C, hC0, fun ρ hρ φ hφ hadm => ?_⟩
  have h := hC Λ ρ hρ φ hφ hadm
  exact ll_fin_odd hε (ll_J_nonneg _ φ) (energy3_nonneg _ _ φ) hC0 (by rw [mul_one]; exact h)

/-! ### Part 3: group (G2) -/

theorem exists_const_group2_le {N k : ℕ} (hN : k + 7 ≤ N) {K : ℝ} (Λ : ℝ) {ε : ℝ}
    (hε : 0 < ε) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ ρ : EuclideanSpace ℝ (Fin n) → ℝ, IsCutoff ρ Λ →
      ∀ φ : EuclideanSpace ℝ (Fin n) → ℝ, ContDiff ℝ ∞ φ →
      (∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 3, jp N φ x ∈ admSet n N K) →
      |∫ x, ρ x ^ 2 * (2 * ∑ I : Fin (k + 3) → Fin n,
          rPot (k + 2) (jp N φ x) I * wSym k (jp N φ x) I)| ≤
        ε * (∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 2, ‖jet (k + 5) φ x‖ ^ 2)
        + ε * energy3 (k + 2) 3 φ + C * (1 + energy3 (k + 1) 3 φ) := by
  have hA := isOpen_admSet (n := n) (N := N) K
  have hlow := admSet_low (n := n) (N := N) K
  obtain ⟨C, hC0, hC⟩ := ll_sum_jet (ι := (Fin (k + 3) → Fin n) × Fin n × Fin n × Fin n × Fin n)
    (by omega : k + 5 ≤ N) Λ
    (fun i q => -2 * (totalD (gg i.2.1 i.2.2.1 i.2.2.2.1 i.2.2.2.2) q
      (EuclideanSpace.single (i.1 0) 1) * wSym k q i.1))
    (fun i => ((isSymbol_const hA (-2 : ℝ)).mul
      ((isSymbol_totalD_gg hA hlow i.2.1 i.2.2.1 i.2.2.2.1 i.2.2.2.2 (i.1 0)).mul
        (isSymbol_wSym_apply hA hlow k (by omega) i.1))).mono (by omega) (by omega))
    (fun i => (Fin.append (Fin.tail i.1)
      ![i.2.2.2.2, i.2.1, i.2.2.1, i.2.2.2.1] : Fin (k + 2 + 4) → Fin n)) hε
  refine ⟨C, hC0, fun ρ hρ φ hφ hadm => ?_⟩
  have h := hC ρ hρ φ hφ hadm
  have e : ∀ x, ρ x ^ 2 * (2 * ∑ I : Fin (k + 3) → Fin n,
          rPot (k + 2) (jp N φ x) I * wSym k (jp N φ x) I) =
      ρ x ^ 2 * ∑ i : (Fin (k + 3) → Fin n) × Fin n × Fin n × Fin n × Fin n,
        (jet (k + 5 + 1) φ x ((Fin.append (Fin.tail i.1)
          ![i.2.2.2.2, i.2.1, i.2.2.1, i.2.2.2.1] : Fin (k + 2 + 4) → Fin n)) *
          (-2 * (totalD (gg i.2.1 i.2.2.1 i.2.2.2.1 i.2.2.2.2) (jp N φ x)
            (EuclideanSpace.single (i.1 0) 1) * wSym k (jp N φ x) i.1))) := by
    intro x
    congr 1
    simp only [Fintype.sum_prod_type, Finset.mul_sum, rPot_apply]
    refine Finset.sum_congr rfl fun I _ => ?_
    simp only [neg_mul, mul_neg, Finset.sum_mul, Finset.mul_sum, ← Finset.sum_neg_distrib]
    refine Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun b _ =>
      Finset.sum_congr rfl fun c _ => Finset.sum_congr rfl fun d _ => ?_
    rw [slotN_jp (k + 2 + 4) (by omega)]
    ring
  have e2 := integral_congr_ae (μ := volume) (Filter.Eventually.of_forall e)
  rw [e2]
  exact ll_fin_jet hε (energy3_nonneg _ _ φ) h

/-! ### Part 4: group (G4) -/

/-- The coordinate basis vector of `Z n m`. -/
noncomputable def ll_e (n m : ℕ) (K' : Fin m → Fin n) : Z n m := WithLp.toLp 2 (Pi.single K' (1 : ℝ))

theorem ll_expand {m : ℕ} (z : Z n m) : z = ∑ K', z K' • ll_e n m K' := by
  ext J
  simp [ll_e, Pi.single_apply]

theorem ll_cL_topPot (k : ℕ) (h2 : k + 2 < N + 1) (q : P n N) (I : Fin (k + 3) → Fin n) :
    ll_cL k h2 q (topPot (k + 2) q) I =
      ∑ K' : Fin (k + 2) → Fin n, topPot (k + 2) q K' * ll_cL k h2 q (ll_e n (k + 2) K') I := by
  conv_lhs => rw [ll_expand (topPot (k + 2) q)]
  simp [map_sum, map_smul]

theorem exists_const_group4_le {N k : ℕ} (hk : 1 ≤ k) (hN : k + 7 ≤ N) {K : ℝ} (Λ : ℝ) {ε : ℝ}
    (hε : 0 < ε) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ ρ : EuclideanSpace ℝ (Fin n) → ℝ, IsCutoff ρ Λ →
      ∀ φ : EuclideanSpace ℝ (Fin n) → ℝ, ContDiff ℝ ∞ φ →
      (∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 3, jp N φ x ∈ admSet n N K) →
      |∫ x, ρ x ^ 2 * (2 * ∑ I : Fin (k + 3) → Fin n,
          cTop k (jp N φ x) (topPot (k + 2) (jp N φ x)) I * wSym k (jp N φ x) I)| ≤
        ε * (∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 2, ‖jet (k + 5) φ x‖ ^ 2)
        + ε * energy3 (k + 2) 3 φ + C * (1 + energy3 (k + 1) 3 φ) := by
  have hA := isOpen_admSet (n := n) (N := N) K
  have hlow := admSet_low (n := n) (N := N) K
  have h2 : k + 2 < N + 1 := by omega
  obtain ⟨C, hC0, hC⟩ := ll_sum_jet
    (ι := (Fin (k + 3) → Fin n) × (Fin (k + 2) → Fin n) × Fin n × Fin n × Fin n × Fin n)
    (by omega : k + 5 ≤ N) Λ
    (fun i q => -2 * (ll_cL k h2 q (ll_e n (k + 2) i.2.1) i.1 *
      (gg i.2.2.1 i.2.2.2.1 i.2.2.2.2.1 i.2.2.2.2.2 q * wSym k q i.1)))
    (fun i => ((isSymbol_const hA (-2 : ℝ)).mul
      ((ll_isSymbol_cL_apply hA hlow k hk (by omega) h2 (ll_e n (k + 2) i.2.1) i.1).mul
        ((isSymbol_gg hA hlow i.2.2.1 i.2.2.2.1 i.2.2.2.2.1 i.2.2.2.2.2).mul
          (isSymbol_wSym_apply hA hlow k (by omega) i.1)))).mono (by omega) (by omega))
    (fun i => (Fin.append i.2.1
      ![i.2.2.2.2.2, i.2.2.1, i.2.2.2.1, i.2.2.2.2.1] : Fin (k + 2 + 4) → Fin n)) hε
  refine ⟨C, hC0, fun ρ hρ φ hφ hadm => ?_⟩
  have h := hC ρ hρ φ hφ hadm
  have e : ∀ x, ρ x ^ 2 * (2 * ∑ I : Fin (k + 3) → Fin n,
          cTop k (jp N φ x) (topPot (k + 2) (jp N φ x)) I * wSym k (jp N φ x) I) =
      ρ x ^ 2 * ∑ i : (Fin (k + 3) → Fin n) × (Fin (k + 2) → Fin n) × Fin n × Fin n × Fin n × Fin n,
        (jet (k + 5 + 1) φ x ((Fin.append i.2.1
          ![i.2.2.2.2.2, i.2.2.1, i.2.2.2.1, i.2.2.2.2.1] : Fin (k + 2 + 4) → Fin n)) *
          (-2 * (ll_cL k h2 (jp N φ x) (ll_e n (k + 2) i.2.1) i.1 *
            (gg i.2.2.1 i.2.2.2.1 i.2.2.2.2.1 i.2.2.2.2.2 (jp N φ x) * wSym k (jp N φ x) i.1)))) := by
    intro x
    congr 1
    simp only [Fintype.sum_prod_type, Finset.mul_sum, Finset.sum_mul, ll_cTop_eq k h2,
      ll_cL_topPot, topPot_apply]
    refine Finset.sum_congr rfl fun I _ => Finset.sum_congr rfl fun K' _ => ?_
    simp only [neg_mul, mul_neg, Finset.sum_mul, Finset.mul_sum, ← Finset.sum_neg_distrib]
    refine Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun b _ =>
      Finset.sum_congr rfl fun c _ => Finset.sum_congr rfl fun d _ => ?_
    rw [slotN_jp (k + 2 + 4) (by omega)]
    ring
  have e2 := integral_congr_ae (μ := volume) (Filter.Eventually.of_forall e)
  rw [e2]
  exact ll_fin_jet hε (energy3_nonneg _ _ φ) h

/-! ### Part 5: group (G3) -/

theorem ll_isSymbol_g3_term {A : Set (P n N)} (hA : IsOpen A)
    (hlow : ∀ m : Fin (N + 1), (m : ℕ) ≤ 3 → ∃ K : ℝ, ∀ p ∈ A, ‖p.2 m‖ ≤ K)
    (k : ℕ) (hN : k + 6 ≤ N) (I : Fin (k + 3) → Fin n) :
    IsSymbol A (2 * k + 3) (k + 3)
      (fun q : P n N => 2 * (Dz (ll_Rm n N k) q I * wSym k q I)) := by
  have hRm := ll_isSymbol_Rm hA hlow k hN
  have hD := (hRm.Dz (by omega) hlow).clm_comp
    (PiLp.proj (𝕜 := ℝ) 2 (fun _ : Fin (k + 2 + 1) → Fin n => ℝ) I)
  have h1 := (hD.mul (isSymbol_wSym_apply hA hlow k (by omega) I)).mono
    (w' := 2 * k + 3) (ℓ' := k + 3) (by omega) (by omega)
  exact ((isSymbol_const hA (2 : ℝ)).mul h1).mono (by omega) (by omega)

/-- One coordinate of group (G3): one integration by parts moves the derivative onto the weight. -/
theorem ll_g3_term {N k : ℕ} (hN : k + 7 ≤ N) {K : ℝ} (Λ : ℝ) (I : Fin (k + 3) → Fin n) {ε : ℝ}
    (hε : 0 < ε) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ ρ : EuclideanSpace ℝ (Fin n) → ℝ, IsCutoff ρ Λ →
      ∀ φ : EuclideanSpace ℝ (Fin n) → ℝ, ContDiff ℝ ∞ φ →
      (∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 3, jp N φ x ∈ admSet n N K) →
      |∫ x, ρ x ^ 2 * (2 * (Dz (ll_Rm n N k) (jp N φ x) I * wSym k (jp N φ x) I))| ≤
        ε * energy3 (k + 2) 3 φ + C := by
  have hA := isOpen_admSet (n := n) (N := N) K
  have hlow := admSet_low (n := n) (N := N) K
  have hRm := ll_isSymbol_Rm hA hlow k (by omega : k + 6 ≤ N)
  have hFm : IsSymbol (admSet n N K) (k + 2) (k + 2)
      (fun q : P n N => ll_Rm n N k q (Fin.tail I)) :=
    hRm.clm_comp (PiLp.proj (𝕜 := ℝ) 2 (fun _ : Fin (k + 2) → Fin n => ℝ) (Fin.tail I))
  have hw : IsSymbol (admSet n N K) k k (fun q : P n N => wSym k q I) :=
    isSymbol_wSym_apply hA hlow k (by omega) I
  have hw2 : IsSymbol (admSet n N K) k k (fun q : P n N => 2 * wSym k q I) :=
    ((isSymbol_const hA (2 : ℝ)).mul hw).mono (by omega) (by omega)
  have hg : IsSymbol (admSet n N K) (2 * k + 2) (k + 2)
      (fun q : P n N => ll_Rm n N k q (Fin.tail I) * (2 * wSym k q I)) :=
    (hFm.mul hw2).mono (by omega) (by omega)
  have hf3 := ll_isSymbol_g3_term hA hlow k (by omega : k + 6 ≤ N) I
  obtain ⟨CR, hCR0, hCR⟩ := ibb_pt hFm (by omega : 3 + min (k + 2) (k + 2) ≤ N)
  obtain ⟨Cf, hCf0, hf⟩ := ibb_f_td (by omega : k + 4 ≤ N) hw2
  obtain ⟨C3, hC30, hC3⟩ := evp_Om_mul3 k (k + 2) (k + 1) 0 (k + 2) (k + 1) 0 (by omega) (by omega)
  have hΛ : 0 ≤ 2 * |Λ| + 1 := by positivity
  obtain ⟨C', hC', hodd⟩ := ll_odd_of_pt (n := n) (k := k) K hε
    (CR * ((2 * |Λ| + 1) * Cf) * C3) (by positivity)
  refine ⟨C', hC', fun ρ hρ φ hφ hadm => ?_⟩
  have hsupp : ∀ x ∈ tsupport ρ, jp N φ x ∈ admSet n N K :=
    fun x hx => hadm x (ibb_ball_sub x (hρ.supp hx))
  have hK : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 3, ‖jet 3 φ x‖ ≤ K :=
    fun x hx => ibb_norm3 (by omega) (hadm x hx)
  obtain ⟨W, hWdef⟩ : ∃ W : EuclideanSpace ℝ (Fin n) → ℝ,
      W = locWeight ρ (fun q : P n N => 2 * wSym k q I) φ := ⟨_, rfl⟩
  obtain ⟨F, hFdef⟩ : ∃ F : EuclideanSpace ℝ (Fin n) → ℝ,
      F = fun y => ll_Rm n N k (jp N φ y) (Fin.tail I) := ⟨_, rfl⟩
  have hW : ContDiff ℝ ∞ W := by
    rw [hWdef]; exact contDiff_locWeight hw2 hφ hρ.smooth hsupp
  have hcW : HasCompactSupport W := by
    rw [hWdef]; exact hasCompactSupport_locWeight (ibb_hasCompactSupport hρ) _ φ
  have hW0 : ∀ y, y ∉ tsupport ρ → fderiv ℝ W y = 0 := by
    intro y hy
    rw [hWdef]; exact fderiv_locWeight_eq_zero _ φ hy
  have hU : IsOpen {y | jp N φ y ∈ admSet n N K} := isOpen_jp_preimage hFm hφ
  have hFd : ∀ x ∈ tsupport ρ, DifferentiableAt ℝ F x := by
    intro x hx
    rw [hFdef]
    exact ((contDiffOn_symbol_jp hFm hφ).contDiffAt (hU.mem_nhds (hsupp x hx))).differentiableAt
      (by simp)
  have hts : tsupport W ⊆ tsupport ρ := by
    refine closure_minimal (fun x hx => ?_) (isClosed_tsupport ρ)
    by_contra hn
    apply hx
    rw [hWdef]
    exact locWeight_eq_zero _ φ hn
  have hid : ∀ x, ρ x ^ 2 * (2 * (Dz (ll_Rm n N k) (jp N φ x) I * wSym k (jp N φ x) I)) =
      fderiv ℝ F x (EuclideanSpace.single (I 0) 1) * W x := by
    intro x
    by_cases hx : x ∈ tsupport ρ
    · have := Dz_jp_apply hRm (by omega : k + 2 + 4 ≤ N) hφ (hsupp x hx) I
      rw [this, hWdef, hFdef]
      simp only [locWeight]
      ring
    · have h0 : ρ x = 0 := image_eq_zero_of_notMem_tsupport hx
      rw [hWdef]
      simp [locWeight, h0]
  have hint1 : Integrable (fun x => fderiv ℝ F x (EuclideanSpace.single (I 0) 1) * W x) := by
    have e : (fun x => fderiv ℝ F x (EuclideanSpace.single (I 0) 1) * W x) =
        locWeight ρ (fun q : P n N => 2 * (Dz (ll_Rm n N k) q I * wSym k q I)) φ :=
      funext fun x => (hid x).symm
    rw [e]
    exact (contDiff_locWeight hf3 hφ hρ.smooth hsupp).continuous.integrable_of_hasCompactSupport
      (hasCompactSupport_locWeight (ibb_hasCompactSupport hρ) _ φ)
  have hcont : Continuous (fun x => F x * fderiv ℝ W x (EuclideanSpace.single (I 0) 1)) := by
    refine continuous_iff_continuousAt.2 fun x => ?_
    by_cases hx : x ∈ tsupport ρ
    · exact (hFd x hx).continuousAt.mul
        (((hW.continuous_fderiv (by simp)).clm_apply continuous_const).continuousAt)
    · have : (fun y => F y * fderiv ℝ W y (EuclideanSpace.single (I 0) 1)) =ᶠ[nhds x]
          fun _ => 0 := by
        filter_upwards [(isClosed_tsupport ρ).isOpen_compl.mem_nhds hx] with y hy
        rw [hW0 y hy]
        simp
      exact continuousAt_const.congr this.symm
  have hint2 : Integrable (fun x => F x * fderiv ℝ W x (EuclideanSpace.single (I 0) 1)) :=
    hcont.integrable_of_hasCompactSupport
      ((hcW.fderiv_apply ℝ (EuclideanSpace.single (I 0) (1 : ℝ))).mul_left)
  have hint3 : Integrable (fun x => F x * W x) := by
    have e : (fun x => F x * W x) =
        locWeight ρ (fun q : P n N => ll_Rm n N k q (Fin.tail I) * (2 * wSym k q I)) φ := by
      funext x
      rw [hWdef, hFdef]
      simp only [locWeight]
      ring
    rw [e]
    exact (contDiff_locWeight hg hφ hρ.smooth hsupp).continuous.integrable_of_hasCompactSupport
      (hasCompactSupport_locWeight (ibb_hasCompactSupport hρ) _ φ)
  have hibp := integral_mul_fderiv_eq_neg_fderiv_mul_of_integrable (μ := volume) (f := F) (g := W)
    (v := EuclideanSpace.single (I 0) 1) hint1 hint2 hint3 (fun x hx => hFd x (hts hx))
    (fun x _ => hW.differentiable (by simp) x)
  have e2 : (∫ x, ρ x ^ 2 * (2 * (Dz (ll_Rm n N k) (jp N φ x) I * wSym k (jp N φ x) I))) =
      -∫ x, F x * fderiv ℝ W x (EuclideanSpace.single (I 0) 1) := by
    rw [hibp, neg_neg]
    exact integral_congr_ae (Filter.Eventually.of_forall hid)
  rw [e2, abs_neg]
  refine hodd φ hφ hK _ hcont (fun x hx => ?_) (fun x hx => ?_)
  · have hxs : x ∉ tsupport ρ := fun h => hx (hρ.supp h)
    simp [hW0 x hxs]
  · by_cases hxs : x ∈ tsupport ρ
    · have h1 := hCR φ x (hsupp x hxs)
      rw [Real.norm_eq_abs] at h1
      have h2 := ibb_wfd (by omega : k + 4 ≤ N) hw2 hf hρ hφ (hsupp x hxs) (I 0)
      have h3 := hC3 (aJet φ x) (aJet_nonneg φ x)
      rw [Om_zero, mul_one] at h3
      have hO1 : 0 ≤ Om (k + 2) (k + 2) (aJet φ x) := evp_Om_nonneg _ _ (aJet_nonneg φ x)
      have hO2 : 0 ≤ Om (k + 1) (k + 1) (aJet φ x) := evp_Om_nonneg _ _ (aJet_nonneg φ x)
      rw [← hWdef] at h2
      rw [abs_mul]
      have hFx : F x = ll_Rm n N k (jp N φ x) (Fin.tail I) := by rw [hFdef]
      rw [hFx]
      calc |ll_Rm n N k (jp N φ x) (Fin.tail I)| *
            |fderiv ℝ W x (EuclideanSpace.single (I 0) 1)|
          ≤ (CR * Om (k + 2) (k + 2) (aJet φ x)) *
              ((2 * |Λ| + 1) * Cf * Om (k + 1) (k + 1) (aJet φ x)) :=
            mul_le_mul h1 h2 (abs_nonneg _) (mul_nonneg hCR0 hO1)
        _ = CR * ((2 * |Λ| + 1) * Cf) *
              (Om (k + 2) (k + 2) (aJet φ x) * Om (k + 1) (k + 1) (aJet φ x)) := by ring
        _ ≤ CR * ((2 * |Λ| + 1) * Cf) * (C3 * Om (2 * k + 3) (k + 2) (aJet φ x)) :=
            mul_le_mul_of_nonneg_left h3 (by positivity)
        _ = _ := by ring
    · rw [hW0 x hxs]
      simp only [_root_.zero_apply, mul_zero, abs_zero]
      exact mul_nonneg (by positivity) (evp_Om_nonneg _ _ (aJet_nonneg φ x))

theorem exists_const_group3_le {N k : ℕ} (hN : k + 7 ≤ N) {K : ℝ} (Λ : ℝ) {ε : ℝ}
    (hε : 0 < ε) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ ρ : EuclideanSpace ℝ (Fin n) → ℝ, IsCutoff ρ Λ →
      ∀ φ : EuclideanSpace ℝ (Fin n) → ℝ, ContDiff ℝ ∞ φ →
      (∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 3, jp N φ x ∈ admSet n N K) →
      |∫ x, ρ x ^ 2 * (2 * ∑ I : Fin (k + 3) → Fin n,
          Dz (fun q : P n N => jetSym potSym (k + 2) q - topPot (k + 2) q) (jp N φ x) I *
            wSym k (jp N φ x) I)| ≤
        ε * (∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 2, ‖jet (k + 5) φ x‖ ^ 2)
        + ε * energy3 (k + 2) 3 φ + C * (1 + energy3 (k + 1) 3 φ) := by
  classical
  have hA := isOpen_admSet (n := n) (N := N) K
  have hlow := admSet_low (n := n) (N := N) K
  set ε' : ℝ := ε / (((Fintype.card (Fin (k + 3) → Fin n) : ℕ) : ℝ) + 1) with hε'
  have hε'0 : 0 < ε' := by positivity
  have hcard : ε' * ((Fintype.card (Fin (k + 3) → Fin n) : ℕ) : ℝ) ≤ ε := by
    rw [hε', div_mul_eq_mul_div, div_le_iff₀ (by positivity)]
    nlinarith [Nat.cast_nonneg (α := ℝ) (Fintype.card (Fin (k + 3) → Fin n))]
  have hex : ∀ I : Fin (k + 3) → Fin n, ∃ C : ℝ, 0 ≤ C ∧ ∀ ξ ∈ {ξ | LLAdm N K Λ ξ},
      |∫ x, ξ.1 x ^ 2 * (2 * (Dz (ll_Rm n N k) (jp N ξ.2 x) I * wSym k (jp N ξ.2 x) I))| ≤
        ε' * energy3 (k + 2) 3 ξ.2 + C * 1 := by
    intro I
    obtain ⟨C, hC0, hC⟩ := ll_g3_term (N := N) (k := k) (K := K) hN Λ I hε'0
    exact ⟨C, hC0, fun ξ hξ => by rw [mul_one]; exact hC ξ.1 hξ.1 ξ.2 hξ.2.1 hξ.2.2⟩
  obtain ⟨C, hC0, hC⟩ := ll_sum_bound (S := {ξ | LLAdm N K Λ ξ})
    (fun I ξ => ∫ x, ξ.1 x ^ 2 * (2 * (Dz (ll_Rm n N k) (jp N ξ.2 x) I * wSym k (jp N ξ.2 x) I)))
    (fun ξ => energy3 (k + 2) 3 ξ.2) (fun _ => 1) (fun ξ _ => energy3_nonneg _ _ _) hcard hex
  refine ⟨C, hC0, fun ρ hρ φ hφ hadm => ?_⟩
  have hsupp : ∀ x ∈ tsupport ρ, jp N φ x ∈ admSet n N K :=
    fun x hx => hadm x (ibb_ball_sub x (hρ.supp hx))
  have hint : ∀ I ∈ (Finset.univ : Finset (Fin (k + 3) → Fin n)),
      Integrable (fun x => ρ x ^ 2 *
        (2 * (Dz (ll_Rm n N k) (jp N φ x) I * wSym k (jp N φ x) I))) := by
    intro I _
    exact (contDiff_locWeight (ll_isSymbol_g3_term hA hlow k (by omega) I) hφ hρ.smooth
      hsupp).continuous.integrable_of_hasCompactSupport
      (hasCompactSupport_locWeight (ibb_hasCompactSupport hρ) _ φ)
  have e : ∫ x, ρ x ^ 2 * (2 * ∑ I : Fin (k + 3) → Fin n,
          Dz (fun q : P n N => jetSym potSym (k + 2) q - topPot (k + 2) q) (jp N φ x) I *
            wSym k (jp N φ x) I) =
      ∑ I : Fin (k + 3) → Fin n, ∫ x, ρ x ^ 2 *
        (2 * (Dz (ll_Rm n N k) (jp N φ x) I * wSym k (jp N φ x) I)) := by
    rw [← integral_finsetSum _ hint]
    congr 1
    funext x
    simp only [Finset.mul_sum]
  rw [e]
  exact ll_fin_odd hε (ll_J_nonneg _ φ) (energy3_nonneg _ _ φ) hC0 (hC (ρ, φ) ⟨hρ, hφ, hadm⟩)

/-! ### Part 6: the four groups together -/

theorem ll_isSymbol_pair {A : Set (P n N)} (hA : IsOpen A)
    (hlow : ∀ m : Fin (N + 1), (m : ℕ) ≤ 3 → ∃ K : ℝ, ∀ p ∈ A, ‖p.2 m‖ ≤ K)
    (k : ℕ) (hN : k + 3 ≤ N) {w ℓ : ℕ} {T : P n N → Z n (k + 3)} (hT : IsSymbol A w ℓ T) :
    IsSymbol A (w + k) (max ℓ k)
      (fun q : P n N => 2 * ∑ I : Fin (k + 3) → Fin n, T q I * wSym k q I) := by
  have hs := IsSymbol.sum hA (Finset.univ : Finset (Fin (k + 3) → Fin n))
    (fun I _ => (hT.clm_comp (PiLp.proj (𝕜 := ℝ) 2 (fun _ : Fin (k + 3) → Fin n => ℝ) I)).mul
      (isSymbol_wSym_apply hA hlow k hN I))
  exact ((isSymbol_const hA (2 : ℝ)).mul hs).mono (by omega) (by omega)

theorem ll_isSymbol_cTop {A : Set (P n N)} (hA : IsOpen A)
    (hlow : ∀ m : Fin (N + 1), (m : ℕ) ≤ 3 → ∃ K : ℝ, ∀ p ∈ A, ‖p.2 m‖ ≤ K)
    (k : ℕ) (hk : 1 ≤ k) (hN : k + 3 ≤ N) {w ℓ : ℕ} {𝔷 : P n N → Z n (k + 2)}
    (h : IsSymbol A w ℓ 𝔷) :
    IsSymbol A w (max (k - 1) ℓ) (fun q : P n N => cTop k q (𝔷 q)) := by
  have h2 : k + 2 < N + 1 := by omega
  have hb := IsSymbol.bilin
    (ContinuousLinearMap.apply ℝ (Z n (k + 3)) :
      Z n (k + 2) →L[ℝ] (Z n (k + 2) →L[ℝ] Z n (k + 3)) →L[ℝ] Z n (k + 3)).flip
    (ll_isSymbol_cL hA hlow k hk hN h2) h
  have hb' := hb.mono (w' := w) (ℓ' := max (k - 1) ℓ) (by omega) le_rfl
  exact hb'.congr fun q _ => ll_cTop_eq k h2 q (𝔷 q)

theorem exists_const_lower_le {N k : ℕ} (hk : 1 ≤ k) (hN : k + 7 ≤ N) {K : ℝ} (Λ : ℝ) {ε : ℝ}
    (hε : 0 < ε) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ ρ : EuclideanSpace ℝ (Fin n) → ℝ, IsCutoff ρ Λ →
      ∀ φ : EuclideanSpace ℝ (Fin n) → ℝ, ContDiff ℝ ∞ φ →
      (∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 3, jp N φ x ∈ admSet n N K) →
      |∫ x, ρ x ^ 2 * (2 * ∑ I : Fin (k + 3) → Fin n,
          (rPot (k + 2) (jp N φ x) I
            + Dz (fun q : P n N => jetSym potSym (k + 2) q - topPot (k + 2) q) (jp N φ x) I
            + cTop k (jp N φ x) (topPot (k + 2) (jp N φ x)) I
            + cTop k (jp N φ x) (jetSym potSym (k + 2) (jp N φ x) - topPot (k + 2) (jp N φ x)) I)
          * wSym k (jp N φ x) I)| ≤
        ε * (∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 2, ‖jet (k + 5) φ x‖ ^ 2)
        + ε * energy3 (k + 2) 3 φ + C * (1 + energy3 (k + 1) 3 φ) := by
  have hA := isOpen_admSet (n := n) (N := N) K
  have hlow := admSet_low (n := n) (N := N) K
  have hε4 : 0 < ε / 4 := by positivity
  obtain ⟨C2, h20, h2⟩ := exists_const_group2_le (n := n) (N := N) (k := k) (K := K)
    (by omega) Λ hε4
  obtain ⟨C3, h30, h3⟩ := exists_const_group3_le (n := n) (N := N) (k := k) (K := K)
    hN Λ hε4
  obtain ⟨C4, h40, h4⟩ := exists_const_group4_le (n := n) (N := N) (k := k) (K := K)
    hk hN Λ hε4
  obtain ⟨C5, h50, h5⟩ := exists_const_group5_le (n := n) (N := N) (k := k) (K := K)
    hk hN Λ hε4
  refine ⟨C2 + C3 + C4 + C5, by positivity, fun ρ hρ φ hφ hadm => ?_⟩
  have hsupp : ∀ x ∈ tsupport ρ, jp N φ x ∈ admSet n N K :=
    fun x hx => hadm x (ibb_ball_sub x (hρ.supp hx))
  have hRm := ll_isSymbol_Rm hA hlow k (by omega : k + 6 ≤ N)
  have s2 := ll_isSymbol_pair hA hlow k (by omega) (isSymbol_rPot hA hlow (k + 2))
  have s3 := ll_isSymbol_pair hA hlow k (by omega) (hRm.Dz (by omega) hlow)
  have s4 := ll_isSymbol_pair hA hlow k (by omega)
    (ll_isSymbol_cTop hA hlow k hk (by omega) (isSymbol_topPot hA hlow (k + 2)))
  have s5 := ll_isSymbol_pair hA hlow k (by omega)
    (ll_isSymbol_cTop hA hlow k hk (by omega) hRm)
  have hI : ∀ {w ℓ : ℕ} {𝔣 : P n N → ℝ}, IsSymbol (admSet n N K) w ℓ 𝔣 →
      Integrable (fun x => ρ x ^ 2 * 𝔣 (jp N φ x)) := fun hs =>
    (contDiff_locWeight hs hφ hρ.smooth hsupp).continuous.integrable_of_hasCompactSupport
      (hasCompactSupport_locWeight (ibb_hasCompactSupport hρ) _ φ)
  have i2 := hI s2
  have i3 : Integrable (fun x => ρ x ^ 2 * (2 * ∑ I : Fin (k + 3) → Fin n,
          Dz (fun q : P n N => jetSym potSym (k + 2) q - topPot (k + 2) q) (jp N φ x) I *
            wSym k (jp N φ x) I)) := hI s3
  have i4 := hI s4
  have i5 : Integrable (fun x => ρ x ^ 2 * (2 * ∑ I : Fin (k + 3) → Fin n,
          cTop k (jp N φ x) (jetSym potSym (k + 2) (jp N φ x) - topPot (k + 2) (jp N φ x)) I *
            wSym k (jp N φ x) I)) := hI s5
  have e : (fun x => ρ x ^ 2 * (2 * ∑ I : Fin (k + 3) → Fin n,
          (rPot (k + 2) (jp N φ x) I
            + Dz (fun q : P n N => jetSym potSym (k + 2) q - topPot (k + 2) q) (jp N φ x) I
            + cTop k (jp N φ x) (topPot (k + 2) (jp N φ x)) I
            + cTop k (jp N φ x) (jetSym potSym (k + 2) (jp N φ x) - topPot (k + 2) (jp N φ x)) I)
          * wSym k (jp N φ x) I)) =
      fun x => ρ x ^ 2 * (2 * ∑ I : Fin (k + 3) → Fin n,
          rPot (k + 2) (jp N φ x) I * wSym k (jp N φ x) I)
        + ρ x ^ 2 * (2 * ∑ I : Fin (k + 3) → Fin n,
          Dz (fun q : P n N => jetSym potSym (k + 2) q - topPot (k + 2) q) (jp N φ x) I *
            wSym k (jp N φ x) I)
        + ρ x ^ 2 * (2 * ∑ I : Fin (k + 3) → Fin n,
          cTop k (jp N φ x) (topPot (k + 2) (jp N φ x)) I * wSym k (jp N φ x) I)
        + ρ x ^ 2 * (2 * ∑ I : Fin (k + 3) → Fin n,
          cTop k (jp N φ x) (jetSym potSym (k + 2) (jp N φ x) - topPot (k + 2) (jp N φ x)) I *
            wSym k (jp N φ x) I) := by
    funext x
    simp only [add_mul, Finset.sum_add_distrib, mul_add]
  rw [e, integral_add ((i2.fun_add i3).fun_add i4) i5, integral_add (i2.fun_add i3) i4, integral_add i2 i3]
  have b2 := h2 ρ hρ φ hφ hadm
  have b3 := h3 ρ hρ φ hφ hadm
  have b4 := h4 ρ hρ φ hφ hadm
  have b5 := h5 ρ hρ φ hφ hadm
  refine (ll_abs4 _ _ _ _).trans ?_
  nlinarith [b2, b3, b4, b5]

/-- info: 'CWFlow.Extension.exists_const_group2_le' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_const_group2_le

/-- info: 'CWFlow.Extension.exists_const_group3_le' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_const_group3_le

/-- info: 'CWFlow.Extension.exists_const_group4_le' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_const_group4_le

/-- info: 'CWFlow.Extension.exists_const_group5_le' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_const_group5_le

/-- info: 'CWFlow.Extension.exists_const_lower_le' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_const_lower_le

end CWFlow.Extension
