import CWFlow.Extension.Symbols
import Mathlib.Analysis.Calculus.FDeriv.Bilinear
import Mathlib.Analysis.Calculus.ContDiff.Basic
import Mathlib.Analysis.Calculus.ContDiff.Operations

/-!
# Symbol classes are closed under continuous bilinear products

If `𝔭` is a symbol of weight `w` and `𝔮` a symbol of weight `w'`, then `B (𝔭, 𝔮)` is a symbol of
weight `w + w'` (and factor weight `max ℓ ℓ'`).  The proof is an induction on the order of
differentiation, using the slot projections `iotaP s` to split the extra argument of the derivative
according to its slot (so that each piece carries a single integer weight).
-/

open scoped ContDiff

namespace CWFlow.Extension

variable {n N : ℕ}

/-! ### Slot projections -/

/-- The projection of a variation onto its `s`-component. -/
noncomputable def iotaP : Slot N → (P n N →L[ℝ] P n N)
  | none => (ContinuousLinearMap.fst ℝ _ _).prod 0
  | some m => (0 : P n N →L[ℝ] EuclideanSpace ℝ (Fin n)).prod
      ((ContinuousLinearMap.single ℝ (fun m : Fin (N + 1) => Z n m) m).comp (slotProjL m))

theorem iotaP_none_apply (δ : P n N) : iotaP (none : Slot N) δ = (δ.1, 0) := rfl

theorem iotaP_some_apply (m : Fin (N + 1)) (δ : P n N) :
    iotaP (some m : Slot N) δ = (0, Pi.single m (δ.2 m)) := rfl

theorem iotaP_sum (δ : P n N) : ∑ s : Slot N, iotaP s δ = δ := by
  rw [Fintype.sum_option]
  rw [iotaP_none_apply]
  simp only [iotaP_some_apply]
  refine Prod.ext ?_ ?_
  · simp [Prod.fst_sum]
  · simp only [Prod.snd_add, Prod.snd_sum]
    simp [Finset.univ_sum_single]

theorem iotaP_idem (s : Slot N) (δ : P n N) : iotaP s (iotaP s δ) = iotaP s δ := by
  cases s with
  | none => simp [iotaP_none_apply]
  | some m => simp [iotaP_some_apply]

theorem comp_iotaP_ne {s t : Slot N} (h : t ≠ s) (δ : P n N) : comp t (iotaP s δ) = 0 := by
  cases s with
  | none =>
    cases t with
    | none => exact absurd rfl h
    | some m' => simp [iotaP_none_apply, comp]
  | some m =>
    cases t with
    | none => simp [iotaP_some_apply, comp]
    | some m' =>
      have : m' ≠ m := fun e => h (by rw [e])
      simp [iotaP_some_apply, comp, Pi.single_eq_of_ne this]

theorem comp_iotaP_self (s : Slot N) (δ : P n N) : comp s (iotaP s δ) = comp s δ := by
  cases s with
  | none => simp [iotaP_none_apply, comp]
  | some m => simp [iotaP_some_apply, comp]

theorem comp_le_norm (s : Slot N) (δ : P n N) : comp s δ ≤ ‖δ‖ := by
  cases s with
  | none => exact norm_fst_le δ
  | some m => exact (norm_le_pi_norm (δ.2) m).trans (norm_snd_le δ)

theorem norm_iotaP_le (s : Slot N) (δ : P n N) : ‖iotaP s δ‖ ≤ comp s δ := by
  cases s with
  | none => simp [iotaP_none_apply, comp, Prod.norm_def]
  | some m =>
    rw [iotaP_some_apply, Prod.norm_def, norm_zero, max_le_iff]
    refine ⟨comp_nonneg _ _, ?_⟩
    refine (pi_norm_le_iff_of_nonneg (comp_nonneg (some m : Slot N) δ)).2 fun m' => ?_
    by_cases h : m' = m
    · subst h
      simp [comp]
    · show ‖(Pi.single m (δ.2 m) : J n N) m'‖ ≤ _
      rw [Pi.single_eq_of_ne h]
      simpa using comp_nonneg (some m : Slot N) δ


/-! ### Truncated bounds with integer weights -/

/-- The weighted sum appearing in the bound. -/
noncomputable def bdS (v : ℤ) (ℓ r : ℕ) (p : P n N) (δ : Fin r → P n N) : ℝ :=
  ∑ s : Fin r → Slot N, OmZ (v - ∑ i, (wt (s i) : ℤ)) ℓ (aOf p) * ∏ i, comp (s i) (δ i)

theorem bdS_nonneg (v : ℤ) (ℓ r : ℕ) (p : P n N) (δ : Fin r → P n N) : 0 ≤ bdS v ℓ r p δ :=
  Finset.sum_nonneg fun s _ =>
    mul_nonneg (OmZ_nonneg (aOf_nonneg p)) (Finset.prod_nonneg fun i _ => comp_nonneg _ _)

theorem bdS_zero (v : ℤ) (ℓ : ℕ) (p : P n N) (δ : Fin 0 → P n N) :
    bdS v ℓ 0 p δ = OmZ v ℓ (aOf p) := by
  simp [bdS]

theorem bdS_succ (v : ℤ) (ℓ r : ℕ) (p : P n N) (δ' : Fin (r + 1) → P n N) :
    bdS v ℓ (r + 1) p δ' = ∑ t : Slot N,
      bdS (v - wt t) ℓ r p (Fin.init δ') * comp t (δ' (Fin.last r)) := by
  unfold bdS
  rw [← Equiv.sum_comp (Fin.snocEquiv (fun _ : Fin (r + 1) => Slot N))]
  rw [Fintype.sum_prod_type]
  refine Finset.sum_congr rfl fun t _ => ?_
  rw [Finset.sum_mul]
  refine Finset.sum_congr rfl fun σ _ => ?_
  simp only [Fin.snocEquiv_apply, Fin.sum_univ_castSucc, Fin.prod_univ_castSucc, Fin.snoc_castSucc,
    Fin.snoc_last, Fin.init]
  have : v - (∑ i : Fin r, (wt (σ i) : ℤ) + (wt t : ℤ)) = v - wt t - ∑ i : Fin r, (wt (σ i) : ℤ) := by
    ring
  rw [this]
  ring

section BU

variable {W : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W] {A : Set (P n N)}

/-- The bounds of all orders `≤ R`, with an integer weight. -/
def BU (R : ℕ) (A : Set (P n N)) (v : ℤ) (ℓ : ℕ) (𝔭 : P n N → W) : Prop :=
  ∀ r ≤ R, ∃ C : ℝ, ∀ p ∈ A, ∀ δ : Fin r → P n N,
    ‖iteratedFDeriv ℝ r 𝔭 p δ‖ ≤ C * bdS v ℓ r p δ

theorem BU.mono {R R' : ℕ} {v : ℤ} {ℓ : ℕ} {𝔭 : P n N → W} (h : BU R' A v ℓ 𝔭) (hR : R ≤ R') :
    BU R A v ℓ 𝔭 := fun r hr => h r (hr.trans hR)

theorem BU.congr {R : ℕ} {v : ℤ} {ℓ : ℕ} {f g : P n N → W} (hA : IsOpen A)
    (h : ∀ q ∈ A, f q = g q) (hb : BU R A v ℓ f) : BU R A v ℓ g := by
  intro r hr
  obtain ⟨C, hC⟩ := hb r hr
  refine ⟨C, fun p hp δ => ?_⟩
  have hev : f =ᶠ[nhds p] g := Filter.eventuallyEq_of_mem (hA.mem_nhds hp) h
  have : iteratedFDeriv ℝ r f p = iteratedFDeriv ℝ r g p :=
    (Filter.EventuallyEq.iteratedFDeriv ℝ hev r).eq_of_nhds
  rw [← this]
  exact hC p hp δ

theorem BU.add {R : ℕ} {v : ℤ} {ℓ : ℕ} {f g : P n N → W} (hA : IsOpen A)
    (hf : ContDiffOn ℝ ∞ f A) (hg : ContDiffOn ℝ ∞ g A) (bf : BU R A v ℓ f) (bg : BU R A v ℓ g) :
    BU R A v ℓ (fun p => f p + g p) := by
  intro r hr
  obtain ⟨C, hC⟩ := bf r hr
  obtain ⟨C', hC'⟩ := bg r hr
  refine ⟨C + C', fun p hp δ => ?_⟩
  have e : iteratedFDeriv ℝ r (fun p => f p + g p) p =
      iteratedFDeriv ℝ r f p + iteratedFDeriv ℝ r g p :=
    iteratedFDeriv_add_apply
      ((hf.contDiffAt (hA.mem_nhds hp)).of_le (by exact_mod_cast le_top))
      ((hg.contDiffAt (hA.mem_nhds hp)).of_le (by exact_mod_cast le_top))
  rw [e, add_apply, add_mul]
  exact (norm_add_le _ _).trans (add_le_add (hC p hp δ) (hC' p hp δ))

/-- The derivative of `𝔭`, restricted to the slot `s` of the new argument. -/
noncomputable def Dsl (s : Slot N) (𝔭 : P n N → W) : P n N → (P n N →L[ℝ] W) :=
  fun q => (fderiv ℝ 𝔭 q).comp (iotaP s)

/-- Precomposition with the slot projection, as a continuous linear map. -/
noncomputable def Phi (s : Slot N) : (P n N →L[ℝ] W) →L[ℝ] (P n N →L[ℝ] W) :=
  (ContinuousLinearMap.compL ℝ (P n N) (P n N) W).flip (iotaP s)

theorem Phi_apply (s : Slot N) (L : P n N →L[ℝ] W) : Phi s L = L.comp (iotaP s) := rfl

theorem contDiffOn_Dsl (hA : IsOpen A) {𝔭 : P n N → W} (h𝔭 : ContDiffOn ℝ ∞ 𝔭 A) (s : Slot N) :
    ContDiffOn ℝ ∞ (Dsl s 𝔭) A := by
  have h1 : ContDiffOn ℝ ∞ (fderiv ℝ 𝔭) A := h𝔭.fderiv_of_isOpen hA (by simp)
  exact (Phi s).contDiff.comp_contDiffOn h1

theorem iFD_Dsl (hA : IsOpen A) {𝔭 : P n N → W} (h𝔭 : ContDiffOn ℝ ∞ 𝔭 A) (s : Slot N)
    {p : P n N} (hp : p ∈ A) (r : ℕ) (δ : Fin r → P n N) :
    iteratedFDeriv ℝ r (Dsl s 𝔭) p δ =
      (iteratedFDeriv ℝ r (fun y => fderiv ℝ 𝔭 y) p δ).comp (iotaP s) := by
  have h1 : ContDiffOn ℝ ∞ (fderiv ℝ 𝔭) A := h𝔭.fderiv_of_isOpen hA (by simp)
  have hd : ContDiffAt ℝ r (fderiv ℝ 𝔭) p :=
    (h1.contDiffAt (hA.mem_nhds hp)).of_le (by exact_mod_cast le_top)
  have e : Dsl s 𝔭 = (Phi s) ∘ (fderiv ℝ 𝔭) := rfl
  rw [e, (Phi s).iteratedFDeriv_comp_left hd le_rfl]
  rfl


theorem bu_deriv (hA : IsOpen A) {𝔭 : P n N → W} (h𝔭 : ContDiffOn ℝ ∞ 𝔭 A) {R : ℕ} {v : ℤ}
    {ℓ : ℕ} (hb : BU (R + 1) A v ℓ 𝔭) (s : Slot N) : BU R A (v - wt s) ℓ (Dsl s 𝔭) := by
  intro r hr
  obtain ⟨C, hC⟩ := hb (r + 1) (by omega)
  refine ⟨|C|, fun p hp δ => ?_⟩
  rw [iFD_Dsl hA h𝔭 s hp]
  refine ContinuousLinearMap.opNorm_le_bound _
    (mul_nonneg (abs_nonneg C) (bdS_nonneg _ _ _ _ _)) fun δ₀ => ?_
  rw [ContinuousLinearMap.comp_apply]
  have h1 := hC p hp (Fin.snoc δ (iotaP s δ₀))
  rw [iteratedFDeriv_succ_apply_right] at h1
  simp only [Fin.init_snoc, Fin.snoc_last] at h1
  rw [bdS_succ, Finset.sum_eq_single s] at h1
  · simp only [Fin.init_snoc, Fin.snoc_last, comp_iotaP_self] at h1
    have h2 : C * (bdS (v - wt s) ℓ r p δ * comp s δ₀) ≤
        |C| * bdS (v - wt s) ℓ r p δ * ‖δ₀‖ := by
      have hS := bdS_nonneg (v - wt s) ℓ r p δ
      calc C * (bdS (v - wt s) ℓ r p δ * comp s δ₀)
          ≤ |C| * (bdS (v - wt s) ℓ r p δ * comp s δ₀) :=
            mul_le_mul_of_nonneg_right (le_abs_self C) (mul_nonneg hS (comp_nonneg _ _))
        _ ≤ |C| * (bdS (v - wt s) ℓ r p δ * ‖δ₀‖) :=
            mul_le_mul_of_nonneg_left (mul_le_mul_of_nonneg_left (comp_le_norm s δ₀) hS)
              (abs_nonneg C)
        _ = _ := by ring
    exact h1.trans h2
  · intro t _ ht
    simp only [Fin.snoc_last]
    rw [comp_iotaP_ne ht, mul_zero]
  · intro h; exact absurd (Finset.mem_univ s) h

theorem bu_succ (hA : IsOpen A) {F : P n N → W} (hF : ContDiffOn ℝ ∞ F A) {R : ℕ} {v : ℤ}
    {ℓ : ℕ} (h0 : BU 0 A v ℓ F) (hs : ∀ s : Slot N, BU R A (v - wt s) ℓ (Dsl s F)) :
    BU (R + 1) A v ℓ F := by
  intro r hr
  rcases r with _ | r
  · exact h0 0 le_rfl
  have hr' : r ≤ R := by omega
  choose C hC using fun s => hs s r hr'
  set K : ℝ := ∑ s : Slot N, |C s| with hK
  refine ⟨K, fun p hp δ' => ?_⟩
  rw [iteratedFDeriv_succ_apply_right, bdS_succ]
  set M := iteratedFDeriv ℝ r (fun y => fderiv ℝ F y) p (Fin.init δ') with hM
  set δ₀ := δ' (Fin.last r) with hδ₀
  have hsum : M δ₀ = ∑ s : Slot N, M (iotaP s δ₀) := by
    conv_lhs => rw [← iotaP_sum δ₀]
    exact map_sum M _ _
  have hterm : ∀ s : Slot N, ‖M (iotaP s δ₀)‖ ≤
      |C s| * (bdS (v - wt s) ℓ r p (Fin.init δ') * comp s δ₀) := by
    intro s
    have hN := iFD_Dsl hA hF s hp r (Fin.init δ')
    rw [← hM] at hN
    have hnn := bdS_nonneg (v - wt s) ℓ r p (Fin.init δ')
    have e1 : M (iotaP s δ₀) = (iteratedFDeriv ℝ r (Dsl s F) p (Fin.init δ')) (iotaP s δ₀) := by
      rw [hN, ContinuousLinearMap.comp_apply, iotaP_idem]
    rw [e1]
    calc ‖(iteratedFDeriv ℝ r (Dsl s F) p (Fin.init δ')) (iotaP s δ₀)‖
        ≤ ‖iteratedFDeriv ℝ r (Dsl s F) p (Fin.init δ')‖ * ‖iotaP s δ₀‖ :=
          ContinuousLinearMap.le_opNorm _ _
      _ ≤ (|C s| * bdS (v - wt s) ℓ r p (Fin.init δ')) * comp s δ₀ := by
          refine mul_le_mul ?_ (norm_iotaP_le s δ₀) (norm_nonneg _)
            (mul_nonneg (abs_nonneg _) hnn)
          exact (hC s p hp (Fin.init δ')).trans
            (mul_le_mul_of_nonneg_right (le_abs_self _) hnn)
      _ = _ := by ring
  rw [hsum]
  refine (norm_sum_le _ _).trans ?_
  rw [Finset.mul_sum]
  refine Finset.sum_le_sum fun s _ => (hterm s).trans ?_
  have hle : |C s| ≤ K :=
    Finset.single_le_sum (f := fun s => |C s|) (fun s _ => abs_nonneg _) (Finset.mem_univ s)
  exact mul_le_mul_of_nonneg_right hle
    (mul_nonneg (bdS_nonneg _ _ _ _ _) (comp_nonneg _ _))

end BU

section Bilin

variable {W₁ W₂ W₃ : Type*} [NormedAddCommGroup W₁] [NormedSpace ℝ W₁]
  [NormedAddCommGroup W₂] [NormedSpace ℝ W₂] [NormedAddCommGroup W₃] [NormedSpace ℝ W₃]
  {A : Set (P n N)}

theorem contDiffOn_bilin (B : W₁ →L[ℝ] W₂ →L[ℝ] W₃) {𝔭 : P n N → W₁} {𝔮 : P n N → W₂}
    (h𝔭 : ContDiffOn ℝ ∞ 𝔭 A) (h𝔮 : ContDiffOn ℝ ∞ 𝔮 A) :
    ContDiffOn ℝ ∞ (fun p => B (𝔭 p) (𝔮 p)) A :=
  (B.contDiff.comp_contDiffOn h𝔭).clm_apply h𝔮

theorem bu_zero (B : W₁ →L[ℝ] W₂ →L[ℝ] W₃) {v v' : ℤ} {ℓ ℓ' : ℕ} {𝔭 : P n N → W₁}
    {𝔮 : P n N → W₂} (b₁ : BU 0 A v ℓ 𝔭) (b₂ : BU 0 A v' ℓ' 𝔮) :
    BU 0 A (v + v') (max ℓ ℓ') (fun p => B (𝔭 p) (𝔮 p)) := by
  intro r hr
  obtain rfl : r = 0 := by omega
  obtain ⟨C₁, h₁⟩ := b₁ 0 le_rfl
  obtain ⟨C₂, h₂⟩ := b₂ 0 le_rfl
  obtain ⟨K, hK0, hK⟩ := exists_OmZ_mul v v'
  refine ⟨‖B‖ * (|C₁| * |C₂|) * K, fun p hp δ => ?_⟩
  have e1 := h₁ p hp δ
  have e2 := h₂ p hp δ
  simp only [iteratedFDeriv_zero_apply, bdS_zero] at e1 e2 ⊢
  have hX := OmZ_nonneg (v := v) (ℓ := ℓ) (aOf_nonneg p)
  have hY := OmZ_nonneg (v := v') (ℓ := ℓ') (aOf_nonneg p)
  have e1' : ‖𝔭 p‖ ≤ |C₁| * OmZ v ℓ (aOf p) :=
    e1.trans (mul_le_mul_of_nonneg_right (le_abs_self _) hX)
  have e2' : ‖𝔮 p‖ ≤ |C₂| * OmZ v' ℓ' (aOf p) :=
    e2.trans (mul_le_mul_of_nonneg_right (le_abs_self _) hY)
  have h3 : ‖B (𝔭 p) (𝔮 p)‖ ≤ ‖B‖ * ‖𝔭 p‖ * ‖𝔮 p‖ := B.le_opNorm₂ _ _
  have h4 : ‖𝔭 p‖ * ‖𝔮 p‖ ≤ (|C₁| * OmZ v ℓ (aOf p)) * (|C₂| * OmZ v' ℓ' (aOf p)) :=
    mul_le_mul e1' e2' (norm_nonneg _) (mul_nonneg (abs_nonneg _) hX)
  have h5 := hK ℓ ℓ' (aOf p) (aOf_nonneg p)
  calc ‖B (𝔭 p) (𝔮 p)‖ ≤ ‖B‖ * (‖𝔭 p‖ * ‖𝔮 p‖) := by rw [← mul_assoc]; exact h3
    _ ≤ ‖B‖ * ((|C₁| * OmZ v ℓ (aOf p)) * (|C₂| * OmZ v' ℓ' (aOf p))) :=
        mul_le_mul_of_nonneg_left h4 (norm_nonneg B)
    _ = ‖B‖ * (|C₁| * |C₂|) * (OmZ v ℓ (aOf p) * OmZ v' ℓ' (aOf p)) := by ring
    _ ≤ ‖B‖ * (|C₁| * |C₂|) * (K * OmZ (v + v') (max ℓ ℓ') (aOf p)) :=
        mul_le_mul_of_nonneg_left h5 (by positivity)
    _ = _ := by ring

end Bilin

universe u1 u2 u3

theorem bu_bilin (R : ℕ) : ∀ {W₁ : Type u1} {W₂ : Type u2} {W₃ : Type u3}
    [NormedAddCommGroup W₁] [NormedSpace ℝ W₁] [NormedAddCommGroup W₂] [NormedSpace ℝ W₂]
    [NormedAddCommGroup W₃] [NormedSpace ℝ W₃] {A : Set (P n N)} (_ : IsOpen A) {v v' : ℤ}
    {ℓ ℓ' : ℕ} {𝔭 : P n N → W₁} {𝔮 : P n N → W₂} (B : W₁ →L[ℝ] W₂ →L[ℝ] W₃),
    ContDiffOn ℝ ∞ 𝔭 A → ContDiffOn ℝ ∞ 𝔮 A → BU R A v ℓ 𝔭 → BU R A v' ℓ' 𝔮 →
    BU R A (v + v') (max ℓ ℓ') (fun p => B (𝔭 p) (𝔮 p)) := by
  induction R with
  | zero =>
    intro W₁ W₂ W₃ _ _ _ _ _ _ A hA v v' ℓ ℓ' 𝔭 𝔮 B h𝔭 h𝔮 b₁ b₂
    exact bu_zero B b₁ b₂
  | succ R ih =>
    intro W₁ W₂ W₃ _ _ _ _ _ _ A hA v v' ℓ ℓ' 𝔭 𝔮 B h𝔭 h𝔮 b₁ b₂
    refine bu_succ hA (contDiffOn_bilin B h𝔭 h𝔮)
      (bu_zero B (b₁.mono (Nat.zero_le _)) (b₂.mono (Nat.zero_le _))) fun s => ?_
    have e1 := ih (W₁ := P n N →L[ℝ] W₁) (W₂ := W₂) (W₃ := P n N →L[ℝ] W₃) hA
      (B.precompL (P n N)) (contDiffOn_Dsl hA h𝔭 s) h𝔮 (bu_deriv hA h𝔭 b₁ s)
      (b₂.mono (Nat.le_succ R))
    have e2 := ih (W₁ := W₁) (W₂ := P n N →L[ℝ] W₂) (W₃ := P n N →L[ℝ] W₃) hA
      (B.precompR (P n N)) h𝔭 (contDiffOn_Dsl hA h𝔮 s) (b₁.mono (Nat.le_succ R))
      (bu_deriv hA h𝔮 b₂ s)
    have w1 : v - (wt s : ℤ) + v' = v + v' - wt s := by ring
    have w2 : v + (v' - (wt s : ℤ)) = v + v' - wt s := by ring
    rw [w1] at e1
    rw [w2] at e2
    have hs1 := contDiffOn_bilin (B.precompL (P n N)) (contDiffOn_Dsl hA h𝔭 s) h𝔮
    have hs2 := contDiffOn_bilin (B.precompR (P n N)) h𝔭 (contDiffOn_Dsl hA h𝔮 s)
    have e3 := BU.add hA hs2 hs1 e2 e1
    refine BU.congr hA (fun q hq => ?_) e3
    have hd𝔭 : DifferentiableAt ℝ 𝔭 q :=
      ((h𝔭.contDiffAt (hA.mem_nhds hq)).of_le (by exact_mod_cast le_top) :
        ContDiffAt ℝ 1 𝔭 q).differentiableAt one_ne_zero
    have hd𝔮 : DifferentiableAt ℝ 𝔮 q :=
      ((h𝔮.contDiffAt (hA.mem_nhds hq)).of_le (by exact_mod_cast le_top) :
        ContDiffAt ℝ 1 𝔮 q).differentiableAt one_ne_zero
    simp only [Dsl]
    rw [ContinuousLinearMap.fderiv_of_bilinear B hd𝔭 hd𝔮]
    refine ContinuousLinearMap.ext fun δ => ?_
    simp

/-- A symbol gives the truncated bounds of every order. -/
theorem IsSymbol.bu {W : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W] {A : Set (P n N)}
    {w ℓ : ℕ} {𝔭 : P n N → W} (h : IsSymbol A w ℓ 𝔭) (R : ℕ) : BU R A (w : ℤ) ℓ 𝔭 :=
  fun r _ => h.bound r

/-- Symbols are closed under continuous bilinear products. -/
theorem IsSymbol.bilin {W₁ W₂ W₃ : Type*} [NormedAddCommGroup W₁] [NormedSpace ℝ W₁]
    [NormedAddCommGroup W₂] [NormedSpace ℝ W₂] [NormedAddCommGroup W₃] [NormedSpace ℝ W₃]
    {A : Set (P n N)} {w w' ℓ ℓ' : ℕ} {𝔭 : P n N → W₁} {𝔮 : P n N → W₂}
    (B : W₁ →L[ℝ] W₂ →L[ℝ] W₃) (h₁ : IsSymbol A w ℓ 𝔭) (h₂ : IsSymbol A w' ℓ' 𝔮) :
    IsSymbol A (w + w') (max ℓ ℓ') (fun p => B (𝔭 p) (𝔮 p)) := by
  refine ⟨h₁.isOpen, contDiffOn_bilin B h₁.smooth h₂.smooth, fun r => ?_⟩
  have key := bu_bilin r h₁.isOpen B h₁.smooth h₂.smooth (h₁.bu r) (h₂.bu r)
  obtain ⟨C, hC⟩ := key r le_rfl
  refine ⟨C, fun p hp δ => ?_⟩
  have := hC p hp δ
  simpa [bdS, Nat.cast_add] using this

theorem IsSymbol.smul {W : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W]
    {A : Set (P n N)} {w w' ℓ ℓ' : ℕ} {𝔭 : P n N → ℝ} {𝔮 : P n N → W}
    (h₁ : IsSymbol A w ℓ 𝔭) (h₂ : IsSymbol A w' ℓ' 𝔮) :
    IsSymbol A (w + w') (max ℓ ℓ') (fun p => 𝔭 p • 𝔮 p) :=
  IsSymbol.bilin (ContinuousLinearMap.lsmul ℝ ℝ) h₁ h₂

theorem IsSymbol.mul {A : Set (P n N)} {w w' ℓ ℓ' : ℕ} {𝔭 𝔮 : P n N → ℝ}
    (h₁ : IsSymbol A w ℓ 𝔭) (h₂ : IsSymbol A w' ℓ' 𝔮) :
    IsSymbol A (w + w') (max ℓ ℓ') (fun p => 𝔭 p * 𝔮 p) :=
  IsSymbol.bilin (ContinuousLinearMap.mul ℝ ℝ) h₁ h₂

/-- info: 'CWFlow.Extension.IsSymbol.bilin' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms IsSymbol.bilin

/-- info: 'CWFlow.Extension.IsSymbol.smul' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms IsSymbol.smul

/-- info: 'CWFlow.Extension.IsSymbol.mul' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms IsSymbol.mul

end CWFlow.Extension
