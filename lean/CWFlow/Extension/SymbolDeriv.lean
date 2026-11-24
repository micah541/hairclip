import CWFlow.Extension.Symbols
import Mathlib.Analysis.Calculus.FDeriv.Prod
import Mathlib.Analysis.Calculus.FDeriv.Pi
import Mathlib.Topology.Algebra.Module.FiniteDimension

/-!
# The total derivative of a function on jet space

`contract e` contracts an `(m+1)`-jet with a vector in its first index; the derivative of the
`m`-jet of a smooth function is `contract` of its `(m+1)`-jet.  `totalD 𝔭 p` is the derivative of
`𝔭` on jet space in the direction `(e, shift e p.2)`, and for a symbol `𝔭` of weight `w` with
`w + 4 ≤ N` it is the derivative of `y ↦ 𝔭 (jp N φ y)` (the missing top slot does not matter).
-/

open scoped ContDiff

namespace CWFlow.Extension

variable {n N : ℕ}

/-! ### Part 1: contraction -/

/-- The linear map underlying `contract`. -/
noncomputable def contractLin {m : ℕ} (e : EuclideanSpace ℝ (Fin n)) :
    Z n (m + 1) →ₗ[ℝ] Z n m where
  toFun z := WithLp.toLp 2 (fun I => ∑ i, e i * z (Fin.cons (α := fun _ => Fin n) i I))
  map_add' z z' := by
    ext I
    simp [mul_add, Finset.sum_add_distrib]
  map_smul' c z := by
    ext I
    simp [Finset.mul_sum, mul_left_comm]

/-- Contraction of an `(m+1)`-jet with a vector in the first index. -/
noncomputable def contract {m : ℕ} (e : EuclideanSpace ℝ (Fin n)) : Z n (m + 1) →L[ℝ] Z n m :=
  LinearMap.toContinuousLinearMap (contractLin e)

theorem contract_apply {m : ℕ} (e : EuclideanSpace ℝ (Fin n)) (z : Z n (m + 1)) (I : Fin m → Fin n) :
    contract e z I = ∑ i, e i * z (Fin.cons (α := fun _ => Fin n) i I) := rfl

theorem contract_add {m : ℕ} (e e' : EuclideanSpace ℝ (Fin n)) (z : Z n (m + 1)) :
    contract (e + e') z = contract e z + contract e' z := by
  ext I
  simp [contract_apply, add_mul, Finset.sum_add_distrib]

theorem contract_smul {m : ℕ} (c : ℝ) (e : EuclideanSpace ℝ (Fin n)) (z : Z n (m + 1)) :
    contract (c • e) z = c • contract e z := by
  ext I
  simp [contract_apply, Finset.mul_sum, mul_assoc]

theorem fderiv_jet_eq_contract {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ) (m : ℕ)
    (x e : EuclideanSpace ℝ (Fin n)) :
    fderiv ℝ (CWFlow.Interp.jet m φ) x e = contract e (CWFlow.Interp.jet (m + 1) φ x) := by
  ext I
  have he : e = ∑ i, e i • EuclideanSpace.single i (1 : ℝ) := by
    ext l
    simp [Pi.single_apply]
  conv_lhs => rw [he]
  rw [map_sum, contract_apply]
  simp only [map_smul, WithLp.ofLp_sum, WithLp.ofLp_smul, Finset.sum_apply, Pi.smul_apply,
    smul_eq_mul]
  refine Finset.sum_congr rfl fun i _ => ?_
  have := CWFlow.Interp.fderiv_jet_apply hφ m x i I
  rw [show (fderiv ℝ (CWFlow.Interp.jet m φ) x (EuclideanSpace.single i 1)).ofLp I = _ from this]

/-! ### Part 2: shift and total derivative -/

/-- The jet-space part of the derivative of `x ↦ jp N φ x` in direction `e`, with the top slot set
to `0`. -/
noncomputable def shift (e : EuclideanSpace ℝ (Fin n)) (ζ : J n N) : J n N := fun m =>
  if h : (m : ℕ) + 1 < N + 1 then contract e (ζ ⟨(m : ℕ) + 1, h⟩) else 0

theorem shift_add_left (e e' : EuclideanSpace ℝ (Fin n)) (ζ : J n N) :
    shift (e + e') ζ = shift e ζ + shift e' ζ := by
  funext m
  unfold shift
  by_cases h : (m : ℕ) + 1 < N + 1
  · simp only [dif_pos h, Pi.add_apply, contract_add]
  · simp only [dif_neg h, Pi.add_apply, add_zero]

theorem shift_smul_left (c : ℝ) (e : EuclideanSpace ℝ (Fin n)) (ζ : J n N) :
    shift (c • e) ζ = c • shift e ζ := by
  funext m
  unfold shift
  by_cases h : (m : ℕ) + 1 < N + 1
  · simp only [dif_pos h, Pi.smul_apply, contract_smul]
  · simp only [dif_neg h, Pi.smul_apply, smul_zero]

/-- `e ↦ (e, shift e ζ)` as a linear map. -/
noncomputable def shiftLin (ζ : J n N) :
    EuclideanSpace ℝ (Fin n) →ₗ[ℝ] EuclideanSpace ℝ (Fin n) × J n N where
  toFun e := (e, shift e ζ)
  map_add' e e' := by
    ext
    · rfl
    · rw [shift_add_left]; rfl
  map_smul' c e := by
    ext
    · rfl
    · rw [shift_smul_left]; rfl

/-- The total derivative of a function on jet space. -/
noncomputable def totalD {W : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W] (𝔭 : P n N → W)
    (p : P n N) : EuclideanSpace ℝ (Fin n) →L[ℝ] W :=
  (fderiv ℝ 𝔭 p).comp (LinearMap.toContinuousLinearMap (shiftLin p.2))

theorem totalD_apply {W : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W] (𝔭 : P n N → W)
    (p : P n N) (e : EuclideanSpace ℝ (Fin n)) :
    totalD 𝔭 p e = fderiv ℝ 𝔭 p (e, shift e p.2) := rfl

/-! ### Part 3: realization -/

theorem jp_snd (N : ℕ) (φ : EuclideanSpace ℝ (Fin n) → ℝ) (x : EuclideanSpace ℝ (Fin n))
    (m : Fin (N + 1)) : (jp N φ x).2 m = CWFlow.Interp.jet (m : ℕ) φ x := rfl

theorem hasFDerivAt_jp {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ)
    (x : EuclideanSpace ℝ (Fin n)) :
    HasFDerivAt (jp N φ)
      ((ContinuousLinearMap.id ℝ (EuclideanSpace ℝ (Fin n))).prod
        (ContinuousLinearMap.pi fun m : Fin (N + 1) =>
          fderiv ℝ (CWFlow.Interp.jet (m : ℕ) φ) x)) x := by
  refine (hasFDerivAt_id x).prodMk ?_
  refine hasFDerivAt_pi'.2 fun m => ?_
  exact ((CWFlow.Interp.contDiff_jet hφ (m : ℕ)).differentiable (by simp) x).hasFDerivAt

theorem hasFDerivAt_symbol_jp {W : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W]
    {A : Set (P n N)} {w ℓ : ℕ} {𝔭 : P n N → W} (h : IsSymbol A w ℓ 𝔭) (hN : w + 4 ≤ N)
    {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ) {x : EuclideanSpace ℝ (Fin n)}
    (hx : jp N φ x ∈ A) :
    HasFDerivAt (fun y => 𝔭 (jp N φ y)) (totalD 𝔭 (jp N φ x)) x := by
  have hd : DifferentiableAt ℝ 𝔭 (jp N φ x) :=
    ((h.smooth.contDiffAt (h.isOpen.mem_nhds hx)).differentiableAt (by simp))
  have hc := hd.hasFDerivAt.comp x (hasFDerivAt_jp (N := N) hφ x)
  refine hc.congr_fderiv ?_
  ext e
  rw [totalD_apply]
  simp only [ContinuousLinearMap.comp_apply]
  rw [← sub_eq_zero, ← map_sub]
  refine h.fderiv_eq_zero _ hx _ ?_
  intro s hs
  cases s with
  | none => simp [comp, wt]
  | some m =>
    have hm : (m : ℕ) + 1 < N + 1 := by
      have : (m : ℕ) - 3 ≤ w := hs
      omega
    simp only [comp, norm_eq_zero]
    have h1 : (fderiv ℝ (CWFlow.Interp.jet (m : ℕ) φ) x e) =
        contract e (CWFlow.Interp.jet ((m : ℕ) + 1) φ x) :=
      fderiv_jet_eq_contract hφ _ x e
    have hm' : (m : ℕ) < N := by omega
    simp [shift, h1, jp, hm']

/-! ### Part 4: pointwise bound along jets -/

theorem aOf_jp (φ : EuclideanSpace ℝ (Fin n) → ℝ) (x : EuclideanSpace ℝ (Fin n)) (j : ℕ) :
    aOf (jp N φ x) j = if 3 + j ≤ N then ‖CWFlow.Interp.jet (3 + j) φ x‖ else 0 := by
  unfold aOf
  by_cases h : 3 + j ≤ N
  · rw [dif_pos (by omega), if_pos h]; rfl
  · rw [dif_neg (by omega), if_neg h]

theorem norm_symbol_jp_le {W : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W]
    {A : Set (P n N)} {w ℓ : ℕ} {𝔭 : P n N → W} (h : IsSymbol A w ℓ 𝔭) :
    ∃ C : ℝ, ∀ (φ : EuclideanSpace ℝ (Fin n) → ℝ) (x : EuclideanSpace ℝ (Fin n)), jp N φ x ∈ A →
      ‖𝔭 (jp N φ x)‖ ≤ C * Om w ℓ (fun j => if 3 + j ≤ N then ‖CWFlow.Interp.jet (3 + j) φ x‖ else 0) := by
  obtain ⟨C, hC⟩ := h.norm_le
  refine ⟨C, fun φ x hx => ?_⟩
  have := hC _ hx
  rwa [show aOf (jp N φ x) = fun j => if 3 + j ≤ N then ‖CWFlow.Interp.jet (3 + j) φ x‖ else 0 from
    funext (aOf_jp φ x)] at this

/-! ### Part 5: smoothness along jets -/

theorem contDiff_jp {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ) :
    ContDiff ℝ ∞ (jp N φ) :=
  contDiff_id.prodMk (contDiff_pi.2 fun m => CWFlow.Interp.contDiff_jet hφ (m : ℕ))

theorem isOpen_jp_preimage {W : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W]
    {A : Set (P n N)} {w ℓ : ℕ} {𝔭 : P n N → W} (h : IsSymbol A w ℓ 𝔭)
    {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ) :
    IsOpen {y | jp N φ y ∈ A} :=
  h.isOpen.preimage (contDiff_jp hφ).continuous

theorem contDiffOn_symbol_jp {W : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W]
    {A : Set (P n N)} {w ℓ : ℕ} {𝔭 : P n N → W} (h : IsSymbol A w ℓ 𝔭)
    {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ) :
    ContDiffOn ℝ ∞ (fun y => 𝔭 (jp N φ y)) {y | jp N φ y ∈ A} :=
  h.smooth.comp (contDiff_jp hφ).contDiffOn (fun _ hy => hy)

/-- info: 'CWFlow.Extension.fderiv_jet_eq_contract' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms fderiv_jet_eq_contract

/-- info: 'CWFlow.Extension.hasFDerivAt_symbol_jp' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms hasFDerivAt_symbol_jp

/-- info: 'CWFlow.Extension.norm_symbol_jp_le' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms norm_symbol_jp_le

/-- info: 'CWFlow.Extension.contDiffOn_symbol_jp' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms contDiffOn_symbol_jp

end CWFlow.Extension
