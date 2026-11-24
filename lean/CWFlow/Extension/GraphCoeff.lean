import CWFlow.Extension.SymbolTotalD
import Mathlib.LinearAlgebra.Matrix.NonsingularInverse
import Mathlib.LinearAlgebra.Matrix.PosDef
import Mathlib.Analysis.Matrix.PosDef
import Mathlib.Algebra.Order.Star.Real
import Mathlib.Analysis.Matrix.Spectrum
import Mathlib.Analysis.SpecialFunctions.Sqrt
import Mathlib.Analysis.Calculus.ContDiff.WithLp

/-!
# Coefficients of a Lagrangian graph

For the graph `x ↦ (x, ∇φ(x))` in `ℝⁿ × ℝⁿ` we write the induced metric `g = 1 + u uᵀ`
(`u = D²φ`), its inverse, the volume density `√det g` and the Christoffel symbols
`Γ^c_{ab} = g^{cd} ζ₃(a,b,f) u_{df}` as explicit globally smooth functions of `(ζ₂, ζ₃)`, and show
that they are symbols of weight `0`.
-/

open scoped ContDiff Matrix

namespace CWFlow.Extension

variable {n N : ℕ}

/-! ### Part 1: definitions -/

/-- The `m`-th jet slot, or `0` if it does not exist. -/
noncomputable def slotOrZero (N m : ℕ) : P n N →L[ℝ] Z n m :=
  if h : m < N + 1 then slotProjL ⟨m, h⟩ else 0

/-- The pair `(ζ₂, ζ₃)` of a point of jet space (`0` for a slot that does not exist). -/
noncomputable def low23 : P n N →L[ℝ] Z n 2 × Z n 3 :=
  (slotOrZero N 2).prod (slotOrZero N 3)

theorem low23_fst (p : P n N) (h : 2 < N + 1) : (low23 p).1 = p.2 ⟨2, h⟩ := by
  simp [low23, slotOrZero, h]
  rfl

theorem low23_snd (p : P n N) (h : 3 < N + 1) : (low23 p).2 = p.2 ⟨3, h⟩ := by
  simp [low23, slotOrZero, h]
  rfl

theorem low23_high (δ : P n N)
    (hδ : δ.1 = 0 ∧ ∀ m : Fin (N + 1), (m : ℕ) ≤ 3 → δ.2 m = 0) : low23 δ = 0 := by
  have key : ∀ m : ℕ, m ≤ 3 → slotOrZero (n := n) N m δ = 0 := by
    intro m hm
    unfold slotOrZero
    by_cases h : m < N + 1
    · rw [dif_pos h, slotProjL_apply]
      exact hδ.2 ⟨m, h⟩ hm
    · rw [dif_neg h]; rfl
  simp only [low23, ContinuousLinearMap.prod_apply]
  rw [key 2 (by norm_num), key 3 (by norm_num)]
  rfl

/-- The metric as a matrix: `g = 1 + u uᵀ`. -/
noncomputable def gMat (u : Z n 2) : Matrix (Fin n) (Fin n) ℝ :=
  fun a b => (if a = b then 1 else 0) + ∑ c, u ![a, c] * u ![b, c]

/-- The metric as a 2-tensor. -/
noncomputable def gFun (u : Z n 2) : Z n 2 := WithLp.toLp 2 (fun I => gMat u (I 0) (I 1))
/-- The inverse metric as a 2-tensor. -/
noncomputable def ginvFun (u : Z n 2) : Z n 2 :=
  WithLp.toLp 2 (fun I => (gMat u)⁻¹ (I 0) (I 1))
/-- The volume density `√det g`. -/
noncomputable def volFun (u : Z n 2) : ℝ := Real.sqrt (gMat u).det
/-- The Christoffel symbols `Γ^c_{ab}` as a 3-tensor with index order `(a, b, c)`. -/
noncomputable def gamFun (q : Z n 2 × Z n 3) : Z n 3 :=
  WithLp.toLp 2 (fun I => ∑ d, ∑ f, (gMat q.1)⁻¹ (I 2) d * q.2 ![I 0, I 1, f] * q.1 ![d, f])

/-! ### Part 2: algebra of the metric -/

/-- The matrix of `u`. -/
def uMat (u : Z n 2) : Matrix (Fin n) (Fin n) ℝ := fun a c => u ![a, c]

theorem gMat_eq (u : Z n 2) : gMat u = 1 + uMat u * (uMat u)ᵀ := by
  ext a b
  simp [gMat, uMat, Matrix.mul_apply, Matrix.one_apply]

theorem gMat_posDef (u : Z n 2) : (gMat u).PosDef := by
  rw [gMat_eq]
  refine Matrix.PosDef.one.add_posSemidef ?_
  simpa using Matrix.posSemidef_self_mul_conjTranspose (uMat u)


theorem det_gMat_pos (u : Z n 2) : 0 < (gMat u).det := (gMat_posDef u).det_pos

theorem one_le_det_gMat (u : Z n 2) : 1 ≤ (gMat u).det := by
  have hP := gMat_posDef u
  have hH := hP.isHermitian
  rw [hH.det_eq_prod_eigenvalues]
  have hle : ∀ i, (1 : ℝ) ≤ hH.eigenvalues i := by
    intro i
    rw [hH.eigenvalues_eq i]
    set v : Fin n → ℝ := ⇑(hH.eigenvectorBasis i) with hv
    have h1 : v ⬝ᵥ v = 1 := by
      have := (orthonormal_iff_ite.1 hH.eigenvectorBasis.orthonormal) i i
      rw [EuclideanSpace.inner_eq_star_dotProduct, if_pos rfl] at this
      simpa [hv] using this
    have h2 : v ⬝ᵥ (gMat u *ᵥ v) = v ⬝ᵥ v + (v ᵥ* uMat u) ⬝ᵥ (v ᵥ* uMat u) := by
      rw [gMat_eq, Matrix.add_mulVec, dotProduct_add, Matrix.one_mulVec, ← Matrix.mulVec_mulVec,
        Matrix.dotProduct_mulVec, Matrix.mulVec_transpose]
    have h3 : 0 ≤ (v ᵥ* uMat u) ⬝ᵥ (v ᵥ* uMat u) :=
      Finset.sum_nonneg fun i _ => mul_self_nonneg _
    simp only [RCLike.re_to_real, star_trivial]
    rw [h2, h1]
    linarith
  calc (1 : ℝ) = ∏ _i : Fin n, (1 : ℝ) := by simp
    _ ≤ _ := ?_
  simp only [RCLike.ofReal_real_eq_id, id]
  exact Finset.prod_le_prod (fun _ _ => zero_le_one) fun i _ => hle i

theorem gMat_mul_inv (u : Z n 2) : gMat u * (gMat u)⁻¹ = 1 :=
  Matrix.mul_nonsing_inv _ (isUnit_iff_ne_zero.2 (det_gMat_pos u).ne')

theorem gMat_inv_mul (u : Z n 2) : (gMat u)⁻¹ * gMat u = 1 :=
  Matrix.nonsing_inv_mul _ (isUnit_iff_ne_zero.2 (det_gMat_pos u).ne')

theorem quad_g (u : Z n 2) (w : Fin n → ℝ) :
    w ⬝ᵥ (gMat u *ᵥ w) = w ⬝ᵥ w + (w ᵥ* uMat u) ⬝ᵥ (w ᵥ* uMat u) := by
  rw [gMat_eq, Matrix.add_mulVec, dotProduct_add, Matrix.one_mulVec]
  congr 1
  rw [← Matrix.mulVec_mulVec, Matrix.dotProduct_mulVec, Matrix.mulVec_transpose]

theorem quad_inv_eq (u : Z n 2) (v : Fin n → ℝ) :
    ∑ a, ∑ b, (gMat u)⁻¹ a b * v a * v b = v ⬝ᵥ ((gMat u)⁻¹ *ᵥ v) := by
  simp only [dotProduct, Matrix.mulVec, Finset.mul_sum]
  refine Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun b _ => by ring

/-- Key estimate: with `w = g⁻¹ v`, `vᵀ g⁻¹ v = wᵀ w + |wᵀ M|²`. -/
theorem quad_inv_ge (u : Z n 2) (v : Fin n → ℝ) :
    ((gMat u)⁻¹ *ᵥ v) ⬝ᵥ ((gMat u)⁻¹ *ᵥ v) ≤ v ⬝ᵥ ((gMat u)⁻¹ *ᵥ v) := by
  set w := (gMat u)⁻¹ *ᵥ v with hw
  have hv : gMat u *ᵥ w = v := by
    rw [hw, Matrix.mulVec_mulVec, gMat_mul_inv, Matrix.one_mulVec]
  have h2 := quad_g u w
  have h3 : 0 ≤ (w ᵥ* uMat u) ⬝ᵥ (w ᵥ* uMat u) :=
    Finset.sum_nonneg fun i _ => mul_self_nonneg _
  rw [hv] at h2
  rw [dotProduct_comm v w]
  linarith

theorem ginvFun_quadratic_nonneg (u : Z n 2) (v : Fin n → ℝ) :
    0 ≤ ∑ a, ∑ b, (gMat u)⁻¹ a b * v a * v b := by
  rw [quad_inv_eq]
  refine le_trans ?_ (quad_inv_ge u v)
  exact Finset.sum_nonneg fun i _ => mul_self_nonneg _

theorem ginvFun_quadratic_le (u : Z n 2) (v : Fin n → ℝ) :
    ∑ a, ∑ b, (gMat u)⁻¹ a b * v a * v b ≤ ∑ a, v a ^ 2 := by
  have hQ := ginvFun_quadratic_nonneg u v
  have hge := quad_inv_ge u v
  rw [quad_inv_eq] at hQ ⊢
  set w := (gMat u)⁻¹ *ᵥ v with hw
  set Q := v ⬝ᵥ w with hQdef
  have hcs := Finset.sum_mul_sq_le_sq_mul_sq Finset.univ v w
  have hV : ∑ a, v a ^ 2 = v ⬝ᵥ v := by simp [dotProduct, sq]
  have hS : ∑ a, w a ^ 2 = w ⬝ᵥ w := by simp [dotProduct, sq]
  have hQ' : ∑ a, v a * w a = Q := rfl
  rw [hV, hS, hQ'] at hcs
  rw [hV]
  have hV0 : 0 ≤ v ⬝ᵥ v := Finset.sum_nonneg fun i _ => mul_self_nonneg _
  by_contra hlt
  push Not at hlt
  nlinarith [mul_le_mul_of_nonneg_left hge hV0]

/-! ### Part 3: global smoothness -/

theorem contDiff_coord {m : ℕ} (I : Fin m → Fin n) : ContDiff ℝ ∞ (fun u : Z n m => u I) :=
  (contDiff_piLp 2).1 contDiff_id I

theorem contDiff_gMat_entry (a b : Fin n) : ContDiff ℝ ∞ (fun u : Z n 2 => gMat u a b) := by
  unfold gMat
  refine contDiff_const.add (ContDiff.sum fun c _ => ?_)
  exact (contDiff_coord _).mul (contDiff_coord _)

theorem contDiff_det_of {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    {f : E → Matrix (Fin n) (Fin n) ℝ} (hf : ∀ i j, ContDiff ℝ ∞ fun x => f x i j) :
    ContDiff ℝ ∞ fun x => (f x).det := by
  simp only [Matrix.det_apply, Units.smul_def, zsmul_eq_mul]
  refine ContDiff.sum fun σ _ => contDiff_const.mul ?_
  exact contDiff_prod fun i _ => hf _ _

theorem contDiff_adjugate_of {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    {f : E → Matrix (Fin n) (Fin n) ℝ} (hf : ∀ i j, ContDiff ℝ ∞ fun x => f x i j) (i j : Fin n) :
    ContDiff ℝ ∞ fun x => (f x).adjugate i j := by
  simp only [Matrix.adjugate_apply]
  refine contDiff_det_of fun k l => ?_
  simp only [Matrix.updateRow_apply]
  by_cases h : k = j
  · simp only [h, if_true]
    exact contDiff_const
  · simp only [h, if_false]
    exact hf _ _

theorem contDiff_det_gMat : ContDiff ℝ ∞ (fun u : Z n 2 => (gMat u).det) :=
  contDiff_det_of contDiff_gMat_entry

theorem contDiff_ginv_entry (a b : Fin n) : ContDiff ℝ ∞ (fun u : Z n 2 => (gMat u)⁻¹ a b) := by
  have : ∀ u : Z n 2, (gMat u)⁻¹ a b = ((gMat u).det)⁻¹ * (gMat u).adjugate a b := by
    intro u
    rw [Matrix.inv_def, Ring.inverse_eq_inv']
    rfl
  simp only [this]
  exact (contDiff_det_gMat.inv fun u => (det_gMat_pos u).ne').mul
    (contDiff_adjugate_of contDiff_gMat_entry a b)

theorem contDiff_gFun : ContDiff ℝ ∞ (gFun : Z n 2 → Z n 2) :=
  contDiff_piLp' (𝕜 := ℝ) 2 fun I => contDiff_gMat_entry _ _

theorem contDiff_ginvFun : ContDiff ℝ ∞ (ginvFun : Z n 2 → Z n 2) :=
  contDiff_piLp' (𝕜 := ℝ) 2 fun I => contDiff_ginv_entry _ _

theorem contDiff_volFun : ContDiff ℝ ∞ (volFun : Z n 2 → ℝ) := by
  refine contDiff_iff_contDiffAt.2 fun u => ?_
  exact (Real.contDiffAt_sqrt (det_gMat_pos u).ne').comp u contDiff_det_gMat.contDiffAt

theorem contDiff_gamFun : ContDiff ℝ ∞ (gamFun : Z n 2 × Z n 3 → Z n 3) := by
  refine contDiff_piLp' (𝕜 := ℝ) 2 fun I => ?_
  show ContDiff ℝ ∞ (fun q : Z n 2 × Z n 3 =>
    ∑ d, ∑ f, (gMat q.1)⁻¹ (I 2) d * q.2 ![I 0, I 1, f] * q.1 ![d, f])
  have hfst : ContDiff ℝ ∞ (Prod.fst : Z n 2 × Z n 3 → Z n 2) := contDiff_fst
  have hsnd : ContDiff ℝ ∞ (Prod.snd : Z n 2 × Z n 3 → Z n 3) := contDiff_snd
  refine ContDiff.sum (s := (Finset.univ : Finset (Fin n))) fun d _ =>
    ContDiff.sum (s := (Finset.univ : Finset (Fin n))) fun f _ => ?_
  have h1 : ContDiff ℝ ∞ (fun q : Z n 2 × Z n 3 => (gMat q.1)⁻¹ (I 2) d) :=
    (contDiff_ginv_entry (I 2) d).comp hfst
  have h2 : ContDiff ℝ ∞ (fun q : Z n 2 × Z n 3 => q.2 ![I 0, I 1, f]) :=
    (contDiff_coord _).comp hsnd
  have h3 : ContDiff ℝ ∞ (fun q : Z n 2 × Z n 3 => q.1 ![d, f]) :=
    (contDiff_coord _).comp hfst
  exact (h1.mul h2).mul h3

/-! ### Part 4: symbols of weight zero -/

section Symbols

variable {A : Set (P n N)}

theorem isSymbol_g (hA : IsOpen A) (hb : Bornology.IsBounded (low23 '' A)) :
    IsSymbol A 0 0 (fun p : P n N => gFun (low23 p).1) :=
  isSymbol_comp_low hA low23 low23_high hb
    (c := fun q : Z n 2 × Z n 3 => gFun q.1) (contDiff_gFun.comp contDiff_fst)

theorem isSymbol_ginv (hA : IsOpen A) (hb : Bornology.IsBounded (low23 '' A)) :
    IsSymbol A 0 0 (fun p : P n N => ginvFun (low23 p).1) :=
  isSymbol_comp_low hA low23 low23_high hb
    (c := fun q : Z n 2 × Z n 3 => ginvFun q.1) (contDiff_ginvFun.comp contDiff_fst)

theorem isSymbol_vol (hA : IsOpen A) (hb : Bornology.IsBounded (low23 '' A)) :
    IsSymbol A 0 0 (fun p : P n N => volFun (low23 p).1) :=
  isSymbol_comp_low hA low23 low23_high hb
    (c := fun q : Z n 2 × Z n 3 => volFun q.1) (contDiff_volFun.comp contDiff_fst)

theorem isSymbol_gam (hA : IsOpen A) (hb : Bornology.IsBounded (low23 '' A)) :
    IsSymbol A 0 0 (fun p : P n N => gamFun (low23 p)) :=
  isSymbol_comp_low hA low23 low23_high hb (c := gamFun) contDiff_gamFun

end Symbols

/-- info: 'CWFlow.Extension.contDiff_ginvFun' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms contDiff_ginvFun

/-- info: 'CWFlow.Extension.ginvFun_quadratic_le' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms ginvFun_quadratic_le

/-- info: 'CWFlow.Extension.one_le_det_gMat' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms one_le_det_gMat

/-- info: 'CWFlow.Extension.contDiff_gFun' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms contDiff_gFun

/-- info: 'CWFlow.Extension.contDiff_volFun' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms contDiff_volFun

/-- info: 'CWFlow.Extension.contDiff_gamFun' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms contDiff_gamFun

/-- info: 'CWFlow.Extension.isSymbol_g' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms isSymbol_g

/-- info: 'CWFlow.Extension.isSymbol_ginv' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms isSymbol_ginv

/-- info: 'CWFlow.Extension.isSymbol_vol' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms isSymbol_vol

/-- info: 'CWFlow.Extension.isSymbol_gam' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms isSymbol_gam

end CWFlow.Extension
