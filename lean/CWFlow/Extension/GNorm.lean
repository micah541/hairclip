import CWFlow.Extension.Symbols
import Mathlib.Analysis.Matrix.Order
import Mathlib.LinearAlgebra.Matrix.PosDef
import Mathlib.LinearAlgebra.SesquilinearForm.Basic
import Mathlib.Analysis.InnerProductSpace.PiL2
import Mathlib.Analysis.Calculus.ContDiff.WithLp
import Mathlib.Analysis.Calculus.ContDiff.FiniteDimension

/-!
# The norm of a tensor with respect to an inverse metric

`gNormSq h T = Σ_{I,J} (Π_i h(I_i, J_i)) T_I T_J`, the squared norm of the tensor `T` with `m`
indices for the metric whose inverse is `h`.  For a symmetric positive semidefinite `h` bounded by
the identity (resp. below by `lam`) it is a positive semidefinite symmetric form, bounded by `‖T‖²`
(resp. below by `lam ^ m ‖T‖²`).  Everything is first proved for plain coordinate functions
`(Fin m → Fin n) → ℝ`, by induction on `m` splitting off the first index; the key analytic input is
the Schur product theorem (`Matrix.PosSemidef.hadamard`).
-/

open scoped ContDiff Matrix

namespace CWFlow.Extension

variable {n : ℕ}

section Plain

variable (H : Fin n → Fin n → ℝ)

/-- The form on plain coordinate functions. -/
noncomputable def pI {m : ℕ} (S T : (Fin m → Fin n) → ℝ) : ℝ :=
  ∑ I : Fin m → Fin n, ∑ J : Fin m → Fin n, (∏ i, H (I i) (J i)) * S I * T J

/-- The slice of a tensor with first index `a`. -/
def slice {m : ℕ} (S : (Fin (m + 1) → Fin n) → ℝ) (a : Fin n) : (Fin m → Fin n) → ℝ :=
  fun I => S (Fin.cons a I)

theorem sum_cons' (j : ℕ) (g : (Fin (j + 1) → Fin n) → ℝ) :
    ∑ J, g J = ∑ i : Fin n, ∑ I : Fin j → Fin n, g (Fin.cons i I) := by
  rw [← (Fin.consEquiv (fun _ : Fin (j + 1) => Fin n)).sum_comp, Fintype.sum_prod_type]
  rfl

theorem pI_succ {m : ℕ} (S T : (Fin (m + 1) → Fin n) → ℝ) :
    pI H S T = ∑ a, ∑ b, H a b * pI H (slice S a) (slice T b) := by
  unfold pI
  simp only [sum_cons' m, Finset.mul_sum]
  refine Finset.sum_congr rfl fun a _ => ?_
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun b _ => ?_
  refine Finset.sum_congr rfl fun I _ => ?_
  refine Finset.sum_congr rfl fun J _ => ?_
  simp only [Fin.prod_univ_succ, Fin.cons_zero, Fin.cons_succ, slice]
  ring

theorem pI_zero (S T : (Fin 0 → Fin n) → ℝ) :
    pI H S T = S default * T default := by
  simp [pI]

theorem pI_symm (hs : ∀ a b, H a b = H b a) {m : ℕ} (S T : (Fin m → Fin n) → ℝ) :
    pI H S T = pI H T S := by
  unfold pI
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun I _ => Finset.sum_congr rfl fun J _ => ?_
  have : ∏ i, H (J i) (I i) = ∏ i, H (I i) (J i) :=
    Finset.prod_congr rfl fun i _ => hs _ _
  rw [this]; ring

/-- The form as a bilinear map. -/
noncomputable def pB (m : ℕ) :
    ((Fin m → Fin n) → ℝ) →ₗ[ℝ] ((Fin m → Fin n) → ℝ) →ₗ[ℝ] ℝ :=
  LinearMap.mk₂ ℝ (pI H)
    (fun S₁ S₂ T => by
      simp only [pI, Pi.add_apply, Finset.sum_add_distrib, mul_add, add_mul])
    (fun c S T => by
      simp only [pI, Pi.smul_apply, smul_eq_mul, Finset.mul_sum]
      refine Finset.sum_congr rfl fun I _ => Finset.sum_congr rfl fun J _ => ?_
      ring)
    (fun S T₁ T₂ => by
      simp only [pI, Pi.add_apply, Finset.sum_add_distrib, mul_add])
    (fun c S T => by
      simp only [pI, Pi.smul_apply, smul_eq_mul, Finset.mul_sum]
      refine Finset.sum_congr rfl fun I _ => Finset.sum_congr rfl fun J _ => ?_
      ring)

theorem pB_apply {m : ℕ} (S T : (Fin m → Fin n) → ℝ) : pB H m S T = pI H S T := rfl

theorem pI_comb {m : ℕ} (x : Fin n → ℝ) (f : Fin n → (Fin m → Fin n) → ℝ) :
    pI H (∑ a, x a • f a) (∑ b, x b • f b) = ∑ a, ∑ b, pI H (f a) (f b) * x a * x b := by
  rw [← pB_apply]
  simp only [map_sum, map_smul, LinearMap.sum_apply, LinearMap.smul_apply, pB_apply,
    smul_eq_mul]
  simp only [Finset.mul_sum]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun b _ => ?_
  ring

theorem psd_of {K : Fin n → Fin n → ℝ} (hs : ∀ a b, K a b = K b a)
    (hn : ∀ v : Fin n → ℝ, 0 ≤ ∑ a, ∑ b, K a b * v a * v b) : (Matrix.of K).PosSemidef := by
  rw [Matrix.posSemidef_iff_dotProduct_mulVec]
  refine ⟨?_, fun x => ?_⟩
  · ext a b
    simp [Matrix.conjTranspose_apply, hs a b]
  · have := hn x
    refine this.trans_eq ?_
    simp only [dotProduct, Matrix.mulVec, Matrix.of_apply, star_trivial, Finset.mul_sum]
    refine Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun b _ => ?_
    ring

theorem sum_mul_nonneg {K G : Fin n → Fin n → ℝ} (hK : (Matrix.of K).PosSemidef)
    (hG : (Matrix.of G).PosSemidef) : 0 ≤ ∑ a, ∑ b, K a b * G a b := by
  have := (hK.hadamard hG).dotProduct_mulVec_nonneg (fun _ => 1)
  simpa [dotProduct, Matrix.mulVec, Matrix.hadamard] using this

theorem gram_psd (hs : ∀ a b, H a b = H b a) {m : ℕ}
    (hpos : ∀ T : (Fin m → Fin n) → ℝ, 0 ≤ pI H T T) (f : Fin n → (Fin m → Fin n) → ℝ) :
    (Matrix.of fun a b => pI H (f a) (f b)).PosSemidef :=
  psd_of (fun a b => pI_symm H hs _ _) (fun v => by
    rw [← pI_comb]; exact hpos _)

theorem pI_nonneg (hs : ∀ a b, H a b = H b a)
    (hn : ∀ v : Fin n → ℝ, 0 ≤ ∑ a, ∑ b, H a b * v a * v b) :
    ∀ (m : ℕ) (T : (Fin m → Fin n) → ℝ), 0 ≤ pI H T T
  | 0, T => by rw [pI_zero]; exact mul_self_nonneg _
  | m + 1, T => by
    rw [pI_succ]
    exact sum_mul_nonneg (psd_of hs hn) (gram_psd H hs (pI_nonneg hs hn m) (slice T))

theorem sum_slice_sq {m : ℕ} (T : (Fin (m + 1) → Fin n) → ℝ) :
    ∑ I, T I ^ 2 = ∑ a, ∑ I, slice T a I ^ 2 := sum_cons' m _

theorem pI_le (hs : ∀ a b, H a b = H b a)
    (hn : ∀ v : Fin n → ℝ, 0 ≤ ∑ a, ∑ b, H a b * v a * v b)
    (hle : ∀ v : Fin n → ℝ, ∑ a, ∑ b, H a b * v a * v b ≤ ∑ a, v a ^ 2) :
    ∀ (m : ℕ) (T : (Fin m → Fin n) → ℝ), pI H T T ≤ ∑ I, T I ^ 2
  | 0, T => by
    rw [pI_zero]; simp [sq]
  | m + 1, T => by
    have hG := gram_psd H hs (pI_nonneg H hs hn m) (slice T)
    have hK : (Matrix.of fun a b => (if a = b then (1 : ℝ) else 0) - H a b).PosSemidef := by
      refine psd_of (fun a b => ?_) (fun v => ?_)
      · rw [hs a b]; simp [eq_comm]
      · have := hle v
        have e : ∑ a, ∑ b, ((if a = b then (1 : ℝ) else 0) - H a b) * v a * v b
            = ∑ a, v a ^ 2 - ∑ a, ∑ b, H a b * v a * v b := by
          simp [sub_mul, Finset.sum_sub_distrib, ite_mul, sq]
        rw [e]; linarith
    have h1 := sum_mul_nonneg hK hG
    have h2 : ∑ a, ∑ b, ((if a = b then (1 : ℝ) else 0) - H a b) *
        pI H (slice T a) (slice T b)
        = ∑ a, pI H (slice T a) (slice T a) - ∑ a, ∑ b, H a b * pI H (slice T a) (slice T b) := by
      simp [sub_mul, Finset.sum_sub_distrib, ite_mul]
    have h3 : ∑ a, pI H (slice T a) (slice T a) ≤ ∑ a, ∑ I, slice T a I ^ 2 :=
      Finset.sum_le_sum fun a _ => pI_le hs hn hle m (slice T a)
    rw [pI_succ, sum_slice_sq]
    linarith

theorem le_pI {lam : ℝ} (hlam : 0 ≤ lam) (hs : ∀ a b, H a b = H b a)
    (hn : ∀ v : Fin n → ℝ, 0 ≤ ∑ a, ∑ b, H a b * v a * v b)
    (hge : ∀ v : Fin n → ℝ, lam * ∑ a, v a ^ 2 ≤ ∑ a, ∑ b, H a b * v a * v b) :
    ∀ (m : ℕ) (T : (Fin m → Fin n) → ℝ), lam ^ m * ∑ I, T I ^ 2 ≤ pI H T T
  | 0, T => by
    rw [pI_zero]; simp [sq]
  | m + 1, T => by
    have hG := gram_psd H hs (pI_nonneg H hs hn m) (slice T)
    have hK : (Matrix.of fun a b => H a b - (if a = b then lam else 0)).PosSemidef := by
      refine psd_of (fun a b => ?_) (fun v => ?_)
      · rw [hs a b]; simp [eq_comm]
      · have := hge v
        have e : ∑ a, ∑ b, (H a b - (if a = b then lam else 0)) * v a * v b
            = ∑ a, ∑ b, H a b * v a * v b - lam * ∑ a, v a ^ 2 := by
          simp [sub_mul, Finset.sum_sub_distrib, ite_mul, sq, Finset.mul_sum, mul_assoc]
        rw [e]; linarith
    have h1 := sum_mul_nonneg hK hG
    have h2 : ∑ a, ∑ b, (H a b - (if a = b then lam else 0)) *
        pI H (slice T a) (slice T b)
        = ∑ a, ∑ b, H a b * pI H (slice T a) (slice T b) - lam * ∑ a, pI H (slice T a) (slice T a) := by
      simp [sub_mul, Finset.sum_sub_distrib, ite_mul, Finset.mul_sum]
    have h3 : ∑ a, lam ^ m * ∑ I, slice T a I ^ 2 ≤ ∑ a, pI H (slice T a) (slice T a) :=
      Finset.sum_le_sum fun a _ => le_pI hlam hs hn hge m (slice T a)
    have h6 : ∑ a, lam ^ m * ∑ I, slice T a I ^ 2 = lam ^ m * ∑ I, T I ^ 2 := by
      rw [← Finset.mul_sum, sum_slice_sq]
    have h5 := mul_le_mul_of_nonneg_left h3 hlam
    rw [h6] at h5
    have h7 : lam ^ (m + 1) * ∑ I, T I ^ 2 = lam * (lam ^ m * ∑ I, T I ^ 2) := by ring
    rw [pI_succ, h7]
    linarith

end Plain

section Tensor

/-- `|T|²_h = Σ_{I,J} (Π_i h(I_i, J_i)) T_I T_J`. -/
noncomputable def gNormSq {m : ℕ} (h : Z n 2) (T : Z n m) : ℝ :=
  ∑ I : Fin m → Fin n, ∑ J : Fin m → Fin n, (∏ i, h ![I i, J i]) * T I * T J

/-- The corresponding bilinear form. -/
noncomputable def gInner {m : ℕ} (h : Z n 2) (S T : Z n m) : ℝ :=
  ∑ I : Fin m → Fin n, ∑ J : Fin m → Fin n, (∏ i, h ![I i, J i]) * S I * T J

/-- The matrix `h` as a plain function. -/
def hh (h : Z n 2) : Fin n → Fin n → ℝ := fun a b => h ![a, b]

theorem gInner_eq_pI {m : ℕ} (h : Z n 2) (S T : Z n m) :
    gInner h S T = pI (hh h) (fun I => S I) (fun I => T I) := rfl

theorem gNormSq_eq_gInner {m : ℕ} (h : Z n 2) (T : Z n m) : gNormSq h T = gInner h T T := rfl

theorem normSq_eq_sum {m : ℕ} (T : Z n m) : ‖T‖ ^ 2 = ∑ I : Fin m → Fin n, (T I) ^ 2 := by
  rw [PiLp.norm_sq_eq_of_L2]
  simp [Real.norm_eq_abs, sq_abs]

theorem gNormSq_nonneg {m : ℕ} {h : Z n 2} (hsymm : ∀ a b, h ![a, b] = h ![b, a])
    (hnn : ∀ v : Fin n → ℝ, 0 ≤ ∑ a, ∑ b, h ![a, b] * v a * v b) (T : Z n m) :
    0 ≤ gNormSq h T :=
  pI_nonneg (hh h) hsymm hnn m _

theorem gNormSq_le {m : ℕ} {h : Z n 2} (hsymm : ∀ a b, h ![a, b] = h ![b, a])
    (hnn : ∀ v : Fin n → ℝ, 0 ≤ ∑ a, ∑ b, h ![a, b] * v a * v b)
    (hle : ∀ v : Fin n → ℝ, ∑ a, ∑ b, h ![a, b] * v a * v b ≤ ∑ a, v a ^ 2) (T : Z n m) :
    gNormSq h T ≤ ‖T‖ ^ 2 := by
  rw [normSq_eq_sum]
  exact pI_le (hh h) hsymm hnn hle m _

theorem le_gNormSq {m : ℕ} {h : Z n 2} {lam : ℝ} (hlam : 0 ≤ lam)
    (hsymm : ∀ a b, h ![a, b] = h ![b, a])
    (hnn : ∀ v : Fin n → ℝ, 0 ≤ ∑ a, ∑ b, h ![a, b] * v a * v b)
    (hge : ∀ v : Fin n → ℝ, lam * ∑ a, v a ^ 2 ≤ ∑ a, ∑ b, h ![a, b] * v a * v b) (T : Z n m) :
    lam ^ m * ‖T‖ ^ 2 ≤ gNormSq h T := by
  rw [normSq_eq_sum]
  exact le_pI (hh h) hlam hsymm hnn hge m _

/-- Cauchy–Schwarz for the form. -/
theorem gInner_sq_le {m : ℕ} {h : Z n 2} (hsymm : ∀ a b, h ![a, b] = h ![b, a])
    (hnn : ∀ v : Fin n → ℝ, 0 ≤ ∑ a, ∑ b, h ![a, b] * v a * v b) (S T : Z n m) :
    gInner h S T ^ 2 ≤ gNormSq h S * gNormSq h T := by
  have hB : (pB (hh h) m).IsSymm := ⟨fun x y => by
    simp only [pB_apply]; exact pI_symm (hh h) hsymm _ _⟩
  have := LinearMap.BilinForm.apply_sq_le_of_symm (pB (hh h) m)
    (fun x => pI_nonneg (hh h) hsymm hnn m x) hB (fun I => S I) (fun I => T I)
  exact this

theorem contDiff_coordG {ι : Type} [Fintype ι] (i : ι) :
    ContDiff ℝ ∞ (fun x : PiLp 2 (fun _ : ι => ℝ) => x i) :=
  contDiff_piLp_apply 2 (𝕜 := ℝ) (E := fun _ : ι => ℝ) (i := i)

theorem contDiff_gNormSq {m : ℕ} :
    ContDiff ℝ ∞ (fun q : Z n 2 × Z n m => gNormSq q.1 q.2) := by
  unfold gNormSq
  refine ContDiff.sum fun I _ => ContDiff.sum fun J _ => ?_
  refine ContDiff.mul (ContDiff.mul ?_ ?_) ?_
  · exact contDiff_prod fun i _ =>
      (contDiff_coordG ![I i, J i]).comp contDiff_fst
  · exact (contDiff_coordG I).comp contDiff_snd
  · exact (contDiff_coordG J).comp contDiff_snd

/-- The elementary bilinear form `S T ↦ S I * T J`. -/
noncomputable def gE {m : ℕ} (I J : Fin m → Fin n) : Z n m →L[ℝ] Z n m →L[ℝ] ℝ :=
  (PiLp.proj 2 (fun _ : Fin m → Fin n => ℝ) I).smulRight
    (PiLp.proj 2 (fun _ : Fin m → Fin n => ℝ) J)

/-- The form as a continuous bilinear map depending smoothly on `h`. -/
noncomputable def gForm {m : ℕ} (h : Z n 2) : Z n m →L[ℝ] Z n m →L[ℝ] ℝ :=
  ∑ I : Fin m → Fin n, ∑ J : Fin m → Fin n, (∏ i, h ![I i, J i]) • gE I J

theorem gForm_apply {m : ℕ} (h : Z n 2) (S T : Z n m) : gForm h S T = gInner h S T := by
  simp [gForm, gE, gInner, mul_assoc]

theorem contDiff_gInner_left {m : ℕ} (S T : Z n m) :
    ContDiff ℝ ∞ (fun h : Z n 2 => gInner h S T) := by
  unfold gInner
  refine ContDiff.sum fun I _ => ContDiff.sum fun J _ => ?_
  exact ((contDiff_prod fun i _ => contDiff_coordG ![I i, J i]).mul contDiff_const).mul
    contDiff_const

theorem contDiff_gForm {m : ℕ} : ContDiff ℝ ∞ (gForm (n := n) (m := m)) := by
  rw [contDiff_clm_apply_iff]
  intro S
  rw [contDiff_clm_apply_iff]
  intro T
  simp_rw [gForm_apply]
  exact contDiff_gInner_left S T

end Tensor

/-- info: 'CWFlow.Extension.gNormSq_nonneg' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms gNormSq_nonneg

/-- info: 'CWFlow.Extension.gNormSq_le' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms gNormSq_le

/-- info: 'CWFlow.Extension.le_gNormSq' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms le_gNormSq

/-- info: 'CWFlow.Extension.contDiff_gNormSq' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms contDiff_gNormSq

end CWFlow.Extension
