import CWFlow.Extension.ImmCoord
import CWFlow.Extension.TensorOps
import CWFlow.Extension.GNorm
import CWFlow.Interp.Jets

/-!
# Pullback of a tensor by a matrix in every slot

Finite-sum algebra for `pullZ A T (I) = Σ_K (Π_l A(K_l, I_l)) T(K)`: elementary properties, the
product rule, the behaviour of `gNormSq` and of `actG` under the pullback.
-/

open scoped ContDiff Matrix

namespace CWFlow.Extension

variable {n : ℕ}

/-- A matrix as a 2-tensor. -/
noncomputable def toZ2 (M : Matrix (Fin n) (Fin n) ℝ) : Z n 2 :=
  WithLp.toLp 2 (fun I => M (I 0) (I 1))

theorem toZ2_apply (M : Matrix (Fin n) (Fin n) ℝ) (a b : Fin n) : toZ2 M ![a, b] = M a b := rfl

/-! ### Part 1 -/

theorem pullZ_one {m : ℕ} (T : Z n m) : pullZ (1 : Matrix (Fin n) (Fin n) ℝ) T = T := by
  ext I
  rw [pullZ_apply]
  have h : ∀ K : Fin m → Fin n,
      (∏ l, (1 : Matrix (Fin n) (Fin n) ℝ) (K l) (I l)) = if K = I then 1 else 0 := by
    intro K
    simp only [Matrix.one_apply]
    rw [Finset.prod_boole]
    simp [funext_iff]
  simp [h]

theorem pullZ_add {m : ℕ} (A : Matrix (Fin n) (Fin n) ℝ) (S T : Z n m) :
    pullZ A (S + T) = pullZ A S + pullZ A T := by
  ext I
  simp [pullZ_apply, mul_add, Finset.sum_add_distrib]

theorem pullZ_sub {m : ℕ} (A : Matrix (Fin n) (Fin n) ℝ) (S T : Z n m) :
    pullZ A (S - T) = pullZ A S - pullZ A T := by
  ext I
  simp [pullZ_apply, mul_sub, Finset.sum_sub_distrib]

theorem pullZ_smul {m : ℕ} (A : Matrix (Fin n) (Fin n) ℝ) (c : ℝ) (T : Z n m) :
    pullZ A (c • T) = c • pullZ A T := by
  ext I
  simp [pullZ_apply, Finset.mul_sum, mul_left_comm]

/-- The continuous linear map `T ↦ pullZ A T`. -/
noncomputable def pullZL {m : ℕ} (A : Matrix (Fin n) (Fin n) ℝ) : Z n m →L[ℝ] Z n m :=
  LinearMap.toContinuousLinearMap
    { toFun := fun T => pullZ A T
      map_add' := pullZ_add A
      map_smul' := fun c T => pullZ_smul A c T }

theorem pullZL_apply {m : ℕ} (A : Matrix (Fin n) (Fin n) ℝ) (T : Z n m) :
    pullZL A T = pullZ A T := rfl

/-- Splitting off the first index. -/
theorem pullZ_cons {m : ℕ} (A : Matrix (Fin n) (Fin n) ℝ) (S : Z n (m + 1)) (a : Fin n)
    (I : Fin m → Fin n) :
    pullZ A S (Fin.cons a I) =
      ∑ a' : Fin n, A a' a * ∑ K : Fin m → Fin n, (∏ l, A (K l) (I l)) * S (Fin.cons a' K) := by
  rw [pullZ_apply, sum_cons' m]
  refine Finset.sum_congr rfl fun a' _ => ?_
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl fun K _ => ?_
  rw [Fin.prod_univ_succ]
  simp only [Fin.cons_zero, Fin.cons_succ]
  ring

/-! ### Part 2: the derivative -/

theorem fderiv_pullZ_apply {m : ℕ} {A : EuclideanSpace ℝ (Fin n) → Matrix (Fin n) (Fin n) ℝ}
    {T : EuclideanSpace ℝ (Fin n) → Z n m} {z : EuclideanSpace ℝ (Fin n)}
    (hA : ∀ j i : Fin n, DifferentiableAt ℝ (fun y => A y j i) z) (hT : DifferentiableAt ℝ T z)
    (v : EuclideanSpace ℝ (Fin n)) (I : Fin m → Fin n) :
    fderiv ℝ (fun y => pullZ (A y) (T y) I) z v =
      (∑ K : Fin m → Fin n, (∑ i : Fin m,
          fderiv ℝ (fun y => A y (K i) (I i)) z v * ∏ l ∈ Finset.univ.erase i, A z (K l) (I l)) *
            T z K)
      + pullZ (A z) (fderiv ℝ T z v) I := by
  have hTK : ∀ K : Fin m → Fin n, DifferentiableAt ℝ (fun y => T y K) z := fun K =>
    (PiLp.proj (𝕜 := ℝ) 2 (fun _ : Fin m → Fin n => ℝ) K).differentiableAt.comp z hT
  have hP : ∀ K : Fin m → Fin n, HasFDerivAt (fun y => ∏ l, A y (K l) (I l))
      (∑ i : Fin m, (∏ l ∈ Finset.univ.erase i, A z (K l) (I l)) •
        fderiv ℝ (fun y => A y (K i) (I i)) z) z := fun K =>
    HasFDerivAt.finsetProd (u := Finset.univ) (g := fun l y => A y (K l) (I l))
      (fun l _ => (hA _ _).hasFDerivAt)
  have hS : HasFDerivAt (fun y => pullZ (A y) (T y) I)
      (∑ K : Fin m → Fin n, ((∏ l, A z (K l) (I l)) • fderiv ℝ (fun y => T y K) z +
        T z K • ∑ i : Fin m, (∏ l ∈ Finset.univ.erase i, A z (K l) (I l)) •
          fderiv ℝ (fun y => A y (K i) (I i)) z)) z := by
    have := HasFDerivAt.fun_sum (u := (Finset.univ : Finset (Fin m → Fin n)))
      (A := fun K y => (∏ l, A y (K l) (I l)) * T y K)
      (A' := fun K => (∏ l, A z (K l) (I l)) • fderiv ℝ (fun y => T y K) z +
        T z K • ∑ i : Fin m, (∏ l ∈ Finset.univ.erase i, A z (K l) (I l)) •
          fderiv ℝ (fun y => A y (K i) (I i)) z)
      (fun K _ => (hP K).mul (hTK K).hasFDerivAt)
    exact this
  rw [hS.fderiv]
  simp only [ContinuousLinearMap.sum_apply, ContinuousLinearMap.add_apply,
    ContinuousLinearMap.smul_apply, smul_eq_mul, pullZ_apply, Finset.sum_add_distrib,
    Finset.mul_sum, Finset.sum_mul]
  rw [add_comm]
  congr 1
  · refine Finset.sum_congr rfl fun K _ => ?_
    refine Finset.sum_congr rfl fun i _ => ?_
    ring
  · refine Finset.sum_congr rfl fun K _ => ?_
    congr 1
    exact (CWFlow.Interp.fderiv_piLp_apply hT v K).symm

/-! ### Part 3: norms -/

theorem pz_sum_comm4 {α β γ δ : Type*} [Fintype α] [Fintype β] [Fintype γ] [Fintype δ]
    (g : α → β → γ → δ → ℝ) :
    ∑ I, ∑ I', ∑ K, ∑ K', g I I' K K' = ∑ K, ∑ K', ∑ I, ∑ I', g I I' K K' := by
  calc ∑ I, ∑ I', ∑ K, ∑ K', g I I' K K' = ∑ I, ∑ K, ∑ I', ∑ K', g I I' K K' :=
        Finset.sum_congr rfl fun I _ => Finset.sum_comm
    _ = ∑ I, ∑ K, ∑ K', ∑ I', g I I' K K' :=
        Finset.sum_congr rfl fun I _ => Finset.sum_congr rfl fun K _ => Finset.sum_comm
    _ = ∑ K, ∑ I, ∑ K', ∑ I', g I I' K K' := Finset.sum_comm
    _ = ∑ K, ∑ K', ∑ I, ∑ I', g I I' K K' :=
        Finset.sum_congr rfl fun K _ => Finset.sum_comm

theorem pz_prod_sum {m : ℕ} (g : Fin m → Fin n → ℝ) :
    ∏ l, ∑ a, g l a = ∑ I : Fin m → Fin n, ∏ l, g l (I l) := by
  rw [Finset.prod_univ_sum]
  simp

theorem pz_key {m : ℕ} (A h : Matrix (Fin n) (Fin n) ℝ) (K K' : Fin m → Fin n) :
    ∑ I : Fin m → Fin n, ∑ I' : Fin m → Fin n,
      (∏ l, h (I l) (I' l)) * ((∏ l, A (K l) (I l)) * (∏ l, A (K' l) (I' l))) =
    ∏ l, (A * h * Aᵀ) (K l) (K' l) := by
  have e1 : ∀ I I' : Fin m → Fin n, (∏ l, h (I l) (I' l)) *
      ((∏ l, A (K l) (I l)) * (∏ l, A (K' l) (I' l))) =
      ∏ l, (A (K l) (I l) * h (I l) (I' l) * A (K' l) (I' l)) := by
    intro I I'
    rw [Finset.prod_mul_distrib, Finset.prod_mul_distrib]
    ring
  simp_rw [e1]
  have e2 : ∀ l, (A * h * Aᵀ) (K l) (K' l) =
      ∑ a, ∑ b, A (K l) a * h a b * A (K' l) b := by
    intro l
    rw [Matrix.mul_apply]
    simp only [Matrix.mul_apply, Matrix.transpose_apply, Finset.sum_mul]
    exact Finset.sum_comm
  simp_rw [e2]
  rw [pz_prod_sum]
  refine Finset.sum_congr rfl fun I _ => ?_
  rw [← pz_prod_sum (fun l b => A (K l) (I l) * h (I l) b * A (K' l) b)]

theorem gNormSq_pullZ {m : ℕ} (A h : Matrix (Fin n) (Fin n) ℝ) (T : Z n m) :
    gNormSq (toZ2 h) (pullZ A T) = gNormSq (toZ2 (A * h * Aᵀ)) T := by
  unfold gNormSq
  simp only [pullZ_apply, toZ2_apply]
  have hr : ∀ K K' : Fin m → Fin n, (∏ l, (A * h * Aᵀ) (K l) (K' l)) * T K * T K' =
      ∑ I : Fin m → Fin n, ∑ I' : Fin m → Fin n,
        (∏ l, h (I l) (I' l)) * ((∏ l, A (K l) (I l)) * (∏ l, A (K' l) (I' l))) *
        T K * T K' := by
    intro K K'
    rw [← pz_key, Finset.sum_mul, Finset.sum_mul]
    refine Finset.sum_congr rfl fun I _ => ?_
    rw [Finset.sum_mul, Finset.sum_mul]
  simp_rw [hr]
  simp only [Finset.mul_sum, Finset.sum_mul]
  rw [pz_sum_comm4, Finset.sum_comm]
  refine Finset.sum_congr rfl fun K _ => Finset.sum_congr rfl fun K' _ => ?_
  refine Finset.sum_congr rfl fun I _ => Finset.sum_congr rfl fun I' _ => ?_
  ring

theorem gNormSq_pullZ_inv {m : ℕ} {A g g' : Matrix (Fin n) (Fin n) ℝ} (hA : IsUnit A.det)
    (hg : IsUnit g.det) (hg' : g' = Aᵀ * g * A) (T : Z n m) :
    gNormSq (toZ2 g'⁻¹) (pullZ A T) = gNormSq (toZ2 g⁻¹) T := by
  rw [gNormSq_pullZ]
  congr 2
  have hAt : IsUnit (Aᵀ).det := by rwa [Matrix.det_transpose]
  have h1 : g'⁻¹ = A⁻¹ * (g⁻¹ * (Aᵀ)⁻¹) := by
    rw [hg', Matrix.mul_inv_rev, Matrix.mul_inv_rev]
  rw [h1]
  calc A * (A⁻¹ * (g⁻¹ * (Aᵀ)⁻¹)) * Aᵀ = (A * A⁻¹) * g⁻¹ * ((Aᵀ)⁻¹ * Aᵀ) := by
        simp only [Matrix.mul_assoc]
    _ = g⁻¹ := by
        rw [Matrix.mul_nonsing_inv A hA, Matrix.nonsing_inv_mul _ hAt]
        simp

theorem sqrt_det_conj {A g g' : Matrix (Fin n) (Fin n) ℝ} (hg' : g' = Aᵀ * g * A) :
    Real.sqrt g'.det = |A.det| * Real.sqrt g.det := by
  have : g'.det = (A.det * A.det) * g.det := by
    rw [hg', Matrix.det_mul, Matrix.det_mul, Matrix.det_transpose]
    ring
  rw [this, Real.sqrt_mul (mul_self_nonneg _), Real.sqrt_mul_self_eq_abs]

/-! ### Part 4: the Christoffel correction -/

/-- The involution `(K, b) ↦ (update K i b, K i)`. -/
def pz_swapEquiv {m : ℕ} (i : Fin m) :
    ((Fin m → Fin n) × Fin n) ≃ ((Fin m → Fin n) × Fin n) where
  toFun p := (Function.update p.1 i p.2, p.1 i)
  invFun p := (Function.update p.1 i p.2, p.1 i)
  left_inv p := by
    obtain ⟨K, b⟩ := p
    simp
  right_inv p := by
    obtain ⟨K, b⟩ := p
    simp

theorem pz_swap {m : ℕ} (i : Fin m) (F : (Fin m → Fin n) → Fin n → ℝ) :
    ∑ K, ∑ b, F K b = ∑ K, ∑ c, F (Function.update K i c) (K i) := by
  rw [← Fintype.sum_prod_type' (fun K b => F K b), ← Fintype.sum_prod_type'
    (fun K c => F (Function.update K i c) (K i))]
  exact ((pz_swapEquiv (n := n) i).sum_comp (fun p => F p.1 p.2)).symm

theorem pz_prod_update {m : ℕ} (A : Matrix (Fin n) (Fin n) ℝ) (I K : Fin m → Fin n) (i : Fin m)
    (c : Fin n) :
    ∏ l, A (K l) (Function.update I i c l) =
      A (K i) c * ∏ l ∈ Finset.univ.erase i, A (K l) (I l) := by
  rw [← Finset.mul_prod_erase Finset.univ _ (Finset.mem_univ i), Function.update_self]
  congr 1
  refine Finset.prod_congr rfl fun l hl => ?_
  rw [Function.update_of_ne (Finset.ne_of_mem_erase hl)]

theorem pz_prod_erase_update {m : ℕ} (A : Matrix (Fin n) (Fin n) ℝ) (I K : Fin m → Fin n)
    (i : Fin m) (c : Fin n) :
    ∏ l ∈ Finset.univ.erase i, A (Function.update K i c l) (I l) =
      ∏ l ∈ Finset.univ.erase i, A (K l) (I l) := by
  refine Finset.prod_congr rfl fun l hl => ?_
  rw [Function.update_of_ne (Finset.ne_of_mem_erase hl)]

theorem pz_prod_split {m : ℕ} (A : Matrix (Fin n) (Fin n) ℝ) (I K : Fin m → Fin n) (i : Fin m) :
    (∏ l ∈ Finset.univ.erase i, A (K l) (I l)) * A (K i) (I i) = ∏ l, A (K l) (I l) :=
  Finset.prod_erase_mul _ _ (Finset.mem_univ i)

theorem pz_gprime {A B : Matrix (Fin n) (Fin n) ℝ} (hBA : B * A = 1) (G G' : Z n 3)
    (D : Fin n → Fin n → Fin n → ℝ)
    (hG' : ∀ a b c : Fin n, G' ![a, b, c] =
      ∑ c' : Fin n, B c c' * ((∑ a' : Fin n, ∑ b' : Fin n, A a' a * A b' b * G ![a', b', c']) +
        D a b c'))
    (a b k : Fin n) :
    ∑ c, A k c * G' ![a, b, c] =
      (∑ a' : Fin n, A a' a * ∑ b' : Fin n, A b' b * G ![a', b', k]) + D a b k := by
  have hAB : A * B = 1 := mul_eq_one_comm.mp hBA
  have h1 : ∀ c, A k c * G' ![a, b, c] = ∑ c', (A k c * B c c') *
      ((∑ a' : Fin n, ∑ b' : Fin n, A a' a * A b' b * G ![a', b', c']) + D a b c') := by
    intro c
    rw [hG', Finset.mul_sum]
    refine Finset.sum_congr rfl fun c' _ => ?_
    ring
  simp_rw [h1]
  rw [Finset.sum_comm]
  simp_rw [← Finset.sum_mul]
  have h2 : ∀ c', ∑ c, A k c * B c c' = if k = c' then 1 else 0 := by
    intro c'
    have := congrFun (congrFun hAB k) c'
    rw [Matrix.mul_apply] at this
    rw [this]
    simp [Matrix.one_apply]
  simp_rw [h2]
  simp only [ite_mul, one_mul, zero_mul, Finset.sum_ite_eq, Finset.mem_univ, if_true]
  congr 1
  refine Finset.sum_congr rfl fun a' _ => ?_
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl fun b' _ => ?_
  ring

theorem actG_pullZ {m : ℕ} {A B : Matrix (Fin n) (Fin n) ℝ} (hBA : B * A = 1) (G G' : Z n 3)
    (D : Fin n → Fin n → Fin n → ℝ)
    (hG' : ∀ a b c : Fin n, G' ![a, b, c] =
      ∑ c' : Fin n, B c c' * ((∑ a' : Fin n, ∑ b' : Fin n, A a' a * A b' b * G ![a', b', c']) +
        D a b c'))
    (T : Z n m) (a : Fin n) (I : Fin m → Fin n) :
    actG G' (pullZ A T) (Fin.cons a I) =
      pullZ A (actG G T) (Fin.cons a I)
      + ∑ K : Fin m → Fin n, (∑ i : Fin m,
          D a (I i) (K i) * ∏ l ∈ Finset.univ.erase i, A (K l) (I l)) * T K := by
  have hL : actG G' (pullZ A T) (Fin.cons a I) =
      ∑ i : Fin m, ∑ K : Fin m → Fin n,
        (∏ l ∈ Finset.univ.erase i, A (K l) (I l)) *
          ((∑ a' : Fin n, A a' a * ∑ b' : Fin n, A b' (I i) * G ![a', b', K i]) + D a (I i) (K i))
          * T K := by
    rw [actG_apply]
    simp only [Fin.cons_zero, Fin.tail_cons, pullZ_apply]
    refine Finset.sum_congr rfl fun i _ => ?_
    calc ∑ c, G' ![a, I i, c] *
          ∑ K : Fin m → Fin n, (∏ l, A (K l) (Function.update I i c l)) * T K
        = ∑ c, ∑ K : Fin m → Fin n, G' ![a, I i, c] *
            (A (K i) c * ∏ l ∈ Finset.univ.erase i, A (K l) (I l)) * T K := by
          refine Finset.sum_congr rfl fun c _ => ?_
          rw [Finset.mul_sum]
          refine Finset.sum_congr rfl fun K _ => ?_
          rw [pz_prod_update]
          ring
      _ = ∑ K : Fin m → Fin n, ∑ c, G' ![a, I i, c] *
            (A (K i) c * ∏ l ∈ Finset.univ.erase i, A (K l) (I l)) * T K := Finset.sum_comm
      _ = _ := by
          refine Finset.sum_congr rfl fun K _ => ?_
          rw [← pz_gprime hBA G G' D hG' a (I i) (K i), Finset.mul_sum, Finset.sum_mul]
          refine Finset.sum_congr rfl fun c _ => ?_
          ring
  have hR : pullZ A (actG G T) (Fin.cons a I) =
      ∑ i : Fin m, ∑ K : Fin m → Fin n,
        (∏ l ∈ Finset.univ.erase i, A (K l) (I l)) *
          (∑ a' : Fin n, A a' a * ∑ b' : Fin n, A b' (I i) * G ![a', b', K i]) * T K := by
    have hi : ∀ (a' : Fin n) (i : Fin m),
        ∑ K : Fin m → Fin n, ∑ b' : Fin n,
          (∏ l ∈ Finset.univ.erase i, A (K l) (I l)) * (A b' (I i) * G ![a', b', K i] * T K) =
        ∑ K : Fin m → Fin n, ∑ c : Fin n,
          (∏ l, A (K l) (I l)) * (G ![a', K i, c] * T (Function.update K i c)) := by
      intro a' i
      rw [pz_swap i (fun K b => (∏ l ∈ Finset.univ.erase i, A (K l) (I l)) *
        (A b (I i) * G ![a', b, K i] * T K))]
      refine Finset.sum_congr rfl fun K _ => Finset.sum_congr rfl fun c _ => ?_
      simp only [Function.update_self]
      rw [pz_prod_erase_update, ← pz_prod_split A I K i]
      ring
    have hM : ∀ i : Fin m, ∑ K : Fin m → Fin n,
        (∏ l ∈ Finset.univ.erase i, A (K l) (I l)) *
          (∑ a' : Fin n, A a' a * ∑ b' : Fin n, A b' (I i) * G ![a', b', K i]) * T K =
        ∑ a' : Fin n, A a' a * ∑ K : Fin m → Fin n, ∑ c : Fin n,
          (∏ l, A (K l) (I l)) * (G ![a', K i, c] * T (Function.update K i c)) := by
      intro i
      simp_rw [← hi]
      simp only [Finset.mul_sum, Finset.sum_mul]
      rw [Finset.sum_comm]
      refine Finset.sum_congr rfl fun a' _ => Finset.sum_congr rfl fun K _ =>
        Finset.sum_congr rfl fun b' _ => ?_
      ring
    simp_rw [hM]
    rw [pullZ_cons]
    simp only [actG_apply, Fin.cons_zero, Fin.tail_cons]
    simp only [Finset.mul_sum]
    calc _ = ∑ a' : Fin n, ∑ i : Fin m, ∑ K : Fin m → Fin n, ∑ c : Fin n,
          A a' a * ((∏ l, A (K l) (I l)) * (G ![a', K i, c] * T (Function.update K i c))) := by
          refine Finset.sum_congr rfl fun a' _ => ?_
          rw [Finset.sum_comm]
      _ = _ := Finset.sum_comm
  rw [hL, hR]
  simp only [mul_add, add_mul, Finset.sum_add_distrib]
  congr 1
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun K _ => ?_
  rw [Finset.sum_mul]
  refine Finset.sum_congr rfl fun i _ => ?_
  ring

/-- info: 'CWFlow.Extension.fderiv_pullZ_apply' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms fderiv_pullZ_apply

/-- info: 'CWFlow.Extension.gNormSq_pullZ' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms gNormSq_pullZ

/-- info: 'CWFlow.Extension.gNormSq_pullZ_inv' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms gNormSq_pullZ_inv

/-- info: 'CWFlow.Extension.actG_pullZ' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms actG_pullZ

/-- info: 'CWFlow.Extension.sqrt_det_conj' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms sqrt_det_conj

end CWFlow.Extension
