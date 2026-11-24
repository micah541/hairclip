import CWFlow.Extension.ImmReparam
import CWFlow.Extension.ImmDiv
import CWFlow.Extension.VolEvol

/-!
# The potential, the mean curvature form and the velocity under a change of coordinates

For a Lagrangian immersion `X` and a smooth map `τ` with invertible Jacobian `A = jacM τ z`:
`potX (X ∘ τ) z = potX X (τ z)`, `muC (X ∘ τ) z c = Σ_{c'} A_{c'c} muC X (τ z) c'` and
`velX (X ∘ τ) z = velX X (τ z)`. This is finite-sum algebra using `A g'⁻¹ Aᵀ = g⁻¹` for
`g' = Aᵀ g A`, plus the chain rule.
-/

open scoped ContDiff Matrix Topology RealInnerProductSpace

namespace CWFlow.Extension

variable {n : ℕ} {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V] {J : V →L[ℝ] V}
  {X : EuclideanSpace ℝ (Fin n) → V} {τ : EuclideanSpace ℝ (Fin n) → EuclideanSpace ℝ (Fin n)}
  {Ω N : Set (EuclideanSpace ℝ (Fin n))} {z : EuclideanSpace ℝ (Fin n)}

/-! ### Part 1: algebra -/

theorem pr_sum4 (g : (Fin 4 → Fin n) → ℝ) :
    ∑ K, g K = ∑ a, ∑ b, ∑ c, ∑ d, g ![a, b, c, d] := by
  rw [sum_cons' 3]
  refine Finset.sum_congr rfl fun a _ => ?_
  rw [sum_cons' 2]
  refine Finset.sum_congr rfl fun b _ => ?_
  rw [sum_cons' 1]
  refine Finset.sum_congr rfl fun c _ => ?_
  rw [sum_cons' 0]
  refine Finset.sum_congr rfl fun d _ => ?_
  rw [Fintype.sum_unique]
  congr 1

/-- `A g'⁻¹ Aᵀ = g⁻¹` entrywise. -/
theorem pr_key {A g g' : Matrix (Fin n) (Fin n) ℝ} (hA : IsUnit A.det) (hg' : g' = Aᵀ * g * A)
    (a' b' : Fin n) : ∑ a, ∑ b, A a' a * g'⁻¹ a b * A b' b = g⁻¹ a' b' := by
  have hAt : IsUnit (Aᵀ).det := by rwa [Matrix.det_transpose]
  have h1 : g'⁻¹ = A⁻¹ * (g⁻¹ * (Aᵀ)⁻¹) := by
    rw [hg', Matrix.mul_inv_rev, Matrix.mul_inv_rev]
  have h2 : A * g'⁻¹ * Aᵀ = g⁻¹ := by
    rw [h1]
    calc A * (A⁻¹ * (g⁻¹ * (Aᵀ)⁻¹)) * Aᵀ = (A * A⁻¹) * g⁻¹ * ((Aᵀ)⁻¹ * Aᵀ) := by
          simp only [Matrix.mul_assoc]
      _ = g⁻¹ := by
          rw [Matrix.mul_nonsing_inv A hA, Matrix.nonsing_inv_mul _ hAt]
          simp
  rw [← h2]
  simp only [Matrix.mul_apply, Matrix.transpose_apply, Finset.sum_mul]
  rw [Finset.sum_comm]

theorem pr_rev4 {M ι κ ρ σ : Type*} [AddCommMonoid M] [Fintype ι] [Fintype κ] [Fintype ρ]
    [Fintype σ] (f : ι → κ → ρ → σ → M) :
    ∑ a, ∑ b, ∑ c, ∑ d, f a b c d = ∑ d, ∑ c, ∑ a, ∑ b, f a b c d := by
  calc ∑ a, ∑ b, ∑ c, ∑ d, f a b c d = ∑ a, ∑ b, ∑ d, ∑ c, f a b c d :=
        Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun b _ => Finset.sum_comm
    _ = ∑ a, ∑ d, ∑ b, ∑ c, f a b c d := Finset.sum_congr rfl fun a _ => Finset.sum_comm
    _ = ∑ d, ∑ a, ∑ b, ∑ c, f a b c d := Finset.sum_comm
    _ = ∑ d, ∑ a, ∑ c, ∑ b, f a b c d :=
        Finset.sum_congr rfl fun d _ => Finset.sum_congr rfl fun a _ => Finset.sum_comm
    _ = ∑ d, ∑ c, ∑ a, ∑ b, f a b c d := Finset.sum_congr rfl fun d _ => Finset.sum_comm

theorem pr_rot4 {M ι κ ρ σ : Type*} [AddCommMonoid M] [Fintype ι] [Fintype κ] [Fintype ρ]
    [Fintype σ] (f : ι → κ → ρ → σ → M) :
    ∑ a, ∑ b, ∑ c, ∑ d, f a b c d = ∑ d, ∑ a, ∑ b, ∑ c, f a b c d := by
  calc ∑ a, ∑ b, ∑ c, ∑ d, f a b c d = ∑ a, ∑ b, ∑ d, ∑ c, f a b c d :=
        Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun b _ => Finset.sum_comm
    _ = ∑ a, ∑ d, ∑ b, ∑ c, f a b c d := Finset.sum_congr rfl fun a _ => Finset.sum_comm
    _ = ∑ d, ∑ a, ∑ b, ∑ c, f a b c d := Finset.sum_comm

theorem pr_rot5 {M ι κ ρ σ υ : Type*} [AddCommMonoid M] [Fintype ι] [Fintype κ] [Fintype ρ]
    [Fintype σ] [Fintype υ] (f : ι → κ → ρ → σ → υ → M) :
    ∑ a, ∑ b, ∑ c, ∑ d, ∑ e, f a b c d e = ∑ e, ∑ a, ∑ b, ∑ c, ∑ d, f a b c d e := by
  calc ∑ a, ∑ b, ∑ c, ∑ d, ∑ e, f a b c d e = ∑ a, ∑ b, ∑ c, ∑ e, ∑ d, f a b c d e :=
        Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun b _ =>
          Finset.sum_congr rfl fun c _ => Finset.sum_comm
    _ = ∑ a, ∑ b, ∑ e, ∑ c, ∑ d, f a b c d e :=
        Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun b _ => Finset.sum_comm
    _ = ∑ a, ∑ e, ∑ b, ∑ c, ∑ d, f a b c d e := Finset.sum_congr rfl fun a _ => Finset.sum_comm
    _ = ∑ e, ∑ a, ∑ b, ∑ c, ∑ d, f a b c d e := Finset.sum_comm

theorem pr_swap22 {M ι κ ρ σ : Type*} [AddCommMonoid M] [Fintype ι] [Fintype κ] [Fintype ρ]
    [Fintype σ] (f : ι → κ → ρ → σ → M) :
    ∑ a, ∑ b, ∑ c, ∑ d, f a b c d = ∑ c, ∑ d, ∑ a, ∑ b, f a b c d := by
  calc ∑ a, ∑ b, ∑ c, ∑ d, f a b c d = ∑ a, ∑ c, ∑ b, ∑ d, f a b c d :=
        Finset.sum_congr rfl fun a _ => Finset.sum_comm
    _ = ∑ c, ∑ a, ∑ b, ∑ d, f a b c d := Finset.sum_comm
    _ = ∑ c, ∑ a, ∑ d, ∑ b, f a b c d :=
        Finset.sum_congr rfl fun c _ => Finset.sum_congr rfl fun a _ => Finset.sum_comm
    _ = ∑ c, ∑ d, ∑ a, ∑ b, f a b c d := Finset.sum_congr rfl fun c _ => Finset.sum_comm

theorem pr_reorder (f : (Fin 4 → Fin n) → ℝ) :
    ∑ a, ∑ b, ∑ c, ∑ d, f ![d, a, b, c] = ∑ I, f I := by
  rw [pr_sum4, pr_rot4 (fun a b c d => f ![d, a, b, c])]

theorem pr_pull3 (A : Matrix (Fin n) (Fin n) ℝ) (T : Z n 3) (a b c : Fin n) :
    pullZ A T ![a, b, c] = ∑ a', ∑ b', ∑ c', A a' a * A b' b * A c' c * T ![a', b', c'] := by
  rw [pullZ_apply, ir_sum3]
  simp [Fin.prod_univ_three]

theorem pr_inner {A g g' : Matrix (Fin n) (Fin n) ℝ} (hA : IsUnit A.det) (hg' : g' = Aᵀ * g * A)
    (K : Fin 4 → Fin n) :
    ∑ I : Fin 4 → Fin n, g'⁻¹ (I 3) (I 0) * g'⁻¹ (I 1) (I 2) * ∏ l, A (K l) (I l) =
      g⁻¹ (K 3) (K 0) * g⁻¹ (K 1) (K 2) := by
  rw [← pr_key hA hg' (K 3) (K 0), ← pr_key hA hg' (K 1) (K 2), pr_sum4]
  simp only [Fin.prod_univ_four, Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two,
    Matrix.cons_val_three, Matrix.head_cons, Matrix.tail_cons, Finset.sum_mul, Finset.mul_sum]
  conv_rhs => rw [pr_rot4]
  refine Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun b _ =>
    Finset.sum_congr rfl fun c _ => Finset.sum_congr rfl fun d _ => ?_
  ring

theorem pr_contractI {A g g' : Matrix (Fin n) (Fin n) ℝ} (hA : IsUnit A.det)
    (hg' : g' = Aᵀ * g * A) (T : Z n 4) :
    ∑ I : Fin 4 → Fin n, g'⁻¹ (I 3) (I 0) * g'⁻¹ (I 1) (I 2) * pullZ A T I =
      ∑ I : Fin 4 → Fin n, g⁻¹ (I 3) (I 0) * g⁻¹ (I 1) (I 2) * T I := by
  calc ∑ I : Fin 4 → Fin n, g'⁻¹ (I 3) (I 0) * g'⁻¹ (I 1) (I 2) * pullZ A T I
      = ∑ I : Fin 4 → Fin n, ∑ K : Fin 4 → Fin n,
          g'⁻¹ (I 3) (I 0) * g'⁻¹ (I 1) (I 2) * ((∏ l, A (K l) (I l)) * T K) := by
        simp only [pullZ_apply, Finset.mul_sum]
    _ = ∑ K : Fin 4 → Fin n, ∑ I : Fin 4 → Fin n,
          g'⁻¹ (I 3) (I 0) * g'⁻¹ (I 1) (I 2) * ((∏ l, A (K l) (I l)) * T K) := Finset.sum_comm
    _ = ∑ K : Fin 4 → Fin n, (∑ I : Fin 4 → Fin n,
          g'⁻¹ (I 3) (I 0) * g'⁻¹ (I 1) (I 2) * ∏ l, A (K l) (I l)) * T K := by
        refine Finset.sum_congr rfl fun K _ => ?_
        rw [Finset.sum_mul]
        exact Finset.sum_congr rfl fun I _ => by ring
    _ = ∑ K : Fin 4 → Fin n, g⁻¹ (K 3) (K 0) * g⁻¹ (K 1) (K 2) * T K := by
        simp_rw [pr_inner hA hg']

theorem pr_contract {A g g' : Matrix (Fin n) (Fin n) ℝ} (hA : IsUnit A.det)
    (hg' : g' = Aᵀ * g * A) (T : Z n 4) :
    ∑ a, ∑ b, ∑ c, ∑ d, g'⁻¹ c d * g'⁻¹ a b * pullZ A T ![d, a, b, c] =
      ∑ a, ∑ b, ∑ c, ∑ d, g⁻¹ c d * g⁻¹ a b * T ![d, a, b, c] := by
  have h1 := pr_reorder (fun I : Fin 4 → Fin n => g'⁻¹ (I 3) (I 0) * g'⁻¹ (I 1) (I 2) *
    pullZ A T I)
  have h2 := pr_reorder (fun I : Fin 4 → Fin n => g⁻¹ (I 3) (I 0) * g⁻¹ (I 1) (I 2) * T I)
  simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two,
    Matrix.cons_val_three, Matrix.head_cons, Matrix.tail_cons] at h1 h2
  rw [h1, h2]
  exact pr_contractI hA hg' T

theorem pr_potX_eq (J : V →L[ℝ] V) (X : EuclideanSpace ℝ (Fin n) → V)
    (x : EuclideanSpace ℝ (Fin n)) :
    potX J X x = -∑ I : Fin 4 → Fin n,
      (gX X x)⁻¹ (I 3) (I 0) * (gX X x)⁻¹ (I 1) (I 2) * SX J X 1 x I := by
  unfold potX
  congr 1
  have h1 := pr_reorder (fun I : Fin 4 → Fin n => (gX X x)⁻¹ (I 3) (I 0) * (gX X x)⁻¹ (I 1) (I 2) *
    SX J X 1 x I)
  simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two,
    Matrix.cons_val_three, Matrix.head_cons, Matrix.tail_cons] at h1
  rw [← h1]
  rfl

theorem potX_comp (hΩ : IsOpen Ω) (hX : ContDiffOn ℝ ∞ X Ω)
    (hinj : ∀ x ∈ Ω, Function.Injective (fderiv ℝ X x))
    (hlag : ∀ x ∈ Ω, ∀ u w, ⟪J (fderiv ℝ X x u), fderiv ℝ X x w⟫ = 0) (hN : IsOpen N)
    (hτ : ContDiffOn ℝ ∞ τ N) (hmaps : ∀ z ∈ N, τ z ∈ Ω)
    (hdet : ∀ z ∈ N, (jacM τ z).det ≠ 0) (hz : z ∈ N) :
    potX J (X ∘ τ) z = potX J X (τ z) := by
  have hA : IsUnit (jacM τ z).det := isUnit_iff_ne_zero.2 (hdet z hz)
  have hg := gX_comp hΩ hX hN hτ hmaps hz
  rw [pr_potX_eq, pr_potX_eq, SX_comp hΩ hX hinj hlag hN hτ hmaps hdet 1 hz]
  congr 1
  exact pr_contractI hA hg _

theorem muC_comp (hΩ : IsOpen Ω) (hX : ContDiffOn ℝ ∞ X Ω)
    (hlag : ∀ x ∈ Ω, ∀ u w, ⟪J (fderiv ℝ X x u), fderiv ℝ X x w⟫ = 0) (hN : IsOpen N)
    (hτ : ContDiffOn ℝ ∞ τ N) (hmaps : ∀ z ∈ N, τ z ∈ Ω)
    (hdet : ∀ z ∈ N, (jacM τ z).det ≠ 0) (hz : z ∈ N) (c : Fin n) :
    muC J (X ∘ τ) z c = ∑ c', jacM τ z c' c * muC J X (τ z) c' := by
  have hA : IsUnit (jacM τ z).det := isUnit_iff_ne_zero.2 (hdet z hz)
  have hk := pr_key hA (gX_comp hΩ hX hN hτ hmaps hz)
  unfold muC
  rw [cX_comp hΩ hX hlag hN hτ hmaps hz]
  simp only [pr_pull3, ← hk, Finset.mul_sum, Finset.sum_mul, mul_neg, Finset.sum_neg_distrib,
    neg_inj]
  rw [pr_rot5]
  refine Finset.sum_congr rfl fun c' _ => ?_
  rw [pr_swap22]
  refine Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun b _ =>
    Finset.sum_congr rfl fun c _ => Finset.sum_congr rfl fun d _ => ?_
  ring

theorem velX_comp (hΩ : IsOpen Ω) (hX : ContDiffOn ℝ ∞ X Ω)
    (hinj : ∀ x ∈ Ω, Function.Injective (fderiv ℝ X x))
    (hlag : ∀ x ∈ Ω, ∀ u w, ⟪J (fderiv ℝ X x u), fderiv ℝ X x w⟫ = 0) (hN : IsOpen N)
    (hτ : ContDiffOn ℝ ∞ τ N) (hmaps : ∀ z ∈ N, τ z ∈ Ω)
    (hdet : ∀ z ∈ N, (jacM τ z).det ≠ 0) (hz : z ∈ N) :
    velX J (X ∘ τ) z = velX J X (τ z) := by
  have hA : IsUnit (jacM τ z).det := isUnit_iff_ne_zero.2 (hdet z hz)
  have hk := pr_key hA (gX_comp hΩ hX hN hτ hmaps hz)
  have hpot : potX J (X ∘ τ) =ᶠ[𝓝 z] potX J X ∘ τ := by
    filter_upwards [hN.mem_nhds hz] with w hw
    exact potX_comp hΩ hX hinj hlag hN hτ hmaps hdet hw
  have hpd : DifferentiableAt ℝ (potX J X) (τ z) :=
    ((vv_contDiffOn_potX (J := J) hΩ hX hinj).contDiffAt
      (hΩ.mem_nhds (hmaps z hz))).differentiableAt (by simp)
  have hfd : ∀ b, fderiv ℝ (potX J (X ∘ τ)) z (eb b) =
      ∑ b', jacM τ z b' b * fderiv ℝ (potX J X) (τ z) (eb b') := by
    intro b
    rw [hpot.fderiv_eq, fderiv_comp z hpd (ir_dτ hN hτ hz)]
    show fderiv ℝ (potX J X) (τ z) (fderiv ℝ τ z (eb b)) = _
    rw [ir_fderiv_tau]
    simp [map_sum, map_smul]
  unfold velX
  simp only [hfd, dX_comp hΩ hX hN hτ hmaps hz, map_sum, map_smul, ← hk, Finset.smul_sum,
    Finset.sum_smul, Finset.mul_sum, Finset.sum_mul, smul_smul]
  rw [pr_swap22]
  refine Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun b _ =>
    Finset.sum_congr rfl fun c _ => Finset.sum_congr rfl fun d _ => ?_
  congr 1
  ring

/-- info: 'CWFlow.Extension.potX_comp' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms potX_comp

/-- info: 'CWFlow.Extension.muC_comp' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms muC_comp

/-- info: 'CWFlow.Extension.velX_comp' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms velX_comp

end CWFlow.Extension
