import CWFlow.Extension.ImmSmooth
import CWFlow.Extension.PullZ

/-!
# Behaviour of the coordinate tensors under a change of coordinates

For a Lagrangian immersion `X` and a smooth map `τ` with invertible Jacobian, the coordinate tensors of
`X ∘ τ` at `z` are the pullbacks by `A = jacM τ z` of those of `X` at `τ z`:
`gX`, `cX`, `SX k` (tensoriality of the covariant derivatives of `C`), and `gamX` transforms with the
inhomogeneous term `D = ∂∂τ`. Consequently the energy densities transform as densities.
-/

open scoped ContDiff Matrix RealInnerProductSpace Topology

namespace CWFlow.Extension

variable {n : ℕ} {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V] {J : V →L[ℝ] V}
  {X : EuclideanSpace ℝ (Fin n) → V} {τ : EuclideanSpace ℝ (Fin n) → EuclideanSpace ℝ (Fin n)}
  {Ω N : Set (EuclideanSpace ℝ (Fin n))} {z : EuclideanSpace ℝ (Fin n)}

/-! ### Part 1: first and second derivatives, metric -/

theorem ir_one_add_one_le : (1 : WithTop ℕ∞) + 1 ≤ ∞ := by
  rw [one_add_one_eq_two]; exact WithTop.coe_le_coe.2 le_top

theorem ir_expand (w : EuclideanSpace ℝ (Fin n)) : w = ∑ j, w j • eb j := by
  ext i
  simp [eb, Pi.single_apply]

theorem ir_lin_expand {F : Type*} [NormedAddCommGroup F] [NormedSpace ℝ F]
    (L : EuclideanSpace ℝ (Fin n) →L[ℝ] F) (v : EuclideanSpace ℝ (Fin n)) :
    L v = ∑ i, v i • L (eb i) := by
  conv_lhs => rw [ir_expand v]
  simp [map_sum, map_smul]

theorem ir_fderiv_tau (τ : EuclideanSpace ℝ (Fin n) → EuclideanSpace ℝ (Fin n))
    (z : EuclideanSpace ℝ (Fin n)) (i : Fin n) :
    fderiv ℝ τ z (eb i) = ∑ j, jacM τ z j i • eb j := by
  conv_lhs => rw [ir_expand (fderiv ℝ τ z (eb i))]
  rfl

theorem ir_dτ (hN : IsOpen N) (hτ : ContDiffOn ℝ ∞ τ N) {w : EuclideanSpace ℝ (Fin n)}
    (hw : w ∈ N) : DifferentiableAt ℝ τ w :=
  (hτ.contDiffAt (hN.mem_nhds hw)).differentiableAt (by simp)

theorem ir_dX (hΩ : IsOpen Ω) (hX : ContDiffOn ℝ ∞ X Ω) {y : EuclideanSpace ℝ (Fin n)}
    (hy : y ∈ Ω) : DifferentiableAt ℝ X y :=
  ((hX y hy).contDiffAt (hΩ.mem_nhds hy)).differentiableAt (by simp)

theorem ir_fcomp_ev (hΩ : IsOpen Ω) (hX : ContDiffOn ℝ ∞ X Ω) (hN : IsOpen N)
    (hτ : ContDiffOn ℝ ∞ τ N) (hmaps : ∀ z ∈ N, τ z ∈ Ω) (hz : z ∈ N) :
    fderiv ℝ (X ∘ τ) =ᶠ[𝓝 z] fun w => (fderiv ℝ X (τ w)).comp (fderiv ℝ τ w) := by
  filter_upwards [hN.mem_nhds hz] with w hw
  exact fderiv_comp w (ir_dX hΩ hX (hmaps w hw)) (ir_dτ hN hτ hw)

theorem ir_fcomp (hΩ : IsOpen Ω) (hX : ContDiffOn ℝ ∞ X Ω) (hN : IsOpen N)
    (hτ : ContDiffOn ℝ ∞ τ N) (hmaps : ∀ z ∈ N, τ z ∈ Ω) (hz : z ∈ N)
    (v : EuclideanSpace ℝ (Fin n)) :
    fderiv ℝ (X ∘ τ) z v = fderiv ℝ X (τ z) (fderiv ℝ τ z v) := by
  rw [(ir_fcomp_ev hΩ hX hN hτ hmaps hz).self_of_nhds]
  rfl

theorem ir_fderiv_fderiv_comp (hΩ : IsOpen Ω) (hX : ContDiffOn ℝ ∞ X Ω) (hN : IsOpen N)
    (hτ : ContDiffOn ℝ ∞ τ N) (hmaps : ∀ z ∈ N, τ z ∈ Ω) (hz : z ∈ N)
    (a c : EuclideanSpace ℝ (Fin n)) :
    fderiv ℝ (fderiv ℝ (X ∘ τ)) z a c =
      fderiv ℝ (fderiv ℝ X) (τ z) (fderiv ℝ τ z a) (fderiv ℝ τ z c) +
        fderiv ℝ X (τ z) (fderiv ℝ (fderiv ℝ τ) z a c) := by
  have hτz := hτ.contDiffAt (hN.mem_nhds hz)
  have hτd := hτz.differentiableAt (by simp)
  have hXz := (hX (τ z) (hmaps z hz)).contDiffAt (hΩ.mem_nhds (hmaps z hz))
  have hDX := (hXz.fderiv_right (m := 1) ir_one_add_one_le).differentiableAt (by simp)
  have hDτ := (hτz.fderiv_right (m := 1) ir_one_add_one_le).differentiableAt (by simp)
  have hA := hDX.hasFDerivAt.comp z hτd.hasFDerivAt
  have h : HasFDerivAt (fun w => (fderiv ℝ X (τ w)).comp (fderiv ℝ τ w)) _ z :=
    hA.clm_comp hDτ.hasFDerivAt
  rw [(ir_fcomp_ev hΩ hX hN hτ hmaps hz).fderiv_eq, h.fderiv]
  simp only [ContinuousLinearMap.add_apply, ContinuousLinearMap.coe_comp', Function.comp_apply,
    ContinuousLinearMap.compL_apply, ContinuousLinearMap.flip_apply]
  abel

theorem ir_fderiv_apply_const {F : Type*} [NormedAddCommGroup F] [NormedSpace ℝ F]
    {f : EuclideanSpace ℝ (Fin n) → EuclideanSpace ℝ (Fin n) →L[ℝ] F}
    (hf : DifferentiableAt ℝ f z) (v u : EuclideanSpace ℝ (Fin n)) :
    fderiv ℝ (fun y => f y u) z v = fderiv ℝ f z v u := by
  have h := (ContinuousLinearMap.apply ℝ F u).hasFDerivAt.comp z hf.hasFDerivAt
  have h2 : fderiv ℝ (fun y => f y u) z = (ContinuousLinearMap.apply ℝ F u).comp (fderiv ℝ f z) :=
    h.fderiv
  rw [h2]
  rfl

/-- `D a b c' = ∂_a ∂_b τ^{c'}` as a component of the second derivative of `τ`. -/
theorem ir_D_apply (hN : IsOpen N) (hτ : ContDiffOn ℝ ∞ τ N) (hz : z ∈ N) (a b c' : Fin n) :
    fderiv ℝ (fun y => jacM τ y c' b) z (eb a) =
      fderiv ℝ (fderiv ℝ τ) z (eb a) (eb b) c' := by
  have hτz := hτ.contDiffAt (hN.mem_nhds hz)
  have hDτ := (hτz.fderiv_right (m := 1) ir_one_add_one_le).differentiableAt (by simp)
  have h1 : DifferentiableAt ℝ (fun y => fderiv ℝ τ y (eb b)) z :=
    (ContinuousLinearMap.apply ℝ (EuclideanSpace ℝ (Fin n)) (eb b)).differentiableAt.comp z hDτ
  have h2 := CWFlow.Interp.fderiv_piLp_apply h1 (eb a) c'
  rw [ir_fderiv_apply_const hDτ] at h2
  exact h2.symm

theorem dX_comp (hΩ : IsOpen Ω) (hX : ContDiffOn ℝ ∞ X Ω) (hN : IsOpen N)
    (hτ : ContDiffOn ℝ ∞ τ N) (hmaps : ∀ z ∈ N, τ z ∈ Ω) (hz : z ∈ N) (i : Fin n) :
    dX (X ∘ τ) z i = ∑ j, jacM τ z j i • dX X (τ z) j := by
  unfold dX
  rw [ir_fcomp hΩ hX hN hτ hmaps hz, ir_fderiv_tau]
  simp [map_sum, map_smul]

theorem ddX_comp (hΩ : IsOpen Ω) (hX : ContDiffOn ℝ ∞ X Ω) (hN : IsOpen N)
    (hτ : ContDiffOn ℝ ∞ τ N) (hmaps : ∀ z ∈ N, τ z ∈ Ω) (hz : z ∈ N) (a b : Fin n) :
    ddX (X ∘ τ) z a b =
      (∑ a', ∑ b', (jacM τ z a' a * jacM τ z b' b) • ddX X (τ z) a' b')
      + ∑ c', fderiv ℝ (fun y => jacM τ y c' b) z (eb a) • dX X (τ z) c' := by
  unfold ddX
  rw [ir_fderiv_fderiv_comp hΩ hX hN hτ hmaps hz, ir_fderiv_tau, ir_fderiv_tau]
  congr 1
  · simp only [map_sum, map_smul, ContinuousLinearMap.sum_apply, ContinuousLinearMap.smul_apply,
      Finset.smul_sum, smul_smul]
    rw [Finset.sum_comm]
    refine Finset.sum_congr rfl fun a' _ => Finset.sum_congr rfl fun b' _ => ?_
    rw [mul_comm]
  · conv_lhs => rw [ir_lin_expand (fderiv ℝ X (τ z)) (fderiv ℝ (fderiv ℝ τ) z (eb a) (eb b))]
    refine Finset.sum_congr rfl fun c' _ => ?_
    rw [ir_D_apply hN hτ hz]
    rfl

theorem gX_comp (hΩ : IsOpen Ω) (hX : ContDiffOn ℝ ∞ X Ω) (hN : IsOpen N)
    (hτ : ContDiffOn ℝ ∞ τ N) (hmaps : ∀ z ∈ N, τ z ∈ Ω) (hz : z ∈ N) :
    gX (X ∘ τ) z = (jacM τ z)ᵀ * gX X (τ z) * jacM τ z := by
  ext a b
  simp only [gX]
  rw [dX_comp hΩ hX hN hτ hmaps hz a, dX_comp hΩ hX hN hτ hmaps hz b]
  simp only [Matrix.mul_apply, Matrix.transpose_apply, sum_inner, inner_sum, inner_smul_left,
    inner_smul_right, conj_trivial, Finset.sum_mul, Finset.mul_sum, gX]
  refine Finset.sum_congr rfl fun j _ => Finset.sum_congr rfl fun k _ => ?_
  ring

theorem ir_inj_tau (hdet : (jacM τ z).det ≠ 0) : Function.Injective (fderiv ℝ τ z) := by
  rw [injective_iff_map_eq_zero]
  intro v hv
  have h : (jacM τ z) *ᵥ (fun i => v i) = 0 := by
    funext j
    have h1 : fderiv ℝ τ z v j = ∑ i, v i * jacM τ z j i := by
      rw [ir_lin_expand (fderiv ℝ τ z) v]
      simp [jacM]
    have h2 := congrArg (fun w : EuclideanSpace ℝ (Fin n) => w j) hv
    simp only [h1] at h2
    simp only [Matrix.mulVec, dotProduct, Pi.zero_apply]
    simpa [mul_comm] using h2
  have := Matrix.eq_zero_of_mulVec_eq_zero hdet h
  ext i
  exact congrFun this i

theorem injective_fderiv_comp (hΩ : IsOpen Ω) (hX : ContDiffOn ℝ ∞ X Ω)
    (hinj : ∀ x ∈ Ω, Function.Injective (fderiv ℝ X x)) (hN : IsOpen N)
    (hτ : ContDiffOn ℝ ∞ τ N) (hmaps : ∀ z ∈ N, τ z ∈ Ω)
    (hdet : ∀ z ∈ N, (jacM τ z).det ≠ 0) (hz : z ∈ N) :
    Function.Injective (fderiv ℝ (X ∘ τ) z) := by
  have : fderiv ℝ (X ∘ τ) z = (fderiv ℝ X (τ z)).comp (fderiv ℝ τ z) := by
    ext v
    exact ir_fcomp hΩ hX hN hτ hmaps hz v
  rw [this]
  exact (hinj _ (hmaps z hz)).comp (ir_inj_tau (hdet z hz))

theorem lagrangian_comp (hΩ : IsOpen Ω) (hX : ContDiffOn ℝ ∞ X Ω)
    (hlag : ∀ x ∈ Ω, ∀ u w, ⟪J (fderiv ℝ X x u), fderiv ℝ X x w⟫ = 0) (hN : IsOpen N)
    (hτ : ContDiffOn ℝ ∞ τ N) (hmaps : ∀ z ∈ N, τ z ∈ Ω) (hz : z ∈ N)
    (u w : EuclideanSpace ℝ (Fin n)) :
    ⟪J (fderiv ℝ (X ∘ τ) z u), fderiv ℝ (X ∘ τ) z w⟫ = 0 := by
  rw [ir_fcomp hΩ hX hN hτ hmaps hz, ir_fcomp hΩ hX hN hτ hmaps hz]
  exact hlag _ (hmaps z hz) _ _

theorem contDiffOn_comp_reparam (hX : ContDiffOn ℝ ∞ X Ω) (hτ : ContDiffOn ℝ ∞ τ N)
    (hmaps : ∀ z ∈ N, τ z ∈ Ω) : ContDiffOn ℝ ∞ (X ∘ τ) N :=
  hX.comp hτ hmaps

theorem volX_comp (hΩ : IsOpen Ω) (hX : ContDiffOn ℝ ∞ X Ω) (hN : IsOpen N)
    (hτ : ContDiffOn ℝ ∞ τ N) (hmaps : ∀ z ∈ N, τ z ∈ Ω) (hz : z ∈ N) :
    volX (X ∘ τ) z = |(jacM τ z).det| * volX X (τ z) := by
  unfold volX
  exact sqrt_det_conj (gX_comp hΩ hX hN hτ hmaps hz)

/-! ### Part 2: the form `C` and the Christoffel symbols -/

theorem ir_sum3 (g : (Fin 3 → Fin n) → ℝ) :
    ∑ K, g K = ∑ a, ∑ b, ∑ c, g ![a, b, c] := by
  rw [sum_cons' 2]
  refine Finset.sum_congr rfl fun a _ => ?_
  rw [sum_cons' 1]
  refine Finset.sum_congr rfl fun b _ => ?_
  rw [sum_cons' 0]
  refine Finset.sum_congr rfl fun c _ => ?_
  rw [Fintype.sum_unique]
  congr 1

theorem cX_comp (hΩ : IsOpen Ω) (hX : ContDiffOn ℝ ∞ X Ω)
    (hlag : ∀ x ∈ Ω, ∀ u w, ⟪J (fderiv ℝ X x u), fderiv ℝ X x w⟫ = 0) (hN : IsOpen N)
    (hτ : ContDiffOn ℝ ∞ τ N) (hmaps : ∀ z ∈ N, τ z ∈ Ω) (hz : z ∈ N) :
    cX J (X ∘ τ) z = pullZ (jacM τ z) (cX J X (τ z)) := by
  have hzero : ∀ c' j, ⟪dX X (τ z) c', J (dX X (τ z) j)⟫ = 0 := by
    intro c' j
    rw [real_inner_comm]
    exact hlag _ (hmaps z hz) _ _
  ext I
  rw [pullZ_apply, ir_sum3]
  show ⟪ddX (X ∘ τ) z (I 0) (I 1), J (dX (X ∘ τ) z (I 2))⟫ = _
  rw [ddX_comp hΩ hX hN hτ hmaps hz, dX_comp hΩ hX hN hτ hmaps hz]
  simp only [map_sum, map_smul, inner_add_left, sum_inner, inner_sum, inner_smul_left,
    inner_smul_right, conj_trivial, hzero, mul_zero, Finset.sum_const_zero, add_zero]
  simp only [Fin.prod_univ_three, cX]
  simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two, Matrix.head_cons,
    Matrix.tail_cons]
  simp only [Finset.mul_sum]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun a' _ => ?_
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun b' _ => Finset.sum_congr rfl fun c' _ => ?_
  ring

theorem ir_ginv_A {A g g' : Matrix (Fin n) (Fin n) ℝ} (hA : IsUnit A.det) (hg' : g' = Aᵀ * g * A) :
    g'⁻¹ * Aᵀ = A⁻¹ * g⁻¹ := by
  have hAt : IsUnit (Aᵀ).det := by rwa [Matrix.det_transpose]
  have h1 : g'⁻¹ = A⁻¹ * (g⁻¹ * (Aᵀ)⁻¹) := by
    rw [hg', Matrix.mul_inv_rev, Matrix.mul_inv_rev]
  rw [h1, Matrix.mul_assoc, Matrix.mul_assoc, Matrix.nonsing_inv_mul _ hAt, Matrix.mul_one]

theorem gamX_comp (hΩ : IsOpen Ω) (hX : ContDiffOn ℝ ∞ X Ω)
    (hinj : ∀ x ∈ Ω, Function.Injective (fderiv ℝ X x)) (hN : IsOpen N)
    (hτ : ContDiffOn ℝ ∞ τ N) (hmaps : ∀ z ∈ N, τ z ∈ Ω)
    (hdet : ∀ z ∈ N, (jacM τ z).det ≠ 0) (hz : z ∈ N) (a b c : Fin n) :
    gamX (X ∘ τ) z ![a, b, c] =
      ∑ c', (jacM τ z)⁻¹ c c' *
        ((∑ a', ∑ b', jacM τ z a' a * jacM τ z b' b * gamX X (τ z) ![a', b', c'])
          + fderiv ℝ (fun y => jacM τ y c' b) z (eb a)) := by
  have hAu : IsUnit (jacM τ z).det := isUnit_iff_ne_zero.2 (hdet z hz)
  have hg' := gX_comp hΩ hX hN hτ hmaps hz
  have hGi := gX_inv_mul (hinj _ (hmaps z hz))
  have hΓ : ∀ a' b' p, gamX X (τ z) ![a', b', p] =
      ∑ e, (gX X (τ z))⁻¹ p e * ⟪ddX X (τ z) a' b', dX X (τ z) e⟫ := fun _ _ _ => rfl
  -- the expansion of `w d' = ⟪∂_a∂_b (X∘τ), ∂_{d'} X⟫`
  have hw : ∀ d', ⟪ddX (X ∘ τ) z a b, dX X (τ z) d'⟫ =
      (∑ a', ∑ b', jacM τ z a' a * jacM τ z b' b * ⟪ddX X (τ z) a' b', dX X (τ z) d'⟫)
        + ∑ e, fderiv ℝ (fun y => jacM τ y e b) z (eb a) * gX X (τ z) e d' := by
    intro d'
    rw [ddX_comp hΩ hX hN hτ hmaps hz]
    simp only [inner_add_left, sum_inner, inner_smul_left, conj_trivial, gX]
  have hkey : ∀ p, ∑ d', (gX X (τ z))⁻¹ p d' * ⟪ddX (X ∘ τ) z a b, dX X (τ z) d'⟫ =
      (∑ a', ∑ b', jacM τ z a' a * jacM τ z b' b * gamX X (τ z) ![a', b', p])
        + fderiv ℝ (fun y => jacM τ y p b) z (eb a) := by
    intro p
    simp only [hw, mul_add, Finset.sum_add_distrib]
    congr 1
    · simp only [hΓ, Finset.mul_sum]
      rw [Finset.sum_comm]
      refine Finset.sum_congr rfl fun a' _ => ?_
      rw [Finset.sum_comm]
      refine Finset.sum_congr rfl fun b' _ => ?_
      refine Finset.sum_congr rfl fun d' _ => ?_
      ring
    · have h1 : ∀ e, ∑ d', (gX X (τ z))⁻¹ p d' * gX X (τ z) e d' = if p = e then 1 else 0 := by
        intro e
        have := congrFun (congrFun hGi p) e
        rw [Matrix.mul_apply, Matrix.one_apply] at this
        rw [← this]
        refine Finset.sum_congr rfl fun d' _ => ?_
        rw [gX_symm]
      have h2 : ∀ e, ∑ d', (gX X (τ z))⁻¹ p d' *
          (fderiv ℝ (fun y => jacM τ y e b) z (eb a) * gX X (τ z) e d') =
          fderiv ℝ (fun y => jacM τ y e b) z (eb a) * (if p = e then 1 else 0) := by
        intro e
        rw [← h1 e, Finset.mul_sum]
        refine Finset.sum_congr rfl fun d' _ => ?_
        ring
      simp only [Finset.mul_sum]
      rw [Finset.sum_comm]
      simp [h2]
  have hAB := ir_ginv_A hAu hg'
  have h2 : ∀ d', ∑ d, (gX (X ∘ τ) z)⁻¹ c d * jacM τ z d' d =
      ∑ p, (jacM τ z)⁻¹ c p * (gX X (τ z))⁻¹ p d' := by
    intro d'
    have := congrFun (congrFun hAB c) d'
    rw [Matrix.mul_apply, Matrix.mul_apply] at this
    simpa [Matrix.transpose_apply] using this
  have h1 : ∀ d, ⟪ddX (X ∘ τ) z a b, dX (X ∘ τ) z d⟫ =
      ∑ d', jacM τ z d' d * ⟪ddX (X ∘ τ) z a b, dX X (τ z) d'⟫ := by
    intro d
    rw [dX_comp hΩ hX hN hτ hmaps hz d]
    simp only [inner_sum, inner_smul_right]
  show ∑ d, (gX (X ∘ τ) z)⁻¹ c d * ⟪ddX (X ∘ τ) z a b, dX (X ∘ τ) z d⟫ = _
  calc ∑ d, (gX (X ∘ τ) z)⁻¹ c d * ⟪ddX (X ∘ τ) z a b, dX (X ∘ τ) z d⟫
      = ∑ d', (∑ d, (gX (X ∘ τ) z)⁻¹ c d * jacM τ z d' d) *
          ⟪ddX (X ∘ τ) z a b, dX X (τ z) d'⟫ := by
        simp only [h1, Finset.mul_sum]
        rw [Finset.sum_comm]
        refine Finset.sum_congr rfl fun d' _ => ?_
        rw [Finset.sum_mul]
        refine Finset.sum_congr rfl fun d _ => ?_
        ring
    _ = ∑ d', ∑ p, (jacM τ z)⁻¹ c p *
          ((gX X (τ z))⁻¹ p d' * ⟪ddX (X ∘ τ) z a b, dX X (τ z) d'⟫) := by
        refine Finset.sum_congr rfl fun d' _ => ?_
        rw [h2, Finset.sum_mul]
        refine Finset.sum_congr rfl fun p _ => ?_
        ring
    _ = ∑ p, (jacM τ z)⁻¹ c p *
          ∑ d', (gX X (τ z))⁻¹ p d' * ⟪ddX (X ∘ τ) z a b, dX X (τ z) d'⟫ := by
        rw [Finset.sum_comm]
        refine Finset.sum_congr rfl fun p _ => ?_
        rw [Finset.mul_sum]
    _ = _ := by
        refine Finset.sum_congr rfl fun p _ => ?_
        rw [hkey]

/-! ### Part 3: the covariant derivatives -/

theorem ir_jac_diff (hN : IsOpen N) (hτ : ContDiffOn ℝ ∞ τ N) (hz : z ∈ N) (j i : Fin n) :
    DifferentiableAt ℝ (fun y => jacM τ y j i) z := by
  have hτz := hτ.contDiffAt (hN.mem_nhds hz)
  have hDτ := (hτz.fderiv_right (m := 1) ir_one_add_one_le).differentiableAt (by simp)
  have h1 : DifferentiableAt ℝ (fun y => fderiv ℝ τ y (eb i)) z :=
    (ContinuousLinearMap.apply ℝ (EuclideanSpace ℝ (Fin n)) (eb i)).differentiableAt.comp z hDτ
  exact (PiLp.proj (𝕜 := ℝ) 2 (fun _ : Fin n => ℝ) j).differentiableAt.comp z h1

theorem ir_pullZ_sum {m : ℕ} (A : Matrix (Fin n) (Fin n) ℝ) (c : Fin n → ℝ) (v : Fin n → Z n m)
    (I : Fin m → Fin n) :
    pullZ A (∑ a', c a' • v a') I = ∑ a', c a' * pullZ A (v a') I := by
  have h : pullZ A (∑ a', c a' • v a') = ∑ a', c a' • pullZ A (v a') := by
    have h := map_sum (pullZL A) (fun a' => c a' • v a') Finset.univ
    simp only [map_smul, pullZL_apply] at h
    exact h
  rw [h]
  simp

theorem ir_pullZ_curry {m : ℕ} (A : Matrix (Fin n) (Fin n) ℝ)
    (L : EuclideanSpace ℝ (Fin n) →L[ℝ] Z n m) (a : Fin n) (I : Fin m → Fin n) :
    pullZ A (curryZ L) (Fin.cons a I) = ∑ a', A a' a * pullZ A (L (eb a')) I := by
  rw [pullZ_cons]
  refine Finset.sum_congr rfl fun a' _ => ?_
  congr 1

theorem SX_comp (hΩ : IsOpen Ω) (hX : ContDiffOn ℝ ∞ X Ω)
    (hinj : ∀ x ∈ Ω, Function.Injective (fderiv ℝ X x))
    (hlag : ∀ x ∈ Ω, ∀ u w, ⟪J (fderiv ℝ X x u), fderiv ℝ X x w⟫ = 0) (hN : IsOpen N)
    (hτ : ContDiffOn ℝ ∞ τ N) (hmaps : ∀ z ∈ N, τ z ∈ Ω)
    (hdet : ∀ z ∈ N, (jacM τ z).det ≠ 0) (k : ℕ) (hz : z ∈ N) :
    SX J (X ∘ τ) k z = pullZ (jacM τ z) (SX J X k (τ z)) := by
  induction k generalizing z with
  | zero => exact cX_comp hΩ hX hlag hN hτ hmaps hz
  | succ k ih =>
    have hX' := contDiffOn_comp_reparam hX hτ hmaps
    have hinj' : ∀ x ∈ N, Function.Injective (fderiv ℝ (X ∘ τ) x) := fun x hx =>
      injective_fderiv_comp hΩ hX hinj hN hτ hmaps hdet hx
    have hBA : (jacM τ z)⁻¹ * jacM τ z = 1 :=
      Matrix.nonsing_inv_mul _ (isUnit_iff_ne_zero.2 (hdet z hz))
    have hev : SX J (X ∘ τ) k =ᶠ[𝓝 z] fun y => pullZ (jacM τ y) (SX J X k (τ y)) := by
      filter_upwards [hN.mem_nhds hz] with y hy using ih hy
    have hdS : DifferentiableAt ℝ (SX J (X ∘ τ) k) z :=
      ((contDiffOn_SX (J := J) hN hX' (hinj' ) k).contDiffAt (hN.mem_nhds hz)).differentiableAt
        (by simp)
    have hF : DifferentiableAt ℝ (fun y => pullZ (jacM τ y) (SX J X k (τ y))) z :=
      hdS.congr_of_eventuallyEq hev.symm
    have hτd := ir_dτ hN hτ hz
    have hSd : DifferentiableAt ℝ (SX J X k) (τ z) :=
      ((contDiffOn_SX (J := J) hΩ hX hinj k).contDiffAt
        (hΩ.mem_nhds (hmaps z hz))).differentiableAt (by simp)
    have hT : DifferentiableAt ℝ (fun y => SX J X k (τ y)) z := hSd.comp z hτd
    suffices h : ∀ (a : Fin n) (I : Fin (k + 3) → Fin n),
        SX J (X ∘ τ) (k + 1) z (Fin.cons a I) =
          pullZ (jacM τ z) (SX J X (k + 1) (τ z)) (Fin.cons a I) by
      ext I
      have := h (I 0) (Fin.tail I)
      rwa [Fin.cons_self_tail] at this
    intro a I
    have hL : SX J (X ∘ τ) (k + 1) z (Fin.cons a I) =
        fderiv ℝ (fun y => pullZ (jacM τ y) (SX J X k (τ y)) I) z (eb a)
          - actG (gamX (X ∘ τ) z) (SX J (X ∘ τ) k z) (Fin.cons a I) := by
      rw [SX_succ, PiLp.sub_apply, curryZ_apply, hev.fderiv_eq]
      simp only [Fin.cons_zero, Fin.tail_cons]
      rw [CWFlow.Interp.fderiv_piLp_apply hF]
    have hR : pullZ (jacM τ z) (SX J X (k + 1) (τ z)) (Fin.cons a I) =
        pullZ (jacM τ z) (curryZ (fderiv ℝ (SX J X k) (τ z))) (Fin.cons a I)
          - pullZ (jacM τ z) (actG (gamX X (τ z)) (SX J X k (τ z))) (Fin.cons a I) := by
      rw [SX_succ, pullZ_sub]
      rfl
    have hP : pullZ (jacM τ z) (fderiv ℝ (fun y => SX J X k (τ y)) z (eb a)) I =
        pullZ (jacM τ z) (curryZ (fderiv ℝ (SX J X k) (τ z))) (Fin.cons a I) := by
      have hc : fderiv ℝ (fun y => SX J X k (τ y)) z =
          (fderiv ℝ (SX J X k) (τ z)).comp (fderiv ℝ τ z) := fderiv_comp z hSd hτd
      rw [ir_pullZ_curry, hc]
      simp only [ContinuousLinearMap.comp_apply]
      rw [ir_fderiv_tau]
      simp only [map_sum, map_smul]
      exact ir_pullZ_sum _ _ (fun a' => fderiv ℝ (SX J X k) (τ z) (eb a')) I
    have hF1 := fderiv_pullZ_apply (A := jacM τ) (T := fun y => SX J X k (τ y))
      (ir_jac_diff hN hτ hz) hT (eb a) I
    have hG1 := actG_pullZ hBA (gamX X (τ z)) (gamX (X ∘ τ) z)
      (fun a b c' => fderiv ℝ (fun y => jacM τ y c' b) z (eb a))
      (gamX_comp hΩ hX hinj hN hτ hmaps hdet hz) (SX J X k (τ z)) a I
    rw [hL, hR, hF1, ih hz, hG1, hP]
    ring

/-! ### Part 4: norms and densities -/

theorem gNormSq_SX_comp (hΩ : IsOpen Ω) (hX : ContDiffOn ℝ ∞ X Ω)
    (hinj : ∀ x ∈ Ω, Function.Injective (fderiv ℝ X x))
    (hlag : ∀ x ∈ Ω, ∀ u w, ⟪J (fderiv ℝ X x u), fderiv ℝ X x w⟫ = 0) (hN : IsOpen N)
    (hτ : ContDiffOn ℝ ∞ τ N) (hmaps : ∀ z ∈ N, τ z ∈ Ω)
    (hdet : ∀ z ∈ N, (jacM τ z).det ≠ 0) (k : ℕ) (hz : z ∈ N) :
    gNormSq (ginvZ (X ∘ τ) z) (SX J (X ∘ τ) k z) =
      gNormSq (ginvZ X (τ z)) (SX J X k (τ z)) := by
  rw [SX_comp hΩ hX hinj hlag hN hτ hmaps hdet k hz]
  exact gNormSq_pullZ_inv (isUnit_iff_ne_zero.2 (hdet z hz))
    (isUnit_iff_ne_zero.2 (gX_det_pos (hinj _ (hmaps z hz))).ne')
    (gX_comp hΩ hX hN hτ hmaps hz) _

theorem eDens_comp (hΩ : IsOpen Ω) (hX : ContDiffOn ℝ ∞ X Ω)
    (hinj : ∀ x ∈ Ω, Function.Injective (fderiv ℝ X x))
    (hlag : ∀ x ∈ Ω, ∀ u w, ⟪J (fderiv ℝ X x u), fderiv ℝ X x w⟫ = 0) (hN : IsOpen N)
    (hτ : ContDiffOn ℝ ∞ τ N) (hmaps : ∀ z ∈ N, τ z ∈ Ω)
    (hdet : ∀ z ∈ N, (jacM τ z).det ≠ 0) (k : ℕ) (hz : z ∈ N) :
    eDens J (X ∘ τ) k z = eDens J X k (τ z) * |(jacM τ z).det| := by
  unfold eDens
  rw [gNormSq_SX_comp hΩ hX hinj hlag hN hτ hmaps hdet k hz, volX_comp hΩ hX hN hτ hmaps hz]
  ring

theorem aNormSq_comp (hΩ : IsOpen Ω) (hX : ContDiffOn ℝ ∞ X Ω)
    (hinj : ∀ x ∈ Ω, Function.Injective (fderiv ℝ X x))
    (hlag : ∀ x ∈ Ω, ∀ u w, ⟪J (fderiv ℝ X x u), fderiv ℝ X x w⟫ = 0) (hN : IsOpen N)
    (hτ : ContDiffOn ℝ ∞ τ N) (hmaps : ∀ z ∈ N, τ z ∈ Ω)
    (hdet : ∀ z ∈ N, (jacM τ z).det ≠ 0) (hz : z ∈ N) :
    aNormSq J (X ∘ τ) z = aNormSq J X (τ z) := by
  unfold aNormSq
  rw [cX_comp hΩ hX hlag hN hτ hmaps hz]
  exact gNormSq_pullZ_inv (isUnit_iff_ne_zero.2 (hdet z hz))
    (isUnit_iff_ne_zero.2 (gX_det_pos (hinj _ (hmaps z hz))).ne')
    (gX_comp hΩ hX hN hτ hmaps hz) _

/-- info: 'CWFlow.Extension.gX_comp' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms gX_comp

/-- info: 'CWFlow.Extension.gamX_comp' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms gamX_comp

/-- info: 'CWFlow.Extension.SX_comp' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms SX_comp

/-- info: 'CWFlow.Extension.eDens_comp' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms eDens_comp

/-- info: 'CWFlow.Extension.aNormSq_comp' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms aNormSq_comp

end CWFlow.Extension
