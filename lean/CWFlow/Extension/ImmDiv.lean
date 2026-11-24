import CWFlow.Extension.ImmSmooth

/-!
# Coordinate identities for immersions: the potential `div(JH)` as a divergence

Basis-free-free (light) identities for `X : ℝⁿ → V`, split off from `ImmPot` so that they can be
used without its heavy imports: Gram matrix helpers (`ip_gi_g`, `ip_isUnit_gX`, ...), the
differentiability helpers `ip_diff_*`, the derivative of the inverse metric (`fderiv_gX_inv`), the
covariant derivative of `C` in coordinates (`SX_one_apply`), the algebra helpers (`ip_mat_alg`,
`ip_alg_i`, `ip_alg_dc`) and `potX_eq_div` (`potX` is `g^{ac}(∂_a μ_c − Γ^e_{ac} μ_e)`).
-/

open scoped ContDiff Matrix InnerProductSpace RealInnerProductSpace Topology

open CWFlow

namespace CWFlow.Extension

variable {n : ℕ} {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V] {J : V →L[ℝ] V}
  {X : EuclideanSpace ℝ (Fin n) → V} {Ω : Set (EuclideanSpace ℝ (Fin n))}
  {x : EuclideanSpace ℝ (Fin n)}

theorem ip_sum_eb (v : EuclideanSpace ℝ (Fin n)) : ∑ a, v a • eb a = v := by
  ext i
  simp [Pi.single_apply]

theorem ip_ext_eb {f g : EuclideanSpace ℝ (Fin n) →L[ℝ] ℝ} (h : ∀ c, f (eb c) = g (eb c)) :
    f = g := by
  refine ContinuousLinearMap.ext fun u => ?_
  have h1 : ∀ k : EuclideanSpace ℝ (Fin n) →L[ℝ] ℝ, k u = ∑ a, u a * k (eb a) := fun k => by
    conv_lhs => rw [← ip_sum_eb u]
    simp [map_sum, map_smul]
  rw [h1 f, h1 g]
  exact Finset.sum_congr rfl fun a _ => by rw [h a]

theorem ip_L_eq_sum (w : EuclideanSpace ℝ (Fin n)) :
    fderiv ℝ X x w = ∑ a, w a • dX X x a := by
  conv_lhs => rw [← ip_sum_eb w]
  simp [map_sum, map_smul, dX]

theorem ip_gX_symm (X : EuclideanSpace ℝ (Fin n) → V) (x : EuclideanSpace ℝ (Fin n)) :
    (gX X x)ᵀ = gX X x := by
  ext a b
  simp [gX, real_inner_comm]

theorem ip_ginv_symm (a b : Fin n) : (gX X x)⁻¹ a b = (gX X x)⁻¹ b a := by
  have h : ((gX X x)⁻¹)ᵀ = (gX X x)⁻¹ := by
    rw [Matrix.transpose_nonsing_inv, ip_gX_symm]
  exact (congrFun (congrFun h a) b).symm

theorem ip_inner_L (w w' : EuclideanSpace ℝ (Fin n)) :
    ⟪fderiv ℝ X x w, fderiv ℝ X x w'⟫ =
      ∑ a, w a * (Matrix.mulVec (gX X x) (fun b => w' b)) a := by
  rw [ip_L_eq_sum w, ip_L_eq_sum w', sum_inner]
  refine Finset.sum_congr rfl fun a _ => ?_
  rw [inner_smul_left, inner_sum]
  simp only [inner_smul_right, Matrix.mulVec, dotProduct, gX, conj_trivial]
  congr 1
  exact Finset.sum_congr rfl fun b _ => by ring

theorem ip_isUnit_gX (hinj : Function.Injective (fderiv ℝ X x)) : IsUnit (gX X x).det := by
  rw [← Matrix.isUnit_iff_isUnit_det, ← Matrix.mulVec_injective_iff_isUnit]
  intro v w hvw
  have h0 : Matrix.mulVec (gX X x) (v - w) = 0 := by
    rw [Matrix.mulVec_sub, hvw, sub_self]
  set ww : EuclideanSpace ℝ (Fin n) := WithLp.toLp 2 (v - w) with hww
  have h1 : ⟪fderiv ℝ X x ww, fderiv ℝ X x ww⟫ = 0 := by
    rw [ip_inner_L]
    have : (fun b => ww b) = v - w := rfl
    rw [this, h0]
    simp
  rw [inner_self_eq_zero] at h1
  have h2 : ww = 0 := hinj (by rw [h1, map_zero])
  have h3 : v - w = 0 := by
    funext i
    have := congrArg (fun u : EuclideanSpace ℝ (Fin n) => u i) h2
    simpa [hww] using this
  exact sub_eq_zero.mp h3

theorem ip_gi_g (hinj : Function.Injective (fderiv ℝ X x)) (b c : Fin n) :
    ∑ a, (gX X x)⁻¹ b a * gX X x a c = if b = c then 1 else 0 := by
  have := congrFun (congrFun (Matrix.nonsing_inv_mul _ (ip_isUnit_gX hinj)) b) c
  rw [Matrix.mul_apply, Matrix.one_apply] at this
  exact this

theorem ip_g_gi (hinj : Function.Injective (fderiv ℝ X x)) (b c : Fin n) :
    ∑ a, gX X x b a * (gX X x)⁻¹ a c = if b = c then 1 else 0 := by
  have := congrFun (congrFun (Matrix.mul_nonsing_inv _ (ip_isUnit_gX hinj)) b) c
  rw [Matrix.mul_apply, Matrix.one_apply] at this
  exact this

theorem ip_inner_J_left (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    (u v : V) : ⟪J u, v⟫ = -⟪u, J v⟫ := by
  calc ⟪J u, v⟫ = ⟪J (J u), J v⟫ := (hJi _ _).symm
    _ = ⟪-u, J v⟫ := by rw [hJJ]
    _ = -⟪u, J v⟫ := inner_neg_left _ _

theorem ip_inner_J_self (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    (v : V) : ⟪J v, v⟫ = 0 := by
  have h := ip_inner_J_left hJi hJJ v v
  rw [real_inner_comm (J v) v] at h
  linarith

theorem ip_inner_L_J_L (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    (hlag : ∀ u w, ⟪J (fderiv ℝ X x u), fderiv ℝ X x w⟫ = 0)
    (t u : EuclideanSpace ℝ (Fin n)) : ⟪fderiv ℝ X x t, J (fderiv ℝ X x u)⟫ = 0 := by
  have h := ip_inner_J_left hJi hJJ (fderiv ℝ X x t) (fderiv ℝ X x u)
  rw [hlag t u] at h
  linarith

theorem ip_contDiffAt (hΩ : IsOpen Ω) (hX : ContDiffOn ℝ ∞ X Ω) (hx : x ∈ Ω) :
    ContDiffAt ℝ ∞ X x :=
  hX.contDiffAt (hΩ.mem_nhds hx)

theorem ip_diff_fderiv (hX : ContDiffAt ℝ ∞ X x) : DifferentiableAt ℝ (fderiv ℝ X) x :=
  (hX.fderiv_right (m := ∞) le_rfl).differentiableAt (by simp)

theorem ip_diff_fderiv2 (hX : ContDiffAt ℝ ∞ X x) :
    DifferentiableAt ℝ (fderiv ℝ (fderiv ℝ X)) x :=
  ((hX.fderiv_right (m := ∞) le_rfl).fderiv_right (m := ∞) le_rfl).differentiableAt (by simp)

theorem ip_hasFDerivAt_dX (hX : ContDiffAt ℝ ∞ X x) (c : Fin n) :
    HasFDerivAt (fun y => dX X y c)
      ((ContinuousLinearMap.apply ℝ V (eb c)).comp (fderiv ℝ (fderiv ℝ X) x)) x :=
  (ContinuousLinearMap.apply ℝ V (eb c)).hasFDerivAt.comp x (ip_diff_fderiv hX).hasFDerivAt

theorem ip_diff_dX (hX : ContDiffAt ℝ ∞ X x) (c : Fin n) :
    DifferentiableAt ℝ (fun y => dX X y c) x :=
  (ip_hasFDerivAt_dX hX c).differentiableAt

theorem ip_fderiv_dX (hX : ContDiffAt ℝ ∞ X x) (a c : Fin n) :
    fderiv ℝ (fun y => dX X y c) x (eb a) = ddX X x a c := by
  rw [(ip_hasFDerivAt_dX hX c).fderiv]
  rfl

theorem ip_diff_ddX (hX : ContDiffAt ℝ ∞ X x) (a b : Fin n) :
    DifferentiableAt ℝ (fun y => ddX X y a b) x := by
  have h := ((ip_diff_fderiv2 hX).clm_apply (differentiableAt_const (eb a))).clm_apply
    (differentiableAt_const (eb b))
  exact h

theorem ip_diff_cX (hX : ContDiffAt ℝ ∞ X x) (I : Fin 3 → Fin n) :
    DifferentiableAt ℝ (fun y => cX J X y I) x :=
  (ip_diff_ddX hX (I 0) (I 1)).inner ℝ (J.differentiableAt.comp x (ip_diff_dX hX (I 2)))

theorem ip_diff_cX_vec (hX : ContDiffAt ℝ ∞ X x) : DifferentiableAt ℝ (cX J X) x :=
  (differentiableAt_piLp 2).2 fun I => ip_diff_cX hX I

theorem ip_gamX_eq (d p e : Fin n) :
    gamX X x ![d, p, e] = ∑ f, (gX X x)⁻¹ e f * ⟪ddX X x d p, dX X x f⟫ := rfl

theorem ip_diffAt_det {M : EuclideanSpace ℝ (Fin n) → Matrix (Fin n) (Fin n) ℝ}
    (h : ∀ i j, DifferentiableAt ℝ (fun y => M y i j) x) :
    DifferentiableAt ℝ (fun y => (M y).det) x := by
  simp_rw [Matrix.det_apply']
  refine DifferentiableAt.fun_sum fun σ _ => ?_
  exact (differentiableAt_const _).mul
    (HasFDerivAt.finsetProd fun i _ => (h _ _).hasFDerivAt).differentiableAt

theorem ip_diffAt_gX (hX : ContDiffAt ℝ ∞ X x) (p q : Fin n) :
    DifferentiableAt ℝ (fun y => gX X y p q) x :=
  (ip_diff_dX hX p).inner ℝ (ip_diff_dX hX q)

theorem ip_diffAt_ginv (hΩ : IsOpen Ω) (hX : ContDiffOn ℝ ∞ X Ω)
    (hinj : ∀ x ∈ Ω, Function.Injective (fderiv ℝ X x)) (hx : x ∈ Ω) (a b : Fin n) :
    DifferentiableAt ℝ (fun y => (gX X y)⁻¹ a b) x := by
  have hXa := ip_contDiffAt hΩ hX hx
  have hG := fun p q => ip_diffAt_gX hXa p q
  have hd : DifferentiableAt ℝ (fun y => (gX X y).det) x := ip_diffAt_det hG
  have hne : (gX X x).det ≠ 0 := (ip_isUnit_gX (hinj x hx)).ne_zero
  have hadj : DifferentiableAt ℝ (fun y => (gX X y).adjugate a b) x := by
    simp_rw [Matrix.adjugate_apply]
    refine ip_diffAt_det fun p q => ?_
    simp only [Matrix.updateRow_apply]
    by_cases hp : p = b
    · simp [hp]
    · simpa [hp] using hG p q
  have e : (fun y => (gX X y)⁻¹ a b) = fun y => ((gX X y).det)⁻¹ * (gX X y).adjugate a b := by
    funext y
    rw [Matrix.inv_def, Ring.inverse_eq_inv']
    rfl
  rw [e]
  exact (hd.inv hne).mul hadj

theorem ip_fderiv_gX (hX : ContDiffAt ℝ ∞ X x) (p q d : Fin n) :
    fderiv ℝ (fun y => gX X y p q) x (eb d) =
      ⟪ddX X x d p, dX X x q⟫ + ⟪dX X x p, ddX X x d q⟫ := by
  have h := fderiv_inner_apply ℝ (ip_diff_dX hX p) (ip_diff_dX hX q) (eb d)
  rw [ip_fderiv_dX hX, ip_fderiv_dX hX] at h
  show fderiv ℝ (fun t => ⟪dX X t p, dX X t q⟫) x (eb d) = _
  rw [h, add_comm]

theorem ip_gam_contract (hinj : Function.Injective (fderiv ℝ X x)) (d p q : Fin n) :
    ∑ e, gamX X x ![d, p, e] * gX X x e q = ⟪ddX X x d p, dX X x q⟫ := by
  simp only [ip_gamX_eq, Finset.sum_mul]
  rw [Finset.sum_comm]
  have h : ∀ f, ∑ e, (gX X x)⁻¹ e f * ⟪ddX X x d p, dX X x f⟫ * gX X x e q =
      ⟪ddX X x d p, dX X x f⟫ * (if f = q then 1 else 0) := by
    intro f
    rw [← ip_gi_g hinj f q, Finset.mul_sum]
    exact Finset.sum_congr rfl fun e _ => by rw [ip_ginv_symm e f]; ring
  simp only [h]
  simp

theorem ip_ginv_deriv_eq (hΩ : IsOpen Ω) (hX : ContDiffOn ℝ ∞ X Ω)
    (hinj : ∀ x ∈ Ω, Function.Injective (fderiv ℝ X x)) (hx : x ∈ Ω) (a q d : Fin n) :
    ∑ p, (fderiv ℝ (fun y => (gX X y)⁻¹ a p) x (eb d) * gX X x p q +
      (gX X x)⁻¹ a p * fderiv ℝ (fun y => gX X y p q) x (eb d)) = 0 := by
  have hXa := ip_contDiffAt hΩ hX hx
  have hev : (fun y => ∑ p, (gX X y)⁻¹ a p * gX X y p q) =ᶠ[𝓝 x]
      fun _ => if a = q then (1 : ℝ) else 0 := by
    filter_upwards [hΩ.mem_nhds hx] with y hy
    exact ip_gi_g (X := X) (x := y) (hinj y hy) a q
  have hd : HasFDerivAt (fun y => ∑ p, (gX X y)⁻¹ a p * gX X y p q)
      (∑ p, ((gX X x)⁻¹ a p • fderiv ℝ (fun y => gX X y p q) x +
        gX X x p q • fderiv ℝ (fun y => (gX X y)⁻¹ a p) x)) x :=
    HasFDerivAt.fun_sum fun p _ =>
      ((ip_diffAt_ginv hΩ hX hinj hx a p).hasFDerivAt).fun_mul (ip_diffAt_gX hXa p q).hasFDerivAt
  have h0 : fderiv ℝ (fun y => ∑ p, (gX X y)⁻¹ a p * gX X y p q) x = 0 := by
    rw [hev.fderiv_eq]
    simp
  have h1 := congrArg (fun L : EuclideanSpace ℝ (Fin n) →L[ℝ] ℝ => L (eb d))
    (hd.fderiv.symm.trans h0)
  simp only [ContinuousLinearMap.sum_apply, ContinuousLinearMap.add_apply,
    ContinuousLinearMap.smul_apply, smul_eq_mul, ContinuousLinearMap.zero_apply] at h1
  refine (Finset.sum_congr rfl fun p _ => ?_).trans h1
  ring

/-- Matrix form of `∂_d g_{pq}`'s inverse: `U = −(g⁻¹ Γ + Γᵀ g⁻¹)`. -/
theorem ip_mat_alg (G Gi U K Γd : Matrix (Fin n) (Fin n) ℝ) (hGs : Gᵀ = G) (hGGi : G * Gi = 1)
    (hGiG : Gi * G = 1) (h1 : U * G + Gi * K = 0) (h2 : K = Γd * G + (Γd * G)ᵀ) :
    U = -(Gi * Γd + Γdᵀ * Gi) := by
  have hU : U = -(Gi * K) * Gi := by
    calc U = U * (G * Gi) := by rw [hGGi, Matrix.mul_one]
      _ = (U * G) * Gi := by rw [Matrix.mul_assoc]
      _ = -(Gi * K) * Gi := by rw [eq_neg_of_add_eq_zero_left h1]
  have e1 : Gi * (Γd * G) * Gi = Gi * Γd := by
    calc Gi * (Γd * G) * Gi = Gi * Γd * (G * Gi) := by simp only [Matrix.mul_assoc]
      _ = Gi * Γd := by rw [hGGi, Matrix.mul_one]
  have e2 : Gi * (G * Γdᵀ) * Gi = Γdᵀ * Gi := by
    calc Gi * (G * Γdᵀ) * Gi = (Gi * G) * Γdᵀ * Gi := by simp only [Matrix.mul_assoc]
      _ = Γdᵀ * Gi := by rw [hGiG, Matrix.one_mul]
  rw [hU, h2, Matrix.neg_mul, Matrix.mul_add, Matrix.add_mul, Matrix.transpose_mul, hGs, e1, e2]

/-- The matrix `(∂_d g^{ab})`. -/
noncomputable def ip_dGinv (X : EuclideanSpace ℝ (Fin n) → V) (x : EuclideanSpace ℝ (Fin n))
    (d : Fin n) : Matrix (Fin n) (Fin n) ℝ :=
  fun a b => fderiv ℝ (fun y => (gX X y)⁻¹ a b) x (eb d)

/-- The matrix `(∂_d g_{pq})`. -/
noncomputable def ip_dG (X : EuclideanSpace ℝ (Fin n) → V) (x : EuclideanSpace ℝ (Fin n))
    (d : Fin n) : Matrix (Fin n) (Fin n) ℝ :=
  fun p q => fderiv ℝ (fun y => gX X y p q) x (eb d)

/-- The matrix `(Γ^e_{dp})_{p e}`. -/
noncomputable def ip_gamM (X : EuclideanSpace ℝ (Fin n) → V) (x : EuclideanSpace ℝ (Fin n))
    (d : Fin n) : Matrix (Fin n) (Fin n) ℝ :=
  fun p e => gamX X x ![d, p, e]

theorem fderiv_gX_inv (hΩ : IsOpen Ω) (hX : ContDiffOn ℝ ∞ X Ω)
    (hinj : ∀ x ∈ Ω, Function.Injective (fderiv ℝ X x)) (hx : x ∈ Ω) (a b d : Fin n) :
    fderiv ℝ (fun y => (gX X y)⁻¹ a b) x (eb d) =
      -∑ e, ((gX X x)⁻¹ a e * gamX X x ![d, e, b] + (gX X x)⁻¹ b e * gamX X x ![d, e, a]) := by
  have hXa := ip_contDiffAt hΩ hX hx
  have hL := hinj x hx
  have hu := ip_isUnit_gX hL
  have hGGi : gX X x * (gX X x)⁻¹ = 1 := Matrix.mul_nonsing_inv _ hu
  have hGiG : (gX X x)⁻¹ * gX X x = 1 := Matrix.nonsing_inv_mul _ hu
  have h1 : ip_dGinv X x d * gX X x + (gX X x)⁻¹ * ip_dG X x d = 0 := by
    ext a q
    have h := ip_ginv_deriv_eq hΩ hX hinj hx a q d
    simp only [Matrix.add_apply, Matrix.mul_apply, Matrix.zero_apply, ip_dGinv, ip_dG]
    rw [← Finset.sum_add_distrib]
    exact h
  have h2 : ip_dG X x d = ip_gamM X x d * gX X x + (ip_gamM X x d * gX X x)ᵀ := by
    ext p q
    simp only [Matrix.add_apply, Matrix.transpose_apply, Matrix.mul_apply, ip_dG, ip_gamM]
    rw [ip_fderiv_gX hXa p q d, ip_gam_contract hL d p q, ip_gam_contract hL d q p]
    congr 1
    exact real_inner_comm _ _
  have hmat := ip_mat_alg (gX X x) (gX X x)⁻¹ (ip_dGinv X x d) (ip_dG X x d) (ip_gamM X x d)
    (ip_gX_symm X x) hGGi hGiG h1 h2
  have h := congrFun (congrFun hmat a) b
  simp only [Matrix.neg_apply, Matrix.add_apply, Matrix.mul_apply, Matrix.transpose_apply,
    ip_gamM] at h
  show ip_dGinv X x d a b = _
  rw [h, neg_inj, ← Finset.sum_add_distrib]
  exact Finset.sum_congr rfl fun e _ => by rw [ip_ginv_symm e b]; ring

theorem ip_upd0 (a b c e : Fin n) : Function.update ![a, b, c] 0 e = ![e, b, c] := by
  funext i
  fin_cases i <;> simp

theorem ip_upd1 (a b c e : Fin n) : Function.update ![a, b, c] 1 e = ![a, e, c] := by
  funext i
  fin_cases i <;> simp

theorem ip_upd2 (a b c e : Fin n) : Function.update ![a, b, c] 2 e = ![a, b, e] := by
  funext i
  fin_cases i <;> simp

theorem ip_tail4 (d a b c : Fin n) : Fin.tail ![d, a, b, c] = ![a, b, c] := by
  funext i
  fin_cases i <;> rfl

theorem SX_one_apply (hΩ : IsOpen Ω) (hX : ContDiffOn ℝ ∞ X Ω)
    (hinj : ∀ x ∈ Ω, Function.Injective (fderiv ℝ X x)) (hx : x ∈ Ω) (a b c d : Fin n) :
    SX J X 1 x ![d, a, b, c] =
      fderiv ℝ (fun y => cX J X y ![a, b, c]) x (eb d)
      - ∑ e, (gamX X x ![d, a, e] * cX J X x ![e, b, c] + gamX X x ![d, b, e] * cX J X x ![a, e, c]
              + gamX X x ![d, c, e] * cX J X x ![a, b, e]) := by
  have hXa := ip_contDiffAt hΩ hX hx
  have h : SX J X 1 x = curryZ (fderiv ℝ (SX J X 0) x) - actG (gamX X x) (SX J X 0 x) :=
    SX_succ J X 0 x
  have e1 : curryZ (fderiv ℝ (cX J X) x) ![d, a, b, c] =
      fderiv ℝ (fun y => cX J X y ![a, b, c]) x (eb d) := by
    rw [curryZ_apply, ip_tail4]
    exact CWFlow.Interp.fderiv_piLp_apply (ip_diff_cX_vec hXa) (eb d) ![a, b, c]
  have e2 : actG (gamX X x) (cX J X x) ![d, a, b, c] =
      ∑ e, (gamX X x ![d, a, e] * cX J X x ![e, b, c] + gamX X x ![d, b, e] * cX J X x ![a, e, c]
              + gamX X x ![d, c, e] * cX J X x ![a, b, e]) := by
    rw [actG_apply, Fin.sum_univ_three]
    simp only [ip_tail4, ip_upd0, ip_upd1, ip_upd2, Finset.sum_add_distrib]
    rfl
  rw [h, SX_zero, PiLp.sub_apply, e1, e2]

theorem ip_ginvZ_apply (c d : Fin n) : ginvZ X x ![c, d] = (gX X x)⁻¹ c d := rfl

theorem ip_sum4 (f : Fin n → Fin n → Fin n → Fin n → ℝ) :
    ∑ a, ∑ b, ∑ c, ∑ d, f a b c d = ∑ d, ∑ c, ∑ a, ∑ b, f a b c d := by
  calc ∑ a, ∑ b, ∑ c, ∑ d, f a b c d = ∑ a, ∑ b, ∑ d, ∑ c, f a b c d :=
        Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun b _ => Finset.sum_comm
    _ = ∑ a, ∑ d, ∑ b, ∑ c, f a b c d := Finset.sum_congr rfl fun a _ => Finset.sum_comm
    _ = ∑ d, ∑ a, ∑ b, ∑ c, f a b c d := Finset.sum_comm
    _ = ∑ d, ∑ a, ∑ c, ∑ b, f a b c d :=
        Finset.sum_congr rfl fun d _ => Finset.sum_congr rfl fun a _ => Finset.sum_comm
    _ = ∑ d, ∑ c, ∑ a, ∑ b, f a b c d := Finset.sum_congr rfl fun d _ => Finset.sum_comm

theorem ip_alg_i (gi : Fin n → Fin n → ℝ) (hs : ∀ a b, gi a b = gi b a)
    (G1 : Fin n → Fin n → ℝ) (C2 : Fin n → Fin n → ℝ) :
    ∑ a, ∑ b, (-∑ e, (gi a e * G1 e b + gi b e * G1 e a)) * C2 a b =
      -(∑ a, ∑ b, gi a b * ∑ e, G1 a e * C2 e b) -
        ∑ a, ∑ b, gi a b * ∑ e, G1 b e * C2 a e := by
  have hA : ∑ a, ∑ b, ∑ e, gi a e * (G1 e b * C2 a b) =
      ∑ a, ∑ b, ∑ e, gi a b * (G1 b e * C2 a e) :=
    Finset.sum_congr rfl fun a _ => Finset.sum_comm
  have hB : ∑ a, ∑ b, ∑ e, gi b e * (G1 e a * C2 a b) =
      ∑ a, ∑ b, ∑ e, gi a b * (G1 a e * C2 e b) := by
    calc ∑ a, ∑ b, ∑ e, gi b e * (G1 e a * C2 a b)
        = ∑ b, ∑ a, ∑ e, gi b e * (G1 e a * C2 a b) := Finset.sum_comm
      _ = ∑ b, ∑ e, ∑ a, gi b e * (G1 e a * C2 a b) :=
        Finset.sum_congr rfl fun b _ => Finset.sum_comm
      _ = ∑ b, ∑ a, ∑ e, gi a b * (G1 a e * C2 e b) :=
        Finset.sum_congr rfl fun b _ => Finset.sum_congr rfl fun e _ =>
          Finset.sum_congr rfl fun a _ => by rw [hs b e]
      _ = ∑ a, ∑ b, ∑ e, gi a b * (G1 a e * C2 e b) := Finset.sum_comm
  have hl : ∑ a, ∑ b, (-∑ e, (gi a e * G1 e b + gi b e * G1 e a)) * C2 a b =
      -(∑ a, ∑ b, ∑ e, gi a e * (G1 e b * C2 a b) +
        ∑ a, ∑ b, ∑ e, gi b e * (G1 e a * C2 a b)) := by
    simp only [neg_mul, Finset.sum_mul, add_mul, Finset.sum_add_distrib, Finset.sum_neg_distrib,
      mul_assoc]
  have hr : -(∑ a, ∑ b, gi a b * ∑ e, G1 a e * C2 e b) -
        ∑ a, ∑ b, gi a b * ∑ e, G1 b e * C2 a e =
      -(∑ a, ∑ b, ∑ e, gi a b * (G1 a e * C2 e b)) -
        ∑ a, ∑ b, ∑ e, gi a b * (G1 b e * C2 a e) := by
    simp only [Finset.mul_sum]
  rw [hl, hr, hA, hB]
  ring

theorem ip_alg_dc (gi : Fin n → Fin n → ℝ) (hs : ∀ a b, gi a b = gi b a)
    (G1 : Fin n → Fin n → ℝ) (Gc : Fin n → ℝ) (C3 : Fin n → Fin n → Fin n → ℝ)
    (DC Dg : Fin n → Fin n → ℝ) (c : Fin n)
    (hD : ∀ a b, Dg a b = -∑ e, (gi a e * G1 e b + gi b e * G1 e a)) (g0 : ℝ) :
    g0 * (-(∑ a, ∑ b, (Dg a b * C3 a b c + gi a b * DC a b)) -
        ∑ e, Gc e * (-∑ a, ∑ b, gi a b * C3 a b e)) =
      -∑ a, ∑ b, g0 * gi a b *
        (DC a b - ∑ e, (G1 a e * C3 e b c + G1 b e * C3 a e c + Gc e * C3 a b e)) := by
  have hP : ∑ a, ∑ b, Dg a b * C3 a b c =
      -(∑ a, ∑ b, gi a b * ∑ e, G1 a e * C3 e b c) -
        ∑ a, ∑ b, gi a b * ∑ e, G1 b e * C3 a e c := by
    simp only [hD]
    exact ip_alg_i gi hs G1 (fun a b => C3 a b c)
  have hS : ∑ e, Gc e * ∑ a, ∑ b, gi a b * C3 a b e =
      ∑ a, ∑ b, gi a b * ∑ e, Gc e * C3 a b e := by
    calc ∑ e, Gc e * ∑ a, ∑ b, gi a b * C3 a b e
        = ∑ e, ∑ a, ∑ b, Gc e * (gi a b * C3 a b e) := by
          simp only [Finset.mul_sum]
      _ = ∑ a, ∑ e, ∑ b, Gc e * (gi a b * C3 a b e) := Finset.sum_comm
      _ = ∑ a, ∑ b, ∑ e, Gc e * (gi a b * C3 a b e) :=
          Finset.sum_congr rfl fun a _ => Finset.sum_comm
      _ = ∑ a, ∑ b, gi a b * ∑ e, Gc e * C3 a b e :=
          Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun b _ => by
            rw [Finset.mul_sum]
            exact Finset.sum_congr rfl fun e _ => by ring
  have e1 : ∑ a, ∑ b, (Dg a b * C3 a b c + gi a b * DC a b) =
      (∑ a, ∑ b, Dg a b * C3 a b c) + ∑ a, ∑ b, gi a b * DC a b := by
    simp only [Finset.sum_add_distrib]
  have e2 : ∑ e, Gc e * (-∑ a, ∑ b, gi a b * C3 a b e) =
      -(∑ e, Gc e * ∑ a, ∑ b, gi a b * C3 a b e) := by
    simp only [mul_neg, Finset.sum_neg_distrib]
  have e3 : ∀ a b, g0 * gi a b *
        (DC a b - ∑ e, (G1 a e * C3 e b c + G1 b e * C3 a e c + Gc e * C3 a b e)) =
      g0 * (gi a b * DC a b - gi a b * ∑ e, G1 a e * C3 e b c -
        gi a b * ∑ e, G1 b e * C3 a e c - gi a b * ∑ e, Gc e * C3 a b e) := by
    intro a b
    rw [Finset.sum_add_distrib, Finset.sum_add_distrib]
    ring
  simp only [e3, ← Finset.mul_sum, Finset.sum_sub_distrib]
  rw [e1, e2, hP, hS]
  ring

/-- `μ` in coordinates, divergence form of `potX`. -/
theorem potX_eq_div (hΩ : IsOpen Ω) (hX : ContDiffOn ℝ ∞ X Ω)
    (hinj : ∀ x ∈ Ω, Function.Injective (fderiv ℝ X x)) (hx : x ∈ Ω) :
    potX J X x = ∑ a, ∑ c, (gX X x)⁻¹ a c *
      (fderiv ℝ (fun y => muC J X y c) x (eb a) - ∑ e, gamX X x ![a, c, e] * muC J X x e) := by
  have hXa := ip_contDiffAt hΩ hX hx
  have hmu : ∀ c d, fderiv ℝ (fun y => muC J X y c) x (eb d) =
      -∑ a, ∑ b, (fderiv ℝ (fun y => (gX X y)⁻¹ a b) x (eb d) * cX J X x ![a, b, c] +
        (gX X x)⁻¹ a b * fderiv ℝ (fun y => cX J X y ![a, b, c]) x (eb d)) := by
    intro c d
    have hd : HasFDerivAt (fun y => muC J X y c)
        (-∑ a, ∑ b, ((gX X x)⁻¹ a b • fderiv ℝ (fun y => cX J X y ![a, b, c]) x +
          cX J X x ![a, b, c] • fderiv ℝ (fun y => (gX X y)⁻¹ a b) x)) x :=
      (HasFDerivAt.fun_sum fun a _ => HasFDerivAt.fun_sum fun b _ =>
        ((ip_diffAt_ginv hΩ hX hinj hx a b).hasFDerivAt).fun_mul
          (ip_diff_cX hXa ![a, b, c]).hasFDerivAt).fun_neg
    rw [hd.fderiv]
    simp only [ContinuousLinearMap.neg_apply, ContinuousLinearMap.sum_apply,
      ContinuousLinearMap.add_apply, ContinuousLinearMap.smul_apply, smul_eq_mul]
    rw [neg_inj]
    exact Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun b _ => by ring
  have hmuX : ∀ e, muC J X x e = -∑ a, ∑ b, (gX X x)⁻¹ a b * cX J X x ![a, b, e] :=
    fun e => rfl
  symm
  unfold potX
  rw [ip_sum4]
  rw [← Finset.sum_neg_distrib]
  refine Finset.sum_congr rfl fun d _ => ?_
  rw [← Finset.sum_neg_distrib]
  refine Finset.sum_congr rfl fun c _ => ?_
  simp only [hmu c d, hmuX, ip_ginvZ_apply, SX_one_apply hΩ hX hinj hx]
  rw [ip_ginv_symm d c]
  exact ip_alg_dc (fun a b => (gX X x)⁻¹ a b) (fun a b => ip_ginv_symm a b)
    (fun e b => gamX X x ![d, e, b]) (fun e => gamX X x ![d, c, e])
    (fun a b e => cX J X x ![a, b, e])
    (fun a b => fderiv ℝ (fun y => cX J X y ![a, b, c]) x (eb d))
    (fun a b => fderiv ℝ (fun y => (gX X y)⁻¹ a b) x (eb d)) c
    (fun a b => fderiv_gX_inv hΩ hX hinj hx a b d) ((gX X x)⁻¹ c d)

/-- info: 'CWFlow.Extension.fderiv_gX_inv' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms fderiv_gX_inv

/-- info: 'CWFlow.Extension.SX_one_apply' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms SX_one_apply

/-- info: 'CWFlow.Extension.potX_eq_div' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms potX_eq_div

end CWFlow.Extension
