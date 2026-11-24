import CWFlow.Extension.ImmDiv
import CWFlow.Scalar.Reparam

/-!
# The potential `div(JH)` in coordinates

For a Lagrangian immersion `X : ℝⁿ → V` (open `Ω`, `X` smooth, `DX` injective, `⟨J ∂_a X, ∂_b X⟩ = 0`):

* Part 1: the basis-free `ginvOf`, `trOf`, `perpOf` in the coordinate frame `eb a`
  (`ginvOf_eq_sum`, `trOf_eq_sum`, `perpOf_eq_sum`);
* Part 2: the mean curvature form `μ_c = ⟨J H, ∂_c X⟩ = −g^{ab} C_{abc}` (`muX`, `muX_eq`);
* Part 3: `div(JH) = g^{ac} (∂_a μ_c − Γ^e_{ac} μ_e)` (`flowPot_eq_div`);
* Part 4: `div(JH) = −g^{cd} g^{ab} (∇_d C)_{abc}` (`potX`, `flowPot_eq_potX`);
* Part 5: the velocity `J ∇ div(JH)` in coordinates (`flowVel_eq_sum`).
-/

open scoped ContDiff Matrix InnerProductSpace RealInnerProductSpace Topology

open CWFlow

namespace CWFlow.Extension

variable {n : ℕ} {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V] {J : V →L[ℝ] V}
  {X : EuclideanSpace ℝ (Fin n) → V} {Ω : Set (EuclideanSpace ℝ (Fin n))}
  {x : EuclideanSpace ℝ (Fin n)}

/-! ### Part 1: the basis-free operations in coordinates -/

theorem ginvOf_eq_sum (hinj : Function.Injective (fderiv ℝ X x))
    (ξ : EuclideanSpace ℝ (Fin n) →L[ℝ] ℝ) :
    ginvOf (fderiv ℝ X x) ξ = ∑ a, ∑ b, ((gX X x)⁻¹ a b * ξ (eb b)) • eb a := by
  apply gOf_injective hinj
  rw [gOf_ginvOf hinj]
  apply ip_ext_eb
  intro c
  rw [gOf_apply]
  have key : ∀ b, ∑ a, (gX X x)⁻¹ a b * gX X x a c = if b = c then 1 else 0 := fun b => by
    rw [← ip_gi_g hinj b c]
    exact Finset.sum_congr rfl fun a _ => by rw [ip_ginv_symm]
  have e : fderiv ℝ X x (∑ a, ∑ b, ((gX X x)⁻¹ a b * ξ (eb b)) • eb a) =
      ∑ a, ∑ b, ((gX X x)⁻¹ a b * ξ (eb b)) • dX X x a := by
    simp [map_sum, map_smul, dX]
  rw [e, sum_inner]
  have e2 : ∀ a, ⟪∑ b, ((gX X x)⁻¹ a b * ξ (eb b)) • dX X x a, fderiv ℝ X x (eb c)⟫ =
      ∑ b, ((gX X x)⁻¹ a b * ξ (eb b)) * gX X x a c := fun a => by
    rw [sum_inner]
    exact Finset.sum_congr rfl fun b _ => by
      rw [inner_smul_left]; rfl
  simp only [e2]
  rw [Finset.sum_comm]
  have e3 : ∀ b, ∑ a, ((gX X x)⁻¹ a b * ξ (eb b)) * gX X x a c =
      ξ (eb b) * (if b = c then 1 else 0) := fun b => by
    rw [← key b, Finset.mul_sum]
    exact Finset.sum_congr rfl fun a _ => by ring
  simp only [e3]
  simp

theorem trOf_eq_sum (hinj : Function.Injective (fderiv ℝ X x)) {W : Type*} [AddCommGroup W]
    [Module ℝ W]
    (B : EuclideanSpace ℝ (Fin n) →ₗ[ℝ] EuclideanSpace ℝ (Fin n) →ₗ[ℝ] W) :
    trOf (fderiv ℝ X x) (fun a c => B a c) = ∑ a, ∑ b, (gX X x)⁻¹ a b • B (eb a) (eb b) := by
  have h := sum_frame_eq (fun a : Fin n => eb a)
    (fun a : Fin n => (EuclideanSpace.proj a : EuclideanSpace ℝ (Fin n) →L[ℝ] ℝ))
    (fun v => ip_sum_eb v) B
    (ginvOf (fderiv ℝ X x) : (EuclideanSpace ℝ (Fin n) →L[ℝ] ℝ) →ₗ[ℝ] EuclideanSpace ℝ (Fin n))
  refine Eq.trans h.symm ?_
  refine Finset.sum_congr rfl fun a _ => ?_
  have e : ginvOf (fderiv ℝ X x) (EuclideanSpace.proj a : EuclideanSpace ℝ (Fin n) →L[ℝ] ℝ) =
      ∑ a', (gX X x)⁻¹ a a' • eb a' := by
    rw [ginvOf_eq_sum hinj]
    refine Finset.sum_congr rfl fun a' _ => ?_
    rw [← Finset.sum_smul]
    congr 1
    have hp : ∀ i, (EuclideanSpace.proj a : EuclideanSpace ℝ (Fin n) →L[ℝ] ℝ) (eb i) =
        if a = i then 1 else 0 := fun i => by
      show (eb i) a = _
      simp [Pi.single_apply]
    simp only [hp, mul_ite, mul_one, mul_zero, Finset.sum_ite_eq, Finset.mem_univ, if_true]
    exact ip_ginv_symm a' a
  show B (eb a) (ginvOf (fderiv ℝ X x) _) = _
  rw [e]
  simp [map_sum, map_smul]

theorem perpOf_eq_sum (hinj : Function.Injective (fderiv ℝ X x)) (v : V) :
    perpOf (fderiv ℝ X x) v =
      v - ∑ a, ∑ b, ((gX X x)⁻¹ a b * ⟪v, dX X x b⟫) • dX X x a := by
  unfold perpOf
  rw [ginvOf_eq_sum hinj]
  congr 1
  simp [map_sum, map_smul, dX]

/-! ### Part 2: the mean curvature form -/

/-- `μ_c = ⟨J H, ∂_c X⟩`. -/
noncomputable def muX (J : V →L[ℝ] V) (X : EuclideanSpace ℝ (Fin n) → V)
    (x : EuclideanSpace ℝ (Fin n)) (c : Fin n) : ℝ := ⟪J (meanCurv X x), dX X x c⟫

theorem meanCurv_eq_sum (hinj : Function.Injective (fderiv ℝ X x)) :
    meanCurv X x = ∑ a, ∑ b, (gX X x)⁻¹ a b • perpOf (fderiv ℝ X x) (ddX X x a b) := by
  exact trOf_eq_sum (X := X) (x := x) hinj
    (((fderiv ℝ (fderiv ℝ X) x : EuclideanSpace ℝ (Fin n) →L[ℝ] EuclideanSpace ℝ (Fin n) →L[ℝ] V
      ).toLinearMap₁₂.compr₂ (perpLin (fderiv ℝ X x))))

theorem ip_inner_perp_L (hinj : Function.Injective (fderiv ℝ X x)) (w : V)
    (u : EuclideanSpace ℝ (Fin n)) : ⟪perpOf (fderiv ℝ X x) w, fderiv ℝ X x u⟫ = 0 := by
  unfold perpOf
  rw [inner_sub_left, inner_ginvOf hinj]
  simp

theorem inner_meanCurv_dX (hinj : Function.Injective (fderiv ℝ X x)) (c : Fin n) :
    ⟪meanCurv X x, dX X x c⟫ = 0 := by
  rw [meanCurv_eq_sum hinj, sum_inner]
  refine Finset.sum_eq_zero fun a _ => ?_
  rw [sum_inner]
  refine Finset.sum_eq_zero fun b _ => ?_
  rw [inner_smul_left]
  have := ip_inner_perp_L hinj (ddX X x a b) (eb c)
  simp [dX, this]

theorem ip_inner_perp_J (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    (hinj : Function.Injective (fderiv ℝ X x))
    (hlag : ∀ u w, ⟪J (fderiv ℝ X x u), fderiv ℝ X x w⟫ = 0) (w : V) (c : Fin n) :
    ⟪perpOf (fderiv ℝ X x) w, J (dX X x c)⟫ = ⟪w, J (dX X x c)⟫ := by
  unfold perpOf
  have h := ip_inner_L_J_L hJi hJJ hlag
    (ginvOf (fderiv ℝ X x) (((innerSL ℝ : V →L[ℝ] V →L[ℝ] ℝ) w).comp (fderiv ℝ X x))) (eb c)
  rw [inner_sub_left]
  unfold dX
  linarith

theorem muX_eq (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    (hΩ : IsOpen Ω) (hX : ContDiffOn ℝ ∞ X Ω)
    (hinj : ∀ x ∈ Ω, Function.Injective (fderiv ℝ X x))
    (hlag : ∀ x ∈ Ω, ∀ u w, ⟪J (fderiv ℝ X x u), fderiv ℝ X x w⟫ = 0) (hx : x ∈ Ω) (c : Fin n) :
    muX J X x c = -∑ a, ∑ b, (gX X x)⁻¹ a b * cX J X x ![a, b, c] := by
  unfold muX
  rw [ip_inner_J_left hJi hJJ, neg_inj, meanCurv_eq_sum (hinj x hx), sum_inner]
  refine Finset.sum_congr rfl fun a _ => ?_
  rw [sum_inner]
  refine Finset.sum_congr rfl fun b _ => ?_
  rw [inner_smul_left, ip_inner_perp_J hJi hJJ (hinj x hx) (hlag x hx)]
  rfl

/-! ### Part 3: the potential as a divergence -/

theorem flowPot_eq_div (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    (hΩ : IsOpen Ω) (hX : ContDiffOn ℝ ∞ X Ω)
    (hinj : ∀ x ∈ Ω, Function.Injective (fderiv ℝ X x)) (hx : x ∈ Ω) :
    flowPot J X x = ∑ a, ∑ c, (gX X x)⁻¹ a c *
      (fderiv ℝ (fun y => muX J X y c) x (eb a) - ∑ e, gamX X x ![a, c, e] * muX J X x e) := by
  have hXa := ip_contDiffAt hΩ hX hx
  have hdJ := differentiableAt_J_meanCurv J hXa (hinj x hx)
  have hL := hinj x hx
  have hflow : flowPot J X x = ∑ a, ∑ c, (gX X x)⁻¹ a c *
      ⟪fderiv ℝ (fun z => J (meanCurv X z)) x (eb a), dX X x c⟫ := by
    let B : EuclideanSpace ℝ (Fin n) →ₗ[ℝ] EuclideanSpace ℝ (Fin n) →ₗ[ℝ] ℝ :=
      LinearMap.mk₂ ℝ (fun a c => ⟪fderiv ℝ (fun z => J (meanCurv X z)) x a, fderiv ℝ X x c⟫)
        (fun a a' c => by simp [inner_add_left]) (fun r a c => by simp [inner_smul_left])
        (fun a c c' => by simp [inner_add_right]) (fun r a c => by simp [inner_smul_right])
    have e : flowPot J X x = trOf (fderiv ℝ X x) (fun a c => B a c) := rfl
    rw [e, trOf_eq_sum hL B]
    rfl
  have hprod : ∀ a c, fderiv ℝ (fun y => muX J X y c) x (eb a) =
      ⟪J (meanCurv X x), ddX X x a c⟫ +
        ⟪fderiv ℝ (fun z => J (meanCurv X z)) x (eb a), dX X x c⟫ := fun a c => by
    have h := fderiv_inner_apply ℝ hdJ (ip_diff_dX hXa c) (eb a)
    rw [ip_fderiv_dX hXa] at h
    exact h
  have hdec : ∀ a c, ⟪J (meanCurv X x), ddX X x a c⟫ =
      ∑ e, gamX X x ![a, c, e] * muX J X x e +
        ⟪J (meanCurv X x), perpOf (fderiv ℝ X x) (ddX X x a c)⟫ := fun a c => by
    have hv : ddX X x a c = perpOf (fderiv ℝ X x) (ddX X x a c) +
        ∑ e, ∑ f, ((gX X x)⁻¹ e f * ⟪ddX X x a c, dX X x f⟫) • dX X x e := by
      rw [perpOf_eq_sum hL]
      abel
    have h1 : ⟪J (meanCurv X x), ddX X x a c⟫ =
        ⟪J (meanCurv X x), perpOf (fderiv ℝ X x) (ddX X x a c)⟫ +
          ∑ e, ∑ f, ((gX X x)⁻¹ e f * ⟪ddX X x a c, dX X x f⟫) * muX J X x e := by
      conv_lhs => rw [hv]
      rw [inner_add_right, inner_sum]
      congr 1
      refine Finset.sum_congr rfl fun e _ => ?_
      rw [inner_sum]
      exact Finset.sum_congr rfl fun f _ => inner_smul_right _ _ _
    rw [h1, add_comm]
    congr 1
    refine Finset.sum_congr rfl fun e _ => ?_
    rw [← Finset.sum_mul, ip_gamX_eq]
  have hperp0 : ∑ a, ∑ c, (gX X x)⁻¹ a c *
      ⟪J (meanCurv X x), perpOf (fderiv ℝ X x) (ddX X x a c)⟫ = 0 := by
    have h : ⟪J (meanCurv X x), ∑ a, ∑ c, (gX X x)⁻¹ a c • perpOf (fderiv ℝ X x) (ddX X x a c)⟫
        = 0 := by
      rw [← meanCurv_eq_sum hL]
      exact ip_inner_J_self hJi hJJ _
    simp only [inner_sum, inner_smul_right] at h
    exact h
  rw [hflow]
  have hterm : ∀ a c, (gX X x)⁻¹ a c * ⟪fderiv ℝ (fun z => J (meanCurv X z)) x (eb a), dX X x c⟫ =
      (gX X x)⁻¹ a c * (fderiv ℝ (fun y => muX J X y c) x (eb a) -
        ∑ e, gamX X x ![a, c, e] * muX J X x e) -
      (gX X x)⁻¹ a c * ⟪J (meanCurv X x), perpOf (fderiv ℝ X x) (ddX X x a c)⟫ := by
    intro a c
    rw [hprod, hdec a c]
    ring
  simp only [hterm, Finset.sum_sub_distrib]
  rw [hperp0, sub_zero]


/-! ### Part 4: the main formula -/

theorem flowPot_eq_potX (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    (hΩ : IsOpen Ω) (hX : ContDiffOn ℝ ∞ X Ω)
    (hinj : ∀ x ∈ Ω, Function.Injective (fderiv ℝ X x))
    (hlag : ∀ x ∈ Ω, ∀ u w, ⟪J (fderiv ℝ X x u), fderiv ℝ X x w⟫ = 0) (hx : x ∈ Ω) :
    flowPot J X x = potX J X x := by
  have hev : ∀ c, (fun y => muX J X y c) =ᶠ[𝓝 x] fun y => muC J X y c := fun c => by
    filter_upwards [hΩ.mem_nhds hx] with y hy
    exact muX_eq hJi hJJ hΩ hX hinj hlag hy c
  have hmuX : ∀ e, muX J X x e = muC J X x e := fun e => muX_eq hJi hJJ hΩ hX hinj hlag hx e
  rw [flowPot_eq_div hJi hJJ hΩ hX hinj hx, potX_eq_div hΩ hX hinj hx]
  simp only [(hev _).fderiv_eq, hmuX]

/-! ### Part 5: the velocity -/

theorem flowVel_eq_sum (hinj : Function.Injective (fderiv ℝ X x)) :
    flowVel J X x =
      ∑ a, ∑ b, ((gX X x)⁻¹ a b * fderiv ℝ (flowPot J X) x (eb b)) • J (dX X x a) := by
  unfold flowVel
  rw [ginvOf_eq_sum hinj]
  simp [map_sum, map_smul, dX]

theorem flowVel_eq_velX (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    (hΩ : IsOpen Ω) (hX : ContDiffOn ℝ ∞ X Ω)
    (hinj : ∀ x ∈ Ω, Function.Injective (fderiv ℝ X x))
    (hlag : ∀ x ∈ Ω, ∀ u w, ⟪J (fderiv ℝ X x u), fderiv ℝ X x w⟫ = 0) (hx : x ∈ Ω) :
    flowVel J X x = velX J X x := by
  have hev : flowPot J X =ᶠ[𝓝 x] potX J X := by
    filter_upwards [hΩ.mem_nhds hx] with y hy
    exact flowPot_eq_potX hJi hJJ hΩ hX hinj hlag hy
  rw [flowVel_eq_sum (hinj x hx), hev.fderiv_eq]
  rfl

/-- info: 'CWFlow.Extension.trOf_eq_sum' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms trOf_eq_sum

/-- info: 'CWFlow.Extension.muX_eq' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms muX_eq

/-- info: 'CWFlow.Extension.flowPot_eq_div' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms flowPot_eq_div

/-- info: 'CWFlow.Extension.flowPot_eq_potX' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms flowPot_eq_potX

/-- info: 'CWFlow.Extension.flowVel_eq_sum' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms flowVel_eq_sum

/-- info: 'CWFlow.Extension.flowVel_eq_velX' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms flowVel_eq_velX

end CWFlow.Extension
