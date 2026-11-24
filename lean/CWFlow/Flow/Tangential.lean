import CWFlow.Scalar.Pointwise
import CWFlow.Scalar.Graph

/-!
# The graph velocity minus the flow velocity is tangential (Chen–Warren, Claim 3.3, pointwise part)

For a Lagrangian immersion with derivative `L : E → V` (`dim V = 2 dim E`), a vector `w` with
`⟨J w, L ·⟩ = 0` is tangent: `V = range L ⊕ J(range L)` orthogonally, and `w ⊥ J(range L)`
(`mem_range_of_inner_J_eq_zero`). Its tangential coordinates are `tanPart L w = g⁻¹⟨w, L ·⟩`.

For a graph `F_u = Ψ(·, Du)` under a symplectic `Ψ` (`DF_u a = DΨ(a, D²u a)`), the vertical vector
`DΨ(0, η)` and `J L g⁻¹δ` have the same pairing `⟨J ·, L a⟩` up to sign: `η a` and `−δ a`. So for
`η = −δ`, the difference `DΨ(0, η) − J L g⁻¹ δ` is tangential (`vertical_sub_flowVel_mem_range`,
`graphFu_vertical_sub_flowVel`). With `η = D(∂ₜu)`, `δ = D div(JH)` (the scalar equation
`∂ₜu = −div(JH)`), the graph velocity `∂ₜF_u` differs from `flowVel = J∇div(JH)` by a tangent vector.
-/

open Filter Set

open scoped InnerProductSpace RealInnerProductSpace Topology ContDiff

noncomputable section

namespace CWFlow

variable {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]

/-- The tangential coordinates `g⁻¹⟨v, L ·⟩` of `v` with respect to `L`. -/
def tanPart (L : E →L[ℝ] V) (v : V) : E :=
  ginvOf L (((innerSL ℝ : V →L[ℝ] V →L[ℝ] ℝ) v).comp L)

omit [FiniteDimensional ℝ E] in
theorem perpOf_eq_sub_tanPart (L : E →L[ℝ] V) (v : V) : perpOf L v = v - L (tanPart L v) := rfl

/-- `v − L (tanPart L v)` is orthogonal to the range of `L`. -/
theorem inner_sub_tanPart {L : E →L[ℝ] V} (hL : Function.Injective L) (v : V) (a : E) :
    ⟪v - L (tanPart L v), L a⟫ = 0 := by
  rw [inner_sub_left, tanPart, inner_ginvOf hL]
  simp

/-! ### Lagrangian splitting -/

section Splitting

variable {J : V →L[ℝ] V}

/-- A vector orthogonal to `range L` and to `J(range L)` vanishes, for `L` injective and Lagrangian with
`dim V = 2 dim E`. -/
theorem eq_zero_of_inner_eq_zero_of_inner_J_eq_zero [FiniteDimensional ℝ V]
    (hJJ : ∀ v, J (J v) = -v) {L : E →L[ℝ] V} (hL : Function.Injective L)
    (hLag : ∀ a b, ⟪J (L a), L b⟫ = 0)
    (hdim : Module.finrank ℝ V = 2 * Module.finrank ℝ E) {p : V}
    (h1 : ∀ a, ⟪p, L a⟫ = 0) (h2 : ∀ a, ⟪p, J (L a)⟫ = 0) : p = 0 := by
  let M : E × E →ₗ[ℝ] V := (L : E →ₗ[ℝ] V).coprod ((J.comp L : E →L[ℝ] V) : E →ₗ[ℝ] V)
  have hM : ∀ x : E × E, M x = L x.1 + J (L x.2) := fun _ => rfl
  have hMinj : Function.Injective M := by
    rw [injective_iff_map_eq_zero]
    rintro ⟨a, b⟩ hab
    rw [hM] at hab
    have ha : L a = 0 := by
      have h := congrArg (fun x => ⟪x, L a⟫) hab
      simp only [inner_add_left, hLag, add_zero, inner_zero_left] at h
      exact inner_self_eq_zero.mp h
    have hb : L b = 0 := by
      rw [ha, zero_add] at hab
      have h := congrArg J hab
      rw [hJJ, map_zero] at h
      exact neg_eq_zero.mp h
    have ha' : a = 0 := hL (by rw [ha, map_zero])
    have hb' : b = 0 := hL (by rw [hb, map_zero])
    rw [ha', hb', Prod.mk_zero_zero]
  have hMsurj : Function.Surjective M :=
    (LinearMap.injective_iff_surjective_of_finrank_eq_finrank
      (by rw [Module.finrank_prod, hdim]; ring)).mp hMinj
  obtain ⟨x, hx⟩ := hMsurj p
  have hpp : ⟪p, p⟫ = 0 := by
    nth_rewrite 2 [← hx]
    rw [hM, inner_add_right, h1, h2, add_zero]
  exact inner_self_eq_zero.mp hpp

/-- **Lagrangian splitting.** For `L` injective and Lagrangian with `dim V = 2 dim E`, a vector `w` with
`⟨J w, L ·⟩ = 0` is tangent: `w = L (tanPart L w)`. -/
theorem mem_range_of_inner_J_eq_zero [FiniteDimensional ℝ V]
    (hJi : ∀ u v, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v, J (J v) = -v)
    {L : E →L[ℝ] V} (hL : Function.Injective L) (hLag : ∀ a b, ⟪J (L a), L b⟫ = 0)
    (hdim : Module.finrank ℝ V = 2 * Module.finrank ℝ E) {w : V}
    (hw : ∀ a, ⟪J w, L a⟫ = 0) : w = L (tanPart L w) := by
  refine (sub_eq_zero.mp (eq_zero_of_inner_eq_zero_of_inner_J_eq_zero hJJ hL hLag hdim
    (inner_sub_tanPart hL w) fun a => ?_))
  have e1 : ⟪w, J (L a)⟫ = 0 := by
    have h := hJi (J w) (L a)
    rw [hJJ, inner_neg_left, hw] at h
    exact neg_eq_zero.mp h
  have e2 : ⟪L (tanPart L w), J (L a)⟫ = 0 := by
    rw [real_inner_comm]
    exact hLag _ _
  rw [inner_sub_left, e1, e2, sub_zero]

end Splitting

/-! ### The pointwise tangential lemma -/

section Graph

variable {J : V →L[ℝ] V} {Ψ : E × (E →L[ℝ] ℝ) → V} {q : E × (E →L[ℝ] ℝ)}

/-- A linear map of the form `a ↦ DΨ(a, H a)` is injective when `Ψ` is symplectic at `q`. -/
theorem injective_of_fderiv_prod
    (hs : ∀ a b, ⟪J (fderiv ℝ Ψ q a), fderiv ℝ Ψ q b⟫ = omegaCanVal a b)
    {H : E →L[ℝ] (E →L[ℝ] ℝ)} {L : E →L[ℝ] V} (hLdef : ∀ a, L a = fderiv ℝ Ψ q (a, H a)) :
    Function.Injective L := by
  intro a b hab
  rw [hLdef, hLdef] at hab
  exact congrArg Prod.fst (injective_fderiv_of_symplectic hs hab)

omit [FiniteDimensional ℝ E] in
/-- A linear map of the form `a ↦ DΨ(a, H a)` with `H` symmetric is Lagrangian when `Ψ` is symplectic
at `q`. -/
theorem lagrangian_of_fderiv_prod
    (hs : ∀ a b, ⟪J (fderiv ℝ Ψ q a), fderiv ℝ Ψ q b⟫ = omegaCanVal a b)
    {H : E →L[ℝ] (E →L[ℝ] ℝ)} (hH : ∀ a b, H a b = H b a)
    {L : E →L[ℝ] V} (hLdef : ∀ a, L a = fderiv ℝ Ψ q (a, H a)) (a b : E) :
    ⟪J (L a), L b⟫ = 0 := by
  rw [hLdef, hLdef, hs]
  simp only [omegaCanVal]
  rw [hH a b, sub_self]

/-- `⟨J (DΨ(0, η) − J L g⁻¹δ), L a⟫ = η a + δ a = 0` for `η = −δ`. -/
theorem inner_J_vertical_sub_eq_zero (hJJ : ∀ v, J (J v) = -v)
    (hs : ∀ a b, ⟪J (fderiv ℝ Ψ q a), fderiv ℝ Ψ q b⟫ = omegaCanVal a b)
    {H : E →L[ℝ] (E →L[ℝ] ℝ)} {L : E →L[ℝ] V} (hLdef : ∀ a, L a = fderiv ℝ Ψ q (a, H a))
    (η δ : E →L[ℝ] ℝ) (hηδ : η = -δ) (a : E) :
    ⟪J (fderiv ℝ Ψ q (0, η) - J (L (ginvOf L δ))), L a⟫ = 0 := by
  rw [map_sub, hJJ, sub_neg_eq_add, inner_add_left,
    inner_ginvOf (injective_of_fderiv_prod hs hLdef), hLdef a, inner_J_fderiv_vertical hs, hηδ]
  simp

/-- **Pointwise tangential lemma**, for `L` given by `L a = DΨ(a, H a)`: with `η = −δ`,
`DΨ(0, η) − J L g⁻¹δ` is tangent to `L`. -/
theorem vertical_sub_eq_tanPart [FiniteDimensional ℝ V]
    (hJi : ∀ u v, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v, J (J v) = -v)
    (hdim : Module.finrank ℝ V = 2 * Module.finrank ℝ E)
    (hs : ∀ a b, ⟪J (fderiv ℝ Ψ q a), fderiv ℝ Ψ q b⟫ = omegaCanVal a b)
    {H : E →L[ℝ] (E →L[ℝ] ℝ)} (hH : ∀ a b, H a b = H b a)
    {L : E →L[ℝ] V} (hLdef : ∀ a, L a = fderiv ℝ Ψ q (a, H a))
    (η δ : E →L[ℝ] ℝ) (hηδ : η = -δ) :
    fderiv ℝ Ψ q (0, η) - J (L (ginvOf L δ)) =
      L (tanPart L (fderiv ℝ Ψ q (0, η) - J (L (ginvOf L δ)))) :=
  mem_range_of_inner_J_eq_zero hJi hJJ (injective_of_fderiv_prod hs hLdef)
    (lagrangian_of_fderiv_prod hs hH hLdef) hdim
    (inner_J_vertical_sub_eq_zero hJJ hs hLdef η δ hηδ)

/-- **Pointwise tangential lemma** for `L = DΨ ∘ (id, H)` (the derivative of a graph immersion,
`fderiv_graphFu_apply`): with `η = −δ`, `v = DΨ(0, η) − J L g⁻¹δ` satisfies `v = L (tanPart L v)`. -/
theorem vertical_sub_flowVel_mem_range [FiniteDimensional ℝ V]
    (hJi : ∀ u v, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v, J (J v) = -v)
    (hdim : Module.finrank ℝ V = 2 * Module.finrank ℝ E)
    (hs : ∀ a b, ⟪J (fderiv ℝ Ψ q a), fderiv ℝ Ψ q b⟫ = omegaCanVal a b)
    {H : E →L[ℝ] (E →L[ℝ] ℝ)} (hH : ∀ a b, H a b = H b a)
    (η δ : E →L[ℝ] ℝ) (hηδ : η = -δ) :
    fderiv ℝ Ψ q (0, η) -
        J (((fderiv ℝ Ψ q).comp ((ContinuousLinearMap.id ℝ E).prod H))
          (ginvOf ((fderiv ℝ Ψ q).comp ((ContinuousLinearMap.id ℝ E).prod H)) δ)) =
      ((fderiv ℝ Ψ q).comp ((ContinuousLinearMap.id ℝ E).prod H))
        (tanPart ((fderiv ℝ Ψ q).comp ((ContinuousLinearMap.id ℝ E).prod H))
          (fderiv ℝ Ψ q (0, η) -
            J (((fderiv ℝ Ψ q).comp ((ContinuousLinearMap.id ℝ E).prod H))
              (ginvOf ((fderiv ℝ Ψ q).comp ((ContinuousLinearMap.id ℝ E).prod H)) δ)))) :=
  vertical_sub_eq_tanPart hJi hJJ hdim hs hH (fun _ => rfl) η δ hηδ

/-! ### On the graph `F_u` -/

variable {u : E → ℝ} {y : E}

/-- **The graph velocity minus the flow velocity is tangential** (pointwise). For `F_u = Ψ(·, Du)` with
`Ψ` symplectic at `(y, Du(y))` and `η = −δ`:
`DΨ(0, η) − J DF_u g⁻¹δ = DF_u (tanPart DF_u (DΨ(0, η) − J DF_u g⁻¹δ))`. -/
theorem graphFu_vertical_sub_tangential [FiniteDimensional ℝ V]
    (hJi : ∀ u v, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v, J (J v) = -v)
    (hdim : Module.finrank ℝ V = 2 * Module.finrank ℝ E)
    (hu : ContDiffAt ℝ 2 u y) (hΨ : DifferentiableAt ℝ Ψ (graphMap u y))
    (hs : ∀ a b, ⟪J (fderiv ℝ Ψ (graphMap u y) a), fderiv ℝ Ψ (graphMap u y) b⟫ = omegaCanVal a b)
    (η δ : E →L[ℝ] ℝ) (hηδ : η = -δ) :
    fderiv ℝ Ψ (graphMap u y) (0, η) -
        J (fderiv ℝ (graphFu Ψ u) y (ginvOf (fderiv ℝ (graphFu Ψ u) y) δ)) =
      fderiv ℝ (graphFu Ψ u) y (tanPart (fderiv ℝ (graphFu Ψ u) y)
        (fderiv ℝ Ψ (graphMap u y) (0, η) -
          J (fderiv ℝ (graphFu Ψ u) y (ginvOf (fderiv ℝ (graphFu Ψ u) y) δ)))) :=
  vertical_sub_eq_tanPart hJi hJJ hdim hs (hu.isSymmSndFDerivAt (by simp))
    (fderiv_graphFu_apply hΨ hu) η δ hηδ

/-- The same with `δ = D div(JH)[F_u](y)`, so that the `J` term is `flowVel J F_u y`: for
`η = −D div(JH)[F_u](y)`, `DΨ(0, η) − flowVel J F_u y` is tangent to `F_u` at `y`. -/
theorem graphFu_vertical_sub_flowVel [FiniteDimensional ℝ V]
    (hJi : ∀ u v, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v, J (J v) = -v)
    (hdim : Module.finrank ℝ V = 2 * Module.finrank ℝ E)
    (hu : ContDiffAt ℝ 2 u y) (hΨ : DifferentiableAt ℝ Ψ (graphMap u y))
    (hs : ∀ a b, ⟪J (fderiv ℝ Ψ (graphMap u y) a), fderiv ℝ Ψ (graphMap u y) b⟫ = omegaCanVal a b)
    (η : E →L[ℝ] ℝ) (hη : η = -fderiv ℝ (flowPot J (graphFu Ψ u)) y) :
    fderiv ℝ Ψ (graphMap u y) (0, η) - flowVel J (graphFu Ψ u) y =
      fderiv ℝ (graphFu Ψ u) y (tanPart (fderiv ℝ (graphFu Ψ u) y)
        (fderiv ℝ Ψ (graphMap u y) (0, η) - flowVel J (graphFu Ψ u) y)) :=
  graphFu_vertical_sub_tangential hJi hJJ hdim hu hΨ hs η _ hη

end Graph

end CWFlow
