import CWFlow.Extension.ImmCoord
import CWFlow.Extension.EvolPointwise

/-!
# The graph bridge

For the graph `X(x) = y + U x + J U ∇φ(x)` of a gradient over a Lagrangian plane, the coordinate tensors
of the immersion `X` defined in `ImmCoord` (metric, Christoffel symbols, Lagrangian second fundamental
form, its covariant derivatives, energy density) are the symbols on jet space (`gMat`, `gamFun`, `tSym`,
`eSym`) evaluated on the jets of `φ`.
-/

open CWFlow.Interp
open scoped ContDiff Matrix RealInnerProductSpace Topology

namespace CWFlow.Extension

variable {n : ℕ} {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]

/-! ### Part 0: algebra of the frame -/

theorem gb_inner_J_left {J : V →L[ℝ] V} (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫)
    (hJJ : ∀ v : V, J (J v) = -v) (u v : V) : ⟪J u, v⟫ = -⟪u, J v⟫ := by
  have h := hJi u (J v)
  rw [hJJ v, inner_neg_right] at h
  rw [← h, neg_neg]

theorem gb_UU {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} (a b : EuclideanSpace ℝ (Fin n)) :
    ⟪U a, U b⟫ = ⟪a, b⟫ := U.inner_map_map a b

theorem gb_JUJU {J : V →L[ℝ] V} (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫)
    {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} (a b : EuclideanSpace ℝ (Fin n)) :
    ⟪J (U a), J (U b)⟫ = ⟪a, b⟫ := by
  rw [hJi, U.inner_map_map]

theorem gb_UJU {J : V →L[ℝ] V} {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} (hU : IsLagFrame J U)
    (a b : EuclideanSpace ℝ (Fin n)) : ⟪U a, J (U b)⟫ = 0 := by
  rw [real_inner_comm]; exact hU.lag b a

theorem gb_JUU {J : V →L[ℝ] V} {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} (hU : IsLagFrame J U)
    (a b : EuclideanSpace ℝ (Fin n)) : ⟪J (U a), U b⟫ = 0 := hU.lag a b

/-! ### Part 1: derivatives of the graph map -/

theorem gb_sum_eb (t : Fin n → ℝ) : ∑ c, t c • eb c = WithLp.toLp 2 t := by
  ext l
  simp [eb, Finset.sum_apply, Pi.single_apply]

theorem gb_jet_one {φ : EuclideanSpace ℝ (Fin n) → ℝ} (x : EuclideanSpace ℝ (Fin n)) (c : Fin n) :
    jet 1 φ x ![c] = fderiv ℝ φ x (eb c) := by
  rw [jet_apply, iteratedFDeriv_one_apply]
  rfl

theorem contDiff_gradE {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ) :
    ContDiff ℝ ∞ (gradE φ) := by
  refine contDiff_piLp' (𝕜 := ℝ) 2 fun c => ?_
  exact (hφ.fderiv_right (m := ∞) (by simp)).clm_apply contDiff_const

theorem contDiff_graphX {J : V →L[ℝ] V} {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} {y : V}
    {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ) : ContDiff ℝ ∞ (graphX J U y φ) := by
  unfold graphX
  exact (contDiff_const.add U.toContinuousLinearMap.contDiff).add
    ((J.comp U.toContinuousLinearMap).contDiff.comp (contDiff_gradE hφ))

theorem gb_gradE_comp {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ)
    (x : EuclideanSpace ℝ (Fin n)) (a : Fin n) :
    fderiv ℝ (gradE φ) x (eb a) = WithLp.toLp 2 (fun c => jet 2 φ x ![a, c]) := by
  ext c
  rw [show (fderiv ℝ (gradE φ) x (eb a)).ofLp c = fderiv ℝ (fun z => gradE φ z c) x (eb a) from
    fderiv_piLp_apply ((contDiff_gradE hφ).differentiable (by simp) x) (eb a) c]
  have h1 : jet 2 φ x ![a, c] = fderiv ℝ (fun z => jet 1 φ z ![c]) x (eb a) :=
    jet_cons_eq_fderiv hφ 1 x a ![c]
  rw [h1]
  congr 1
  congr 1
  funext z
  rw [gb_jet_one]
  rfl

theorem gb_fderiv_graphX {J : V →L[ℝ] V} {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} {y : V}
    {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ) (x : EuclideanSpace ℝ (Fin n)) :
    fderiv ℝ (graphX J U y φ) x =
      U.toContinuousLinearMap + (J.comp U.toContinuousLinearMap).comp (fderiv ℝ (gradE φ) x) := by
  have hg : DifferentiableAt ℝ (gradE φ) x := (contDiff_gradE hφ).differentiable (by simp) x
  have h1 : HasFDerivAt (fun z => (J.comp U.toContinuousLinearMap) (gradE φ z))
      ((J.comp U.toContinuousLinearMap).comp (fderiv ℝ (gradE φ) x)) x :=
    (J.comp U.toContinuousLinearMap).hasFDerivAt.comp x hg.hasFDerivAt
  have h2 : HasFDerivAt (fun z => U z) U.toContinuousLinearMap x :=
    U.toContinuousLinearMap.hasFDerivAt
  have h3 := (h2.const_add y).add h1
  exact (show graphX J U y φ = fun z => (y + U z) + (J.comp U.toContinuousLinearMap) (gradE φ z) from
    by funext z; simp [graphX, add_assoc]) ▸ h3.fderiv

theorem dX_graphX {J : V →L[ℝ] V} {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} {y : V}
    {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ) (x : EuclideanSpace ℝ (Fin n))
    (a : Fin n) :
    dX (graphX J U y φ) x a = U (eb a) + J (U (∑ c, jet 2 φ x ![a, c] • eb c)) := by
  unfold dX
  rw [gb_fderiv_graphX hφ, gb_sum_eb, ← gb_gradE_comp hφ]
  rfl

theorem ddX_graphX {J : V →L[ℝ] V} {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} {y : V}
    {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ) (x : EuclideanSpace ℝ (Fin n))
    (a b : Fin n) :
    ddX (graphX J U y φ) x a b = J (U (∑ c, jet 3 φ x ![a, b, c] • eb c)) := by
  have hd : DifferentiableAt ℝ (fderiv ℝ (graphX J U y φ)) x :=
    ((contDiff_graphX (J := J) (U := U) (y := y) hφ).fderiv_right (m := ∞) (by simp)).differentiable
      (by simp) x
  have h1 : fderiv ℝ (fun z => fderiv ℝ (graphX J U y φ) z (eb b)) x (eb a) =
      ddX (graphX J U y φ) x a b := by
    rw [fderiv_clm_apply hd (differentiableAt_const _)]
    simp [ddX]
  rw [← h1]
  have h2 : (fun z => fderiv ℝ (graphX J U y φ) z (eb b)) =
      fun z => U (eb b) + (J.comp U.toContinuousLinearMap)
        (WithLp.toLp 2 (fun c => jet 2 φ z ![b, c])) := by
    funext z
    have := dX_graphX (J := J) (U := U) (y := y) hφ z b
    rw [gb_sum_eb] at this
    exact this
  rw [h2]
  have hj : DifferentiableAt ℝ (fun z => (WithLp.toLp 2 (fun c => jet 2 φ z ![b, c]) :
      EuclideanSpace ℝ (Fin n))) x := by
    refine (differentiableAt_piLp 2).2 fun c => ?_
    have := (contDiff_jet hφ 2).differentiable (by simp) x
    exact ((differentiableAt_piLp 2).1 this ![b, c])
  have h3 : fderiv ℝ (fun z => (J.comp U.toContinuousLinearMap)
      (WithLp.toLp 2 (fun c => jet 2 φ z ![b, c]) : EuclideanSpace ℝ (Fin n))) x =
      (J.comp U.toContinuousLinearMap).comp (fderiv ℝ (fun z => (WithLp.toLp 2
        (fun c => jet 2 φ z ![b, c]) : EuclideanSpace ℝ (Fin n))) x) :=
    ((J.comp U.toContinuousLinearMap).hasFDerivAt.comp x hj.hasFDerivAt).fderiv
  rw [fderiv_const_add, h3]
  simp only [ContinuousLinearMap.comp_apply]
  rw [gb_sum_eb]
  refine congrArg J (congrArg U ?_)
  ext c
  rw [show (fderiv ℝ (fun z => (WithLp.toLp 2 (fun c => jet 2 φ z ![b, c]) :
      EuclideanSpace ℝ (Fin n))) x (eb a)).ofLp c =
      fderiv ℝ (fun z => jet 2 φ z ![b, c]) x (eb a) from
    fderiv_piLp_apply hj (eb a) c]
  exact (jet_cons_eq_fderiv hφ 2 x a ![b, c]).symm


/-! ### Part 2: the metric, the Christoffel symbols, the form `C` -/

theorem gb_inner_toLp (f g : Fin n → ℝ) :
    ⟪(WithLp.toLp 2 f : EuclideanSpace ℝ (Fin n)), WithLp.toLp 2 g⟫ = ∑ i, f i * g i := by
  simp [PiLp.inner_apply, mul_comm]

theorem gb_inner_toLp_eb (g : Fin n → ℝ) (c : Fin n) :
    ⟪(WithLp.toLp 2 g : EuclideanSpace ℝ (Fin n)), eb c⟫ = g c := by
  simp [eb, EuclideanSpace.inner_single_right]

theorem gb_inner_eb_eb (a b : Fin n) : ⟪eb a, eb b⟫ = if a = b then (1 : ℝ) else 0 := by
  simp [eb, EuclideanSpace.inner_single_left]

theorem gb_jet2_symm {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ)
    (x : EuclideanSpace ℝ (Fin n)) (a c : Fin n) : jet 2 φ x ![a, c] = jet 2 φ x ![c, a] :=
  jet_swap_first hφ 0 x a c ![]

/-- `∂_a X = U e_a + J U (Σ_c u_{ac} e_c)`, with the vector written as a `toLp`. -/
theorem gb_dX' {J : V →L[ℝ] V} {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} {y : V}
    {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ) (x : EuclideanSpace ℝ (Fin n))
    (a : Fin n) :
    dX (graphX J U y φ) x a =
      U (eb a) + J (U (WithLp.toLp 2 (fun c => jet 2 φ x ![a, c]))) := by
  rw [dX_graphX hφ, gb_sum_eb]

theorem gb_ddX' {J : V →L[ℝ] V} {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} {y : V}
    {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ) (x : EuclideanSpace ℝ (Fin n))
    (a b : Fin n) :
    ddX (graphX J U y φ) x a b =
      J (U (WithLp.toLp 2 (fun c => jet 3 φ x ![a, b, c]))) := by
  rw [ddX_graphX hφ, gb_sum_eb]

theorem gb_inner_dX {J : V →L[ℝ] V} (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫)
    {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} (hU : IsLagFrame J U) {y : V}
    {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ) (x : EuclideanSpace ℝ (Fin n))
    (a b : Fin n) :
    ⟪dX (graphX J U y φ) x a, dX (graphX J U y φ) x b⟫ =
      (if a = b then 1 else 0) + ∑ c, jet 2 φ x ![a, c] * jet 2 φ x ![b, c] := by
  rw [gb_dX' hφ, gb_dX' hφ, inner_add_left, inner_add_right, inner_add_right, gb_UU,
    gb_UJU hU, gb_JUU hU, gb_JUJU hJi, gb_inner_eb_eb, gb_inner_toLp]
  ring

theorem gX_graphX {J : V →L[ℝ] V} (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫)
    {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} (hU : IsLagFrame J U) {y : V}
    {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ) (x : EuclideanSpace ℝ (Fin n)) :
    gX (graphX J U y φ) x = gMat (jet 2 φ x) := by
  ext a b
  exact gb_inner_dX hJi hU hφ x a b

theorem ginvZ_graphX {J : V →L[ℝ] V} (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫)
    {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} (hU : IsLagFrame J U) {y : V}
    {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ) (x : EuclideanSpace ℝ (Fin n)) :
    ginvZ (graphX J U y φ) x = ginvFun (jet 2 φ x) := by
  unfold ginvZ ginvFun
  rw [gX_graphX hJi hU hφ]

theorem volX_graphX {J : V →L[ℝ] V} (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫)
    {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} (hU : IsLagFrame J U) {y : V}
    {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ) (x : EuclideanSpace ℝ (Fin n)) :
    volX (graphX J U y φ) x = volFun (jet 2 φ x) := by
  unfold volX volFun
  rw [gX_graphX hJi hU hφ]

theorem gb_inner_ddX_dX {J : V →L[ℝ] V} (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫)
    {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} (hU : IsLagFrame J U) {y : V}
    {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ) (x : EuclideanSpace ℝ (Fin n))
    (a b d : Fin n) :
    ⟪ddX (graphX J U y φ) x a b, dX (graphX J U y φ) x d⟫ =
      ∑ f, jet 3 φ x ![a, b, f] * jet 2 φ x ![d, f] := by
  rw [gb_ddX' hφ, gb_dX' hφ, inner_add_right, gb_JUU hU, gb_JUJU hJi, gb_inner_toLp, zero_add]

theorem gamX_graphX {J : V →L[ℝ] V} (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫)
    {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} (hU : IsLagFrame J U) {y : V}
    {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ) (x : EuclideanSpace ℝ (Fin n)) :
    gamX (graphX J U y φ) x = gamFun (jet 2 φ x, jet 3 φ x) := by
  unfold gamX gamFun
  ext I
  simp only [gb_inner_ddX_dX hJi hU hφ, gX_graphX hJi hU hφ, Finset.mul_sum]
  refine Finset.sum_congr rfl fun d _ => Finset.sum_congr rfl fun f _ => ?_
  ring

theorem cX_graphX {J : V →L[ℝ] V} (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫)
    (hJJ : ∀ v : V, J (J v) = -v)
    {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} (hU : IsLagFrame J U) {y : V}
    {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ) (x : EuclideanSpace ℝ (Fin n)) :
    cX J (graphX J U y φ) x = jet 3 φ x := by
  ext I
  show ⟪ddX (graphX J U y φ) x (I 0) (I 1), J (dX (graphX J U y φ) x (I 2))⟫ = _
  rw [gb_ddX' hφ, gb_dX' hφ, map_add, hJJ, inner_add_right, gb_JUJU hJi, inner_neg_right,
    gb_JUU hU, gb_inner_toLp_eb, neg_zero, add_zero]
  have hI : I = ![I 0, I 1, I 2] := by
    funext i; fin_cases i <;> rfl
  rw [hI]
  rfl

theorem gb_fderiv_graphX_apply {J : V →L[ℝ] V} {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} {y : V}
    {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ) (x u : EuclideanSpace ℝ (Fin n)) :
    fderiv ℝ (graphX J U y φ) x u = U u + J (U (fderiv ℝ (gradE φ) x u)) := by
  rw [gb_fderiv_graphX hφ]
  rfl

theorem gb_hess_symm {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ)
    (x u w : EuclideanSpace ℝ (Fin n)) :
    ⟪fderiv ℝ (gradE φ) x u, w⟫ = ⟪u, fderiv ℝ (gradE φ) x w⟫ := by
  have key : ∀ a b : Fin n, ⟪fderiv ℝ (gradE φ) x (eb a), eb b⟫ =
      ⟪eb a, fderiv ℝ (gradE φ) x (eb b)⟫ := by
    intro a b
    rw [gb_gradE_comp hφ, gb_gradE_comp hφ, gb_inner_toLp_eb, real_inner_comm, gb_inner_toLp_eb]
    exact gb_jet2_symm hφ x a b
  have hu : u = ∑ i, u i • eb i := by rw [gb_sum_eb]
  have hw : w = ∑ i, w i • eb i := by rw [gb_sum_eb]
  conv_lhs => rw [hu, hw]
  conv_rhs => rw [hu, hw]
  simp only [map_sum, map_smul, sum_inner, inner_sum, inner_smul_left, inner_smul_right, key]

theorem lagrangian_graphX {J : V →L[ℝ] V} (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫)
    (hJJ : ∀ v : V, J (J v) = -v)
    {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} (hU : IsLagFrame J U) {y : V}
    {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ) (x u w : EuclideanSpace ℝ (Fin n)) :
    ⟪J (fderiv ℝ (graphX J U y φ) x u), fderiv ℝ (graphX J U y φ) x w⟫ = 0 := by
  rw [gb_fderiv_graphX_apply hφ, gb_fderiv_graphX_apply hφ]
  simp only [map_add, hJJ, inner_add_left, inner_add_right, inner_neg_left, gb_JUU hU,
    gb_JUJU hJi, gb_UJU hU, gb_UU]
  rw [gb_hess_symm hφ x u w]
  ring

theorem injective_fderiv_graphX {J : V →L[ℝ] V} (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫)
    (hJJ : ∀ v : V, J (J v) = -v)
    {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} (hU : IsLagFrame J U) {y : V}
    {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ) (x : EuclideanSpace ℝ (Fin n)) :
    Function.Injective (fderiv ℝ (graphX J U y φ) x) := by
  rw [injective_iff_map_eq_zero]
  intro v hv
  have h : ⟪fderiv ℝ (graphX J U y φ) x v, U v⟫ = ⟪v, v⟫ := by
    rw [gb_fderiv_graphX_apply hφ, inner_add_left, gb_UU, gb_JUU hU, add_zero]
  rw [hv, inner_zero_left] at h
  exact inner_self_eq_zero.1 h.symm

/-! ### Part 3: the covariant derivatives and the energy density -/

theorem low23_jp {N : ℕ} (hN : 3 ≤ N) (φ : EuclideanSpace ℝ (Fin n) → ℝ)
    (x : EuclideanSpace ℝ (Fin n)) : low23 (jp N φ x) = (jet 2 φ x, jet 3 φ x) := by
  refine Prod.ext ?_ ?_
  · rw [low23_fst _ (by omega)]; rfl
  · rw [low23_snd _ (by omega)]; rfl

theorem gb_exists_adm {N : ℕ} (φ : EuclideanSpace ℝ (Fin n) → ℝ) (x : EuclideanSpace ℝ (Fin n)) :
    ∃ K : ℝ, jp N φ x ∈ admSet n N K := by
  refine ⟨1 + ∑ m : Fin (N + 1), ‖(jp N φ x).2 m‖, fun m _ => ?_⟩
  have := Finset.single_le_sum (f := fun m : Fin (N + 1) => ‖(jp N φ x).2 m‖)
    (fun m _ => norm_nonneg _) (Finset.mem_univ m)
  linarith

theorem gb_SX_fun {J : V →L[ℝ] V} (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫)
    (hJJ : ∀ v : V, J (J v) = -v)
    {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} (hU : IsLagFrame J U) {y : V}
    {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ) {N : ℕ} (k : ℕ) (hN : k + 3 ≤ N) :
    SX J (graphX J U y φ) k = fun z => tSym k (jp N φ z) := by
  induction k with
  | zero =>
    funext z
    rw [SX_zero, cX_graphX hJi hJJ hU hφ, tSym_zero_jp (by omega)]
  | succ k ih =>
    have ih := ih (by omega)
    funext x
    obtain ⟨K, hK⟩ := gb_exists_adm (N := N) φ x
    have hA := isOpen_admSet (n := n) (N := N) K
    have hlow := admSet_low (n := n) (N := N) K
    have hT := isSymbol_tSym hA hlow k (by omega)
    have hdiff : DifferentiableAt ℝ (fun z => tSym k (jp N φ z)) x := by
      have h1 := contDiffOn_symbol_jp (N := N) hT hφ
      have h2 : IsOpen {y | jp N φ y ∈ admSet n N K} := isOpen_jp_preimage hT hφ
      exact (h1.contDiffAt (h2.mem_nhds hK)).differentiableAt (by simp)
    ext I
    rw [SX_succ, PiLp.sub_apply, curryZ_apply, actG_apply,
      tSym_succ_jp_apply hA hlow k (by omega) hφ hK I, ih,
      fderiv_piLp_apply hdiff, gamX_graphX hJi hU hφ]
    congr 1
    have e1 : gamFun (jet 2 φ x, jet 3 φ x) = gamSym (jp N φ x) := by
      unfold gamSym; rw [low23_jp (by omega)]
    rw [e1]

theorem SX_graphX {J : V →L[ℝ] V} (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫)
    (hJJ : ∀ v : V, J (J v) = -v)
    {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} (hU : IsLagFrame J U) {y : V}
    {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ) {N k : ℕ} (hN : k + 3 ≤ N)
    (x : EuclideanSpace ℝ (Fin n)) :
    SX J (graphX J U y φ) k x = tSym k (jp N φ x) :=
  congrFun (gb_SX_fun hJi hJJ hU hφ k hN) x

theorem eDens_graphX {J : V →L[ℝ] V} (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫)
    (hJJ : ∀ v : V, J (J v) = -v)
    {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} (hU : IsLagFrame J U) {y : V}
    {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ) {N k : ℕ} (hN : k + 3 ≤ N)
    (x : EuclideanSpace ℝ (Fin n)) :
    eDens J (graphX J U y φ) k x = eSym k (jp N φ x) := by
  unfold eDens eSym
  rw [ginvZ_graphX hJi hU hφ, volX_graphX hJi hU hφ, SX_graphX hJi hJJ hU hφ hN,
    low23_jp (by omega)]

theorem aNormSq_graphX {J : V →L[ℝ] V} (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫)
    (hJJ : ∀ v : V, J (J v) = -v)
    {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} (hU : IsLagFrame J U) {y : V}
    {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ) (x : EuclideanSpace ℝ (Fin n)) :
    aNormSq J (graphX J U y φ) x = gNormSq (ginvFun (jet 2 φ x)) (jet 3 φ x) := by
  unfold aNormSq
  rw [ginvZ_graphX hJi hU hφ, cX_graphX hJi hJJ hU hφ]

/-- info: 'CWFlow.Extension.gX_graphX' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms gX_graphX

/-- info: 'CWFlow.Extension.SX_graphX' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms SX_graphX

/-- info: 'CWFlow.Extension.eDens_graphX' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms eDens_graphX

/-- info: 'CWFlow.Extension.lagrangian_graphX' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms lagrangian_graphX

end CWFlow.Extension
