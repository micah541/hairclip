import CWFlow.Extension.EnergySymbol

/-!
# The geometry of an immersion in coordinates

For a map `X : ℝⁿ → V` into a real inner product space with a linear map `J : V → V`, the coordinate
expressions of the tensors used in the extension theorem:

* `dX X x a = ∂_a X`, `ddX X x a b = ∂_a ∂_b X`;
* `gX X x`: the induced metric `g_{ab} = ⟨∂_a X, ∂_b X⟩` as a matrix; `ginvZ`, `volX = √det g`;
* `gamX X x`: the Christoffel symbols `Γ^c_{ab} = g^{cd} ⟨∂_a ∂_b X, ∂_d X⟩`, index order `(a, b, c)`;
* `cX J X x`: the Lagrangian second fundamental form `C_{abc} = ⟨∂_a ∂_b X, J ∂_c X⟩`;
* `SX J X k = ∇^k C`, a tensor with `k + 3` indices, the derivative indices first:
  `∇S = DS − Γ · S` with the operations `curryZ`, `actG` of `TensorOps`;
* `eDens J X k = |∇^k C|²_g √det g`, the energy density, and `aNormSq J X = |C|²_g` (`= |A|²_g` for a
  Lagrangian immersion in a space of dimension `2n`).

These are definitions only. They are the same operations as the symbols `gMat`, `gamFun`, `tSym`, `eSym`
on jet space, which they become for a graph `x ↦ y + U x + J U ∇φ(x)` over a Lagrangian plane.
-/

open scoped ContDiff Matrix RealInnerProductSpace

namespace CWFlow.Extension

variable {n : ℕ} {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]

/-- The `a`-th coordinate vector of `ℝⁿ`. -/
noncomputable abbrev eb (a : Fin n) : EuclideanSpace ℝ (Fin n) := EuclideanSpace.single a (1 : ℝ)

/-- `∂_a X`. -/
noncomputable def dX (X : EuclideanSpace ℝ (Fin n) → V) (x : EuclideanSpace ℝ (Fin n)) (a : Fin n) : V :=
  fderiv ℝ X x (eb a)

/-- `∂_a ∂_b X`. -/
noncomputable def ddX (X : EuclideanSpace ℝ (Fin n) → V) (x : EuclideanSpace ℝ (Fin n)) (a b : Fin n) :
    V :=
  fderiv ℝ (fderiv ℝ X) x (eb a) (eb b)

/-- The induced metric `g_{ab} = ⟨∂_a X, ∂_b X⟩`. -/
noncomputable def gX (X : EuclideanSpace ℝ (Fin n) → V) (x : EuclideanSpace ℝ (Fin n)) :
    Matrix (Fin n) (Fin n) ℝ :=
  fun a b => ⟪dX X x a, dX X x b⟫

/-- The inverse metric as a 2-tensor. -/
noncomputable def ginvZ (X : EuclideanSpace ℝ (Fin n) → V) (x : EuclideanSpace ℝ (Fin n)) : Z n 2 :=
  WithLp.toLp 2 (fun I => (gX X x)⁻¹ (I 0) (I 1))

/-- The volume density `√det g`. -/
noncomputable def volX (X : EuclideanSpace ℝ (Fin n) → V) (x : EuclideanSpace ℝ (Fin n)) : ℝ :=
  Real.sqrt (gX X x).det

/-- The Christoffel symbols `Γ^c_{ab} = g^{cd} ⟨∂_a ∂_b X, ∂_d X⟩`, with index order `(a, b, c)`. -/
noncomputable def gamX (X : EuclideanSpace ℝ (Fin n) → V) (x : EuclideanSpace ℝ (Fin n)) : Z n 3 :=
  WithLp.toLp 2 (fun I => ∑ d, (gX X x)⁻¹ (I 2) d * ⟪ddX X x (I 0) (I 1), dX X x d⟫)

/-- The Lagrangian second fundamental form `C_{abc} = ⟨∂_a ∂_b X, J ∂_c X⟩`. -/
noncomputable def cX (J : V →L[ℝ] V) (X : EuclideanSpace ℝ (Fin n) → V)
    (x : EuclideanSpace ℝ (Fin n)) : Z n 3 :=
  WithLp.toLp 2 (fun I => ⟪ddX X x (I 0) (I 1), J (dX X x (I 2))⟫)

/-- `∇^k C`, with `k + 3` indices, the derivative indices first. -/
noncomputable def SX (J : V →L[ℝ] V) (X : EuclideanSpace ℝ (Fin n) → V) :
    (k : ℕ) → EuclideanSpace ℝ (Fin n) → Z n (k + 3)
  | 0 => cX J X
  | k + 1 => fun x => curryZ (fderiv ℝ (SX J X k) x) - actG (gamX X x) (SX J X k x)

theorem SX_zero (J : V →L[ℝ] V) (X : EuclideanSpace ℝ (Fin n) → V) : SX J X 0 = cX J X := rfl

theorem SX_succ (J : V →L[ℝ] V) (X : EuclideanSpace ℝ (Fin n) → V) (k : ℕ)
    (x : EuclideanSpace ℝ (Fin n)) :
    SX J X (k + 1) x = curryZ (fderiv ℝ (SX J X k) x) - actG (gamX X x) (SX J X k x) := rfl

/-- The energy density `|∇^k C|²_g √det g`. -/
noncomputable def eDens (J : V →L[ℝ] V) (X : EuclideanSpace ℝ (Fin n) → V) (k : ℕ)
    (x : EuclideanSpace ℝ (Fin n)) : ℝ :=
  gNormSq (ginvZ X x) (SX J X k x) * volX X x

/-- `|C|²_g`. -/
noncomputable def aNormSq (J : V →L[ℝ] V) (X : EuclideanSpace ℝ (Fin n) → V)
    (x : EuclideanSpace ℝ (Fin n)) : ℝ :=
  gNormSq (ginvZ X x) (cX J X x)

/-- `div(JH)` in coordinates: `−g^{cd} g^{ab} (∇_d C)_{abc}`. For a Lagrangian immersion this is
`CWFlow.flowPot J X` (`flowPot_eq_potX` in `ImmPot`). -/
noncomputable def potX (J : V →L[ℝ] V) (X : EuclideanSpace ℝ (Fin n) → V)
    (x : EuclideanSpace ℝ (Fin n)) : ℝ :=
  -∑ a, ∑ b, ∑ c, ∑ d, ginvZ X x ![c, d] * ginvZ X x ![a, b] * SX J X 1 x ![d, a, b, c]

/-- The mean curvature form in coordinates, `μ_c = −g^{ab} C_{abc}` (`= ⟨JH, ∂_c X⟩` for a Lagrangian
immersion). -/
noncomputable def muC (J : V →L[ℝ] V) (X : EuclideanSpace ℝ (Fin n) → V)
    (x : EuclideanSpace ℝ (Fin n)) (c : Fin n) : ℝ :=
  -∑ a, ∑ b, (gX X x)⁻¹ a b * cX J X x ![a, b, c]

/-- The velocity `J ∇ div(JH) = g^{ab} ∂_b(potX) J ∂_a X` in coordinates. -/
noncomputable def velX (J : V →L[ℝ] V) (X : EuclideanSpace ℝ (Fin n) → V)
    (x : EuclideanSpace ℝ (Fin n)) : V :=
  ∑ a, ∑ b, ((gX X x)⁻¹ a b * fderiv ℝ (potX J X) x (eb b)) • J (dX X x a)

/-- The Jacobian matrix `A_{ji} = ∂_i τ^j` of a map of `ℝⁿ`. -/
noncomputable def jacM (τ : EuclideanSpace ℝ (Fin n) → EuclideanSpace ℝ (Fin n))
    (z : EuclideanSpace ℝ (Fin n)) : Matrix (Fin n) (Fin n) ℝ :=
  fun j i => fderiv ℝ τ z (eb i) j

/-- The pullback of a tensor by a matrix in every slot:
`(pullZ A T)(I) = Σ_J (Π_l A(J_l, I_l)) T(J)`. -/
noncomputable def pullZ {m : ℕ} (A : Matrix (Fin n) (Fin n) ℝ) (T : Z n m) : Z n m :=
  WithLp.toLp 2 (fun I => ∑ K : Fin m → Fin n, (∏ l, A (K l) (I l)) * T K)

theorem pullZ_apply {m : ℕ} (A : Matrix (Fin n) (Fin n) ℝ) (T : Z n m) (I : Fin m → Fin n) :
    pullZ A T I = ∑ K : Fin m → Fin n, (∏ l, A (K l) (I l)) * T K := rfl

/-- The gradient of a function on `ℝⁿ` as a vector. -/
noncomputable def gradE (φ : EuclideanSpace ℝ (Fin n) → ℝ) (x : EuclideanSpace ℝ (Fin n)) :
    EuclideanSpace ℝ (Fin n) :=
  WithLp.toLp 2 (fun c => fderiv ℝ φ x (eb c))

/-- A Lagrangian frame: an isometric embedding of `ℝⁿ` onto a Lagrangian plane of `V`. -/
structure IsLagFrame (J : V →L[ℝ] V) (U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V) : Prop where
  lag : ∀ a b : EuclideanSpace ℝ (Fin n), ⟪J (U a), U b⟫ = 0

/-- The graph of `∇φ` over the Lagrangian plane `y + U(ℝⁿ)`: `x ↦ y + U x + J U ∇φ(x)`. -/
noncomputable def graphX (J : V →L[ℝ] V) (U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V) (y : V)
    (φ : EuclideanSpace ℝ (Fin n) → ℝ) (x : EuclideanSpace ℝ (Fin n)) : V :=
  y + U x + J (U (gradE φ x))

end CWFlow.Extension
