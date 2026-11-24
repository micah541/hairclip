import CWFlow.Weinstein.Interface

/-!
# Graphs of exact forms under a Weinstein map (Chen–Warren, Prop. 3.2, the pointwise part)

Let `Ψ : E × E* → V` pull `ω₀` back to `ω_can` on an open set `U` (as for a local Weinstein map).
For `u : E → ℝ`, put `F_u(y) = Ψ(y, Du(y))` (`graphFu`).

* `F_u` is Lagrangian (`graphFu_lagrangian`): the graph of the closed form `du` is `ω_can`-isotropic,
  `ω_can((a, D²u a), (b, D²u b)) = D²u(a, b) − D²u(b, a) = 0`.
* `F_u` is an immersion (`injective_fderiv_graphFu`): `DΨ` is injective because `ω_can` is
  nondegenerate (`injective_fderiv_of_symplectic`).
* **The Hamiltonian potential.** For a time-dependent `u(t, y)`, the velocity
  `∂ₜF_{u(t)} = DΨ(0, D∂ₜu)` (`hasDerivAt_graphFu_time`) is vertical, and
  `⟨J ∂ₜF, DF_u a⟩ = D(∂ₜu)(a)` (`inner_J_time_fderiv_graphFu`). So the normal velocity of `F_{u(t)}`
  is `−J∇_g(∂ₜu)`, and `F_{u(t)}` moves by the flow `∂ₜF = J∇div(JH)` (normal part) iff
  `∂ₜu = −div(JH)[F_u]` up to a function of `t`. This is the scalar equation.

These facts use only `Ψ*ω₀ = ω_can`, so they hold for any target (almost Kähler included) once stated
with `pullback2`; here they are stated for `ω₀`.
-/

open Filter Set

open scoped InnerProductSpace RealInnerProductSpace Topology ContDiff

noncomputable section

namespace CWFlow

variable {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V] {J : V →L[ℝ] V}
variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]

/-- `y ↦ (y, Du(y))`, the graph of `du`. -/
def graphMap (u : E → ℝ) (y : E) : E × (E →L[ℝ] ℝ) := (y, fderiv ℝ u y)

/-- `F_u = Ψ ∘ (id, Du)`. -/
def graphFu (Ψ : E × (E →L[ℝ] ℝ) → V) (u : E → ℝ) : E → V := Ψ ∘ graphMap u

theorem hasFDerivAt_graphMap {u : E → ℝ} {y : E} (hu : ContDiffAt ℝ 2 u y) :
    HasFDerivAt (graphMap u)
      ((ContinuousLinearMap.id ℝ E).prod (fderiv ℝ (fderiv ℝ u) y)) y := by
  have hD := (hu.fderiv_right (m := 1) (by norm_num)).differentiableAt (by norm_num)
  exact (hasFDerivAt_id y).prodMk hD.hasFDerivAt

/-- `DF_u(a) = DΨ(a, D²u(a))`. -/
theorem fderiv_graphFu_apply {Ψ : E × (E →L[ℝ] ℝ) → V} {u : E → ℝ} {y : E}
    (hΨ : DifferentiableAt ℝ Ψ (graphMap u y)) (hu : ContDiffAt ℝ 2 u y) (a : E) :
    fderiv ℝ (graphFu Ψ u) y a = fderiv ℝ Ψ (graphMap u y) (a, fderiv ℝ (fderiv ℝ u) y a) := by
  rw [graphFu, fderiv_comp y hΨ (hasFDerivAt_graphMap hu).differentiableAt,
    (hasFDerivAt_graphMap hu).fderiv]
  rfl

variable {Ψ : E × (E →L[ℝ] ℝ) → V} {q : E × (E →L[ℝ] ℝ)}

/-- A map pulling `ω₀` back to the nondegenerate `ω_can` has injective derivative. -/
theorem injective_fderiv_of_symplectic
    (hs : ∀ a b, ⟪J (fderiv ℝ Ψ q a), fderiv ℝ Ψ q b⟫ = omegaCanVal a b) :
    Function.Injective (fderiv ℝ Ψ q) := by
  rw [injective_iff_map_eq_zero]
  intro v hv
  have h : ∀ w, omegaCanVal v w = 0 := fun w => by rw [← hs, hv, map_zero, inner_zero_left]
  have h2 : v.2 = 0 := by
    ext e
    simpa [omegaCanVal] using h (e, 0)
  have h1 : v.1 = 0 := by
    refine SeparatingDual.eq_zero_of_forall_dual_eq_zero (R := ℝ) fun f => ?_
    simpa [omegaCanVal, h2] using h (0, f)
  exact Prod.ext h1 h2

/-- The vertical pairing: `⟨J DΨ(0, ζ), DΨ(w)⟩ = ζ(w₁)`. -/
theorem inner_J_fderiv_vertical
    (hs : ∀ a b, ⟪J (fderiv ℝ Ψ q a), fderiv ℝ Ψ q b⟫ = omegaCanVal a b)
    (ζ : E →L[ℝ] ℝ) (w : E × (E →L[ℝ] ℝ)) :
    ⟪J (fderiv ℝ Ψ q (0, ζ)), fderiv ℝ Ψ q w⟫ = ζ w.1 := by
  rw [hs]
  simp [omegaCanVal]

variable {u : E → ℝ} {y : E}

/-- `F_u` is Lagrangian. -/
theorem graphFu_lagrangian (hΨ : DifferentiableAt ℝ Ψ (graphMap u y)) (hu : ContDiffAt ℝ 2 u y)
    (hs : ∀ a b, ⟪J (fderiv ℝ Ψ (graphMap u y) a), fderiv ℝ Ψ (graphMap u y) b⟫ = omegaCanVal a b)
    (a b : E) : ⟪J (fderiv ℝ (graphFu Ψ u) y a), fderiv ℝ (graphFu Ψ u) y b⟫ = 0 := by
  rw [fderiv_graphFu_apply hΨ hu, fderiv_graphFu_apply hΨ hu, hs]
  simp only [omegaCanVal]
  rw [hu.isSymmSndFDerivAt (by simp) a b, sub_self]

/-- `F_u` is an immersion. -/
theorem injective_fderiv_graphFu (hΨ : DifferentiableAt ℝ Ψ (graphMap u y))
    (hu : ContDiffAt ℝ 2 u y)
    (hs : ∀ a b, ⟪J (fderiv ℝ Ψ (graphMap u y) a), fderiv ℝ Ψ (graphMap u y) b⟫ = omegaCanVal a b) :
    Function.Injective (fderiv ℝ (graphFu Ψ u) y) := by
  intro a b hab
  rw [funext (fderiv_graphFu_apply hΨ hu)] at hab
  have := injective_fderiv_of_symplectic hs hab
  exact congrArg Prod.fst this

/-! ### Time-dependent graphs: the Hamiltonian potential -/

/-- `∂ₜ D_y w(t, ·) = D_y(∂ₜ w)`, written as `D²w(t, y)(1, 0) ∘ inr`. -/
theorem hasDerivAt_slice_fderiv {w : ℝ × E → ℝ} (hw : ContDiff ℝ 2 w) (t : ℝ) (y : E) :
    HasDerivAt (fun s => fderiv ℝ (fun z => w (s, z)) y)
      ((fderiv ℝ (fderiv ℝ w) (t, y) (1, 0)).comp (ContinuousLinearMap.inr ℝ ℝ E)) t := by
  have hwd : ∀ p, DifferentiableAt ℝ w p := fun p => hw.differentiable (by simp) p
  have hD := (hw.fderiv_right (m := 1) (by norm_num)).differentiable one_ne_zero
  have hcurve := (hD (t, y)).hasFDerivAt.comp_hasDerivAt t
    ((hasDerivAt_id t).prodMk (hasDerivAt_const t y))
  have hslice : ∀ s, fderiv ℝ (fun z => w (s, z)) y =
      (fderiv ℝ w (s, y)).comp (ContinuousLinearMap.inr ℝ ℝ E) := fun s => by
    have h : HasFDerivAt (fun z => w (s, z))
        ((fderiv ℝ w (s, y)).comp (ContinuousLinearMap.inr ℝ ℝ E)) y :=
      (hwd (s, y)).hasFDerivAt.comp y (hasFDerivAt_prodMk_right s y)
    exact h.fderiv
  simp only [hslice]
  have hpre := ((ContinuousLinearMap.compL ℝ E (ℝ × E) ℝ).flip
    (ContinuousLinearMap.inr ℝ ℝ E)).hasFDerivAt.comp_hasDerivAt t hcurve
  refine hpre.congr_deriv (ContinuousLinearMap.ext fun a => ?_)
  simp

/-- The spatial derivative of `∂ₜw` is `D²w(t, y)(1, 0) ∘ inr`. -/
theorem fderiv_time_deriv_apply {w : ℝ × E → ℝ} (hw : ContDiff ℝ 2 w) (t : ℝ) (y a : E) :
    fderiv ℝ (fun z => fderiv ℝ w (t, z) (1, 0)) y a = fderiv ℝ (fderiv ℝ w) (t, y) (1, 0) (0, a) := by
  have hD := (hw.fderiv_right (m := 1) (by norm_num)).differentiable one_ne_zero
  have h : HasFDerivAt (fun z => fderiv ℝ w (t, z))
      ((fderiv ℝ (fderiv ℝ w) (t, y)).comp (ContinuousLinearMap.inr ℝ ℝ E)) y :=
    (hD (t, y)).hasFDerivAt.comp y (hasFDerivAt_prodMk_right t y)
  have h2 := h.clm_apply (hasFDerivAt_const ((1 : ℝ), (0 : E)) y)
  rw [h2.fderiv]
  simp only [ContinuousLinearMap.add_apply, ContinuousLinearMap.coe_comp', Function.comp_apply,
    ContinuousLinearMap.inr_apply, ContinuousLinearMap.zero_apply, map_zero,
    ContinuousLinearMap.flip_apply, zero_add]
  have hsym := (hw.contDiffAt (x := (t, y))).isSymmSndFDerivAt (n := 2) (by simp)
    ((0 : ℝ), a) ((1 : ℝ), (0 : E))
  simp [hsym]

/-- The velocity of a time-dependent graph is vertical: `∂ₜF_{u(t)} = DΨ(0, D_y ∂ₜw)`. -/
theorem hasDerivAt_graphFu_time {w : ℝ × E → ℝ} (hw : ContDiff ℝ 2 w) {t : ℝ} {y : E}
    (hΨ : DifferentiableAt ℝ Ψ (graphMap (fun z => w (t, z)) y)) :
    HasDerivAt (fun s => graphFu Ψ (fun z => w (s, z)) y)
      (fderiv ℝ Ψ (graphMap (fun z => w (t, z)) y)
        (0, (fderiv ℝ (fderiv ℝ w) (t, y) (1, 0)).comp (ContinuousLinearMap.inr ℝ ℝ E))) t := by
  have hg : HasDerivAt (fun s => graphMap (fun z => w (s, z)) y)
      ((0 : E), (fderiv ℝ (fderiv ℝ w) (t, y) (1, 0)).comp (ContinuousLinearMap.inr ℝ ℝ E)) t :=
    (hasDerivAt_const t y).prodMk (hasDerivAt_slice_fderiv hw t y)
  exact hΨ.hasFDerivAt.comp_hasDerivAt t hg

/-- **The Hamiltonian potential of a time-dependent graph is `∂ₜu`:**
`⟨J ∂ₜF_{u(t)}, DF_{u(t)} a⟩ = D(∂ₜu)(a)`. -/
theorem inner_J_time_fderiv_graphFu {w : ℝ × E → ℝ} (hw : ContDiff ℝ 2 w) {t : ℝ} {y : E}
    (hΨ : DifferentiableAt ℝ Ψ (graphMap (fun z => w (t, z)) y))
    (hs : ∀ a b, ⟪J (fderiv ℝ Ψ (graphMap (fun z => w (t, z)) y) a),
      fderiv ℝ Ψ (graphMap (fun z => w (t, z)) y) b⟫ = omegaCanVal a b) (a : E) :
    ⟪J (fderiv ℝ Ψ (graphMap (fun z => w (t, z)) y)
        (0, (fderiv ℝ (fderiv ℝ w) (t, y) (1, 0)).comp (ContinuousLinearMap.inr ℝ ℝ E))),
      fderiv ℝ (graphFu Ψ (fun z => w (t, z))) y a⟫ =
      fderiv ℝ (fun z => fderiv ℝ w (t, z) (1, 0)) y a := by
  have hu : ContDiffAt ℝ 2 (fun z => w (t, z)) y :=
    (hw.comp (contDiff_const.prodMk contDiff_id)).contDiffAt
  rw [fderiv_graphFu_apply hΨ hu, inner_J_fderiv_vertical hs, fderiv_time_deriv_apply hw]
  rfl

end CWFlow
