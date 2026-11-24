import CWFlow.Extension.GraphBridge
import CWFlow.Extension.ImmSmooth
import CWFlow.Interp.BallPotential
import Mathlib.Analysis.Calculus.InverseFunctionTheorem.ContDiff

/-!
# Local graph representation of a Lagrangian immersion

Horizontal and vertical coordinates `horC`, `verC` of a vector of `V` with respect to a Lagrangian frame
`U` (so `v = y₀ + U (horC v) + J U (verC v)`), and the representation of a smooth Lagrangian immersion
whose horizontal projection is a diffeomorphism from an open set onto a ball as the graph of a gradient
(`exists_graph_of_proj`).
-/

open Metric
open scoped ContDiff Topology RealInnerProductSpace

namespace CWFlow.Extension

variable {n : ℕ} {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
  {J : V →L[ℝ] V} {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} {y₀ : V}

/-! ### Part 1: horizontal and vertical coordinates -/

/-- The horizontal coordinate `Uᵀ (v − y₀)`. -/
noncomputable def horC (U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V) (y₀ v : V) :
    EuclideanSpace ℝ (Fin n) :=
  WithLp.toLp 2 (fun a => ⟪U (eb a), v - y₀⟫)

/-- The vertical coordinate `(JU)ᵀ (v − y₀)`. -/
noncomputable def verC (J : V →L[ℝ] V) (U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V) (y₀ v : V) :
    EuclideanSpace ℝ (Fin n) :=
  WithLp.toLp 2 (fun a => ⟪J (U (eb a)), v - y₀⟫)

/-- The linear part of `horC`. -/
noncomputable def horCL (U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V) : V →L[ℝ] EuclideanSpace ℝ (Fin n) :=
  (EuclideanSpace.equiv (Fin n) ℝ).symm.toContinuousLinearMap.comp
    (ContinuousLinearMap.pi fun a => innerSL ℝ (U (eb a)))

/-- The linear part of `verC`. -/
noncomputable def verCL (J : V →L[ℝ] V) (U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V) :
    V →L[ℝ] EuclideanSpace ℝ (Fin n) :=
  (EuclideanSpace.equiv (Fin n) ℝ).symm.toContinuousLinearMap.comp
    (ContinuousLinearMap.pi fun a => innerSL ℝ (J (U (eb a))))

theorem horCL_apply (w : V) : horCL U w = WithLp.toLp 2 (fun a => ⟪U (eb a), w⟫) := rfl

theorem verCL_apply (w : V) : verCL J U w = WithLp.toLp 2 (fun a => ⟪J (U (eb a)), w⟫) := rfl

theorem horC_eq (v : V) : horC U y₀ v = horCL U (v - y₀) := rfl

theorem verC_eq (v : V) : verC J U y₀ v = verCL J U (v - y₀) := rfl

theorem hasFDerivAt_horC (v : V) : HasFDerivAt (horC U y₀) (horCL U) v := by
  have h := (horCL U).hasFDerivAt (x := v - y₀) |>.comp v ((hasFDerivAt_id v).sub_const y₀)
  exact h.congr_fderiv (by simp)

theorem hasFDerivAt_verC (v : V) : HasFDerivAt (verC J U y₀) (verCL J U) v := by
  have h := (verCL J U).hasFDerivAt (x := v - y₀) |>.comp v ((hasFDerivAt_id v).sub_const y₀)
  exact h.congr_fderiv (by simp)

theorem horC_graphX (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    (hU : IsLagFrame J U) (φ : EuclideanSpace ℝ (Fin n) → ℝ) (x : EuclideanSpace ℝ (Fin n)) :
    horC U y₀ (graphX J U y₀ φ x) = x := by
  ext a
  have : ⟪U (eb a), graphX J U y₀ φ x - y₀⟫ = x a := by
    have h : graphX J U y₀ φ x - y₀ = U x + J (U (gradE φ x)) := by
      simp [graphX, add_assoc]
    rw [h, inner_add_right, gb_UU, gb_UJU hU, add_zero]
    simp [eb, EuclideanSpace.inner_single_left]
  exact this

theorem verC_graphX (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    (hU : IsLagFrame J U) (φ : EuclideanSpace ℝ (Fin n) → ℝ) (x : EuclideanSpace ℝ (Fin n)) :
    verC J U y₀ (graphX J U y₀ φ x) = gradE φ x := by
  ext a
  have : ⟪J (U (eb a)), graphX J U y₀ φ x - y₀⟫ = gradE φ x a := by
    have h : graphX J U y₀ φ x - y₀ = U x + J (U (gradE φ x)) := by
      simp [graphX, add_assoc]
    rw [h, inner_add_right, gb_JUU hU, gb_JUJU hJi, zero_add]
    simp [eb, EuclideanSpace.inner_single_left]
  exact this

theorem gl_orthonormal (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hU : IsLagFrame J U) :
    Orthonormal ℝ (Sum.elim (fun a : Fin n => U (eb a)) (fun a : Fin n => J (U (eb a)))) := by
  rw [orthonormal_iff_ite]
  rintro (i | i) (j | j)
  · simp [gb_UU, gb_inner_eb_eb]
  · simp [gb_UJU hU]
  · simp [gb_JUU hU]
  · simp [gb_JUJU hJi, gb_inner_eb_eb]

theorem decomp_horC_verC (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    (hU : IsLagFrame J U) [FiniteDimensional ℝ V] (hdim : Module.finrank ℝ V = 2 * n) (v : V) :
    v = y₀ + U (horC U y₀ v) + J (U (verC J U y₀ v)) := by
  have hon := gl_orthonormal hJi hU
  have hcard : Fintype.card (Fin n ⊕ Fin n) = Module.finrank ℝ V := by
    rw [hdim]; simp [two_mul]
  have hsp : ⊤ ≤ Submodule.span ℝ (Set.range
      (Sum.elim (fun a : Fin n => U (eb a)) (fun a : Fin n => J (U (eb a))))) :=
    (hon.linearIndependent.span_eq_top_of_card_eq_finrank' hcard).ge
  let b := OrthonormalBasis.mk hon hsp
  have h := b.sum_repr' (v - y₀)
  rw [Fintype.sum_sum_type] at h
  simp only [b, OrthonormalBasis.coe_mk, Sum.elim_inl, Sum.elim_inr] at h
  have e1 : U (horC U y₀ v) = ∑ a, ⟪U (eb a), v - y₀⟫ • U (eb a) := by
    have hh : horC U y₀ v = ∑ a, ⟪U (eb a), v - y₀⟫ • eb a := (gb_sum_eb _).symm
    rw [hh]; simp [map_sum, map_smul]
  have e2 : J (U (verC J U y₀ v)) = ∑ a, ⟪J (U (eb a)), v - y₀⟫ • J (U (eb a)) := by
    have hh : verC J U y₀ v = ∑ a, ⟪J (U (eb a)), v - y₀⟫ • eb a := (gb_sum_eb _).symm
    rw [hh]; simp [map_sum, map_smul]
  rw [e1, e2, add_assoc, h]
  abel

/-! ### Part 2: graph representation -/

local notation "𝔼" => EuclideanSpace ℝ (Fin n)

theorem gl_exists_equiv {L : 𝔼 →L[ℝ] 𝔼} (hL : Function.Injective L) :
    ∃ f' : 𝔼 ≃L[ℝ] 𝔼, (f' : 𝔼 →L[ℝ] 𝔼) = L := by
  refine ⟨(LinearEquiv.ofInjectiveEndo (L : 𝔼 →ₗ[ℝ] 𝔼) hL).toContinuousLinearEquiv, ?_⟩
  ext w
  rfl

theorem exists_graph_of_proj (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    (hU : IsLagFrame J U) [FiniteDimensional ℝ V] (hdim : Module.finrank ℝ V = 2 * n)
    {X : 𝔼 → V} {W : Set 𝔼} (hW : IsOpen W)
    (hX : ContDiffOn ℝ ∞ X W)
    (hlag : ∀ y ∈ W, ∀ u w, ⟪J (fderiv ℝ X y u), fderiv ℝ X y w⟫ = 0)
    (hpinj : ∀ y ∈ W, Function.Injective (horCL U ∘L fderiv ℝ X y))
    (hinjOn : Set.InjOn (fun y => horC U y₀ (X y)) W)
    {c : 𝔼} {R : ℝ} (hball : ball c R ⊆ (fun y => horC U y₀ (X y)) '' W) :
    ∃ (s : 𝔼 → 𝔼) (φ : 𝔼 → ℝ),
      ContDiffOn ℝ ∞ s (ball c R) ∧ Set.MapsTo s (ball c R) W ∧
      (∀ x ∈ ball c R, horC U y₀ (X (s x)) = x) ∧
      (∀ y ∈ W, horC U y₀ (X y) ∈ ball c R → s (horC U y₀ (X y)) = y) ∧
      ContDiffOn ℝ ∞ φ (ball c R) ∧ φ c = 0 ∧
      ∀ x ∈ ball c R, X (s x) = graphX J U y₀ φ x := by
  classical
  set χ : 𝔼 → 𝔼 := fun y => horC U y₀ (X y) with hχdef
  have hhor : ContDiff ℝ ∞ (horC U y₀) := by
    have : ContDiff ℝ ∞ (fun v : V => horCL U (v - y₀)) :=
      (horCL U).contDiff.comp (contDiff_id.sub contDiff_const)
    exact this
  have hver : ContDiff ℝ ∞ (verC J U y₀) := by
    have : ContDiff ℝ ∞ (fun v : V => verCL J U (v - y₀)) :=
      (verCL J U).contDiff.comp (contDiff_id.sub contDiff_const)
    exact this
  have hXat : ∀ y ∈ W, ContDiffAt ℝ ∞ X y := fun y hy => hX.contDiffAt (hW.mem_nhds hy)
  have hχat : ∀ y ∈ W, ContDiffAt ℝ ∞ χ y := fun y hy => hhor.contDiffAt.comp y (hXat y hy)
  have hχder : ∀ y ∈ W, HasFDerivAt χ (horCL U ∘L fderiv ℝ X y) y := fun y hy =>
    (hasFDerivAt_horC (X y)).comp y ((hXat y hy).differentiableAt (by simp)).hasFDerivAt
  set s : 𝔼 → 𝔼 := Function.invFunOn χ W with hs
  have hsW : ∀ x ∈ ball c R, s x ∈ W ∧ χ (s x) = x := fun x hx => by
    obtain ⟨y, hy, hyx⟩ := hball hx
    exact ⟨Function.invFunOn_mem ⟨y, hy, hyx⟩, Function.invFunOn_eq ⟨y, hy, hyx⟩⟩
  have hleft : ∀ y ∈ W, s (χ y) = y := fun y hy => hinjOn.leftInvOn_invFunOn hy
  have hsat : ∀ x ∈ ball c R, ContDiffAt ℝ ∞ s x := by
    intro x hx
    obtain ⟨hyW, hyχ⟩ := hsW x hx
    obtain ⟨f', hf'⟩ := gl_exists_equiv (hpinj (s x) hyW)
    have hd : HasFDerivAt χ (f' : 𝔼 →L[ℝ] 𝔼) (s x) := by rw [hf']; exact hχder _ hyW
    have hn : (∞ : WithTop ℕ∞) ≠ 0 := by simp
    have hcd := hχat _ hyW
    have hg : ContDiffAt ℝ ∞ (hcd.localInverse hd hn) (χ (s x)) := hcd.to_localInverse hd hn
    rw [hyχ] at hg
    have hgx : hcd.localInverse hd hn x = s x := by
      have := hcd.localInverse_apply_image hd hn
      rwa [hyχ] at this
    have hev1 : ∀ᶠ x' in 𝓝 x, hcd.localInverse hd hn x' ∈ W :=
      hg.continuousAt.eventually_mem (by rw [hgx]; exact hW.mem_nhds hyW)
    have hev2 : ∀ᶠ x' in 𝓝 x, χ (hcd.localInverse hd hn x') = x' := by
      have := (hcd.hasStrictFDerivAt' hd hn).eventually_right_inverse
      rwa [hyχ] at this
    have hev3 : ∀ᶠ x' in 𝓝 x, x' ∈ ball c R := isOpen_ball.mem_nhds hx
    refine hg.congr_of_eventuallyEq ?_
    filter_upwards [hev1, hev2, hev3] with x' h1 h2 h3
    obtain ⟨h4, h5⟩ := hsW x' h3
    exact hinjOn h4 h1 (h5.trans h2.symm)
  have hsC : ContDiffOn ℝ ∞ s (ball c R) := fun x hx => (hsat x hx).contDiffWithinAt
  have hmaps : Set.MapsTo s (ball c R) W := fun x hx => (hsW x hx).1
  set v : 𝔼 → 𝔼 := fun x => verC J U y₀ (X (s x)) with hv
  have hvC : ContDiffOn ℝ ∞ v (ball c R) := hver.comp_contDiffOn (hX.comp hsC hmaps)
  have hdec : ∀ x ∈ ball c R, X (s x) = y₀ + U x + J (U (v x)) := by
    intro x hx
    have h := decomp_horC_verC (y₀ := y₀) hJi hJJ hU hdim (X (s x))
    have h2 : horC U y₀ (X (s x)) = x := (hsW x hx).2
    rw [h2] at h
    exact h
  have hsym : ∀ x ∈ ball c R, ∀ a b : 𝔼, ⟪fderiv ℝ v x a, b⟫ = ⟪a, fderiv ℝ v x b⟫ := by
    intro x hx a b
    have hxW := (hsW x hx).1
    have hsd : DifferentiableAt ℝ s x := (hsat x hx).differentiableAt (by simp)
    have hvd : DifferentiableAt ℝ v x :=
      (hvC.contDiffAt (isOpen_ball.mem_nhds hx)).differentiableAt (by simp)
    have hXd : DifferentiableAt ℝ X (s x) := (hXat _ hxW).differentiableAt (by simp)
    have h1 : HasFDerivAt (fun z => X (s z)) ((fderiv ℝ X (s x)).comp (fderiv ℝ s x)) x :=
      hXd.hasFDerivAt.comp x hsd.hasFDerivAt
    have h2 : HasFDerivAt (fun z => y₀ + U z + J (U (v z)))
        (U.toContinuousLinearMap + (J.comp U.toContinuousLinearMap).comp (fderiv ℝ v x)) x :=
      ((U.toContinuousLinearMap.hasFDerivAt (x := x)).const_add y₀).add
        ((J.comp U.toContinuousLinearMap).hasFDerivAt.comp x hvd.hasFDerivAt)
    have h3 : HasFDerivAt (fun z => X (s z))
        (U.toContinuousLinearMap + (J.comp U.toContinuousLinearMap).comp (fderiv ℝ v x)) x :=
      h2.congr_of_eventuallyEq (by
        filter_upwards [isOpen_ball.mem_nhds hx] with z hz using hdec z hz)
    have h4 := h1.unique h3
    have hl := hlag (s x) hxW (fderiv ℝ s x a) (fderiv ℝ s x b)
    have e : ∀ w, fderiv ℝ X (s x) (fderiv ℝ s x w) = U w + J (U (fderiv ℝ v x w)) := fun w => by
      have := congrArg (fun L => L w) h4
      simpa using this
    rw [e, e] at hl
    simp only [map_add, hJJ, inner_add_left, inner_add_right, inner_neg_left, gb_JUU hU,
      gb_JUJU hJi, gb_UJU hU, gb_UU, neg_zero, add_zero, zero_add] at hl
    linarith
  obtain ⟨φ, hφC, hφ0, hφd⟩ := CWFlow.Interp.exists_potential_ball hvC hsym
  have hgrad : ∀ x ∈ ball c R, gradE φ x = v x := by
    intro x hx
    ext b
    show fderiv ℝ φ x (eb b) = v x b
    rw [hφd x hx]
    simp [eb, EuclideanSpace.inner_single_right]
  refine ⟨s, φ, hsC, hmaps, fun x hx => (hsW x hx).2, fun y hy _ => hleft y hy, hφC, hφ0, ?_⟩
  intro x hx
  rw [hdec x hx]
  simp [graphX, hgrad x hx]

/-- info: 'CWFlow.Extension.decomp_horC_verC' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms decomp_horC_verC

/-- info: 'CWFlow.Extension.exists_graph_of_proj' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_graph_of_proj

end CWFlow.Extension
