import CWFlow.Extension.GraphRadius
import CWFlow.Extension.JetOfSymbol

/-!
# Graph data at every point of a Lagrangian immersion

* `exists_lagFrame`: every Lagrangian subspace given as the image of an injective linear map has an orthonormal
  Lagrangian frame `U` (an isometric embedding of `ℝⁿ` with the same range) with `verCL J U ∘ L = 0`.
* `graph_at`: the graph radius theorem applies at every point `q`, with the frame produced by the tangent plane,
  and its conclusion is strengthened by `‖gradE φ x‖ ≤ c₀ ‖x‖`, `‖verC‖ ≤ c₀ ‖x‖` and a Lipschitz bound for `F ∘ σ`.
-/

open Metric Set
open CWFlow.Interp
open scoped ContDiff Manifold Topology RealInnerProductSpace

set_option linter.unusedSectionVars false

namespace CWFlow.Extension

/-! ### Part 1: Lagrangian frames -/

theorem exists_lagFrame {n : ℕ} {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
    [FiniteDimensional ℝ V] {J : V →L[ℝ] V} (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫)
    (hJJ : ∀ v : V, J (J v) = -v)
    {L : EuclideanSpace ℝ (Fin n) →L[ℝ] V} (hL : Function.Injective L)
    (hlag : ∀ u w, ⟪J (L u), L w⟫ = 0) :
    ∃ U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V, IsLagFrame J U ∧
      (∀ w, ∃ a, U a = L w) ∧ (∀ a, ∃ w, L w = U a) ∧ verCL J U ∘L L = 0 := by
  classical
  set P : Submodule ℝ V := LinearMap.range (L : EuclideanSpace ℝ (Fin n) →ₗ[ℝ] V) with hP
  have hfin : Module.finrank ℝ P = n := by
    rw [hP, LinearMap.finrank_range_of_inj hL, finrank_euclideanSpace_fin]
  let b : OrthonormalBasis (Fin n) ℝ P :=
    (stdOrthonormalBasis ℝ P).reindex (finCongr hfin)
  let e : EuclideanSpace ℝ (Fin n) ≃ₗᵢ[ℝ] P := b.repr.symm
  let U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V := P.subtypeₗᵢ.comp e.toLinearIsometry
  have hUa : ∀ a, U a = (e a : V) := fun a => rfl
  have h1 : ∀ a, ∃ w, L w = U a := fun a => by
    obtain ⟨w, hw⟩ := (e a).2
    exact ⟨w, by rw [hUa]; exact hw⟩
  have h2 : ∀ w, ∃ a, U a = L w := fun w => by
    refine ⟨e.symm ⟨L w, ⟨w, rfl⟩⟩, ?_⟩
    rw [hUa]; simp
  have hlag' : ∀ a b, ⟪J (U a), U b⟫ = 0 := fun a b => by
    obtain ⟨w, hw⟩ := h1 a
    obtain ⟨w', hw'⟩ := h1 b
    rw [← hw, ← hw']; exact hlag w w'
  refine ⟨U, ⟨hlag'⟩, h2, h1, ?_⟩
  ext w : 1
  rw [ContinuousLinearMap.comp_apply, verCL_apply]
  ext a
  obtain ⟨w', hw'⟩ := h1 (eb a)
  simp only [WithLp.ofLp_zero, Pi.zero_apply]
  show ⟪J (U (eb a)), L w⟫ = 0
  rw [← hw']; exact hlag w' w

/-! ### Part 2: graph data at every point -/

variable {n : ℕ} {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
  {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
  {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
  (ι : EuclideanSpace ℝ (Fin n) ≃L[ℝ] E)
  {J : V →L[ℝ] V} {F : M → V}

/-- `gradE φ` is differentiable at a point of smoothness, with derivative bounded by the Hessian. -/
theorem ga_gradE_at {φ : EuclideanSpace ℝ (Fin n) → ℝ} {x : EuclideanSpace ℝ (Fin n)}
    (hφ : ContDiffAt ℝ ∞ φ x) :
    DifferentiableAt ℝ (gradE φ) x ∧ ∀ v, ‖fderiv ℝ (gradE φ) x v‖ ≤ ‖jet 2 φ x‖ * ‖v‖ := by
  have hdf : DifferentiableAt ℝ (fderiv ℝ φ) x :=
    (hφ.fderiv_right (m := ∞) (by simp)).differentiableAt (by simp)
  have hdg : DifferentiableAt ℝ (gradE φ) x := by
    refine (differentiableAt_piLp 2).2 fun c => ?_
    exact hdf.clm_apply (differentiableAt_const _)
  refine ⟨hdg, fun v => ?_⟩
  have hcomp : ∀ a, fderiv ℝ (gradE φ) x (eb a) =
      WithLp.toLp 2 (fun c => jet 2 φ x ![a, c]) := by
    intro a
    ext c
    rw [show (fderiv ℝ (gradE φ) x (eb a)).ofLp c = fderiv ℝ (fun z => gradE φ z c) x (eb a) from
      fderiv_piLp_apply hdg (eb a) c]
    have h1 : jet 2 φ x ![a, c] = fderiv ℝ (fun z => jet 1 φ z ![c]) x (eb a) :=
      jet_cons_eq_fderiv_of_contDiffAt hφ 1 a ![c]
    rw [h1]
    congr 1
    congr 1
    funext z
    rw [gb_jet_one]
    rfl
  have h1 : fderiv ℝ (gradE φ) x v =
      WithLp.toLp 2 (fun b => ∑ a, v a * jet 2 φ x ![a, b]) := by
    rw [ir_lin_expand (fderiv ℝ (gradE φ) x) v]
    ext b
    simp [hcomp]
  rw [h1]
  exact sl_frob (jet 2 φ x) v

theorem ga_lip {φ : EuclideanSpace ℝ (Fin n) → ℝ} {c₀ : ℝ} (hc₀ : 0 < c₀)
    (hφ : ContDiffOn ℝ ∞ φ (ball 0 6))
    (hj : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6, ‖jet 2 φ x‖ ≤ c₀) :
    ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6, ∀ x' ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6,
      ‖gradE φ x - gradE φ x'‖ ≤ c₀ * ‖x - x'‖ := by
  intro x hx x' hx'
  have hat : ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6, ContDiffAt ℝ ∞ φ y := fun y hy =>
    hφ.contDiffAt (isOpen_ball.mem_nhds hy)
  exact Convex.norm_image_sub_le_of_norm_fderiv_le (f := gradE φ) (C := c₀)
    (fun y hy => (ga_gradE_at (hat y hy)).1)
    (fun y hy => ContinuousLinearMap.opNorm_le_bound _ hc₀.le fun v =>
      ((ga_gradE_at (hat y hy)).2 v).trans
        (mul_le_mul_of_nonneg_right (hj y hy) (norm_nonneg _)))
    (convex_ball _ _) hx' hx

theorem graph_at [T2Space M] [CompactSpace M]
    (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    [FiniteDimensional ℝ V] (hdim : Module.finrank ℝ V = 2 * n) {F : M → V} (hF : IsLagImm ι J F)
    {c₀ c₁ : ℝ} (hc₀ : 0 < c₀) (hc₀1 : 6 * c₀ < 1) (hc₁ : 0 ≤ c₁)
    (hsmall : 6 * (c₁ * (1 + c₀ ^ 2) ^ ((3 : ℝ) / 2)) < c₀)
    (hA : ∀ z : M, ∀ y ∈ chartDomE ι z, aNormSq J (F ∘ chartInv ι z) y ≤ c₁ ^ 2) (q : M) :
    ∃ (U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V) (σ : EuclideanSpace ℝ (Fin n) → M)
      (φ : EuclideanSpace ℝ (Fin n) → ℝ),
      IsLagFrame J U ∧ IsParam ι (ball 0 6) σ ∧ σ 0 = q ∧ ContDiffOn ℝ ∞ φ (ball 0 6) ∧ φ 0 = 0 ∧
      (∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6,
        F (σ x) = graphX J U (F q) φ x ∧ horC U (F q) (F (σ x)) = x ∧
        ‖jet 2 φ x‖ ≤ c₀ ∧ ‖gradE φ x‖ ≤ c₀ * ‖x‖ ∧ ‖verC J U (F q) (F (σ x))‖ ≤ c₀ * ‖x‖) ∧
      (∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6, ∀ x' ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6,
        ‖F (σ x) - F (σ x')‖ ≤ (1 + c₀) * ‖x - x'‖) ∧
      ∀ ρ, 0 < ρ → ρ ≤ 6 →
        σ '' ball 0 ρ = connectedComponentIn
          {p | ‖horC U (F q) (F p)‖ < ρ ∧ ‖verC J U (F q) (F p)‖ < 1} q := by
  classical
  have hmem := im_chartFwd_mem ι (mem_extChartAt_source (I := 𝓘(ℝ, E)) q)
  obtain ⟨U, hU, -, -, hv⟩ := exists_lagFrame hJi hJJ (grf_dFq_inj ι hF q)
    (fun u w => hF.lag q _ hmem u w)
  obtain ⟨σ, φ, hP, hσ0, hφ, hφ0, hb, hcomp⟩ :=
    graph_radius ι hJi hJJ hU hdim hF hc₀ hc₀1 hc₁ hsmall hA (p := q) (y₀ := F q) rfl hv
  have hLip := ga_lip hc₀ hφ (fun x hx => (hb x hx).2.2.1)
  have hg0 : gradE φ 0 = 0 := by
    have h := verC_graphX (y₀ := F q) hJi hJJ hU φ 0
    rw [← (hb 0 (mem_ball_self (by norm_num))).1, hσ0, verC_eq, sub_self, map_zero] at h
    exact h.symm
  have hgrad : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6, ‖gradE φ x‖ ≤ c₀ * ‖x‖ := by
    intro x hx
    have := hLip x hx 0 (mem_ball_self (by norm_num))
    simpa [hg0] using this
  refine ⟨U, σ, φ, hU, hP, hσ0, hφ, hφ0, fun x hx => ?_, ?_, hcomp⟩
  · obtain ⟨h1, h2, h3, -⟩ := hb x hx
    refine ⟨h1, h2, h3, hgrad x hx, ?_⟩
    rw [h1, verC_graphX hJi hJJ hU]
    exact hgrad x hx
  · intro x hx x' hx'
    have h1 := (hb x hx).1
    have h2 := (hb x' hx').1
    have hd : F (σ x) - F (σ x') = U (x - x') + J (U (gradE φ x - gradE φ x')) := by
      rw [h1, h2]; simp only [graphX, map_sub]; abel
    have hJn : ∀ v : V, ‖J v‖ = ‖v‖ := fun v => by
      rw [norm_eq_sqrt_real_inner, norm_eq_sqrt_real_inner, hJi]
    rw [hd]
    calc ‖U (x - x') + J (U (gradE φ x - gradE φ x'))‖
        ≤ ‖U (x - x')‖ + ‖J (U (gradE φ x - gradE φ x'))‖ := norm_add_le _ _
      _ = ‖x - x'‖ + ‖gradE φ x - gradE φ x'‖ := by rw [hJn, U.norm_map, U.norm_map]
      _ ≤ ‖x - x'‖ + c₀ * ‖x - x'‖ := by gcongr; exact hLip x hx x' hx'
      _ = (1 + c₀) * ‖x - x'‖ := by ring

/-- info: 'CWFlow.Extension.exists_lagFrame' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_lagFrame

/-- info: 'CWFlow.Extension.graph_at' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms graph_at

end CWFlow.Extension
