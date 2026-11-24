import CWFlow.Extension.TimeGraphUnif
import CWFlow.Extension.ChartGraphEvol
import CWFlow.Extension.ChartBounds
import CWFlow.Extension.VelBound

/-!
# The flow in a graph chart, on the whole persistence interval

A flow `Φ` of Lagrangian immersions (`IsFlowM`) together with a graph chart `(p, U, σ, φ₀)` of `Φ t₀` whose
horizontal coordinate `χ_t = horC ∘ Φ t ∘ σ` stays `C¹`-close to the identity and whose slope stays small on
`(t₀, T)` (`hpers`, the conclusion of `exists_persist`) is described on all of `(t₀, T)` by one graph function
`φ(t, ·)` over the ball of radius `5` (`exists_piece_flow`):

* `∂ₜφ = potSym (jp 4 φ_t)`;
* the bounds `|∇φ| ≤ 1/2`, `|D²φ| ≤ 2c₀`, `|D³φ| ≤ 2c₁`;
* `σ ∘ ξ_t` is a parametrization (`IsParam`) of the piece of `M` over the ball of radius `5`;
* the horizontal coordinate obeys `∂ₜ χ = vSym (jp 5 φ_t ∘ χ_t)`.

The proof applies `exists_time_graph_unif` to `Y (t, y) = Φ t (σ y)` and then converts, slice by slice
(`PfSlice`), the potential and the velocity of `Y_t` to the jet expressions of `φ_t`, as in `cge_core`.
-/

open Metric Set
open CWFlow.Interp
open scoped ContDiff Manifold Topology RealInnerProductSpace

namespace CWFlow.Extension

section Slice

variable {n : ℕ} {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
  {J : V →L[ℝ] V} {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} {y₀ : V}
  {X : EuclideanSpace ℝ (Fin n) → V} {ξ : EuclideanSpace ℝ (Fin n) → EuclideanSpace ℝ (Fin n)}
  {φ : EuclideanSpace ℝ (Fin n) → ℝ}

local notation "𝔼" => EuclideanSpace ℝ (Fin n)

/-- The data of one time slice: a Lagrangian immersion `X` of the ball of radius `6`, whose horizontal map is
`C¹`-close to the identity on the ball of radius `11/2`, with inverse `ξ` on the ball of radius `5` and graph
function `φ` there. -/
structure PfSlice (J : V →L[ℝ] V) (U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V) (y₀ : V)
    (X : EuclideanSpace ℝ (Fin n) → V) (ξ : EuclideanSpace ℝ (Fin n) → EuclideanSpace ℝ (Fin n))
    (φ : EuclideanSpace ℝ (Fin n) → ℝ) : Prop where
  hX : ContDiffOn ℝ ∞ X (ball (0 : EuclideanSpace ℝ (Fin n)) 6)
  hinj : ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6, Function.Injective (fderiv ℝ X y)
  hlag : ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6, ∀ u w, ⟪J (fderiv ℝ X y u), fderiv ℝ X y w⟫ = 0
  hξ : ContDiffOn ℝ ∞ ξ (ball (0 : EuclideanSpace ℝ (Fin n)) 5)
  hφ : ContDiffOn ℝ ∞ φ (ball (0 : EuclideanSpace ℝ (Fin n)) 5)
  pos : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 5,
    ξ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) (11 / 2) ∧ horC U y₀ (X (ξ x)) = x ∧
      X (ξ x) = graphX J U y₀ φ x
  left : ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) (11 / 2),
    horC U y₀ (X y) ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 5 → ξ (horC U y₀ (X y)) = y
  bd : ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) (11 / 2),
    ‖fderiv ℝ (fun y' => horC U y₀ (X y')) y - ContinuousLinearMap.id ℝ _‖ ≤ 1 / 8

theorem pf_inj_of_close {L : 𝔼 →L[ℝ] 𝔼} (h : ‖L - ContinuousLinearMap.id ℝ 𝔼‖ ≤ 1 / 8) :
    Function.Injective L := by
  rw [injective_iff_map_eq_zero]
  intro v hv
  have h1 : ‖v‖ ≤ ‖(L - ContinuousLinearMap.id ℝ 𝔼) v‖ := by
    have : (L - ContinuousLinearMap.id ℝ 𝔼) v = -v := by simp [hv]
    rw [this, norm_neg]
  have h2 := (L - ContinuousLinearMap.id ℝ 𝔼).le_opNorm v
  have h3 : ‖v‖ ≤ 1 / 8 * ‖v‖ := h1.trans (h2.trans (by gcongr))
  have : ‖v‖ = 0 := by linarith [norm_nonneg v]
  exact norm_eq_zero.1 this

theorem pf_ball56 : ball (0 : 𝔼) (11 / 2) ⊆ ball (0 : 𝔼) 6 := ball_subset_ball (by norm_num)

theorem pf_chi_diff (h : PfSlice J U y₀ X ξ φ) {y : 𝔼} (hy : y ∈ ball (0 : 𝔼) 6) :
    DifferentiableAt ℝ (fun y' => horC U y₀ (X y')) y := by
  have hYd : DifferentiableAt ℝ X y := ir_dX isOpen_ball h.hX hy
  exact DifferentiableAt.comp (g := horC U y₀) (f := X) y
    (hasFDerivAt_horC (U := U) (y₀ := y₀) (X y)).differentiableAt hYd

theorem pf_fderiv_chi (h : PfSlice J U y₀ X ξ φ) {y : 𝔼} (hy : y ∈ ball (0 : 𝔼) 6) :
    fderiv ℝ (fun y' => horC U y₀ (X y')) y = (horCL U).comp (fderiv ℝ X y) := by
  have hYd : DifferentiableAt ℝ X y := ir_dX isOpen_ball h.hX hy
  have := (hasFDerivAt_horC (U := U) (y₀ := y₀) (X y)).comp y hYd.hasFDerivAt
  exact this.fderiv

theorem pf_ξ_det (h : PfSlice J U y₀ X ξ φ) {x : 𝔼} (hx : x ∈ ball (0 : 𝔼) 5) :
    (jacM ξ x).det ≠ 0 := by
  have hτd : DifferentiableAt ℝ ξ x := ir_dτ isOpen_ball h.hξ hx
  have hmem := (h.pos x hx).1
  refine im_det_ne_zero_of_left_inv hτd (pf_chi_diff h (pf_ball56 hmem)) ?_
  filter_upwards [isOpen_ball.mem_nhds hx] with w hw
  exact (h.pos w hw).2.1

theorem pf_chi_deriv (h : PfSlice J U y₀ X ξ φ) {x : 𝔼} (hx : x ∈ ball (0 : 𝔼) 5) :
    (fderiv ℝ (fun y' => horC U y₀ (X y')) (ξ x)).comp (fderiv ℝ ξ x) =
      ContinuousLinearMap.id ℝ 𝔼 := by
  have hτd : DifferentiableAt ℝ ξ x := ir_dτ isOpen_ball h.hξ hx
  have hmem := (h.pos x hx).1
  have h1 := fderiv_comp x (pf_chi_diff h (pf_ball56 hmem)) hτd
  have hev : ((fun y' => horC U y₀ (X y')) ∘ ξ) =ᶠ[𝓝 x] id := by
    filter_upwards [isOpen_ball.mem_nhds hx] with w hw
    exact (h.pos w hw).2.1
  rw [← h1, hev.fderiv_eq]
  exact fderiv_id

theorem pf_graph_germ (h : PfSlice J U y₀ X ξ φ) {x : 𝔼} (hx : x ∈ ball (0 : 𝔼) 5) {ψ : 𝔼 → ℝ}
    (hev : φ =ᶠ[𝓝 x] ψ) : graphX J U y₀ ψ =ᶠ[𝓝 x] X ∘ ξ := by
  have hev' : ∀ᶠ y in 𝓝 x, φ =ᶠ[𝓝 y] ψ := hev.eventually_nhds
  filter_upwards [hev', isOpen_ball.mem_nhds hx] with y hy hy5
  rw [Function.comp_apply, (h.pos y hy5).2.2]
  exact (cge_graphX_congr hy).symm

theorem pf_pot (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    (hU : IsLagFrame J U) (h : PfSlice J U y₀ X ξ φ) {x : 𝔼} (hx : x ∈ ball (0 : 𝔼) 5) :
    potX J X (ξ x) = potSym (jp 4 φ x) := by
  obtain ⟨ψ, hψ, hev⟩ := cb_exists_global_eventually isOpen_ball h.hφ hx
  rw [cb_jp_congr 4 hev, ← potX_graphX hJi hJJ hU y₀ hψ (le_refl 4) x,
    vm_potX_congr (pf_graph_germ h hx hev)]
  have hmaps : ∀ z ∈ ball (0 : 𝔼) 5, ξ z ∈ ball (0 : 𝔼) 6 := fun z hz => pf_ball56 (h.pos z hz).1
  exact (potX_comp (J := J) isOpen_ball h.hX h.hinj h.hlag isOpen_ball h.hξ hmaps
    (fun z hz => pf_ξ_det h hz) hx).symm

theorem pf_jet2 (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    (hU : IsLagFrame J U) (h : PfSlice J U y₀ X ξ φ) {x : 𝔼} (hx : x ∈ ball (0 : 𝔼) 5) :
    jet 2 φ x = slopeZ J U (fderiv ℝ X (ξ x)) := by
  obtain ⟨ψ, hψ, hev⟩ := cb_exists_global_eventually isOpen_ball h.hφ hx
  have hgerm := pf_graph_germ h hx hev
  have hmem := (h.pos x hx).1
  have hmem6 : ξ x ∈ ball (0 : 𝔼) 6 := pf_ball56 hmem
  rw [jet_congr_of_eventuallyEq hev 2, ← slopeZ_graphX hJi hJJ hU hψ y₀ x, hgerm.fderiv_eq]
  have hτd : DifferentiableAt ℝ ξ x := ir_dτ isOpen_ball h.hξ hx
  have hXd : DifferentiableAt ℝ X (ξ x) := ir_dX isOpen_ball h.hX hmem6
  rw [fderiv_comp x hXd hτd]
  have hleft : Function.LeftInverse (fderiv ℝ (fun y' => horC U y₀ (X y')) (ξ x))
      (fderiv ℝ ξ x) := fun v => by
    have := congrArg (fun L : 𝔼 →L[ℝ] 𝔼 => L v) (pf_chi_deriv h hx)
    simpa using this
  obtain ⟨T, hT⟩ := gl_exists_equiv hleft.injective
  have hL : Function.Injective (horCL U ∘L fderiv ℝ X (ξ x)) := by
    rw [← pf_fderiv_chi h hmem6]
    exact pf_inj_of_close (h.bd _ hmem)
  rw [← hT]
  exact slopeZ_comp T hL

theorem pf_vel (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    (hU : IsLagFrame J U) (h : PfSlice J U y₀ X ξ φ) {y : 𝔼} (hy : y ∈ ball (0 : 𝔼) (11 / 2))
    (hy5 : horC U y₀ (X y) ∈ ball (0 : 𝔼) 5) (c : Fin n) :
    horCL U (velX J X y) c = vSym (jp 5 φ (horC U y₀ (X y))) c := by
  have hy6 : y ∈ ball (0 : 𝔼) 6 := pf_ball56 hy
  obtain ⟨ψ, hψ, hev⟩ := cb_exists_global_eventually isOpen_ball h.hφ hy5
  have hχc : ContinuousAt (fun z => horC U y₀ (X z)) y := (pf_chi_diff h hy6).continuousAt
  have hgraph : ∀ᶠ z in 𝓝 y, graphX J U y₀ φ (horC U y₀ (X z)) =
      graphX J U y₀ ψ (horC U y₀ (X z)) := by
    have hev' : ∀ᶠ x' in 𝓝 (horC U y₀ (X y)), graphX J U y₀ φ x' = graphX J U y₀ ψ x' := by
      filter_upwards [hev.eventually_nhds] with x' hx' using cge_graphX_congr hx'
    exact hχc.tendsto.eventually hev'
  have h1 : ∀ᶠ z in 𝓝 y, z ∈ ball (0 : 𝔼) (11 / 2) := isOpen_ball.mem_nhds hy
  have h2 : ∀ᶠ z in 𝓝 y, horC U y₀ (X z) ∈ ball (0 : 𝔼) 5 :=
    hχc.preimage_mem_nhds (isOpen_ball.mem_nhds hy5)
  have heq : ∀ᶠ z in 𝓝 y, X z = graphX J U y₀ ψ (horC U y₀ (X z)) := by
    filter_upwards [hgraph, h1, h2] with z hz hz1 hz2
    rw [← hz]
    have := (h.pos _ hz2).2.2
    rw [h.left z hz1 hz2] at this
    exact this
  obtain ⟨hvel, -, -, -⟩ := vb_reparam hJi hJJ hU y₀ hψ (X := X) isOpen_ball h.hX h.hinj hy6 heq
  have e : velX J X y = velX J (graphX J U y₀ ψ) (horC U y₀ (X y)) := hvel.eq_of_nhds
  rw [e, cb_jp_congr 5 hev]
  have hadm : ∃ K, jp 5 ψ (horC U y₀ (X y)) ∈ admSet n 5 K := by
    refine ⟨‖jp 5 ψ (horC U y₀ (X y))‖ + 1, fun m _ => ?_⟩
    have := (norm_le_pi_norm (jp 5 ψ (horC U y₀ (X y))).2 m).trans
      (norm_snd_le (jp 5 ψ (horC U y₀ (X y))))
    linarith
  exact horCL_velX_graphX hJi hJJ hU y₀ hψ (N := 5) (by norm_num) _ hadm c

end Slice

section Param

variable {n : ℕ} {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
  {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
  {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]

local notation "𝔼" => EuclideanSpace ℝ (Fin n)

theorem pf_isParam (ι : EuclideanSpace ℝ (Fin n) ≃L[ℝ] E) {J : V →L[ℝ] V} {F : M → V}
    (hF : IsLagImm ι J F) {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} {y₀ : V}
    {σ : EuclideanSpace ℝ (Fin n) → M} (hσ : IsParam ι (ball (0 : 𝔼) 6) σ)
    {ξ : EuclideanSpace ℝ (Fin n) → EuclideanSpace ℝ (Fin n)}
    (hξ : ContDiffOn ℝ ∞ ξ (ball (0 : 𝔼) 5))
    (hpos : ∀ x ∈ ball (0 : 𝔼) 5,
      ξ x ∈ ball (0 : 𝔼) (11 / 2) ∧ horC U y₀ (F (σ (ξ x))) = x)
    (hleft : ∀ y ∈ ball (0 : 𝔼) (11 / 2), horC U y₀ (F (σ y)) ∈ ball (0 : 𝔼) 5 →
      ξ (horC U y₀ (F (σ y))) = y) :
    IsParam ι (ball (0 : 𝔼) 5) (fun x => σ (ξ x)) := by
  have h56 : ball (0 : 𝔼) (11 / 2) ⊆ ball (0 : 𝔼) 6 := ball_subset_ball (by norm_num)
  have hmaps : ∀ x ∈ ball (0 : 𝔼) 5, ξ x ∈ ball (0 : 𝔼) 6 := fun x hx => h56 (hpos x hx).1
  have hξc : ContinuousOn ξ (ball (0 : 𝔼) 5) := hξ.continuousOn
  have hcont : ContinuousOn (fun x => σ (ξ x)) (ball (0 : 𝔼) 5) := hσ.cont.comp hξc hmaps
  have hts : ∀ z : M, ∀ x ∈ ball (0 : 𝔼) 5 ∩ (fun x => σ (ξ x)) ⁻¹' (extChartAt 𝓘(ℝ, E) z).source,
      ContDiffAt ℝ ∞ (chartFwd ι z ∘ fun x => σ (ξ x)) x := by
    intro z x hx
    obtain ⟨hxB, hxS⟩ := hx
    have hW := im_isOpen_W ι hσ z
    have hmem : ξ x ∈ ball (0 : 𝔼) 6 ∩ σ ⁻¹' (extChartAt 𝓘(ℝ, E) z).source := ⟨hmaps x hxB, hxS⟩
    have h1 : ContDiffAt ℝ ∞ (chartFwd ι z ∘ σ) (ξ x) :=
      (hσ.trans_smooth z).contDiffAt (hW.mem_nhds hmem)
    have h2 : ContDiffAt ℝ ∞ ξ x := hξ.contDiffAt (isOpen_ball.mem_nhds hxB)
    exact h1.comp x h2
  refine
  { isOpen := isOpen_ball
    injOn := ?_
    isOpen_image := ?_
    cont := hcont
    inv := ⟨fun q => horC U y₀ (F q),
      (grl_continuous_horC.comp (grl_continuous_F ι hF)).continuousOn,
      fun x hx => (hpos x hx).2⟩
    trans_smooth := fun z x hx => (hts z x hx).contDiffWithinAt
    trans_det := ?_ }
  · intro a ha b hb hab
    have h1 : ξ a = ξ b := hσ.injOn (hmaps a ha) (hmaps b hb) hab
    calc a = horC U y₀ (F (σ (ξ a))) := (hpos a ha).2.symm
      _ = horC U y₀ (F (σ (ξ b))) := by rw [h1]
      _ = b := (hpos b hb).2
  · obtain ⟨κ, hκc, hκ⟩ := hσ.inv
    have hχ : ContinuousOn (fun y => horC U y₀ (F (σ y))) (ball (0 : 𝔼) 6) :=
      (grl_continuous_horC.comp (grl_continuous_F ι hF)).comp_continuousOn hσ.cont
    have hSo : IsOpen (ball (0 : 𝔼) (11 / 2) ∩
        (fun y => horC U y₀ (F (σ y))) ⁻¹' ball (0 : 𝔼) 5) :=
      (hχ.mono h56).isOpen_inter_preimage isOpen_ball isOpen_ball
    have hopen := hκc.isOpen_inter_preimage hσ.isOpen_image hSo
    have himg : (fun x => σ (ξ x)) '' ball (0 : 𝔼) 5 = σ '' ball (0 : 𝔼) 6 ∩
        κ ⁻¹' (ball (0 : 𝔼) (11 / 2) ∩ (fun y => horC U y₀ (F (σ y))) ⁻¹' ball (0 : 𝔼) 5) := by
      ext q
      constructor
      · rintro ⟨x, hx, rfl⟩
        refine ⟨⟨ξ x, hmaps x hx, rfl⟩, ?_⟩
        show κ (σ (ξ x)) ∈ ball (0 : 𝔼) (11 / 2) ∩
          (fun y => horC U y₀ (F (σ y))) ⁻¹' ball (0 : 𝔼) 5
        rw [hκ _ (hmaps x hx)]
        exact ⟨(hpos x hx).1, by
          show horC U y₀ (F (σ (ξ x))) ∈ ball (0 : 𝔼) 5
          rw [(hpos x hx).2]; exact hx⟩
      · rintro ⟨⟨y, hy, rfl⟩, hS⟩
        have hS' : κ (σ y) ∈ ball (0 : 𝔼) (11 / 2) ∩
            (fun y => horC U y₀ (F (σ y))) ⁻¹' ball (0 : 𝔼) 5 := hS
        rw [hκ y hy] at hS'
        obtain ⟨hy', hy5⟩ := hS'
        refine ⟨horC U y₀ (F (σ y)), hy5, ?_⟩
        show σ (ξ (horC U y₀ (F (σ y)))) = σ y
        rw [hleft y hy' hy5]
    rw [himg]
    exact hopen
  · intro z x hx
    obtain ⟨hxB, hxS⟩ := hx
    have hτ : DifferentiableAt ℝ (chartFwd ι z ∘ fun x => σ (ξ x)) x :=
      (hts z x ⟨hxB, hxS⟩).differentiableAt (by simp)
    have hτx : (chartFwd ι z ∘ fun x => σ (ξ x)) x ∈ chartDomE ι z := im_chartFwd_mem ι hxS
    have hFd : DifferentiableAt ℝ (F ∘ chartInv ι z) ((chartFwd ι z ∘ fun x => σ (ξ x)) x) :=
      ((hF.smooth z).contDiffAt ((im_isOpen_chartDomE ι z).mem_nhds hτx)).differentiableAt
        (by simp)
    have hΦ : DifferentiableAt ℝ (horC U y₀ ∘ (F ∘ chartInv ι z))
        ((chartFwd ι z ∘ fun x => σ (ξ x)) x) :=
      (hasFDerivAt_horC ((F ∘ chartInv ι z) ((chartFwd ι z ∘ fun x => σ (ξ x)) x))).differentiableAt.comp
        _ hFd
    refine im_det_ne_zero_of_left_inv hτ hΦ ?_
    have h1 : ∀ᶠ w in 𝓝 x, w ∈ ball (0 : 𝔼) 5 := isOpen_ball.mem_nhds hxB
    have h2 : ∀ᶠ w in 𝓝 x, σ (ξ w) ∈ (extChartAt 𝓘(ℝ, E) z).source :=
      (hcont.continuousAt (isOpen_ball.mem_nhds hxB)).preimage_mem_nhds
        ((isOpen_extChartAt_source z).mem_nhds hxS)
    filter_upwards [h1, h2] with w hw1 hw2
    simp only [Function.comp_apply, id]
    rw [im_chartInv_chartFwd ι hw2]
    exact (hpos w hw1).2

end Param

section Main


theorem exists_piece_flow {n : ℕ} {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    [FiniteDimensional ℝ E] {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
    {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V] [FiniteDimensional ℝ V]
    (ι : EuclideanSpace ℝ (Fin n) ≃L[ℝ] E) {J : V →L[ℝ] V}
    (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    (hdim : Module.finrank ℝ V = 2 * n)
    {c₀ c₁ : ℝ} (hc₀ : 0 ≤ c₀) (hc₀1 : c₀ ≤ 1 / 100) (hc₁ : 0 ≤ c₁)
    {I : Set ℝ} {Φ : ℝ → M → V} (hΦ : IsFlowM ι J I Φ) {t₀ T : ℝ} (hT : t₀ < T) (hI : Ioo t₀ T ⊆ I)
    (hA : ∀ t ∈ Ioo t₀ T, ∀ z : M, ∀ y ∈ chartDomE ι z, aNormSq J (Φ t ∘ chartInv ι z) y ≤ c₁ ^ 2)
    {p : M} {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} {σ : EuclideanSpace ℝ (Fin n) → M}
    {φ₀ : EuclideanSpace ℝ (Fin n) → ℝ} (hch : IsGraphChart ι J (Φ t₀) c₀ p U σ φ₀)
    (hpers : ∀ t ∈ Ioo t₀ T, ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) (11 / 2),
      ‖fderiv ℝ (fun y' => horC U (Φ t₀ p) (Φ t (σ y'))) y - ContinuousLinearMap.id ℝ _‖ ≤ 1 / 8 ∧
      ‖horC U (Φ t₀ p) (Φ t (σ y)) - y‖ ≤ 1 / 8 ∧
      ‖verC J U (Φ t₀ p) (Φ t (σ y))‖ ≤ 1 / 2 ∧
      ‖slopeZ J U (fderiv ℝ (fun y' => Φ t (σ y')) y)‖ ≤ 2 * c₀) :
    ∃ (ξ : ℝ × EuclideanSpace ℝ (Fin n) → EuclideanSpace ℝ (Fin n))
      (φ : ℝ × EuclideanSpace ℝ (Fin n) → ℝ),
      ContDiffOn ℝ ∞ φ (Ioo t₀ T ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) 5) ∧
      ∀ t ∈ Ioo t₀ T,
        IsParam ι (ball (0 : EuclideanSpace ℝ (Fin n)) 5) (fun x => σ (ξ (t, x))) ∧
        (∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 5,
          ξ (t, x) ∈ ball (0 : EuclideanSpace ℝ (Fin n)) (11 / 2) ∧
          horC U (Φ t₀ p) (Φ t (σ (ξ (t, x)))) = x ∧
          Φ t (σ (ξ (t, x))) = graphX J U (Φ t₀ p) (fun x' => φ (t, x')) x ∧
          HasDerivAt (fun s => φ (s, x)) (potSym (jp 4 (fun x' => φ (t, x')) x)) t ∧
          ‖gradE (fun x' => φ (t, x')) x‖ ≤ 1 / 2 ∧
          ‖jet 2 (fun x' => φ (t, x')) x‖ ≤ 2 * c₀ ∧
          ‖jet 3 (fun x' => φ (t, x')) x‖ ≤ 2 * c₁) ∧
        (∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) (11 / 2),
          horC U (Φ t₀ p) (Φ t (σ y)) ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 5 →
          ξ (t, horC U (Φ t₀ p) (Φ t (σ y))) = y ∧
          ∀ c : Fin n, HasDerivAt (fun s => horC U (Φ t₀ p) (Φ s (σ y)) c)
            (vSym (jp 5 (fun x' => φ (t, x')) (horC U (Φ t₀ p) (Φ t (σ y)))) c) t) := by
  classical
  have hY : ContDiffOn ℝ ∞ (fun q : ℝ × (EuclideanSpace ℝ (Fin n)) => Φ q.1 (σ q.2))
      (Ioo t₀ T ×ˢ ball (0 : (EuclideanSpace ℝ (Fin n))) 6) :=
    (cge_contDiffOn_family ι hΦ.smooth hch.param hΦ.isOpen).mono (prod_mono hI subset_rfl)
  have himm : ∀ t ∈ Ioo t₀ T, ∀ y ∈ ball (0 : (EuclideanSpace ℝ (Fin n))) 6,
      Function.Injective (fderiv ℝ (Φ t ∘ σ) y) ∧
        ∀ u w, ⟪J (fderiv ℝ (Φ t ∘ σ) y u), fderiv ℝ (Φ t ∘ σ) y w⟫ = 0 := fun t ht y hy =>
    ⟨((hΦ.imm t (hI ht)).comp_param ι hch.param).2.1 y hy,
      ((hΦ.imm t (hI ht)).comp_param ι hch.param).2.2 y hy⟩
  obtain ⟨ξ, φ, hξ, hφ, hper⟩ := exists_time_graph_unif hJi hJJ hch.frame hdim (Φ t₀ p) hT
    (Y := fun q : ℝ × (EuclideanSpace ℝ (Fin n)) => Φ q.1 (σ q.2)) hY
    (fun t ht y hy => (himm t ht y hy).1) (fun t ht y hy => (himm t ht y hy).2)
    (fun t ht y hy => cge_hasDerivAt_family ι hΦ hch.param (hI ht) hy)
    (fun t ht y hy => ⟨(hpers t ht y hy).1, (hpers t ht y hy).2.1⟩)
  refine ⟨ξ, φ, hφ, fun t ht => ?_⟩
  obtain ⟨-, hpos, hleft⟩ := hper t ht
  have hS : PfSlice J U (Φ t₀ p) (fun y => Φ t (σ y)) (fun x => ξ (t, x)) (fun x => φ (t, x)) :=
    { hX := tg_slice hY ht
      hinj := fun y hy => (himm t ht y hy).1
      hlag := fun y hy => (himm t ht y hy).2
      hξ := tg_slice hξ ht
      hφ := tg_slice hφ ht
      pos := fun x hx => ⟨(hpos x hx).1, (hpos x hx).2.1, (hpos x hx).2.2.1⟩
      left := hleft
      bd := fun y hy => (hpers t ht y hy).1 }
  have hpar : IsParam ι (ball (0 : (EuclideanSpace ℝ (Fin n))) 5) (fun x => σ (ξ (t, x))) :=
    pf_isParam ι (hΦ.imm t (hI ht)) hch.param hS.hξ
      (fun x hx => ⟨(hpos x hx).1, (hpos x hx).2.1⟩) hleft
  refine ⟨hpar, fun x hx => ?_, fun y hy hy5 => ⟨hleft y hy hy5, fun c => ?_⟩⟩
  · obtain ⟨hm, hh, hg, hd⟩ := hpos x hx
    have hj2 : ‖jet 2 (fun x' => φ (t, x')) x‖ ≤ 2 * c₀ := by
      rw [pf_jet2 hJi hJJ hch.frame hS hx]
      exact (hpers t ht _ hm).2.2.2
    refine ⟨hm, hh, hg, hd.congr_deriv (pf_pot hJi hJJ hch.frame hS hx), ?_, hj2, ?_⟩
    · have e : gradE (fun x' => φ (t, x')) x = verC J U (Φ t₀ p) (Φ t (σ (ξ (t, x)))) := by
        rw [hg]
        exact (verC_graphX hJi hJJ hch.frame _ x).symm
      rw [e]
      exact (hpers t ht _ hm).2.2.1
    · have h1 := cb_aNormSq_param_le ι (hΦ.imm t (hI ht)) (hA t ht) hpar hx
      have h2 := cb_aNormSq_graph hJi hJJ hch.frame isOpen_ball hS.hφ
        (F := Φ t) (σ := fun x => σ (ξ (t, x))) (y₀ := Φ t₀ p)
        (fun x' hx' => (hpos x' hx').2.2.1) hx
      rw [h2] at h1
      have h3 := norm_jet3_sq_le (fun x' => φ (t, x')) x
      have hs0 : 0 ≤ ‖jet 2 (fun x' => φ (t, x')) x‖ := norm_nonneg _
      have hs1 : ‖jet 2 (fun x' => φ (t, x')) x‖ ^ 2 ≤ (2 * c₀) ^ 2 := by gcongr
      have hs2 : (1 + ‖jet 2 (fun x' => φ (t, x')) x‖ ^ 2) ^ 3 ≤ 4 := by
        have : ‖jet 2 (fun x' => φ (t, x')) x‖ ^ 2 ≤ 1 / 2500 := by nlinarith
        have h4 : (1 + ‖jet 2 (fun x' => φ (t, x')) x‖ ^ 2) ^ 3 ≤ (1 + 1 / 2500) ^ 3 := by
          gcongr
        linarith [h4, (by norm_num : ((1 : ℝ) + 1 / 2500) ^ 3 ≤ 4)]
      have h5 : ‖jet 3 (fun x' => φ (t, x')) x‖ ^ 2 ≤ (2 * c₁) ^ 2 := by
        calc ‖jet 3 (fun x' => φ (t, x')) x‖ ^ 2
            ≤ (1 + ‖jet 2 (fun x' => φ (t, x')) x‖ ^ 2) ^ 3 *
              gNormSq (ginvFun (jet 2 (fun x' => φ (t, x')) x))
                (jet 3 (fun x' => φ (t, x')) x) := h3
          _ ≤ (1 + ‖jet 2 (fun x' => φ (t, x')) x‖ ^ 2) ^ 3 * c₁ ^ 2 :=
              mul_le_mul_of_nonneg_left h1 (by positivity)
          _ ≤ 4 * c₁ ^ 2 := by gcongr
          _ = (2 * c₁) ^ 2 := by ring
      exact (pow_le_pow_iff_left₀ (norm_nonneg _) (by positivity) two_ne_zero).1 h5
  · have hy6 : y ∈ ball (0 : (EuclideanSpace ℝ (Fin n))) 6 := pf_ball56 hy
    have hv : HasDerivAt (fun s => Φ s (σ y)) (velX J (fun y' => Φ t (σ y')) y) t :=
      cge_hasDerivAt_family ι hΦ hch.param (hI ht) hy6
    have h2 := (hasFDerivAt_horC (U := U) (y₀ := Φ t₀ p) (Φ t (σ y))).comp_hasDerivAt t hv
    have h3 := (EuclideanSpace.proj c : (EuclideanSpace ℝ (Fin n)) →L[ℝ] ℝ).hasFDerivAt.comp_hasDerivAt t h2
    exact h3.congr_deriv (pf_vel hJi hJJ hch.frame hS hy hy5 c)

end Main

/-- info: 'CWFlow.Extension.exists_piece_flow' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_piece_flow

end CWFlow.Extension
