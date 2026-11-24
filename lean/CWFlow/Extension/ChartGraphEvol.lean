import CWFlow.Extension.TimeGraphFlow
import CWFlow.Extension.GraphChart
import CWFlow.Extension.GaugeTerm
import CWFlow.Extension.VolM

/-!
# A flow of Lagrangian immersions in a graph chart is a graph evolution

A flow `Φ` of Lagrangian immersions of a compact manifold `M` (`IsFlowM`) is read in a graph chart taken at
time `t₀`. For `t` near `t₀` it is described by a globally smooth graph function `Ψ (t, x)` solving the scalar
equation `∂ₜΨ = potSym (jp N Ψ)` on the ball of radius `3` with bounded low jets (`GraphEvol N 2 3`), together
with the reparametrization `ξ` (`exists_chart_graphEvol`).
-/

open Metric Set
open CWFlow.Interp
open scoped ContDiff Manifold Topology RealInnerProductSpace

namespace CWFlow.Extension

section Lemmas

variable {n : ℕ}

local notation "𝔼" => EuclideanSpace ℝ (Fin n)

/-! ### Part 1: bounds on the low jets of the chart function -/

theorem cge_norm_jet_one (f : 𝔼 → ℝ) (x : 𝔼) : ‖jet 1 f x‖ = ‖gradE f x‖ := by
  have h : ‖jet 1 f x‖ ^ 2 = ‖gradE f x‖ ^ 2 := by
    rw [norm_jet_sq, EuclideanSpace.norm_sq_eq]
    refine Fintype.sum_equiv (Equiv.funUnique (Fin 1) (Fin n)) _ _ fun I => ?_
    rw [iteratedFDeriv_one_apply]
    rfl
  exact (sq_eq_sq₀ (norm_nonneg _) (norm_nonneg _)).mp h

theorem cge_norm_fderiv_le (f : 𝔼 → ℝ) (x : 𝔼) : ‖fderiv ℝ f x‖ ≤ ‖gradE f x‖ := by
  rw [← cge_norm_jet_one, ← norm_iteratedFDeriv_one (𝕜 := ℝ) (x := x) f]
  exact norm_iteratedFDeriv_le_norm_jet 1 f x

theorem cge_jet_bounds {c₀ K₃ : ℝ} (hc₀ : 0 ≤ c₀) (hc₀1 : c₀ ≤ 1 / 100) (hK₃ : K₃ ≤ 1)
    {φ₀ : 𝔼 → ℝ} (hφ : ContDiffOn ℝ ∞ φ₀ (ball (0 : 𝔼) 6)) (h0 : φ₀ 0 = 0)
    (hgr : ∀ x ∈ ball (0 : 𝔼) 6, ‖jet 2 φ₀ x‖ ≤ c₀ ∧ ‖gradE φ₀ x‖ ≤ c₀ * ‖x‖)
    (hK : ∀ x ∈ ball (0 : 𝔼) 6, ‖jet 3 φ₀ x‖ ≤ K₃) :
    ∀ x ∈ closedBall (0 : 𝔼) 3, ∀ m : ℕ, m ≤ 3 → ‖jet m φ₀ x‖ < 2 := by
  intro x hx m hm
  have hx6 : x ∈ ball (0 : 𝔼) 6 := closedBall_subset_ball (by norm_num) hx
  have hx3 : ‖x‖ ≤ 3 := by simpa using hx
  interval_cases m
  · rw [norm_jet_zero, Real.norm_eq_abs]
    have hb : ∀ y ∈ closedBall (0 : 𝔼) 3, ‖fderiv ℝ φ₀ y‖ ≤ 3 / 100 := by
      intro y hy
      have hy6 : y ∈ ball (0 : 𝔼) 6 := closedBall_subset_ball (by norm_num) hy
      have hy3 : ‖y‖ ≤ 3 := by simpa using hy
      calc ‖fderiv ℝ φ₀ y‖ ≤ ‖gradE φ₀ y‖ := cge_norm_fderiv_le _ _
        _ ≤ c₀ * ‖y‖ := (hgr y hy6).2
        _ ≤ 1 / 100 * 3 := by nlinarith [norm_nonneg y]
        _ = 3 / 100 := by norm_num
    have hd : ∀ y ∈ closedBall (0 : 𝔼) 3, DifferentiableAt ℝ φ₀ y := fun y hy =>
      ((hφ.contDiffAt (isOpen_ball.mem_nhds
        (closedBall_subset_ball (by norm_num) hy))).differentiableAt (by simp))
    have := (convex_closedBall (0 : 𝔼) 3).norm_image_sub_le_of_norm_fderiv_le hd hb
      (mem_closedBall_self (by norm_num)) hx
    rw [h0, sub_zero, sub_zero, Real.norm_eq_abs] at this
    nlinarith [norm_nonneg x]
  · rw [cge_norm_jet_one]
    have := (hgr x hx6).2
    nlinarith [norm_nonneg x]
  · have := (hgr x hx6).1
    linarith
  · have := hK x hx6
    linarith

/-! ### Part 2: the family `Y (t, y) = Φ t (σ y)` -/

section Family

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
  {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
  {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]

theorem cge_contDiffOn_family (ι : EuclideanSpace ℝ (Fin n) ≃L[ℝ] E) {I : Set ℝ}
    {Φ : ℝ → M → V}
    (hs : ∀ z : M, ContDiffOn ℝ ∞ (fun q : ℝ × 𝔼 => Φ q.1 (chartInv ι z q.2))
      (I ×ˢ chartDomE ι z))
    {B : Set 𝔼} {σ : 𝔼 → M} (hσ : IsParam ι B σ) (hI : IsOpen I) :
    ContDiffOn ℝ ∞ (fun q : ℝ × 𝔼 => Φ q.1 (σ q.2)) (I ×ˢ B) := by
  intro q hq
  obtain ⟨t, y⟩ := q
  set z := σ y with hz
  have hW := im_isOpen_W ι hσ z
  have hO : IsOpen (I ×ˢ (B ∩ σ ⁻¹' (extChartAt 𝓘(ℝ, E) z).source)) := hI.prod hW
  have hmem : (t, y) ∈ I ×ˢ (B ∩ σ ⁻¹' (extChartAt 𝓘(ℝ, E) z).source) :=
    ⟨hq.1, im_mem_W (σ := σ) hq.2⟩
  have hH : ContDiffOn ℝ ∞ (fun q : ℝ × 𝔼 => ((q.1, chartFwd ι z (σ q.2)) : ℝ × 𝔼))
      (I ×ˢ (B ∩ σ ⁻¹' (extChartAt 𝓘(ℝ, E) z).source)) :=
    contDiffOn_fst.prodMk ((hσ.trans_smooth z).comp contDiffOn_snd (fun q hq => hq.2))
  have hG := (hs z).comp hH (fun q hq => ⟨hq.1, im_chartFwd_mem ι hq.2.2⟩)
  have h2 : ContDiffOn ℝ ∞ (fun q : ℝ × 𝔼 => Φ q.1 (σ q.2))
      (I ×ˢ (B ∩ σ ⁻¹' (extChartAt 𝓘(ℝ, E) z).source)) :=
    hG.congr fun q hq => by
      simp only [Function.comp_apply]
      rw [im_chartInv_chartFwd ι hq.2.2]
  exact ((h2.contDiffAt (hO.mem_nhds hmem))).contDiffWithinAt

theorem cge_velX_congr {J : V →L[ℝ] V} {X X' : 𝔼 → V} {x : 𝔼} (h : X =ᶠ[𝓝 x] X') :
    velX J X x = velX J X' x := by
  have hp : potX J X =ᶠ[𝓝 x] potX J X' := by
    filter_upwards [h.eventually_nhds] with y hy using vm_potX_congr hy
  unfold velX
  rw [gX_congr h, hp.fderiv_eq]
  simp only [dX_congr h]

theorem cge_hasDerivAt_family (ι : EuclideanSpace ℝ (Fin n) ≃L[ℝ] E) {J : V →L[ℝ] V} {I : Set ℝ}
    {Φ : ℝ → M → V} (hΦ : IsFlowM ι J I Φ) {B : Set 𝔼} {σ : 𝔼 → M} (hσ : IsParam ι B σ)
    {t : ℝ} (ht : t ∈ I) {y : 𝔼} (hy : y ∈ B) :
    HasDerivAt (fun s => Φ s (σ y)) (velX J (fun y' => Φ t (σ y')) y) t := by
  set z := σ y with hz
  have hyW := im_mem_W (E := E) (σ := σ) hy
  have hW := im_isOpen_W ι hσ z
  have hw : chartFwd ι z (σ y) ∈ chartDomE ι z := im_chartFwd_mem ι (mem_extChartAt_source z)
  have h1 := hΦ.vel t ht z _ hw
  have hinv : chartInv ι z (chartFwd ι z (σ y)) = σ y :=
    im_chartInv_chartFwd ι (mem_extChartAt_source z)
  rw [hinv] at h1
  have himm := hΦ.imm t ht
  have h2 := velX_comp (J := J) (X := Φ t ∘ chartInv ι z) (τ := chartFwd ι z ∘ σ)
    (Ω := chartDomE ι z) (N := B ∩ σ ⁻¹' (extChartAt 𝓘(ℝ, E) z).source)
    (im_isOpen_chartDomE ι z) (himm.smooth z) (himm.inj z) (himm.lag z) hW (hσ.trans_smooth z)
    (fun w hw => im_chartFwd_mem ι hw.2) (hσ.trans_det z) hyW
  have h3 : (fun y' => Φ t (σ y')) =ᶠ[𝓝 y] ((Φ t ∘ chartInv ι z) ∘ (chartFwd ι z ∘ σ)) := by
    filter_upwards [hW.mem_nhds hyW] with w hw
    simp only [Function.comp_apply]
    rw [im_chartInv_chartFwd ι hw.2]
  rw [cge_velX_congr h3, h2]
  exact h1

end Family

/-! ### Part 3: globalization by a cutoff -/

theorem cge_exists_cutoff : ∃ ψ : 𝔼 → ℝ, ContDiff ℝ ∞ ψ ∧ (∀ x ∈ closedBall (0 : 𝔼) 4, ψ x = 1) ∧
    tsupport ψ ⊆ closedBall (0 : 𝔼) (9 / 2) := by
  let b : ContDiffBump (0 : 𝔼) := ⟨4, 9 / 2, by norm_num, by norm_num⟩
  refine ⟨b, b.contDiff, fun x hx => b.one_of_mem_closedBall hx, ?_⟩
  rw [b.tsupport_eq]

open Classical in
/-- The globalization of `φ'` (defined on the ball of radius `5`) by the cutoff `ψ`. -/
noncomputable def cge_glue (ψ : 𝔼 → ℝ) (φ' : ℝ × 𝔼 → ℝ) (q : ℝ × 𝔼) : ℝ :=
  if q.2 ∈ ball (0 : 𝔼) 5 then ψ q.2 * φ' q else 0

theorem cge_glue_eq {ψ : 𝔼 → ℝ} (hψ1 : ∀ x ∈ closedBall (0 : 𝔼) 4, ψ x = 1) {φ' : ℝ × 𝔼 → ℝ}
    {q : ℝ × 𝔼} (hq : q.2 ∈ closedBall (0 : 𝔼) 4) : cge_glue ψ φ' q = φ' q := by
  have h5 : q.2 ∈ ball (0 : 𝔼) 5 := closedBall_subset_ball (by norm_num) hq
  simp [cge_glue, h5, hψ1 _ hq]

theorem cge_contDiffOn_glue {T : Set ℝ} (hT : IsOpen T) {ψ : 𝔼 → ℝ} (hψ : ContDiff ℝ ∞ ψ)
    (hsupp : tsupport ψ ⊆ closedBall (0 : 𝔼) (9 / 2)) {φ' : ℝ × 𝔼 → ℝ}
    (hφ' : ContDiffOn ℝ ∞ φ' (T ×ˢ ball (0 : 𝔼) 5)) :
    ContDiffOn ℝ ∞ (cge_glue ψ φ') (T ×ˢ (univ : Set 𝔼)) := by
  intro q hq
  by_cases hx : q.2 ∈ ball (0 : 𝔼) 5
  · have hO : IsOpen (T ×ˢ ball (0 : 𝔼) 5) := hT.prod isOpen_ball
    have h1 : ContDiffOn ℝ ∞ (fun q : ℝ × 𝔼 => ψ q.2 * φ' q) (T ×ˢ ball (0 : 𝔼) 5) :=
      (hψ.comp contDiff_snd).contDiffOn.mul hφ'
    have h2 : ContDiffOn ℝ ∞ (cge_glue ψ φ') (T ×ˢ ball (0 : 𝔼) 5) :=
      h1.congr fun q' hq' => by simp [cge_glue, hq'.2]
    exact (h2.contDiffAt (hO.mem_nhds ⟨hq.1, hx⟩)).contDiffWithinAt
  · have hnot : q.2 ∉ closedBall (0 : 𝔼) (9 / 2) := by
      intro h
      apply hx
      exact closedBall_subset_ball (by norm_num) h
    have hopen : IsOpen {q' : ℝ × 𝔼 | q'.2 ∉ closedBall (0 : 𝔼) (9 / 2)} :=
      (isClosed_closedBall.isOpen_compl).preimage continuous_snd
    have hev : cge_glue ψ φ' =ᶠ[𝓝 q] fun _ => (0 : ℝ) := by
      filter_upwards [hopen.mem_nhds hnot] with q' hq'
      have : ψ q'.2 = 0 := image_eq_zero_of_notMem_tsupport fun h => hq' (hsupp h)
      simp [cge_glue, this]
    exact ((contDiffAt_const (c := (0 : ℝ))).congr_of_eventuallyEq hev).contDiffWithinAt

theorem cge_graphX_congr {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
    {J : V →L[ℝ] V} {U : 𝔼 →ₗᵢ[ℝ] V} {y : V} {f g : 𝔼 → ℝ} {x : 𝔼} (h : f =ᶠ[𝓝 x] g) :
    graphX J U y f x = graphX J U y g x := by
  unfold graphX gradE
  rw [h.fderiv_eq]

theorem cge_graphX_sub_const {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
    {J : V →L[ℝ] V} {U : 𝔼 →ₗᵢ[ℝ] V} {y : V} (f : 𝔼 → ℝ) (c : ℝ) (x : 𝔼) :
    graphX J U y (fun x' => f x' - c) x = graphX J U y f x := by
  unfold graphX gradE
  simp only [fderiv_sub_const]

/-! ### Part 4: the core -/

theorem cge_core {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V] [FiniteDimensional ℝ V]
    {J : V →L[ℝ] V} (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    (hdim : Module.finrank ℝ V = 2 * n) (N : ℕ) (hN : 4 ≤ N)
    {c₀ K₃ : ℝ} (hc₀ : 0 ≤ c₀) (hc₀1 : c₀ ≤ 1 / 100) (hK₃ : K₃ ≤ 1)
    {U : 𝔼 →ₗᵢ[ℝ] V} (hU : IsLagFrame J U) (y₀ : V) {a b t₀ : ℝ} (ht₀ : t₀ ∈ Ioo a b)
    {Y : ℝ × 𝔼 → V} (hY : ContDiffOn ℝ ∞ Y (Ioo a b ×ˢ ball (0 : 𝔼) 6))
    (hinj : ∀ t ∈ Ioo a b, ∀ y ∈ ball (0 : 𝔼) 6,
      Function.Injective (fderiv ℝ (fun y' => Y (t, y')) y))
    (hlag : ∀ t ∈ Ioo a b, ∀ y ∈ ball (0 : 𝔼) 6, ∀ u w,
      ⟪J (fderiv ℝ (fun y' => Y (t, y')) y u), fderiv ℝ (fun y' => Y (t, y')) y w⟫ = 0)
    (hvel : ∀ t ∈ Ioo a b, ∀ y ∈ ball (0 : 𝔼) 6,
      HasDerivAt (fun s => Y (s, y)) (velX J (fun y' => Y (t, y')) y) t)
    (hgraph : ∀ y ∈ ball (0 : 𝔼) 6, horC U y₀ (Y (t₀, y)) = y)
    {φ₀ : 𝔼 → ℝ} (hφ₀ : ContDiffOn ℝ ∞ φ₀ (ball (0 : 𝔼) 6)) (h0 : φ₀ 0 = 0)
    (hgr : ∀ x ∈ ball (0 : 𝔼) 6, Y (t₀, x) = graphX J U y₀ φ₀ x ∧ ‖jet 2 φ₀ x‖ ≤ c₀ ∧
      ‖gradE φ₀ x‖ ≤ c₀ * ‖x‖)
    (hK : ∀ x ∈ ball (0 : 𝔼) 6, ‖jet 3 φ₀ x‖ ≤ K₃) :
    ∃ δ > 0, Ioo (t₀ - δ) (t₀ + δ) ⊆ Ioo a b ∧
      ∃ (ξ : ℝ × 𝔼 → 𝔼) (Ψ : ℝ × 𝔼 → ℝ),
        GraphEvol N 2 3 (Ioo (t₀ - δ) (t₀ + δ)) Ψ ∧
        ContDiffOn ℝ ∞ ξ (Ioo (t₀ - δ) (t₀ + δ) ×ˢ ball (0 : 𝔼) 5) ∧
        (∀ x ∈ ball (0 : 𝔼) 4, Ψ (t₀, x) = φ₀ x) ∧
        (∀ x ∈ ball (0 : 𝔼) 5, ξ (t₀, x) = x) ∧
        ∀ t ∈ Ioo (t₀ - δ) (t₀ + δ),
          InjOn (fun y => horC U y₀ (Y (t, y))) (ball (0 : 𝔼) (11 / 2)) ∧
          (∀ y ∈ ball (0 : 𝔼) (11 / 2),
            Function.Injective (fderiv ℝ (fun y' => horC U y₀ (Y (t, y'))) y) ∧
            ‖fderiv ℝ (fun y' => horC U y₀ (Y (t, y'))) y - ContinuousLinearMap.id ℝ _‖ ≤ 1 / 4) ∧
          (∀ x ∈ ball (0 : 𝔼) 5,
            ξ (t, x) ∈ ball (0 : 𝔼) (11 / 2) ∧ horC U y₀ (Y (t, ξ (t, x))) = x) ∧
          (∀ x ∈ ball (0 : 𝔼) 4,
            Y (t, ξ (t, x)) = graphX J U y₀ (fun x' => Ψ (t, x')) x) ∧
          (∀ y ∈ ball (0 : 𝔼) (11 / 2),
            horC U y₀ (Y (t, y)) ∈ ball (0 : 𝔼) 5 → ξ (t, horC U y₀ (Y (t, y))) = y) := by
  classical
  obtain ⟨δ₁, hδ₁, hsub, ξ, φ, hξ, hφ, hper⟩ := exists_time_graph_flow hJi hJJ hU hdim y₀ ht₀ hY
    hinj hlag hvel hgraph (r := 5) (R' := 11 / 2) (by norm_num) (by norm_num) (by norm_num)
  have hT₁o : IsOpen (Ioo (t₀ - δ₁) (t₀ + δ₁)) := isOpen_Ioo
  have ht₀T : t₀ ∈ Ioo (t₀ - δ₁) (t₀ + δ₁) := ⟨by linarith, by linarith⟩
  have hb45 : ball (0 : 𝔼) 4 ⊆ ball (0 : 𝔼) 5 := ball_subset_ball (by norm_num)
  have hb56 : ball (0 : 𝔼) 5 ⊆ ball (0 : 𝔼) 6 := ball_subset_ball (by norm_num)
  -- the relation at `t₀`
  have hξ0 : ∀ x ∈ ball (0 : 𝔼) 5, ξ (t₀, x) = x := by
    intro x hx
    obtain ⟨hmem, hh, -, -⟩ := (hper t₀ ht₀T).2.2.1 x hx
    have := hgraph (ξ (t₀, x)) (ball_subset_ball (by norm_num) hmem)
    exact (hh.symm.trans this).symm
  have hgrad0 : ∀ x ∈ ball (0 : 𝔼) 5, gradE (fun x' => φ (t₀, x')) x = gradE φ₀ x := by
    intro x hx
    have h1 := ((hper t₀ ht₀T).2.2.1 x hx).2.2.1
    rw [hξ0 x hx] at h1
    have h2 := (hgr x (hb56 hx)).1
    calc gradE (fun x' => φ (t₀, x')) x
        = verC J U y₀ (graphX J U y₀ (fun x' => φ (t₀, x')) x) :=
          (verC_graphX hJi hJJ hU _ x).symm
      _ = verC J U y₀ (Y (t₀, x)) := by rw [h1]
      _ = verC J U y₀ (graphX J U y₀ φ₀ x) := by rw [h2]
      _ = gradE φ₀ x := verC_graphX hJi hJJ hU _ x
  have hφs : ContDiffOn ℝ ∞ (fun x' => φ (t₀, x')) (ball (0 : 𝔼) 5) := tg_slice hφ ht₀T
  obtain ⟨c, hc⟩ : ∃ c : ℝ, c = φ (t₀, 0) := ⟨_, rfl⟩
  have heq : ∀ x ∈ ball (0 : 𝔼) 5, φ (t₀, x) = φ₀ x + c := by
    have := eqOn_ball_of_fderiv_eq (x₀ := (0 : 𝔼)) (R := 5) (φ := fun x' => φ (t₀, x'))
      (ψ := fun x => φ₀ x + c) (hφs.differentiableOn (by simp))
      (((hφ₀.mono hb56).differentiableOn (by simp)).add_const c)
      (fun x hx => by
        rw [fderiv_add_const]
        exact gr_fderiv_eq_of_gradE (hgrad0 x hx))
      (by show φ (t₀, 0) = φ₀ 0 + c; rw [h0, zero_add, hc])
    exact fun x hx => this hx
  -- the globalized graph function
  obtain ⟨ψ, hψ, hψ1, hψs⟩ := cge_exists_cutoff (n := n)
  have hφ'S : ContDiffOn ℝ ∞ (fun q : ℝ × 𝔼 => φ q - c)
      (Ioo (t₀ - δ₁) (t₀ + δ₁) ×ˢ ball (0 : 𝔼) 5) := hφ.sub contDiffOn_const
  have hΨ : ContDiffOn ℝ ∞ (cge_glue ψ fun q : ℝ × 𝔼 => φ q - c)
      (Ioo (t₀ - δ₁) (t₀ + δ₁) ×ˢ (univ : Set 𝔼)) := cge_contDiffOn_glue hT₁o hψ hψs hφ'S
  obtain ⟨Ψ, hΨdef⟩ : ∃ Ψ : ℝ × 𝔼 → ℝ, Ψ = cge_glue ψ fun q : ℝ × 𝔼 => φ q - c := ⟨_, rfl⟩
  rw [← hΨdef] at hΨ
  have hΨeq : ∀ q : ℝ × 𝔼, q.2 ∈ closedBall (0 : 𝔼) 4 → Ψ q = φ q - c := fun q hq => by
    rw [hΨdef]; exact cge_glue_eq hψ1 hq
  have hΨ0 : ∀ x ∈ ball (0 : 𝔼) 4, Ψ (t₀, x) = φ₀ x := by
    intro x hx
    rw [hΨeq (t₀, x) (ball_subset_closedBall hx), heq x (hb45 hx)]
    ring
  have hΨslice : ∀ t ∈ Ioo (t₀ - δ₁) (t₀ + δ₁), ContDiff ℝ ∞ (fun y => Ψ (t, y)) := fun t ht =>
    contDiff_slice hT₁o hΨ ht
  have hgerm : ∀ t ∈ Ioo (t₀ - δ₁) (t₀ + δ₁), ∀ x ∈ ball (0 : 𝔼) 4,
      (fun y => Ψ (t, y)) =ᶠ[𝓝 x] (fun y => φ (t, y) - c) := by
    intro t ht x hx
    filter_upwards [isOpen_ball.mem_nhds hx] with y hy using
      hΨeq (t, y) (ball_subset_closedBall hy)
  -- the graph identity
  have hgrid : ∀ t ∈ Ioo (t₀ - δ₁) (t₀ + δ₁), ∀ x ∈ ball (0 : 𝔼) 4,
      Y (t, ξ (t, x)) = graphX J U y₀ (fun x' => Ψ (t, x')) x := by
    intro t ht x hx
    rw [((hper t ht).2.2.1 x (hb45 hx)).2.2.1, cge_graphX_congr (hgerm t ht x hx),
      cge_graphX_sub_const (fun x' => φ (t, x')) c x]
  -- the time derivative
  have hdt : ∀ t ∈ Ioo (t₀ - δ₁) (t₀ + δ₁), ∀ x ∈ ball (0 : 𝔼) 4,
      dt Ψ (t, x) = potX J (fun y => Y (t, y)) (ξ (t, x)) := by
    intro t ht x hx
    have hdiff : DifferentiableAt ℝ Ψ (t, x) :=
      (hΨ.contDiffAt ((hT₁o.prod isOpen_univ).mem_nhds ⟨ht, mem_univ x⟩)).differentiableAt
        (by simp)
    have h1 : HasDerivAt (fun s => Ψ (s, x)) (dt Ψ (t, x)) t :=
      tf_hasDerivAt_slice_t hdiff.hasFDerivAt
    have h2 : HasDerivAt (fun s => Ψ (s, x)) (potX J (fun y => Y (t, y)) (ξ (t, x))) t := by
      refine (((hper t ht).2.2.1 x (hb45 hx)).2.2.2.sub_const c).congr_of_eventuallyEq ?_
      filter_upwards [hT₁o.mem_nhds ht] with s hs
      exact hΨeq (s, x) (ball_subset_closedBall hx)
    exact h1.unique h2
  -- the equation
  have heqn : ∀ t ∈ Ioo (t₀ - δ₁) (t₀ + δ₁), ∀ x ∈ ball (0 : 𝔼) 3,
      dt Ψ (t, x) = potSym (jp N (fun y => Ψ (t, y)) x) := by
    intro t ht x hx
    have hx4 : x ∈ ball (0 : 𝔼) 4 := ball_subset_ball (by norm_num) hx
    have hx5 : x ∈ ball (0 : 𝔼) 5 := hb45 hx4
    have htab : t ∈ Ioo a b := hsub ht
    rw [hdt t ht x hx4, ← potX_graphX hJi hJJ hU y₀ (hΨslice t ht) hN x]
    have hgermG : (graphX J U y₀ fun y => Ψ (t, y)) =ᶠ[𝓝 x]
        ((fun y => Y (t, y)) ∘ fun z => ξ (t, z)) := by
      filter_upwards [isOpen_ball.mem_nhds hx4] with y hy
      exact (hgrid t ht y hy).symm
    rw [vm_potX_congr hgermG]
    have hX : ContDiffOn ℝ ∞ (fun y => Y (t, y)) (ball (0 : 𝔼) 6) := tg_slice hY htab
    have hτ : ContDiffOn ℝ ∞ (fun z => ξ (t, z)) (ball (0 : 𝔼) 5) := tg_slice hξ ht
    have hmaps : ∀ z ∈ ball (0 : 𝔼) 5, ξ (t, z) ∈ ball (0 : 𝔼) 6 := fun z hz =>
      ball_subset_ball (by norm_num) ((hper t ht).2.2.1 z hz).1
    refine (potX_comp (J := J) isOpen_ball hX (hinj t htab) (hlag t htab) isOpen_ball hτ hmaps
      ?_ hx5).symm
    intro z hz
    have hτd : DifferentiableAt ℝ (fun z => ξ (t, z)) z := ir_dτ isOpen_ball hτ hz
    have hYd : DifferentiableAt ℝ (fun y => Y (t, y)) (ξ (t, z)) :=
      ir_dX isOpen_ball hX (hmaps z hz)
    have hτ' : DifferentiableAt ℝ (fun y => horC U y₀ (Y (t, y))) (ξ (t, z)) :=
      DifferentiableAt.comp (g := horC U y₀) (f := fun y => Y (t, y)) (ξ (t, z))
        (hasFDerivAt_horC (U := U) (y₀ := y₀) (Y (t, ξ (t, z)))).differentiableAt hYd
    refine im_det_ne_zero_of_left_inv hτd hτ' ?_
    filter_upwards [isOpen_ball.mem_nhds hz] with w hw
    exact ((hper t ht).2.2.1 w hw).2.1
  -- admissibility near `t₀`
  have hAopen := isOpen_symbol_slice_dom (N := N) (isOpen_admSet (n := n) (N := N) 2) hT₁o hΨ
  have hW0 : ({t₀} : Set ℝ) ×ˢ closedBall (0 : 𝔼) 3 ⊆
      {q : ℝ × 𝔼 | q ∈ Ioo (t₀ - δ₁) (t₀ + δ₁) ×ˢ (univ : Set 𝔼) ∧
        jp N (fun y => Ψ (q.1, y)) q.2 ∈ admSet n N 2} := by
    rintro ⟨t, x⟩ ⟨ht, hx⟩
    have ht' : t = t₀ := ht
    subst ht'
    refine ⟨⟨ht₀T, mem_univ _⟩, ?_⟩
    intro m hm
    have hx' : x ∈ closedBall (0 : 𝔼) 3 := hx
    have hx4 : x ∈ ball (0 : 𝔼) 4 := closedBall_subset_ball (by norm_num) hx'
    have hev : (fun y => Ψ (t, y)) =ᶠ[𝓝 x] φ₀ := by
      filter_upwards [isOpen_ball.mem_nhds hx4] with y hy using hΨ0 y hy
    show ‖jet (m : ℕ) (fun y => Ψ (t, y)) x‖ < 2
    rw [jet_congr_of_eventuallyEq hev]
    exact cge_jet_bounds hc₀ hc₀1 hK₃ hφ₀ h0 (fun x hx => ⟨(hgr x hx).2.1, (hgr x hx).2.2⟩) hK x
      hx' m hm
  obtain ⟨u, v, hu, hv, htu, hSv, huv⟩ := generalized_tube_lemma (isCompact_singleton (x := t₀))
    (isCompact_closedBall (0 : 𝔼) 3) hAopen hW0
  obtain ⟨ε, hε, hεu⟩ := Metric.isOpen_iff.1 hu t₀ (htu rfl)
  have hδpos : 0 < min ε δ₁ := lt_min hε hδ₁
  have hTu : Ioo (t₀ - min ε δ₁) (t₀ + min ε δ₁) ⊆ u := fun t ht =>
    hεu (by
      rw [Real.ball_eq_Ioo]
      exact ⟨by linarith [min_le_left ε δ₁, ht.1], by linarith [min_le_left ε δ₁, ht.2]⟩)
  have hTT₁ : Ioo (t₀ - min ε δ₁) (t₀ + min ε δ₁) ⊆ Ioo (t₀ - δ₁) (t₀ + δ₁) :=
    Ioo_subset_Ioo (by linarith [min_le_right ε δ₁]) (by linarith [min_le_right ε δ₁])
  refine ⟨min ε δ₁, hδpos, hTT₁.trans hsub, ξ, Ψ, ⟨isOpen_Ioo, hΨ.mono (prod_mono hTT₁ subset_rfl),
    fun t ht x hx => ?_, fun t ht x hx => heqn t (hTT₁ ht) x hx⟩, hξ.mono (prod_mono hTT₁ subset_rfl),
    hΨ0, hξ0, fun t ht => ?_⟩
  · have hmem : (t, x) ∈ u ×ˢ v := ⟨hTu ht, hSv (ball_subset_closedBall hx)⟩
    exact (huv hmem).2
  · obtain ⟨hinj', hder, hpos, hleft⟩ := hper t (hTT₁ ht)
    exact ⟨hinj', hder, fun x hx => ⟨(hpos x hx).1, (hpos x hx).2.1⟩,
      fun x hx => hgrid t (hTT₁ ht) x hx, hleft⟩

end Lemmas

/-! ### Part 5: the theorem -/

theorem exists_chart_graphEvol {n : ℕ} {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    [FiniteDimensional ℝ E]
    {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
    {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V] [FiniteDimensional ℝ V]
    (ι : EuclideanSpace ℝ (Fin n) ≃L[ℝ] E) {J : V →L[ℝ] V}
    (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    (hdim : Module.finrank ℝ V = 2 * n) (N : ℕ) (hN : 4 ≤ N)
    {c₀ K₃ : ℝ} (hc₀ : 0 ≤ c₀) (hc₀1 : c₀ ≤ 1 / 100) (hK₃ : K₃ ≤ 1)
    {I : Set ℝ} {Φ : ℝ → M → V} (hΦ : IsFlowM ι J I Φ) {t₀ : ℝ} (ht₀ : t₀ ∈ I)
    {p : M} {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} {σ : EuclideanSpace ℝ (Fin n) → M}
    {φ₀ : EuclideanSpace ℝ (Fin n) → ℝ} (hch : IsGraphChart ι J (Φ t₀) c₀ p U σ φ₀)
    (hK : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6, ‖jet 3 φ₀ x‖ ≤ K₃) :
    ∃ δ > 0, Ioo (t₀ - δ) (t₀ + δ) ⊆ I ∧
      ∃ (ξ : ℝ × EuclideanSpace ℝ (Fin n) → EuclideanSpace ℝ (Fin n))
        (Ψ : ℝ × EuclideanSpace ℝ (Fin n) → ℝ),
        GraphEvol N 2 3 (Ioo (t₀ - δ) (t₀ + δ)) Ψ ∧
        ContDiffOn ℝ ∞ ξ (Ioo (t₀ - δ) (t₀ + δ) ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) 5) ∧
        (∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 4, Ψ (t₀, x) = φ₀ x) ∧
        (∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 5, ξ (t₀, x) = x) ∧
        ∀ t ∈ Ioo (t₀ - δ) (t₀ + δ),
          InjOn (fun y => horC U (Φ t₀ p) (Φ t (σ y))) (ball (0 : EuclideanSpace ℝ (Fin n)) (11 / 2)) ∧
          (∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) (11 / 2),
            Function.Injective (fderiv ℝ (fun y' => horC U (Φ t₀ p) (Φ t (σ y'))) y) ∧
            ‖fderiv ℝ (fun y' => horC U (Φ t₀ p) (Φ t (σ y'))) y - ContinuousLinearMap.id ℝ _‖ ≤ 1 / 4) ∧
          (∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 5,
            ξ (t, x) ∈ ball (0 : EuclideanSpace ℝ (Fin n)) (11 / 2) ∧
            horC U (Φ t₀ p) (Φ t (σ (ξ (t, x)))) = x) ∧
          (∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 4,
            Φ t (σ (ξ (t, x))) = graphX J U (Φ t₀ p) (fun x' => Ψ (t, x')) x) ∧
          (∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) (11 / 2),
            horC U (Φ t₀ p) (Φ t (σ y)) ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 5 →
              ξ (t, horC U (Φ t₀ p) (Φ t (σ y))) = y) := by
  obtain ⟨ε, hε, hεI⟩ := Metric.isOpen_iff.1 hΦ.isOpen t₀ ht₀
  have hab : Ioo (t₀ - ε) (t₀ + ε) ⊆ I := by
    rw [← Real.ball_eq_Ioo]
    exact hεI
  have ht₀ab : t₀ ∈ Ioo (t₀ - ε) (t₀ + ε) := ⟨by linarith, by linarith⟩
  have hY : ContDiffOn ℝ ∞ (fun q : ℝ × EuclideanSpace ℝ (Fin n) => Φ q.1 (σ q.2))
      (Ioo (t₀ - ε) (t₀ + ε) ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) 6) :=
    (cge_contDiffOn_family ι hΦ.smooth hch.param hΦ.isOpen).mono (prod_mono hab subset_rfl)
  have himm : ∀ t ∈ Ioo (t₀ - ε) (t₀ + ε), ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6,
      Function.Injective (fderiv ℝ (Φ t ∘ σ) y) ∧
        ∀ u w, ⟪J (fderiv ℝ (Φ t ∘ σ) y u), fderiv ℝ (Φ t ∘ σ) y w⟫ = 0 := fun t ht y hy =>
    ⟨((hΦ.imm t (hab ht)).comp_param ι hch.param).2.1 y hy,
      ((hΦ.imm t (hab ht)).comp_param ι hch.param).2.2 y hy⟩
  obtain ⟨δ, hδ, hsub, ξ, Ψ, h1, h2, h3, h4, h5⟩ := cge_core hJi hJJ hdim N hN hc₀ hc₀1 hK₃
    hch.frame (Φ t₀ p) ht₀ab (Y := fun q : ℝ × EuclideanSpace ℝ (Fin n) => Φ q.1 (σ q.2)) hY
    (fun t ht y hy => (himm t ht y hy).1) (fun t ht y hy => (himm t ht y hy).2)
    (fun t ht y hy => cge_hasDerivAt_family ι hΦ hch.param (hab ht) hy)
    (fun y hy => (hch.graph y hy).2.1) hch.smooth hch.zero
    (fun x hx => ⟨(hch.graph x hx).1, (hch.graph x hx).2.2.1, (hch.graph x hx).2.2.2⟩) hK
  exact ⟨δ, hδ, hsub.trans hab, ξ, Ψ, h1, h2, h3, h4, h5⟩

/-- info: 'CWFlow.Extension.exists_chart_graphEvol' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_chart_graphEvol

end CWFlow.Extension
