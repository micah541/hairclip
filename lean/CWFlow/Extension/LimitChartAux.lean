import CWFlow.Extension.GraphTower
import CWFlow.Extension.GaugeTerm
import CWFlow.Extension.SliceSmooth
import CWFlow.Extension.ChartBounds
import CWFlow.Extension.TimeGraph
import CWFlow.Extension.IntM
import CWFlow.Interp.TowerBound
import CWFlow.Interp.OdeTower
import CWFlow.Interp.BallInterp

/-!
# Auxiliary lemmas for the limit chart

Generic ingredients for `exists_limit_chart`: restriction of a parametrization, constants do not
change the potential, the velocity or the higher jets of a graph function, the shift-and-bound
step for the graph function (`lc_shift_bounds`) and the transfer of derivative bounds from the
ODE for the horizontal coordinate to the graph map (`lc_ode_comp`).
-/

open Metric Set
open CWFlow.Interp
open scoped ContDiff Manifold Topology RealInnerProductSpace

namespace CWFlow.Extension

section Param

theorem lc_isParam_mono {n : ℕ} {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    [FiniteDimensional ℝ E] {M : Type*} [TopologicalSpace M] [ChartedSpace E M]
    [IsManifold 𝓘(ℝ, E) ∞ M] (ι : EuclideanSpace ℝ (Fin n) ≃L[ℝ] E)
    {B B' : Set (EuclideanSpace ℝ (Fin n))} {σ : EuclideanSpace ℝ (Fin n) → M}
    (h : IsParam ι B σ) (hB' : IsOpen B') (hsub : B' ⊆ B) : IsParam ι B' σ := by
  obtain ⟨κ, hκc, hκ⟩ := h.inv
  refine ⟨hB', h.injOn.mono hsub, ?_, h.cont.mono hsub,
    ⟨κ, hκc.mono (Set.image_mono hsub), fun y hy => hκ y (hsub hy)⟩, fun z => ?_, fun z y hy => ?_⟩
  · have e : σ '' B' = σ '' B ∩ κ ⁻¹' B' := by
      ext m
      constructor
      · rintro ⟨y, hy, rfl⟩
        refine ⟨⟨y, hsub hy, rfl⟩, ?_⟩
        show κ (σ y) ∈ B'
        rw [hκ y (hsub hy)]
        exact hy
      · rintro ⟨⟨y, hy, rfl⟩, hm⟩
        refine ⟨y, ?_, rfl⟩
        have : κ (σ y) ∈ B' := hm
        rwa [hκ y hy] at this
    rw [e]
    exact hκc.isOpen_inter_preimage h.isOpen_image hB'
  · exact (h.trans_smooth z).mono (inter_subset_inter_left _ hsub)
  · exact h.trans_det z y ⟨hsub hy.1, hy.2⟩

end Param

section Const

variable {n : ℕ} {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]

theorem lc_gradE_sub_const (φ : EuclideanSpace ℝ (Fin n) → ℝ) (c : ℝ) :
    gradE (fun y => φ y - c) = gradE φ := by
  funext x
  unfold gradE
  simp only [fderiv_sub_const]

theorem lc_graphX_sub_const (J : V →L[ℝ] V) (U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V) (y₀ : V)
    (φ : EuclideanSpace ℝ (Fin n) → ℝ) (c : ℝ) :
    graphX J U y₀ (fun y => φ y - c) = graphX J U y₀ φ := by
  funext x
  unfold graphX
  rw [lc_gradE_sub_const]

theorem lc_jp_congr {N : ℕ} {φ ψ : EuclideanSpace ℝ (Fin n) → ℝ} {x : EuclideanSpace ℝ (Fin n)}
    (h : φ =ᶠ[𝓝 x] ψ) : jp N φ x = jp N ψ x :=
  Prod.ext rfl (funext fun m => jet_congr_of_eventuallyEq h _)

theorem lc_exists_adm {N : ℕ} (p : P n N) : ∃ K, p ∈ admSet n N K := by
  refine ⟨1 + ∑ m : Fin (N + 1), ‖p.2 m‖, fun m _ => ?_⟩
  have := Finset.single_le_sum (f := fun m : Fin (N + 1) => ‖p.2 m‖) (fun m _ => norm_nonneg _)
    (Finset.mem_univ m)
  linarith

theorem lc_potSym_sub_const_global {J : V →L[ℝ] V} (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫)
    (hJJ : ∀ v : V, J (J v) = -v) {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} (hU : IsLagFrame J U)
    {ψ : EuclideanSpace ℝ (Fin n) → ℝ} (hψ : ContDiff ℝ ∞ ψ) (c : ℝ) {N : ℕ} (hN : 4 ≤ N)
    (x : EuclideanSpace ℝ (Fin n)) :
    potSym (jp N (fun y => ψ y - c) x) = potSym (jp N ψ x) := by
  rw [← potX_graphX hJi hJJ hU 0 (hψ.sub contDiff_const) hN x,
    ← potX_graphX hJi hJJ hU 0 hψ hN x, lc_graphX_sub_const]

theorem lc_vSym_sub_const_global {J : V →L[ℝ] V} (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫)
    (hJJ : ∀ v : V, J (J v) = -v) {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} (hU : IsLagFrame J U)
    {ψ : EuclideanSpace ℝ (Fin n) → ℝ} (hψ : ContDiff ℝ ∞ ψ) (c : ℝ) {N : ℕ} (hN : 5 ≤ N)
    (x : EuclideanSpace ℝ (Fin n)) (a : Fin n) :
    vSym (jp N (fun y => ψ y - c) x) a = vSym (jp N ψ x) a := by
  rw [← horCL_velX_graphX hJi hJJ hU 0 (hψ.sub contDiff_const) hN x (lc_exists_adm _) a,
    ← horCL_velX_graphX hJi hJJ hU 0 hψ hN x (lc_exists_adm _) a, lc_graphX_sub_const]

theorem lc_potSym_sub_const {J : V →L[ℝ] V} (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫)
    (hJJ : ∀ v : V, J (J v) = -v) {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} (hU : IsLagFrame J U)
    {s : Set (EuclideanSpace ℝ (Fin n))} (hs : IsOpen s) {φ : EuclideanSpace ℝ (Fin n) → ℝ}
    (hφ : ContDiffOn ℝ ∞ φ s) {x : EuclideanSpace ℝ (Fin n)} (hx : x ∈ s) (c : ℝ) {N : ℕ}
    (hN : 4 ≤ N) : potSym (jp N (fun y => φ y - c) x) = potSym (jp N φ x) := by
  obtain ⟨ψ, hψ, hev⟩ := cb_exists_global_eventually hs hφ hx
  have hev' : (fun y => φ y - c) =ᶠ[𝓝 x] fun y => ψ y - c := by
    filter_upwards [hev] with y hy
    simp only [hy]
  rw [lc_jp_congr hev', lc_jp_congr hev]
  exact lc_potSym_sub_const_global hJi hJJ hU hψ c hN x

theorem lc_vSym_sub_const {J : V →L[ℝ] V} (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫)
    (hJJ : ∀ v : V, J (J v) = -v) {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} (hU : IsLagFrame J U)
    {s : Set (EuclideanSpace ℝ (Fin n))} (hs : IsOpen s) {φ : EuclideanSpace ℝ (Fin n) → ℝ}
    (hφ : ContDiffOn ℝ ∞ φ s) {x : EuclideanSpace ℝ (Fin n)} (hx : x ∈ s) (c : ℝ) {N : ℕ}
    (hN : 5 ≤ N) (a : Fin n) : vSym (jp N (fun y => φ y - c) x) a = vSym (jp N φ x) a := by
  obtain ⟨ψ, hψ, hev⟩ := cb_exists_global_eventually hs hφ hx
  have hev' : (fun y => φ y - c) =ᶠ[𝓝 x] fun y => ψ y - c := by
    filter_upwards [hev] with y hy
    simp only [hy]
  rw [lc_jp_congr hev', lc_jp_congr hev]
  exact lc_vSym_sub_const_global hJi hJJ hU hψ c hN x a

theorem lc_jet_sub_const_global {ψ : EuclideanSpace ℝ (Fin n) → ℝ} (hψ : ContDiff ℝ ∞ ψ) (c : ℝ)
    {m : ℕ} (hm : 1 ≤ m) (x : EuclideanSpace ℝ (Fin n)) :
    jet m (fun y => ψ y - c) x = jet m ψ x := by
  have key : iteratedFDeriv ℝ m (fun y => ψ y - c) x = iteratedFDeriv ℝ m ψ x := by
    have e : (fun y => ψ y - c) = ψ - fun _ => c := rfl
    rw [e, iteratedFDeriv_sub_apply (hψ.contDiffAt.of_le (by exact_mod_cast le_top))
      (contDiffAt_const.of_le (by exact_mod_cast le_top)),
      iteratedFDeriv_const_of_ne (by omega)]
    simp
  unfold jet
  rw [key]

theorem lc_jet_sub_const {s : Set (EuclideanSpace ℝ (Fin n))} (hs : IsOpen s)
    {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiffOn ℝ ∞ φ s) {x : EuclideanSpace ℝ (Fin n)}
    (hx : x ∈ s) (c : ℝ) {m : ℕ} (hm : 1 ≤ m) :
    jet m (fun y => φ y - c) x = jet m φ x := by
  obtain ⟨ψ, hψ, hev⟩ := cb_exists_global_eventually hs hφ hx
  have hev' : (fun y => φ y - c) =ᶠ[𝓝 x] fun y => ψ y - c := by
    filter_upwards [hev] with y hy
    simp only [hy]
  rw [jet_congr_of_eventuallyEq hev', jet_congr_of_eventuallyEq hev]
  exact lc_jet_sub_const_global hψ c hm x

end Const

section SmoothBdd

theorem lc_smoothBdd_congr {X Y : Type*} [NormedAddCommGroup X] [NormedSpace ℝ X]
    [NormedAddCommGroup Y] [NormedSpace ℝ Y] {S : Set X} (hS : IsOpen S) {f g : X → Y}
    (h : EqOn f g S) (hf : SmoothBdd S f) : SmoothBdd S g := by
  refine ⟨hf.1.congr (fun x hx => (h hx).symm), fun m => ?_⟩
  obtain ⟨C, hC⟩ := hf.2 m
  refine ⟨C, fun x hx => ?_⟩
  have hev : f =ᶠ[𝓝 x] g := by
    filter_upwards [hS.mem_nhds hx] with y hy using h hy
  rw [← (Filter.EventuallyEq.iteratedFDeriv ℝ hev m).eq_of_nhds]
  exact hC x hx

theorem lc_smoothBdd_comp {X Y Z : Type*} [NormedAddCommGroup X] [NormedSpace ℝ X]
    [NormedAddCommGroup Y] [NormedSpace ℝ Y] [NormedAddCommGroup Z] [NormedSpace ℝ Z]
    {S : Set X} {S' : Set Y} (hS : IsOpen S) (hS' : IsOpen S') {g : Y → Z} {f : X → Y}
    (hg : SmoothBdd S' g) (hf : SmoothBdd S f) (hmaps : MapsTo f S S') :
    SmoothBdd S (fun x => g (f x)) := by
  refine ⟨hg.1.comp hf.1 hmaps, fun m => ?_⟩
  choose Cg hCg using hg.2
  choose Cf hCf using hf.2
  set C : ℝ := ∑ i ∈ Finset.range (m + 1), |Cg i| with hC
  set D : ℝ := 1 + ∑ i ∈ Finset.range (m + 1), |Cf i| with hD
  have hD1 : 1 ≤ D := by
    have : 0 ≤ ∑ i ∈ Finset.range (m + 1), |Cf i| := Finset.sum_nonneg fun i _ => abs_nonneg _
    linarith
  refine ⟨m.factorial * C * D ^ m, fun x hx => ?_⟩
  rw [← iteratedFDerivWithin_of_isOpen m hS hx]
  have := norm_iteratedFDerivWithin_comp_le (g := g) (f := f) (n := m) (s := S) (t := S') (x := x)
    (N := ∞) hg.1 hf.1 (tw_le_infty m) hS'.uniqueDiffOn hS.uniqueDiffOn hmaps hx (C := C) (D := D)
    (fun i hi => by
      rw [iteratedFDerivWithin_of_isOpen i hS' (hmaps hx)]
      refine (hCg i _ (hmaps hx)).trans ?_
      refine (le_abs_self _).trans ?_
      exact Finset.single_le_sum (f := fun i => |Cg i|) (fun i _ => abs_nonneg _)
        (Finset.mem_range.2 (by omega)))
    (fun i hi1 him => by
      rw [iteratedFDerivWithin_of_isOpen i hS hx]
      refine (hCf i x hx).trans ?_
      refine (le_abs_self _).trans ?_
      have h1 : |Cf i| ≤ D := by
        have := Finset.single_le_sum (f := fun i => |Cf i|) (fun i _ => abs_nonneg _)
          (Finset.mem_range.2 (by omega : i < m + 1))
        have h' : |Cf i| ≤ ∑ i ∈ Finset.range (m + 1), |Cf i| := this
        linarith
      exact h1.trans (le_self_pow₀ hD1 (by omega)))
  exact this

end SmoothBdd

section Bounds

variable {n : ℕ}

theorem lc_exists_bound_jet {W : Type*} [NormedAddCommGroup W] [InnerProductSpace ℝ W]
    [FiniteDimensional ℝ W] {O : Set (EuclideanSpace ℝ (Fin n))} (hO : IsOpen O) {R : ℝ}
    (hR : 0 < R) (hRO : closedBall (0 : EuclideanSpace ℝ (Fin n)) R ⊆ O)
    {f : EuclideanSpace ℝ (Fin n) → W} (hf : ContDiffOn ℝ ∞ f O) (j : ℕ) :
    ∃ C : ℝ, ∀ y ∈ closedBall (0 : EuclideanSpace ℝ (Fin n)) R, ‖jet j f y‖ ≤ C := by
  obtain ⟨R', hR', g, hg, -, hgf⟩ := exists_global_of_contDiffOn hO hR hRO hf
  obtain ⟨C, hC⟩ := (isCompact_closedBall (0 : EuclideanSpace ℝ (Fin n)) R).exists_bound_of_continuousOn
    (f := jet j g) (contDiff_jet hg j).continuous.continuousOn
  refine ⟨C, fun y hy => ?_⟩
  have hev : f =ᶠ[𝓝 y] g := by
    filter_upwards [isOpen_ball.mem_nhds (closedBall_subset_ball hR' hy)] with z hz using
      (hgf z hz).symm
  rw [jet_congr_of_eventuallyEq hev]
  exact hC y hy

theorem lc_hasDerivAt_euclid {v : ℝ → EuclideanSpace ℝ (Fin n)} {w : EuclideanSpace ℝ (Fin n)}
    {t : ℝ} (h : ∀ c, HasDerivAt (fun s => v s c) (w c) t) : HasDerivAt v w t := by
  have hg : HasDerivAt (fun s => fun c => v s c) (fun c => w c) t := hasDerivAt_pi.2 h
  exact (PiLp.hasFDerivAt_toLp (𝕜 := ℝ) 2 (fun c : Fin n => v t c)).comp_hasDerivAt t hg

theorem lc_gradE_slice {F : ℝ × EuclideanSpace ℝ (Fin n) → ℝ} {t : ℝ}
    {x : EuclideanSpace ℝ (Fin n)} (hd : DifferentiableAt ℝ F (t, x)) :
    gradE (fun y => F (t, y)) x =
      WithLp.toLp 2 (fun c => fderiv ℝ F (t, x) (0, EuclideanSpace.single c (1 : ℝ))) := by
  have h2 : HasFDerivAt (fun y => F (t, y))
      ((fderiv ℝ F (t, x)).comp (ContinuousLinearMap.inr ℝ ℝ (EuclideanSpace ℝ (Fin n)))) x :=
    hd.hasFDerivAt.comp x (hasFDerivAt_prodMk_right t x)
  unfold gradE
  rw [h2.fderiv]
  rfl

end Bounds

section Shift

variable {n : ℕ}

theorem lc_contDiffOn_vSym {a b r : ℝ} {F : ℝ × EuclideanSpace ℝ (Fin n) → ℝ}
    (hF : ContDiffOn ℝ ∞ F (Ioo a b ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) r))
    (hlow : ∀ q ∈ Ioo a b ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) r, ∀ m ≤ 3,
      ‖jet m (fun y => F (q.1, y)) q.2‖ < 2) (c : Fin n) :
    ContDiffOn ℝ ∞ (fun q : ℝ × EuclideanSpace ℝ (Fin n) =>
      vSym (jp 5 (fun y => F (q.1, y)) q.2) c)
      (Ioo a b ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) r) := by
  intro q hq
  obtain ⟨G, ρ, hρ, hsub, hG, hGF⟩ := gtw_loc hF hq.2
  have hsym := isSymbol_vSym (isOpen_admSet (n := n) (N := 5) 2) (admSet_low 2) (le_refl 5) c
  have hGsm := contDiffOn_symbol_slice hsym isOpen_Ioo hG
  have hopen := isOpen_symbol_slice_dom (N := 5) (isOpen_admSet (n := n) (N := 5) 2) isOpen_Ioo hG
  have hjp : ∀ q' ∈ Ioo a b ×ˢ ball q.2 ρ,
      jp 5 (fun y => F (q'.1, y)) q'.2 = jp 5 (fun y => G (q'.1, y)) q'.2 := fun q' hq' =>
    gtw_jp_congr 5 isOpen_ball (fun y hy => (hGF q'.1 hq'.1 y hy).symm) hq'.2
  have hqD : q ∈ {q' : ℝ × EuclideanSpace ℝ (Fin n) |
      q' ∈ Ioo a b ×ˢ (univ : Set (EuclideanSpace ℝ (Fin n))) ∧
        jp 5 (fun y => G (q'.1, y)) q'.2 ∈ admSet n 5 2} := by
    refine ⟨⟨hq.1, mem_univ _⟩, ?_⟩
    rw [← hjp q ⟨hq.1, mem_ball_self hρ⟩]
    exact fun m hm => hlow q hq (m : ℕ) hm
  have hat : ContDiffAt ℝ ∞ (fun q' : ℝ × EuclideanSpace ℝ (Fin n) =>
      vSym (jp 5 (fun y => G (q'.1, y)) q'.2) c) q := hGsm.contDiffAt (hopen.mem_nhds hqD)
  have hev : (fun q' : ℝ × EuclideanSpace ℝ (Fin n) => vSym (jp 5 (fun y => F (q'.1, y)) q'.2) c)
      =ᶠ[𝓝 q] fun q' => vSym (jp 5 (fun y => G (q'.1, y)) q'.2) c := by
    filter_upwards [(isOpen_Ioo.prod (isOpen_ball (x := q.2) (ε := ρ))).mem_nhds
      (show q ∈ Ioo a b ×ˢ ball q.2 ρ from ⟨hq.1, mem_ball_self hρ⟩)] with q' hq'
    rw [hjp q' hq']
  exact (hat.congr_of_eventuallyEq hev).contDiffWithinAt

theorem lc_smoothBdd_V {a b r : ℝ} {F : ℝ × EuclideanSpace ℝ (Fin n) → ℝ}
    (hF : ContDiffOn ℝ ∞ F (Ioo a b ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) r))
    (hlow : ∀ q ∈ Ioo a b ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) r, ∀ m ≤ 3,
      ‖jet m (fun y => F (q.1, y)) q.2‖ < 2)
    (hVb : ∀ (k : ℕ) (c : Fin n), ∃ C : ℝ,
      ∀ q ∈ Ioo a b ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) r,
      ‖iteratedFDeriv ℝ k (fun q : ℝ × EuclideanSpace ℝ (Fin n) =>
        vSym (jp 5 (fun y => F (q.1, y)) q.2) c) q‖ ≤ C) :
    SmoothBdd (Ioo a b ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) r)
      (fun q : ℝ × EuclideanSpace ℝ (Fin n) =>
        (WithLp.toLp 2 (fun c : Fin n => vSym (jp 5 (fun y => F (q.1, y)) q.2) c) :
          EuclideanSpace ℝ (Fin n))) :=
  ot_smoothBdd_piLp (isOpen_Ioo.prod isOpen_ball) fun I =>
    ⟨lc_contDiffOn_vSym hF hlow I, fun k => hVb k I⟩

theorem lc_smoothBdd_graph {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V] {a b r : ℝ}
    {F : ℝ × EuclideanSpace ℝ (Fin n) → ℝ}
    (hF : ContDiffOn ℝ ∞ F (Ioo a b ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) r))
    (hb : ∀ k : ℕ, ∃ C : ℝ, ∀ q ∈ Ioo a b ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) r,
      ‖iteratedFDeriv ℝ k F q‖ ≤ C)
    (J : V →L[ℝ] V) (U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V) (y₀ : V) :
    SmoothBdd (Ioo a b ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) r)
      (fun q : ℝ × EuclideanSpace ℝ (Fin n) =>
        y₀ + U q.2 + J (U (gradE (fun y => F (q.1, y)) q.2))) := by
  set S := Ioo a b ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) r with hSdef
  have hS : IsOpen S := isOpen_Ioo.prod isOpen_ball
  have hSb : Bornology.IsBounded S := (Metric.isBounded_Ioo a b).prod Metric.isBounded_ball
  have hD : SmoothBdd S (fderiv ℝ F) := tw_smoothBdd_fderiv hS ⟨hF, hb⟩
  have hc : ∀ c : Fin n, SmoothBdd S (fun q => fderiv ℝ F q
      ((0 : ℝ), EuclideanSpace.single c (1 : ℝ))) := fun c =>
    tw_smoothBdd_clm_left hS (ContinuousLinearMap.apply ℝ ℝ
      ((0 : ℝ), EuclideanSpace.single c (1 : ℝ))) hD
  have hv : SmoothBdd S (fun q : ℝ × EuclideanSpace ℝ (Fin n) =>
      (WithLp.toLp 2 (fun c : Fin n => fderiv ℝ F q
        ((0 : ℝ), EuclideanSpace.single c (1 : ℝ))) : EuclideanSpace ℝ (Fin n))) :=
    ot_smoothBdd_piLp hS (fun I => hc I)
  have h1 : SmoothBdd S (fun q : ℝ × EuclideanSpace ℝ (Fin n) => U q.2) :=
    gtw_smoothBdd_clm hSb (U.toContinuousLinearMap.comp
      (ContinuousLinearMap.snd ℝ ℝ (EuclideanSpace ℝ (Fin n))))
  have h2 : SmoothBdd S (fun q : ℝ × EuclideanSpace ℝ (Fin n) =>
      J (U (WithLp.toLp 2 (fun c : Fin n => fderiv ℝ F q
        ((0 : ℝ), EuclideanSpace.single c (1 : ℝ))) : EuclideanSpace ℝ (Fin n)))) :=
    tw_smoothBdd_clm_left hS (J.comp U.toContinuousLinearMap) hv
  have h3 := tw_smoothBdd_add hS (tw_smoothBdd_add hS (tw_smoothBdd_const S y₀) h1) h2
  refine lc_smoothBdd_congr hS (fun q hq => ?_) h3
  obtain ⟨t, x⟩ := q
  have hd : DifferentiableAt ℝ F (t, x) :=
    (hF.contDiffAt (hS.mem_nhds hq)).differentiableAt (by simp)
  show y₀ + U x + J (U (WithLp.toLp 2 (fun c : Fin n => fderiv ℝ F (t, x)
      ((0 : ℝ), EuclideanSpace.single c (1 : ℝ))))) =
    y₀ + U x + J (U (gradE (fun y => F (t, y)) x))
  rw [lc_gradE_slice hd]

theorem lc_ode_comp {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V] {J : V →L[ℝ] V}
    {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} {y₀ : V} {t₂ t₁ T : ℝ} (h21 : t₂ < t₁) (h1T : t₁ < T)
    {F : ℝ × EuclideanSpace ℝ (Fin n) → ℝ}
    (hF : ContDiffOn ℝ ∞ F (Ioo t₂ T ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) 3))
    (hlow : ∀ q ∈ Ioo t₂ T ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) 3, ∀ m ≤ 3,
      ‖jet m (fun y => F (q.1, y)) q.2‖ < 2)
    (hGb : ∀ k : ℕ, ∃ C : ℝ, ∀ q ∈ Ioo t₂ T ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) 3,
      ‖iteratedFDeriv ℝ k F q‖ ≤ C)
    (hVb : ∀ (k : ℕ) (c : Fin n), ∃ C : ℝ,
      ∀ q ∈ Ioo t₂ T ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) 3,
      ‖iteratedFDeriv ℝ k (fun q : ℝ × EuclideanSpace ℝ (Fin n) =>
        vSym (jp 5 (fun y => F (q.1, y)) q.2) c) q‖ ≤ C)
    {χ : ℝ × EuclideanSpace ℝ (Fin n) → EuclideanSpace ℝ (Fin n)}
    (hχ : ContDiffOn ℝ ∞ χ (Ioo t₂ T ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) 2))
    (hmaps : ∀ q ∈ Ioo t₂ T ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) 2,
      (q.1, χ q) ∈ Ioo t₂ T ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) 3)
    (hode : ∀ q ∈ Ioo t₂ T ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) 2, ∀ c : Fin n,
      HasDerivAt (fun s => χ (s, q.2) c) (vSym (jp 5 (fun y => F (q.1, y)) (χ q)) c) q.1)
    (hinit : ∀ j : ℕ, ∃ C : ℝ, ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 2,
      ‖jet j (fun y' => χ (t₁, y')) y‖ ≤ C) :
    SmoothBdd (Ioo t₁ T ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) 2)
      (fun q : ℝ × EuclideanSpace ℝ (Fin n) =>
        y₀ + U (χ q) + J (U (gradE (fun y => F (q.1, y)) (χ q)))) := by
  have hV := lc_smoothBdd_V hF hlow hVb
  have hS3 : IsOpen (Ioo t₂ T ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) 3) :=
    isOpen_Ioo.prod isOpen_ball
  have hχb := ode_iteratedFDeriv_bounded (n := n) (t₁ := t₁) (b := T) h1T
    (B := ball (0 : EuclideanSpace ℝ (Fin n)) 2) isOpen_ball hS3 hV h21 hχ hmaps
    (fun q hq => lc_hasDerivAt_euclid (fun c => hode q hq c)) hinit
  have hSc : IsOpen (Ioo t₁ T ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) 2) :=
    isOpen_Ioo.prod isOpen_ball
  have hScb : Bornology.IsBounded (Ioo t₁ T ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) 2) :=
    (Metric.isBounded_Ioo t₁ T).prod Metric.isBounded_ball
  have hχS : SmoothBdd (Ioo t₁ T ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) 2) χ :=
    ⟨hχ.mono (prod_mono (Ioo_subset_Ioo_left h21.le) subset_rfl), hχb⟩
  have hf : SmoothBdd (Ioo t₁ T ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) 2)
      (fun q : ℝ × EuclideanSpace ℝ (Fin n) => (q.1, χ q)) := by
    have e1 := gtw_smoothBdd_clm hScb
      ((ContinuousLinearMap.inl ℝ ℝ (EuclideanSpace ℝ (Fin n))).comp
        (ContinuousLinearMap.fst ℝ ℝ (EuclideanSpace ℝ (Fin n))))
    have e2 := tw_smoothBdd_clm_left hSc
      (ContinuousLinearMap.inr ℝ ℝ (EuclideanSpace ℝ (Fin n))) hχS
    refine lc_smoothBdd_congr hSc (fun q _ => ?_) (tw_smoothBdd_add hSc e1 e2)
    simp
  have hG := lc_smoothBdd_graph hF hGb J U y₀
  exact lc_smoothBdd_comp hSc hS3 hG hf
    (fun q hq => hmaps q ⟨⟨by linarith [hq.1.1], hq.1.2⟩, hq.2⟩)

theorem lc_shift_bounds {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V] {J : V →L[ℝ] V}
    (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} (hU : IsLagFrame J U) {t₀ T : ℝ} (hT : t₀ < T)
    {φ : ℝ × EuclideanSpace ℝ (Fin n) → ℝ}
    (hφ : ContDiffOn ℝ ∞ φ (Ioo t₀ T ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) 5))
    (hder : ∀ t ∈ Ioo t₀ T, ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 5,
      HasDerivAt (fun s => φ (s, x)) (potSym (jp 4 (fun y => φ (t, y)) x)) t)
    (hg : ∀ t ∈ Ioo t₀ T, ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 5,
      ‖gradE (fun y => φ (t, y)) x‖ ≤ 1 / 2)
    (h2 : ∀ t ∈ Ioo t₀ T, ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 5,
      ‖jet 2 (fun y => φ (t, y)) x‖ < 2)
    (h3 : ∀ t ∈ Ioo t₀ T, ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 5,
      ‖jet 3 (fun y => φ (t, y)) x‖ < 2)
    (hCm : ∀ m : ℕ, ∃ C : ℝ, 1 ≤ m → ∀ t ∈ Ioo t₀ T,
      ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 3, ‖jet m (fun y => φ (t, y)) x‖ ≤ C) :
    ∃ t₂ ∈ Ioo t₀ T, ∃ c₂ : ℝ,
      (∀ q ∈ Ioo t₂ T ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) 3, ∀ m ≤ 3,
        ‖jet m (fun y => φ (q.1, y) - c₂) q.2‖ < 2) ∧
      (∀ k : ℕ, ∃ C : ℝ, ∀ q ∈ Ioo t₂ T ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) 3,
        ‖iteratedFDeriv ℝ k (fun q : ℝ × EuclideanSpace ℝ (Fin n) => φ q - c₂) q‖ ≤ C) ∧
      (∀ (k : ℕ) (c : Fin n), ∃ C : ℝ,
        ∀ q ∈ Ioo t₂ T ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) 3,
        ‖iteratedFDeriv ℝ k (fun q : ℝ × EuclideanSpace ℝ (Fin n) =>
          vSym (jp 5 (fun y => φ (q.1, y) - c₂) q.2) c) q‖ ≤ C) := by
  classical
  have hB35 : ball (0 : EuclideanSpace ℝ (Fin n)) 3 ⊆ ball (0 : EuclideanSpace ℝ (Fin n)) 5 :=
    ball_subset_ball (by norm_num)
  have hsl : ∀ t ∈ Ioo t₀ T, ContDiffOn ℝ ∞ (fun y => φ (t, y))
      (ball (0 : EuclideanSpace ℝ (Fin n)) 5) := fun t ht => tg_slice hφ ht
  choose Cm hCm' using hCm
  set Cj : ℕ → ℝ := fun m => max (Cm m) 0 with hCj
  -- low jets of the shifted slices
  have hlowc : ∀ t ∈ Ioo t₀ T, ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 5, ∀ c : ℝ,
      |φ (t, x) - c| < 2 → ∀ m : ℕ, m ≤ 3 → ‖jet m (fun y => φ (t, y) - c) x‖ < 2 := by
    intro t ht x hx c hc m hm
    interval_cases m
    · rw [norm_jet_zero, Real.norm_eq_abs]
      exact hc
    · rw [lc_jet_sub_const isOpen_ball (hsl t ht) hx c le_rfl, cb_norm_jet_one]
      linarith [hg t ht x hx]
    · rw [lc_jet_sub_const isOpen_ball (hsl t ht) hx c (by norm_num)]
      exact h2 t ht x hx
    · rw [lc_jet_sub_const isOpen_ball (hsl t ht) hx c (by norm_num)]
      exact h3 t ht x hx
  -- the drift of the graph function is bounded
  have hsym := isSymbol_potSym (isOpen_admSet (n := n) (N := 4) 2) (admSet_low 2) (le_refl 4)
  have hsb := gtw_smoothBdd_symbol hsym (fun p hp => hp.1) (gtw_bdd_Sp n 4 3 Cj)
  obtain ⟨Cp, hCp⟩ := hsb.2 0
  have hdrift : ∀ t ∈ Ioo t₀ T, ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 3,
      |potSym (jp 4 (fun y => φ (t, y)) x)| ≤ Cp := by
    intro t ht x hx
    have hx5 := hB35 hx
    have hc0 : |φ (t, x) - φ (t, x)| < 2 := by simp
    have hmem : jp 4 (fun y => φ (t, y) - φ (t, x)) x ∈ gtwSp n 4 3 Cj := by
      refine ⟨fun m hm => hlowc t ht x hx5 _ hc0 (m : ℕ) hm, ?_, fun m => ?_⟩
      · exact mem_ball_zero_iff.1 hx
      · show ‖jet (m : ℕ) (fun y => φ (t, y) - φ (t, x)) x‖ < Cj m + 1
        by_cases hm0 : (m : ℕ) = 0
        · rw [hm0, norm_jet_zero]
          have h0 : 0 ≤ Cj 0 := le_max_right _ _
          simp only [sub_self, norm_zero]
          linarith
        · rw [lc_jet_sub_const isOpen_ball (hsl t ht) hx5 _ (by omega)]
          have h1 := hCm' (m : ℕ) (by omega) t ht x hx
          have h2' : Cm m ≤ Cj m := le_max_left _ _
          linarith
    have := hCp _ hmem
    rw [norm_iteratedFDeriv_zero, Real.norm_eq_abs] at this
    rw [← lc_potSym_sub_const hJi hJJ hU isOpen_ball (hsl t ht) hx5 (φ (t, x)) (le_refl 4)]
    exact this
  -- the choice of `t₂`
  set M : ℝ := |Cp| with hM
  have hM0 : 0 ≤ M := abs_nonneg _
  obtain ⟨δ, hδ0, hδ1, hδ2⟩ : ∃ δ : ℝ, 0 < δ ∧ δ ≤ (T - t₀) / 2 ∧ M * δ ≤ 1 / 4 := by
    refine ⟨min ((T - t₀) / 2) (1 / (4 * (M + 1))), lt_min (by linarith) (by positivity),
      min_le_left _ _, ?_⟩
    calc M * min ((T - t₀) / 2) (1 / (4 * (M + 1))) ≤ M * (1 / (4 * (M + 1))) :=
          mul_le_mul_of_nonneg_left (min_le_right _ _) hM0
      _ ≤ 1 / 4 := by
          rw [mul_one_div, div_le_div_iff₀ (by positivity) (by norm_num)]
          linarith
  set t₂ : ℝ := T - δ with ht₂def
  have ht₂ : t₂ ∈ Ioo t₀ T := ⟨by linarith, by linarith⟩
  have hsub : ∀ t ∈ Ioo t₂ T, t ∈ Ioo t₀ T := fun t ht => ⟨by linarith [ht.1, ht₂.1], ht.2⟩
  have htime : ∀ t ∈ Ioo t₂ T, |φ (t, 0) - φ (t₂, 0)| ≤ 1 / 4 := by
    intro t ht
    have h := (convex_Ioo t₀ T).norm_image_sub_le_of_norm_hasDerivWithin_le
      (f := fun s => φ (s, 0)) (f' := fun s => potSym (jp 4 (fun y => φ (s, y)) 0)) (C := M)
      (fun s hs => (hder s hs 0 (mem_ball_self (by norm_num))).hasDerivWithinAt)
      (fun s hs => by
        rw [Real.norm_eq_abs]
        exact (hdrift s hs 0 (mem_ball_self (by norm_num))).trans (le_abs_self _))
      ht₂ (hsub t ht)
    rw [Real.norm_eq_abs, Real.norm_eq_abs] at h
    calc |φ (t, 0) - φ (t₂, 0)| ≤ M * |t - t₂| := h
      _ ≤ M * δ := by
          refine mul_le_mul_of_nonneg_left ?_ hM0
          rw [abs_of_pos (by linarith [ht.1])]
          linarith [ht.2]
      _ ≤ 1 / 4 := hδ2
  have hspace : ∀ t ∈ Ioo t₀ T, ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 3,
      |φ (t, x) - φ (t, 0)| < 3 / 2 := by
    intro t ht x hx
    have h := (convex_ball (0 : EuclideanSpace ℝ (Fin n)) 3).norm_image_sub_le_of_norm_fderiv_le
      (f := fun y => φ (t, y)) (C := 1 / 2)
      (fun y hy => ((hsl t ht).contDiffAt (isOpen_ball.mem_nhds (hB35 hy))).differentiableAt
        (by simp))
      (fun y hy => (cb_norm_fderiv_le _ _).trans (hg t ht y (hB35 hy)))
      (mem_ball_self (by norm_num)) hx
    simp only [Real.norm_eq_abs, sub_zero] at h
    have hxn : ‖x‖ < 3 := by simpa using hx
    linarith
  set c₂ : ℝ := φ (t₂, 0) with hc₂
  have hsm : ∀ t ∈ Ioo t₂ T, ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 3, |φ (t, x) - c₂| < 2 := by
    intro t ht x hx
    have e1 := hspace t (hsub t ht) x hx
    have e2 := htime t ht
    calc |φ (t, x) - c₂| = |(φ (t, x) - φ (t, 0)) + (φ (t, 0) - φ (t₂, 0))| := by
          congr 1; ring
      _ ≤ |φ (t, x) - φ (t, 0)| + |φ (t, 0) - φ (t₂, 0)| := abs_add_le _ _
      _ < 2 := by linarith
  have hS3 : Ioo t₂ T ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) 3 ⊆
      Ioo t₀ T ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) 5 :=
    prod_mono (Ioo_subset_Ioo_left ht₂.1.le) hB35
  have hF : ContDiffOn ℝ ∞ (fun q : ℝ × EuclideanSpace ℝ (Fin n) => φ q - c₂)
      (Ioo t₂ T ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) 3) :=
    (hφ.sub contDiffOn_const).mono hS3
  have heq : ∀ q ∈ Ioo t₂ T ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) 3,
      HasDerivAt (fun s => φ (s, q.2) - c₂)
        (potSym (jp 4 (fun y => φ (q.1, y) - c₂) q.2)) q.1 := by
    intro q hq
    have hqO := hS3 hq
    have := (hder q.1 hqO.1 q.2 hqO.2).sub_const c₂
    exact this.congr_deriv
      (lc_potSym_sub_const hJi hJJ hU isOpen_ball (hsl q.1 hqO.1) hqO.2 c₂ (le_refl 4)).symm
  have hlowS : ∀ q ∈ Ioo t₂ T ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) 3, ∀ m ≤ 3,
      ‖jet m (fun y => φ (q.1, y) - c₂) q.2‖ < 2 := fun q hq m hm =>
    hlowc q.1 (hS3 hq).1 q.2 (hB35 hq.2) c₂ (hsm q.1 hq.1 q.2 hq.2) m hm
  have hbd : ∀ m : ℕ, ∃ C : ℝ, ∀ q ∈ Ioo t₂ T ×ˢ ball (0 : EuclideanSpace ℝ (Fin n)) 3,
      ‖jet m (fun y => φ (q.1, y) - c₂) q.2‖ ≤ C := by
    intro m
    by_cases hm : 1 ≤ m
    · refine ⟨Cm m, fun q hq => ?_⟩
      rw [lc_jet_sub_const isOpen_ball (hsl q.1 (hS3 hq).1) (hB35 hq.2) c₂ hm]
      exact hCm' m hm q.1 (hS3 hq).1 q.2 hq.2
    · have : m = 0 := by omega
      subst this
      exact ⟨2, fun q hq => (hlowS q hq 0 (by norm_num)).le⟩
  obtain ⟨hGb, hVb⟩ := graph_iteratedFDeriv_bounded (F := fun q : ℝ × EuclideanSpace ℝ (Fin n) =>
    φ q - c₂) hF heq hlowS hbd
  exact ⟨t₂, ht₂, c₂, hlowS, hGb, hVb⟩

end Shift

end CWFlow.Extension
