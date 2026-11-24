import CWFlow.Extension.GraphAt
import CWFlow.Extension.FlowM
import CWFlow.Extension.GraphBridge
import CWFlow.Extension.GraphLocal
import CWFlow.Extension.EnergySymbol
import CWFlow.Extension.ImmSmooth
import CWFlow.Interp.BallInterp
import Mathlib.Analysis.Calculus.BumpFunction.Basic

/-!
# An adapted cover by graph charts

For a Lagrangian immersion `F : M → V` of a compact manifold with small second fundamental form, there are
finitely many graph charts `σ j` (with frames `U j` and potentials `φ j`, centered at `p j`) such that the
unit balls `σ j '' ball 0 1` cover `M`, and the number of charts is bounded by the volume of `F`.
-/

open Metric Set MeasureTheory
open CWFlow.Interp
open scoped ContDiff Manifold Topology RealInnerProductSpace

set_option linter.unusedSectionVars false

namespace CWFlow.Extension

variable {n : ℕ} {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
  {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
  {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
  (ι : EuclideanSpace ℝ (Fin n) ≃L[ℝ] E)
  {J : V →L[ℝ] V} {F : M → V}

/-! ### Part 1: `horC` and `verC` are 1-Lipschitz -/

theorem cv_bessel (w : Fin n → V) (hw : Orthonormal ℝ w) (v : V) :
    ‖(WithLp.toLp 2 (fun a => ⟪w a, v⟫) : EuclideanSpace ℝ (Fin n))‖ ≤ ‖v‖ := by
  rw [EuclideanSpace.norm_eq]
  refine Real.sqrt_le_iff.2 ⟨norm_nonneg _, ?_⟩
  simpa using hw.sum_inner_products_le (s := Finset.univ) v

theorem cv_orthonormal_U (U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V) :
    Orthonormal ℝ (fun a : Fin n => U (eb a)) := by
  rw [orthonormal_iff_ite]
  intro a b
  rw [gb_UU, gb_inner_eb_eb]

theorem cv_orthonormal_JU (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫)
    (U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V) :
    Orthonormal ℝ (fun a : Fin n => J (U (eb a))) := by
  rw [orthonormal_iff_ite]
  intro a b
  rw [gb_JUJU hJi, gb_inner_eb_eb]

theorem cv_horC_lip (U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V) (y₀ v v' : V) :
    ‖horC U y₀ v - horC U y₀ v'‖ ≤ ‖v - v'‖ := by
  have : horC U y₀ v - horC U y₀ v' =
      (WithLp.toLp 2 (fun a => ⟪U (eb a), v - v'⟫) : EuclideanSpace ℝ (Fin n)) := by
    ext a
    simp [horC, ← inner_sub_right]
  rw [this]
  exact cv_bessel _ (cv_orthonormal_U U) _

theorem cv_verC_lip (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫)
    (U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V) (y₀ v v' : V) :
    ‖verC J U y₀ v - verC J U y₀ v'‖ ≤ ‖v - v'‖ := by
  have : verC J U y₀ v - verC J U y₀ v' =
      (WithLp.toLp 2 (fun a => ⟪J (U (eb a)), v - v'⟫) : EuclideanSpace ℝ (Fin n)) := by
    ext a
    simp [verC, ← inner_sub_right]
  rw [this]
  exact cv_bessel _ (cv_orthonormal_JU hJi U) _

/-! ### Part 2: smoothness of `F` and of the cutoffs -/

theorem cv_contMDiff (hF : IsLagImm ι J F) : ContMDiff 𝓘(ℝ, E) 𝓘(ℝ, V) ∞ F := by
  intro p
  have h1 : ContDiffOn ℝ ∞ (F ∘ (extChartAt 𝓘(ℝ, E) p).symm) (extChartAt 𝓘(ℝ, E) p).target := by
    have := (hF.smooth p).comp ι.symm.contDiff.contDiffOn
      (s := (extChartAt 𝓘(ℝ, E) p).target) (fun e he => by simpa [chartDomE] using he)
    refine this.congr fun e he => ?_
    simp [chartInv]
  have h2 : ContMDiffOn 𝓘(ℝ, E) 𝓘(ℝ, V) ∞ (F ∘ (extChartAt 𝓘(ℝ, E) p).symm)
      (extChartAt 𝓘(ℝ, E) p).target := h1.contMDiffOn
  have hc : ContMDiffOn 𝓘(ℝ, E) 𝓘(ℝ, E) ∞ (extChartAt 𝓘(ℝ, E) p) (extChartAt 𝓘(ℝ, E) p).source := by
    rw [extChartAt_source]; exact contMDiffOn_extChartAt
  have h3 := h2.comp hc
    (fun x hx => (extChartAt 𝓘(ℝ, E) p).map_source hx)
  have h4 : ContMDiffOn 𝓘(ℝ, E) 𝓘(ℝ, V) ∞ F (extChartAt 𝓘(ℝ, E) p).source :=
    h3.congr fun x hx => by
      show F x = F ((extChartAt 𝓘(ℝ, E) p).symm (extChartAt 𝓘(ℝ, E) p x))
      rw [(extChartAt 𝓘(ℝ, E) p).left_inv hx]
  exact (h4.contMDiffAt ((isOpen_extChartAt_source p).mem_nhds (mem_extChartAt_source p)))

theorem cv_cutoff : ∃ η : EuclideanSpace ℝ (Fin n) → ℝ, ContDiff ℝ ∞ η ∧ (∀ x, 0 ≤ η x) ∧ (∀ x, η x ≤ 1) ∧
    (∀ x, ‖x‖ ≤ 1 / 16 → η x = 1) ∧ (∀ x, 1 / 8 ≤ ‖x‖ → η x = 0) := by
  let b : ContDiffBump (0 : EuclideanSpace ℝ (Fin n)) := ⟨1 / 16, 1 / 8, by norm_num, by norm_num⟩
  refine ⟨b, b.contDiff, b.nonneg', fun _ => b.le_one, fun x hx => ?_, fun x hx => ?_⟩
  · exact b.one_of_mem_closedBall (by simpa [b] using hx)
  · exact b.zero_of_le_dist (by simpa [b] using hx)

theorem cv_contDiff_horC (U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V) (y₀ : V) :
    ContDiff ℝ ∞ (horC U y₀) := by
  have : horC U y₀ = fun v => horCL U (v - y₀) := rfl
  rw [this]
  exact (horCL U).contDiff.comp (contDiff_id.sub contDiff_const)

open Classical in
/-- The cutoff on `M` attached to a graph chart. -/
noncomputable def cvChi (η : EuclideanSpace ℝ (Fin n) → ℝ) (U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V)
    (y₀ : V) (σ : EuclideanSpace ℝ (Fin n) → M) (F : M → V) : M → ℝ :=
  fun p => if p ∈ σ '' ball (0 : EuclideanSpace ℝ (Fin n)) 6 then η (horC U y₀ (F p)) else 0

theorem cv_chi_vanish {η : EuclideanSpace ℝ (Fin n) → ℝ} (hη0 : ∀ x, 1 / 8 ≤ ‖x‖ → η x = 0)
    {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} {y₀ : V} {σ : EuclideanSpace ℝ (Fin n) → M}
    (hhor : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6, horC U y₀ (F (σ x)) = x) :
    ∀ p, p ∉ σ '' ball (0 : EuclideanSpace ℝ (Fin n)) (1 / 8) → cvChi η U y₀ σ F p = 0 := by
  intro p hp
  unfold cvChi
  split_ifs with h
  · obtain ⟨x, hx, rfl⟩ := h
    rw [hhor x hx]
    refine hη0 x ?_
    by_contra hlt
    exact hp ⟨x, by simpa using not_le.1 hlt, rfl⟩
  · rfl

theorem cv_chi_tsupport [T2Space M] {η : EuclideanSpace ℝ (Fin n) → ℝ} (hη0 : ∀ x, 1 / 8 ≤ ‖x‖ → η x = 0)
    {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} {y₀ : V} {σ : EuclideanSpace ℝ (Fin n) → M}
    (hσ : IsParam ι (ball 0 6) σ)
    (hhor : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6, horC U y₀ (F (σ x)) = x) :
    tsupport (cvChi η U y₀ σ F) ⊆ σ '' closedBall (0 : EuclideanSpace ℝ (Fin n)) (1 / 8) := by
  have hK : IsCompact (σ '' closedBall (0 : EuclideanSpace ℝ (Fin n)) (1 / 8)) :=
    (isCompact_closedBall _ _).image_of_continuousOn
      (hσ.cont.mono (closedBall_subset_ball (by norm_num)))
  refine closure_minimal (fun p hp => ?_) hK.isClosed
  by_contra hn
  exact hp (cv_chi_vanish hη0 hhor p (fun h => hn (image_mono ball_subset_closedBall h)))

theorem cv_chi_contMDiff [T2Space M] {η : EuclideanSpace ℝ (Fin n) → ℝ} (hη : ContDiff ℝ ∞ η)
    (hη0 : ∀ x, 1 / 8 ≤ ‖x‖ → η x = 0) (hFs : ContMDiff 𝓘(ℝ, E) 𝓘(ℝ, V) ∞ F)
    {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} {y₀ : V} {σ : EuclideanSpace ℝ (Fin n) → M}
    (hσ : IsParam ι (ball 0 6) σ)
    (hhor : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6, horC U y₀ (F (σ x)) = x) :
    ContMDiff 𝓘(ℝ, E) 𝓘(ℝ) ∞ (cvChi η U y₀ σ F) := by
  have hg : ContMDiff 𝓘(ℝ, E) 𝓘(ℝ) ∞ (fun p => η (horC U y₀ (F p))) :=
    (hη.comp (cv_contDiff_horC U y₀)).contMDiff.comp hFs
  have hK : IsCompact (σ '' closedBall (0 : EuclideanSpace ℝ (Fin n)) (1 / 8)) :=
    (isCompact_closedBall _ _).image_of_continuousOn
      (hσ.cont.mono (closedBall_subset_ball (by norm_num)))
  have hvan := cv_chi_vanish (η := η) hη0 (U := U) (y₀ := y₀) (F := F) hhor
  intro p
  by_cases hp : p ∈ σ '' ball (0 : EuclideanSpace ℝ (Fin n)) 6
  · refine (hg p).congr_of_eventuallyEq ?_
    filter_upwards [hσ.isOpen_image.mem_nhds hp] with q hq
    simp [cvChi, hq]
  · have : cvChi η U y₀ σ F =ᶠ[𝓝 p] fun _ => 0 := by
      filter_upwards [hK.isClosed.isOpen_compl.mem_nhds (fun hm => hp
        (image_mono (closedBall_subset_ball (by norm_num : (1 / 8 : ℝ) < 6)) hm))]
        with q hq using hvan q (fun h => hq (image_mono ball_subset_closedBall h))
    exact (contMDiffAt_const).congr_of_eventuallyEq this

/-! ### Part 3: the volume lower bound -/

theorem cv_one_le_volX (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V}
    (hU : IsLagFrame J U) {y₀ : V} {σ : EuclideanSpace ℝ (Fin n) → M}
    {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiffOn ℝ ∞ φ (ball 0 6))
    (hgr : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6, F (σ x) = graphX J U y₀ φ x)
    {x : EuclideanSpace ℝ (Fin n)} (hx : x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6) :
    1 ≤ volX (F ∘ σ) x := by
  have hxn : ‖x‖ < 6 := by simpa using hx
  set r : ℝ := (6 - ‖x‖) / 3 with hr
  have hr0 : 0 < r := by rw [hr]; linarith
  let b : ContDiffBump x := ⟨r, 2 * r, hr0, by linarith⟩
  have hsub : closedBall x (2 * r) ⊆ ball (0 : EuclideanSpace ℝ (Fin n)) 6 := by
    intro z hz
    rw [mem_closedBall, dist_eq_norm] at hz
    rw [mem_ball, dist_zero_right]
    calc ‖z‖ = ‖(z - x) + x‖ := by rw [sub_add_cancel]
      _ ≤ ‖z - x‖ + ‖x‖ := norm_add_le _ _
      _ < 6 := by linarith
  have hψ : ContDiff ℝ ∞ (fun z => b z * φ z) := by
    rw [contDiff_iff_contDiffAt]
    intro z
    by_cases hz : z ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6
    · exact b.contDiffAt.mul (hφ.contDiffAt (isOpen_ball.mem_nhds hz))
    · have hz' : z ∈ (closedBall x (2 * r))ᶜ := fun h => hz (hsub h)
      have hev : (fun z => b z * φ z) =ᶠ[𝓝 z] fun _ => 0 := by
        filter_upwards [isClosed_closedBall.isOpen_compl.mem_nhds hz'] with w hw
        have : b w = 0 := b.zero_of_le_dist (by
          rw [mem_compl_iff, mem_closedBall, not_le] at hw; exact hw.le)
        simp [this]
      exact contDiffAt_const.congr_of_eventuallyEq hev
  have hψφ : ∀ w ∈ ball x r, b w * φ w = φ w := by
    intro w hw
    have : b w = 1 := b.one_of_mem_closedBall (mem_closedBall.2 (mem_ball.1 hw).le)
    simp [this]
  have hev : ∀ z ∈ ball x r, (fun w => b w * φ w) =ᶠ[𝓝 z] φ := by
    intro z hz
    filter_upwards [isOpen_ball.mem_nhds hz] with w hw using hψφ w hw
  have hgrad : ∀ z ∈ ball x r, gradE (fun w => b w * φ w) z = gradE φ z := by
    intro z hz
    unfold gradE
    rw [(hev z hz).fderiv_eq]
  have hFe : F ∘ σ =ᶠ[𝓝 x] graphX J U y₀ (fun w => b w * φ w) := by
    filter_upwards [isOpen_ball.mem_nhds (mem_ball_self hr0)] with z hz
    have hz6 : z ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6 :=
      hsub (ball_subset_closedBall.trans (closedBall_subset_closedBall (by linarith)) hz)
    simp only [Function.comp]
    rw [hgr z hz6]
    simp only [graphX, hgrad z hz]
  calc (1 : ℝ) ≤ volFun (jet 2 (fun w => b w * φ w) x) := one_le_volFun _
    _ = volX (graphX J U y₀ (fun w => b w * φ w)) x := (volX_graphX hJi hU hψ x).symm
    _ = volX (F ∘ σ) x := (volX_congr hFe).symm

theorem cv_vol_lower [T2Space M] [CompactSpace M] (A : IntAtlas (E := E) M)
    (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hF : IsLagImm ι J F)
    {η : EuclideanSpace ℝ (Fin n) → ℝ} (hη : ContDiff ℝ ∞ η) (hη0 : ∀ x, 0 ≤ η x)
    (hη1 : ∀ x, ‖x‖ ≤ 1 / 16 → η x = 1) (hη2 : ∀ x, 1 / 8 ≤ ‖x‖ → η x = 0)
    {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} (hU : IsLagFrame J U) {y₀ : V}
    {σ : EuclideanSpace ℝ (Fin n) → M} (hσ : IsParam ι (ball 0 6) σ)
    {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiffOn ℝ ∞ φ (ball 0 6))
    (hgr : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6,
      F (σ x) = graphX J U y₀ φ x ∧ horC U y₀ (F (σ x)) = x) :
    (volume (ball (0 : EuclideanSpace ℝ (Fin n)) (1 / 16))).toReal ≤
      intM ι A F (cvChi η U y₀ σ F) (densOf volX) := by
  have hhor := fun x hx => (hgr x hx).2
  have hsm := cv_chi_contMDiff ι hη hη2 (cv_contMDiff ι hF) hσ hhor
  have hsupp : tsupport (cvChi η U y₀ σ F) ⊆ σ '' ball 0 6 :=
    (cv_chi_tsupport ι hη2 hσ hhor).trans
      (image_mono (closedBall_subset_ball (by norm_num)))
  rw [intM_eq_param ι A hF hsm (isNatDens_densOf_volX J) hσ hsupp]
  have hint := im_integrableOn_param ι hF hsm (isNatDens_densOf_volX J) hσ hsupp
  have hfun : ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6,
      densOf volX (F ∘ σ) (cvChi η U y₀ σ F ∘ σ) y = η y * volX (F ∘ σ) y := by
    intro y hy
    have : cvChi η U y₀ σ F (σ y) = η y := by
      simp only [cvChi, if_pos (mem_image_of_mem σ hy), hhor y hy]
    simp only [densOf, Function.comp_apply, this]
  have hge : ∀ y ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6, 1 ≤ volX (F ∘ σ) y := fun y hy =>
    cv_one_le_volX hJi hU hφ (fun x hx => (hgr x hx).1) hy
  have hB : ball (0 : EuclideanSpace ℝ (Fin n)) (1 / 16) ⊆ ball 0 6 := ball_subset_ball (by norm_num)
  calc (volume (ball (0 : EuclideanSpace ℝ (Fin n)) (1 / 16))).toReal
      = ∫ y in ball (0 : EuclideanSpace ℝ (Fin n)) (1 / 16), (1 : ℝ) := by
        simp [Measure.real]
    _ ≤ ∫ y in ball (0 : EuclideanSpace ℝ (Fin n)) (1 / 16),
        densOf volX (F ∘ σ) (cvChi η U y₀ σ F ∘ σ) y := by
        refine setIntegral_mono_on (integrableOn_const measure_ball_lt_top.ne) (hint.mono_set hB)
          measurableSet_ball fun y hy => ?_
        rw [hfun y (hB hy), hη1 y (by simpa using (mem_ball_zero_iff.1 hy).le), one_mul]
        exact hge y (hB hy)
    _ ≤ ∫ y in ball (0 : EuclideanSpace ℝ (Fin n)) 6,
        densOf volX (F ∘ σ) (cvChi η U y₀ σ F ∘ σ) y := by
        refine setIntegral_mono_set hint ?_ hB.eventuallyLE
        refine (ae_restrict_iff' measurableSet_ball).2 (ae_of_all _ fun y hy => ?_)
        rw [hfun y hy]
        exact mul_nonneg (hη0 y) (zero_le_one.trans (hge y hy))

theorem cv_chi_bounds {η : EuclideanSpace ℝ (Fin n) → ℝ} (hη0 : ∀ x, 0 ≤ η x) (hη1 : ∀ x, η x ≤ 1)
    (U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V) (y₀ : V) (σ : EuclideanSpace ℝ (Fin n) → M) (p : M) :
    0 ≤ cvChi η U y₀ σ F p ∧ cvChi η U y₀ σ F p ≤ 1 := by
  unfold cvChi
  split_ifs
  · exact ⟨hη0 _, hη1 _⟩
  · exact ⟨le_rfl, zero_le_one⟩

/-! ### Part 4: counting -/

theorem cv_contMDiff_sum {α : Type*} (s : Finset α) (h : α → M → ℝ)
    (hh : ∀ a ∈ s, ContMDiff 𝓘(ℝ, E) 𝓘(ℝ) ∞ (h a)) :
    ContMDiff 𝓘(ℝ, E) 𝓘(ℝ) ∞ (fun p => ∑ a ∈ s, h a p) := by
  classical
  induction s using Finset.induction_on with
  | empty => simpa using contMDiff_const
  | insert a s ha ih =>
    have h' : ContMDiff 𝓘(ℝ, E) 𝓘(ℝ) ∞ (fun p => h a p + ∑ b ∈ s, h b p) :=
      (hh a (Finset.mem_insert_self a s)).add (ih fun b hb => hh b (Finset.mem_insert_of_mem hb))
    simpa [Finset.sum_insert ha] using h'

theorem cv_volX_nonneg (X : EuclideanSpace ℝ (Fin n) → V) (y : EuclideanSpace ℝ (Fin n)) :
    0 ≤ volX X y := Real.sqrt_nonneg _

theorem cv_intM_mono [T2Space M] [CompactSpace M] (A : IntAtlas (E := E) M) (hF : IsLagImm ι J F)
    {h₁ h₂ : M → ℝ} (hs₁ : ContMDiff 𝓘(ℝ, E) 𝓘(ℝ) ∞ h₁) (hs₂ : ContMDiff 𝓘(ℝ, E) 𝓘(ℝ) ∞ h₂)
    (hle : ∀ p, h₁ p ≤ h₂ p) :
    intM ι A F h₁ (densOf volX) ≤ intM ι A F h₂ (densOf volX) := by
  have hd := isNatDens_densOf_volX (n := n) J
  have hsum := intM_densOf_sum ι A hF hd (Finset.univ : Finset (Fin 2)) ![h₁, fun p => h₂ p - h₁ p]
    (by
      intro a _
      fin_cases a
      · exact hs₁
      · exact hs₂.sub hs₁)
  have heq : (fun p => ∑ a ∈ (Finset.univ : Finset (Fin 2)), (![h₁, fun p => h₂ p - h₁ p] : Fin 2 → M → ℝ) a p)
      = h₂ := by
    funext p
    simp [Fin.sum_univ_two]
  rw [heq, Fin.sum_univ_two] at hsum
  have hnn : 0 ≤ intM ι A F (fun p => h₂ p - h₁ p) (densOf volX) :=
    intM_densOf_nonneg ι A (fun p => sub_nonneg.2 (hle p)) cv_volX_nonneg
  simp only [Matrix.cons_val_zero, Matrix.cons_val_one] at hsum
  linarith

theorem cv_count [T2Space M] [CompactSpace M] (A : IntAtlas (E := E) M) (hF : IsLagImm ι J F)
    (χ : M → M → ℝ) (hχ : ∀ q, ContMDiff 𝓘(ℝ, E) 𝓘(ℝ) ∞ (χ q)) {v : ℝ}
    (hlow : ∀ q, v ≤ intM ι A F (χ q) (densOf volX)) (S : Finset M)
    (hS : ∀ p, ∑ q ∈ S, χ q p ≤ 1) : (S.card : ℝ) * v ≤ volM ι A F := by
  have hd := isNatDens_densOf_volX (n := n) J
  have h1 : intM ι A F (fun p => ∑ q ∈ S, χ q p) (densOf volX) =
      ∑ q ∈ S, intM ι A F (χ q) (densOf volX) :=
    intM_densOf_sum ι A hF hd S χ (fun q _ => hχ q)
  have h2 : intM ι A F (fun p => ∑ q ∈ S, χ q p) (densOf volX) ≤ volM ι A F :=
    cv_intM_mono ι A hF (cv_contMDiff_sum S χ fun q _ => hχ q) contMDiff_const hS
  calc (S.card : ℝ) * v = ∑ _q ∈ S, v := by simp
    _ ≤ ∑ q ∈ S, intM ι A F (χ q) (densOf volX) := Finset.sum_le_sum fun q _ => hlow q
    _ = _ := h1.symm
    _ ≤ volM ι A F := h2

/-! ### Part 5: nested charts -/

theorem cv_nested (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) {c₀ : ℝ} (hc₀ : 0 < c₀) (hc₀1 : 6 * c₀ < 1)
    {p : M} {Up : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} {σq σp : EuclideanSpace ℝ (Fin n) → M}
    (hqc : ContinuousOn σq (ball 0 6))
    (hqL : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6, ∀ x' ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6,
      ‖F (σq x) - F (σq x')‖ ≤ (1 + c₀) * ‖x - x'‖)
    (hph : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6, horC Up (F p) (F (σp x)) = x)
    (hpv : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6, ‖verC J Up (F p) (F (σp x))‖ ≤ c₀ * ‖x‖)
    (hpc : σp '' ball 0 1 = connectedComponentIn
      {r | ‖horC Up (F p) (F r)‖ < 1 ∧ ‖verC J Up (F p) (F r)‖ < 1} p)
    {m : M} (hm₁ : m ∈ σq '' ball (0 : EuclideanSpace ℝ (Fin n)) (1 / 8))
    (hm₂ : m ∈ σp '' ball (0 : EuclideanSpace ℝ (Fin n)) (1 / 8)) :
    σq '' ball (0 : EuclideanSpace ℝ (Fin n)) (1 / 8) ⊆ σp '' ball 0 1 := by
  obtain ⟨x₁, hx₁, rfl⟩ := hm₁
  obtain ⟨x₀, hx₀, hx₀m⟩ := hm₂
  have h6 : ball (0 : EuclideanSpace ℝ (Fin n)) (1 / 8) ⊆ ball 0 6 := ball_subset_ball (by norm_num)
  have hx₀n : ‖x₀‖ < 1 / 8 := by simpa using hx₀
  have hx₁n : ‖x₁‖ < 1 / 8 := by simpa using hx₁
  have hhm : horC Up (F p) (F (σq x₁)) = x₀ := by rw [← hx₀m]; exact hph x₀ (h6 hx₀)
  have hvm : ‖verC J Up (F p) (F (σq x₁))‖ ≤ c₀ * ‖x₀‖ := by
    rw [← hx₀m]; exact hpv x₀ (h6 hx₀)
  have hsub : σq '' ball (0 : EuclideanSpace ℝ (Fin n)) (1 / 8) ⊆
      {r | ‖horC Up (F p) (F r)‖ < 1 ∧ ‖verC J Up (F p) (F r)‖ < 1} := by
    rintro _ ⟨x, hx, rfl⟩
    have hxn : ‖x‖ < 1 / 8 := by simpa using hx
    have hd : ‖F (σq x) - F (σq x₁)‖ < (1 + c₀) * (1 / 4) := by
      calc ‖F (σq x) - F (σq x₁)‖ ≤ (1 + c₀) * ‖x - x₁‖ := hqL x (h6 hx) x₁ (h6 hx₁)
        _ < (1 + c₀) * (1 / 4) := by
          refine mul_lt_mul_of_pos_left ?_ (by linarith)
          calc ‖x - x₁‖ ≤ ‖x‖ + ‖x₁‖ := norm_sub_le _ _
            _ < 1 / 4 := by linarith
    refine ⟨?_, ?_⟩
    · have h1 := cv_horC_lip Up (F p) (F (σq x)) (F (σq x₁))
      have h2 := norm_le_insert' (horC Up (F p) (F (σq x))) (horC Up (F p) (F (σq x₁)))
      have h3 : ‖horC Up (F p) (F (σq x₁))‖ = ‖x₀‖ := by rw [hhm]
      linarith
    · have h1 := cv_verC_lip hJi Up (F p) (F (σq x)) (F (σq x₁))
      have h2 := norm_le_insert' (verC J Up (F p) (F (σq x))) (verC J Up (F p) (F (σq x₁)))
      nlinarith
  have hpre : IsPreconnected (σq '' ball (0 : EuclideanSpace ℝ (Fin n)) (1 / 8)) :=
    ((convex_ball _ _).isPreconnected).image _ (hqc.mono h6)
  have hmS : σq x₁ ∈ connectedComponentIn
      {r | ‖horC Up (F p) (F r)‖ < 1 ∧ ‖verC J Up (F p) (F r)‖ < 1} p := by
    rw [← hpc]
    exact ⟨x₀, ball_subset_ball (by norm_num) hx₀, hx₀m⟩
  have := hpre.subset_connectedComponentIn (mem_image_of_mem σq hx₁) hsub
  rw [hpc, connectedComponentIn_eq hmS]
  exact this

/-! ### Part 6: the cover -/

theorem exists_cover [T2Space M] [CompactSpace M] (A : IntAtlas (E := E) M)
    (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    [FiniteDimensional ℝ V] (hdim : Module.finrank ℝ V = 2 * n) {F : M → V} (hF : IsLagImm ι J F)
    {c₀ c₁ : ℝ} (hc₀ : 0 < c₀) (hc₀1 : 6 * c₀ < 1) (hc₁ : 0 ≤ c₁)
    (hsmall : 6 * (c₁ * (1 + c₀ ^ 2) ^ ((3 : ℝ) / 2)) < c₀)
    (hA : ∀ z : M, ∀ y ∈ chartDomE ι z, aNormSq J (F ∘ chartInv ι z) y ≤ c₁ ^ 2) :
    ∃ (N : ℕ) (p : Fin N → M) (U : Fin N → (EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V))
      (σ : Fin N → EuclideanSpace ℝ (Fin n) → M) (φ : Fin N → EuclideanSpace ℝ (Fin n) → ℝ),
      (N : ℝ) * (volume (ball (0 : EuclideanSpace ℝ (Fin n)) (1 / 16))).toReal ≤ volM ι A F ∧
      (∀ j, IsLagFrame J (U j) ∧ IsParam ι (ball 0 6) (σ j) ∧ σ j 0 = p j ∧
        ContDiffOn ℝ ∞ (φ j) (ball 0 6) ∧ φ j 0 = 0 ∧
        (∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6,
          F (σ j x) = graphX J (U j) (F (p j)) (φ j) x ∧ horC (U j) (F (p j)) (F (σ j x)) = x ∧
          ‖jet 2 (φ j) x‖ ≤ c₀ ∧ ‖gradE (φ j) x‖ ≤ c₀ * ‖x‖) ∧
        ∀ ρ, 0 < ρ → ρ ≤ 6 → σ j '' ball 0 ρ = connectedComponentIn
          {q | ‖horC (U j) (F (p j)) (F q)‖ < ρ ∧ ‖verC J (U j) (F (p j)) (F q)‖ < 1} (p j)) ∧
      ∀ q : M, ∃ j, q ∈ σ j '' ball 0 1 := by
  classical
  have hg := fun q => graph_at ι hJi hJJ hdim hF hc₀ hc₀1 hc₁ hsmall hA q
  choose U σ φ hU hP hσ0 hφ hφ0 hb hLip hcomp using hg
  obtain ⟨η, hη, hη0, hηle, hη1, hη2⟩ := cv_cutoff (n := n)
  set v : ℝ := (volume (ball (0 : EuclideanSpace ℝ (Fin n)) (1 / 16))).toReal with hv
  have hvpos : 0 < v := ENNReal.toReal_pos
    (Metric.measure_ball_pos volume _ (by norm_num)).ne' measure_ball_lt_top.ne
  have hhor : ∀ q, ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6,
      horC (U q) (F q) (F (σ q x)) = x := fun q x hx => (hb q x hx).2.1
  let χ : M → M → ℝ := fun q => cvChi η (U q) (F q) (σ q) F
  have hχ : ∀ q, ContMDiff 𝓘(ℝ, E) 𝓘(ℝ) ∞ (χ q) := fun q =>
    cv_chi_contMDiff ι hη hη2 (cv_contMDiff ι hF) (hP q) (hhor q)
  have hlow : ∀ q, v ≤ intM ι A F (χ q) (densOf volX) := fun q =>
    cv_vol_lower ι A hJi hF hη hη0 hη1 hη2 (hU q) (hP q) (hφ q)
      (fun x hx => ⟨(hb q x hx).1, (hb q x hx).2.1⟩)
  let adm : Finset M → Prop := fun S => ∀ a ∈ S, ∀ b ∈ S, a ≠ b →
    Disjoint (σ a '' ball (0 : EuclideanSpace ℝ (Fin n)) (1 / 8)) (σ b '' ball 0 (1 / 8))
  have hsum : ∀ S, adm S → ∀ p, ∑ q ∈ S, χ q p ≤ 1 := by
    intro S hS p
    by_cases h : ∃ q₀ ∈ S, p ∈ σ q₀ '' ball (0 : EuclideanSpace ℝ (Fin n)) (1 / 8)
    · obtain ⟨q₀, hq₀, hp⟩ := h
      rw [Finset.sum_eq_single_of_mem q₀ hq₀ (fun q hq hne =>
        cv_chi_vanish (F := F) hη2 (hhor q) p
          (fun hpq => Set.disjoint_left.1 (hS q hq q₀ hq₀ hne) hpq hp))]
      exact (cv_chi_bounds (F := F) hη0 hηle (U q₀) (F q₀) (σ q₀) p).2
    · rw [Finset.sum_eq_zero (fun q hq => cv_chi_vanish (F := F) hη2 (hhor q) p
        (fun hpq => h ⟨q, hq, hpq⟩))]
      exact zero_le_one
  have hcard : ∀ S, adm S → (S.card : ℝ) * v ≤ volM ι A F := fun S hS =>
    cv_count ι A hF χ hχ hlow S (hsum S hS)
  have hK : ∀ S, adm S → S.card ≤ ⌊volM ι A F / v⌋₊ := fun S hS =>
    Nat.le_floor ((le_div_iff₀ hvpos).2 (hcard S hS))
  let T : Set ℕ := {k | ∃ S : Finset M, adm S ∧ S.card = k}
  have hne : T.Nonempty := ⟨0, ∅, by intro a ha; simp at ha, rfl⟩
  have hbdd : BddAbove T := ⟨⌊volM ι A F / v⌋₊, fun k ⟨S, hS, hk⟩ => hk ▸ hK S hS⟩
  obtain ⟨S, hS, hSc⟩ := Nat.sSup_mem hne hbdd
  have hq0 : ∀ q, q ∈ σ q '' ball (0 : EuclideanSpace ℝ (Fin n)) (1 / 8) := fun q =>
    ⟨0, mem_ball_self (by norm_num), hσ0 q⟩
  have hmax : ∀ q, ∃ p ∈ S, (σ q '' ball (0 : EuclideanSpace ℝ (Fin n)) (1 / 8) ∩
      σ p '' ball (0 : EuclideanSpace ℝ (Fin n)) (1 / 8)).Nonempty := by
    intro q
    by_contra hcon
    have hdis : ∀ p ∈ S, Disjoint (σ q '' ball (0 : EuclideanSpace ℝ (Fin n)) (1 / 8))
        (σ p '' ball (0 : EuclideanSpace ℝ (Fin n)) (1 / 8)) := fun p hp =>
      Set.disjoint_left.2 fun m hm1 hm2 => hcon ⟨p, hp, ⟨m, hm1, hm2⟩⟩
    have hqS : q ∉ S := fun hq => Set.disjoint_left.1 (hdis q hq) (hq0 q) (hq0 q)
    have hadm : adm (insert q S) := by
      intro a ha b hb hab
      rcases Finset.mem_insert.1 ha with rfl | ha'
      · rcases Finset.mem_insert.1 hb with rfl | hb'
        · exact absurd rfl hab
        · exact hdis b hb'
      · rcases Finset.mem_insert.1 hb with rfl | hb'
        · exact (hdis a ha').symm
        · exact hS a ha' b hb' hab
    have hle : (insert q S).card ≤ sSup T := le_csSup hbdd ⟨insert q S, hadm, rfl⟩
    rw [Finset.card_insert_of_notMem hqS, ← hSc] at hle
    omega
  have hcov : ∀ q, ∃ p ∈ S, q ∈ σ p '' ball (0 : EuclideanSpace ℝ (Fin n)) 1 := by
    intro q
    obtain ⟨p, hp, m, hm1, hm2⟩ := hmax q
    refine ⟨p, hp, ?_⟩
    have := cv_nested hJi hc₀ hc₀1 (p := p) (hP q).cont
      (fun x hx x' hx' => hLip q x hx x' hx') (fun x hx => (hb p x hx).2.1)
      (fun x hx => (hb p x hx).2.2.2.2) (hcomp p 1 one_pos (by norm_num)) hm1 hm2
    exact this (hq0 q)
  refine ⟨S.card, fun j => ((S.equivFin.symm j : S) : M),
    fun j => U ((S.equivFin.symm j : S) : M), fun j => σ ((S.equivFin.symm j : S) : M),
    fun j => φ ((S.equivFin.symm j : S) : M), hcard S hS, fun j => ?_, fun q => ?_⟩
  · refine ⟨hU _, hP _, hσ0 _, hφ _, hφ0 _, fun x hx => ?_, hcomp _⟩
    obtain ⟨h1, h2, h3, h4, -⟩ := hb _ x hx
    exact ⟨h1, h2, h3, h4⟩
  · obtain ⟨p, hp, hq⟩ := hcov q
    refine ⟨S.equivFin ⟨p, hp⟩, ?_⟩
    simpa using hq

/-- info: 'CWFlow.Extension.exists_cover' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_cover

end CWFlow.Extension
