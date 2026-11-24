import Mathlib.Analysis.Calculus.FDeriv.Extend
import Mathlib.Analysis.Calculus.ContDiff.Defs
import Mathlib.Analysis.Calculus.ContDiff.FTaylorSeries
import Mathlib.Analysis.Calculus.MeanValue
import Mathlib.Analysis.Convex.Topology
import Mathlib.Tactic

/-!
# Smooth extension to the closure of a convex open set

A smooth function on an open convex set `U` of a normed space all of whose iterated derivatives
are bounded extends to a smooth function on `closure U`.  The extension is built derivative by
derivative: each iterated derivative is Lipschitz on `U` (mean value inequality), hence extends
continuously to the closure, and the extended derivatives form a Taylor series on `closure U`
(`hasFDerivWithinAt_closure_of_tendsto_fderiv`).
-/

open scoped ContDiff Topology
open Filter Set

namespace CWFlow.Interp

/-- A Lipschitz map on a subset `U` into a complete space extends continuously to `closure U`. -/
theorem exists_continuousOn_closure_of_lipschitz {E G : Type*} [NormedAddCommGroup E]
    [NormedAddCommGroup G] [CompleteSpace G] {U : Set E} {g : E → G} {K : ℝ}
    (hK : ∀ x ∈ U, ∀ y ∈ U, ‖g x - g y‖ ≤ K * ‖x - y‖) :
    ∃ h : E → G, ContinuousOn h (closure U) ∧ EqOn h g U := by
  classical
  have hlim : ∀ x ∈ closure U, ∃ c, Tendsto g (𝓝[U] x) (𝓝 c) := by
    intro x hx
    have : (𝓝[U] x).NeBot := mem_closure_iff_nhdsWithin_neBot.mp hx
    rw [← cauchy_map_iff_exists_tendsto, Metric.cauchy_iff]
    refine ⟨inferInstance, fun ε hε => ?_⟩
    set δ : ℝ := ε / (4 * (|K| + 1)) with hδ
    have hKp : 0 < |K| + 1 := by positivity
    have hδ0 : 0 < δ := by positivity
    refine ⟨g '' (U ∩ Metric.ball x δ), image_mem_map ?_, ?_⟩
    · exact inter_mem_nhdsWithin U (Metric.ball_mem_nhds x hδ0)
    · rintro _ ⟨a, ⟨haU, hax⟩, rfl⟩ _ ⟨b, ⟨hbU, hbx⟩, rfl⟩
      rw [dist_eq_norm]
      refine lt_of_le_of_lt (hK a haU b hbU) ?_
      have h1 : ‖a - b‖ < 2 * δ := by
        have := dist_triangle_right a b x
        rw [dist_eq_norm] at this
        rw [Metric.mem_ball] at hax hbx
        linarith
      calc K * ‖a - b‖ ≤ |K| * ‖a - b‖ := mul_le_mul_of_nonneg_right (le_abs_self K) (norm_nonneg _)
        _ ≤ (|K| + 1) * (2 * δ) := by
          apply mul_le_mul (by linarith) h1.le (norm_nonneg _) hKp.le
        _ < ε := by
          rw [hδ]
          field_simp
          linarith
  choose! c hc using hlim
  have hcg : EqOn c g U := by
    intro x hx
    have : (𝓝[U] x).NeBot := mem_closure_iff_nhdsWithin_neBot.mp (subset_closure hx)
    have hg : Tendsto g (𝓝[U] x) (𝓝 (g x)) := by
      rw [tendsto_iff_norm_sub_tendsto_zero]
      have hs : Tendsto (fun y => ‖K‖ * ‖y - x‖) (𝓝[U] x) (𝓝 0) := by
        have : Tendsto (fun y : E => ‖K‖ * ‖y - x‖) (𝓝 x) (𝓝 (‖K‖ * ‖x - x‖)) :=
          ((continuous_const.mul (continuous_norm.comp (continuous_id.sub continuous_const)))
            |>.tendsto x)
        simpa using this.mono_left nhdsWithin_le_nhds
      refine squeeze_zero' (Eventually.of_forall fun _ => norm_nonneg _) ?_ hs
      filter_upwards [self_mem_nhdsWithin] with y hy
      have := hK y hy x hx
      calc ‖g y - g x‖ ≤ K * ‖y - x‖ := this
        _ ≤ ‖K‖ * ‖y - x‖ := mul_le_mul_of_nonneg_right (le_abs_self K) (norm_nonneg _)
    exact tendsto_nhds_unique (hc x (subset_closure hx)) hg
  refine ⟨c, ?_, hcg⟩
  have hlip : ∀ x ∈ closure U, ∀ y ∈ closure U, dist (c x) (c y) ≤ (‖K‖₊ : ℝ) * dist x y := by
    intro x hx y hy
    have : (𝓝[U] x).NeBot := mem_closure_iff_nhdsWithin_neBot.mp hx
    have : (𝓝[U] y).NeBot := mem_closure_iff_nhdsWithin_neBot.mp hy
    set l : Filter (E × E) := 𝓝[U] x ×ˢ 𝓝[U] y with hl
    have : l.NeBot := Filter.prod_neBot.mpr ⟨inferInstance, inferInstance⟩
    have h1 : Tendsto (fun p : E × E => p.1) l (𝓝 x) :=
      (tendsto_fst : Tendsto Prod.fst l (𝓝[U] x)).mono_right nhdsWithin_le_nhds
    have h2 : Tendsto (fun p : E × E => p.2) l (𝓝 y) :=
      (tendsto_snd : Tendsto Prod.snd l (𝓝[U] y)).mono_right nhdsWithin_le_nhds
    have t1 : Tendsto (fun p : E × E => ‖g p.1 - g p.2‖) l (𝓝 ‖c x - c y‖) :=
      (((hc x hx).comp tendsto_fst).sub ((hc y hy).comp tendsto_snd)).norm
    have t2 : Tendsto (fun p : E × E => ‖K‖ * ‖p.1 - p.2‖) l (𝓝 (‖K‖ * ‖x - y‖)) :=
      tendsto_const_nhds.mul (h1.sub h2).norm
    have ev : ∀ᶠ p : E × E in l, ‖g p.1 - g p.2‖ ≤ ‖K‖ * ‖p.1 - p.2‖ := by
      have e1 : ∀ᶠ p : E × E in l, p.1 ∈ U := tendsto_fst.eventually self_mem_nhdsWithin
      have e2 : ∀ᶠ p : E × E in l, p.2 ∈ U := tendsto_snd.eventually self_mem_nhdsWithin
      filter_upwards [e1, e2] with p hp1 hp2
      exact (hK _ hp1 _ hp2).trans
        (mul_le_mul_of_nonneg_right (le_abs_self K) (norm_nonneg _))
    have := le_of_tendsto_of_tendsto t1 t2 ev
    simpa [dist_eq_norm] using this
  exact (LipschitzOnWith.of_dist_le_mul hlip).continuousOn

/-- A smooth function on an open convex set with all iterated derivatives bounded extends to a
smooth function on the closure. -/
theorem exists_contDiffOn_closure_of_bounded {E F : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    [NormedAddCommGroup F] [NormedSpace ℝ F] [CompleteSpace F]
    {U : Set E} (hU : IsOpen U) (hconv : Convex ℝ U) {f : E → F} (hf : ContDiffOn ℝ ∞ f U)
    (hb : ∀ k : ℕ, ∃ C : ℝ, ∀ x ∈ U, ‖iteratedFDeriv ℝ k f x‖ ≤ C) :
    ∃ g : E → F, ContDiffOn ℝ ∞ g (closure U) ∧ Set.EqOn g f U := by
  classical
  choose C hC using hb
  have hdiff : ∀ k : ℕ, ∀ x ∈ U, DifferentiableAt ℝ (iteratedFDeriv ℝ k f) x := fun k x hx =>
    (hf.contDiffAt (hU.mem_nhds hx)).differentiableAt_iteratedFDeriv
      (WithTop.coe_lt_coe.2 (ENat.natCast_lt_top k))
  have hlip : ∀ k : ℕ, ∀ x ∈ U, ∀ y ∈ U,
      ‖iteratedFDeriv ℝ k f x - iteratedFDeriv ℝ k f y‖ ≤ C (k + 1) * ‖x - y‖ := by
    intro k x hx y hy
    exact hconv.norm_image_sub_le_of_norm_fderiv_le (hdiff k)
      (fun z hz => by rw [norm_fderiv_iteratedFDeriv]; exact hC _ z hz) hy hx
  choose P hPc hPe using fun k : ℕ => exists_continuousOn_closure_of_lipschitz (hlip k)
  have hder : ∀ (k : ℕ), ∀ x ∈ closure U,
      HasFDerivWithinAt (P k) ((P (k + 1) x).curryLeft) (closure U) x := by
    intro k x hx
    have hx' : ContinuousWithinAt (P (k + 1)) U x := (hPc (k + 1) x hx).mono subset_closure
    refine hasFDerivWithinAt_closure_of_tendsto_fderiv ?_ hconv hU
      (fun y hy => (hPc k y hy).mono subset_closure) ?_
    · intro y hy
      exact ((hdiff k y hy).differentiableWithinAt).congr (hPe k) (hPe k hy)
    · have h1 : Tendsto (fun y => (continuousMultilinearCurryLeftEquiv ℝ
          (fun _ : Fin (k + 1) => E) F) (P (k + 1) y)) (𝓝[U] x)
          (𝓝 ((continuousMultilinearCurryLeftEquiv ℝ
            (fun _ : Fin (k + 1) => E) F) (P (k + 1) x))) :=
        ((continuousMultilinearCurryLeftEquiv ℝ (fun _ : Fin (k + 1) => E) F).continuous
          |>.tendsto _).comp hx'
      refine h1.congr' ?_
      filter_upwards [self_mem_nhdsWithin] with y hy
      have e : P k =ᶠ[𝓝 y] iteratedFDeriv ℝ k f :=
        Filter.eventuallyEq_of_mem (hU.mem_nhds hy) (hPe k)
      rw [e.fderiv_eq, fderiv_iteratedFDeriv, Function.comp_apply, hPe (k + 1) hy]
  have hT : HasFTaylorSeriesUpToOn ((⊤ : ℕ∞) : ℕ∞ω) (fun x => (P 0 x).curry0)
      (fun x m => P m x) (closure U) :=
    ⟨fun x _ => rfl, fun m _ x hx => hder m x hx, fun m _ => hPc m⟩
  refine ⟨fun x => (P 0 x).curry0, hT.contDiffOn, fun x hx => ?_⟩
  simp [hPe 0 hx]


/-- info: 'CWFlow.Interp.exists_contDiffOn_closure_of_bounded' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_contDiffOn_closure_of_bounded

end CWFlow.Interp
