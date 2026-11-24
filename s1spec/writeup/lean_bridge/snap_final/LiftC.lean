import CWFlow.Curve.CircleChart
import CWFlow.Extension.IntM

/-!
# Joint smoothness of the lifted family from chart readings (SCOPE W5, part of (T4))

If `(t, y) ↦ Φ t (chartInv ι1 z y)` is smooth on `[0, T'] × chartDomE ι1 z` for every chart `z` (the conclusion of the
repository's `Extension.fc_contDiffOn_chart_of_isFlow`), then `(t, θ) ↦ Φ t (e^{iθ})` is smooth near every
`(t₀, θ₀)` with `0 < t₀ < T'`.
-/

open Set Function Filter
open scoped Manifold Topology ContDiff

noncomputable section

namespace CWFlow.Curve

variable {Φ : ℝ → Circle → ℂ}

/-- Joint smoothness of the lift at `(t₀, θ₀)` from smoothness of the chart readings. -/
theorem contDiffAt_lift_of_chart {T' : ℝ}
    (hsm : ∀ z : Circle, ContDiffOn ℝ ∞ (fun q : ℝ × E1 => Φ q.1 (Extension.chartInv ι1 z q.2))
      (Icc 0 T' ×ˢ Extension.chartDomE ι1 z)) {t₀ : ℝ} (ht₀ : t₀ ∈ Ioo 0 T') (θ₀ : ℝ) :
    ContDiffAt ℝ ∞ (fun q : ℝ × ℝ => Φ q.1 (Circle.exp q.2)) (t₀, θ₀) := by
  set z := Circle.exp θ₀
  have hG := hsm z
  have hw : Aθ θ₀ ∈ liftDom z := exists_liftDom θ₀
  have hT := liftTau_mem_target hw
  have hm : ContDiffAt ℝ ∞ (fun q : ℝ × ℝ => ((q.1, liftTau z (Aθ q.2)) : ℝ × E1)) (t₀, θ₀) := by
    have h1 : ContDiffAt ℝ ∞ (liftTau z) (Aθ θ₀) :=
      (contDiffOn_liftTau z).contDiffAt ((isOpen_liftDom z).mem_nhds hw)
    have h3 : ContDiffAt ℝ ∞ (fun s : ℝ => liftTau z (Aθ s)) θ₀ :=
      ContDiffAt.comp (g := liftTau z) (f := fun s : ℝ => Aθ s) θ₀ h1 Aθ.contDiff.contDiffAt
    have h2 : ContDiffAt ℝ ∞ (fun q : ℝ × ℝ => liftTau z (Aθ q.2)) (t₀, θ₀) :=
      ContDiffAt.comp (g := fun s : ℝ => liftTau z (Aθ s)) (f := fun q : ℝ × ℝ => q.2) (t₀, θ₀) h3
        contDiffAt_snd
    exact contDiffAt_fst.prodMk h2
  have hmem : Icc 0 T' ×ˢ Extension.chartDomE ι1 z ∈ 𝓝 ((t₀, liftTau z (Aθ θ₀)) : ℝ × E1) :=
    prod_mem_nhds (Icc_mem_nhds ht₀.1 ht₀.2) ((isOpen_extChartAt_target z).mem_nhds hT)
  have hc := (hG.contDiffAt hmem).comp (t₀, θ₀) hm
  refine hc.congr_of_eventuallyEq ?_
  have hopen : IsOpen {q : ℝ × ℝ | Aθ q.2 ∈ liftDom z} :=
    (isOpen_liftDom z).preimage (Aθ.continuous.comp continuous_snd)
  filter_upwards [hopen.mem_nhds (show Aθ θ₀ ∈ liftDom z from hw)] with q hq
  have hq' : Aθ q.2 ∈ liftDom z := hq
  show Φ q.1 (Circle.exp q.2) = Φ q.1 (Extension.chartInv ι1 z (liftTau z (Aθ q.2)))
  have e : Extension.chartInv ι1 z (liftTau z (Aθ q.2)) = Circle.exp q.2 := by
    have := chartSymm_liftTau hq'
    rw [Aθ_symm_Aθ] at this
    exact this
  rw [e]

end CWFlow.Curve
