import CWFlow.Flow.GaugeEvolution
import CWFlow.Flow.FlowS
import CWFlow.Flow.Potential
import CWFlow.Scalar.PotFamily

/-!
# The potential of the gauge section

`flowPotAt E J Φ s m` is `div(JH)` of the immersion `Φ_s` at `m ∈ M` (read in the chart at `m`), and
`gaugePot E J Φ ψ t x = −∫₀ᵗ flowPotAt E J Φ s (ψ_s x) ds`.

If `Φ_t(ψ_t(φ_α⁻¹ y)) = Ψ(y, ξ(t, y))` on `[0, T] × B` with `ξ(0, ·) = 0` (the gauge), then in the chart at `α`
(`gaugePot_chart`): `gaugePot` is jointly smooth on `[0, T] × B`, its spatial derivative is `ξ`, and its time
derivative is `−flowPotAt E J Φ t (ψ_t ·)`. So `Φ_t ∘ ψ_t` is the graph of the differential of `gaugePot t`.
For velocity `κ • flowVel` (`gaugePot_chartS`) the same holds for `κ * gaugePot`, with time derivative
`−κ flowPotAt …`.
-/

open Set Filter
open scoped InnerProductSpace RealInnerProductSpace Topology ContDiff Manifold

noncomputable section

set_option linter.unusedSectionVars false

namespace CWFlow

namespace Bridge

universe u

variable {E : Type u} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
  [NeZero (Module.finrank ℝ E)] [CompleteSpace E]
variable {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
  [CompactSpace M] [T2Space M] [SigmaCompactSpace M]
variable {V : Type u} [NormedAddCommGroup V] [InnerProductSpace ℝ V] [FiniteDimensional ℝ V]

variable (E) in
/-- `div(JH)` of the immersion `Φ_s` at the point `m`. -/
def flowPotAt (J : V →L[ℝ] V) (Φ : ℝ → M → V) (s : ℝ) (m : M) : ℝ :=
  flowPot J (chartRep E (Φ s) m) (chart E m m)

variable (E) in
/-- The potential `−∫₀ᵗ div(JH)[Φ_s](ψ_s x) ds`. -/
def gaugePot (J : V →L[ℝ] V) (Φ : ℝ → M → V) (ψ : ℝ → M → M) (t : ℝ) (x : M) : ℝ :=
  -∫ s in (0 : ℝ)..t, flowPotAt E J Φ s (ψ s x)

theorem gaugePot_zero (J : V →L[ℝ] V) (Φ : ℝ → M → V) (ψ : ℝ → M → M) (x : M) :
    gaugePot E J Φ ψ 0 x = 0 := by
  simp [gaugePot]

section PotS

variable {J : V →L[ℝ] V}
variable {T : ℝ} (hT : 0 < T) {κ : ℝ} {Φ : ℝ → M → V} (hΦ : IsFlowS E J κ T Φ) {ψ : ℝ → M → M}
  (hψ : ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, E) ∞ (fun q : ℝ × M => ψ q.1 q.2)
    (Icc 0 T ×ˢ univ))
  {Ψ : E × (E →L[ℝ] ℝ) → V} {Q : Set (E × (E →L[ℝ] ℝ))} (hQ : IsOpen Q) (hΨ : ContDiffOn ℝ ∞ Ψ Q)
  (hs : ∀ q ∈ Q, ∀ a b, ⟪J (fderiv ℝ Ψ q a), fderiv ℝ Ψ q b⟫ = omegaCanVal a b)
  {α : M} {B : Set E} (hB : IsOpen B) (hBU : B ⊆ (chart E α).target)
  {ξ : ℝ × E → (E →L[ℝ] ℝ)} (hξ : ContDiffOn ℝ ∞ ξ (Icc 0 T ×ˢ B))
  (hG : ∀ t ∈ Icc (0 : ℝ) T, ∀ y ∈ B, (y, ξ (t, y)) ∈ Q ∧
    Φ t (ψ t ((chart E α).symm y)) = Ψ (y, ξ (t, y)))
include hQ hΨ hs hB hξ hG

/-- The derivative of the section `Γ_t`. -/
theorem fderiv_gaugeSec {t : ℝ} (ht : t ∈ Icc (0 : ℝ) T) {z : E} (hz : z ∈ B) (b : E) :
    fderiv ℝ (gaugeSec Ψ ξ t) z b =
      fderiv ℝ Ψ (z, ξ (t, z)) (b, fderiv ℝ (fun z => ξ (t, z)) z b) := by
  have hξt : ContDiffOn ℝ ∞ (fun z => ξ (t, z)) B :=
    hξ.comp (contDiff_const.prodMk contDiff_id).contDiffOn (fun z hz => ⟨ht, hz⟩)
  have hΨz : HasFDerivAt Ψ (fderiv ℝ Ψ (z, ξ (t, z))) (z, ξ (t, z)) :=
    ((hΨ.contDiffAt (hQ.mem_nhds (hG t ht z hz).1)).differentiableAt (by simp)).hasFDerivAt
  have hξz := ((hξt.contDiffAt (hB.mem_nhds hz)).differentiableAt (by simp)).hasFDerivAt
  have hpair : HasFDerivAt (fun z => (z, ξ (t, z)))
      ((ContinuousLinearMap.id ℝ E).prod (fderiv ℝ (fun z => ξ (t, z)) z)) z :=
    (hasFDerivAt_id z).prodMk hξz
  have h : HasFDerivAt (gaugeSec Ψ ξ t) _ z :=
    HasFDerivAt.comp (g := Ψ) (f := fun z => (z, ξ (t, z))) z hΨz hpair
  rw [h.fderiv]
  rfl

theorem injective_fderiv_gaugeSec {t : ℝ} (ht : t ∈ Icc (0 : ℝ) T) {z : E} (hz : z ∈ B) :
    Function.Injective (fderiv ℝ (gaugeSec Ψ ξ t) z) :=
  injective_of_fderiv_prod (hs _ (hG t ht z hz).1) (H := fderiv ℝ (fun z => ξ (t, z)) z)
    (fderiv_gaugeSec hQ hΨ hs hB hξ hG ht hz)

include hT in
/-- `(s, y) ↦ div(JH)[Γ_s](y)` is jointly smooth. -/
theorem contDiffOn_flowPot_gaugeSec :
    ContDiffOn ℝ ∞ (fun p : ℝ × E => flowPot J (gaugeSec Ψ ξ p.1) p.2) (Icc 0 T ×ˢ B) := by
  have hw : ContDiffOn ℝ ∞ (fun p : ℝ × E => Ψ (p.2, ξ p)) (Icc 0 T ×ˢ B) :=
    ContDiffOn.comp (g := Ψ) (f := fun p : ℝ × E => (p.2, ξ p)) hΨ (contDiffOn_snd.prodMk hξ)
      (fun p hp => (hG p.1 hp.1 p.2 hp.2).1)
  have hinj : ∀ p ∈ Icc (0 : ℝ) T ×ˢ B, Function.Injective
      (fderiv ℝ (fun z => (fun p : ℝ × E => Ψ (p.2, ξ p)) (p.1, z)) p.2) := fun p hp =>
    injective_fderiv_gaugeSec hQ hΨ hs hB hξ hG hp.1 hp.2
  have h := contDiffOn_flowPot_family J hT hB hw hinj
  exact h

include hΦ hψ hBU

/-- `div(JH)` of the section is `div(JH)` of the flow at the image point. -/
theorem flowPot_gaugeSec {t : ℝ} (ht : t ∈ Icc (0 : ℝ) T) {y : E} (hy : y ∈ B) :
    flowPot J (gaugeSec Ψ ξ t) y = flowPotAt E J Φ t (ψ t ((chart E α).symm y)) := by
  obtain ⟨p, hp⟩ : ∃ p, p = ψ t ((chart E α).symm y) := ⟨_, rfl⟩
  have hpsrc : p ∈ (chart E p).source := mem_extChartAt_source (I := 𝓘(ℝ, E)) p
  have hβ : ψ t ((chart E α).symm y) ∈ (chart E p).source := by rw [← hp]; exact hpsrc
  have hyt : y ∈ (chart E α).target := hBU hy
  have hψt := contMDiff_slice hψ ht
  have hNo : IsOpen (B ∩ diffeoDom E (ψ t) α p) := hB.inter (isOpen_diffeoDom' hψt.continuous α p)
  have hτ : ContDiffOn ℝ ∞ (diffeoRep E (ψ t) α p) (B ∩ diffeoDom E (ψ t) α p) :=
    (contDiffOn_diffeoRep' hψt α p).mono inter_subset_right
  have heqz : ∀ z ∈ B ∩ diffeoDom E (ψ t) α p,
      gaugeSec Ψ ξ t z = chartRep E (Φ t) p (diffeoRep E (ψ t) α p z) := by
    intro z hz
    have h := chartRep_comp_diffeo (Φ t) hz.2
    rw [← h]
    show Ψ (z, ξ (t, z)) = Φ t (ψ t ((chart E α).symm z))
    exact (hG t ht z hz.1).2.symm
  have hev2 : ∀ z ∈ B ∩ diffeoDom E (ψ t) α p,
      gaugeSec Ψ ξ t =ᶠ[𝓝 z] chartRep E (Φ t) p ∘ diffeoRep E (ψ t) α p := fun z hz => by
    filter_upwards [hNo.mem_nhds hz] with w hw using heqz w hw
  have hρ : ∀ z ∈ B ∩ diffeoDom E (ψ t) α p, DifferentiableAt ℝ (diffeoRep E (ψ t) α p) z :=
    fun z hz => (hτ.contDiffAt (hNo.mem_nhds hz)).differentiableAt (by simp)
  have hXc : ∀ z ∈ B ∩ diffeoDom E (ψ t) α p,
      ContDiffAt ℝ ∞ (chartRep E (Φ t) p) (diffeoRep E (ψ t) α p z) := fun z hz =>
    ((hΦ.imm t ht).smooth p).contDiffAt ((isOpen_chart_target (E := E) p).mem_nhds
      ((chart E p).map_source hz.2.2))
  have hτinv : ∀ z ∈ B ∩ diffeoDom E (ψ t) α p,
      ∃ e : E ≃L[ℝ] E, (e : E →L[ℝ] E) = fderiv ℝ (diffeoRep E (ψ t) α p) z := by
    intro z hz
    have hcomp := fderiv_comp z ((hXc z hz).differentiableAt (by simp)) (hρ z hz)
    have hinj : Function.Injective (fderiv ℝ (diffeoRep E (ψ t) α p) z) := by
      intro b b' hab
      apply injective_fderiv_gaugeSec hQ hΨ hs hB hξ hG ht hz.1
      rw [(hev2 z hz).fderiv_eq, hcomp, ContinuousLinearMap.comp_apply,
        ContinuousLinearMap.comp_apply, hab]
    exact ⟨(LinearMap.linearEquivOfInjective (fderiv ℝ (diffeoRep E (ψ t) α p) z : E →ₗ[ℝ] E) hinj
      rfl).toContinuousLinearEquiv, ContinuousLinearMap.ext fun _ => rfl⟩
  have hyN : y ∈ B ∩ diffeoDom E (ψ t) α p := ⟨hy, hyt, hβ⟩
  have hτy : diffeoRep E (ψ t) α p y = chart E p p := by
    show chart E p (ψ t ((chart E α).symm y)) = chart E p p
    rw [← hp]
  have h := flowPot_comp hNo hτ hτinv hXc
    (fun z hz => (hΦ.imm t ht).inj p _ ((chart E p).map_source hz.2.2)) J hyN
  rw [hτy] at h
  rw [flowPot_congr_of_eventuallyEq (hev2 y hyN), h, ← hp]
  rfl

include hT

/-- **The potential in a chart**, for velocity `κ • flowVel` (`∂ₜξ = −κ D flowPot`). `κ * gaugePot` is jointly
smooth on `[0, T] × B`, its spatial derivative is the gauge section `ξ`, and its time derivative is
`−κ div(JH)`. -/
theorem gaugePot_chartS (hJJ : ∀ v, J (J v) = -v) (hξ0 : ∀ y ∈ B, ξ (0, y) = 0) :
    ContDiffOn ℝ ∞ (fun q : ℝ × E => κ * gaugePot E J Φ ψ q.1 ((chart E α).symm q.2)) (Icc 0 T ×ˢ B) ∧
    ∀ t ∈ Icc (0 : ℝ) T, ∀ y ∈ B,
      HasFDerivAt (fun z => κ * gaugePot E J Φ ψ t ((chart E α).symm z)) (ξ (t, y)) y ∧
      HasDerivWithinAt (fun s => κ * gaugePot E J Φ ψ s ((chart E α).symm y))
        (-(κ * flowPotAt E J Φ t (ψ t ((chart E α).symm y)))) (Icc 0 T) t := by
  have hθ0 := contDiffOn_flowPot_gaugeSec hT hQ hΨ hs hB hξ hG
  have hθ : ContDiffOn ℝ ∞ (fun p : ℝ × E => κ * flowPot J (gaugeSec Ψ ξ p.1) p.2)
      (Icc 0 T ×ˢ B) := contDiffOn_const.mul hθ0
  have heq : ∀ p ∈ Icc (0 : ℝ) T ×ˢ B, κ * gaugePot E J Φ ψ p.1 ((chart E α).symm p.2) =
      timePot (fun p : ℝ × E => κ * flowPot J (gaugeSec Ψ ξ p.1) p.2) p := by
    intro p hp
    show κ * -∫ s in (0 : ℝ)..p.1, flowPotAt E J Φ s (ψ s ((chart E α).symm p.2)) =
      -∫ s in (0 : ℝ)..p.1, κ * flowPot J (gaugeSec Ψ ξ s) p.2
    rw [intervalIntegral.integral_const_mul, mul_neg]
    congr 2
    refine intervalIntegral.integral_congr fun s hsI => ?_
    rw [uIcc_of_le hp.1.1] at hsI
    exact (flowPot_gaugeSec hΦ hψ hQ hΨ hs hB hBU hξ hG ⟨hsI.1, hsI.2.trans hp.1.2⟩ hp.2).symm
  have hθd : ∀ p ∈ Icc (0 : ℝ) T ×ˢ B,
      DifferentiableAt ℝ (fun z => κ * flowPot J (gaugeSec Ψ ξ p.1) z) p.2 := fun p hp =>
    ((ContDiffOn.comp (g := fun p : ℝ × E => κ * flowPot J (gaugeSec Ψ ξ p.1) p.2)
      (f := fun z : E => ((p.1, z) : ℝ × E)) hθ (contDiff_const.prodMk contDiff_id).contDiffOn
      (fun z hz => ⟨hp.1, hz⟩)).contDiffAt (hB.mem_nhds hp.2)).differentiableAt (by simp)
  have hD := (contDiffOn_slice_fderiv hT hB hθ).continuousOn
  have hξd : ∀ t ∈ Icc (0 : ℝ) T, ∀ y ∈ B, HasDerivWithinAt (fun s => ξ (s, y))
      (-fderiv ℝ (fun z => κ * flowPot J (gaugeSec Ψ ξ t) z) y) (Icc 0 T) t := by
    intro t ht y hy
    have h := hasDerivWithinAt_time hT hB hξ ht hy
    have hdiff : DifferentiableAt ℝ (fun z => flowPot J (gaugeSec Ψ ξ t) z) y :=
      ((ContDiffOn.comp (g := fun p : ℝ × E => flowPot J (gaugeSec Ψ ξ p.1) p.2)
        (f := fun z : E => ((t, z) : ℝ × E)) hθ0 (contDiff_const.prodMk contDiff_id).contDiffOn
        (fun z hz => ⟨ht, hz⟩)).contDiffAt (hB.mem_nhds hy)).differentiableAt (by simp)
    have hdt : dtW ξ T (t, y) = -fderiv ℝ (fun z => κ * flowPot J (gaugeSec Ψ ξ t) z) y := by
      refine ContinuousLinearMap.ext fun a => ?_
      rw [dtW_gauge_eqS hJJ hT hΦ hψ hQ hΨ hs hB hBU hξ hG ht hy a, fderiv_const_mul hdiff]
      simp
    rwa [hdt] at h
  have h2 : ∀ p ∈ Icc (0 : ℝ) T ×ˢ B, HasFDerivAt
      (fun z => timePot (fun p : ℝ × E => κ * flowPot J (gaugeSec Ψ ξ p.1) p.2) (p.1, z))
      (ξ p) p.2 :=
    fun p hp => hasFDerivAt_timePot hB hθ.continuousOn hθd hD hξ0 hξd hp.1 hp.2
  have h1 : ∀ p ∈ Icc (0 : ℝ) T ×ˢ B, HasDerivWithinAt
      (fun s => timePot (fun p : ℝ × E => κ * flowPot J (gaugeSec Ψ ξ p.1) p.2) (s, p.2))
      (-(κ * flowPot J (gaugeSec Ψ ξ p.1) p.2)) (Icc 0 T) p.1 :=
    fun p hp => hasDerivWithinAt_timePot hθ.continuousOn hp.1 hp.2
  refine ⟨(contDiffOn_of_partial hT hB hθ hξ h2 h1).congr heq, fun t ht y hy => ⟨?_, ?_⟩⟩
  · refine (h2 (t, y) ⟨ht, hy⟩).congr_of_eventuallyEq ?_
    filter_upwards [hB.mem_nhds hy] with z hz using heq (t, z) ⟨ht, hz⟩
  · have h := (h1 (t, y) ⟨ht, hy⟩).congr (fun s hsI => heq (s, y) ⟨hsI, hy⟩) (heq (t, y) ⟨ht, hy⟩)
    rw [flowPot_gaugeSec hΦ hψ hQ hΨ hs hB hBU hξ hG ht hy] at h
    exact h

end PotS

/-- **The potential in a chart.** `gaugePot` is jointly smooth on `[0, T] × B`, its spatial derivative is the
gauge section `ξ`, and its time derivative is `−div(JH)`. -/
theorem gaugePot_chart {J : V →L[ℝ] V}
    {T : ℝ} (hT : 0 < T) {Φ : ℝ → M → V} (hΦ : IsFlow E J T Φ) {ψ : ℝ → M → M}
    (hψ : ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, E) ∞ (fun q : ℝ × M => ψ q.1 q.2)
      (Icc 0 T ×ˢ univ))
    {Ψ : E × (E →L[ℝ] ℝ) → V} {Q : Set (E × (E →L[ℝ] ℝ))} (hQ : IsOpen Q) (hΨ : ContDiffOn ℝ ∞ Ψ Q)
    (hs : ∀ q ∈ Q, ∀ a b, ⟪J (fderiv ℝ Ψ q a), fderiv ℝ Ψ q b⟫ = omegaCanVal a b)
    {α : M} {B : Set E} (hB : IsOpen B) (hBU : B ⊆ (chart E α).target)
    {ξ : ℝ × E → (E →L[ℝ] ℝ)} (hξ : ContDiffOn ℝ ∞ ξ (Icc 0 T ×ˢ B))
    (hG : ∀ t ∈ Icc (0 : ℝ) T, ∀ y ∈ B, (y, ξ (t, y)) ∈ Q ∧
      Φ t (ψ t ((chart E α).symm y)) = Ψ (y, ξ (t, y)))
    (hJJ : ∀ v, J (J v) = -v) (hξ0 : ∀ y ∈ B, ξ (0, y) = 0) :
    ContDiffOn ℝ ∞ (fun q : ℝ × E => gaugePot E J Φ ψ q.1 ((chart E α).symm q.2)) (Icc 0 T ×ˢ B) ∧
    ∀ t ∈ Icc (0 : ℝ) T, ∀ y ∈ B,
      HasFDerivAt (fun z => gaugePot E J Φ ψ t ((chart E α).symm z)) (ξ (t, y)) y ∧
      HasDerivWithinAt (fun s => gaugePot E J Φ ψ s ((chart E α).symm y))
        (-flowPotAt E J Φ t (ψ t ((chart E α).symm y))) (Icc 0 T) t := by
  obtain ⟨h1, h2⟩ := gaugePot_chartS hT hΦ.toS hψ hQ hΨ hs hB hBU hξ hG hJJ hξ0
  simp only [one_mul] at h1 h2
  exact ⟨h1, h2⟩

end Bridge

end CWFlow
