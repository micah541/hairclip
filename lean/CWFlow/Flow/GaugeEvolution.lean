import CWFlow.Flow.GraphUnique
import CWFlow.Flow.FlowS

/-!
# The evolution of the gauge section

Let `Φ` be a flow and suppose that, in the chart at `α` and for `y` in an open set `B`,
`Φ_t(ψ_t(φ_α⁻¹ y)) = Ψ(y, ξ(t, y))` for a symplectic map `Ψ` (cotangent coordinates) and a jointly smooth
family of maps `ψ_t`. Then (`dtW_gauge_eq`)

  `∂ₜξ = −D(flowPot[Γ_t])`,  `Γ_t = Ψ(·, ξ(t, ·))`.

More generally (`dtW_gauge_eqS`), if `Φ` has velocity `κ • flowVel` (`IsFlowS`), then `∂ₜξ = −κ D(flowPot[Γ_t])`.

Proof: `∂ₜΓ = DΨ(0, ∂ₜξ)` is vertical; on the other hand `Γ_t = Φ_t ∘ ψ_t`, so `∂ₜΓ` is the flow velocity of
`Φ_t` (which is the flow velocity `J DΓ g⁻¹ D flowPot` of `Γ_t`, by invariance) plus a vector tangent to `Φ_t`.
Pairing with `J ·` against `DΓ_t a`: the vertical vector gives `∂ₜξ(a)`, the flow velocity `−D flowPot(a)`,
and the tangent vector `0` (Lagrangian).
-/

open Set Filter
open scoped InnerProductSpace RealInnerProductSpace Topology ContDiff Manifold

noncomputable section

set_option linter.unusedSectionVars false

namespace CWFlow

section Mono

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
variable {M : Type*} [TopologicalSpace M] [ChartedSpace E M]
variable {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]

/-- A flow on `[0, T]` is a flow on `[0, δ]`, `δ ≤ T`. -/
theorem IsFlow.mono {J : V →L[ℝ] V} {T δ : ℝ} {Φ : ℝ → M → V} (h : IsFlow E J T Φ) (hδ : δ ≤ T) :
    IsFlow E J δ Φ where
  smooth := h.smooth.mono (prod_mono (Icc_subset_Icc_right hδ) subset_rfl)
  imm := fun t ht => h.imm t ⟨ht.1, ht.2.trans hδ⟩
  vel := fun t ht x z hx => (h.vel t ⟨ht.1, ht.2.trans hδ⟩ x z hx).mono (Icc_subset_Icc_right hδ)

end Mono

namespace Bridge

universe u

variable {E : Type u} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
  [NeZero (Module.finrank ℝ E)] [CompleteSpace E]
variable {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
  [CompactSpace M] [T2Space M] [SigmaCompactSpace M]
variable {V : Type u} [NormedAddCommGroup V] [InnerProductSpace ℝ V] [FiniteDimensional ℝ V]
variable {J : V →L[ℝ] V}

/-- The section `Γ_t = Ψ(·, ξ(t, ·))`. -/
abbrev gaugeSec (Ψ : E × (E →L[ℝ] ℝ) → V) (ξ : ℝ × E → (E →L[ℝ] ℝ)) (t : ℝ) (z : E) : V :=
  Ψ (z, ξ (t, z))

/-- **Evolution of the gauge section** for velocity `κ • flowVel`: `∂ₜξ = −κ D(flowPot[Γ_t])`. -/
theorem dtW_gauge_eqS {κ : ℝ} (hJJ : ∀ v, J (J v) = -v)
    {T : ℝ} (hT : 0 < T) {Φ : ℝ → M → V} (hΦ : IsFlowS E J κ T Φ) {ψ : ℝ → M → M}
    (hψ : ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, E) ∞ (fun q : ℝ × M => ψ q.1 q.2)
      (Icc 0 T ×ˢ univ))
    {Ψ : E × (E →L[ℝ] ℝ) → V} {Q : Set (E × (E →L[ℝ] ℝ))} (hQ : IsOpen Q) (hΨ : ContDiffOn ℝ ∞ Ψ Q)
    (hs : ∀ q ∈ Q, ∀ a b, ⟪J (fderiv ℝ Ψ q a), fderiv ℝ Ψ q b⟫ = omegaCanVal a b)
    {α : M} {B : Set E} (hB : IsOpen B) (hBU : B ⊆ (chart E α).target)
    {ξ : ℝ × E → (E →L[ℝ] ℝ)} (hξ : ContDiffOn ℝ ∞ ξ (Icc 0 T ×ˢ B))
    (hG : ∀ t ∈ Icc (0 : ℝ) T, ∀ y ∈ B, (y, ξ (t, y)) ∈ Q ∧
      Φ t (ψ t ((chart E α).symm y)) = Ψ (y, ξ (t, y)))
    {t : ℝ} (ht : t ∈ Icc (0 : ℝ) T) {y : E} (hy : y ∈ B) (a : E) :
    dtW ξ T (t, y) a = -(κ * fderiv ℝ (flowPot J (fun z => Ψ (z, ξ (t, z)))) y a) := by
  obtain ⟨p, hp⟩ : ∃ p, p = ψ t ((chart E α).symm y) := ⟨_, rfl⟩
  have hpsrc : p ∈ (chart E p).source := mem_extChartAt_source (I := 𝓘(ℝ, E)) p
  have hβ : ψ t ((chart E α).symm y) ∈ (chart E p).source := by rw [← hp]; exact hpsrc
  have hyt : y ∈ (chart E α).target := hBU hy
  have hc : chart E p p ∈ (chart E p).target := (chart E p).map_source hpsrc
  have hceq : psiRep E ψ α p (t, y) = chart E p p := by
    show chart E p (ψ t ((chart E α).symm y)) = chart E p p
    rw [← hp]
  have hSsub : Icc (0 : ℝ) T ×ˢ B ⊆ Icc (0 : ℝ) T ×ˢ (chart E α).target := prod_mono subset_rfl hBU
  have hUD : UniqueDiffWithinAt ℝ (Icc (0 : ℝ) T ×ˢ B) (t, y) :=
    ((uniqueDiffOn_Icc hT).prod hB.uniqueDiffOn) _ ⟨ht, hy⟩
  have hq₀ : (y, ξ (t, y)) ∈ Q := (hG t ht y hy).1
  -- the derivative of `Γ(s, z) = Ψ(z, ξ(s, z))`
  have hΨd : HasFDerivAt Ψ (fderiv ℝ Ψ (y, ξ (t, y))) (y, ξ (t, y)) :=
    ((hΨ.contDiffAt (hQ.mem_nhds hq₀)).differentiableAt (by simp)).hasFDerivAt
  have hξd := hasFDerivWithinAt_prod hT hB hξ ht hy
  have hin : HasFDerivWithinAt (fun q : ℝ × E => (q.2, ξ q))
      ((ContinuousLinearMap.snd ℝ ℝ E).prod (fderivWithin ℝ ξ (Icc 0 T ×ˢ B) (t, y)))
      (Icc 0 T ×ˢ B) (t, y) := hasFDerivWithinAt_snd.prodMk hξd
  have hA := HasFDerivAt.comp_hasFDerivWithinAt (g := Ψ) (f := fun q : ℝ × E => (q.2, ξ q)) (t, y)
    hΨd hin
  -- the same map through the flow
  have hP := (contDiffWithinAt_psiRep hψ ht hyt hβ).mono hSsub
  have hPd : HasFDerivWithinAt (psiRep E ψ α p)
      (fderivWithin ℝ (psiRep E ψ α p) (Icc 0 T ×ˢ B) (t, y)) (Icc 0 T ×ˢ B) (t, y) :=
    (hP.differentiableWithinAt (by simp)).hasFDerivWithinAt
  have hN : {q : ℝ × E | ψ q.1 ((chart E α).symm q.2) ∈ (chart E p).source} ∈
      𝓝[Icc 0 T ×ˢ B] (t, y) :=
    nhdsWithin_mono _ hSsub
      ((contMDiffOn_psiPre hψ α (t, y) ⟨ht, hyt⟩).continuousWithinAt.preimage_mem_nhdsWithin
        ((isOpen_extChartAt_source (I := 𝓘(ℝ, E)) p).mem_nhds hβ))
  have hΦrep := contDiffOn_flowRep hΦ.smooth p
  have hcp : psiRep E ψ α p (t, y) ∈ (chart E p).target := by rw [hceq]; exact hc
  have hΦd := hasFDerivWithinAt_prod hT (isOpen_chart_target (E := E) p) hΦrep ht hcp
  have hκ : HasFDerivWithinAt (fun q : ℝ × E => (q.1, psiRep E ψ α p q))
      ((ContinuousLinearMap.fst ℝ ℝ E).prod
        (fderivWithin ℝ (psiRep E ψ α p) (Icc 0 T ×ˢ B) (t, y)))
      (Icc 0 T ×ˢ B) (t, y) := hasFDerivWithinAt_fst.prodMk hPd
  have hmaps : MapsTo (fun q : ℝ × E => (q.1, psiRep E ψ α p q))
      ((Icc 0 T ×ˢ B) ∩ {q : ℝ × E | ψ q.1 ((chart E α).symm q.2) ∈ (chart E p).source})
      (Icc 0 T ×ˢ (chart E p).target) := fun q hq => ⟨hq.1.1, (chart E p).map_source hq.2⟩
  have h1 := HasFDerivWithinAt.comp (g := flowRep E Φ p)
    (f := fun q : ℝ × E => (q.1, psiRep E ψ α p q)) (t, y) hΦd (hκ.mono inter_subset_left) hmaps
  have h2 := h1.mono_of_mem_nhdsWithin (inter_mem self_mem_nhdsWithin hN)
  have hev : (Ψ ∘ fun q : ℝ × E => (q.2, ξ q)) =ᶠ[𝓝[Icc 0 T ×ˢ B] (t, y)]
      flowRep E Φ p ∘ fun q : ℝ × E => (q.1, psiRep E ψ α p q) := by
    filter_upwards [inter_mem self_mem_nhdsWithin hN] with q hq
    simp only [flowRep, psiRep, Function.comp_apply]
    rw [(chart E p).left_inv hq.2]
    exact (hG q.1 hq.1.1 q.2 hq.1.2).2.symm
  have h3 := h2.congr_of_eventuallyEq hev (hev.self_of_nhdsWithin ⟨ht, hy⟩)
  have hEq := hUD.eq hA h3
  have hEq1 : fderiv ℝ Ψ (y, ξ (t, y)) (0, fderivWithin ℝ ξ (Icc 0 T ×ˢ B) (t, y) (1, 0)) =
      fderivWithin ℝ (flowRep E Φ p) (Icc 0 T ×ˢ (chart E p).target) (t, psiRep E ψ α p (t, y))
        (1, fderivWithin ℝ (psiRep E ψ α p) (Icc 0 T ×ˢ B) (t, y) (1, 0)) := by
    have h := congrArg (fun L : ℝ × E →L[ℝ] V => L (1, 0)) hEq
    exact h
  rw [fderivWithin_prod_apply hT hB hξ ht hy,
    fderivWithin_prod_apply hT (isOpen_chart_target (E := E) p) hΦrep ht hcp, one_smul, one_smul,
    hceq] at hEq1
  simp only [map_zero, add_zero] at hEq1
  have hdt : dtW (flowRep E Φ p) T (t, chart E p p) =
      κ • flowVel J (chartRep E (Φ t) p) (chart E p p) := by
    have hfun : (fun s => flowRep E Φ p (s, chart E p p)) = fun s => Φ s p := by
      funext s
      simp only [flowRep]
      rw [(chart E p).left_inv hpsrc]
    show derivWithin (fun s => flowRep E Φ p (s, chart E p p)) (Icc 0 T) t = _
    rw [hfun]
    exact (hΦ.vel t ht p p hpsrc).derivWithin (uniqueDiffOn_Icc hT t ht)
  -- the slice `Γ_t` and its derivative
  have hξt : ContDiffOn ℝ ∞ (fun z => ξ (t, z)) B :=
    hξ.comp (contDiff_const.prodMk contDiff_id).contDiffOn (fun z hz => ⟨ht, hz⟩)
  have hΓ : ∀ z ∈ B, HasFDerivAt (gaugeSec Ψ ξ t)
      ((fderiv ℝ Ψ (z, ξ (t, z))).comp
        ((ContinuousLinearMap.id ℝ E).prod (fderiv ℝ (fun z => ξ (t, z)) z))) z := by
    intro z hz
    have hΨz : HasFDerivAt Ψ (fderiv ℝ Ψ (z, ξ (t, z))) (z, ξ (t, z)) :=
      ((hΨ.contDiffAt (hQ.mem_nhds (hG t ht z hz).1)).differentiableAt (by simp)).hasFDerivAt
    have hξz := ((hξt.contDiffAt (hB.mem_nhds hz)).differentiableAt (by simp)).hasFDerivAt
    have hpair : HasFDerivAt (fun z => (z, ξ (t, z)))
        ((ContinuousLinearMap.id ℝ E).prod (fderiv ℝ (fun z => ξ (t, z)) z)) z :=
      (hasFDerivAt_id z).prodMk hξz
    exact HasFDerivAt.comp (g := Ψ) (f := fun z => (z, ξ (t, z))) z hΨz hpair
  have hLz : ∀ z ∈ B, ∀ b, fderiv ℝ (gaugeSec Ψ ξ t) z b =
      fderiv ℝ Ψ (z, ξ (t, z)) (b, fderiv ℝ (fun z => ξ (t, z)) z b) := by
    intro z hz b
    rw [(hΓ z hz).fderiv]
    rfl
  have hinjz : ∀ z ∈ B, Function.Injective (fderiv ℝ (gaugeSec Ψ ξ t) z) := fun z hz =>
    injective_of_fderiv_prod (hs _ (hG t ht z hz).1) (H := fderiv ℝ (fun z => ξ (t, z)) z)
      (hLz z hz)
  -- `Γ_t = Φ_t ∘ ψ_t` in charts
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
  have hXd : ∀ z ∈ B ∩ diffeoDom E (ψ t) α p,
      DifferentiableAt ℝ (chartRep E (Φ t) p) (diffeoRep E (ψ t) α p z) := fun z hz =>
    (((hΦ.imm t ht).smooth p).contDiffAt ((isOpen_chart_target (E := E) p).mem_nhds
      ((chart E p).map_source hz.2.2))).differentiableAt (by simp)
  have hτinv : ∀ z ∈ B ∩ diffeoDom E (ψ t) α p,
      ∃ e : E ≃L[ℝ] E, (e : E →L[ℝ] E) = fderiv ℝ (diffeoRep E (ψ t) α p) z := by
    intro z hz
    have hcomp := fderiv_comp z (hXd z hz) (hρ z hz)
    have hinj : Function.Injective (fderiv ℝ (diffeoRep E (ψ t) α p) z) := by
      intro b b' hab
      apply hinjz z hz.1
      rw [(hev2 z hz).fderiv_eq, hcomp, ContinuousLinearMap.comp_apply,
        ContinuousLinearMap.comp_apply, hab]
    exact ⟨(LinearMap.linearEquivOfInjective (fderiv ℝ (diffeoRep E (ψ t) α p) z : E →ₗ[ℝ] E) hinj
      rfl).toContinuousLinearEquiv, ContinuousLinearMap.ext fun _ => rfl⟩
  have hyN : y ∈ B ∩ diffeoDom E (ψ t) α p := ⟨hy, hyt, hβ⟩
  have hτy : diffeoRep E (ψ t) α p y = chart E p p := by
    show chart E p (ψ t ((chart E α).symm y)) = chart E p p
    rw [← hp]
  have hfv : flowVel J (gaugeSec Ψ ξ t) y = flowVel J (chartRep E (Φ t) p) (chart E p p) := by
    have h := flowVel_comp' (isOpen_chart_target (E := E) p) ((hΦ.imm t ht).smooth p)
      ((hΦ.imm t ht).inj p) J hNo hτ hτinv (fun w hw => (chart E p).map_source hw.2.2) hyN
    rw [hτy] at h
    rw [flowVel_congr_of_eventuallyEq (hev2 y hyN), h]
  have hLa : fderiv ℝ (gaugeSec Ψ ξ t) y a =
      fderiv ℝ (chartRep E (Φ t) p) (chart E p p) (fderiv ℝ (diffeoRep E (ψ t) α p) y a) := by
    rw [(hev2 y hyN).fderiv_eq, fderiv_comp y (hXd y hyN) (hρ y hyN),
      ContinuousLinearMap.comp_apply, hτy]
  -- pair with `J ·` against `DΓ_t a`
  rw [hdt, ← hfv] at hEq1
  have hX' : (fun z => flowRep E Φ p (t, z)) = chartRep E (Φ t) p := rfl
  rw [hX'] at hEq1
  have key := congrArg (fun v : V => ⟪J v, fderiv ℝ (gaugeSec Ψ ξ t) y a⟫) hEq1
  have hl : ⟪J (fderiv ℝ Ψ (y, ξ (t, y)) (0, dtW ξ T (t, y))), fderiv ℝ (gaugeSec Ψ ξ t) y a⟫ =
      dtW ξ T (t, y) a := by
    rw [hLz y hy a]
    exact inner_J_fderiv_vertical (hs _ hq₀) (dtW ξ T (t, y)) (a, fderiv ℝ (fun z => ξ (t, z)) y a)
  have hr1 : ⟪J (flowVel J (gaugeSec Ψ ξ t) y), fderiv ℝ (gaugeSec Ψ ξ t) y a⟫ =
      -fderiv ℝ (flowPot J (gaugeSec Ψ ξ t)) y a := by
    show ⟪J (J (fderiv ℝ (gaugeSec Ψ ξ t) y (ginvOf (fderiv ℝ (gaugeSec Ψ ξ t) y)
      (fderiv ℝ (flowPot J (gaugeSec Ψ ξ t)) y)))), fderiv ℝ (gaugeSec Ψ ξ t) y a⟫ = _
    rw [hJJ, inner_neg_left, inner_ginvOf (hinjz y hy)]
  have hr2 : ∀ v : E, ⟪J (fderiv ℝ (chartRep E (Φ t) p) (chart E p p) v),
      fderiv ℝ (gaugeSec Ψ ξ t) y a⟫ = 0 := fun v => by
    rw [hLa]
    exact (hΦ.imm t ht).lag p _ hc v _
  rw [hl, map_add, inner_add_left, hr2, add_zero, map_smul, real_inner_smul_left, hr1] at key
  rw [key]
  ring

/-- **Evolution of the gauge section**: `∂ₜξ = −D(flowPot[Γ_t])`. -/
theorem dtW_gauge_eq (hJJ : ∀ v, J (J v) = -v)
    {T : ℝ} (hT : 0 < T) {Φ : ℝ → M → V} (hΦ : IsFlow E J T Φ) {ψ : ℝ → M → M}
    (hψ : ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, E) ∞ (fun q : ℝ × M => ψ q.1 q.2)
      (Icc 0 T ×ˢ univ))
    {Ψ : E × (E →L[ℝ] ℝ) → V} {Q : Set (E × (E →L[ℝ] ℝ))} (hQ : IsOpen Q) (hΨ : ContDiffOn ℝ ∞ Ψ Q)
    (hs : ∀ q ∈ Q, ∀ a b, ⟪J (fderiv ℝ Ψ q a), fderiv ℝ Ψ q b⟫ = omegaCanVal a b)
    {α : M} {B : Set E} (hB : IsOpen B) (hBU : B ⊆ (chart E α).target)
    {ξ : ℝ × E → (E →L[ℝ] ℝ)} (hξ : ContDiffOn ℝ ∞ ξ (Icc 0 T ×ˢ B))
    (hG : ∀ t ∈ Icc (0 : ℝ) T, ∀ y ∈ B, (y, ξ (t, y)) ∈ Q ∧
      Φ t (ψ t ((chart E α).symm y)) = Ψ (y, ξ (t, y)))
    {t : ℝ} (ht : t ∈ Icc (0 : ℝ) T) {y : E} (hy : y ∈ B) (a : E) :
    dtW ξ T (t, y) a = -fderiv ℝ (flowPot J (fun z => Ψ (z, ξ (t, z)))) y a := by
  have h := dtW_gauge_eqS (κ := 1) hJJ hT hΦ.toS hψ hQ hΨ hs hB hBU hξ hG ht hy a
  rw [h, one_mul]

end Bridge

end CWFlow
