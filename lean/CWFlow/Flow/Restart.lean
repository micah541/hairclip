import CWFlow.Flow.NormalTransport
import CWFlow.Flow.RestartScalar
import CWFlow.Flow.UniqueMain
import CWFlow.Flow.GaugeGraph
import CWFlow.Flow.GraphImm
import CWFlow.Flow.InducedMetric
import CWFlow.Bridge.Existence

/-!
# Restarting a flow beyond its time of existence

`isFlow_restart`: a flow `Φ` of Lagrangian immersions of a compact manifold on `[0, T]` extends to a flow on
`[0, T + ε]` for some `ε > 0`, with the same values on `[0, T]` (Chen-Warren, blueprint step 5, "restart without
matching derivatives").

Proof.
* The reversed path `Γ_s = Φ_{T-s}` is an `IsFlowS (-1)` (`rst_reverse`), so by the gauge theorem
  (`exists_gauge_graphS`) for a short time `Γ_s ∘ ψ_s = graphImm (w_s)` over `Γ_0 = Φ_T`, with `w` jointly smooth,
  `w_0 = 0`; the flow `χ'` of `tanFieldN` of `G'_s = graphImm (w_s)` inverts `ψ` (`psi_comp_chiN`), so
  `G'_s ∘ χ'_s = Γ_s` (`rst_transport`).
* The scalar data `w_s` tend to the datum `0` as `s → 0⁺`; by continuous dependence for the Mantegazza-Martinazzi
  equation, for small `s` the scalar equation has a solution from `w_s` on a time interval `[0, T₁]` with `T₁`
  independent of `s`, hence a flow `Θ''` from `G'_s` (`rst_scalar`); `Θ_t = Θ''_t ∘ χ'_s` is a flow from `Φ_{T-s}`.
* By uniqueness `Θ_t = Φ_{T-s+t}` on `[0, s]`, so `Φ` followed by `Θ` is a flow on `[0, T + T₁ - s]` (`rst_glue`).
-/

open Set Filter Bundle
open DifferentialGeometry DifferentialGeometry.Tensor0SBundle ShortTime ShortTime.HFlow CWFlow.JetRec
open DifferentialGeometry.Geometry.Curvature DifferentialGeometry.Geometry.Connection
open DifferentialGeometry.Integral.L2
open DifferentialGeometry.Analysis.Sobolev DifferentialGeometry.Analysis.Spectral
open DifferentialGeometry.Analysis.Parabolic.TensorSpectral DifferentialGeometry.Analysis.Elliptic
open scoped InnerProductSpace RealInnerProductSpace Topology ContDiff Manifold
open CWFlow.Bridge

noncomputable section

set_option linter.unusedSectionVars false

namespace CWFlow

universe u

variable {E : Type u} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
  [NeZero (Module.finrank ℝ E)]
variable {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
  [CompactSpace M] [T2Space M] [SigmaCompactSpace M]
variable {V : Type u} [NormedAddCommGroup V] [InnerProductSpace ℝ V] [FiniteDimensional ℝ V]

/-! ### Time reversal and time shift -/

/-- The time reversal of a flow on `[0, T]` is an `IsFlowS (-1)`. -/
theorem rst_reverse {J : V →L[ℝ] V} {T : ℝ} {Φ : ℝ → M → V} (h : IsFlow E J T Φ) :
    IsFlowS E J (-1) T (fun s => Φ (T - s)) where
  smooth := by
    have h1 : ContMDiff (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, ℝ) ∞ (fun q : ℝ × M => T - q.1) :=
      ContMDiff.comp (g := fun s : ℝ => T - s) (f := fun q : ℝ × M => q.1)
        (contDiff_const.sub contDiff_id).contMDiff contMDiff_fst
    have h2 : ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) ∞
        (fun q : ℝ × M => ((T - q.1, q.2) : ℝ × M)) (Icc 0 T ×ˢ univ) :=
      h1.contMDiffOn.prodMk contMDiffOn_snd
    exact ContMDiffOn.comp (g := fun q : ℝ × M => Φ q.1 q.2)
      (f := fun q : ℝ × M => ((T - q.1, q.2) : ℝ × M)) h.smooth h2
      (fun q hq => ⟨⟨by linarith [hq.1.2], by linarith [hq.1.1]⟩, mem_univ _⟩)
  imm := fun t ht => h.imm (T - t) ⟨by linarith [ht.2], by linarith [ht.1]⟩
  vel := fun t ht x z hx => by
    have h1 : HasDerivWithinAt (fun s : ℝ => T - s) (-1) (Icc 0 T) t :=
      (hasDerivWithinAt_id t _).const_sub T
    exact HasDerivWithinAt.scomp t
      (h.vel (T - t) ⟨by linarith [ht.2], by linarith [ht.1]⟩ x z hx) h1
      (fun s hs => ⟨by linarith [hs.2], by linarith [hs.1]⟩)

/-- A flow on `[0, T₁]`, restarted to the time interval `[a, a + T₁]`. -/
theorem rst_shiftR {J : V →L[ℝ] V} {T₁ a : ℝ} {Θ : ℝ → M → V} (h : IsFlow E J T₁ Θ) :
    ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, V) ∞ (fun q : ℝ × M => Θ (q.1 - a) q.2)
      (Icc a (a + T₁) ×ˢ univ) ∧
    ∀ t ∈ Icc a (a + T₁), ∀ x z : M, x ∈ (chart E z).source →
      HasDerivWithinAt (fun s => Θ (s - a) x)
        (flowVel J (chartRep E (Θ (t - a)) z) (chart E z x)) (Icc a (a + T₁)) t := by
  refine ⟨?_, fun t ht x z hx => ?_⟩
  · have h1 : ContMDiff (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, ℝ) ∞ (fun q : ℝ × M => q.1 - a) :=
      ContMDiff.comp (g := fun s : ℝ => s - a) (f := fun q : ℝ × M => q.1)
        (contDiff_id.sub contDiff_const).contMDiff contMDiff_fst
    have h2 : ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) ∞
        (fun q : ℝ × M => ((q.1 - a, q.2) : ℝ × M)) (Icc a (a + T₁) ×ˢ univ) :=
      h1.contMDiffOn.prodMk contMDiffOn_snd
    exact ContMDiffOn.comp (g := fun q : ℝ × M => Θ q.1 q.2)
      (f := fun q : ℝ × M => ((q.1 - a, q.2) : ℝ × M)) h.smooth h2
      (fun q hq => ⟨⟨by linarith [hq.1.1], by linarith [hq.1.2]⟩, mem_univ _⟩)
  · have h1 : HasDerivWithinAt (fun s : ℝ => s - a) 1 (Icc a (a + T₁)) t :=
      (hasDerivWithinAt_id t _).sub_const a
    have hm : MapsTo (fun s : ℝ => s - a) (Icc a (a + T₁)) (Icc 0 T₁) :=
      fun s hs => ⟨by linarith [hs.1], by linarith [hs.2]⟩
    have h2 := HasDerivWithinAt.scomp t
      (h.vel (t - a) ⟨by linarith [ht.1], by linarith [ht.2]⟩ x z hx) h1 hm
    rw [one_smul] at h2
    exact h2

/-! ### Gluing -/

/-- Local nature of joint smoothness in time. -/
theorem rst_local_smooth {S S' : Set ℝ} {Ψ Φ' : ℝ → M → V} (hag : ∀ s ∈ S', Ψ s = Φ' s) {t : ℝ}
    (ht : t ∈ S') (hn : S' ∈ 𝓝[S] t) (x : M)
    (h : ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, V) ∞ (fun q : ℝ × M => Φ' q.1 q.2) (S' ×ˢ univ)) :
    ContMDiffWithinAt (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, V) ∞ (fun q : ℝ × M => Ψ q.1 q.2)
      (S ×ˢ univ) (t, x) := by
  have h1 := (h (t, x) ⟨ht, mem_univ _⟩).congr
    (f₁ := fun q : ℝ × M => Ψ q.1 q.2) (fun q hq => by simp only [hag q.1 hq.1])
    (by simp only [hag t ht])
  refine h1.mono_of_mem_nhdsWithin ?_
  rw [nhdsWithin_prod_eq]
  exact prod_mem_prod hn (by rw [nhdsWithin_univ]; exact univ_mem)

/-- Local nature of the time derivative. -/
theorem rst_local_vel {S S' : Set ℝ} {Ψ Φ' : ℝ → M → V} {x : M} {t : ℝ} {v : V}
    (hag : ∀ s ∈ S', Ψ s x = Φ' s x) (ht : t ∈ S') (hn : S' ∈ 𝓝[S] t)
    (h : HasDerivWithinAt (fun s => Φ' s x) v S' t) : HasDerivWithinAt (fun s => Ψ s x) v S t :=
  (h.congr (fun s hs => hag s hs) (hag t ht)).mono_of_mem_nhdsWithin hn

/-- **Gluing**: `Φ` on `[0, a + s]` followed by `Θ(· - a)` is a flow on `[0, a + T₁]`, if `Θ` is a flow on
`[0, T₁]` agreeing with `Φ(a + ·)` on `[0, s]`, `0 < s ≤ T₁`. -/
theorem rst_glue {J : V →L[ℝ] V} {a s T₁ : ℝ} (ha : 0 ≤ a) (hs : 0 < s) (hsT : s ≤ T₁)
    {Φ Θ : ℝ → M → V} (hΦ : IsFlow E J (a + s) Φ) (hΘ : IsFlow E J T₁ Θ)
    (hag : ∀ t ∈ Icc (0 : ℝ) s, Θ t = Φ (a + t)) :
    ∃ Ψ : ℝ → M → V, IsFlow E J (a + T₁) Ψ ∧ ∀ t ∈ Icc (0 : ℝ) (a + s), Ψ t = Φ t := by
  obtain ⟨Ψ, hL, hR'⟩ : ∃ Ψ : ℝ → M → V, (∀ t, t ≤ a + s → Ψ t = Φ t) ∧
      ∀ t, a + s < t → Ψ t = Θ (t - a) :=
    ⟨fun t => if t ≤ a + s then Φ t else Θ (t - a), fun t ht => if_pos ht,
      fun t ht => if_neg (not_le.2 ht)⟩
  have hR : ∀ t ∈ Icc a (a + T₁), Ψ t = Θ (t - a) := by
    intro t ht
    by_cases h : t ≤ a + s
    · rw [hL t h, hag (t - a) ⟨by linarith [ht.1], by linarith⟩, add_sub_cancel]
    · exact hR' t (not_le.1 h)
  obtain ⟨hsm, hvel⟩ := rst_shiftR (a := a) hΘ
  refine ⟨Ψ, ⟨?_, ?_, ?_⟩, fun t ht => hL t ht.2⟩
  · rintro ⟨t, x⟩ ⟨ht, -⟩
    by_cases h : t < a + s
    · exact rst_local_smooth (S' := Icc 0 (a + s)) (fun r hr => hL r hr.2) ⟨ht.1, h.le⟩
        (mem_nhdsWithin.2 ⟨Iio (a + s), isOpen_Iio, h, fun r hr => ⟨hr.2.1, hr.1.le⟩⟩) x hΦ.smooth
    · have h' : a < t := by linarith [not_lt.1 h]
      exact rst_local_smooth (S' := Icc a (a + T₁)) (Φ' := fun r => Θ (r - a)) hR ⟨h'.le, ht.2⟩
        (mem_nhdsWithin.2 ⟨Ioi a, isOpen_Ioi, h', fun r hr => ⟨hr.1.le, hr.2.2⟩⟩) x hsm
  · intro t ht
    by_cases h : t ≤ a + s
    · rw [hL t h]
      exact hΦ.imm t ⟨ht.1, h⟩
    · rw [hR t ⟨by linarith [not_le.1 h], ht.2⟩]
      exact hΘ.imm (t - a) ⟨by linarith [not_le.1 h], by linarith [ht.2]⟩
  · intro t ht x z hx
    by_cases h : t < a + s
    · rw [hL t h.le]
      exact rst_local_vel (S' := Icc 0 (a + s)) (fun r hr => by rw [hL r hr.2]) ⟨ht.1, h.le⟩
        (mem_nhdsWithin.2 ⟨Iio (a + s), isOpen_Iio, h, fun r hr => ⟨hr.2.1, hr.1.le⟩⟩)
        (hΦ.vel t ⟨ht.1, h.le⟩ x z hx)
    · have h' : a < t := by linarith [not_lt.1 h]
      rw [hR t ⟨h'.le, ht.2⟩]
      exact rst_local_vel (S' := Icc a (a + T₁)) (Φ' := fun r => Θ (r - a))
        (fun r hr => by rw [hR r hr]) ⟨h'.le, ht.2⟩
        (mem_nhdsWithin.2 ⟨Ioi a, isOpen_Ioi, h', fun r hr => ⟨hr.1.le, hr.2.2⟩⟩)
        (hvel t ⟨h'.le, ht.2⟩ x z hx)

/-! ### Transport: the reversed flow is a graph, up to reparametrization -/

/-- The slice of a jointly smooth scalar family has smooth chart readings. -/
theorem rst_contDiffOn_chartRep {f : M → ℝ} (hf : ContMDiff 𝓘(ℝ, E) 𝓘(ℝ, ℝ) ∞ f) (β : M) :
    ContDiffOn ℝ ∞ (chartRep E f β) (chart E β).target :=
  contMDiffOn_iff_contDiffOn.mp (hf.comp_contMDiffOn (contMDiffOn_extChartAt_symm (I := 𝓘(ℝ, E)) β))

/-- **Transport for the reversed flow.** For an `IsFlowS (-1)` `Γ`, for a short time `δ` there are a jointly
smooth family `w` of functions with `w₀ = 0` and diffeomorphisms `χ'_s` with
`graphImm (Γ 0) (w_s) ∘ χ'_s = Γ_s`. -/
theorem rst_transport [CompleteSpace E] {J : V →L[ℝ] V} (hJi : ∀ u v, ⟪J u, J v⟫ = ⟪u, v⟫)
    (hJJ : ∀ v, J (J v) = -v) (hdim : Module.finrank ℝ V = 2 * Module.finrank ℝ E)
    {T : ℝ} (hT : 0 < T) {Γ : ℝ → M → V} (hΓ : IsFlowS E J (-1) T Γ) :
    ∃ δ > 0, δ ≤ T ∧ ∃ w : ℝ → M → ℝ,
      ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, ℝ) ∞ (fun q : ℝ × M => w q.1 q.2) (Icc 0 δ ×ˢ univ) ∧
      (∀ x, w 0 x = 0) ∧
      ∃ χ' : ℝ → (M ≃ₘ⟮𝓘(ℝ, E), 𝓘(ℝ, E)⟯ M), ∀ s ∈ Icc (0 : ℝ) δ, ∀ x : M,
        graphImm E J (Γ 0) (w s) (χ' s x) = Γ s x := by
  have hF : IsLagImmersion E J (Γ 0) := hΓ.imm 0 ⟨le_rfl, hT.le⟩
  obtain ⟨δ, hδ, hδT, ψ, hψ, hψ0, hwsm, -, hcov⟩ :=
    Bridge.exists_gauge_graphS (κ := -1) hJi hJJ hF hdim hT hΓ rfl
  obtain ⟨w, hw⟩ : ∃ w : ℝ → M → ℝ, ∀ t x, w t x = -1 * gaugePot E J Γ ψ t x :=
    ⟨fun t x => -1 * gaugePot E J Γ ψ t x, fun _ _ => rfl⟩
  have hwsm' : ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, ℝ) ∞ (fun q : ℝ × M => w q.1 q.2)
      (Icc 0 δ ×ˢ univ) := hwsm.congr (fun q _ => hw q.1 q.2)
  have hwt : ∀ t, w t = fun x => -1 * gaugePot E J Γ ψ t x := fun t => funext (hw t)
  have hcov' : ∀ t ∈ Icc (0 : ℝ) δ, ∀ x : M, (∃ α, Covers E J (Γ 0) (w t) α x) ∧
      Γ t (ψ t x) = graphImm E J (Γ 0) (w t) x := by
    intro t ht x
    rw [hwt t]
    exact hcov t ht x
  obtain ⟨G', hG'⟩ : ∃ G' : ℝ → M → V, ∀ t, G' t = graphImm E J (Γ 0) (w t) := ⟨_, fun _ => rfl⟩
  have hG'eq : ∀ t ∈ Icc (0 : ℝ) δ, ∀ x : M, Γ t (ψ t x) = G' t x := by
    intro t ht x
    rw [hG' t]
    exact (hcov' t ht x).2
  have himm' : ∀ t ∈ Icc (0 : ℝ) δ, IsLagImmersion E J (G' t) := by
    intro t ht
    rw [hG' t]
    exact isLagImmersion_graphImm (rst_contDiffOn_chartRep (contMDiff_slice_fn hwsm' ht))
      hJi hJJ hF (fun x => (hcov' t ht x).1)
  have hG'sm : ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, V) ∞ (fun q : ℝ × M => G' q.1 q.2)
      (Icc 0 δ ×ˢ univ) := by
    have h1 : ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) ∞
        (fun q : ℝ × M => ((q.1, ψ q.1 q.2) : ℝ × M)) (Icc 0 δ ×ˢ univ) :=
      contMDiffOn_fst.prodMk hψ
    refine (ContMDiffOn.comp (g := fun q : ℝ × M => Γ q.1 q.2) hΓ.smooth h1 ?_).congr ?_
    · exact fun q hq => ⟨⟨hq.1.1, hq.1.2.trans hδT⟩, mem_univ _⟩
    · intro q hq
      exact (hG'eq q.1 hq.1 q.2).symm
  have hsmF : ∀ α : M, ContDiffOn ℝ ∞ (flowRep E G' α) (Icc 0 δ ×ˢ (chart E α).target) :=
    fun α => contDiffOn_flowRep hG'sm α
  obtain ⟨Y, hYsm, hYeq⟩ := DifferentialGeometry.Analysis.ODE.seeley_time_extend (I := 𝓘(ℝ, E))
    (tanFieldN E G' δ) δ hδ (contMDiffOn_tanFieldN hδ hsmF himm')
  obtain ⟨χ', h0, -, -, hode⟩ := exists_reparam_flow Y hYsm hδ.le
  have hχ : ∀ t ∈ Icc (0 : ℝ) δ, ∀ x, HasMFDerivAt 𝓘(ℝ, ℝ) 𝓘(ℝ, E) (fun s => χ' s x) t
      ((1 : ℝ →L[ℝ] ℝ).smulRight (tanFieldN E G' δ t (χ' t x))) := by
    intro t ht x
    have h := hode t (Icc_subset_Ioo_reparam ht) x
    rw [hYeq t ht] at h
    exact h
  refine ⟨δ, hδ, hδT, w, hwsm', fun x => by rw [hw, gaugePot_zero, mul_zero], χ', fun s hs x => ?_⟩
  have hψχ := psi_comp_chiN hJi hJJ hδ hsmF himm' (hΓ.mono hδT) hψ hψ0 hG'eq
    (χ := fun t x => χ' t x) h0 hχ hs x
  rw [← hG' s, ← hG'eq s hs, hψχ]

/-! ### The scalar side: restart from the datum `w_s` -/

/-- **Restart on the scalar side.** If `w` is jointly smooth on `[0, δ] × M` with `w₀ = 0` over the Lagrangian
immersion `F`, then for some `s ∈ (0, δ]` and `T₁ > s` there is a flow of Lagrangian immersions on `[0, T₁]`
starting from `graphImm F (w_s)`. -/
theorem rst_scalar [CompleteSpace E] {J : V →L[ℝ] V} (hJi : ∀ u v, ⟪J u, J v⟫ = ⟪u, v⟫)
    (hJJ : ∀ v, J (J v) = -v) (hdim : Module.finrank ℝ V = 2 * Module.finrank ℝ E) {F : M → V}
    (hF : IsLagImmersion E J F) {δ : ℝ} (hδ : 0 < δ) {w : ℝ → M → ℝ}
    (hwsm : ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, ℝ) ∞ (fun q : ℝ × M => w q.1 q.2)
      (Icc 0 δ ×ˢ univ))
    (hw0 : ∀ x, w 0 x = 0) :
    ∃ T₁ s : ℝ, 0 < s ∧ s ≤ δ ∧ s < T₁ ∧ ∃ Θ : ℝ → M → V, IsFlow E J T₁ Θ ∧
      Θ 0 = graphImm E J F (w s) := by
  obtain ⟨g, -⟩ := exists_inducedMetric hF
  obtain ⟨T₁, hT₁, u0, hu0⟩ := Bridge.exists_scalar_solution g hJi hJJ hF
  -- the data `δ / (k + 2)`
  have hseq0 : Tendsto (fun k : ℕ => δ / ((k : ℝ) + 2)) atTop (𝓝 0) :=
    tendsto_const_nhds.div_atTop (tendsto_atTop_add_const_right _ _ tendsto_natCast_atTop_atTop)
  have hpos : ∀ k : ℕ, 0 < δ / ((k : ℝ) + 2) := fun k => by positivity
  have hseq : Tendsto (fun k : ℕ => δ / ((k : ℝ) + 2)) atTop (𝓝[>] 0) :=
    tendsto_nhdsWithin_iff.2 ⟨hseq0, Eventually.of_forall hpos⟩
  have hle : ∀ k : ℕ, δ / ((k : ℝ) + 2) ≤ δ := fun k => by
    rw [div_le_iff₀ (by positivity)]
    nlinarith [(Nat.cast_nonneg k : (0 : ℝ) ≤ k)]
  have hconv : ∀ j : ℕ, ∀ ε > 0, ∀ᶠ k : ℕ in atTop, ∀ x : M,
      fiberNorm g 0 (0 + j) x ((iteratedCovGrad g 0 0 j
        (scalarFam g w (δ / ((k : ℝ) + 2)) - (0 : SmoothCcTensor g 0 0))).toSection x) ≤ ε := by
    intro j ε hε
    have h := tendsto_scalarFam_zero g hδ hwsm hw0 j ε hε
    filter_upwards [hseq.eventually h] with k hk
    simpa only [sub_zero] using hk
  obtain ⟨wv, hev, -⟩ := (ShortTime.mantegazza_martinazzi (cwData g hJi hJJ hF) (by norm_num) g).2.2
    0 T₁ u0 hu0 (fun k : ℕ => scalarFam g w (δ / ((k : ℝ) + 2))) hconv
  obtain ⟨k, hk1, hk2⟩ := (hev.and (hseq0.eventually (gt_mem_nhds hT₁))).exists
  have hs0 := hpos k
  have hsδ := hle k
  have hsl : ContMDiff 𝓘(ℝ, E) 𝓘(ℝ, ℝ) ∞ (w (δ / ((k : ℝ) + 2))) :=
    contMDiff_slice_fn hwsm ⟨hs0.le, hsδ⟩
  obtain ⟨χ'', h0'', hflow⟩ := exists_flow_of_solution g hJi hJJ hF hdim hk1
  refine ⟨T₁, δ / ((k : ℝ) + 2), hs0, hsδ, hk2, _, hflow, ?_⟩
  funext x
  show graphImm E J F (fnOf g (wv k 0)) (χ'' 0 x) = _
  rw [h0'', hk1.initial, fnOf_scalarFam g hsl]

/-! ### The restart theorem -/

/-- **Restart.** A flow of Lagrangian immersions of a compact manifold on `[0, T]` extends to a flow on
`[0, T + ε]` for some `ε > 0`. -/
theorem isFlow_restart {J : V →L[ℝ] V} (hJi : ∀ u v, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v, J (J v) = -v)
    (hdim : Module.finrank ℝ V = 2 * Module.finrank ℝ E) {T : ℝ} (hT : 0 < T) {Φ : ℝ → M → V}
    (hΦ : IsFlow E J T Φ) :
    ∃ ε > 0, ∃ Ψ : ℝ → M → V, IsFlow E J (T + ε) Ψ ∧ ∀ t ∈ Icc (0 : ℝ) T, Ψ t = Φ t := by
  haveI : CompleteSpace E := FiniteDimensional.complete ℝ E
  have hΓ : IsFlowS E J (-1) T (fun s => Φ (T - s)) := rst_reverse hΦ
  obtain ⟨δ, hδ, hδT, w, hwsm, hw0, χ', hχ'⟩ := rst_transport hJi hJJ hdim hT hΓ
  have hF : IsLagImmersion E J ((fun s => Φ (T - s)) 0) := hΓ.imm 0 ⟨le_rfl, hT.le⟩
  obtain ⟨T₁, s, hs, hsδ, hsT₁, Θ'', hΘ'', hΘ''0⟩ := rst_scalar hJi hJJ hdim hF hδ hwsm hw0
  have hsT : s ≤ T := hsδ.trans hδT
  -- the flow from `Φ (T - s)`
  have hΘ : IsFlow E J T₁ (fun t x => Θ'' t (χ' s x)) := hΘ''.comp_diffeo (χ' s)
  have hΘ0 : ∀ x, Θ'' 0 (χ' s x) = Φ (T - s) x := by
    intro x
    rw [hΘ''0]
    exact hχ' s ⟨hs.le, hsδ⟩ x
  -- uniqueness on `[0, s]`
  have hsh : IsFlow E J s (fun t => Φ (T - s + t)) := by
    have h := hΦ.shift (a := T - s) (by linarith)
    rwa [show T - (T - s) = s by ring] at h
  have hΘs : IsFlow E J s (fun t x => Θ'' t (χ' s x)) := hΘ.mono hsT₁.le
  have hun := flow_unique_of_immersion hJi hJJ hdim hΘs hsh (by
    funext x
    show Θ'' 0 (χ' s x) = Φ (T - s + 0) x
    rw [add_zero]
    exact hΘ0 x)
  have hΦ' : IsFlow E J (T - s + s) Φ := by rwa [sub_add_cancel]
  obtain ⟨Ψ, hΨ, hΨΦ⟩ := rst_glue (a := T - s) (by linarith) hs hsT₁.le hΦ' hΘ
    (fun t ht => hun t ht)
  refine ⟨T₁ - s, by linarith, Ψ, ?_, ?_⟩
  · rw [show T + (T₁ - s) = T - s + T₁ by ring]
    exact hΨ
  · rw [sub_add_cancel] at hΨΦ
    exact hΨΦ

/-- info: 'CWFlow.isFlow_restart' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms isFlow_restart

end CWFlow
