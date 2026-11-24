import CWFlow.Curve.Lift
import CWFlow.Flow.MaximalFlow
import CWFlow.Flow.Restart

/-!
# The `L²` blow-up criterion for curve diffusion flow (SCOPE W6)

* `cdf_extends_of_bounded_energy` (T2): a curve diffusion flow on `[0, T)` with `∫ κ² ds` bounded on `[t₀, T)`
  extends past `T`;
* `cdf_exists_maximal` (T1): the maximal flow from an immersed closed curve; if `Tmax < ∞` then
  `∫ κ² ds` is unbounded on `[0, Tmax)`.

Route (SCOPE §2.2 (A)): (T4) the lift is `IsCDF`; (T5) `κ²` is bounded on `[t₀, T)`; (T3) so `|A|²` is bounded
in every chart for the shifted flow `Φ(t₀ + ·)`; `flow_extends_of_bounded_A` extends it; `rst_glue` glues it to
`Φ`, and uniqueness identifies the glued flow with `Φ` on `[0, T)`.
-/

open Set Function Filter
open scoped Manifold Topology ContDiff RealInnerProductSpace ENNReal

noncomputable section

namespace CWFlow.Curve

open CWFlow

/-- **(T2) Extension form of the criterion** (replaces DKS (D1) as used in C-1). -/
theorem cdf_extends_of_bounded_energy {T : ℝ} (hT : 0 < T) {Φ : ℝ → Circle → ℂ}
    (hΦ : ∀ T' ∈ Ioo 0 T, IsFlow E1 Jc T' Φ) {K t₀ : ℝ} (ht₀ : t₀ ∈ Ioo 0 T)
    (hK : ∀ t ∈ Ico t₀ T, bendingEnergy (lift (Φ t)) ≤ K) :
    ∃ ε > 0, ∃ Ψ : ℝ → Circle → ℂ, IsFlow E1 Jc (T + ε) Ψ ∧ ∀ t ∈ Ico 0 T, Ψ t = Φ t := by
  -- (T4) + (T5): `κ²` is bounded on `[t₀, T)`
  have hcdf := isCDF_of_forall_isFlow hΦ
  obtain ⟨C, hC⟩ := hcdf.curv_sq_le ht₀ (K := K) (fun t ht => hK t ht)
  -- the shifted flow `Φ(· + t₀)`
  set Φs : ℝ → Circle → ℂ := fun t => Φ (t + t₀) with hΦs
  have hΦs_flow : ∀ T' ∈ Ioo 0 (T - t₀), IsFlow E1 Jc T' Φs := by
    intro T' hT'
    have h := IsFlow.shift_gen (hΦ (T' + t₀) ⟨by linarith [hT'.1, ht₀.1], by linarith [hT'.2]⟩)
      ht₀.1.le
    have e : T' + t₀ - t₀ = T' := by ring
    rw [e] at h
    exact h
  have hTt : 0 < T - t₀ := by linarith [ht₀.2]
  -- (T3): `|A|² ≤ max C 0` in every chart, for the shifted flow
  have hA : ∀ t ∈ Ico 0 (T - t₀), ∀ z : Circle, ∀ y ∈ Extension.chartDomE ι1 z,
      Extension.aNormSq Jc (Φs t ∘ Extension.chartInv ι1 z) y ≤ (√(max C 0)) ^ 2 := by
    intro t ht z y hy
    obtain ⟨T', hT', htT'⟩ : ∃ T' ∈ Ioo 0 (T - t₀), t ≤ T' :=
      ⟨(t + (T - t₀)) / 2, ⟨by linarith [ht.1], by linarith [ht.2]⟩, by linarith [ht.2]⟩
    have himm := (hΦs_flow T' hT').imm t ⟨ht.1, htT'⟩
    obtain ⟨θ, -, he⟩ := aNormSq_eq_curv_sq himm z hy
    rw [he, Real.sq_sqrt (le_max_right _ _)]
    have := hC (t + t₀) ⟨by linarith [ht.1], by linarith [ht.2]⟩ θ
    exact this.trans (le_max_left _ _)
  -- extend the shifted flow
  obtain ⟨ε, hε, Θ, hΘ, hΘeq⟩ := flow_extends_of_bounded_A ι1 Jc_isom Jc_JJ
    Complex.finrank_real_complex hTt hΦs_flow hA
  -- glue `Θ` to `Φ` at `t₀` with overlap `s`
  set s := (T - t₀) / 2 with hs
  have hs0 : 0 < s := by positivity
  have hsT : s ≤ T - t₀ + ε := by linarith
  have hΦa : IsFlow E1 Jc (t₀ + s) Φ := hΦ _ ⟨by linarith [ht₀.1], by linarith⟩
  have hag : ∀ t ∈ Icc (0 : ℝ) s, Θ t = Φ (t₀ + t) := by
    intro t ht
    rw [hΘeq t ⟨ht.1, by linarith [ht.2]⟩]
    show Φ (t + t₀) = Φ (t₀ + t)
    rw [add_comm]
  obtain ⟨Ψ, hΨ, hΨeq⟩ := rst_glue ht₀.1.le hs0 hsT hΦa hΘ hag
  have hTe : t₀ + (T - t₀ + ε) = T + ε := by ring
  rw [hTe] at hΨ
  refine ⟨ε, hε, Ψ, hΨ, fun t ht => ?_⟩
  -- uniqueness on `[0, T'']`, `t ≤ T'' < T`
  set T'' := (t + T) / 2 with hT''
  have hT''T : T'' < T := by linarith [ht.2]
  have hT''0 : 0 < T'' := by linarith [ht.1, ht.2]
  have h0 : Ψ 0 = Φ 0 := hΨeq 0 ⟨le_rfl, by linarith [ht₀.1]⟩
  exact cdf_unique (hΨ.mono (by linarith)) (hΦ T'' ⟨hT''0, hT''T⟩) h0 t
    ⟨ht.1, by linarith [ht.2]⟩

/-- **(T1) The `L²` blow-up criterion** (replaces DKS Thm 3.1 as used in Lemma "blow-up criterion, the form
needed"): the maximal curve diffusion flow from an immersed closed curve; if `Tmax < ∞` then `∫ κ² ds` is
unbounded on `[0, Tmax)`. -/
theorem cdf_exists_maximal {F₀ : Circle → ℂ} (hF₀ : IsLagImmersion E1 Jc F₀) :
    ∃ Tmax : ℝ≥0∞, 0 < Tmax ∧ ∃ Φ : ℝ → Circle → ℂ, Φ 0 = F₀ ∧
      (∀ T : ℝ, 0 < T → ENNReal.ofReal T < Tmax → IsFlow E1 Jc T Φ) ∧
      (∀ (T : ℝ) (Ψ : ℝ → Circle → ℂ), 0 < T → IsFlow E1 Jc T Ψ → Ψ 0 = F₀ →
        ENNReal.ofReal T < Tmax ∧ ∀ t ∈ Icc (0 : ℝ) T, Ψ t = Φ t) ∧
      (Tmax ≠ ⊤ → ∀ K : ℝ, ∃ t ∈ Ico (0 : ℝ) Tmax.toReal, K < bendingEnergy (lift (Φ t))) := by
  obtain ⟨Tmax, hpos, Φ, h0, h1, h2, -⟩ := cdf_exists_maximal_A hF₀
  refine ⟨Tmax, hpos, Φ, h0, h1, h2, fun hfin K => ?_⟩
  by_contra hcon
  push_neg at hcon
  have heq : Tmax = ENNReal.ofReal Tmax.toReal := (ENNReal.ofReal_toReal hfin).symm
  have hTpos : 0 < Tmax.toReal := ENNReal.toReal_pos hpos.ne' hfin
  have hΦ : ∀ T' ∈ Ioo 0 Tmax.toReal, IsFlow E1 Jc T' Φ := fun T' hT' =>
    h1 T' hT'.1 (by rw [heq]; exact (ENNReal.ofReal_lt_ofReal_iff hTpos).2 hT'.2)
  obtain ⟨ε, hε, Ψ, hΨ, hΨeq⟩ := cdf_extends_of_bounded_energy hTpos hΦ (K := K)
    (t₀ := Tmax.toReal / 2) ⟨by positivity, by linarith⟩
    (fun t ht => hcon t ⟨by linarith [ht.1, hTpos], ht.2⟩)
  have hΨ0 : Ψ 0 = F₀ := by rw [hΨeq 0 ⟨le_rfl, hTpos⟩, h0]
  have hlt := (h2 (Tmax.toReal + ε) Ψ (by linarith) hΨ hΨ0).1
  have hle : Tmax ≤ ENNReal.ofReal (Tmax.toReal + ε) :=
    heq.le.trans (ENNReal.ofReal_le_ofReal (by linarith))
  exact absurd hlt (not_lt.2 hle)

/-! ## Axiom audit -/

/-- info: 'CWFlow.Curve.cdf_exists_maximal' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms cdf_exists_maximal
/-- info: 'CWFlow.Curve.cdf_extends_of_bounded_energy' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms cdf_extends_of_bounded_energy
/-- info: 'CWFlow.Curve.aNormSq_eq_curv_sq' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms aNormSq_eq_curv_sq
/-- info: 'CWFlow.Curve.IsFlow.lift_isCDF' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms IsFlow.lift_isCDF
/-- info: 'CWFlow.Curve.IsCDF.curv_sq_le' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms IsCDF.curv_sq_le

end CWFlow.Curve
