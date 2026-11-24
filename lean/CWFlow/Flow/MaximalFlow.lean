import CWFlow.Flow.ExtensionTheorem
import CWFlow.Flow.Main

/-!
# The maximal flow (Chen–Warren, Cor. "maximal flow")

For a Lagrangian immersion `F₀` of a closed manifold into a Hermitian vector space there is a maximal existence
time `Tmax ∈ (0, ∞]` and a family `Φ` which is an `IsFlow` on `[0, T]` for every `T < Tmax`, such that any flow
from `F₀` on `[0, T]` has `T < Tmax` and agrees with `Φ`. If `Tmax < ∞`, the second fundamental form is
unbounded on `[0, Tmax)`.

`Tmax` is the supremum of the times of existence of flows from `F₀`. It is not attained, because a flow on `[0, T]`
restarts (`isFlow_restart`); flows are glued through uniqueness (`flow_unique_of_immersion`); and a bound on `|A|²`
would allow an extension past `Tmax` (`flow_extends_of_bounded_A`).
-/

open Set
open scoped ContDiff Manifold Topology RealInnerProductSpace ENNReal

set_option linter.unusedSectionVars false

namespace CWFlow

universe u

variable {E : Type u} [NormedAddCommGroup E] [NormedSpace ℝ E] [FiniteDimensional ℝ E]
  [NeZero (Module.finrank ℝ E)]
variable {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
  [CompactSpace M] [T2Space M] [SigmaCompactSpace M]
variable {V : Type u} [NormedAddCommGroup V] [InnerProductSpace ℝ V] [FiniteDimensional ℝ V]

/-- Being a flow on `[0, T]` depends only on the values on `[0, T]`. -/
theorem mf_isFlow_congr {J : V →L[ℝ] V} {T : ℝ} {Φ Ψ : ℝ → M → V} (h : IsFlow E J T Ψ)
    (heq : ∀ t ∈ Icc (0 : ℝ) T, Φ t = Ψ t) : IsFlow E J T Φ where
  smooth := h.smooth.congr (fun q hq => by
    show Φ q.1 q.2 = Ψ q.1 q.2
    rw [heq q.1 hq.1])
  imm := fun t ht => by rw [heq t ht]; exact h.imm t ht
  vel := fun t ht x z hx => by
    have e : chartRep E (Φ t) z = chartRep E (Ψ t) z := by rw [heq t ht]
    rw [e]
    exact (h.vel t ht x z hx).congr (fun s hs => by rw [heq s hs]) (by rw [heq t ht])

/-- **The maximal flow** (Chen–Warren, Cor. "maximal flow", in `ℂⁿ`). -/
theorem exists_maximal_flow {n : ℕ} (ι : EuclideanSpace ℝ (Fin n) ≃L[ℝ] E) {J : V →L[ℝ] V}
    (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    (hdim : Module.finrank ℝ V = 2 * n) {F₀ : M → V} (hF₀ : IsLagImmersion E J F₀) :
    ∃ Tmax : ℝ≥0∞, 0 < Tmax ∧ ∃ Φ : ℝ → M → V, Φ 0 = F₀ ∧
      (∀ T : ℝ, 0 < T → ENNReal.ofReal T < Tmax → IsFlow E J T Φ) ∧
      (∀ (T : ℝ) (Ψ : ℝ → M → V), 0 < T → IsFlow E J T Ψ → Ψ 0 = F₀ →
        ENNReal.ofReal T < Tmax ∧ ∀ t ∈ Icc (0 : ℝ) T, Ψ t = Φ t) ∧
      (Tmax ≠ ⊤ → ∀ K : ℝ, ¬ ∀ t ∈ Ico (0 : ℝ) Tmax.toReal, ∀ z : M, ∀ y ∈ Extension.chartDomE ι z,
        Extension.aNormSq J (Φ t ∘ Extension.chartInv ι z) y ≤ K ^ 2) := by
  classical
  have hdim' : Module.finrank ℝ V = 2 * Module.finrank ℝ E := by
    rw [hdim, ← ι.toLinearEquiv.finrank_eq, finrank_euclideanSpace_fin]
  obtain ⟨T₀, hT₀, Φ₀, h00, hs, hi, hv⟩ := exists_flow_of_immersion hJi hJJ hF₀ hdim'
  set S : Set ℝ := {T | 0 < T ∧ ∃ Ψ : ℝ → M → V, IsFlow E J T Ψ ∧ Ψ 0 = F₀} with hSdef
  have hS0 : T₀ ∈ S := ⟨hT₀, Φ₀, ⟨hs, hi, hv⟩, h00⟩
  have hext : ∀ T ∈ S, ∃ T' ∈ S, T < T' := by
    rintro T ⟨hT, Ψ, hΨ, h0⟩
    obtain ⟨ε, hε, Ψ', hΨ', hagree⟩ := isFlow_restart hJi hJJ hdim' hT hΨ
    exact ⟨T + ε, ⟨by linarith, Ψ', hΨ', by rw [hagree 0 ⟨le_rfl, hT.le⟩, h0]⟩, by linarith⟩
  have hch : ∀ T : ℝ, T ∈ S → ∃ Ψ : ℝ → M → V, IsFlow E J T Ψ ∧ Ψ 0 = F₀ := fun T hT => hT.2
  choose! fl hfl using hch
  have uniq : ∀ T ∈ S, ∀ T' ∈ S, ∀ t : ℝ, 0 ≤ t → t ≤ T → t ≤ T' → fl T t = fl T' t := by
    intro T hT T' hT' t ht0 htT htT'
    exact flow_unique_of_immersion hJi hJJ hdim' ((hfl T hT).1.mono (min_le_left T T'))
      ((hfl T' hT').1.mono (min_le_right T T')) (by rw [(hfl T hT).2, (hfl T' hT').2]) t
      ⟨ht0, le_min htT htT'⟩
  obtain ⟨Φ, hΦ0, hΦeq⟩ : ∃ Φ : ℝ → M → V, Φ 0 = F₀ ∧ ∀ T ∈ S, ∀ t ∈ Icc (0 : ℝ) T, Φ t = fl T t := by
    refine ⟨fun t => if h : ∃ T, T ∈ S ∧ t ≤ T then fl h.choose t else F₀, ?_, ?_⟩
    · have h : ∃ T, T ∈ S ∧ (0 : ℝ) ≤ T := ⟨T₀, hS0, hT₀.le⟩
      dsimp only
      rw [dif_pos h]
      exact (hfl _ (Classical.choose_spec h).1).2
    · intro T hT t ht
      have h : ∃ T, T ∈ S ∧ t ≤ T := ⟨T, hT, ht.2⟩
      dsimp only
      rw [dif_pos h]
      exact uniq _ (Classical.choose_spec h).1 T hT t ht.1 (Classical.choose_spec h).2 ht.2
  set Tmax : ℝ≥0∞ := sSup (ENNReal.ofReal '' S) with hTm
  have hpos : 0 < Tmax :=
    lt_of_lt_of_le (ENNReal.ofReal_pos.2 hT₀) (le_sSup ⟨T₀, hS0, rfl⟩)
  have hlt : ∀ T ∈ S, ENNReal.ofReal T < Tmax := by
    intro T hT
    obtain ⟨T', hT', hTT'⟩ := hext T hT
    exact lt_of_lt_of_le ((ENNReal.ofReal_lt_ofReal_iff hT'.1).2 hTT') (le_sSup ⟨T', hT', rfl⟩)
  have hfl1 : ∀ T : ℝ, 0 < T → ENNReal.ofReal T < Tmax → IsFlow E J T Φ := by
    intro T hT hTlt
    obtain ⟨_, ⟨T', hT', rfl⟩, hlt'⟩ := lt_sSup_iff.1 hTlt
    have hTT' : T < T' := (ENNReal.ofReal_lt_ofReal_iff hT'.1).1 hlt'
    exact mf_isFlow_congr ((hfl T' hT').1.mono hTT'.le)
      (fun t ht => hΦeq T' hT' t ⟨ht.1, ht.2.trans hTT'.le⟩)
  refine ⟨Tmax, hpos, Φ, hΦ0, hfl1, ?_, ?_⟩
  · intro T Ψ hT hΨ h0
    have hTS : T ∈ S := ⟨hT, Ψ, hΨ, h0⟩
    refine ⟨hlt T hTS, fun t ht => ?_⟩
    exact (flow_unique_of_immersion hJi hJJ hdim' hΨ (hfl T hTS).1
      (h0.trans (hfl T hTS).2.symm) t ht).trans (hΦeq T hTS t ht).symm
  · intro hfin K hK
    have heq : Tmax = ENNReal.ofReal Tmax.toReal := (ENNReal.ofReal_toReal hfin).symm
    have hTpos : 0 < Tmax.toReal := ENNReal.toReal_pos hpos.ne' hfin
    obtain ⟨ε, hε, Ψ, hΨ, hagree⟩ := flow_extends_of_bounded_A ι hJi hJJ hdim hTpos (Φ := Φ)
      (fun T' hT' => hfl1 T' hT'.1 (by
        rw [heq]; exact (ENNReal.ofReal_lt_ofReal_iff hTpos).2 hT'.2)) hK
    have hmem : Tmax.toReal + ε ∈ S :=
      ⟨by linarith, Ψ, hΨ, by rw [hagree 0 ⟨le_rfl, hTpos⟩, hΦ0]⟩
    have h1 := hlt _ hmem
    have h2 : Tmax ≤ ENNReal.ofReal (Tmax.toReal + ε) := by
      exact heq.le.trans (ENNReal.ofReal_le_ofReal (by linarith))
    exact absurd h1 (not_lt.2 h2)

/-- info: 'CWFlow.exists_maximal_flow' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_maximal_flow

end CWFlow
