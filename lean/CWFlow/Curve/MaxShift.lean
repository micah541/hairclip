import CWFlow.Curve.Criterion
import CWFlow.Flow.Restart

/-!
# The time-shifted maximal flow; endpoint gluing (review W7, item S3b)

* `cdf_maximal_shift`: if `Φ` is the maximal curve diffusion flow from `Φ 0` with maximal time `Tmax` (the
  clauses of `cdf_exists_maximal`, T1) and `0 ≤ t₁`, `ofReal t₁ < Tmax`, then `t ↦ Φ (t₁ + t)` is the maximal
  flow from `Φ t₁`, with maximal time `Tmax - ofReal t₁`: it is an `IsFlow` on every `[0, T]` with
  `ofReal T < Tmax - ofReal t₁`; every `IsFlow` from `Φ t₁` on `[0, T]` has `ofReal T < Tmax - ofReal t₁` and
  equals `Φ (t₁ + ·)` there; and `∫ κ² ds` is unbounded on `[t₁, Tmax)` if `Tmax < ∞`.
  (C-1 Theorem (v).)  Proof: `IsFlow.shift_gen`, `rst_glue` (overlap form), `cdf_unique`, and T1 applied
  to `Φ t₁`.
* `isFlow_glue_endpoint`: flows on `[0, a]` and `[0, b]` matching at the seam (`Θ 0 = Φ a`) glue to a flow on
  `[0, a + b]` (U-1 Lemma "gluing", endpoint form).  Proof: `isFlow_restart` gives an overlap, then
  `rst_glue`, then uniqueness for the right piece.
-/

open Set Function Filter
open scoped Manifold Topology ContDiff ENNReal

noncomputable section

namespace CWFlow.Curve

open CWFlow

/-- A real time strictly between `t₁` and `Tmax`. -/
theorem exists_real_btwn_of_ofReal_lt {t₁ : ℝ} {Tmax : ℝ≥0∞} (h : ENNReal.ofReal t₁ < Tmax) :
    ∃ r : ℝ, t₁ < r ∧ ENNReal.ofReal r < Tmax := by
  obtain ⟨r, hr0, h1, h2⟩ := ENNReal.lt_iff_exists_real_btwn.1 h
  exact ⟨r, (ENNReal.ofReal_lt_ofReal_iff'.1 h1).1, h2⟩

/-- `ofReal T < Tmax - ofReal t₁ ↔ ofReal (t₁ + T) < Tmax` for `0 ≤ t₁`, `0 ≤ T`. -/
theorem ofReal_lt_sub_iff {t₁ T : ℝ} {Tmax : ℝ≥0∞} (ht₁ : 0 ≤ t₁) (hT : 0 ≤ T) :
    ENNReal.ofReal T < Tmax - ENNReal.ofReal t₁ ↔ ENNReal.ofReal (t₁ + T) < Tmax := by
  rw [lt_tsub_iff_left, ENNReal.ofReal_add ht₁ hT]

/-- **S3b (i).** The time-shifted maximal flow is the maximal flow from `Φ t₁`, with maximal time
`Tmax - ofReal t₁`. Hypotheses `h1`, `h2` are the clauses of T1 (`cdf_exists_maximal`) for `F₀ = Φ 0`. -/
theorem cdf_maximal_shift {Tmax : ℝ≥0∞} {Φ : ℝ → Circle → ℂ}
    (h1 : ∀ T : ℝ, 0 < T → ENNReal.ofReal T < Tmax → IsFlow E1 Jc T Φ)
    (h2 : ∀ (T : ℝ) (Ψ : ℝ → Circle → ℂ), 0 < T → IsFlow E1 Jc T Ψ → Ψ 0 = Φ 0 →
      ENNReal.ofReal T < Tmax ∧ ∀ t ∈ Icc (0 : ℝ) T, Ψ t = Φ t)
    {t₁ : ℝ} (ht₁ : 0 ≤ t₁) (ht₁T : ENNReal.ofReal t₁ < Tmax) :
    0 < Tmax - ENNReal.ofReal t₁ ∧
      (∀ T : ℝ, 0 < T → ENNReal.ofReal T < Tmax - ENNReal.ofReal t₁ →
        IsFlow E1 Jc T (fun t => Φ (t₁ + t))) ∧
      (∀ (T : ℝ) (Ψ : ℝ → Circle → ℂ), 0 < T → IsFlow E1 Jc T Ψ → Ψ 0 = Φ t₁ →
        ENNReal.ofReal T < Tmax - ENNReal.ofReal t₁ ∧ ∀ t ∈ Icc (0 : ℝ) T, Ψ t = Φ (t₁ + t)) ∧
      (Tmax - ENNReal.ofReal t₁ ≠ ⊤ → ∀ K : ℝ, ∃ t ∈ Ico (0 : ℝ) (Tmax - ENNReal.ofReal t₁).toReal,
        K < bendingEnergy (lift (Φ (t₁ + t)))) := by
  -- (a) the shifted map is a flow on `[0, T]` whenever `ofReal (t₁ + T) < Tmax`
  have hA : ∀ T : ℝ, 0 < T → ENNReal.ofReal (t₁ + T) < Tmax →
      IsFlow E1 Jc T (fun t => Φ (t₁ + t)) := by
    intro T hT hlt
    have h := IsFlow.shift_gen (h1 (t₁ + T) (by linarith) hlt) ht₁
    have e1 : t₁ + T - t₁ = T := by ring
    have e2 : (fun t => Φ (t + t₁)) = fun t => Φ (t₁ + t) := funext fun t => by rw [add_comm]
    rw [e1, e2] at h
    exact h
  -- (b) every flow from `Φ t₁` lives inside `[t₁, Tmax)` and is the shifted map
  have hB : ∀ (T : ℝ) (Ψ : ℝ → Circle → ℂ), 0 < T → IsFlow E1 Jc T Ψ → Ψ 0 = Φ t₁ →
      ENNReal.ofReal (t₁ + T) < Tmax ∧ ∀ t ∈ Icc (0 : ℝ) T, Ψ t = Φ (t₁ + t) := by
    intro T Ψ hT hΨ hΨ0
    obtain ⟨r, hr, hrT⟩ := exists_real_btwn_of_ofReal_lt ht₁T
    set s := min (r - t₁) T with hs_def
    have hs0 : 0 < s := lt_min (by linarith) hT
    have hsT : s ≤ T := min_le_right _ _
    have hsr : t₁ + s ≤ r := by linarith [min_le_left (r - t₁) T]
    have hΦr : IsFlow E1 Jc r Φ := h1 r (by linarith) hrT
    have hΦs : IsFlow E1 Jc (t₁ + s) Φ := hΦr.mono hsr
    have hsh : IsFlow E1 Jc s (fun t => Φ (t₁ + t)) := by
      have h := IsFlow.shift_gen hΦs ht₁
      have e1 : t₁ + s - t₁ = s := by ring
      have e2 : (fun t => Φ (t + t₁)) = fun t => Φ (t₁ + t) := funext fun t => by rw [add_comm]
      rw [e1, e2] at h
      exact h
    have hag : ∀ t ∈ Icc (0 : ℝ) s, Ψ t = Φ (t₁ + t) :=
      cdf_unique (hΨ.mono hsT) hsh (by show Ψ 0 = Φ (t₁ + 0); rw [hΨ0, add_zero])
    obtain ⟨Ψ', hΨ', hΨ'eq⟩ := rst_glue ht₁ hs0 hsT hΦs hΨ hag
    have hΨ'0 : Ψ' 0 = Φ 0 := hΨ'eq 0 ⟨le_rfl, by linarith⟩
    have hlt : ENNReal.ofReal (t₁ + T) < Tmax := (h2 (t₁ + T) Ψ' (by linarith) hΨ' hΨ'0).1
    exact ⟨hlt, cdf_unique hΨ (hA T hT hlt) (by show Ψ 0 = Φ (t₁ + 0); rw [hΨ0, add_zero])⟩
  -- positivity of the remaining time
  have hpos : 0 < Tmax - ENNReal.ofReal t₁ := tsub_pos_of_lt ht₁T
  -- (c) blow-up: T1 for `Φ t₁`, and its maximal time is `Tmax - ofReal t₁`
  have himm : IsLagImmersion E1 Jc (Φ t₁) := by
    obtain ⟨r, hr, hrT⟩ := exists_real_btwn_of_ofReal_lt ht₁T
    exact (h1 r (by linarith) hrT).imm t₁ ⟨ht₁, hr.le⟩
  obtain ⟨Tm', hTm'pos, Φ', hΦ'0, hΦ'1, hΦ'2, hΦ'3⟩ := cdf_exists_maximal himm
  -- `ofReal T < Tm' ↔ ofReal T < Tmax - ofReal t₁` for `T > 0`
  have hiff : ∀ T : ℝ, 0 < T → (ENNReal.ofReal T < Tm' ↔ ENNReal.ofReal T < Tmax - ENNReal.ofReal t₁) := by
    intro T hT
    rw [ofReal_lt_sub_iff ht₁ hT.le]
    constructor
    · intro h
      exact (hB T Φ' hT (hΦ'1 T hT h) hΦ'0).1
    · intro h
      exact (hΦ'2 T (fun t => Φ (t₁ + t)) hT (hA T hT h) (by show Φ (t₁ + 0) = Φ t₁; rw [add_zero])).1
  have heqT : Tm' = Tmax - ENNReal.ofReal t₁ := by
    refine le_antisymm (not_lt.1 fun hc => ?_) (not_lt.1 fun hc => ?_)
    · obtain ⟨r, -, hr1, hr2⟩ := ENNReal.lt_iff_exists_real_btwn.1 hc
      have hrpos : 0 < r := ENNReal.ofReal_pos.1 (lt_of_le_of_lt (by simp) hr1)
      exact absurd ((hiff r hrpos).1 hr2) (not_lt.2 hr1.le)
    · obtain ⟨r, -, hr1, hr2⟩ := ENNReal.lt_iff_exists_real_btwn.1 hc
      have hrpos : 0 < r := ENNReal.ofReal_pos.1 (lt_of_le_of_lt (by simp) hr1)
      exact absurd ((hiff r hrpos).2 hr2) (not_lt.2 hr1.le)
  refine ⟨hpos, fun T hT h => hA T hT ((ofReal_lt_sub_iff ht₁ hT.le).1 h),
    fun T Ψ hT hΨ hΨ0 => ?_, fun hfin K => ?_⟩
  · obtain ⟨hlt, heq⟩ := hB T Ψ hT hΨ hΨ0
    exact ⟨(ofReal_lt_sub_iff ht₁ hT.le).2 hlt, heq⟩
  · rw [← heqT] at hfin ⊢
    obtain ⟨t, ht, hK⟩ := hΦ'3 hfin K
    refine ⟨t, ht, ?_⟩
    -- `Φ' t = Φ (t₁ + t)`
    have hTpos : 0 < Tm'.toReal := ENNReal.toReal_pos hTm'pos.ne' hfin
    set T := (t + Tm'.toReal) / 2 with hT_def
    have hT0 : 0 < T := by linarith [ht.1, ht.2]
    have hTlt : ENNReal.ofReal T < Tm' := by
      rw [← ENNReal.ofReal_toReal hfin]
      exact (ENNReal.ofReal_lt_ofReal_iff hTpos).2 (by linarith [ht.2])
    have h := (hB T Φ' hT0 (hΦ'1 T hT0 hTlt) hΦ'0).2 t ⟨ht.1, by linarith [ht.2]⟩
    rw [← h]
    exact hK

/-- **S3b (ii). Endpoint gluing** (U-1 Lemma "gluing"): flows on `[0, a]` and `[0, b]` with
`Θ 0 = Φ a` glue to a flow `Ψ` on `[0, a + b]`, `Ψ = Φ` on `[0, a]` and `Ψ (a + t) = Θ t` on `[0, b]`. -/
theorem isFlow_glue_endpoint {a b : ℝ} (ha : 0 < a) (hb : 0 < b) {Φ Θ : ℝ → Circle → ℂ}
    (hΦ : IsFlow E1 Jc a Φ) (hΘ : IsFlow E1 Jc b Θ) (hseam : Θ 0 = Φ a) :
    ∃ Ψ : ℝ → Circle → ℂ, IsFlow E1 Jc (a + b) Ψ ∧ (∀ t ∈ Icc (0 : ℝ) a, Ψ t = Φ t) ∧
      ∀ t ∈ Icc (0 : ℝ) b, Ψ (a + t) = Θ t := by
  -- extend `Φ` a little past `a`
  obtain ⟨ε, hε, Φe, hΦe, hΦeq⟩ := isFlow_restart Jc_isom Jc_JJ hdim1 ha hΦ
  set s := min ε b with hs_def
  have hs0 : 0 < s := lt_min hε hb
  have hsb : s ≤ b := min_le_right _ _
  have hΦes : IsFlow E1 Jc (a + s) Φe := hΦe.mono (by linarith [min_le_left ε b])
  -- `Θ = Φe (a + ·)` on `[0, s]` by uniqueness
  have hsh : IsFlow E1 Jc s (fun t => Φe (a + t)) := by
    have h := IsFlow.shift_gen hΦes ha.le
    have e1 : a + s - a = s := by ring
    have e2 : (fun t => Φe (t + a)) = fun t => Φe (a + t) := funext fun t => by rw [add_comm]
    rw [e1, e2] at h
    exact h
  have hag : ∀ t ∈ Icc (0 : ℝ) s, Θ t = Φe (a + t) := by
    refine cdf_unique (hΘ.mono hsb) hsh ?_
    rw [hseam, add_zero, hΦeq a ⟨ha.le, le_rfl⟩]
  obtain ⟨Ψ, hΨ, hΨeq⟩ := rst_glue ha.le hs0 hsb hΦes hΘ hag
  have hleft : ∀ t ∈ Icc (0 : ℝ) a, Ψ t = Φ t := fun t ht => by
    rw [hΨeq t ⟨ht.1, by linarith [ht.2]⟩, hΦeq t ht]
  refine ⟨Ψ, hΨ, hleft, ?_⟩
  -- right piece: `Ψ (a + ·)` and `Θ` are flows on `[0, b]` from `Φ a`
  have hshΨ : IsFlow E1 Jc b (fun t => Ψ (a + t)) := by
    have h := IsFlow.shift_gen hΨ ha.le
    have e1 : a + b - a = b := by ring
    have e2 : (fun t => Ψ (t + a)) = fun t => Ψ (a + t) := funext fun t => by rw [add_comm]
    rw [e1, e2] at h
    exact h
  refine cdf_unique hshΨ hΘ ?_
  show Ψ (a + 0) = Θ 0
  rw [add_zero, hseam, hleft a ⟨ha.le, le_rfl⟩]

/-! ## Axiom audit -/

/-- info: 'CWFlow.Curve.cdf_maximal_shift' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms cdf_maximal_shift
/-- info: 'CWFlow.Curve.isFlow_glue_endpoint' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms isFlow_glue_endpoint

end CWFlow.Curve
