import CWFlow.Moser.Field
import CWFlow.Moser.Flow
import DifferentialGeometry.Analysis.Calculus.Cutoff.Basic

/-!
# Moser's relative theorem in a chart: the canonical time-one map

Let `λ, β` be 1-forms on a finite-dimensional normed space `W`, smooth on an open `O`. Put `Ω_t = dλ + t·dβ` and
`X_t = −(ι_· Ω_t)⁻¹ β` (`Moser/Field.lean`). The good set `moserGood` consists of the `(t, q)` with `q ∈ O` and
`ι_· Ω_t(q)` invertible. A *Moser curve* from `q` is a solution `γ : [0, 1] → W` of `γ' = X_t(γ)` with
`γ(0) = q`, whose graph lies in the good set.

* Moser curves are unique (`moserCurve_unique`).
* Their starting points form an open set `moserDom` (`isOpen_moserDom`).
* The time-one map `moserEnd q = γ(1)` is smooth there (`contDiffOn_moserEnd`).
* `moserEnd*(dλ + dβ) = dλ` (`moserEnd_symplectic`).
* A zero `z` of `β` whose constant curve is good lies in `moserDom` and is fixed (`mem_moserDom_of_zero`,
  `moserEnd_of_zero`).

The map `moserEnd` is defined by the uncut ODE and involves no choices. This is what makes it natural under
changes of chart, and hence what makes it glue (`Weinstein/Glue.lean`). The proofs work with cut-off flows: a
cut-off `f = χ X` with `χ = 1` on a compact neighbourhood of a compact part of the good set (`exists_cutoff`) has
a global smooth flow (`Moser/Flow.lean`). Moser curves near a given one are trajectories of that flow
(`exists_local_flow`), and `moser_flow` gives the Moser identity along them.
-/

open ContinuousAlternatingMap Filter Set Function

open scoped Topology ContDiff

noncomputable section

namespace CWFlow

open DifferentialGeometry.Analysis

variable {W : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W]

/-- The Lie derivative depends only on the germ of the form. -/
theorem lieAlt_congr {Ω₁ Ω₂ : W → W [⋀^Fin 2]→L[ℝ] ℝ} {X : W → W} {p : W} (h : Ω₁ =ᶠ[𝓝 p] Ω₂)
    (a b : W) : lieAlt Ω₁ X p a b = lieAlt Ω₂ X p a b := by
  rw [lieAlt, lieAlt, h.fderiv_eq, h.eq_of_nhds]

/-- The Moser field vanishes wherever `β` does. -/
theorem moserField_of_zero {lam β : W → W →L[ℝ] ℝ} {z : W} (hz : β z = 0) (t : ℝ) :
    moserField lam β t z = 0 := by
  have h : toAlt1 β z = 0 := by
    ext v
    simp [hz]
  rw [moserField, h, _root_.map_zero, neg_zero]

theorem Icc_mem_nhdsWithin_Ici_of_mem {t : ℝ} (ht : t ∈ Ico (0 : ℝ) 1) : Icc (0 : ℝ) 1 ∈ 𝓝[Ici t] t :=
  mem_nhdsWithin.2 ⟨Iio 1, isOpen_Iio, ht.2, fun _ hy => ⟨ht.1.trans hy.2, hy.1.le⟩⟩

variable [FiniteDimensional ℝ W]

/-! ### Moser curves -/

section Defs

variable (lam β : W → W →L[ℝ] ℝ) (O : Set W)

/-- The good set: `q ∈ O` and `ι_· Ω_t(q)` invertible. -/
def moserGood : Set (ℝ × W) :=
  {p | p.2 ∈ O ∧ ∃ e : W ≃L[ℝ] (W [⋀^Fin 1]→L[ℝ] ℝ), (e : _ →L[ℝ] _) = iotaT lam β p.1 p.2}

/-- A Moser curve from `q`: `γ(0) = q`, `γ' = X_t(γ)` on `[0, 1]`, with graph in the good set. -/
structure IsMoserCurve (q : W) (γ : ℝ → W) : Prop where
  start : γ 0 = q
  good : ∀ t ∈ Icc (0 : ℝ) 1, (t, γ t) ∈ moserGood lam β O
  deriv : ∀ t ∈ Icc (0 : ℝ) 1, HasDerivWithinAt γ (moserField lam β t (γ t)) (Icc 0 1) t

/-- The starting points of Moser curves. -/
def moserDom : Set W := {q | ∃ γ, IsMoserCurve lam β O q γ}

open Classical in
/-- The time-one map `q ↦ γ(1)` (and `q ↦ q` off `moserDom`). -/
def moserEnd (q : W) : W := if h : ∃ γ, IsMoserCurve lam β O q γ then h.choose 1 else q

end Defs

variable {lam β : W → W →L[ℝ] ℝ} {O : Set W}

theorem IsMoserCurve.continuousOn {q : W} {γ : ℝ → W} (h : IsMoserCurve lam β O q γ) :
    ContinuousOn γ (Icc 0 1) :=
  fun t ht => (h.deriv t ht).continuousWithinAt

theorem IsMoserCurve.mem {q : W} {γ : ℝ → W} (h : IsMoserCurve lam β O q γ) {t : ℝ}
    (ht : t ∈ Icc (0 : ℝ) 1) : γ t ∈ O :=
  (h.good t ht).1

theorem moserDom_subset : moserDom lam β O ⊆ O := fun _ ⟨_, h⟩ => h.start ▸ h.mem ⟨le_rfl, zero_le_one⟩

theorem IsMoserCurve.isCompact_graph {q : W} {γ : ℝ → W} (h : IsMoserCurve lam β O q γ) :
    IsCompact ((fun t => (t, γ t)) '' Icc (0 : ℝ) 1) :=
  isCompact_Icc.image_of_continuousOn (continuousOn_id.prodMk h.continuousOn)

theorem IsMoserCurve.graph_subset {q : W} {γ : ℝ → W} (h : IsMoserCurve lam β O q γ) :
    (fun t => (t, γ t)) '' Icc (0 : ℝ) 1 ⊆ moserGood lam β O :=
  image_subset_iff.2 fun t ht => h.good t ht

/-- Two integral curves on `[0, 1]` of a smooth compactly supported field, with the same start, agree. -/
theorem eqOn_of_cutoff {f : ℝ → W → W} (hf : ContDiff ℝ ∞ (uncurry f))
    (hc : HasCompactSupport (uncurry f)) {γ γ' : ℝ → W} (hγ : ContinuousOn γ (Icc 0 1))
    (hγd : ∀ t ∈ Ico (0 : ℝ) 1, HasDerivWithinAt γ (f t (γ t)) (Ici t) t)
    (hγ' : ContinuousOn γ' (Icc 0 1))
    (hγ'd : ∀ t ∈ Ico (0 : ℝ) 1, HasDerivWithinAt γ' (f t (γ' t)) (Ici t) t) (h0 : γ 0 = γ' 0) :
    EqOn γ γ' (Icc 0 1) := by
  obtain ⟨C, hC⟩ := lipschitzWith_slice hf hc
  exact ODE_solution_unique_of_mem_Icc_right (v := f) (s := fun _ => univ) (K := C)
    (fun t _ => (hC t).lipschitzOnWith) hγ hγd (fun _ _ => mem_univ _) hγ' hγ'd
    (fun _ _ => mem_univ _) h0

/-- A trajectory of a cut-off field that stays where the cut-off is `1` is a Moser curve. -/
theorem isMoserCurve_psi {f : ℝ → W → W} (hf : ContDiff ℝ ∞ (uncurry f))
    (hc : HasCompactSupport (uncurry f)) {L : Set (ℝ × W)} (hLG : L ⊆ moserGood lam β O)
    (hfX : ∀ t y, (t, y) ∈ L → f t y = moserField lam β t y) {q : W}
    (hq : ∀ s ∈ Icc (0 : ℝ) 1, (s, psi hf hc s q) ∈ L) :
    IsMoserCurve lam β O q (fun s => psi hf hc s q) where
  start := psi_zero hf hc q
  good s hs := hLG (hq s hs)
  deriv s hs := by
    rw [← hfX s _ (hq s hs)]
    exact (hasDerivAt_psi hf hc q s).hasDerivWithinAt

/-! ### Results needing smoothness of the data -/

section Smooth

variable (hO : IsOpen O) (hlam : ContDiffOn ℝ ∞ lam O) (hβ : ContDiffOn ℝ ∞ β O)
include hO hlam hβ

theorem isOpen_moserGood : IsOpen (moserGood lam β O) := by
  rw [isOpen_iff_mem_nhds]
  rintro ⟨t, q⟩ ⟨hq, hinv⟩
  filter_upwards [eventually_iotaT_invertible (hlam.contDiffAt (hO.mem_nhds hq))
    (hβ.contDiffAt (hO.mem_nhds hq)) hinv,
    continuous_snd.continuousAt.preimage_mem_nhds (hO.mem_nhds hq)] with p h1 h2
  exact ⟨h2, h1⟩

theorem contDiffOn_moserField :
    ContDiffOn ℝ ∞ (fun p : ℝ × W => moserField lam β p.1 p.2) (moserGood lam β O) := by
  rintro ⟨t, q⟩ ⟨hq, hinv⟩
  exact (contDiffAt_moserField (hlam.contDiffAt (hO.mem_nhds hq)) (hβ.contDiffAt (hO.mem_nhds hq))
    hinv).contDiffWithinAt

/-- A cut-off of the Moser field, equal to it on a compact neighbourhood `L` of a given compact part `C` of the
good set, and vanishing at the zeros of `β`. -/
theorem exists_cutoff {C : Set (ℝ × W)} (hC : IsCompact C) (hCG : C ⊆ moserGood lam β O) :
    ∃ (f : ℝ → W → W) (hf : ContDiff ℝ ∞ (uncurry f)) (hc : HasCompactSupport (uncurry f))
      (L : Set (ℝ × W)), C ⊆ interior L ∧ L ⊆ moserGood lam β O ∧
      (∀ t y, (t, y) ∈ L → f t y = moserField lam β t y) ∧ (∀ t z, β z = 0 → f t z = 0) := by
  have hG := isOpen_moserGood hO hlam hβ
  obtain ⟨L, hL, hCL, hLG⟩ := exists_compact_between hC hG hCG
  obtain ⟨M, hM, hLM, hMG⟩ := exists_compact_between hL hG hLG
  obtain ⟨χ, hχ, hχ1, hχM, -⟩ := exists_bump_one_on hL isOpen_interior hLM
  have hχc : HasCompactSupport χ :=
    hM.of_isClosed_subset (isClosed_tsupport χ) (hχM.trans interior_subset)
  refine ⟨fun t q => χ (t, q) • moserField lam β t q,
    contDiff_cutoff_smul hG hχ (hχM.trans (interior_subset.trans hMG))
      (contDiffOn_moserField hO hlam hβ),
    hχc.smul_right (f' := fun p : ℝ × W => moserField lam β p.1 p.2), L, hCL, hLG,
    fun t y hy => ?_, fun t z hz => ?_⟩
  · have h1 : χ (t, y) = 1 := hχ1 hy
    simp only [h1, one_smul]
  · simp only [moserField_of_zero hz, smul_zero]

/-- **Uniqueness of Moser curves.** -/
theorem moserCurve_unique {q : W} {γ γ' : ℝ → W} (h : IsMoserCurve lam β O q γ)
    (h' : IsMoserCurve lam β O q γ') : EqOn γ γ' (Icc 0 1) := by
  obtain ⟨f, hf, hc, L, hCL, -, hfX, -⟩ := exists_cutoff hO hlam hβ
    (h.isCompact_graph.union h'.isCompact_graph) (union_subset h.graph_subset h'.graph_subset)
  have key : ∀ {δ : ℝ → W}, IsMoserCurve lam β O q δ →
      (fun t => (t, δ t)) '' Icc (0 : ℝ) 1 ⊆ (fun t => (t, γ t)) '' Icc (0 : ℝ) 1 ∪
        (fun t => (t, γ' t)) '' Icc (0 : ℝ) 1 →
      ∀ t ∈ Ico (0 : ℝ) 1, HasDerivWithinAt δ (f t (δ t)) (Ici t) t := by
    intro δ hδ hsub t ht
    have hL : (t, δ t) ∈ L :=
      interior_subset (hCL (hsub ⟨t, Ico_subset_Icc_self ht, rfl⟩))
    rw [hfX t _ hL]
    exact (hδ.deriv t (Ico_subset_Icc_self ht)).mono_of_mem_nhdsWithin
      (Icc_mem_nhdsWithin_Ici_of_mem ht)
  exact eqOn_of_cutoff hf hc h.continuousOn (key h subset_union_left) h'.continuousOn
    (key h' subset_union_right) (h.start.trans h'.start.symm)

/-- Near a starting point of a Moser curve, the Moser curves are the trajectories of one cut-off flow, and they
stay in the interior of the region where the cut-off is `1`. -/
theorem exists_local_flow {q₀ : W} (h₀ : q₀ ∈ moserDom lam β O) :
    ∃ (f : ℝ → W → W) (hf : ContDiff ℝ ∞ (uncurry f)) (hc : HasCompactSupport (uncurry f))
      (L : Set (ℝ × W)) (V : Set W), L ⊆ moserGood lam β O ∧
      (∀ t y, (t, y) ∈ L → f t y = moserField lam β t y) ∧ IsOpen V ∧ q₀ ∈ V ∧
      ∀ q ∈ V, ∀ s ∈ Icc (0 : ℝ) 1, (s, psi hf hc s q) ∈ interior L := by
  obtain ⟨γ, hγ⟩ := h₀
  obtain ⟨f, hf, hc, L, hCL, hLG, hfX, -⟩ :=
    exists_cutoff hO hlam hβ hγ.isCompact_graph hγ.graph_subset
  have hγL : ∀ t ∈ Icc (0 : ℝ) 1, (t, γ t) ∈ interior L := fun t ht => hCL ⟨t, ht, rfl⟩
  have heq : EqOn (fun s => psi hf hc s q₀) γ (Icc 0 1) := by
    refine eqOn_of_cutoff hf hc (fun t _ => (hasDerivAt_psi hf hc q₀ t).continuousAt.continuousWithinAt)
      (fun t _ => (hasDerivAt_psi hf hc q₀ t).hasDerivWithinAt) hγ.continuousOn (fun t ht => ?_)
      (by rw [psi_zero, hγ.start])
    rw [hfX t _ (interior_subset (hγL t (Ico_subset_Icc_self ht)))]
    exact (hγ.deriv t (Ico_subset_Icc_self ht)).mono_of_mem_nhdsWithin
      (Icc_mem_nhdsWithin_Ici_of_mem ht)
  have hcont : Continuous (fun p : ℝ × W => (p.1, psi hf hc p.1 p.2)) :=
    continuous_fst.prodMk (contDiff_psi hf hc).continuous
  have hn : IsOpen {p : ℝ × W | (p.1, psi hf hc p.1 p.2) ∈ interior L} :=
    isOpen_interior.preimage hcont
  have hp : Icc (0 : ℝ) 1 ×ˢ {q₀} ⊆ {p : ℝ × W | (p.1, psi hf hc p.1 p.2) ∈ interior L} := by
    rintro ⟨s, q⟩ ⟨hs, hq⟩
    rw [mem_singleton_iff] at hq
    subst hq
    simp only [mem_setOf_eq]
    rw [show psi hf hc s q = γ s from heq hs]
    exact hγL s hs
  obtain ⟨u, V, -, hV, hsu, hqV, huV⟩ := generalized_tube_lemma isCompact_Icc isCompact_singleton hn hp
  exact ⟨f, hf, hc, L, V, hLG, hfX, hV, hqV rfl, fun q hq s hs =>
    huV (mk_mem_prod (hsu hs) hq)⟩

theorem isOpen_moserDom : IsOpen (moserDom lam β O) := by
  rw [isOpen_iff_mem_nhds]
  intro q₀ h₀
  obtain ⟨f, hf, hc, L, V, hLG, hfX, hV, hq₀, hVL⟩ := exists_local_flow hO hlam hβ h₀
  filter_upwards [hV.mem_nhds hq₀] with q hq
  exact ⟨_, isMoserCurve_psi hf hc hLG hfX fun s hs => interior_subset (hVL q hq s hs)⟩

theorem moserEnd_eq {q : W} {γ : ℝ → W} (h : IsMoserCurve lam β O q γ) : moserEnd lam β O q = γ 1 := by
  have hex : ∃ γ, IsMoserCurve lam β O q γ := ⟨γ, h⟩
  rw [moserEnd, dif_pos hex]
  exact moserCurve_unique hO hlam hβ hex.choose_spec h ⟨zero_le_one, le_rfl⟩

theorem moserEnd_mem {q : W} (hq : q ∈ moserDom lam β O) : moserEnd lam β O q ∈ O := by
  obtain ⟨γ, hγ⟩ := hq
  rw [moserEnd_eq hO hlam hβ hγ]
  exact hγ.mem ⟨zero_le_one, le_rfl⟩

/-- Near each point of `moserDom`, `moserEnd` is the time-one map of a cut-off flow. -/
theorem moserEnd_eventuallyEq {q₀ : W} (h₀ : q₀ ∈ moserDom lam β O) :
    ∃ (f : ℝ → W → W) (hf : ContDiff ℝ ∞ (uncurry f)) (hc : HasCompactSupport (uncurry f))
      (L : Set (ℝ × W)), L ⊆ moserGood lam β O ∧
      (∀ t y, (t, y) ∈ L → f t y = moserField lam β t y) ∧
      (∀ s ∈ Icc (0 : ℝ) 1, (s, psi hf hc s q₀) ∈ interior L) ∧
      moserEnd lam β O =ᶠ[𝓝 q₀] psi hf hc 1 := by
  obtain ⟨f, hf, hc, L, V, hLG, hfX, hV, hq₀, hVL⟩ := exists_local_flow hO hlam hβ h₀
  refine ⟨f, hf, hc, L, hLG, hfX, hVL q₀ hq₀, ?_⟩
  filter_upwards [hV.mem_nhds hq₀] with q hq
  exact moserEnd_eq hO hlam hβ
    (isMoserCurve_psi hf hc hLG hfX fun s hs => interior_subset (hVL q hq s hs))

theorem contDiffOn_moserEnd : ContDiffOn ℝ ∞ (moserEnd lam β O) (moserDom lam β O) := by
  intro q₀ h₀
  obtain ⟨f, hf, hc, -, -, -, -, heq⟩ := moserEnd_eventuallyEq hO hlam hβ h₀
  have hψ : ContDiff ℝ ∞ (psi hf hc 1) := (contDiff_psi hf hc).comp (contDiff_const.prodMk contDiff_id)
  exact (hψ.contDiffAt.congr_of_eventuallyEq heq).contDiffWithinAt

/-- The Moser identity for a cut-off flow, along a trajectory staying where the cut-off is `1`. -/
theorem moser_identity_cutoff {f : ℝ → W → W} (hf : ContDiff ℝ ∞ (uncurry f))
    (hc : HasCompactSupport (uncurry f)) {L : Set (ℝ × W)} (hLG : L ⊆ moserGood lam β O)
    (hfX : ∀ t y, (t, y) ∈ L → f t y = moserField lam β t y) {x : W}
    (hx : ∀ s ∈ Icc (0 : ℝ) 1, (s, psi hf hc s x) ∈ interior L) (v w : W) :
    omegaT lam β 1 (psi hf hc 1 x) ![fderiv ℝ (psi hf hc 1) x v, fderiv ℝ (psi hf hc 1) x w] =
      extDeriv (toAlt1 lam) x ![v, w] := by
  have hlamA : ∀ q ∈ O, ContDiffAt ℝ ∞ lam q := fun q hq => hlam.contDiffAt (hO.mem_nhds hq)
  have hβA : ∀ q ∈ O, ContDiffAt ℝ ∞ β q := fun q hq => hβ.contDiffAt (hO.mem_nhds hq)
  have hdiff : ∀ q ∈ O, DifferentiableAt ℝ lam q ∧ DifferentiableAt ℝ β q := fun q hq =>
    ⟨(hlamA q hq).differentiableAt (by simp), (hβA q hq).differentiableAt (by simp)⟩
  have hfs : ∀ s, ContDiff ℝ ∞ (f s) := fun s => hf.comp (contDiff_const.prodMk contDiff_id)
  have hslice : ∀ s y, (s, y) ∈ interior L → ∀ᶠ y' in 𝓝 y, (s, y') ∈ interior L := fun s y hy =>
    (continuous_const.prodMk continuous_id).continuousAt.preimage_mem_nhds (isOpen_interior.mem_nhds hy)
  have hmemO : ∀ s ∈ Icc (0 : ℝ) 1, psi hf hc s x ∈ O := fun s hs =>
    (hLG (interior_subset (hx s hs))).1
  have hcurve : ∀ s ∈ Icc (0 : ℝ) 1, HasDerivAt (fun r => omegaT lam β r (psi hf hc r x))
      (extDeriv (toAlt1 β) (psi hf hc s x) +
        fderiv ℝ (omegaT lam β s) (psi hf hc s x) (f s (psi hf hc s x))) s := fun s hs => by
    set y := psi hf hc s x with hy
    have hyO : y ∈ O := hmemO s hs
    have hA := ((contDiffAt_extDeriv_toAlt1 (hlamA y hyO)).differentiableAt (by simp)).hasFDerivAt
    have hB := ((contDiffAt_extDeriv_toAlt1 (hβA y hyO)).differentiableAt (by simp)).hasFDerivAt
    have hψ := hasDerivAt_psi hf hc x s
    have h1 := hA.comp_hasDerivAt s hψ
    have h2 := (hasDerivAt_id' s).smul (hB.comp_hasDerivAt s hψ)
    have hΩ : HasFDerivAt (omegaT lam β s) _ y := hA.add (hB.const_smul s)
    refine (h1.add h2).congr_deriv ?_
    rw [hΩ.fderiv]
    simp only [ContinuousLinearMap.add_apply, ContinuousLinearMap.coe_smul', Pi.smul_apply,
      Function.comp_apply, one_smul]
    abel
  have hmoser : ∀ s ∈ Icc (0 : ℝ) 1, ∀ a b, extDeriv (toAlt1 β) (psi hf hc s x) ![a, b] +
      lieAlt (omegaT lam β s) (f s) (psi hf hc s x) a b = 0 := fun s hs a b => by
    set y := psi hf hc s x with hy
    have hyL : (s, y) ∈ interior L := hx s hs
    have hyO : y ∈ O := hmemO s hs
    have heqΩ : omegaT lam β s =ᶠ[𝓝 y] extDeriv (toAlt1 (fun q => lam q + s • β q)) := by
      filter_upwards [hO.mem_nhds hyO] with q hq
      exact omegaT_eq (hdiff q hq).1 (hdiff q hq).2 s
    rw [lieAlt_congr heqΩ]
    refine moser_condition_exact ?_ (differentiableAt_toAlt1 (hdiff y hyO).2)
      ((hfs s).differentiable (by simp) y) ?_ a b
    · exact contDiffAt_toAlt1 (((hlamA y hyO).add ((hβA y hyO).const_smul s)).of_le
        (WithTop.coe_le_coe.2 le_top))
    · filter_upwards [hslice s y hyL] with q hq
      intro b
      have hqG := hLG (interior_subset hq)
      rw [hfX s q (interior_subset hq), ← omegaT_eq (hdiff q hqG.1).1 (hdiff q hqG.1).2 s,
        moserField_spec hqG.2 b]
      simp
  have h := moser_flow hf hc (Ω := omegaT lam β) (Ωdot := fun _ => extDeriv (toAlt1 β)) x hcurve
    hmoser v w
  rw [h, omegaT, zero_smul, add_zero]

/-- **Moser's identity for the canonical time-one map:** `moserEnd*(dλ + dβ) = dλ` on `moserDom`. -/
theorem moserEnd_symplectic {q : W} (hq : q ∈ moserDom lam β O) (a b : W) :
    omegaT lam β 1 (moserEnd lam β O q)
        ![fderiv ℝ (moserEnd lam β O) q a, fderiv ℝ (moserEnd lam β O) q b] =
      extDeriv (toAlt1 lam) q ![a, b] := by
  obtain ⟨f, hf, hc, L, hLG, hfX, hqL, heq⟩ := moserEnd_eventuallyEq hO hlam hβ hq
  rw [heq.fderiv_eq, heq.eq_of_nhds]
  exact moser_identity_cutoff hO hlam hβ hf hc hLG hfX hqL a b

/-- A zero of `β` along which the constant curve is good is a starting point of a Moser curve. -/
theorem mem_moserDom_of_zero {z : W} (hz : β z = 0)
    (hgood : ∀ t ∈ Icc (0 : ℝ) 1, (t, z) ∈ moserGood lam β O) : z ∈ moserDom lam β O :=
  ⟨fun _ => z, rfl, hgood, fun t _ => by
    rw [moserField_of_zero hz]
    exact hasDerivWithinAt_const t _ z⟩

/-- Zeros of `β` are fixed by the time-one map. -/
theorem moserEnd_of_zero {z : W} (hz : β z = 0) (hzD : z ∈ moserDom lam β O) :
    moserEnd lam β O z = z := by
  obtain ⟨γ, hγ⟩ := hzD
  obtain ⟨f, hf, hc, L, hCL, -, hfX, hf0⟩ :=
    exists_cutoff hO hlam hβ hγ.isCompact_graph hγ.graph_subset
  have h := eqOn_of_cutoff hf hc hγ.continuousOn (fun t ht => ?_) continuousOn_const
    (fun t _ => ?_) hγ.start (⟨zero_le_one, le_rfl⟩ : (1 : ℝ) ∈ Icc (0 : ℝ) 1)
  · rw [moserEnd_eq hO hlam hβ hγ]
    exact h
  · rw [hfX t _ (interior_subset (hCL ⟨t, Ico_subset_Icc_self ht, rfl⟩))]
    exact (hγ.deriv t (Ico_subset_Icc_self ht)).mono_of_mem_nhdsWithin
      (Icc_mem_nhdsWithin_Ici_of_mem ht)
  · rw [hf0 t z hz]
    exact hasDerivWithinAt_const t _ z

end Smooth

end CWFlow
