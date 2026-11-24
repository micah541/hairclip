import CWFlow.Extension.SymbolTime
import Mathlib.Analysis.Calculus.ParametricIntegral
import Mathlib.MeasureTheory.Measure.Haar.InnerProductSpace

/-!
# Smoothness of slices of a time-dependent function, and differentiation under the integral

For `Φ` jointly smooth on `T ×ˢ univ` (`T ⊆ ℝ` open):
* the slices `Φ (t, ·)` and `dt Φ (t, ·)` are smooth (`contDiff_slice`, `contDiff_dt_slice`);
* the jets of the slices, and a symbol evaluated on them, are jointly smooth in `(t, x)`
  (`contDiffOn_jp_slice`, `contDiffOn_symbol_slice`), on an open domain (`isOpen_symbol_slice_dom`);
* the time derivative passes under an integral against a compactly supported weight
  (`hasDerivAt_integral_mul_slice`).
-/

open scoped ContDiff
open CWFlow.Interp

namespace CWFlow.Extension

variable {n N : ℕ}

/-! ### Part 1: slices -/

theorem contDiff_slice {Φ : ℝ × EuclideanSpace ℝ (Fin n) → ℝ} {T : Set ℝ} (hT : IsOpen T)
    (hΦ : ContDiffOn ℝ ∞ Φ (T ×ˢ (Set.univ : Set (EuclideanSpace ℝ (Fin n))))) {t : ℝ}
    (ht : t ∈ T) : ContDiff ℝ ∞ (fun y => Φ (t, y)) := by
  have h : ContDiffOn ℝ ∞ (fun y : EuclideanSpace ℝ (Fin n) => Φ (t, y)) Set.univ :=
    hΦ.comp (contDiff_const.prodMk contDiff_id).contDiffOn
      (fun y _ => ⟨ht, Set.mem_univ _⟩)
  exact contDiffOn_univ.1 h

theorem contDiff_dt_slice {Φ : ℝ × EuclideanSpace ℝ (Fin n) → ℝ} {T : Set ℝ} (hT : IsOpen T)
    (hΦ : ContDiffOn ℝ ∞ Φ (T ×ˢ (Set.univ : Set (EuclideanSpace ℝ (Fin n))))) {t : ℝ}
    (ht : t ∈ T) : ContDiff ℝ ∞ (fun y => dt Φ (t, y)) := by
  have hO : IsOpen (T ×ˢ (Set.univ : Set (EuclideanSpace ℝ (Fin n)))) := hT.prod isOpen_univ
  have h : ContDiffOn ℝ ∞ (fun y : EuclideanSpace ℝ (Fin n) => dt Φ (t, y)) Set.univ :=
    (contDiffOn_dt hO hΦ).comp (contDiff_const.prodMk contDiff_id).contDiffOn
      (fun y _ => ⟨ht, Set.mem_univ _⟩)
  exact contDiffOn_univ.1 h

/-! ### Part 2: jets and symbols along slices -/

theorem contDiffOn_jp_slice {Φ : ℝ × EuclideanSpace ℝ (Fin n) → ℝ} {T : Set ℝ} (hT : IsOpen T)
    (hΦ : ContDiffOn ℝ ∞ Φ (T ×ˢ (Set.univ : Set (EuclideanSpace ℝ (Fin n))))) :
    ContDiffOn ℝ ∞
      (fun q : ℝ × EuclideanSpace ℝ (Fin n) => jp N (fun y => Φ (q.1, y)) q.2)
      (T ×ˢ Set.univ) := by
  have hO : IsOpen (T ×ˢ (Set.univ : Set (EuclideanSpace ℝ (Fin n)))) := hT.prod isOpen_univ
  refine (contDiff_snd.contDiffOn).prodMk ?_
  refine contDiffOn_pi.2 fun m => ?_
  refine contDiffOn_piLp' (𝕜 := ℝ) 2 (fun I => ?_)
  set e : Fin (m : ℕ) → ℝ × EuclideanSpace ℝ (Fin n) :=
    fun k => (0, EuclideanSpace.single (I k) 1) with he
  have hDm : ContDiffOn ℝ ∞ (iteratedFDeriv ℝ (m : ℕ) Φ) (T ×ˢ Set.univ) := by
    intro q hq
    exact ((hΦ.contDiffAt (hO.mem_nhds hq)).iteratedFDeriv_right
      (infty_add_natCast_le (m : ℕ))).contDiffWithinAt
  have hev : ContDiffOn ℝ ∞ (fun q => iteratedFDeriv ℝ (m : ℕ) Φ q e) (T ×ˢ Set.univ) :=
    (ContinuousMultilinearMap.apply ℝ (fun _ : Fin (m : ℕ) => ℝ × EuclideanSpace ℝ (Fin n)) ℝ
      e).contDiff.comp_contDiffOn hDm
  refine hev.congr fun q hq => ?_
  show jet (m : ℕ) (fun y => Φ (q.1, y)) q.2 I = _
  rw [jet_apply, iteratedFDeriv_slice hO hΦ (m : ℕ) (show (q.1, q.2) ∈ _ from hq)]

theorem contDiffOn_symbol_slice {W : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W]
    {A : Set (P n N)} {w ℓ : ℕ} {𝔭 : P n N → W} (h : IsSymbol A w ℓ 𝔭)
    {Φ : ℝ × EuclideanSpace ℝ (Fin n) → ℝ} {T : Set ℝ} (hT : IsOpen T)
    (hΦ : ContDiffOn ℝ ∞ Φ (T ×ˢ (Set.univ : Set (EuclideanSpace ℝ (Fin n))))) :
    ContDiffOn ℝ ∞
      (fun q : ℝ × EuclideanSpace ℝ (Fin n) => 𝔭 (jp N (fun y => Φ (q.1, y)) q.2))
      {q | q ∈ T ×ˢ (Set.univ : Set (EuclideanSpace ℝ (Fin n))) ∧
        jp N (fun y => Φ (q.1, y)) q.2 ∈ A} :=
  h.smooth.comp ((contDiffOn_jp_slice hT hΦ).mono (fun q hq => hq.1)) (fun q hq => hq.2)

theorem isOpen_symbol_slice_dom {A : Set (P n N)} (hA : IsOpen A)
    {Φ : ℝ × EuclideanSpace ℝ (Fin n) → ℝ} {T : Set ℝ} (hT : IsOpen T)
    (hΦ : ContDiffOn ℝ ∞ Φ (T ×ˢ (Set.univ : Set (EuclideanSpace ℝ (Fin n))))) :
    IsOpen {q : ℝ × EuclideanSpace ℝ (Fin n) |
      q ∈ T ×ˢ (Set.univ : Set (EuclideanSpace ℝ (Fin n))) ∧
        jp N (fun y => Φ (q.1, y)) q.2 ∈ A} := by
  have hO : IsOpen (T ×ˢ (Set.univ : Set (EuclideanSpace ℝ (Fin n)))) := hT.prod isOpen_univ
  have := (contDiffOn_jp_slice (N := N) hT hΦ).continuousOn.isOpen_inter_preimage hO hA
  simpa [Set.inter_def] using this

/-! ### Part 3: differentiation under the integral -/

theorem slc_continuous_mul {G : ℝ × EuclideanSpace ℝ (Fin n) → ℝ}
    {O : Set (ℝ × EuclideanSpace ℝ (Fin n))} (hO : IsOpen O) (hG : ContinuousOn G O)
    {ρ : EuclideanSpace ℝ (Fin n) → ℝ} (hρ : Continuous ρ) {s : ℝ}
    (hs : ∀ x ∈ tsupport ρ, (s, x) ∈ O) : Continuous fun x => ρ x * G (s, x) := by
  rw [continuous_iff_continuousAt]
  intro x
  by_cases hx : x ∈ tsupport ρ
  · have h1 : ContinuousAt G (s, x) := hG.continuousAt (hO.mem_nhds (hs x hx))
    have h2 : ContinuousAt (fun y : EuclideanSpace ℝ (Fin n) => (s, y)) x :=
      (continuous_const.prodMk continuous_id).continuousAt
    exact hρ.continuousAt.mul (h1.comp h2)
  · have hev : (fun y => ρ y * G (s, y)) =ᶠ[nhds x] fun _ => (0 : ℝ) := by
      filter_upwards [(isClosed_tsupport ρ).isOpen_compl.mem_nhds hx] with y hy
      rw [image_eq_zero_of_notMem_tsupport hy, zero_mul]
    exact continuousAt_const.congr hev.symm

theorem hasDerivAt_integral_mul_slice {F : ℝ × EuclideanSpace ℝ (Fin n) → ℝ}
    {O : Set (ℝ × EuclideanSpace ℝ (Fin n))} (hO : IsOpen O) (hF : ContDiffOn ℝ ∞ F O)
    {ρ : EuclideanSpace ℝ (Fin n) → ℝ} (hρ : Continuous ρ) (hρc : HasCompactSupport ρ)
    {t : ℝ} (ht : ∀ x ∈ tsupport ρ, (t, x) ∈ O) :
    HasDerivAt (fun s : ℝ => ∫ x, ρ x * F (s, x))
      (∫ x, ρ x * fderiv ℝ F (t, x) (1, 0)) t := by
  set S : Set (EuclideanSpace ℝ (Fin n)) := tsupport ρ with hS
  have hSc : IsCompact S := hρc
  -- tube lemma
  obtain ⟨u, v, hu, hv, htu, hSv, huv⟩ := generalized_tube_lemma (isCompact_singleton (x := t))
    hSc hO (fun p hp => by
      obtain ⟨h1, h2⟩ := hp
      rw [Set.mem_singleton_iff] at h1
      have := ht p.2 h2
      rw [← h1] at this
      exact this)
  obtain ⟨ε, hε, hεu⟩ := Metric.isOpen_iff.1 hu t (htu rfl)
  set δ : ℝ := ε / 2 with hδ
  have hδ0 : 0 < δ := by positivity
  have hsub : ∀ s ∈ Metric.closedBall t δ, ∀ x ∈ S, (s, x) ∈ O := by
    intro s hs x hx
    refine huv ⟨hεu ?_, hSv hx⟩
    rw [Metric.mem_ball]
    rw [Metric.mem_closedBall] at hs
    linarith
  -- bound
  have hD : ContDiffOn ℝ ∞ (fun q => fderiv ℝ F q (1, 0)) O :=
    (hF.fderiv_of_isOpen hO (m := ∞) (by simp)).clm_apply contDiffOn_const
  have hDc : ContinuousOn (fun q => fderiv ℝ F q (1, 0)) O := hD.continuousOn
  have hK : IsCompact (Metric.closedBall t δ ×ˢ S) := (isCompact_closedBall t δ).prod hSc
  obtain ⟨M, hM⟩ := hK.exists_bound_of_continuousOn
    (hDc.mono (fun q hq => hsub q.1 hq.1 q.2 hq.2))
  have hcont : ∀ s ∈ Metric.ball t δ, Continuous fun x => ρ x * F (s, x) := fun s hs =>
    slc_continuous_mul hO hF.continuousOn hρ
      (fun x hx => hsub s (Metric.ball_subset_closedBall hs) x hx)
  have hcont' : ∀ s ∈ Metric.ball t δ,
      Continuous fun x => ρ x * fderiv ℝ F (s, x) (1, 0) := fun s hs =>
    slc_continuous_mul hO hDc hρ (fun x hx => hsub s (Metric.ball_subset_closedBall hs) x hx)
  have htb : t ∈ Metric.ball t δ := Metric.mem_ball_self hδ0
  have hcs : ∀ (g : EuclideanSpace ℝ (Fin n) → ℝ), True →
      HasCompactSupport (fun x => ρ x * g x) := fun g _ => hρc.mul_right
  refine (hasDerivAt_integral_of_dominated_loc_of_deriv_le
    (F := fun s x => ρ x * F (s, x)) (F' := fun s x => ρ x * fderiv ℝ F (s, x) (1, 0))
    (bound := fun x => ‖ρ x‖ * M) (s := Metric.ball t δ) (x₀ := t)
    (Metric.ball_mem_nhds t hδ0) ?_ ?_ ?_ ?_ ?_ ?_).2
  · filter_upwards [Metric.ball_mem_nhds t hδ0] with s hs
    exact (hcont s hs).aestronglyMeasurable
  · exact (hcont t htb).integrable_of_hasCompactSupport (hcs (fun x => F (t, x)) trivial)
  · exact (hcont' t htb).aestronglyMeasurable
  · refine Filter.Eventually.of_forall fun a s hs => ?_
    by_cases ha : a ∈ S
    · rw [norm_mul]
      exact mul_le_mul_of_nonneg_left
        (hM (s, a) ⟨Metric.ball_subset_closedBall hs, ha⟩) (norm_nonneg _)
    · rw [image_eq_zero_of_notMem_tsupport ha]
      simp
  · exact ((hρ.norm.integrable_of_hasCompactSupport hρc.norm).mul_const M)
  · refine Filter.Eventually.of_forall fun a s hs => ?_
    by_cases ha : a ∈ S
    · have hsO : (s, a) ∈ O := hsub s (Metric.ball_subset_closedBall hs) a ha
      have hd : DifferentiableAt ℝ F (s, a) :=
        (hF.contDiffAt (hO.mem_nhds hsO)).differentiableAt (by simp)
      have hc : HasDerivAt (fun r : ℝ => (r, a)) ((1 : ℝ), (0 : EuclideanSpace ℝ (Fin n))) s :=
        (hasDerivAt_id s).prodMk (hasDerivAt_const s a)
      exact (hd.hasFDerivAt.comp_hasDerivAt s hc).const_mul (ρ a)
    · have h0 := image_eq_zero_of_notMem_tsupport ha
      have : (fun r : ℝ => ρ a * F (r, a)) = fun _ => 0 := by
        funext r; rw [h0, zero_mul]
      rw [this, h0, zero_mul]
      exact hasDerivAt_const s 0

/-- info: 'CWFlow.Extension.contDiff_slice' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms contDiff_slice

/-- info: 'CWFlow.Extension.contDiffOn_symbol_slice' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms contDiffOn_symbol_slice

/-- info: 'CWFlow.Extension.hasDerivAt_integral_mul_slice' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms hasDerivAt_integral_mul_slice

end CWFlow.Extension
