import CWFlow.Interp.Jets
import Mathlib.Analysis.Calculus.FDeriv.Symmetric
import Mathlib.Analysis.Calculus.FDeriv.Prod
import Mathlib.Analysis.Calculus.Deriv.Prod

/-!
# Spatial jets commute with the time derivative

For `Φ` smooth on an open set `Ω ⊆ ℝ × ℝⁿ` the `m`-jet in `x` of the slice `Φ (s, ·)` is differentiable in
`s`, with derivative the `m`-jet of the slice of `dt Φ = ∂ₜΦ` (`hasDerivAt_jet_slice`). The key inputs are
the slice formula `iteratedFDeriv_slice` and the symmetry `iteratedFDeriv_cons_eq_snoc` of higher
derivatives of a smooth (not nec. analytic) function.
-/

set_option linter.unusedSectionVars false

open scoped ContDiff
open CWFlow.Interp

namespace CWFlow.Extension

section Symm

variable {F : Type*} [NormedAddCommGroup F] [NormedSpace ℝ F]
  {W : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W]

theorem iteratedFDeriv_succ_succ_eq {f : F → W} {q : F} (k : ℕ)
    (a b : F) (w : Fin k → F) :
    iteratedFDeriv ℝ (k + 2) f q (Fin.cons a (Fin.cons b w)) =
      fderiv ℝ (fderiv ℝ (iteratedFDeriv ℝ k f)) q a b w := by
  have h : fderiv ℝ (fderiv ℝ (iteratedFDeriv ℝ k f)) q = ((continuousMultilinearCurryLeftEquiv ℝ
      (fun _ : Fin (k + 1) => F) W : _ ≃ₗᵢ[ℝ] _) : _ →L[ℝ] _).comp
        (fderiv ℝ (iteratedFDeriv ℝ (k + 1) f) q) := by
    rw [fderiv_iteratedFDeriv]
    exact ContinuousLinearEquiv.comp_fderiv (𝕜 := ℝ) (E := ContinuousMultilinearMap ℝ (fun _ : Fin (k + 1) => F) W) (F := F →L[ℝ] ContinuousMultilinearMap ℝ (fun _ : Fin k => F) W) (iso := (continuousMultilinearCurryLeftEquiv ℝ (fun _ : Fin (k + 1) => F) W).toContinuousLinearEquiv) (f := iteratedFDeriv ℝ (k + 1) f)
  rw [h, iteratedFDeriv_succ_apply_left]
  rfl

/-- Moving the first argument of the `(m+1)`-st derivative to the last position. -/
theorem iteratedFDeriv_cons_eq_snoc {U : Set F} (hU : IsOpen U) {f : F → W}
    (hf : ∀ y ∈ U, ContDiffAt ℝ ∞ f y) (m : ℕ) :
    ∀ q ∈ U, ∀ (v : F) (w : Fin m → F),
      iteratedFDeriv ℝ (m + 1) f q (Fin.cons v w) = iteratedFDeriv ℝ (m + 1) f q (Fin.snoc w v) := by
  induction m with
  | zero =>
    intro q _ v w
    congr 1
    funext i
    fin_cases i
    simp [Fin.snoc]
  | succ m ih =>
    intro q hq v w
    obtain ⟨a, w', rfl⟩ : ∃ a w', w = Fin.cons a w' := ⟨w 0, Fin.tail w, (Fin.cons_self_tail w).symm⟩
    have hq' : ContDiffAt ℝ ∞ (iteratedFDeriv ℝ m f) q :=
      (hf q hq).iteratedFDeriv_right (i := m) (infty_add_natCast_le _)
    have hs := hq'.isSymmSndFDerivAt (n := ∞) (by simp; exact WithTop.coe_le_coe.2 le_top)
    have hD : DifferentiableAt ℝ (iteratedFDeriv ℝ (m + 1) f) q :=
      ((hf q hq).iteratedFDeriv_right (i := m + 1) (infty_add_natCast_le _)).differentiableAt (by simp)
    have step1 : iteratedFDeriv ℝ (m + 1 + 1) f q (Fin.cons v (Fin.cons a w')) =
        iteratedFDeriv ℝ (m + 1 + 1) f q (Fin.cons a (Fin.cons v w')) := by
      rw [iteratedFDeriv_succ_succ_eq, iteratedFDeriv_succ_succ_eq, hs v a]
    have step2 : fderiv ℝ (iteratedFDeriv ℝ (m + 1) f) q a (Fin.cons v w') =
        fderiv ℝ (iteratedFDeriv ℝ (m + 1) f) q a (Fin.snoc w' v) := by
      have e1 : fderiv ℝ (fun y => iteratedFDeriv ℝ (m + 1) f y (Fin.cons v w')) q a =
          fderiv ℝ (iteratedFDeriv ℝ (m + 1) f) q a (Fin.cons v w') := by
        rw [fderiv_continuousMultilinear_apply_const hD]
        rfl
      have e2 : fderiv ℝ (fun y => iteratedFDeriv ℝ (m + 1) f y (Fin.snoc w' v)) q a =
          fderiv ℝ (iteratedFDeriv ℝ (m + 1) f) q a (Fin.snoc w' v) := by
        rw [fderiv_continuousMultilinear_apply_const hD]
        rfl
      rw [← e1, ← e2]
      have hev : (fun y => iteratedFDeriv ℝ (m + 1) f y (Fin.cons v w')) =ᶠ[nhds q]
          fun y => iteratedFDeriv ℝ (m + 1) f y (Fin.snoc w' v) := by
        filter_upwards [hU.mem_nhds hq] with y hy
        exact ih y hy v w'
      rw [hev.fderiv_eq]
    rw [step1, iteratedFDeriv_succ_apply_left]
    simp only [Fin.cons_zero, Fin.tail_cons]
    rw [step2, ← Fin.cons_snoc_eq_snoc_cons, iteratedFDeriv_succ_apply_left]
    simp

theorem natCast_le_infty (m : ℕ) : (m : WithTop ℕ∞) ≤ ∞ :=
  WithTop.coe_le_coe.2 le_top

end Symm

section Slice

variable {n : ℕ} {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V] [FiniteDimensional ℝ V]

local notation "E" => EuclideanSpace ℝ (Fin n)

/-- The `m`-th derivative of a slice `y ↦ Φ (s, y)` is the `m`-th derivative of `Φ` evaluated on
vectors tangent to the slice. -/
theorem iteratedFDeriv_slice {Φ : ℝ × E → V} {Ω : Set (ℝ × E)}
    (hΩ : IsOpen Ω) (hΦ : ContDiffOn ℝ ∞ Φ Ω) (m : ℕ) {s : ℝ} {x : E} (h : (s, x) ∈ Ω)
    (v : Fin m → E) :
    iteratedFDeriv ℝ m (fun y => Φ (s, y)) x v =
      iteratedFDeriv ℝ m Φ (s, x) (fun k => (0, v k)) := by
  set a : ℝ × E := (s, 0) with ha
  set Ψ : ℝ × E → V := fun p => Φ (a + p) with hΨdef
  set Ω' : Set (ℝ × E) := (fun p => a + p) ⁻¹' Ω with hΩ'def
  have hΩ' : IsOpen Ω' := hΩ.preimage (continuous_const.add continuous_id)
  have hΨ : ContDiffOn ℝ ∞ Ψ Ω' :=
    hΦ.comp (contDiff_const.add contDiff_id).contDiffOn (fun p hp => hp)
  set ι : E →L[ℝ] ℝ × E := ContinuousLinearMap.inr ℝ ℝ E with hι
  have hslice : (fun y => Φ (s, y)) = Ψ ∘ ι := by
    funext y; simp [Ψ, a, ι]
  have hιx : ι x ∈ Ω' := by simpa [Ω', a, ι] using h
  have h1 := ι.iteratedFDerivWithin_comp_right hΨ hΩ'.uniqueDiffOn
    (hΩ'.preimage ι.continuous).uniqueDiffOn hιx (natCast_le_infty m)
  rw [iteratedFDerivWithin_of_isOpen m (hΩ'.preimage ι.continuous) (by simpa using hιx),
    iteratedFDerivWithin_of_isOpen m hΩ' hιx] at h1
  rw [hslice, h1]
  have h2 : iteratedFDeriv ℝ m Ψ (ι x) = iteratedFDeriv ℝ m Φ (a + ι x) :=
    iteratedFDeriv_comp_add_left m a (ι x)
  have h3 : a + ι x = (s, x) := by simp [a, ι]
  rw [h2, h3]
  simp [ι]

/-- The time derivative of `Φ`, as a function of `(t, x)`. -/
noncomputable def dt (Φ : ℝ × E → V) (q : ℝ × E) : V :=
  fderiv ℝ Φ q (1, 0)

theorem contDiffOn_dt {Φ : ℝ × E → V} {Ω : Set (ℝ × E)}
    (hΩ : IsOpen Ω) (hΦ : ContDiffOn ℝ ∞ Φ Ω) : ContDiffOn ℝ ∞ (dt Φ) Ω :=
  (hΦ.fderiv_of_isOpen hΩ (m := ∞) (by simp)).clm_apply contDiffOn_const

/-- The iterated derivatives of `dt Φ` are those of `Φ` with an extra time-direction argument. -/
theorem iteratedFDeriv_dt {Φ : ℝ × E → V} {Ω : Set (ℝ × E)}
    (hΩ : IsOpen Ω) (hΦ : ContDiffOn ℝ ∞ Φ Ω) (m : ℕ) {q : ℝ × E} (hq : q ∈ Ω)
    (w : Fin m → ℝ × E) :
    iteratedFDeriv ℝ m (dt Φ) q w = iteratedFDeriv ℝ (m + 1) Φ q (Fin.snoc w (1, 0)) := by
  have hf : ContDiffAt ℝ ∞ (fderiv ℝ Φ) q :=
    (hΦ.fderiv_of_isOpen hΩ (m := ∞) (by simp)).contDiffAt (hΩ.mem_nhds hq)
  set ev : (ℝ × E →L[ℝ] V) →L[ℝ] V := ContinuousLinearMap.apply ℝ V ((1 : ℝ), (0 : E)) with hev
  have h := ev.iteratedFDeriv_comp_left hf (i := m) (natCast_le_infty m)
  have hdt : dt Φ = ev ∘ fderiv ℝ Φ := by
    funext p; simp [dt, ev]
  rw [hdt, h, iteratedFDeriv_succ_apply_right]
  simp [ev]

/-- One coordinate of the jet: derivative in the time variable. -/
theorem hasDerivAt_jet_slice_apply {Φ : ℝ × E → V}
    {Ω : Set (ℝ × E)} (hΩ : IsOpen Ω) (hΦ : ContDiffOn ℝ ∞ Φ Ω)
    (m : ℕ) {t : ℝ} {x : E} (htx : (t, x) ∈ Ω) (I : Fin m → Fin n) :
    HasDerivAt (fun s : ℝ => jet m (fun y => Φ (s, y)) x I)
      (jet m (fun y => dt Φ (t, y)) x I) t := by
  set e : Fin m → ℝ × E := fun k => (0, EuclideanSpace.single (I k) 1) with he
  have hΦt : ContDiffAt ℝ ∞ Φ (t, x) := hΦ.contDiffAt (hΩ.mem_nhds htx)
  have hDm : ContDiffAt ℝ ∞ (iteratedFDeriv ℝ m Φ) (t, x) :=
    hΦt.iteratedFDeriv_right (infty_add_natCast_le m)
  have hd : DifferentiableAt ℝ (iteratedFDeriv ℝ m Φ) (t, x) := hDm.differentiableAt (by simp)
  have hl := hd.hasFDerivAt.continuousMultilinear_apply_const e
  have hs : HasDerivAt (fun s : ℝ => (s, x)) ((1 : ℝ), (0 : E)) t :=
    (hasDerivAt_id t).prodMk (hasDerivAt_const t x)
  have h := hl.comp_hasDerivAt t hs
  have hev : (fun s : ℝ => jet m (fun y => Φ (s, y)) x I) =ᶠ[nhds t]
      ((fun p => iteratedFDeriv ℝ m Φ p e) ∘ fun s : ℝ => (s, x)) := by
    have hopen : IsOpen ((fun s : ℝ => (s, x)) ⁻¹' Ω) :=
      hΩ.preimage (continuous_id.prodMk continuous_const)
    filter_upwards [hopen.mem_nhds htx] with s hs
    rw [jet_apply, iteratedFDeriv_slice hΩ hΦ m hs]
    rfl
  refine (h.congr_of_eventuallyEq hev).congr_deriv ?_
  rw [jet_apply, iteratedFDeriv_slice hΩ (contDiffOn_dt hΩ hΦ) m htx, iteratedFDeriv_dt hΩ hΦ m htx,
    ← iteratedFDeriv_cons_eq_snoc hΩ (fun y hy => hΦ.contDiffAt (hΩ.mem_nhds hy)) m _ htx,
    iteratedFDeriv_succ_apply_left]
  rfl

/-- For `Φ` smooth on the open set `Ω ⊆ ℝ × ℝⁿ`, the `m`-jet in `x` of the slice `Φ (s, ·)` is
differentiable in `s`, with derivative the `m`-jet of the slice of `∂ₜΦ`. -/
theorem hasDerivAt_jet_slice {Φ : ℝ × E → V}
    {Ω : Set (ℝ × E)} (hΩ : IsOpen Ω) (hΦ : ContDiffOn ℝ ∞ Φ Ω)
    (m : ℕ) {t : ℝ} {x : E} (htx : (t, x) ∈ Ω) :
    HasDerivAt (fun s : ℝ => jet m (fun y => Φ (s, y)) x) (jet m (fun y => dt Φ (t, y)) x) t := by
  have hg : HasDerivAt (fun s : ℝ => fun I : Fin m → Fin n => jet m (fun y => Φ (s, y)) x I)
      (fun I : Fin m → Fin n => jet m (fun y => dt Φ (t, y)) x I) t :=
    hasDerivAt_pi.2 fun I => hasDerivAt_jet_slice_apply hΩ hΦ m htx I
  exact (PiLp.hasFDerivAt_toLp (𝕜 := ℝ) 2 (fun I : Fin m → Fin n =>
    jet m (fun y => Φ (t, y)) x I)).comp_hasDerivAt t hg

end Slice

/-- info: 'CWFlow.Extension.hasDerivAt_jet_slice' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms hasDerivAt_jet_slice

/-- info: 'CWFlow.Extension.contDiffOn_dt' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms contDiffOn_dt

end CWFlow.Extension
