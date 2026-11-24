import CWFlow.Moser.Forms
import DifferentialGeometry.Analysis.ODE.Flow.CompactSupport

/-!
# Flows of compactly supported time-dependent fields on a normed space

For `f : ℝ → W → W` smooth with compact support in `(t, x)`, the autonomized field `(1, f)` on `ℝ × W` is
constant outside a compact set, hence complete (DifferentialGeometry,
`exists_globalIntegralCurve_of_hasCompactSupport_sub_const`), and its flow is jointly smooth
(`contMDiff_curveAt`). The time-dependent flow `ψ t x` (from time `0`) is the `W`-component started at `(0, x)`:
`ψ 0 x = x`, `∂ₜ ψ t x = f t (ψ t x)`, and `(t, x) ↦ ψ t x` is smooth.
-/

open Filter Set Function

open scoped Topology ContDiff Manifold

noncomputable section

namespace CWFlow

open DifferentialGeometry.Analysis.ODE

variable {W : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W] [FiniteDimensional ℝ W]

/-- The autonomized field `(t, x) ↦ (1, f t x)`. -/
def autField (f : ℝ → W → W) (p : ℝ × W) : ℝ × W := (1, f p.1 p.2)

variable {f : ℝ → W → W}

theorem contDiff_autField (hf : ContDiff ℝ ∞ (uncurry f)) : ContDiff ℝ ∞ (autField f) :=
  contDiff_const.prodMk hf

theorem autField_complete (hf : ContDiff ℝ ∞ (uncurry f)) (hc : HasCompactSupport (uncurry f))
    (p : ℝ × W) :
    ∃ γ : ℝ → ℝ × W, γ 0 = p ∧ IsMIntegralCurve (I := 𝓘(ℝ, ℝ × W)) γ (autField f) := by
  have hsupp : HasCompactSupport (fun q => autField f q - ((1 : ℝ), (0 : W))) := by
    have e : (fun q => autField f q - ((1 : ℝ), (0 : W))) =
        (fun y : W => ((0 : ℝ), y)) ∘ uncurry f := by
      funext q
      simp [autField]
      rfl
    rw [e]
    exact hc.comp_left (by simp)
  exact exists_globalIntegralCurve_of_hasCompactSupport_sub_const (autField f)
    ((contDiff_autField hf).of_le (by simp)) _ hsupp p

/-- The flow of the autonomized field. -/
def flowAut (hf : ContDiff ℝ ∞ (uncurry f)) (hc : HasCompactSupport (uncurry f)) (t : ℝ)
    (p : ℝ × W) : ℝ × W :=
  curveAt (I := 𝓘(ℝ, ℝ × W)) (autField f) (autField_complete hf hc) p t

variable (hf : ContDiff ℝ ∞ (uncurry f)) (hc : HasCompactSupport (uncurry f))

theorem flowAut_zero (p : ℝ × W) : flowAut hf hc 0 p = p :=
  curveAt_zero (I := 𝓘(ℝ, ℝ × W)) _ _ p

theorem hasDerivAt_flowAut (p : ℝ × W) (t : ℝ) :
    HasDerivAt (fun s => flowAut hf hc s p) (autField f (flowAut hf hc t p)) t := by
  have h := curveAt_integralCurve (I := 𝓘(ℝ, ℝ × W)) (autField f) (autField_complete hf hc) p t
  rw [hasDerivAt_iff_hasFDerivAt]
  exact hasMFDerivAt_iff_hasFDerivAt.1 h

theorem contDiff_flowAut :
    ContDiff ℝ ∞ (fun q : ℝ × (ℝ × W) => flowAut hf hc q.1 q.2) := by
  have h := contMDiff_curveAt (I := 𝓘(ℝ, ℝ × W)) (autField f)
    (contMDiff_vectorSpace_iff_contDiff.2 (contDiff_autField hf)) (autField_complete hf hc)
  rw [← contMDiff_iff_contDiff, modelWithCornersSelf_prod, ← chartedSpaceSelf_prod]
  exact h

/-- The time component of the flow is `s + t`. -/
theorem flowAut_fst (s : ℝ) (x : W) (t : ℝ) : (flowAut hf hc t (s, x)).1 = s + t := by
  have hd : ∀ r, HasDerivAt (fun r => (flowAut hf hc r (s, x)).1 - r) 0 r := fun r => by
    have h := ((ContinuousLinearMap.fst ℝ ℝ W).hasFDerivAt.comp_hasDerivAt r
      (hasDerivAt_flowAut hf hc (s, x) r)).sub (hasDerivAt_id r)
    exact h.congr_deriv (by simp [autField])
  have := is_const_of_deriv_eq_zero (fun r => (hd r).differentiableAt) (fun r => (hd r).deriv) t 0
  simp only [flowAut_zero, sub_zero] at this
  linarith

/-- The time-dependent flow from time `0`. -/
def psi (t : ℝ) (x : W) : W := (flowAut hf hc t (0, x)).2

theorem psi_zero (x : W) : psi hf hc 0 x = x := by
  simp [psi, flowAut_zero]

theorem hasDerivAt_psi (x : W) (t : ℝ) :
    HasDerivAt (fun s => psi hf hc s x) (f t (psi hf hc t x)) t := by
  have h := (ContinuousLinearMap.snd ℝ ℝ W).hasFDerivAt.comp_hasDerivAt t
    (hasDerivAt_flowAut hf hc (0, x) t)
  have e : (flowAut hf hc t (0, x)).1 = t := by rw [flowAut_fst]; ring
  exact h.congr_deriv (by simp [autField, psi, e])

theorem contDiff_psi : ContDiff ℝ ∞ (fun q : ℝ × W => psi hf hc q.1 q.2) :=
  contDiff_snd.comp ((contDiff_flowAut hf hc).comp
    (contDiff_fst.prodMk (contDiff_const.prodMk contDiff_snd)))

/-! ### The variational equation, from symmetry of second derivatives -/

omit [FiniteDimensional ℝ W] in
/-- `D_y f t = D(uncurry f)(t, y) ∘ inr`. -/
theorem fderiv_slice_apply {g : ℝ × W → W} {t : ℝ} {y : W} (hg : DifferentiableAt ℝ g (t, y))
    (w : W) : fderiv ℝ (fun z => g (t, z)) y w = fderiv ℝ g (t, y) (0, w) := by
  have h : HasFDerivAt (fun z => g (t, z)) ((fderiv ℝ g (t, y)).comp (ContinuousLinearMap.inr ℝ ℝ W))
      y := hg.hasFDerivAt.comp y (hasFDerivAt_prodMk_right t y)
  rw [h.fderiv]
  rfl

omit [FiniteDimensional ℝ W] in
/-- **Variational equation.** If `Ψ : ℝ × W → W` is `C²` and `∂ₜΨ = f(t, Ψ)` with `f` `C¹`, then the spatial
derivative `Dψ_t = D_x Ψ(t, ·)` satisfies `∂ₜ Dψ_t = D f_t(Ψ) ∘ Dψ_t`. -/
theorem hasDerivAt_spatial_fderiv {Ψ : ℝ × W → W} (hΨ : ContDiff ℝ 2 Ψ)
    (hf1 : ContDiff ℝ 1 (uncurry f))
    (hdt : ∀ t x, HasDerivAt (fun s => Ψ (s, x)) (f t (Ψ (t, x))) t) (t : ℝ) (x : W) :
    HasDerivAt (fun s => fderiv ℝ (fun y => Ψ (s, y)) x)
      ((fderiv ℝ (f t) (Ψ (t, x))).comp (fderiv ℝ (fun y => Ψ (t, y)) x)) t := by
  have hΨd : ∀ q, DifferentiableAt ℝ Ψ q := fun q => hΨ.differentiable (by simp) q
  have hD := (hΨ.fderiv_right (m := 1) (by norm_num)).differentiable one_ne_zero
  -- `s ↦ DΨ(s, x)` and its derivative
  have hcurve := (hD (t, x)).hasFDerivAt.comp_hasDerivAt t
    ((hasDerivAt_id t).prodMk (hasDerivAt_const t x))
  have hslice : ∀ s, fderiv ℝ (fun y => Ψ (s, y)) x =
      (fderiv ℝ Ψ (s, x)).comp (ContinuousLinearMap.inr ℝ ℝ W) := fun s =>
    ContinuousLinearMap.ext fun w => fderiv_slice_apply (hΨd (s, x)) w
  simp only [hslice]
  have hpre := ((ContinuousLinearMap.compL ℝ W (ℝ × W) W).flip
    (ContinuousLinearMap.inr ℝ ℝ W)).hasFDerivAt.comp_hasDerivAt t hcurve
  refine hpre.congr_deriv (ContinuousLinearMap.ext fun w => ?_)
  -- `∂ₜΨ = f(t, Ψ)` as functions of `q = (t, x)`
  have hdt' : ∀ q : ℝ × W, fderiv ℝ Ψ q (1, 0) = f q.1 (Ψ q) := fun q => by
    have h1 := (hΨd q).hasFDerivAt.comp_hasDerivAt q.1
      ((hasDerivAt_id q.1).prodMk (hasDerivAt_const q.1 q.2))
    have h2 := hdt q.1 q.2
    have := h1.unique h2
    simpa using this
  have hsym := (hΨ.contDiffAt (x := (t, x))).isSymmSndFDerivAt (n := 2) (by simp) ((1 : ℝ), (0 : W)) ((0 : ℝ), w)
  -- derivative of `q ↦ f q.1 (Ψ q)` in the direction `(0, w)`
  have hfq : HasFDerivAt (fun q : ℝ × W => f q.1 (Ψ q))
      ((fderiv ℝ (uncurry f) (t, Ψ (t, x))).comp
        ((ContinuousLinearMap.fst ℝ ℝ W).prod (fderiv ℝ Ψ (t, x)))) (t, x) := by
    have h := ((hf1.differentiable one_ne_zero) (t, Ψ (t, x))).hasFDerivAt.comp (t, x)
      ((hasFDerivAt_fst).prodMk (hΨd (t, x)).hasFDerivAt)
    exact h
  have hev : (fun q => fderiv ℝ Ψ q (1, 0)) = fun q : ℝ × W => f q.1 (Ψ q) := funext hdt'
  have hD2 : fderiv ℝ (fderiv ℝ Ψ) (t, x) ((0 : ℝ), w) ((1 : ℝ), (0 : W)) =
      fderiv ℝ (fun q => fderiv ℝ Ψ q (1, 0)) (t, x) (0, w) := by
    rw [fderiv_clm_apply (hD (t, x)) (differentiableAt_const _)]
    simp
  rw [hev, hfq.fderiv] at hD2
  have hfs : fderiv ℝ (f t) (Ψ (t, x)) (fderiv ℝ Ψ (t, x) (0, w)) =
      fderiv ℝ (uncurry f) (t, Ψ (t, x)) (0, fderiv ℝ Ψ (t, x) (0, w)) := by
    have := fderiv_slice_apply (g := uncurry f) (t := t) (y := Ψ (t, x))
      ((hf1.differentiable one_ne_zero) _) (fderiv ℝ Ψ (t, x) (0, w))
    simpa [uncurry] using this
  simp only [ContinuousLinearMap.comp_apply, ContinuousLinearMap.compL_apply,
    ContinuousLinearMap.flip_apply, ContinuousLinearMap.inr_apply]
  rw [hsym, hD2, hfs]
  simp

theorem hasDerivAt_fderiv_psi (x : W) (t : ℝ) :
    HasDerivAt (fun s => fderiv ℝ (psi hf hc s) x)
      ((fderiv ℝ (f t) (psi hf hc t x)).comp (fderiv ℝ (psi hf hc t) x)) t :=
  hasDerivAt_spatial_fderiv (Ψ := fun q => psi hf hc q.1 q.2) ((contDiff_psi hf hc).of_le (WithTop.coe_le_coe.2 le_top))
    (hf.of_le (WithTop.coe_le_coe.2 le_top)) (fun t x => hasDerivAt_psi hf hc x t) t x

/-! ### Moser's identity along the flow -/

theorem fderiv_psi_zero (x : W) : fderiv ℝ (psi hf hc 0) x = ContinuousLinearMap.id ℝ W := by
  have h : psi hf hc 0 = id := funext (psi_zero hf hc)
  rw [h, fderiv_id]

/-- **Moser along the flow.** If the Moser condition `Ω̇ + L_f Ω = 0` holds along the trajectory of `x` for
`s ∈ [0, 1]`, then `(ψ₁*Ω₁)(x) = Ω₀(x)`. -/
theorem moser_flow {Ω Ωdot : ℝ → W → W [⋀^Fin 2]→L[ℝ] ℝ} (x : W)
    (hcurve : ∀ s ∈ Icc (0 : ℝ) 1, HasDerivAt (fun r => Ω r (psi hf hc r x))
      (Ωdot s (psi hf hc s x) + fderiv ℝ (Ω s) (psi hf hc s x) (f s (psi hf hc s x))) s)
    (hmoser : ∀ s ∈ Icc (0 : ℝ) 1, ∀ a b,
      Ωdot s (psi hf hc s x) ![a, b] + lieAlt (Ω s) (f s) (psi hf hc s x) a b = 0)
    (v w : W) :
    Ω 1 (psi hf hc 1 x) ![fderiv ℝ (psi hf hc 1) x v, fderiv ℝ (psi hf hc 1) x w] = Ω 0 x ![v, w] := by
  set g : ℝ → ℝ := fun r =>
    Ω r (psi hf hc r x) ![fderiv ℝ (psi hf hc r) x v, fderiv ℝ (psi hf hc r) x w] with hg
  have hd : ∀ s ∈ Icc (0 : ℝ) 1, HasDerivAt g 0 s := fun s hs => by
    have h := hasDerivAt_pullback_alt (X := f) (ψ := fun r => psi hf hc r x)
      (Dψ := fun r => fderiv ℝ (psi hf hc r) x) (hcurve s hs) (hasDerivAt_fderiv_psi hf hc x s) v w
    rwa [hmoser s hs] at h
  have hconst := constant_of_has_deriv_right_zero (f := g) (a := 0) (b := 1)
    (fun s hs => (hd s hs).continuousAt.continuousWithinAt)
    (fun s hs => (hd s (Ico_subset_Icc_self hs)).hasDerivWithinAt) 1 ⟨zero_le_one, le_rfl⟩
  rw [hg] at hconst
  simp only at hconst
  rw [hconst, psi_zero, fderiv_psi_zero]
  rfl

/-! ### Fixed points and stability near them -/

include hf hc in
theorem lipschitzWith_slice : ∃ C, ∀ t, LipschitzWith C (f t) := by
  obtain ⟨C, hC⟩ := ContDiff.lipschitzWith_of_hasCompactSupport hc hf (by simp)
  refine ⟨C, fun t => ?_⟩
  have h := hC.comp (LipschitzWith.prodMk_left t)
  rw [mul_one] at h
  exact h

/-- A point where the field vanishes for all time is fixed by the flow. -/
theorem psi_fixed {k : W} (hk : ∀ s, f s k = 0) {s : ℝ} (hs : 0 ≤ s) : psi hf hc s k = k := by
  obtain ⟨C, hC⟩ := lipschitzWith_slice hf hc
  have h := ODE_solution_unique_of_mem_Icc_right (v := f) (s := fun _ => Set.univ) (K := C)
    (f := fun r => psi hf hc r k) (g := fun _ => k) (a := 0) (b := s)
    (fun t _ => (hC t).lipschitzOnWith)
    (fun r _ => (hasDerivAt_psi hf hc k r).continuousAt.continuousWithinAt)
    (fun r _ => (hasDerivAt_psi hf hc k r).hasDerivWithinAt) (fun _ _ => Set.mem_univ _)
    continuousOn_const (fun r _ => by simpa [hk r] using hasDerivWithinAt_const r (Ici r) k)
    (fun _ _ => Set.mem_univ _) (psi_zero hf hc k)
  exact h ⟨hs, le_rfl⟩

/-- Points near a compact set `K` of fixed points stay in any open `U ⊇ K` for `s ∈ [0, 1]`. -/
theorem exists_nhds_psi_mem {K : Set W} (hK : IsCompact K) (hfix : ∀ k ∈ K, ∀ s, f s k = 0)
    {U : Set W} (hU : IsOpen U) (hKU : K ⊆ U) :
    ∃ V : Set W, IsOpen V ∧ K ⊆ V ∧ ∀ x ∈ V, ∀ s ∈ Icc (0 : ℝ) 1, psi hf hc s x ∈ U := by
  have hcont : Continuous (fun q : ℝ × W => psi hf hc q.1 q.2) := (contDiff_psi hf hc).continuous
  have hn : IsOpen {q : ℝ × W | psi hf hc q.1 q.2 ∈ U} := hU.preimage hcont
  have hp : Icc (0 : ℝ) 1 ×ˢ K ⊆ {q : ℝ × W | psi hf hc q.1 q.2 ∈ U} := by
    rintro ⟨s, k⟩ ⟨hs, hk⟩
    simp only [Set.mem_setOf_eq]
    rw [psi_fixed hf hc (hfix k hk) hs.1]
    exact hKU hk
  obtain ⟨u, v, -, hv, hsu, hKv, huv⟩ := generalized_tube_lemma isCompact_Icc hK hn hp
  exact ⟨v, hv, hKv, fun x hx s hs => huv (Set.mk_mem_prod (hsu hs) hx)⟩

end CWFlow
