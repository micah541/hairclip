import CWFlow.Moser.Theorem

/-!
# Naturality of the canonical Moser map

Let `T : W → W` be smooth on an open `N`, with `DT` invertible there. Suppose `T` pulls the data `(λ₂, β₂)` back
to `(λ₁, β₁)` on `N`: `λ₁ = T*λ₂` and `β₁ = T*β₂` (`pb1form`). Then the following hold.

* `T*Ω₂,t = Ω₁,t` (`omegaT_natural`); this uses `d T* = T* d` (`dOne_pb1form`, from the symmetry of `D²T`).
* `T` maps the good set to the good set (`good_natural`).
* `T` intertwines the Moser fields: `X₂(T q) = DT X₁(q)` (`moserField_natural`).
* `T` maps Moser curves staying in `N` to Moser curves (`moserCurve_natural`), so
  `moserEnd₂ (T q) = T (moserEnd₁ q)` (`moserEnd_natural`).

This is what makes the chart-level Weinstein maps glue.
-/

open ContinuousAlternatingMap Filter Set Function

open scoped Topology ContDiff

noncomputable section

namespace CWFlow

variable {W : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W]

/-- Pullback of a field of linear functionals: `(T*β)(q) = β(T q) ∘ DT(q)`. -/
def pb1form (T : W → W) (β : W → W →L[ℝ] ℝ) (q : W) : W →L[ℝ] ℝ := (β (T q)).comp (fderiv ℝ T q)

@[simp] theorem pb1form_apply (T : W → W) (β : W → W →L[ℝ] ℝ) (q a : W) :
    pb1form T β q a = β (T q) (fderiv ℝ T q a) := rfl

/-- `d(T*β) = T*(dβ)`. -/
theorem dOne_pb1form {β : W → W →L[ℝ] ℝ} {T : W → W} {q : W} (hβ : DifferentiableAt ℝ β (T q))
    (hT : ContDiffAt ℝ 2 T q) (a b : W) :
    dOne (pb1form T β) q a b = dOne β (T q) (fderiv ℝ T q a) (fderiv ℝ T q b) := by
  have hTd : DifferentiableAt ℝ T q := hT.differentiableAt (by norm_num)
  have hD : DifferentiableAt ℝ (fderiv ℝ T) q :=
    (hT.fderiv_right (m := 1) (by norm_num)).differentiableAt (by norm_num)
  have h : HasFDerivAt (fun p => (β (T p)).comp (fderiv ℝ T p)) _ q :=
    (hβ.hasFDerivAt.comp q hTd.hasFDerivAt).clm_comp hD.hasFDerivAt
  have e : pb1form T β = fun p => (β (T p)).comp (fderiv ℝ T p) := rfl
  simp only [dOne]
  rw [e, h.fderiv]
  simp only [ContinuousLinearMap.add_apply, ContinuousLinearMap.coe_comp', Function.comp_apply,
    ContinuousLinearMap.compL_apply, ContinuousLinearMap.flip_apply]
  rw [hT.isSymmSndFDerivAt (by simp) a b]
  ring

theorem omegaT_apply {lam β : W → W →L[ℝ] ℝ} {q : W} (hlam : DifferentiableAt ℝ lam q)
    (hβ : DifferentiableAt ℝ β q) (t : ℝ) (a b : W) :
    omegaT lam β t q ![a, b] = dOne lam q a b + t * dOne β q a b := by
  simp only [omegaT, ContinuousAlternatingMap.add_apply, ContinuousAlternatingMap.smul_apply,
    smul_eq_mul]
  rw [extDeriv_toAlt1 hlam, extDeriv_toAlt1 hβ]

theorem alt1_ext {f g : W [⋀^Fin 1]→L[ℝ] ℝ} (h : ∀ b, f ![b] = g ![b]) : f = g := by
  ext w
  have hw : w = ![w 0] := by funext i; fin_cases i; rfl
  rw [hw, h]

variable [FiniteDimensional ℝ W]

section Natural

variable {lam₁ β₁ lam₂ β₂ : W → W →L[ℝ] ℝ} {O₁ O₂ N : Set W} {T : W → W}
variable (hO₂ : IsOpen O₂) (hlam₂ : ContDiffOn ℝ ∞ lam₂ O₂) (hβ₂ : ContDiffOn ℝ ∞ β₂ O₂)
  (hN : IsOpen N) (hTN : MapsTo T N O₂) (hT : ContDiffOn ℝ ∞ T N)
  (hTinv : ∀ q ∈ N, ∃ e : W ≃L[ℝ] W, (e : W →L[ℝ] W) = fderiv ℝ T q)
  (hlamT : ∀ q ∈ N, lam₁ q = pb1form T lam₂ q) (hβT : ∀ q ∈ N, β₁ q = pb1form T β₂ q)
include hO₂ hlam₂ hβ₂ hN hTN hT hlamT hβT

omit [FiniteDimensional ℝ W] in
/-- `T*Ω₂,t = Ω₁,t` on `N`. -/
theorem omegaT_natural {q : W} (hq : q ∈ N) (t : ℝ) (a b : W) :
    omegaT lam₂ β₂ t (T q) ![fderiv ℝ T q a, fderiv ℝ T q b] = omegaT lam₁ β₁ t q ![a, b] := by
  have hq₂ : T q ∈ O₂ := hTN hq
  have hl₂ := (hlam₂.contDiffAt (hO₂.mem_nhds hq₂)).differentiableAt (by simp)
  have hb₂ := (hβ₂.contDiffAt (hO₂.mem_nhds hq₂)).differentiableAt (by simp)
  have hT2 : ContDiffAt ℝ 2 T q :=
    (hT.contDiffAt (hN.mem_nhds hq)).of_le (WithTop.coe_le_coe.2 le_top)
  have hl₁e : lam₁ =ᶠ[𝓝 q] pb1form T lam₂ := by
    filter_upwards [hN.mem_nhds hq] with p hp using hlamT p hp
  have hb₁e : β₁ =ᶠ[𝓝 q] pb1form T β₂ := by
    filter_upwards [hN.mem_nhds hq] with p hp using hβT p hp
  have hT1 : DifferentiableAt ℝ T q := hT2.differentiableAt (by norm_num)
  have hD : DifferentiableAt ℝ (fderiv ℝ T) q :=
    (hT2.fderiv_right (m := 1) (by norm_num)).differentiableAt (by norm_num)
  have hl₁ : DifferentiableAt ℝ lam₁ q :=
    ((hl₂.hasFDerivAt.comp q hT1.hasFDerivAt).clm_comp hD.hasFDerivAt).differentiableAt.congr_of_eventuallyEq
      hl₁e
  have hb₁ : DifferentiableAt ℝ β₁ q :=
    ((hb₂.hasFDerivAt.comp q hT1.hasFDerivAt).clm_comp hD.hasFDerivAt).differentiableAt.congr_of_eventuallyEq
      hb₁e
  have e1 : dOne lam₁ q a b = dOne lam₂ (T q) (fderiv ℝ T q a) (fderiv ℝ T q b) := by
    rw [dOne, hl₁e.fderiv_eq, ← dOne, dOne_pb1form hl₂ hT2]
  have e2 : dOne β₁ q a b = dOne β₂ (T q) (fderiv ℝ T q a) (fderiv ℝ T q b) := by
    rw [dOne, hb₁e.fderiv_eq, ← dOne, dOne_pb1form hb₂ hT2]
  rw [omegaT_apply hl₂ hb₂, omegaT_apply hl₁ hb₁, e1, e2]

include hTinv

/-- `T` maps the good set to the good set. -/
theorem good_natural {t : ℝ} {q : W} (hq : q ∈ N) (hg : (t, q) ∈ moserGood lam₁ β₁ O₁) :
    (t, T q) ∈ moserGood lam₂ β₂ O₂ := by
  obtain ⟨e, he⟩ := hTinv q hq
  obtain ⟨e₁, he₁⟩ := hg.2
  have hinj : Function.Injective (omegaT lam₂ β₂ t (T q)).curryLeft := by
    rw [injective_iff_map_eq_zero]
    intro v hv
    have hv' : v = fderiv ℝ T q (e.symm v) := by rw [← he]; simp
    have hz : iotaT lam₁ β₁ t q (e.symm v) = 0 := alt1_ext fun b => by
      have h := congrArg (fun f => f ![fderiv ℝ T q b]) hv
      simp only [curryLeft_apply_apply] at h
      rw [hv'] at h
      have h2 : iotaT lam₁ β₁ t q (e.symm v) ![b] = omegaT lam₁ β₁ t q ![e.symm v, b] := rfl
      rw [h2, ← omegaT_natural hO₂ hlam₂ hβ₂ hN hTN hT hlamT hβT hq]
      simpa using h
    rw [← he₁] at hz
    have h0 : e.symm v = 0 := e₁.injective (hz.trans (_root_.map_zero e₁).symm)
    rw [hv', h0, _root_.map_zero]
  exact ⟨hTN hq, curryEquiv _ hinj, by rw [coe_curryEquiv]; rfl⟩

/-- `T` intertwines the Moser fields. -/
theorem moserField_natural {t : ℝ} {q : W} (hq : q ∈ N) (hg : (t, q) ∈ moserGood lam₁ β₁ O₁) :
    moserField lam₂ β₂ t (T q) = fderiv ℝ T q (moserField lam₁ β₁ t q) := by
  obtain ⟨e, he⟩ := hTinv q hq
  obtain ⟨e₂, he₂⟩ := (good_natural hO₂ hlam₂ hβ₂ hN hTN hT hTinv hlamT hβT hq hg).2
  have key : iotaT lam₂ β₂ t (T q) (moserField lam₂ β₂ t (T q)) =
      iotaT lam₂ β₂ t (T q) (fderiv ℝ T q (moserField lam₁ β₁ t q)) := alt1_ext fun c => by
    have hc : c = fderiv ℝ T q (e.symm c) := by rw [← he]; simp
    have h1 : iotaT lam₂ β₂ t (T q) (moserField lam₂ β₂ t (T q)) ![c] =
        omegaT lam₂ β₂ t (T q) ![moserField lam₂ β₂ t (T q), c] := rfl
    have h2 : iotaT lam₂ β₂ t (T q) (fderiv ℝ T q (moserField lam₁ β₁ t q)) ![c] =
        omegaT lam₂ β₂ t (T q) ![fderiv ℝ T q (moserField lam₁ β₁ t q), c] := rfl
    rw [h1, h2, moserField_spec ⟨e₂, he₂⟩ c, hc,
      omegaT_natural hO₂ hlam₂ hβ₂ hN hTN hT hlamT hβT hq, moserField_spec hg.2, hβT q hq,
      pb1form_apply]
  rw [← he₂] at key
  exact e₂.injective key

/-- `T` maps Moser curves staying in `N` to Moser curves. -/
theorem moserCurve_natural {q : W} {γ : ℝ → W} (h : IsMoserCurve lam₁ β₁ O₁ q γ)
    (hγN : ∀ t ∈ Icc (0 : ℝ) 1, γ t ∈ N) : IsMoserCurve lam₂ β₂ O₂ (T q) (T ∘ γ) where
  start := by rw [Function.comp_apply, h.start]
  good t ht := good_natural hO₂ hlam₂ hβ₂ hN hTN hT hTinv hlamT hβT (hγN t ht) (h.good t ht)
  deriv t ht := by
    have hTd : DifferentiableAt ℝ T (γ t) :=
      (hT.contDiffAt (hN.mem_nhds (hγN t ht))).differentiableAt (by simp)
    rw [Function.comp_apply, moserField_natural hO₂ hlam₂ hβ₂ hN hTN hT hTinv hlamT hβT (hγN t ht)
      (h.good t ht)]
    exact hTd.hasFDerivAt.comp_hasDerivWithinAt t (h.deriv t ht)

/-- **Naturality of the canonical time-one map.** -/
theorem moserEnd_natural (hO₁ : IsOpen O₁) (hlam₁ : ContDiffOn ℝ ∞ lam₁ O₁)
    (hβ₁ : ContDiffOn ℝ ∞ β₁ O₁) {q : W} {γ : ℝ → W} (h : IsMoserCurve lam₁ β₁ O₁ q γ)
    (hγN : ∀ t ∈ Icc (0 : ℝ) 1, γ t ∈ N) :
    moserEnd lam₂ β₂ O₂ (T q) = T (moserEnd lam₁ β₁ O₁ q) := by
  rw [moserEnd_eq hO₂ hlam₂ hβ₂ (moserCurve_natural hO₂ hlam₂ hβ₂ hN hTN hT hTinv hlamT hβT h hγN),
    moserEnd_eq hO₁ hlam₁ hβ₁ h, Function.comp_apply]

end Natural

end CWFlow
