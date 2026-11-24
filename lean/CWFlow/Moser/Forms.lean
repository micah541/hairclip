import Mathlib

/-!
# Moser's argument for 2-forms on a normed space

Forms are Mathlib's `W → W [⋀^Fin k]→L[ℝ] ℝ` with `extDeriv`. For a 2-form `Ω` and a vector field `X`,
`lieAlt Ω X p a b = DΩ[X](a, b) + Ω(DX a, b) + Ω(a, DX b)` is `(L_X Ω)(a, b)`.

* `moser_condition`: if `Ω_t = ω + t·dβ` (`ω` constant, `β ∈ C²`) and `ι_{X} Ω_t = −β` near `p`, then
  `dβ + L_X Ω_t = 0` at `p`. (Cartan's formula, used only through scalar identities: no derivative of a
  form-valued interior product is needed.)
* `hasDerivAt_pullback`: along a flow `ψ_s` of `X_s` with variational equation for `Dψ_s`,
  `d/ds Ω_s(ψ_s)(Dψ_s v, Dψ_s w) = (Ω̇_s + L_{X_s} Ω_s)(Dψ_s v, Dψ_s w)`; `pullback_const` concludes.

Note: type ascriptions `DifferentiableAt ℝ f p` for form-valued `f` select topological-vector-space
instances that do not unify with the normed ones used by the calculus lemmas; such facts are obtained
from terms, never ascribed.
-/

open ContinuousAlternatingMap Filter

open scoped Topology

noncomputable section

namespace CWFlow

variable {W : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W]

/-- `(L_X Ω)(a, b) = DΩ[X](a, b) + Ω(DX a, b) + Ω(a, DX b)`. -/
def lieAlt (Ω : W → W [⋀^Fin 2]→L[ℝ] ℝ) (X : W → W) (p a b : W) : ℝ :=
  fderiv ℝ Ω p (X p) ![a, b] + Ω p ![fderiv ℝ X p a, b] + Ω p ![a, fderiv ℝ X p b]

theorem alt2_swap (f : W [⋀^Fin 2]→L[ℝ] ℝ) (a b : W) : f ![a, b] = -f ![b, a] := by
  have h := f.map_swap ![b, a] (show (0 : Fin 2) ≠ 1 by decide)
  have e : (![b, a] ∘ Equiv.swap (0 : Fin 2) 1) = ![a, b] := by
    ext i
    fin_cases i <;> rfl
  rw [e] at h
  exact h

/-- `dβ(a, b) = Dβ[a](b) − Dβ[b](a)` for a 1-form. -/
theorem extDeriv_one_apply {β : W → W [⋀^Fin 1]→L[ℝ] ℝ} {p : W} (hβ : DifferentiableAt ℝ β p)
    (a b : W) : extDeriv β p ![a, b] = fderiv ℝ β p a ![b] - fderiv ℝ β p b ![a] := by
  rw [extDeriv_apply hβ, Fin.sum_univ_two]
  simp only [fderiv_continuousAlternatingMap_apply_const_apply hβ]
  have e0 : Fin.removeNth 0 ![a, b] = ![b] := by ext i; fin_cases i; rfl
  have e1 : Fin.removeNth 1 ![a, b] = ![a] := by ext i; fin_cases i; rfl
  rw [e0, e1]
  simp
  ring

/-- `dΩ(x, a, b) = DΩ[x](a, b) − DΩ[a](x, b) + DΩ[b](x, a)` for a 2-form. -/
theorem extDeriv_two_apply {Ω : W → W [⋀^Fin 2]→L[ℝ] ℝ} {p : W} (hΩ : DifferentiableAt ℝ Ω p)
    (x a b : W) :
    extDeriv Ω p ![x, a, b] =
      fderiv ℝ Ω p x ![a, b] - fderiv ℝ Ω p a ![x, b] + fderiv ℝ Ω p b ![x, a] := by
  rw [extDeriv_apply hΩ, Fin.sum_univ_three]
  simp only [fderiv_continuousAlternatingMap_apply_const_apply hΩ]
  have e0 : Fin.removeNth 0 ![x, a, b] = ![a, b] := by ext i; fin_cases i <;> rfl
  have e1 : Fin.removeNth 1 ![x, a, b] = ![x, b] := by ext i; fin_cases i <;> rfl
  have e2 : Fin.removeNth 2 ![x, a, b] = ![x, a] := by ext i; fin_cases i <;> rfl
  rw [e0, e1, e2]
  simp
  ring

/-- `D(q ↦ Ω q (X q, b))[a] = DΩ[a](X, b) + Ω(DX a, b)`. -/
theorem fderiv_iota_alt {Ω : W → W [⋀^Fin 2]→L[ℝ] ℝ} {X : W → W} {p : W}
    (hΩ : DifferentiableAt ℝ Ω p) (hX : DifferentiableAt ℝ X p) (a b : W) :
    fderiv ℝ (fun q => Ω q ![X q, b]) p a =
      fderiv ℝ Ω p a ![X p, b] + Ω p ![fderiv ℝ X p a, b] := by
  have hg : ∀ i : Fin 2, HasFDerivAt (fun q => (![X q, b] : Fin 2 → W) i)
      ((![fderiv ℝ X p, 0] : Fin 2 → W →L[ℝ] W) i) p := by
    intro i
    fin_cases i
    · exact hX.hasFDerivAt
    · exact hasFDerivAt_const b p
  have h := hΩ.hasFDerivAt.continuousAlternatingMap_apply hg
  have e : (fun q => Ω q ![X q, b]) = fun q => Ω q (fun i => (![X q, b] : Fin 2 → W) i) := rfl
  rw [e, h.fderiv]
  have u0 : Function.update ![X p, b] 0 (fderiv ℝ X p a) = ![fderiv ℝ X p a, b] := by
    ext i; fin_cases i <;> rfl
  have u1 : Function.update ![X p, b] 1 (0 : W) = ![X p, 0] := by
    ext i; fin_cases i <;> rfl
  simp [Fin.sum_univ_two, u0, u1]

/-- **The Moser condition.** `Ω_t = ω + t·dβ`, `ι_X Ω_t = −β` near `p` ⟹ `dβ + L_X Ω_t = 0` at `p`. -/
theorem moser_condition (ω : W [⋀^Fin 2]→L[ℝ] ℝ) {β : W → W [⋀^Fin 1]→L[ℝ] ℝ} {X : W → W}
    {p : W} {t : ℝ} (hβ : ContDiffAt ℝ 2 β p) (hX : DifferentiableAt ℝ X p)
    (hι : ∀ᶠ q in 𝓝 p, ∀ b, (ω + t • extDeriv β q) ![X q, b] = -β q ![b]) (a b : W) :
    extDeriv β p ![a, b] + lieAlt (fun q => ω + t • extDeriv β q) X p a b = 0 := by
  have hβd := hβ.differentiableAt (by norm_num)
  have hd := (hβ.fderiv_right (m := 1) (by norm_num)).differentiableAt (by norm_num)
  have hdβ : DifferentiableAt ℝ (extDeriv β) p := by
    exact alternatizeUncurryFinCLM _ _ _ |>.hasFDerivAt.comp p hd.hasFDerivAt |>.differentiableAt
  have hsm : HasFDerivAt (t • extDeriv β) (t • fderiv ℝ (extDeriv β) p) p :=
    hdβ.hasFDerivAt.const_smul t
  have hΩ' : HasFDerivAt (fun q => ω + t • extDeriv β q) (t • fderiv ℝ (extDeriv β) p) p :=
    hsm.const_add ω
  have hΩd := hΩ'.differentiableAt
  -- closedness of `Ω_t`
  have hclosed : extDeriv (fun q => ω + t • extDeriv β q) p = 0 := by
    show extDeriv ((fun _ => ω) + t • extDeriv β) p = 0
    rw [extDeriv_add (differentiableAt_const ω) hsm.differentiableAt, extDeriv_smul,
      extDeriv_extDeriv_apply hβ (by simp)]
    simp only [extDeriv, fderiv_const, Pi.zero_apply, smul_zero, add_zero]
    ext v
    simp [alternatizeUncurryFin_apply]
  have hc : extDeriv (fun q => ω + t • extDeriv β q) p ![X p, a, b] = 0 := by rw [hclosed]; rfl
  rw [extDeriv_two_apply hΩd] at hc
  -- the scalar identities `D(Ω(X, c))[d] = −Dβ[d](c)`
  have hs : ∀ c d : W, fderiv ℝ (fun q => ω + t • extDeriv β q) p d ![X p, c] +
      (ω + t • extDeriv β p) ![fderiv ℝ X p d, c] = -fderiv ℝ β p d ![c] := fun c d => by
    rw [← fderiv_iota_alt hΩd hX d c]
    have hev : (fun q => (ω + t • extDeriv β q) ![X q, c]) =ᶠ[𝓝 p] fun q => -β q ![c] :=
      hι.mono fun q hq => hq c
    rw [hev.fderiv_eq, fderiv_fun_neg, ContinuousLinearMap.neg_apply,
      fderiv_continuousAlternatingMap_apply_const_apply hβd]
  have h1 := hs b a
  have h2 := hs a b
  have hsw := alt2_swap (ω + t • extDeriv β p) a (fderiv ℝ X p b)
  rw [lieAlt, extDeriv_one_apply hβd, hsw]
  linarith

/-- **The Moser condition, exact form.** If `Ω = dγ` (`γ ∈ C²`) and `ι_X Ω = −β` near `p` (`β ∈ C¹`), then
`dβ + L_X Ω = 0` at `p`. (For `γ = λ + tβ` this is the Moser condition for `Ω_t = dλ + t·dβ`.) -/
theorem moser_condition_exact {γ β : W → W [⋀^Fin 1]→L[ℝ] ℝ} {X : W → W} {p : W}
    (hγ : ContDiffAt ℝ 2 γ p) (hβd : DifferentiableAt ℝ β p) (hX : DifferentiableAt ℝ X p)
    (hι : ∀ᶠ q in 𝓝 p, ∀ b, extDeriv γ q ![X q, b] = -β q ![b]) (a b : W) :
    extDeriv β p ![a, b] + lieAlt (extDeriv γ) X p a b = 0 := by
  have hd := (hγ.fderiv_right (m := 1) (by norm_num)).differentiableAt (by norm_num)
  have hΩd : DifferentiableAt ℝ (extDeriv γ) p := by
    exact alternatizeUncurryFinCLM _ _ _ |>.hasFDerivAt.comp p hd.hasFDerivAt |>.differentiableAt
  have hc : extDeriv (extDeriv γ) p ![X p, a, b] = 0 := by
    rw [extDeriv_extDeriv_apply hγ (by simp)]; rfl
  rw [extDeriv_two_apply hΩd] at hc
  have hs : ∀ c d : W, fderiv ℝ (extDeriv γ) p d ![X p, c] + extDeriv γ p ![fderiv ℝ X p d, c] =
      -fderiv ℝ β p d ![c] := fun c d => by
    rw [← fderiv_iota_alt hΩd hX d c]
    have hev : (fun q => extDeriv γ q ![X q, c]) =ᶠ[𝓝 p] fun q => -β q ![c] :=
      hι.mono fun q hq => hq c
    rw [hev.fderiv_eq, fderiv_fun_neg, ContinuousLinearMap.neg_apply,
      fderiv_continuousAlternatingMap_apply_const_apply hβd]
  have h1 := hs b a
  have h2 := hs a b
  have hsw := alt2_swap (extDeriv γ p) a (fderiv ℝ X p b)
  rw [lieAlt, extDeriv_one_apply hβd, hsw]
  linarith

/-- **The pullback along a flow.** If `s ↦ Ω_s(ψ_s)` has derivative `Ω̇ + DΩ[X]` and `Dψ_s` solves the
variational equation, then `d/ds Ω_s(ψ_s)(Dψ_s v, Dψ_s w) = (Ω̇ + L_X Ω)(Dψ v, Dψ w)`. -/
theorem hasDerivAt_pullback_alt {Ω Ωdot : ℝ → W → W [⋀^Fin 2]→L[ℝ] ℝ} {X : ℝ → W → W}
    {ψ : ℝ → W} {Dψ : ℝ → W →L[ℝ] W} {s : ℝ}
    (hcurve : HasDerivAt (fun r => Ω r (ψ r)) (Ωdot s (ψ s) + fderiv ℝ (Ω s) (ψ s) (X s (ψ s))) s)
    (hDψ : HasDerivAt Dψ ((fderiv ℝ (X s) (ψ s)).comp (Dψ s)) s) (v w : W) :
    HasDerivAt (fun r => Ω r (ψ r) ![Dψ r v, Dψ r w])
      (Ωdot s (ψ s) ![Dψ s v, Dψ s w] + lieAlt (Ω s) (X s) (ψ s) (Dψ s v) (Dψ s w)) s := by
  have hv : HasDerivAt (fun r => Dψ r v) (fderiv ℝ (X s) (ψ s) (Dψ s v)) s := by
    simpa using hDψ.clm_apply (hasDerivAt_const s v)
  have hw : HasDerivAt (fun r => Dψ r w) (fderiv ℝ (X s) (ψ s) (Dψ s w)) s := by
    simpa using hDψ.clm_apply (hasDerivAt_const s w)
  have hg : ∀ i : Fin 2, HasFDerivAt (fun r => (![Dψ r v, Dψ r w] : Fin 2 → W) i)
      ((1 : ℝ →L[ℝ] ℝ).smulRight
        ((![fderiv ℝ (X s) (ψ s) (Dψ s v), fderiv ℝ (X s) (ψ s) (Dψ s w)] : Fin 2 → W) i)) s := by
    intro i
    fin_cases i
    · exact hv.hasFDerivAt
    · exact hw.hasFDerivAt
  have h := (hcurve.hasFDerivAt.continuousAlternatingMap_apply hg).hasDerivAt
  have e : (fun r => Ω r (ψ r) ![Dψ r v, Dψ r w]) =
      fun r => Ω r (ψ r) (fun i => (![Dψ r v, Dψ r w] : Fin 2 → W) i) := rfl
  rw [e]
  refine h.congr_deriv ?_
  have u0 : Function.update ![Dψ s v, Dψ s w] 0 (fderiv ℝ (X s) (ψ s) (Dψ s v)) =
      ![fderiv ℝ (X s) (ψ s) (Dψ s v), Dψ s w] := by ext i; fin_cases i <;> rfl
  have u1 : Function.update ![Dψ s v, Dψ s w] 1 (fderiv ℝ (X s) (ψ s) (Dψ s w)) =
      ![Dψ s v, fderiv ℝ (X s) (ψ s) (Dψ s w)] := by ext i; fin_cases i <;> rfl
  simp [lieAlt, Fin.sum_univ_two, u0, u1]
  ring

/-- **Moser.** If `Ω̇ + L_X Ω = 0` along the flow for all `s`, then `ψ_s*Ω_s` is constant in `s`. -/
theorem pullback_const_alt {Ω Ωdot : ℝ → W → W [⋀^Fin 2]→L[ℝ] ℝ} {X : ℝ → W → W}
    {ψ : ℝ → W} {Dψ : ℝ → W →L[ℝ] W}
    (hcurve : ∀ s, HasDerivAt (fun r => Ω r (ψ r))
      (Ωdot s (ψ s) + fderiv ℝ (Ω s) (ψ s) (X s (ψ s))) s)
    (hDψ : ∀ s, HasDerivAt Dψ ((fderiv ℝ (X s) (ψ s)).comp (Dψ s)) s)
    (hmoser : ∀ s a b, Ωdot s (ψ s) ![a, b] + lieAlt (Ω s) (X s) (ψ s) a b = 0) (v w : W) (s : ℝ) :
    Ω s (ψ s) ![Dψ s v, Dψ s w] = Ω 0 (ψ 0) ![Dψ 0 v, Dψ 0 w] := by
  have h : ∀ r, HasDerivAt (fun r => Ω r (ψ r) ![Dψ r v, Dψ r w]) 0 r := fun r => by
    have := hasDerivAt_pullback_alt (hcurve r) (hDψ r) v w
    rwa [hmoser] at this
  exact is_const_of_deriv_eq_zero (fun r => (h r).differentiableAt) (fun r => (h r).deriv) s 0

end CWFlow
