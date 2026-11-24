import CWFlow.Weinstein.Bridge

/-!
# The Moser vector field

On a finite-dimensional normed space `W`, for 1-forms `λ, β` (fields of functionals) put
`Ω_t(q) = dλ(q) + t·dβ(q)` (`omegaT`) and `X_t(q) = −(ι_· Ω_t(q))⁻¹ β(q)` (`moserField`). Where `ι_· Ω_t(q)` is
invertible, `Ω_t(q)(X_t(q), b) = −β q b` (`moserField_spec`); invertibility is open in `(t, q)`, and `X_t(q)` is
smooth in `(t, q)` there. At a point where `dβ = 0` and `ι_· dλ` is injective, `ι_· Ω_t` is invertible for
every `t` (`iotaT_invertible`). `omegaT_eq` identifies `Ω_t = d(λ + tβ)`, the form required by
`moser_condition_exact`.
-/

open ContinuousAlternatingMap Filter

open scoped Topology ContDiff

noncomputable section

namespace CWFlow

variable {W : Type*} [NormedAddCommGroup W] [NormedSpace ℝ W] [FiniteDimensional ℝ W]

instance finiteDimensional_alt1 : FiniteDimensional ℝ (W [⋀^Fin 1]→L[ℝ] ℝ) :=
  LinearEquiv.finiteDimensional
    (ofSubsingletonLIE (𝕜 := ℝ) (E := W) (F := ℝ) (0 : Fin 1)).toLinearEquiv

theorem finrank_alt1 : Module.finrank ℝ (W [⋀^Fin 1]→L[ℝ] ℝ) = Module.finrank ℝ W := by
  rw [← (ofSubsingletonLIE (𝕜 := ℝ) (E := W) (F := ℝ) (0 : Fin 1)).toLinearEquiv.finrank_eq,
    ← (LinearMap.toContinuousLinearMap : Module.Dual ℝ W ≃ₗ[ℝ] (W →L[ℝ] ℝ)).finrank_eq,
    Subspace.dual_finrank_eq]

/-- A 2-form with injective `v ↦ ι_v ω` gives an isomorphism `W ≃L W*`. -/
def curryEquiv (σ : W [⋀^Fin (1 + 1)]→L[ℝ] ℝ) (hinj : Function.Injective σ.curryLeft) :
    W ≃L[ℝ] (W [⋀^Fin 1]→L[ℝ] ℝ) :=
  (LinearMap.linearEquivOfInjective (σ.curryLeft : W →ₗ[ℝ] _) hinj
    (finrank_alt1 (W := W)).symm).toContinuousLinearEquiv

theorem coe_curryEquiv (σ : W [⋀^Fin (1 + 1)]→L[ℝ] ℝ) (hinj : Function.Injective σ.curryLeft) :
    (curryEquiv σ hinj : W →L[ℝ] (W [⋀^Fin 1]→L[ℝ] ℝ)) =
      curryLeftLI (𝕜 := ℝ) (E := W) (F := ℝ) (n := 1) σ :=
  ContinuousLinearMap.ext fun _ => rfl

/-- `Ω_t(q) = dλ(q) + t·dβ(q)`. -/
def omegaT (lam β : W → W →L[ℝ] ℝ) (t : ℝ) (q : W) : W [⋀^Fin (1 + 1)]→L[ℝ] ℝ :=
  extDeriv (toAlt1 lam) q + t • extDeriv (toAlt1 β) q

/-- `v ↦ ι_v Ω_t(q)`. -/
def iotaT (lam β : W → W →L[ℝ] ℝ) (t : ℝ) (q : W) : W →L[ℝ] (W [⋀^Fin 1]→L[ℝ] ℝ) :=
  curryLeftLI (𝕜 := ℝ) (E := W) (F := ℝ) (n := 1) (omegaT lam β t q)

/-- The Moser field `X_t(q) = −(ι_· Ω_t(q))⁻¹ β(q)`. -/
def moserField (lam β : W → W →L[ℝ] ℝ) (t : ℝ) (q : W) : W :=
  -(ContinuousLinearMap.inverse (iotaT lam β t q)) (toAlt1 β q)

theorem toAlt1_add_smul (lam β : W → W →L[ℝ] ℝ) (t : ℝ) :
    toAlt1 (fun q => lam q + t • β q) = toAlt1 lam + t • toAlt1 β := by
  funext q
  ext v
  simp

/-- `Ω_t = d(λ + tβ)`. -/
theorem omegaT_eq {lam β : W → W →L[ℝ] ℝ} {q : W} (hlam : DifferentiableAt ℝ lam q)
    (hβ : DifferentiableAt ℝ β q) (t : ℝ) :
    omegaT lam β t q = extDeriv (toAlt1 (fun q => lam q + t • β q)) q := by
  rw [toAlt1_add_smul, extDeriv_add (differentiableAt_toAlt1 hlam)
    ((differentiableAt_toAlt1 hβ).const_smul t), extDeriv_smul, omegaT]

/-- Where `ι_· Ω_t(q)` is invertible, `Ω_t(q)(X_t(q), b) = −β q b`. -/
theorem moserField_spec {lam β : W → W →L[ℝ] ℝ} {t : ℝ} {q : W}
    (hinv : ∃ e : W ≃L[ℝ] (W [⋀^Fin 1]→L[ℝ] ℝ), (e : _ →L[ℝ] _) = iotaT lam β t q)
    (b : W) : omegaT lam β t q ![moserField lam β t q, b] = -β q b := by
  obtain ⟨e, he⟩ := hinv
  have h1 : omegaT lam β t q ![moserField lam β t q, b] =
      iotaT lam β t q (moserField lam β t q) ![b] := rfl
  rw [h1, moserField, ← he, ContinuousLinearMap.inverse_equiv, map_neg]
  simp

/-- At a point where `dβ = 0` and `ι_· dλ` is injective, `ι_· Ω_t` is invertible for every `t`. -/
theorem iotaT_invertible {lam β : W → W →L[ℝ] ℝ} {q : W} (t : ℝ)
    (h0 : extDeriv (toAlt1 β) q = 0)
    (hinj : Function.Injective (extDeriv (toAlt1 lam) q).curryLeft) :
    ∃ e : W ≃L[ℝ] (W [⋀^Fin 1]→L[ℝ] ℝ), (e : _ →L[ℝ] _) = iotaT lam β t q :=
  ⟨curryEquiv _ hinj, by rw [coe_curryEquiv, iotaT, omegaT, h0, smul_zero, add_zero]⟩

theorem contDiffAt_extDeriv_toAlt1 {β : W → W →L[ℝ] ℝ} {q₀ : W} (hβ : ContDiffAt ℝ ∞ β q₀) :
    ContDiffAt ℝ ∞ (extDeriv (toAlt1 β)) q₀ := by
  have hf := (contDiffAt_toAlt1 hβ).fderiv_right (m := ∞) le_rfl
  exact ((alternatizeUncurryFinCLM ℝ W ℝ).contDiff.contDiffAt).comp q₀ hf

theorem contDiffAt_omegaT {lam β : W → W →L[ℝ] ℝ} {t₀ : ℝ} {q₀ : W}
    (hlam : ContDiffAt ℝ ∞ lam q₀) (hβ : ContDiffAt ℝ ∞ β q₀) :
    ContDiffAt ℝ ∞ (fun p : ℝ × W => omegaT lam β p.1 p.2) (t₀, q₀) := by
  have h1 : ContDiffAt ℝ ∞ (fun p : ℝ × W => extDeriv (toAlt1 lam) p.2) (t₀, q₀) := by
    exact ContDiffAt.comp (g := extDeriv (toAlt1 lam)) (f := Prod.snd) (t₀, q₀)
      (contDiffAt_extDeriv_toAlt1 hlam) contDiffAt_snd
  have h2 : ContDiffAt ℝ ∞ (fun p : ℝ × W => extDeriv (toAlt1 β) p.2) (t₀, q₀) := by
    exact ContDiffAt.comp (g := extDeriv (toAlt1 β)) (f := Prod.snd) (t₀, q₀)
      (contDiffAt_extDeriv_toAlt1 hβ) contDiffAt_snd
  exact h1.add (contDiffAt_fst.smul h2)

theorem contDiffAt_iotaT {lam β : W → W →L[ℝ] ℝ} {t₀ : ℝ} {q₀ : W}
    (hlam : ContDiffAt ℝ ∞ lam q₀) (hβ : ContDiffAt ℝ ∞ β q₀) :
    ContDiffAt ℝ ∞ (fun p : ℝ × W => iotaT lam β p.1 p.2) (t₀, q₀) := by
  have hc := LinearIsometry.contDiff (𝕜 := ℝ) (E := W [⋀^Fin (1 + 1)]→L[ℝ] ℝ)
    (F := W →L[ℝ] W [⋀^Fin 1]→L[ℝ] ℝ) (n := ∞) (curryLeftLI (𝕜 := ℝ) (E := W) (F := ℝ) (n := 1))
  exact ContDiffAt.comp (t₀, q₀) hc.contDiffAt (contDiffAt_omegaT hlam hβ)

/-- The Moser field is smooth in `(t, q)` near a point where `ι_· Ω_t(q)` is invertible. -/
theorem contDiffAt_moserField {lam β : W → W →L[ℝ] ℝ} {t₀ : ℝ} {q₀ : W}
    (hlam : ContDiffAt ℝ ∞ lam q₀) (hβ : ContDiffAt ℝ ∞ β q₀)
    (hinv : ∃ e : W ≃L[ℝ] (W [⋀^Fin 1]→L[ℝ] ℝ), (e : _ →L[ℝ] _) = iotaT lam β t₀ q₀) :
    ContDiffAt ℝ ∞ (fun p : ℝ × W => moserField lam β p.1 p.2) (t₀, q₀) := by
  obtain ⟨e, he⟩ := hinv
  have hI : ContDiffAt ℝ ∞ (fun p : ℝ × W => ContinuousLinearMap.inverse (iotaT lam β p.1 p.2))
      (t₀, q₀) := by
    have h0 := contDiffAt_map_inverse (𝕜 := ℝ) (n := ∞) e
    rw [he] at h0
    exact h0.comp (t₀, q₀) (contDiffAt_iotaT hlam hβ)
  have hb : ContDiffAt ℝ ∞ (fun p : ℝ × W => toAlt1 β p.2) (t₀, q₀) := by
    exact ContDiffAt.comp (g := toAlt1 β) (f := Prod.snd) (t₀, q₀) (contDiffAt_toAlt1 hβ)
      contDiffAt_snd
  exact (hI.clm_apply hb).neg

/-- Invertibility of `ι_· Ω_t(q)` is an open condition in `(t, q)`. -/
theorem eventually_iotaT_invertible {lam β : W → W →L[ℝ] ℝ} {t₀ : ℝ} {q₀ : W}
    (hlam : ContDiffAt ℝ ∞ lam q₀) (hβ : ContDiffAt ℝ ∞ β q₀)
    (hinv : ∃ e : W ≃L[ℝ] (W [⋀^Fin 1]→L[ℝ] ℝ), (e : _ →L[ℝ] _) = iotaT lam β t₀ q₀) :
    ∀ᶠ p : ℝ × W in 𝓝 (t₀, q₀),
      ∃ e : W ≃L[ℝ] (W [⋀^Fin 1]→L[ℝ] ℝ), (e : _ →L[ℝ] _) = iotaT lam β p.1 p.2 := by
  have hΛ := (contDiffAt_iotaT (t₀ := t₀) hlam hβ).continuousAt
  obtain ⟨e, he⟩ := hinv
  have hmem : iotaT lam β t₀ q₀ ∈ Set.range ((↑) : (W ≃L[ℝ] (W [⋀^Fin 1]→L[ℝ] ℝ)) → _) := ⟨e, he⟩
  filter_upwards [hΛ.preimage_mem_nhds (ContinuousLinearEquiv.isOpen.mem_nhds hmem)] with p hp
  exact hp

end CWFlow
