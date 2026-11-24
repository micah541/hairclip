import DifferentialGeometry.Analysis.ODE.TimeDependentFlow.DiffeomorphismFamily.CenteredInterval

/-!
# Reparametrization flows on a closed manifold (Chen–Warren, Claim 3.3 ingredients)

* `hasDerivAt_comp_curve`: the chain rule for `s ↦ F (s, χ s)` with `F : ℝ × M → V` and a curve
  `χ` in `M` with velocity `v` at `t`:
  `d/ds F (s, χ s) = ∂ₜF (t, χ t) + d(F t)_{χ t} v`.
* `exists_reparam_flow`: a smooth time-dependent vector field `X` on a closed manifold has a flow
  `χ t` of diffeomorphisms with `χ 0 = id`, jointly smooth (with jointly smooth inverses) on
  `Ioo (-1) (T + 1) × M ⊇ Icc 0 T × M`, and `∂ₜ χ t x = X t (χ t x)`. A thin wrapper around
  DifferentialGeometry's `exists_diffeomorph_flow_on_centered_interval`.
* `exists_reparam_flow_hasDerivAt`: the two combined,
  `d/ds F (s, χ s x) = ∂ₜF (t, χ t x) + d(F t)_{χ t x} (X t (χ t x))`.
-/

open Set Function Bundle

open scoped Manifold Topology ContDiff

noncomputable section

namespace CWFlow

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
  {M : Type*} [TopologicalSpace M] [ChartedSpace E M]
  {V : Type*} [NormedAddCommGroup V] [NormedSpace ℝ V]

/-- Chain rule for `s ↦ F (s, χ s)`: the time derivative plus the spatial derivative applied to the
velocity `v` of `χ` at `t`. -/
theorem hasDerivAt_comp_curve {F : ℝ × M → V} {χ : ℝ → M} {t : ℝ}
    (hF : MDifferentiableAt (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, V) F (t, χ t))
    {v : TangentSpace 𝓘(ℝ, E) (χ t)}
    (hχ : HasMFDerivAt 𝓘(ℝ, ℝ) 𝓘(ℝ, E) χ t ((1 : ℝ →L[ℝ] ℝ).smulRight v)) :
    HasDerivAt (fun s => F (s, χ s))
      (mfderiv 𝓘(ℝ, ℝ) 𝓘(ℝ, V) (fun s => F (s, χ t)) t (1 : ℝ) +
        mfderiv 𝓘(ℝ, E) 𝓘(ℝ, V) (fun y => F (t, y)) (χ t) v) t := by
  have hFd := hF.hasMFDerivAt
  -- the full derivative along the curve
  have hcurve := (hasMFDerivAt_id (I := 𝓘(ℝ, ℝ)) t).prodMk hχ
  have hG := hFd.comp t hcurve
  -- re-type the derivative with the normed-space instances of `V` (not those of `TangentSpace`)
  have hGf' : HasFDerivAt (fun s => F (s, χ s)) (_ : ℝ →L[ℝ] V) t :=
    hasMFDerivAt_iff_hasFDerivAt.mp hG
  have hGf := hGf'.hasDerivAt
  -- the partial derivative in time
  have htime := hFd.comp t ((hasMFDerivAt_id (I := 𝓘(ℝ, ℝ)) t).prodMk
    (hasMFDerivAt_const (I := 𝓘(ℝ, ℝ)) (I' := 𝓘(ℝ, E)) (χ t) t))
  have htime' : mfderiv 𝓘(ℝ, ℝ) 𝓘(ℝ, V) (fun s => F (s, χ t)) t = _ := htime.mfderiv
  -- the partial derivative in space
  have hspace := hFd.comp (χ t) ((hasMFDerivAt_const (I := 𝓘(ℝ, E)) (I' := 𝓘(ℝ, ℝ)) t (χ t)).prodMk
    (hasMFDerivAt_id (I := 𝓘(ℝ, E)) (χ t)))
  have hspace' : mfderiv 𝓘(ℝ, E) 𝓘(ℝ, V) (fun y => F (t, y)) (χ t) = _ := hspace.mfderiv
  refine hGf.congr_deriv ?_
  rw [htime', hspace']
  -- the remaining identity `D (1, v) = D (1, 0) + D (0, v)`, stated with the normed-space instances
  have key : ∀ D : ℝ × E →L[ℝ] V,
      D ((1 : ℝ), ((1 : ℝ →L[ℝ] ℝ) (1 : ℝ)) • (v : E)) = D ((1 : ℝ), 0) + D (0, (v : E)) := by
    intro D
    rw [← map_add]
    congr 1
    refine Prod.ext ?_ ?_
    · change (1 : ℝ) = 1 + 0
      rw [add_zero]
    · change (1 : ℝ) • (v : E) = 0 + (v : E)
      rw [zero_add, one_smul]
  exact key (mfderiv (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, V) F (t, χ t))

variable [FiniteDimensional ℝ E] [IsManifold 𝓘(ℝ, E) ∞ M] [CompactSpace M] [T2Space M]

/-- The reparametrization flow of a smooth time-dependent vector field `X` on a closed manifold:
a family of diffeomorphisms `χ t`, `χ 0 = id`, jointly smooth (with jointly smooth inverses) for
`t ∈ Ioo (-1) (T + 1)` (an open interval containing `Icc 0 T`), with `∂ₜ χ t x = X t (χ t x)`. -/
theorem exists_reparam_flow (X : ℝ → ∀ x : M, TangentSpace 𝓘(ℝ, E) x)
    (hX : ContMDiff (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) (𝓘(ℝ, E).prod 𝓘(ℝ, E)) ∞
      (fun q : ℝ × M => (TotalSpace.mk' E q.2 (X q.1 q.2) : TangentBundle 𝓘(ℝ, E) M)))
    {T : ℝ} (hT : 0 ≤ T) :
    ∃ χ : ℝ → (M ≃ₘ⟮𝓘(ℝ, E), 𝓘(ℝ, E)⟯ M),
      (∀ x, χ 0 x = x) ∧
      ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, E) ∞ (fun q : ℝ × M => χ q.1 q.2)
        (Ioo (-1) (T + 1) ×ˢ univ) ∧
      ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, E) ∞ (fun q : ℝ × M => (χ q.1).symm q.2)
        (Ioo (-1) (T + 1) ×ˢ univ) ∧
      (∀ t ∈ Ioo (-1) (T + 1), ∀ x,
        HasMFDerivAt 𝓘(ℝ, ℝ) 𝓘(ℝ, E) (fun s => χ s x) t
          ((1 : ℝ →L[ℝ] ℝ).smulRight (X t (χ t x)))) := by
  obtain ⟨Φ, h0, hsm, hsm', hode⟩ :=
    DifferentialGeometry.Analysis.ODE.exists_diffeomorph_flow_on_centered_interval
      (I := 𝓘(ℝ, E)) X hX 0 (T + 1) (by linarith)
  have hsub : Ioo (-1) (T + 1) ⊆ Ioo (0 - (T + 1)) (0 + (T + 1)) :=
    Ioo_subset_Ioo (by linarith) (by linarith)
  exact ⟨Φ, h0, hsm.mono (prod_mono hsub subset_rfl), hsm'.mono (prod_mono hsub subset_rfl),
    fun t ht x => hode t (hsub ht) x⟩

theorem Icc_subset_Ioo_reparam {T : ℝ} : Icc 0 T ⊆ Ioo (-1) (T + 1) :=
  Icc_subset_Ioo (by norm_num) (by linarith)

/-- The reparametrization flow together with the chain rule along it: for `F : ℝ × M → V` smooth,
`d/ds F (s, χ s x) = ∂ₜF (t, χ t x) + d(F t)_{χ t x} (X t (χ t x))`. -/
theorem exists_reparam_flow_hasDerivAt (X : ℝ → ∀ x : M, TangentSpace 𝓘(ℝ, E) x)
    (hX : ContMDiff (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) (𝓘(ℝ, E).prod 𝓘(ℝ, E)) ∞
      (fun q : ℝ × M => (TotalSpace.mk' E q.2 (X q.1 q.2) : TangentBundle 𝓘(ℝ, E) M)))
    {T : ℝ} (hT : 0 ≤ T) (F : ℝ × M → V)
    (hF : ContMDiff (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, V) ∞ F) :
    ∃ χ : ℝ → (M ≃ₘ⟮𝓘(ℝ, E), 𝓘(ℝ, E)⟯ M),
      (∀ x, χ 0 x = x) ∧
      ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, E) ∞ (fun q : ℝ × M => χ q.1 q.2)
        (Ioo (-1) (T + 1) ×ˢ univ) ∧
      ContMDiffOn (𝓘(ℝ, ℝ).prod 𝓘(ℝ, E)) 𝓘(ℝ, E) ∞ (fun q : ℝ × M => (χ q.1).symm q.2)
        (Ioo (-1) (T + 1) ×ˢ univ) ∧
      (∀ t ∈ Ioo (-1) (T + 1), ∀ x,
        HasMFDerivAt 𝓘(ℝ, ℝ) 𝓘(ℝ, E) (fun s => χ s x) t
          ((1 : ℝ →L[ℝ] ℝ).smulRight (X t (χ t x)))) ∧
      (∀ t ∈ Ioo (-1) (T + 1), ∀ x,
        HasDerivAt (fun s => F (s, χ s x))
          (mfderiv 𝓘(ℝ, ℝ) 𝓘(ℝ, V) (fun s => F (s, χ t x)) t (1 : ℝ) +
            mfderiv 𝓘(ℝ, E) 𝓘(ℝ, V) (fun y => F (t, y)) (χ t x) (X t (χ t x))) t) := by
  obtain ⟨χ, h0, hsm, hsm', hode⟩ := exists_reparam_flow X hX hT
  exact ⟨χ, h0, hsm, hsm', hode, fun t ht x =>
    hasDerivAt_comp_curve (χ := fun s => χ s x) (hF.mdifferentiableAt (by simp)) (hode t ht x)⟩

end CWFlow
