import CWFlow.Extension.SymbolLow

/-!
# Localized weights

A cutoff `ρ` times a symbol evaluated on the jets of `φ`, supported where the jets are admissible,
is a globally smooth compactly supported weight `ρ² · 𝔣(jets of φ)`, with explicit first and second
partial derivatives.
-/

open scoped ContDiff

namespace CWFlow.Extension

variable {n N : ℕ}

/-- The localized weight `ρ² · 𝔣(jets of φ)`. -/
noncomputable def locWeight (ρ : EuclideanSpace ℝ (Fin n) → ℝ) (𝔣 : P n N → ℝ)
    (φ : EuclideanSpace ℝ (Fin n) → ℝ) (x : EuclideanSpace ℝ (Fin n)) : ℝ :=
  ρ x ^ 2 * 𝔣 (jp N φ x)

theorem lw_eventuallyEq_zero {ρ : EuclideanSpace ℝ (Fin n) → ℝ} (𝔣 : P n N → ℝ)
    (φ : EuclideanSpace ℝ (Fin n) → ℝ) {x : EuclideanSpace ℝ (Fin n)} (hx : x ∉ tsupport ρ) :
    locWeight ρ 𝔣 φ =ᶠ[nhds x] fun _ => 0 := by
  have hn : (tsupport ρ)ᶜ ∈ nhds x := (isClosed_tsupport ρ).isOpen_compl.mem_nhds hx
  filter_upwards [hn] with y hy
  have : ρ y = 0 := image_eq_zero_of_notMem_tsupport hy
  simp [locWeight, this]

theorem locWeight_eq_zero {ρ : EuclideanSpace ℝ (Fin n) → ℝ} (𝔣 : P n N → ℝ)
    (φ : EuclideanSpace ℝ (Fin n) → ℝ) {x : EuclideanSpace ℝ (Fin n)} (hx : x ∉ tsupport ρ) :
    locWeight ρ 𝔣 φ x = 0 := by
  have : ρ x = 0 := image_eq_zero_of_notMem_tsupport hx
  simp [locWeight, this]

theorem contDiff_locWeight {A : Set (P n N)} {w ℓ : ℕ} {𝔣 : P n N → ℝ} (h : IsSymbol A w ℓ 𝔣)
    {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ) {ρ : EuclideanSpace ℝ (Fin n) → ℝ}
    (hρ : ContDiff ℝ ∞ ρ) (hsupp : ∀ x ∈ tsupport ρ, jp N φ x ∈ A) :
    ContDiff ℝ ∞ (locWeight ρ 𝔣 φ) := by
  rw [contDiff_iff_contDiffAt]
  intro x
  by_cases hx : x ∈ tsupport ρ
  · have hU := isOpen_jp_preimage h hφ
    have hon : ContDiffOn ℝ ∞ (locWeight ρ 𝔣 φ) {y | jp N φ y ∈ A} :=
      (hρ.pow 2).contDiffOn.mul (contDiffOn_symbol_jp h hφ)
    exact hon.contDiffAt (hU.mem_nhds (hsupp x hx))
  · exact (contDiffAt_const (c := (0 : ℝ))).congr_of_eventuallyEq (lw_eventuallyEq_zero 𝔣 φ hx)

theorem hasCompactSupport_locWeight {ρ : EuclideanSpace ℝ (Fin n) → ℝ} (hc : HasCompactSupport ρ)
    (𝔣 : P n N → ℝ) (φ : EuclideanSpace ℝ (Fin n) → ℝ) : HasCompactSupport (locWeight ρ 𝔣 φ) := by
  refine hc.mono' ?_
  intro x hx
  by_contra hn
  exact hx (locWeight_eq_zero 𝔣 φ hn)

/-- First partial derivatives, at points where the jets are admissible. -/
theorem fderiv_locWeight {A : Set (P n N)} {w ℓ : ℕ} {𝔣 : P n N → ℝ} (h : IsSymbol A w ℓ 𝔣)
    (hN : w + 4 ≤ N) {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ)
    {ρ : EuclideanSpace ℝ (Fin n) → ℝ} (hρ : ContDiff ℝ ∞ ρ) {x : EuclideanSpace ℝ (Fin n)}
    (hx : jp N φ x ∈ A) (e : EuclideanSpace ℝ (Fin n)) :
    fderiv ℝ (locWeight ρ 𝔣 φ) x e =
      2 * ρ x * fderiv ℝ ρ x e * 𝔣 (jp N φ x) + ρ x ^ 2 * totalD 𝔣 (jp N φ x) e := by
  have hρd : HasFDerivAt ρ (fderiv ℝ ρ x) x :=
    (hρ.differentiable (by simp)).differentiableAt.hasFDerivAt
  have hg := hasFDerivAt_symbol_jp h hN hφ hx
  have hH : HasFDerivAt (fun y => ρ y * ρ y * 𝔣 (jp N φ y)) _ x := (hρd.mul hρd).mul hg
  have key : locWeight ρ 𝔣 φ = fun y => ρ y * ρ y * 𝔣 (jp N φ y) := by
    funext y; simp [locWeight, sq]
  rw [key, hH.fderiv]
  simp
  ring

/-- The first partial derivative vanishes off the support of the cutoff. -/
theorem fderiv_locWeight_eq_zero {ρ : EuclideanSpace ℝ (Fin n) → ℝ} (𝔣 : P n N → ℝ)
    (φ : EuclideanSpace ℝ (Fin n) → ℝ) {x : EuclideanSpace ℝ (Fin n)} (hx : x ∉ tsupport ρ) :
    fderiv ℝ (locWeight ρ 𝔣 φ) x = 0 := by
  rw [(lw_eventuallyEq_zero 𝔣 φ hx).fderiv_eq]
  simp

/-- Second partial derivatives, at points where the jets are admissible; `hlow` is the standing bound on the low jets. -/
theorem fderiv_fderiv_locWeight {A : Set (P n N)} {w ℓ : ℕ} {𝔣 : P n N → ℝ} (h : IsSymbol A w ℓ 𝔣)
    (hlow : ∀ m : Fin (N + 1), (m : ℕ) ≤ 3 → ∃ K : ℝ, ∀ p ∈ A, ‖p.2 m‖ ≤ K)
    (hN : w + 5 ≤ N) {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ)
    {ρ : EuclideanSpace ℝ (Fin n) → ℝ} (hρ : ContDiff ℝ ∞ ρ) {x : EuclideanSpace ℝ (Fin n)}
    (hx : jp N φ x ∈ A) (e e' : EuclideanSpace ℝ (Fin n)) :
    fderiv ℝ (fun y => fderiv ℝ (locWeight ρ 𝔣 φ) y e') x e =
      fderiv ℝ (fun y => 2 * ρ y * fderiv ℝ ρ y e') x e * 𝔣 (jp N φ x)
      + 2 * ρ x * fderiv ℝ ρ x e' * totalD 𝔣 (jp N φ x) e
      + 2 * ρ x * fderiv ℝ ρ x e * totalD 𝔣 (jp N φ x) e'
      + ρ x ^ 2 * totalD (fun q => totalD 𝔣 q e') (jp N φ x) e := by
  have hU := isOpen_jp_preimage h hφ
  have hG : IsSymbol A (w + 1) (max ℓ w + 1) (fun q => totalD 𝔣 q e') :=
    (h.totalD (by omega) hlow).clm_comp (ContinuousLinearMap.apply ℝ ℝ e')
  have hρd : HasFDerivAt ρ (fderiv ℝ ρ x) x :=
    (hρ.differentiable (by simp)).differentiableAt.hasFDerivAt
  have hdρ : ContDiff ℝ ∞ (fderiv ℝ ρ) := hρ.fderiv_right (by simp)
  have hF1 : Differentiable ℝ (fun y => 2 * ρ y * fderiv ℝ ρ y e') :=
    ((contDiff_const.mul hρ).mul (hdρ.clm_apply contDiff_const)).differentiable (by simp)
  have hF1x := (hF1 x).hasFDerivAt
  have hg := hasFDerivAt_symbol_jp h (by omega) hφ hx
  have hg2 := hasFDerivAt_symbol_jp hG (by omega) hφ hx
  have hH : HasFDerivAt (fun y => 2 * ρ y * fderiv ℝ ρ y e' * 𝔣 (jp N φ y)
      + ρ y * ρ y * totalD 𝔣 (jp N φ y) e') _ x := (hF1x.mul hg).add ((hρd.mul hρd).mul hg2)
  have hev : (fun y => fderiv ℝ (locWeight ρ 𝔣 φ) y e') =ᶠ[nhds x]
      fun y => 2 * ρ y * fderiv ℝ ρ y e' * 𝔣 (jp N φ y) + ρ y * ρ y * totalD 𝔣 (jp N φ y) e' := by
    filter_upwards [hU.mem_nhds hx] with y hy
    rw [fderiv_locWeight h (by omega) hφ hρ hy e']
    ring
  rw [hev.fderiv_eq, hH.fderiv]
  simp
  ring

/-- info: 'CWFlow.Extension.contDiff_locWeight' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms contDiff_locWeight

/-- info: 'CWFlow.Extension.fderiv_locWeight' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms fderiv_locWeight

/-- info: 'CWFlow.Extension.fderiv_fderiv_locWeight' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms fderiv_fderiv_locWeight

end CWFlow.Extension
