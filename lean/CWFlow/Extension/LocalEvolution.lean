import CWFlow.Extension.LeadingMain
import CWFlow.Extension.LeadingLower

/-!
# The local evolution inequality

Assembly of the pointwise time derivative of the energy density (`EvolPointwise`), the main term estimate
(`LeadingMain`) and the lower order estimate (`LeadingLower`) into a differential inequality for the cut-off
energy of a graph function solving the scalar equation on the ball of radius `3`.
-/

open CWFlow.Interp Metric MeasureTheory
open scoped ContDiff

namespace CWFlow.Extension

variable {n N : ℕ}

/-! ### Additivity of `cTop` -/

theorem le_cTop_add (k : ℕ) (q : P n N) (z₁ z₂ : Z n (k + 2)) :
    cTop k q (z₁ + z₂) = cTop k q z₁ + cTop k q z₂ := by
  unfold cTop
  rw [← map_add]
  congr 1
  refine Prod.ext (by simp) ?_
  funext m
  by_cases h : (m : ℕ) = k + 2
  · obtain ⟨m, hm⟩ := m
    simp only at h
    subst h
    simp
  · simp [h]

/-! ### Decomposition of the main expression -/

theorem le_decomp {k : ℕ} (hN : k + 7 ≤ N) {K : ℝ} {q : P n N} (hq : q ∈ admSet n N K) :
    2 * qSym (k + 3) q (jetSym potSym (k + 3) q + cTop k q (jetSym potSym (k + 2) q)) (tSym k q)
        * volFun (low23 q).1
      = 2 * ∑ I : Fin (k + 3) → Fin n, topPot (k + 3) q I * wSym k q I
        + 2 * ∑ I : Fin (k + 3) → Fin n,
          (rPot (k + 2) q I
            + Dz (fun r : P n N => jetSym potSym (k + 2) r - topPot (k + 2) r) q I
            + cTop k q (topPot (k + 2) q) I
            + cTop k q (jetSym potSym (k + 2) q - topPot (k + 2) q) I) * wSym k q I := by
  have hA := isOpen_admSet (n := n) (N := N) K
  have hlow := admSet_low (n := n) (N := N) K
  have hnhds : admSet n N K ∈ nhds q := hA.mem_nhds hq
  have hd1 : DifferentiableAt ℝ (topPot (k + 2) : P n N → Z n (k + 2)) q :=
    ((isSymbol_topPot hA hlow (k + 2)).smooth.contDiffAt hnhds).differentiableAt (by simp)
  have hd2 : DifferentiableAt ℝ
      (fun r : P n N => jetSym potSym (k + 2) r - topPot (k + 2) r) q :=
    ((isSymbol_jetSym_potSym_sub hA hlow (k + 2) (by omega)).smooth.contDiffAt
      hnhds).differentiableAt (by simp)
  have hfun : (jetSym potSym (k + 2) : P n N → Z n (k + 2)) =
      fun r => topPot (k + 2) r + (jetSym potSym (k + 2) r - topPot (k + 2) r) := by
    funext r; abel
  have hD1 : jetSym potSym (k + 3) q = topPot (k + 3) q + rPot (k + 2) q +
      Dz (fun r : P n N => jetSym potSym (k + 2) r - topPot (k + 2) r) q := by
    have e1 : jetSym potSym (k + 3) q = Dz (jetSym potSym (k + 2)) q := rfl
    have e2 := Dz_topPot_sub (n := n) (N := N) (k + 2) (by omega) q hd1
    rw [e1, hfun, Dz_add hd1 hd2]
    rw [← hfun]
    have : Dz (topPot (k + 2)) q = topPot (k + 3) q + rPot (k + 2) q := by
      rw [← e2]; abel
    rw [this]
  have hC : cTop k q (jetSym potSym (k + 2) q) = cTop k q (topPot (k + 2) q) +
      cTop k q (jetSym potSym (k + 2) q - topPot (k + 2) q) := by
    rw [← le_cTop_add]; congr 1; abel
  set S2 : Z n (k + 3) := rPot (k + 2) q
    + Dz (fun r : P n N => jetSym potSym (k + 2) r - topPot (k + 2) r) q
    + cTop k q (topPot (k + 2) q)
    + cTop k q (jetSym potSym (k + 2) q - topPot (k + 2) q) with hS2
  have hS : jetSym potSym (k + 3) q + cTop k q (jetSym potSym (k + 2) q) =
      topPot (k + 3) q + S2 := by
    rw [hD1, hC, hS2]; abel
  have h1 := qSym_mul_vol_eq (n := n) (N := N) k q (topPot (k + 3) q)
  have h2 := qSym_mul_vol_eq (n := n) (N := N) k q S2
  rw [hS, map_add, add_apply,
    show ∀ a b V : ℝ, 2 * (a + b) * V = 2 * (a * V) + 2 * (b * V) from fun a b V => by ring, h1, h2]
  have h3 : ∀ I : Fin (k + 3) → Fin n, S2 I = rPot (k + 2) q I
      + Dz (fun r : P n N => jetSym potSym (k + 2) r - topPot (k + 2) r) q I
      + cTop k q (topPot (k + 2) q) I
      + cTop k q (jetSym potSym (k + 2) q - topPot (k + 2) q) I := fun I => by
    rw [hS2]; simp
  simp only [h3]

/-! ### Differentiation under the integral -/

theorem le_tsupport_sq (ρ : EuclideanSpace ℝ (Fin n) → ℝ) :
    tsupport (fun x => ρ x ^ 2) ⊆ tsupport ρ := by
  apply closure_mono
  intro x hx h0
  exact hx (by simp [h0])

theorem le_compactSupport_sq {ρ : EuclideanSpace ℝ (Fin n) → ℝ} {Λ : ℝ} (hρ : IsCutoff ρ Λ) :
    HasCompactSupport (fun x => ρ x ^ 2) :=
  IsCompact.of_isClosed_subset (ibb_hasCompactSupport hρ) (isClosed_tsupport _)
    (le_tsupport_sq ρ)

theorem le_hasDerivAt_integral {k : ℕ} (hN : k + 3 ≤ N) {K : ℝ} {T : Set ℝ}
    {Φ : ℝ × EuclideanSpace ℝ (Fin n) → ℝ} (h : GraphEvol N K 3 T Φ)
    {ρ : EuclideanSpace ℝ (Fin n) → ℝ} {Λ : ℝ} (hρ : IsCutoff ρ Λ) {t : ℝ} (ht : t ∈ T) :
    HasDerivAt (fun s : ℝ => ∫ x, ρ x ^ 2 * eSym k (jp N (fun y => Φ (s, y)) x))
      (∫ x, ρ x ^ 2 * fderiv ℝ
        (fun q : ℝ × EuclideanSpace ℝ (Fin n) => eSym k (jp N (fun y => Φ (q.1, y)) q.2))
        (t, x) (1, 0)) t := by
  have hA := isOpen_admSet (n := n) (N := N) K
  have hlow := admSet_low (n := n) (N := N) K
  have hO := isOpen_symbol_slice_dom (N := N) hA h.isOpen h.smooth
  have hF := contDiffOn_symbol_slice (isSymbol_eSym hA hlow k hN) h.isOpen h.smooth
  have := hasDerivAt_integral_mul_slice hO hF (ρ := fun x => ρ x ^ 2)
    (hρ.smooth.continuous.pow 2) (le_compactSupport_sq hρ) (t := t) (fun x hx => by
      have hx' : x ∈ tsupport ρ := le_tsupport_sq ρ hx
      exact ⟨⟨ht, Set.mem_univ _⟩, h.adm t ht x (ibb_ball_sub x (hρ.supp hx'))⟩)
  exact this

theorem le_integrable_deriv {k : ℕ} (hN : k + 3 ≤ N) {K : ℝ} {T : Set ℝ}
    {Φ : ℝ × EuclideanSpace ℝ (Fin n) → ℝ} (h : GraphEvol N K 3 T Φ)
    {ρ : EuclideanSpace ℝ (Fin n) → ℝ} {Λ : ℝ} (hρ : IsCutoff ρ Λ) {t : ℝ} (ht : t ∈ T) :
    Integrable (fun x => ρ x ^ 2 * fderiv ℝ
        (fun q : ℝ × EuclideanSpace ℝ (Fin n) => eSym k (jp N (fun y => Φ (q.1, y)) q.2))
        (t, x) (1, 0)) := by
  have hA := isOpen_admSet (n := n) (N := N) K
  have hlow := admSet_low (n := n) (N := N) K
  have hO := isOpen_symbol_slice_dom (N := N) hA h.isOpen h.smooth
  have hF := contDiffOn_symbol_slice (isSymbol_eSym hA hlow k hN) h.isOpen h.smooth
  have hGc := ((hF.fderiv_of_isOpen hO (m := ∞) (by simp)).clm_apply
    (contDiffOn_const (c := ((1 : ℝ), (0 : EuclideanSpace ℝ (Fin n)))))).continuousOn
  have hc := slc_continuous_mul hO hGc (ρ := fun x => ρ x ^ 2) (hρ.smooth.continuous.pow 2)
    (s := t) (fun x hx => by
      have hx' : x ∈ tsupport ρ := le_tsupport_sq ρ hx
      exact ⟨⟨ht, Set.mem_univ _⟩, h.adm t ht x (ibb_ball_sub x (hρ.supp hx'))⟩)
  exact hc.integrable_of_hasCompactSupport (le_compactSupport_sq hρ).mul_right

theorem le_fderiv_eq {k : ℕ} (hN : k + 3 ≤ N) {K : ℝ} {T : Set ℝ}
    {Φ : ℝ × EuclideanSpace ℝ (Fin n) → ℝ} (h : GraphEvol N K 3 T Φ) {t : ℝ} (ht : t ∈ T)
    {x : EuclideanSpace ℝ (Fin n)} (hx : x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 3) {E : ℝ}
    (hE : HasDerivAt (fun s : ℝ => eSym k (jp N (fun y => Φ (s, y)) x)) E t) :
    fderiv ℝ (fun q : ℝ × EuclideanSpace ℝ (Fin n) => eSym k (jp N (fun y => Φ (q.1, y)) q.2))
      (t, x) (1, 0) = E := by
  have hA := isOpen_admSet (n := n) (N := N) K
  have hlow := admSet_low (n := n) (N := N) K
  have hO := isOpen_symbol_slice_dom (N := N) hA h.isOpen h.smooth
  have hF := contDiffOn_symbol_slice (isSymbol_eSym hA hlow k hN) h.isOpen h.smooth
  have hmem : (t, x) ∈ {q : ℝ × EuclideanSpace ℝ (Fin n) |
      q ∈ T ×ˢ (Set.univ : Set (EuclideanSpace ℝ (Fin n))) ∧
        jp N (fun y => Φ (q.1, y)) q.2 ∈ admSet n N K} :=
    ⟨⟨ht, Set.mem_univ _⟩, h.adm t ht x hx⟩
  have hd : DifferentiableAt ℝ (fun q : ℝ × EuclideanSpace ℝ (Fin n) =>
      eSym k (jp N (fun y => Φ (q.1, y)) q.2)) (t, x) :=
    (hF.contDiffAt (hO.mem_nhds hmem)).differentiableAt (by simp)
  have hc : HasDerivAt (fun r : ℝ => (r, x)) ((1 : ℝ), (0 : EuclideanSpace ℝ (Fin n))) t :=
    (hasDerivAt_id t).prodMk (hasDerivAt_const t x)
  have h3 := hd.hasFDerivAt.comp_hasDerivAt t hc
  exact h3.unique hE

/-! ### Integrability of the two pieces -/

theorem le_integrable_P {k : ℕ} (hN : k + 7 ≤ N) {K : ℝ} {ρ : EuclideanSpace ℝ (Fin n) → ℝ}
    {Λ : ℝ} (hρ : IsCutoff ρ Λ) {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ)
    (hsupp : ∀ x ∈ tsupport ρ, jp N φ x ∈ admSet n N K) :
    Integrable (fun x => ρ x ^ 2 * (2 * ∑ I : Fin (k + 3) → Fin n,
      topPot (k + 3) (jp N φ x) I * wSym k (jp N φ x) I)) := by
  have hA := isOpen_admSet (n := n) (N := N) K
  have hlow := admSet_low (n := n) (N := N) K
  have s := ll_isSymbol_pair hA hlow k (by omega) (isSymbol_topPot hA hlow (k + 3))
  exact (contDiff_locWeight s hφ hρ.smooth hsupp).continuous.integrable_of_hasCompactSupport
    (hasCompactSupport_locWeight (ibb_hasCompactSupport hρ) _ φ)

theorem le_integrable_Q {k : ℕ} (hk : 1 ≤ k) (hN : k + 7 ≤ N) {K : ℝ}
    {ρ : EuclideanSpace ℝ (Fin n) → ℝ} {Λ : ℝ} (hρ : IsCutoff ρ Λ)
    {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ)
    (hsupp : ∀ x ∈ tsupport ρ, jp N φ x ∈ admSet n N K) :
    Integrable (fun x => ρ x ^ 2 * (2 * ∑ I : Fin (k + 3) → Fin n,
          (rPot (k + 2) (jp N φ x) I
            + Dz (fun q : P n N => jetSym potSym (k + 2) q - topPot (k + 2) q) (jp N φ x) I
            + cTop k (jp N φ x) (topPot (k + 2) (jp N φ x)) I
            + cTop k (jp N φ x) (jetSym potSym (k + 2) (jp N φ x) - topPot (k + 2) (jp N φ x)) I)
          * wSym k (jp N φ x) I)) := by
  have hA := isOpen_admSet (n := n) (N := N) K
  have hlow := admSet_low (n := n) (N := N) K
  have hRm := ll_isSymbol_Rm hA hlow k (by omega : k + 6 ≤ N)
  have s2 := ll_isSymbol_pair hA hlow k (by omega) (isSymbol_rPot hA hlow (k + 2))
  have s3 := ll_isSymbol_pair hA hlow k (by omega) (hRm.Dz (by omega) hlow)
  have s4 := ll_isSymbol_pair hA hlow k (by omega)
    (ll_isSymbol_cTop hA hlow k hk (by omega) (isSymbol_topPot hA hlow (k + 2)))
  have s5 := ll_isSymbol_pair hA hlow k (by omega)
    (ll_isSymbol_cTop hA hlow k hk (by omega) hRm)
  have hI : ∀ {w ℓ : ℕ} {𝔣 : P n N → ℝ}, IsSymbol (admSet n N K) w ℓ 𝔣 →
      Integrable (fun x => ρ x ^ 2 * 𝔣 (jp N φ x)) := fun hs =>
    (contDiff_locWeight hs hφ hρ.smooth hsupp).continuous.integrable_of_hasCompactSupport
      (hasCompactSupport_locWeight (ibb_hasCompactSupport hρ) _ φ)
  have i2 := hI s2
  have i3 : Integrable (fun x => ρ x ^ 2 * (2 * ∑ I : Fin (k + 3) → Fin n,
          Dz (fun q : P n N => jetSym potSym (k + 2) q - topPot (k + 2) q) (jp N φ x) I *
            wSym k (jp N φ x) I)) := hI s3
  have i4 := hI s4
  have i5 : Integrable (fun x => ρ x ^ 2 * (2 * ∑ I : Fin (k + 3) → Fin n,
          cTop k (jp N φ x) (jetSym potSym (k + 2) (jp N φ x) - topPot (k + 2) (jp N φ x)) I *
            wSym k (jp N φ x) I)) := hI s5
  have e : (fun x => ρ x ^ 2 * (2 * ∑ I : Fin (k + 3) → Fin n,
          (rPot (k + 2) (jp N φ x) I
            + Dz (fun q : P n N => jetSym potSym (k + 2) q - topPot (k + 2) q) (jp N φ x) I
            + cTop k (jp N φ x) (topPot (k + 2) (jp N φ x)) I
            + cTop k (jp N φ x) (jetSym potSym (k + 2) (jp N φ x) - topPot (k + 2) (jp N φ x)) I)
          * wSym k (jp N φ x) I)) =
      fun x => ρ x ^ 2 * (2 * ∑ I : Fin (k + 3) → Fin n,
          rPot (k + 2) (jp N φ x) I * wSym k (jp N φ x) I)
        + ρ x ^ 2 * (2 * ∑ I : Fin (k + 3) → Fin n,
          Dz (fun q : P n N => jetSym potSym (k + 2) q - topPot (k + 2) q) (jp N φ x) I *
            wSym k (jp N φ x) I)
        + ρ x ^ 2 * (2 * ∑ I : Fin (k + 3) → Fin n,
          cTop k (jp N φ x) (topPot (k + 2) (jp N φ x)) I * wSym k (jp N φ x) I)
        + ρ x ^ 2 * (2 * ∑ I : Fin (k + 3) → Fin n,
          cTop k (jp N φ x) (jetSym potSym (k + 2) (jp N φ x) - topPot (k + 2) (jp N φ x)) I *
            wSym k (jp N φ x) I) := by
    funext x
    simp only [add_mul, Finset.sum_add_distrib, mul_add]
  rw [e]
  exact ((i2.fun_add i3).fun_add i4).fun_add i5

/-! ### The local evolution inequality -/

/-- The local evolution inequality: for a graph function solving the scalar equation on the ball of radius 3,
the cut-off energy of order `k` has time derivative bounded above by minus a multiple of the cut-off energy two
orders higher, plus small multiples of higher energies on larger balls, plus lower order energy. -/
theorem exists_const_local_evolution {N k : ℕ} (hk : 1 ≤ k) (hN : k + 7 ≤ N) {K : ℝ} (hK : 0 ≤ K)
    (Λ : ℝ) {ε : ℝ} (hε : 0 < ε) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ (T : Set ℝ) (Φ : ℝ × EuclideanSpace ℝ (Fin n) → ℝ), GraphEvol N K 3 T Φ →
      ∀ ρ : EuclideanSpace ℝ (Fin n) → ℝ, IsCutoff ρ Λ → ∀ t ∈ T,
      ∃ D : ℝ,
        HasDerivAt (fun s : ℝ => ∫ x, ρ x ^ 2 * eSym k (jp N (fun y => Φ (s, y)) x)) D t ∧
        D ≤ -(2 * ((1 + K ^ 2)⁻¹) ^ (k + 5)) *
              (∫ x, ρ x ^ 2 * ‖jet (k + 5) (fun y => Φ (t, y)) x‖ ^ 2)
            + ε * (∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 2,
                ‖jet (k + 5) (fun y => Φ (t, y)) x‖ ^ 2)
            + ε * energy3 (k + 2) 3 (fun y => Φ (t, y))
            + C * (1 + energy3 (k + 1) 3 (fun y => Φ (t, y))) := by
  classical
  have hε2 : 0 < ε / 2 := by positivity
  obtain ⟨C₀, hC₀0, hC₀⟩ := exists_const_dt_eSym_sub_le (n := n) N k hk hN K
  obtain ⟨C₁, hC₁0, hC₁⟩ := exists_const_integral_Om_odd (n := n) (m := k + 2) (by omega)
    (ρ := 2) (R := 3) (by norm_num) (by norm_num) K (ε := ε / (2 * (C₀ + 1))) (by positivity)
  obtain ⟨C₂, hC₂0, hC₂⟩ := exists_const_main_le (n := n) (N := N) (k := k) hN hK Λ hε2
  obtain ⟨C₃, hC₃0, hC₃⟩ := exists_const_lower_le (n := n) (N := N) (k := k) hk hN (K := K) Λ hε2
  refine ⟨C₀ * C₁ + C₂ + C₃, by positivity, fun T Φ h ρ hρ t ht => ?_⟩
  have hφ : ContDiff ℝ ∞ (fun y : EuclideanSpace ℝ (Fin n) => Φ (t, y)) :=
    contDiff_slice h.isOpen h.smooth ht
  have hadm : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 3,
      jp N (fun y => Φ (t, y)) x ∈ admSet n N K := fun x hx => h.adm t ht x hx
  have hsupp : ∀ x ∈ tsupport ρ, jp N (fun y => Φ (t, y)) x ∈ admSet n N K :=
    fun x hx => hadm x (ibb_ball_sub x (hρ.supp hx))
  refine ⟨_, le_hasDerivAt_integral (by omega) h hρ ht, ?_⟩
  -- pointwise bound of the derivative
  have hpt : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 3,
      fderiv ℝ (fun q : ℝ × EuclideanSpace ℝ (Fin n) =>
          eSym k (jp N (fun y => Φ (q.1, y)) q.2)) (t, x) (1, 0) ≤
        (2 * ∑ I : Fin (k + 3) → Fin n, topPot (k + 3) (jp N (fun y => Φ (t, y)) x) I *
            wSym k (jp N (fun y => Φ (t, y)) x) I
          + 2 * ∑ I : Fin (k + 3) → Fin n,
          (rPot (k + 2) (jp N (fun y => Φ (t, y)) x) I
            + Dz (fun q : P n N => jetSym potSym (k + 2) q - topPot (k + 2) q)
                (jp N (fun y => Φ (t, y)) x) I
            + cTop k (jp N (fun y => Φ (t, y)) x) (topPot (k + 2) (jp N (fun y => Φ (t, y)) x)) I
            + cTop k (jp N (fun y => Φ (t, y)) x) (jetSym potSym (k + 2)
                (jp N (fun y => Φ (t, y)) x) - topPot (k + 2) (jp N (fun y => Φ (t, y)) x)) I)
          * wSym k (jp N (fun y => Φ (t, y)) x) I)
        + C₀ * Om (2 * k + 3) (k + 2) (aJet (fun y => Φ (t, y)) x) := by
    intro x hx
    obtain ⟨h1, h2⟩ := hC₀ 3 T Φ h t ht x hx
    rw [le_fderiv_eq (by omega) h ht hx h1]
    rw [h.jet_dt_eq ht hx (k + 3) (by omega), h.jet_dt_eq ht hx (k + 2) (by omega),
      le_decomp hN (hadm x hx)] at h2
    have := (abs_le.1 h2).2
    linarith
  set φ : EuclideanSpace ℝ (Fin n) → ℝ := fun y => Φ (t, y) with hφdef
  have hOm0 : ∀ x, 0 ≤ C₀ * Om (2 * k + 3) (k + 2) (aJet φ x) := fun x =>
    mul_nonneg hC₀0 (evp_Om_nonneg _ _ (aJet_nonneg φ x))
  have hbound : ∀ x, ρ x ^ 2 * fderiv ℝ (fun q : ℝ × EuclideanSpace ℝ (Fin n) =>
          eSym k (jp N (fun y => Φ (q.1, y)) q.2)) (t, x) (1, 0) ≤
      (ρ x ^ 2 * (2 * ∑ I : Fin (k + 3) → Fin n,
          topPot (k + 3) (jp N φ x) I * wSym k (jp N φ x) I)
        + ρ x ^ 2 * (2 * ∑ I : Fin (k + 3) → Fin n,
          (rPot (k + 2) (jp N φ x) I
            + Dz (fun q : P n N => jetSym potSym (k + 2) q - topPot (k + 2) q) (jp N φ x) I
            + cTop k (jp N φ x) (topPot (k + 2) (jp N φ x)) I
            + cTop k (jp N φ x) (jetSym potSym (k + 2) (jp N φ x) - topPot (k + 2) (jp N φ x)) I)
          * wSym k (jp N φ x) I))
      + (ball (0 : EuclideanSpace ℝ (Fin n)) 2).indicator
          (fun y => C₀ * Om (2 * k + 3) (k + 2) (aJet φ y)) x := by
    intro x
    by_cases hx : x ∈ tsupport ρ
    · have hx2 : x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 2 := hρ.supp hx
      have hx3 := ibb_ball_sub x hx2
      rw [Set.indicator_of_mem hx2]
      have h1 := hpt x hx3
      have hsq0 : 0 ≤ ρ x ^ 2 := sq_nonneg _
      have hsq1 : ρ x ^ 2 ≤ 1 := by
        have := hρ.le_one x
        rw [← sq_abs]
        exact pow_le_one₀ (abs_nonneg _) this
      have h2 := mul_le_mul_of_nonneg_left h1 hsq0
      have h3 := mul_le_of_le_one_left (hOm0 x) hsq1
      nlinarith [h2, h3]
    · have h0 : ρ x ^ 2 = 0 := by rw [image_eq_zero_of_notMem_tsupport hx]; norm_num
      rw [h0]
      have hi : 0 ≤ (ball (0 : EuclideanSpace ℝ (Fin n)) 2).indicator
          (fun y => C₀ * Om (2 * k + 3) (k + 2) (aJet φ y)) x :=
        Set.indicator_nonneg (fun y _ => hOm0 y) x
      simpa using hi
  have iP := le_integrable_P (k := k) hN hρ hφ hsupp
  have iQ := le_integrable_Q hk hN hρ hφ hsupp
  have iInd : Integrable ((ball (0 : EuclideanSpace ℝ (Fin n)) 2).indicator
      (fun y => C₀ * Om (2 * k + 3) (k + 2) (aJet φ y))) :=
    (integrable_indicator_iff measurableSet_ball).2
      (integrableOn_ball_of_continuous (continuous_const.mul (continuous_Om hφ _ _)) 2)
  have hmono := integral_mono (le_integrable_deriv (by omega) h hρ ht) ((iP.fun_add iQ).fun_add iInd)
    hbound
  rw [integral_add (iP.fun_add iQ) iInd, integral_add iP iQ,
    integral_indicator measurableSet_ball, integral_const_mul] at hmono
  have bP := hC₂ ρ hρ φ hφ hadm
  have bQ := hC₃ ρ hρ φ hφ hadm
  have hK3 : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 3, ‖jet 3 φ x‖ ≤ K :=
    fun x hx => ibb_norm3 (by omega) (hadm x hx)
  have bO : ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 2, Om (2 * k + 3) (k + 2) (aJet φ x) ≤
      ε / (2 * (C₀ + 1)) * energy3 (k + 2) 3 φ + C₁ := by
    have := hC₁ φ hφ hK3
    have e : 2 * (k + 2) - 1 = 2 * k + 3 := by omega
    rwa [e] at this
  have hq : C₀ * (ε / (2 * (C₀ + 1))) ≤ ε / 2 := by
    calc C₀ * (ε / (2 * (C₀ + 1))) = ε / 2 * (C₀ / (C₀ + 1)) := by field_simp
      _ ≤ ε / 2 * 1 :=
        mul_le_mul_of_nonneg_left (div_le_one_of_le₀ (by linarith) (by positivity)) hε2.le
      _ = ε / 2 := mul_one _
  have hE2 := energy3_nonneg (n := n) (k + 2) 3 φ
  have hE1 := energy3_nonneg (n := n) (k + 1) 3 φ
  have hm1 := mul_le_mul_of_nonneg_left bO hC₀0
  have hm2 := mul_le_mul_of_nonneg_right hq hE2
  have hm3 := mul_nonneg (mul_nonneg hC₀0 hC₁0) hE1
  have hQ' := (le_abs_self _).trans bQ
  nlinarith [hmono, bP, hQ', hm1, hm2, hm3]

end CWFlow.Extension

/-- info: 'CWFlow.Extension.exists_const_local_evolution' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms CWFlow.Extension.exists_const_local_evolution
