import CWFlow.Extension.GraphLocal
import CWFlow.Extension.IBPBounds

/-!
# The potential, the horizontal velocity and the gauge term of a graph

For the graph `x ↦ y₀ + U x + J U ∇φ(x)` over a Lagrangian plane:

* `potX_graphX`: the potential `div(JH)` is the jet-space symbol `potSym` on the jets of `φ`;
* `vSym`: the horizontal component of the velocity `J∇div(JH)` is a jet-space symbol of weight `(2, 2)`
  (`isSymbol_vSym`, `horCL_velX_graphX`);
* `exists_const_gauge_term`: an integral estimate for the gauge term
  `∫ eSym k · ⟨∇(ρ²), V⟩`.
-/

open Metric Set MeasureTheory
open CWFlow.Interp
open scoped ContDiff Topology RealInnerProductSpace

namespace CWFlow.Extension

/-! ### Part 1: the potential of a graph -/

theorem potX_graphX {n : ℕ} {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
    {J : V →L[ℝ] V} (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} (hU : IsLagFrame J U) (y₀ : V)
    {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ) {N : ℕ} (hN : 4 ≤ N) (x) :
    potX J (graphX J U y₀ φ) x = potSym (jp N φ x) := by
  unfold potX potSym
  rw [ginvZ_graphX hJi hU hφ, SX_graphX hJi hJJ hU hφ (N := N) (k := 1) (by omega),
    low23_jp (by omega)]

/-! ### Part 2: the horizontal velocity as a symbol -/

/-- The horizontal velocity of the graph flow: `V_c = −Σ_{a,b} g^{ab} ∂_b(pot) u_{ac}`. -/
noncomputable def vSym {n N : ℕ} (p : P n N) (c : Fin n) : ℝ :=
  -∑ a, ∑ b, ginvFun (low23 p).1 ![a, b] * totalD (potSym (N := N)) p (EuclideanSpace.single b 1) *
    (low23 p).1 ![a, c]

theorem gt_isSymbol_low23_fst {n N : ℕ} {A : Set (P n N)} (hA : IsOpen A)
    (hlow : ∀ m : Fin (N + 1), (m : ℕ) ≤ 3 → ∃ K : ℝ, ∀ p ∈ A, ‖p.2 m‖ ≤ K) (a c : Fin n) :
    IsSymbol A 0 0 (fun p : P n N => (low23 p).1 ![a, c]) := by
  have h : IsSymbol A 0 0 (fun p : P n N => (low23 p).1) :=
    isSymbol_comp_low hA low23 low23_high (bdd_low23 hlow)
      (c := fun q : Z n 2 × Z n 3 => q.1) contDiff_fst
  exact h.clm_comp (PiLp.proj (𝕜 := ℝ) 2 (fun _ : Fin 2 → Fin n => ℝ) ![a, c])

theorem isSymbol_vSym {n N : ℕ} {A : Set (P n N)} (hA : IsOpen A)
    (hlow : ∀ m : Fin (N + 1), (m : ℕ) ≤ 3 → ∃ K : ℝ, ∀ p ∈ A, ‖p.2 m‖ ≤ K) (hN : 5 ≤ N) (c : Fin n) :
    IsSymbol A 2 2 (fun p => vSym p c) := by
  have hG := isSymbol_ginv hA (bdd_low23 hlow)
  have hg : ∀ a b : Fin n, IsSymbol A 0 0 (fun p : P n N => ginvFun (low23 p).1 ![a, b]) :=
    fun a b => hG.clm_comp (PiLp.proj (𝕜 := ℝ) 2 (fun _ : Fin 2 → Fin n => ℝ) ![a, b])
  have hP := (isSymbol_potSym hA hlow (by omega : 4 ≤ N)).totalD (by omega) hlow
  have hP' : ∀ b : Fin n, IsSymbol A 2 2
      (fun p : P n N => totalD (potSym (N := N)) p (EuclideanSpace.single b 1)) := fun b =>
    (hP.clm_comp (ContinuousLinearMap.apply ℝ ℝ
      (EuclideanSpace.single b (1 : ℝ) : EuclideanSpace ℝ (Fin n)))).mono (by omega) (by omega)
  have hterm : ∀ a b : Fin n, IsSymbol A 2 2 (fun p : P n N =>
      ginvFun (low23 p).1 ![a, b] * totalD (potSym (N := N)) p (EuclideanSpace.single b 1) *
        (low23 p).1 ![a, c]) := fun a b =>
    ((((hg a b).mul (hP' b)).mul (gt_isSymbol_low23_fst hA hlow a c)).mono (by omega) (by omega))
  have hs := IsSymbol.sum hA (Finset.univ : Finset (Fin n)) fun a _ =>
    IsSymbol.sum hA (Finset.univ : Finset (Fin n)) fun b _ => hterm a b
  exact hs.neg

theorem gt_inner_JdX {n : ℕ} {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
    {J : V →L[ℝ] V} (hJJ : ∀ v : V, J (J v) = -v)
    {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} (hU : IsLagFrame J U) (y₀ : V)
    {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ) (x : EuclideanSpace ℝ (Fin n))
    (a c : Fin n) :
    ⟪U (eb c), J (dX (graphX J U y₀ φ) x a)⟫ = -jet 2 φ x ![a, c] := by
  rw [gb_dX' hφ, map_add, hJJ, inner_add_right, gb_UJU hU, inner_neg_right, gb_UU,
    real_inner_comm, gb_inner_toLp_eb, zero_add]

theorem horCL_velX_graphX {n : ℕ} {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
    {J : V →L[ℝ] V} (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} (hU : IsLagFrame J U) (y₀ : V)
    {φ : EuclideanSpace ℝ (Fin n) → ℝ} (hφ : ContDiff ℝ ∞ φ) {N : ℕ} (hN : 5 ≤ N) (x)
    (hadm : ∃ K, jp N φ x ∈ admSet n N K) (c : Fin n) :
    horCL U (velX J (graphX J U y₀ φ) x) c = vSym (jp N φ x) c := by
  obtain ⟨K, hK⟩ := hadm
  have hpot : potX J (graphX J U y₀ φ) = fun y => potSym (jp N φ y) :=
    funext fun y => potX_graphX hJi hJJ hU y₀ hφ (by omega) y
  have hd : HasFDerivAt (fun y => potSym (jp N φ y)) (totalD (potSym (N := N)) (jp N φ x)) x :=
    hasFDerivAt_symbol_jp (isSymbol_potSym (isOpen_admSet K) (admSet_low K) (by omega))
      (by omega) hφ hK
  have hb : ∀ b : Fin n, fderiv ℝ (potX J (graphX J U y₀ φ)) x (eb b) =
      totalD (potSym (N := N)) (jp N φ x) (EuclideanSpace.single b 1) := by
    intro b
    rw [hpot, hd.fderiv]
  have e : horCL U (velX J (graphX J U y₀ φ) x) c = ⟪U (eb c), velX J (graphX J U y₀ φ) x⟫ := rfl
  rw [e, velX]
  simp only [inner_sum, inner_smul_right, gt_inner_JdX hJJ hU y₀ hφ, hb,
    gX_graphX hJi hU hφ]
  unfold vSym
  rw [low23_jp (by omega)]
  simp only [mul_neg, Finset.sum_neg_distrib]
  congr 1

/-! ### Part 3: the gauge term -/

theorem gt_dsq {n : ℕ} {ρ : EuclideanSpace ℝ (Fin n) → ℝ} (hρ : ContDiff ℝ ∞ ρ)
    (x : EuclideanSpace ℝ (Fin n)) (c : Fin n) :
    fderiv ℝ (fun y => ρ y ^ 2) x (EuclideanSpace.single c 1) =
      2 * ρ x * fderiv ℝ ρ x (EuclideanSpace.single c 1) := by
  have hd : HasFDerivAt ρ (fderiv ℝ ρ x) x := ((hρ.differentiable (by simp)) x).hasFDerivAt
  have h : HasFDerivAt (fun y => ρ y * ρ y) (ρ x • fderiv ℝ ρ x + ρ x • fderiv ℝ ρ x) x :=
    hd.mul hd
  rw [show (fun y => ρ y ^ 2) = fun y => ρ y * ρ y from funext fun y => sq _, h.fderiv]
  simp
  ring

theorem exists_const_gauge_term {n N k : ℕ} (hk : 1 ≤ k) (hN : k + 7 ≤ N) {K : ℝ} (Λ : ℝ) {ε : ℝ}
    (hε : 0 < ε) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ ρ : EuclideanSpace ℝ (Fin n) → ℝ, IsCutoff ρ Λ →
      ∀ φ : EuclideanSpace ℝ (Fin n) → ℝ, ContDiff ℝ ∞ φ →
      (∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 3, jp N φ x ∈ admSet n N K) →
      |∫ x, eSym k (jp N φ x) * ∑ c, fderiv ℝ (fun y => ρ y ^ 2) x (EuclideanSpace.single c 1) *
          vSym (jp N φ x) c| ≤ ε * energy3 (k + 2) 3 φ + C := by
  have hA := isOpen_admSet (n := n) (N := N) K
  have hlow := admSet_low (n := n) (N := N) K
  have hS : ∀ c : Fin n, IsSymbol (admSet n N K) (2 * k + 3) (k + 2)
      (fun p : P n N => eSym k p * vSym p c) := fun c =>
    (((isSymbol_eSym hA hlow k (by omega)).mul (isSymbol_vSym hA hlow (by omega) c)).mono
      (by omega) (by omega))
  have hex : ∀ c : Fin n, ∃ C : ℝ, 0 ≤ C ∧ ∀ (φ : EuclideanSpace ℝ (Fin n) → ℝ)
      (x : EuclideanSpace ℝ (Fin n)), jp N φ x ∈ admSet n N K →
      ‖eSym k (jp N φ x) * vSym (jp N φ x) c‖ ≤ C * Om (2 * k + 3) (k + 2) (aJet φ x) :=
    fun c => ibb_pt (hS c) (by omega)
  choose C hC0 hC using hex
  obtain ⟨M, hMdef⟩ : ∃ M : ℝ, M = 2 * |Λ| * ∑ c, C c := ⟨_, rfl⟩
  have hM0 : 0 ≤ M := by
    rw [hMdef]
    exact mul_nonneg (by positivity) (Finset.sum_nonneg fun c _ => hC0 c)
  obtain ⟨C', hC', hint⟩ := exists_const_integral_Om_odd (n := n) (m := k + 2) (by omega)
    (ρ := 2) (R := 3) (by norm_num) (by norm_num) K (ε := ε / (M + 1)) (div_pos hε (by linarith))
  refine ⟨M * C', mul_nonneg hM0 hC', fun ρ hρ φ hφ hadm => ?_⟩
  have hK : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 3, ‖jet 3 φ x‖ ≤ K :=
    fun x hx => ibb_norm3 (by omega) (hadm x hx)
  have hom : ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 2, Om (2 * k + 3) (k + 2) (aJet φ x) ≤
      ε / (M + 1) * energy3 (k + 2) 3 φ + C' := by
    have := hint φ hφ hK
    have e : 2 * (k + 2) - 1 = 2 * k + 3 := by omega
    rwa [e] at this
  set h : EuclideanSpace ℝ (Fin n) → ℝ := fun x => eSym k (jp N φ x) *
    ∑ c, fderiv ℝ (fun y => ρ y ^ 2) x (EuclideanSpace.single c 1) * vSym (jp N φ x) c with hh
  have hb : ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 2,
      |h x| ≤ M * Om (2 * k + 3) (k + 2) (aJet φ x) := by
    intro x hx
    have hxa := hadm x (ibb_ball_sub x hx)
    have hρ1 := hρ.le_one x
    calc |h x| = |∑ c, 2 * ρ x * fderiv ℝ ρ x (EuclideanSpace.single c 1) *
          (eSym k (jp N φ x) * vSym (jp N φ x) c)| := by
          congr 1
          simp only [hh]
          rw [Finset.mul_sum]
          refine Finset.sum_congr rfl fun c _ => ?_
          rw [gt_dsq hρ.smooth]
          ring
      _ ≤ ∑ c, |2 * ρ x * fderiv ℝ ρ x (EuclideanSpace.single c 1) *
          (eSym k (jp N φ x) * vSym (jp N φ x) c)| := Finset.abs_sum_le_sum_abs _ _
      _ ≤ ∑ c, 2 * |Λ| * (C c * Om (2 * k + 3) (k + 2) (aJet φ x)) :=
          Finset.sum_le_sum fun c _ =>
            ibb_term hρ1 (ibb_d1 hρ x c)
              (by simpa [Real.norm_eq_abs] using hC c φ x hxa)
      _ = M * Om (2 * k + 3) (k + 2) (aJet φ x) := by
          rw [← Finset.mul_sum, ← Finset.sum_mul, hMdef]
          ring
  have hz : ∀ x, x ∉ ball (0 : EuclideanSpace ℝ (Fin n)) 2 → h x = 0 := by
    intro x hx
    have hxs : x ∉ tsupport ρ := fun h' => hx (hρ.supp h')
    have h0 : ρ x = 0 := image_eq_zero_of_notMem_tsupport hxs
    simp [hh, gt_dsq hρ.smooth, h0]
  have hS' : closedBall (0 : EuclideanSpace ℝ (Fin n)) 2 ⊆ {y | jp N φ y ∈ admSet n N K} :=
    fun y hy => hadm y (closedBall_subset_ball (by norm_num) hy)
  have hcont : ContinuousOn h (closedBall (0 : EuclideanSpace ℝ (Fin n)) 2) := by
    have he := (contDiffOn_symbol_jp (isSymbol_eSym hA hlow k (by omega)) hφ).continuousOn.mono hS'
    have hv : ∀ c : Fin n, ContinuousOn (fun y => vSym (jp N φ y) c)
        (closedBall (0 : EuclideanSpace ℝ (Fin n)) 2) := fun c =>
      (contDiffOn_symbol_jp (isSymbol_vSym hA hlow (by omega) c) hφ).continuousOn.mono hS'
    have hd : ∀ c : Fin n, Continuous (fun y => fderiv ℝ (fun y => ρ y ^ 2) y
        (EuclideanSpace.single c 1)) := fun c =>
      ((hρ.smooth.pow 2).continuous_fderiv (by simp)).clm_apply continuous_const
    exact he.mul (continuousOn_finsetSum _ fun c _ => (hd c).continuousOn.mul (hv c))
  have hint_h : IntegrableOn h (ball (0 : EuclideanSpace ℝ (Fin n)) 2) :=
    (hcont.integrableOn_compact (isCompact_closedBall 0 2)).mono_set ball_subset_closedBall
  have hq : M * (ε / (M + 1)) ≤ ε := by
    calc M * (ε / (M + 1)) = ε * (M / (M + 1)) := by ring
      _ ≤ ε * 1 := mul_le_mul_of_nonneg_left (div_le_one_of_le₀ (by linarith) (by linarith)) hε.le
      _ = ε := mul_one _
  rw [← setIntegral_eq_integral_of_forall_compl_eq_zero hz]
  calc |∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 2, h x|
      ≤ ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 2, |h x| := by
        simpa [Real.norm_eq_abs] using norm_integral_le_integral_norm (μ := volume.restrict
          (ball (0 : EuclideanSpace ℝ (Fin n)) 2)) h
    _ ≤ ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 2, M * Om (2 * k + 3) (k + 2) (aJet φ x) :=
        setIntegral_mono_on hint_h.abs
          ((integrableOn_ball_of_continuous (continuous_Om hφ _ _) 2).const_mul _)
          measurableSet_ball hb
    _ = M * ∫ x in ball (0 : EuclideanSpace ℝ (Fin n)) 2, Om (2 * k + 3) (k + 2) (aJet φ x) :=
        integral_const_mul _ _
    _ ≤ M * (ε / (M + 1) * energy3 (k + 2) 3 φ + C') := mul_le_mul_of_nonneg_left hom hM0
    _ ≤ ε * energy3 (k + 2) 3 φ + M * C' := by
        nlinarith [mul_le_mul_of_nonneg_right hq (energy3_nonneg (n := n) (k + 2) 3 φ)]

/-- info: 'CWFlow.Extension.potX_graphX' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms potX_graphX

/-- info: 'CWFlow.Extension.isSymbol_vSym' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms isSymbol_vSym

/-- info: 'CWFlow.Extension.horCL_velX_graphX' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms horCL_velX_graphX

/-- info: 'CWFlow.Extension.exists_const_gauge_term' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_const_gauge_term

end CWFlow.Extension
