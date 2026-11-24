import CWFlow.Extension.ChartSup
import CWFlow.Extension.GraphAt
import CWFlow.Extension.GaugeTerm
import CWFlow.Extension.PotReparam

/-!
# A uniform bound for the velocity and its derivative

For a Lagrangian immersion `F : M → V` of a compact manifold with small second fundamental form and bounded
energies, the velocity `velX J (F ∘ σ')` read in any parametrization `σ'` is bounded, and its derivative is
bounded by `C ‖D(F ∘ σ')‖` (`exists_const_vel_bound`).  The proof reads everything in a graph chart centred at
the point: for a graph `x ↦ y₀ + U x + J U ∇ψ(x)` the velocity is `vbV J U (jp N ψ x)` with `vbV` a globally
defined (on the admissible set) smooth function of the jets, and reparametrization invariance of the velocity
(`velX_comp`) transfers this to `F ∘ σ'`.
-/

open Metric Set
open CWFlow.Interp
open scoped ContDiff Manifold Topology RealInnerProductSpace

namespace CWFlow.Extension

/-! ### Part 1: the velocity of a graph as a function of jets -/

/-- The coefficients `g^{ab} ∂_b(pot)` as a function on jet space. -/
noncomputable def vbS {n N : ℕ} (p : P n N) (a b : Fin n) : ℝ :=
  ginvFun (low23 p).1 ![a, b] * totalD (potSym (N := N)) p (EuclideanSpace.single b 1)

/-- The vector `Σ_{a,b} S_{ab} e_a` (horizontal part of the velocity). -/
noncomputable def vbα {n N : ℕ} (p : P n N) : EuclideanSpace ℝ (Fin n) :=
  ∑ a, ∑ b, vbS p a b • eb a

/-- The vector `Σ_{a,b} S_{ab} Σ_c u_{ac} e_c` (vertical part of the velocity). -/
noncomputable def vbβ {n N : ℕ} (p : P n N) : EuclideanSpace ℝ (Fin n) :=
  ∑ a, ∑ b, vbS p a b • (WithLp.toLp 2 (fun c => (low23 p).1 ![a, c]) : EuclideanSpace ℝ (Fin n))

/-- The velocity of the graph over the frame `U`, as a function of the jets. -/
noncomputable def vbV {n N : ℕ} {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
    (J : V →L[ℝ] V) (U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V) (p : P n N) : V :=
  J (U (vbα p)) - U (vbβ p)

theorem vb_velX_graphX {n N : ℕ} {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
    {J : V →L[ℝ] V} (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} (hU : IsLagFrame J U) (y₀ : V)
    {ψ : EuclideanSpace ℝ (Fin n) → ℝ} (hψ : ContDiff ℝ ∞ ψ) (hN : 5 ≤ N)
    {x : EuclideanSpace ℝ (Fin n)} {K : ℝ} (hK : jp N ψ x ∈ admSet n N K) :
    velX J (graphX J U y₀ ψ) x = vbV J U (jp N ψ x) := by
  have hpot : potX J (graphX J U y₀ ψ) = fun y => potSym (jp N ψ y) :=
    funext fun y => potX_graphX hJi hJJ hU y₀ hψ (by omega) y
  have hd : HasFDerivAt (fun y => potSym (jp N ψ y)) (totalD (potSym (N := N)) (jp N ψ x)) x :=
    hasFDerivAt_symbol_jp (isSymbol_potSym (isOpen_admSet K) (admSet_low K) (by omega))
      (by omega) hψ hK
  have hb : ∀ b : Fin n, fderiv ℝ (potX J (graphX J U y₀ ψ)) x (eb b) =
      totalD (potSym (N := N)) (jp N ψ x) (EuclideanSpace.single b 1) := by
    intro b
    rw [hpot, hd.fderiv]
  unfold velX vbV vbα vbβ vbS
  simp only [hb, gX_graphX hJi hU hψ, gb_dX' hψ, low23_jp (by omega : 3 ≤ N)]
  have hterm : ∀ a b : Fin n,
      ((gMat (jet 2 ψ x))⁻¹ a b *
          totalD (potSym (N := N)) (jp N ψ x) (EuclideanSpace.single b 1)) •
        J (U (eb a) + J (U (WithLp.toLp 2 (fun c => jet 2 ψ x ![a, c])))) =
      (ginvFun (jet 2 ψ x) ![a, b] *
          totalD (potSym (N := N)) (jp N ψ x) (EuclideanSpace.single b 1)) • J (U (eb a)) -
      (ginvFun (jet 2 ψ x) ![a, b] *
          totalD (potSym (N := N)) (jp N ψ x) (EuclideanSpace.single b 1)) •
        U (WithLp.toLp 2 (fun c => jet 2 ψ x ![a, c])) := by
    intro a b
    rw [map_add, hJJ, smul_add, smul_neg, ← sub_eq_add_neg]
    rfl
  simp only [hterm, Finset.sum_sub_distrib, map_sum, map_smul]

/-! ### Part 1b: smoothness and bounds on compact sets of jet space -/

theorem vb_smooth {n N : ℕ} (hN : 5 ≤ N) (K : ℝ) :
    ContDiffOn ℝ ∞ (vbα (n := n) (N := N)) (admSet n N K) ∧
      ContDiffOn ℝ ∞ (vbβ (n := n) (N := N)) (admSet n N K) := by
  have hg : ∀ a b : Fin n, ContDiff ℝ ∞ (fun p : P n N => ginvFun (low23 p).1 ![a, b]) := fun a b =>
    (contDiff_coord ![a, b]).comp
      (contDiff_ginvFun.comp (contDiff_fst.comp (low23 (n := n) (N := N)).contDiff))
  have hu : ∀ a c : Fin n, ContDiff ℝ ∞ (fun p : P n N => (low23 p).1 ![a, c]) := fun a c =>
    (contDiff_coord ![a, c]).comp (contDiff_fst.comp (low23 (n := n) (N := N)).contDiff)
  have hP := (isSymbol_potSym (isOpen_admSet (n := n) (N := N) K) (admSet_low (n := n) (N := N) K)
    (by omega : 4 ≤ N)).totalD (by omega) (admSet_low (n := n) (N := N) K)
  have hS : ∀ a b : Fin n, ContDiffOn ℝ ∞ (fun p : P n N => vbS p a b) (admSet n N K) :=
    fun a b => by
      unfold vbS
      exact (hg a b).contDiffOn.mul (hP.smooth.clm_apply contDiffOn_const)
  refine ⟨?_, ?_⟩
  · unfold vbα
    exact ContDiffOn.sum fun a _ => ContDiffOn.sum fun b _ => (hS a b).smul contDiffOn_const
  · unfold vbβ
    refine ContDiffOn.sum fun a _ => ContDiffOn.sum fun b _ => (hS a b).smul ?_
    exact (contDiff_piLp' (𝕜 := ℝ) 2 fun c => hu a c).contDiffOn

theorem vb_bound_center {n N : ℕ} (hN : 5 ≤ N) (R : ℝ) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ p : P n N, ‖p‖ ≤ R →
      ‖vbα p‖ ≤ C ∧ ‖vbβ p‖ ≤ C ∧ ‖fderiv ℝ (vbα (n := n) (N := N)) p‖ ≤ C ∧
        ‖fderiv ℝ (vbβ (n := n) (N := N)) p‖ ≤ C := by
  obtain ⟨h1, h2⟩ := vb_smooth (n := n) hN (|R| + 1)
  have hA : IsOpen (admSet n N (|R| + 1)) := isOpen_admSet _
  have hcpt : IsCompact (closedBall (0 : P n N) R) := isCompact_closedBall _ _
  have hsub : closedBall (0 : P n N) R ⊆ admSet n N (|R| + 1) := by
    intro p hp m _
    have hpR : ‖p‖ ≤ R := by simpa using hp
    have h3 : ‖p.2 m‖ ≤ ‖p‖ := (norm_le_pi_norm p.2 m).trans (le_max_right _ _)
    have := le_abs_self R
    linarith
  obtain ⟨C1, hC1⟩ := hcpt.exists_bound_of_continuousOn (h1.continuousOn.mono hsub)
  obtain ⟨C2, hC2⟩ := hcpt.exists_bound_of_continuousOn (h2.continuousOn.mono hsub)
  obtain ⟨C3, hC3⟩ := hcpt.exists_bound_of_continuousOn
    ((h1.continuousOn_fderiv_of_isOpen hA (by simp)).mono hsub)
  obtain ⟨C4, hC4⟩ := hcpt.exists_bound_of_continuousOn
    ((h2.continuousOn_fderiv_of_isOpen hA (by simp)).mono hsub)
  refine ⟨|C1| + |C2| + |C3| + |C4|, by positivity, fun p hp => ?_⟩
  have hpm : p ∈ closedBall (0 : P n N) R := by simpa using hp
  have a1 := le_abs_self C1
  have a2 := le_abs_self C2
  have a3 := le_abs_self C3
  have a4 := le_abs_self C4
  have b1 := abs_nonneg C1
  have b2 := abs_nonneg C2
  have b3 := abs_nonneg C3
  have b4 := abs_nonneg C4
  have c1 := hC1 p hpm
  have c2 := hC2 p hpm
  have c3 := hC3 p hpm
  have c4 := hC4 p hpm
  refine ⟨?_, ?_, ?_, ?_⟩ <;> linarith

theorem vb_fderiv_jp_le (n N : ℕ) :
    ∃ c : ℝ, 1 ≤ c ∧ ∀ {ψ : EuclideanSpace ℝ (Fin n) → ℝ}, ContDiff ℝ ∞ ψ →
      ∀ (x : EuclideanSpace ℝ (Fin n)) (r : ℝ), 0 ≤ r → (∀ m ≤ N + 1, ‖jet m ψ x‖ ≤ r) →
      ∀ w, ‖fderiv ℝ (jp N ψ) x w‖ ≤ c * (1 + r) * ‖w‖ := by
  set c : ℝ := 1 + ∑ m : Fin (N + 1), ‖contractCLM n (m : ℕ)‖ with hc
  have hs : 0 ≤ ∑ m : Fin (N + 1), ‖contractCLM n (m : ℕ)‖ :=
    Finset.sum_nonneg fun m _ => ContinuousLinearMap.opNorm_nonneg _
  refine ⟨c, by rw [hc]; linarith, ?_⟩
  intro ψ hψ x r hr hjet w
  have hc1 : 1 ≤ c := by rw [hc]; linarith
  rw [(hasFDerivAt_jp (N := N) hψ x).fderiv]
  simp only [ContinuousLinearMap.prod_apply, ContinuousLinearMap.id_apply]
  rw [Prod.norm_def]
  refine max_le ?_ ?_
  · have : 1 * 1 * ‖w‖ ≤ c * (1 + r) * ‖w‖ := by
      gcongr
      · linarith
    linarith
  · refine (pi_norm_le_iff_of_nonneg (by positivity)).2 fun m => ?_
    show ‖fderiv ℝ (jet (m : ℕ) ψ) x w‖ ≤ _
    rw [fderiv_jet_eq_contract hψ]
    have h1 : ‖contract w (jet ((m : ℕ) + 1) ψ x)‖ ≤
        ‖contractCLM n (m : ℕ)‖ * ‖w‖ * ‖jet ((m : ℕ) + 1) ψ x‖ := by
      have := (contractCLM n (m : ℕ) w).le_opNorm (jet ((m : ℕ) + 1) ψ x)
      have h2 := (contractCLM n (m : ℕ)).le_opNorm w
      calc ‖contract w (jet ((m : ℕ) + 1) ψ x)‖
          = ‖contractCLM n (m : ℕ) w (jet ((m : ℕ) + 1) ψ x)‖ := rfl
        _ ≤ ‖contractCLM n (m : ℕ) w‖ * ‖jet ((m : ℕ) + 1) ψ x‖ := this
        _ ≤ ‖contractCLM n (m : ℕ)‖ * ‖w‖ * ‖jet ((m : ℕ) + 1) ψ x‖ := by gcongr
    have h3 : ‖contractCLM n (m : ℕ)‖ ≤ ∑ m : Fin (N + 1), ‖contractCLM n (m : ℕ)‖ :=
      Finset.single_le_sum (f := fun m : Fin (N + 1) => ‖contractCLM n (m : ℕ)‖)
        (fun m _ => ContinuousLinearMap.opNorm_nonneg _) (Finset.mem_univ m)
    have h4 : ‖jet ((m : ℕ) + 1) ψ x‖ ≤ r := hjet _ (by have := m.2; omega)
    calc ‖contract w (jet ((m : ℕ) + 1) ψ x)‖
        ≤ ‖contractCLM n (m : ℕ)‖ * ‖w‖ * ‖jet ((m : ℕ) + 1) ψ x‖ := h1
      _ ≤ (c * ‖w‖) * r := by
          gcongr
          linarith
      _ ≤ c * (1 + r) * ‖w‖ := by
          have := norm_nonneg w
          nlinarith [mul_nonneg (by linarith : (0:ℝ) ≤ c) this]

/-! ### Part 2: dependence on the germ -/

theorem vb_potX_congr {n : ℕ} {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
    {J : V →L[ℝ] V} {X X' : EuclideanSpace ℝ (Fin n) → V} {x : EuclideanSpace ℝ (Fin n)}
    (h : X =ᶠ[𝓝 x] X') : potX J X x = potX J X' x := by
  unfold potX
  rw [ginvZ_congr h, SX_congr h 1]

theorem vb_velX_congr {n : ℕ} {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
    {J : V →L[ℝ] V} {X X' : EuclideanSpace ℝ (Fin n) → V} {x : EuclideanSpace ℝ (Fin n)}
    (h : X =ᶠ[𝓝 x] X') : velX J X x = velX J X' x := by
  have hp : potX J X =ᶠ[𝓝 x] potX J X' := by
    filter_upwards [h.eventually_nhds] with y hy using vb_potX_congr hy
  unfold velX
  rw [gX_congr h, hp.fderiv_eq]
  simp only [dX_congr h]

/-! ### Part 3: reparametrization onto a graph -/

theorem vb_norm_le_fderiv_graphX {n : ℕ} {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
    {J : V →L[ℝ] V} {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} (hU : IsLagFrame J U) (y₀ : V)
    {ψ : EuclideanSpace ℝ (Fin n) → ℝ} (hψ : ContDiff ℝ ∞ ψ) (x w : EuclideanSpace ℝ (Fin n)) :
    ‖w‖ ≤ ‖fderiv ℝ (graphX J U y₀ ψ) x w‖ := by
  have h : ⟪fderiv ℝ (graphX J U y₀ ψ) x w, U w⟫ = ⟪w, w⟫ := by
    rw [gb_fderiv_graphX_apply hψ, inner_add_left, gb_UU, gb_JUU hU, add_zero]
  have h2 := real_inner_le_norm (fderiv ℝ (graphX J U y₀ ψ) x w) (U w)
  rw [h, real_inner_self_eq_norm_mul_norm, U.norm_map] at h2
  by_cases hw : w = 0
  · simp [hw]
  · exact le_of_mul_le_mul_right h2 (norm_pos_iff.2 hw)

theorem vb_reparam {n : ℕ} {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
    {J : V →L[ℝ] V} (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} (hU : IsLagFrame J U) (y₀ : V)
    {ψ : EuclideanSpace ℝ (Fin n) → ℝ} (hψ : ContDiff ℝ ∞ ψ)
    {X : EuclideanSpace ℝ (Fin n) → V} {B : Set (EuclideanSpace ℝ (Fin n))} (hB : IsOpen B)
    (hX : ContDiffOn ℝ ∞ X B) (hinj : ∀ z ∈ B, Function.Injective (fderiv ℝ X z))
    {y : EuclideanSpace ℝ (Fin n)} (hy : y ∈ B)
    (heq : ∀ᶠ z in 𝓝 y, X z = graphX J U y₀ ψ (horC U y₀ (X z))) :
    (velX J X =ᶠ[𝓝 y] fun z => velX J (graphX J U y₀ ψ) (horC U y₀ (X z))) ∧
    DifferentiableAt ℝ (fun z => horC U y₀ (X z)) y ∧
    (∀ v, ‖fderiv ℝ (fun z => horC U y₀ (X z)) y v‖ ≤ ‖fderiv ℝ X y v‖) ∧
    fderiv ℝ X y = (fderiv ℝ (graphX J U y₀ ψ) (horC U y₀ (X y))).comp
      (fderiv ℝ (fun z => horC U y₀ (X z)) y) := by
  set G := graphX J U y₀ ψ with hGdef
  set τ : EuclideanSpace ℝ (Fin n) → EuclideanSpace ℝ (Fin n) := fun z => horC U y₀ (X z) with hτdef
  have hG : ContDiff ℝ ∞ G := contDiff_graphX hψ
  have hhor : ContDiff ℝ ∞ (horC U y₀) := by
    have : ContDiff ℝ ∞ (fun v : V => horCL U (v - y₀)) :=
      (horCL U).contDiff.comp (contDiff_id.sub contDiff_const)
    exact this
  have hτ : ContDiffOn ℝ ∞ τ B := hhor.comp_contDiffOn hX
  have hτd : ∀ z ∈ B, DifferentiableAt ℝ τ z := fun z hz =>
    (hτ.contDiffAt (hB.mem_nhds hz)).differentiableAt (by simp)
  obtain ⟨t, ht, hto, hyt⟩ := _root_.eventually_nhds_iff.1 heq
  have hNo : IsOpen (t ∩ B) := hto.inter hB
  have hyN : y ∈ t ∩ B := ⟨hyt, hy⟩
  have hXτ : ∀ z ∈ t ∩ B, X =ᶠ[𝓝 z] G ∘ τ := fun z hz => by
    filter_upwards [hNo.mem_nhds hz] with w hw using ht w hw.1
  have hder : ∀ z ∈ t ∩ B, fderiv ℝ X z = (fderiv ℝ G (τ z)).comp (fderiv ℝ τ z) := fun z hz => by
    rw [(hXτ z hz).fderiv_eq]
    exact fderiv_comp z (hG.differentiable (by simp) (τ z)) (hτd z hz.2)
  have hinjτ : ∀ z ∈ t ∩ B, Function.Injective (fderiv ℝ τ z) := fun z hz => by
    have h1 := hinj z hz.2
    rw [hder z hz] at h1
    exact Function.Injective.of_comp (f := ⇑(fderiv ℝ G (τ z))) h1
  have hdet : ∀ z ∈ t ∩ B, (jacM τ z).det ≠ 0 := fun z hz => by
    rw [← im_det_eq_jacM, ContinuousLinearMap.det]
    have h := LinearEquiv.isUnit_det'
      (LinearEquiv.ofInjectiveEndo (fderiv ℝ τ z : EuclideanSpace ℝ (Fin n) →ₗ[ℝ] EuclideanSpace ℝ (Fin n))
        (hinjτ z hz))
    have e : ((LinearEquiv.ofInjectiveEndo
        (fderiv ℝ τ z : EuclideanSpace ℝ (Fin n) →ₗ[ℝ] EuclideanSpace ℝ (Fin n)) (hinjτ z hz) :
        EuclideanSpace ℝ (Fin n) ≃ₗ[ℝ] EuclideanSpace ℝ (Fin n)) :
        EuclideanSpace ℝ (Fin n) →ₗ[ℝ] EuclideanSpace ℝ (Fin n)) =
        (fderiv ℝ τ z : EuclideanSpace ℝ (Fin n) →ₗ[ℝ] EuclideanSpace ℝ (Fin n)) :=
      LinearMap.ext fun v => by
        rw [LinearEquiv.coe_coe, LinearEquiv.coe_ofInjectiveEndo]
    rw [e] at h
    exact h.ne_zero
  have hvel : ∀ z ∈ t ∩ B, velX J X z = velX J G (τ z) := fun z hz => by
    rw [vb_velX_congr (hXτ z hz)]
    exact velX_comp (J := J) (X := G) (τ := τ) (Ω := univ) (N := t ∩ B) isOpen_univ
      hG.contDiffOn (fun x _ => injective_fderiv_graphX hJi hJJ hU hψ x)
      (fun x _ u w => lagrangian_graphX hJi hJJ hU hψ x u w) hNo (hτ.mono inter_subset_right)
      (fun _ _ => mem_univ _) hdet hz
  refine ⟨?_, hτd y hy, fun v => ?_, hder y hyN⟩
  · filter_upwards [hNo.mem_nhds hyN] with z hz using hvel z hz
  · rw [hder y hyN, ContinuousLinearMap.comp_apply]
    exact vb_norm_le_fderiv_graphX hU y₀ hψ _ _

/-! ### Part 4: the velocity of a graph with bounded jets -/

theorem vb_norm_J {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V] {J : V →L[ℝ] V}
    (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (v : V) : ‖J v‖ = ‖v‖ := by
  rw [norm_eq_sqrt_real_inner, norm_eq_sqrt_real_inner, hJi]

theorem vb_graph_bound {n N : ℕ} {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
    {J : V →L[ℝ] V} (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    (hN : 5 ≤ N) (R : ℝ) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V}, IsLagFrame J U → ∀ (y₀ : V)
      {ψ : EuclideanSpace ℝ (Fin n) → ℝ}, ContDiff ℝ ∞ ψ →
      ∀ {x : EuclideanSpace ℝ (Fin n)}, ‖x‖ ≤ R → (∀ m ≤ N + 1, ‖jet m ψ x‖ ≤ R) →
      ‖velX J (graphX J U y₀ ψ) x‖ ≤ C ∧
        DifferentiableAt ℝ (velX J (graphX J U y₀ ψ)) x ∧
        ∀ w, ‖fderiv ℝ (velX J (graphX J U y₀ ψ)) x w‖ ≤ C * ‖w‖ := by
  obtain ⟨C0, hC0, hb⟩ := vb_bound_center (n := n) (N := N) hN R
  obtain ⟨c, hc1, hjp⟩ := vb_fderiv_jp_le n N
  refine ⟨2 * C0 * (1 + c * (1 + |R|)), by positivity, ?_⟩
  intro U hU y₀ ψ hψ x hxR hjet
  have hR0 : 0 ≤ R := (norm_nonneg _).trans (hjet 0 (by omega))
  have hRabs : |R| = R := abs_of_nonneg hR0
  have hpR : ‖jp N ψ x‖ ≤ R := by
    rw [Prod.norm_def]
    refine max_le hxR ((pi_norm_le_iff_of_nonneg hR0).2 fun m => ?_)
    exact hjet m (by have := m.2; omega)
  have hpK : jp N ψ x ∈ admSet n N (|R| + 1) := by
    intro m _
    have h3 : ‖(jp N ψ x).2 m‖ ≤ ‖jp N ψ x‖ :=
      (norm_le_pi_norm (jp N ψ x).2 m).trans (le_max_right _ _)
    have := le_abs_self R
    linarith
  obtain ⟨hbα, hbβ, hbα', hbβ'⟩ := hb _ hpR
  obtain ⟨hα, hβ⟩ := vb_smooth (n := n) hN (|R| + 1)
  have hAopen : IsOpen (admSet n N (|R| + 1)) := isOpen_admSet _
  have hαd : DifferentiableAt ℝ (vbα (n := n) (N := N)) (jp N ψ x) :=
    (hα.contDiffAt (hAopen.mem_nhds hpK)).differentiableAt (by simp)
  have hβd : DifferentiableAt ℝ (vbβ (n := n) (N := N)) (jp N ψ x) :=
    (hβ.contDiffAt (hAopen.mem_nhds hpK)).differentiableAt (by simp)
  set G := graphX J U y₀ ψ with hGdef
  have hev : velX J G =ᶠ[𝓝 x] fun y => vbV J U (jp N ψ y) := by
    have hopen : IsOpen {y | jp N ψ y ∈ admSet n N (|R| + 1)} :=
      hAopen.preimage (contDiff_jp hψ).continuous
    filter_upwards [hopen.mem_nhds hpK] with y hy
    exact vb_velX_graphX hJi hJJ hU y₀ hψ hN hy
  have hvα : HasFDerivAt (fun p : P n N => J (U (vbα p)))
      ((J.comp U.toContinuousLinearMap).comp (fderiv ℝ (vbα (n := n) (N := N)) (jp N ψ x)))
      (jp N ψ x) :=
    (J.comp U.toContinuousLinearMap).hasFDerivAt.comp (jp N ψ x) hαd.hasFDerivAt
  have hvβ : HasFDerivAt (fun p : P n N => U (vbβ p))
      (U.toContinuousLinearMap.comp (fderiv ℝ (vbβ (n := n) (N := N)) (jp N ψ x)))
      (jp N ψ x) :=
    U.toContinuousLinearMap.hasFDerivAt.comp (jp N ψ x) hβd.hasFDerivAt
  have hvV := (hvα.sub hvβ).comp x (hasFDerivAt_jp (N := N) hψ x)
  have hvx : HasFDerivAt (velX J G) _ x := hvV.congr_of_eventuallyEq hev
  refine ⟨?_, hvx.differentiableAt, fun w => ?_⟩
  · rw [vb_velX_graphX hJi hJJ hU y₀ hψ hN hpK]
    unfold vbV
    calc ‖J (U (vbα (jp N ψ x))) - U (vbβ (jp N ψ x))‖
        ≤ ‖J (U (vbα (jp N ψ x)))‖ + ‖U (vbβ (jp N ψ x))‖ := norm_sub_le _ _
      _ = ‖vbα (jp N ψ x)‖ + ‖vbβ (jp N ψ x)‖ := by rw [vb_norm_J hJi, U.norm_map, U.norm_map]
      _ ≤ 2 * C0 * (1 + c * (1 + |R|)) := by
          have h1 : 0 ≤ c * (1 + |R|) := by positivity
          nlinarith
  · rw [hvx.fderiv]
    simp only [ContinuousLinearMap.comp_apply, sub_apply]
    rw [← (hasFDerivAt_jp (N := N) hψ x).fderiv]
    set δ := fderiv ℝ (jp N ψ) x w with hδ
    have hδb : ‖δ‖ ≤ c * (1 + R) * ‖w‖ := hjp hψ x R hR0 hjet w
    have e1 : ‖(J.comp U.toContinuousLinearMap) (fderiv ℝ (vbα (n := n) (N := N)) (jp N ψ x) δ)‖ ≤
        C0 * ‖δ‖ := by
      show ‖J (U (fderiv ℝ (vbα (n := n) (N := N)) (jp N ψ x) δ))‖ ≤ _
      rw [vb_norm_J hJi, U.norm_map]
      exact (ContinuousLinearMap.le_opNorm _ _).trans (mul_le_mul_of_nonneg_right hbα' (norm_nonneg _))
    have e2 : ‖U.toContinuousLinearMap (fderiv ℝ (vbβ (n := n) (N := N)) (jp N ψ x) δ)‖ ≤
        C0 * ‖δ‖ := by
      show ‖U (fderiv ℝ (vbβ (n := n) (N := N)) (jp N ψ x) δ)‖ ≤ _
      rw [U.norm_map]
      exact (ContinuousLinearMap.le_opNorm _ _).trans (mul_le_mul_of_nonneg_right hbβ' (norm_nonneg _))
    calc _ ≤ _ := norm_sub_le _ _
      _ ≤ C0 * ‖δ‖ + C0 * ‖δ‖ := add_le_add e1 e2
      _ ≤ C0 * (c * (1 + R) * ‖w‖) + C0 * (c * (1 + R) * ‖w‖) := by gcongr
      _ ≤ 2 * C0 * (1 + c * (1 + |R|)) * ‖w‖ := by
          rw [hRabs]
          have h1 : 0 ≤ c * (1 + R) * ‖w‖ := by positivity
          have h2 : 0 ≤ C0 * (c * (1 + R) * ‖w‖) := mul_nonneg hC0 h1
          have h3 : 0 ≤ C0 * ‖w‖ := mul_nonneg hC0 (norm_nonneg _)
          nlinarith

theorem vb_contDiff_horC {n : ℕ} {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
    {U : EuclideanSpace ℝ (Fin n) →ₗᵢ[ℝ] V} (y₀ : V) : ContDiff ℝ ∞ (horC U y₀) := by
  have : ContDiff ℝ ∞ (fun v : V => horCL U (v - y₀)) :=
    (horCL U).contDiff.comp (contDiff_id.sub contDiff_const)
  exact this

/-! ### Part 5: the main theorem -/

theorem exists_const_vel_bound {n : ℕ} {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    [FiniteDimensional ℝ E] {M : Type*} [TopologicalSpace M] [ChartedSpace E M] [IsManifold 𝓘(ℝ, E) ∞ M]
    [T2Space M] [CompactSpace M]
    {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V] [FiniteDimensional ℝ V]
    (ι : EuclideanSpace ℝ (Fin n) ≃L[ℝ] E) (A : IntAtlas (E := E) M) {J : V →L[ℝ] V}
    (hJi : ∀ u v : V, ⟪J u, J v⟫ = ⟪u, v⟫) (hJJ : ∀ v : V, J (J v) = -v)
    (hdim : Module.finrank ℝ V = 2 * n)
    {c₀ c₁ : ℝ} (hc₀ : 0 < c₀) (hc₀1 : c₀ ≤ 1 / 100) (hc₁ : 0 ≤ c₁)
    (hsmall : 6 * (c₁ * (1 + c₀ ^ 2) ^ ((3 : ℝ) / 2)) < c₀)
    (hK₃ : (1 + c₀ ^ 2) ^ ((3 : ℝ) / 2) * c₁ ≤ 1) (B : ℕ → ℝ) :
    ∃ C : ℝ, ∀ {F : M → V}, IsLagImm ι J F →
      (∀ z : M, ∀ y ∈ chartDomE ι z, aNormSq J (F ∘ chartInv ι z) y ≤ c₁ ^ 2) →
      (∀ l, energyM ι A J l F ≤ B l) →
      ∀ {B' : Set (EuclideanSpace ℝ (Fin n))} {σ' : EuclideanSpace ℝ (Fin n) → M}, IsParam ι B' σ' →
      ∀ y ∈ B', ‖velX J (F ∘ σ') y‖ ≤ C ∧
        ∀ v, ‖fderiv ℝ (velX J (F ∘ σ')) y v‖ ≤ C * ‖fderiv ℝ (F ∘ σ') y v‖ := by
  classical
  choose C hC using fun m : ℕ => exists_const_chart_sup ι A hJi hJJ hc₀.le hc₀1 hK₃ m B
  set R : ℝ := 1 + ∑ m ∈ Finset.range 7, |C m| with hRdef
  have hR0 : 0 ≤ R := by rw [hRdef]; positivity
  obtain ⟨Cg, hCg0, hCg⟩ := vb_graph_bound (n := n) (N := 5) hJi hJJ le_rfl R
  refine ⟨Cg, ?_⟩
  intro F hF hA hE B' σ' hσ' y hy
  set q : M := σ' y with hq
  obtain ⟨U, σ, φ, hU, hP, hσ0, hφ, hφ0, hb, -, hcomp⟩ :=
    graph_at ι hJi hJJ hdim hF hc₀ (by linarith) hc₁ hsmall hA q
  have hch : IsGraphChart ι J F c₀ q U σ φ :=
    ⟨hU, hP, hσ0, hφ, hφ0, fun x hx => ⟨(hb x hx).1, (hb x hx).2.1, (hb x hx).2.2.1,
      (hb x hx).2.2.2.1⟩, hcomp⟩
  have hK := graphChart_jet3_le ι hJi hJJ hF hc₀.le hc₁ hA hch
  have hjφ : ∀ m, ∀ x ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 2, ‖jet m φ x‖ ≤ C m :=
    fun m x hx => hC m hF (fun l _ => hE l) hch hK x hx
  have h0ball : (0 : EuclideanSpace ℝ (Fin n)) ∈ ball (0 : EuclideanSpace ℝ (Fin n)) 6 :=
    mem_ball_self (by norm_num)
  obtain ⟨ψ, hψ, hev⟩ := cb_exists_global_eventually isOpen_ball hφ h0ball
  have hjψ : ∀ m ≤ 5 + 1, ‖jet m ψ 0‖ ≤ R := by
    intro m hm
    rw [← jet_congr_of_eventuallyEq hev m]
    have h1 : ‖jet m φ 0‖ ≤ C m := hjφ m 0 (mem_ball_self (by norm_num))
    have h2 : |C m| ≤ ∑ m ∈ Finset.range 7, |C m| :=
      Finset.single_le_sum (f := fun m => |C m|) (fun m _ => abs_nonneg _)
        (Finset.mem_range.2 (by omega))
    have h3 := le_abs_self (C m)
    rw [hRdef]
    linarith
  set G := graphX J U (F q) ψ with hGdef
  -- the reparametrization
  obtain ⟨hFσ, hFinj, -⟩ := IsLagImm.comp_param ι hF hσ'
  have hgraph : ∀ᶠ x in 𝓝 (0 : EuclideanSpace ℝ (Fin n)),
      graphX J U (F q) φ x = graphX J U (F q) ψ x := by
    have hev' : ∀ᶠ x in 𝓝 (0 : EuclideanSpace ℝ (Fin n)), φ =ᶠ[𝓝 x] ψ := hev.eventually_nhds
    filter_upwards [hev'] with x hx
    unfold graphX gradE
    rw [hx.fderiv_eq]
  have hcont : ContinuousAt (fun z => F (σ' z)) y :=
    (hFσ.continuousOn.continuousAt (hσ'.isOpen.mem_nhds hy))
  have hτc : ContinuousAt (fun z => horC U (F q) (F (σ' z))) y :=
    (vb_contDiff_horC (U := U) (F q)).continuous.continuousAt.comp hcont
  have hτ0 : horC U (F q) ((F ∘ σ') y) = 0 := by
    show horC U (F q) (F q) = 0
    rw [horC_eq, sub_self, map_zero]
  have himg : ∀ᶠ z in 𝓝 y, σ' z ∈ σ '' ball (0 : EuclideanSpace ℝ (Fin n)) 6 := by
    refine (hσ'.cont.continuousAt (hσ'.isOpen.mem_nhds hy)).eventually_mem ?_
    exact hP.isOpen_image.mem_nhds ⟨0, h0ball, hσ0⟩
  have hg2 : ∀ᶠ z in 𝓝 y, graphX J U (F q) φ (horC U (F q) (F (σ' z))) =
      graphX J U (F q) ψ (horC U (F q) (F (σ' z))) := by
    have := hτc.tendsto
    have h0 : horC U (F q) (F (σ' y)) = 0 := hτ0
    rw [h0] at this
    exact this.eventually hgraph
  have heq : ∀ᶠ z in 𝓝 y, F (σ' z) = G (horC U (F q) (F (σ' z))) := by
    filter_upwards [himg, hg2] with z hz1 hz2
    obtain ⟨x, hx, hxz⟩ := hz1
    rw [← hxz] at hz2 ⊢
    rw [(hb x hx).2.1] at hz2 ⊢
    rw [(hb x hx).1]
    exact hz2
  obtain ⟨hvel, hτd, hτn, -⟩ := vb_reparam hJi hJJ hU (F q) hψ (X := F ∘ σ') hσ'.isOpen hFσ
    hFinj hy heq
  obtain ⟨hv0, hvd, hvD⟩ := hCg hU (F q) hψ (x := 0) (by simpa using hR0) hjψ
  have key : ∀ w, ‖fderiv ℝ (velX J G) ((fun z => horC U (F q) ((F ∘ σ') z)) y) w‖ ≤
      Cg * ‖w‖ := by
    intro w
    rw [show (fun z => horC U (F q) ((F ∘ σ') z)) y = 0 from hτ0]
    exact hvD w
  have hvd' : DifferentiableAt ℝ (velX J G) ((fun z => horC U (F q) ((F ∘ σ') z)) y) := by
    rw [show (fun z => horC U (F q) ((F ∘ σ') z)) y = 0 from hτ0]
    exact hvd
  refine ⟨?_, fun v => ?_⟩
  · have e : velX J (F ∘ σ') y = velX J G ((fun z => horC U (F q) ((F ∘ σ') z)) y) :=
      hvel.eq_of_nhds
    rw [e, show (fun z => horC U (F q) ((F ∘ σ') z)) y = 0 from hτ0]
    exact hv0
  · have hcomp := fderiv_comp y hvd' hτd
    have hfd : fderiv ℝ (velX J (F ∘ σ')) y = (fderiv ℝ (velX J G)
        ((fun z => horC U (F q) ((F ∘ σ') z)) y)).comp
        (fderiv ℝ (fun z => horC U (F q) ((F ∘ σ') z)) y) := by
      rw [hvel.fderiv_eq]
      exact hcomp
    rw [hfd, ContinuousLinearMap.comp_apply]
    exact (key _).trans (mul_le_mul_of_nonneg_left (hτn v) hCg0)

/-- info: 'CWFlow.Extension.exists_const_vel_bound' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in #print axioms exists_const_vel_bound

end CWFlow.Extension
