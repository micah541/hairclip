# T3 design note: one validated P1 slab in the sine class (final, after prototyping)

## 0. Equation (checked against hls.py)
Unknowns: th(α) = Σ_{n≥1} u_n sin nα, with θ = −π/2 + α + th. The clock is τ′ (dτ′ = (2π/L)⁴dt), and λ_n = n⁴ − n²:

  u_n′ = −λ_n u_n + N_n(u),  N = −∂(v² + v³/3) + (1+v)·∂⁻¹(w² − ⟨w²⟩),  v = th_α,  w = th_αα.

This is the task's N after one integration by parts: ∂⁻¹(th_ααα th_α − mean) = v v_α − ∂⁻¹(w² − ⟨w²⟩), using (v v_α)(0) = 0.
- hls.py integrates −th_αααα − Gθ_α. That equals −(∂⁴+∂²)th + N_task + th_αα(0)(1+th_α), and the last term vanishes for odd th.
- Checked in floats at t = 2: the IBP form agrees with the original to 5e-12 (on |N| ≈ 7e3), and the hls form agrees with ours to roundoff.
- ∂⁻¹ is the zero-mean antiderivative. On odd functions it vanishes at 0, so T_v(0) = 0.
- The class is invariant.
- Slaved quantities: (ln L)′ = −½Σn⁴u_n², t′ = (L/2π)⁴, X0′ = (L/2π)(Σn³u_n, 0).

## 1. Method: (a) — self-consistent bounds + Lohner block + polynomial tail, with a Chebyshev-in-time base
Per slab s ∈ [0,h], with h·θ_max⁴ = c ≈ 2 (θ_max = max|θ_α|). All rigorous quantities are in arb at 128 bits.
1. **Base ū.** Chebyshev–Lobatto polynomials of degree K = 16 in s, for M sine modes.
   - Float collocation Newton. GMRES with a dense block for modes ≤ M_f and a per-mode (2/h·D + λ_n)⁻¹ above.
   - Then 3 Newton polishes with the collocation residual computed in arb (≈1e-9 → 1e-23 → 1e-31).
   - The initial value is imposed on modes with λ_n h ≤ 8K². Stiffer modes satisfy the ODE at s = 0 (slaved start; well conditioned only for λh ≫ 2K²). The mismatch x̄ − ū(0) goes into z(0), where the Duhamel bounds damp it exactly.
   - The centre x̄ is carried as exact 128-bit arb midpoints. Rounding it to doubles would re-inject a 1e-16 inconsistency.
   - The rigorous residual r = ū′ + Λū − N(ū) has degree 3K. It is evaluated at 3K+1 nodes, transformed exactly, and re-expanded on each sub-interval for local sup bounds.
2. **Error split.** e = u − ū = P(s)ρ + z(s), with |ρ_i| ≤ w_i.
   - Generator columns G = [Q; S Q] are full M-vectors: Q is orthonormal on the **slow block** n ≤ M_a with λ_{M_a}h ≤ 8 (M_a ≈ 1.41θ_max at c = 2), and S Q = P_r(h)P_a(h)⁻¹Q is slaved.
   - Why this block: cond P_a(h) ≲ e⁸, so the slaved high parts are accurate and consistent. A block containing stiff modes (λh ≳ 30) makes the near-null directions inconsistent. Their defect then feeds back into ρ: a measured instability, ×20–1000 per slab.
   - P is the float collocation solution of P′ = (−Λ + J(s))P, Chebyshev degree K_P = 16. That degree resolves e^{−8s/h}; K_P = 8 left 1e-3 relative defects on the boundary columns.
   - The defect D = P′ + ΛP − JP is computed in arb exactly. J's coefficient functions are truncated at time degree K_J = 6, so D has degree 2K_J + K_P, and an arb majorant covers the remainder.
3. **z bounds.** z′ = −Λz + Jz + Q(e) − r − Dρ, with Q = N(ū+e) − N(ū) − Je.
   - Duhamel with exact e^{−λ_n s}: componentwise bounds for n ≤ N_c = 4M, and a tail T n^{−q} beyond.
   - Ten sub-intervals: [0, h/512, h/64] plus 8 uniform ones.
   - On each, a float fixed-point iteration Z ← RHS(Z)/(1−μ), μ = 3%, gives a uniform relative margin, with Q fed back after a failed round. This is followed by one strict arb check Z > max(a, e^{−λΔ}a + φ(Δ)[ĴZ + Q̂ + r̂ + d̂]) (Zgliczyński's first-exit argument).
   - **Intra-slab absorption** at each sub-interval end s_i: z_a = P_a(s_i)ζ, ρ ← ρ + P_a(s_i)⁻¹ζ, z_r ← z_r − S_iζ.
   - Every bound uses product matrices: |S_i| with S_i = P_r P_a⁻¹, and |B P_a(s_i)⁻¹|. Bounding |P_r||P_a⁻¹| instead loses a factor cond(P_a) and closes a feedback loop.
   - The crude bound therefore never acts over more than one sub-interval. Over a whole slab it amplifies by ~1e4 at c = 2.
4. **Lohner step at s = h.**
   - P_a(h)·diag(w) = QR with column sorting.
   - B = Q⁻¹P_a(h) in arb; w′ = |B|w + Σ_i|B P_a(s_i)⁻¹|ζ_i + |Q⁻¹|z_a(h).
   - Box′ = |P_r(h) − G′_rB|(w and ζ_i terms) + z_r(h) + |G′_r||Q⁻¹|z_a(h) + rad ū(h).
   - Tail′ comes from the tail inequality.

**Why not (b).**
- A radii-polynomial slab needs an approximate space-time inverse of dimension (K+1)M ≈ 6·10³ (t = 2) to 1.7·10⁴ (Γ0). Its Z₁ bound costs O(((K+1)M)³) ≈ 10¹¹–10¹³ arb operations per slab, i.e. hours in python-flint.
- Its output is a norm ball. Chaining balls multiplies the slab norms, and the majorant propagator has ℓ¹ norm ≈ 1.5–2 per slab at t = 2, giving 10⁷⁰⁺ over the run unless a variational block is added, which is (a).
- In (a) the stiff part is never inverted. The arb work is O((K_J+K_P)·M_a·M^{1.35}) for the defect, plus O(K·M log M) for the residual.

**Datum requirement (found here).** The base polynomial cannot represent an unresolved initial layer in modes with 8 < λh < 8K²: their Lagrange tail pollutes the whole slab. The incoming datum must lie on the slow manifold of the M-mode sine Galerkin flow.
- Chained slabs satisfy this automatically: ū(h) satisfies the ODE at a collocation node.
- The hls snapshots do not. At t = 2 their stiff modes are off by 1e-5–1e-2 relative, and exact Γ0 has an ODE residual of 1.1e4 at s = 0.
- P0/T2 must therefore define Γ0 consistently, e.g. Γ0 := ū(h) of an arb-polished preliminary slab after our own ETDRK4 relaxation, then re-close b_2 and recertify embeddedness.

## 2. Norms
- |e|_j = Σn^j|e_n|, with j = 3 the main norm (C³ in α). We also report j = 0..3.
- Physical C² of the curve ≲ (2π/L)|e|_1 + |δL| terms.
- Tail: |e_n| ≤ T n^{−q} for n > N_c, with **q = 6**.
- N has order 3, but in the IBP form its convolutions need Σn²|e_n|; |e|_3 and Q need q > 4, so q = 6 leaves a margin.
- The tail map contracts by C_J/((N_c+1)²−1), with C_J ≈ ‖ĉ‖₁ ≈ 620 at t = 2. Measured: T drops by ≥ 1e-6 per slab.

## 3. Explicit bounds (componentwise majorants of even, nonnegative sequences)
Forms: σ for sine series, γ for cosine series. Products: γ*γ, γ*σ, −σ*σ/4. Derivatives: γ = kσ/2, σ = −2kγ. ∂⁻¹ of a mean-free cosine series: σ = 2γ/k, with the mean dropped.
- **|a*b| ≤ |a|*|b|.** This is sharper than the ℓ¹_{(1+n)^j} algebra.
- **J.** With E_m = m y_m/2, H_m = m² y_m, Y_m = (Σ^{w̄}*H)_m/m:
  |Je|_n ≤ 2n(Γ^{c̄}*E)_n + (E*Σ^{W̄})_n + ((δ₀+Γ^{v̄})*Y)_n, where c̄ = 2v̄ + v̄² and Γ, Σ are slab sups (Chebyshev coefficient sums).
- **Q.** −∂((1+v̄)ε² + ε³/3) + εY_e + (1+v̄)Z₂ + εZ₂, with |Z₂|_m ≤ (H*H)_m/(2m).
- **Far tail.** For n > R: (a*b)_n ≤ 2^p n^{−p}(ā‖b‖₁ + b̄‖a‖₁), with ā = sup_{|m|>R/2} a_m|m|^p.
  - Base sequences have support ≤ 2M = N_c/2, so they contribute ā = 0.
  - Inside n ≤ N_c, y is carried explicitly to N_c + 2M, and far contributions are added as sup|a_far|·‖b‖₁.
- **Semigroup.** sup_{s≤Δ}|z_n| ≤ max(a_n, e^{−λΔ}a_n + φ_n f_n), with φ_n = −expm1(−λ_nΔ)/λ_n. The end value is e^{−λΔ}a_n + φ_n f_n.
- **Block.** As in §1.4. All absolute values are taken of products, never products of absolute values of ill-conditioned factors.

## 4. Slab output (chainable)
- **Endpoint, same format as the input:** (x̄′ exact 128-bit, G′ (M × M_a), w′, box′ for n ≤ N_c, tail T′).
- **Over the whole slab:**
  - the base ū (arb Chebyshev coefficients);
  - sup_s|e_n| ≤ Σ_k|P_k|w + max_i Z^{(i)}_n, hence sup|e|_j;
  - the enclosure of ∫⟨u_αα²⟩dτ′ = −Δln L (exact Chebyshev integral plus error).
- t follows from L by quadrature, and k_max from ū plus |e|_1.

## 5. Production additions (T3 production step)
- **Γ0 (option (c), t2/gamma0_final.json).** The T2 curve is a float snapshot at t = 0.001, before the self-crossing at t ≈ 0.0028. Its fast modes are off the slow manifold by δ_n ≈ 1e-10 (n = 186–350), 3e-12, 9e-14, 1e-16 and 1e-20 up to n = 1000, so λ_nδ_n ≈ 0.3–1. A polynomial slab cannot carry e^{−λ_n s}δ_n once λ_n h ≫ 8: with Δt = t/2 the residual is 2e-3.
  - *Slaving.* Modes 186..1000 are replaced by the values at s = 0 of a free-start Lobatto collocation over [0, 1.8t] in τ′. The ODE is imposed at all K + 1 nodes for n ≥ 186 and the initial value below; the solution is arb-polished. Free start is well conditioned only for λh ≳ 26, and λ_186h_s = 29.
  - *Second pass.* Modes ≥ 250 are then replaced by the s = 0 values of the first validated slab's own base (h₁ = t/2, free start for λh₁ > 26). The slab therefore starts with Γ0 − ū(0) at double rounding, 9e-23.
    - The one-pass mismatch (5.5e-18 at n = 250) sat in the box. The a priori J majorant coupled it into the block edge, J_blk ≈ 1.6e7·max Z_rest at n = 185, and it was absorbed as η = 2.5e-19: ρ 8e-24 → 1.3e-20 in slab 1.
    - With the second pass, ρ goes 8.14e-24 → 8.26e-24.
  - b_2 is then re-closed (arb Newton plus Krawczyk) and embeddedness re-certified (gamma0.py final, G0SNAP = gamma0s_src.npz). Modes ≤ 185 and > 1000 are the T2 values.
  - Rejected: (a) graded slabs from the unslaved curve need h₁ ≈ 1e-3t, then ≈ 9 doublings at M = 1000, about +1.3 h. (b) A later snapshot has the same float inconsistency.
- **Graded initial slabs.**
  - Δt_k = ½t_k with K = 24 (Chebyshev ratio 13 at Γ0, 10 required), while tκ_max⁴ < 4, i.e. t ≲ 0.03: 9 slabs.
  - Slab 1 uses a slaved start for λh > 26. The slaved Γ0 agrees with it to 5.5e-18, and the difference goes to the initial box.
  - M_a = #{λh ≤ 8}: 185 → 167 → 151 → … → 81.
  - The intra-slab absorption cap for slab k+1 is 1.5·η_k/max w_k (slab 1: preset), since a cap that is too small costs a full redo.
- **Adaptive h, every slab.** h ≤ c/θ_max⁴, shrunk until the float base's Chebyshev decay ratio r (degrees 1–6) satisfies r^K ≥ 1e24. Near the singularity the time scale is T − t, not 1/κ⁴. With K = 16 the rule halves the c = 2 length (swing); with K = 24 it accepts c = 2 at k ≈ 13 (r = 15.1). Use K = 24 from t = 4 on.
- **Block resizing at the Lohner step**, to a′ = #{λ h_next ≤ 8}.
  - *Growth.* Passenger columns ride the slab with ρ = 0. They enter P but not the defect bound. The seed for new mode m is e_m + (Λ_r − λ_m − J_rr)⁻¹J_{r,m}, shifted by its own decay rate. The unshifted quasi-static seed left a ~1e-4 inconsistency, which gave a measured feedback ×4–7 per slab.
  - *Drop.* A wide QR on rows 1..a′. Slaving uses the least-squares map G′_r = P_r(h)P_a(h)⁺Q. The annihilated direction goes to the box through Rmis = P_r(h) − G′_rB.
- **M resizing.** New modes enter with centre 0 (their error stays in the box). Dropped modes go to the box as |x̄_n| + |G_n|w + box_n; N_c = 4M follows, and the tail absorbs the overflow.
  - Rule: M += 80 while the top 40 modes of x̄ exceed 1e-30. In the swing this gave 480 → 800 by θ = 41.
  - Near T the truncation residual N_n(ū), n > M, is what sets M. At k ≈ 13, M = 1000 gives residual 3e-6 (box |·|_3 9e-8), while M = 1500 gives 2e-13. Hence M ≈ 1300 at k = 12.
- **Quadratures, exact Chebyshev algebra in arb.**
  - G_base(s) = ∫₀ˢ ½Σn⁴ū_n², with |G − G_base| ≤ L_err.
  - L′ = L·e^{−(G_base(h) ± L_err)}.
  - Δt = (L/2π)⁴[∫₀ʰ Σ_{j≤14}(−4G_base)^j/j! ± h·R₁₄]·e^{±4L_err}.
  - ΔX0 = (L/2π)[∫₀ʰ e^{−G_base}Σn³ū_n ± h(sup|p̄|(R + e^{L_err} − 1) + |e|_3)], with X0 = (x0, 0) on the axis.
  - k_max at the slab end: a grid of 2048 cells on [0, π] with a third-order Taylor bound, then branch and bound with direct arb evaluations, ± |e(h)|_1. On the t = 2 state the enclosure has width 9e-14 and lies at the tip α = π.

## 6. Lemma (first exit and Galerkin limit), self-contained
**Setting.**
- Sequences u = (u_n)_{n≥1} are sine coefficients. N_n(u) is given by the convolution formulas of §0. For |u_m| ≤ B_m with Σ m³B_m < ∞, every sum in N_n converges absolutely.
- On I = [0, Δ] fix: the polynomials ū, P (modes ≤ M), ρ with |ρ_i| ≤ w_i, r = ū′ + Λū − N(ū), D = P′ + ΛP − JP, and Q(s, e) = N(ū+e) − N(ū) − J(s)e.
- Let Y = (Z_1, …, Z_{N_c}, T(N_c+1)^{−q}, …) with q = 6, and W = {z : |z_n| ≤ Y_n ∀n}. Majorants Ĵ, Q̂ (§3) satisfy, for z ∈ W and s ∈ I: |(J z)_n| ≤ Ĵ(Y)_n and |Q_n(s, Pρ + z)| ≤ Q̂(Y)_n. The forcing bounds are |r_n| ≤ r̂_n and |(Dρ)_n| ≤ d̂_n.
- Put F_n = Ĵ(Y)_n + Q̂(Y)_n + r̂_n + d̂_n and g_n(t) = e^{−λ_n t}a_n + φ_n(t)F_n, where φ_n(t) = (1 − e^{−λ_n t})/λ_n (φ = t if λ_n = 0) and |z_n(0)| ≤ a_n.
- **(H) Isolation:** Y_n > max(a_n, g_n(Δ)) for every n. This is checked in arb: componentwise for n ≤ N_c, and for n > N_c through the strict tail inequality on T.

**Claim.** Let u ∈ C⁰(I; H⁵) ∩ C¹(I; L²) solve u_s = −(∂⁴+∂²)u + N(u) in the sine class, with z(0) := u(0) − ū(0) − P(0)ρ satisfying |z_n(0)| ≤ a_n. Then z(s) ∈ W for all s ∈ I, and |z_n(Δ)| ≤ g_n(Δ).

**Proof.**
1. *Galerkin system.* For N ≥ M, let u^N solve u^N_n′ = −λ_n u^N_n + N_n(u^N) for n ≤ N, with u^N_m ≡ 0 for m > N and u^N(0) = Π_N u(0). This is a polynomial ODE on ℝ^N, so a unique local solution exists and continues as long as it stays bounded.
   Set z^N = u^N − ū − Pρ, which has finitely many modes. Since N_n(ū + Pρ + z^N) = N_n(ū) + (J z^N)_n + Q_n(s, Pρ + z^N), and since r_n and (Dρ)_n involve only ū and P, every n ≤ N satisfies z^N_n′ = −λ_n z^N_n + F^N_n with F^N_n = (J z^N)_n + Q_n − r_n − (Dρ)_n.
2. *First exit.* We have |z^N_n(0)| ≤ a_n < Y_n, and only finitely many coordinates move. Let s* = sup{s ≤ Δ : z^N([0,s]) ⊂ W}; then s* > 0.
   On [0, s*], z^N ∈ W, so |F^N_n| ≤ F_n, and Duhamel gives |z^N_n(s*)| ≤ g_n(s*). Here g_n′(t) = e^{−λ_n t}(F_n − λ_n a_n) has constant sign, so g_n(s*) ≤ max(g_n(0), g_n(Δ)) < Y_n by (H), for all n ≤ N.
   If s* < Δ, continuity of the finitely many z^N_n contradicts maximality. Hence u^N exists on I, z^N(I) ⊂ W, and |z^N_n(Δ)| ≤ g_n(Δ).
3. *Compactness.* Uniformly in N and s: |u^N_n| ≤ B_n := sup_I|ū_n| + sup_I|(Pρ)_n| + Y_n, where B_n = Y_n = T n^{−6} for n > max(M, N_c).
   Let Ñ_n := sup{|N_n(v)| : |v_m| ≤ B_m}. It is finite because the convolution sums in N_n are dominated by sums of m^k B_m with k ≤ 3, and these converge since q = 6 > 4.
   Hence |∂_s u^N_n| ≤ λ_nB_n + Ñ_n. For each n the family (u^N_n)_N is bounded and equi-Lipschitz on I. Arzelà–Ascoli and a diagonal argument give a subsequence N_k with u^{N_k}_n → u^∞_n uniformly on I, for every n.
4. *The limit solves the mild system.* u^N_n(s) = e^{−λ_n s}u_n(0) + ∫₀ˢ e^{−λ_n(s−σ)}N_n(u^N(σ))dσ for N ≥ n.
   Each convolution series in N_n converges termwise and is dominated by the B-majorant series, so N_n(u^{N_k}(σ)) → N_n(u^∞(σ)). Also |N_n(u^N)| ≤ Ñ_n. Dominated convergence in σ then shows that u^∞ satisfies the mild equations, and the closed bounds pass to the limit: z^∞(I) ⊂ W and |z^∞_n(Δ)| ≤ g_n(Δ).
5. *Regularity.* For n > max(M, N_c): |u^∞_n| ≤ T n^{−6} and |N_n(u^∞)| ≤ Ñ_n ≤ C n^{−4} (N costs two derivatives in the IBP form).
   So u^∞_n is C¹ with |u^∞_n′| ≤ C n^{−2}. Since Σn^{10}n^{−12} < ∞ and Σn^{−4} < ∞, u^∞ ∈ C⁰(I; H⁵) ∩ C¹(I; L²), and it solves the PDE in L².
6. *Uniqueness.* Let δ = u − u^∞ and A = max over I of ‖v‖_∞, ‖w‖_∞, ‖W‖_∞ for both solutions (finite, since H⁵ ⊂ C³). Then
   ½(d/ds)‖δ‖² = −‖δ_αα‖² + ‖δ_α‖² + ⟨N(u) − N(u^∞), δ⟩,
   where N(u) − N(u^∞) = −∂(δ_α X) + δ_α W(u) + (1 + v^∞)∂⁻¹((w + w^∞)δ_αα − mean), with ‖X‖_∞ ≤ C_A.
   - Integrating by parts, ⟨−∂(δ_αX), δ⟩ ≤ C_A‖δ_α‖².
   - The second term is ≤ C_A‖δ_α‖‖δ‖.
   - ‖∂⁻¹(…)‖_∞ ≤ C_A‖δ_αα‖_{L¹}, so the third term is ≤ C_A‖δ_αα‖‖δ‖.
   With ‖δ_α‖² ≤ ‖δ‖‖δ_αα‖ and Young: (d/ds)‖δ‖² ≤ C‖δ‖². Since δ(0) = 0, Gronwall gives δ ≡ 0. So u = u^∞ satisfies the claim. ∎

**Remarks.**
- (i) Steps 1–5 also give *existence* in C⁰H⁵ ∩ C¹L² on I, from data in H⁵, and Step 6 gives uniqueness. Chaining sub-intervals and slabs therefore constructs the CDF solution from Γ0 on all of P1, with no external existence theorem.
- (ii) The intra-slab absorption, the Lohner step and M/block resizing are exact changes of representation, so they need no analysis.
- (iii) (H) is verified with a 3% uniform margin (μ). By the monotonicity of g_n, the sup over I of |z_n| is attained at s = 0 or s = Δ.
