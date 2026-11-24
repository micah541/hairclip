# Lemma F: far-field energies of the P2 error (step 4, T8)

Status. §2 (Lemma F) and §4 (Lemma J) are proved here for the **linear** error equation, with the modulation and residual forcing included. The energies are half-curve energies; the far-end boundary term is §3.1 (Lemma B and hypothesis (SE)). Lemma N3′ (§6.1) bounds the nonlinear terms of the exact graph operator and the modulation terms 𝒬 at every level j ≤ 3, with x-dependent majorants on the certified far boxes and without level 4 (`p2/t8se/n3sym.py`, `n3qsym.py`, `n3wcheck.py`, `n3far.py` → `camp.py --n3`). The algebraic core, the level-j identity (I_j), is machine-checked: `lfsym.py` (κ̃-coefficients, j ≤ 4) and `cforms.py` (generic coefficients β, b₂, c₀, j ≤ 3) run an Euler-operator test, and negative controls fail as they should. Every constant is an explicit sup or integral of base quantities. Certified (arb) values on the polished base, slabs 5–187: `cconst.py`, `merge.py`, `camp.py` → `constants/` (T8_REPORT §3). Float values: `fb.py`, `weights.py`, `lfnum.py`, `jcas.py`, `amp.py`. Inputs taken from T5/T7 without being re-derived are tagged **[T7]**.

Notation: x is base arclength from the tip, y := x − X0 with X0 = 2.2490776515 (as in Lemma S), and c := ¼ − ã.

## 1. Setting

**1.1 The far-field equation.** T5 §2.1 (signs confirmed in T1 `signs.py`) gives
e_τ = 𝓛e − Σᵢxᵢgᵢ − (1−ϑ)r + 𝒩 + 𝒬, with 𝓛e = −e'''' + b₂e'' + b₁e' + b₀e and
b₂ = −κ̃², b₁ = −5κ̃κ̃' − c q̃ − b̃t̃₁ + ϖ, b₀ = c − 3κ̃'² − 4κ̃κ̃''.
Inserting ϖ = −κ̃κ̃' + Q + c(q̃ − x) + b̃t̃₁ (T1 §3) gives

 b₁ = −c y + β, β := Q − cX0 − 6κ̃κ̃' + A, Q := ∫₀ˣκ̃'², b₀ = c + c₀, c₀ := −3κ̃'² − 4κ̃κ̃''.

Here A is the residual part of the transport (T6 §5: ϖ = ϖ* + A, A' = κ̃r, A(0) = 0); the float study of T8_REPORT §2 used A = 0. The terms −cq̃ and −b̃t̃₁ cancel exactly in b₁, so b₁ involves neither the closure q̃ nor the turning t̃₁.

On the far field y ≥ y_a ≥ 14.001 we have x ≥ 16.25, so ζ_∞ = 1 and v := e' is T5's E_I variable. Differentiating once:

 **v_τ = −v'''' − c y v' + (βv + b₂v' + c₀e)' + f**, f := −x₁g₁' − x₂g₂' − x₃(gˢ₃ + gʳ)' − r' + (𝒩+𝒬)'.  (1.1)

Here g₃ = Ṽ/n = gˢ₃ + gʳ + r/n with gˢ₃ := (cp̃ + b̃ñ)/n and gʳ := −κ̃''/n, by T5's definition r := Ṽ + κ̃'' − cp̃ − b̃ñ. So x₃g₃ + (1−ϑ)r = x₃(gˢ₃ + gʳ) + r: the clock factor drops out of the residual forcing. The derivation of (1.1) is ∂ₓ(𝓛e) = −v'''' + (b₂v')' + (−cyv + βv)' + cv + (c₀e)' = −v'''' − cyv' + (βv + b₂v' + c₀e)'.

We write w_m := ∂ₓᵐv for m ≥ 0 and w₋₁ := e. Lemma F treats f without the nonlinear part (𝒩+𝒬)'; §6.1 adds 𝒩 and 𝒬 at levels 0–3.

**1.2 Parities.** The base is exactly symmetric about both ends (T1: κ̃' = κ̃''' = 0 at x = 0 and x = ℓ̃). We view the half curve as one half of the closed doubled curve x ∈ [0, 2ℓ̃]. About the far end x = ℓ̃:
- κ̃, b₂, c₀ and e are even, and w_m has parity (−1)^{m+1};
- β^{(m)} has parity (−1)^{m+1} for m ≥ 1;
- b₁ + ℓ̃_τ is odd up to the constant D := b₁(ℓ̃) + ℓ̃_τ. This holds because κ̃κ̃', Q − Q(ℓ̃) and y − y_ℓ are odd, with y_ℓ := ℓ̃ − X0. With T1's length equation ℓ̃_τ = cℓ̃ − Q(ℓ̃) + ε_ℓ and b₁(ℓ̃) = Q(ℓ̃) − cℓ̃ + A(ℓ̃) one gets **D = A(ℓ̃) + ε_ℓ**, |D| ≤ 6.9e-13 (T6: sup|A(ℓ̃)| ≤ 6.7e-13, sup|ε_ℓ| ≤ 2.3e-14). T6's scalar D_end = ϖ(ℓ̃) + ℓ̃_τ equals D + cq̃(ℓ̃) + b̃t̃₁(ℓ̃); it is not the transport coefficient (the closure terms cancel in b₁, §1.1).

**1.3 Weights.** For each level j ∈ {0..4} take a ramp [y_a,j, y_b,j] and set
 T_j(x,τ) := 28 η_j(y)² Y_ℓ(x,τ)^{2j−2} Ω(Y_ℓ/Y_h(τ)), where:
- η_j = S(t), S(u) = u⁴(35 − 84u + 70u² − 20u³), t = ln(y/y_a)/ln(y_b/y_a), clipped to [0,1]. The T5 design is y_a, y_b = 14, 28 for all j (the Lemma S ramp).
- Y_ℓ := y_ℓ Φ(y/y_ℓ) is a far-end flattening. Φ(s) = s for s ≤ ½, and Φ' = 1 − ψ(2s−1) with ψ(u) = u⁵(21/4 − 49u/8 + 15u²/8). Then Φ' is odd about s = 1, so Y_ℓ is even about ℓ̃ on the doubled curve, and Y_ℓ = y on the half of the curve nearer the tip.
- Ω ≡ 1 for the *global* energies. The *localized* energies of Remark 6.3 use an outer factor that follows the characteristics.

The properties used are: (W1) T_j ∈ W^{4,∞}, since η² is C³ at y_b with T'''' bounded; (W2) T_j = 0 for y ≤ y_a and T_j ≍ (y − y_a)⁸ there, so T'²/T and T''²/T are bounded; (W3) T_j is even about ℓ̃. (The jets are checked against finite differences in `weights.py`; the only mismatch is at the C³ point y = y_b.)

Regions: the ramp R_j := {y_a,j < y < y_b,j} and the core C_j := {y ≥ y_b,j}, the latter including the flattened far end.

**1.4 Energies and data.** The energies are the half-curve integrals F_j := ∫_{x*}^{ℓ̃(τ)} T_j w_j², the dissipations D_j := ∫T_j w_{j+2}², and N_{j+1}² := ∫T_j w_{j+1}² (same domain). The ramp data are the E_I norms of T5 restricted to the ramp, V_m := (∫_{R} W (∂ᵐv)²)^{1/2} ≤ T5's V_m, with W := (1+x²)^{−2}. The cost of a stronger weight is ∫T_j w_j² vs ∫W w_j², i.e. y^{2j+2} on the ramp. This conversion factor drives the amplification (T8_REPORT §2).

**1.5 Regularity.** κ̃ ∈ C⁴ with κ̃^{(5)} ∈ L^∞ (T1), so e ∈ W^{7,∞} and v ∈ W^{6,∞} on the far field (T5 §2.1). The identity (I_j) below involves w_{≤ j+2} and κ̃^{(≤5)} only. It is proved for smooth data and extends by density, because every term is continuous in the norms that are available. The time derivative is taken a.e.; F_j is absolutely continuous.

## 2. Statement

**Lemma F.** Let j ∈ {0,…,4}, let the weights satisfy (W1)–(W3), and let the ramps satisfy y_a,i < y_b,j for i < j, so that T_i > 0 on C_j. Let F_j be the half-curve energies of §1.4 and assume the symmetric-end hypothesis (SE) of §3.1. Along the linear part of (1.1), for a.e. τ:

 F_j' ≤ −λ_j F_j − θ_j D_j + 2√F_j g_j + h_j,  (F)

 g_j = Σ_{i<j} ½c_ji √F_i + |x|γ_j + ½a_j Ȳ,
 h_j = Σ_{i<j} d_ji F_i + 𝓘_j + e_j Ȳ² + δ_s⁻¹(𝓡_j + |x|γʳ_j)² + 𝓡̂_j²,

with the ramp data

 𝓘_j = k_jj V_j² + Σ_{i<j}(k_ji V_i² + r_ji V_j V_i) + re_j V_j Ȳ  (standard ramp),
 𝓘_j = kD_j D_{j−2} + Σ_{i<j} sD_ji (D_{j−2}F_i)^{1/2} + …  (nested ramp: R_j ⊂ C_{j−2}).

The constants are defined as follows. Fix a = 0.01, ν = 10⁻³, δ_P ∈ (0, 1.5], δ_s ∈ (0, ½] and δ₆ = 0.1. Let n_j be the number of (j+1,·)-pairs in (I_j), plus 1 for the forcing, recall P_j, q_j, p_ji from §3.2, and let ε_j := |D|/m_j be the far-end constant of Lemma B (§3.1).

- q̂_j := ∂_τT_j + q_j + ½P_j'' + P_j²/(4δ_P T_j) + (a n_j/4 + ε_j/2)(T_j'')⁺ + ν|T_j''|.
- ΔF_j := a n_j/8 + ε_j(2/y_ℓ + 5/4), σ_j := c + sup_{C_j} q̂_j/T_j + ΔF_j and λ_j := c − σ_j.
- θ_j := 2 − δ_P − a n_j/2 − δ₆[j=4] − δ_s − ε_j. (The design opt4n has δ_P = 1.5 at j = 3; it needs δ_s < ½ there, and the certified run uses δ_s = ¼.)
- k_jj := sup_{R_j} (q̂_j + (λ_j + ΔF_j) T_j)⁺/W. In the nested case, kD_j := sup_{R_j} (q̂_j + (λ_j + ΔF_j)T_j)⁺/T_{j−2}. (Earlier versions omitted ΔF_jT_j here; it is the ramp share of the F_j-part of the Young steps, §3.4.)
- c_ji := sup_{C_j} |p_ji|/√(T_jT_i). The e-coupling adds 2A₂ to c_j0; A₂ would suffice.
- r_ji := sup_{R_j}|p_ji|/W. In the nested case, sD_ji := sup_{R_j}|p_ji|/√(T_{j−2}T_i).
- (j+1,0)-pairs p: d_j0 += sup_{C_j} p²/(T_jT_0)/(2a) and k_j0 += sup_{R_j} p²/(T_jW)/(2a).
- (j,e)-pairs p: a_j := ε̂₀(∫_{C_j}p²/T_j)^{1/2}, A₂ := (∫_{C_j}(p²/T_j)H²)^{1/2}, and re_j := (sup_{R_j}|p|/W · ∫_{R_j}|p|ê²)^{1/2}.
- (j+1,e)- and (j+2,e)-pairs p: e_j += [2ε̂₀²∫_{C_j}p²/T_j + ∫_{R_j}(p²/T_j)ê²]/κ and d_j0 += 2∫_{C_j}(p²/T_j)H²/κ, with κ = 2a for (j+1,e) and κ = 4δ₆ for (j+2,e).
- γ_j := (Σᵢ∫T_j(∂^{j+1}gᵢˢ)²)^{1/2} and γʳ_j := (∫T_j(∂^{j−1}gʳ)²)^{1/2}, where gˢ := (p̃/‖χp̃‖_ρ, ñ/‖χñ‖_ρ, (cp̃+b̃ñ)/n) and gʳ := −κ̃''/n with n = ‖χṼ‖_ρ (T5 normalisation). For j = 0, ∂^{−1}gʳ := −κ̃'/n.
- 𝓡_j := (∫T_j(∂^{j−1}r)²)^{1/2}, which is T5's R^F_j up to the Y_ℓ flattening, with ∂^{−1}r := −∫_x^{ℓ̃}r.
- 𝓡̂_j² := (8/a)∫(T_j'²/T_j)ρ_j² + ν⁻¹∫|T_j''|ρ_j², with ρ_j := ∂^{j−1}(r + x₃gʳ).

The e-bound uses three functions:
- ê(x) := K_e + Ĩ(x)/√μ_e, with Ĩ(x)² := ∫_{x*}^x W^{−1} and |e(x*)| ≤ K_eȲ **[T7]**;
- ε̂₀ := ê(x₁), where x₁ is the start of C₀;
- H(x)² := ∫_{x₁}^x T₀^{−1} on C₀, with H := 0 below x₁.

**Validity.** All sups are finite by (W2). If λ_j > 0 and θ_j > 0, then (F) holds with these constants (any λ ≤ λ_j may be used in (F) and in k_jj). In the nested case, R_j ⊂ C_{j−2} and T_i > 0 on R_j are required.

**Hand-off coverage.** The window energies of Lemma S (cutoff S(ln(y/14)/ln 2) on [14,28]) are dominated by F_j whenever η_j ≥ S(ln(y/14)/ln 2), for instance for every log-ramp [14, y_b] with y_b ≤ 28.

## 3. Proof

**3.1 Half-curve energies and the far-end boundary term.** T_j vanishes to order 8 at y_a (W2), so nothing happens at x*. At the moving far end,

 F_j' = ∫_{x*}^{ℓ̃}(∂_τT_j)w_j² + 2∫_{x*}^{ℓ̃}T_j w_j ∂_τw_j + ℓ̃_τT_j(ℓ̃)w_j(ℓ̃)².

Insert (I_j); its total-derivative part 𝔉_j leaves 𝔉_j(ℓ̃). The flux of the transport term 2T_jw_j·b₁w_{j+1} = (T_jb₁w_j²)' − (T_jb₁)'w_j² is T_jb₁w_j². Hence

 F_j' = ∫_{x*}^{ℓ̃}[(∂_τT_j)w_j² + density of (I_j)] + 𝔅^D_j + 𝔅^S_j,
 𝔅^D_j := D·T_j(ℓ̃)w_j(ℓ̃)², D = b₁(ℓ̃) + ℓ̃_τ = A(ℓ̃) + ε_ℓ (§1.2),

and 𝔅^S_j is the sum of the fluxes at ℓ̃ of the fourth-order, b₂-, c₀-, forcing and X_j-terms. 𝔅^S_j contains no b₁; it is a bilinear form in the boundary jets (e, w₀, …, w_{j+3})(ℓ̃) with coefficients from T_j^{(≤3)}(ℓ̃) (T_j' = T_j''' = 0 there by (W3)) and the coefficient jets at ℓ̃.

**Lemma B (transport defect).** With m_j := (Φ(½)/Φ(1))^{2j−2} for j ≥ 2, m₀ = m₁ = 1 (Φ(½) = ½, Φ(1) = 0.8828125) and ε_j := |D|/m_j,

 |𝔅^D_j| ≤ ε_j[(2/y_ℓ + 5/4)F_j + D_j + ½∫(T_j'')⁺w_j²].

*Proof.* P := {y_ℓ/2 ≤ y ≤ y_ℓ} has length d = y_ℓ/2 and lies in C_j. On P, T_j ≥ m_jT_j(ℓ̃): T₀ = 28Y_ℓ⁻² decreases, T₁ = 28 is constant, T_j increases for j ≥ 2. For f = w_j, the trace inequality f(ℓ̃)² ≤ d⁻¹∫_P f² + 2‖f‖_{L²(P)}‖f'‖_{L²(P)} (average f(ℓ̃)² = f(z)² + ∫_z^ℓ̃(f²)' over z ∈ P) gives T_j(ℓ̃)w_j(ℓ̃)² ≤ m_j⁻¹(F_j/d + 2√F_jN_{j+1}) ≤ m_j⁻¹((2/y_ℓ + 1)F_j + N²_{j+1}). Finally N²_{j+1} ≤ D_j + ¼F_j + ½∫(T_j'')⁺w_j² (S3). ∎

Lemma B is what the constants ΔF_j, θ_j and q̂_j of §2 contain. With |D| ≤ 6.9e-13, ε_j ≤ 2.1e-11: negligible, but included in the certified run.

**(SE) Symmetric end.** 𝔅^S_j ≤ 0.

(SE) holds with 𝔅^S_j = 0 when the error, the base coefficients and the residual are symmetric about the far-end normal line. This is the closed-curve argument: write 2F_j as an integral over the doubled far field [x*, 2ℓ̃ − x*], on which (1.1) holds with the even/odd continuations of the coefficients (§1.2). Under x ↦ 2ℓ̃ − x every density is even by §1.2 and (W3), and the transport/∂_τT pair is mirror-symmetric (b₁(2ℓ̃−x) = −b₁(x) − 2ℓ̃_τ, ∂_τ[T(2ℓ̃−x)] = (∂_τT)(2ℓ̃−x) + 2ℓ̃_τT'(2ℓ̃−x)), so no flux appears at ℓ̃. The evenness of T_j is what makes the fourth-order fluxes vanish; with the unflattened power weight, 2T'w_jw_j'' + T'''w_j² would survive at ℓ̃ for odd j, and a positive term 2|T'|w_j'² for j = 0.

*Remark (far end off the axis).* The base is symmetric about its own far-end normal line N(τ), since κ̃ is exactly even about ℓ̃. The true curve is symmetric about the axis. By T6, N is off the axis by |q̃(ℓ̃)| ≤ 4.6e-9 and tilted by |t̃₁(ℓ̃)| ≤ 2.4e-13. The offset grows like ℓ̃ (q̃(ℓ̃)/ℓ̃ ≈ 1.3e-14 on slabs ≥ 100), and the turning defect θ̃(ℓ̃) + 3π/2 ≤ 1.65e-14 is constant. This is the signature of a fixed physical rigid rotation of about 1.6e-14 rad. So (SE) fails only through this rigid mode; by §1.1 the offset does *not* enter the transport term. There are two ways to restore (SE):
- (a) Compare the true curve with the base rotated by that angle. Rotations about the origin are symmetries of the rescaled flow up to the axial drift, which changes r by at most |b̃|·1.7e-14 ≤ 1.1e-15 (T6 would enclose this). This is a choice in P2's representation, for T5/T9; T4 must check the start tolerance, since the tip moves by about 4e-14.
- (b) Bound 𝔅^S_j directly, as (jets of the rigid mode at ℓ̃, of size 1.7e-14·ℓ̃(19.4/ℓ̃)^m) × (traces of w_{≤j+3} at ℓ̃). This needs pointwise far-end control of w_{j+3}, which the level-j energies do not give.

The mode's own energies are ≤ 1e-20 relative to F_j, so the effect is nil; but (SE) is an input of Lemma F, not a consequence.

**3.2 The identity (I_j)** (machine-checked: `lfsym.py`, output `lf_forms.pkl`). With T = T_j,

 2T w_j ∂ʲ[v_τ − f] ≐ −2T w_{j+2}² + P_j w_{j+1}² + q_j w_j² + Σ_{i=−1}^{j−1} p_ji w_j w_i + X_j,

where ≐ means equality modulo a total x-derivative, and
- P_j = 4T'' + 2κ̃²T. The second summand is the anti-diffusive −2b₂T term.
- q_j = −T'''' + (cy − β)T' + (1−2j)cT + (2j+1)β'T + (Tb₂)'' − (j+1)(Tb₂')' + j(j+1)Tb₂'' + 2c₀T.
- p_ji = 2C(j+1,j+1−i)Tβ^{(j+1−i)}·[i ≥ 0, j+1−i ≥ 2] + 2C(j+1,j+2−i)Tb₂^{(j+2−i)}·[1 ≤ i ≤ j−1] + 2C(j+1,j−i)Tc₀^{(j−i)}·[i ≤ j−1].
- X_j keeps every coefficient derivative at κ̃^{(≤5)}. It replaces three high-order terms by integrations by parts onto (Tw_jw_i):
  - j=3: the c₀^{(4)}-term 2Tc₀^{(4)}w₃e becomes −2c₀'''(T'w₃e + Tw₄e + Tw₃w₀).
  - j=4, β^{(5)} and c₀^{(4)} terms: these become −(2β^{(4)} + 10c₀''')(T'w₄w₀ + Tw₅w₀ + Tw₄w₁).
  - j=4, c₀^{(5)}-term: this becomes 2c₀'''(T''w₄e + 2T'w₅e + 2T'w₄w₀ + Tw₆e + 2Tw₅w₀ + Tw₄w₁).

Derivation:
- From −2Tw_jw_j'''', three integrations by parts give −2Tw_{j+2}² + 4T''w_{j+1}² − T''''w_j².
- From ∂ʲ(−cyw₁) = −cyw_{j+1} − jcw_j: the term c(yT' + (1−2j)T)w_j².
- The β-, b₂- and c₀-terms are expanded by Leibniz. The β and b₂ terms with m = 0, 1 (and m = 2 for b₂) are folded into q_j and P_j by one or two integrations by parts.

The Euler test checks the whole identity for arbitrary T, κ̃ and Q with Q' = κ̃'². It reports κ̃-derivatives ≤ 5 in all coefficients. Pairs present: j=0: (0,e); j=1: (1,e),(1,0); j=2: (2,e),(2,0),(2,1); j=3: (3,e),(3,0..2),(4,e); j=4: (4,e),(4,0..3),(5,e),(5,0),(6,e).

**3.3 Bounds.**
- **(S1)** Exactly, ∫P w_{j+1}² = −∫P w_j w_{j+2} + ½∫P''w_j². Then |∫Pw_jw_{j+2}| ≤ δ_P D_j + (4δ_P)⁻¹∫(P²/T)w_j².
- **(S2) (j,i), 0 ≤ i < j.** On C_j the weights are regular: ∫|p||w_jw_i| ≤ c_ji√(F_jF_i). On R_j use the ramp data: ≤ r_ji V_jV_i. In the nested case use instead ≤ (∫_{R_j}T_{j−2}w_j²)^{1/2}(∫(p²/T_{j−2})w_i²)^{1/2} ≤ sD_ji√(D_{j−2}F_i).
- **(S3) (j+1,i), i ∈ {0,e}.** First N_{j+1}² = −∫Tw_jw_{j+2} + ½∫T''w_j² ≤ D_j + ¼F_j + ½∫(T'')⁺w_j². Then |∫pw_{j+1}w_i| ≤ N_{j+1}M_i ≤ (a/2)N_{j+1}² + M_i²/(2a), where M_i² = ∫(p²/T)w_i² is split over C_j (giving F₀) and R_j (giving V₀ or ê).
- **(S4) (j+2,e).** |∫pw_{j+2}e| ≤ δ₆D_j + (4δ₆)⁻¹∫(p²/T)e².
- **(S5) e.** Write e(x) = e(x*) + ∫_{x*}^x v. Then |e| ≤ ê(x)Ȳ on R₀, using μ_eV₀² ≤ Λ_e ≤ Ȳ². On C₀, |e| ≤ ε̂₀Ȳ + H(x)√F₀, by Cauchy–Schwarz with T₀ and because ê is increasing. Hence ∫_{C_j}|p||w_j||e| ≤ √F_j(a_jȲ + A₂√F₀) by Minkowski, and ∫_{R_j}|p||w_j||e| ≤ re_jV_jȲ. For j ≥ 1 the A₂-part goes into c_j0. For j = 0 it is a self-coupling, and λ₀ is lowered by 2A₂ (≤ 10⁻⁸; the certified run does this).
- **(S6) Forcing.**
  - Residual and rough modulation: two integrations by parts give −2∫(Tw_j)''ρ_j = −2∫(T''w_j + 2T'w_{j+1} + Tw_{j+2})ρ_j (their fluxes at ℓ̃ belong to 𝔅^S_j). With ‖ρ_j‖_T ≤ 𝓡_j + |x|γʳ_j, this is ≤ δ_sD_j + δ_s⁻¹(𝓡_j + |x|γʳ_j)² + (a/2)N_{j+1}² + 𝓡̂_j² + ν∫|T''|w_j². (For j=0, ∂^{−1}r = S₁ − S₁(ℓ̃) vanishes at ℓ̃.)
  - Smooth modulation: Cauchy–Schwarz in x and in i gives ≤ 2|x|√F_jγ_j. For this smooth forcing the two-integrations route is far worse: it pays the y⁴ weight mismatch between D_j and F_j.

**3.4 Collecting.**
- All w_j²-terms are ∫q̂_jw_j² + ΔF_jF_j, where ΔF_j is the F_j part of the N_{j+1}² bounds in (S3), (S6) and Lemma B. On C_j they are ≤ (sup_{C_j}q̂_j/T_j + ΔF_j)∫_{C_j}T_jw_j² = −λ_j∫_{C_j}T_jw_j². On R_j, write ∫_{R_j}(q̂_j + ΔF_jT_j)w_j² = −λ_j∫_{R_j}T_jw_j² + ∫_{R_j}(q̂_j + (λ_j + ΔF_j)T_j)w_j²; the last integral is ≤ k_jjV_j². In the nested case it is ≤ kD_jD_{j−2}, using T_{j−2} ≡ T_{j−2}|_{core} on R_j.
- The D_j-coefficient is −2 + δ_P + an_j/2 + δ₆[j=4] + δ_s + ε_j = −θ_j.

Together with 𝔅^D_j (Lemma B) and 𝔅^S_j ≤ 0 (SE), this proves (F). ∎

## 4. Lemma J (time-weighted ramp data)

**Hypotheses.** T5's closure is run with a reserved fraction κ_𝔇 of its dissipation surplus. Concretely:
- ½Λ_e' ≤ −ϱΛ_e − κ_𝔇𝔇 + a_cY + Ξ'², where 𝔇 = (θ₀/4)B² + (ε/8)D₄² + (μ_e/8)V₂² + (μ_eε_I/8)V₄² and Ξ'² = Ξ²/(1−κ_𝔇).
- Ȳ' = −ϱȲ + Φ with Φ = a_c + √ϱΞ', √Λ_e ≤ Ȳ ≤ 1 and Ȳ ≥ Ξ'/√ϱ **[T7]**.

As written in T5 (Y), the surplus is consumed by the residual Young steps, so κ_𝔇 = 0 there; we use κ_𝔇 = ½.

We also use:
- (I0) and (I2) of T5 §2.4 with |x| ≤ C_xA;
- the interpolations W₁² ≤ AB + c₁A², W₃² ≤ BD₄ + c₁B², V₁² ≤ V₀V₂ + c₁W V₀², V₃² ≤ V₂V₄ + c₁W V₂² **[T7/S3]**.

Let 0 < λ < 2 min(μ₁, μ₂), ω(s) := e^{−λ(τ−s)}, and J_X := ∫_{τ0}^τ ωX².

**Conclusion.**
1. ∫ω𝔇 ≤ 𝔅 := κ_𝔇⁻¹[½ + (ϱ̄ − λ/2)⁺/λ].
2. J_A ≤ 1/λ, J_B ≤ 4𝔅/θ₀, J_{D4} ≤ 8𝔅/ε, J_{W1} ≤ √(J_AJ_B) + c₁J_A, J_{W3} ≤ √(J_BJ_{D4}) + c₁J_B.
3. (μ₁ − λ/2)J_{V0} + ½J_{V2} ≤ K + C_xγ₀√(J_AJ_{V0}), with K := ½ω(τ0)/μ_e + a_fJ_{W1} + b_fJ_B.
4. J_{V4} ≤ 2ω(τ0)/(μ_eε_I) + M(λ)(c₃√J_{D4} + c₂√J_{W3} + c₁ᴵ√J_B + c₀ᴵ√J_{W1} + C_xγ₂√J_A)², with M(λ) := min_κ 4[(1+Ψ₁/κ)² + (Ψ₁κ + 2Ψ₁√c₁W + 2c₁W)²/(4(μ₂ − λ/2))].
5. J_{V1} ≤ √(J_{V0}J_{V2}) + c₁W J_{V0} and J_{V3} ≤ √(J_{V2}J_{V4}) + c₁W J_{V2}. The Λ-routes J_{V2} ≤ 8𝔅/μ_e and J_{V4} ≤ 8𝔅/(μ_eε_I) are also valid.
6. Residual terms. Suppose the majorant satisfies Ȳ ≥ Ȳ_* (the certified run uses Ȳ_* = 10⁻¹²; the budget is about 2.3·10⁻¹²). The term 𝓡₁(V₂ + 2Ψ₁V₁ + Ψ₂V₀) of (I0) is bounded by δ(V₀² + V₂²) + 𝓡₁²C_res/(4δ), where C_res = (1+Ψ₁)² + (Ψ₁(1 + 2√c₁W) + Ψ₂)² and Ψ₂ := 2c₁W (this uses V₁ ≤ ½(V₀ + V₂) + √c₁W V₀). Hence item 3 holds with μ₁ → μ₁ − δ, ½J_{V2} → (½ − δ)J_{V2}, and K += C_res(𝓡̄₁/Ȳ_*)²/(4δλ). In (I2), 𝓡ᴵ_s enters X additively; it adds 𝓡̄ᴵ_s/(Ȳ_*√λ) to √∫ωX². The run uses δ = 10⁻³.

**Proof.**
1. Multiply the closure by ω and integrate, using ½∫ωΛ' = ½[ωΛ] − (λ/2)∫ωΛ. Since Y ≤ Ȳ and Ξ'² ≤ √ϱΞ'Ȳ, we have a_cY + Ξ'² ≤ ΦȲ. Write ∫ωΦȲ = ∫ω(½(Ȳ²)' + ϱȲ²) = ½Ȳ(τ)² − ½ω(τ0)Ȳ(τ0)² + ∫ω(ϱ − λ/2)Ȳ². This gives κ_𝔇∫ω𝔇 ≤ ½Ȳ(τ)² + ∫ω(ϱ − λ/2)(Ȳ² − Λ_e). The last integral is ≤ (ϱ̄ − λ/2)⁺/λ.
2. The pointwise bound Λ_e ≤ Ȳ² ≤ 1 gives J_A, and item 1 gives J_B and J_{D4}. The remaining bounds follow from the interpolations and Cauchy–Schwarz in time.
3. Multiply (I0) by ω.
4. Multiply (I2) by ω. Use V₃ ≤ √(V₂V₄) + √c₁W V₂ ≤ V₄/(2κ) + (κ/2 + √c₁W)V₂, then Young: XαV₄ ≤ ¼V₄² + α²X² and Xβ̂V₂ ≤ (μ₂ − λ/2)V₂² + β̂²X²/(4(μ₂ − λ/2)). Finally apply Minkowski in time to √∫ωX². ∎

For the T7 inputs **[T7]** (S1 reference values / T5 defaults): μ_e = 0.02, ε = 10⁻⁵, ε_I = 3·10⁻⁵, θ₀ = 0.02, ϱ̄ = 0.03, c₁ = 0.5804, a_f = 0.0696, b_f = 0.0394, μ₁ = 0.1552, μ₂ = 0.6555, c^I_{0..3} = (0.425, 0.424, 0.526, 0.442), c₁W = 0.8477, Ψ₁ = 0.39, C_x = 4, γ₀ = 0.05, γ₂ = 0.01, κ_𝔇 = ½. With λ = 0.244 this gives, per unit Ȳ²:
- J₀..J₄ = (205, 255, 31.7, 1.24·10⁴, 4.82·10⁶);
- 73% of J₄ is the c₃²J_{D4} = c₃²·8𝔅/ε term, i.e. the 1/√ε loss.

## 5. Corollary C (comparison and hand-off)

**C1 (supersolution).** Suppose F' ≤ −λF + 2g√F + h with λ, g, h ≥ 0. Define H := e^{−∫λ}F(τ0) + ∫e^{−∫_s^τλ}h and G := ∫e^{−½∫_s^τλ}g. Then √F ≤ √H + G.

Proof. Z := (√H + G)² satisfies Z' = −λZ + 2g√Z + h + Gh/√H ≥ −λZ + 2g√Z + h. On any interval where F > Z > 0 we get (F − Z)' ≤ (−λ + g/√Z)(F − Z). Grönwall from the last time with F ≤ Z gives F ≤ Z; add δ > 0 to F(τ0) and let δ → 0. ∎

**C2 (triangular system).** Apply C1 level by level, with √F_i ≤ √H_i + G_i for i < j inside g_j and h_j, and with ∫ωh_j bounded through Lemma J:
- ∫ω𝓘_j ≤ Σ k_jm J_m + Σ r_ji √(J_jJ_i) + re_j √(J_j/λ);
- in the nested case, ∫ωD_{j−2} ≤ θ_{j−2}⁻¹[H_{j−2}^{(λ)} + 2 sup(g_{j−2}√F_{j−2})/λ]. Here ω and H_{j−2}^{(λ)} use the rate λ = min(λ_j, λ_{j−2}), since the term (λ − λ_{j−2})⁺∫ωF_{j−2} must vanish.

`amp.py` implements this with per-slab sups (sup over θ = 0, ½, 1) and the running sup of k_jm.

**C3 (windows, item 3).** z_j = (y∂_y)ʲv = Σ_m S(j,m) y^m w_m (Stirling), so √E_j^{[y,1.5y]}[err] ≤ κ_ch Σ_m S(j,m)√F_m. For E₃ this is √F₃ + 3√F₂ + √F₁. Here κ_ch is T9's chart factor (≤ 1.02).

**C4 (pointwise level 3, item 4).** For I = [y, 1.5y]:
 sup_I|w₃|² ≤ ∫_Iw₃²/|I| + 2‖w₃‖_{L²(I)}‖w₄‖_{L²(I)},
 with ‖w₃‖² ≤ F₃/(28y⁴) and either ‖w₄‖² ≤ F₄/(28y⁶) or ‖w₄‖² ≤ D₂(τ*)/(28y²).

**C5 (good time).** Integrating (F) at j = 2 over a window W = [τ_a, τ_b] with ω ≡ 1 gives θ₂∫_W D₂ ≤ F₂(τ_a) + ∫_W(2g₂√F₂ + h₂ + inflow). So there is τ* ∈ W with D₂(τ*) ≤ |W|⁻¹∫_W D₂. Lemma O does not make this selection: its H4, ∫(h₀⁽⁵⁾)² ≤ M² per cylinder, is a hypothesis on the data at t₁ = t(τ*), so step 4 must choose τ* (review A1). The bound is quadratic in Ȳ and is dominated by F₂(τ_a)/(θ₂|W|): with the hand-off sub-window W′ (|W′| = 1.46·10⁻⁴) it is D₂(τ*) ≤ 1.25·10⁻⁸ at Ȳ = 10⁻¹², and with one whole slab (|W| = 0.15) it is 1.26·10⁻¹¹ (`camp.py --window 180-180`). On a cylinder support S (y ≥ y₀ ≥ 22) the localized dissipation obeys ∫_S w₄² ≤ D_S/(28y₀²) with D_S := ∫_S T₂w₄² ≤ D₂. A localized energy ∫T₂ψ²w₂² (ψ = 1 on S) does not improve the leading term: its value at τ_a is bounded only by F₂(τ_a), and its inflow (the ramp, or transport through an inner cut-off) is that of F₂. So the window length, not the localization, decides. `p2/t8se/h4loc.py` evaluates H4 per cylinder (error part L⁷·1.001²·D₂(τ*)/(28y₀²), base part from the T6 hulls of κ̃^{(≤3)}). With W′, H4 fails by ×16–279 on the cone cylinders (r_c = 200–380) and by ×4.95 on the far cylinder centred at 376. Even the Ȳ-free part, 2.6·10⁻¹⁰ (the review's 2.9·10⁻⁹ also contains the part linear in Ȳ), fails by ×5.8 at r_c = 380. With one slab (slab 179 or 180), the worst cylinder (cone r_c = 370–380) has ∫(h₀⁽⁵⁾)² ≤ 0.31·M²; the far cylinders have ≤ 0.0052·M²; and the Ȳ-free part is ≤ 2.7·10⁻¹³, i.e. ≤ 0.6% of M². Item 4 then uses the same τ*, and the binding requirement becomes item 3 (Ȳ ≤ 1.79·10⁻¹¹). H1–H4 must then hold over the whole window slab. By Lemma N3′, level 4 is not needed for the nonlinear bootstrap either.

## 6. Remarks

**6.1 Nonlinear and modulation terms, levels 0–3 (Lemma N3′): exact operator, x-dependent majorants, no level 4.**

This replaces the earlier Lemma N3. That lemma expanded ∂⁴(𝒩+𝒬) fully, which produced the trilinear classes w₃²w₄ and w₃³, and it used T5's box with placeholder majorants a₀ = 20η_F, a₁ = 40(η_F/x + |κ̃|). On base v2, T7's Hessian sups grow with x (x²n₂ ≤ 6.2·10¹⁴ at slab 180), so the placeholders were not justified. N3′ uses the exact operator, has no trilinear classes, covers every level j ≤ 3 and the modulation part 𝒬, and is evaluated on the certified far boxes.

*Operator.* Let α := 1 − κ̃e and J² := α² + e'². At the base modulation and with the clock fixed (T7 `nl7.py`), the equation is e_τ = F(e) with

 F(e) = (J/α)[−κ_ss + c(−q̃e' + (p̃+e)α)/J + b̃(−t̃₁e' + ñα)/J] − Ṽ + e'(ϖ − eθ̃_τ)/α.

F = A₄e'''' + A₃e''' + A₂ with A₄ = −J⁻⁴ and A₃, A₂ functions of z := (e, e', e''). The nonlinear remainder is

 𝒩 := F(e) − F(0) − DF(0)e = Ãe'''' + B̃e''' + C̃,  Ã := 1 − J⁻⁴, B̃ := A₃(z) − A₃(0), C̃ := A₂(z) − A₂(0) − ∇A₂(0)·z.

*Cancellation: no coefficient grows with x.* Insert the base identities

 ϖ = −κ̃κ̃' + Q + c(q̃ − x) + b̃t̃₁ + A,  θ̃_τ = −κ̃''' − κ̃²κ̃' + κ̃Q − cxκ̃ + ρ_θ with ρ_θ = r' + κ̃A,  Ṽ = −κ̃'' + cp̃ + b̃ñ + r.

The quadratic part of the tangential term is ee'(κ̃ϖ − θ̃_τ)/α = ee'(cκ̃q̃ + b̃κ̃t̃₁ + κ̃''' − r')/α. Its first two parts cancel the quadratic parts −cq̃κ̃ee'/α and −b̃t̃₁κ̃ee'/α of the normal velocity. The remaining c- and b̃-terms are linear in e. Hence

 𝒩 = 𝒩_K + (κ̃''' − r')ee'/α,  where 𝒩_K is the nonlinear part of −(J/α)κ_ss.

So 𝒩 depends on none of c, b̃, p̃, q̃, ñ, t̃₁, x, Q and A. Its coefficients are functions of κ̃, …, κ̃''' and r' only, and none of them grows with x. `n3sym.py` checks A₄ = −J⁻⁴, F(0) = −r and the independence of c, …, Q symbolically; A drops out of the generated forms. T7's growth of n₂ comes from the interval hull of the cancelling sum, with q̃, x and Q taken as independent intervals; it is not a property of 𝒩.

*Modulation part.* For the terms 𝒬, with the true clock and modulation (T5 §2.1: c → (1−ϑ)c − x₁, b̃ → (1−ϑ)b̃ − x₂, κ_ss → (1−ϑ)κ_ss), the exact split is

 𝒩 + 𝒬 = (1−ϑ)𝒩 + ϑe'''' + ϑκ̃²e'' − ϑc₀e + ϑP + g·e' − (ϑc + x₁)e,

 P := (κ̃''' − r')ee'/α,  g := U/α + 5ϑκ̃κ̃',  U := (x₁ + ϑc)q̃ + (x₂ + ϑb̃)t̃₁,  c₀ = −3κ̃'² − 4κ̃κ̃''.

`n3qsym.py` checks this split symbolically. The first five terms have the structure of 𝒩, with Ã → (1−ϑ)Ã + ϑ and C̃ → (1−ϑ)C̃ + ϑκ̃²e'' − ϑc₀e + ϑP, and no coefficient grows with x; they go through the weak form below. The last two terms carry the velocity U ≈ (x₁ + ϑc)x, which does grow with x; they go through the strong form. The modulation bounds come from T5's (M): |x| ≤ m̄(1 − m̄ē(A+W₁))⁻¹(ℓ̄A + c_φℛ₂ + n_φη₁(A+W₁+B) + strip). The linear part is T7's C_xΛ^{1/2}, and the residual part x_r := m̄c_φℛ₂ is added explicitly per slab (T7's m̄; c_φ is the 2-norm of (c_φ,j)_{j=1..3}, the largest over the T7 certificates, since the largest component understates it ×1.30 (T7); T6's ℛ₂; x_r ≤ 2.6·10⁻¹⁶ on slabs 2–180), in camp and here (review B6). The nonlinear and E_ij parts are O(η₁Ȳ) relative. T7 Lemma 4c(M) bounds them on the box by x_NL ≤ 1.95·10⁻¹⁵. This fits in the slack C_x(Ȳ − Y_cert) ≥ 2.44·10⁻¹² (ratio ≤ 8.0·10⁻⁴; t7/constants.json meta.modulation_nl) that the certified (Y) level leaves for every Ȳ ≥ 10⁻¹². That level is 1.33·10⁻¹³ (p2/ycomp.py, SLAB0 = 2), with Y(τ₀) ≤ 6.21·10⁻¹³ at the switch. So |x| ≤ C_xȲ + x_r and |ϑ| = |x₃|/n ≤ (C_xȲ + x_r)/n, with the per-slab lower bound of n = ‖χṼ‖_ρ (cconst's nV).

*Weak form, level 3.* Recall w_m = ∂^{m+1}e. At level 3, 𝒩 contributes S := 2∫T₃w₃∂⁴𝒩 to F₃'. Two integrations by parts give S = 2∫(T₃w₃)''∂²𝒩, with

 ∂²𝒩 = Ãw₅ + a₄w₄ + a₃w₃ + a♭,

 a₄ := 2Ã' + B̃,  a₃ := Ã'' + 2B̃' + (∂_{e''}A₃)e''' + ∂_{e''}C̃,  a♭ := [B̃'']♭e''' + [C̃'']♭.

Here ' is the total x-derivative along the solution, and [·]♭ drops the chain-rule term that contains e''''. Since Ã depends on (e, e') and A₃, A₂ depend on z, the coefficients a₄, a₃ and those inside a♭ involve only e, …, e''', the boxed jets. So ∂²𝒩 is linear in the free jets w₃, w₄, w₅, and **there are no trilinear classes**. Write T := T₃ and ϕ := T'(6Ã' + 4B̃) − 2T''Ã. One more integration by parts of the w₄w₅- and w₄²-terms gives

 S = 2∫TÃw₅² + ∫Ψ₄₅w₄w₅ + ∫Ψ₃₅w₃w₅ + ∫Ψ₃₄w₃w₄ + ∫Ψ₃₃w₃² + 2∫(Tw₅ + 2T'w₄ + T''w₃)a♭,

 Ψ₄₅ = 2Ta₄,  Ψ₃₅ = 2Ta₃ + 4T''Ã − T'(6Ã' + 4B̃),  Ψ₃₄ = T'(4a₃ − 6Ã'' − 4B̃') − 2T''B̃ + 2T'''Ã,  Ψ₃₃ = 2T''a₃.

There are no boundary terms. At y_a, T, …, T''' vanish (W2). At ℓ̃ the closed-curve argument of (SE) applies, because 𝒩 is even when e is. S involves w_{≤5}, κ̃^{(≤5)} and r^{(≤3)} only, so the identity extends by density as in §1.5.

*Weak form, levels j ≤ 2.* For j = 1, 2 the contribution 2∫T_jw_j∂^{j+1}𝒩 equals 2∫(T_jw_j)''∂^{j−1}𝒩, with ∂^{j−1}𝒩 = Ãw_{j+2} + α₁w_{j+1} + α₀w_j + α♭ and

 j = 2:  α₁ = Ã' + B̃,  α₀ = B̃' + ∂_{e''}C̃,  α♭ = ∂_XC̃ + (∂_eC̃)e' + (∂_{e'}C̃)e'';
 j = 1:  α₁ = B̃,  α₀ = ∫₀¹∂_{e''}C̃(e, e', se'')ds,  α♭ = C̃(e, e', 0).

With ϕ := 4T'α₁ − 2T''Ã − 2T'Ã', the same steps give S_j = 2∫TÃw_{j+2}² + Σ∫Ψ_{ab}w_aw_b + 2∫(Tw_{j+2} + 2T'w_{j+1} + T''w_j)α♭, where

 Ψ_{j+1,j+2} = 2Tα₁,  Ψ_{j,j+2} = 2Tα₀ + 4T''Ã − T'(4α₁ − 2Ã'),  Ψ_{j,j+1} = T'(4α₀ − 4α₁' + 2Ã'') − 2T''α₁ + 4T''Ã' + 2T'''Ã,  Ψ_{jj} = 2T''α₀

(j = 3 is the case above). At level 0, 2∫T₀w₀∂𝒩 = −2∫(T₀w₀)'𝒩, and one more integration by parts of the Ãw₃-terms gives

 S₀ = 2∫TÃw₂² + ∫2T(Ã' − B̃)w₁w₂ + ∫[4T''Ã + 2T'(2Ã' − B̃)]w₀w₂ + ∫[2T'''Ã + 4T''Ã' + 2T'Ã'']w₀w₁ − 2∫(Tw₁ + T'w₀)C̃,

with |C̃| ≤ γ_e|e| + γ₀|w₀| + γ₁|w₁|, γᵢ := ½ΣⱼZⱼ sup|∂ᵢ∂ⱼA₂|. `n3wcheck.py` checks every one of these identities (the difference of the integrands is an exact derivative; Euler-operator test).

*Strong form for g·e' and (ϑc + x₁)e.* At level j, 2∫T_jw_j∂^{j+1}(gw₀) = ∫[(2j+1)g'T − gT']w_j² + 2Σ_{k=2}^{j+1} C(j+1,k)∫Tw_j g^{(k)}w_{j+1−k}. The flux Tgw_j² vanishes at ℓ̃ because g(ℓ̃) = 0 (q̃(ℓ̃) = t̃₁(ℓ̃) = κ̃'(ℓ̃) = 0 on base v2), and T = 0 at y_a. The x-growth of U is absorbed: gτ₁ is bounded because τ₁ = O(1/Y_ℓ). In g^{(j+1)} the top jet e^{(j+1)} = w_j occurs once, with coefficient ∂_eg = Uκ̃/α²; that part costs 2 sup|∂_eg·e'| in σ_j. Hence

 Δσ_j ∋ sup_C[(2j+1)|g'| + |g||τ₁| + 2|∂_eg||e'|] + 2(x̄ + ϑ̄|c|)  (at j = 0 the ∂_eg-term is part of g'),  Δk_jj ∋ sup_R[((2j+1)|g'| + 2|∂_eg||e'|)T + |g||T'|]/W,
 Δc_{j,j+1−k} = 2C(j+1,k) sup_{C∪R}|g^{(k)}♭|Y_ℓ^{k−1}  (k = 2, …, j+1; ♭: the top jet removed),

where x̄ := C_xȲ + x_r and ϑ̄ := x̄/n. The Δc_{jm} enter g_j as ½Δc_{jm}√F_m. The majorants of g^{(k)} use the Frenet jets of q̃, t̃₁ (q̃' = 1 + κ̃p̃, p̃' = −κ̃q̃, t̃₁' = κ̃ñ, ñ' = −κ̃t̃₁) and (1/α)^{(m)}, m ≤ 4, from `n3qforms.py`.

*Boxes.* The bootstrap hypothesis 𝔹′ is the following, for the slab and x ≥ x_a:
- y ≥ 28: |e| ≤ 1.7·E(x) with E := η₁ + √F₀·y^{3/2}/√84 and √F₀ ≤ A₀Ȳ + B₀ (camp's split, per slab), and |e'| ≤ ζ₁, x|e''| ≤ ζ₂, x²|e'''| ≤ ζ₃;
- y < 28 (x ≤ 34.25): |e|, |e'| ≤ η₁ = 10⁻⁵ (T5's tip box, T5 §2.5; T7 Lemma 4b certifies it: c_e1Ȳ = 2.8·10⁻⁸ at Ȳ = 1.7·10⁻¹², and `n3boxes.py` checks this), and x|e''|, x²|e'''| ≤ η_F = 10⁻² (the e'' part is the assembly's Lemma 4a).

The defaults are ζ₁ = 1.5·10⁻³, ζ₂ = 3·10⁻⁴, ζ₃ = η_F, with Ȳ = 1.7·10⁻¹² in E and in the modulation bound. `n3boxes.py` checks all four for the final camp run, at Ȳ = 1.7·10⁻¹² (every smaller Ȳ is then covered, since the bounds are monotone in Ȳ): |e'| ≤ 5.47·10⁻⁴, x|e''| ≤ 1.20·10⁻⁴, x²|e'''| ≤ 9.33·10⁻³, and √F₀ (with N3′) ≤ 1.00·√F₀ (used). At Ȳ = 10⁻¹² the margins are ×2.76, ×3.97, ×1.69 and ×1.70 (review B2: the first version assumed 7·10⁻⁴, 1.65·10⁻⁴ and factor 1.25, so its true margins were ×1.29, ×2.19, ×1.76 and ×1.25). The |e| box is a fixed point: n3far is run with the camp output without N3′, then once more with the camp output of that n3 file. (E bounds |e| because |e| ≤ |e(28)| + H₂₈√F₀, with |e(28)| ≤ η₁ and H₂₈(y)² := ∫₂₈^y T₀⁻¹ ≤ (y³ − 28³)/84.)

*Majorants per x-range.* Cut each slab's far field into the ramp R, [x_b, 30.25], [30.25, 34.25] and the dyadic ranges [2ᵏ·34.25, 2ᵏ⁺¹·34.25] up to ℓ̃. On each range, the base data are arb hulls of κ̃, …, κ̃^{(5)}, r', r'', r''' and p̃, q̃, ñ, t̃₁ over the range and the slab. The box enters through (Z₀, Z₁, Z₂, Z₃) := (1.25E(x_hi), ζ₁, ζ₂/x_lo, ζ₃/x_lo²) for y ≥ 28, and (η₁, η₁, η_F/x_lo, η_F/x_lo²) below. The box is convex and contains 0, so mean-value forms with sups over the box give (i, j = 0, 1, 2 index z):
- Ã⁺, |Ã|, |Ã'| and |Ã''| directly;
- |B̃| ≤ Σᵢ Zᵢ sup|∂ᵢA₃|, |B̃'| ≤ Σᵢ Zᵢ sup|∂_X∂ᵢA₃| + Σᵢ Zᵢ₊₁ sup|∂ᵢA₃|, |∂_{e''}C̃| ≤ Σⱼ Zⱼ sup|∂ⱼ∂_{e''}A₂|; then |a₄| ≤ 2|Ã'| + |B̃| and |a₃| ≤ |Ã''| + 2|B̃'| + Z₃ sup|∂_{e''}A₃| + |∂_{e''}C̃|;
- |a♭| ≤ β♭_e|e| + Σ_{m≤2} β♭_m|w_m|. In [C̃'']♭ = ∂_X²C̃ + 2Σᵢ(∂_X∂ᵢC̃)zᵢ' + Σᵢⱼ(∂ᵢ∂ⱼA₂)zᵢ'zⱼ' + Σ_{i≤1}(∂ᵢC̃)zᵢ'' and in [B̃'']♭e''', every monomial is one jet times boxed factors. Keep the jet and box the rest. The α♭ of the levels j ≤ 2 and the γ's of level 0 are treated the same way.

For 𝒬, all majorants of 𝒩 are multiplied by 1 + ϑ̄, Ã⁺ and |Ã| gain ϑ̄(1 + |Ã|) and ϑ̄, |∂_{e''}C̃| gains ϑ̄κ̃², the Hessian sups of A₂ gain ϑ̄ times those of P, and the linear terms ϑκ̃²e'' and −ϑc₀e add their own derivatives to the α₀ and β♭ of each level. ∂_X differentiates through the base only (κ̃^{(k)} → κ̃^{(k+1)}, r^{(k)} → r^{(k+1)}). `n3forms.py` (generated by `n3sym.py`) and `n3qforms.py` (generated by `n3qsym.py`) give every derivative needed as an arb function.

*Constants.* On the core (y ≥ y_b), with τ_m := T^{(m)}/T, set

 m₃₅ := 2|a₃| + 4|τ₂||Ã| + |τ₁|(6|Ã'| + 4|B̃|),  m₃₄ := |τ₁|(2|Ã''| + 4|B̃'| + 4Z₃ sup|∂_{e''}A₃| + 4|∂_{e''}C̃|) + 2|τ₂||B̃| + 2|τ₃||Ã|,

 m₄₅ := 2|a₄|,  m₃₃ := 2|τ₂||a₃|,

so that |Ψ_{ab}| ≤ m_{ab}T (sups over the core pieces; m₄₅ also over R₃). On the ramp set M₃₅ := sup_R Ψ₃₅²/(TW), and M₃₄ likewise, using the same bounds with T^{(m)} in place of Tτ_m (T^{(m)2}/T is bounded by (W2)). Also set M₃₃ := sup_R |Ψ₃₃|/W, M_f := sup_R T''²/(TW), M_T := sup_R |T''|/W and τ̄₂ := sup_C |τ₂|. Young's inequality is used with the weights δ_D = 10, δ_R = 10⁻², δ_A = 0.05, δ_B = 3·10⁻³, δ_C = 1, δ_f = 10⁻³ and δ₃₄ := (4M₃₄/M_T)^{1/2}, together with N₄² ≤ D₃ + ¼F₃ + ½∫(T'')⁺w₃² (S3). With c_N := ½m₃₄ + m₄₅/(2δ_D) + ½δ₃₄ + δ_B this gives

 Δθ₃ = 2 sup Ã⁺ + ½(m₃₅ + m₄₅)δ_D + δ_R + δ_A + c_N,

 Δσ₃ = m₃₅/(2δ_D) + ½m₃₄ + m₃₃ + δ_C τ̄₂² + c_N(¼ + ½τ̄₂),

 Δk₃₃ = M₃₅/(4δ_R) + M₃₄/(2δ₃₄) + M₃₃ + δ_f M_f + ½c_N M_T.

The a♭ terms are paired as follows, with ‖a‖²_T := ∫Ta² and ‖a‖²_{(1)} := ∫(T'²/T)a² (the script doubles the 1/δ_A and 4/δ_B terms, which is spare):

 2∫Tw₅a♭ ≤ δ_A D₃ + ‖a♭‖²_T/δ_A,  4∫T'w₄a♭ ≤ δ_B N₄² + 4‖a♭‖²_{(1)}/δ_B,

 2∫T''w₃a♭ ≤ δ_C τ̄₂²F₃ + δ_f M_f V₃² + ‖a♭‖²_{T,C}/δ_C + ‖a♭‖²_{T,R}/δ_f.

Three facts are used: (Σ of four terms)² ≤ 4Σ(·)²; T₃/T_m = Y_ℓ^{6−2m} on the core and y^{6−2m} on the ramp (one ramp for all levels); and on y ≥ 28, e² ≤ 2η₁² + 2H₂₈²F₀. Together they add to h₃ the terms Σ_{m≤2}(Δd₃ₘF_m + Δk₃ₘV_m²) + hf₃, where

 Δd₃ₘ = 4[f_C sup_C β♭_m²Y_ℓ^{6−2m} + f₁ sup_C τ₁²β♭_m²Y_ℓ^{6−2m} + f_R sup_R β♭_m²y^{6−2m}],  plus 8[f_C∫_C Tβ♭_e²H₂₈² + f₁∫_C Tτ₁²β♭_e²H₂₈²] for m = 0,

 Δk₃ₘ = 4f₁ sup_R β♭_m²T'²/(TW),

 hf₃ = 8η₁²[f_C∫_C Tβ♭_e² + f₁∫_C Tτ₁²β♭_e²] + 4η₁²[f_R∫_R Tβ♭_e² + f₁∫_R(T'²/T)β♭_e²],

 f_C = 2/δ_A + 1/δ_C,  f_R = 2/δ_A + 1/δ_f,  f₁ = 8/δ_B.

Levels j = 1, 2 use the same formulas with (j, j+1, j+2) in place of (3, 4, 5), T_j/T_m = Y_ℓ^{2(j−m)}, j + 1 terms in α♭ instead of 4, and the weights δ_D = 1, δ_R = δ_A = 10⁻³, δ_B = 10⁻⁴; θ₂ enters item 4 through C5, so the D-side weights are kept small there. Level 0 pairs the three Ψ's as above (with N₁, √F₀, √D₀ and the ramp datum V₀), and the C̃-term as 2∫T|w₁||C̃| ≤ 2γ₁N₁² + 2γ₀√F₀N₁ + 2N₁(E_abs + E_F√F₀), and so on for T'w₀C̃; here E_abs² := 2η₁²∫_C Tγ_e² + η₁²∫_R Tγ_e² and E_F² := 2∫_{y≥28}Tγ_e²H₂₈². Every product of two norms is split by Young (weights 10⁻⁴ on N₁² for the absolute and ramp T'-terms), and N₁² ≤ D₀ + (¼ + ½τ̄₂)F₀ + ½M_T V₀² (S3).

**Lemma N3′.** Assume 𝔹′, (SE), |x| ≤ x̄ and |ϑ| ≤ x̄/n. Then (F) at each j ≤ 3 holds for the equation with 𝒩 + 𝒬 included, provided θ_j, λ_j, k_jj, k_jm and c_jm are replaced by θ_j − Δθ_j, λ_j − Δσ_j, k_jj + Δk_jj, k_jm + Δk_jm and c_jm + Δc_jm, and h_j gains Σ_{m<j}Δd_jmF_m + hf_j. All of these are per-slab arb bounds.

*Evaluation (base v2, slabs 2–180; certified).* The inputs are the SE-corrected T6 data that `p2/se_corr.py` writes (p2/t6/jets and residuals.json), with its outward rounding. Run `cconst.py` (TOPLOSS = 0.0196) and `merge.py` to get constants_fin.json, and `camp.py` without N3′ to get amp_fin.json. Then run `n3far.py --slabs 2-180 --amp constants/amp_fin.json --out constants/n3_it0.json` (about 30 s, one process) and `camp.py ... --n3 constants/n3_it0.json`. Then run n3far once more with `--amp` set to that camp output (→ n3.json), camp again, and `n3boxes.py`. Every output JSON has an `args` block with the argv, the env and the input paths. The cconst slab files also record sha256 digests of their jets, of residuals.json and of the T7 inputs. A cached slab is reused only if these match, and merge.py checks that they agree. camp enters each Δd_jmF_m through √F_m ≤ A_mȲ + B_m, with an adaptive split of the cross term.
- sup Ã⁺ = 0.0161 at slab 180 on the 1.7-box (sup|κ̃e| ≤ 4.1·10⁻³; the certified value is 2.4·10⁻³). The top loss is 2 sup Ã⁺ ≈ 8 sup|κ̃e|. The assembly's TOPLOSS = 2 sup|A₄ − A₄(0)| ≤ 0.0196, on the unscaled box, is the same quantity.
- Level 3: Δθ₃ ≤ 0.096, which takes θ₃ from 0.220 to ≥ 0.124; Δσ₃ ≤ 2.9·10⁻³; Δk₃₃ ≤ 1.7·10⁹ (+0.15% of k₃₃); Δd₃ₘ ≤ (1.6·10⁵, 74, 0.044); Δk₃ₘ ≤ (9.5·10³, 8.0·10⁵, 4.0·10⁷); Δc₃ₘ ≤ (0.044, 2.8·10⁻³, 9.3·10⁻⁵); hf₃ ≤ 9.1·10⁻¹².
- Levels 0–2: Δθ_j ≤ 0.035 (0.032 of it is the top term); Δσ_j ≤ 2.0·10⁻³, 1.3·10⁻³, 1.3·10⁻³; Δk_jj ≤ 7.3, 8.3·10³, 4.6·10⁶ (≤ 0.05% of k_jj); Δd₁₀ ≤ 1.4, Δd₂ₘ ≤ (2.5·10³, 0.18); Δc_jm ≤ 7.1·10⁻⁴; hf_j ≤ 3·10⁻¹³.
- 𝒬: x̄ ≤ 1.33·10⁻¹¹ and ϑ̄ ≤ 5.8·10⁻⁸ (at Ȳ = 1.7·10⁻¹²; slab 180, n ≥ 2.3·10⁻⁴). Its strong-form share of Δσ_j is ≤ 1.7·10⁻⁶; it is the only source of the Δc_jm.
- Requirement (pipeline of STEP4_RESULT, current KCH, review B4 and B6 applied). On W′ (SUBW): Ȳ ≤ 4.99·10⁻¹² (item 4), against 5.11·10⁻¹² without N3′. On a one-slab window (`--window 180-180`): Ȳ ≤ 1.79·10⁻¹¹ (item 3, E3), against 1.80·10⁻¹¹ without N3′.
- 𝔹′ closes, as checked above (thinnest margin ×1.69, x²|e'''|). N3′ is certified for Ȳ ≤ 1.7·10⁻¹² (the run's Ȳ); the x²|e'''| box itself closes up to Ȳ ≈ 1.83·10⁻¹².
- *Review B1: the residual of the corrected base, computed directly.* `p2/t6c/slabc.py` runs T6's engine on the (SE)-corrected base itself, so se_corr's majorants are no longer used. a₁ = −δ is exact. a₂ is the root of Im Ŵ(ℓ) = 0; on every slab, by Banach, it is unique in a ball of radius ≤ 7.6·10⁻²⁶ (contraction ≤ 2·10⁻¹⁴). a₁′ and a₂′ come from T6's θ_τ = ρ_θ + Φ − Φ(0) and from its velocity identity. The geometry uses e^{iθ} from θ itself, at time degree 20. The two window elements are cut at u = 0.6 and 0.8. On slab 11 a knot crosses u = 0.6, so that slab is also cut in time. The corrected far end encloses 0: |Im Ŵ(ℓ)| ≤ 2.2·10⁻²⁰ and |θ̂(ℓ) + 3π/2| ≤ 1.9·10⁻²⁵. With the bump set to zero, the run reproduces T6. Maxima over slabs 2–180, direct against se_corr: R1 4.6·10⁻¹⁵ vs 2.6·10⁻¹⁴; 𝓡₀ 3.5·10⁻⁹ vs 8.7·10⁻⁹; 𝓡₁ 4.4·10⁻¹¹ vs 4.5·10⁻¹⁰; 𝓡₂ 2.5·10⁻¹⁰ vs 6.9·10⁻⁹; 𝓡₃ 4.9·10⁻⁹ vs 1.7·10⁻⁷. se_corr left out that S₁(ℓ) moves; this shifts ∂⁻¹r over the whole far field. The direct level-0 hats are up to ×3 per slab, but their maxima are below se_corr's. Chain on these data (`constants/*_b1*`): B_j = 7.6·10⁻⁷, 1.8·10⁻⁹, 5.1·10⁻⁸, 7.8·10⁻⁶. The far-box margins at Ȳ = 10⁻¹² are ×6.43, ×4.69, ×1.99 and ×1.70. The requirement is Ȳ ≤ 5.15·10⁻¹² on W′ and ≤ 1.79·10⁻¹¹ on one slab. H4 on one slab is ≤ 0.23·M². ϱ_c is still 0.04711. The (Y) level is 3.2·10⁻¹⁴.

*Remarks.*
- The class terms (A), (B) of the old N3 in `cconst.py` (κ_N, `NLFAC`) are superseded; keeping them is conservative (θ₃ −5·10⁻⁵ and k₃₃ +2.7·10⁵, approximately).
- With `--n3`, cconst's TOPLOSS counts the top term a second time; it may be set to 0.
- `n3far.py --noQ --levels 3` reproduces the level-3-only evaluation of the first version of N3′ (with that version's boxes).

*Conclusion.* N3′ uses F_j, D_j, N_{j+1} (through S3), the lower F_m and the ramp data V₀–V₃. **Level 4 is not used.** Together with C5, level 4 is dropped altogether.

**6.2 Variants.**
- Nested ramps (R_j ⊂ C_{j−2}) replace W-data at level j by D_{j−2}. This is good for level 4 (A₄ falls 40×) and bad for level 3.
- Narrow or steep ramps lose through (S1): the P²/T term at the steep foot of S dominates the transport flux.
- With δ_P = 1.5 (most of D_j used) the best level-3 ramp is [14, 22].

**6.3 Localization.** A transport-following outer factor Ω(Y_ℓ/Y_h), with Y_h' = c_hY_h and c_h < c, removes the bend from the hand-off energies. However:
- its 4th-order commutators at small Y_h cost more than the bend contributes;
- in the global T5 design the bend contributes ≤ 3% of √F₃, through modulation (γ_j on the bend), and < 10⁻⁵ through c_ji.

It is therefore not recommended.

**6.4 Not covered.** These remain outside the lemma:
- the T7 constants: Lemma J inputs, C_x and the modulation constants are taken from t7/constants.json, some still pending or provisional there; μ_e, ϱ_c and K_e are placeholders;
- T9's chart factors and the base's own window energies (T1);
- hypothesis (SE) at the far end (§3.1, Remark), and T4's F_j(τ₀);
- the bootstrap boxes of 6.1 at Ȳ above 1.7·10⁻¹² (they close up to Ȳ ≈ 1.83·10⁻¹²), and the tip box η₁ of 6.1 (review A3).
