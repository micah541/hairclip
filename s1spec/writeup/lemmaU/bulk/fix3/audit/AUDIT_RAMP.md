# G0 ramp audit: the E_I inequality as the verifier codes it (levels 0 and 2)

Fix round 3, unit G0. Nov 16, 2026. Opus. Read-only with respect to the proof tree: nothing was edited or rerun. The only new files are in this directory. The only computation is `g0_ramp_cert.py` (arb, 200 bits, nice 19, 85 MB, 1.3 s; a second run is byte-identical apart from the timing field).

Labels:
- **[cert]**: my arb computation in `g0_ramp_cert.json`, or an exact decimal comparison;
- **[read]**: verified by reading the cited file:line;
- **[hand]**: my derivation by hand, not machine-checked;
- **[est]**: estimate.

Paths: `V/` = `/home/claude/s1spec/step3nl_v33/verifier/`, `S/` = `/home/claude/s1spec/step3nl_v33/`.

## Verdict: KILL. A verifier rerun is needed, and I stopped there.

The kill condition of §8 holds in its strongest form.
- **The rounding reserve of the inner-ramp couplings is negative.** The verifier types four certified upper bounds rounded *down* ([cert], §4):
  - CF0 (c3): relative −9.65e-5;
  - CF2 (c1): −6.53e-5;
  - CF3 (c0): −4.05e-5;
  - AF (a_f): −5.4e-7.
  So the inner-ramp linear commutator itself is not covered as coded, and no uncovered inner-ramp term can be absorbed by rounding.
- **There are uncovered inner-ramp terms**, so the kill check triggers even without the rounding defect:
  1. **Level 2, quasilinear commutator** −d001·[∂⁴, ζ]g on [10.25, 16.25]. Coefficients relative to CF are (7.74, 7.74, 7.74, 7.73)e-5 [cert]. Together with the truncation, the required increase of c3 is ≥ 1.74e-4 relative.
  2. **Level 2, other nonlinear terms on the inner ramp.** The code bounds ζ𝒩′ there with v-norms V_j and with the x ≥ 16.25 sup `supv1`. Neither is valid where ζ < 1. Uncovered are the lower-order commutator pieces and, for example, d110·w″·w⁗, which is ≈ 0.66·0.297·sup|w″|_ramp·D4·‖D2v‖ [est]. At sup|w″| ~ 1e-2 that is ~0.4% of the CF0 coupling [est].
  3. **Level 0, nonlinear monomials on the inner ramp.** The quadratic and cubic monomials of 𝒩 (sups d200 0.0133, d110 0.297, d300 3.02 on [10.25, 16.25], box 1e-2 [read]) are not in term (a), which carries the linear monomials only.
  4. **Level 2, modulation commutator** −(a q + b cosθ)ζ′g on the inner ramp. It is ≤ 1.271·|(a,b)|·W1·‖D2v‖ [cert], i.e. 6.6e-7 × CF3 at ā = 2.19103e-7. It is the smallest of the four. **ROUTE §4's suspicion is confirmed, and the term is about 3× smaller than estimated there.**
- **There are verifier-side gaps that are not ramp terms**, found while writing out the inequality on [34.25, 60]. Two of them are not small [hand/est]:
  - **(N1) Top-order cross terms of 1 − J⁻⁴ are missing.**
    - At level 0 the code charges the w⁗ term of 𝒩 as d001·V2² (`rhs_iv.py:147`). The integration by parts behind that form leaves (sup|c′ − ψ′c| + Ψ1·sup|c|)·V1V2 + O(ψ′²c)·V0V2. With c = d001 ≤ 2.27e-2 that is ≥ 0.059·V1V2 + ~0.04·V0V2.
    - The code has 1.67e-3·V1V2. ASSEMBLY.md:78 still lists "0.025 V1V2", but it was lost in `inner_closure_float2.py:92` and in `rhs_iv.py`.
    - At level 2 the code has d001·V4² (`rhs_iv.py:178`) where the pairing gives d001·V4·‖D2v‖. The missing part is d001·V4·(2.61V3 + 1.7V2) ≈ 0.059·V3V4 + 0.039·V2V4.
    - Worst case, about 10–20% of the E_I dissipation at the box coefficient d001 = 2.27e-2 [est]. The tail rate is ρ = 2.69e-3, so a rerun is needed to know whether it still closes.
  - **(N2) The β-, ψ1- and ψ2-couplings of the v-equation are absent from I0 and I2.** These are |β|·‖ζψ̂3′‖_W = 0.26|β| and the analogous terms. ASSEMBLY.md:80 lists them "[C norms]"; `rhs_iv.py` has no β at all. Size up to ≈ 4.6e-4 in rate at α = 1.5e-4, and ≈ 9e-5 at α → 0, against ρ = 2.69e-3 [est].
  - Smaller omissions are listed in §2.3.

**Consequence.** The verifier's E_I inequality (tag fin12f4) is not an upper bound for ½·d/dτ(V0² + ε_I·V2²) as written. Fixing it needs a written derivation and a verifier rerun with a new tag, followed by re-pinning downstream. Per ROUTE §11 Q3 this is Micah's decision. Nothing on that path was started.

**The G2 flux-lemma list is still well defined** (§5), because the outer-ramp terms are the same whether or not the verifier is rerun.

## 1. The inequality as coded

Variables: v = ζw′ = ζg, V_j = ‖v^{(j)}‖_W; ζ = 1 − S7((16.25 − x)/6) on [10.25, 16.25], 1 on [16.25, 57], S7((60 − x)/3) on [57, 60]; ρ-norms A, W1, B, W3, D4 of w̃. Shorthand: Vq = V1 + 1.67·V0 (`rhs_iv.py:142`); D2v = V4 + 2.61·V3 + 1.7·V2 (`:176`). The (a,b) bound is CAB·Am = |(a,b)| (Euclidean), with Am = Amx + Amc (`:143-145`, `constants.py:73`).

**Level 0** (multiplied by μ). The positive part I0p is at `rhs_iv.py:147-149`, the additions at `:150-174`, the negative part at `:189`:
```
-(1-θ)(μ1 V0² + ½V2²)
+ AF W1² + BF B²                                              [inner-ramp linear commutator]
+ d001 V2² + (d010 V2 + d100 V1 + A1 V0)·Vq                   [𝒩, x ≥ 16.25, TQ constants.py:28]
+ (d200 s V1 + d110 s V2 + d300 s² V1)·Vq,  s = supv1 = max(sI2, S1V/32)   (:100-102, :141)
+ |(a,b)|·0.2889·|α|·V0 + |(a,b)|·(38.5 + 0.66)·V0² + |(a,b)|·2.5·W1²     (CABQ, CAB3985, CAB25: :45-46)
+ FnW(α)·Vq                                                    [forcing, L group, :150]
+ 0.66·Σ_j (c0_j + dw_j S_ramp)·(D4, W3, B, W1, A)_j·Vq        [term (a): 𝒩 on the inner ramp, linear monomials, :152-161]
+ 9.9e-12·S_all·V0                                             [c0′w, :164-168]
+ FLUX0                                                        [outer-ramp flux, L group, :170-172]
```

**Level 2** (multiplied by μ·ε_I), at `rhs_iv.py:176-187` and `:189`:
```
-(1-θ)(½V4² + μ2 V2²)
+ (CF0 D4 + CF1 W3 + CF2 B + CF3 W1)·D2v                      [inner-ramp linear commutator]
+ d001 V4² + (d010 V3 + d110 s V3 + 2 d200 s V3 + d110 s (3V3 + 1.31V2) + d300 s² V2)·D2v
+ |(a,b)|·39.16·V2² + 8.19·|(a,b)|·V2V3 + |(a,b)|·0.3·|α|·D2v   (62·Am = (62/7.57)|(a,b)|)
+ |α|³·D2v                                                     [forcing, L group, :181]
+ FLUX2                                                        [outer-ramp flux, :183-185]
```
No other terms exist: `grep` finds no β, no c0′, no ψ1, ψ2 and no Gårding e-terms in `rhs_iv.py`.

## 2. What the inequality must contain, term by term

The v-equation (ASSEMBLY.md:70, BOOTSTRAP.md:11):
v_τ = L1v − [L1,ζ]g + ζc0′w + ζΦ_F′ + ζ𝒩′ + ζM′, with M = −βψ̂3 − 2αβĥ2 − aψ1 − bψ2 + aA1(u) + bB1(u).
- A1(u) = −(u − u′q/α_g) and B1(u) = u′cosθ/α_g (`agent_int/v2/int_v56.py:4-6`).
- Their w-parts in the g-equation are a(q g′ + (q′ − 1)g) and b(cosθ g′ − k sinθ g), up to O(k) corrections of 1/α_g.
- Pairings: level 0 uses ⟨v, ·⟩_W. Level 2 uses ⟨v″, (·)″⟩_W = ⟨D2v, ·⟩_W, with D2v = (Wv″)″/W.

Columns: **IR** = inner ramp [10.25, 16.25]; **Z1** = ζ ≡ 1 on [16.25, 57], which contains [34.25, 57]; **OR** = outer ramp [57, 60]. "cov" = covered by the cited constant; "G2" = belongs to the flux lemma; "✗" = uncovered by the verifier.

### 2.1 Level 0

| # | Term; where ζ is commuted or integrated by parts | IR | Z1 | OR |
|---|---|---|---|---|
| 0.1 | ⟨v, L1v⟩ ≤ −μ1V0² − ½V2²: no commutation (v ∈ H⁴, compact support) | cov MU1 = 0.1552 (`constants.py:19`; lint_cert_p, XA = 10.25) | cov | cov |
| 0.2 | −⟨v, [L1,ζ]g⟩. Inner: identity = ∫A_z g² + B_z g′² (`agent_coup/verify_identity.py`) | cov in form AF W1² + BF B² (`constants.py:23`). **AF typed 6.9555e-2 < raw 6.955503763e-2 ✗ (−5.4e-7 rel)**; BF ok | — | G2 (F1: K0 ≤ 3.9038e-21, done) |
| 0.3 | ⟨v, ζc0′w⟩ ≤ ‖ζc0′‖_W·sup\|w\|·V0, direct | cov C0P_NORM 9.9e-12 ≥ 9.875e-12 (v2 V5) | cov | cov |
| 0.4 | ζΦ_F′, one IBP: −∫(Wζv)′Φ_F ≤ ‖Φ_F‖_W(V1 + κV0), κ = sup(\|ζ′\| + ψ′ζ) | κ = 0.5340 [cert] | 1.2827 [cert] | 1.5883 [cert] ≤ 1.67: cov (FNW `constants.py:31`; C⁴ values 1.243e-8α + 4.03e-8α² are below FNW for all α ≤ 2.8 [cert arithmetic]) |
| 0.5a | ζ𝒩′, one IBP: −∫(Wζv)′𝒩, linear monomials | cov term (a): factor 0.66 ≥ sup(W/ρ)^½ = 0.65580 [cert]; κ_in ≤ 1.67 | cov: d010 V2, d100 V1, A1 V0 (·Vq) | **G2** (𝒩 is in g, g′, g″, g‴ on [57, 60], not in V_j) |
| 0.5b | same, quadratic and cubic monomials d200 w″², d110 w″w‴, d300 w″³ | **✗ not in term (a)** (`pb_ramp.py` KEYS are linear only) | cov (s·V1, s·V2, s²·V1) | G2 |
| 0.5c | same, top order d001 w⁗ = d001 v‴ | in term (a) (linear) | **✗ (N1)**: coded d001V2²; needs + (sup\|c′−ψ′c\| + Ψ1 sup\|c\|)V1V2 + sup\|(Wψ′c)′/W\|·V0V2 ≈ 0.059V1V2 + 0.04V0V2 [hand/est] | G2 |
| 0.5d | A0·w (w itself) | in term (a) | ✗ minor: A0 ≤ 9.2e-9 on [16.25, 60] (`r60_T3`), not coded | G2 |
| 0.6a | a-transport, main part a(q v′ + (q′−1)v): IBP, coefficient ½ψ′q + ½q′ − 1 | sup 37.179, inf 1.046 [cert] | same | same: cov by 38.5 (`CAB3985`) |
| 0.6b | b-transport, main part b(cosθ v′ − k sinθ v): coefficient ≤ ½ψ′ + 2\|k\| | 0.6528 [cert]; with 0.6a, Euclidean 37.185 ≤ 39.16 cov | cov | cov |
| 0.6c | transport commutators −aqζ′g − b cosθ ζ′g, paired with ζg | ≤ \|a\|·0.4336 + \|b\|·0.0726 ≤ 2.5\|(a,b)\| [cert: 0.0726 ≥ sup ζ′(W/ρ)^½ · sup(W/ρ)^½]: cov CAB25 | — | **G2** (−a∫Wqζζ′g², −b∫W cosθ ζζ′g²; sign of a unknown; ROUTE §4) |
| 0.6d | aA1(u0)′ + bB1(u0)′, direct norm with ζ | cov: \|(a,b)\|·√(0.2731² + 0.01589²) = 0.2736 ≤ 0.2889 (CABQ) | cov | cov |
| 0.6e | −βψ̂3′ − 2αβĥ2′ − aψ1′ − bψ2′, direct (norms 0.2596, 1.796, 1.709e-7, 1.683e-8, v2 V5) | **✗ (N2)** not coded anywhere | ✗ | ✗ (inside the same norm) |
| 0.6f | O(k) parts of A1, B1 through 1/α_g | not audited (k ≤ 1.7e-4 at σ = 8) | — | — |

### 2.2 Level 2

‖D2v‖_W ≤ V4 + 2Ψ1V3 + sup|ψ′² − ψ″|·V2. With 2Ψ1 ≤ 2.6097 and sup|ψ′² − ψ″| ≤ 1.69544 on [10.25, 60] [cert], `D2v` (`:176`) is covered.

| # | Term | IR | Z1 | OR |
|---|---|---|---|---|
| 2.1 | ⟨v″, (L1v)″⟩ ≤ −½V4² − μ2V2² + e-terms (\|e1\| ≤ 5.66e-6, \|e0\| ≤ 5.11e-6) | e-terms ✗ minor (ASSEMBLY:86 lists 5.7e-6(V3V2 + V2²); not coded) | same | same |
| 2.2 | −⟨D2v, [L1,ζ]g⟩ | form cov: c_j = sup\|coef_j\|(W/ρ)^½ (`run_coup.py:71-72`). **CF0, CF2, CF3 typed below raw ✗** (§4) | — | **G2** (F1: P2 pointwise; T3, O4 need the L¹ channel) |
| 2.3 | ⟨D2v, ζc0′w⟩ ≤ 9.9e-12·S_all·‖D2v‖ | ✗ minor, not coded | ✗ minor | ✗ minor |
| 2.4 | ⟨D2v, ζΦ_F′⟩ direct: ‖ζΦ_F′‖_W ≤ 2.70e-8\|α\| + 1.156e-7α² + 2.603e-3\|α\|³ (v3 F1a/F1b) | ✗: coded 1·\|α\|³ only (`:181`), not an upper bound for \|α\| ≲ 1.6e-4, i.e. on the whole run (STEP3_REVIEW.md:63, never fixed). Absolute size tiny | same | same |
| 2.5a | top order: ζ·d001·w⁽⁵⁾ = d001(v⁗ − [∂⁴,ζ]g) | **✗** commutator: d001·(4ζ′, 6ζ″, 4ζ‴, ζ⁗)(W/ρ)^½ ≤ (3.42, 4.07, 3.28, 3.28)e-5 on (D4, W3, B, W1) [cert], = 7.74e-5 × CF | **✗ (N1)**: coded d001V4²; needs d001·V4·‖D2v‖, missing ≈ 0.059V3V4 + 0.039V2V4 [cert coefficients × d001] | **G2** (the commutator part, with d001 at the core box) |
| 2.5b | Dd001·w⁗ + d010 w⁗, Dd010 w‴ + d100 w‴, … (Dd = total x-derivative, including ∂_{U1}d·(u0″ + w″)) | ✗ (ζw⁗ ≠ v‴ on the ramp; needs ρ-norms ×0.66) | partly ✗: d010 V3 coded; Dd001 V3 (≈ 4U1(u0″ + s)·V3 [est]), d100 V2 and the A0, A1 terms missing | G2 |
| 2.5c | quadratic and cubic: (d200 w″²)′, (d110 w″w‴)′, (d300 w″³)′ | ✗ (`supv1` is a sup over x ≥ 16.25 and V_j ≠ ρ-norms there; e.g. d110 w″w⁗ ≈ 0.66·0.297·s_ramp·D4 [est]) | partly ✗: 2d200 s w‴ coded with V3 instead of V2 (costs ×15.35 by Hardy); 3d300 s²V2 coded with coefficient 1; d110 terms cov (GN identity) | G2 |
| 2.6a | a-transport main part, ⟨v″, (qv′ + (q′−1)v)″⟩_W: coefficient of ∫Wv″² is ½ψ′q + (5/2)q′ − 1; q″, q‴ terms tiny | sup 39.1788 [cert]; with b (0.6529) Euclidean 39.1842 **> 39.16** | same | same |
| — | ↳ the deficit 0.0242\|(a,b)\|V2² is absorbed by the coded 8.19\|(a,b)\|V2V3 via Hardy–Poincaré V2 ≤ (2/inf ψ′)V3 = 15.35V3 (v″(10.25) = 0; inf ψ′ = 0.13031 [cert]), which needs 0.371 ≤ 8.19: cov, but the lemma is unwritten | | | |
| 2.6b | transport commutators −(aq + b cosθ)ζ′g paired with D2v | **✗** ≤ 1.2709\|(a,b)\|·W1·‖D2v‖ [cert] (= 6.55e-7 × CF3 at ā) | — | **G2** |
| 2.6c | aA1(u0)′ + bB1(u0)′ direct | cov 0.2736 ≤ 0.3 (CAB03) | cov | cov |
| 2.6d | β, ψ1, ψ2 terms direct | ✗ (N2) | ✗ | ✗ |

### 2.3 Summary of verifier-side items, by size
1. **Probably material; a rerun is needed to know.**
   - (N1) top-order cross terms, levels 0 and 2;
   - (N2) β-coupling in I0 and I2;
   - level-0 and level-2 nonlinear quadratic and cubic monomials on the inner ramp (0.5b, 2.5c), at ~100× the coded term-(a) W3 coefficient [est, s_ramp ~ 1e-2].
2. **Small but uncovered.**
   - CF/AF truncation;
   - the quasilinear commutator 2.5a on the inner ramp;
   - the modulation commutator 2.6b on the inner ramp;
   - Dd001 terms;
   - the 2d200 and 3d300 index and coefficient errors.
3. **Tiny.** c0′w at level 2; FnW′ α, α² parts; A0w; ψ1, ψ2; Gårding e-terms.
4. **Strip constants** (not a ramp item, recorded for the rerun): SB, SA, IA, DS1, DS2 are also typed below their raw values (−4.6e-6 … −2.4e-5 relative, §4).

## 3. The form behind the constants [read + cert]

| Constant | Where | Form | Certified need | Status |
|---|---|---|---|---|
| 38.5 | CAB3985, levels 0 and 2 | ½·sup\|1 − qψ′\| with the "F" value 77 (ASSEMBLY:81); a-transport main part | level 0: 37.179; level 2: 39.179 (the extra 3/2 comes from (qv′)″ ∋ 2q′v″ and the −w part of A1) | level 0 cov; level 2 needs the Hardy step |
| 0.66 (in CAB3985) | same | ½·sup\|(W cosθ)′/W\| ≤ ½·sup√(ψ′² + k²); v3 F2 gives 1.306/2 = 0.653 | 0.6528 / 0.6529 | cov |
| 0.66 (term (a), `ADD_I0_RAMP[0]`) | `constants.py:35` | sup over the inner ramp of (W/ρ)^½ | 0.65580 | cov (reserve 0.64%) |
| 2.5 | CAB25 | sup over the inner ramp of qζ′W/ρ (v3 F2: 0.4357) | 0.4336 for a | cov (×5.7) |
| 62 | `c62`, level 2 V2V3 | ASSEMBLY:90 "2\|a\|V2V3" with the old \|(a,b)\| ≤ 31.2A, so 2·31.2 = 62.4 rounded to 62 [inferred]. Now 62·Am = 8.19\|(a,b)\|. My derivation has no V2V3 term from transport | — | surplus; absorbs the 0.0242 deficit (0.371 needed) |
| 0.3 | CAB03 | ≥ ‖ζA1(u0)′‖_W/α + ‖ζB1(u0)′‖_W/α = 0.2731 + 0.01589 (ASSEMBLY labels it ‴; v3 ‖ζA1(u0)‴‖/α = 7.2e-5) | 0.2736 (Euclidean) | cov |
| 0.273 + 0.0159 | CABQ, level 0 | the same norms; 0.273 is typed below 0.2731, but cov by Euclidean \|(a,b)\| | 0.2736 ≤ 0.2889 | cov |
| 1.67 | Vq | sup(\|ζ′\| + ψ′ζ) over supp ζ | 1.5883 (outer ramp) | cov |
| 2.61, 1.7 | D2v | 2Ψ1; sup\|ψ′² − ψ″\| | 2.6097; 1.69544 | cov (thin: 6e-4, 3e-3) |
| CF | level 2 | c_j = sup\|coef_j([L1,ζ])\|(W/ρ)^½; essentially C(4,k)\|ζ^{(k)}\|(W/ρ)^½ (my looser bounds 0.4429, 0.5267, 0.4242, 0.4251) | raw `coupling_constants.json` | **typed below raw** |

## 4. Rounding reserves [cert: exact decimal comparison against `agent_coup/coupling_constants.json` P = 2 rigorous]

| Constant | Raw certified upper bound | Typed (`constants.py`) | Relative reserve |
|---|---|---|---|
| CF0 = c3 | 0.442142686616 | 0.4421 | **−9.65e-5** |
| CF1 = c2 | 0.525899483214 | 0.5259 | +9.8e-7 |
| CF2 = c1 | 0.423527657147 | 0.4235 | **−6.53e-5** |
| CF3 = c0 | 0.425017231881 | 0.4250 | **−4.05e-5** |
| AF = a_f | 0.069555037634 | 6.9555e-2 | **−5.4e-7** |
| BF = b_f | 0.039358915436 | 3.9359e-2 | +2.15e-6 |
| SB, SA, IA | — | — | −4.6e-6, −2.1e-6, −5.7e-7 |
| DS0, DS1, DS2, DS3 (d3, d2, d1, d0) | — | — | +1.9e-5, −2.4e-5, −7.6e-6, +3.2e-6 |

The verifier rounds every typed literal ±1 ulp (`rhs_iv.py:14-16`). That makes the double contain the decimal, not the certified value.

**Kill-check arithmetic.** The required relative increase of the inner-ramp couplings is 9.65e-5 (truncation) + 7.74e-5 (quasilinear commutator) + 6.6e-7 (modulation, W1 only), plus the uncovered nonlinear pieces (2.5b, 2.5c). The available reserve is negative. So the kill condition holds.

The only nearby reserves are in modulation terms (CAB25 ×5.7, c62), and they multiply different monomials (W1², V2V3). They cannot absorb D4·D2v, W3·D2v or B·D2v terms.

## 5. Final term list for the G2 flux lemma (outer ramp [57, 60])

This holds unchanged after the verifier is fixed. In each item, everything is expressed through g = w′ on [57, 60], so through (A)'s z_j, and carries W ≤ W(57) ≤ 1.44e-18.

**Level 0:**
1. the linear commutator −⟨v, [L1,ζ]g⟩_W. Already done in F1 (K0), drift included.
2. the transport commutators −a∫Wqζζ′g² and −b∫W cosθ ζζ′g², with |a|, |b| ≤ 2.19103e-7. Use them as a ball in the commutator operator, as ROUTE §4 says; do not reuse the drift sign.
3. the full nonlinear part ∫_{57}^{60}(Wζv)′𝒩. That is the quasilinear (1 − J⁻⁴)g‴, the lower linear monomials, d200 g′², d110 g′g″, d300 g′³ and A0 w, with the core box of (A) as coefficient source.
   - The verifier's V-norm terms do not cover [57, 60]; there they count something else.
   - So the whole [57, 60] part goes to G2, not a difference.

**Level 2:**
1. the linear commutator −⟨D2v, [L1,ζ]g⟩_W. F1 has P2, T3 and O4; the T3 and O4 parts go through the L¹ channel and the trace route of ROUTE §2.2.
2. the quasilinear commutator −⟨D2v, (1 − J⁻⁴)[∂⁴,ζ]g⟩_W. Its relative size against item 1 is ≤ sup|1 − J⁻⁴| on the core box, ≈ 2.3e-2 at the verifier box. Note that the d001·v⁗ part pairs with ‖D2v‖ (not only V4) on [57, 60] too; that is the same (N1) defect, here localized.
3. the remaining nonlinear terms ⟨D2v, ζ𝒩′⟩_W on [57, 60] other than d001·v⁗: Dd001·g‴, d010·g‴, d100·g″, (d200 g′²)′, (d110 g′g″)′, (d300 g′³)′, and the A0, A1 terms.
   - Each one can be bounded by ‖D2v‖_{W,[57,60]} times a g-norm on the ramp.
   - Pointwise or time-integrated, by order: g‴ goes through the L¹ channel as T3 does.
4. the transport commutators −⟨D2v, (aq + b cosθ)ζ′g⟩_W.

**Not in G2** (covered in the verifier on [57, 60] once its own defects are fixed): forcing (κ = 1.5883 ≤ 1.67 at level 0, direct at level 2), c0′w, A1(u0)′ and B1(u0)′, and the transport main parts (0.6a, 0.6b, 2.6a).

## 6. What a rerun would have to change (input for Micah's decision; nothing done)

1. **Write the E_I inequality out analytically first.** No document derives it: CERT_PLAN says "bookkeeping … to be written out analytically", and ASSEMBLY.md and the code differ in several places. A sympy Euler-operator check of each IBP, as in F1, would make it auditable.
2. **Constants.**
   - Outward-round CF, AF and BF from `coupling_constants.json`, and the strip constants SB, SA, IA, DS.
   - Add the inner-ramp level-2 terms in ρ-norms: the quasilinear commutator, the lower-order and quadratic 𝒩′ terms with ×0.66, and the modulation commutator 1.271|(a,b)|W1.
   - Add the quadratic and cubic monomials to term (a).
3. **Missing terms.**
   - Fix (N1) at levels 0 and 2, including the Dd001 terms, and correct 2d200·V3 → V2 and d300 → 3d300.
   - Add (N2): β, ψ1, ψ2 at both levels, using the closure's β bound.
   - Add c0′w at level 2, the full ‖ζΦ_F′‖, the Gårding e-terms and A0·w.
4. **Optional but probably decisive for margins.** Recertify the far T3 coefficients at the actual α ≤ 1.5e-4. d001 ≈ 2u′² is dominated by the α ≤ 1e-3 box: |u0′| ≤ 0.085 there against ≈ 0.013 at α = 1.5e-4 [est]. This would shrink the (N1) terms by ~10× [est].
5. **Then** rerun at both α points with a new tag. Re-derive mod_sups and the J's; re-pin (A) → U → zone3 and every consumer of the fin12f4 outputs (V33_VERIFIER.md list). The compute is minutes: 215 s per α point in fin12f4.

## 7. Honest limits
- The "needed" side of §2 is a hand derivation [hand]. I did not machine-check the IBP identities. The (N1) coefficients are lower bounds of what any IBP route must pay (the 2Ψ1·c cross terms), not a full certified expression.
- I did not evaluate the verifier's rate with the missing terms. That would be a rerun, which the kill rule forbids. Rate impacts are [est] from Young's inequality against the coded dissipation. The verifier may well close after the fixes; this audit does not decide that.
- sup|w″| on the inner ramp along the trajectory was not extracted. It is not in the stored outputs; s_ramp ~ 1e-2 is [est].
- The O(k) corrections of A1 and B1 (1/α_g) and the strip-side modulation terms were not audited.
- I read the code as it is at the sha256 values below. ASSEMBLY.md and VERIFIER.md are not proofs, and I did not treat them as such.

## Files
- `g0_ramp_cert.py`: arb sups (κ, transport coefficients, Hardy constant, inner-ramp commutator coefficients, D2v constants, rounding reserves). It has an args block with input sha256s.
- `g0_ramp_cert.json`: the output; a rerun is byte-identical apart from `runtime_s`.

Inputs read (sha256 prefixes):
- `V/constants.py` 3f74fdd6008e0a8c
- `V/rhs_iv.py` 4165b774b7083041
- `V/closure.py` d1a1ac116840b73b
- `S/ASSEMBLY.md` 551e72ebe183366f
- `S/agent_coup/coupling_constants.json` b46b0d6cd19bb624
- `S/agent_int/v2/int_constants_v2.json` a47108ec703ae34b
- `S/agent_int/v3/int_constants_v3.json` b9f43492a6faf682
- `S/agent_int/planB/planB_ramp.json` c8fc5bce02b14cf4
- `S/agent_int/planB/r60_T3_boxw1_0.02.json` dcf2b21314d4cd8f
