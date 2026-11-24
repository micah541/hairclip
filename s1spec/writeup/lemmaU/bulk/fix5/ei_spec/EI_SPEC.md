# EI_SPEC: derivation of the E_I differential inequality (spec for verifier v4)

Fix round 5, unit "ei_spec". Nov 19, 2026. Fable (design). Read-only on the proof tree; everything new is in this directory.
Nothing here is closed. This is a specification, not a proof and not an implementation: verifier v4 was **not** started.

Labels:
- **[sym]**: exact sympy identity, checked by `ei_identities.py` -> `ei_identities.json` (all residuals 0);
- **[cert]**: certified number quoted from an existing output (file named);
- **[read]**: verified by reading the cited file:line;
- **[hand]**: my derivation by hand, not machine-checked;
- **[est]**: estimate, not a bound;
- **[probe]**: rigorous interval branch-and-bound evaluation (the fin12f4 `bb.py`) of an inequality whose *new* constants are [est]. So the result is [est].

`GI` = `/home/claude/s1spec/writeup/lemmaU/bulk/fix5/g0_indep/G0_INDEP.md`, the independent check (Fable, Nov 19; verdict FIX-AND-RERUN). It finished while this spec was being written; I read it afterwards and folded it in (§0 item 8, §4.6, §4.9, §6, §7, §8). Items taken from it are marked [GI].

Paths: `V/` = `/home/claude/s1spec/step3nl_v33/verifier/`, `S/` = `/home/claude/s1spec/step3nl_v33/`, `G0` = `/home/claude/s1spec/writeup/lemmaU/bulk/fix3/audit/g0_ramp_cert.json`. Line numbers of `rhs_iv.py`, `constants.py`, `closure.py`, `unc.py` are those of tag fin12f4.

## 0. Summary

1. **Feasibility: the corrected inequality very probably still closes at the current boxes, with no lever pulled** [probe]. With every term of this spec added (new constants [est], far-field coefficients at the present alpha <= 1e-3 box):
   - tail anchor rate: -2.6930e-3 -> **-2.6872e-3** (rho moves by 0.2%);
   - first march cell: +3.03441e-2 -> +3.03515e-2;
   - with all new terms multiplied by 3: -2.628e-3 and +3.0483e-2;
   - full march at alpha_1 = 1.5e-4: closes in 10 steps, peak x1.00002, J_m x1.0005, sup|(a,b)| x1.00001, sup|beta| unchanged (§8);
   - **the sensitive output is J, not rho**: in the x3 run the unconditional pass drops from delta_t = 0.2 to 0.01 and the J constants jump x6 (§8). At x1 delta_t stays 0.2.
   Reason: the maximising shapes are A-dominated (inner energy); E_I has its own dissipation mu_1 = 0.155 against a net rate of 2.7e-3, so E_I-side terms of a few percent of the E_I dissipation do not reach the maximiser. The audit's "10-20% of the E_I dissipation" for (N1) is correct in size and harmless in effect.
2. **Two design constraints matter more than the sizes** (§5.3):
   - every new term must respect the three degree conditions of `tail_cert.py` (radial degree <= 2 in (D4, V4); x-valuation >= 2 in H; joint (alpha, x)-valuation >= 2 in L). Two natural bounds violate them (the inner-ramp monomial w''·w'''' at level 2, and the beta-coupling with the constant strip term). The spec gives forms that pass; the probe passes `tail_cert` with them.
   - `unc.py` must majorise every L term; it does not today (§11, SF-3).
3. **Scheme.** The quasilinear top order is kept globally as (c v''')' with c = 1 - J^-4, so there is no boundary term at x = 16.25 or x = 57. Each ramp then carries the commutator c·[d^4, zeta]g + c'·[d^3, zeta]g. On the inner ramp it is bounded with rho-norms. On the outer ramp it goes to the flux lemma, exactly as ROUTE_DECISION §3 asks.
4. **Outer ramp.** §4.10 defines the two functionals K0^OR, K2^OR that the flux lemma v2 (unit G2) must bound. They are consistent with AUDIT §5 and ROUTE §1; three tiny items are added (c0'w, A0·w and the O(k) modulation remainders on [57, 60]).
5. **Levers** (§9) are not needed for closure [probe]. The cheapest useful one is a recertification of the far coefficients at the alpha-level 2.5e-4, which costs no downstream change.
6. **Coding plan** (§10): new directory `step3nl_v4/`, tag `fin13v4`, about 26-30 agent-hours plus one Fable review. Compute: under 1 h in total, one process, < 800 MB.
7. **Side findings outside the E_I block** (§11): three items in the same verifier that v4 should fix or justify. I did not audit the E/E2 blocks.

8. **Independent check [GI] folded in.** It agrees with this spec on the scheme (global top-order term, ramp commutators) and on every missing term. Four points needed changes or explicit answers here:
   - (a) *Level-2 inner-ramp nonlinear block* (GI 2-10): GI says it cannot be bounded in rho-norms as the code is organised, because it needs w^(5) and the obvious bound has radial degree 2.25. Both obstacles are removed by the reformulation of §4.3, §4.4, §4.6: w^(5) never appears on the ramp (it is inside the global (c v''')', and the remainder is the commutator c K4 + c' K3, of order <= w''''), and the two degree-2.25 monomials are bounded with one factor converted to v through the cutoff (S1z, S2z). So it is a reformulation, as GI says, and this spec is that reformulation. It still needs new certified inputs (the derivative coefficients of stage C2).
   - (b) *beta strip constant* (GI 2.1, 4.4.1): §4.9 uses the state form, which costs nothing in floor or slopes but needs the modulation lemma restated (U-3). The fallback (constant kept) costs F0 x7.4 and J-slopes x8.4 [cert GI] and a change of the tail lemma. Details in §4.9.
   - (c) *More wrong-side constants and one wrong norm* [GI 1.3, 0-14, 0-1]: DS1, DS2, SB, SA, IA, d010, d110, d300, FN[0], FN[1], C3, FNW[0]; `rhs_iv.py:144` uses FnW where the modulation lemma has ||chi Phi_F||_rho = Fn; MU1 = 0.1552 is above the L-int certificate it cites and valid only through `v2_T4_P2_0.25.json`. All are in §6 and §7 now. My first draft had copied the FnW error into ABc; corrected.
   - (d) *Feasibility evidence*: GI's own trial rerun and my probe are independent implementations and agree (§8).

## 1. Setting and notation

### 1.1 Objects [read: S/ASSEMBLY.md:8-24, S/BOOTSTRAP.md:7-14, V/constants.py]
- x = sigma - sigma_tip on the half line (even class). u = alpha·psi3 + alpha^2·h2 + w = u0 + w. (psi3, h2 are the C^4 versions of `agent_int/v3/c4fun.py`; see U-6 in §11.)
- g = w'. zeta = 0 on x <= 10.25, 1 - S7((16.25 - x)/6) on **IR** = [10.25, 16.25], 1 on **Z1** = [16.25, 57], S7((60 - x)/3) on **OR** = [57, 60], 0 beyond. S7(t) = t^4(35 - 84t + 70t^2 - 20t^3), so zeta is C^3 with zeta^(j) = 0 (j = 1, 2, 3) at 10.25, 16.25, 57, 60, and zeta'''' is bounded with jumps there.
- v = zeta·g. W = e^{-psi}, psi' = (1 - s)·4x/(1 + x^2) + s·phi', s = S7((x - 30.25)/8), phi = (1 + x^2)^{2/3}/4, rho = e^{-phi}.
- V_j = ||v^(j)||_W (L^2(W dx) on [10.25, 60]). Free variables V0, V2, V4; V1 <= (V0(V2 + c1W·V0))^{1/2}, V3 <= (V2(V4 + c1W·V2))^{1/2}, c1W = 0.8477 (`rhs_iv.py:94`).
- rho-norms of w~ = chi·w (chi = 1 for sigma <= 29, so w~ = w on IR): A, W1, B, W3, D4 = ||w~^(j)||_rho, j = 0..4.
- E_I = V0^2 + eps_I·V2^2, eps_I = 3e-5; Lambda = E + mu·E_I, mu = 0.08.
- Profile ODE: k' = m, m' = p/4, p' = -kq, q' = 1 + kp. theta' = k.

### 1.2 Derived sup quantities used below
Each is an interval function of the verifier variables. "rad" is its degree in (D4, V4) (radial mode of `tail_cert.py`); "val" its x-valuation.

| Name | Meaning | Bound in the verifier's variables | rad | val | Code |
|---|---|---|---|---|---|
| S0r, S1r, S2r | sup over IR of \|w\|, \|w'\|, \|w''\| | `mx(mid[j], regS["12<x<=16.25"][j])`, j = 0, 1, 2 (K_rho Sobolev, regions [8, 12] and [12, 16.25]) | 0, 0, 1/4 | 1 | exists (`rhs_iv.py:96-99, 155`); S2r new as a name |
| sI1r, sI2r | sup over [16.25, 34.25] of \|v\|, \|v'\| (= \|w'\|, \|w''\| there) | `gm(V0, KW(PSI1W·V0 + 2V1))`, `gm(V1, KW(PSI1W·V1 + 2V2))` **before** the max with S0V, S1V/32 | 0 | 1 | exists (`:100`), but only the max is kept (`:102`) |
| sZ | sup over Z1 of \|g'\| | max(sI2r, S1V/32) | 0 | 0 | = `supv1` (`:141`) |
| sOR | sup over OR of \|v'\| | S1V/54.75 + (35/48)·S0V = 0.0221 (v' = zeta g' + zeta' g; \|zeta'\| <= 35/48 on OR [exact]; y >= 54.75 on OR). [GI 2-15] quotes "up to 0.030", i.e. with S1V/32; either is above the coded 0.01875 | 0 | 0 | **new** |
| S1z | sup over IR of \|zeta g'\| | (2·KWIR·V1·V2)^{1/2} + (35/96)·S1r, KWIR = 1/W(16.25) (integration from the left end, see below) | 0 | 1 | **new** |
| S2z | sup over IR of \|zeta g''\| | (2·KWIR·V2·V3)^{1/2} + 2(35/96)·S2r + Z2i·S1r, Z2i = sup_IR\|zeta''\| | 1/4 | 1 | **new** |
| sg1 | sup over [10.25, 60] of \|g'\| | max(sZ, S2r) (on OR \|g'\| <= S1V/54.75 <= sZ) | 1/4 | 0 | **new** as a name |
| sG | sup over [10.25, 60] of \|v'\| | max(sI2r, S1V/32, sOR, S2r + (35/96)·S1r) | 1/4 | 0 | **new** (code uses `supv1`, which misses sOR: 0.0221 > 0.01875 [hand]) |
| supv2 | sup over [16.25, 34.25] of \|v''\| | `gm(V2, KW(PSI1W·V2 + 2V3))` | 1/4 | 1 | exists (`:177`) |
| sw16 | \|w(16.25)\| | `sr[0]` (K_rho(16.25) Sobolev) | 0 | 1 | exists (`:99`) |
| S_all34 | sup over [0, 34.25] of \|w\| | max over regions of `regS[.][0]`, far one = sw16 + Phi(34.25)^{1/2}·V0 | 0 | 1 | exists as `S_all` (`:166-167`) but is used as if it covered [0, 60] |

Sobolev rule behind the V-rows: for f with f·W^{1/2} -> 0 at the right end of supp, f(x)^2 W(x) = int_x^inf W(psi' f^2 - 2 f f') <= sup(psi')·||f||_W^2 + 2||f||_W·||f'||_W. So sup over [a, b] of f^2 <= W(b)^{-1}(...) [hand; this is how KW = 1/W(right end) arises: 1/W(20) = 401^2 = 1.608e5 against the coded 1.609e5]. Here sup(psi') is over the whole support (PSI1W = 1.305), not over the region.
For S1z, S2z integrate from the left instead: f = v' or v'' vanishes at 10.25, and f(x)^2 W(x) = int_{10.25}^x W(2 f f' - psi' f^2) <= 2||f||_W·||f'||_W because psi' >= 0; W is decreasing, so sup over IR of f^2 <= 2·KWIR·||f||_W·||f'||_W [hand].

### 1.3 Standing hypotheses (boxes). v4 must **assert** each one at every march step and in the tail
| # | Hypothesis | Where it comes from | fin12f4 status |
|---|---|---|---|
| B-a | alpha <= ramp/far table level (2.5e-4 at the production points) | barrier | `ramp_level()` raises otherwise [read `constants.py:49-52`] |
| B-b | S0r, S1r <= eta_IR (IR jet box). v4: eta_IR = 1e-3 | E via K_rho | logged, not asserted [read `closure.py:259-260`]; max 7.474e-4 [cert log] |
| B-c | \|w'\| <= 0.02 on [16.25, 60]: sI1r <= 0.02 and S0V <= 0.02 | E_I; (A) interface (I2) | logged, **not asserted** [read `closure.py:261-262`] |
| B-d | \|w\| <= 0.02 on [16.25, 26], <= 0.5 on [26, 60] (U0-box of the far coefficients, `r60_t3.py:20`) | on [16.25, 34.25]: sw16 + Phi(34.25)^{1/2}V0 <= 5.91e-4 [cert log `far_w`]; on [34.25, 60]: \|w(34.25)\| + 25.75·S0V <= 0.395 [hand] | first part logged; **second part checked nowhere** |
| B-e | S0V = sup\|w'\| <= 0.0153, S1V = sup y\|w''\| <= 0.6 on [34.25, 60] | bootstrap input from (A), interface (I2) | env of `run_v33.sh` |
| B-f | modulation validity box: alpha <= 5e-4, A, W1 <= 1e-3, strip sups | agent_mod/ibp | post-check in `mod_sups.py` docstring only |

## 2. The v-equation

### 2.1 The w-equation [read S/BOOTSTRAP.md:8-11, S/ASSEMBLY.md:70, `agent_int/v2/int_v56.py:4-6`, `agent_int/planB/pb_nlin.py:15-20`]
u_tau = F(u) - a(p + u - u'q/al_g) + b(sin(theta) + u'cos(theta)/al_g), al_g = 1 - k·u, and alpha_tau = lamt·alpha + beta. Hence

  w_tau = L w + Phi_F + N(w) + M,
- L = -d^4 - k^2 d^2 - (5km + q/4) d + c0, c0 = 1/4 - 3m^2 - kp (`graph_eq.py`);
- Phi_F = F(u0) - lamt·(alpha·psi3 + 2 alpha^2·h2): the slow-manifold residual; depends on x and alpha only;
- N(w) = F(u0 + w) - F(u0) - L w;
- M = -beta·psi3 - 2 alpha beta·h2 - a·psi1 - b·psi2 + a·A1(u) + b·B1(u), psi1 = p, psi2 = -sin(theta), A1(u) = -u + u'·q/al_g(u), B1(u) = u'·cos(theta)/al_g(u).

[hand] I did not re-derive the graph equation or the sign conventions of (a, b); I take them from the files above. Only absolute values of a, b, beta are used below, so their signs do not matter.

### 2.2 Structure of N [read `S/graph_struct.py`, `graph_struct.pkl`, `pb_nlin.py:15-27`]
F is a polynomial in (u'', u''', u''''), with monomials (0,0,1), (0,1,0), (1,0,0), (1,1,0), (2,0,0), (3,0,0), (0,0,0) only. With U0 = u0 + w, U1 = u0' + w', a_j = u0^(j), J^2 = (1 - k U0)^2 + U1^2:

  N = c·w'''' + N_low,  **c = 1 - J^-4** (exact: the u'''' coefficient of F is -J^-4 [read pkl]),
  N_low = d010·w''' + d100·w'' + d200·w''^2 + d110·w'' w''' + d300·w''^3 + d000,

where d_tj0 = d_tj0(x; U0, U1, a2, a3, a4) (exact polynomial coefficients, no box on w'', w''', w''''), d110 = 10 J^-6 U1, d300 = 3 J^-6 (1 - 6 J^-2 U1^2), and d000 = G(x, U0, U1) - G(x, u0, a1) - c0·w + (5km + q/4)·w' with |d000| <= A0|w| + A1|w'| (mean value over the box; `pb_nlin.py:25-27`).

### 2.3 Differentiating and cutting off
[sym I9] (L w)' = L1 g + c0'·w, with L1 = -d^4 - k^2 d^2 - (7km + q/4) d + e0, e0 = -8m^2 - (5/2)kp. So

  g_tau = L1 g + c0' w + Phi_F' + (c g''')' + N_low' + M'.

Multiply by zeta and commute zeta only where needed:
- zeta·L1 g = L1 v - [L1, zeta] g;
- [sym I2] zeta·(c g''')' = (c v''')' - c·K4 - c'·K3, with K3 = [d^3, zeta]g = 3 zeta' g'' + 3 zeta'' g' + zeta''' g and K4 = [d^4, zeta]g = 4 zeta' g''' + 6 zeta'' g'' + 4 zeta''' g' + zeta'''' g;
- [sym I11 and its analogue for B1] A1(u) = A1(u0) - w + g·Q_u + E_A·w, Q_u = q/al_g(u), E_A = u0' q k/(al_g(u) al_g(u0)); B1(u) = B1(u0) + g·C_u + E_B·w, C_u = cos(theta)/al_g(u), E_B = u0' k cos(theta)/(al_g(u) al_g(u0)). Write Q_u = q + dQ, C_u = cos(theta) + dC, with dQ = q k u/(1 - k u), dC = cos(theta) k u/(1 - k u). Then

  zeta·(A1(u))' = zeta·A1(u0)' + [q v' + (q' - 1) v] + (dQ·v)' - zeta'·Q_u·g + zeta·(E_A w)',
  zeta·(B1(u))' = zeta·B1(u0)' + [cos(theta) v' - k sin(theta) v] + (dC·v)' - zeta'·C_u·g + zeta·(E_B w)'.

**The v-equation** (exact [hand, with the [sym] pieces above]):

  v_tau = L1 v                                                   (T1)
        - [L1, zeta] g                                           (T2)
        + (c v''')'                                              (T3)
        - c·K4 - c'·K3                                           (T4)
        + zeta·N_low'                                            (T5)
        + zeta·c0'·w                                             (T6)
        + zeta·Phi_F'                                            (T7)
        + a[q v' + (q' - 1)v] + b[cos(theta) v' - k sin(theta) v]   (T8)
        - zeta'(a Q_u + b C_u) g                                 (T9)
        + a·zeta·A1(u0)' + b·zeta·B1(u0)'                        (T10)
        - beta·zeta·(psi3' + 2 alpha h2') - a·zeta·psi1' - b·zeta·psi2'   (T11)
        + a[(dQ v)' + zeta (E_A w)'] + b[(dC v)' + zeta (E_B w)']        (T12)

T2, T4, T9 are supported on IR and OR. This is ASSEMBLY.md:70 with T3 + T4 + T5 = zeta·N' written out and M' expanded.

## 3. Energy identities

- Level 0: (1/2) d/dtau V0^2 = <v, v_tau>_W.
- Level 2: (1/2) d/dtau V2^2 = <v'', v_tau''>_W = <D2v, v_tau>_W, where D2v := (W v'')''/W = v'''' - 2 psi' v''' + (psi'^2 - psi'') v'' [sym I3, I4]. The two boundary terms vanish because v_tau = zeta·g_tau and (v_tau)' vanish at 10.25 and 60.
- ||D2v||_W <= V4 + 2 PSI1·V3 + PSI2·V2 =: **DV**, with 2 PSI1 <= 2.6097 and PSI2 = sup|psi'^2 - psi''| <= 1.69544 on [10.25, 60] [cert G0]. Code: `D2v = V4 + 2.61 V3 + 1.7 V2` (`rhs_iv.py:176`), covered.
- Level-0 integration by parts used for T5 and T7 [sym I10]: <v, zeta f'>_W = -int (W zeta v)' f, and (W zeta v)' = W(zeta v' + (zeta' - psi' zeta) v). No boundary term (zeta = 0 at both ends).
- ||zeta v' + (zeta' - psi' zeta) v||_{W, zone} <= V1 + KAP_zone·V0, KAP_zone = sup_zone(|zeta'| + psi' zeta): IR 0.5340, Z1 1.2827, OR 1.5883 [cert G0]. Code: `Vq = V1 + 1.67 V0`, covered.

**Regularity (not written anywhere; needs a short lemma, U-1 in §11).** u0 is C^4 with bounded piecewise fifth derivative; zeta is C^3 with bounded zeta''''. So v is C^3 with v'''' in L^inf, and v_tau in H^2-type classes. The identities above then hold for a.e. tau, with Lambda absolutely continuous, provided the solution is smooth for tau > tau_1. FLUX_REPORT §1 uses the same form ("H^4-continuous").

**Zone rule.** No integration by parts is done on a sub-zone. Integrals of terms that are not expressed through v are split as int_IR + int_Z1 + int_OR after the (global) identities above, and each piece is bounded by Cauchy-Schwarz on its zone. On Z1, g^(j) = v^(j), so ||1_Z1 g^(j)||_W <= V_j. On IR, ||1_IR f||_W <= RWR·||f||_rho with RWR = sup_IR (W/rho)^{1/2} <= 0.65580 [cert G0] and g^(j) = w~^(j+1). The OR pieces are collected in §4.10.

## 4. Term by term

Format: identity; bound at level 0 (multiplies mu); bound at level 2 (multiplies mu·eps_I); constants; code status. "H" / "L" is the group of `rhs_pn` (x-degree >= 2 / <= 1).

### 4.1 T1: linear part
- Level 0: <v, L1 v>_W <= -MU1·V0^2 - (1/2) V2^2, MU1 = 0.1552. Code `:189`. **Provenance to correct** [GI 0-1]: the L-int certificate it cites (`agent_lint/runs/cert_PX2_XA10.25_A230.25.json`) gives 0.155167 < 0.1552; the typed value is valid only through the independent `agent_int/v2/v2_T4_P2_0.25.json` (0.15552). v4: read MU1 from the v2_T4 JSON, rounded **down** (it is a lower bound), and cite that file. I did not recheck the two numbers myself.
- Level 2: (L1 v)'' = M1 v'' + e1 v' + e0 v [read `agent_int/v2/int_garding.py:1-10`], so
  <v'', (L1 v)''>_W <= -(1/2) V4^2 - MU2·V2^2 + **E1·V1 V2 + E0·V0 V2**, MU2 = 0.6555, E1 = sup|e1| <= 5.66e-6, E0 = sup|e0| <= 5.11e-6 [cert `v2_T4_P2_0.25.json`].
  Code: the e-terms are **absent** (`:176-189`). ASSEMBLY.md:86 lists them with the wrong monomials ("V3V2 + V2^2"). **New** (H).

### 4.2 T2: linear commutator
- IR, level 0: -<v, [L1, zeta]g>_IR = int (A_z g^2 + B_z g'^2) <= AF·W1^2 + BF·B^2 [cert identity `agent_coup/verify_identity.py`; raw a_f 0.069555037634, b_f 0.039358915436 in `coupling_constants.json`]. Code `:147`, with AF typed **below** raw. **Corrected: outward rounding.**
- IR, level 2: |<D2v, [L1, zeta]g>_IR| <= (CF0·D4 + CF1·W3 + CF2·B + CF3·W1)·DV, CF_j = sup_IR |coef_j([L1, zeta])|·(W/rho)^{1/2} [cert `run_coup.py:71-72`; raw c3..c0 = 0.442142686616, 0.525899483214, 0.423527657147, 0.425017231881]. Code `:178`, CF0, CF2, CF3 typed below raw. **Corrected: outward rounding.**
- OR, both levels: -> K0^OR, K2^OR (§4.10; F1 done for the pointwise parts).

### 4.3 T3: quasilinear top order (the audit's N1)
c = 1 - J^-4 on the whole support. Write c' = c_x + c_U0·U1 + c_U1·(a2 + g'), so |c'| <= CA + CB·|g'| pointwise, CA := sup|c_x + c_U0 U1 + c_U1 a2|, CB := sup|c_U1| over the box (per zone). Put CBAR := sup|c| (this is the existing "d001") and

  **C1 := CA + CB·sg1** (sup of |c'| on the support; sg1 = sup|g'| from §1.2. Note: sg1, not sG; sG = sup|v'| is used only in the GN identity of §4.5).

- Level 0 [sym I1]: <v, (c v''')'>_W = int W c v''^2 + int W (c' - 2 psi' c) v' v'' - int W ((psi' c)' - psi'^2 c) v v''
  <= CBAR·V2^2 + (C1 + 2 PSI1·CBAR)·V1 V2 + (PSI1·C1 + PSI2·CBAR)·V0 V2.
  Code: `d001 * V2V2` only (`:147`). **Corrected** (two new monomials, H).
- Level 2: <D2v, (c v''')'>_W = int W D2v (c v'''' + c' v''') <= **CBAR·V4·DV + C1·V3·DV**.
  Code: `d001 * V4 * V4` (`:178`). **Corrected** (H).
- Zone-wise sups are allowed and recommended: replace each product of sups by the max over the zones {IR, [16.25, 34.25], [34.25, 60]} of the zone-wise products (psi' is largest where c is largest, so the gain is small; it matters only for the box sizes).
- Optional sharpening at level 2 (not in the baseline): int W c v''''(-2 psi' v''') = int (W c psi')' v'''^2, which replaces 2 PSI1·CBAR·V3V4 by sup|(c psi')' - psi'^2 c|·V3^2. Not needed [probe].
- Sizes at the present far box (alpha <= 1e-3, |w'| <= 0.02): CBAR = 2.27e-2 [cert `r60_T3_boxw1_0.02.json`]; CB ~ 4 sup|U1| ~ 0.42, CA ~ CB·sup|u0''| ~ 6.5e-4, C1 ~ 0.01 [est]. At the alpha-level 2.5e-4: CBAR ~ 3.5e-3, CB ~ 0.17 [est].

### 4.4 T4: quasilinear commutator -c K4 - c' K3
- IR, both levels, direct norm. Define (new constants, sup over IR, with (W/rho)^{1/2} inside the sup):
  k4 = sup|(4 zeta', 6 zeta'', 4 zeta''', zeta'''')|(W/rho)^{1/2} <= (0.44291, 0.52671, 0.42416, 0.42506) [cert G0 `binom_zeta_k_sqrtWrho`];
  k3 = sup|(3 zeta', 3 zeta'', zeta''')|(W/rho)^{1/2} [**new**; <= (3/4) k4_1, (1/2) k4_2, (1/4) k4_3 trivially].
  RQ := CBAR_IR·(k4_1 D4 + k4_2 W3 + k4_3 B + k4_4 W1) + (CA_IR + CB_IR·S2r)·(k3_1 W3 + k3_2 B + k3_3 W1).
  Level 0: <= V0·RQ. Level 2: <= DV·RQ.
  Code: level 0 has `0.66 * d001_ramp * D4 * Vq` inside term (a) (`:152-161`), which covers the D4 part only and belongs to the old decomposition; level 2 has nothing. **New** (H). Size: CBAR_IR <= 6.93e-5 at (2.5e-4, eta 1e-3) [cert `planB_ramp.json`], i.e. 7.7e-5 of CF [cert G0].
- OR -> §4.10.

### 4.5 T5 on Z1: lower-order nonlinear terms
Level 0: -int_Z1 W (v' - psi' v) N_low, with on Z1: w'' = v', w''' = v'', w' = v. Bound (H):

  (V1 + KAP_Z1·V0)·[ d010·V2 + d100·V1 + A1·V0 + d200·sZ·V1 + d110·sZ·V2 + d300·sZ^2·V1 + (hA0·sw16 + HA0·V0) ].

- The last bracket bounds ||1_Z1 A0 w||_W by |w(x)| <= |w(16.25)| + Phi(x)^{1/2} V0 on Z1 (Phi(x) = int_{16.25}^x 1/W; valid because g = v on Z1): hA0 := (int_Z1 W A0(x)^2 dx)^{1/2}, HA0 := (int_Z1 W A0(x)^2 Phi(x) dx)^{1/2}, A0(x) the box-sup function. **New** (tiny: A0 <= 9.2e-9 on [16.25, 20], <= 3.3e-11 beyond [cert r60]).
- Code `:147-148`: all but the A0 term exist, with `Vq` (1.67 >= KAP). **Unchanged except + A0 term.**

Level 2: <D2v, 1_Z1 N_low'>_W <= DV·||1_Z1 N_low'||_W. Total derivative, with Dd := d^a + d^b·g' (d^a := sup|d_x + d_U0·U1 + d_U1·a2 + sum_j d_{a_j} a_{j+1}|, d^b := sup|d_U1|; the fifth jet a5 = u0^(5) enters only through D(d000), because d000 contains a4):

  N_low' = d010 g''' + (Dd010 + d100) g'' + Dd100 g' + 2 d200 g' g'' + Dd200 g'^2 + d110 (g''^2 + g' g''') + Dd110 g' g'' + 3 d300 g'^2 g'' + Dd300 g'^3 + D(d000),
  D(d000) = E1c·g' + E0c(w, g), |E1c| <= A1 (it is d_U1 G(U) + (5km + q/4)), |E0c| <= B0·|w| + B1·|g| (mean value in (w, g) over the box; E0c(0, 0) = 0 [hand]).

Bound of ||1_Z1 N_low'||_W, monomial by monomial (highest derivative in L^2(W), the others in sup):

| Monomial | Bound | In code (`:179`)? |
|---|---|---|
| d010 g''' | d010·V3 | yes |
| (Dd010 + d100) g'' | (d010^a + d100)·V2 + d010^b·sZ·V2 | **no** |
| Dd100 g' | d100^a·V1 + d100^b·sZ·V1 | **no** |
| 2 d200 g' g'' | 2 d200·sZ·V2 | coded as 2 d200·s·**V3** (wrong index; V2 <= 15.35 V3 only by Hardy) |
| Dd200 g'^2 | d200^a·sZ·V1 + d200^b·sZ^2·V1 | **no** |
| d110 g''^2 | d110·sG·(3 V3 + PSI1·V2) [sym I8: ||v''^2||_W <= sup|v'|·(3V3 + PSI1 V2), global sup sG] | yes, but with `supv1` instead of sG |
| d110 g' g''' | d110·sZ·V3 | yes |
| Dd110 g' g'' | d110^a·sZ·V2 + d110^b·sZ^2·V2 | **no** |
| 3 d300 g'^2 g'' | 3 d300·sZ^2·V2 | coded with coefficient 1 |
| Dd300 g'^3 | d300^a·sZ^2·V1 + d300^b·sZ^3·V1 | **no** |
| E1c g' | A1·V1 | **no** |
| E0c | B1·V0 + hB0·sw16 + HB0·V0 | **no** |

**Corrected and extended** (H). The derivative coefficients need a new certified stage ("t3d", §6). u0^(5) is available: the C^4 functions of `c4fun.py` have exact jets of any order on the pieces and on the extension [read `agent_int/v3/NOTES_v3.md`].
Sizes [est] from the pkl structure: d110^b = 10, d300^b ~ 36 sup|U1| ~ 3.8, d200 ~ 9 sup|u0''|, d010 ~ 10 sup|U1| sup|u0''|; all others small.

### 4.6 T5 on IR: lower-order nonlinear terms in rho-norms
Level 0: |int_IR (W zeta v)' N_low| <= (V1 + KAP_IR·V0)·RWR·[ d010·W3 + d100·B + A1·W1 + A0·A + d200·S2r·B + d110·S2r·W3 + d300·S2r^2·B ], coefficients = IR sups.
- Code: term (a) (`:152-161`) has 0.66·(d001 D4 + d010 W3 + d100 B + A1 W1 + A0 A)·Vq. In this spec the D4 entry moves to §4.4 and **the three nonlinear monomials are new** (IR sups 0.01331, 0.2968, 3.024 at the 1e-3/1e-2 box [cert `planB_T3_0.001_0.01_0.01_1.json`]).
- **The table interpolation c0 + dw·S_ramp of `constants.py:39-48` is not justified** [hand]: a chord through eta = 0 and eta = 1e-2 bounds the sup over the eta-box only if that sup is convex in eta, which is not proved (for A1 the chord value at 1e-3 equals the certified value to 4 digits). Fix: use the certified table row eta = 1e-3 as a constant and assert B-b (S_ramp <= 7.474e-4 [cert]).

Level 2: |<D2v, 1_IR zeta N_low'>_W| <= DV·RWR·R5, with (g^(j) = w^(j+1); zeta <= 1 except where marked)

  R5 := d010·D4 + (d010^a + d100)·W3 + d010^b·S2r·W3 + d100^a·B + d100^b·S2r·B + 2 d200·S2r·W3 + d200^a·S2r·B + d200^b·S2r^2·B
        + d110·(**S1z·D4** + **S2z·W3**) + d110^a·S2r·W3 + d110^b·S2r^2·W3 + 3 d300·S2r^2·W3 + d300^a·S2r^2·B + d300^b·S2r^3·B + A1·B + B1·W1 + B0·A.

- **Why S1z, S2z.** The monomial d110·(g' g''' + g''^2) = d110·(w'' w'''' + w'''^2) bounded as S2r·D4 has radial degree 1/4 + 1 and, times DV, 2.25 > 2. That breaks the radial lemma (POS(r)/r^2 nonincreasing), whatever the size of the coefficient. Using the cutoff, zeta g' = v' - zeta' g, so sup_IR|zeta g'| <= S1z, which is bounded by the energy (radial degree 0); then d110·S1z·D4·DV has degree 2. Likewise zeta g'' = v'' - 2 zeta' g' - zeta'' g gives S2z (degree 1/4) times W3 (degree 1/2): 1.75. [hand; the probe's `tail_cert` passes with these forms.]
- Code: nothing on IR at level 2 beyond CF. **New** (H). Size: d010_IR <= 3.2e-5 against CF0 = 0.442; d110_IR·S1z ~ 0.3·1e-3 [est].
- **Answer to [GI] 2-10 ("missing, and not derivable as the code is organised").** GI names two obstacles. (i) "zeta D' needs w^(5), which no energy controls": true for zeta·N' taken as a whole; in this spec the only fifth-order piece, zeta c g'''', is never bounded on the ramp. It is part of (c v''')' (§4.3, global, paired as CBAR·V4·DV), and what remains on IR is c K4 + c' K3 (§4.4) and N_low' (here), both of order <= g''' = w''''. So D4 is the highest rho-norm needed, and it is available (E2 dissipation). (ii) "sup|w''|·D4·||D2v|| has radial degree 2.25": avoided by S1z, S2z above. GI proposes the same cure ("one factor converted to v", degree 1.75); the form here is slightly different (the factor converted is the sup factor, so the D4 term has degree exactly 2 and the W3 term 1.75), and both pass `tail_cert` (GI's `g0all_tail.out`; my probe). Every other monomial of R5 has degree <= 2 without help (list in §5.3).
- What remains open for this block is not the formulation but the inputs: d^a, d^b, B0, B1 on IR have no certified values anywhere in the tree [GI 2-5, 2-9 say the same for the far field]. They are stage C2.

OR -> §4.10.

### 4.7 T6: zeta c0' w
- Level 0: <= V0·C0W, level 2: <= DV·C0W, with
  **C0W := C0P·S_all34 + C0Hb·V0**, C0P = ||zeta c0'||_{W,[10.25,57]} <= 9.9e-12 [cert v2 V5: 9.875e-12 for the whole support], C0Hb := (int_{34.25}^{57} W c0'^2 (Phi(x) - Phi(34.25)) dx)^{1/2} [**new**, est < 1e-30]. (On [34.25, 57]: |w(x)| <= |w(34.25)| + (Phi(x) - Phi(34.25))^{1/2} V0, since g = v there.)
- Code: level 0 `c0pn * S_all * V0` (`:164-168`): S_all is a sup over [0, 34.25] only (PHI34), so [34.25, 60] is **not covered** as coded. Level 2: absent. **Corrected / new** (H).
- OR piece (w on [57, 60] is not controlled by V0) -> §4.10. Size ~ 1e-100.

### 4.8 T7: forcing
- Level 0: <= FNW(alpha)·(V1 + KAP·V0), KAP = 1.5883 <= 1.67. FNW for the C^4 functions: 1.243e-8·al + 4.03e-8·al^2 + 0.1823·al^3 [cert `int_constants_v3.json`]. Code `:150` uses (1.09e-7, 6.5e-9, 0.1823): an upper bound for al <= 2.8 [cert arithmetic, AUDIT 0.4], but **type the C^4 triple** in v4. (L)
- Level 2: <= FNW2(alpha)·DV, FNW2 = 2.70e-8·al + 1.156e-7·al^2 + 2.603e-3·al^3 [cert v3 F1a/F1b]. Code `:181`: `al^3 * D2v`, not an upper bound for al < 1.6e-4. **Corrected** (L).
- The R3 correction (fix4_r3) concerns (A)'s forcing r0, not this term [read R3_REPORT §2].

### 4.9 Modulation
Bounds for the parameters, as interval functions (existing: `Amx`, `Amc`, `:143-145`):
- |(a, b)| <= AB := ABx + ABc, ABx = (CAB0 + MAB·al)·A + QAB·Q, ABc = PHIAB·**Fn(al)** + strip, Q = (A + W1 + B)^2. The modulation lemma has ||chi Phi_F||_rho here, i.e. Fn (`constants.py:30`), as `mod_sups.py:25` uses. `rhs_iv.py:144` uses FnW (the W-norm): too small by 7.3x at small alpha [GI 0-14]. **Corrected.** Size: 0.2% of |(a, b)| [GI].
- |beta| <= BT := BTx + BTc, BTx = (BETA0 + BETA1·al)·A + QB·Q + **BTs**, BTc = BETAF·Fn(al), where **BTs is the strip part as a linear form in the state**: BTs := sb0·S0s + sb1·sI1r + sb2·sI2r + sb3·supv2, S0s = sw16 + Phi(34.25)^{1/2} V0. Coefficients from `agent_mod/ibp/mod_ibp.json` key `strip.beta`: linear_alpha_S0 4.33e-8 (times al, to be confirmed), commutator_S0_S3 = (1.54e-9, 1.30e-9, 7.75e-10, 1.89e-10), plus the quadratic and cubic coefficients times the box 0.02 [read]. **This restatement of the modulation lemma is not written anywhere** (U-3).
  **Why not the constant STRIP_BETA, and what each choice costs.** beta multiplies ||zeta psi3'||_W·V0 with no factor alpha, so a constant in beta gives a non-decaying term 0.2596·1.06e-9·V0 = 2.75e-10·V0 of x-degree 1 and alpha-degree 0. Three ways to treat it:
  1. *State form (this spec).* BTs is linear in the state, the term is quadratic (H), no floor, no change of F0 or of the J-slopes. In size it is negligible: coefficients ~1e-9 times the Sobolev factor KW^{1/2} ~ 1.5e3, i.e. ~1e-5 of the E_I dissipation [est]. Cost: the modulation lemma must be restated with its strip part as a linear form (U-3; text and a small script, no rerun of `agent_mod`), and S_3 must not be needed as a box (SF-2).
  2. *Constant kept, absorbed by Young against mu·MU1·V0^2.* Floor contribution 9.76e-21 per unit tau, i.e. **F0 x7.4, sqrt(L_F) 9.87e-10 -> 2.87e-9, J-slopes x8.4** [cert GI 2.1]. Closure is not at risk, but the slopes feed (A) (H-J-W window form, kappa_J; C6 of ROUTE §5 has a reserve of only ~6e-7 for the J constants, and the slopes enter the window bounds), alpha_F triples (alpha-condition margin 53.7 -> ~18 [est]), and the tail lemma needs a new slot: `closure.py:143-147` has none for a degree-1 floor when QUAD = 1 [GI 4.4.2], and the Young step must be done inside the cell bound (the same construction as `unc.Phi`).
  3. *Constant kept, through the existing f1 slot of `tail()`* (f1·sqrt(Lambda) against rho instead of mu·MU1): floor sqrt(L_F) ~ 2 f1/rho = 5.8e-8 with f1 = sqrt(mu)·2.75e-10 [hand], above the peak. Not usable; and `closure.py:195` asserts f1 = 0.
  Recommendation: option 1. If U-3 cannot be proved as stated, option 2 is the fallback, at the cost above plus a larger rerun of the (A) chain's window data.
  The 0.657·Fn(al) part of beta (BTc) is a genuine forcing in L: 2.73e-11·V0 at al = 1.5e-4, which doubles the I0 forcing (coded 1.67·FnW = 2.84e-11) [cert GI 2.1]. It is in both my probe and GI's rerun.

Terms:

| Term | Level 0 | Level 2 | Constants | Code |
|---|---|---|---|---|
| T8 main [sym I5, I6, I7] | AB·T0·V0^2 | AB·(T2·V2^2 + 2 Q2B·V1 V2 + Q2B·V0·(V3 + PSI1·V2)) | T0 = sup\|(psi' q/2 + q'/2 - 1, psi' cos/2 - k sin/2)\|_2 <= 37.185; T2 = sup\|(psi' q/2 + 5q'/2 - 1, psi' cos/2 - 5 k sin/2)\|_2 <= 39.1842 [cert G0, computed there with the slightly larger b-forms psi'/2 + 2\|k\| and its level-2 analogue; recompute with these exact forms]; Q2B = sup\|(q'', (cos theta)'')\|_2, q'' = mp - k^2 q [**new**, tiny] | level 0 `CAB3985` 39.16 >= 37.185: ok. Level 2: 39.16 **< 39.1842**, compensated only by the unexplained `c62·Am·V2·V3` through an unwritten Hardy step [AUDIT 2.6a]. **Corrected**: type T2 honestly; delete c62 |
| T9 on IR | AB·TC0·W1^2 | AB·TC2·W1·DV | TC0 = sup_IR\|(q, cos) zeta zeta' W/rho\|_2 <= 0.4336 + 0.0726; TC2 = sup_IR\|(q, cos)\| \|zeta'\| (W/rho)^{1/2} <= 1.2709 [cert G0], with Q_u, C_u instead of q, cos: factor (1 + 2e-7) [est] | level 0 `CAB25` (2.5): ok. Level 2: **absent**. **New** |
| T9 on OR | -> K0^OR | -> K2^OR | | absent |
| T10 | AB·al·NAB·V0 | AB·al·NAB·DV | NAB = \|(0.2731, 0.01589)\|_2 = 0.2736 [cert v2 V5 / v3 F3] | ok (`CABQ` 0.2889, `CAB03` 0.3). ABc part in L, ABx part in H, as now |
| T11 beta | BT·(NP3 + 2 al·NH2)·V0 | BT·(NP3 + 2 al·NH2)·DV | NP3 = \|\|zeta psi3'\|\|_W <= 0.2596, NH2 <= 1.796 [cert v2 V5] | **absent (N2)**. BTx part in H, BTc part in L |
| T11 psi1, psi2 | ABx·NP12·V0 + [ABc·NP12·V0] | same with DV | NP12 = \|(1.709e-7, 1.683e-8)\|_2 [cert v2 V5]; [GI] uses the sum 1.88e-7 | **absent**. The ABc part is an x-degree-1 term with a constant strip piece: replace the constant strip in AB by its state form (as for beta), or carry the strip·NP12 product (5.5e-9·1.72e-7 ~ 1e-15) as a floor f1; the first is cleaner |
| T12 (dQ v)', (dC v)' | AB·(DQ0·V0·V1 + DQ1·V0^2) | AB·(DQ0·V1 + DQ1·V0)·DV | DQ0 = sup\|(dQ, dC)\|_2, DQ1 = sup\|(dQ', dC')\|_2 over the box; \|k q\| <= 1.3e-3 on IR, \|U0\| <= 1e-3 [**new**, est 1e-6] | **absent (AUDIT 0.6f)** |
| T12 zeta (E_A w)', zeta (E_B w)' | AB·V0·(EA0·V0 + RWR·EA1_IR·A + hE·sw16 + HE·V0) | same with DV | EA0 = sup\|(E_A, E_B)\|_2, EA1 = sup of the x-derivative (it contains u'' through al_g'; bound with sG), hE, HE Hardy integrals as in §4.5 [**new**, tiny] | **absent** |

Notes.
- T8 level 2: [sym I6] gives exactly W[(psi' Q/2 + 5Q'/2 - 1) v''^2 + 2 Q'' v' v'' - Q'' v v''' + psi' Q'' v v''] up to a total derivative, with Q = q (profile only, so Q'' needs no bound on w). The u-dependence of Q_u is entirely in T12.
- Euclidean pairing: |a f + b h| <= |(a, b)|·(f^2 + h^2)^{1/2} pointwise, then sup. cos(theta) and sin(theta) may be bounded by 1 (as in v2 V5).
- The norms NP3, NH2, NAB, NP12 were certified for the C^3 functions; for the C^4 ones the change is ~1e-7 relative [read NOTES_v3]. v4 should recompute them with `c4fun` (U-6).

### 4.10 Outer ramp: the interface with the flux lemma (unit G2 of ROUTE_DECISION)
Everything on [57, 60] that is not expressed through v goes into two functionals. With g, w and the coefficient boxes of (A) on [57, 60]:

  **K0^OR** := int_57^60 W v·( -[L1, zeta]g - c K4 - c' K3 - zeta'(a Q_u + b C_u) g + zeta c0' w + zeta (a E_A' + b E_B') w ) dx - int_57^60 (W zeta v)' N_low dx,

  **K2^OR** := int_57^60 W D2v·( -[L1, zeta]g - c K4 - c' K3 - zeta'(a Q_u + b C_u) g + zeta N_low' + zeta c0' w + zeta (a E_A' + b E_B') w ) dx.

Then the verifier's inequality is: (1/2) d/dtau V0^2 <= [§4.1-4.9, level 0] + K0^OR and (1/2) d/dtau V2^2 <= [§4.1-4.9, level 2] + K2^OR.

What stays in the verifier on OR (global V-norm terms, no flux): T1, T3, T7, T8, T10, T11, the (dQ v)' and E_A v parts of T12, and the coefficient sups CBAR, C1 (which use B-c, B-e on OR), sOR inside sG.

Consumption, as ROUTE §2.1, §5 decided (unchanged by this spec):
- pointwise parts of K0^OR + eps_I·K2^OR <= FLUX0 + eps_I·FLUX2 (the verifier adds mu·(FLUX0 + eps_I·FLUX2) to L, `:170-172, 183-185`; floor F0, `closure.py:144`);
- the g''' pieces (T3 of F1, the (3,·) monomials of N_low' and of c K4) and the endpoint term go through the time-integrated lemma (G3) with the level-1 and level-2 exports of (A).

Relation to AUDIT §5: item lists agree. Differences: (i) the piece d001·v'''' is **not** in K2^OR (it is inside T3, global); what remains is the commutator c K4 + c' K3, so the "N1 on [57, 60]" remark of AUDIT §5 level 2 item 2 is taken care of in the verifier; (ii) three formal items are added: zeta c0' w, A0·w (inside N_low) and the E_A', E_B' terms on OR. All three carry k, m, p at sigma >= 54.75 (<= 1e-43 [cert hstail bounds quoted in FLUX_REPORT §2]).

**FLUX0, FLUX2 stay inputs.** v4 does not need G2 to exist; it needs the same statement as fin12f4 ("K^OR <= FLUX pointwise" as hypothesis (I3)), now with K^OR defined as above.

## 5. The assembled inequality

### 5.1 Level 0 (multiplied by mu)
```
NEG0 = (1-theta)(MU1 V0^2 + 1/2 V2^2)
I0p (H) =
   AF W1^2 + BF B^2                                                        [4.2]
 + CBAR V2^2 + (C1 + 2 PSI1 CBAR) V1 V2 + (PSI1 C1 + PSI2 CBAR) V0 V2      [4.3]
 + V0 * RQ                                                                 [4.4]
 + (V1 + KAP_Z1 V0) * (d010 V2 + d100 V1 + A1 V0 + d200 sZ V1 + d110 sZ V2 + d300 sZ^2 V1 + hA0 sw16 + HA0 V0)   [4.5]
 + (V1 + KAP_IR V0) * RWR * (d010i W3 + d100i B + A1i W1 + A0i A + d200i S2r B + d110i S2r W3 + d300i S2r^2 B)   [4.6]
 + V0 * C0W                                                                [4.7]
 + AB T0 V0^2 + AB TC0 W1^2 + ABx al NAB V0                                [4.9: T8, T9, T10]
 + BTx (NP3 + 2 al NH2) V0 + ABx NP12 V0                                   [4.9: T11]
 + AB (DQ0 V0 V1 + DQ1 V0^2) + AB V0 (EA0 V0 + RWR EA1i A + hE sw16 + HE V0)   [4.9: T12]
L0 =
   FNW(al) (V1 + KAP V0)                                                   [4.8]
 + ABc al NAB V0 + BTc (NP3 + 2 al NH2) V0                                 [4.9]
 + FLUX0                                                                   [4.10]
```
### 5.2 Level 2 (multiplied by mu·eps_I)
```
NEG2 = (1-theta)(1/2 V4^2 + MU2 V2^2)
I2p (H) =
   E1 V1 V2 + E0 V0 V2                                                     [4.1]
 + (CF0 D4 + CF1 W3 + CF2 B + CF3 W1) DV                                   [4.2]
 + CBAR V4 DV + C1 V3 DV                                                   [4.3]
 + RQ DV                                                                   [4.4]
 + (table of 4.5) DV                                                       [4.5]
 + RWR R5 DV                                                               [4.6]
 + C0W DV                                                                  [4.7]
 + AB (T2 V2^2 + 2 Q2B V1 V2 + Q2B V0 (V3 + PSI1 V2)) + AB TC2 W1 DV + ABx al NAB DV    [4.9]
 + BTx (NP3 + 2 al NH2) DV + ABx NP12 DV + AB (DQ0 V1 + DQ1 V0) DV + AB (EA0 V0 + RWR EA1i A + hE sw16 + HE V0) DV
L2 =
   FNW2(al) DV + ABc al NAB DV + BTc (NP3 + 2 al NH2) DV + FLUX2
```
with DV = V4 + 2 PSI1 V3 + PSI2 V2. The E, E2 blocks of `rhs_pn` (`:106-139`) are outside this spec (but see §11).

### 5.3 Structural conditions (must be re-certified by `tail_cert.py` on the v4 code)
1. Radial: every H and L term has degree <= 2 in (D4, V4). Degrees: W3, V3: 1/2; S2r, S2z, supv2, sG: 1/4; S1z, sZ, sI1r, sI2r, S0r, S1r, sw16: 0; DV: 1. Critical terms: CF0·D4·DV, d010i·D4·DV, d110i·S1z·D4·DV, CBAR·V4·DV (all exactly 2); 3 d300i·S2r^2·W3·DV and d110^b·S2r^2·W3·DV (exactly 2).
2. s-mode: every H term has x-valuation >= 2. Terms that contain sZ or sG (valuation 0, because of the constants S0V, S1V) must carry two further state factors; all terms above do. The strip sups in BTs must be the raw sI1r, sI2r (valuation 1), not the maxima with S0V, S1V (my first probe failed `tail_cert` exactly here).
3. Joint mode (floors off): every L term has joint (alpha, x)-valuation >= 2. BTc and ABc·al and FNW, FNW2 carry a factor alpha; no constant-times-V term is allowed (f1 must stay 0).
4. Monotonicity: all quantities enter with nonnegative coefficients; NEG is unchanged.

[probe] The probe implementation of §5.1-5.2 passes `tail_cert` in all three modes: radial H (0, 2), L (0, 1); s H (2, 5), L (0, 1); joint, floors off, H (2, 5), L (2, 5).

## 6. Constants: what is needed, how to certify, what exists

All new sups: python-flint arb, 200 bits, box subdivision as in the script named, **upper** endpoint, then `math.nextafter` up when converted to a double; inputs with sha256 in an args block. No typed decimal may be below its raw value: v4 reads the raw JSON and rounds up (§10, stage C5).

| Constant | Definition | Interval / box | Method (script to copy) | Status |
|---|---|---|---|---|
| MU1, MU2 | Garding rates | [10.25, 60] | - | exists [cert] |
| E1, E0 | sup\|e1\|, sup\|e0\| | [10.25, 60] | `int_garding.py` output | exists [cert], unused |
| PSI1, PSI2, KAP_IR/Z1/OR, RWR, k4_j | §3, §4.4 | IR, supp | `g0_ramp_cert.py` | exists [cert G0] |
| k3_j, Z2i, KWIR | §4.4, §1.2 | IR | extend `g0_ramp_cert.py` | **new**, seconds |
| AF, BF, CF0-3 (and SB, SA, IA, DS0-3 for the E blocks) | raw in `coupling_constants.json` | - | round up from raw | exists [cert]; typed below raw: AF, CF0, CF2, CF3, DS1, DS2, SB, SA, IA |
| TQ d010, d110, d300, d100; FN[0], FN[1], C3; FNW[0] | typed 1.67e-3, 1.086, 3.314, 1.1e-4; 7.93e-7, 2.769e-6, 12.104; 1.09e-7 | - | round up from raw (1.67114e-3, 1.086422, 3.3140088, 1.129e-4; 7.93250e-7, 2.769321e-6, 12.104385; 1.092e-7) | typed below raw by 3e-6 .. 2.6e-2 [cert GI 1.3]. FN and C3 belong to the E blocks and the alpha-equation (the dominant forcing), not to E_I; they are fixed by the same stage C5 |
| MU1 | Garding rate, a **lower** bound | [10.25, 60] | read from `v2_T4_P2_0.25.json`, round down | provenance wrong [GI] |
| CBAR, d010, d100, d200, d110, d300, A0, A1 | coefficient sups of N | zones IR (box eta 1e-3), [16.25, 34.25], [34.25, 60] (box B-c, B-d); alpha-levels 2.5e-4 and 1e-3 | `pb_nlin.py` + `r60_t3.py` + `pb_ramp.py` | far: exists at 1e-3 only; IR: linear ones at 3 levels, nonlinear ones at 1e-3/1e-2 only. **Extend** |
| CA, CB | sup\|c_x + c_U0 U1 + c_U1 a2\|, sup\|c_U1\| | same zones, levels | new sympy derivative of the (0,0,1) coefficient, then as above | **new** |
| d^a, d^b for d010, d100, d200, d110, d300; B0, B1 | §4.5 | same | "t3d": sympy total derivative of `Fpert` with rules U0' = U1, a_j' = a_{j+1}, profile ODE; needs a5 from `c4fun` | **new** (main new certification, est 4-5 agent-h) |
| hA0, HA0, hB0, HB0, hE, HE, C0Hb | weighted integrals with Phi | Z1 | sum of len x enclosure, as `coup_common.py` does for Phi(34.25) | **new**, tiny |
| C0P | \|\|zeta c0'\|\|_W | [10.25, 57] | `int_v56.py` | exists as 9.9e-12 (whole support) |
| FNW, FNW2 | forcing norms | supp | - | exists [cert v3]; retype |
| T0, T2, TC0, TC2 | transport sups | supp / IR | `g0_ramp_cert.py` | exists [cert G0] |
| Q2B, DQ0, DQ1, EA0, EA1 | §4.9 | supp, box | extend `int_v56.py` | **new**, tiny |
| NP3, NH2, NAB, NP12 | W-norms | supp | `int_v56.py` with `c4fun` | exists for C^3; recompute for C^4 |
| CAB0, MAB, QAB, PHIAB, BETA0, BETA1, QB, BETAF | modulation lemma | - | - | exists [cert]; strip part to be restated as a linear form (U-3) |

## 7. Old code terms -> new spec terms

| # | fin12f4 code (`rhs_iv.py`) | Spec | Status |
|---|---|---|---|
| 1 | `:189` NEG, level 0 and 2 | §4.1 | unchanged |
| 2 | `:147` `AF W1W1 + BF BB` | §4.2 | **corrected** (AF rounding) |
| 3 | `:147` `d001 V2V2` | §4.3: CBAR V2^2 + (C1 + 2 PSI1 CBAR) V1V2 + (PSI1 C1 + PSI2 CBAR) V0V2 | **corrected** (N1) |
| 4 | `:147` `(d010 V2 + d100 V1 + A1 V0) Vq` | §4.5 level 0 | unchanged; **+ new** A0 term |
| 5 | `:148` `(d200 s V1 + d110 s V2 + d300 s^2 V1) Vq` | §4.5 level 0 | unchanged |
| 6 | `:149` `CABQ Amx al V0` (+ `:151` Amc in L) | §4.9 T10 | unchanged |
| 7 | `:149` `CAB3985 Am V0V0` | §4.9 T8 level 0 (T0) | unchanged (retype T0) |
| 8 | `:149` `CAB25 Am W1W1` | §4.9 T9 level 0 | unchanged |
| 9 | `:150` `FnW (V1 + 1.67 V0)` in L | §4.8 | unchanged (retype C^4 triple) |
| 10 | `:152-161` term (a), D4 entry | §4.4: V0·RQ | **replaced** (new decomposition) |
| 11 | `:152-161` term (a), W3, B, W1, A entries | §4.6 level 0 | **corrected** (table chord -> row eta 1e-3) |
| 12 | - | §4.6 level 0 nonlinear IR monomials | **new** |
| 13 | `:164-168` `c0pn S_all V0` | §4.7 | **corrected** (range of S_all; + C0Hb) |
| 14 | `:170-174` FLUX0 | §4.10 | unchanged as an input; K0^OR redefined |
| 15 | `:178` `(CF·(D4, W3, B, W1)) D2v` | §4.2 | **corrected** (rounding) |
| 16 | `:178` `d001 V4 V4` | §4.3: CBAR V4 DV + C1 V3 DV | **corrected** (N1) |
| 17 | `:179` `d010 V3 D2v` | §4.5 table row 1 | unchanged |
| 18 | `:179` `d110 supv1 V3 D2v` | §4.5 row d110 g' g''' | unchanged |
| 19 | `:179` `2 d200 supv1 V3 D2v` | §4.5: 2 d200 sZ V2 | **corrected** (index) |
| 20 | `:179` `d110 supv1 (3V3 + 1.31V2) D2v` | §4.5: d110 sG (3V3 + PSI1 V2) | **corrected** (global sup sG) |
| 21 | `:179` `d300 supv1^2 V2 D2v` | §4.5: 3 d300 sZ^2 V2 | **corrected** (factor 3) |
| 22 | `:179` else-branch `d110 supv2 V2` (GN = 0) | - | **deleted** (not valid beyond 34.25) |
| 23 | - | §4.5 rows with d^a, d^b, E1c, E0c | **new** |
| 24 | - | §4.4 RQ·DV; §4.6 RWR·R5·DV | **new** |
| 25 | `:180` `CAB3985 Am V2V2` | §4.9 T8 level 2: T2 (39.1842) + Q2B terms | **corrected** |
| 26 | `:180` `c62 Am V2 V3` | - | **deleted** (no derivation; was silently covering #25) |
| 27 | `:180` `CAB03 Amx al D2v` (+ `:182`) | §4.9 T10 level 2 | unchanged |
| 28 | `:181` `al^3 D2v` in L | §4.8: FNW2(al) DV | **corrected** |
| 29 | `:183-187` FLUX2 | §4.10 | unchanged as an input; K2^OR redefined |
| 30 | - | §4.1 e-terms | **new** |
| 31 | - | §4.7 level 2 | **new** |
| 32 | - | §4.9 T9 level 2 (TC2) | **new** |
| 33 | - | §4.9 T11 beta, psi1, psi2, both levels (N2) | **new** |
| 34 | - | §4.9 T12, both levels | **new** |
| 35 | `constants.py:22-25` typed SB, SA, IA, AF, BF, CF, DS | raw values rounded up | **corrected** (AF, CF0, CF2, CF3, DS1, DS2, SB, SA, IA were below raw) |
| 36 | `constants.py:39-48` chord (secant) table | row eta = 1e-3 + assert | **corrected** ([GI 0-7]: the secant happens to hold at eta = 1e-3 with 0.01% to spare for d010 and A1; it is still not a proof) |
| 38 | `rhs_iv.py:144` `am_phi * FnW(al)` | §4.9: PHIAB·Fn(al) | **corrected** [GI 0-14] |
| 39 | `constants.py:16, 28, 30, 31` typed C3, TQ d010/d110/d300/d100, FN[0], FN[1], FNW[0] | raw values rounded up; far TQ replaced by the zone x level tables of C2 | **corrected** [GI 1.3] |
| 40 | `constants.py:19` MU1 with the L-int citation | value from `v2_T4_P2_0.25.json`, rounded down | **corrected** (provenance) [GI 0-1] |
| 37 | `unc.py:19-28` `lk()` hand-written L majorant | generated from the same L list as `rhs_pn` | **corrected** (SF-3) |

## 8. Feasibility

Probe: a copy of the fin12f4 verifier in the session scratchpad with §5 implemented (file kept here as `probe/rhs_probe.py`, diff against `V/rhs_iv.py` in `probe/rhs_probe.diff`). Same environment as `run_v33.sh`. New constants are the [est] values of `PROBE_CONSTS` in that file; existing constants as typed (CF raised by 2e-3 relative to cover rounding, RQ and R5). One process, nice 19.

**Baseline check** [probe = cert here, no new constant]: the copy reproduces fin12f4 exactly: tail anchor U = -0.0026930418625762047, first cell U = 0.03034405538749633.

**Cells** (alpha <= a_hi, Lambda cell as in fin12f4):

| Variant | tail anchor rate_s (a* = 3.685e-5) | first march cell (alpha <= 1.6e-4) |
|---|---|---|
| fin12f4 | -2.69304e-3 | +3.034406e-2 |
| corrected, far coefficients at the 1e-3 box | **-2.68724e-3** | +3.035151e-2 |
| same, all new terms x3 | -2.62781e-3 | +3.048312e-2 |
| corrected, far coefficients [est] at the 2.5e-4 box | -2.68859e-3 | +3.034873e-2 |
| same, all new terms x3 | -2.68686e-3 | +3.037025e-2 |

**Full march at alpha_1 = 1.5e-4, sqrt(Lambda_1) = 5e-8** (corrected, 1e-3 box): both runs close in 10 steps with the tail lemma; `tail_cert` PASS. One process, nice 19; 368 s and 74 MB (x1), 1214 s and 150 MB (x3; it shared the CPU with another job for part of the time).

| quantity | fin12f4 | corrected (x1) | ratio | all new terms x3 | ratio |
|---|---|---|---|---|---|
| peak sqrt(Lambda) | 5.307458e-8 | 5.307569e-8 | 1.000021 | 5.310132e-8 | 1.00050 |
| tail rho | 2.693042e-3 | 2.687525e-3 | 0.99795 | 2.621352e-3 | 0.9734 |
| sqrt(L_F) | 9.871704e-10 | 9.881832e-10 | 1.00103 | 1.000578e-9 | 1.0136 |
| alpha_F | 6.866519e-7 | 6.873414e-7 | 1.00100 | 6.956293e-7 | 1.0131 |
| J0 const / slope | 6.027237e-12 / 3.454743e-18 | 6.030381e-12 / 3.457438e-18 | 1.00052 / 1.00078 | **3.592047e-11 / 2.964025e-17** | **5.96 / 8.58** |
| J4 const | 6.236181e-8 | 6.239434e-8 | 1.00052 | 3.716571e-7 | 5.96 |
| V1 delta_t | 0.2 | 0.2 | | **0.01** | |
| s_I1, s_I2 | 1.160338e-3, 1.545618e-2 | 1.160362e-3, 1.545650e-2 | 1.000021 | 1.160922e-3, 1.546397e-2 | 1.00050 |
| M1_ramp, S_ramp | 0.432773, 7.473592e-4 | 0.432782, 7.473749e-4 | 1.000021 | 0.432991, 7.477357e-4 | 1.00050 |
| sup\|(a,b)\|, sup\|beta\| (`mod_sups`, 6 digits up) | 2.19103e-7, 3.33811e-9 | 2.19106e-7, 3.33811e-9 | 1.000014, 1 | 2.19175e-7, 3.33811e-9 | 1.00033, 1 |

Caveats of the march: the J values come from `unc.Phi`, whose hand-written L majorant does not contain the probe's new L terms (the beta forcing), so the probe's J are slightly low (SF-3); sup|beta| is dominated at step 0 by terms that do not depend on the peak.

**The sensitive output is J, not rho.** In the x3 run rho loses 2.7% but the J constants jump x6, because the unconditional pass falls from delta_t = 0.2 to 0.01 on its grid (`closure.py:186-190`): the tail anchor condition "H - (1 - delta)(NEG - theta·NEGV) <= 0" has little room at delta = 0.2 (r1* = -0.00268 with theta = 0.05). J_m and the slopes are inputs of (A) (H-J-W). So the real feasibility question for the downstream chain is whether delta_t = 0.2 survives with certified constants. At x1 it does [probe], and [GI]'s run agrees (delta_t = 0.2). The margin is somewhere between x1 and x3 on the new terms; I did not locate it. Levers if it fails: L1-L3 of §9 (they shrink the new terms), or a finer delta grid (0.15, 0.1: J grows like 1/delta, so x1.3-x2 instead of x6).

**Independent evidence [GI 2.2, 2.3].** GI patched its own copy of the verifier (different term forms and coefficients, written before it saw the audit or this spec) and ran the full march at alpha_1 = 1.5e-4 with everything on: 10 steps, `tail_cert` PASS, **rho = 2.688223e-3 (-0.18%)**, peak 5.307553e-8 (x1.000018), J_0..J_4 x1.000465, J0 slope x1.0007, sqrt(L_F) x1.0009, s_I1, s_I2, M1_ramp, S_ramp x1.000018 [cert as bounds of its patched right-hand side; coefficients hand/est]. My probe gives rho -0.21% and the same peak to 5 digits. Two independent implementations agreeing to 3e-5 of rho is good evidence that the size estimate is right. It is not a certificate.

Neither run covers:
- alpha_1 = 1.25e-4 (the remark point);
- the outer ramp (K^OR is an input in both; the flux lemma v2 is open);
- the beta strip constant as an absolute term (GI left it out; my probe has it in state form with [est] coefficients). If option 2 of §4.9 had to be used, F0 x7.4 and J-slopes x8.4 [cert GI];
- FN[0], FN[1], C3 rounded outward (about 1e-6 in the rate [est GI]);
- the corrected forcing part of |(a, b)| (Fn instead of FnW at `rhs_iv.py:144`): both runs still have FnW there; 0.2% of |(a, b)| [GI];
- certified values of any new constant.

Reading:
- The corrected E_I inequality closes with rho essentially unchanged [est]. The tail needs rate_s <= -0.002 (`closure.py:127`); the margin to that threshold is x1.34 as before. If it were missed, the march simply continues to smaller alpha (L5 of §9).
- The maximiser stays A-dominated (shape w ~ 1e-3..3e-2 in the probe's `best`), which is why E_I-side terms hardly register.
- What is **not** covered by the probe: the true values of the new constants (CA, CB, d^a, d^b, B0, B1, k3, the Hardy integrals, Q2B, DQ, EA). The x3 runs bound the effect of a factor 3 error on all of them at once.
- The remark point alpha_1 = 1.25e-4 was not probed (same structure, smaller alpha).

## 9. Levers, if the certified constants turn out worse than estimated

In order of invasiveness. Downstream cost is relative to the unavoidable baseline of §10.4.

| Lever | Effect | Downstream cost |
|---|---|---|
| L1. Recertify the far coefficients at alpha-level 2.5e-4 (the march never exceeds 1.6e-4) and by zone | CBAR 2.27e-2 -> ~3.5e-3, d010 /10, d110 /2.6, d200 /4 [est] | none (one more certification table, minutes) |
| L2. Level-2 cross term by parts: 2 PSI1 CBAR V3V4 -> sup\|(c psi')' - psi'^2 c\| V3^2 | removes the largest N1 piece at level 2 | none (one more sup) |
| L3. Tighter far box on [16.25, 34.25] (\|w'\| <= 2e-3 instead of 0.02; certified sI1 <= 1.16e-3) | CBAR, CB on that zone /5 [est] | none, but a new asserted hypothesis with margin x1.7 |
| L4. Young weights | none needed: the branch-and-bound takes the sup over shapes, so no splitting constant is chosen by hand | - |
| L5. More march steps before the tail (automatic) | absorbs a rate loss up to ~8e-3 at alpha -> 0 [est from d rate/d alpha ~ 220] | peak sqrt(Lambda) up by ~e^{delta·5.5}; all consumers of the peak move |
| L6. theta, delta of the unconditional pass; finer delta_t grid (add 0.15, 0.1) | J constants only; protects against the x6 jump seen in the x3 probe | (A) chain inputs J_m |
| L7. mu (window [2.906e-3, 0.1135]) or eps_I | changes Lambda itself | **step-4 hand-off stages (H1) and both kits**; avoid |
| L8. S0V, S1V (inputs from (A)) or the hand-off level 5e-8 | sZ, sOR, boxes | (A) margins; step 4; avoid |

Recommendation: implement L1 in v4 from the start (it is part of the "t3d" stage anyway, since the derivative coefficients must be computed per level), keep L2, L3 in reserve, and do not touch L7, L8.

## 10. Coding plan for verifier v4

Rules: new directory `/home/claude/s1spec/step3nl_v4/`, new tag **fin13v4**; `step3nl_v33` and `step3nl` stay byte-identical. Every changed script is backed up next to it as `*_pre_v4.py` (byte copy of the v33 file, checked by `cmp`). One process, `nice -n 19`, < 800 MB. Every stage writes an args block (argv, env, input paths with sha256). SHA256SUMS at the end. Ask Micah before stage P (it re-pins everything downstream).

### 10.1 Stages
| Stage | Content | Copy from | Model, time | Compute |
|---|---|---|---|---|
| S0 setup | `cp -a step3nl_v33 step3nl_v4`; remove old run outputs from the copy's `verifier/`; make all read paths explicit (the scripts read `/home/claude/s1spec/step3nl/graph_struct.pkl`, `step3/hstail.py`, `mc_prof.txt`: read-only, pin by sha256) | - | Sonnet, 1 h | - |
| S1 identities | extend `ei_identities.py`: symbolic N_low' with the total-derivative rules, cross-checked at random float points against finite differences; the v-equation residual (T1..T12 sum minus zeta·g_tau) = 0 symbolically | this directory | Opus, 2-3 h | seconds |
| C1 cut_sups | k3, Z2i, KWIR and re-export of the G0 sups | `fix3/audit/g0_ramp_cert.py` | Sonnet, 1 h | seconds, 85 MB |
| C2 t3d | zone x level tables of CBAR, CA, CB, d, d^a, d^b, A0, A1, B0, B1 (IR box eta 1e-3; far boxes B-c, B-d); a5 from `c4fun` | `agent_int/planB/pb_nlin.py`, `r60_t3.py`, `pb_ramp.py`, `agent_int/v3/c4fun.py` | Opus, 5-6 h | est 5-15 min, < 500 MB (check: `planB/r60.log`, `ramp.log` give the v33 timings) |
| C3 mod_x | Q2B, DQ0, DQ1, EA0, EA1; NP3, NH2, NAB, NP12 for the C^4 functions | `agent_int/v2/int_v56.py`, `agent_int/v3/int_f23.py` | Opus, 2 h | 1-2 min |
| C4 hardy | hA0, HA0, hB0, HB0, hE, HE, C0Hb | `agent_coup/coup_common.py` (Phi tabulation) | Sonnet, 1.5 h | 1-2 min |
| C5 constants | `constants_v4.json` built from the raw JSONs by one script that rounds **up** and records (raw, typed, source sha256); `constants.py` reads it and contains no typed certified literal; a test asserts typed >= raw for every entry | `verifier/constants.py` | Sonnet, 2 h | seconds |
| R rhs | `rhs_iv.py` per §5 (E/E2 blocks untouched except the items of §11 that Micah approves); L terms as an explicit list consumed by both `rhs_pn` and `unc.lk`; `closure.py`: asserts of §1.3, strip state form in the E_I terms, args block; `mod_sups.py` unchanged unless U-3 is adopted there too | `verifier/rhs_iv.py`, `unc.py`, `closure.py`, `tail_cert.py` | Opus, 5-6 h | - |
| T tests | (i) regression: with env `NEW=0` and the old typed constants the v4 code reproduces the fin12f4 cell values bit for bit (two cells are quoted in §8); (ii) `tail_cert.py` PASS; (iii) `check_vs_float`-type point test of every new term against an independent float evaluation; (iv) `sanity.py`, `sanity_unc.py` | - | Opus/Sonnet, 2-3 h | 5 min |
| P production | `tail_cert`; `closure.py 1.5e-4`; `closure.py 1.25e-4`; `mod_sups.py`; then the same again into a second directory and `cmp` (byte-identical apart from the timing fields `secs`, `runtime_s`) | `run_v33.sh` -> `run_v4.sh` | Sonnet, 1 h | est 2 x (6-8 min per alpha point) = 25-35 min; RAM to be read from the probe (§8) |
| D record | `V4_VERIFIER.md`: fin12f4 vs fin13v4 table as in `V33_VERIFIER.md`; SHA256SUMS; list of consumers | `V33_VERIFIER.md` | Sonnet, 1 h | - |
| W write-up | `ei_lemma_v4.tex`: this spec as a lemma with proof (statement = §5 + §4.10, hypotheses = §1.3), plus the regularity lemma U-1 | this file | Opus, 3 h | - |
| V review | independent review from the scripts and the tex, not from reports | - | Fable, one pass | - |

Total: about 26-30 agent-hours plus the review (sum of the rows: 28.5). Compute under 1 h in total.

### 10.2 Certification method per constant
See §6, column "Method". Common rules: arb at 200 bits; box widths as in the copied script; every sup is the max of upper endpoints; integrals are sums of length x upper enclosure; conversion to double by `float(x.upper())` then `nextafter(+inf)`; each JSON entry stores the raw upper bound, and the verifier's interval constant is built from the raw value (`_c()` of `rhs_iv.py:13-16` must stop widening a typed decimal by +-1 ulp around a value that is below the certified one).

### 10.3 Checks
- byte-identical double run (two directories), args blocks, SHA256SUMS;
- `tail_cert` PASS on the production code path (it executes the verbatim `rhs_pn`);
- asserts of §1.3 active in the production runs, with the margins logged;
- f1 = 0 (`closure.py:195`);
- the regression test T(i).

### 10.4 Downstream, after v4 (each step needs Micah's go; none is part of v4)
1. **(A) chain t2 -> t4 -> zone3**, new tag (v36): consumes peak, J_m and slopes, abar, betabar, sI1, sI2, M1_ramp, S_ramp, rho, L_F, alpha_F from fin13v4. **Reapply the R3 patch set `fix4_r3/patches/` (r3_assembly_arb, r3_assembly_rep, r3_u34_sweeps, and the path patches)**; R3_REPORT says so explicitly. Add ROUTE's level-2 export r_2 = 1/1000 in the same rerun (units G1, F2').
2. **c_dU.** c_dU = 1.96498e-8 is only x1.018 below the typed bound 2e-8 (R3_REPORT). It scales with betabar + DL·alpha_1. From the probe, betabar does not move at 6 digits (3.33811e-9) and abar moves by 1.4e-5 relative (x1 run; 3.3e-4 in the x3 run) [probe]. If v4 raises betabar by more than 1.7%, raise the typed bound (it is used by `u1_ghat`'s M1z check and zone3's assertion) and recheck those two. Conversely, U-3 (strip part of beta as a state form in `mod_sups.py` too) would *lower* betabar, since STRIP_BETA = 1.06e-9 is a third of betabar = 3.34e-9 [read]; that is a lever for c_dU, not needed for E_I.
3. Flux lemma v2 (G2) and the L^1 consumer (G3) on the v36 interface, with K^OR of §4.10.
4. Step 4: new manifest (v34 file; the closed manifest and manifest_v33 stay), re-pin of the stages that pin `A_interface.json` / U_cert / zone3; new driver that reads `EXCEPTION_v32_*.json` (CLAUDE.md, Nov 15).
5. Kits: kit_step3 v4 (adds C1-C5, R, P, and later flux_cert_v2, l1_check), kit_step4b v4.
6. Texts: `step3_v33_main.tex` successor (§ verifier), numbers file, erratum entries for every row marked corrected/new/deleted in §7.

## 11. What is not yet justified, and side findings

Not yet justified (inside this spec):
- **U-1** Regularity lemma for the two energy identities (§3). Not written anywhere.
- **U-2** The whole "needed" side is [hand] except the identities marked [sym]. In particular the expansion of N_low' (§4.5) and the mean-value form of D(d000) must be machine-checked in stage S1 before any constant is computed.
- **U-3** The strip part of the modulation bounds as a linear form in strip sups (§4.9). I read the coefficients in `mod_ibp.json`; I did not check the lemma's derivation, nor that `linear_alpha_S0` carries a factor alpha.
- **U-4** Hypothesis B-d on [34.25, 60] (|w| <= 0.5) by |w(34.25)| + 25.75·S0V <= 0.395 [hand]; must be asserted in v4.
- **U-5** The O(k) terms T12 are bounded, not computed; all their constants are [est].
- **U-6** C^3 versus C^4 slow-manifold functions: the level-2 energy needs the C^4 ones (NOTES_v3 "KEY FINDING"); several typed norms are for the C^3 ones. The difference is ~1e-7 relative [read], but v4 should use one set.
- **U-7** The probe and [GI]'s rerun are experiments with estimated constants. They do not certify closure, and they share the gaps listed in §8.
- **U-8** [GI 4.4.2] If any degree-1, alpha-degree-0 term survives (the strip constant of |(a, b)| times NP12, ~1e-15), the tail lemma's floor formula has no slot for it with QUAD = 1. This spec avoids all such terms (state forms); v4's `tail_cert` must confirm L has joint valuation >= 2 with floors off and that the only floors are FLUX0, FLUX2.

Side findings in the same verifier, outside E_I (I did not audit the E and E2 blocks; these are things I ran into) [read unless marked]:
- **SF-1** (also noted in [GI] §3 item 2: "the stakes are higher in E0") The E0 block has no term for the strip part of beta: the beta-loop there is Pi3 (linear), MOD_ALPHA (alpha-linear), `qb752` (quadratic) and 5.72·Fn (forcing), and `CLIN` = 0 (`constants.py:66`, `rhs_iv.py:107-129`). As a constant, 7.52·STRIP_BETA·A would be a floor far above the peak [hand: f1 = 8e-9]; as a state form it is small (coefficient ~8.6e-4 on A·sqrt(V1V2), about 3e-4 of the available dissipation [est]). It should be added in state form.
- **SF-2** The strip hypothesis "S_i <= 0.02" of the modulation lemma includes S_3 = sup|w'''| on sigma in [29, 32] (`commutator_S0_S3`). S_0, S_1, S_2 are bounded by the energy (5.9e-4, 1.16e-3, 1.55e-2 [cert log]); S_3 is not bounded pointwise in tau by Lambda (it needs V3, i.e. V4). I found no check of it. Suggested fix: integrate the strip commutator pairing by parts so that only S_0 (or L^2 norms) appears; otherwise treat c·S_3 with supv2 inside the alpha-equation, which needs a time-integrated argument.
- **SF-3** `unc.lk()` (`unc.py:19-28`) majorises L by a hand-written list that omits `CABQ·Amc·al·V0` (`rhs_iv.py:151`) and `CAB03·Amc·al·D2v` (`:182`), so Phi(alpha) is not an upper bound of sup(L - delta N) as coded. Size tiny [est], but the J constants are certified numbers consumed by (A).
- **SF-4** The box hypotheses of §1.3 are logged but not asserted in `closure.py`.

## Files in this directory
- `EI_SPEC.md`: this file.
- `ei_identities.py`, `ei_identities.json`: exact sympy checks I1-I11 (all residuals 0).
- `probe/rhs_probe.py`: the probe's `rhs_iv.py` (feasibility experiment only; **not** verifier v4, constants [est]).
- `probe/rhs_probe.diff`: its diff against `V/rhs_iv.py`.
- `probe/`: logs and rate tables of the probe runs.
- `SHA256SUMS`.
