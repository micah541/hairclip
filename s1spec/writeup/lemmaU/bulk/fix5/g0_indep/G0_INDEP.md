# G0_INDEP: independent check of the step-3 tip-zone verifier's E_I inequality (tag fin12f4)

Reviewer: Claude (Fable), Nov 19, 2026. Read-only on the proof tree; everything here lives in `writeup/lemmaU/bulk/fix5/g0_indep/`.
Sections 0-3 were written to disk **before** `fix3/audit/` was opened; that version is kept as `G0_INDEP_sections0-3_pre_audit.md` (sha256 c03b0d8e...). After the audit was read, three things changed in this file: one retraction in row 0-15, sec. 2.3 (results of a run that was already under way), and sec. 4 (comparison).

Labels: **[cert]** = computed here in python-flint `arb` or by the verifier's own outward-rounded B&B; **[hand]** = derived by hand from the equation (one identity machine-checked, see below); **[est]** = estimate.

## 0. Verdict in three lines

- **The coded inequality is not a valid upper bound for (1/2) d/dtau (V0^2 + eps_I V2^2).** Eighteen terms of the v-equation are absent (11) or under-covered (7) in `verifier/rhs_iv.py`, 17 typed constants are rounded the wrong way, and no written derivation of the inequality exists (ASSEMBLY.md sec. 5-6 is a term list with [F]/[A] tags, and it disagrees with the code in at least 7 places).
- **Verdict: FIX-AND-RERUN.** Every gap I found is a missing estimate with room, not a failed one. With all of them switched on in a patched copy (coefficients [hand]/[est], B&B [cert]) the tail rate moves from -0.0026930 to -0.0026881 and the first march cell from +0.0271016 to +0.0271067. The worst shape is purely inner (A-dominant); E_I enters it only at second order.
- **Not fatal, but not cosmetic either.** The largest single omission is 24 % of the level-0 E_I dissipation with the boxes the code types (3 % with the boxes actually reached). The inner ramp at level 2 has never been derived at all and needs a reformulation, not just a constant.

## 0.1 What was checked, and how

Object: `step3nl_v33/verifier/rhs_iv.py` (242 lines), `constants.py`, `closure.py`, `unc.py`, `../inner_closure_float2.py`, and the scripts that produce the typed constants (`agent_coup/`, `agent_int/v2`, `v3`, `planB`, `c4switch`, `agent_lint/`, `agent_sob/`, `agent_mod/ibp`).
Production environment (`run_v33.sh`): `ADD=abcd FIX=1 MOD=ibp QUAD=1 CMODE=abs RAMP=table GN=1`, alpha1 = 1.5e-4, sqrt(Lambda1) = 5e-8.

The equation (from `graph_struct.py`, `agent_int/planB/pb_nlin.py:16-24`, BOOTSTRAP.md "Unknowns and equations"):

    w_tau = L w + D(w) + Phi_F - beta psi3 - 2 alpha beta h2 - a psi1 - b psi2 + a A1(u) + b B1(u),     u = u0 + w, u0 = alpha psi3 + alpha^2 h2
    D = F(u0+w) - F(u0) - L w = c w'''' + e w''' + f w'' + d200 w''^2 + d110 w''w''' + d300 w''^3 + R0,   c = 1 - J^-4,  |R0| <= A0|w| + A1|w'|
    A1(u) = -u + u' q/alpha_g,   B1(u) = u' cos(theta)/alpha_g,   psi1 = p,  psi2 = -sin(theta)

With g = w', v = zeta g, L1 = -d^4 - k^2 d^2 - (7km + q/4) d - (8m^2 + 5kp/2):

    v_tau = L1 v - [L1,zeta] g + zeta c0' w + zeta D' + zeta Phi_F' + zeta M',     M = -beta psi3 - 2 alpha beta h2 - a psi1 - b psi2 + a A1(u) + b B1(u)

Level 0 pairs this with W v. Level 2 pairs its second derivative with W v'', which the code does as <., D2v>_W with D2v = (W v'')''/W, ||D2v|| <= V4 + 2.61 V3 + 1.7 V2 (`rhs_iv.py:176`).

I did **not** audit: the inner E0/E2 inequality (except where its constants sit in the same lines), the modulation lemma, the outer-ramp flux lemma (term (c), known open as F1), Lemma S inputs S0V/S1V, the B&B and tail machinery. Line numbers of `constants.py` and `rhs_iv.py` refer to the proof-tree files.

Machine check: the level-0 top-order identity used in item 0-3 below is verified symbolically (`check_identity.py`, Euler operator = 0).

## 1. Term inventory

Status: **OK** covered; **SLACK** not coded but dominated by a coded term; **LOW** coded with too small a constant or the wrong variable; **MISSING** absent; **EXT** delegated to another lemma.

### 1.1 Level 0: (1/2) d/dtau V0^2

| # | term | needed bound [hand] | code | status |
|---|---|---|---|---|
| 0-1 | <L1 v, v>_W | -mu1 V0^2 - (1/2) V2^2 | `rhs_iv.py:189` (NEG: mu1t, halft), `constants.py:19` | OK. mu1 = 0.1552 is above the L-int certificate it cites (0.155167, `agent_lint/runs/cert_PX2_XA10.25_A230.25.json`) but below the independent 0.15552 of `agent_int/v2/v2_T4_P2_0.25.json`. Valid, wrong provenance. |
| 0-2 | linear commutator, inner ramp: -<[L1,zeta]g, v>_W = int A_z g^2 + B_z g'^2 | a_f W1^2 + b_f B^2 | `rhs_iv.py:147` (AF, BF) | OK in form (identity verified in `agent_coup/verify_identity.py`; drift term is inside A_z, `coup_common.py:78-79`). AF typed 6.9555e-2 < certified 6.95550376e-2 (LOW by 5e-7). |
| 0-3 | quasilinear top order where zeta = 1: int W v (c v''')' = int W c v''^2 + int W (c' - 2 psi' c) v'v'' - int (W psi' c)' v v'' | d001 V2^2 + sup\|c' - 2psi'c\| V1V2 + sup\|psi''c + psi'c' - psi'^2 c\| V0V2 | `rhs_iv.py:147`: `d001*V2V2` only | **MISSING** (both cross terms). ASSEMBLY.md:78 has "0.025 V1V2"; the code has nothing. |
| 0-4 | quasilinear lower order where zeta = 1: -int W (D - c w'''') (v' - psi' v) | (d010 V2 + d100 V1 + A1 V0 + A0\|\|w\|\|_W + d200 s V1 + d110 s V2 + d300 s^2 V1)(V1 + 1.305 V0), s = sup\|w''\| | `rhs_iv.py:147-148` | OK except A0\|\|w\|\|_W, which is absent: SLACK (A0 <= 9.3e-9 beyond 16.25; the typed A1 = 7.1e-5 is the ramp value and exceeds the far value 1.3e-7 by 500x). d010, d110, d300 typed below their certified sups (LOW, 7e-4, 4e-4, 3e-6). |
| 0-5 | quasilinear, inner ramp, linear deviation: -int W D (zeta v' + (zeta' - psi' zeta) v) with D in rho-norms | 0.66 (c_r D4 + e_r W3 + f_r B + A1_r W1 + A0_r A)(V1 + 0.76 V0) if no second integration by parts is done on the ramp | `rhs_iv.py:152-161` (term (a), table) | OK **only if** item 0-3 is restricted to zeta = 1, which produces boundary terms at x = 16.25. Done globally (the only clean way), 0-3 leaves the commutator c (3 zeta' w''' + 3 zeta'' w'' + zeta''' w')(v' - psi' v) + c zeta' w'''' v on the ramp: the W3, B, W1 coefficients must be (1.10 c_r + e_r, 0.63 c_r + f_r, 0.24 c_r + A1_r), i.e. 3.0x, 2.7x, 2.1x the coded ones. **LOW.** |
| 0-6 | quasilinear, inner ramp, monomials d200 w''^2, d110 w''w''', d300 w''^3 | 0.66 (d200_r S2 B + d110_r S2 W3 + d300_r S2^2 B)(V1 + 0.76 V0), S2 = sup_ramp\|w''\| (rho-Sobolev) | none. `planB/pb_ramp.py:11` tabulates only the five linear coefficients | **MISSING.** |
| 0-7 | S-dependence of the ramp table | coefficient at the certified box | `constants.py:48`: secant between eta = 0 and eta = 1e-2 | LOW in principle: a secant bounds a sup only if the sup is convex in eta. Checked against the eta = 1e-3 rows: holds with 0.01 % to spare for d010 and A1. Use the eta = 1e-3 row (S_ramp <= 7.5e-4). |
| 0-8 | zeta c0' w | \|\|zeta c0'\|\|_W sup\|w\| V0 | `rhs_iv.py:164-168` | OK on [10.25, 34.25]. On (34.25, 60] `S_all` is not a bound for sup\|w\|, but sup\|c0'\| <= 5.6e-35 there (`int_constants_v2.json`): SLACK. |
| 0-9 | forcing zeta Phi_F' | \|\|Phi_F\|\|_W (V1 + 1.67 V0) | `rhs_iv.py:150`, `constants.py:31` | OK as a sum, wrong term by term. FNW = (1.09e-7, 6.5e-9, 0.1823) are the C^3-profile values; the decomposition uses the C^4 profiles (FN in line 30 is from `c4switch`), for which the certified values are (1.243e-8, 4.03e-8, 0.1823) (`v3/int_constants_v3.json`). The alpha^2 coefficient is 6.2x too small, the alpha coefficient 8.8x too large; the sum is valid for alpha <= 1. Also 1.09e-7 < certified 1.092e-7. |
| 0-10 | -beta zeta psi3' - 2 alpha beta zeta h2' | (0.2596 + 3.592\|alpha\|) \|beta\| V0, \|beta\| <= (0.00244 + 66.53\|alpha\|) A + 5.53 Q + 0.657 Fn + 1.06e-9 | none | **MISSING.** ASSEMBLY.md:80 lists it; neither `rhs_iv.py` nor `inner_closure_float2.py:91-94` has it. |
| 0-11 | -a zeta psi1' - b zeta psi2' | 1.88e-7 \|(a,b)\| V0 | none | **MISSING** (ASSEMBLY.md:80 lists it). |
| 0-12 | a zeta A1(u0)' + b zeta B1(u0)' | (0.2731 + 0.0159)\|alpha\|\|(a,b)\| V0 | `rhs_iv.py:149,151` (CABQ) | OK. |
| 0-13 | a zeta (q w'/alpha_g - w)' + b zeta (w' cos theta/alpha_g)' | \|a\| (1/2) sup\|1 - q psi'\| V0^2 + 0.653\|b\| V0^2 + (0.436\|a\| + 0.16\|b\|) W1^2 | `rhs_iv.py:149` (CAB3985 = 38.5 + 0.66, CAB25 = 2.5) | OK (certified sup\|1 - q psi'\| = 74.37). |
| 0-14 | forcing part of \|(a,b)\| | 3.8404 \|\|chi Phi_F\|\|_rho = 3.8404 Fn(alpha) | `rhs_iv.py:144`: `am_phi * FnW(al)` | **LOW**: uses the W-norm residual FnW instead of the rho-norm Fn (7.3x too small at small alpha, 66x in the alpha^3 term). `mod_sups.py:25` uses Fn, so the code contradicts itself. |
| 0-15 | everything on the outer ramp [57, 60] (linear commutator, D with zeta < 1, a q zeta zeta' g^2, forcing with sup zeta' = 0.73) | absolute constant | `rhs_iv.py:170-172` (FLUX0) | EXT (flux lemma, open as F1). *Corrected after reading AUDIT_RAMP (sec. 4):* my pre-audit text said 1.67 fails on [57,60] because sup\|zeta'\| + sup psi' is about 2.0. The quantity needed is sup(\|zeta'\| + psi' zeta) <= 1.5833 [cert, `check_kappa_outer.py`], so 1.67 is valid there and the forcing needs no outer-ramp exception. |

### 1.2 Level 2: (eps_I/2) d/dtau V2^2

| # | term | needed bound [hand] | code | status |
|---|---|---|---|---|
| 2-1 | <(L1 v)'', v''>_W | -(1/2)V4^2 - mu2 V2^2 + e1 V1V2 + e0 V0V2 | `rhs_iv.py:189` | Dissipation OK (mu2 >= 0.65552). The e-terms (5.67e-6, 5.11e-6) are **MISSING** (ASSEMBLY.md:86 has them). |
| 2-2 | linear commutator, inner ramp | (c3 D4 + c2 W3 + c1 B + c0 W1) \|\|D2v\|\| | `rhs_iv.py:178`, `constants.py:24` | OK in form. CF typed (0.4421, 0.5259, 0.4235, 0.4250) vs certified (0.442143, 0.525899, 0.423528, 0.425017): three of four rounded **down** (9.7e-5, 6.5e-5, 4.1e-5). |
| 2-3 | top order, zeta = 1: c v'''' paired with D2v | d001 V4 (V4 + 2.61 V3 + 1.7 V2) | `rhs_iv.py:178`: `d001*(V4*V4)` | **MISSING** cross terms d001 V4 (2.61 V3 + 1.7 V2). |
| 2-4 | (c' + e) v''' | (sup\|c_U1 u0''\| + sup\|c_U1\| s + d010) V3 \|\|D2v\|\|, c_U1 = 4 J^-6 u' | `rhs_iv.py:179`: `d010*V3` | **LOW**: c' absent (linear part 6.5e-4, nonlinear part 0.42 s, against the coded d110 s = 1.09 s). |
| 2-5 | (e' + f) v'' + f' v' + R0' | (d100 + \|e'\|) V2 + \|f'\| V1 + derivative-level A-terms | none | **MISSING.** No certified bounds for e', f', R0' exist anywhere in the tree. |
| 2-6 | 2 d200 w''w''' | 2 d200 s V2 | `rhs_iv.py:179`: `two*d200*supv1*V3` | **LOW**: wrong variable (V3 is not >= V2). |
| 2-7 | d110 (w'''^2 + w''w'''') | d110 s (3V3 + 1.31V2) + d110 s V3 | `rhs_iv.py:179` | OK. |
| 2-8 | 3 d300 w''^2 w''' | 3 d300 s^2 V2 | `rhs_iv.py:179`: `d300*supv1^2*V2` | **LOW**: factor 3. |
| 2-9 | d200' w''^2 + d110' w''w''' + d300' w''^3 | sups of coefficient derivatives | none | **MISSING.** |
| 2-10 | inner ramp, zeta D' with zeta < 1 | needs zeta w^(5) = v'''' - 4zeta'w'''' - 6zeta''w''' - 4zeta'''w'' - zeta''''w'; then rho-norm terms 0.66 [(1.46 c_r + e_r) D4 + ...] \|\|D2v\|\| and monomials | none (the far-field V-terms are applied as if zeta = 1) | **MISSING, and not derivable as the code is organised.** A direct rho-norm bound of zeta D' needs w^(5), which no energy controls; bounding zeta w''w'''' by sup\|w''\| D4 \|\|D2v\|\| grows like r^2.25 in the radial variable (D4, V4) and breaks the radial lemma (`tail_cert.py`, mode radial needs degree <= 2). It must be written with one factor converted to v (sup\|w''\| (V3 + rho-terms) \|\|D2v\|\|, degree 1.75). |
| 2-11 | forcing (zeta Phi_F')'' | (2.701e-8\|alpha\| + 1.156e-7 alpha^2 + 2.603e-3\|alpha\|^3) \|\|D2v\|\| (`v3/int_constants_v3.json`, F1a/F1b) | `rhs_iv.py:181`: `al^3 * D2v` | **LOW** for alpha < 1.64e-4, i.e. on the whole run: by 1.20x at alpha = 1.5e-4, 19.9x at the tail anchor, 5.7e4x at alpha_F [cert]. ASSEMBLY.md:89 tags it "[F] not yet certified"; it was certified in v3 and never wired in. |
| 2-12 | beta, psi1, psi2 | (0.2596 + 3.592\|alpha\|)\|beta\| \|\|D2v\|\| + 1.88e-7\|(a,b)\| \|\|D2v\|\| | none | **MISSING.** |
| 2-13 | a, b transport | \|a\| (1/2) sup\|q psi' + 3\| V2^2 = 39.19\|a\| V2^2, plus 0.653\|b\| V2^2 | `rhs_iv.py:180`: (38.5 + 0.66)\|(a,b)\| V2^2 + 8.19\|(a,b)\| V2V3 | OK through the Euclidean bound and the V2V3 slack; the constant 38.5 + 0.66 = 39.16 is the level-0 one and is below 39.19. ASSEMBLY.md:90 tags these "[F constants]". |
| 2-14 | zeta c0' w at level 2 | 9.9e-12 sup\|w\| \|\|D2v\|\| | none | MISSING (negligible). |
| 2-15 | outer ramp | absolute constant | `rhs_iv.py:183-185` (FLUX2) | EXT. GN's sup\|v'\| is also wrong there (v' = zeta w'' + zeta' w', up to 0.030 against the coded 0.01875). |

### 1.3 Typed constants: direction of rounding

All of these are meant as upper bounds. Ratio typed/certified - 1 [cert, `g0_numbers.json`]:

| constant (`constants.py` line) | typed | certified source value | rel. error |
|---|---|---|---|
| AF (23) | 6.9555e-2 | 6.955503763e-2 | -5.4e-7 |
| CF[0], CF[2], CF[3] (24) | 0.4421, 0.4235, 0.4250 | 0.4421427, 0.4235277, 0.4250172 | -9.7e-5, -6.5e-5, -4.1e-5 |
| DS[1], DS[2] (25) | 1.5418e-2, 2.3648e-2 | 1.54183748e-2, 2.36481789e-2 | -2.4e-5, -7.6e-6 |
| SB, SA, IA (22) | 3.5955e-5, 1.3372e-4, 9.5997e-12 | 3.59551663e-5, 1.33720281e-4, 9.59970542e-12 | -4.6e-6, -2.1e-6, -5.7e-7 |
| TQ d010, d110, d300 (28) | 1.67e-3, 1.086, 3.314 | 1.67114e-3, 1.086422, 3.3140088 | -6.8e-4, -3.9e-4, -2.7e-6 |
| TQ d100 (28) | 1.1e-4 | 1.129e-4 (ramp region only) | -2.6e-2 (harmless: far value is 2.5e-5) |
| FN[0], FN[1], C3 (30, 16) | 7.93e-7, 2.769e-6, 12.104 | 7.93250e-7, 2.769321e-6, 12.104385 | -3.2e-4, -1.2e-4, -3.2e-5 |
| FNW[0] (31) | 1.09e-7 | 1.092e-7 | -1.8e-3 |

Rounded correctly: BF, CF[1], DS[0], DS[3], PHI34, d001, d200, A0, QPSI, BETA0, BETA1, BETAF, QB, CAB0, MAB, QAB, PHIAB, MU2, the Sobolev constants (`sob_constants_P2.py` rounds outward).

FN is not an E_I constant, but it is the dominant forcing of the whole inequality (5.72 Fn A, and beta in the alpha-equation), and it is too small by 3.2e-4.

The TQ sups were computed for the box alpha <= 1e-3, |w'| <= 0.02; the run has alpha <= 1.5e-4. The true sups on the run are therefore smaller than the typed numbers, so none of these down-roundings can make the inequality false in fact. They do make the typed numbers not certified.

### 1.4 Where ASSEMBLY.md and the code disagree

| ASSEMBLY.md | code |
|---|---|
| :78 "0.025 V1V2" (top order) | absent |
| :80 \|beta\| 0.26 V0 + 2\|alpha\|\|beta\| 1.8 V0 | absent |
| :80 \|a\| 1.7e-7 V0 + \|b\| 1.7e-8 V0 | absent |
| :82 inner-ramp couplings "[A, small]" | linear part only |
| :86 5.7e-6 (V3V2 + V2^2) | absent |
| :89 \|\|zeta Phi_F'\|\|_W \|\|D2v\|\| "[F]" | alpha^3 \|\|D2v\|\| with coefficient 1.0 |
| :90 2\|a\| V2V3 | 62/CAB \|(a,b)\| V2V3 (larger; harmless) |

## 2. Sizes of the uncovered terms

### 2.1 By hand, arithmetic in arb (`g0_numbers.py`, `g0_numbers.json`)

| item | size | relative to |
|---|---|---|
| 0-3 top-order cross terms, typed boxes (d001 = 0.0227, sup\|w''\| = 0.01875): 0.0678 V1V2 + 0.0497 V0V2 [hand coefficients] | sup over shapes of ratio = **0.239** [cert] | level-0 dissipation mu1 V0^2 + V2^2/2 |
| same, boxes actually reached (alpha <= 1.5e-4, \|w'\| <= 0.0153; d001 about 2.2e-3 [est]) | 0.030 [cert arithmetic, est input] | same |
| 2-3 level-2 cross terms, d001 = 0.0227 | **0.110** [cert] | V4^2/2 + mu2 V2^2 |
| same, d001 about 2.2e-3 [est] | 0.011 | same |
| 0-10 beta coupling: cross term kappa a_ v_ in rate units (A = s a_, V0 = s v_/sqrt(mu)) | kappa = 9.14e-4 at alpha = 1.5e-4, 3.59e-4 at the tail anchor, 1.79e-4 at alpha = 0 [cert] | kappa/2 = **17 %, 6.7 %, 3.3 % of rho = 2.693e-3** |
| same, effect on the top eigenvalue (the V0 direction is damped at 0.155, the A direction at about 0.015) | kappa^2/(4 * 0.140) = 1.5e-6, 2.3e-7, 5.7e-8 [cert] | 0.06 % of rho |
| 0-10 degree-0 part, forcing: 0.2596 * 0.657 Fn(alpha) V0 | 2.73e-11 at alpha = 1.5e-4 [cert] | coded I0 forcing 1.67 FnW = 2.84e-11: the missing piece **doubles** the I0 forcing |
| 0-10 degree-0 part, strip constant of the code's own beta bound (STRIP_BETA = 1.06e-9, `constants.py:76`) | 2.75e-10 V0, not decaying in alpha. Young against mu mu1 V0^2: 9.76e-21 per unit tau [cert] | **7.4x the flux floor F0 = 1.31e-21**; sqrt(L_F) 9.87e-10 -> 2.87e-9, J-slopes x8.4 |
| 0-5 ramp commutator of the top-order term | coded W3, B, W1 ramp coefficients low by 3.0x, 2.7x, 2.1x [hand] | VERIFIER.md states 5-10x room at alpha = 1e-3 |
| 0-6 ramp monomials | d110_r S2 = 0.30 * 5e-3 = 1.6e-3 on W3 [est] against the coded 2.6e-5 to 3.2e-5 | **about 50x the coded W3 ramp coefficient** |
| 2-10 level-2 ramp, linear part | about 6.6e-5 D4 \|\|D2v\|\| [est] | 1.5e-4 of CF[0] = 0.4421, same order as the down-rounding of CF |
| 2-11 level-2 forcing | needed/coded = 1.20, 19.9, 5.7e4 at alpha = 1.5e-4, 3.68e-5, 6.9e-7 [cert]. Absolute size after Young: 2.1e-29 per unit tau at alpha = 1.5e-4 [cert] | rho * Lambda = 7e-18: irrelevant by 11 orders |
| 0-11, 2-12 psi1', psi2' | 1.88e-7 \|(a,b)\| V0 | 1.7e4 times smaller than 0-10 |
| 0-14 forcing part of \|(a,b)\| | 3.84 (Fn - FnW) = 3.9e-10 at alpha = 1.5e-4 [hand] | 0.2 % of \|(a,b)\| = 2.2e-7 |
| 2-1 e-terms | 5.7e-6 | 1e-5 of mu2 |
| 1.3 down-rounded constants | 5e-7 to 1.8e-3 | the constants themselves |

### 2.2 By rerunning the verifier's own B&B with the terms added

A copy of the verifier (`vp/`, originals hashed in `vp_original_sha256.txt`, patch in `apply_patch.py` and `apply_patch_constants.py`) has the missing terms behind an environment switch `G0FIX`:

- `x` item 0-3 with the typed-box coefficients (0.0600 + 0.42 s) V1V2 + (0.0394 + 0.549 s) V0V2;
- `b` items 0-10 and 2-12 (degree >= 1 part in H, the 0.657 Fn part in L; the strip constant is not in the run, see 2.1);
- `p` items 0-11, 2-12 (the part multiplying the strip constant of \|(a,b)\| is a degree-1 floor term and is switched off with the floors; it is 1e-31 per unit tau after Young [est]);
- `m` item 0-6 with the alpha <= 1e-3 ramp sups (d200, d110, d300) = (0.01331, 0.2968, 3.024);
- `q` items 2-1 (e), 2-3, 2-4, 2-5 (as 2 d100 V2 + 2 A1 V1 + 2 A1 V0, [est]), 2-6, 2-8;
- `f` item 2-11;
- `s` item 2-10 as an [est] model: 6x the ramp table times \|\|D2v\|\|, plus the monomials in the degree-1.75 form;
- `c` 13 of the 17 constants of sec. 1.3 rounded outward (not FN[0], FN[1], C3, which would add about 1e-6 to the rate [est], and not d100);
- `RAMPMULT=4` item 0-5.

Certified B&B bound U of the rate (tolerance 5e-4 above the point lower bound lb), same environment as fin12f4 (`cells.jsonl`):

| G0FIX | tail anchor (alpha <= 3.6849e-5, Ly in [2.8065e-15, 1.01x], floors off): U / lb | first march cell (alpha <= 1.5e-4, Ly about 2.5e-15): U / lb |
|---|---|---|
| none (reproduces fin12f4) | -0.0026930 / -0.0031889 | +0.0271016 / +0.0266021 |
| x | -0.0026891 / -0.0031889 | +0.0271020 / +0.0266021 |
| b | -0.0026891 / -0.0031884 | +0.0271054 / +0.0266055 |
| p | -0.0026930 / -0.0031889 | +0.0271016 / +0.0266021 |
| m | -0.0026930 / -0.0031889 | +0.0271018 / +0.0266021 |
| q | -0.0026898 / -0.0031889 | +0.0271020 / +0.0266021 |
| f | -0.0026930 / -0.0031889 | +0.0271016 / +0.0266021 |
| s | -0.0026929 / -0.0031889 | +0.0271001 / +0.0266021 |
| c | -0.0026928 / -0.0031889 | +0.0271020 / +0.0266021 |
| all | -0.0026886 / -0.0031883 | +0.0271055 / +0.0266066 |
| none, RAMPMULT=4 | -0.0026888 / -0.0031888 | not run |
| all, RAMPMULT=4 | **-0.0026881 / -0.0031876** | **+0.0271067 / +0.0266068** |

U carries the B&B tolerance, so differences between rows below about 5e-6 in U are search noise; lb (a point value of the patched right-hand side at its own maximiser) is the cleaner indicator and moves by at most 1.3e-6 at the anchor and 4.7e-6 at the first cell.

So with every item of sec. 1 added, the certified rate changes by 4.9e-6 at the tail anchor (0.18 % of rho) and by 5.1e-6 at the first cell. The baseline row reproduces the fin12f4 log to all printed digits.

These are [cert] as B&B bounds of the patched right-hand side. The patched right-hand side itself is **not** a certified inequality: its added coefficients are [hand] or [est] (items 2-5, 2-9, 2-10 in particular have no certified inputs).

### 2.3 Full march + tail + unconditional pass with everything on

`run_full.sh`: `G0FIX=bcfmpqsx RAMPMULT=4`, otherwise the fin12f4 environment, alpha1 = 1.5e-4, sqrt(Lambda1) = 5e-8. 335 s, one process, nice 19, 55 MB.
`tail_cert.py` on the patched right-hand side: **PASS** (radial, s, joint), `g0all_tail.out`. The march closes in 10 steps and the tail lemma applies (`g0all_1.5e-4.out`, `vp/rate_table_g0all_a0.00015.json`, comparison in `g0all_vs_fin12f4.json`):

| quantity | fin12f4 | all terms added | ratio |
|---|---|---|---|
| steps | 10 | 10 | |
| peak sqrt(Lambda) | 5.307458e-8 | 5.307553e-8 | 1.000018 |
| tail rho | 2.693042e-3 | 2.688223e-3 | 0.998211 |
| tail anchor a*, L* | 3.684882e-5, 2.806455e-15 | 3.684882e-5, 2.806563e-15 | 1, 1.000038 |
| floor sqrt(L_F) | 9.871704e-10 | 9.880549e-10 | 1.000896 |
| alpha_F | 6.866519e-7 | 6.872539e-7 | 1.000877 |
| J0, J1, J2, J4 (constants) | 6.027237e-12, 8.467211e-12, 1.870781e-12, 6.236181e-8 | 6.030042e-12, 8.471151e-12, 1.871651e-12, 6.239083e-8 | 1.000465 each |
| J0 slope | 3.454743e-18 | 3.457146e-18 | 1.000696 |
| s_I1, s_I2, M1_ramp, S_ramp | 1.160338e-3, 1.545618e-2, 0.432773, 7.473592e-4 | 1.160358e-3, 1.545645e-2, 0.4327807, 7.473725e-4 | 1.000018 each |
| V1 delta_t | 0.2 | 0.2 | |

Reading: [cert] as statements about the patched right-hand side; **not** a certificate of step 3, because the added coefficients are [hand]/[est] and four things are not in the run:

- the strip constant of beta as a floor term in E_I (sec. 2.1: floor up to x7.4 if kept absolute);
- FN[0], FN[1], C3 outward (about 1e-6 in the rate [est]);
- the two items of AUDIT_RAMP I had missed (sec. 4.3: 1.271 |(a,b)| W1 ||D2v||, and the 39.16 -> 39.19 constant), both below 1e-6 of the terms they sit next to;
- the outer ramp (flux lemma).

The alpha1 = 1.25e-4 point was not run.

## 3. Verdict

**FIX-AND-RERUN.**

- Not PASS: the coded right-hand side omits terms that are present in the equation (sec. 1), so `rate_table_fin12f4_*.json` does not certify the E_I inequality, and through it neither Lambda <= barrier nor the J_0..J_4 integrals that feed Lemma S(A) / Lemma U.
- Not KILL in the sense of "plausibly fatal": every omitted term can be bounded with quantities the scheme already controls, none changes a sign or a structural condition once written in the right form, and the certified rate with all of them added is the same to 0.2 %. The reason is structural: the binding direction of the Lyapunov inequality is the inner A-mode (damping 0.0152 against 0.155 for V0), and E_I errors enter it only through cross terms at second order.
- What "fix" means is more than retyping constants. A written derivation of the level-0 and level-2 inequalities has to exist first; the level-2 inner-ramp block needs a new formulation and new certified inputs (item 2-10, and derivative-level coefficient sups for 2-5, 2-9).

The items that matter most:

1. **Top-order cross terms of 1 - J^-4 after integration by parts** (0-3, 2-3). Largest omission relative to the E_I dissipation: 24 % at level 0 and 11 % at level 2 with the boxes the code types. Fix: add them, and recompute d001 and its derivative for the alpha-box of the run (expected about 3 % and 1 %).
2. **beta-coupling absent from E_I** (0-10, 2-12), with psi1', psi2' (0-11). Cross term up to 17 % of rho in rate units, effect 0.06 % of rho. Its degree-0 part doubles the I0 forcing, and the code's own strip constant in the beta bound would raise the floor F0 by 7.4x unless the strip part of beta is made state-dependent. This decision also affects E0 (beta's strip part is absent there too; out of scope here, but it is the same line of reasoning and the stakes are higher in E0).
3. **Inner ramp** (0-5, 0-6, 0-7, 2-10). Level 0: commutator of the quasilinear top-order term under-covered by 2-3x, quadratic/cubic monomials absent (50x the coded W3 coefficient), secant extrapolation of the table. Level 2: the whole nonlinear block is absent and cannot be bounded in rho-norms as it stands.
4. **Level-2 far block is wrong in five places** (2-4, 2-5, 2-6, 2-8, 2-9, plus the e-terms of 2-1) and **the level-2 forcing is coded as alpha^3** although the certified alpha- and alpha^2-norms exist (2-11).
5. **Seventeen typed constants are rounded down** (sec. 1.3), FNW refers to the C^3 profiles, and the forcing part of |(a,b)| uses the wrong norm (0-14). Each is harmless in size; together they mean no number in `constants.py` can be quoted as certified without rechecking it against its source json.

Downstream: anything pinned on fin12f4 (rho, peak sqrt(Lambda), J_0..J_4 and their slopes, s_I1, s_I2, M1_ramp, sup|(a,b)|, sup|beta|) should be expected to move by well under 1 %, except the floor and the J-slopes, which move by up to 8.4x if the strip constant is kept as an absolute term in beta. None of that is certified until the corrected verifier exists.

## 4. Comparison with AUDIT_RAMP (`fix3/audit/AUDIT_RAMP.md`, Opus, Nov 16)

Read after sec. 0-3 were on disk. The two reviews agree on the substance: the coded inequality is not a bound, and the fix is a written derivation plus a rerun.

### 4.1 Items I confirm (found independently, same file:line)

| AUDIT_RAMP | here | remark |
|---|---|---|
| (N1) top-order cross terms, levels 0 and 2 | 0-3, 2-3 | Same coefficients (0.059 V1V2 + 0.04 V0V2; d001 V4 (2.61 V3 + 1.7 V2)). I add the c' part (0.42 s) and a certified ratio: 24 % / 11 % of the dissipation, where the audit estimated 10-20 %. Identity now machine-checked (`check_identity.py`). |
| (N2) beta, psi1, psi2 couplings absent | 0-10, 0-11, 2-12 | Same norms. See 4.2 for the size. |
| level-0 monomials missing on the inner ramp | 0-6 | Same sups (0.0133, 0.297, 3.02). |
| level-2 inner ramp: quasilinear commutator, lower-order and quadratic terms | 2-10 | Audit certifies the commutator at 7.74e-5 of CF; my [est] was 1.5e-4 including the e_r part. Same order. |
| CF0, CF2, CF3, AF typed below certified; SB, SA, IA, DS1, DS2 too | 1.3 | Identical numbers. |
| level-2 forcing coded as alpha^3 | 2-11 | Identical. |
| Dd001 v''' missing; 2 d200 with V3 instead of V2; d300 without the factor 3; d100, A-terms missing at level 2 | 2-4, 2-5, 2-6, 2-8 | Identical. |
| Gårding e-terms, c0'w at level 2, A0 w | 2-1, 2-14, 0-4 | Identical. |
| FNW covers the C^4 residuals only as a sum | 0-9 | Identical. |
| no written derivation; ASSEMBLY.md and code disagree | 0, 1.4 | Identical. |

### 4.2 Items I dispute or qualify

1. **"KILL" versus FIX-AND-RERUN.** The audit's KILL is a process verdict under its own stop rule (negative rounding reserve, uncovered terms), and it says so: "The verifier may well close after the fixes; this audit does not decide that." It could not run the verifier. I did, on a copy: with every item of both reviews added, the certified rate is unchanged to 0.2 % (sec. 2.2, 2.3). So the certificate is dead, the approach is not. If "KILL" is read as "plausibly fatal", I dispute it.
2. **(N1) and (N2) are "probably material".** Not for the rate. (N1) is large relative to the E_I dissipation but the E_I directions are damped 10x faster than the binding inner mode; the certified effect at the tail anchor is 3.9e-6. (N2): the audit's 4.6e-4 is the half cross-coefficient kappa/2; the effect on the top eigenvalue is kappa^2/(4 * 0.14) = 1.5e-6 [cert], and the B&B shows 3.9e-6 in U and 5e-7 in lb.
3. **Ramp monomials "~100x the coded W3 coefficient".** I get about 50x at sup_ramp|w''| = 5e-3 [est]. Either way the certified effect on the rate is below 1e-7 (row `m`).
4. **Row 0.5c, inner ramp: "in term (a)".** Incomplete. Term (a) covers the top-order term on the ramp only if the second integration by parts of (N1) stops at x = 16.25, which leaves pointwise boundary terms there. Done globally it leaves the commutator c (3 zeta' w''' + 3 zeta'' w'' + zeta''' w'), so the W3, B, W1 coefficients of term (a) are 2-3x too small (my 0-5).
5. **My own error, corrected.** I claimed the constant 1.67 fails on the outer ramp. The audit's kappa_out = 1.5883 is right (I get 1.5833 [cert]); 1.67 holds. Row 0-15 is corrected.

### 4.3 Items AUDIT_RAMP has and I missed

1. **Level-2 transport commutator on the inner ramp** (its 2.6b): -(a q + b cos theta) zeta' g paired with D2v, <= 1.271 |(a,b)| W1 ||D2v|| [its cert]. I had the level-0 version only. Size 6.6e-7 of CF3: negligible, but it is a term.
2. **Level-2 a-transport constant** (its 2.6a): certified need 39.1842 (Euclidean with b) against the coded 39.16; the deficit is absorbed by the V2V3 term only through a Hardy-Poincaré step V2 <= 15.35 V3 that is not written anywhere. I had noted 39.19 > 39.16 and waved at "V2V3 slack" without saying that a lemma is needed.
3. **CABQ**: 0.273 is typed below 0.2731; covered because |(a,b)| is Euclidean (0.2736 <= 0.2889). I had marked it plainly OK.
4. **The G2 list** for the outer ramp (its sec. 5). I did not audit the outer ramp.

### 4.4 Items I have and AUDIT_RAMP missed

1. **beta's degree-0 parts in E_I.** The 0.657 Fn(alpha) part of the beta bound doubles the I0 forcing. The strip constant STRIP_BETA = 1.06e-9 that the code itself uses for beta gives a non-decaying forcing 2.75e-10 V0, i.e. 7.4x the flux floor F0 [cert]: sqrt(L_F) 9.87e-10 -> 2.87e-9 and the J-slopes x8.4, unless the strip part of beta is made state-dependent. (The same question exists for beta in E0, outside this review.)
2. **Floor structure.** The strip constant of |(a,b)| times psi1', psi2' is a degree-1 term; `tail_cert.py` fails in mode joint if it is put in L without the floor flag (I hit this). Size 1e-31, but the tail lemma's floor formula (`closure.py:143-147`) has no slot for degree-1 floors when QUAD=1.
3. **Level-2 inner ramp cannot be bounded in the obvious way.** zeta w''w'''' <= sup|w''| D4 ||D2v|| has radial degree 2.25 and breaks the radial lemma; zeta D' needs w^(5). The fix has to convert one factor to v (2-10). The audit's sec. 6 item 2 ("add the inner-ramp level-2 terms in rho-norms ... with x0.66") would run into this.
4. **Forcing part of |(a,b)| uses FnW instead of Fn** (`rhs_iv.py:144` against `mod_sups.py:25`): 7.3x to 66x too small (0-14).
5. **More down-rounded constants:** FN[0], FN[1], C3 (the dominant forcing of the whole inequality, -3.2e-4), FNW[0] (-1.8e-3), TQ d010, d110, d300 (-6.8e-4, -3.9e-4, -2.7e-6).
6. **mu1 provenance:** 0.1552 exceeds the L-int certificate cited for it (0.155167); it is valid only by `v2_T4_P2_0.25.json` (0.15552). The audit lists it as covered by lint_cert_p.
7. **Ramp table secant in S** (0-7): certified only at eta in {0, 1e-3, 1e-2}; holds at eta = 1e-3 with 0.01 % to spare.
8. **Coefficient-derivative terms at level 2** (2-9, and e', f', R0' in 2-5): no certified inputs exist for them.
9. **c0'w on (34.25, 60]:** `S_all` is not a sup of |w| there (harmless, c0' <= 5.6e-35).
10. **The rerun evidence itself** (sec. 2.2, 2.3).

## 5. Files (all in `writeup/lemmaU/bulk/fix5/g0_indep/`)

| file | content |
|---|---|
| `G0_INDEP.md` | this report |
| `G0_INDEP_sections0-3_pre_audit.md` | sec. 0-3 as written before `fix3/audit/` was opened (sha256 c03b0d8e...) |
| `g0_numbers.py`, `g0_numbers.json` | arb sizes of sec. 2.1 and the rounding table of sec. 1.3 |
| `check_identity.py` | sympy Euler-operator check of the level-0 top-order identity |
| `check_kappa_outer.py` | arb: sup over [57,60] of \|zeta'\| + psi' zeta = 1.5833 |
| `vp/` | copy of 8 verifier files; `rhs_iv.py`, `constants.py`, `unc.py` patched (switch `G0FIX`; with `G0FIX` empty and `RAMPMULT=1` it reproduces fin12f4). Reads its inputs from `step3nl_v33/`. |
| `vp_original_sha256.txt` | hashes of the 8 files at copy time (equal to the proof tree now) |
| `apply_patch.py`, `apply_patch_constants.py` | the patch, as applied |
| `cell.py`, `run_cells.sh`, `cells.jsonl` | per-switch B&B bounds of sec. 2.2 |
| `run_full.sh`, `g0all_tail.out`, `g0all_1.5e-4.out`, `g0all_status.txt`, `g0all_vs_fin12f4.json`, `vp/rate_table_g0all_a0.00015.json`, `vp/verify_log_g0all_a0.00015.txt` | the full run of sec. 2.3 |

Inputs reviewed (sha256 prefixes): `verifier/rhs_iv.py` 4165b774b7083041, `constants.py` 3f74fdd6008e0a8c, `closure.py` d1a1ac116840b73b, `unc.py` 2a4c4eab40b5231a, `bb.py` 76bc1e24d787f53b, `ival.py` 74b65bc2505a30cd, `tail_cert.py` b346e6369e25f6bb, `ASSEMBLY.md` 551e72ebe183366f, `inner_closure_float2.py` 1b99b45e8e8ae293.

Nothing in the proof tree was written, edited or rerun in place.
