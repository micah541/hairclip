# Rigorous verifier for the FLOAT closure (inner_closure_float2.py)

**Method.** `rhs_iv.py` re-implements the float rhs with outward-rounded interval arithmetic (`ival.py`, no libm; angles by rational parametrization). Constants live only in `constants.py`. rhs = H + L - NEG (H: POS terms of x-degree>=2, L: lower-degree terms, NEG: explicit quadratics), so a Ly cell costs no dependency loss. `bb.py` is a branch-and-bound over the 5 shape variables (lam, theta, D4, theta_I, V4) with mean-value (AD) bounds; D4/V4 tails are cut per box by the radial lemma (POS(r)/r^2 nonincreasing, POS(r*)<=r*^2/2). alpha and Ly are cells. `closure.py` marches a rigorous upper solution: alpha' <= lamt*alpha+beta (exact exponential step with bb >= beta(a_up,L_up)), Ly' <= 2*r*Ly with r the certified B&B bound over alpha in [0,a_up], Ly in the step enclosure (alpha-monotone). Finally a scale-invariance tail lemma (anchor with rate_s<=-rho; floors dropped) gives alpha, Ly -> floor sqrt(L_F)~2.4e-8, alpha_F~8e-7 (b,c floors make the limit a floor, not 0).

**Monotonicity check.** P-class arithmetic has no subtraction and all constants are asserted >=0, so W1,W3,V1,V3,S[...],supv1,supv2 and alpha enter only with nonnegative coefficients (replacing them by upper bounds is valid). Only NEG has negative terms and it involves only free variables. Numerically checked in `check_vs_float.py` (interval rhs vs float at 300 points, relative width 2e-12; alpha-monotonicity; radial lemma 0/200 violations). `sanity.py`: random/float-grid rates never exceed U; scale lemma 0/4000 violations.

**Added terms** (`ADD` env in constants.py): (a) I0 ramp coupling 0.66*(1.85e-2 D4+...)(V1+1.67V0); (b) 9.9e-12*1e-2*V0 in I0; (c) ramp flux 1e-17 in I0 and I2; (d) strip cubic = 2.6e-3 x [16.25,31.25] cubic terms with E_I sup-bounds.

**Results.**
| alpha1 | ADD | closes | rigorous max sqrt(Ly) | float | cap 7.5e-6 margin | runtime |
|---|---|---|---|---|---|---|
| 1e-3 | bcd | yes | 2.014e-6 | 1.86e-6 | 3.7x | 203 s |
| 2e-3 | bcd | yes | 4.874e-6 | 4.14e-6 | 1.54x (2.05x vs 1e-5) | 471 s |

Tail lemma anchors: (a*=6.8e-5, L*=4.0e-12, rho=2.7e-3) and (a*=5.3e-5, L*=2.4e-11, rho=2.9e-3), alpha-condition OK. Rate table: `rate_table.json`; log: `verify_log.txt`. (Run in parallel, 2 processes; runtimes include contention.)

**Term (a) does NOT close.** With (a) on, the certified lower bound (`cem_probe.py`, point evaluations) of the max rate is +5.57 at alpha=0, scale-invariant (Ly=1e-16..1e-12), at shape lam~0.99, theta_I~4e-3, u~0.124, q~4: the D4 coefficient 1.85e-2*0.66 couples D4 to V1+1.67 V0 and D4's only damping is -(eps/2)D4^2, eps=1e-5. Scan of the D4 coefficient (AR1), alpha=0, Ly=1e-12: 1.85e-2 -> +5.69; 5e-3 -> +0.186; 2e-3 -> -0.011 (baseline value). So the threshold is between 2e-3 and 5e-3 (about 4-9x smaller). This is a sensitivity statement, not a fix; no constants were tuned. Terms b, c, d do not change rate lb at this precision.

**Caveats.** (1) float grid underestimates the true max (alpha=0,Ly=1e-12: float -0.0222 vs certified lb -0.0110; alpha=2e-3: 0.63 vs 0.713), so the rigorous sqrt(Ly) peak is higher than float. (2) Tail lemma relies on joint (alpha,x)-degree>=2 scaling (structural argument + 0/4000 numeric test) and the radial lemma (argued + numeric). (3) B&B leaves accepted at tol 5e-4/rtol 0.5%; the bound U is certified, lb only informational. (4) alpha levels rounded up to 1e-5 multiples. (5) ADD="" closure run not done; ADD="" only used for float comparison. (6) Not formally machine-verified; float rounding handled by relative inflation 5e-16 per op rather than directed rounding modes.

---
## Follow-up (R = 60 constants): certified ramp coefficients, time integrals, hand-off scan   (ADD=abcd, RAMP=table)

**1. Term (a) with certified coefficients.** `constants.py` now reads `../agent_int/planB/planB_ramp.json` (key `table`; the key `RAMP_linear_deviation` lives in `planB_constants.json`, same numbers rounded). Coefficients of D4, W3, B, W1, A (= w'''', w''', w'', w', w) in the I0 term `0.66*(...)*(V1+1.67 V0)`: c(alpha<=lev, S) = c0(lev) + [(c(lev,eta=1e-2) - c0(lev))/1e-2] * S, lev = smallest tabulated alpha in {2.5e-4, 5e-4, 1e-3} >= the cell's alpha bound (alpha <= 2.5e-4 uses the 2.5e-4 values; a cell above 1e-3 is refused, so the alpha1=2e-3 run is outside the table). S_ramp = max(sup|w|, sup|w'|) on [10.25,16.25] from the KREG bounds on [8,12] and [12,16.25] (the certified box is |w|=|w'|<=eta, so both are included; they coincide numerically), as an interval function of the shape variables, so it enters monotonically with nonnegative coefficients. First cells of alpha1=1e-3 have a_up = 1e-3(1+~1e-9): treated as level 1e-3 plus 1e-6 absolute slack in the rate (documented in closure.py).
Effect on the rate: negligible (alpha=1e-3, Ly=1.7e-12: U = 0.17683 without (a), 0.17685 with it). Sensitivity: multiplying all table coefficients by 2, 3, 5, 10, 30 gives U = 0.17688, 0.17692, 0.17736, 0.7268, 11.09, i.e. the real coefficients have a factor of roughly 5-10 of room (fixed 1.85e-2 was ~10-25x too large).
**S_ramp along the trajectory (rigorous worst case over ALL shapes with Lyap<=Ly_up; S <= 14.08*sqrt(Ly)):** it does NOT stay below 1e-2 for large hand-off. Max over the trajectory (upper bounds): see table. Where it exceeds 1e-2 the linear-in-S model is extrapolated beyond the certified box (flagged "S>1e-2" below). Rate-wise it is irrelevant (see the sensitivity), but it is not certified there.

**2. Time integrals** (columns "int"): from Lyap' <= r_theta*Ly - theta*NEGV on each step cell (theta=0.3, NEGV = mu MU1 V0^2 + (mu/2+mu epsI MU2) V2^2; r_theta certified by B&B with the V-dissipation reduced by theta), telescoped, plus the tail lemma with the same theta (tail anchor rate_s_theta equal to the theta=0 one). V1 = the code's V1 = sqrt(V0 (V2 + C1W V0)), bounded by C1W int V0^2 + Cauchy-Schwarz (affine majorant). **Because the (b),(c) floor forcing (c0' w term, flux 1e-17) does not decay, int_{tau1}^infty is not finite from these bounds; the bound is affine in T:** int_{tau1}^T V_i^2 <= const_i + slope_i (T - tau_t) with slopes V0: 2.15e-16, V1: 3.02e-16, V2: 6.67e-17 per unit tau (all runs), tau_t = tail start (column). If b,c are dropped the integrals are finite (slope 0).

**3. Hand-off scan** (tau measured from the hand-off, Lyap(tau1)<=Lambda1; alpha1 = |alpha(tau1)|). All runs close (barrier + tail), peak sqrt(Ly) far below the 7.5e-6 cap:

| alpha1 | sqrt(L1) | peak sqrt(Ly) | margin to 7.5e-6 | max S_ramp | S<=1e-2? | tau_t | int V0^2 | int V1^2 | int V2^2 | runtime |
|---|---|---|---|---|---|---|---|---|---|---|
| 1e-3 | 1e-6 | 1.569e-6 | 4.8x | 2.21e-2 | no | 10.6 | 4.69e-10 | 6.58e-10 | 1.45e-10 | 374 s |
| 1e-3 | 5e-7 | 8.514e-7 | 8.8x | 1.20e-2 | no | 10.4 | 1.32e-10 | 1.86e-10 | 4.11e-11 | 362 s |
| 1e-3 | 3.5e-7 | 6.393e-7 | 11.7x | 9.0e-3 | yes | 10.1 | 7.21e-11 | 1.01e-10 | 2.24e-11 | 150 s |
| 1e-3 | 2.5e-7 | 4.954e-7 | 15.1x | 6.98e-3 | yes | 10.6 | 4.19e-11 | 5.89e-11 | 1.30e-11 | 163 s |
| 5e-4 | 1e-6 | 1.157e-6 | 6.5x | 1.63e-2 | no | 7.4 | 3.14e-10 | 4.42e-10 | 9.76e-11 | 177 s |
| 5e-4 | 5e-7 | 5.885e-7 | 12.7x | 8.3e-3 | yes | 7.7 | 8.02e-11 | 1.13e-10 | 2.49e-11 | 71 s |

(plus affine slopes above.) **Largest hand-off that closes: sqrt(L1)=1e-6 at alpha1=1e-3 and also at 5e-4 (margins 4.8x / 6.5x to 7.5e-6); the largest with S_ramp <= 1e-2 (everything inside the certified box) is 3.5e-7 at alpha1=1e-3 (not yet tried between 3.5e-7 and 5e-7; the S<=1e-2 threshold is sqrt(Ly)<=7.1e-7) and 5e-7 at alpha1=5e-4.** At alpha1=1e-3 the peak is mostly alpha-forced growth (peak/hand-off ratio 1.6 -> 2.0), not hand-off driven.

**Caveats.** (i) S_ramp is a worst-case-over-shapes bound (KREG constants are crude); for larger hand-off the coefficient model is an extrapolation (rate insensitive: factor 5 of room). (ii) alpha1=2e-3 cannot be run with the table (alpha>1e-3). (iii) theta=0.3 integral constants are not optimized; they cost one extra B&B pass per step. (iv) Far-field coefficients on [30,60] unchanged. Files: `verify_log_followup.txt`, `rate_table_followup.json`, `integ.py`, closure.py env `LY1SQ`, `THETA`, `RAMP`, `RAMPMULT`.

---
## Follow-up 2: (b),(c) state-dependent (QUAD=1), recertified far T3 coefficients, far-field slope sups   (alpha1=1e-3, ADD=abcd, THETA=0.3)

**(b)** now `mu*9.9e-12*S_all*V0` in I0 (degree 2, in H), S_all = rigorous sup|w| over ALL regions (KREG bounds for 0<x<=16.25, E_I/far bound sr0+sqrt(PHI34) V0 beyond), replacing the 1e-2 box placeholder.
**(c)** now `CFLUX*(V0^2+V1^2+V2^2)` in I0 and `CFLUX*(V2^2+V3^2+V4^2)` in I2 (V1,V3 = the code's interpolation bounds). **CFLUX is an input, not derived here**: I have no certified commutator/ramp-flux constant (the Lemma O constant 3.2e-18 comes with crude G's, not with state norms), so CFLUX=1e-17 is only the old placeholder magnitude. Sensitivity (B&B cells alpha<=1e-3,Ly~1.45e-13 and the tail-anchor cell alpha<=7e-5,Ly~2.4e-13): CFLUX = 1e-17, 1e-3, 1e-2, 3e-2 change nothing (U = 0.3171 / -0.00268); at 1e-1 the tail cell becomes positive (+0.042), so closure tolerates any CFLUX <= ~3e-2 (so "c tiny" is safe by >10 orders). If the true commutator constant is below that, term (c) can be absorbed.
With both quadratic there are no degree<2 terms: the tail floor is 0 (sqrt(L_F)=0, Lyap -> 0) and **the integrals are finite** (slopes = 0).
**Far T3 coefficients** (r60_T3_boxw1_0.02.json, max over 16.25..60): d001 2.27e-2, d010 1.67e-3, d200 0.0153, d110 1.086, d300 3.314; d100 1.1e-4, A0, A1 unchanged (constants.TQ). Valid only if sup|w'|<=0.02 on x>=16.25.

| sqrt(L1) | alpha1 | peak sqrt(Ly) | margin to 7.5e-6 | max S_ramp | far sup|w'| (E_I/KW, sI1) | far sup|w'| =max(E0,E_I) | far sup|w| | tail rho | J0=int V0^2 | J1=int V1^2 | J2=int V2^2 | runtime |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 3.5e-7 | 1e-3 | 6.393e-7 | 11.7x | 9.0e-3 | 1.398e-2 | 1.398e-2 | 7.11e-3 | 2.0e-3 | 7.211e-11 | 1.013e-10 | 2.238e-11 | 165 s |
| 2.5e-7 | 1e-3 | 4.954e-7 | 15.1x | 6.98e-3 | 1.083e-2 | 1.083e-2 | 5.51e-3 | 3.0e-3 | 4.189e-11 | 5.885e-11 | 1.300e-11 | 180 s |

J's are rigorous upper bounds for int_{tau1}^infty (finite; affine slopes 0). **Item 2:** the E_I-based sup|w'| on x>=16.25 (worst case over all shapes with Lyap<=Ly_up; it is for x>=16.25, hence also [16.25,30]) peaks at 1.40e-2 (3.5e-7) and 1.08e-2 (2.5e-7): both <= 0.02, so the 0.02-box coefficients apply (headroom 1.43x, 1.85x); it is above 1e-2 for the 3.5e-7 case, i.e. the old 1e-2 far box would NOT have applied there. The far sup|w'| bound is 21.9*sqrt(Ly), so 0.02 holds for sqrt(Ly)<=9.1e-7. Terms (a) ramp S_ramp stays <=1e-2 for both. Files: `verify_log_followup2.txt`, `rate_table_abcdQ_L*_a0.001.json`.

---
## Final pass (rigorous rounding, tail lemma, leaf certificates)

**1. Rounding.** `ival.py` no longer inflates by 5e-16 (that scheme was also sound, but by an argument about two successive roundings): every `+ * / sqrt` is an IEEE binary64 correctly rounded operation, and each bound is moved one double outward with `np.nextafter` (`fd` floors at 0; exact zeros stay 0). Gradient (AD) arrays use the same rule. Constants: decimal literals via `float()` (correctly rounded) +- 1 ulp; every constant that used to be formed by float arithmetic (`1/32-PI3`, `SB+SA`, `2*IA`, `CAB*QPSI`, `CAB*(38.5+0.66)`, `MU*EPSI`, `MU1(1-theta)`, ...) and every `1/sqrt(.)` is now computed in 200-bit `arb` from decimal strings and converted outward. The mean-value (MVF) penalty in `bb.mvf_rate_ub` is rebuilt with outward rounding at each step (and uses `g^2/(4n)` unless `g >= 2nw(1+1e-12)`, where the endpoint formula is exact). Underflow: nonzero intermediates are >> 1e-250; a final absolute pad of 1e-150 in `sub_signed` covers any denormal loss. `circ()` (rational angle map) keeps its 1e-15 absolute pad (<= 6 roundings of values in [0,1], error < 7e-16).
Rerun of the two production cases (alpha1=1e-3, ADD=abcd, QUAD=1, theta=0.3; `verify_log_final*.txt`, `rate_table_final_*.json`): **numbers are unchanged to display precision.** Peak sqrt(Ly): 6.392908e-07 and 4.953730e-07 (both identical to the previous digits); all step rate bounds differ by <= 7e-15; tail rho 0.00204 / 0.00303 (same), tau_t 10.05 / 10.55; J0,J1,J2 identical to 5 digits (3.5e-7: 7.21e-11, 1.013e-10, 2.238e-11; 2.5e-7: 4.19e-11, 5.885e-11, 1.300e-11); far sup|w'| (E_I) 1.3976e-2 / 1.0830e-2 (changes at the 1e-16 level, now slightly smaller since the old 5e-16 pads are gone); S_ramp 9.002e-3 / 6.976e-3. `check_vs_float.py` (interval rhs vs float2, now with the pre-recertification T3 maxima) still encloses the float rhs at 300 points. Runtime rose from ~170 s to ~17 min per case when both cases run concurrently (`np.nextafter` is slower than a multiply); cell cost ~6-20 s.

**2. Tail lemma.** `TAIL_LEMMA.tex/.pdf` (4 pages). Structural facts are now *machine-certified* by `tail_cert.py` (log `tail_cert_log.txt`): the verbatim `rhs_pn` code is executed on degree-tracking objects (valuation `nu`: sum=min, product=add, sqrt=half, max=min; top degree `p`: sum=max, ...; closure proved in Lemma 1), giving: radial (D4,V4 only) `H+L in F_2` (POS(r)/r^2 nonincreasing); fixed alpha, x scaled: `H in G_2`, `L in F_2` (so H/s^2 up, L/s^2 down, which is the Ly-cell rule); jointly in (alpha,x): all summands of H and L have valuation >= 2 (with QUAD; with the old constant (b),(c) terms the joint check FAILS, which is exactly the old floor). The random 0/200 test is no longer used. The remaining numerical inequality of the radial lemma, `POS(q*) <= q*^2/2`, is evaluated in interval arithmetic in every B&B cell/box (`find_qstar`, `qcut`). The tail contraction (Lemma 4) is a first-touch argument whose alpha-condition is checked in 120-bit arb; I corrected its formulation: the decay `tot <= -rho*Ly` is only used at `Ly = ell` (not for all `Ly <= ell`).

**3. Leaves.** Confirmed (proof in TAIL_LEMMA.tex, section on leaves): the returned `U` is a rigorous upper bound of the supremum for *any* `tol, rtol, minw`: a box is only dropped when its (outward-rounded) upper bound is `<= lb <= U`, an accepted leaf contributes its own upper bound, and the loop ends when no box is left; `lb` need not even be a valid lower bound for soundness (it is one: a point evaluation). Tolerances only control tightness: `U <= sup + max(tol, rtol|lb|)` unless a leaf was accepted by the minimal-width rule. `leaf_audit.py` (log `leaf_audit_log.txt`) on five representative cells: minimal-width leaves = 0, `U-lb <= max(tol,rtol|lb|)` in all cases (e.g. tail-anchor cell alpha<=7e-5: U=-0.00268, lb=-0.00318).

**4. CFLUX** remains an input (default 1e-17); closure holds for any CFLUX <= ~3e-2.

**New caveat found in this pass (affects the J_i, not the closure).** The first-touch barrier gives only *upper* bounds `Ly <= L(tau)`, `alpha <= a(tau)`. The integral estimates `theta*int NEG_V <= Ly1 + sum max(r_theta,0) Lup D` (and the tail version) use the cell rate bound on the actual trajectory, which needs the actual `Ly` to stay *inside* the cell (`Ly >= Ldn`), resp. `alpha <= kappa*sqrt(Ly)` in the tail; a state with `Ly << L_barrier` has an uncontrolled rate (forcing `L/Ly` grows). So the sharp J_i reported earlier (J0 7.2e-11, ...) are **conditional on that tube hypothesis**. Fully rigorous bounds using only `Ly <= L(tau)` (`crude_int.py`: `mu V0^2 <= Ly`, `mu epsI V2^2 <= Ly`, `Lambda_int = sum Lup_k D_k + L*/rho`):
| sqrt(L1) | Lambda_int | int V0^2 | int V1^2 | int V2^2 |
|---|---|---|---|---|
| 3.5e-7 | 2.03e-10 | 2.53e-9 | 4.65e-7 | 8.45e-5 |
| 2.5e-7 | 8.29e-11 | 1.04e-9 | 1.90e-7 | 3.45e-5 |
(finite, but the V2/V1 bounds lose the factor 1/epsI = 3.3e4 and are weak). A sharp unconditional version needs an alpha-additive (Young) treatment of the degree-1 forcing in the B&B (extra supremum of L^2/(4|A|) per cell) -- not done; tell me if the downstream lemma needs it.

## Fix plan (V1, F2, F3, F4, production)

Environment for all runs: `THETA=0.05 ADD=abcd FIX=1 S0V=0.0153 S1V=3.0 DELTA=0.2 MOD=ibp`, alpha1 = 5e-4. Logs `verify_log_fix5_L{1e-7,2e-7}_a0.0005.txt`, tables `rate_table_fix5_*.json` (march) and `rate_table_unc5_*.json` (march + unconditional pass, produced by `unc_run.py`), sups in `E_I_sups.json`.

**V1 (unconditional J's).** Method in TAIL_LEMMA.tex (V1 section): one extra B&B pass per step on H-(1-delta)N (L dropped) plus the absolute term Phi(alpha). theta = 0.05 is the largest usable value: at the tail anchor the condition H-(1-delta)(NEG-theta NEGV) <= 0 fails for theta >= 0.1 (rate +0.0105 at theta=0.1, -0.0015 at theta=0.05, delta=0.02) with the modulation constants of F3, which is why theta dropped from 0.3. `tail()` now also requires a grid delta_t with r1* <= 0 before accepting an anchor (delta_t = 0.02 in both runs). Results (T-tau_t = time after the tail start):

| sqrt(Lambda1) | J0 (V0^2) | J1 (V1^2) | J2 (V2^2) | J4 (V4^2) | slopes per unit (T-tau_t), V0/V1/V2/V4 |
|---|---|---|---|---|---|
| 1e-7 | 8.9156e-11 | 1.2525e-10 | 2.7673e-11 | 9.2247e-07 | 2.911e-15 / 4.089e-15 / 9.035e-16 / 3.012e-11 |
| 2e-7 | 1.6862e-10 | 2.3688e-10 | 5.2337e-11 | 1.7447e-06 | 1.326e-15 / 1.863e-15 / 4.117e-16 / 1.372e-11 |

These are larger than the theta=0.3 values (approx 1.9e-11 at 2e-7) by the factor 6 in theta and the new modulation constants. The earlier theta=0.3 tables (`*fix3*`, `*unc3*`) failed the tail condition with the IBP constants and are superseded.

**Production march** (alpha1 = 5e-4, no margin hidden: every cell is a rigorous B&B bound):
- sqrt(Lambda1)=1e-7: closes; 21 steps, peak sqrt(Ly) = 1.366e-7, tail start tau = 8.50, tail rho = 2.39e-3 (alpha-condition margin 106), floor sqrt(L_F) = 2.02e-8, alpha_F <= 8.41e-6.
- sqrt(Lambda1)=2e-7: closes; 18 steps, peak 2.514e-7, tail start tau = 8.38, rho = 2.98e-3 (margin 59), floor sqrt(L_F) = 1.81e-8, alpha_F <= 4.23e-6.
- S_ramp = sup(|w|,|w'|) on [10.25,16.25], worst case over all shapes: 1.924e-3 (1e-7), 3.540e-3 (2e-7); certified box 1e-2.
- E_I sups (x >= 16.25): 1e-7: sI1 = 2.987e-3, sI2 = 3.978e-2, ramp M1 = 28 sI2 = 1.114; 2e-7: sI1 = 5.496e-3, sI2 = 7.321e-2, M1 = 2.05 (exceeds 1.2). M0 = sup|w'| = S_ramp as above. sup|w'''| not produced (not needed).
- Modulation validity box: A_max <= 1.37e-7 (2.51e-7), W1_max <= 1.72e-6 (3.16e-6), alpha <= 5e-4 along the whole trajectory (post-check over all steps; `validity_box_ok` true in `E_I_sups.json`), so A, W1 <= 1e-3 holds with a margin of more than 300.

**F2.** S0_v = 0.0153, S1_v = 3.0 (= 0.094*32, placeholder for S1_v/y = 0.094) enter as coefficient inputs on [34.25,60] (`S0V`, `S1V`). GN identity used as `sup|v'|(3V3+1.31V2)`. The weighted GN and the zeta' ramp GN versions are withdrawn/not implemented.

**F3.** `mod_ibp.json` constants: BETA1=66.5, BETAF=0.657, q_beta=5.53 (QB), CAB=3.708+622.4|alpha|, q_ab=51.3, Phi_F coefficient 3.84 in (a,b). Strip terms (< 1e-8 per unit S) are NOT included in the code; they are negligible against the ~1e-3 margins but are an open item to add if a fully literal treatment is needed. The E0 alpha-term (MOD_ALPHA) is the exact certified 6.19|alpha| scaled by 66.5/63.5: this is the "built from the beta bound" reading, a conservative over-estimate of a positive term; if the term is the exact alpha-level form the rescale is unnecessary but harmless. The case is not decided in the code (kept scaled).

**F4.** Term (c) as absolute forcing: I0 = 6.1e-18, I2 = 1.55e-17 (placeholders from `lemmaS_cflux.json`), entering L as `MU*I0 + MU*epsI*I2` and, in the tail, as F0; tail lemma extended to L^1_loc forcing with floor L_F = 2 F0/rho (TAIL_LEMMA.tex).

**Caveats.** Placeholders: S1V, F4 flux constants. Unrun: `sanity_unc.py`. The conditional `integrals()` pass was not rerun for fix5 (superseded by V1).

### Fix plan, update 2 (certified F4 constants, S1V = 0.6, sqrt(Lambda1) = 5e-8)
Env: `THETA=0.05 S1V=0.6 FLUX0=4.1e-21 FLUX2=2.0e-20` (certified K0, K2 per unit tau from lemmaS_consts.json; defaults in constants.py now), rest as above. Tables `rate_table_unc6_L{5e-8,1e-7}.json`, sups `E_I_sups_fix6.json`, logs `verify_log_fix6_*`, `/tmp/y_*.out`. The 2e-7 point was not rerun. J = X/(theta n_k), X = Ly1/2 + sum r1+ L D + sum Phi D + L*/2 + Phi(a*)/rho; slope = Phi(a*) L_F/L* + F0 (ordering in J-parts: Lambda1-proportional = Ly1/2 + r1+ part; forcing = Phi parts; tail anchor = L*/2).
Results in the final message of the session; superseded fix5 numbers for the F4 and S1V placeholders.

### Fix plan, update 3 (assembly inconsistencies; 5e-8 production record)
Environment of every fix5/fix6/fix7 run: `THETA=0.05 ADD=abcd FIX=1 S0V=0.0153 S1V=0.6 (fix5: 3.0) DELTA=0.2 MOD=ibp`, F4 flux 4.1e-21 / 2.0e-20 (fix6, fix7). **S0V actually used = 0.0153** (exported in the launching shell; the constants.py default 0.011 was NOT used; 0.0153 >= Lemma S's 0.0088). The logs do not echo S0V.

**Strip terms, now literal.** constants.py `STRIP_BETA = 0.02*5.3e-8`, `STRIP_AB = 0.02*2.5e-7` (sums of all strip coefficients in mod_ibp.json per unit S, rounded up: beta 5.256e-8, a 1.29e-7, b 2.43e-7; note these are ~1e-7 per unit S, not < 1e-8, but S <= 0.02 makes them 1e-9 / 5e-9 absolute). Added to beta() and the tail alpha-condition (closure.py) and to Am in the (a,b) terms (rhs_iv.py). Rerun `fix7_L5e-8` / `rate_table_unc7_L5e-8.json`: nothing moves (J0 2.6422e-11 -> 2.6424e-11; peak sqrt(Ly) 7.928e-8, rho 2.17e-3, tail margin 161, delta_t = 0.2 unchanged). tail_cert structural verdicts identical with STRIP=0 and 1.

**Certified modulation sups at 5e-8** (`mod_sups_5e-8.json`, max over all 27 march steps of the mod_ibp bounds with A <= far A_max, Q <= 102637 Ly, alpha <= alpha_up, strip included; the tail has smaller alpha, Ly so is covered): sup|(a,b)| <= 3.36e-7 (hence sup|a|, sup|b| <= 3.36e-7), sup|beta| <= 5.8e-9. This replaces the estimate <~ 1e-5.

**5e-8 production results (alpha1 = 5e-4, theta 0.05, delta_t 0.2).** 27 steps, peak sqrt(Ly) 7.93e-8, tail start tau 8.99, rho 2.17e-3 (alpha-condition margin 161), floor sqrt(L_F) 5.49e-10, alpha_F 3.47e-7, S_ramp 1.12e-3. Far-field sups: sI1 1.733e-3, sI2 2.310e-2, ramp M1 = 28 sI2 = 0.646, A_max 7.93e-8, W1_max 9.98e-7 (box 1e-3 holds).

| sqrt(Lambda1) | J0 | J1 | J2 | J4 | slopes V0/V1/V2/V4 |
|---|---|---|---|---|---|
| 5e-8 (fix7) | 2.6424e-11 | 3.7121e-11 | 8.2016e-12 | 2.7340e-07 | 8.827e-19 / 1.240e-18 / 2.740e-19 / 9.133e-15 |
| 1e-7 (fix6, no strip) | 4.6389e-11 | 6.5169e-11 | 1.4399e-11 | 4.7998e-07 | 6.708e-19 / 9.423e-19 / 2.082e-19 / 6.940e-15 |

J0 breakdown at 5e-8 (shares of X): Ly1/2 7.6%, positive-r1 march 4.3% (Lambda1-proportional, 11.9%); march Phi 56.2% + tail Phi/rho 12.8% (forcing, 69%); L*/2 19.1% (tail anchor). At 1e-7: 17.4% + 8.2%; 33.5% + 8.7%; 32.3%.

### Fix plan, update 4 (tail_cert fix, MOD_ALPHA = 25, anchor eta = 0.01)
**Bug found and fixed (reported to coordinator):** in the MOD=ibp variant the x-degree-0 part of the (a,b) bound (3.84 FnW(alpha) + strip) multiplied alpha*V0 and alpha*D2v inside H, i.e. x-degree-1 terms in H (valuation 1 in s mode); H/s^2 is then not nondecreasing, so the cell bound H/s_hi^2 was not guaranteed. They now go to L (rhs_iv.py, `Amx`/`Amc`). tail_cert.py's s-mode check had flagged this (H lo = 1); all runs before fix8 had it (effect numerically tiny; the numbers below supersede fix5-fix7). tail_cert.py now ends with `TAIL_CERT PASS` (radial, s, joint on the floors-off part; floor terms are L-only; TAIL_LEMMA.tex section "Which conditions hold for which part").

**Run (5e-8, alpha1 = 5e-4, theta 0.05, env as before + MODA=25 ETA=0.01, tag fix9m25):** E0 alpha-term MOD_ALPHA = 25 (reviewer's worst case; replaces 6.19*66.5/63.5) and tail anchor cell Ly in [L*, 1.01 L*] (eta = 0.01; also used for the V1 delta_t check and the r1* pass). **It closes**; certified treatment = MOD_ALPHA 25.
Results: 29 steps, peak sqrt(Ly) 8.191e-08, tail start tau 9.69, rho 2.884e-03 (alpha-condition margin 130), floor sqrt(L_F) 4.770e-10, alpha_F 2.438e-07, S_ramp 1.153e-03; sI1 1.791e-03, sI2 2.385e-02, ramp M1 0.668, A_max 8.191e-08, W1_max 1.031e-06. J0 = 2.5423e-11, J1 = 3.5714e-11, J2 = 7.8909e-12, J4 = 2.6304e-07; slopes 6.995e-19 / 9.827e-19 / 2.171e-19 / 7.238e-15 (`E_I_sups_fix9_m25.json`, `rate_table_unc9_m25_L5e-8.json`). The 1e-7 point was not rerun with these settings. The largest closing MOD_ALPHA was not searched (25 closes).

### Fix plan, update 5 (final configuration at smaller hand-off alpha1; sqrt(Lambda1) = 5e-8)
Final configuration: MOD=ibp (MODA=25, strip in), anchor eta = 0.01, theta 0.05, S0V 0.0153, S1V 0.6, certified F4 (4.1e-21 / 2.0e-20), delta 0.2. Modulation constants (certified for |alpha| <= 5e-4) are used unchanged for smaller alpha (conservative); ramp table row 2.5e-4 used (all steps; conservative since alpha <= 1.25e-4 / 1.5e-4). Tables: `rate_table_fin10_a0.000125.json`, `rate_table_fin10_a0.00015.json`, summary `final_alpha_points_5e-8.json`.

| alpha1 | steps | peak sqrt(Ly) | tail tau | rho (margin) | floor sqrt(L_F) / alpha_F | S_ramp | sI1 / sI2 / M1 | sup\|a\|,\|b\| / sup\|beta\| |
|---|---|---|---|---|---|---|---|---|
| 1.25e-4 | 9 | 5.2135e-8 | 4.88 | 2.861e-3 (175) | 4.789e-10 / 3.295e-7 | 7.341e-4 | 1.140e-3 / 1.5183e-2 / 0.4251 | 2.146e-7 / 3.20e-9 |
| 1.5e-4 | 10 | 5.3074e-8 | 5.48 | 2.693e-3 (177) | 4.936e-10 / 3.433e-7 | 7.474e-4 | 1.160e-3 / 1.5456e-2 / 0.4328 | 2.183e-7 / 3.34e-9 |

J (const + slope (T-tau_t)), delta_t = 0.2:

| alpha1 | J0 | J1 | J2 | J4 | slopes V0 / V1 / V2 / V4 |
|---|---|---|---|---|---|
| 1.25e-4 | 5.6765e-12 | 7.9745e-12 | 1.7619e-12 | 5.8733e-08 | 8.365e-19 / 1.175e-18 / 2.596e-19 / 8.655e-15 |
| 1.5e-4 | 6.0272e-12 | 8.4672e-12 | 1.8708e-12 | 6.2362e-08 | 8.637e-19 / 1.213e-18 / 2.681e-19 / 8.936e-15 |

A_max 5.21e-8 / 5.31e-8, W1_max 6.6e-7 / 6.7e-7: validity box holds.
