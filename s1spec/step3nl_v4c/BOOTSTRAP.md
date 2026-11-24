# Step 3 nonlinear bootstrap — closed system of estimates (design, Nov 6 night)

Supersedes the informal "three-region scheme". Everything below is in the even class, x = sigma - sigma_tip,
rho = e^{-phi}, phi = (1+x^2)^{2/3}/4, chi = chi_32(sigma) (1 for sigma <= 29, C^3 smoothstep to 0 at 32).

## Unknowns and equations
u = alpha psi3^ + w on the graph chart, w~ = chi w.  Graph equation (exact, graph_struct.pkl):
  u_tau = F(u) - a (p + u - u'q/alpha_g) + b (sin th + u' cos th/alpha_g),
  F(u) = -J^-4 u'''' + c3(u,u',u'') u''' + F2(u,u',u'') + (1/4)(p + u - u'q/alpha_g),   c3 = 10 J^-6 u'u'' + b01.
  N0(u) = F(u) - L u = -(J^-4 - 1) u'''' + c3 u''' + [F2 + drift - L3 u],  L3 = L + d^4 (no u''' term in L).
  w_tau = L w + alpha r^ - beta psi3^ - a psi1 - b psi2 + N0 + a A1(u) + b B1(u),   beta = alpha_tau - lamt alpha.
Constraints <w~, f_j>_rho = 0 (f1 = chi p, f2 = chi psi2, f3 = chi F - a eta1 - b eta3) fix (a, b, beta):
  sum_i x_i (G_ij + E_ij(u)) = <w~, L* f_j> - <[L,chi] w, f_j> + alpha <chi r^, f_j> + <chi N0, f_j>,
  x = (a, b, beta), G_ij = <chi psi_i, f_j> (psi_3 := psi3^), G_13, G_23 ~ 4e-7 (f3 ~ adjoint mode), G_33 >= 0.0311.

## Energies
Inner (rho):         E = ||w~||^2 + eps ||w~''||^2,  eps ~ 1e-5 (theta-margin absorbs E2's Gaarding growth 14.34; the binding
                     constraint is the linear modulation coupling at E2: eps 421 ||w~||^2 <~ omega0/3).
Slow mode:           Lyapunov functional  Lyap = kappa alpha^2 + E  (triangular lemma, lyap_psi3b: kappa ~ 7.5e4,
                     rate omega' ~ 0.026 at the L^2 level; to be redone with the eps E2 term and the theta-form).
Intermediate (W_I):  v = zeta_I w' (zeta_I: 0 on x <= 6, 1 on x >= 12, cut off again before the chart edge R_c),
                     E_I = ||v||_W^2 + eps_I ||v''||_W^2  (L-int certified: rate >= 0.2435, H^2 dissipation 1/2).
Outer:               original variables, lifespan lemma (O) on balls of radius r/4 for the remaining time.

## Decomposition (revised): u = alpha psi3^ + alpha^2 h2^ + w (second-order slow manifold; forcing O(alpha^3))

## How each energy is closed (what it needs pointwise)
E0 level: pair N0 with rho chi^2 w.  The w'''' term:  int g w w'''' = int g w''^2 + 2 int g' w'w'' + int g'' w w'',
  g = rho chi^2 (1 - J^-4).  g' contains u''; its w''-part gives  2 int rho chi^2 d_{u'}(1-J^-4) w' w''^2  (needs sup|w'|).
  g'' contains u'''; its w'''-part is integrated by parts once more (needs sup|w|, sup|w'|, sup|u''|-coefficients).
  The w''' term (c3 w''') likewise: one integration by parts, then the u'''-part of (c3)' once more.
  Unavoidable cubic terms: int rho chi^2 c w (w'')^3 (from F2 and from the c3 step): need sup|w''| ->
  sup|w''|^2 <= C (||w''|| ||w'''|| + ...) <= C' ||w''||^{3/2} ||w''''||^{1/2}  (weighted interpolation), absorbed by Young
  into eps ||w~''''||^2 (E2 dissipation).
E2 level (simpler than E0): pair N0(u) - N0(u0) with w'''' directly (two integrations by parts on the test side):
  ||chi(N0(u) - N0(u0))|| <= sup|1 - J^-4| ||w''''|| + C(alpha, eta)(||w'''|| + ||w''|| + ...) ; sup|w''| by interpolation.
  (old note:) same scheme one level up (top term -J^-4 |w''''|^2; coefficients need sup|u'''| -> interpolation into
  ||w''''||, ||w'''''|| as above).
So the system closes with pointwise bootstrap bounds on (w, w') only, plus interpolation:
  x <= 12: from E (rho-Sobolev, loss <= rho(12)^-1/2 ~ 32);
  12 <= x <= R_c: from E_I (sup|w'| ~ W_I^-1/2 E_I^1/2, loss ~ x^3; w = w(12) + int w').

## Constants to certify (status)
 [x] theta-certificate (1/9, 1/32) for <L w~, w~>;   [x] psi3^ residual 1.28e-6;   [x] L-int (rate 0.2435)
 [x] <psi3^, f_j>, ||psi3^||, Gram data (psi3ip, lyap_consts2);   [x] <chi Q2(psi3^), f_j>, ||chi Q2|| (q2_ip)
 [ ] (Sob)  weighted Sobolev/interpolation constants on [0, 34] (rho) and [12, R_c] (W_I), explicit, interval arithmetic
 [~] (Jet)  E0 alpha-linear coefficients DONE (alin_E0.py); cubic coefficients float only (nlin_E0_float.py);
            sup bounds over the chart of the coefficient functions: 1 - J^-4, d_u, d_u' of it, c3, d_{u''} c3, F2-coefficients,
            evaluated on jet boxes  u^(j) in alpha psi3^(j) + [-eta_j, eta_j]  (interval arithmetic per x-chunk)
 [ ] (M)    modulation: (G + E(u))^-1 bound, |a|, |b| <= C_ab (||w~|| + sup-terms), |beta| <= (eps1 ||w~|| + ...)/G33
 [ ] (Comm) [L, chi] w pairings on sigma in [29, 32] (rho <= 2.6e-11) via sup|w^(j)|, j <= 3, from E_I
 [x] second-order slow manifold: h2^ certified (h2res 1.57e-6), third-order forcing <= 12.34|alpha|^3 (n3_forcing)
 [x] E2-level Gaarding constants (garding_cert.py);  [x] E2 nonlinear coefficients on jet boxes (nlin_E2.py)
 [ ] ||D_rho f_j|| (certified), weighted Gagliardo-Nirenberg constants, E0 cubic coefficients (rigorous, jet boxes)
 [ ] (N-in) E0 + eps E2 inequality; (H-in) E2 part; (N-int) E_I inequality incl. the zeta_I commutator at 6..12
 [ ] (O)    lifespan lemma, explicit c;  R_f >= 4 (c (1 - 4 a*))^{-1/4}
 [ ] assembly: forward invariance of {Lyap <= L*, E_I <= ..., outer OK}; basin alpha0.

## Size check (float, from earlier runs)
 forcing: alpha^2 chi Q2 (norm 22) direct, + beta-loop 0.31 alpha^2 * n3 (7.52) -> ~ 2.3 alpha^2: NOT binding in the Lyapunov
 frame (C_N ~ 1e-3).  Binding: chart/quasilinear smallness  sup|u'| <= ~0.1 on x <= 34  (|alpha| sup|psi3^'| = 39|alpha| at
 x = 30)  and the Sobolev losses for w  =>  alpha0 ~ 1e-3 .. 2.5e-3.

## Far-field coupling (Nov 7) -- avoid pointwise bounds by using rho << W
Pointwise bounds from the intermediate energy lose K_W^1/2 ~ x^P (P = 3: 3.3e4 at x = 32), and the inner cutoff
commutator <[L,chi]w, w~>_rho, if bounded with sup norms, gives a loop gain ~ 1.8e-9 CFAR^2 (CFAR ~ 4e3 => 0.03 > omega).
Fix: on the strip sigma in [29,32] (and on [20, 31] in general) rho/W <= 0.055 (rho(31)^2/W ~ 4e-13), so integrate by
parts and pair against the intermediate DISSIPATION ||v''||_W = ||w'''||_W instead, e.g.
   |int rho chi chi' w''' w| <= ||v''||_W (sup(rho/W) sup chi'^2)^1/2 ||w~|| = 0.17 ||v''||_W ||w~||  (absorbed: Young).
Same for the zeta_I commutator feeding E_I: move derivatives off w'''' onto (W zeta' v) and pair with ||v''||_W and the
inner dissipation ||w''||_rho:  int W zeta' v'' w'' <= ||v''||_W (sup_ramp W/rho)^1/2 sup|zeta'| ||w''||_rho, W/rho <= 4e-3
on [12, 18]  =>  flux <= (1/8)||v''||^2 + 8.8e-4 ||w''||^2, absorbed by E's (1/9)||w''||^2 if the Lyapunov weight
mu on E_I is <= ~60.  Lyapunov functional:  Lyap = kappa alpha^2 + E + mu E_I (+ eps_I E_I2).
Far-region cubic terms: split [12,20] (rho-norms + sup from E_I, K_W = 6.4e7) and [20,31] (W-norms, rho/W <= 80).
Float closure (inner_closure_float.py, CFAR = 300 placeholder): decay rate 0.021 at alpha = 0; alpha-linear growth
0.46 (alpha/1e-3) (alin 0.19 + modulation 0.27) -> E grows by ~e^3.6 for alpha1 = 1e-3; cubic terms cap sqrt(E) at ~7e-6.
Improvements available: pointwise eigenvalue bound of the alpha-linear form instead of sum of sups; exact rank-2
modulation form on the constrained space; finer x-regions for the cubic sups.

## Outer lemma (O) -- reduced requirements (Nov 7)
* The intermediate energy's outer cutoff zeta_out on y in [R-3, R] enters only through commutator terms weighted by
  W(R) (rho-type there: W(60) ~ 2.4e-23 with the agent_sob normalisation). Hence only CRUDE pointwise bounds are needed in
  the annulus [R-3, R] (graph over S1 with |u'| <= 1/2 and |u^(j)| <= C, j <= 4), not smallness.
* In original variables a fixed physical point at radius r moves outward in y like y' = y/4 (drift); it leaves the
  intermediate region at the time t_r with R lambda(t_r) = r, having been inside (y <= R - 6, E_I-controlled) at t_r'
  with (R - 6) lambda(t_r') = r. Elapsed time t_r - t_r' = lambda(t_r')^4 * O(4 log(R/(R-6))) << (r/4)^4.
  Remaining time after hand-off: T - t1 ~ lambda1^4/(1 - 4a*); for the far curve at distance d0 the time ratio is
  (lambda1/d0)^4 ~ 1e-7.  So (O) needs only SHORT-TIME stability with explicit constants:
  (O1) annulus: scale-invariant C^4 bounds at t_r' propagate over a time << (r/4)^4 with a factor 2;
  (O2) far curve (incl. the lens and the rest of the hairpin): C^k bounds at t1 (from step 4) propagate over
       [t1, T) (ratio ~1e-7); then |V| <= V_max gives no folding into B(x(t), R lambda(t)) since V_max lambda1^4 << d0.
  Both are localized energy estimates for CDF over very short times -- crude constants suffice.

## Far-field coupling: DECIDED DESIGN (Nov 7, evening) -- replaces the CFAR placeholder
Problem: inner cutoff strip (chi_32, sigma in [29,32]) and the intermediate inner ramp couple E and E_I both ways;
the product (W/rho)_ramp * (rho/W)_strip = (W(x_r)/W(x_s)) (rho(x_s)/rho(x_r)) decides feasibility -> SMALL P wins.
Integrate the commutators by parts down to f^2, f'^2 (design/coupling_design*.py, strip_sign.py):
  -<[L,chi]w, chi w>_rho = int A w^2 + B w'^2,  A = -2(rho chi chi')''' + 3(rho chi chi'')'' - 2(rho chi chi''')' + rho chi chi''''
                                                 + (q/4) rho chi chi'  [drift: GOOD sign],   B = 6 rho (chi'^2 - phi' chi chi') >= 0.
  float: A/rho in [-8.2, +3.6] (A+ only on sigma ~ [29.8, 30.6]), B/rho <= 5.76.
  Same identity for the intermediate ramp with (W, zeta, g = w'): A_z, B_z; the drift term (q/4) W zeta zeta' >= 0 there (BAD).
Bounds:  strip  int B w'^2 <= sB ||v||^2_W ;  int A+ w^2 <= 2 int A+ K_rho(xh) (..E..) + 2 (int A+ Phi) ||v||^2_W (Hardy from xh,
         Phi = int_xh^x W^-1);  ramp  flux <= a_f ||w~'||^2 + b_f ||w~''||^2  with ||w~'||^2 <= ||w~|| ||w~''|| + c1 ||w~||^2.
  => window  mu_min = (sB/mu1_strip + sA/mu1_min)/thI  <=  mu  <=  mu_max = thE min_t (om0 + t^2/9)/(a_f (t + c1) + b_f t^2).
Float results (thE = thI = 1/3, actual blended W):
   P = 1.5, ramp [10.25,16.25]: mu in [1.4e-4, 7.5e-3] (ratio 53);  P = 2: [3.2e-3, 0.114] (ratio 36);  P = 3, [12,18]: ratio 7.
   w-level E_I (energy of w itself) FAILS (drift flux into E_I at its inner ramp vs om0); P = 3 v-level fails via the Hardy
   term for A+.  => CHOICE: v = zeta w', W = exp(-phi_I), phi_I' = (1-s) 2P x/(1+x^2) + s phi', s = S((x-30.25)/8),
   P = 2 (backup 1.5 / 1.75), zeta: 0 for x <= 10.25 (sigma = 8), 1 for x >= 16.25; outer ramp [R-3, R], R = 60.
   L-int (lint_cert_p.py, XA = 10.25, sigma >= 8 so no CAPD-box dependence): mu1 >= 0.1552 (P=2), 0.1454 (1.75), 0.1241 (1.5);
   blend >= 0.3357 / 0.2785 / 0.2211;  rho part 1.0378; far tail lint_tail (mu1 >= 0.003 x^4/3).
E_I = ||v||^2_W + eps_I2 ||v''||^2_W  (level 2 needed for sup|w''| on [12, 34.25] in E's E2 terms; its ramp flux pairs
  ||w~''''||_rho (E2 dissipation, weight eps/2) with ||v''''||_W  =>  eps_I2 ~ eps/(mu (W/rho) c^2) ~ 1e-5..1e-4).
Lyapunov  Lyap = kappa alpha^2 + E + mu E_I.  Pointwise far bounds: sup|w'| <~ (K_W (..))^1/2 E_I^1/2, K_W(34.25) = 2.3e6 (P=2).
Constants still to certify for this design: (C2) level-2 Gaarding of L1 in W; (C3) sB, sA, sA0, a_f, b_f + level-2 ramp flux;
  (C4) W-forcing/residual norms; (C5) far nonlinear coefficients; (C6) ||zeta psi_j'||_W, ||zeta (A1, B1)(u0)'||_W, ||zeta psi3^'||_W;
  (C7) outer flux (W(60) ~ 3e-20 for P = 2).
