# Certification plan: slow-mode coordinates (spectral-complement coercivity)

Goal: certify Re<Lv,v>_rho <= -omega ||v||^2 for even v perpendicular to {f1*, f2*, f3*} in L^2(rho), with
omega ~ 0.07 (float target: -0.0850; margin ~0.015). Here fi* = chi*(w_i/rho) - a_i*eta1 - b_i*eta3 are cut-off
adjoint eigenfunctions with tip bumps making them C^3-even, exactly as f3 in the fourth supplement.

Float evidence (adjcut.py): the abscissa on {chi psi_i*}^perp is -0.0850 for cutoffs on sigma in [29,32] and
[19,22], and -0.0848 for [12,15]. The exact spectral complement also gives -0.0850.

Adjoint data (step3/refine_wstar_all.py -> step3/wstar_all.json; frame alpha*wa - wb at sigma = 32):
- lambda = 1:     alpha = -0.230997556221866  (lambda recovered as 0.999999999999979)
- lambda = 1/4:   alpha = +7.084730070047681  (lambda recovered as 0.249999999999630)
- lambda = lamt:  alpha = +0.212944118491896  (lamt = -0.256153222457747, as certified)

Steps:
1. gen_hs3.py: replace the single adjoint generator "w" (aux rows ws0..ws3, parameter lamt) by three generators
   wa, wb, wc. Each needs its own aux rows and its own lambda parameter, with initial data from wstar_all.json
   (lt_frame(lambda), alpha). Update gexpr/fexpr, aux_for/aux_fun/aux_init and params().
2. Change omega: LAM = 0.07, or 0.075. Redo the Morse count (count_zeros.py: n(A_even) = 3 should persist,
   since the next H_sym eigenvalue is -0.2388) and the tail data (hstail.py -> tail_mid.pkl) at the new omega.
   Check the positivity-at-infinity margin with mu shifted by -(omega - 3/64).
3. run8.py: jobs Z for {wa, wb, wc, e1, e3}, plus Q for all pairs. The p, C, S generators are no longer needed
   for the constraints, but keep them if the modulation directions psi1, psi2 enter elsewhere.
4. gram_post8.py: build the f_i* from the new generators, with bump coefficients from the tip values. Then the
   3x3 Gram matrix G must be negative definite, with a certified Schur complement.
5. Lyapunov lemma: the defect of {f*}^perp from L-invariance, via residual norms as before. This gives the decay
   rate on the complement for the slow-mode trapping theorem (fifth supplement, section 7).
6. Then certify psi3 itself (eigenfunction, growth |y|^2.02) and h2, h3 (linear ODE problems with polynomial
   growth at infinity).

## Progress (Nov 6, 11:40)
- gen_hs4.py (adds adjoint generators wa: lambda = 1, wb: lambda = 1/4, inline exact eigenvalues), run9.py: all 24
  validated runs RESULT OK (1-3 s each), still at the OLD omega = 3/64.
- gram_post9.py: G2 (wa, wb block) negative definite, g11 = -18.37, det = 6.0 +- 0.69; Schur complement
  in [-2.18e-3, +7.2e-4], so NOT yet certified (midpoint about -7e-4; the old constraints gave -4.6e-4 in float).
  Cause: f_wa and f_wb are nearly parallel in the A^{-1} metric (det = 118.3 - 111.8), so the cancellation is
  lost when G entries are bounded as intervals before the Schur step.
- Fixes, in order: (1) compute det(G2) and S in shared-symbol affine arithmetic (keep the Aff entries, bound only
  at the end); (2) or orthogonalize, e.g. use wb - (g12/g11) wa, at the affine level; (3) tighten the wb runs
  (alpha = 7.08 inflates its enclosures). Then move omega to 0.07: Morse count, tail data, reruns.
- (12:58) gram_post10.py: Schur test in shared-symbol affine arithmetic gives no gain (S in [-2.7e-3, 1.1e-3]),
  because the uncertainty is in per-run remainders, not shared symbols.
- hsrun4 (BSUB env: stage-B step subdivision) with BSUB=4 leaves the raw integrals unchanged
  (relative radius ~1e-4 for I(wb,wb) and I(wa,wa)). So the step size is not the bottleneck.
- Tip-bump corrections are small (wa: |a| <= 1.3e-4; wb: 2.5e-5; w: 1.3e-6).
- So the g22 = 6.44 +- 9.4e-3 width must come from the boundary/ratio terms of the half-line formula
  (rho1, rho3 = mean-value ratios over W13-type denominators) for the wb generator. Next: print the radii of the
  rho1/rho3 and y1/y3 parts per entry; then rescale wb, or represent it relative to wa, so these terms are well
  conditioned.

## Decision (13:08): use the ORIGINAL certified constraints {f1, f2, f3} (rate 3/64); park the adjoint-constraint certification
Diagnosis: for the wb (lambda = 1/4) constraint all three parts of its Gram entry are about +-2e-2, the main term
included. The width comes from the S1 profile itself (tip position known to about 1e-4), amplified about 100x
because f_wa and f_wb are nearly parallel. Fallback if constants are ever too tight: tighten step 1
(sigma_tip to about 1e-8).

Slow-mode trapping with the original constraints (orthogonal modulation u perp_rho f1, f2):
- Let psi3^ be an explicit approximate eigenfunction with a certified residual r^ = (L - lamt) psi3^ in L^2(rho).
  Put phi3 = (I - Pi) psi3^, where Pi projects onto span{psi1, psi2} along {f1,f2}^perp. Then
  calL phi3 = lamt phi3 + (I - Pi) r^.
- Set c = <u, f3> / <phi3, f3> and v = u - c phi3. Then v lies in {f1, f2, f3}^perp, where the certified coercivity
  -3/64 holds.
- c' = lamt c + [<v, calL^* f3 - lamt f3> + <(I-Pi) r^, f3> c] / <phi3, f3> + nonlinear.
  The coupling is the certified adjoint residual (r3 ~ 5e-5, Pi-terms; see lyap_consts2.json).
- v' = P calL v + P (I-Pi) r^ c + P N(u), with P the projection along phi3 onto {f3}^perp.
- Cost compared with spectral modulation: a, b are linear in v (v is about 1e-3 of c psi3 for the hairpin), and the
  rate is 3/64 instead of 0.085.
- Next concrete items: (1) choose psi3^ (fit from the finite-difference eigenvector, smooth, growth |y|^2.02) and
  certify ||(L - lamt) psi3^||_rho and <phi3, f3> by interval arithmetic; (2) redo the Lyapunov lemma with phi3 from
  psi3^ (the g-term becomes (I-Pi) r^, which is small); (3) nonlinear jet-box bounds for N(c phi3 + v).
- (13:40) psi3^ DONE (untrusted): multi-domain Chebyshev collocation psi3md.py (21 elements on [0,36], degree 14,
  C^3 matching, even BCs at 0, free end at 36). lambda = -0.2561532205 (certified -0.2561532224577, so within 2e-9).
  Float residual ||(L - lambda) psi3^||_rho = 1.4e-7 on [0, 34.25], max weighted 5.5e-7. Coefficients in
  psi3md_14.json. psi3^(2,5,10) = -0.45014, 5.00063, 40.7448.
  Next: interval residual with certified profile hulls (step3/mon.txt MB/MC hulls, or hsrun3 MONALL runs);
  then <phi3, f3>, Pi-terms, Lyapunov lemma with phi3 = (I - Pi) psi3^.
- (13:45) Certified residual for psi3^: interval evaluation on monitored boxes is too crude. Profile hulls on boxes
  of width 0.01 vary by ~0.01, so the residual bound would be ~1e-2. Instead, use the CAPD Taylor method, as for
  the Gram integrals (build_A2 / A2_INT): per collocation element, one hsrun3 run with an extra integral row
  I' = rho(x) * r(x)^2, where r = (L - lamt) psi3^ is evaluated from the profile STATE (p, q, k, m in stage C) and
  psi3^ written as an explicit polynomial in the time variable (x = sigma - st; sigma = 8 - c tau in stage C).
  Layout: MC lines = [p q k m st ...] lo/hi pairs over rescaled-time boxes. Stage B (sigma in [8,32]) carries the
  profile as tail-theorem chunk parameters (gen_hs2.chunks: 24 chunks, 4 remainder amplitudes each, ~1e-4 at
  sigma = 8 and ~1e-20 at 32), so the r^2 row there uses the tail formula plus remainders.
- (13:50) DESIGN for a certified psi3^, mirroring f3 (replaces the polynomial psi3^ for certification):
  * Cutting psi3 off at finite sigma is NOT viable: [L,chi]psi3 ~ (q/4) chi' psi3 with psi3 ~ x^2.02, so the
    weighted residual is ~1e-2 at x = 32, ~3e-4 at 40 and ~1e-6 only near 50. Hence: no cutoff.
  * psi3^ := U on [sigma_tip, 32] and g on [32, inf), minus tip bumps a' eta1 + b' eta3:
    - U = TRUE solution of (L - lamt)U = 0 (profile from the certified tail/stage-C representation) with data
      J1 + beta*J2 at sigma = 32. J1 = exact jet of sigma^A (A = 1 - 4 lamt). J2 = exact jet of exp(S),
      S' = -(sigma/4)^(1/3) (WKB, use its exact derivatives). New generator rows (gen_hs4-style patch):
      u0..u3 with u'''' = -k^2 u'' - (5km + q/4)u' + (1/4 - 3m^2 - kp - lamt)u, NO e^phi scaling (|U| ~ 1e3 at 32).
    - g = sigma^A + beta*exp(S(sigma)) explicitly for sigma >= 32. Its jet at 32 equals J1 + beta J2, so psi3^ is C^3.
    - beta from float shooting (inward growth ~e^24, so beta ~ 1e-10). Tip defects u'(tip), u'''(tip) ~ 1e-9 are
      removed by bumps with a', b' ~ 1e-9.
  * Residual = bump terms (~1e-8) + explicit (L - lamt) g on [32, inf) (sigma^(A-4) and WKB errors, by hand).
    Zero elsewhere, because U is an exact solution for the true profile.
  * Inner products needed (<psi3^, f1>, <psi3^, f2>, <psi3^, f3>, norms): CAPD integral rows as in build_A2.
- (13:55) Float data for the certified psi3^ (psi3beta.py -> psi3beta.json):
  * Leading-order jet sigma^A at 32 excites the decaying mode (inward growth ~e^24): the tip conditions came out
    incompatible (u'''/u = -2.5). FIX: the 4-term polynomial branch sigma^A * sum_j a_j sigma^(-4j), with
    a_j = a_{j-1}(A-4j+4)(A-4j+3)(A-4j+2)(A-4j+1)/j: a = 1, -0.049801, -2.894766, -2878.93 (asymptotic; at
    sigma = 32 the terms are 4.7e-8, 2.6e-12, 2.5e-15).
  * With it: beta = 5.24e-15 (nearly pure polynomial branch), u'(tip) = 0 imposed, u'''(tip)/u = 1.0e-6,
    u(tip) = 1.5507, so psi3^ = U/1.5507.
  * Next: gen_hs5 patch with generator rows u30..u33 (L-rows, unscaled, data J1 + beta J2 as tight intervals); CAPD
    integral rows for <psi3^, f_i>_rho (i = 1, 2, 3), ||psi3^||, and tip values (for the bumps). Plus the explicit
    tail sigma >= 32 (g = series + beta e^S, residual from the truncated a_4 term) by hand.
- (14:00) gen_hs5.py (generator rows u0..u3 for (L - lamt)u = 0, unscaled; integral rows B0..B6, C0..C6 =
  <ws0 u>, <p u>, <C u>, <S u>, <u^2>, <e1 u>, <e3 u> weighted): the g jet at 32 encloses to 1e-56,
  (1115.18, 70.557, 2.2592, 0.0017377). The first CAPD run FAILED with "sqr overflow" after 2 s. Suspects:
  (i) the weight expression R in coeffs() is the scaled weight meant for e^{phi}-scaled generators, so u0^2*R may
  overflow at the start; check c['R'] in stage B against rho; (ii) unscaled u rows of size ~1e3 against the row
  scaling ZSHIFT/QSHIFT. Fix: scale the u rows like the adjoint rows (u~ = u * e^{-phi} or by a constant 1e-3) and
  use the matching weight.
- (14:05) Cause of the gen_hs5 overflow: the inward integration of (L - lamt)u = 0 is exponentially dichotomic.
  The decaying-at-infinity mode grows inward (rate (sigma/4)^(1/3) ~ 1.8-2), and the enclosure width grows
  ~4x per unit (wrapping, plus the chunk remainders times |u| ~ 1e3). Monitored hulls reach +-5e5 by sigma = 24.
  Direct CAPD integration of psi3 is therefore the wrong tool.
  NEW ROUTE: keep the explicit multi-domain polynomial psi3^ (residual 1.4e-7) and certify its residual with CAPD
  integral rows, with psi3^ entering as an EXPLICIT function of the stage time (no dichotomy):
   * stage B (sigma in [8,32], 24 profile chunks of length 1): refit psi3^ per chunk [8+j, 9+j] by degree-16
     Chebyshev and pass the coefficients as per-chunk parameters (the driver already sets npc parameters per chunk);
   * stage C (sigma in [tip, 8]): one global Chebyshev of degree ~60 on [tip, 8] in the rescaled time, inlined;
   * integral row r^2 rho with r = (L - lamt) psi3^ built from the profile (state in C; tail formula in B);
   * beyond sigma = 32: explicit series/WKB extension, residual by hand.
  Inner products <psi3^, f_i> are computed the same way, as rows psi3^ * f_i * rho.
- (14:20) CERTIFIED (interval arithmetic; CAPD profile hulls on 51,246 stage-C boxes, MONALL = 2e-4; psi3res_C.py):
  ||(L - lamt) psi3^||_{L^2(rho), sigma in [tip, 8]} <= 5.38e-3, max weighted |r| <= 4.7e-3 near x = 0.47.
  The floor is the S1 profile enclosure at the tip (m in [-1.7e-4, 4e-5], q width 1.2e-3, k width 7e-5), not the
  boxes: second-order Taylor forms change nothing. Fine for trapping, since the residual forces v ~ (5e-3/(3/64))|c|.
  Tightening it needs a tighter step-1 certificate.
  TODO: sigma in [8, 32] and the tail (expected <= 2e-4: profile ~ ray, k, m, p <= 2e-4 at sigma = 8 and decaying),
  using the tail-theorem bounds of the first supplement. Then inner products <psi3^, f_i> by interval quadrature on
  monitored hulls of the C, S, ws rows.
- (Nov 6) psi3^ CERTIFIED as an explicit C^3 function (psi3fun.py): psi3md2_14.json pieces + degree-7 Hermite
  corrections (exact, enclosed in arb) removing the interface jet mismatches (1e-9..2.5e-8, chebfit), imposing
  u'(0) = u'''(0) = 0, and matching the explicit extension c1 g(s) + c2 E(s) at X = 28 (s = x + stip);
  R g = -4 a_4 s^(A-16) exactly (telescoping), E = WKB exponential (residual factor ~1 but c2 = 1.7e-8).
  Residual (psi3res2.py, HALF-LINE norms as in the linear certificate):
     stage C (CAPD hulls, 51,246 boxes): ||r||^2 <= 1.4479e-5 ; stage B (ray + chunk bounds, Taylor form of the
     residual itself -- the first version had a dependency blow-up to 6.6e-3): 2.03e-10 (st-uncertainty term);
     tail sigma >= 32: 5.1e-19.   TOTAL ||(L - lamt) psi3^||_rho <= 3.8051e-3 (full line 5.38e-3), ||psi3^|| = 7.6996.
  Inner products (psi3ip.py, hulls of outA5.txt, half line): <psi3^, f1> in [-0.8012, -0.7648],
     <psi3^, f2> in [-2.5370, -2.4827], <psi3^, f3> in [0.031107, 0.033604] (float 0.032356).
  Linear triangular lemma with phi3 = (I - Pi) psi3^/N' (lyap_psi3.py): N' >= 0.031105, ||phi3|| <= 247.6,
     ||g|| <= 0.1223 (was 67.7 with phi3 ~ f3), G0 = 4.19 (dominated by eps1 ||phi3||^2), omega' up to 0.0290;
     omega' = 0.0261 with kappa = 7.5e4, C = 10.64  (was 0.0213, C = 14.8).
  Remark: ||phi3|| = 247.6 >> 1/||f3|| = 51: psi3^ and f3 (adjoint mode) make an angle with cos = 0.216 (non-normality).
- (Nov 6, late) KEY OBSTRUCTION found for the nonlinear energy: the certified coercivity has no H^2 margin, and the
  convexity trick only gives theta = (omega - omega')/(2 mu_-) ~ 0.0035 (mu_- = 3.34 at the tip, mu_float.py).
  The quasilinear term (J^{-4} - 1) u'''' ~ 6.2 alpha w'''' needs theta ~ 0.1.
  => RE-CERTIFIED the linear estimate with an H^2 margin (same constraints f1, f2, f3):
     THEOREM (theta-form, step3/): for even v in D(L), v perp f1, f2, f3:
         Re<Lv,v>_rho + (1/9) ||v''||_rho^2 <= -(1/32) ||v||_rho^2.
     Pipeline = original with env FORM_MULT=9/8 (=1/(1-theta), exact in CAPD) and OMEGA_Q=1/32:
       positivity at infinity (theta-form): mu >= 0.001 x^(4/3) on x >= 30  [positivity_cert_T1.json];
       tail theorem (hstail, lam = -1/32, theta-form): eta = 0.0395, Lip 3.1e-3, transversality margin 1.0;
       Morse count (out2P_T1.txt): 2 zeros of W01 + [W01(0) W13(0) < 0] = 3 = n(A_even);
       34 CAPD runs (run8.py T1: outT1*.txt), gram_post8: Schur complement in [-4.93e-4, -1.73e-4] -> G3 < 0, n(G) = 3.
     [certificate_step3_theta_T1.json].  Float scan (float_theta.py): (theta, omega) = (0.1, 1/32) Schur -4.2e-4,
     (1/9, 3/128) -5.4e-4, (0.2, 1/64) -3.1e-4: margins of this size are certifiable.
  Insight for the bootstrap: alpha(tau) decays like e^{-0.256 tau}, so every perturbation PROPORTIONAL TO alpha is
  time-integrable (growth factor exp(C alpha0/0.256)); the energy needs omega_eff > 0 only for the terms quadratic in w.
  Cutoff issue: with the ORIGINAL constraints, the inner energy must use the SAME cutoff chi_32 (then <chi psi_i, w~> = 0
  exactly); a cutoff at R ~ 8 produces a coupling 30 ||w|| <(1-chi_R) psi_2, w~> ~ 3.4 ||w|| ||w||_ann (fatal).
  So: w~ = chi_32 (u - alpha psi3^), constraints <w~, f_j> = 0 (j = 1, 2, 3) fix (a, b, alpha_tau).
- (Nov 6, night) BOOTSTRAP ARCHITECTURE (decided after float scans; supersedes the vague 'three-region scheme'):
  Float facts: (i) rho_beta = exp(-b(1+x^2)^(2/3)) with b < 0.2 is NOT coercive on {f}^perp (beta_scan.py: +0.14 at b = 0.1);
  (ii) the hybrid weight W_h (rho for sigma <= 8, smooth blend to (1+x^2)^(-p), p = 3.5, back to rho-like for sigma >= 36)
  keeps the constrained abscissa (-0.0714 / -0.0475 at theta = 1/9; hybrid2.py) and cuts the Sobolev loss at x = 34
  from 9.6e5 to 1.25e3 (360 at 25, 670 at 30).
  Why plain rho fails: pointwise control on 8 <= |x| <= 34 is needed (region-II coefficients; the constraint-cutoff
  commutator at sigma in [29,32]); rho-Sobolev loses 1e4..1e6 there.  The modulation b (|b| <= 31||w~||, since f2 is not an
  adjoint mode) drives a far-field translation; through W-Sobolev this looks 100x larger than it is.
  SCHEME:
   1. Inner: rho-energies E0 = ||w~||^2, E2 = ||w~''||^2, w~ = chi_32 (u - alpha psi3^), constraints f1, f2, f3
      (3-parameter modulation a, b, beta = alpha_tau - lamt alpha), theta-certificate (1/9, 1/32), Pi3-loss 0.0165.
      E0 + eps E2 closes linearly: the theta-margin -(2/9)E2 absorbs E2's Gaarding growth (max(-mu2) = 14.3 at the tip,
      |e1| <= 16.2, |e0| <= 10.7; garding2.py) for eps < ~4e-4.
   2. Intermediate: W_h-energy for v = zeta_I w' (derivative kills the translation mode: psi2' = 0 on rays) and its H^2
      level; positivity of the W_h-form on [6, oo) (L-int, positivity.py-type certification); no constraints needed.
      w on [8, 34] = w(8) + int w'; w(8) from the inner energy (rho-loss 7.5).
   3. Far field |y| >= R_f ~ 25-30: frozen in original variables (explicit lifespan/local-energy estimate; time ratio
      T_rem/(r/4)^4 <= (5/R_f)^4/(1 - 4a*) ~ 1e-3), released with data from the intermediate energy (W_h-loss 360-670).
   4. Couplings: inner <- intermediate (coefficients on 8..34, commutator at 31..34 via the localized pairing
      rho(31) C sup|w|^2); intermediate <- inner (boundary commutator at 6..9; beta, a forcing).
  Expected basin: alpha0 ~ 1e-3 (k_max ~ 3300 for c3k = -3.3); r^ floor (step-1 profile) and the O(alpha^2) forcing
  (46 alpha0^2) are the next bottlenecks (fixes: tighter step-1 tail at S = 30 + CAPD to 8; second-order h2).
- (Nov 6, night) LEMMA L-int CERTIFIED (lint_cert.py, lint_tail.py; ball arithmetic, step-1 tail bounds, stip interval):
  weight W_I: phi' = 6x/(1+x^2) (i.e. (1+x^2)^-3) on x <= 30.25, C^3 smoothstep blend to (1/3)x(1+x^2)^(-1/3) on
  [30.25, 38.25], rho-type beyond.  Operator on v = w':  L1 = -d^4 - k^2 d^2 - (7km + q/4) d - 8m^2 - (5/2)kp.
  For v in C_c^oo((12, oo)):  Re<L1 v, v>_{W_I} <= -1/2 ||v''||^2 - int W_I mu1 v^2  with
     mu1 >= 0.2435 on [12, 30.25] (min at 12), >= 0.5625 on the blend, >= 1.030 on [38.25, 68.25],
     mu1 >= 0.003 x^(4/3) on [68.25, oo)  (smallest 0.0093).
  So the intermediate (slope) energy decays at rate >= 0.2435 with H^2 dissipation 1/2. (Float: the polynomial weight is
  NOT positive below x ~ 11: the anti-diffusion 2 W''/W ~ 2p^2/x^2 wins; lint_float.py.)
- (Nov 6, evening) STEP-1 TIGHTENING DONE (s1_tight/; same eps-normalisation Z_0(8) = eps as step 1).
  Tail at S = 8: (k,k',k'') = eps*Yhat(S) + P(S)(0, d, conj d), Yhat = EXACT decaying solution of k''' + (sigma/4)k = 0
  (Laplace integral K_j = (1/2 pi i) int t^j e^{sigma t + t^4} dt, arb quadrature + tail bound; laplace.py),
  |d| <= rn1 = 4.37e-13 (Lyapunov-Perron for the CUBIC remainder only, step-1 constants; tail2.py, L = 0.3014),
  q = S + 4kk' + 4J, J = eps^2 Jhat + J', Jhat = cn^2 int_8^oo K1^2 = 0.631426835(1) (double contour integral +
  explicit remainder; jhat_rig.py), |J'| <= 3.6e-16.  CAPD input gen_input2.py (eps-column carries dq/deps).
  Results (s1proof unchanged):  M(1.45e-4) in [0.04339455686, 0.04339455704] (was width 1.7e-5),
     eps* in [1.5346198892965e-4, 1.5346198899299e-4]  (width 6.3e-14, rel 4e-10; was 8.3e-9)   [bracket2.json]
     tau_tip in [10.249077650575, 10.249077651520] -> sigma_tip in [-2.2490776515195, -2.2490776505754] (was +-5.4e-5),
     D in [2.50080627654, 2.50080627668], k(0) in [1.53982063762, 1.53982063822]; geometry.py: ALL CHECKS PASSED.
     |sigma_tip - STIP(psi3^)| <= 5.5e-10 (psi3^ was built on an accurate float profile).
  Downstream: tight profile boxes prof_tight.txt (s1prof.cpp: 51,246 boxes, centre enclosure + hull; centre widths
  <= 2.7e-9).  psi3res3.py (second-order Taylor form of the residual per box, profile derivatives from the ODE):
     ||(L - lamt) psi3^||_rho (half line) <= 1.2772e-6  (C 5.7e-7, B 1.14e-6 [ray + chunk-bound profile terms], T 0)
     -- was 3.8051e-3.  Floor now = psi3^'s own fit (|R(xc)| <= 2.1e-6 at piece interfaces) and the stage-B chunk bounds.
  lyap_psi3b.json: g = 4.1e-5 (was 0.122), G0 = 4.072 (eps1 ||phi3||^2 dominates), omega' max 0.02909 (unchanged).
  => the r^ forcing (was 0.022|alpha| effective) is gone.  Remaining basin bottleneck: O(alpha^2) forcing
     (||chi Q2(psi3^)||^2 <= 535, i.e. norm <= 23.1; q2_ip.json) -> second-order correction h2.
- (Nov 6, night) SECOND-ORDER SLOW MANIFOLD CERTIFIED (needed: the energy method amplifies the alpha^2 forcing ~125x
  relative to the true response alpha^2 h2, which would shrink the basin to alpha0 ~ 1.6e-4; see BOOTSTRAP.md).
  * Divergence-form E0 identity (nlin_E0.py, checked numerically in nlin_E0_check.py to 1e-8):
      int Psi w F_curv(u) = int T_c(u; w),  T_c = -[Phi J^-1 w'' + B1 w' + B0 w] kappa,  Phi = Psi/alpha_g,
    T_c/rho depends on (phi', phi'', k, m, p, u, u', u'', w, w', w'') only and has degree <= 2 in (u'', w''): the E0
    estimate needs pointwise bounds on (w, w') only (no sup|w''|; the (w'')^3 terms of the expanded form cancel).
    Float sizes (nlin_E0_float.py): alpha-linear coefficients a_ij <= 54|alpha| near the tip; cubic coefficients <= 180.
  * Q2(psi3^) is a 12-term polynomial vanishing on the ray (q2_derivs.pkl).
  * h2^ (h2md.py n = 14, h2fun.py = psi3fun construction with 2 lamt, A2 = 1 - 8 lamt = 3.049): solves
    (L - 2 lamt) h2 = -Q2(psi3^), even, admissible; h2^(0) = 0.2874, ||h2^||_rho^2 in [170.74, 171.83] (half line).
    h2res.py: ||(L - 2 lamt) h2^ + Q2(psi3^)||_rho <= 1.5719e-6 (C 1.57e-6; Bt = sigma in [8,20] with the EXACT tight tail
    profile 6.6e-8 [ray + chunk bounds fail here: 1.3e-2, the Q2 ~ k psi3^2 forcing cancels against (L-L_ray) h2];
    B 3e-10; T 4e-12).
  * Third-order forcing (n3_forcing.py, arb power series in alpha at the exact point 0 + interval c5 remainder):
      ||chi (N0(alpha psi3^ + alpha^2 h2^) - alpha^2 Q2(psi3^))||_rho <= 12.34 |alpha|^3  (|alpha| <= 3e-3),
                                                                    <= 16.69 |alpha|^3  (|alpha| <= 1e-2).
    (An interval expansion point A = [-a, a] overestimates c3 240x by dependency; never use it.)
  => total w-forcing  alpha r^ + alpha^2 r^_h + [N0(u0) - alpha^2 Q2] <= 1.3e-6|alpha| + 1.6e-6 alpha^2 + 12.34|alpha|^3.
- (Nov 6, night) INNER-ENERGY CONSTANTS, first batch (all rigorous, arb on CAPD/chunk boxes):
  * alin_E0.py (alpha <= 3e-3, u0 = alpha psi3^ + alpha^2 h2^): sup |a_ij|/|alpha| of the E0 alpha-linear form
      x <= 12:   a00 56.5, a01 58.7, a02 31.2, a10 30.2, a11 39.1, a12 0.53, a20 6.0, a21 5.3, a22 6.3
      12..31.25: a1j <= 2.9, a20 16.4, a21 31.7, a22 15.2 (these grow ~ alpha psi3^'^2: smaller for alpha <= 1e-3)
    => ||w||, ||w'|| parts (~100|alpha|) exceed omega0 = 1/32 - 0.016: handled by time-integrability (alpha decays at
       ~0.256 while ||w~|| << 100|alpha|): growth factor exp(2 C int|alpha|) ~ e^{1.8} for alpha1 = 1e-3; the w'' parts
       (<= 37|alpha|) are absorbed by the theta-margin (1/9)||w~''||^2.
  * garding_cert.py: E2-level Gaarding (Mt = level-2 operator): -mu2 <= 14.343 (x<=1), 2.841 (1..3), -0.180 (3..12),
    -0.289 (12..34.25); |e1| <= 16.28, 13.86, 2.17, 3e-8; |e0| <= 10.67, 5.25, 0.46, 2e-8 (same regions).
    With eps = 3e-4 in E = ||w~||^2 + eps ||w~''||^2 the linear part decays at omega_E ~ 0.0141 (float bookkeeping).
  * Cubic terms (float coefficients; E0 divergence form): need sup|w| <~ 3e-4 and sup|w'| <~ 4e-4 near the tip
    (monomials w w''^2 (31) and w'^2 w'' (12.4) are the binding ones)  <=>  sqrt(E) <~ 1.3e-5: requires the hand-off
    deviation from the second-order slow manifold to be <~ 1e-6 in rho-norm (alpha1^3 h3 ~ 5e-7 at alpha1 = 1e-3).
- (Nov 6, night) E2-LEVEL FINDINGS.
  * nlin_E2.py (alpha <= 1e-3, |w|,|w'| <= 1e-4, rigorous): D = F(u0+w) - F(u0) - Lw = sum d_tjl w''^t w'''^j w''''^l + d_000;
      x<=12: |1 - J^-4| (w'''') <= 6.9e-3, w''' 1.2e-2, w'' 4.2e-2, w''^2 4.75, w''w''' 0.136, w''^3 3.05; |d_000| <= 0.54|w| + 0.62|w'|
      12..34.25: 4.3e-3, 6.9e-4, 2.9e-5, 0.014, 0.46, 3.06; d_000: 1.3e-6|w| + 5.0e-3|w'|.
    Pairing: <(chi D)'', w~''>_rho = <chi D, w~'''' - 2 phi' w~''' + (phi'^2 - phi'') w~''>  (only ||chi D|| needed).
  * NEW OBSTRUCTION: the modulation forces E2 linearly:  a <(chi psi1)'', w~''> = a <D_rho f1, w~>,
    D_rho g = g'''' - 2 phi' g''' + (phi'^2 - phi'') g'',  ||D_rho f1|| = 11.8, ||D_rho f2|| <= 6.5 (float, via the ODE),
    |(a, b)| <= Minv Lfn ||w~|| = 31.2 ||w~||  =>  eps * 421 ||w~||^2 must stay below ~ omega0/3  =>  eps <~ 1e-5.
    Consequence: ||w~''|| <= sqrt(E/eps) ~ 316 sqrt(E): pointwise bounds are lossy.  Remedy: handle the E0 cubic terms by
    weighted Gagliardo-Nirenberg + the E0 dissipation (superlinear in E), use the pointwise bounds only to keep the
    jet-box coefficients near their u0 values (|w|, |w'| <~ 1e-2).
- (Nov 6, night) dmod_norms.py (rigorous): ||D_rho f1|| <= 11.854, ||D_rho f2|| <= 8.693 (|A cos + B sin| <= sqrt(A^2+B^2),
  any frame angle), cutoff strip negligible  =>  E2 modulation coupling <= Minv Lfn sqrt(...) = 458.4 ||w~||^2
  =>  eps <= ~1.07e-5 for a loss <= omega0/3.
  Interpolation instead of sup bounds for the E0 cubic terms (e.g. w'^2 w'', coefficient 12.4 at the tip):
    |int_{x<=2} rho c w'^2 w''| <= c (2.1)^1/2 (phi' ||w'||^2 + 2 ||w'|| ||w''||)^1/2 ||w'|| ||w''||  <~ 2c E^1/2 eps^-1/8 ||w''||^2,
  so the requirement relaxes from sqrt(E) <~ 4e-6 (pointwise route) to sqrt(E) <~ 1e-4 near the tip.
  NEXT: rigorous regional cubic coefficients (E0, jet boxes); inner closure in float with all certified constants;
        intermediate energy (W_I) nonlinear terms + zeta_I commutator; outer lemma; assembly.
- (Nov 6, night) cubic_E0.py (rigorous; alpha <= 1e-3, |w|,|w'| <= 1e-2 in the jet box; R3 <= (1/2) sum C |w|^a|w'|^b|w''|^c):
    x<=2:     w w'^2 140, w^2 w' 109, w^2 w'' 100, w^3 42, w w' w'' 18.7, w w''^2 17.5, w'^2 w'' 7.9, w'^3 2.6;
              quartic+: w^3 w'' 146, w^2 w''^2 117, w^3 w''^2 136, w^2 w' w'' 95, w w'^2 w'' 80, w w'^2 w''^2 33
    2<x<=12:  w^2 w' 29.5, w^2 w'' 13.5, w w'^2 17.5, w w' w'' 14.2, w w''^2 8.2, w'^2 w'' 3.6 (quartic <= 45)
    12..31.25: all cubic <= 0.9 except w'^2 w''^2 8.2, w'^3 w'' 8.6, w w'^2 w'' 4.5, w w' w''^2 4.3 (quartic).
  (BUG FIXED on the way: python-flint ball**n is NaN for balls containing 0 and max(0, nan) silently drops it; all
   scripts now multiply explicitly and assert non-NaN.)
- (Nov 7) FLOAT INNER CLOSURE (inner_closure_float.py; certified constants, far-field pointwise bounds modelled as
  CFAR sqrt(E), CFAR = 300 placeholder; eps = 1e-5; worst case over the shape (||w~||, ||w~''||, ||w~''''||)):
    alpha = 0: decay rate 0.021 (linear 0.026 - E2 coupling 0.005); cubic terms cap sqrt(E) at ~7e-6 (rate 0 at E = 1e-10).
    alpha-linear growth 0.46 (alpha/1e-3) (form 0.19 + modulation 0.27), time-integrable.
    Comparison ODE: alpha1 = 1e-3, sqrt(E1) = 1e-6 -> max sqrt(E) 5.0e-6 (margin 1.4x), decays afterwards;
                    alpha1 = 5e-4, sqrt(E1) = 1e-6 -> max sqrt(E) 1.9e-6 (margin 3.6x).
  agent_sob/ (sonnet agent, rigorous): c1 = 0.5804, Phi1(X) = 1.083, K = 1/min rho: [0,1] 1.488, [1,2] 2.078, [2,4] 5.222,
    [4,8] 56.92, [8,12] 992.4, [12,20] 8.0e5, [20,34.25] 1.2e12; W-weight (R = 40/60/80): c1W 0.645/0.848/1.029,
    Psi1 1.140/1.305/1.437, K_W [12,20] 6.4e7, [20,30] 7.3e8, [30,40] 9.5e11; chi commutator
    [L,chi]w = -4chi'w''' - 6chi''w'' - (4chi''' + 2k^2chi')w' - (chi'''' + k^2chi'' + (5km+q/4)chi')w,
    C3..C0 = 2.917, 5.010, 7.779, 16.21; rho_max(strip) 2.02e-11, len 3.
  alin_E0_spec.py (rigorous, alpha <= 1e-3): spectral radius of the symmetrised (w, w') block / |alpha|:
    x<=2 59.3 (sum of |a_ij| was 139.3), 2..12 2.37, far 0.61; |a02+a20| 35.9/2.77/3.91, |a12+a21| 5.6/2.36/7.90,
    |a22| 6.2/1.41/3.68.  (python-flint: abs() of a 0-straddling ball is a 0-centred ball, NOT >= 0 -- take
    arb(abs(z).upper()) before sqrt.)
- (Nov 7, evening) MODULATION FORMS CERTIFIED (agent_mod/, rigorous Gram 8x8 on 51k boxes, Taylor models; on {f1,f2}^perp):
    alpha-level  sup|Q1|/||w||^2 <= 13.549  (was 268);   E2-level  sup Q2/||w||^2 <= 43.327 (was 458.4), inf in [-0.189,-0.184].
  (my float check masked x < 0.3 and was wrong -- the agent caught it; cos(theta) = St C - Ct S in the Lyapunov frame.)
  Float inner closure with these: alpha1 = 1e-3, sqrt(E1) = 1e-6 -> max sqrt(E) 1.83e-6 (was 5.0e-6); alpha1 = 2e-3 -> 5.8e-6.
- (Nov 7, evening) agent_lint (L-int variants): ORIGINAL lint_cert.py does NOT reproduce 0.2435 (fixed grid, smoothstep
  dependency near B2 -> -1.40); lint_cert_p.py (adaptive bisection, NaN guards) is the valid certificate. Variants table in
  agent_lint/LINT_VARIANTS.md; new runs XA = 10.25: P = 2 0.1552, 1.75 0.1454, 1.5 0.1241 (blend 0.336/0.279/0.221).
  w-level variant (lint_w/, operator L): P = 3 fails at XA = 12 (-0.0065), ok from XA = 13 (0.069), 14 (0.122), 16 (0.191).
- (Nov 7, evening) agent_int (T1-T4, P = 3, R = 40/60, alpha <= 1e-3, jet box 1e-2): C3W = 8.83e-3; W-residuals of psi3^, h2^
  2.6e-10, 2.4e-10; far coefficients |1 - J^-4| <= 1.85e-2, w''' 1.5e-3, w''^2 0.015, w''w''' 0.98, w''^3 3.26, A0 2e-6,
  A1 8.6e-6 (nlin_E2's A1 = 5e-3 was a q/4 dependency artefact); sup|u0'| 0.085 at R = 60 (thin chart margin);
  level-2 Gaarding -mu2 <= -0.744.  (P = 3 design is superseded; redo for P = 2 / 1.5.)
- (Nov 7, evening) FAR-FIELD COUPLING DESIGN decided (BOOTSTRAP.md, last section): v-level E_I with P = 2, ramp [10.25, 16.25].
- (Nov 7, evening) agent_int/v2 (new design, P = 2 / 1.5, x in [10.25, R], sigma >= 8.0009): C3W = 0.1823 / 0.6885 (w-level
  forcing ||N0(u0) - a^2 Q2||_W; the v-equation only needs it via integration by parts); W-residuals psi3^ 1.09e-7 (loose on
  [10.25, 12], ray + chunk bounds) / h2^ 6.5e-9; level-2 Gaarding -mu2 <= -0.6555 (= mu1 + 1/2), |e1|,|e0| <= 5.7e-6;
  ||zeta psi3^'|| 0.26, ||zeta h2^'|| 1.8, ||zeta A1(u0)'||/a 0.273, ||zeta B1(u0)'||/a 0.016, ||zeta psi1'|| 1.7e-7, ||zeta psi2'|| 1.7e-8,
  ||zeta c0'|| 1e-11; W(57) 1.3e-18, W(60) 2.8e-20 (P = 2).  V3 far coefficients as T3; on [10.25, 16.25]: w'' 1.1e-4, A0 1.7e-5, A1 7.1e-5.
- (Nov 7, evening) COUPLED FLOAT CLOSURE inner_closure_float2.py (Lyap = E + mu E_I, worst case over the E/E_I split and both shapes;
  coupling constants still [DES] floats, everything else certified): mu = 0.08, epsI = 3e-5, alpha-linear form by regional spectral
  radii (ALIN=spec): alpha1 = 1e-3, ||w~||_1 = 1e-6, ||v||_1 = 3e-6 -> max sqrt(Lyap) 1.94e-6; alpha1 = 2e-3 -> 4.65e-6; 3e-3 fails.
  Nonlinear cap sqrt(Lyap) ~ 7e-6 (inner cubic terms for smooth shapes; level-2 intermediate cubic terms for rough shapes).
  mu * epsI must stay <~ 2e-6 (level-2 ramp flux vs the inner eps-dissipation).
- (Nov 7, evening) agent_coup (rigorous coupling constants, both P; identity -<[L,chi]w,chi w> = int A w^2 + B w'^2 verified symbolically):
  P = 2: sB 3.5955e-5, sA 1.3372e-4, IA 9.60e-12, K_rho(16.25) 3.0208e4, Phi(34.25) 9.8547e6, a_f 6.956e-2, b_f 3.936e-2,
  level-2 ramp c3..c0 = 0.442/0.526/0.424/0.425, level-2 strip d3..d0 = 6.84e-3/1.54e-2/2.36e-2/4.56e-2;
  mu-window [2.906e-3, 0.1135] (ratio >= 39).  P = 1.5: window [1.356e-4, 7.48e-3] (>= 55).  Outer ramp sup|A_zo| 3e-18 (R = 60).
  Coupled float closure with all these: alpha1 = 1e-3 -> max sqrt(Lyap) 1.95e-6; 2e-3 -> 4.68e-6.
  STILL FLOAT/MODELLED in inner_closure_float2: the far Sobolev constants for the P = 2 weight (K_W, Psi1, c1W are P = 3 / design
  values) and the bookkeeping of the E_I inequality itself (to be written out analytically).
- (Nov 7, late) agent_int/v3: the certified psi3^, h2^ are only C^3 (u'''' jumps 1e-7..4e-6 at the 20 piece interfaces), so derivative-level
  quantities (level-2 intermediate energy: Phi_F', w^(5)) need C^4 versions; c4fun.py builds psi3~, h2~ (degree-9 Hermite corrections,
  jets 0-4 matched, jumps < 1e-55, ~1e-7 relative change). For these: C3W' (derivative-level forcing, with zeta) 2.60e-3; derivative
  residual norms 2.7e-8 / 1.16e-7; level-0 W-residuals 1.24e-8 / 4.0e-8; sup|1 - q psi'| <= 74.4 (R = 60); sup|(W cos)'/W| <= 1.31;
  ramp sup q zeta' W/rho = 0.436; ||zeta A1(u0)'''||/a 7.2e-5, B1 2.1e-6, ||zeta psi3~'''|| 2.6e-5, ||zeta h2~'''|| 1.1e-2.
  TODO: switch the whole pipeline to psi3~, h2~ (rerun psi3res3/h2res/inner products; perturbation bounds ~1e-7 relative suffice).
- (Nov 7, late) MODULATION LEMMA (agent_mod/lemma, rigorous): leading |(a,b)| <= 3.708 A (was 31.2 A: w~ perp f1, f2 projects L*f);
  with nonlinear feedback c_A = 34.7 / 8.98 / 7.57 for jet box 1e-2 / 1e-3 / 1e-4 (driven by A0 = 0.54 of nlin_E2 -- sup over alpha <= 1e-3,
  likely partly a dependency artefact); beta: (0.00244 + 63.5|alpha| + 13.4 W1) A + 0.626 ||chi Phi_F|| (+ nonlinear); G33 >= 0.031107,
  ||G^-1||_inf <= 259, Neumann factor 0.080.  NOTE: the old closure's beta forcing 0.63 Fn/G33 double-counted 1/G33 (fixed: 0.626 Fn).
  NEW DANGER: the alpha-part of beta (via E_13 = -alpha <chi A1(psi3^), f3>) times n3 = 7.52 is an alpha-linear form in E0 of crude size
  300-480 |alpha| (float closure: <= ~150 needed at alpha1 = 1e-3; 2e-3 fails for >= 84). Exact rank-form computation requested (agent_mod/bloop).
- (Nov 7, late) HAND-OFF NUMERICS RESOLVED (agent_handoff2/oddtest): the "unstable mode growing like k" was an ODD frame artefact -- the
  numerical curve is symmetric to 1e-9 lambda but about the axis y = +3.477e-7 (tip offset), so the y = 0 frame sees a perpendicular
  translation (cos theta1, eigenvalue 1/4, odd class). With the axis through the tip: r2/|alpha|^3 = 32-35 (k = 400, 800). In the even class
  the remainder is ~32 |alpha|^3: 2.6e-7 at |alpha| = 2e-3 (k ~ 1650), 3.2e-8 at 1e-3 (k ~ 3300).  Design study launched: arms along the
  S1 asymptotic ray (reduce alpha*k_max ~ -3.3, i.e. the arm curvature near the singular point) to hand off at much smaller k_max.
- (Nov 9) C^4 SWITCH (agent_int/c4switch, rigorous): perturbations psi3~ - psi3^, h2~ - h2^ (zero for x >= 28): rho-norms of d^(j), j = 0..4:
  psi3 1.9e-10, 1.1e-9, 8.9e-9, 9.5e-8, 1.2e-6; h2 7.3e-10, 3.0e-9, 2.2e-8, 2.3e-7, 3.0e-6 (W-norms much smaller; sups on [0,34.25] <= 3.9e-6 / 1.5e-5).
  New rho-residuals: ||(L - lamt) psi3~|| <= 7.93e-7 (was 1.277e-6), ||(L - 2 lamt) h2~ + Q2(psi3~)|| <= 2.769e-6 (was 1.572e-6);
  third-order forcing C3 = 12.340 (|alpha| <= 3e-3), 12.104 (|alpha| <= 1e-3).  Inner products / Gram / Lyapunov constants: carry over
  with the S1 error terms (relative changes <= ~1e-6).  => adopt u0 = alpha psi3~ + alpha^2 h2~ everywhere.
- (Nov 9) CIGAR TEST (step4/cigar/cigar.py, float): stadium aspect 25:1 rounds to a circle without self-crossing; aspect 120:1 (r = 0.05,
  N >= 8192 enforced -- default coarsening under-resolved it) retracts from both ends, rims meet, waist self-crosses briefly (two axis
  crossings at x = +-0.38, t ~ 0.13), heals by t ~ 0.18, then a circle; max kappa*r only decreases (1.05 -> 0.16). Valley behind a
  retracting end deepens slowly (min half-width/r 0.94 -> 0.55 over a 340x time increase). => the slot + chamber (pre-built pocket) is what
  step 4 uses; long arms are not essential (keyhole slot in a disk is a candidate starting curve).
- (Nov 9) BETA-LOOP RESOLVED (agent_mod/bloop, rigorous): e13 = <gA, f3> in [-0.02635, -0.02499] (crude 0.347), e23 = -<gB, f3> in
  [0.003445, 0.003617], H33 = 31.0 +- 0.7, E33/alpha in [0.0302, 0.0338]; the beta-loop alpha-form largely CANCELS Q1:
  combined Qfull/alpha on {f1,f2,f3}^perp in [-6.19, 1.90] (Q1 alone [-13.47, 11.02], loop alone [-11.16, 9.15]).  So the alpha-level
  modulation coefficient is 6.19 (was 13.55 + loop).  beta's alpha-part is ~3 |alpha| A (not 63.5).  Float closure: alpha1 = 1e-3 -> max
  sqrt(Lyap) 1.86e-6; 2e-3 -> 4.14e-6; 3e-3 still fails (alpha^3 forcing + alpha-linear growth).  Open: alpha = 0 loop form (item 3; the
  bound eps1 n3/G33 = 0.01607 stands).
- (Nov 9) LEMMA O (writeup/lemmaO/, v1 + v2, energy method, explicit double-precision constants): proof structure fine, but the
  thresholds are far too small: binding quantity pointwise h_xx; t* ~ 1e-7 (X = 2.4) .. 1.6e-6 (X = 5) for M ~ 0.03-0.3, while the
  hand-over at y_h = 30 needs t ~ 4e-5 (cylinder half-width < r). Misses by 100-500x (6-25x at X = 5). Pairing forcing with energy (more
  regularity) gains <= 28%; nested cutoffs do not help. PDE sanity check: true deviations are 1e-5..5e-4 of the bounds (method is crude,
  not the statement). => crude Lemma O is usable only where the time ratio is <~ 1e-7, i.e. rescaled distance y >~ 2.4 / (1e-7)^(1/4) ~ 135
  (and cylinders small enough that the bent cone's turning across them is < 1/4).
  OUTER-REGION PLAN B (cheapest analytically): extend E_I's polynomial region with SCALE-ADAPTED level weights (level-2 weight eps_I (1 +
  (y/60)^4) W, so pointwise control is scale-invariant ~ y^P sqrt(E_I)), P: 2 -> 1 beyond ~40 (rate (2P-1)/8 >= 1/8), out to R ~ 250 in the
  S1-ray chart (slope of u0 ~ 0.37 at y = 250 if |alpha| <= 1e-3), and use Lemma O v2 only beyond (time ratio (2.4/250)^4 ~ 8.5e-9 << 1e-7;
  hand-over M ~ 0.12; far curve with small cylinders, ratios <= 1e-9). Cost: hand-off at |alpha1| <= 1e-3 (k* ~ 3300 for the current clip;
  depends on the shape study) and a recertification of the far-field constants for the new weight. Plan C (keep alpha1 = 2e-3): chart over the
  slow-manifold curve (bent cone) beyond ~60 -- new derivation. Decision pending the shape study.
- (Nov 9) SHAPE STUDY (step4/design/): within the clip family alpha*k ~ 78/arm, but the late length L grows with the arm, so uniform-mesh
  N(k*) ~ 1.3-2.6e5 for every arm length 24-48; best: arm 24, chamber 0.4, W 0.85 (|alpha k| = 2.86, k* ~ 1540 at |alpha| = 2e-3);
  arm 48: |alpha k| = 1.84. Disk-with-slot and arms <= 5 do NOT singularize (relax to a circle) -- the long arms matter. S1 ray at 106.7 deg
  from the axis. Even remainder r2/|alpha|^3 ~ 45. A tip-adapted (non-uniform) mesh is the only big lever (N cut by ~L/10; then k* enters
  only logarithmically).
  DECISION: hand-off at |alpha1| ~ 1e-3 (k* ~ 2900 for the best clip, ~1800 for arm 48), outer region by PLAN B (E_I extended to R = 260 with
  P: 2 -> 1 beyond 40 and scale-adapted level-2 weight; Lemma O beyond ~200); step 4 must use a tip-adapted mesh.
- (Nov 9) LEMMA O v3 (writeup/lemmaO/lemmaO_v3.{tex,pdf}, lemmaO_v3_consts.json, lemmaO_v3_interval.py; python-flint balls, float match 1e-6):
  near the S1 cone, K2 = K3 = 1.5M, M = 0.12: X = 3: t* >= 3.70e-7 vs needed 1.2e-7 (margin 3.08); X = 2.4: t* >= 1.224e-7 vs 5e-8 (2.45).
  Binding term: pointwise h_xx. M-sweep: targets hold for M <= 0.3. Far curve: certified t* >= 1.0e-9 for X = 2 with M <= 3.69 and X = 2.4 with M <= 4.10
  (K = 1.5M); NOT for M ~ 6-10 (short 3-56x) -> far cylinders must be chosen with M <= 3.7 (radius >= ~0.1 physical at k* = 2900 gives ratio <= 1e-9;
  larger radius where the curvature is low). Hand-off checklist H1-H9 (lemmaO_v3.tex): one graph over [-X, X] in the tangent frame; |h0_x| <= 0.12
  (cone) / 0.11 (far); |h0^(j)| <= M, j = 2..4; int (h0^(5))^2 <= M^2 (good time); ||h0_xxx||_L2(-X,X) <= 0.14 (X = 3) / 0.1365 (X = 2.4) / 2.12 (far,
  X = 2) -- NOT implied by the pointwise bounds, check on the data; a priori bounds must close at all neighbouring centers (overlap with the E_I region:
  Lemma O from ~200, E_I to 260, OK); other arm at distance >= 0.1 cylinder units.
- (Nov 9) RIGOROUS VERIFIER (verifier/, VERIFIER.md, rate_table.json, verify_log.txt, run_all.sh; R = 60 constants): with added terms (b) 9.9e-12 V0
  in I0, (c) ramp flux 1e-17, (d) strip cubic with E_I sups: CLOSES. alpha1 = 1e-3: max sqrt(Lyap) 2.014e-6 (float 1.86e-6; margin 3.7x to 7.5e-6);
  alpha1 = 2e-3: 4.874e-6 (1.54x). Float grid underestimates the max rate by up to ~2x near alpha = 0 (-0.0222 vs certified -0.0110).
  Term (a), the I0 inner-ramp coupling 0.66*(d001*D4 + d010*W3 + ...)*(V1 + 1.67 V0), does NOT close with the constants as specified
  (rate +5.6 even at alpha = 0): D4's only damping is (eps/2) D4^2, eps = 1e-5; threshold for the D4 coefficient ~2e-3 (closes) .. 5e-3 (fails).
  DIAGNOSIS: the spec used the [10.25,60] sup d001 = 1.85e-2 at the max alpha, as an alpha-independent constant. The true coefficient is the
  linear deviation (1 - J^-4) restricted to the RAMP [10.25,16.25]: ~ -4 k u0 + 2 u0'^2 + w-part, i.e. O(|alpha|) and 0 at alpha = 0
  (T3_per_region [12,20]: 2.46e-3 at |alpha| <= 1e-3 with w-box 0.01). FIX: ramp-restricted, alpha-dependent coefficients (asked of the Plan B
  agent), then rerun the verifier once with the Plan B constants. Open caveats for the final pass: (i) tail lemma (radial scaling, degree >= 2)
  argued, not machine-proved -> write it up; (ii) rounding by 5e-16 relative inflation, not directed rounding -> switch to arb or nextafter.
- (Nov 9) PLAN B CONSTANTS DONE (agent_int/planB/, planB_constants.json, NOTES_planB.md): MU1 >= 0.11515 (dip at 60; -> 1/8 from below),
  MU2 (level-2 weight) >= 0.128, c1W 0.0936, KW [20,40]..[160,260] = 2.6e6, 1.6e7, 6.4e7, 1.7e8; C3W 7.55 at alpha <= 1e-3 (0.18 at <= 3e-4);
  T3 on [250,260] large (w'''' 0.69, w''^3 341); W(260)/W(10.25) = 6.7e-5 (not 9e-9). Inner-ramp add-on [10.25,16.25]: w'''' coeff
  7.7e-4 / 8.5e-4 / 1.76e-3 at w-box 0 / 1e-3 / 1e-2 (alpha = 1e-3), fit 0.19|a| + 578 a^2 -> term (a) is fixable.
- (Nov 9) STRUCTURAL PROBLEM FOUND (supersedes the Plan B decision): energies of E_I type cannot supply the POINTWISE data the outer
  region needs. Lemma O hypotheses at hand-over y_h need scale-invariant bounds (|w''| <~ M X/y_h ~ 2e-3 at y_h = 200, higher jets
  ~1e-5..); from E_I at the closure level (E_I ~ 5e-11, eps_I = 3e-5, K_W = 1/min W ~ 1e8) one gets sup|w''| ~ 0.1-1 even at a
  Chebyshev good time. Same defect in the R = 60 design (hand-over at y_h = 30: needs |w''| <~ 0.02 and t* >= 3e-4; cylinder
  [22.5, 37.5] reaches into the rho-blend where E_I gives nothing). Any decaying weight with a FIXED outer edge forces per-point
  hand-over at fixed y_h (needs high-jet pointwise data -> fails); an EXPANDING edge (fixed in physical space) with Lemma O applied
  only once at t1 is fine at the edge but the interior far points (r < r_O) need uniform-in-time scale-invariant slope control,
  which a global weighted energy only gives marginally (K(y(tau)) grows exactly as fast as E decays, by scale invariance).
  KEY FACT: in shell variables xi = y/y_k the slope v = u' obeys v_tau = -(xi/4) V_xi - y_k^-4 (a V'' + b V'^2)'' : the far field is
  pure transport up to y^-4 corrections (forcing of the alpha-series is also y^-4 suppressed; its far-field inaccuracy is just
  transported hand-off data).
  NEW PLAN ("Lemma S"): keep the certified R = 60 E_I; add scale-invariant dyadic-shell energies (levels 0, 1 for v = w', cutoffs of
  O(1) relative width) on y >= 24 up to a physically fixed outer edge inside the Lemma O region; inflow from E_I only through
  TIME-INTEGRATED L^2 quantities (int V0^2, int V1^2 dtau); all 4th-order/nonlinear/commutator terms carry y^-4. Rough float
  estimate at sqrt(Lyap) = 2e-6: shell slope |w'| <~ 0.05, scale-invariant y|w''| <~ 0.25 -> parabolicity and T3 coefficients fine.
  Lemma O (crude, certified) used ONCE at t1 for y >= y_O ~ 150-200 with step-4 data (margin 2.3-7). alpha1 = 1e-3 can stay.
  Plan B's R = 260 constants are not needed (sunk), except the inner-ramp add-on.
- (Nov 9) LEMMA S first pass (writeup/lemmaS/lemmaS_consts.py): transport form verified symbolically (S1 tail ~1e-14 at y >= 24; 4th-order
  and b-terms carry e^-4s). BUT the u0/modulation forcing (beta(psi3' + 2 alpha h2'), (a y + b cos th) U_y) is NOT y^-4 suppressed and is
  amplified along characteristics (U ~ y^1.02) while the closure certifies only ~0.011 decay -> deviation shells cannot be uniform in time.
  FIX (sent): split. (A) deviation shells v = w' on [24,60] only (transit 3.7 tau-units -> forcing accumulates over a bounded window);
  (B) full-slope shells f = u' on [48, R(tau)] (a priori bounds only; modulation = transport-speed terms costing int|a|, int|b|).
  Note: the weak certified decay (0.011) makes the inflow integrals ~45x Lyap1/mu; hand-off level may need to drop to ~5e-7.
- (Nov 9) LEMMA S split version (writeup/lemmaS/: lemmaS.tex/pdf, lemmaS_consts.json, lemmaS_interval.py, lemmaS_symbolic.py):
  transport reduction residual 0 (y^-4 <= 3.0e-6 at 24). (A) inner ramp at y in [19,24] (x in [21.25,26.25]).
  S0_v^2 <= E0/l + 2 sqrt(E0 E1), S1_v^2 <= E1/l + 2 sqrt(E1 E2), l = ln 2.5; inflow E_j^in = (1/4) K_j J_j, J_j = int V_j^2 dtau,
  K = 1.72e5, 8.56e7, 4.29e10 (W-weighted; unweighted ramp versions 0.478, 235, 1.17e5). Beta forcing harmless (~1e-5 in v; dominant
  63.5 alpha A). (B) needs moving windows of fixed s-length; cost exp(C 48^-4 + c int|a| + int|b|/48).
  With MY reference J's: S0_v = 0.051 (sqrt Lyap1 = 1e-6) / 0.026 (5e-7) OK; S1_v floor 0.43 / 0.21. NOTE my reference J2 = 2e-6 was far
  too pessimistic: J2 = int V2^2 is the level-0 dissipation, ~ V0^2(tau1) + coupling ~ 1e-10..1e-9 -> E2^in ~ O(1), S1_v ~ 0.36 (1e-6) /
  0.18 (5e-7). Waiting for the verifier's rigorous J0, J1, J2. Caveat: rate constants C_j hand-derived, Prop. 1 only sketched (conclusions
  insensitive: e^{C eps} <= 1.35 for C <= 1e5) -> must be written properly later.
- (Nov 9) VERIFIER FOLLOW-UP (VERIFIER.md, verify_log_followup.txt, rate_table_followup.json, integ.py): term (a) with certified ramp
  coefficients closes in every case (rate insensitive; ~5-10x room). Hand-off scan (alpha1, sqrt Lyap1 -> peak, margin to 7.5e-6,
  max S_ramp): (1e-3, 1e-6) 1.57e-6, 4.8x, 2.2e-2 [box violated]; (1e-3, 5e-7) 8.5e-7, 8.8x, 1.2e-2 [violated]; (1e-3, 3.5e-7) 6.4e-7,
  11.7x, 9.0e-3 OK; (1e-3, 2.5e-7) 5.0e-7, 15x, 7.0e-3 OK; (5e-4, 5e-7) 5.9e-7, 12.7x, 8.3e-3 OK. J0, J1, J2 at (1e-3, 3.5e-7):
  7.21e-11, 1.01e-10, 2.24e-11 (+ linear growth ~2e-16/tau from constant placeholders (b),(c): being recast as quadratic terms).
  Lemma S with these: S0_v = 0.0129, S1_v = 0.219 (sup|w''|(24) <= 9.1e-3).
  DECISION: hand-off target alpha1 = 1e-3, sqrt Lyap1 = 3.5e-7 (backup 2.5e-7). T3 on [16.25,60] being recertified with slope box 0.02.
- (Nov 9) T3 RECERTIFIED (agent_int/planB/r60_T3_boxw1_0.02.json, r60_t3.py): R = 60 regions, alpha <= 1e-3, |w'| <= 0.02, |w| <= 0.02
  / 0.5: max d001 2.27e-2, d010 1.67e-3, d200 0.0153, d110 1.086, d300 3.314 (d100 2.54e-5 far). |w|-dependence nil; exact polynomial in
  w'', w''', w''''.
- (Nov 9) VERIFIER FOLLOW-UP 2 (verify_log_followup2.txt, rate_table_abcdQ_L*_a0.001.json): (b) recast as mu 9.9e-12 S_all V0, (c) as
  CFLUX (V0^2+V1^2+V2^2) [CFLUX an INPUT; safe up to 3e-2; certification asked of the Lemma S agent]. New T3 in. Both hand-offs close:
  sqrt Lyap1 = 3.5e-7 -> peak 6.39e-7 (11.7x), S_ramp 9.0e-3, E_I-based far sup|w'| 1.40e-2 (<= 0.02 box), J = 7.21e-11, 1.013e-10,
  2.24e-11 (finite, to infinity); 2.5e-7 -> 4.95e-7 (15.1x), 6.98e-3, 1.08e-2, J = 4.19e-11, 5.89e-11, 1.30e-11. Lyap -> 0 (no floor).
  STEP 3 REMAINING: CFLUX cert; Lemma S Prop. 1 write-up; final rigor pass (true outward rounding, tail lemma tex, B&B leaf semantics);
  assembly write-up.
- (Nov 9) LEMMA S pass 3 (writeup/lemmaS/, lemmaS_assembled.json, lemmaS_cflux.json, prop1 scripts): Prop. 1 identities machine-verified
  (sympy, levels 0-2), core rate constants certified (arb; eps_in*G <= 0.04). (B) works: moving cosh windows, no cutoffs, Theta_B = 1.007,
  S0_f <= 0.60, S1_f <= 0.62 (int A <= 1.2e-4, int|a| <= 8.8e-4). CFLUX: no multiplicative bound in cut norms; ABSOLUTE additive
  bounds |K0| <= 6e-18, |K2| <= 1.5e-17 (W(57) <= 1.4e-18 thanks to the rho-blend), 8e-21 / 2e-20 with Lemma S sups -> use additive form.
  (A) FAILS as written: inner-ramp [19,24] (0.23 in s) commutators with S^4 cutoff have derivative sups 1e3-1e9, Young squares them ->
  S0_v ~ 1.5, S1_v ~ 33. Without ramp terms: S0_v = 0.026 (with gamma = 1 damping weight), S1_v = 0.54 (my 0.013/0.22 omitted the
  ramp weight and conversion). Suggested fixes: Cauchy-Schwarz with ramp data (needs J3 = int V3^2, i.e. E_I level-2 dissipation),
  wider/smoother ramp (y in [15,30], 0.69 in s). Also: Prop. 1 uses a jet box with V-derivatives up to 7th order (source unclear ->
  possible circularity). Waiting for the independent (Fable) review before redesigning.
- (Nov 9) INDEPENDENT (FABLE) REVIEW -> STEP3_REVIEW.md. Step 3 NOT nearly done. Ranked: (1) Lemma S(A) does not close [fatal as
  stated]; (2) no pointwise bounds on x in [34.25,60] (verifier used K_W of [30,34.25] for the whole far region; nothing gives sup|w'''|)
  [serious]; (3) verifier's |(a,b)| <= 7.57A needs jet box 1e-4 (actual 1e-2 -> 34.7A) and drops the W1,B,W3,D4,sup|w''| terms of the
  modulation lemma; beta's D4 term not absorbable -> fix by integrating <chi N, f_j> by parts [serious]; (4) term (c) has no multiplicative
  form -> absolute forcing -> floor, so Lambda -> 0 not proved [serious]; (5) Lemma O near the cone: time ratio used half-width r (reaches
  the tip); bent-cone slope 0.238 vs limit 0.2437 at y_O = 150, y_O = 200 fails; H6 still points to the Plan B energy [serious];
  (6) certified tail decay is 2.04e-3 (Lambda), not 0.011 [minor]; (7) jet box 1e-2 violated on [16.25,32] (1.40e-2) for the inner cubic /
  modulation coefficients [minor]. Solid: transport structure and the (A)/(B) split idea, inner closure architecture and margins, Lemma O v3
  as a statement, modulation linear algebra, conclusion logic (conditional). Fix plan requested from the reviewer -> STEP3_FIXPLAN.md.
- (Nov 9) VERIFIER FINAL PASS (VERIFIER.md "Final pass", TAIL_LEMMA.tex/pdf, tail_cert.py, leaf_audit.py, verify_log_final*.txt):
  rounding fully outward (nextafter + 200-bit arb constants); production numbers unchanged (peak 6.392908e-7 / 4.953730e-7); tail lemma
  structural facts machine-certified (degree tracking), decay only "at first touch"; leaf acceptance semantics proved (U rigorous for
  any tol). CAVEAT: sharp J's valid only if Ly stays in each rate cell (unproved); rigorous J's (Ly <= barrier only) at 3.5e-7:
  int V0^2 <= 2.5e-9, int V1^2 <= 4.6e-7, int V2^2 <= 8.4e-5 (1/eps_I loss). Sharp unconditional J's need an alpha-additive treatment
  of linear forcing in the B&B. Forwarded to the fix-plan designer.
- (Nov 9, late) FIX PLAN ADOPTED (STEP3_FIXPLAN.md, by the reviewer): hand-off alpha1 = 5e-4, sqrt Lyap1 <= 2e-7 (fallback 3.5e-7); theorem
  version W (trapping with a floor; convergence = later small-to-medium upgrade C); far field option (ii) modified (shells (A) core
  y in [28,60], ramp [14,28], divergence form, Cauchy-Schwarz ramp commutators, only V, V_s pointwise); F2 pointwise data on [34.25,60]
  via S0_v, S1_v coefficients + identity removing sup|w'''|; F3 modulation by IBP (no D4/W3; box |w|,|w'| <= 0.02); F4 term (c) absolute
  forcing from shell energies; F5 Lemma O geometry y_O = 200, half-width 0.75 r, R1 = 175, (B) windows fixed in physical space.
  Dispatched: verifier agent V1 (unconditional J0,J1,J2,J4) + F2 + F3 insertion + F4 + production runs; Lemma S agent F1 + F4 constants
  + F5; new agent (sonnet) F3 modulation IBP lemma -> agent_mod/ibp/mod_ibp.json. Estimated ~60% of step 3 done.
- (Nov 9, late) LEMMA S F1/F4/F5 (writeup/lemmaS/lemmaS_F1.json, _F4.json, _F5_geom.json; design J's, alpha1 = 5e-4, gamma = 1):
  S0_v = 0.0152-0.0154 (box 0.02 OK), S1_v = 1.27-1.35 (limit 3 OK), E2 = 125-155 (plan ~100). Quasilinear level 2 needs the eps(14)
  split of the core rate at ramp points (G_22 ~ 2e5-4e6) -> load-bearing, to be stated as a lemma. F4: K0 <= 1.8e-20, K2 <= 5.5e-19
  (N3 term still placeholder). F5 geometry (true coefficient 1.3056): slope 0.112/0.2437, curvature 0.037/0.12, time margin 2.3.
  Next: (B) prop + region table, rigorous N3, rerun with unconditional J's from V1.
- (Nov 9, late) LEMMA S CORRECTION: quasilinear level-2 runs of (A) were INVALID (gamma'_2 = -0.09 / -6.4); only the linear model
  closes level 2 (E2 ~ 125, S1_v ~ 1.27). Levels 0/1 valid (S0_v = 0.0153). Core level 2 needs G_22 <~ 1.2e5, honest box gives 7.7e5
  (cubic cross terms via Y_3 ~ 1e3, driven by M0, M1). (B) closes at level 2 with the honest box (y_E = 48-50, eps_B = 50^-4), but its
  entry data use (A)'s E2. N3 now rigorous (negligible). Region table and O<->B bootstrap statements written (not proved). run_F1.py CLI
  ready. Asked: GN/L^4 treatment of the cubic terms (or gamma = 1.5 / good-unknown energy), else a list of downstream level-2 uses.
- (Nov 9, late) LEVEL-2 ATTACK (lemmaS_consts.json "level2_attack"): (A) core level 2 INVALID for gamma = 1..3 (G_22 7.7e5-1.1e6);
  even with all jets >= 2 removed (best case for GN) G_22 = 2.5e5 at M1 = 1.4 (needs M1 <~ 1.2; bootstrap needs >= 1.32). Level 2 is
  stable for y >= 44.6 (honest box) / 33.5 (no high jets) since eps = y^-4. DECISION: nested design: levels 0-1 on (A) core [28,60];
  level 2 only from y >= 45 (ramp [45,50], inflow from (A) level-1 dissipation) into (B); E_I's S1_v coefficient on [34.25,60]
  replaced by a weighted GN bound in S0_v and E_I norms (verifier). H6 K2/K3: windows present at tau1 from hand-off data; later
  windows from the nested level 2.
- (Nov 10) NESTED LEVEL 2 (lemmaS_nested.py, lemmaS_GN.py): closes THINLY (gamma'_2 = 0.0616, E2 = 151, S1_v[45,60] = 1.340,
  M1 bootstrap 1.385 <= 1.4: 0.7% margin). CIRCULAR: sup|z1| on [28,45] assumed by levels 0-1, never derived. Jet boxes M2-M4 = 3
  assumed in (A),(B). GN bound for E_I needs pointwise ||w'''||_L2 -> useless in the rho-blend region. Asked the designer (Fable) for a
  decisive circularity-free far-field design (options: sqrt Lyap1 = 1e-7; shift the rho-blend to ~48; parabolic-smoothing lemma for
  jets; effective-time formulation) with a dependency graph and >= 2x margins -> STEP3_FIXPLAN2.md.
- (Nov 10) FIXPLAN2 (designer): energy far field is the right tool; two implementation errors: (1) ramp commutators taken at too high
  order (identity int Theta w d^2(a w'') = int a Theta w''^2 + 2 int a Theta'' w w'' - int (a Theta'')' w w' - int a_s Theta' w'^2 has no
  Theta' w' w'' term -> only order-j ramp data, E2 ~150 -> ~1-5, no J3/J4); (2) M1 must be two-zone: ramp [14,28] read off E_I pointwise
  inside the closure (~1.2, derived), core [28,60] bootstrap box 0.75 (re-derived ~0.30). M2-M4 boxes not needed (divergence form, L^4 of
  v_ss). Design: alpha1 = 5e-4, sqrt Lyap1 <= 1e-7, single shell system levels 0-2 on core [28,60] + ramp [14,28], J0-J2 only; (B), Lemma O
  as before. Margins >= 2.1x (float, scaled). Hand-off now also needs E2(tau1) <= 0.2 on y in [14,60] (step 4). Dispatched T1-T3 (Lemma S
  agent), T4 (verifier: unconditional J's + certified pointwise E_I sups at (5e-4, 1e-7)). Reserve: parabolic-smoothing lemma.
- (Nov 10) MODULATION LEMMA BY IBP (agent_mod/ibp/: MOD_IBP.md, mod_ibp.json, ibp_*.py): residual-0 symbolic IBP; no D4/W3; w''^3
  cancels exactly; float check 1e-11. Certified (outA5 hull boxes, C^4 jets, box |w|,|w'| <= 0.02, A, W1 <= 1e-3), Q = (A+W1+B)^2:
  |beta| <= (0.00244 + 66.5|alpha|)A + 5.53 Q + 0.657 ||chi Phi_F|| + strip;  |(a,b)| <= (3.708 + 622.4|alpha|)A + 51.3 Q + 3.84 ||chi
  Phi_F|| + strip (strip coeff < 1e-8/S). q_ab includes the E-matrix W1*A feedback dropped by the verifier before. Caveat: simplex B&B
  for q in doubles (+1e-12) -> being redone in intervals. Forwarded to the verifier.
- (Nov 10) Modulation q's now interval-certified (arb simplex B&B, simplex.py): q_beta 5.526, q_ab 51.34 (alpha <= 5e-4); 7.623 / 57.50 (1e-3). Unchanged.
- (Nov 10) VERIFIER V1 partial: unconditional J's at alpha1 = 5e-4 (old F3 placeholders): sqrt Lyap1 = 2e-7: J0..J2, J4 = 1.95e-11,
  2.73e-11, 6.04e-12, 2.01e-7; 3.5e-7: 5.07e-11, 7.13e-11, 1.57e-11, 5.25e-7 (~1.4x conditional). Certified modulation constants swapped
  in; F2 inputs S0V = 0.0153, S1V = 3.0; F4 absolute flux placeholders 6.1e-18 / 1.55e-17. Production marches at (5e-4, 1e-7) and
  (5e-4, 2e-7) running (slow: ~35 s/cell under load).
- (Nov 10) LEMMA S per FIXPLAN2 (lemmaS_F2.py/json, lemmaS_T1.py): ALL LEVELS CLOSE (design J's /3.5, sqrt Lyap 1.45e-7, E(tau1) <=
  (2e-7, 1e-4, 0.2), core box (0.11, 0.75), ramp (0.05, 2.5)): E = 4.2e-7, 2.7e-4, 1.15; S0_v = 0.0069 (2.9x), core M0 0.050 (2.2x),
  core M1 0.32 (2.3x), G_22 core 3.9e4 (3.1x); gamma' = 0.246, 0.247, 0.187. T1 identity: sign of (a Theta'')' w w' is +, checked
  (Euler operator). Ramp commutators by repeated exact IBP (Young too lossy). Open: v_ss^2 L^4 term not yet in (GN-idealised); ramp-foot
  cut t >= 1e-3; inputs are design values; F4 and (B) entry to redo. Asked for these; then rerun with the verifier's J's + ramp M's.
- (Nov 10) LEMMA S: ramp-foot cut removed (harmless). L^4 identity int Theta z2^4 = -3 int Theta z1 z2^2 z3 - int Theta' z1 z2^3 (exact,
  kappa = gamma+4) implemented by analogy (M2-M4 -> 3 M1), not yet term by term. Core M1 box 0.75: gamma'_2 = 0.021 (1.06x on G_22);
  box 0.60: gamma'_2 = 0.128 (2x), M1 margin 1.45x (S1_v 0.37 + 0.045). E2 = 3.4 dominated by the ramp share (G_ramp,22 = 1.4e7 ~ M1_ramp^4
  with input 2.5; certified ramp M1 ~1.2 would cut it ~20x). F4: K0 <= 1.3e-21, K2 <= 1.8e-20. (B) valid (E_B = 4.2e-7, 3e-4, 3.5).
  H6: K2 exceeded for y <~ 75 with E2(tau1) = 0.2 -> hand-off requirement E2(tau1) <= 0.03 on [50,150]; K3 needs (B) level 3 from hand-off
  data (strip populated only by parcels present at tau1). Asked: term-by-term L^4, (B) level 3, explicit hand-off list.
- (Nov 10) LEMMA S term-by-term L^4 (lemmaS_L4terms.py): core M1 box 0.6 valid (G_22 9.3e4, 1.6x; M1 margin 1.45x), 0.5 valid (2.5x /
  1.2x), 0.75 INVALID. Ramp still by analogy (M_eff = 3 M1_ramp). H6: K2 needs E2(tau1) <= 0.03 (1.1x); K3 needs E2(tau1) <~ 3e-4 and
  E3(tau1) <~ 5e-4 on [50,150]; (B) level 3 (lemmaS_B3.py, sympy-checked) closes if strip M1 <= 0.3 with assumed |z2|,|z3| <= 1 (to
  re-derive). Hand-off list (9 items) in lemmaS.tex.
- (Nov 10) VERIFIER production (/tmp/q1.out, q2.out): (5e-4, 1e-7): tail rho 2.39e-3, alpha-cond margin 106, floor 2.0e-8, S_ramp 1.9e-3,
  A_max 1.37e-7, W1_max 1.72e-6, far sup|w'| 3.0e-3, sup|w''| 4.0e-2, M1_ramp <= 1.114. (5e-4, 2e-7): S_ramp 3.5e-3, sI1 5.5e-3, M1_ramp 2.05.
- (Nov 10) VERIFIER FIX PLAN DONE (VERIFIER.md "Fix plan", E_I_sups.json, rate_table_unc5_*, verify_log_fix5_*): both points close at
  alpha1 = 5e-4 with certified modulation constants. Tail needs theta = 0.05 (J's ~6x the theta = 0.3 values). Unconditional
  J (+ slope*(T - tau_t)): 1e-7: J0..J2 = 8.92e-11, 1.25e-10, 2.77e-11, J4 9.2e-7, slopes 2.9e-15/4.1e-15/9.0e-16; 2e-7: 1.69e-10,
  2.37e-10, 5.23e-11. Production 1e-7: peak 1.366e-7, rho 2.39e-3, alpha-margin 106, floor 2.0e-8, S_ramp 1.92e-3, sI1 2.99e-3,
  M1_ramp 1.114; 2e-7: peak 2.51e-7, rho 2.98e-3, floor 1.8e-8, S_ramp 3.54e-3, M1_ramp 2.05. Validity box (A, W1 <= 1e-3) OK
  (post-check). Placeholders: S1V = 3.0, F4 flux. Strip modulation terms (<1e-8/S) omitted. Forwarded to Lemma S (damped-shell form for
  floor slopes).
- (Nov 10) LEMMA S WITH CERTIFIED INPUTS (lemmaS_P1_M*.json, lemmaS_P2_M*.json): sqrt Lyap1 = 1e-7 CLOSES (core M1 box 0.6) but all
  margins < 2x: S0_v 0.0147 (1.36x), M0 core 0.058 (1.9x), M1 core 0.477 (1.26x), G_22 9.3e4 (1.6x); E = 1.85e-6, 1.31e-3, 1.39;
  gamma' = 0.246, 0.245, 0.099. (B) valid; F4 K0 <= 4.1e-21, K2 <= 2.0e-20. sqrt Lyap1 = 2e-7 FAILS (S0_v 0.0203, M1 core 1.04).
  H6 strip: K2 <= 0.102 (1.8x), K3 <= 0.120 (1.5x) (hand-off requirement). Limits: inflow from certified J0, J1 (theta = 0.05 tail).
  Next: ramp L^4 term-by-term (eta^4 weight); verifier production at sqrt Lyap1 = 5e-8; projected margins.
- (Nov 10) LEMMA S ramp L^4: term-by-term blocked by the cutoff foot (eta'/eta unbounded); needs J3 (available: J3 <= sqrt(J2 J4) + c
  J2) -> asked. Bracket at 1e-7, box 0.6: ramp M_eff 0 / 3.6 / 7.2 / 14.4 -> M1 margin 1.46 / 1.26 / 1.06 / 0.78. Projection to 5e-8:
  J's have a large Lambda1-independent part (~70% of J0), so 5e-8 gives only ~1.4-1.5x (S0_v 1.5-1.6x, M0 ~1.95x). Lever = forcing/floor
  part of J. Verifier: 5e-8 run with certified F4 constants (K0 4.1e-21, K2 2.0e-20) and S1V = 0.6, plus a J breakdown.
- (Nov 10) LEMMA S CLOSES at sqrt Lyap1 = 1e-7 (R5_* in lemmaS_consts.json): ramp L^4 term-by-term via interpolation + dissipation
  (foot term homogeneous, pure J2 data; J3/J4 unused). Core box 0.6: min margin 1.22x (M1); G_22 1.6x; S0_v 1.36x; M0 1.9x.
  E = 1.85e-6, 1.32e-3, 1.59 (E2: hand-off 0.2, inflow 0.28, ramp 1.11). (B) valid; F4 K0 4.1e-21, K2 2.1e-20. Remaining assumptions:
  sup z2 = 3 in the {2,3,3} class (4% of G_22); J_m <-> z_k conversion vs the verifier's definition; F4 horizon. Asked to close them
  and consolidate the statement + hand-off list.
- (Nov 10) VERIFIER fix6 (E_I_sups_fix6.json, rate_table_unc6_*): certified F4 constants + S1V = 0.6 -> floor 2e-8 -> 5e-10, tail
  delta_t = 0.2, J's ~halved. 1e-7: peak 1.366e-7, rho 2.40e-3, J0..J2 = 4.64e-11, 6.52e-11, 1.44e-11 (slopes ~1e-18), S_ramp 1.92e-3,
  ramp M1 1.114. 5e-8: peak 7.93e-8, rho 2.17e-3, J = 2.64e-11, 3.71e-11, 8.20e-12, ramp M1 0.646. J breakdown: forcing int 2 Phi(alpha)
  = 42% (1e-7) / 69% (5e-8) of the constant part; Lambda1 part 26% / 12%; tail anchor 32% / 19%. Forwarded to Lemma S for the final
  margin table (adopt 1e-7 if min margin >= 1.5x).
- (Nov 10) LEMMA S FINAL RUNS with fix6 J's (writeup/lemmaS/lemmaS_N1_*, N5_*): 1e-7: S0_v 0.0111 (1.81x), S1_v 0.341; best core box
  0.575 -> M1 1.49x, G_22 1.50x, gamma'_2 0.088, M0 2.04x. 5e-8: S0_v 0.0088 (2.27x), S1_v 0.235; box 0.52 -> M1 1.86x, G_22 1.93x,
  gamma'_2 0.124, M0 2.12x. DECISION: hand-off alpha1 = 5e-4, sqrt Lyap1 = 5e-8 (core box 0.52); min margin ~1.86x. Agent was interrupted
  by a server restart; resumed to finish the closing items and the consolidated statement + hand-off list.
- (Nov 10) LEMMA S CONSOLIDATED at (5e-4, 5e-8), core box 0.52: S0_v 0.0088 (2.27x), M0 0.0521 (2.11x), M1 0.280 (1.86x), G_22 7.8e4
  (1.93x), gamma' 0.246/0.246/0.124; E = 6.9e-7, 4.5e-4, 0.350. {2,3,3} class on the core by IBP (derived). J conversion checked (no
  mismatch). F4 uniform in T (K0 1.5e-21, K2 6.5e-21). 6-item hand-off list in lemmaS.pdf. GAPS (used for ALL tau, so bootstrap, not
  hand-off): ramp sup|z2| <= 30 on [14,28] (sensitivity mild); strip sups |z2|,|z3| <= 1, M1 <= 0.3 for (B) level 3. Asked: Lemma P
  (crude parabolic smoothing) for the ramp z2, and (B) up to level 4 (hand-off data) for the strip sups.
- (Nov 10) LEMMA S: strip [50,150] CLOSED (level 3 + Sobolev, no level 4: y|f'| <= 0.168 (1.8x), |z2| <~ 0.2 (2.5x)). Ramp {2,3,3}
  by IBP: sup z2 now only multiplies small cutoff-free data; tolerance sup|z2| <= 1e4 -> M1 margin 1.79x (1e3 -> 1.85x). Remaining gap:
  Lemma P (sup|z2| <= 1e4 on [14,28] for all tau); level-3 localized energy gives ~6e2 at B = 1e3 (uncertified). Asked to finish with
  J3 = sqrt(J2 J4) + c1W J2 (~1.5e-9 at 5e-8) and the first window from hand-off data.
- (Nov 10) LEMMA P: level-3 route blocked (cutoff foot eta ~ t^3 only controls t^3|z2|; need t^{4/3}). Green-kernel / mild-solution
  route outlined (||d^3 G_t||_L1 = c3 t^-3/4, c3 ~ 0.6; drift removed in the characteristic frame; Duhamel on time-1 windows; input
  only sup|v|): expected sup|z2| ~ 10-100 vs tolerance 1e4. Asked to carry it out with certified kernel constants. H-Z2 is the ONLY
  unproved far-field hypothesis.
- (Nov 10) H-Z2 CLOSED without Lemma P: foot term covered by cutoff-free g = (eps omega)^{1/2} z2 on the whole ramp, Sobolev in s +
  Hoelder in time (inputs J2^free 5e-7, pointwise F_max 1.58e-4 from Lambda, J3 = sqrt(J2 J4) + c1W J2 = 1.5e-9); data term 3.1e-3.
  FINAL (lemmaS_P2_M0.52.json): S0_v 0.0088 (2.27x), M1 core 0.279 (1.86x, min), gamma' 0.246/0.246/0.124. Hand-off list: 5 items.
  Remaining: grid maxima -> arb (asked). STEP 3 FAR FIELD CLOSED (modulo that). Dispatched: step-3 assembly write-up (sonnet) ->
  writeup/step3_proof/; resumed step-4 feasibility part 2 (macroscopic phase).
- (Nov 10) LEMMA S FINAL (tag P3; one-page summary at the top of lemmaS.pdf, 14 pp): grid sups already arb-cell enclosures; C_eta now
  arb (15.84); numbers unchanged (min margin M1 core 1.86x, G_22 1.93x, S0_v 2.27x, M0 2.11x). FINAL-RIGOR-PASS LIST (deferred, small):
  (1) Lemma S assembly arithmetic in arb/outward rounding (now doubles); (2) L^4 class weights from sympy floats -> exact/arb;
  (3) sliver x in [16.2491,16.25] (zeta < 1) neglected at < 1e-29 assuming |z_j| <= 1e3 -> replace by an L^2 bound;
  (4) verifier: modulation strip terms (< 1e-8/S) omitted -> include; validity box post-check -> in-cell.
- (Nov 10) STEP-3 ASSEMBLY WRITTEN (writeup/step3_proof/step3_proof.tex/pdf, 15 pp; INCONSISTENCIES.md, 13 items). Theorem version W
  with hand-off H1-H6, conclusions W1-W5; components [C] except Lemma S [A] (now arb-cell sups); bootstrap B1-B4 via D1-D6 with
  strict-improvement factors (min 1.74 S0_v vs verifier input, 1.86 core M1). Key inconsistencies sent out: S0V value used in the 5e-8
  run; hand-off E2 0.2 vs 3e-4 on [50,60]; Lemma O time ratio 1.2e-7 vs 1.6e-7 and far-cylinder margin 1.001; ramp M1 0.646 vs box 0.7
  (+ y|U'| 0.044 -> 1.01); flux linearity for the bootstrap; J2^free provenance; sup|a|, sup|b| not recorded; K3 level unspecified;
  5e-8 results missing from VERIFIER.md. Lemma S/O agent and verifier agent resolving.
- (Nov 10) INCONSISTENCIES 5-9 RESOLVED (Lemma S/O agent; INCONSISTENCIES.md): ramp box bounds full y|f'| (0.646 + 0.0203 + 0.01 = 0.68)
  -> box 0.8 (1.18x); rerun R8: min margin 1.83x (M1 core), S0_v 2.27x, G_22 1.93x. Hand-off: E2[14,60] <= 0.2; E1,E2,E3 <= 3e-5 on
  [50,150] and [60,350] (K2 <= 0.080, K3 <= 0.046). J2^free = 2(51.6 J1e + 3.02e4 J2e) from verifier J's. F4 linear in E (proved).
  Lemma O: needed ratio 1.6e-7, margin 2.31 (y_O 200, half-width 0.75r, X 3); far cylinders X = 2, L = 113 lambda1, M = 1, K = 1.5:
  t* >= 1.232e-8 vs 6.13e-9 (2.008x); step-4 far-curve jets (lambda1 units): |h''| <= 8.9e-3, |h'''| <= 7e-5, |h''''| <= 6.2e-7,
  slope <= 0.21.
- (Nov 10) VERIFIER resolutions: S0V = 0.0153 used in the 5e-8 run; sup|a|, sup|b| <= 3.36e-7, sup|beta| <= 5.8e-9 (mod_sups_5e-8.json);
  strip terms included (coefficients ~1e-7/S, not < 1e-8; absolute ~1e-9..5e-9) -> rerun fix7: nothing moves (J0 2.6424e-11, peak
  7.93e-8, rho 2.17e-3). VERIFIER.md "Fix plan, update 3". Loose end: tail_cert.py default prints "ALL STRUCTURAL CONDITIONS FAIL"
  (checks the homogeneous structure including absolute terms) -> asked to make it check exactly the floors-off part + floor argument.
  Assembly tex being updated with all resolutions; then Fable final review.
- (Nov 10) step3_proof.tex/pdf UPDATED (14 pp): H1-H6 with single strip requirement, Lemma O far-cylinder jets in H5; W1 with certified
  modulation sups; Lemma S (A),(B) tagged [C] (arb-cell sups; assembly arithmetic doubles noted). Strict-improvement table min 1.83 (core
  M1). Section 6 caveats (9): Lemma S assembly arithmetic; floor's (1+eta) factor; E0 alpha-term reading; margins < 2x incl. underived 0.01
  allowance; Lemma O not self-contained + tiny far jets; step-4 feasibility of hand-off, |f| <= 1/2 in (B) not tabulated; 3/64 vs
  omega = 1/32 link (!); citations / sign conventions; upgrade C. FINAL FABLE REVIEW requested -> STEP3_REVIEW2.md.
- (Nov 10) FINAL FABLE REVIEW (STEP3_REVIEW2.md): bootstrap SOUND and non-circular; tail logic sound (tail_cert log ends PASS); NOT closed
  as written because of statement-level defects: (1) embeddedness: theorem assumed/concluded embedded, but S1 is non-embedded -> immersed,
  Lemma O per arc, drop H8; (2) hand-off H3/H4 clash with the bent cone (y|U'| = 0.27 at 350; E1(U) on [50,150] = 5.6e-3 vs 3e-5) ->
  state for the deviation, K2/K3 must include the U-part (~0.117 vs M = 0.12); (3) MOD_ALPHA = 6.19*66.5/63.5 is an assumption -> rerun
  with worst case 25; (4) floor anchor eta = 1e-9 -> certify on [L*, 1.01 L*]; (5) ramp 0.01 allowance -> one interval sup; (6) K 5.535 vs
  1.5 inconsistency, lemmaO_v3 not self-contained, |f| <= 1/2 untabulated; (7) 3/64 vs 1/32 cosmetic (theta-form certified separately).
  Dispatched: verifier (3, 4), Lemma S (2, 5, |f|), assembly (1, 6-K, 7).
- (Nov 10) step3_proof.tex: immersed throughout (H6 -> Regularity, W4/5.3 per arc), Lemma O per arc with H8 dropped, remark that embeddedness is a step-4 initial-curve property; B4 K = 1.5; theta-form certificate cited, 3/64 caveat removed. Waiting: verifier MOD_ALPHA = 25 + anchor; Lemma S deviation thresholds.
- (Nov 10) VERIFIER fix9m25: 5e-8 closes with MOD_ALPHA = 25 (worst case, certified treatment) and anchor eta = 0.01; BUG FIXED: x-degree-0
  part of the (a,b) bound was inside H (breaks H/s^2 monotone) -> moved to L; tail_cert PASS; all pre-fix8 runs superseded. Numbers: peak
  8.19e-8, rho 2.88e-3, alpha-margin 130, floor 4.77e-10, S_ramp 1.15e-3, ramp M1 0.668, J0..J2 = 2.5423e-11, 3.5714e-11, 7.8909e-12,
  J4 2.6304e-7.
- (Nov 10) LEMMA S review-2 fixes: hand-off now for the DEVIATION v = f - U; H4 "y|f'| <= 0.12" withdrawn; (B) per-window boxes (|f| <=
  U(1.5y)+0.1, y|f'| <= 1.05U(1.5y)+0.1; lemmaS_Bbox.py), all gamma' >= 0.188; new hand-off: deviation E0 <= 2e-7, E1..E3 <= 3e-6 on
  windows [50,253], E4 <= 3e-6 at tau1; Lemma O jets incl. U: K2 0.044 (2.7x), K3 0.014 (8.4x), K4 0.041 (2.9x). (B) |f| <= 0.292 (1.71x),
  box margins 1.27x / 1.58x. Ramp allowance certified (alpha^2 sup y|h2''| = 1.96e-4; ramp y|f'| <= 0.6665 old / 0.6885 with fix9m25 ->
  box 0.8). Final rerun with fix9m25 inputs requested.
- (Nov 10) LEMMA S FINAL with fix9m25 (lemmaS_F9_M0.52.json): S0_v 0.0087 (2.30x), M0 0.0520 (2.12x), M1 core 0.281 (1.85x), G_22 7.8e4
  (1.93x), gamma' 0.246/0.246/0.124; ramp y|f'| 0.6885 vs box 0.8 (1.16x); ramp M0 0.0217 (2.3x); (B) windows gamma'_3 >= 0.188, |f| <=
  0.292 (1.71x), box margins 1.27-1.58x; F4 K0 1.46e-21, K2 6.5e-21; Lemma O jets 2.70x / 8.4x / 2.9x. STEP 3 CLOSED modulo the final-rigor
  pass (Lemma S assembly in double, sympy-float class weights, zeta<1 sliver) + literature/conventions + optional upgrade C.
  step3_proof.tex final update requested.
- (Nov 10) step3_proof.pdf FINAL (14 pp, status box: closed modulo final rigor pass). Add to the rigor-pass list: regenerate modulation sups (mod_sups) for fix9m25 (only sup|a| < 1/4 used downstream); restore the first-cell check value.

## Nov 10 night: hand-off moved to k ~ 2.3e4 (step-4 milestone 1)
- **Milestone 1 (float) result:**
  - Passes: tip, sqrt(Lambda), inner shells, chart extent, Lemma O jets (direct K4 = 7.2e-4 at k = 6081, i.e. 170x under M = 0.12).
  - Fails: only the outer-window E3/E4 of the deviation. These follow a pure cubic, v = -40 (y/k)^3, so E_j ~ k^-5.9.
- **Crossings by window:**
  - E4 needs k >= 1.6e4 on the windows to [169,253], and k >= 2.2-2.3e4 for [253,380];
  - E3 needs k >= 1.55e4 on [253,380];
  - E2 needs k >= 1.07e4 on [253,380].
- **Decision:** alpha1 -> 1.25e-4 (alternate 1.5e-4), sqrt(Lambda1) = 5e-8. Verifier rerun requested at both points.
- **Proposal for the Lemma S rerun:** E4 is used ONLY for the K4 jet at tau1 (lemmaS.tex finding 6: "hand-off statement, not propagated"). Replace "E4 <= 3e-6 on all windows" by a direct hand-off bound K4 <= 0.041 (incl. U) on the cone cylinder.
  - Step 4 (T8) certifies jets anyway.
  - It also saves one derivative level in step 4.
  - E3 becomes binding; it has margin ~7 at k = 2.2e4.
- **Also to check:** at alpha1 = 1.25e-4, U(380) ~ 0.072 (vs 0.29). The (B) per-window entry energies (3e-6) can likely be raised. Ask Lemma S for the largest e with all margins >= 2.

## Nov 10 night (2): step 3 re-closed at alpha1 = 1.5e-4 (Lemma S v3)
- **Verifier** (final configuration), certified at both alpha1 = 1.25e-4 and 1.5e-4:
  - peak sqrt(Lambda) 5.21e-8 / 5.31e-8; rho 2.86e-3 / 2.69e-3; floor 4.8e-10 / 4.9e-10;
  - J ~ 4.5x smaller than at fix9m25;
  - modulation sups regenerated: (a,b) <= 2.18e-7, beta <= 3.34e-9.
- **Lemma S v3**, run by me because the agent could not be resumed:
  - core box (0.11, 0.52), ramp box 0.9;
  - certified a/beta sups wired in (new --abar/--beta options; env ABAR/BETA/AL/SQ/RT);
  - min margin 2.05 (was 1.16).
- **Rigor finding:** the legacy beta formula (0.00244 + 63.5 alpha) sqrt(Lambda) was not an upper bound for |beta|. The forcing it feeds is negligible (<= 3e-11).
- **Hand-off changes:**
  - E4 dropped, replaced by a direct K4 <= 0.06 at tau1;
  - window entry raised to 1e-5.
- **Adopted point:** alpha1 = 1.5e-4. Step 4 can hand off at k1 ~ 2e4 with a 10% alpha margin.
- **step3_proof.tex** updated by its agent (15 pp).
- **Rigor pass still open:**
  - assembly double -> arb;
  - sympy-float class weights;
  - zeta < 1 sliver a-priori bound;
  - first-cell check;
  - literature.

## Nov 10 night (3): Lemma S rigor pass closed (writeup/lemmaS/rigor/RIGOR_PASS.md)
- **Closed:** assembly in arb; exact class weights (deviation <= 4.6e-15); unconditional zeta<1 sliver bound (<= 8.3e-8 relative); ramp-commutator tables recomputed in arb (stored doubles differed by <= 2.3e-6 per component, <= 1.4e-7 in any quoted quantity).
- **Page-1 decimals** are rounded outward, in lemmaS.tex, lemmaS_consts.json and step3_proof.tex. Minimum margin 2.05.
- **LaTeX fix:** unescaped underscore in lemmaS.tex rebuilt; 18 pp, no errors.
- **Step-3 rigor items still open:** the verifier first-cell check; literature.

## Nov 10 (15:07 PT): independent reproduction on the user's Windows 11 machine
- **Setup:** kit v1 on Python 3.13.7 with python-flint 0.9.0; i7-10700, 2017 s wall.
- **check.py:** 85 PASS, 0 FAIL, 6 WARN. Identical to the Linux run, with the same 6 known WARNs; kit v2 is fixing them.
- **Reproduction quality:**
  - verifier tables at both alpha points, and their J's, bit-identical;
  - IBP pairings, the Lemma S rigor assembly (arb) and Lemma O exact;
  - float parts within 2e-9 relative.
- **Record:** kit_step3/user_runs/windows11_v1_nov10.txt.

## Nov 10 (16:10 PT): v3.1 closed
- **Kit v2** (91 PASS / 0 FAIL / 0 WARN) built and sent to the user.
- **Pipeline fixes:** constants rounded outward; mod_sups.py added (|(a,b)| <= 2.19103e-7 at alpha1 = 1.5e-4; the old 2.183e-7 did not dominate the formula); Lemma S rerun with safe inputs (59 inequalities hold, min margin 2.05); Lemma O decimals rounded outward; far-curve |h''| threshold 8.9e-3 -> 8.8e-3.
- **step3_proof.tex** updated to v3.1 and rebuilt with 0 errors. Details in CHANGES_v3.1.md.
- **Nov 10, 16:51 PT:** the user's Windows 11 run of kit v2 gives 91 PASS, 0 FAIL, 0 WARN, 0 SKIP. Every comparison is exact, including the float parts, since both runs used PYTHONHASHSEED=0. Record: kit_step3/user_runs/windows11_v2_nov10_check.txt.
