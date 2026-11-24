# Modulation lemma by integration by parts (F3)

Norms: A=||w~||, W1=||w~'||, B=||w~''|| in L2(rho) (half line); Q := (A+W1+B)^2. Jet box eta = 0.02 for |w|,|w'| on the core (sigma<=29). Strip S_i = sup_{sigma in [29,32]}|w^(i)|.

## Derivation (ibp_sym.py, exact sympy, residuals 0 -> verify_symbolic.json)
N(w) = F(a+w)-F(a)-L0 w, a = jets of u0 = alpha psi3 + alpha^2 h2 (C^4), L0 = linearisation at u=0 (checked = DF(0)).
F = c4(U0,U1) u4 + c3(U0,U1,u2) u3 + G(U0,U1,u2) (quasilinear, from graph_struct.pkl), Psi_j = chi f_j rho. Split N = Lin + R:
* Lin = DF(a)w - L0 w = sum dC_i w_i:  int Psi Lin = int w Lam, Lam = sum (-1)^i Dx^i(dC_i Psi). The a5 terms cancel, so u0 in C^4 suffices.
* R = D4 w4 + D3 w3 + Rrest (D4=c4(U)-c4(a), D3=c3(U)-c3(a)). Psi D4 w4 = Dx(Psi D4 w3) - Dx(Psi D4) w3; Theta = Psi D3 - Dx(Psi D4) is affine in w2, P = int Theta dw2, Theta w3 = Dx P - Pex.
  Then int Psi R = int QR, QR = Psi Rrest - Pex, free of w3,w4, QR = O(|w|^2), uses only Psi, Psi', Psi'' and w,w',w''. Degree in w'' is 2: **the w''^3 term cancels exactly**, so there is no sup|w''| B^2 term.
Boundary terms vanish (w, u0, Psi even; chi has 4 vanishing derivatives at sigma=32).
Float cross-check (float_check.py, 3 seeds, random even w, u0, Psi, ODE profile): direct vs IBP agree to rel 1e-11 (abs diff <=4e-11 on O(1-40) values): float_check_*.json.

## Bounds (ibp_run.py, ibp_final.py, interval arithmetic on the 3425 outA5 hull boxes, C4 jets, alpha ball [-alpha1,alpha1])
Pairing x = H R, H=(G^-1)^T as in the old lemma. Per component i in {a,b,beta}:
1. Linear: |<w~,L*f_j>| <= c_j A with c = (13.198, 6.793, 6.645e-5); alpha-part |alpha <w,(d_alpha L)* chi f_j>| <= m_j|alpha|A with m = (230.8, 117.5, 0.139) at alpha1=5e-4 (mean value in alpha over the ball, i.e. exact in alpha). Combined via H: lead = (3.616, 0.867, 0.00244), |(a,b)|: 3.708 (unchanged); new IBP m_i = (76.9, 62.9, 4.48), m_ab = 98.6.
2. Quadratic: |int rho G_ab w_a w_b| <= sum sup|Gamma_ab| N_a N_b; q = max over the simplex (rigorous branch and bound, simplex.py).
3. Cubic: Lagrange remainder in (w0,w1) over the eta-box (8x8 cells to cut dependency): cubic <= eta[(M3/3)(A+W1)^2 + M2 (A+W1)B + 2 M1 B^2].
4. Feedback of E(alpha,w) (old Neumann step, nu=0.0466, hypotheses |alpha|<=alpha1, A,W1<=1e-3, strip |w|,|w'|<=0.02) folded into m and q. It contributes the W1*A couplings (50.9, 90.3, 13.0) which the previous verifier dropped (old "13.4 W1 A"); they are part of q_total.

## Final inequalities, alpha1 = 5e-4, eta = 0.02 (final_5e-4.json, key `final.core`)
|beta|  <= (0.002440 + 66.5|alpha|) A + 5.53 Q + 0.657||chi Phi_F|| + strip
|(a,b)| <= (3.708 + 622.4|alpha|) A + 51.3 Q + 3.84||chi Phi_F|| + strip
(components: a: 3.616+332|alpha|, q 33.0; b: 0.867+520.5|alpha|, q 39.5.)

Mapping to the verifier: BETA0=0.00244, BETA1=m_total(beta)=66.5, q_beta=q_total(beta)=5.53 (IBP part alone: 4.78 incl. cubic, 2.31 quadratic); CAB=3.708, m_ab=622.4, q_ab=51.3 (IBP part alone: 27.9 incl. cubic, 20.1 quadratic; remaining ~23 is the E-matrix W1*A coupling). m_total is dominated by the old E^alpha feedback (255/457/62); the IBP adds 77/63/4.5. Placeholders: q_beta=50 is safe (5.5); q_ab=30 is NOT safe if the E-coupling is included (51.3), fine if only the IBP part (27.9) is meant.
Variant `near` (final.near): far region x>16.25 handled with sups instead (coefficients in run_*.json `far_L1`, ~1e-5): beta q_total 5.26, b IBP part 15.0 (a unchanged, tip-dominated).

Strip terms (negligible, <=1e-8 per unit S): see `strip` in the json: linear |alpha| L1 S0, commutator coefficients of S0..S3, quadratic/cubic L1 coefficients of S0,S1,S2.

Caveats: hull boxes (ws3 up to 20% wide) inflate far-region sups; far-region (x~31) values dominate sup|Gamma| for (a,b) via H_i3 f3; cubic term not small relative to quadratic for beta/b in the core variant (ratio ~1) but comparable to eta*M ; E-feedback bound is crude.
Files: ibp_sym.py, ibp_codegen.py -> ibp_gen.py, ibp_boxes.py, ibp_consts.py, ibp_run.py, simplex.py, ibp_final.py, float_check.py, make_json.py.

## alpha1 = 1e-3 (final_1e-3.json)
|beta| <= (0.00244 + 69.7|alpha|)A + 7.62 Q; |(a,b)| <= (3.708 + 661.1|alpha|)A + 57.5 Q (IBP part alone: beta 7.0, ab 39.4). Both q grow because the far-region Gamma is linear in alpha.

## Rigor update: simplex bounds in ball arithmetic
`simplex.py` now evaluates the cell bounds (value at the upper corner; forms are nonnegative and increasing) in arb, with exact dyadic cell endpoints, and reports the upper end of the largest ball; floats are used only for pruning. Rerun of `ibp_final.py` for both alpha1 gives q values identical to the earlier ones at the printed precision (4 decimals): alpha1=5e-4: q_beta 5.526, q_ab 51.34 (a 32.96, b 39.46); alpha1=1e-3: q_beta 7.623, q_ab 57.50. `mod_ibp.json` regenerated.
