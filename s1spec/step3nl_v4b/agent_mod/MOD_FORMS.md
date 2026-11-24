# Rigorous modulation-form bounds on {f1,f2}^perp (L^2(rho), half line)

**Results** (`mod_forms.json`):

| quantity | rigorous | float (exact ODE derivatives) | old crude bound |
|---|---|---|---|
| R1 = sup abs(Q1)/norm(w)^2 | **<= 13.5491** (witness >= 13.525) | 13.546 (eigs of Q1: -13.546, 11.326) | 268 |
| L2max = sup Q2/norm(w)^2 | **<= 43.3271** (witness >= 43.318) | 43.324 | 458.4 |
| L2min = inf Q2/norm(w)^2 | **in [-0.18936, -0.18428]** | -0.1887 | |

(Q1 separately: sup <= 11.3299, inf >= -13.5491.)

**Float reference in the task is off.** `mod_form2_float.py` masks x in (0.3, 30.5), which drops the tip region where most of L*f and D f lives. With the mask removed (`float_nomask_fd.py`) or with exact ODE derivatives (`float_gram.py`) the float eigenvalues are Q1 [-13.5, 11.3], Q2 [-0.19, 43.3], not [-10.15, 9.42] and [-1.03, 23.96]. The rigorous norms agree with the certified ones: norm(Lf1,Lf2) = 15.067 (lyap_consts2 crude upper bound 15.59), norm(Df1) = 11.775 (dmod_norms 11.854), norm(Df2) = 6.135 (dmod_norms 8.69), 1/lambda_min(M) <= 2.00039 (lyap 2.00046).

**Method.**
1. `mod_gram2.py`: rigorous 8x8 Gram of V = (L*f1, L*f2, gA, gB, D f1, D f2) and F = (f1, f2), with L*, D_rho, chi, psi3^ as in the task. Functions are evaluated as truncated Taylor series in x - c (arb_series). The state series come from the ODE p'=-kq, q'=1+kp, k'=m, m'=p/4, (cos th)'=-k sin th, (sin th)'=k cos th (psi2 = -sin th, cos th = +St*C - Ct*S in the Lyapunov frame; checked against the float theta in `check_frame.py`; this is the opposite sign to the one quoted in the task). Per box: int H = len*H(c) + H'(c)*int(x-c) + (1/2 H''(hull))*int(x-c)^2, with H = rho u_i u_j, H(c) from the tight centre enclosure and H'' from the hull. Boxes: 51246 tight-profile boxes `prof_tight.txt` for x <= tau_tip, with theta from X5 (= theta - theta(sigma=8)); for sigma in [8,32], the 2400 chunk-bound rows of outA5.txt ('MB') split into 4 sub-boxes each, with the C^3 smoothstep chi series (chi strip included). Gram radii are about 1e-6 relative (the plain box-hull version `mod_gram.py` had about 7%, which is useless because L*f1 and L*f2 have correlation -0.998).
2. `mod_forms.py`: Pi-projected Gram = Schur complement G_VV - G_VF Minv G_FV (interval Minv, 2x2). With ell = (<w,Pi Lf1>, <w,Pi Lf2>, <w,Pi u3>, <w,Pi u4>), Q = ell^T C' ell, C' = T^T C T (T encodes ha, hb through Minv; C as specified). The image of the unit ball of F^perp is the ellipsoid G^(1/2) B, so sup Q = lambda_max(C' G) = lambda_max(R C' R^T) with G = R^T R.
3. Rigour against the interval Gram: G is nearly singular (smallest normalised eigenvalue 2e-5 for Q1, 7e-6 for Q2), so interval Cholesky of [G] is not used. Instead, for every true G in [G], G <= G+ := mid(G) + diag(row sums of radii) (diagonal dominance). Then E(G) lies inside E(G+) and sup <= lambda_max(C' G+). G+ is an exact matrix, positive definite by arb Cholesky. lambda_max(R+ C' R+^T) is bounded by proving c*I - mid(B) - diag(row radii) positive definite (arb Cholesky). The inf is bounded the same way using -C'. Inner witnesses are Rayleigh quotients of explicit vectors w = sum beta_a Pi v_a, evaluated in ball arithmetic. Nothing is regularised beyond this majorant, whose cost is about 1e-4 relative.

**Checks.** Rigorous Gram midpoints agree with the independent float Gram (`float_gram.py`, exact ODE derivatives, trapezoid h=0.002) to 1e-5 relative or better. Group 3 / 1 sub-box and group 1 / 4 sub-boxes give the same bounds to 1e-7. Forms on the unconstrained space are not computed.

**Runtime.** Gram 125 s (group 1, 4 sub-boxes; 39 s for group 3), forms under 2 s, single core.

**Caveats.**
- The 8x8 Gram is for the f-version with only f1, f2 constrained. The f3 stretch was NOT done.
- psi3^ (jets of psi3md2_14), the tight profile and the outA5 hulls are trusted as certified inputs; the tip-region/far-region interface (x = tau_tip) uses the true tip location consistently (box lengths exact in tau; the last box overhangs x<0 by <1e-9 and is corrected).
- gA, gB use cos(theta) with the sign above. Flipping only cos(theta) would change gB's sign relative to hb and give a different Q1.
- Files: `mod_gram.py` (crude version + shared symbolics), `mod_gram2.py`, `mod_forms.py`, `make_json.py`, `check_frame.py`, `float_gram.py`, `float_nomask_fd.py`; data `gram8_g1s4.json`, `gram8_g3.json`, `forms_*.json`, `mod_forms.json`.
