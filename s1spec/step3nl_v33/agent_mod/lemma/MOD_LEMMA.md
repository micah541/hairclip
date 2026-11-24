# Modulation lemma (a, b, beta) for u = alpha psi3^ + alpha^2 h2^ + w

All norms are half-line L^2(rho). A = ||w~||, W1 = ||w~'||, B = ||w~''||, W3 = ||w~'''||, D4 = ||w~''''||, S_i = sup over sigma in [29,32] of |w^(i)| (raw w, i = 0..4), S2 also written sup|w''| in the nonlinear terms, Phi = ||chi Phi_F||.

## Setting and hypotheses
- |alpha| <= 1e-3, |w|, |w'| <= eta on the whole support (main: eta = 1e-2; tables also for 1e-3 and 1e-4). Then alpha_g = 1 - k u >= 0.98306 (sup|k| = 1.5403, sup|k u0| = 1.54e-3).
- Constraints <w~, f_j> = 0, j = 1,2,3, with f1 = chi p, f2 = chi psi2, f3 = chi F - a_c eta1 - b_c eta3.
- x = (a, b, beta) solves  x^T (G + E) = R  (R_j as in the task). Equivalently x = H R - H E^T x, with H = (G^-1)^T.

## Statement (numbers for eta = 1e-2; `mod_lemma.json`, key `main`)
Write x_i for a, b, beta, or the vector (a,b). Then
  |x_i| <= T_i + F_i,
  T_i = c_iA A + c_iW1 W1 + c_iB B + c_iW3 W3 + c_iD4 D4 + c_iS2B (S2 B) + c_iS2W3 (S2 W3) + c_iS2S2B (S2^2 B) + c_iPhi Phi + (strip terms, each coefficient <= 9e-9),
  F_i = [ |alpha| n0_i + W1 n1_i + |alpha| A n2_i + S1 n3_i ] * (sum_q tau_q q) / (1 - nu),
with nu = 0.0798 (so 1/(1-nu) = 1.087), valid whenever W1 <= 1e-3, A <= 1e-3, S_i <= 1e3.

| coefficient | a | b | beta | (a,b) |
|---|---|---|---|---|
| c_A (eta=1e-2) | 19.4 | 27.5 | 2.70 | 34.7 |
| c_A (eta=1e-3) | 6.31 | 5.40 | 0.462 | 8.98 |
| c_A (eta=1e-4) | 5.59 | 4.19 | 0.339 | 7.57 |
| leading part of c_A (no nonlinearity) | 3.616 | 0.867 | 0.00244 | 3.708 |
| c_W1 (1e-2 / 1e-3) | 17.0 / 3.01 | 28.6 / 5.07 | 2.90 / 0.514 | 33.2 / 5.90 |
| c_B (1e-2 / 1e-3) | 1.04 / 0.212 | 1.75 / 0.357 | 0.177 / 0.0362 | 2.03 / 0.415 |
| c_W3 (1e-2 / 1e-3) | 0.495 / 0.0799 | 0.833 / 0.134 | 0.0844 / 0.0137 | 0.968 / 0.156 |
| c_D4 (1e-2 / 1e-3) | 0.274 / 0.0469 | 0.461 / 0.0792 | 0.0471 / 0.00827 | 0.537 / 0.0921 |
| c_S2B (1e-2 / 1e-3) | 21.6 / 17.5 | 36.3 / 29.5 | 3.68 / 2.99 | 42.2 / 34.4 |
| c_S2W3 (1e-2 / 1e-3) | 1.02 / 0.669 | 1.76 / 1.16 | 0.209 / 0.143 | 2.04 / 1.34 |
| c_S2S2B (1e-2 / 1e-3) | 14.1 / 12.2 | 23.9 / 20.7 | 2.59 / 2.26 | 27.8 / 24.0 |
| c_Phi | 1.735 | 3.141 | 0.626 | 3.649 |
| n0 (alpha) | 32.2 | 57.8 | 7.84 | |
| n1 (W1) | 6.43 | 11.4 | 1.65 | |
| n2 (alpha A) | 3.06 | 4.64 | 0.280 | |

tau_q (eta = 1e-2): A 40.1, W1 35.0, B 2.14, W3 1.02, D4 0.565, S2B 44.5, S2W3 2.11, S2S2B 29.1, Phi 5.23, strip terms about 1e-8 or smaller. For eta = 1e-3: tau_A 13.0, tau_W1 6.22 (all tau in `mod_lemma.json`).

**Linear corollary** (nonlinear, forcing and strip terms off; this is the closed form of "c0 A + c1 |alpha| (A + W1)"):
  |a| <= (3.616 + 261 |alpha| + 52.1 W1) A,
  |b| <= (0.867 + 468 |alpha| + 92.6 W1) A,
  |(a,b)| <= 3.708 A + sqrt(a_fb^2 + b_fb^2) with a_fb, b_fb the same feedback parts,
  |beta| <= (0.00244 + 63.5 |alpha| + 13.4 W1) A.
The 0.00244 is eps1/G33-type. For comparison, the old crude bound was |(a,b)| <= 31.2 A.

## Proof outline and where each number comes from
1. **G** (`G` in the json). G_11, G_12, G_22 = Gram of (f1, f2) (rigorous 1e-6, `../gram8_g1s4.json`). G_31, G_32 = <chi psi3^, f_j> from `psi3ip.json` (+- strip difference <= 1.1e-8). G_33 in [0.031107, 0.033604]. |G_13| <= 4.01e-7, |G_23| <= 1.88e-7 from lyap_consts2 (`psif3`) plus the strip difference. ||f3|| <= 0.0194765 (lyap_copy.py rerun, `lyap_extra.json`).
2. **H = (G^-1)^T** by the block formula with Schur complement s = G_33 - r M^-1 c (ball Gaussian elimination blows up; the block formula was verified G*G^-1 contains I). H_11..H_22 = M^-1 to 1e-4; H_31, H_32 = O(1e-5); H_13 in [70.9, 78.7], H_23 in [133, 148], H_33 in [29.8, 32.1]; ||G^-1||_inf <= 259.
3. **Leading term.** w~ is orthogonal to f1, f2, so <w~, L*f_j> = <w~, Pi L*f_j> for j = 1,2, with the 2x2 Schur-complement Gram of Pi L*f (13.2^2, 6.79^2, ...). With h_i = (H_i1, H_i2): |sum_j H_ij <w~, L*f_j>| <= sqrt(h_i^T G_pi h_i) A, and the vector norm uses lambda_max(C G_pi C^T). The j = 3 part is |R_3^lin| <= eps1 A (eps1 = 6.645e-5, lyap_consts2), times |H_i3|.
4. **Forcing.** <chi Phi_F, f_j>: |.| <= ||sum_j H_ij f_j|| Phi, from the Gram of (f1, f2, f3).
5. **Commutator** [L, chi] w = -4chi' w''' - 6chi'' w'' - (4chi''' + 2k^2 chi') w' - (chi'''' + k^2 chi'' + (5km + q/4) chi') w (checked symbolically). Pairing with f_j needs only the strip sigma in [29,32], where rho <= 2.6e-11: coefficients kappa_{j,i} = int_strip rho |f_j| |c_i| <= 1.0e-28 (j=1), 7.7e-11 (j=2), 4.8e-11 (j=3).
6. **Nonlinear remainder.** Nw = D = sum d_tjl w''^t w'''^j w''''^l + d_000 (nlin_E2.py, rerun here for eta = 1e-2, 1e-3, alpha <= 1e-3, regions x <= 12 and x > 12). Used: |<chi D, f_j>| <= sum_reg ||f_j||_reg [d001 D4 + d010 W3 + d100 B + A0 A + A1 W1 + d200 S2 B + d110 S2 W3 + d300 S2^2 B]_reg, plus an strip piece int_strip rho|f_j| times the same in S_i. Region norms of f_j: ||f1|| = 1.84 (x <= 12; 5e-6 beyond), ||f2|| = 1.19 (0.034 beyond), ||f3|| = 0.0193 (0.0035 beyond).
7. **E.** Since <chi w, f_j> = <w~, f_j> = 0, E_1j = <chi u0, f_j> - <chi q T(u), f_j> with T(u) = u'/alpha_g, and E_2j = -<chi cos(theta) T(u), f_j>.
   - alpha-part: |E_1j|, |E_2j| <= |alpha| e_{1,2;j} with L^1 hull integrals (e1 = 4.23, 10.8, 0.347; e2 = 1.02, 0.786, 0.0163, j = 1,2,3; alpha <= 1e-3 in the h2 terms).
   - w-part: T(u) - T(u0) = w'/alpha_g + u0' k w/(alpha_g alpha_g0), so |E^w_1j| <= ||q f_j/alpha_g|| W1 + |alpha| ||q f_j k (psi3' + alpha h2')/alpha_g^2|| A (+ strip terms in S1, S0). The L^2 norms ||q f_j||/alpha_g are 1.94, 2.29, 0.0935.
   - E_3j = 2 alpha <chi h2, f_j>, |E_3j| <= 2|alpha| ||h2|| ||f_j|| with ||h2|| <= 13.11.
8. **Neumann.** |x| <= t + N|x|, N_ik = sum_j |H_ij||E_kj|. A weight vector d > 0 with N_ref d <= nu d, nu = 0.0798 (checked in ball arithmetic at alpha = 1e-3, W1 = A = 1e-3, strip sups 1e3), gives |x| <= tau d/(1-nu) with tau = max_k t_k/d_k, hence the form of F_i above. `neumann_test.py` Monte-Carlo checks the algebra (ratio bound/actual stays below 1).

## Caveats
- **D4, W3 and S2 = sup|w''| appear.** The nonlinear pairing was bounded by Cauchy-Schwarz with the certified expanded coefficients, not by integration by parts. That would need the divergence-form structure and would still leave w''-type cubic terms. The D4, W3 coefficients are small (<= 0.54), but the S2*B-type coefficients are 17 to 42 because the d200 (w''^2) coefficient is about 4.8 at the tip independent of eta and alpha.
- **A0, A1 are sups over |alpha| <= 1e-3**, not proportional to alpha. They dominate c_A beyond the leading part: at eta = 1e-4 tip A0 = 0.54 for alpha = 1e-3 and 0.158 for alpha = 1e-4. Rerunning `nlin_E2.py AMAX eta eta` with smaller AMAX sharpens this directly. The eta = 1e-2 jet box mostly inflates c_A, c_W1 by factors 3 to 5 versus eta = 1e-3.
- **E^alpha is crude** (L^1 integrals of absolute values over the 0.01 hull boxes), giving nu = 0.08. This is harmless (nu is far below 1), but n0, n1 could be reduced.
- **Trusted inputs.** G_3j from psi3ip (+-2%), eps1 and psif3 from lyap_consts2 (not re-derived), outA5 hulls, nlin_E2 coefficients. The first-order claim <w~, L*f3> <= eps1 A is taken as in ASSEMBLY.
- **Projection.** The leading term only projects off f1, f2 (not f3), a negligible loss.
- **Region norms** of f_j use hull boxes (about 0.5% above the accurate Gram); strip and tail sums are included. The strip piece is slightly double counted in the region norms.
- **No dynamical validation.** I did not check the modulation identity against a nonlinear flow run; only the commutator formula (symbolic), G*G^-1 = I, and the Neumann algebra (Monte-Carlo) were tested.
- **Frame.** Frame conventions as in `../MOD_FORMS.md` (cos(theta) = St*C - Ct*S).

Scripts: `lem_integrals.py` (hull integrals, `lem_integrals.json`), `lemma_consts.py` (algebra; args: nlin json, tag), `make_lemma_json.py` (-> `mod_lemma.json`), `neumann_test.py`, `lyap_copy.py` (copy of lyap_consts2 dumping `lyap_extra.json`). Runtimes: a few seconds each, except `nlin_E2.py` (about 1-2 min per eta).
