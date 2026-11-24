# Step 4 float feasibility (clip 24/0.4/0.05/0.85). All numbers non-rigorous. Scripts and logs in this directory.

Files: `common.py` (Fourier eval), `f1.py f1sum.py` (F1), `run_f2.py worker.py perturb.py basex.py chart.py f2an.py f2tab.py` (F2), `f3geom.py f3fem.py f3fd.py f3val.py f3run.py f3batch.sh` (F3); logs in `logs/` (`f1_summary.txt`, `f2_table.md`, `f3_*.log/json`).

## Pipeline check (done first)
Even-part chart (both strands, same frame, axis through the tip, joint fit of (x_c, lam, alpha, beta2); `chart.py`): r2/|alpha|^3 = 39.1 (k=100), 40.1 (400), 43.5 (806), 46.6 (1600). This matches ~45|alpha|^3.
The single-strand axis-through-tip fit gives 203|alpha|^3 at k=1600 instead. The snapshots carry a physical odd asymmetry of ~4.3e-9 near the tip (cos-coefficient 6e-10 in theta from the initial data; closure gap 7.8e-7 at the far end). It is constant in physical units, so it grows like k in lam units. All F2 numbers use the even part.

## F1 (resolution). DOF = sum(p_i+1), Chebyshev degree <=128 reference, geometric partition, width-1 elements on [0,4]
Rule: all trailing coefficients of theta < tol (coefficients < 3e-15 are rounding noise). Snapshots k=10, 30 are from my base run; 100, 400, 1600 are the existing design snapshots. sigma_end = 159, 470, 1567, 6288, 25164. Extrapolation: linear in ln k, fitted on k >= 100.
| tol 1e-13: r | k=10 | 30 | 100 | 400 | 1600 | k=2900 |
|---|---|---|---|---|---|---|
| 0.25 | 250 | 296 | 358 | 408 | 462 | 484 |
| 0.5 | 200 | 219 | 281 | 321 | 351 | 368 |
| 1.0 | 171 | 183 | 247 | 271 | 295 | 305 |
Same rule, other tolerances at k=1600 (and k=2900): tol 1e-12: 411/307/259 (431/320/269) for r=0.25/0.5/1. tol 1e-11: 375/272/221. tol 1e-10: 338/245/193. See `logs/f1_summary.txt`.
- Max degree p at tol 1e-13: 19-20 (r=0.25), 28 (r=0.5), 35 (r=1). Element counts at k=1600: 44/26/17.
- DOF grows only like ln k: +25 to +37 per e-fold. All values are far below 2.6e5 uniform. Larger r gives fewer DOF, at the price of p ~ 35.
- sup error of theta'''' (element-wise, tail bound) when truncating at tol 1e-13 is 2.4e-6 (r=0.25), 5e-6 (r=0.5), 9e-5 (r=1). At tol 1e-12 it is 2.5e-5 to 1.3e-4. sup|theta''''| itself is ~1.
  So theta'''' to 1e-12 is not resolved at the 1e-13 level of theta. Extrapolating the geometric coefficient decay below the noise floor gives roughly 160-410 DOF for theta'''' error <= 1e-8 to 1e-12 (rows in the summary). That estimate is noisy (non-monotone in k), so treat it as ~250-400 +-30%.
- Floors: Chebyshev evaluation noise ~1e-15 (theta ~ O(1)). The snapshots are spectrally filtered. Fourier coefficients in the last 10% of modes are up to 2e-12 (k=100), 5e-13 (400), 1.2e-13 (1600), and 5e-13 (k=10).
  So the 1e-13 rows are at or below the data floor for k=10, 100, 400, and only the 1e-12 rows are fully meaningful there. Time-stepping error of the snapshots (quoted 1e-11 to 1e-12) is not measured here.

## F2 (sensitivity modulo T, translation, alpha)
Setup: base run from scratch, snapshots at k=10, 30, 100, 400, 1600 (kmax(t=0)=42.1; it dips to 1.8 at t~4, then rises, so the early levels are reached on the rise at t~4.598). Perturbations: even normal displacement e(s) (Gaussian pair, width 2*lam_j `tip`, or 0.05 `chamber`, peak = delta*lam_j resp. delta*0.05).
- The compensator is chosen so that int kappa e = 0, and dtheta = e' + kappa*int_0^s kappa e is added on the equal-arclength grid, so L and closure are preserved. Checked: closure unchanged, normal displacement matches e to 1e-6, symmetry to 1e-13.
- delta = 1e-6 and 2e-6. Resumed runs reproduce the base bit-for-bit (control from k=400: difference exactly 0).
- Comparison at each k: fit (x_c, lam, alpha, beta2) in both runs, dw = w_pert - w_base on the common chart. The kmax mismatch between runs (~1e-4 relative, which gives dw ~ 1e-7, comparable to the signal) is removed by transporting the base along ln lam with a derivative from extra base snapshots at 1.003 k (`basex.py`).
- Norms: half line, rho (ρ) and W=(1+x²)^-2 on the chart x <= 25.7 (k=100) or 27.7. The chart does not reach 260, so the W numbers are a lower-range version.
- A_rho = ||dw(k)||_rho / ||e/lam_j||_rho (raw size, column "size/δ" = 0.2-2.5). A_mod = ||dw(k)||/||dw(k_j)|| (size taken modulo symmetries at k_j; not available at k_j=10, 30 for the tip chart). dalpha/δ = (α_p - α_b)/δ. The last column is the linearity ratio, (δ=2e-6 result)/(δ=1e-6 result), ideal 1.
Full table: `logs/f2_table.md`. Condensed (delta=1e-6):
| k_j | shape | A_rho(100) | A_rho(400) | A_rho(1600) | A_mod(1600) | A_W(1600) | dα/δ at 1600 | linearity |
|---|---|---|---|---|---|---|---|---|
| 10 | tip | 1.3e-2 | 3.8e-4 | 1.2e-4 | 9e-4 | 3.5e-4 | 2e-4 | floor-limited |
| 10 | chamber | 4.4e-2 | 1.9e-3 | 3.3e-4 | 5.6e-4 | 2.0e-4 | 3.6e-4 | floor-limited |
| 30 | tip | 3.8e-3 | 7.0e-4 | 2.9e-5 | 1.2e-4 | 5.7e-5 | 1.5e-4 | 0.6-1.0 |
| 30 | chamber | 3.3e-3 | 7.2e-4 | 7.6e-5 | 2.9e-4 | 6.1e-5 | -2.1e-3 | 0.2-1.1 |
| 100 | tip | 0.32 (k_j itself) | 1.8e-3 | 1.0e-4 | 3.3e-4 | 6.9e-4 | 2.3e-3 | 1.00 |
| 100 | chamber | 0.29 | 2.5e-3 | 1.5e-4 | 5.2e-4 | 2.2e-3 | 2.2e-2 | 1.00 |
| 400 | tip | - | 0.33 | 1.8e-3 | 5.5e-3 | 1.7e-2 | 1.5e-2 | 1.00 |
| 400 | chamber | - | 0.069 | 4.8e-4 | 6.9e-3 | 4.1e-3 | 2.0e-3 | 1.00 |
- Result: no amplification at all. Every A is <= 0.33, and the largest is the immediate response at k_j. The response decays by ~1e-2 per factor 4 in k (consistent with the stable spectrum, e^{-0.7*5.5}=0.02), so A(k_j -> 1600) <= 7e-3 for k_j >= 100 and <= 3e-4 for the early k_j. This is far below the 1e4 threshold.
- Floor: dw noise ~1e-10 (absolute) from the transport correction and fit. Early-k_j runs (δ·size ~ 1e-6) fall to that floor at k >= 400, so those linearity ratios are not meaningful and the numbers are upper bounds.
- Anomaly: `p10_tip_1` at k=30 gives dw/δ = 100, a failed fit of the poorly formed chart (xe 9.9) in that run only (`p10_tip_2` gives 0.34). The k=30 chart is unreliable in general (end 9.9).
- Δα/δ is O(1e-4 - 0.25): alpha shifts are not amplified. The T-shift and translation are absorbed by the λ, x_c fit and not reported.

## F3 (tail dimension)
Linearization (arclength sigma, N = iT, kappa = theta_s, T' = kappa N): delta kappa = e'' + kappa^2 e,  delta(Y.N) = e - q e' (q = Y.T, p = Y.N), so
L e = -e'''' - kappa^2 e'' - (5 kappa kappa' + q/4) e' + (1/4 - 3 kappa'^2 - 4 kappa kappa'') e,   e_s = L e + tau e', tau = int_0^sigma kappa V_Gamma, V_Gamma = -kappa'' + p/4.
- Sign: the e' coefficient from the tilt is -q/4 (not +). It is the Evans operator of step 2 verbatim.
- Transport: for an arclength-labelled non-stationary base the correction is tau = int kappa V_Gamma. For a quasi-static far curve this turns -q/4 into ~ -sigma/4 (checked analytically on a circle: tau = -sigma/4). For S1 tau = 0.
- If instead the base is labelled by the full flow Y_s = V N + (q/4) T (Lagrangian), tau = +q/4 and the drift cancels. That is a different labelling and was not used.

Validation:
1. Finite differences against the nonlinear V of the displaced curve, using exact displaced-curve formulas and 8th-order sigma differences: error scales as eps^2 and is 1.5e-5 (rel., eps=1e-3), 3.8e-6 (eps=5e-4), for S1, and the same for numerical curves at k=100 and 1600 (`f3fd.py`).
2. Galerkin spectrum of the (non-symmetrized) S1 operator (Hermite C^1 FEM, even class, Sigma=30, h=0.05 and 0.1, same to 1e-5): 1.00000, 0.25000, -0.25615, -0.69989, then -1.2145, -1.638. This reproduces the target 1, 1/4, -0.2562, -0.6999. Sigma=40 gives the same.

Counts (Sigma=35, free boundary, `logs/f3_*.log`). n = #eig of sym_rho above -c; M = smallest M such that the max eigenvalue of sym_rho on the complement of the first M S1-sym eigenvectors is <= -c (M = n in every case).
| Gamma | top sym eigenvalues | n,M (c=.05) | (c=.1) | (c=.25) | (c=.5) | (c=1) |
|---|---|---|---|---|---|---|
| S1 | 1.728, 0.2205, 0.0188, -0.239, -0.500 | 3 | 3 | 4 | 5 | 6 |
| k=1600 | 1.729, 0.2205, 0.019, -0.2385 | 3 | 3 | 4 | 5 | 6 |
| 400 | 1.732, 0.2207, 0.020, -0.237 | 3 | 3 | 4 | 5 | 6 |
| 100 | 1.744, 0.2217, 0.0236, -0.232 | 3 | 3 | 4 | 5 | 6 |
| 30 | 1.779, 0.2266, 0.036, -0.215 | 3 | 3 | 4 | 5 | 7 |
| 10 | 1.857, 0.2501, 0.074, -0.171 | 3 | 3 | 4 | 5 | 7 |
- Robust to: fitted frame (k >= 100, same counts); Sigma = 28 and 45 (same); omitting tau (same for k >= 400; k=100 gives 7 at c=1; k=30 gives n=11 at c=1 and 3/3/4/5 for c <= 0.5).
- Frame for k=10, 30: tip-anchored (lam = 1.5398/kappa_tip, x_c = tip + D_S1*lam); there is no S1 chart at those k.
- Margin: the 4th eigenvalue is -0.17 to -0.24, so c=0.1 holds for all k, but c can't exceed ~0.17 at k=10.
- H^2 variant (eps = 1e-5): gave numbers identical to the L^2 case to 4 digits. I did not verify that the code path actually changes the inner product, so it is unconfirmed. Treat the H^2 level as not done.

## Verdict against STEP4_DESIGN §5
1. F1: 300-500 DOF at k=1600 (r=0.25), ~480 at k=2900 (r=0.25), ~370 (r=0.5), ~305 (r=1), versus 2.6e5 uniform. DOF grow like ln k. Resolving theta'''' to 1e-12 needs roughly 1.5x more and cannot be measured below the ~1e-15 theta noise or the 1e-12 to 1e-13 spectral-filter floor of the snapshots.
2. F2: amplification modulo (T, translation, alpha) is <= 0.33 immediately and <= 7e-3 to k=1600 (rho and W), linear in delta, with no growth. This is far below 1e4. Early-k_j results are floor-limited (1e-10).
3. F3: tail dimension is 3 (c=0.05, 0.1) for S1 and for all k=10..1600; 4 at c=0.25; <= 7 at c=1. This is far below 60, and the "3 plus a few" check on S1 passes. The Lohner core could be as small as ~4-8 modes with the coercive tail beyond.
4. Decision rule: F2 <~ 1e4 and F3 <~ 60 hold throughout, including the early phase, so build the hybrid (Lohner core + coercive tail + weighted energies); no larger early block is indicated.
5. Caveats: single clip, one trajectory, float only; W-norm only on chart x <= 28 (not 260); k=10,30 use a tip-anchored frame and an unreliable chart; the H^2 variant is unverified; time-stepping error of snapshots was not independently measured.
