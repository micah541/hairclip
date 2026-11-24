# Hand-off feasibility measurement (non-rigorous): size of w = u - alpha psi3 - beta2 h2

Target: keyhole run k3_L24_ch0.3 (arms 24, slot 0.05, chamber 0.3). Norms are HALF-line (x >= 0), weight rho, on [0, min(28, chart_end-0.3)].
Check: ||psi3||_rho = 7.697, ||h2||_rho = 13.067, h2(0) = 0.2874 (match the stated values).

## Method (handoff2.py)
* u = normal graph over S1 of the rescaled curve. I replaced measure.graph by a cubic-Hermite intersection (exact tangents): the linear
  interpolation in measure.graph has a ~1e-4 systematic error, far above the 1e-6 target.
* measure.py integrates theta1 by the trapezoid rule (error up to 4e-6, floor ~1e-5 in u). I used a spline antiderivative instead. This changed r2 by less than 5%.
* (x_c, lam) and (alpha, beta2) are fitted jointly by variable projection (least_squares on (x_c, lam), linear LS for alpha, beta2), starting from the measure.py modulation.
  r1 = min ||u - alpha psi3||, r2 = ||u - alpha psi3 - beta2 h2||. The "C" columns have beta2 := alpha^2 and (x_c, lam, alpha) all refitted.
* Runs: run2.py (copy of the keyhole_run/run_hairpin logic, adjustable cdt, tail_tol, N0, resume from a snapshot). Extension (extA) resumed from snap_k3_L24_ch0.3_400
  with cdt 1e-2, tail_tol 1e-10 up to kmax = 6434. Resolution check: fineB (cdt 5e-3, tail_tol 1e-12, same N0 = 16384, from t = 0), continued as extB to 3200.
  repA is a from-scratch repeat of the original settings (it reproduces the existing snapshots bit-for-bit).

## Table (handoff2_k3_L24_ch0.3.txt; full-precision lines there)
| kmax | alpha | alpha*kmax | beta2/alpha^2 | r0 | r1 | r2 | r2/abs(alpha)^3 |
|---|---|---|---|---|---|---|---|
| 100 | -3.31e-2 | -3.31 | -0.19 | 2.4e-1 | 1.7e-3 | 1.41e-3 | 39 |
| 201 | -1.66e-2 | -3.34 | +0.40 | 1.2e-1 | 5.2e-4 | 1.82e-4 | 40 |
| 281 | -1.18e-2 | -3.33 | +0.58 | 8.1e-2 | 3.7e-4 | 8.8e-5 | 53 |
| 403 | -8.21e-3 | -3.31 | +0.77 | 5.6e-2 | 2.4e-4 | 6.7e-5 | 122 |
| 562 | -5.86e-3 | -3.29 | +1.03 | 4.0e-2 | 1.7e-4 | 7.8e-5 | 389 |
| 804 | -4.07e-3 | -3.28 | +1.65 | 2.8e-2 | 1.6e-4 | 1.07e-4 | 1.6e3 |
| 1124 | -2.91e-3 | -3.28 | +3.1 | 2.0e-2 | 1.9e-4 | 1.47e-4 | 6e3 |
| 1608 | -2.05e-3 | -3.30 | +7.4 | 1.4e-2 | 2.5e-4 | 2.10e-4 | 2.4e4 |
| 2248 | -1.50e-3 | -3.37 | +18 | 9.7e-3 | 3.4e-4 | 2.94e-4 | 8.7e4 |
| 3217 | -1.11e-3 | -3.56 | +46 | 6.9e-3 | 4.9e-4 | 4.20e-4 | 3.1e5 |
| 4496 | -8.7e-4 | -3.93 | +102 | 5.1e-3 | 6.8e-4 | 5.87e-4 | 8.8e5 |
| 6434 | -7.4e-4 | -4.76 | +203 | 3.9e-3 | 9.7e-4 | 8.40e-4 | 2.1e6 |

Chart end is 28 (the cap) for kmax >= 140 (22.4 at 100, 12.3 at 50). beta2/alpha^2 is not close to +1 at any single k (it passes through 1 near k = 560, but only because the unstable component below contaminates it).

## Main finding: w is not O(alpha^3). It contains a growing single mode, w ~ k_max^1.
* For kmax >= 560, r2 grows exactly like kmax (r2 ~ 1.3e-7 kmax), and the profile w(x)/kmax is the same function at kmax = 1120, 2240, 4480 (ratio 4.00 per doubling in kmax*w).
  It is NOT in span{p1, sin theta1} (symmetry modes; 1e-6 of its norm). I read it as the S1 unstable mode (eigenvalue +1/4 in -log(T-t)); its identity is not checked against an eigenfunction.
* It is physical, not numerical: the finer run (extB, cdt/2, tail_tol/100) shows the same growth (r2 = 3.5e-4 vs 4.2e-4 at kmax = 3200; the ~15% difference is a different unstable amplitude from the initial layer).
* All three chambers (0.2, 0.3, 0.4) have the SAME sign of this mode, with amplitude in units of ch0.3: 1.35, 1.00, 2.21 (unstable_bracket.txt). It is non-monotone in the chamber radius with a minimum
  ~1 near 0.27, so there is no sign change to bracket in this one-parameter family. Another parameter (arm, slot, ...) is needed to tune it away.
* Removing the fitted component c*phi (phi = w at kmax 1608) from w leaves a clean stable remainder: ||w_stable|| = 31 to 32 |alpha|^3 for kmax = 100, 200, 280, 400, 800
  (32.3, 29.4, 30.4, 31.6, 31.8; unstable_removed_k3_L24_ch0.3.txt). It agrees between coarse and fine runs (31.6 vs 31.6 at 400; 31.8 vs 32.2 at 800: numfloor_stable.txt).
  So w_stable ~ 3.2e-8 at |alpha| = 1e-3 and ~4e-9 at 5e-4.

## alpha*kmax and the hand-off points
* alpha*kmax = -3.31 (400), -3.28 (800 to 1100), -3.30 (1600) before the unstable mode contaminates it (drifts to -3.56 at 3200, -4.8 at 6400). Asymptotic (clean) value: about -3.3.
  Other chambers: -3.65 (ch0.2), -2.99 to -3.1 (ch0.4).
* |alpha| = 1e-3 at kmax ~ 3300, and |alpha| = 5e-4 at kmax ~ 6600 (using -3.3; measured: alpha = 1.11e-3 at 3217, 7.4e-4 at 6434, but contaminated).
* r2 as the run actually evolves: ~4.3e-4 at kmax 3300 and ~8.6e-4 at 6600 (linear extrapolation r2 ~ 1.3e-7 kmax; measured 4.2e-4 at 3217, 8.4e-4 at 6434).
  r2 never goes below 6e-5 (kmax ~ 400, |alpha| ~ 8e-3): r2 reaches its minimum there, because the stable part is falling (31 alpha^3) while the unstable part grows.
  r2/|alpha|^3 reaches 3e5 to 2e6, not O(1).

## Numerical floor
* Spatial/temporal resolution: halving cdt and tightening tail_tol by 100 changes alpha by ~0.5%, r2_stable by <3%, and the unstable amplitude by ~15% (initial-layer timing, tunable).
* The floor for the stable part is below ~1e-7 at kmax >= 800 (r2_stable = 2.2e-6 at 800 itself is a 31 alpha^3 signal, reproduced by both runs). Chart floor: S1 geometry (spline theta1) ~1e-7; Hermite interpolation ~1e-8.
  I did not demonstrate a floor below 1e-6 at |alpha| ~ 1e-3 directly, because the unstable mode dominates there.
* Caveat: runs from kmax 400 on are resumes from snapshots (the first step ramps from a small dt); extA and extB differ in unstable amplitude as above.

## Verdict
Taken literally, ||w||_rho ~ 1e-6 at |alpha| ~ 5e-4 to 1e-3 is NOT realistic for the untuned target hairpin: the actual r2 there is ~4e-4 to 9e-4 (about 400 to 900 times too large), because the trajectory
leaves S1 along an unstable mode (amplitude ~1.3e-7 kmax in rho-norm). It is feasible only if that unstable coefficient is tuned (shooting in a second parameter) down to ~1e-3 of its current value
(a_u/kmax below ~3e-10 for 1e-6 at kmax 3300). Once it is tuned out, the stable remainder is only about 31|alpha|^3 (3e-8 at 1e-3), so the 1e-6 requirement would then be easily met and the accuracy
is limited by how well the unstable component can be cancelled. Note the constant 31 means w_stable/alpha^3 is large, so the h2 truncation is fine but a third-order term would be needed for better than ~3e-8.
(All of this is a non-rigorous numerical estimate. The identification of the growing mode as the unstable eigenmode is by scaling and shape only.)

## ADDENDUM (oddtest, supersedes the Verdict above)
The "unstable mode" was a frame artefact: the tip sits at y = +3.477e-7 (initial-curve asymmetry, closure defect ~1.7e-7; the solver preserves the symmetry to 1e-13 in theta).
The odd part (u_up-u_low)/2 matches cos(theta1) to 0.3-3% (oddtest/odd_k3_L24_ch0.3_all.txt). With the symmetry axis through the tip (now the default in handoff2.py; old frame in old_frame/),
r2 = 32-46 |alpha|^3 throughout (table in handoff2_k3_L24_ch0.3.txt): r2 = 3.6e-7 at k=1608 (alpha=-1.99e-3), 1.3e-7 at 2248, 5.2e-8 at 3217. alpha*k = -3.30 (k=400) -> -3.20 (1608) -> -3.14 (3217).
Hand-off criterion (|alpha|<=2e-3, r2<=3e-7) holds from k* ~ 1700. The symmetrized run (3) was not needed. Remaining odd contamination (rotation) can still inflate r2 if the initial closure defect is large (see design/evenfit.py).
