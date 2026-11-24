# FEAS2 report: step 4 macro phase (float only)

Setup: clip run, T = 4.5988639 (varies ~1e-5 with the step sequence, see Floors). HLS/ETDRK4, cdt = 1e-2. Rescaled tip frame, lambda = 1.5398/kmax. 2 processes max. Scripts and logs are in feas2/ and feas2/logs.

## F2b sensitivity (delta = 1e-6, 2e-6; modulo T-shift, translation, alpha)
A_rel = dw/delta, where dw is the residual w-difference after the chart fit, normalised per delta. The full table, including the W-norm and the physical amplification, is in logs/f2tab.out. Range over the shapes (peak/fixed/far), by injection time t_j:

| t_j | A_rho k=100 | k=400 | k=1600 | worst shape | dalpha/delta k=1600 |
|---|---|---|---|---|---|
| 0 | 1.5e-3 .. 2.4e-2 | 3.6e-5 .. 4.2e-4 | 3.7e-6 .. 3.3e-5 | peak (100), fixed (400, 1600) | 4e-5 (peak) |
| 0.3 | 1.2e-2 .. 1.4e-1 | 2.1e-4 .. 5.3e-3 | 6.4e-5 .. 4.8e-4 | fixed | 1.2e-3 (peak) |
| 1.5 | 2.1e-2 .. 1.0e-1 | 8.2e-4 .. 4.2e-3 | 2.9e-5 .. 1.2e-3 | fixed | 2e-4 |
| 3.0 | 3.2e-2 .. 4.0e-1 | 2.4e-4 .. 5.5e-3 | 6.8e-5 .. 5.4e-4 | fixed | 3.8e-3 (peak) |
| 4.5 | 8.5e-3 .. 6.1e-2 | 5.6e-4 .. 4.5e-3 | 1.1e-4 .. 1.1e-3 | fixed | 1.3e-3 (peak) |

- The W-norm values are within about 2x of the rho values (see logs/f2tab.out).
- No growth with k in any case: A drops by roughly 10^2 from k=100 to k=1600. This is the expected stable-manifold contraction after the T-shift and alpha are removed.
- The residual is largest for localized fixed-width bumps (w = 0.05) injected at t >= 0.3, and largest at t = 3.0 (A = 0.40 at k=100). That is the moment the 5th direction is nearly marginal (F3b).
- Linearity check (d2/d1)/2 is 0.3-1.7 at most entries; some far/k=1600 entries reach 2-4, so the nonlinear and floor noise there is about 1e-6 relative. The k=1600 numbers have the weakest linearity.
- Crossing events (bracketed by a 1.5e-3 scan; refine to 1e-3 needs bisection, not done):

| event | t bracket | ncross | kmax |
|---|---|---|---|
| first crossing | 0.00281 .. 0.00313 | 0 -> 2 | 3.2 |
| second pair appears | 0.01563 .. 0.01594 | 2 -> 4 | 3.4 |
| extra pair disappears | 2.64250 .. 2.64281 | 4 -> 2 | 2.12 |
| end of scan (t = 4.5985) | - | 2 persist | 11.1 |

The 2 remaining crossings persist through the whole rapid phase, so "embedded" holds only in the sense of the final loop geometry.

## F3b tail count of sym_rho(L_Gamma + transport)
Variant A is the actual frame (lambda = 1.5398/kmax, rho tied to mu = 1/4). Eigenvalues are of the form on the full half-line, not on the S1 complement.

| t | kmax | top eig (3rd..6th) | n(>-0.05) | n(>-0.1) | consistent-mu variant E n(>-0.05) |
|---|---|---|---|---|---|
| 0.5 | 2.70 | .189, -.048, -.291, -.539 | 4 | 4 | 15 |
| 2 | 2.20 | .186, -.052, -.276, -.439 | 3 | 4 | 15 |
| 3 | 1.99 | .278, .091, .022, -.367 | 5 | 5 | 16 |
| 4 | 1.84 | .275, .101, -.218, -.440 | 4 | 4 | 6 |
| 4.5 | 2.77 | .196, -.007, -.266, -.537 | 4 | 4 | 4 |
| 4.59 | 5.03 | .120, -.118, -.357, -.607 | 3 | 3 | 3 |

- Frame-A counts: 4 (5 near t = 3, 3 at t = 4.59). S1 itself has 3 (the 3 cutoff directions).
- The extra directions are the 4th and 5th modes. Their rho-weight sits at s ~ 9-12 (s50/s90 between 8.9 and 14.7), which is in the slot-wall zone and beyond. The axis crossings are at s = 2.5-4 and 5-10 (t >= 2.6), so those modes live around and beyond the second crossing and outside the chamber.
- Variant E (frame-consistent mu, wc) gives 10-16 non-coercive directions for t <= 3.5: yes, many more. Most of that is the rho-vs-mu mismatch, since mu goes from -0.09 to 0.2 against the design value 1/4. Variant E is 4 only for t >= 4.25 (mu approaches 1/4).
- T-clock variant F (fixed centre) gives 16, 11, 10 at t = 0.5, 2, 3, then 4, 4, 3. So the 3-direction count is not robust to the frame choice before t about 4.25.
- Fixed-frame rho is not the right weight for the macro phase. A frame-adapted weight is needed there, and the 4/5 count must be redone with it.

## P1 float Galerkin + tail certificate (omega = 0.05, Sigma = 100, FEM h = 0.1)
| case | ncon | N closes | margin at N_close | N for margin >= 0.05 | limiter |
|---|---|---|---|---|---|
| k = 1600 | 3 | 3 | 0.144 | 3 | tail nu_4 - omega |
| k = 100 | 3 | 6 | 0.18 | 6 | tail vs sup-norm penalty |
| k = 10 | 3 | 29 | 0.02-0.08 | 29-30 | tail |
| t = 4.59 (k = 5) | 3 / 4 | 44-47 / 43 | 0.004-0.049 / 0.048-0.15 | 45-69 / 43-44 | cross-term/Galerkin block (3) / tail (4) |
| t = 4.5 | 4 | 54-57 | 0.001-0.054 | 54-60 | cross-term block |
| t = 2.0 | 4 | 76 (cap 80) | 0.07-0.29 | 76 | cross-term block |
| t = 0.5 | 4 | 76 (cap 80) | 0.05-0.25 | 76 | cross-term block |

- The two entries per row are Rc taper 20 and 10 (range). Margins are float values on a bounded nu grid.
- Cost in the macro phase is driven by the sup-norm coefficient difference to S1 (d0 about 4, d2 about 0.7 for t <= 2). The certificate is expensive but closes.

## P2 spacing and certificate count
Budget theta <= m/(omega+m) = 0.2857 for margin 0.02.

| phase | sum theta | sum ds | certificates (sum theta / b; greedy) |
|---|---|---|---|
| layer t < 0.01 | 1.4e7 | 0.002 | not coverable (see below) |
| macro 0.01 < t < 4.55 | 12 | 4.54 | 42; 39 |
| approach k 3.3 .. 7 | 0.58 | 3.0 | 2 |
| rapid k 7 .. 1810 | 0.60 | 22.2 | 2.1; 3 |
| rapid k 1810 .. 6600 (new) | 0.0013 | 4.4 | within the last k >= 1600 certificate |

Allowed Delta s for margin 0.02 (full 40 rows: logs/p2_rows.json):

| k (or t) | Delta s allowed |
|---|---|
| t = 0.01..0.02 | 0.007 |
| t = 0.05..0.3 | 1.4-3.5 |
| t = 0.3..1.5 (ds is 0.06-0.15 steps) | 0.06-0.11 |
| t = 3.0..3.5 | 0.017 (worst) |
| t = 4.0..4.58 | 0.2-1.5 |
| k = 7 / 10 / 14 / 20 | 1.2 / 2.1 / 3.2 / 4.5 |
| k = 40 / 160 / 640 / 1810 | 9.4 / 41 / 175 / 510 |

- theta/ds is about 1.2 k^-1.05 for k > 30 (a fit).
- Extension to k about 6600 (Step 3 hand-off): theta for 1810->2560->3620->5120->6600 is 5.4e-4, 3.8e-4, 2.7e-4, 1.4e-4, which is far under the 0.2857 budget. One certificate with N = 3 covers all k >= 1600 (coercivity of the complement stays 0.239).
- Total from t = 0.01 to k = 6600: about 42-45 certificates (39-42 macro + 2 approach + 3 rapid), about 48 once the layer is excluded and a safety factor is included.
- Layer t < 0.01: the theta values are huge (coercivity at t = 0 is -5.5e5, the state is the original clip curve with kmax = 42). A dedicated initial-layer argument is needed; this study does not cover it.

## Floors
- T varies by about 1e-5 with the step sequence (cdt 1e-2 vs 5e-3: dT = -8.5e-6 at k = 100 and 400). After the T-shift is removed the chart residuals are dw about 1e-8 (rho, k = 100) and 5e-10 (k = 400). The cdt5e-3 test at k = 1600 was not run ("missing 1600"). Everything with delta = 1e-6 at k = 1600 is therefore near the floor, which explains the linearity scatter.
- Time floor at k >= 2560: the state variable t repeats (4.598865035383) while kmax grows, because dt = cdt * lambda^4 falls below 1e-12 against t about 4.6. Step 4 must integrate in tau = T - t (or in s) and not in t.
- FEM (h = 0.1) and Sigma: Sigma = 35 gave spurious eigenvalues; Sigma = 100 is needed for the dense S1 half-line spectrum above nu about 3.4. N_close values depend on this.
- F2b: fixed_2 runs were skipped (the fixed shape is only checked at delta = 1e-6). Nothing here is validated, float only.
- Shared machine: all jobs ran with at most 2 processes. t0_far_2 was killed once and rerun.

## Verdict
1. Cost: macro phase about 42 certificates (N = 44-76 plus 3-4 constraints); the approach and rapid phases to k = 6600 only 5 more (N 47 -> 29 -> 6 -> 3). The hand-off extension to k about 6600 costs nothing extra. The initial layer t < 0.01 is not coverable by this scheme.
2. Sensitivity is benign: after removing T-shift, translation and alpha, amplification falls about 100x from k = 100 to 1600 for all injection times, with the worst case at t = 3.0 (0.40 at k = 100). Crossings: 0 -> 2 at t about 0.003, 2 -> 4 at 0.016, 4 -> 2 at 2.6427; 2 persist.
3. The "3 modulated directions" plan survives only for t >= about 4.55 (k >= 5). Before that, 4 directions are needed (5 near t = 3, where a mode crosses the threshold), with rho-weight in the slot-wall zone s about 9-12.
4. A frame-consistent count (E, F) gives 10-16 non-coercive directions for t <= 3.5, mostly weight/mu mismatch; stated plainly, the fixed S1 weight does not support shadowing in the macro phase, so a time-dependent weight is needed.
5. Floors: T uncertainty about 1e-5, t cannot be used as a clock beyond k about 2500, and the k = 1600 cdt check and all validation remain open; everything here is float.
