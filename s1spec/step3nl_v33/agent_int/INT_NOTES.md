# Intermediate region [12, R], weight W = exp(-psi), R = 40 and 60  (agent_int/)

All numbers: python-flint ball arithmetic, upper bounds rounded outward (int_constants.json, 4 digits).
|alpha| <= 1e-3, jet box |w|,|w'| <= 1e-2, half-line norms ||f||_W^2 = int_12^R W f^2, W = (1+x^2)^-3 for x <= 30.25.
Common code: int_common.py (psi = integral of psi' by acb.integral, psi(30.25) = 20.4603, psi(60) = 52.1018 [matches agent_sob K_W];
W_up(x_lo) = exp(-psi(x_lo)) per box since W decreasing; psi', psi'', psi''', psi'''' on balls with exact clipping at 30.25 / 38.25;
sigma-boxes whose x-range starts at >= 12 rigorously and ends at R + 1e-9; profile = HS.profile_bounds(sigma_lo) on every box, all sigma).

| key | script | bounds |
|---|---|---|
| T1_C3W (R=40 / 60) = 8.834e-3 / 8.834e-3 | int_n3.py | ||N0(a psi3 + a^2 h2) - a^2 Q2(psi3)||_{W,[12,R]} <= C3W |a|^3 (series at exact alpha=0 for c3,c4; interval A only for c5) |
| T1_sup_g3over6 = 10.00 / 10.32 | int_n3.py | sup_[12,R] \|g'''/6\| (c3 at alpha = 0); per region 9.78, 9.78, 10.00, 10.32 |
| T2_psi3_res_W = 2.619e-10 / 2.619e-10 | int_res.py | ||(L - lamt) psi3||_{W,[12,R]} |
| T2_h2_res_W = 2.351e-10 / 2.351e-10 | int_res.py | ||(L - 2lamt) h2 + Q2(psi3)||_{W,[12,R]} (sigma<=20: tight-tail Taylor form; beyond: ray+chunk+Q2) |
| T3_sup (see json, per region + cumulative) | int_nlin.py | sup \|1-J^-4\| (w'''' coeff) 8.436e-3 / 1.847e-2; w''' 9.82e-4 / 1.503e-3; w'' 3.23e-5 / 3.23e-5; w''^2 1.403e-2 / 1.502e-2; w''w''' 0.6543 / 0.9768; w''^3 3.115 / 3.255; A0 2.009e-6; A1 8.593e-6 |
| T3 chart | int_nlin.py | sup\|u0'\| 0.05461 / 0.08501; sup\|u0''\| 1.502e-3 / 1.539e-3; sup\|u0\| 1.014 / 2.409 |
| T4 Gaarding | int_garding.py | sup(-mu2) = -0.7440 (both R; i.e. mu2 >= 0.744), sup\|e1\| 6.307e-8, sup\|e0\| 6.056e-8 (all from [12,20]; beyond 1e-18 and smaller) |

Per-region values: int_constants.json (T1_per_region_*, T2_*_per_region_norm, T3_per_region, T4_per_region).
Regions: [12,20], [20,30.25], [30.25,40], [40,60]; T2 regions additionally split at x = 28 (end of the polynomial pieces / start of the extension).

## T2 per-region contributions (norms, W)
psi3: [12,20] 2.17e-10, [20,28] 1.47e-10, [28,30.25] 6.5e-13, [30.25,40] 3.6e-13, [40,60] 7.2e-15.
h2:   [12,20] 1.04e-10, [20,28] 2.11e-10, [28,30.25] 3.3e-12, [30.25,40] 3.8e-12, [40,60] 9.7e-14.
Dominant terms: piece-fit error of the polynomial elements (ray residual) near x = 12 and in the last element just below x = 28
(pointwise weighted max 4.8e-10 at x = 28.0 for h2); for x >= 28 the explicit-extension residual is negligible and the bound is
governed by the tip-uncertainty term dST/4 |u'| (dST = |ST - STIP| <= 5.5e-10).

## Surprises / remarks
1. W-weighted residuals are ~40x smaller than the rho ones restricted to x >= 12 (rho-parity run, int_T2_rho_1.json: 9.9e-9 / 8.6e-9);
   the extension residuals do NOT matter. BUT the unweighted rigorous pointwise bound grows (tip term: h2 2.2e-7, psi3 1.1e-8 at x = 59);
   it is only the factor sqrt(W) ~ x^-3 that hides it. For a slower-decaying weight the dST term would need a sharper tip enclosure.
2. T1: the pointwise coefficient |g'''/6| is essentially flat (9.78 -> 10.32 on [12,60]); the tail [40,60] adds only 1.4e-10 to ||.||^2.
3. T3: A1 in nlin_E2.py (5e-3) is a DEPENDENCY ARTEFACT: dU1 contains -q/(4 alpha_g) + q/4 which cancels identically on the ray, but
   as balls it leaves width(q)/4 (box 0.02 -> 5e-3, 0.05 -> 1.25e-2). Fixed here by splitting dU1 = dU1|_{q=0} + q * c_q with
   c_q = k U0 (...)/(4 (U0 k - 1)(...)^3) factored exactly (identity checked numerically, affinity in q asserted). A1: 8.6e-6 on [12,20],
   5.3e-8 beyond. The old E2 constants A1 = 5e-3 can be improved by the same trick (not done: existing files untouched).
4. T3 chart: sup|u0'| = 0.085 on [40,60] (0.055 on [12,40]) at alpha = 1e-3; h2 part is 1.7e-3 of it. R = 60 leaves little room
   for alpha > 1e-3 (|u0'| scales linearly in alpha psi3', psi3' ~ 83 at x = 60).
5. T4: mu2 = mu1 + 1/2 up to profile terms (q' = 1 enters c0t as 2 D(a1) = -1/2), so the H^2 level has a dissipation MARGIN, no Gaarding
   growth; mu1 check reproduces lint_cert (min 0.2435 at x = 12, 0.5625 on the blend, 1.1 on [40,60]). e1, e0 are pure profile terms
   (7k^2 q - 37mp/2 ; 5k^3p/2 + 5k^2/2 + 23kmq/2 - 13p^2/8): negligible.
6. Pitfall hit: S^(j)((x-30.25)/8) in expanded form gives +-1.4 on a width-0.02 ball near t = 1 (cancellation), which made mu1 look like -12 on
   [38,38.25] at width 0.02; int_common evaluates S^(j) by exact Taylor shift at the box centre (fixed; widths <= 0.01 used anyway).
7. T2 depends mildly on box width (WS=1: 2.85e-10 / 3.81e-10; 0.5: 2.69e-10 / 2.75e-10; 0.25: 2.62e-10 / 2.35e-10): values are conservative, not converged limits.

## Validation (float, mpmath 40 digits)
T1 float_check_T1.py: finite-difference c3 vs arb series at x = 12, 20, 35, 59.9: agree to <= 5e-6 relative (9.05034 / 9.45785 / 9.89538 / 10.3173).
T2 float_check_T2.py, _T2b.py: pointwise float residuals (ray; and full residual with tight-tail profile for h2) are below the rigorous box bounds at all
   tested x; float W-norm of the h2 residual on [12,20] 7.6e-11 vs rigorous 1.0e-10 (region [12,20], final run); psi3 ray-only float W-norm 2.05e-10 on [12,40] vs rigorous 2.6e-10.
T3 float_check_T3.py: Taylor coefficients of D at x = 45: 1-J^-4 = 7.677e-3, w''' 9.28e-4, w'' 2.10e-5, w''^2 1.313e-2, w''^3 2.897, w''w''' 0.614 (all <= region sups).
T4 float_check_T4.py: float mu1, mu2 at x = 12..59 lie inside the rigorous ball values; mu2 - mu1 = 0.5.
rho-parity: int_res.py ... rho reproduces the order of the stage-B parts of the original rho-certificates for x >= 12 (those totals include x < 12).

## Reproduce
python3 int_n3.py 1e-3 0.25 ; python3 int_res.py 0.25 W ; python3 int_nlin.py 1e-3 1e-2 1e-2 0.5 ; python3 int_garding.py 0.25 ; python3 collate.py
(int_nlin.py spends ~70 s in sympy at import; int_res.py ~2 min.)
