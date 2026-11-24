# v2 intermediate constants: new design, P = 2 (primary) / 1.5 (backup), region [10.25, R]

Scripts here (copies adapted from ../): int_common.py (weight parameter P, psi rigorous by acb.integral, W_up/W_lo, Taylor-shift S^(j)),
int_n3.py (V1), int_res.py (V2), int_nlin.py (V3, q/4 factoring kept), int_garding.py (V4, `[width] [P]`), int_v56.py (V5, V6), collate_v2.py -> int_constants_v2.json.
|alpha|<=1e-3, jet box 1e-2, half-line W-norms, all sigma >= 8.0009 (x >= 10.25). Outward-rounded, NaN-asserted. Per-region data in the json.
V1, V2 are plain W-norms (no zeta). V5 norms include zeta^2 (box sup of zeta; inner ramp S((x-10.25)/6), outer ramp S((R-x)/3)); zeta-free norms also in json.

| | P=2 R=40 | P=2 R=60 | P=1.5 R=40 | P=1.5 R=60 |
|---|---|---|---|---|
| V1 C3W | 0.1823 | 0.1823 | 0.6885 | 0.6885 |
| V1 sup g'''/6 | 12.26 | 12.26 | 12.26 | 12.26 |
| V2 psi3 residual | 1.092e-7 | 1.092e-7 | 3.541e-7 | 3.541e-7 |
| V2 h2 residual | 6.469e-9 | 6.469e-9 | 3.136e-8 | 3.136e-8 |
| V4 sup(-mu2) / sup(-mu1) | -0.6555 / -0.1555 | same | -0.6242 / -0.1242 | same |
| V4 sup\|e1\| / \|e0\| | 5.66e-6 / 5.11e-6 | same | same | same |
| V5 psi1' / psi2' | 1.7e-7 / 1.7e-8 | same | 6.0e-7 / 5.9e-8 | same |
| V5 psi3^' / h2^' | 0.2595 / 1.795 | 0.2596 / 1.796 | 1.204 / 8.99 | 1.204 / 9.00 |
| V5 A1'/alpha, B1'/alpha | 0.273, 0.0159 | same | 1.269, 0.0684 | 1.270, 0.0684 |
| V5 zeta c0' | 9.9e-12 | same | 3.4e-11 | same |
| V6 W(R-3), W(R), int W | 3.59e-8, 1.23e-9, 3.11e-8 | 1.34e-18, 2.76e-20, 1.02e-18 | 1.23e-6, 4.22e-8, 1.065e-6 | 4.59e-17, 9.46e-19, 3.49e-17 |

V3 (P-independent), sup over [10.25, R=40 / 60]: w'''' coefficient (1-J^-4) 8.44e-3 / 1.85e-2; w''' 9.8e-4 / 1.5e-3; w'' 1.13e-4 (on [10.25,16.25], 2.5e-5 beyond); w''^2 0.014 / 0.015;
w''w''' 0.654 / 0.977; w''^3 3.12 / 3.26; A0 1.69e-5; A1 7.10e-5 (7.5e-8 beyond 16.25); sup|u0'| 0.0196 on [10.25,16.25], 0.0250 on [16.25,20], 0.0546 to 40, 0.0850 to 60;
sup|u0''| 1.54e-3; sup|u0'''| 5.25e-6 (3.5e-6 beyond 16.25). Per-region values in the json. V5 sup per region: |kq| 1.3e-3 on [10.25,16.25]; A1'/alpha 20 -> 92; B1'/alpha 1.4 -> 1.5.

V4 mu1 check vs lint_cert_p.py: P=2 min 0.1555 at x=10.25 (cert 0.1552), blend min 0.3358 (cert 0.3357); P=1.5 0.1243 (0.1241), blend 0.2211 (0.2211).
mu2 = mu1 + 1/2 up to profile terms again: H^2 level has a dissipation margin (mu2 >= 0.624 for P=1.5, 0.656 for P=2).

V2 per-region (P=2, R=60): psi3 [10.25,16.25] 1.09e-7, [16.25,20] 1.7e-9, [20,28] 3.7e-9, [28,30.25] 1.9e-11, [30.25,40] 1.2e-11, [40,60] 2.5e-13;
h2 3.0e-9, 7.5e-10, 5.7e-9, 9.5e-11, 1.2e-10, 3.3e-12.

## Caveats
1. psi3 residual is dominated by x in [10.25, 12] and uses the crude ray + chunk bounds there (weighted parts at x = 10.25: ray 9e-8, profile terms 1.5e-7).
   A tight-tail Taylor stage for psi3 (as used for h2) would likely cut it by 10-100x; not done. h2 already uses the tight tail for sigma <= 20.
2. V1 is the same quantity as T1 (||N0(u0) - alpha^2 Q2||_W), NOT its x-derivative: the derivative-level forcing needs u^(5) of psi3^/h2^ (ext_jet supports n<=5 only) and was not computed.
3. psi2' = -k cos(theta): |cos theta| <= 1 was used (no frame data); B1 uses |cos|,|sin| <= 1, theta' = k. ratio A1'/alpha, B1'/alpha: sup over |alpha| <= 1e-3 via u0^(j) = A * (psi3^(j) + A h2^(j)).
4. Profile data are the step-1 tail bounds profile_bounds(sigma_lo) per box (kb ~ 1.7e-4 at sigma = 8); norms involving k, m, p (psi1', psi2', c0') are upper bounds, not sharp.
5. Tip ST enclosure used; box x-ranges start at >= 10.25 and end at R + 1e-9. T2 values depend mildly on box width (conservative); V1 at WS=0.25, V3 at 0.5, V4 at 0.25, V5 at 1.
6. Float checks (v2/float_check_*.py): V1 c3 at x = 10.3, 11, 16 agree to <= 3e-7 relative; h2 full residual float <= rigorous box bound at x = 10.3, 11, 14, float W-norm on [10.3,20] 1.4e-9 <= rigorous 6.5e-9 total.
