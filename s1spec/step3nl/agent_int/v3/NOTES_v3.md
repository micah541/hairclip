# v3: derivative-level constants (P = 2, region [10.25, R], R = 40 / 60), |alpha| <= 1e-3

## KEY FINDING: the certified psi3^, h2^ are only C^3 -> derivative-level quantities are infinite in L^2
At each of the 20 interfaces psi3^ and h2^ have a jump in u'''' (1e-7 .. 4e-6; checked). Hence (L psi3^)', d/dx N0(u0) contain delta masses, and
the level-2 energy (needs v'''' = w^(5) in L^2) is not even defined for w = u - u0. Fix used here (what I did): `c4fun.py` builds C^4 versions
psi3~, h2~ = same Chebyshev pieces and same explicit extension, but degree-9 Hermite corrections removing the jet mismatches of orders 0..4
(and u'(0) = u'''(0) = 0); jets are exact arb enclosures of ANY order on the polynomial pieces, and of any order on the extension (Bell recursion for (e^S)^(d), verified
against ext_jet to 1e-55). C^4 verified: max jump of jets 0..4 at all interfaces and at x = 28 is < 1e-55. psi3~ - psi3^ is ~1e-7 relative (jets 0..4).
All F-constants below are for psi3~, h2~ (u^(5) is piecewise bounded; no Cauchy-strip estimate needed). If psi3~/h2~ replace psi3^/h2^, the level-0 constants change to the
values below (level-0 for C^4): ||(L-lamt)psi3~||_W 1.243e-8 (was 1.09e-7: now with the tight-tail Taylor stage for psi3 too), ||(L-2lamt)h2~+Q2||_W 4.03e-8 (was 6.5e-9: the jet-4
corrections add ~1e-6 near interfaces, mostly the last element before x = 28), C3W (level 0) 0.1823 (unchanged).
Every other earlier constant (T3/V3 coefficients, mu's, V5 norms) is insensitive (change ~1e-7 relative) but formally refers to C^3 functions.

## Results (identical for R = 40 and 60 to the digits shown unless stated; zeta = box sup)
| | value |
|---|---|
| F1a  C3W' = ||zeta d/dx[N0(u0) - a^2 Q2]||_W / |a|^3 | 2.603e-3 (nozeta: 4.50e-2) ; per region [10.25,16.25] 2.47e-3, [16.25,20] 5.96e-4, [20,30.25] 5.26e-4, [30.25,40] 1.90e-4, [40,60] 3.6e-6 |
| F1a  sup of pointwise \|c3'\| | 6.43 (at x ~ 10.25, where zeta = 0), 0.064 / 0.117 / 0.026 / 0.020 on the other regions |
| F1b  ||zeta ((L-lamt)psi3~)'||_W | 2.70e-8 (per region 4.8e-9, 1.7e-9, 2.65e-8 [20,28], 2.1e-11, 5.8e-13, 6.6e-15) |
| F1b  ||zeta ((L-2lamt)h2~+Q2)'||_W | 1.156e-7 (1.39e-8, 5.0e-9, 1.15e-7 [20,28], 3.1e-11, 8.4e-12, 1.8e-13) |
| F2  sup \|1 - q psi'\| | 42.04 (R=40), 74.37 (R=60) |
| F2  sup \|(W cos th)'/W\| <= sup sqrt(psi'^2 + k^2) | 1.140 (R=40), 1.306 (R=60) |
| F2  sup_ramp q zeta' W/rho | 0.4357 (sup zeta' 0.3665, sup W/rho 0.4313 on the ramp) |
| F3  ||zeta A1(u0)'''||_W / alpha | 7.23e-5 (nozeta 1.46e-4) |
| F3  ||zeta B1(u0)'''||_W / alpha | 2.07e-6 (nozeta 2.33e-5) |
| F3  ||zeta psi3~'''||_W , ||zeta h2~'''||_W | 2.63e-5 , 1.12e-2 |
Pointwise sups of the F3 quantities per region are in int_constants_v3.json (A1'''/alpha <= 0.0123, B1'''/alpha <= 0.0039, psi3''' <= 4.3e-3, h2''' <= 1.09).
Level-1 residual weighted pointwise max 8e-8 (psi) / 4e-7 (h2) at x ~ 28 (end of the polynomial pieces).

## Method summary
F1a: dual numbers over arb_series in alpha (forward-mode d/dx: u_j' = u_{j+1}, profile ' from the ODE), series at the exact point 0 for c3', c4', interval A only for c5'
 (same remainder structure as T1); sanity: c0' = c1' = 0, (c2)' = DQ2 (checked symbolically-numerically).
F1b: sympy total derivatives R, R', R'', R''' of the residual (Q2 from q2_derivs.pkl). sigma <= 20: tight-tail profile at the box centre; level 0: R + R'h + R''(box)h^2/2,
 level 1: R' + R''(xc)h + R'''(box)h^2/2. sigma >= 20: ray residual (Taylor form on pieces; closed forms R g = -4 a4 s^(A-16) and the WKB term, both differentiated
 exactly, on the extension) + tip term dST/4|u''| + differentiated profile terms + |DQ2| (h2).
F3: Taylor series in the x-increment (order 3) of u0 = A*yh, profile ODE jets, al_g, cos(theta) = c0 cos Theta - s0 sin Theta, |c0|^2+|s0|^2 = 1.
F2: Pj_box (Taylor-shift S^(j)); W/rho monotone on the ramp certified per box (phi' - psi' > 0).

## Caveats
* The psi3~, h2~ are new functions (not the certified C^3 objects); their derivation (jets, constraints inner products <psi3, f_j>, G, lyapunov constants) was not re-run: effect ~1e-7 relative.
* h2~ has jets of orders 4, 5 that are NOT small near interfaces (jump of u^(5) 1e-5..1e-3, harmless in L^2 but it is why the [20,28] residual region dominates).
* F2's (W cos th)'/W uses |cos|,|sin| <= 1 (sqrt(psi'^2 + k^2) bound, no frame data); B1 uses theta' = k.
* Box widths: F1a 0.25 (x 0.005..0.0125), F1b 0.25, F3/F2 0.5. Residual bounds are conservative (not converged limits). The residual stage Bt has a floor from the step-1 J-bound in q (1.8e-8 at sigma = 8).
* Float checks (float_check_*.py): F1a d/dx c3 by finite differences vs dual series, rel. diff <= 1e-5 (x = 12, 15, 25.5, 35); F1b symbolic R' vs 300-bit finite differences agree to 6 digits for x >= 14
  (at x = 10.5 FD is polluted by the midpoint of the q-enclosure; the symbolic value lies below the rigorous box bound everywhere tested); F3: A1', B1' reproduce the v2 V5 constants (0.2728, 0.01588),
  A1''' vs the ray formula 2 yh3 + sigma yh4 (float <= rigorous).
