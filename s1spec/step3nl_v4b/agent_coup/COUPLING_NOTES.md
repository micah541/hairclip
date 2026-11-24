# Rigorous far-field coupling constants (agent_coup)

Files: `coup_common.py` (ball-arithmetic library), `run_coup.py PX out.json` (all of Q1-Q6 for one P), `verify_identity.py` (+`.log`),
`float_ref.py` (independent float reference), `consolidate.py`, `coup_P2.json`, `coup_P1.5.json`, `coupling_constants.json`.

**Identities.** Verified symbolically (`verify_identity.py`): the Euler operator in w of
(-rho chi w [L,chi]w) - (A w^2 + B w'^2) is identically 0, i.e. the difference is a total derivative, for generic rho, chi, k, and
drift coefficient. The stated A, B are correct as written. The same formulas with 7km, W, zeta (and zeta_out) are used for the ramps.
The compiled polynomial forms A/Wt, B/Wt (Wt = rho or W, using P_j = phi_I^(j), ODE rules k'=m, m'=p/4, p'=-kq, q'=1+kp) match the
generic A, B to double rounding at random points. A, B vanish outside [29,32] (resp. the ramp), because every term carries a derivative of the cutoff.

**Method.** Pieces are float intervals with exact endpoints, split at 30.25 and 38.25 where needed. Smoothstep derivatives use Taylor-shifted
S^(j) (stable). phi_I is tabulated at 2^-9 nodes on [30.25,38.25] by rigorous `acb.integral`, with a mean-value correction in between
(closed forms elsewhere). The tip is the interval hull of the given enclosure: strip pieces are in sigma (x = sigma - tip as a ball), ramp
and outer pieces are in x (sigma = x + tip). k, m, p are +/- radii and q = sigma + [-qs, qs], from `hstail.profile_bounds` (the step-1 tail
bound, valid for sigma >= 8; the ramp starts at sigma = 8.0009). sup = branch and bound (all upper bounds from piece enclosures);
integrals = sum of length * enclosure with adaptive bisection (rel. 1e-3); Phi(x) = int_{16.25}^x 1/W is bounded per 2^-9 segment by
e^{phi_I(g)}(e^{Mh}-1)/M, M = sup phi_I'. mu_max is a rigorous lower bound: the largest lambda such that
(1/9 - lambda b)t^2 - lambda a t + (om0 - lambda a c1) >= 0 for all t >= 0 (checked by discriminant in arb, with a, b the upper bounds). Runtime about 50 s per P (75k pieces).

**Results** (upper bounds; float = independent numpy/sympy or the design scripts). Rigorous/float agree to 0.03-0.07% everywhere, rigorous above float.

| quantity | P=2 rig | P=2 float | P=1.5 rig | P=1.5 float |
|---|---|---|---|---|
| sB | 3.5955e-5 | 3.5930e-5 | 1.1101e-6 | 1.1094e-6 |
| IA | 9.5997e-12 | 9.5971e-12 | same | same |
| sA | 1.3372e-4 | 1.3368e-4 | 4.9844e-6 | 4.9831e-6 |
| K_rho(16.25) | 3.0208e4 | 3.0208e4 | same | same |
| Phi(34.25) | 9.8547e6 | 9.8547e6 | 3.4672e5 | 3.4672e5 |
| a_f / b_f | 6.9555e-2 / 3.9359e-2 | 6.9509e-2 / 3.9333e-2 | 1.0548 / 0.62782 | 1.0541 / 0.62740 |
| drift sup | 8.9508e-2 | 8.9457e-2 | 1.2929 | 1.2922 |
| c3,c2,c1,c0 | .4421 .5259 .4235 .4250 | .4419 .5256 .4233 .4250 | 1.6356 2.0253 1.6777 1.7149 | 1.6346 2.0241 1.6766 1.7148 |
| d3,d2,d1,d0 | 6.841e-3 1.5418e-2 2.3648e-2 4.5627e-2 | 6.837e-3 1.541e-2 2.363e-2 4.563e-2 | 1.1992e-3 2.7258e-3 4.2091e-3 8.1602e-3 | 1.1985e-3 2.7243e-3 4.2067e-3 8.1601e-3 |
| mu window | [2.9061e-3, 0.113493] ratio >= 39.05 | [2.9052e-3, 0.113569] | [1.3556e-4, 7.4792e-3] ratio >= 55.17 | [1.3551e-4, 7.4843e-3] |

Outer ramp (R=60 / R=40) and W(R-3), W(R) are in `coupling_constants.json` (P=2: sup|A_zo| 3.15e-18 / 7.90e-8, sup|B_zo| 1.81e-18 / 5.37e-8,
W(60) 2.7637e-20, W(40) 1.2315e-9).

**Caveats.** (1) mu1_min, mu1_strip are the quoted decimals (0.1552/0.3357, 0.1241/0.2211), not recomputed here. (2) Q5 assumes the L1
form (7km) with weight W and cutoff zeta_out; coef_j are those of [L1, zeta_out]. (3) Final float sums over pieces carry a 1e-9 relative
cushion rather than directed rounding; sup values are outward-rounded. (4) Strip pieces use the clamped smoothstep only inside [29,32],
chi'''' jumps at the endpoints are handled by splitting there. (5) Phi uses the upper endpoint of the piece (valid since Phi increases); the
reported "lower witness" for sA is not a rigorous lower bound. (6) The tip ball is hulled (about 1e-9 wide), profile bounds are taken at
sigma_lo floored to 1/16 (monotone, so valid).
