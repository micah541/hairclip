# T2 report: certified initial curve Gamma0 (source `p2_t0.001.npz`, t = 0.001)
Rerun: `python3 -I gamma0.py build` (121 s total, 1 core) and `python3 -I gamma0.py sanity` (1 s). Outputs: `gamma0_coeffs.json`, `results.json`, `build.log`.

**Convention.** alpha = 0 is already an axis crossing (no shift). theta = -pi/2 + phi, phi = alpha + sum_{n<=M} b_n sin(n alpha). This is the opposite orientation to the pi/2 convention in the task: theta(0) = -pi/2, and the half-curve alpha in [0,pi] lies in the LOWER half-plane (y < 0). X(alpha) = X0 + (L0/2pi) int e^{i theta} with X0 = (x0, 0), so Im X = -(L0/2pi) int cos(phi). Closure condition C(b) = int_0^pi cos(phi) = 0. b_n are the exact doubles `-2 Im(th_hat_n)/N`; L0 = 103.65251688569134 (exact double).

**Truncation.** M = 3762 (not several hundred). The criterion sum_{n>M} n^4|b_n| <= 1e-12 gives tail 9.96e-13. The spectrum decays slowly: the tail is 3e-8 at M = 900 and 0.7 at M = 500. Modes above about 1000 are at the 1e-24 float noise level, but they are kept. Cutting at about 900 would cost about 3e-8 in the n^4-weighted tail.
- Discarded even (cosine) part: max |a_n| = 7.3e-11.
- Discarded constant deviation of mean(theta - alpha) from -pi/2: 8.2e-11.

**Closure.** Closing coefficient j = 2, not b_1: dC/db_1 = -0.014 (bad), dC/db_2 = -1.997.
- Validated quadrature (Taylor models of degree 24 per arc with rigorous remainder; 512 leaves; per-arc radius < 1e-22) gives C(m) = 2.0172e-16 +/- 7e-24 at m = -1.002456557732618 (exact double).
- Krawczyk on X = m +/- 2.0202e-16 gives K(X) strictly inside int X, with |1 - C'(X)/c'| <= 3.0e-12. So there is a unique closing b_2 in X, with width 4.04e-16. (c' = -1.99704.)
- Gamma0 is the member of the one-parameter family with b_2 = the root b*.

**Embeddedness.** All quantities use b_2 in the closing interval.
- 701 arcs.
- Every arc has total turning <= 0.04997 (< pi/2, so it is a graph over a direction); the max adjacent-pair turning is 0.0960 (< pi). Each arc and each adjacent pair is therefore injective.
- Non-adjacent arc boxes (X enclosures, outward-rounded to doubles; exact double comparison) are pairwise disjoint after 3 refinement rounds. Min positive L-infinity box separation is 8.9e-3.
- y < 0 on the middle arcs 49..655: certified box upper bound < 0 on each (first one -0.8746).
- End zones, local arguments:
  - Left zone (arcs 0..48): phi enclosure inside (-1.5, 1.5), so cos(phi) > 0 and y strictly decreases from y(0) = 0.
  - Right zone (arcs 656..700): phi enclosure inside (1.6, 4.7), so cos(phi) < 0. With y(pi) = 0 (this is exactly the closure condition C(b*) = 0), y(alpha) = (L0/2pi) int_alpha^pi cos(phi) < 0.
- The enclosure of X(pi) is -0.8227220 +/- 1.7e-16 in Re and +/- 1e-14 in Im, consistent with y(pi) = 0.
- Hence the half-curve is injective and meets the axis only at its two distinct endpoints. With the mirror symmetry, Gamma0 is embedded.
- Slot half-gap (min |y| on the inner-edge plateau, alpha in [1.75, 2.9]): certified lower bound 0.0248220, float value 0.0248223. The certified value is a box bound, so it is lower by about 3e-7.

**Sanity (float).** Ran hls.py (cdt 1e-2, tol 1e-10) from Gamma0 for dt = 0.009 and compared with `p2_t0.01`: max|d theta| = 1.8e-8, dL = 3.3e-8 (L = 103.1028), dX0 = 4e-10, k_max 3.37275273 vs 3.37275269.

**Caveats / not fully rigorous.**
- The code is not independently reviewed. Rigor rests on python-flint/arb balls and on my Taylor-remainder argument. The remainder uses the Lagrange/integral form with box enclosure, an even K and a positive weight s^K.
- The disjointness test rounds arb bounds outward to doubles with `nextafter`.
- The y < 0 argument in the right zone relies on y(pi) = 0 for b*. The boxes themselves only show |Im X(pi)| < 1e-14.
- "Close to the t = 0.001 state" is measured only by the discarded parts above and by the sanity run. It is not certified.
- M is large, so downstream engines will likely want to re-truncate.
