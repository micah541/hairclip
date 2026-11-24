# STEP 4, milestone 1: float go/no-go on the step-3 hand-off (clip 24:0.4:0.05:0.85)

All numbers are FLOAT (double precision), no rigorous code. Everything is in `step4/m1/`. k = k_max = max curvature of the curve (rescaled: lambda*k_max = K0 = 1.5398 holds to 1e-5 in every snapshot). Time t is not used anywhere.

## 0. Verdict

**NO-GO for k_max in [6000, 10000] against the thresholds of lemmaS.pdf (page-1 summary).** The range of k where ALL hold, within what was run (k <= 10132), is **empty**.

What passes and what does not (k = 10132 is the last point; the two ends of the requested range are k = 6081 and 10132):

* **Tip / modulation (H1): passes.** |alpha| <= 5e-4 from k >= 5.6e3 (alpha*k_max = -2.85 at 2e3, -2.77 at 6e3, -2.74 at 1e4). sqrt(Lambda) <= 5e-8 from k ~ 4.7e3; it is 2.3e-8 at k = 6081 (margin 2.2) and 7.7e-9 at k = 10132 (margin 6.5). The deviation w tracks A = ||w~||_rho ~ 214 ... 250 |alpha|^3.
* **Shells on [14,60] (H2): pass** with margin >= 2 from k = 2020 (E0 = 9.9e-8 vs 2e-7 there; 1.4e-10 at 6081).
* **S1 chart (H4): passes.** The curve is a graph over S1 with |slope| <= 1/2 up to y = 0.117 k_max (714 at 6081, 1190 at 10132), above 380 from k >= 3.3e3. (STEP4_ARCH risk R1, extent ~ k^0.75 and ending near 80, does not occur; the "28" and "16" there were the cap XT = 28 of the old chart code.)
* **Lemma O jets: pass.** Cone cylinder (centre r = 200, L = 50, X = 3): slope <= 0.12 from k >= 5.2e3, K2,K3,K4 and ||h'''|| pass with margins 3.6/340/170/220 at k = 6081 (K2 = 0.034, K3 = 3.6e-4, K4 = 7.2e-4, N3 = 6.5e-4). Far cylinders (X = 2, L = 113 lambda, centres r >= 376): slope <= 0.21 from k >= 4.2e3, and the three jets pass by 13, 210, 650 at 6081 (physical kappa at r >= 226 is 2.6 at k = 1e4, as predicted in STEP4_ARCH).
* **Outer windows (H3), levels 3 and 4: FAIL, and only these.** At k = 10132: E3 on [169,253] = 5.1e-6 (1.7x over 3e-6), E4 on [112,169] = 6.2e-6 (2.1x over), E4 on [169,253] = 4.4e-5 (14.8x over). If "all windows [y,1.5y], y in [50,253]" is read literally (this reaches [253,380]): E2 4.3e-6 (1.4x), E3 3.7e-5 (12x), E4 3.2e-4 (106x). At the ARCH candidate k = 6600 (interpolated): E3[169,253] = 6.4e-5 (21x), E4[169,253] = 5.5e-4 (183x), E4[112,169] = 7.7e-5 (26x), E3[112,169] = 8.9e-6 (3.0x), E4[75,112] = 1.0e-5 (3.4x); windows [50,75] pass.

**Why, and how it scales.** In the far field the deviation from the slow manifold is, to 5%, a pure cubic:

    v(y) = f - U = -c3 (y/k_max)^3,   c3 = 38.1 (k = 2e3) ... 39.8 (k = 8e3) ... 40.2 (y = 112, k = 8e3),   for y in [50, 350]

(values -v/(y/k)^3 in `res_*.json`, key `v_samples`; table in section 4). This is the next-order term of the slow manifold (about +1.9 alpha^3 y^3 with alpha*k ~ -2.74); U contains only the first two orders. Since z_j = (y d/dy)^j v = 3^j v, every level is 9x the previous, and with omega = 28/y

    E_j[ya,yb] = 28 * 9^j * c3^2 / (5 k^6) * (yb^5 - ya^5),

i.e. **E_j ~ k^(-5.9)** (fitted exponents 5.87-5.95 over k = 5e3...1e4 for all windows, section 5) and it does not improve faster. The closed form (c3 = 40) and the power-law fits agree to 1-4% on the crossing points. The k at which each remaining item reaches its threshold (extrapolation beyond k = 10132, not run):

* E3 on [169,253]: k >= 1.11e4. E4 on [112,169]: k >= 1.15e4.
* **E4 on [169,253]: k >= 1.6e4** (alpha ~ 1.7e-4). This is the first k at which every item in the task's list holds.
* All sliding windows y <= 253 (through [253,380]): E4 needs k >= 2.2-2.3e4 (E3 1.55e4, E2 1.07e4).
* Tip items and jets then hold with large margin (sqrt(Lambda) ~ 250 alpha^3 ~ 1e-9).

So the true flow meets the list only at k_1 ~ 1.6e4 (2.3e4 for all sliding windows), about 2.4-3.4 times the k ~ 6600 planned. Per ARCH risk R1 this means raising k_1 (cost: +1.3 to +1.7 units of ln k in the float representation; the float run itself needs N = 2^21 uniform at 1.6e4 by the rule N ~ 100 k). A cheaper alternative, which I did not test: the proof side could subtract the cubic term (include alpha^3 h3 with h3' ~ +1.9 y^3 in the far field in U), after which the deviation is O(alpha^4); that is a step-3 change.

## 1. What was run

Solver `step4/hls.py` unchanged (ETDRK4, HLS tangent-angle formulation, symmetric by construction up to the initial closure defect), driver `cont.py` (copy of the run3/run2 logic, resumes from a snapshot, snapshots when k_max >= level, ckpt every 1000 steps). cdt = 1e-2, tail_tol = 1e-10 as in the design runs, except where stated. Start: `design/snap_clip_24_0.4_0.05_0.85_1600.npz` (k = 1613.02, N = 2^17).

| run (log) | settings | snapshots (k, N, tail) |
|---|---|---|
| main (`run_main.log`), then main2 (`run_main2.log`) | cdt 1e-2, tol 1e-10, Nmax 2^19. The machine restarted once after k = 3036; main2 resumed from that snapshot. | 2020.0 (2^18, 5e-17), 3036.3 (2^18, 6e-13) \| 4045.5 (2^18, 6e-11), 5061.0 (2^19, 7e-15), 6080.7 (2^19, 3e-13), 7018.5 (2^19, 3.6e-12), 8090.7 (2^19, 3e-11), 10121.9 (2^19, **4.8e-10 > tol, cap reached**) |
| big (`run_big.log`) | from 7018.5, Nmax 2^20 | 8098.3 (2^19, 3.0e-11), 10131.6 (2^20, 3.6e-15) |
| res20 | from 7018.5, tol 1e-13 (N = 2^20 throughout) | 8098.3 (2^20, 1.4e-17) |
| hr | from 3036.3, tol 1e-13, Nmax 2^20 | 4045.5 (2^19), 5061.0 (2^19), 6080.7 (2^20), 7018.5 (2^20) |
| fine (`run_fine.log`) | from 1613.0, cdt 5e-3, tol 1e-12 | 2018.6 (2^18), 3002.7 (2^18), 4036.1 (2^19) |

Wall: main+main2 ~ 1050 s for 3036 -> 10122 (3000 steps); big 600 s; hr 1100 s; fine 1100 s. Never more than 2 processes. N = 2^19 was NOT enough: at k >= 8000 the tail criterion 1e-10 under-resolves the inner deviation (section 7), so the final 8098 and 10132 snapshots are from N = 2^20 runs. The tables use: 2020 and 3036 from main; 4045.5, 5061.0, 6080.7, 7018.5 from hr (N >= 2^19, tol 1e-13); 8098.3 from res20; 10131.6 from big. The main-run values at 4046-7018 differ from hr by <= 2% in A and <= 7e-4 in any shell energy (section 7).

## 2. Definitions I used (check these against your reading of the lemma)

* **Frame.** For each snapshot, (x_c, lambda) fitted (item 3 below); Gamma = (gamma - x_c)/lambda; y = sigma = x - 2.2490776515. Even part: the lower strand is mirrored about the horizontal line through the tip (`feas/chart.py` convention) and u, u' are averaged at the same S1 abscissa, which removes the odd part (rotation / 3.9e-7 tip offset from the initial closure defect: its rho-norm is 2.7e-6 at k = 2e3 and 1.3e-5 at 1e4, i.e. >> A, so the symmetrization is essential for Lambda; in the far field the odd slope is <= 1.2e-8).
* **Slope f, deviation v, U.** f = u' of the normal graph over S1. For sigma < 22 the foot point of the normal and u' = (1 - k1 u) tan(phi), phi = angle(curve tangent, S1 tangent) from theta Hermite-interpolated at the foot; for sigma >= 22 S1 is a straight ray (p < 1e-13) and f = tan(angle to the ray) as a function of y = tangential coordinate (cubic spline). U = alpha psi3~' + alpha^2 h2~' from `agent_int/v3/c4fun.py` (C^4 versions; arb midpoints). Beyond x = 28 c4fun's own explicit far-field extension (power law 1.3057 y^1.0246 plus the exponentially small term) is used, which is what step 3 uses; no extra extension was needed. v = f - U.
* **Modulation.** (x_c, lambda, alpha) by least squares on the even u over x <= 28, rho-weighted half-line L2, residual u - alpha psi3~ - alpha^2 h2~. This reproduces the constraints <w,f1> = <w,f2> = 0 (translation and scaling modes; cutoff differences are ~ rho(28) = 6e-10) but replaces <w~,f3> = 0 (f3 = chi F - a eta1 - b eta3, the approximate ADJOINT eigenfunction, which is only available as arb hulls) by <w, psi3~ + 2 alpha h2~> = 0. See caveat C1.
* **Lambda = E + mu E_I** exactly as in ASSEMBLY.md section 0: half-line, trapezoid on the 0.01 grid; E = ||chi32 w||_rho^2 + 1e-5 ||(chi32 w)''||_rho^2, E_I = ||zeta w'||_W^2 + 3e-5 ||(zeta w')''||_W^2, mu = 0.08, zeta rising on [10.25,16.25] and falling on [57,60] (smoothstep S), W = (1+x^2)^-2 for x <= 30.25 (blend beyond). w', w'' by 4th-order differences of the exact slope w' (w itself only enters A).
* **Shell energies.** E_j[ya,yb] = int_{ln ya}^{ln yb} omega z_j^2 ds, z_j = d^j v/ds^j (s = ln y), omega = 28/y, sharp windows (no cutoff theta, which can only lower them), exact end points. Levels j = 0..4. z_j from local degree-12 polynomial least squares on a uniform s-grid (spacing 0.002, half-width 0.4); five other stencils agree to <= 1e-5 relative, and y^3 is reproduced to 8e-6 (end-point rule). E4 is applied to the same windows as E1-E3 (the summary says "E4 <= 3e-6 at tau1"; I read it as the same windows).
* **Windows.** The four requested ([50,75], [75,112], [112,169], [169,253]) and the sliding family [y, 1.5y], y in [50,253] (maximum over y; it always sits at y = 253, i.e. the window [253,380]).
* **Chart extent.** First y with |f| > 0.5 along the strand (the q-range of the graph, i.e. where tangential coordinate stops increasing, is ~2.1x larger and is in the json under `chart_ext`).
* **Lemma O jets** (all from theta via spectral derivatives in alpha; in rescaled units the highest wavenumber is ~13 per unit length, so noise is ~1e-13). For a centre at arc position s_c: psi = phi - phi_c; h_x = tan psi, h_xx = kappa sec^3, h_xxx = sec^4 (kappa_s + 3 kappa^2 tan psi), h_xxxx = sec [4 sec^4 tan(psi) kappa (kappa_s + 3 kappa^2 tan psi) + sec^4 (kappa_ss + 6 kappa kappa_s tan psi + 3 kappa^3 sec^2)], sup over |x| <= X L. Cone: centre r = |Gamma| = 200 (and sup over centres 200-226 and 200-400), L = r_c/4, X = 3, cylinder units h^(j) -> L^(j-1) h^(j), N3 = ||h'''||_{L2(-3,3)} = L^(3/2) ||h_xxx||; thresholds 0.12 for slope and j = 2,3,4 and 0.14 for N3. Far: X = 2, L = 113 (lambda units), centres r_c >= 376 (so that the window stays at r >= 150) every ~5-50 lambda along the whole strand to the far end, thresholds 0.21 (slope), 8.9e-3, 7.0e-5, 6.2e-7, N3 = L^(3/2)||h_xxx||_{L2(|x|<=2L)} <= 1.0; and the pointwise kappa, kappa_s, kappa_ss for r >= 226 on the whole closed curve (the reading in STEP4_ARCH H5). Not checked: Lemma O H1 (sup|h0|), H4 (fifth derivative), H8 (other arms), H9; the time ratios are by construction.

## 3. Modulation fit and Lambda

Full-precision lines are in `res_<tag>.json`. N = base grid size (the snapshot's).

TABLE_MOD

x_c is the real part in the lab frame; it converges to 19.19894876348 (+ 3.925e-7 i), so the singular point is determined to 1e-11. A/|alpha|^3 is ~ 210-250 up to k = 8e3, i.e. w = O(alpha^3) as step 3 assumes (the earlier agent_handoff2 number "r2/|alpha|^3 = 32-46" is for a free h2 coefficient; with beta2 = alpha^2 imposed it is ~ 5x larger, which is what U = alpha psi3' + alpha^2 h2' implies). At k = 10132, A = 7.7e-9 is above 250 |alpha|^3 = 4.9e-9: an unexplained excess of ~6e-9 (section 7). The 2^19, tol 1e-10 snapshot at 10122 gave A = 5.6e-7 (under-resolved, discarded).

## 4. Every quantity against its threshold

Format: value, then `OK xM` (M = threshold/value) or `FAIL /M` (M = value/threshold). Columns are k_max. H1 = hand-off item (1), H2 = (2), H3 = windows, H3s = sliding windows, H4 = chart, H5c = cone cylinders, H5f = far cylinders and pointwise jets.

TABLE_MASTER

Cubic law check (-v/(y/k)^3, `v_samples` in the json): k = 2020: 39.5 (y = 50), 38.6 (112), 38.0 (169), 38.1 (253), 40.2 (350); k = 4046: 40.0, 39.7, 39.3, 38.6, 38.0; k = 8098: 40.2, 40.3, 40.2, 39.8, 39.2; k = 10132: 40.2 (y = 253: 6.24e-4 / (253/10132)^3 = 40.1).

## 5. Where it holds, scaling, extrapolation

First k (from the fits/closed form) at which each group holds by itself: |alpha| 5.6e3; sqrt(Lambda) 4.7e3; [14,60] 2.0e3; chart 3.3e3; cone slope 5.2e3; far slope 4.2e3; windows E0,E1 (to 253) 5.8e3; E2 7.7e3; **E3 1.11e4; E4 1.6e4**; sliding windows E2 1.07e4, E3 1.55e4, **E4 2.2e4**.

Exponents and extrapolated crossings for the windows that can fail (E ~ k^-p, linear fit of log E on log k over the five points k = 5061, 6081, 7018, 8098, 10132; closed form E_j = 28*9^j c3^2 (yb^5 - ya^5)/(5 k^6) with c3 = 40):

TABLE_SCALE

Largest window start y* such that all windows [y,1.5y] with y <= y* satisfy the level (first violation, 0.01 grid in ln y; `scan.py`): E0, E1, E2, E3, E4 at k = 5061: 142, 159, 102, 66, 50; k = 6081: 178, 198, 128, 82, 53; k = 7018: 213, 237, 151, 97, 63; k = 8098: 252, >253, 181, 115, 74; k = 10132: >253, >253, 237, 151, 97. (E4 would pass all windows up to y ~ 97, E3 up to ~151, at k = 10132; at k = 6081 up to 53 and 82.) So a shorter strip (y <= ~100) or no level 4 would be reachable around k = 8-10e3, a full strip to 253 is not.

## 6. Notes on items that pass

* Cone cylinders: the slope of the cone cylinder is 0.075 + (r_c/200 - 1)*... grows ~ linearly in the centre; sup over centres 200-226: slope 0.112 at k = 6081, 0.067 at 10132; over 200-400: 0.196 and 0.118 (so if the cone family has to reach centre 400 the slope threshold 0.12 needs k ~ 1e4). K2 at r = 200 is 0.034 at 6081, in line with the Lemma S derived 0.044 at alpha = 5e-4; K3 and K4 are 40x and 60x below the derived 0.014, 0.041, because the derived values are bounds from E3, E4 = 3e-6 whereas the true jets are tiny (the E3, E4 thresholds, not the jets, are what fails).
* Far cylinders: the maximum sits at the smallest centre (r = 376) and falls like 1/k; physical kappa at r = 226 is kappa_r/lambda = 3.96e-4/1.52e-4 = 2.6 at k = 1e4, kappa_s and kappa_ss are 1300x and 3400x below the thresholds.
* Consistency of the two chart methods in the overlap y in [14,22] (normal-graph slope vs ray slope): 1.8e-8 (the ray ignores p(sigma) ~ 1e-7...1e-9); the normal-graph value is used there.

## 7. Numerical floors (honest list)

1. **Interpolation of the strand (inner part).** The chart code I inherited (`feas/chart.py:graph`, `agent_handoff2/handoff2.py:hermite_graph`) has an inert Newton step (the derivative `evd` is divided by the cell width, so dt is smaller by that factor): the foot parameter stays at the linear-secant guess and u has an h^2 error (1e-8, with a period-3 saw-tooth in u' of 1e-6 near the tip). I fixed this in `prof.py:graph2` and take the slope from theta (cubic Hermite with exact kappa) instead of the tangent of the position Hermite cubic. With the fix the error scales as h^4: at k = 6081, upsampling 8 -> 16 changes A by 7.6% (1.7e-9); I use upsampling 16 (ds/lambda ~ 0.02), residual estimated <= 1.5e-10 (0.7% of A). Numbers from the older agents' tables therefore carry a ~1e-8 floor in u.
2. **Spatial/temporal resolution of the solve.** (a) cdt/2, tol/100 (fine) vs main at the same k: after correcting for the tiny k offsets (dk/k = -7e-4, -1.1e-2, -2.3e-3), far-field and [14,60] shell energies agree to <= 1e-4 relative and agree with the k^-5.9 law; A agrees to 0.03% (2000), 0.06% (3000), 1.5% (4000). (b) hr (tol 1e-13, N 2^19-2^20) vs main (tol 1e-10, N 2^18-2^19): A changes by -1.3% (k 4046), -1.3% (5061), -1.5% (6081), -1.9% (7018), i.e. 1e-9, 5e-10, 3.4e-10, 2.9e-10; shell energies by <= 7e-4 ([14,60] E0) and <= 3e-5 in [169,253]. (c) **N = 2^19 with tail 3e-11 at k = 8098 gives A = 1.63e-8 versus 9.81e-9 at N = 2^20** (a +6.5e-9 error, +66%); the 2^19 snapshot at 10122 (tail 4.8e-10) gives A = 5.6e-7. So the tail test at 1e-10 is not a sufficient criterion for the inner deviation once A < 3e-8; use <= 1e-12 (the hr/res20/big snapshots do). Far-field energies at 8098: 2^19 vs 2^20 differ by <= 3e-5 for [169,253] and 5.5e-4 (E0) to 3.9e-3 (E4) on [14,60].
3. **Unexplained excess of A at the largest k.** A(10132) = 7.7e-9 vs 250|alpha|^3 = 4.9e-9, A/|alpha|^3 jumps from 250 (8098) to 390. N is converged (tail 4e-15), cdt was 1e-2 only (fine-vs-main at k <= 4000 says time-step effects are <= 1.5%, but untested at 1e4). Candidates: the S1 profile / theta1 spline (handoff2 quotes ~1e-7 in theta1, which would show as A ~ 1e-7, so it is at most 7.7e-9 here), position round-off (1e-14 absolute / lambda 1.5e-4 = 7e-11), or a genuine O(alpha^4) term. Taken as a floor of ~6e-9 on A; the margin to 5e-8 is still 6.5.
4. **Derivative levels.** Stencil gains (200 points, degree 12): 1.5e3 for j = 3, 5.8e4 for j = 4. With data noise <= 1e-12 (theta round-off and spline error) the z_4 noise is <= 6e-8 against the level that matters (|z_4| ~ 8e-3 for E4 = 3e-6); the 5-stencil test gives 1e-5 relative on the energies.
5. **Double precision / t.** Tracking k, not t, throughout. Positions come from an FFT antiderivative of exp(i theta): absolute error ~1e-14 on coordinates of size ~ 20 -> 7e-11 in rescaled units, below 1% of A for k <= 8e3 and ~1% at 1e4. The far-field slope uses only theta (round-off 1e-16).
6. **Restarts.** Every resume restarts the step-size ramp (~150 small steps); the main2 and big runs, both started from different snapshots, agree at 8091 / 8098 up to the k offset.
7. **Not independent of fit choices:** alpha is determined to ~ 1e-6 relative (up8 vs up16: 8e-7); an alpha error dalpha shifts v at y = 253 by dalpha*psi3'(253) = 380 dalpha ~ 4e-10 -> negligible against v ~ 6e-4.

## 8. Caveats

* **C1 (modulation).** Lambda is evaluated at the least-squares modulation, not at <w~,f3>_rho = 0 with the certified adjoint function f3 (which I could not evaluate: only arb hull rows exist). The LS fit makes w orthogonal to psi3~ + 2 alpha h2~, so A_f3 <= A / cos(gamma) where gamma is the angle between f3 and psi3~ in L2(rho). The margin on sqrt(Lambda) is 1.3 at k = 5061, 2.2 at 6081, 3.4 at 7018, 6.5 at 10132, so the verdict "tip passes" survives gamma < 38 deg at 5061 and < 63 deg at 6081. A computation of the adjoint eigenfunction would settle this; it does not change the NO-GO, which comes from E3/E4.
* **C2 (E4 convention).** If "E4 <= 3e-6 at tau1" is meant on [50,60] only, or on the strip to 150 only, the outcome changes: E4 on [75,112] passes at k >= 8.6e3 and E4 on [112,169] needs 1.15e4. Please confirm which windows.
* **C3 (extrapolation).** k = 1.6e4 / 2.3e4 are extrapolations of a clean law (c3 = 38-40, p = 5.9, 8 snapshots), not runs. The law is the leading far-field behaviour; I did not test it beyond k = 10132 (c3 drifts upward by ~5% over k = 2e3 -> 1e4, which would lower the k_cross by ~1%).
* Cone/far family boundaries are my choice (centres 200-226 for the cone, >= 376 for far); the sup over cone centres up to 400 is also in the json (`cone_sup`).
* One odd detail: 1/3 of the budget went to the restart and the N = 2^20 reruns. The main-run tables at 4046-7018 and 8091/10122 (N = 2^19) are in `res_main_*.json` for comparison; they are not used.

## 9. Files (all in `step4/m1/`)

`cont.py` (resumable run), `geo.py` (S1 on x <= 64 + C4 psi3~/h2~ cache `c4_cache.npz`), `prof.py` (symmetrized fit, graph2, ray profile, stage 1), `en.py` (Lambda, v, shell energies), `jets.py` (Lemma O), `drive.py SNAP UP TAG` -> `p_TAG.npz`, `res_TAG.json`; `make_report.py`, `tabs.py`, `scaling_md.py`, `scan.py` (tables); snapshots `snap_*.npz` (main, main2, big, res20, hr, fine), logs `run_*.log`. Old stage-1 results from before the Newton fix were discarded.
