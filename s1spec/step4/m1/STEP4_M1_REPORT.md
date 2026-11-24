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

    v(y) = f - U = -c3 (y/k_max)^3,   c3 = 38-40 (k = 2e3) ... 40-40.6 (k = 8e3-1e4),   for y in [50, 350]

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

| tag | k_max | N | lambda | lambda*k_max | x_c (lab, real) | alpha | alpha*k_max | alpha/5e-4 | A=\|w~\|_rho | A/\|alpha\|^3 | sqrt(mu E_I) | sqrt(Lambda) | sup \|w'\| on [16.25,57] | S1-chart: first y with \|f\|>0.5 | odd-part norm (rho, x<=28) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| main_1600 | 1613.0 | 131072 | 9.5454e-04 | 1.53970 | 19.19894876294 | -1.7748e-03 | -2.8627 | 3.55 | 1.176e-06 | 210 | 2.24e-07 | 1.20e-06 | 1.53e-03 | 189 | 2.1e-06 |
| main_2000 | 2020.0 | 262144 | 7.6223e-04 | 1.53973 | 19.19894876326 | -1.4093e-03 | -2.8468 | 2.82 | 5.996e-07 | 214 | 1.14e-07 | 6.10e-07 | 7.87e-04 | 237 | 2.7e-06 |
| main_3000 | 3036.3 | 262144 | 5.0712e-04 | 1.53976 | 19.19894876344 | -9.2822e-04 | -2.8183 | 1.86 | 1.768e-07 | 221 | 3.38e-08 | 1.80e-07 | 2.34e-04 | 357 | 4.0e-06 |
| hr_4000 | 4045.5 | 524288 | 3.8061e-04 | 1.53977 | 19.19894876347 | -6.9175e-04 | -2.7985 | 1.38 | 7.496e-08 | 226 | 1.43e-08 | 7.63e-08 | 9.93e-05 | 475 | 5.3e-06 |
| hr_5000 | 5061.0 | 524288 | 3.0424e-04 | 1.53978 | 19.19894876348 | -5.4991e-04 | -2.7831 | 1.10 | 3.836e-08 | 231 | 7.31e-09 | 3.90e-08 | 5.08e-05 | 594 | 6.6e-06 |
| hr_6000 | 6080.7 | 1048576 | 2.5323e-04 | 1.53979 | 19.19894876348 | -4.5563e-04 | -2.7706 | 0.91 | 2.221e-08 | 235 | 4.22e-09 | 2.26e-08 | 2.93e-05 | 714 | 8.0e-06 |
| hr_7000 | 7018.5 | 1048576 | 2.1939e-04 | 1.53979 | 19.19894876348 | -3.9336e-04 | -2.7608 | 0.79 | 1.449e-08 | 238 | 2.75e-09 | 1.47e-08 | 1.91e-05 | 824 | 9.2e-06 |
| res20_8000 | 8098.3 | 1048576 | 1.9014e-04 | 1.53980 | 19.19894876348 | -3.3971e-04 | -2.7511 | 0.68 | 9.808e-09 | 250 | 1.80e-09 | 9.97e-09 | 1.25e-05 | 951 | 1.1e-05 |
| big_10000 | 10131.6 | 1048576 | 1.5198e-04 | 1.53980 | 19.19894876348 | -2.7004e-04 | -2.7360 | 0.54 | 7.681e-09 | 390 | 9.52e-10 | 7.74e-09 | 6.41e-06 | 1190 | 1.3e-05 |


x_c is the real part in the lab frame; it converges to 19.19894876348 (+ 3.925e-7 i), so the singular point is determined to 1e-11. A/|alpha|^3 is ~ 210-250 up to k = 8e3, i.e. w = O(alpha^3) as step 3 assumes (the earlier agent_handoff2 number "r2/|alpha|^3 = 32-46" is for a free h2 coefficient; with beta2 = alpha^2 imposed it is ~ 5x larger, which is what U = alpha psi3' + alpha^2 h2' implies). At k = 10132, A = 7.7e-9 is above 250 |alpha|^3 = 4.9e-9: an unexplained excess of ~6e-9 (section 7). The 2^19, tol 1e-10 snapshot at 10122 gave A = 5.6e-7 (under-resolved, discarded).

## 4. Every quantity against its threshold

Format: value, then `OK xM` (M = threshold/value) or `FAIL /M` (M = value/threshold). Columns are k_max. H1 = hand-off item (1), H2 = (2), H3 = windows, H3s = sliding windows, H4 = chart, H5c = cone cylinders, H5f = far cylinders and pointwise jets.

| quantity | k=2020 | k=3036 | k=4046 | k=5061 | k=6081 | k=7018 | k=8098 | k=10132 |
|---|---|---|---|---|---|---|---|---|
| |alpha|  (thr 5e-4) | 1.41e-03 FAIL /2.82 | 9.28e-04 FAIL /1.86 | 6.92e-04 FAIL /1.38 | 5.50e-04 FAIL /1.10 | 4.56e-04 OK x1.1 | 3.93e-04 OK x1.3 | 3.40e-04 OK x1.5 | 2.70e-04 OK x1.9 |
| sqrt(Lambda)  (thr 5e-8) | 6.10e-07 FAIL /12.21 | 1.80e-07 FAIL /3.60 | 7.63e-08 FAIL /1.53 | 3.90e-08 OK x1.3 | 2.26e-08 OK x2.2 | 1.47e-08 OK x3.4 | 9.97e-09 OK x5.0 | 7.74e-09 OK x6.5 |
|   of which sqrt(E) = (A^2+eps B^2)^.5 | 6.00e-07 | 1.77e-07 | 7.50e-08 | 3.84e-08 | 2.22e-08 | 1.45e-08 | 9.81e-09 | 7.68e-09 |
|   of which sqrt(mu E_I) | 1.14e-07 | 3.38e-08 | 1.43e-08 | 7.31e-09 | 4.22e-09 | 2.75e-09 | 1.80e-09 | 9.52e-10 |
| E0 on [14,60]  (thr 2e-07) | 9.96e-08 OK x2.0 | 8.78e-09 OK x22.8 | 1.58e-09 OK x126.4 | 4.14e-10 OK x482.6 | 1.38e-10 OK x1447.7 | 5.85e-11 OK x3416.9 | 2.49e-11 OK x8032.9 | 6.62e-12 OK x30197.4 |
| E1 on [14,60]  (thr 0.0001) | 8.99e-07 OK x111.2 | 7.99e-08 OK x1252.3 | 1.44e-08 OK x6928.6 | 3.79e-09 OK x26400.4 | 1.26e-09 OK x79118.5 | 5.36e-10 OK x186655.7 | 2.28e-10 OK x439358.1 | 5.99e-11 OK x1670803.9 |
| E2 on [14,60]  (thr 0.2) | 7.95e-06 OK x25160.2 | 7.13e-07 OK x280695.8 | 1.29e-07 OK x1545222.7 | 3.41e-08 OK x5870354.0 | 1.14e-08 OK x17560796.7 | 4.83e-09 OK x41385409.7 | 2.05e-09 OK x97386851.6 | 5.38e-10 OK x371538376.7 |
| E0 on [50, 75]  (thr 2e-07) | 2.64e-07 FAIL /1.32 | 2.34e-08 OK x8.5 | 4.23e-09 OK x47.3 | 1.11e-09 OK x180.3 | 3.70e-10 OK x540.2 | 1.57e-10 OK x1274.2 | 6.67e-11 OK x2997.7 | 1.76e-11 OK x11360.2 |
| E1 on [50, 75]  (thr 3e-06) | 2.34e-06 OK x1.3 | 2.10e-07 OK x14.3 | 3.81e-08 OK x78.8 | 1.00e-08 OK x299.4 | 3.35e-09 OK x895.8 | 1.42e-09 OK x2111.0 | 6.04e-10 OK x4966.6 | 1.59e-10 OK x18924.1 |
| E2 on [50, 75]  (thr 3e-06) | 2.05e-05 FAIL /6.84 | 1.85e-06 OK x1.6 | 3.39e-07 OK x8.9 | 8.95e-08 OK x33.5 | 3.00e-08 OK x100.0 | 1.27e-08 OK x235.3 | 5.42e-09 OK x553.0 | 1.42e-09 OK x2108.4 |
| E3 on [50, 75]  (thr 3e-06) | 1.78e-04 FAIL /59.27 | 1.62e-05 FAIL /5.41 | 2.99e-06 OK x1.0 | 7.93e-07 OK x3.8 | 2.67e-07 OK x11.2 | 1.14e-07 OK x26.4 | 4.85e-08 OK x61.9 | 1.27e-08 OK x235.5 |
| E4 on [50, 75]  (thr 3e-06) | 1.53e-03 FAIL /511.28 | 1.41e-04 FAIL /46.84 | 2.61e-05 FAIL /8.69 | 6.97e-06 FAIL /2.32 | 2.36e-06 OK x1.3 | 1.01e-06 OK x3.0 | 4.31e-07 OK x7.0 | 1.14e-07 OK x26.4 |
| E0 on [75, 112]  (thr 2e-07) | 1.90e-06 FAIL /9.52 | 1.71e-07 OK x1.2 | 3.11e-08 OK x6.4 | 8.19e-09 OK x24.4 | 2.74e-09 OK x73.0 | 1.16e-09 OK x172.0 | 4.95e-10 OK x404.3 | 1.30e-10 OK x1537.5 |
| E1 on [75, 112]  (thr 3e-06) | 1.67e-05 FAIL /5.56 | 1.51e-06 OK x2.0 | 2.76e-07 OK x10.9 | 7.31e-08 OK x41.1 | 2.45e-08 OK x122.4 | 1.04e-08 OK x287.6 | 4.44e-09 OK x675.4 | 1.17e-09 OK x2571.0 |
| E2 on [75, 112]  (thr 3e-06) | 1.45e-04 FAIL /48.37 | 1.32e-05 FAIL /4.39 | 2.43e-06 OK x1.2 | 6.46e-07 OK x4.6 | 2.18e-07 OK x13.8 | 9.28e-08 OK x32.3 | 3.96e-08 OK x75.7 | 1.04e-08 OK x287.6 |
| E3 on [75, 112]  (thr 3e-06) | 1.26e-03 FAIL /420.81 | 1.14e-04 FAIL /38.08 | 2.12e-05 FAIL /7.06 | 5.67e-06 FAIL /1.89 | 1.92e-06 OK x1.6 | 8.21e-07 OK x3.7 | 3.51e-07 OK x8.5 | 9.29e-08 OK x32.3 |
| E4 on [75, 112]  (thr 3e-06) | 1.11e-02 FAIL /3714.17 | 9.85e-04 FAIL /328.37 | 1.83e-04 FAIL /60.99 | 4.92e-05 FAIL /16.41 | 1.67e-05 FAIL /5.58 | 7.20e-06 FAIL /2.40 | 3.09e-06 FAIL /1.03 | 8.22e-07 OK x3.6 |
| E0 on [112, 169]  (thr 2e-07) | 1.45e-05 FAIL /72.47 | 1.31e-06 FAIL /6.53 | 2.39e-07 FAIL /1.20 | 6.34e-08 OK x3.2 | 2.13e-08 OK x9.4 | 9.08e-09 OK x22.0 | 3.87e-09 OK x51.7 | 1.02e-09 OK x196.4 |
| E1 on [112, 169]  (thr 3e-06) | 1.27e-04 FAIL /42.31 | 1.14e-05 FAIL /3.81 | 2.10e-06 OK x1.4 | 5.60e-07 OK x5.4 | 1.89e-07 OK x15.9 | 8.07e-08 OK x37.2 | 3.45e-08 OK x87.0 | 9.09e-09 OK x330.0 |
| E2 on [112, 169]  (thr 3e-06) | 1.12e-03 FAIL /374.11 | 9.92e-05 FAIL /33.08 | 1.83e-05 FAIL /6.11 | 4.90e-06 FAIL /1.63 | 1.66e-06 OK x1.8 | 7.12e-07 OK x4.2 | 3.05e-07 OK x9.8 | 8.08e-08 OK x37.1 |
| E3 on [112, 169]  (thr 3e-06) | 1.03e-02 FAIL /3417.59 | 8.63e-04 FAIL /287.72 | 1.59e-04 FAIL /52.86 | 4.26e-05 FAIL /14.18 | 1.45e-05 FAIL /4.83 | 6.23e-06 FAIL /2.08 | 2.68e-06 OK x1.1 | 7.13e-07 OK x4.2 |
| E4 on [112, 169]  (thr 3e-06) | 1.01e-01 FAIL /33651.82 | 7.62e-03 FAIL /2540.09 | 1.37e-03 FAIL /456.93 | 3.67e-04 FAIL /122.28 | 1.25e-04 FAIL /41.69 | 5.40e-05 FAIL /17.99 | 2.33e-05 FAIL /7.77 | 6.24e-06 FAIL /2.08 |
| E0 on [169, 253]  (thr 2e-07) | 1.07e-04 FAIL /532.72 | 9.41e-06 FAIL /47.07 | 1.73e-06 FAIL /8.65 | 4.61e-07 FAIL /2.30 | 1.55e-07 OK x1.3 | 6.65e-08 OK x3.0 | 2.84e-08 OK x7.0 | 7.51e-09 OK x26.6 |
| E1 on [169, 253]  (thr 3e-06) | 9.69e-04 FAIL /322.94 | 8.23e-05 FAIL /27.44 | 1.51e-05 FAIL /5.03 | 4.03e-06 FAIL /1.34 | 1.36e-06 OK x2.2 | 5.85e-07 OK x5.1 | 2.51e-07 OK x12.0 | 6.66e-08 OK x45.1 |
| E2 on [169, 253]  (thr 3e-06) | 9.38e-03 FAIL /3125.35 | 7.27e-04 FAIL /242.48 | 1.31e-04 FAIL /43.72 | 3.50e-05 FAIL /11.66 | 1.19e-05 FAIL /3.96 | 5.11e-06 FAIL /1.70 | 2.20e-06 OK x1.4 | 5.86e-07 OK x5.1 |
| E3 on [169, 253]  (thr 3e-06) | 1.03e-01 FAIL /34412.17 | 6.64e-03 FAIL /2213.35 | 1.15e-03 FAIL /382.97 | 3.03e-04 FAIL /101.09 | 1.03e-04 FAIL /34.28 | 4.43e-05 FAIL /14.76 | 1.91e-05 FAIL /6.37 | 5.12e-06 FAIL /1.71 |
| E4 on [169, 253]  (thr 3e-06) | 1.43e+00 FAIL /477041.89 | 6.53e-02 FAIL /21763.75 | 1.03e-02 FAIL /3444.40 | 2.65e-03 FAIL /882.56 | 8.89e-04 FAIL /296.26 | 3.82e-04 FAIL /127.24 | 1.65e-04 FAIL /54.93 | 4.43e-05 FAIL /14.78 |
| E0 sliding max, windows [y,1.5y], y in [50,253]  (thr 2e-07) | 8.85e-04 FAIL /4426.19 | 6.98e-05 FAIL /349.17 | 1.26e-05 FAIL /63.00 | 3.35e-06 FAIL /16.76 | 1.13e-06 FAIL /5.67 | 4.86e-07 FAIL /2.43 | 2.09e-07 FAIL /1.04 | 5.54e-08 OK x3.6 |
| E1 sliding max, windows [y,1.5y], y in [50,253]  (thr 3e-06) | 9.59e-03 FAIL /3195.49 | 6.34e-04 FAIL /211.41 | 1.11e-04 FAIL /36.91 | 2.92e-05 FAIL /9.74 | 9.88e-06 FAIL /3.29 | 4.24e-06 FAIL /1.41 | 1.83e-06 OK x1.6 | 4.87e-07 OK x6.2 |
| E2 sliding max, windows [y,1.5y], y in [50,253]  (thr 3e-06) | 1.32e-01 FAIL /44163.39 | 6.13e-03 FAIL /2042.58 | 9.92e-04 FAIL /330.80 | 2.56e-04 FAIL /85.33 | 8.60e-05 FAIL /28.65 | 3.69e-05 FAIL /12.28 | 1.59e-05 FAIL /5.29 | 4.25e-06 FAIL /1.42 |
| E3 sliding max, windows [y,1.5y], y in [50,253]  (thr 3e-06) | 2.89e+00 FAIL /963770.42 | 6.73e-02 FAIL /22424.01 | 9.36e-03 FAIL /3118.80 | 2.29e-03 FAIL /763.03 | 7.53e-04 FAIL /250.89 | 3.20e-04 FAIL /106.70 | 1.37e-04 FAIL /45.79 | 3.68e-05 FAIL /12.26 |
| E4 sliding max, windows [y,1.5y], y in [50,253]  (thr 3e-06) | 1.16e+02 FAIL /38564676.98 | 9.28e-01 FAIL /309321.96 | 9.81e-02 FAIL /32684.12 | 2.16e-02 FAIL /7190.77 | 6.76e-03 FAIL /2254.87 | 2.82e-03 FAIL /938.62 | 1.19e-03 FAIL /397.86 | 3.17e-04 FAIL /105.75 |
| chart: first y with |f|>0.5  (need >= 380) | 237 FAIL | 357 FAIL | 475 OK | 594 OK | 714 OK | 824 OK | 951 OK | 1190 OK |
| cone r_c=200, slope  (thr 0.12) | 2.95e-01 FAIL /2.46 | 1.98e-01 FAIL /1.65 | 1.49e-01 FAIL /1.24 | 1.19e-01 OK x1.0 | 9.93e-02 OK x1.2 | 8.61e-02 OK x1.4 | 7.46e-02 OK x1.6 | 5.96e-02 OK x2.0 |
| cone r_c=200, K2  (thr 0.12) | 1.12e-01 OK x1.1 | 7.03e-02 OK x1.7 | 5.15e-02 OK x2.3 | 4.06e-02 OK x3.0 | 3.35e-02 OK x3.6 | 2.89e-02 OK x4.1 | 2.50e-02 OK x4.8 | 1.99e-02 OK x6.0 |
| cone r_c=200, K3  (thr 0.12) | 1.29e-02 OK x9.3 | 3.35e-03 OK x35.8 | 1.32e-03 OK x91.0 | 6.44e-04 OK x186.4 | 3.56e-04 OK x337.5 | 2.24e-04 OK x536.3 | 1.81e-04 OK x663.7 | 2.30e-04 OK x522.5 |
| cone r_c=200, K4  (thr 0.12) | 3.97e-03 OK x30.2 | 7.69e-04 OK x156.0 | 8.57e-04 OK x140.1 | 8.08e-04 OK x148.4 | 7.24e-04 OK x165.8 | 6.45e-04 OK x186.0 | 5.73e-04 OK x209.5 | 4.64e-04 OK x258.6 |
| cone r_c=200, N3  (thr 0.14) | 1.44e-02 OK x9.7 | 4.73e-03 OK x29.6 | 2.13e-03 OK x65.6 | 1.12e-03 OK x124.5 | 6.47e-04 OK x216.5 | 4.16e-04 OK x336.4 | 2.70e-04 OK x518.3 | 1.65e-04 OK x849.2 |
| cone r_c in [200,226], slope  (thr 0.12) | 3.32e-01 FAIL /2.77 | 2.22e-01 FAIL /1.85 | 1.67e-01 FAIL /1.39 | 1.34e-01 FAIL /1.12 | 1.12e-01 OK x1.1 | 9.69e-02 OK x1.2 | 8.41e-02 OK x1.4 | 6.72e-02 OK x1.8 |
| cone r_c in [200,226], K2  (thr 0.12) | 1.30e-01 FAIL /1.08 | 8.03e-02 OK x1.5 | 5.85e-02 OK x2.1 | 4.60e-02 OK x2.6 | 3.80e-02 OK x3.2 | 3.27e-02 OK x3.7 | 2.82e-02 OK x4.2 | 2.25e-02 OK x5.3 |
| cone r_c in [200,226], K3  (thr 0.12) | 1.90e-02 OK x6.3 | 4.97e-03 OK x24.1 | 1.94e-03 OK x61.9 | 9.42e-04 OK x127.4 | 5.23e-04 OK x229.4 | 3.29e-04 OK x364.8 | 2.07e-04 OK x579.6 | 2.30e-04 OK x521.3 |
| cone r_c in [200,226], K4  (thr 0.12) | 6.94e-03 OK x17.3 | 1.20e-03 OK x99.9 | 8.57e-04 OK x140.1 | 8.57e-04 OK x140.0 | 7.84e-04 OK x153.0 | 7.10e-04 OK x169.1 | 6.32e-04 OK x189.8 | 5.19e-04 OK x231.2 |
| cone r_c in [200,226], N3  (thr 0.14) | 1.99e-02 OK x7.0 | 6.55e-03 OK x21.4 | 2.98e-03 OK x47.0 | 1.59e-03 OK x88.3 | 9.29e-04 OK x150.8 | 6.01e-04 OK x233.1 | 3.87e-04 OK x362.1 | 2.06e-04 OK x680.6 |
| far cyl. (X=2,L=113, r_c>=376), slope  (thr 0.21) | 4.17e-01 FAIL /1.98 | 2.85e-01 FAIL /1.36 | 2.17e-01 FAIL /1.03 | 1.75e-01 OK x1.2 | 1.47e-01 OK x1.4 | 1.28e-01 OK x1.6 | 1.11e-01 OK x1.9 | 8.93e-02 OK x2.4 |
| far cyl. (X=2,L=113, r_c>=376), h2  (thr 0.0089) | 2.36e-03 OK x3.8 | 1.44e-03 OK x6.2 | 1.05e-03 OK x8.5 | 8.23e-04 OK x10.8 | 6.79e-04 OK x13.1 | 5.86e-04 OK x15.2 | 5.05e-04 OK x17.6 | 4.02e-04 OK x22.1 |
| far cyl. (X=2,L=113, r_c>=376), h3  (thr 7e-05) | 7.77e-06 OK x9.0 | 2.30e-06 OK x30.5 | 1.00e-06 OK x70.0 | 5.36e-07 OK x130.6 | 3.28e-07 OK x213.2 | 2.26e-07 OK x309.7 | 1.57e-07 OK x445.6 | 9.02e-08 OK x776.3 |
| far cyl. (X=2,L=113, r_c>=376), h4  (thr 6.2e-07) | 6.87e-08 OK x9.0 | 1.22e-08 OK x50.7 | 4.01e-09 OK x154.6 | 1.79e-09 OK x346.9 | 9.49e-10 OK x653.5 | 5.86e-10 OK x1058.0 | 3.66e-10 OK x1694.6 | 2.33e-10 OK x2666.6 |
| far cyl. (X=2,L=113, r_c>=376), N3  (thr 1) | 7.25e-02 OK x13.8 | 2.55e-02 OK x39.3 | 1.24e-02 OK x80.8 | 7.11e-03 OK x140.7 | 4.60e-03 OK x217.6 | 3.30e-03 OK x302.6 | 2.40e-03 OK x416.9 | 1.47e-03 OK x678.9 |
| pointwise kappa at r>=226  (thr 0.0089) | 1.75e-03 OK x5.1 | 1.23e-03 OK x7.2 | 9.52e-04 OK x9.4 | 7.72e-04 OK x11.5 | 6.49e-04 OK x13.7 | 5.66e-04 OK x15.7 | 4.92e-04 OK x18.1 | 3.96e-04 OK x22.5 |
| pointwise kappa_s at r>=226  (thr 7e-05) | 1.30e-06 OK x53.9 | 5.94e-07 OK x117.8 | 3.35e-07 OK x209.1 | 2.14e-07 OK x327.2 | 1.48e-07 OK x472.3 | 1.11e-07 OK x629.2 | 8.36e-08 OK x837.7 | 5.34e-08 OK x1311.2 |
| pointwise kappa_ss at r>=226  (thr 6.2e-07) | 1.27e-09 OK x489.7 | 3.73e-10 OK x1663.0 | 2.22e-10 OK x2791.1 | 2.57e-10 OK x2416.6 | 2.50e-10 OK x2483.0 | 2.34e-10 OK x2649.7 | 2.14e-10 OK x2895.9 | 1.80e-10 OK x3436.3 |


Cubic law check (-v/(y/k)^3, `v_samples` in the json): k = 2020: 39.5 (y = 50), 38.6 (112), 38.0 (169), 38.1 (253), 40.2 (350); k = 4046: 40.0, 39.7, 39.3, 38.6, 38.0; k = 8098: 40.2, 40.3, 40.2, 39.8, 39.2; k = 10132: 40.6 (y = 50), 40.5 (112), 40.4 (169), 40.1 (253), 39.6 (350).

## 5. Where it holds, scaling, extrapolation

First k (from the fits/closed form) at which each group holds by itself: |alpha| 5.6e3; sqrt(Lambda) 4.7e3; [14,60] 2.0e3; chart 3.3e3; cone slope 5.2e3; far slope 4.2e3; windows E0,E1 (to 253) 5.8e3; E2 7.7e3; **E3 1.11e4; E4 1.6e4**; sliding windows E2 1.07e4, E3 1.55e4, **E4 2.2e4**.

Exponents and extrapolated crossings for the windows that can fail (E ~ k^-p, linear fit of log E on log k over the five points k = 5061, 6081, 7018, 8098, 10132; closed form E_j = 28*9^j c3^2 (yb^5 - ya^5)/(5 k^6) with c3 = 40):

| window | level | exponent p (E ~ k^-p, fit over k>=5000) | E at k=10132 | threshold | fitted k_cross | closed-form k_cross (c3=40) |
|---|---|---|---|---|---|---|
| [112, 169] | E0 | 5.95 | 1.02e-09 | 2e-07 | 4174 | 4187 |
| [112, 169] | E1 | 5.94 | 9.09e-09 | 3e-06 | 3816 | 3845 |
| [112, 169] | E2 | 5.92 | 8.08e-08 | 3e-06 | 5501 | 5546 |
| [112, 169] | E3 | 5.89 | 7.13e-07 | 3e-06 | 7942 | 7998 |
| [112, 169] | E4 | 5.87 | 6.24e-06 | 3e-06 | 11482 | 11535 |
| [169, 253] | E0 | 5.93 | 7.51e-09 | 2e-07 | 5827 | 5854 |
| [169, 253] | E1 | 5.91 | 6.66e-08 | 3e-06 | 5321 | 5377 |
| [169, 253] | E2 | 5.89 | 5.86e-07 | 3e-06 | 7681 | 7754 |
| [169, 253] | E3 | 5.88 | 5.12e-06 | 3e-06 | 11095 | 11184 |
| [169, 253] | E4 | 5.89 | 4.43e-05 | 3e-06 | 15994 | 16130 |
| sliding max (last window [253, 380]) | E0 | 5.91 | 5.54e-08 | 2e-07 | 8155 | 8220 |
| sliding max (last window [253, 380]) | E1 | 5.90 | 4.87e-07 | 3e-06 | 7444 | 7549 |
| sliding max (last window [253, 380]) | E2 | 5.90 | 4.25e-06 | 3e-06 | 10740 | 10888 |
| sliding max (last window [253, 380]) | E3 | 5.95 | 3.68e-05 | 3e-06 | 15416 | 15703 |
| sliding max (last window [253, 380]) | E4 | 6.07 | 3.17e-04 | 3e-06 | 21747 | 22648 |


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
* Provenance: the percentages quoted in section 7 items 1-3 (up8/up16, fine-vs-main, hr-vs-main) are from my run notes; the raw res_*.json/table_*.md files are the authoritative numbers.

## 9. Files (all in `step4/m1/`)

`cont.py` (resumable run), `geo.py` (S1 on x <= 64 + C4 psi3~/h2~ cache `c4_cache.npz`), `prof.py` (symmetrized fit, graph2, ray profile, stage 1), `en.py` (Lambda, v, shell energies), `jets.py` (Lemma O), `drive.py SNAP UP TAG` -> `p_TAG.npz`, `res_TAG.json`; `make_report.py`, `tabs.py`, `scaling_md.py`, `scan.py` (tables); snapshots `snap_*.npz` (main, main2, big, res20, hr, fine), logs `run_*.log`. Old stage-1 results from before the Newton fix were discarded.
