# Plan B (R = 260) far-field recertification. Scripts: pb_common, pb_lint, pb_sob, pb_n3, pb_res, pb_nlin, pb_mod, pb_ramp, collate_planB; float_mu*.py (coarse float scans), logs *.log. Result: planB_constants.json
All interval numbers: python-flint balls, outward rounded, NaN asserts, C^4 functions psi3~, h2~ (c4fun.py), |alpha|<=1e-3, jet box 1e-2. Boxes: x from 10.25 (sigma>=8.0009) to 260.
Float scans (float_mu.py, float_mu2.py) done first: no show-stopper.

## Certified results (R=260, P as specified)
* MU1 (a): global min 0.11515 at x=60 (profile: 0.1545 at 10.25, 0.31-0.34 on [20,40], dip to 0.115 at 55-60 where P drops, then rises 0.1166 (70) -> 0.1227 (250), -> 1/8). per-region in json.
* MU2 (b): plain W: >= 0.6151 (= mu1+1/2). Level-2 weight (1+(x/60)^4)W: >= 0.12813 (min at x=260; 0.1867 at 100, 0.1318 at 200). e1<=5.7e-6, e0<=5.1e-6 (as before; zero beyond x~20).
* KW, PSI1W, C1W (d): c1W = 0.09362 (all R), Psi1 = sup psi' from 10.25: 0.3872 (decreasing; from 60: 0.03334; from 250: 0.008), KW: 7.03e4, 1.61e5, 2.56e6 ([20,40]), 1.59e7 ([40,80]), 6.4e7 ([80,160]), 1.68e8 ([160,260]).
  Level-2 weight W2 = W*(1+(x/60)^4) is INCREASING beyond x~60 (psi2' < 0): I2 uses (psi2')_+ and KW2 = 1/min W2 (4.5e6 over [10.25,260]); c1W2 = 0.0935.
* FnW (c): ||(L-lt)psi3~||_W = 7.19e-9, ||(L-2lt)h2~+Q2||_W = 2.13e-8 (derivative level with zeta: 1.80e-8, 7.94e-8); C3W (alpha^3 coeff) 7.549 at alpha<=1e-3 (0.181 at 3e-4; 0.189 at 5e-4; 0.353 at 7e-4);
  derivative level C3W' (zeta) 0.1911 (0.0023 at 3e-4).
* T3 coefficients per region incl. [250,260]: json T3_per_region; MOD (V5/F2/F3 analogues + ||zeta q psi3''||_W = 0.518, c0' norm 9.9e-12): json MOD_*.
* RAMP_linear_deviation (add-on, [10.25,16.25]): w'''' coefficient 7.7e-4 (alpha 1e-3, w-box 0) ; 8.5e-4 (eta 1e-3); 1.76e-3 (eta 1e-2) -> below 2e-3. Table + a,b split in json
  (1-J^-4 ~ 2u0'^2 dominates: quadratic in alpha, b~580; linear part a~0.19).

## Worse than 2x vs R=60 design
* C3W level 0: 0.18 -> 7.55; C3W': 0.0026 -> 0.19. Cause: u0' = alpha psi3~' reaches 0.37-0.42 on [200,260] (sup|u0'| 0.085 at 60); the alpha-series remainder (c4, c5 ~ x^4..x^6) dominates for alpha_max=1e-3; identical to the old value (0.18) for alpha_max <= 3e-4. The actual forcing is ~1e-8; only the C3 |alpha|^3 form is lossy. Almost all of C3W' sits on [200,260].
* T3 w''^3 3.3 -> 341, w'' sq 0.015 -> 1.86, w''w''' 0.98 -> 13, w'''' coeff 0.0185 -> 0.69, w''' 0.0015 -> 0.024 on [250,260] (|u0'| ~ 0.42).  On [10.25,60] unchanged.
* ||zeta psi3~'||_W 0.26 -> 0.51, ||zeta h2~'||_W 1.8 -> 24.6, A1' 0.27 -> 0.60 (W decays only like x^-2).  sup A1'/alpha 542, h2' 3.6e4 at x~250.
* K_W at the outer end: 3.6e19 (x=60) -> 1.7e8 at 260 (better).  sup|1-q psi'|: 77 -> 2.78 (better); sup sqrt(psi'^2+k^2) 1.31 -> 0.387.
## (e) Outer ramp [250,260] -- NOT as negligible as stated
W(260)/W(10.25) = 6.7e-5 (not 9e-9; W(260) itself is 5.9e-9, W(250)=6.4e-9, int_250^260 W <= 6.2e-8). Coefficient sups on [250,260] are the large ones above (w''^3 341, w''^2 1.86), and the C3W' ramp part is 0.107. Flux = sqrt(W) 8e-5 x coefficients: the cubic coefficient (341) x sqrt(W)... not tiny; needs sup|w''| control there. d001 0.69 on the ramp also means 1-J^-4 is O(1) there (u0'~0.42).
## Caveats
Residual/T3 boxes: widths 0.005-0.2 (single resolution, conservative). Tail beyond 260 not certified (zeta=0). mu1 asymptote reached only ~0.1227 at 250 (< 1/8). Level-2 e1,e0 bounds computed with profile tail bounds only. 
