# Weighted interpolation / Sobolev constants and cutoff commutator

Everything below is produced by `sob_constants.py` (arb balls, cells of width 2^-12, so all breakpoints 30.25, 38.25 are exact grid points; `sob_constants.json`). Constants are rounded outward to 4 significant digits. `float_check.py` is an independent numpy/scipy float check (agrees to all printed digits).

Notation: rho = e^-phi, phi=(1+x^2)^(2/3)/4, X=34.25. Facts used (certified on cells): phi'>0 on (0,oo), phi''>0 (min phi'' >= 0.01054 on [0,X]), so phi' increasing, rho decreasing, Phi1(x)=sup_[x,X] phi' = phi'(X). psi'>0 on [6,80] (min >= 0.195), so W decreasing.

## (I1) ||f'||^2 <= ||f|| ||f''|| + c1 ||f||^2,  c1 = (1/2) sup (phi'^2 - phi'')_+
Proof: int rho f'^2 = -int rho f f'' + int phi' rho f f' (boundary terms vanish: f'(0)=0, f(X)=0); int phi' rho f f' = (1/2) int phi' rho (f^2)' = (1/2) int (phi'^2 - phi'') rho f^2. Then Cauchy-Schwarz. Also valid for f'' (even): ||f'''||^2 <= ||f''|| ||f''''|| + c1 ||f''||^2.
**c1 = 0.5804** (sup (phi'^2-phi'')_+ = 1.161, attained at x=X).

## (I2) rho g^2(x) <= Phi1(x) ||g||^2_[x,X] + 2 ||g||_[x,X] ||g'||_[x,X]   (g = f or f')
Proof: rho g^2(x) = -int_x^X (rho g^2)' = int_x^X (phi' g^2 - 2 g g') rho; bound phi' <= Phi1(x), Cauchy-Schwarz.
Hence sup_region g^2 <= K (Phi1(X) ||g||^2 + 2||g|| ||g'||), K = 1/min rho = exp(phi(right end)).
**Phi1(X) = phi'(X) = 1.083.**

| region | [0,1] | [1,2] | [2,4] | [4,8] | [8,12] | [12,20] | [20,34.25] |
|---|---|---|---|---|---|---|---|
| K | 1.488 | 2.078 | 5.222 | 56.92 | 992.4 | 8.019e5 | 1.212e12 |

## (I3) Same for W = e^-psi on [a,R], a=6, v supported in (a,R)
(I1_W) ||v'||_W^2 <= ||v||_W ||v''||_W + c1W ||v||_W^2, c1W = (1/2) sup_[a,R] (psi'^2-psi'')_+ (same proof, both boundary terms vanish).
(I2_W) W v^2(x) <= Psi1(x) ||v||^2_[x,R] + 2||v||_[x,R]||v'||_[x,R], Psi1(x) = sup_[x,R] psi' (integrate from x to R).
psi' = (1-s)A + sB, A=6x/(1+x^2), B=(1/3)x(1+x^2)^(-1/3); psi'' = (1-s)A' + sB' + s'(B-A), A'=6(1-x^2)/(1+x^2)^2, B'=phi'', s' = S'(t)/8, S'(t)=140 t^3 (1-t)^3, t=(x-30.25)/8.
**Normalisation (K_W depends on it): W = (1+x^2)^-3 exactly for x <= 30.25**, psi(x)=3 log(1+x^2) there, psi continued by integrating psi' (rigorous acb integral on [30.25,38.25], closed form beyond: psi = psi(38.25)+phi(x)-phi(38.25)). psi(38.25)-psi(30.25) = 5.1416.

| R | 40 | 60 | 80 |
|---|---|---|---|
| c1W | 0.6448 | 0.8477 | 1.029 |
| Psi1 = sup_[x,R] psi' (same for every region start in 6..60, attained at x=R) | 1.140 | 1.305 | 1.437 |

| region | [12,20] | [20,30] | [30,40] | [40,60] | [60,80] |
|---|---|---|---|---|---|
| K_W = exp(psi(right)) | 6.449e7 | 7.315e8 | 9.519e11 | 4.242e22 | 3.549e34 |

(Use region [l,r] only if r <= R.)  Note psi' is not monotone (A decreasing, B increasing); sup over [x,R] was taken on cells, and equals max(psi'(x), psi'(R)) here.

## (I4) int_0^b rho |f'|^2 |f''| <= sqrt(K_b) (Phi1(X)||f'||^2 + 2||f'|| ||f''||)^(1/2) ||f'|| ||f''||,  K_b = e^{phi(b)}
Proof: <= sup_[0,b]|f'| int rho |f'||f''| <= sup|f'| ||f'|| ||f''||; sup_[0,b] |f'|^2 <= K_b(...) by (I2) with g=f'.
**K_2 = 2.078, K_12 = 992.4** (Phi1(X)=1.083).

## Part (b): cutoff commutator
chi(sigma)=S((32-sigma)/3); d/dx = d/dsigma, so chi^(j) = (-1/3)^j S^(j)(u).
L = -d^4 - k^2 d^2 - (5km+q/4) d + c0.  Zeroth-order term commutes. Using Leibniz:

[L,chi] w = -4chi' w''' - 6chi'' w'' - (4chi''' + 2k^2 chi') w' - (chi'''' + k^2 chi'' + (5km+q/4) chi') w.

sup_[0,1]|S^(j)| (cells, 2^19, vs exact critical points): 2.18757 (2.1875=35/16), 7.51427 (7.51319), 52.5031 (52.5), 840.016 (840).
sup|chi^(j)|, j=1..4 <= **0.7292, 0.8350, 1.945, 10.38**.
With |k|,|m|,|p| <= 1e-20, |q| <= 32+1e-20 (sigma<=32), |5km+q/4| <= 8.001:
**C3 = 4 sup|chi'| = 2.917; C2 = 6 sup|chi''| = 5.010; C1 = 4 sup|chi'''| + 2k^2 sup|chi'| = 7.779; C0 = sup|chi''''| + k^2 sup|chi''| + 8.001 sup|chi'| = 16.21.**
Strip: supp chi' = sigma in [29,32], so x in [29-tip_hi, 32-tip_lo] = [31.2490776505, 34.2490776516] inside [31.249,34.250]. Length: 3 exactly (support of chi'); 3.001 for the stated strip.
rho_max = rho(31.249) <= **2.024e-11** (same with tight endpoint). rho_max*len <= 6.073e-11 (len=3.001; about 6.07e-11 with len=3).
|<[L,chi]w,g>_rho| <= rho_max * len * (C0 sup|w| + C1 sup|w'| + C2 sup|w''| + C3 sup|w'''|) * sup|g|.
