# L-int variants (rigorous, python-flint; lint_cert_p.py)

A2=30.25, B2=38.25 for all (blend never failed, so 34.25/42.25 not needed). Far tail: see below.

| PX | XA | mu1 poly [XA,A2] (at x) | mu1 blend | mu1 rho [38.25,68.25] | runtime s |
|---|---|---|---|---|---|
| 3.0 | 12 | 0.2435 (12.00) | 0.5623 (30.25) | 1.0378 (38.25) | 31.1 |
| 2.5 | 12 | 0.2553 (12.00) | 0.4496 (30.25) | 1.0378 (38.25) | 32.5 |
| 2.5 | 10 | 0.1011 (10.00) | 0.4496 (30.25) | 1.0378 (38.25) | 31.8 |
| 2.0 | 12 | 0.2210 (12.00) | 0.3357 (30.25) | 1.0378 (38.25) | 36.1 |
| 2.0 | 10 | 0.1415 (10.00) | 0.3357 (30.25) | 1.0378 (38.25) | 47.8 |
| 1.5 | 12 | 0.1555 (12.00) | 0.2211 (30.25) | 1.0378 (38.25) | 48.3 |
| 1.5 | 10 | 0.1178 (10.00) | 0.2211 (30.25) | 1.0378 (38.25) | 49.8 |
| 1.5 | 8 | 0.0338 (8.00) | 0.2211 (30.25) | 1.0378 (38.25) | 47.9 |

Skipped (float mu1<0 on [XA,A2]): (3,6), (3,8), (3,10), (2.5,6), (2.5,8), (2,6), (2,8), (1.5,6). All minima are at x = XA for the poly region and at A2 (blend start) for the blend; rho min at B2.

## Findings / caveats
1. The supplied lint_cert.py does NOT reproduce 0.2435 as written: its fixed 0.005 grid gives enclosures of width ~3 near B2 (dependency problem in the smoothstep polynomial near t=1), so it prints mu1 >= -1.40 (x=38.245); its 'crude tail' block returns NaN (X up to 1e12, P_j balls NaN). lint_cert.json has different keys (poly/blend/rho_68), i.e. came from an older script version. In lint_cert_p.py pieces are adaptively bisected (until enclosure width <= 0.01, depth <= 14); the maths is unchanged. With PX=3, XA=12 this reproduces poly 0.2435 @12, blend 0.562, rho 1.038 (original json: 0.2435 / 0.5625 / 1.030).
2. For XA < 10.2491 (sigma = x + STIP < 8) the step-1 tail bounds (sigma >= 8) do not apply; the clamp max(.,8) in the original would be INVALID. For those pieces k,m,p,q are the hulls of the CAPD boxes in mc_prof.txt (x = 8 - tau - STIP; mapping checked numerically: q ~ 8-tau, p' = -kq, m' = p/4). This affects (PX,XA) = (1.5, 8) and all XA=10 rows (sigma_min = 7.75); their rigour rests on mc_prof.txt boxes being validated enclosures.
3. Far tail: lint_tail.py depends only on the rho-type weight (x >= 68.25 > B2) and on no PX/XA/A2/B2; it was re-run unchanged: mu1 >= 0.003 x^(4/3), ok=True (smallest 0.0093). It applies to every variant (also to B2=42.25). The crude tail block in lint_cert is NaN and should be ignored.
4. NaN guard: every piece lower bound is tested for NaN (float(x)!=float(x)); no NaN in any certified run. Powers are explicit products.

## Recommendation
- Minimising x^PX with mu1 >= ~0.1: (PX, XA) = (1.5, 10), mu1 = 0.118 (blend 0.221).
- (1.5, 8) gives only 0.034; (2, 10) gives 0.1415, (2.5, 10) 0.101 -- PX=1.5 costs ~0.02 vs PX=2 in rate but loses a factor x^0.5 pointwise.
- If a larger margin is needed: (2, 10), or keep XA=12 (PX=1.5 certified 0.1555; PX=2.5 0.2553).