# Beta-loop, exact alpha-linear coefficients (rigorous)

**Definitions and sign check.** x^T(G+E) = R gives x = H R - H E^T x, H = (G^-1)^T. To first order in alpha, beta = -alpha H33 (e13 a + e23 b) (H31, H32 ~ 1e-5, contribution below 1e-4 relative; the E33 feedback factor 1/(1 + H33 E33) = 1 +- 2e-3 at alpha = 1e-3). Here E_13 = -<chi A1(u), f3> gives e13 = <gA, f3>. E_23 = -<chi B1(u), f3> gives e23 = -<gB, f3>, with gA = chi(psi - q psi'), gB = chi psi' cos(theta). The energy term -beta<chi psi3, w~> = + alpha H33 (e13<w,ha> + e23<w,hb>)<w,g0>, g0 = chi psi3^. This matches your sign. a = <w,ha>, b = <w,hb> to leading order.

**Pairings with f3** (`bloop_pair.py`: tight Taylor-model integrals of rho*v, times the hull of F from outA5, bump a_c eta1 + b_c eta3 added by L1 bound; radius 2-5%, all from the F hull):
- <gA,f3> in [-0.026351, -0.024995]
- <gB,f3> in [-0.003617, -0.003445]
- <g0,f3> = G33 in [0.031717, 0.032994] (psi3ip: [0.031107, 0.033604])
- <chi h2,f3> in [0.015106, 0.016884], so E33/alpha = 2<chi h2,f3> in [0.0302, 0.0338]
- e13 = -0.02567, e23 = +0.003531 (midpoint F, the point-F values, non-rigorous). These are 0.347 -> 0.026 and 0.0163 -> 0.0035 compared to the crude L1 bounds.
- H33 = 1/G33 in 31.0 +- 0.7. So H33 e13 in -0.80 +- 0.043 and H33 e23 in 0.109 +- 0.006.
- alpha^2 remainders (L1, alpha <= 1e-3): |E13 - alpha e13| <= 0.308 alpha^2, |E23 - alpha e23| <= 0.088 alpha^2.

**Forms on the constrained space** (sup, inf of Q/alpha per ||w||^2, rigorous; `bloop.json`). The coefficients are enclosed in the interval boxes above; the sup over the box is taken at vertices, valid because lambda_max is convex in the linear parameters:

| form | {f1,f2}^perp | {f1,f2,f3}^perp |
|---|---|---|
| Q1 alone | [-13.55, 11.33] | [-13.47, 11.02] |
| beta-loop piece alone | [-11.56, 9.30] | [-11.16, 9.15] |
| **Qfull = Q1 + loop** | **[-6.01, 1.72]** | **[-6.19, 1.90]** |

So |Qfull| <= 6.19 |alpha| ||w||^2 on {f1,f2,f3}^perp (R1 = 13.55 -> 6.2). Float eigenvalues of the midpoint form are about [-5.74, 1.34], so the rigorous bound is within about 8% of the heuristic. The loop piece is a large term (size 9-11) that largely cancels against Q1; the crude 300-480 estimate is far off.

**Method.** The 8x8 Gram (Lf1, Lf2, gA, gB, g0, f1, f2, f3) comes from `gram7_g3.json` (mod_gram2-type rigorous, 5 vectors plus f1, f2), the f3 row from `pair_f3.json`, <f1,f3> and <f2,f3> from lyap_consts2 (psif3), and ||f3|| from lyap_extra. Schur-complement projections (first off f1, f2, then f3), the majorant G+ = mid + diag(rad), and arb Cholesky bounds as in `../mod_forms.py`.

**Caveats.**
- The f3 pairings inherit the 2-5% hull uncertainty of F (outA5), so the form bounds carry roughly +-0.2.
- The intervals on H33 e13 and H33 e23 are treated as independent although they share F; the true uncertainty is smaller (see the midpoint F heuristic: sup 1.58, inf -5.97).
- Item (3), the alpha=0 loop form <w,(L*-lamt)f3><w,g0>/G33 on {f1,f2,f3}^perp, was NOT done: it needs derivative jets of F and its Gram with (L*-lamt)f3. The existing bound eps1 n3/G33 = 0.01607 stands.
- The w-dependent parts of E (|alpha| A, W1 terms), the +H13-type j=3 contribution to a,b (<= 0.005 A) and the h2-part of E33 times beta0 are NOT in Qfull. The first are covered by the lemma.
- The same sign conventions as MOD_FORMS (cos theta = St C - Ct S).

Scripts: `bloop_pair.py`, `bloop_rem.py`, `bloop_forms.py`, `bg2.py` (mod_gram2 variant), `mg_local.py`, `bloop_hull.py` (crude first pass). Runtime about 1.5 min total.
