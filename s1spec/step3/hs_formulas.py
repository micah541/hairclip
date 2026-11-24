"""Closed forms for the conical-end analysis of the Jacobi system of H_sym (weight phi = beta(1+x^2)^{2/3}).

Scaled coordinates  y = D(x) Y,  D = diag(1, P1, P1^2, P1^3),  P_j = phi^{(j)}(x).  On the ray (k = m = p = 0,
q = sigma = x + stip) the scaled matrix is  P1 * N(a, c) - D^{-1} D',
    N = [[0,1,0,0],[0,0,1,0],[0,-a,1,-1],[-c,0,0,1]],  a = A0/P1^2,  c = (V0 - lam)/P1^4,
    D^{-1}D' = diag(0, r, 2r, 3r),  r = P2/P1.
In the chart u = x^{-1/3} all of  P1*u, r, a, c, da/dx, dc/dx, 1/P1, 1/P1^3  are given by expressions
analytic in u on [0, u0] (checked below against direct evaluation in x).
"""
import sympy as sp

x, u, stip, lam = sp.symbols('x u stip lam', real=True)
beta = sp.Rational(1, 4)
phi = beta*(1 + x**2)**sp.Rational(2, 3)
P = [sp.diff(phi, x, j) for j in range(0, 6)]          # P[1..5]
sigma = x + stip
A0 = 2*(P[1]**2 - P[2])
r4 = P[1]**4 - 6*P[1]**2*P[2] + 3*P[2]**2 + 4*P[1]*P[3] - P[4]
V0 = -r4/2 - sigma*P[1]/8 + sp.Rational(3, 8)
a_x = A0/P[1]**2
c_x = (V0 - lam)/P[1]**4
r_x = P[2]/P[1]
da_x = sp.diff(a_x, x)
dc_x = sp.diff(c_x, x)

def to_u(e):
    """substitute x = u^{-3}, write powers of (1+u^-6) as u^{..}(1+u^6)^{..}"""
    e = e.subs(x, u**-3)
    e = sp.powsimp(sp.powdenest(e, force=True), force=True)
    return e

w6 = sp.symbols('w6', positive=True)   # stands for (1+u^6)

def u_form(expr_x):
    # rewrite (1 + x^2)^alpha = x^{2 alpha} (1 + x^{-2})^alpha with x^{-2} = u^6
    s = sp.symbols('s', positive=True)
    e = expr_x.subs(1 + x**2, s)
    e = e.subs(s, x**2*w6)
    e = e.subs(x, u**-3)
    e = sp.simplify(sp.powsimp(sp.expand_power_base(e, force=True), force=True))
    return e

QUANT = dict(P1u=P[1]*x**sp.Rational(-1, 3), r=r_x, a=a_x, c=c_x, dadx=da_x, dcdx=dc_x,
             invP1=1/P[1], invP13=1/P[1]**3)

if __name__ == "__main__":
    import mpmath as mp
    out = {}
    for name, e in QUANT.items():
        eu = u_form(e)
        out[name] = eu
        # numeric check at x = 50
        xv = 50.0; uv = xv**(-1/3)
        d = float(e.subs({x: xv, stip: -2.249, lam: -0.05}))
        f = float(eu.subs({u: uv, w6: 1 + uv**6, stip: -2.249, lam: -0.05}))
        print(f"{name:6s}: direct {d:+.15e}  u-form {f:+.15e}  rel.diff {abs(d-f)/max(1e-300,abs(d)):.1e}")
    import pickle
    pickle.dump({k: sp.srepr(v) for k, v in out.items()}, open('hs_uforms.pkl', 'wb'))
    for k, v in out.items():
        print(k, '=', sp.sstr(v)[:300])
