#!/usr/bin/env python3
"""Sanity check (NOT part of the proof; floating point, mpmath 30 digits): for a test function g, the reduced forms of
flux_reduce.py integrate to the same K0, K2 as the defining integrals  int W v Hc  and  int (W v'')'' Hc  over [57, 60]."""
import pickle, sympy as sp, mpmath as mp, sys, os, json
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import flux_reduce as FR
mp.mp.dps = 30
x = FR.x
u = (60 - x)/3
zeta = u**4*(35 - 84*u + 70*u**2 - 20*u**3)
Wx = sp.exp(-((1 + x**2)**sp.Rational(2, 3)/4 - sp.Rational(1, 4)*(1 + 57**2)**sp.Rational(2, 3)))   # W up to a constant factor
Cx = (x - sp.Rational(2249, 1000))/4
Kx = sp.Rational(1, 10**3)*sp.sin(x)**2          # a nonzero 'k^2' to exercise the K terms
res = {}
for gname, gx in (('poly-exp', sp.exp(sp.Rational(1, 7)*x)*sp.cos(x) + (x - 58)**3/10), ('trig', sp.sin(3*x) + sp.Rational(1, 3)*sp.cos(x/2))):
    forms = pickle.load(open(os.path.join(HERE, 'flux_forms.pkl'), 'rb'))
    subsmap = lambda e: e.subs({FR.Wf: Wx, FR.Zf: zeta, FR.Cf: Cx, FR.Kf: Kx}).doit()
    Hc = sum(c*sp.diff(gx, x, j) for j, c in enumerate([0]*0) ) if False else None
    v = zeta*gx
    Hcx = 4*sp.diff(zeta, x)*sp.diff(gx, x, 3) + 6*sp.diff(zeta, x, 2)*sp.diff(gx, x, 2) + 4*sp.diff(zeta, x, 3)*sp.diff(gx, x) + sp.diff(zeta, x, 4)*gx \
          + Kx*(2*sp.diff(zeta, x)*sp.diff(gx, x) + sp.diff(zeta, x, 2)*gx) + Cx*sp.diff(zeta, x)*gx
    K0def = Wx*v*Hcx; K2def = sp.diff(Wx*sp.diff(v, x, 2), x, 2)*Hcx
    for lev, defn in ((0, K0def), (2, K2def)):
        red = 0
        for key, srep in forms[lev].items():
            a, b = map(int, key.split(','))
            c = eval(srep, {**sp.__dict__, 'x': x})
            red += subsmap(c)*sp.diff(gx, x, a)*sp.diff(gx, x, b)
        f1 = sp.lambdify(x, defn, 'mpmath'); f2 = sp.lambdify(x, red, 'mpmath')
        I1 = mp.quad(f1, [57, 58, 59, 60]); I2 = mp.quad(f2, [57, 58, 59, 60])
        res['%s/K%d' % (gname, lev)] = (mp.nstr(I1, 20), mp.nstr(I2, 20), mp.nstr(abs(I1 - I2)/max(abs(I1), mp.mpf('1e-30')), 5))
        print(gname, 'K%d' % lev, 'definition', mp.nstr(I1, 15), 'reduced', mp.nstr(I2, 15), 'rel diff', mp.nstr(abs(I1 - I2)/abs(I1), 3), flush=True)
json.dump(res, open(os.path.join(HERE, 'sanity_quad.json'), 'w'), indent=1)
