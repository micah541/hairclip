#!/usr/bin/env python3
import pickle, json, math, sys, re
import sympy as sp
from flint import arb, ctx
ctx.prec = 100
sys.path.insert(0, '.')
import lemmaS_cflux as C
x = C.x
forms = pickle.load(open('cflux_forms.pkl', 'rb'))
Xmin = 57; Xmax = 60; X0 = 2.2490776515
# ---- arb-friendly lambdify (fractional powers -> arb symbols) ----
def arb_func(expr, extra=()):
    expr = sp.sympify(expr)
    reps = {}
    def rep(e):
        k = (e.exp.p, e.exp.q); nm = 'R%d_%d' % (k[0] if k[0] >= 0 else 10**6 - k[0], k[1]) + ('n' if k[0] < 0 else '')
        reps[nm] = arb(k[0])/arb(k[1]); return sp.Pow(e.base, sp.Symbol(nm))
    expr = expr.replace(lambda e: e.is_Pow and e.exp.is_Rational and not e.exp.is_Integer, rep)
    names = sorted(reps)
    f = sp.lambdify([x] + [sp.Symbol(n) for n in names], expr, modules=[{'exp': lambda z: z.exp(), 'sqrt': lambda z: z.sqrt()}])
    vals = [reps[n] for n in names]
    return lambda xv: arb(f(xv, *vals))

# ---- ingredients ----
zeta = C.zeta_expr
zders = {m: arb_func(sp.diff(zeta, x, m)) if m else arb_func(zeta) for m in range(0, 8)}
wr = {n: arb_func(C.wratio(n)) for n in range(0, 8)}
ydef = lambda xv: xv - arb(str(X0))
# ---- hypothesis box ----
MBOX = dict(M0=0.25, M1=0.4, M2=1.4, M3=5.8, M4=29.4, M5=150.0)   # |V^(j)_x| <= Mj / y^j   (from s-jets 0.4,1,2,4; see tex)
UMAX = 0.11
def abounds(box=MBOX, ymin=57 - X0):
    res, Fa, Fb = C.A_bounds(box, ymin, UMAX, rmax=4)
    return res, Fa, Fb

def sup_over_ramp(coef, Ab, nsub=300):
    """sup_{x in [57,60]} |coef|/w , coef = polynomial in zt^(m), w^(n), A_i^(r)."""
    coef = sp.expand(coef)
    zs = {sp.Derivative(C.zt, (x, m)) if m else C.zt: m for m in range(0, 8)}
    ws = {sp.Derivative(C.w, (x, n)) if n else C.w: n for n in range(0, 8)}
    syms = {}
    gens = []
    for atom in list(coef.atoms(sp.Derivative)) + list(coef.atoms(sp.Function)):
        gens.append(atom)
    gens = list(dict.fromkeys(gens))
    P = sp.Poly(coef, *gens)
    best = arb(0)
    for i in range(nsub):
        lo = arb(Xmin) + arb(Xmax - Xmin)*i/nsub; hi = arb(Xmin) + arb(Xmax - Xmin)*(i + 1)/nsub
        ball = arb((lo + hi)/2, (hi - lo)/2); yb = ydef(ball)
        tot = arb(0)
        for mon, c in P.terms():
            t = arb(abs(float(c))) if c.is_Number else None
            t = arb(sp.Rational(abs(c)).p)/arb(sp.Rational(abs(c)).q)
            for gsym, p in zip(gens, mon):
                if p == 0: continue
                if gsym in zs:
                    f = abs(zders[zs[gsym]](ball)); fac = arb(f.upper())
                elif gsym in ws:
                    n = ws[gsym]; f = abs(wr[n](ball)); fac = arb(f.upper())       # w^(n)/w
                else:   # A_i^(r)
                    nm = gsym.func.__name__ if isinstance(gsym, sp.Function) else gsym.expr.func.__name__
                    r = 0 if isinstance(gsym, sp.Function) else gsym.derivative_count
                    i_ = int(nm[1:]); fac = arb(Ab[(i_, r)])/(yb.lower())**r if r else arb(Ab[(i_, r)])
                t *= fac**p
            tot += t
        best = max(best, arb(tot.upper()), key=lambda z: float(z.upper()))
    return float(best.upper())

def psi_lower(xend, n=3000):
    """rigorous LOWER bound of psi(xend), W=e^{-psi}, psi=2ln(1+x^2) for x<=30.25, psi'=(1-s)4x/(1+x^2)+s*phi' beyond (s=S((x-30.25)/8))"""
    a = arb(30.25); psi0 = 2*(1 + a*a).log()
    S = lambda t: t**3*(10 - 15*t + 6*t*t)
    tot = arb(0); L = arb(xend) - a
    def dpsi(xb):
        t = (xb - a)/8
        s = arb(1) if float(t.lower()) >= 1 else S(t)
        if float(t.lower()) < 1 and float(t.upper()) > 1: s = arb(0, 1) + arb(0.5)
        phip = xb/(3*(1 + xb*xb)**(arb(1)/3))
        return (1 - s)*4*xb/(1 + xb*xb) + s*phip
    for i in range(n):
        lo = a + L*i/n; hi = a + L*(i + 1)/n
        b = arb((lo + hi)/2, (hi - lo)/2); d = dpsi(b)
        tot += arb(float(d.lower()))*(hi - lo)
    return psi0 + tot

if __name__ == '__main__':
    Ab, Fa, Fb = abounds()
    out = dict(box=MBOX, Fa=Fa, Fb=Fb, Abar={f'{k[0]},{k[1]}': v for k, v in Ab.items()})
    # coefficient sups
    G = {}
    for name in ('K0', 'K2'):
        G[name] = {}
        for (a_, b_), c in forms[name].items():
            G[name][f'{a_},{b_}'] = sup_over_ramp(c, Ab)
            print(name, (a_, b_), G[name][f'{a_},{b_}'], flush=True)
    out['Gamma'] = G
    pl = psi_lower(57)
    W57 = float(arb(-pl).exp().upper())
    out['W57_upper'] = W57; out['intW_ramp_upper'] = 3*W57
    out['psi57_lower'] = float(pl.lower())
    json.dump(out, open('cflux_raw.json', 'w'), indent=1)
    print('W(57) <=', W57)
