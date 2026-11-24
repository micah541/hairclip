#!/usr/bin/env python3
"""Port of lemmaS_F2.ramp_level with NO doubles: exact Rational table coefficients, arb Abr / Ghat_ramp / YC_ramp (from weights_exact), L = ln 2 as an arb ball,
derivative polynomials of eta^2 exact.  Writes rigor/ramp_tables.json (arb upper endpoints of g['ramp'][k][m], j=0,1,2).  Same adaptive refinement as sup_on (the
refinement heuristic only chooses cells; every returned number is the maximum of arb upper endpoints over a cover of [0,1]).
Rerun: cd lemmaS && PYTHONHASHSEED=0 python3 rigor/ramp_arb.py   (needs weights_exact.json)"""
import os, sys, json, math
HERE = os.path.dirname(os.path.abspath(__file__)); PAR = os.path.dirname(HERE)
os.chdir(PAR); sys.path.insert(0, PAR); sys.path.insert(0, HERE)
import sympy as sp
from flint import arb, ctx
import arbio, weights_exact as WX
ctx.prec = 200
A = WX.A; Q_ = sp.Rational
STIR = {0: {0: 1}, 1: {1: 1}, 2: {1: 1, 2: 1}, 3: {1: 1, 2: 3, 3: 1}}
LN2 = arb(2).log(); t_ = sp.Symbol('t'); Ls = sp.Symbol('Ls'); gam = 1
Sx = t_**3*(10 - 15*t_ + 6*t_**2); e2s = Sx**2
E2F = [sp.lambdify(t_, sp.diff(e2s, t_, l), modules=[{}]) for l in range(8)]
def e2(t, l): return E2F[l](t)/LN2**l
Dq = {}
for m in (1, 2, 3, 4):
    R = sum(sp.binomial(m, l)*sp.diff(e2s, t_, l)/Ls**l*(-gam)**(m - l) for l in range(1, m + 1)); f = sp.expand(R); Dq[m] = []
    for q_ in range(4):
        Dq[m].append(sp.lambdify((t_, Ls), f, modules=[{}])); f = sp.expand(sp.diff(f, t_)/Ls - gam*f)
def yy(t): return 14*(t*LN2).exp()
X0 = arb(22490776515)/10**10
om = lambda t: (-(arb(gam))*((yy(t)/28).log())).exp()
eps = lambda t: yy(t)**-4
Wt = lambda t: (1 + (yy(t) + X0)**2)**-2
rho = lambda t: abs((e2(t, 2) - 2*(gam + 4)*e2(t, 1))*om(t)*eps(t))/2
def expand_ab(a, b, d):
    if a == b: return {(a, d): 1}
    if b == a + 1: return {(a, d + 1): arb(1)/2}
    out = {}
    for key, v in list(expand_ab(a + 1, b - 1, d).items()) + list(expand_ab(a, b - 1, d + 1).items()): out[key] = out.get(key, 0) + v
    return out
def ramp_level(j, T, Abr, Ghr, YCr, tmin=0.0):
    n = j + 2; pairs = [(a, b, m) for (a, b), m in T[j].items() if (a, b) != (n, n)]
    def cd(mons, d, t):
        tot = arb(0)
        for (m, i, r), co in mons.items():
            if m == 0: continue
            for p_ in range(d + 1):
                for q_ in range(d + 1 - p_):
                    u = d - p_ - q_; mult = math.factorial(d)//(math.factorial(p_)*math.factorial(q_)*math.factorial(u))
                    tot += A(co)*mult*4**p_*Abr[(i, r + u)]*abs(Dq[m][q_](t, LN2))*om(t)
        return tot*eps(t)
    wts = []
    for (a, b, mons) in pairs:
        if b <= j: wts.append((mons, {(a, 0): 1} if a == b else {(a, 0): arb(1)/2, (b, 0): arb(1)/2}))
        else: wts.append((mons, expand_ab(a, b, 0)))
    def q(k, t):
        tot = arb(0)
        for mons, w in wts:
            for (kk, d), v in w.items():
                if kk == k: tot += arb(v)*cd(mons, d, t)
        tot += Ghr[k]*eps(t)*e2(t, 0)*om(t)
        if k == j: tot += YCr[j + 1]*rho(t)
        return tot
    def sup_on(f, lo=tmin, n_=200, tol=0.02, depth=12):
        def ev(a_, b_): return abs(f(arb(a_).union(arb(b_)))).upper()
        def pt(a_, b_): return float(abs(f(arb((a_ + b_)/2))).mid())
        cells = [(lo + (1 - lo)*i/n_, lo + (1 - lo)*(i + 1)/n_) for i in range(n_)]
        lb = max(pt(a_, b_) for a_, b_ in cells); best = arb(0); stack = [(a_, b_, 0) for a_, b_ in cells]
        while stack:
            a_, b_, dd = stack.pop(); u = ev(a_, b_)
            if float(u) <= lb*(1 + tol) or dd >= depth:
                if u > best: best = u
                continue
            m_ = (a_ + b_)/2; lb = max(lb, pt(a_, b_)); stack += [(a_, m_, dd + 1), (m_, b_, dd + 1)]
        return best
    out = {}
    for k in range(j + 1): out[k] = {m: sup_on(lambda t, m=m, k=k: q(k, t)*yy(t)**(2*m - 1)/Wt(t)) for m in STIR[k]}
    return out
if __name__ == '__main__':
    WE = arbio.load(json.load(open(os.path.join(HERE, 'weights_exact.json')))['res']); T = WX.build_tables()
    Abr = WX.ab_exact(Q_(5, 100), Q_(9, 10), 6); res = {}
    for j in (0, 1, 2):
        Ghr = {int(k): arbio.up(v) for k, v in WE['ramp'][str(j)]['Ghat'].items()}; YCr = {int(k): arbio.up(v) for k, v in WE['ramp'][str(j)]['YC'].items()}
        res[j] = ramp_level(j, T, Abr, Ghr, YCr); print('level', j, {k: {m: '%.8g' % float(v) for m, v in d.items()} for k, d in res[j].items()}, flush=True)
    json.dump(arbio.dump(res), open(os.path.join(HERE, 'ramp_tables.json'), 'w'))
