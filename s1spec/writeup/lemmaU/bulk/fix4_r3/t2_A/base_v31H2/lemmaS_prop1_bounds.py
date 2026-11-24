#!/usr/bin/env python3
"""Certified constants for Prop. 1 (energy inequalities), from the machine-verified forms prop1_forms.pkl.
Triangular system E_0 <- E_1 <- E_2 (level j couples only to k<=j).  Weight omega = vartheta(s) e^{-gam (s-s0)}, vartheta = S(t)^4 on the inner
ramp s in [ln19, ln24], 1 beyond; eps = e^{-4s}.  See lemmaS.tex Sec. 3."""
import pickle, json, math, sys, itertools
import sympy as sp
from flint import arb, ctx
ctx.prec = 100
sys.path.insert(0, '.')
from lemmaS_interval import sup_arb
import lemmaS_prop1_sym as P
from sympy.functions.combinatorial.numbers import bell as bellpoly
s = P.s
res = pickle.load(open('prop1_forms.pkl', 'rb'))

# ---------- (n-1,n) pair -> diagonal, re-verified ----------
def diag_reduce(Q):
    terms = dict(Q)
    while True:
        off = [k for k in terms if k[0] != k[1] and terms[k] != 0]
        if not off: return {k: sp.expand(v) for k, v in terms.items() if v != 0}
        a, b = max(off, key=lambda k: (k[1], k[1] - k[0])); c = terms.pop((a, b))
        def add(k, v): k = tuple(sorted(k)); terms[k] = terms.get(k, 0) + v
        if b == a + 1: add((a, a), -sp.Rational(1, 2)*P.D(c))
        else: add((a + 1, b - 1), -c); add((a, b - 1), -P.D(c))

def prepare(j):
    Q, Qtr = res[j]
    Qd = diag_reduce(Q)
    ok = P.is_total_derivative(sp.expand(P.to_u(P.form_expr(Q)) - P.to_u(P.form_expr(Qd))))
    return Qd, Qtr, ok

def table(Q):
    """c_ab / e^{-4s} as {(a,b): {(m,i,r): coef}}  for  om^(m) * A_i^(r)"""
    out = {}
    for (a, b), c in Q.items():
        e = sp.expand(sp.simplify(c.subs(P.ep, sp.exp(-4*s)).doit()*sp.exp(4*s)))
        tab = {}
        P_ = sp.Poly(e, *(list(e.atoms(sp.Derivative)) + [f for f in e.atoms(sp.Function) if not isinstance(f, sp.Derivative) and f.func.__name__ in ('om', 'A0', 'A1', 'A2')]))
        for mon, co in P_.terms():
            m = i = r = None; m = 0
            for gsym, p in zip(P_.gens, mon):
                if p == 0: continue
                assert p == 1
                if isinstance(gsym, sp.Derivative): nm, k = gsym.expr.func.__name__, gsym.derivative_count
                else: nm, k = gsym.func.__name__, 0
                if nm == 'om': m = k
                else: i, r = int(nm[1:]), k
            tab[(m, i, r)] = tab.get((m, i, r), 0) + float(abs(co))
        out[(a, b)] = tab
    return out

# ---------- hypothesis box, A bounds ----------
BOX = dict(M0=0.6, M1=0.8, M2=1.2, M3=2.0, M4=3.0, M5=4.0, M6=5.0, M7=6.0, Umax=0.5, p=1.0246)
def sup_rat(fk, M0, N=600):
    p = sp.Symbol('p'); f = sp.lambdify(p, fk, modules=[{}]); best = 0.0
    for i in range(N):
        lo = arb(M0)*(2*i - N)/N; hi = arb(M0)*(2*i + 2 - N)/N
        v = abs(f(arb((lo + hi)/2, (hi - lo)/2)))
        best = max(best, float(arb(v).upper()))
    return best
def A_bounds_s(B=BOX, rmax=6):
    p = sp.Symbol('p'); a = (1 + p**2)**-2; b = -3*p*(1 + p**2)**-3
    Fa = [sup_rat(sp.diff(a, p, k), B['M0']) for k in range(rmax + 3)]
    Fb = [sup_rat(sp.diff(b, p, k), B['M0']) for k in range(rmax + 3)]
    Vj = [B['M%d' % i] for i in range(8)]; pp = B['p']; Um = B['Umax']
    def faa(F, r):
        if r == 0: return F[0]
        return sum(F[k]*float(bellpoly(r, k, list(Vj[1:r - k + 2]))) for k in range(1, r + 1))
    Ab = {}
    for r in range(rmax + 1):
        Ab[(2, r)] = faa(Fa, r)
        Ab[(1, r)] = faa(Fa, r) + sum(math.comb(r, q)*faa(Fb, q)*(Vj[r - q + 1] + pp**(r - q + 1)*Um) for q in range(r + 1))
        Ab[(0, r)] = sum(math.comb(r, q)*(faa(Fa[1:], q)*(pp**2 - pp)*pp**(r - q)*Um + faa(Fb[1:], q)*(2*pp)**(r - q)*pp**2*Um**2) for q in range(r + 1))
    return Ab, Fa, Fb

# ---------- ramp quantities ----------
t = sp.Symbol('t')
Sx = t**3*(10 - 15*t + 6*t**2)
def theta_sup(r):
    e = sp.cancel(sp.diff(Sx**4, t, r)/Sx**2)
    f = sp.lambdify(t, e, modules=[{}])
    return float(sup_arb(lambda tt: abs(f(tt)), 0, 1).upper())
SA, SB = math.log(19), math.log(24); DS = SB - SA

def constants(gam=1.0, B=BOX, La=19.0, verbose=False):
    DS = math.log(24.0/La)
    Ab, Fa, Fb = A_bounds_s(B)
    amin = (1 + B['M0']**2)**-2
    Th = [1.0] + [theta_sup(r)/DS**r for r in (1, 2, 3, 4)]
    OmR = [0.0] + [sum(math.comb(m, r)*Th[r]*gam**(m - r) for r in range(1, m + 1)) for m in (1, 2, 3, 4)]
    OmC = [gam**m for m in range(5)]
    kk = gam + 4
    S2 = Sx**2
    e = sp.cancel(sp.diff(S2, t, 2)/DS**2 - 2*kk*sp.diff(S2, t)/DS + kk**2*S2)
    f = sp.lambdify(t, e, modules=[{}]); D2 = float(sup_arb(lambda tt: abs(f(tt)), 0, 1).upper())
    S4 = Sx**4
    f1 = sp.lambdify(t, sp.diff(S4, t)/DS, modules=[{}]); f2 = sp.lambdify(t, sp.diff(S4, t, 2)/DS**2, modules=[{}])
    rho2 = float(sup_arb(lambda tt: abs(f2(tt)), 0, 1).upper()) + 2*kk*float(sup_arb(lambda tt: abs(f1(tt)), 0, 1).upper())
    out = dict(amin=amin, OmR=OmR, Theta=Th, D2=D2, rho2=rho2, Ab={f'{k[0]},{k[1]}': v for k, v in Ab.items()})
    G = {}
    for j in (0, 1, 2):
        Q, Qtr, ok = prepare(j); assert ok
        tab = table(Q); n = j + 2
        DC = {k: sum(co*OmC[m]*Ab[(i, r)] for (m, i, r), co in tab.get((k, k), {}).items()) for k in range(n + 1)}
        DR = {k: sum(co*OmR[m]*Ab[(i, r)] for (m, i, r), co in tab.get((k, k), {}).items()) for k in range(n + 1)}
        Bc, Br = 1.5*amin, 0.4*amin
        GC = {k: DC[k] for k in range(j + 1)}; GR = {k: DR[k] for k in range(j + 1)}
        GC[j] += DC[j + 1]**2/(4*Bc) + DC[j + 1]*kk**2/2
        GR[j] += DR[j + 1]**2/(4*Br) + DR[j + 1]*D2/2 + DC[j + 1]*rho2/2
        out[f'GhatC_{j}'] = {str(k): v for k, v in GC.items()}
        out[f'GhatR_{j}'] = {str(k): v for k, v in GR.items()}
        out[f'DC_{j}'] = {str(k): v for k, v in DC.items()}; out[f'DR_{j}'] = {str(k): v for k, v in DR.items()}
        G[j] = (GC, GR)
    return out, G

if __name__ == '__main__':
    out, G = constants()
    json.dump(out, open('prop1_consts_raw.json', 'w'), indent=1, default=str)
    print(json.dumps({k: out[k] for k in ('amin','OmR','D2','rho2','GhatC_0','GhatC_1','GhatC_2','GhatR_0','GhatR_1','GhatR_2')}, indent=1))
