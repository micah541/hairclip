#!/usr/bin/env python3
"""Symbolic derivation + machine verification of the energy identities (Prop. 1).
Level-j dissipative form:  T_j = -2 om z_j d^j( eps L G ),  G = A2 z2 + A1 z1 + A0 z0,  L = d^2 - 5 d + 6,  z_k = u^(k).
Step 1: P_j = -2 (-1)^j d^j(om z_j) eps L G   (IBP j times)           -> check T_j - P_j = d/ds(...)  [Euler operator = 0]
Step 2: Q_j = reduction of P_j by IBP to indices <= j+2              -> check P_j - Q_j = d/ds(...)
Transport: T^tr_j = -2 om z_j d^j(c z1) -> Q^tr_j (indices <= j)     -> same check.
Output: pickle of coefficient tables  c_ab = e^{-4s} * sum coef * om^(m) * A_i^(r)."""
import sympy as sp, pickle, sys, itertools
from sympy.calculus.euler import euler_equations
s = sp.Symbol('s')
u = sp.Function('u')(s)
om = sp.Function('om')(s); ep = sp.Function('ep')(s); cc = sp.Function('c')(s)
A = {i: sp.Function('A%d' % i)(s) for i in (0, 1, 2)}
Z = sp.symbols('z0:12')
def D(e, k=1): return sp.diff(e, s, k)
def zfun(k): return u if k == 0 else sp.Derivative(u, (s, k))
def to_u(expr): return expr.subs({Z[k]: zfun(k) for k in range(11, -1, -1)})
def from_u(expr):
    for k in range(11, 0, -1): expr = expr.subs(sp.Derivative(u, (s, k)), Z[k])
    return expr.subs(u, Z[0])
def Lop(f): return D(f, 2) - 5*D(f) + 6*f
Gexpr = A[2]*D(u, 2) + A[1]*D(u) + A[0]*u

def is_total_derivative(expr):
    """expr quadratic in u; check Euler operator wrt u vanishes (=> total s-derivative)"""
    if sp.expand(expr) == 0: return True
    eqs = euler_equations(expr, [u], s)
    if not eqs: return True
    return sp.simplify(sp.expand(eqs[0].lhs)) == 0

def collect_terms(expr):
    e = sp.expand(from_u(sp.expand(expr)))
    P = sp.Poly(e, *Z)
    out = {}
    for mon, coef in P.terms():
        idx = [k for k, p in enumerate(mon) for _ in range(p)]
        assert len(idx) == 2, (mon,)
        a, b = sorted(idx)
        out[(a, b)] = out.get((a, b), 0) + coef
    return out

def reduce_terms(terms, n):
    """IBP until all indices <= n.  returns dict (a,b)->c, a<=b<=n"""
    terms = dict(terms)
    while True:
        bad = [k for k in terms if k[1] > n and terms[k] != 0]
        if not bad: return {k: sp.expand(v) for k, v in terms.items() if v != 0}
        a, b = max(bad, key=lambda k: k[1]); c = terms.pop((a, b))
        def add(k, v):
            k = tuple(sorted(k)); terms[k] = terms.get(k, 0) + v
        if b == a: raise RuntimeError('diag index > n')
        if b == a + 1: add((a, a), -sp.Rational(1, 2)*D(c))
        else:
            add((a + 1, b - 1), -c); add((a, b - 1), -D(c))

def form_expr(terms):
    return sum(c*Z[a]*Z[b] for (a, b), c in terms.items())

def level(j, verbose=True):
    n = j + 2
    T = -2*om*zfun(j)*D(ep*Lop(Gexpr), j) if j > 0 else -2*om*u*ep*Lop(Gexpr)
    T = sp.expand(T)
    P = sp.expand(-2*(-1)**j*D(om*zfun(j), j)*ep*Lop(Gexpr))
    ok1 = is_total_derivative(sp.expand(T - P))
    Q = reduce_terms(collect_terms(P), n)
    ok2 = is_total_derivative(sp.expand(P - to_u(form_expr(Q))))
    # transport
    Ttr = sp.expand(-2*om*zfun(j)*D(cc*D(u), j))
    Qtr = reduce_terms(collect_terms(Ttr), j)
    ok3 = is_total_derivative(sp.expand(Ttr - to_u(form_expr(Qtr))))
    if verbose: print('level', j, 'T=P:', ok1, ' P=Q:', ok2, ' transport:', ok3, ' #terms', len(Q), len(Qtr), flush=True)
    return Q, Qtr, (ok1, ok2, ok3)

if __name__ == '__main__':
    res = {}
    for j in (0, 1, 2):
        Q, Qtr, ok = level(j)
        assert all(ok)
        res[j] = (Q, Qtr)
    pickle.dump(res, open('prop1_forms.pkl', 'wb'))
    # top coefficient check
    for j in (0, 1, 2):
        Q = res[j][0]; n = j + 2
        top = sp.simplify(Q[(n, n)].subs(ep, sp.exp(-4*s)))
        print('level', j, 'top diag coeff:', top)
