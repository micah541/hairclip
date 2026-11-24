"""Lemma F, step 1: the exact level-j energy identity, machine-checked.

v = e', w_m = d^m v (w_{-1} := e), y = x - X0.  Linear far-field equation (T5 coefficients, b1 = -c y + beta):
    v_tau = -v'''' - c y v' + (beta v + b2 v' + c0 e)',
    beta = Q - c X0 - 6 k k',  b2 = -k^2,  c0 = -3 k'^2 - 4 k k'',  Q' = k'^2,  k = kappa~, c = 1/4 - a~.
Claim (pointwise in x, modulo a total x-derivative):
    2 T w_j d^j[v_tau]  =  -2 T w_{j+2}^2 + P w_{j+1}^2 + q w_j^2 + sum_{i=-1}^{j-1} p_{ji} w_j w_i + X_j
with P = 4T'' + 2k^2 T, q, p_{ji} below, and X_j the terms produced by the integrations by parts that keep every coefficient
derivative at kappa^(<=5) (base is C^4 with kappa^(5) in L^inf): j = 3: c0^(4) e-term; j = 4: beta^(5), c0^(4), c0^(5).
The test: the Euler-Lagrange expression (w.r.t. e) of LHS - RHS vanishes identically (T, k, Q arbitrary, Q' = k'^2).
usage: python3 lfsym.py [jmax]   (prints OK per level; ~1-3 min for j = 4)"""
import sys, pickle, time
import sympy as sp
from math import comb

x, c, X0 = sp.symbols('x c X0', real=True)
T = sp.Function('T')(x); k = sp.Function('k')(x); Q = sp.Function('Q')(x)
y = x - X0
NJ = 16
E = sp.symbols('E0:%d' % NJ)                 # jet variables: E[n] = e^(n)

def Dx(f):
    """total x-derivative on (T, k, Q, jets of e), with Q' = k'^2"""
    out = sp.diff(f, x).subs(sp.Derivative(Q, x), sp.diff(k, x) ** 2)
    for n in range(NJ - 1):
        if f.has(E[n]): out += E[n + 1] * sp.diff(f, E[n])
    return out

def Dn(f, n):
    for _ in range(n): f = Dx(f)
    return f

def w(m):           # w_m = e^(m+1); w_{-1} = e
    return E[m + 1]

beta = Q - c * X0 - 6 * k * sp.diff(k, x)
b2 = -k ** 2
c0 = -3 * sp.diff(k, x) ** 2 - 4 * k * sp.diff(k, x, 2)

def coef_d(h, m):  # m-th derivative of a base coefficient (Q' -> k'^2)
    return Dn(h, m)

def lhs(j):
    v = E[1]
    rhs_v = -E[5] - c * y * E[2] + Dx(beta * v + b2 * E[2] + c0 * E[0])     # v_tau
    return 2 * T * w(j) * Dn(rhs_v, j)

def parts(j):
    """RHS pieces as sympy expressions in T, k, Q (coefficient functions), and the list of bilinear pairs"""
    Tp = [Dn(T, m) for m in range(5)]
    P = 4 * Tp[2] + 2 * k ** 2 * T
    q = (-Tp[4] + c * (y * Tp[1] + (1 - 2 * j) * T) - Dx(T * beta) + 2 * (j + 1) * T * coef_d(beta, 1)
         + Dn(T * b2, 2) - (j + 1) * Dx(T * coef_d(b2, 1)) + j * (j + 1) * T * coef_d(b2, 2) + 2 * T * c0)
    pairs = {}          # (a, b) -> coefficient, meaning coefficient * w_a * w_b  (a >= b), a,b in {-1..j+2}
    def add(a, b, cf):
        key = (max(a, b), min(a, b)); pairs[key] = pairs.get(key, 0) + cf
    special = []        # (description, coefficient-order) of terms that are integrated by parts
    for i in range(-1, j):
        # beta: m = j+1-i >= 2, i >= 0
        m = j + 1 - i
        if i >= 0 and m >= 2:
            if m <= 4: add(j, i, 2 * comb(j + 1, m) * T * coef_d(beta, m))
            else:      # beta^(5) (j=4, i=0): -2C int beta^(4) (T' w_j w_i + T w_{j+1} w_i + T w_j w_{i+1})
                cf = 2 * comb(j + 1, m); h4 = coef_d(beta, m - 1)
                add(j, i, -cf * h4 * Tp[1]); add(j + 1, i, -cf * h4 * T); add(j, i + 1, -cf * h4 * T); special.append(('beta', m))
        # b2: m = j+2-i in [3, j+1]
        m = j + 2 - i
        if 3 <= m <= j + 1: add(j, i, 2 * comb(j + 1, m) * T * coef_d(b2, m))
        # c0: m = j - i in [1, j+1]
        m = j - i
        if 1 <= m <= j + 1:
            cf = 2 * comb(j + 1, m)
            if m <= 3: add(j, i, cf * T * coef_d(c0, m))
            elif m == 4 or (m == 5 and False):
                h3 = coef_d(c0, 3)
                add(j, i, -cf * h3 * Tp[1]); add(j + 1, i, -cf * h3 * T); add(j, i + 1, -cf * h3 * T); special.append(('c0', m))
            elif m == 5:   # two IBPs: cf int c0^(3) (T w_j w_i)''
                h3 = coef_d(c0, 3)
                add(j, i, cf * h3 * Tp[2]); add(j + 1, i, 2 * cf * h3 * Tp[1]); add(j, i + 1, 2 * cf * h3 * Tp[1])
                add(j + 2, i, cf * h3 * T); add(j + 1, i + 1, 2 * cf * h3 * T); add(j, i + 2, cf * h3 * T); special.append(('c0', m))
    return dict(P=P, q=q, pairs=pairs, special=special)

def rhs(j, pt):
    out = -2 * T * w(j + 2) ** 2 + pt['P'] * w(j + 1) ** 2 + pt['q'] * w(j) ** 2
    for (a, b), cf in pt['pairs'].items(): out += cf * w(a) * w(b)
    return out

def euler(L, nmax):
    tot = 0
    for n in range(nmax + 1):
        dL = sp.diff(L, E[n])
        if dL != 0: tot += (-1) ** n * Dn(dL, n)
    return sp.expand(tot)

def max_order(expr, f):
    """highest derivative order of function f appearing in expr"""
    mo = 0 if expr.has(f) else -1
    for d in expr.atoms(sp.Derivative):
        if d.expr == f: mo = max(mo, sum(cnt for _, cnt in d.variable_count))
    return mo

if __name__ == '__main__':
    jmax = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    store = {}
    for j in range(jmax + 1):
        t0 = time.time()
        pt = parts(j)
        diff = sp.expand(lhs(j) - rhs(j, pt))
        Eu = euler(diff, j + 8)
        ok = (Eu == 0)
        # regularity: highest kappa derivative in the RHS coefficients
        allc = [pt['P'], pt['q']] + list(pt['pairs'].values())
        kmax = max(max_order(sp.expand(cf), k) for cf in allc); Qmax = max(max_order(sp.expand(cf), Q) for cf in allc)
        print(f"level {j}: Euler(LHS-RHS) == 0: {ok};  pairs {sorted(pt['pairs'].keys())};  IBP'd {pt['special']};  max kappa-derivative {kmax}, Q-derivative {Qmax}  ({time.time()-t0:.0f}s)", flush=True)
        assert ok
        store[j] = pt
    pickle.dump(store, open('lf_forms.pkl', 'wb'))
    print('saved lf_forms.pkl')
