#!/usr/bin/env python3
"""CFLUX add-on: zeta-commutators of the 4th-order part on the E_I outer ramp x in [57,60] (zeta: 1 -> 0, quintic smoothstep).
v = zeta g,  g = w' (deviation slope, uncut).  v_tau = zeta*g_tau => v_tau = (L1 v) + H,  H = D^2 X(zeta g) - zeta D^2 X g,  X f = A2 f'' + A1 f' + A0 f.
 K0 = int W zeta g H,   K2 = int W (zeta g)'' H''          (the commutator contributions to (1/2)d/dtau V0^2 and V2^2).
Step 1 (sympy): reduce K0, K2 by integration by parts to bilinear forms in g,..,g^(3) with ramp-supported coefficients; machine check (Euler operator).
Step 2 (arb): certified sup over the ramp of each coefficient / W, with W = e^{-psi}, psi' = phi' on x>=38.25, phi=(1+x^2)^{2/3}/4,
          A-derivative bounds from a jet box (default M0..M4x) via Bell polynomials.
Step 3: absolute bounds using N_j^2 = int_ramp W (g^(j))^2 <= int_ramp W dx * (Mg_j / y^j)^2 and the outflow term."""
import sympy as sp, itertools, json, math, sys
from sympy.calculus.euler import euler_equations
from sympy.functions.combinatorial.numbers import bell as bellpoly
from flint import arb, ctx
ctx.prec = 100
x = sp.Symbol('x')
g = sp.Function('g')(x); zt = sp.Function('zt')(x); w = sp.Function('w')(x)
A = {i: sp.Function('A%d' % i)(x) for i in (0, 1, 2)}
G = sp.symbols('g0:12')
D = lambda e, k=1: sp.diff(e, x, k)
def Xop(f): return A[2]*D(f, 2) + A[1]*D(f) + A[0]*f
H = sp.expand(D(Xop(zt*g), 2) - zt*D(Xop(g), 2))
def to_sym(e):
    e = sp.expand(e)
    for k in range(11, 0, -1): e = e.subs(sp.Derivative(g, (x, k)), G[k])
    return e.subs(g, G[0])
def to_fun(e):
    for k in range(11, -1, -1): e = e.subs(G[k], D(g, k) if k else g)
    return e
def terms_of(e):
    P = sp.Poly(to_sym(e), *G); out = {}
    for mon, c in P.terms():
        idx = [k for k, p in enumerate(mon) for _ in range(p)]; assert len(idx) == 2
        a, b = sorted(idx); out[(a, b)] = out.get((a, b), 0) + c
    return out
def reduce(terms, n):
    terms = dict(terms)
    while True:
        bad = [k for k in terms if k[1] > n and terms[k] != 0]
        if not bad: return {k: sp.expand(v) for k, v in terms.items() if v != 0}
        a, b = max(bad, key=lambda k: k[1]); c = terms.pop((a, b))
        def add(k, v): k = tuple(sorted(k)); terms[k] = terms.get(k, 0) + v
        if b == a + 1: add((a, a), -sp.Rational(1, 2)*D(c))
        else: add((a+1, b-1), -c); add((a, b-1), -D(c))
def is_td(e):
    e = sp.expand(e)
    if e == 0: return True
    eqs = euler_equations(e, [g], x)
    return (not eqs) or sp.simplify(sp.expand(eqs[0].lhs)) == 0
def build():
    out = {}
    for name, integrand, n in (('K0', w*zt*g*H, 2), ('K2', w*D(zt*g, 2)*D(H, 2), 3)):
        integrand = sp.expand(integrand)
        Q = reduce(terms_of(integrand), n)
        chk = is_td(sp.expand(integrand - to_fun(sum(c*G[a]*G[b] for (a, b), c in Q.items()))))
        print(name, 'reduced to indices <=', n, ' #pairs', len(Q), ' total-derivative check:', chk, flush=True)
        assert chk
        out[name] = Q
    return out

# ---------- numerics ----------
phi = (1 + x**2)**sp.Rational(2, 3)/4
Wexpr = sp.exp(-phi)
def wratio(n):   # w^(n)/w as explicit function of x  (psi' = phi' on the ramp)
    return sp.simplify(sp.diff(Wexpr, x, n)/Wexpr)
Sx = lambda t: t**3*(10 - 15*t + 6*t*t)
zeta_expr = 1 - Sx((x - 57)/3)       # quintic smoothstep 1 -> 0 on [57,60]  (ASSUMED shape)
def lam(e): return sp.lambdify(x, e, 'mpmath')

def abs_poly_bound(expr, subs_abs):
    """expr polynomial in symbols; replace symbols by their (nonneg) bounds after taking abs of coefficients"""
    expr = sp.expand(expr); P = sp.Poly(expr, *list(subs_abs.keys())); tot = 0
    for mon, c in P.terms():
        t = abs(c)
        for sym, p in zip(P.gens, mon): t *= subs_abs[sym]**p
        tot += t
    return tot

def A_bounds(M, y, Umax, rmax=4):
    """|A_i^(r)| <= Abar_{i,r}/y^r via Bell polynomials; M = dict M0..M4 (jet box of V in x-scaled units), Umax sup|U|.
    coefficients a(p)=(1+p^2)^-2, b(p)=-3p(1+p^2)^-3."""
    p = sp.Symbol('p')
    a = (1 + p**2)**-2; b = -3*p*(1 + p**2)**-3
    M0 = M['M0']
    def supd(f, k):   # sup |f^(k)| on |p|<=M0 (arb branch&bound)
        from lemmaS_interval import sup_arb
        fk = sp.diff(f, p, k); fl = sp.lambdify(p, fk, modules=[{'__builtins__': {}}, 'math'])
        # evaluate with arb via string substitution
        return fk
    # simple rigorous bounds: derivatives of a,b are rational; evaluate sup via dense interval partition
    def sup_rat(fk):
        N = 400; best = arb(0); f = sp.lambdify(p, fk, modules=[{'Pow': lambda a_, b_: a_**b_}])
        for i in range(N):
            lo = arb(M0)*(2*i - N)/N; hi = arb(M0)*(2*i + 2 - N)/N
            ball = arb((lo + hi)/2, (hi - lo)/2)
            v = f(ball); u = abs(v)
            best = max(best, arb(u.upper()), key=lambda z: float(z.upper()))
        return float(best.upper())
    Fa = [sup_rat(sp.diff(a, p, k)) for k in range(0, rmax + 2)]
    Fb = [sup_rat(sp.diff(b, p, k)) for k in range(0, rmax + 2)]
    Vj = [M['M%d' % i]/y**i for i in range(0, rmax + 2)]
    def faa(F, r):     # |D^r F(V)| <= sum_k F_k B_{r,k}(|V'|,|V''|,...)
        if r == 0: return F[0]
        return sum(F[k]*float(bellpoly(r, k, list(Vj[1:r - k + 2]))) for k in range(1, r + 1))
    U = lambda j: Umax*(1.0246**j if j == 0 else math.prod(1.0246 - i for i in range(j)))/y**j
    res = {}
    for r in range(0, rmax + 1):
        res[(2, r)] = faa(Fa, r)
        # A1 = b(V)(V_x+U_x): Leibniz
        res[(1, r)] = sum(math.comb(r, q)*faa(Fb, q)*(Vj[r - q + 1] + abs(U(r - q + 1))) for q in range(0, r + 1))
        # A0 = Da*(U_xx) + Db*U_x^2 (Da,Db = mean-value quotients; derivatives bounded like a',b' compositions)
        Da = lambda q: faa(Fa[1:], q); Db = lambda q: faa(Fb[1:], q)
        res[(0, r)] = sum(math.comb(r, q)*(Da(q)*abs(U(r - q + 2)) + Db(q)*sum(math.comb(r - q, m)*abs(U(m + 1)*U(r - q - m + 1)) for m in range(0, r - q + 1))) for q in range(0, r + 1))
    return res, Fa, Fb

if __name__ == '__main__':
    sys.path.insert(0, '.')
    forms = build()
    import pickle; pickle.dump(forms, open('cflux_forms.pkl', 'wb'))
