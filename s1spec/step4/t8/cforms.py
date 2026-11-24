"""T8 certified: the level-j identity (I_j) of LEMMA_F.md with GENERIC coefficient functions beta, b2, c0 (so that the
residual part A of the transport, beta -> beta + A, and any other exact coefficient can be inserted), machine-checked
by the Euler-operator test (as lfsym.py), and exported as polynomials for interval evaluation.

    v_tau = -v'''' - c y v' + (beta v + b2 v' + c0 e)',     y = x - X0,  w_m = d^m v,  w_{-1} = e
    2 T w_j d^j[v_tau] == -2 T w_{j+2}^2 + P w_{j+1}^2 + q w_j^2 + sum p_ab w_a w_b     (mod d/dx)

Export (forms.json): per level j, polynomials in the symbols
    T0..T4 (T^(m)), B0..B4 (beta^(m)), D0..D5 (b2^(m)), C0..C3 (c0^(m)), c, y
with exact rational coefficients, as lists [coef_num, coef_den, {symbol: power}].
usage: python cforms.py [jmax]      (jmax = 3 default; ~20 s)"""
import os, sys, json, time
import sympy as sp
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
x, c, X0 = sp.symbols('x c X0', real=True)
T = sp.Function('T')(x); be = sp.Function('beta')(x); b2 = sp.Function('b2')(x); c0 = sp.Function('c0')(x)
y = x - X0
NJ = 14
E = sp.symbols('E0:%d' % NJ)


def Dx(f):
    out = sp.diff(f, x)
    for n in range(NJ - 1):
        if f.has(E[n]): out += E[n + 1] * sp.diff(f, E[n])
    return out


def Dn(f, n):
    for _ in range(n): f = Dx(f)
    return f


def w(m): return E[m + 1]


def lhs(j):
    rhs_v = -E[5] - c * y * E[2] + Dx(be * E[1] + b2 * E[2] + c0 * E[0])
    return 2 * T * w(j) * Dn(rhs_v, j)


def parts(j):
    Tp = [sp.diff(T, x, m) for m in range(5)]
    d = lambda h, m: sp.diff(h, x, m)
    P = 4 * Tp[2] + 2 * (-b2) * T                       # 4T'' - 2 b2 T  (= 4T'' + 2 kappa^2 T)
    q = (-Tp[4] + c * (y * Tp[1] + (1 - 2 * j) * T) - sp.diff(T * be, x) + 2 * (j + 1) * T * d(be, 1)
         + sp.diff(T * b2, x, 2) - (j + 1) * sp.diff(T * d(b2, 1), x) + j * (j + 1) * T * d(b2, 2) + 2 * T * c0)
    pairs = {}
    def add(a, b, cf):
        key = (max(a, b), min(a, b)); pairs[key] = pairs.get(key, 0) + cf
    for i in range(-1, j):
        m = j + 1 - i
        if i >= 0 and m >= 2:
            assert m <= 4
            add(j, i, 2 * comb(j + 1, m) * T * d(be, m))
        m = j + 2 - i
        if 3 <= m <= j + 1: add(j, i, 2 * comb(j + 1, m) * T * d(b2, m))
        m = j - i
        if 1 <= m <= j + 1:
            cf = 2 * comb(j + 1, m)
            if m <= 3: add(j, i, cf * T * d(c0, m))
            elif m == 4:      # j = 3, i = -1: 2 T c0'''' w3 e -> -2 c0''' (T' w3 e + T w4 e + T w3 w0)
                h3 = d(c0, 3)
                add(j, i, -cf * h3 * Tp[1]); add(j + 1, i, -cf * h3 * T); add(j, i + 1, -cf * h3 * T)
            else: raise NotImplementedError('j >= 4 not exported (level 4 dropped: LEMMA_F C5, N3)')
    return dict(P=P, q=q, pairs=pairs)


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


SYM = {}
def _symbols():
    if SYM: return SYM
    for m in range(5): SYM[('T', m)] = sp.Symbol('T%d' % m)
    for m in range(5): SYM[('B', m)] = sp.Symbol('B%d' % m)
    for m in range(6): SYM[('D', m)] = sp.Symbol('D%d' % m)
    for m in range(4): SYM[('C', m)] = sp.Symbol('C%d' % m)
    return SYM


def to_poly(expr):
    """replace derivatives by symbols; y stays as symbol 'y' (x - X0)"""
    S = _symbols(); sub = {}
    for (nm, f) in (('T', T), ('B', be), ('D', b2), ('C', c0)):
        for m in range(6, 0, -1):
            if (nm, m) in S: sub[sp.Derivative(f, (x, m))] = S[(nm, m)]
        sub[f] = S[(nm, 0)]
    ys = sp.Symbol('y')
    e = sp.expand(sp.expand(expr).subs(sub).subs(x, ys + X0))
    assert not e.has(x) and not e.has(X0), e
    P = sp.Poly(e, *sorted(e.free_symbols, key=str)) if e.free_symbols else None
    out = []
    if P is None:
        if e != 0: out.append([int(sp.Rational(e).p), int(sp.Rational(e).q), {}])
        return out
    gens = [str(g) for g in P.gens]
    for mon, cf in P.terms():
        cf = sp.Rational(cf)
        out.append([int(cf.p), int(cf.q), {gens[k]: int(p) for k, p in enumerate(mon) if p}])
    return out


def main():
    jmax = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    out = {}
    for j in range(jmax + 1):
        t0 = time.time(); pt = parts(j)
        Eu = euler(sp.expand(lhs(j) - rhs(j, pt)), j + 8)
        ok = (Eu == 0)
        print(f'level {j}: Euler(LHS - RHS) == 0 for generic T, beta, b2, c0: {ok}; pairs {sorted(pt["pairs"])} ({time.time()-t0:.0f}s)', flush=True)
        assert ok
        # negative control: drop one pair -> must fail
        if j >= 1:
            k0 = sorted(pt['pairs'])[0]; bad = dict(pt); bad['pairs'] = {k: v for k, v in pt['pairs'].items() if k != k0}
            assert euler(sp.expand(lhs(j) - rhs(j, bad)), j + 8) != 0
        out[str(j)] = dict(P=to_poly(pt['P']), q=to_poly(pt['q']),
                           pairs={f'{a},{b}': to_poly(cf) for (a, b), cf in pt['pairs'].items()})
    with open(os.path.join(HERE, 'forms.json'), 'w') as f: json.dump(out, f, indent=0)
    print('saved forms.json (negative controls failed as they should)')


if __name__ == '__main__':
    main()
