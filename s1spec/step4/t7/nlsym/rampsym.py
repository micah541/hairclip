"""Symbolic structure of the intrinsic nonlinear remainder N (T8 n3sym form) and of its x-derivative, as polynomials in the
jets e'', e''', e'''', e^(5) with coefficients in (e, e', k..k^(4), r', r'').  Output: nlsym/rampforms.json (sympy srepr of
each coefficient) for nlterms7.py.  usage: python3 -I nlsym/rampsym.py"""
import os, json, time
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
E = sp.symbols('e0:6'); K = sp.symbols('k0:6'); r1, r2, r3 = sp.symbols('r1 r2 r3')


def D(f):
    return (sum(sp.diff(f, E[i]) * E[i + 1] for i in range(5)) + sum(sp.diff(f, K[i]) * K[i + 1] for i in range(5))
            + sp.diff(f, r1) * r2 + sp.diff(f, r2) * r3)


def operator():
    """T7 nl7 operator with the base identities and the drift-free representative (c = b = p = q = n = t1 = x = Q = A = 0;
    T8 n3sym checks that N does not depend on them): w = -k k', theta_tau = -k''' - k^2 k' + r', V = -k'' + r"""
    al = 1 - K[0] * E[0]; J = sp.sqrt(al ** 2 + E[1] ** 2)
    kap = (al * E[2] - E[1] * D(al) + K[0] * (al ** 2 + E[1] ** 2)) / J ** 3
    kss = D(D(kap) / J) / J
    w = -K[0] * K[1]; tht = -K[3] - K[0] ** 2 * K[1] + r1; V = -K[2]
    return (J / al) * (-kss) - V + E[1] * (w - E[0] * tht) / al


def main():
    t0 = time.time()
    F = operator()
    zero = {v: 0 for v in E[:5]}
    N = F - F.subs(zero) - sum(sp.diff(F, v).subs(zero) * v for v in E[:5])
    print('N built', round(time.time() - t0), 's', flush=True)
    jets = (E[2], E[3], E[4])
    def poly(expr, gens):
        num, den = sp.fraction(sp.together(expr))
        P = sp.Poly(sp.expand(num), *gens)
        return {m: sp.simplify(cf / den) if False else cf / den for m, cf in P.terms()}, den
    PN, denN = poly(N, jets)
    print('N monomials (e2,e3,e4):', sorted(PN.keys()), round(time.time() - t0), 's', flush=True)
    DN = D(N)
    PD, denD = poly(DN, (E[2], E[3], E[4], E[5]))
    print('DN monomials (e2,e3,e4,e5):', sorted(PD.keys()), round(time.time() - t0), 's', flush=True)
    out = dict(N={str(m): sp.srepr(c) for m, c in PN.items()}, DN={str(m): sp.srepr(c) for m, c in PD.items()},
               note='coefficients of monomials e2^a e3^b e4^c (N) and e2^a e3^b e4^c e5^d (DN); symbols e0, e1, k0..k5, r1, r2')
    json.dump(out, open(os.path.join(HERE, 'rampforms.json'), 'w'))
    print('done', round(time.time() - t0), 's')


if __name__ == '__main__':
    main()
