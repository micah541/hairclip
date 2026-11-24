"""Rigorous upper bound of max_{n in simplex} sqrt(sum_f Qf(n)^2), Qf = g00 n0^2+g01 n0n1+g11 n1^2+g02 n0n2+g12 n1n2+g22 n2^2, n0+n1+n2=1, n>=0, g>=0.
Q is increasing in each n_i, so the bound on a cell [lo0,hi0]x[lo1,hi1] (n2 <= 1-lo0-lo1) is the value at the upper corner.  All bounds are evaluated in
arb ball arithmetic (cell endpoints are dyadic, hence exact); the returned value is the upper end of the largest ball over all leaves.  Candidate lower
bounds (only used to prune, never reported) are plain floats."""
from math import sqrt
from flint import arb
def _Ff(forms, n0, n1, n2):
    s = 0.0
    for g in forms:
        v = g['00']*n0*n0 + g['01']*n0*n1 + g['11']*n1*n1 + g['02']*n0*n2 + g['12']*n1*n2 + g['22']*n2*n2; s += v*v
    return sqrt(s)
def _Fa(fa, n0, n1, n2):
    s = arb(0)
    for g in fa:
        v = g['00']*n0*n0 + g['01']*n0*n1 + g['11']*n1*n1 + g['02']*n0*n2 + g['12']*n1*n2 + g['22']*n2*n2; s += v*v
    return s.sqrt()
def simplex_max(forms, tol=2e-3, maxcells=2000000):
    fa = [{k: arb(float(v)) for k, v in g.items()} for g in forms]
    def ub(lo0, hi0, lo1, hi1):
        n2 = arb(1) - arb(lo0) - arb(lo1)
        n2 = arb(1) if n2.upper() > 1 else n2
        return float(_Fa(fa, arb(hi0), arb(hi1), n2).upper())
    best_lo = 0.0; cells = [(0.0, 1.0, 0.0, 1.0)]; n = 0; leaves = []
    while cells and n < maxcells:
        lo0, hi0, lo1, hi1 = cells.pop(); n += 1
        if lo0 + lo1 > 1: continue
        u = ub(lo0, hi0, lo1, hi1)
        for c0 in (lo0, (lo0+hi0)/2, hi0):
            for c1 in (lo1, (lo1+hi1)/2, hi1):
                if c0 + c1 <= 1: best_lo = max(best_lo, _Ff(forms, c0, c1, 1-c0-c1)*(1-1e-12))
        if u <= best_lo*(1+tol): leaves.append(u); continue
        m0, m1 = (lo0+hi0)/2, (lo1+hi1)/2
        cells += [(lo0, m0, lo1, m1), (m0, hi0, lo1, m1), (lo0, m0, m1, hi1), (m0, hi0, m1, hi1)]
    rest = [ub(*c) for c in cells if c[0] + c[2] <= 1]
    assert not cells, "maxcells reached (result still a valid bound)" if False else True
    return max(leaves + rest + [0.0]), best_lo
