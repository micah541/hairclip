#!/usr/bin/env python3
"""Term-by-term jet bookkeeping.  A_i^(r) = sum over monomials in solution jets Z_j (= z_j, j>=2) with sup-coefficients (V_s -> M1, V -> M0, U_j explicit).
Returns Aterms[(i,r)] = {Zmultiset(tuple): bound}."""
import sympy as sp, math, itertools, sys
sys.path.insert(0, '.')
import lemmaS_prop1_bounds as PB
from sympy.functions.combinatorial.numbers import bell as _bell
def Aterms(M0, M1, Umax=0.05, p=1.0246, rmax=6):
    pp = sp.Symbol('p'); a = (1 + pp**2)**-2; b = -3*pp*(1 + pp**2)**-3
    Fa = [PB.sup_rat(sp.diff(a, pp, k), M0) for k in range(rmax + 3)]; Fb = [PB.sup_rat(sp.diff(b, pp, k), M0) for k in range(rmax + 3)]
    Z = sp.symbols('Z0:10')
    def Vj(j):                      # symbolic bound for V_j
        if j == 1: return sp.Float(M1)
        return sp.Float(p**j*Umax) + Z[j]
    def Vb(j): return sp.Float(M1) if j == 1 else sp.Float(p**j*Umax) + Z[j]
    def faa(F, r):
        if r == 0: return sp.Float(F[0])
        return sum(sp.Float(F[k])*_bell(r, k, [Vj(j) for j in range(1, r - k + 2)]) for k in range(1, r + 1))
    def collect(expr):
        expr = sp.expand(expr); out = {}
        for term in sp.Add.make_args(expr):
            c, mon = term.as_independent(*Z); key = []
            for f in sp.Mul.make_args(mon):
                if f == 1: continue
                base, ex = f.as_base_exp(); key += [int(str(base)[1:])]*int(ex)
            key = tuple(sorted(key)); out[key] = out.get(key, 0) + abs(float(c))
        return out
    def merge(*ds):
        o = {}
        for d in ds:
            for k, v in d.items(): o[k] = o.get(k, 0) + v
        return o
    scale = lambda d, s: {k: v*s for k, v in d.items()}
    At = {}
    for r in range(rmax + 1):
        At[(2, r)] = collect(faa(Fa, r))
        t1 = merge(*[scale(collect(faa(Fb, q)*(Vj(r - q + 1) + sp.Float(p**(r - q + 1)*Umax))), math.comb(r, q)) for q in range(r + 1)])
        At[(1, r)] = merge(At[(2, r)], t1)
        t0 = merge(*[merge(scale(collect(faa(Fa[1:], q)), math.comb(r, q)*(p**2 - p)*p**(r - q)*Umax), scale(collect(faa(Fb[1:], q)), math.comb(r, q)*(2*p)**(r - q)*p**2*Umax**2)) for q in range(r + 1)])
        At[(0, r)] = t0
    return At
if __name__ == '__main__':
    import lemmaS_F1 as F
    T = F.build_tables(); At = Aterms(0.11, 0.75)
    for j in (2,):
        gam = 1.0; W = {}
        for (a, b), mons in T[j].items():
            if (a, b) == (j + 2, j + 2): continue
            for (m, i, r), co in mons.items():
                for zk, v in At[(i, r)].items():
                    key = (a, b, zk); W[key] = W.get(key, 0) + abs(float(co))*gam**m*v
        H = {}
        for (a, b, zk), v in W.items():
            Hset = tuple(sorted([x for x in (a, b) if x >= 2] + list(zk))); H.setdefault((Hset, tuple(x for x in (a, b) if x < 2)), 0); H[(Hset, tuple(x for x in (a, b) if x < 2))] += v
        for k, v in sorted(H.items(), key=lambda kv: -kv[1])[:40]: print(k, '%.3g' % v)

# ---------------- L4 / interpolation bounds for the extra jet factors ----------------
KAP_DEFAULT = 5.0
def classify(T, j, M0, M1, gam, M0z=0.02, Umax=0.05, At=None):
    """items: ('P',x,y,w) = w*N_x*N_y ; ('L4',r,w) = w*int|z_2|^2|z_r| ; ('S33',w) = w*int|z_2| z_3^2 (sup z_2 assumed).  Weights include sups of z_0 (M0z), z_1 (M1)."""
    At = At or Aterms(M0, M1, Umax); n = j + 2; sup = {0: M0z, 1: M1}; items = []
    for (a, b), mons in T[j].items():
        if (a, b) == (n, n): continue
        W = {}
        for (m, i, r), co in mons.items():
            for zk, v in At[(i, r)].items():
                if zk: W[zk] = W.get(zk, 0) + abs(float(co))*gam**m*v
        for X, w in W.items():
            lows = [x for x in (a, b) if x < 2]; R = sorted([x for x in (a, b) if x >= 2] + list(X)); Ml = 1.0
            for l in lows: Ml *= sup[l]
            if len(R) == 1:
                if len(lows) != 2: raise ValueError(('linear term', (a, b), X, w))
                lo = sorted(lows, key=lambda l: sup[l]); items.append(('P', lo[1], R[0], w*sup[lo[0]]))
            elif len(R) == 2: items.append(('P', R[0], R[1], w*Ml))
            elif len(R) == 3 and R.count(2) >= 2:
                r = list(R); r.remove(2); r.remove(2); items.append(('L4', r[0], w*Ml))
            elif len(R) == 3 and R.count(3) >= 2 and 2 in R: items.append(('S33', w*Ml))
            else: raise ValueError(('unhandled', (a, b), X, R, w))
    return items
def l4_extra_YC(T, j, M0, M1, gam, kap=KAP_DEFAULT, M0z=0.02, M2star=3.0, Umax=0.05, fallback=None, s33='assumed', dfrac=0.1):
    """core aggregation of classify() into additive Y_k weights (quadratic Young); see classify for the cases"""
    n = j + 2; amin = (1 + M0**2)**-2
    add = {k: 0.0 for k in range(n + 1)}; log = []; dissip = 0.0
    def put(x, y, wt, tag):
        nonlocal dissip
        if x > y: x, y = y, x
        if y > n: raise ValueError(('index beyond dissipation', x, y, tag))
        if y == n and x == n: dissip += wt
        elif y == n:
            d = dfrac*amin; add[x] += wt**2/(4*d); dissip += d
        elif x == y: add[x] += wt
        else: add[x] += wt/2; add[y] += wt/2
        log.append((tag, (x, y), wt))
    try:
        items = classify(T, j, M0, M1, gam, M0z, Umax)
    except ValueError:
        if fallback is None: raise
        return _old_fallback(T, j, M0, M1, gam, kap, M0z, M2star, Umax, fallback)
    for it in items:
        if it[0] == 'P': put(it[1], it[2], it[3], it)
        elif it[0] == 'L4': put(3, it[1], 3*M1*it[2], it); put(2, it[1], kap*M1*it[2], it)
        elif s33 == 'ibp' and n >= 4:
            # int Th|z2| z3^2 = int Th z3 phi', phi = z2|z2|/2 ; IBP: <= 1/2 int (kap Th|z3| + Th|z4|) z2^2 <= 1/2 ||z2^2||_Th (kap N3 + N4),  ||z2^2||_Th <= M1(3N3 + kap N2)
            w = it[1]; put(3, 3, 1.5*kap*M1*w, it); put(3, 4, 1.5*M1*w, it); put(2, 3, 0.5*kap**2*M1*w, it); put(2, 4, 0.5*kap*M1*w, it)
        else: put(3, 3, it[1]*M2star, it)
    return add, log, dissip, 0.5*amin

def _old_fallback_impl(T, j, M0, M1, gam, kap=KAP_DEFAULT, M0z=0.02, M2star=3.0, Umax=0.05, fallback=None):
    """returns (YCadd dict, log).  Pair (a,b) with extra jets X (z_k, k>=2): sup for z_0 (M0z), z_1 (M1); |R|=2: L2 x L2;
    |R|=3 with two z_2: ||z_2^2||_Theta <= M1(3 N3 + kap N2) (exact identity int Th z2^4 = -3 int Th z1 z2^2 z3 - int Th' z1 z2^3, |Th'|=kap Th);
    {2,3,3}: sup route with M2star (flagged)."""
    At = Aterms(M0, M1, Umax); n = j + 2; amin = (1 + M0**2)**-2
    add = {k: 0.0 for k in range(n + 1)}; log = []; dissip = 0.0; sup = {0: M0z, 1: M1}
    def put(x, y, wt, tag):
        nonlocal dissip
        if x > y: x, y = y, x
        if y > n: raise ValueError(('index beyond dissipation', x, y, tag))
        if y == n and x == n: dissip += wt
        elif y == n:
            d = 0.1*amin; add[x] += wt**2/(4*d); dissip += d
        elif x == y: add[x] += wt
        else: add[x] += wt/2; add[y] += wt/2
        log.append((tag, (x, y), wt))
    for (a, b), mons in T[j].items():
        if (a, b) == (n, n): continue
        W = {}
        for (m, i, r), co in mons.items():
            for zk, v in At[(i, r)].items():
                if zk: W[zk] = W.get(zk, 0) + abs(float(co))*gam**m*v
        for X, w in W.items():
            lows = [x for x in (a, b) if x < 2]; R = sorted([x for x in (a, b) if x >= 2] + list(X)); Ml = 1.0
            for l in lows: Ml *= sup[l]
            if len(R) == 1:
                if len(lows) != 2: raise ValueError(('linear term', (a, b), X, w))
                lo = sorted(lows, key=lambda l: sup[l]); put(lo[1], R[0], w*sup[lo[0]], ((a, b), X, 'one low sup')); continue
            if len(R) == 2: put(R[0], R[1], w*Ml, ((a, b), X))
            elif len(R) == 3 and R.count(2) >= 2:
                r = [x for x in R]; r.remove(2); r.remove(2); r = r[0]
                put(3, r, 3*M1*w*Ml, ((a, b), X, 'L4a')); put(2, r, kap*M1*w*Ml, ((a, b), X, 'L4b'))
            elif len(R) == 3 and R.count(3) >= 2 and 2 in R: put(3, 3, w*Ml*M2star, ((a, b), X, 'sup z2 (assumed)'))
            elif fallback is not None and len(R) >= 3:       # assumed sup boxes for the lowest-order high factors (flagged)
                rest = R[-2:]; ww = w*Ml
                for k in R[:-2]: ww *= fallback[k]
                put(rest[0], rest[1], ww, ((a, b), X, 'fallback sup', tuple(R[:-2])))
            else: raise ValueError(('unhandled', (a, b), X, R, w))
    return add, log, dissip, 0.5*amin

def _old_fallback(T, j, M0, M1, gam, kap, M0z, M2star, Umax, fallback):
    return _old_fallback_impl(T, j, M0, M1, gam, kap, M0z, M2star, Umax, fallback)
def ghat_from_YC(YC, j, amin, kap=KAP_DEFAULT):
    Bc = 1.5*amin; G = {k: YC[k] for k in range(j + 1)}
    G[j] += YC[j + 1]*(YC[j + 1]/(4*Bc) + kap**2/2); return G
