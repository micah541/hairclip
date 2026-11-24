"""Lemma U, Phase 1, task P2: certified (arb) library for the inner-ramp constants of the Lemma-U energy.

Weight   Theta = W(sigma) D(s),  W = 1/(e^{-mu(sigma+a)} + 1 + e^{nu(sigma - ln1.5 - a)}),  D = 1 + g e^{-4s},  sigma = s - s_rho,  s = ln y.
Ramp     y in [28, 56] (rescaled), eta = S(t), S the quintic smoothstep, t = (s - ln28)/ln2;  centres s_rho >= ln 56, so sigma <= s - ln56 <= 0 on the ramp.
Params   (mu, nu) = (4, 4), a in [0, 0.3] (interval), g in [G_LO, G_HI] (interval).

Every number returned is an arb upper endpoint of an enclosure over a cover of the domain (branch and bound; the refinement heuristic only chooses cells).
Exact inputs: table coefficients (sympy Rationals from lemmaS/rigor/tables_exact.pkl, read-only), Faa di Bruno/Bell polynomials (exact), box values as Rationals.
Ports (with the weight factor gam^m replaced by an honest Gamma_m = sup|Theta^(m)/Theta| and kappa by kappa~) of lemmaS/rigor/weights_exact.py:
  build_At_sym, eval_At, classify_x, l4_extra_YC_x, budget_ok, ghat_x, ramp (foot block).  Nothing under writeup/lemmaS* is imported or modified."""
import sys, os, math, pickle, heapq
sys.dont_write_bytecode = True
import sympy as sp
from sympy.functions.combinatorial.numbers import bell as bellpoly
from flint import arb, ctx
ctx.prec = 128
Q_ = sp.Rational
TABLES = '/home/claude/s1spec/writeup/lemmaS/rigor/tables_exact.pkl'
PQ = Q_(10246, 10000)
MU, NU = 4, 4
A_LO, A_HI = Q_(0), Q_(3, 10)
G_LO, G_HI = Q_(300000), Q_(600000)

def A(q):
    q = Q_(q); return arb(int(q.p))/arb(int(q.q))
def ipow(x, k):
    """x^k for integer k >= 0 by multiplication (flint 0.9: arb(0, r)**1 returns nan)"""
    r = arb(1)
    for _ in range(k): r = r*x
    return r
def hull(lo, hi):
    return arb(lo).union(arb(hi))
LN2 = arb(2).log(); LN15 = (arb(3)/2).log(); LN28 = arb(28).log(); LN56 = arb(56).log()

def load_tables():
    T = pickle.load(open(TABLES, 'rb'))
    for j in T:
        for mons in T[j].values():
            assert all(co >= 0 for co in mons.values())          # |coef| tables
    return T

# ------------------------------------------------------------------ A-coefficient bounds (port of weights_exact, U and p as parameters)
_pp = sp.Symbol('p'); _a = (1 + _pp**2)**-2; _b = -3*_pp*(1 + _pp**2)**-3
_lam = {}; _supc = {}
def _fun(kind, k):
    if (kind, k) not in _lam: _lam[(kind, k)] = sp.lambdify(_pp, sp.diff(_a if kind == 'a' else _b, _pp, k), modules=[{}])
    return _lam[(kind, k)]
def sup_rat_arb(kind, k, M0q, N=600):
    key = (kind, k, M0q)
    if key in _supc: return _supc[key]
    f = _fun(kind, k); M0 = A(M0q); best = arb(0)
    for i in range(N):
        lo = M0*(2*i - N)/N; hi = M0*(2*i + 2 - N)/N
        u = abs(f(lo.union(hi))).upper()
        if u > best: best = u
    _supc[key] = best; return best

def build_At_sym(UQ, rmax=6):
    """symbolic Aterms with U = UQ (Rational): keys (i,r) -> {Z-multiset: [coefficient expressions in fa*, fb*, m1]}"""
    FA = sp.symbols('fa0:%d' % (rmax + 3)); FB = sp.symbols('fb0:%d' % (rmax + 3)); m1 = sp.Symbol('m1'); Z = sp.symbols('Z0:10')
    def Vj(j): return m1 if j == 1 else PQ**j*UQ + Z[j]
    def faa(F_, r):
        if r == 0: return F_[0]
        return sum(F_[k]*bellpoly(r, k, [Vj(j) for j in range(1, r - k + 2)]) for k in range(1, r + 1))
    def collect(expr, scale=1):
        expr = sp.expand(expr*scale); out = {}
        for term in sp.Add.make_args(expr):
            c, mon = term.as_independent(*Z); key = []
            for f in sp.Mul.make_args(mon):
                if f == 1: continue
                base, ex = f.as_base_exp(); key += [int(str(base)[1:])]*int(ex)
            out.setdefault(tuple(sorted(key)), []).append(c)
        return out
    def merge(*ds):
        o = {}
        for d in ds:
            for k, v in d.items(): o.setdefault(k, []).extend(v)
        return o
    At = {}
    for r in range(rmax + 1):
        At[(2, r)] = collect(faa(FA, r))
        t1 = merge(*[collect(faa(FB, q)*(Vj(r - q + 1) + PQ**(r - q + 1)*UQ), sp.binomial(r, q)) for q in range(r + 1)])
        At[(1, r)] = merge(At[(2, r)], t1)
        At[(0, r)] = merge(*[merge(collect(faa(FA[1:], q), sp.binomial(r, q)*(PQ**2 - PQ)*PQ**(r - q)*UQ), collect(faa(FB[1:], q), sp.binomial(r, q)*(2*PQ)**(r - q)*PQ**2*UQ**2)) for q in range(r + 1)])
    return At, FA, FB, m1

def eval_At(Asym, M0q, M1q, rmax=6):
    At, FA, FB, m1 = Asym
    vals = {m1: A(M1q)}
    for k in range(rmax + 3): vals[FA[k]] = sup_rat_arb('a', k, M0q); vals[FB[k]] = sup_rat_arb('b', k, M0q)
    out = {}
    for key, d in At.items():
        o = {}
        for zk, lst in d.items():
            tot = arb(0)
            for c in lst:
                coef, rest = c.as_coeff_Mul(); v = A(abs(coef))
                for f in sp.Mul.make_args(rest):
                    if f == 1: continue
                    base, ex = f.as_base_exp(); v = v*vals[base]**int(ex)
                tot += v
            o[zk] = tot
        out[key] = o
    return out
def ab_from_At(At):
    """Z-free part (honest: includes the p^j U parts of the jets j >= 2)"""
    return {k: d.get((), arb(0)) for k, d in At.items()}
def ab_exact(M0q, M1q, UQ, rmax=4):
    """port of weights_exact.ab_exact (box with M2.. = 0, U = UQ), used only for the comparison with feas_P / (A)"""
    V = [Q_(M0q), Q_(M1q)] + [Q_(0)]*6
    Fa = [sup_rat_arb('a', k, M0q) for k in range(rmax + 3)]; Fb = [sup_rat_arb('b', k, M0q) for k in range(rmax + 3)]
    def faa_num(F_, r):
        if r == 0: return F_[0]
        return sum((F_[k]*A(bellpoly(r, k, list(V[1:r - k + 2]))) for k in range(1, r + 1)), arb(0))
    pp, Um = PQ, Q_(UQ); Ab = {}
    for r in range(rmax + 1):
        Ab[(2, r)] = faa_num(Fa, r)
        Ab[(1, r)] = Ab[(2, r)] + sum((A(math.comb(r, q))*faa_num(Fb, q)*A(V[r - q + 1] + pp**(r - q + 1)*Um) for q in range(r + 1)), arb(0))
        Ab[(0, r)] = sum((A(math.comb(r, q))*(faa_num(Fa[1:], q)*A((pp**2 - pp)*pp**(r - q)*Um) + faa_num(Fb[1:], q)*A((2*pp)**(r - q)*pp**2*Um**2)) for q in range(r + 1)), arb(0))
    return Ab

# ------------------------------------------------------------------ core YC and L4 classes with a general weight (Gm[m] >= sup|Theta^(m)/Theta|)
def YC_gen(tab, j, Ab, Gm):
    n = j + 2; YC = {k: arb(0) for k in range(n + 1)}
    for (a, b), mons in tab.items():
        if (a, b) == (n, n): continue
        gC = sum((A(co)*Gm[m]*Ab[(i, r)] for (m, i, r), co in mons.items()), arb(0))
        if a == b: YC[a] += gC
        else: YC[a] += gC/2; YC[b] += gC/2
    return YC

def classify_gen(T, j, M0zq, M1q, At, Gm):
    n = j + 2; sup = {0: A(M0zq), 1: A(M1q)}; supq = {0: Q_(M0zq), 1: Q_(M1q)}; items = []
    for (a, b), mons in T[j].items():
        if (a, b) == (n, n): continue
        W = {}
        for (m, i, r), co in mons.items():
            for zk, v in At[(i, r)].items():
                if zk: W[zk] = W.get(zk, arb(0)) + A(co)*Gm[m]*v
        for X, w in W.items():
            lows = [x for x in (a, b) if x < 2]; R = sorted([x for x in (a, b) if x >= 2] + list(X)); Ml = arb(1)
            for l in lows: Ml = Ml*sup[l]
            if len(R) == 1:
                if len(lows) != 2: raise ValueError(('linear term', (a, b), X))
                lo = sorted(lows, key=lambda l: supq[l]); items.append(('P', lo[1], R[0], w*sup[lo[0]]))
            elif len(R) == 2: items.append(('P', R[0], R[1], w*Ml))
            elif len(R) == 3 and R.count(2) >= 2:
                r = list(R); r.remove(2); r.remove(2); items.append(('L4', r[0], w*Ml))
            elif len(R) == 3 and R.count(3) >= 2 and 2 in R: items.append(('S33', w*Ml))
            else: raise ValueError(('unhandled', (a, b), X, R))
    return items

def l4_extra_gen(T, j, M0q, M0zq, M1q, At, Gm, kap, dfrac):
    """port of l4_extra_YC_x with s33 = 'ibp', weights |coef| Gamma_m, kappa -> kappa~ (arb)"""
    n = j + 2; amin = (1 + A(M0q)**2)**-2
    add = {k: arb(0) for k in range(n + 1)}; dissip = [arb(0)]; ndp = [0]; M1 = A(M1q)
    def put(x, y, wt):
        if x > y: x, y = y, x
        if y > n: raise ValueError(('index beyond dissipation', x, y))
        if y == n and x == n: dissip[0] += wt
        elif y == n:
            d = A(dfrac)*amin; add[x] += wt**2/(4*d); ndp[0] += 1
        elif x == y: add[x] += wt
        else: add[x] += wt/2; add[y] += wt/2
    items = classify_gen(T, j, M0zq, M1q, At, Gm)
    for it in items:
        if it[0] == 'P': put(it[1], it[2], it[3])
        elif it[0] == 'L4': put(3, it[1], 3*M1*it[2]); put(2, it[1], kap*M1*it[2])
        elif n >= 4:
            w = it[1]; put(3, 3, kap*3/2*M1*w); put(3, 4, M1*w*3/2); put(2, 3, kap**2/2*M1*w); put(2, 4, kap/2*M1*w)
        else: raise ValueError('S33 below level 2')
    return add, dissip[0], ndp[0], items

def budget(nn, ndp, dfrac, M0q, extra=Q_(0)):
    """dissipation used = nn + (ndp*dfrac + extra)*amin <= amin/2 ; returns (ok, slack/amin lower bound)"""
    amin = (1 + A(M0q)**2)**-2; fr = ndp*Q_(dfrac) + extra
    sl = (A(Q_(1, 2) - fr)) - nn/amin
    return bool(sl.lower() >= 0), float(sl.lower())

def ghat_gen(YC, j, M0q, kap2):
    amin = (1 + A(M0q)**2)**-2; Bc = A(Q_(3, 2))*amin; G = {k: YC[k] for k in range(j + 1)}
    G[j] = G[j] + YC[j + 1]*(YC[j + 1]/(4*Bc) + kap2/2); return G

# ------------------------------------------------------------------ weight W: R_k = W^(k)/W as polynomials in p, q (sympy, exact integer coefficients)
_p, _q = sp.symbols('p q')
def R_polys(mu, nu, kmax=4):
    R = {1: mu*_p - nu*_q}; dp = _p*(-mu + mu*_p - nu*_q); dq = _q*(nu + mu*_p - nu*_q)
    for k in range(1, kmax):
        R[k + 1] = sp.expand(sp.diff(R[k], _p)*dp + sp.diff(R[k], _q)*dq + R[1]*R[k])
    return {k: sp.Poly(v, _p, _q) for k, v in R.items()}
RP = R_polys(MU, NU)
class TP2:
    """centred (Taylor) form of an integer polynomial P(p, q): P(p0 + dp, q0 + dq) = sum c_ij dp^i dq^j, c_ij = D^(i,j)P(p0,q0)/(i! j!), exact expansion"""
    def __init__(self, P):
        e = P.as_expr(); deg = P.total_degree(); self.terms = []
        for i in range(deg + 1):
            for j in range(deg + 1 - i):
                d = sp.Poly(sp.diff(e, _p, i, _q, j)/(sp.factorial(i)*sp.factorial(j)), _p, _q)
                if not d.is_zero: self.terms.append((i, j, [(a, b, Q_(c)) for (a, b), c in d.terms()]))
    def __call__(self, p, q):
        p0 = arb(p.mid()); q0 = arb(q.mid()); dp = p - p0; dq = q - q0; tot = arb(0)
        for i, j, cs in self.terms:
            c = sum((A(cc)*ipow(p0, a)*ipow(q0, b) for a, b, cc in cs), arb(0)); tot += c*ipow(dp, i)*ipow(dq, j)
        return tot
RT2 = {k: TP2(v) for k, v in RP.items()}
def Rk(p, q, kmax=3):
    return [arb(1)] + [RT2[k](p, q) for k in range(1, kmax + 1)]

C15 = LN15
def pq_from_sigma(sig, a):
    """sig, a: arb balls.  p = P/Den, q = Q/Den, Den = P + 1 + Q, P = e^{-mu(sig+a)}, Q = e^{nu(sig-c-a)} (single-occurrence forms)"""
    p = 1/(1 + (MU*(sig + a)).exp() + ((MU + NU)*sig - NU*C15 + (MU - NU)*a).exp())
    q = 1/(1 + (-NU*(sig - C15 - a)).exp() + (-(MU + NU)*sig + NU*C15 - (MU - NU)*a).exp())
    return p, q
def pq_lower_tail(sigT, a):
    """sigma in (-inf, sigT]: p in [p_lo, 1], q in [0, Q(sigT)]"""
    X = (MU*(sigT + a)).exp() + ((MU + NU)*sigT - NU*C15 + (MU - NU)*a).exp()
    plo = (1/(1 + X)).lower(); qhi = ((NU*(sigT - C15 - a)).exp()).upper()
    return hull(plo, 1), hull(0, qhi)
def pq_upper_tail(sigT, a):
    """sigma in [sigT, inf): q in [q_lo, 1], p in [0, P/Q(sigT)]"""
    X = (-NU*(sigT - C15 - a)).exp() + (-(MU + NU)*sigT + NU*C15 - (MU - NU)*a).exp()
    qlo = (1/(1 + X)).lower(); phi = ((-(MU + NU)*sigT + NU*C15 - (MU - NU)*a).exp()).upper()
    return hull(0, phi), hull(qlo, 1)

def theta_s(s, g):
    """theta = g e^{-4s}/(1 + g e^{-4s}) = 1/(1 + e^{4s}/g)  (single occurrence of s and g)"""
    return 1/(1 + (4*s).exp()/g)

# ------------------------------------------------------------------ cutoff eta = S(t): derivatives in s of eta and eta^2 (exact integer polynomials in t)
_t = sp.Symbol('t'); _S = _t**3*(10 - 15*_t + 6*_t**2)
E2P = [sp.Poly(sp.diff(_S**2, _t, l), _t) for l in range(5)]
EP = [sp.Poly(sp.diff(_S, _t, l), _t) for l in range(4)]
class TP1:
    """centred (Taylor) form of an integer polynomial in t"""
    def __init__(self, P):
        e = P.as_expr(); deg = P.degree(); self.terms = []
        for k in range(deg + 1):
            d = sp.Poly(sp.diff(e, _t, k)/sp.factorial(k), _t)
            if not d.is_zero: self.terms.append((k, [(a, Q_(c)) for (a,), c in d.terms()]))
    def __call__(self, t):
        t0 = arb(t.mid()); h = t - t0; tot = arb(0)
        for k, cs in self.terms:
            c = sum((A(cc)*ipow(t0, a) for a, cc in cs), arb(0)); tot += c*ipow(h, k)
        return tot
_TPC = {}
def _tp(P, t):
    key = id(P)
    if key not in _TPC: _TPC[key] = TP1(P)
    return _TPC[key](t)
def e2(t, l): return _tp(E2P[l], t)/ipow(LN2, l)
def eta(t, l): return _tp(EP[l], t)/ipow(LN2, l)

# ------------------------------------------------------------------ branch and bound (certified sup of |f| over a cover of a box domain)
def bb_sup(f, boxes, tol=0.002, maxeval=60000, mindiam=None, smart=True):
    """f(box) -> arb (enclosure of the target on the box); f(box, point=True) -> arb at the box centre (lower-bound heuristic).
    boxes: list of tuples of (lo, hi) float pairs; a component may carry a tail flag handled inside f.  Returns (upper bound, best lower, nevals, nboxes left).
    The result is the max of arb upper endpoints over a cover of the union of the boxes: certified."""
    def split_smart(b):
        """halve the dimension whose collapse to its midpoint lowers the enclosure most (costs one evaluation per finite dimension;
        only chooses the cell split, every bound is still an enclosure)"""
        best = None
        for i, (lo, hi) in enumerate(b):
            if isinstance(hi, float) and math.isinf(hi): continue
            m = (lo + hi)/2
            if not (lo < m < hi): continue
            c = list(b); c[i] = (m, m); u = ub(tuple(c))
            if best is None or u < best[0]: best = (u, i)
        if best is None: return None, 0
        k = best[1]; lo, hi = b[k]; m = (lo + hi)/2
        b1 = list(b); b2 = list(b); b1[k] = (lo, m); b2[k] = (m, hi); return (tuple(b1), tuple(b2)), sum(1 for (lo, hi) in b if not (isinstance(hi, float) and math.isinf(hi)))
    def split(b):
        widths = [(hi - lo) if not (isinstance(hi, float) and math.isinf(hi)) else -1 for (lo, hi) in b]
        # scaled widths: f may declare scales via f.scales
        sc = getattr(f, 'scales', [1.0]*len(b))
        k = max(range(len(b)), key=lambda i: widths[i]*sc[i])
        lo, hi = b[k]; m = (lo + hi)/2
        if not (lo < m < hi): return None
        b1 = list(b); b2 = list(b); b1[k] = (lo, m); b2[k] = (m, hi); return tuple(b1), tuple(b2)
    def ub(b):
        v = abs(f(b)).upper(); u = float(v)
        if u != u: u = math.inf
        if arb(u) < v: u = math.nextafter(u, math.inf)
        return u
    def lb(b):
        try:
            v = f(b, point=True)
        except TypeError:
            return 0.0
        x = abs(v).lower(); return float(x) if float(x) == float(x) else 0.0
    heap = []; best_lb = 0.0; nev = 0; done_max = 0.0
    for b in boxes:
        u = ub(b); nev += 1; best_lb = max(best_lb, lb(b)); heapq.heappush(heap, (-u, b))
    while heap:
        u, b = heapq.heappop(heap); u = -u
        if u <= max(best_lb*(1 + tol), done_max) or nev >= maxeval:
            return max(u, done_max), best_lb, nev, len(heap)
        if smart:
            ch, extra = split_smart(b); nev += extra
        else: ch = split(b)
        if ch is None: done_max = max(done_max, u); continue
        for c in ch:
            uc = ub(c); nev += 1; best_lb = max(best_lb, lb(c)); heapq.heappush(heap, (-uc, c))
    return done_max, best_lb, nev, 0
