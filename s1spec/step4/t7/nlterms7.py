"""T7 Lemma 4b (pointwise tip box from Lambda_e) and Lemma 4c (the nonlinear terms of (Y), with Q and the nonlinear and
E_ij parts of (M)), certified in arb.  T7_REPORT 4b, 4c.  Review step4/STEP4_REVIEW.md A2, A3; T8 note on (M) (B6).

Lemma 4b (base-free): for every time, with Lambda = Lambda_e = A^2 + eps B^2 + mu_e (V0^2 + eps_I V2^2),
   |e(x)| <= c_e(x) Lambda^(1/2),  |e'(x)| <= c_e1(x) Lambda^(1/2)  on [0, XE] (XE = 34.249 = end of supp chi_32), with
   (i)   x <= X1 (u = e):  rho u^2 <= Phi1 A^2 + 2 A W1 <= (Phi1 + 2s) Lambda,  s^2 = 1/(2 sqrt eps) + c1;
   (ii)  x <= X1:          rho u'^2 <= Phi1 W1^2 + 2 W1 B <= 2 C_A Lambda,  C_A = (m31 eps^-3/2 + c1/(4 eps))^1/2 + Phi1/2 (1/(2 sqrt eps) + c1);
   (iii) x in [16.25, 34.25] (zeta = 1, v = e'): f = psi_c (1+x^2)^-1 v with psi_c = 1 on [10.25, x0],
         cos^2(pi (x - x0)/(2 l_c)) on [x0, x0 + l_c], 0 beyond (l_c = 100; W = (1+x^2)^-2 there since ell~/2 >= 205),
         sup f^2 <= K_line (|f|^2 + eps_I |f''|^2), K_line = eps_I^(-1/4)/(2 sqrt 2), |f''| <= a V2 + b V0 (V1 <= kV1 V0 + V2/2):
         |e'(x0)| <= (1 + x0^2) (K_line c_om / mu_e)^(1/2) Lambda^(1/2);
   (iv)  x >= 16.25: |e(x)| <= |e(16.25)| + H(x) V0, H^2 = int_16.25^x (1+s^2)^2 ds (W >= (1+x^2)^-2 also on the flattened part).
Lemma 4c: see T7_REPORT 4c.  Output: nl_terms.json (per certificate interval: closure terms for closure_cert.py --nl, the
pieces, the (M) bound; plus tip_box).
usage: python3 -I nlterms7.py [--work run_v2se] [--Ybox 1.7e-12] [--n3 ../p2/t8se/constants/n3.json] [--ids 1-25]
       then closure_cert.py --work W --mu 0.036 --nl W/nl_terms.json --out closure_cert_nl.json; merge.py publishes both."""
import os, sys, json, math, time, argparse
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); S4 = os.path.dirname(HERE)
for p in (HERE, os.path.join(HERE, 'nlsym'), os.path.join(S4, 'p2', 't8se'), os.path.join(S4, 't6')):
    if p not in sys.path: sys.path.insert(0, p)
from flint import arb, ctx
ctx.prec = 128
import nlforms_gen as G
from td2 import TD2
import n3qforms
import reader as T6R

EPS, EPSI, C1, PHI2 = arb('1e-5'), arb('3e-5'), arb('0.5804'), arb('1.3')
X0R, XS = arb('10.25'), arb('16.25')          # zeta ramp [10.25, 16.25]; zeta = 1 for x >= 16.25
XM = arb('30.25')                             # middle region M = [16.25, 30.25]; T8's core box from 30.249 on
ZERO = arb(0)
def up(v): return math.nextafter(float(arb(v).upper()), math.inf)
def U(v): return arb(repr(float(v)))          # stored floats are certified bounds: take them exactly
def ab(v): return abs(arb(v))
def sq(v): return arb(up(v)).sqrt()             # sqrt of an upper bound (monotone): an upper bound
def ub(v): return arb(up(v))                    # exact nonnegative upper bound
def rho(x): return (-(1 + x * x) ** (arb(2) / 3) / 4).exp()
def Wf(x): return 1 / (1 + x * x) ** 2
def hull(vals):
    lo = min(float(arb(v).lower()) for v in vals); hi = max(float(arb(v).upper()) for v in vals)
    return arb((lo + hi) / 2, (hi - lo) / 2) + arb(0, 1e-300)


# ============================================================================================ Lemma 4b
def lemma4b(C, MU, c1W, slabs_ell_min):
    import certdata as CD
    XE, X1 = CD.XE, CD.X1
    PHI1 = XE / (3 * (1 + XE * XE) ** (arb(1) / 3))                 # sup phi' on [0, XE] (phi' increasing)
    s = (1 / (2 * EPS.sqrt()) + C1).sqrt()
    Ke = lambda x: ((PHI1 + 2 * s) / rho(x)).sqrt()
    m31 = 3 * arb(3).sqrt() / 16
    CA = (m31 * EPS ** (-arb(3) / 2) + C1 / (4 * EPS)).sqrt() + PHI1 / 2 * (1 / (2 * EPS.sqrt()) + C1)
    E1tip = lambda x: (2 * CA / rho(x)).sqrt()
    Kline = EPSI ** (-arb(1) / 4) / (2 * arb(2).sqrt())
    lc = arb(100); dl = arb('1e-4'); pi = arb.pi()
    s1 = 2 * X0R / (1 + X0R ** 2); s2 = (6 * X0R ** 2 - 2) / (1 + X0R ** 2) ** 2
    kV1 = arb(1) / 2 + c1W.sqrt()
    a_ = 1 + s1 + pi / (2 * lc)
    b_ = 2 * s1 * kV1 + s2 + (pi / lc) * (kV1 + s1) + pi ** 2 / (2 * lc ** 2)
    c_om = max(1 + EPSI * b_ ** 2 * (1 + 1 / dl), a_ ** 2 * (1 + dl), key=lambda z: float(z.upper()))
    assert slabs_ell_min >= 2 * (34.25 + 100), slabs_ell_min           # support inside [10.25, ell~/2]: no flattening
    E1I = lambda x: (1 + x * x) * (Kline * c_om / MU).sqrt()
    H2 = lambda x: (x + 2 * x ** 3 / 3 + x ** 5 / 5) - (XS + 2 * XS ** 3 / 3 + XS ** 5 / 5)
    KS = Ke(XS)
    ce = lambda x: Ke(x) if float(x.upper()) <= 16.25 else (KS ** 2 + H2(x) / MU).sqrt()
    ce1 = lambda x: E1tip(x) if float(x.upper()) <= 16.25 else (E1I(x) if float(x.lower()) > float(X1.upper()) else min(E1tip(x), E1I(x), key=lambda z: float(z.upper())))
    sup_e = ce(XE); sup_e1 = max(E1tip(XS), E1I(XE), key=lambda z: float(z.upper()))
    pts = [arb(0), X0R, XS, arb('24.25'), arb(28), XM, X1, XE]
    tab = [dict(x=float(p.mid()), c_e=up(ce(p)), c_e1=up(ce1(p))) for p in pts]
    out = dict(c_e=up(sup_e), c_e1=up(sup_e1), table=tab, K_e_16=up(KS), C_A=up(CA), K_line=up(Kline), c_omega=up(c_om), a=up(a_), b=up(b_),
               sigma1=up(s1), sigma2=up(s2), Phi1=up(PHI1), l_c=100, delta=1e-4,
               statement='|e(x)| <= c_e(x) Lambda_e^1/2, |e\'(x)| <= c_e1(x) Lambda_e^1/2 on [0, 34.249]; sup: c_e, c_e1 (Lemma 4b)',
               far='|e(x)| <= (K_e(16.25)^2 + H(x)^2/mu_e)^1/2 Lambda_e^1/2 for all x in [16.25, ell~], H^2 = int_16.25^x (1+s^2)^2')
    return out, dict(KS=KS, MU=MU, H2=H2, ce_sup=sup_e, ce1_sup=sup_e1, PHI1=PHI1, XE=XE, X1=X1, ce=ce)


# ============================================================================================ polynomial bookkeeping
def mmul(k1, k2):
    d = dict(k1)
    for a, e in k2: d[a] = d.get(a, Fr(0)) + e
    return tuple(sorted((a, e) for a, e in d.items() if e != 0))


class Pl:
    """sum of monomials c * prod atom^exponent with c >= 0 (arb upper bounds), exponents Fractions"""
    __slots__ = ('d',)
    def __init__(s, d=None): s.d = d if d is not None else {}
    @staticmethod
    def a(name, c=1, e=1): return Pl({((name, Fr(e)),): arb(c)})
    @staticmethod
    def c(v): return Pl({(): arb(v)})
    def __add__(s, o):
        if not isinstance(o, Pl): o = Pl.c(o)
        d = dict(s.d)
        for k, v in o.d.items(): d[k] = d.get(k, ZERO) + v
        return Pl(d)
    __radd__ = __add__
    def __mul__(s, o):
        if not isinstance(o, Pl): o = arb(o); return Pl({k: v * o for k, v in s.d.items()})
        d = {}
        for k1, v1 in s.d.items():
            for k2, v2 in o.d.items():
                k = mmul(k1, k2); d[k] = d.get(k, ZERO) + v1 * v2
        return Pl(d)
    __rmul__ = __mul__
    def subs(s, name, P):
        out = Pl()
        for k, v in s.d.items():
            e = dict(k).get(name)
            if e is None: out = out + Pl({k: v}); continue
            assert e.denominator == 1, (name, e)
            term = Pl({tuple((a, x) for a, x in k if a != name): v})
            for _ in range(int(e)): term = term * P
            out = out + term
        return out
    def deg(s): return {k: sum(e for _, e in k) for k in s.d}


def sqrtprod(P1, P2):
    """sqrt(P1 P2) <= sum_ij sqrt(c_i c_j m_i m_j) (nonnegative terms)"""
    d = {}
    for k1, v1 in P1.d.items():
        for k2, v2 in P2.d.items():
            k = tuple(sorted((a, e / 2) for a, e in mmul(k1, k2))); d[k] = d.get(k, ZERO) + sq(v1 * v2)
    return Pl(d)


CL10 = {'A', 'B', 'D4', 'V0', 'V2', 'V4', 'W1', 'W3', 'V1', 'V3'}
BASIC = ('A', 'B', 'V0', 'V2', 'SL', 'D4', 'V4')
BOXABLE = ('A', 'B', 'V0', 'V2', 'SL')


class Reducer:
    """turns a Pl into closure terms (coef, a, b).  Quadratic monomials in the 10 closure atoms pass unchanged (closure_cert
    optimises their AM-GM split); sqrt(Lambda) X uses AM-GM with t = (cap_SL/cap_X)^1/2; every other monomial is expanded into
    the basic atoms (S, S3, S3z and the sqrt atoms by their bounds), its excess degree is removed by the box x_i <= box_i on
    the box-able atoms (greedily, smallest box_i cap_i^1/2 first: A and sqrt(Lambda) before V0, B, V2), and the remaining
    degree-2 monomial is split by the weighted AM-GM  prod x_i^a_i <= sum_i (a_i/2) lam_i x_i^2,  lam_i = cap_i/prod cap^(a/2)
    (any positive weights are valid; cap_i only steers the split towards atoms with closure room).  Dissipation degree
    (exponents of D4, V4) <= 2 and total degree >= 2 are asserted."""
    CAP = dict(A=arb('0.01'), SL=arb('0.01'), B=arb('0.0125'), V0=arb('0.0027'), V2=arb('0.009'), D4=arb('2.5e-6'), V4=arb('2.7e-7'))
    def __init__(s, MU, c1W, Ybox, expand):
        s.MU = MU; s.Y = arb(Ybox); s.expand = expand; c = s.CAP
        s.box = dict(A=s.Y, SL=s.Y, B=s.Y / EPS.sqrt(), V0=s.Y / MU.sqrt(), V2=s.Y / (MU * EPSI).sqrt())
        s.cap = dict(c); s.cap.update(W1=(c['A'] * c['B']).sqrt(), W3=(c['B'] * c['D4']).sqrt(), V1=(c['V0'] * c['V2']).sqrt(),
                                      V3=(c['V2'] * c['V4']).sqrt())
        s.order = sorted(BOXABLE, key=lambda a: float((s.box[a] * s.cap[a].sqrt()).mid()))
        s.terms = {}; s.maxdiss = Fr(0); s.n_reduced = 0
    def emit(s, c, a, b):
        key = tuple(sorted((a, b))); s.terms[key] = s.terms.get(key, ZERO) + c
    def lam(s, c):                                   # c * Lambda
        s.emit(c, 'A', 'A'); s.emit(c * EPS, 'B', 'B'); s.emit(c * s.MU, 'V0', 'V0'); s.emit(c * s.MU * EPSI, 'V2', 'V2')
    def add(s, P):
        for k, c in P.d.items():
            if float(c.upper()) == 0: continue
            d = dict(k); tot = sum(d.values())
            if tot == 2 and all(e.denominator == 1 for e in d.values()) and set(d) <= CL10:
                at = [a for a, e in d.items() for _ in range(int(e))]; s.emit(c, at[0], at[1]); continue
            if tot == 2 and set(d) <= CL10 | {'SL'} and d.get('SL') == 1:
                X = [a for a in d if a != 'SL'][0]; t = (s.cap['SL'] / s.cap[X]).sqrt()
                s.lam(c * t / 2); s.emit(c / (2 * t), X, X); continue
            if tot == 2 and d == {'SL': Fr(2)}: s.lam(c); continue
            Q = Pl({k: c})
            for name in ('S3z', 'S3', 'S', 'W1', 'W3', 'V1', 'V3'):
                if name in d: Q = Q.subs(name, s.expand[name])
            for k2, c2 in Q.d.items(): s.mono(c2, dict(k2))
    @staticmethod
    def pw(x, e): return (x.log() * arb(e.numerator) / arb(e.denominator)).exp()
    def mono(s, c, d):
        assert set(d) <= set(BASIC), d
        p = sum(d.get(a, Fr(0)) for a in ('D4', 'V4')); q = sum(d.get(a, Fr(0)) for a in BOXABLE); tot = p + q
        assert p <= 2, ('dissipation degree > 2', d)
        assert tot >= 2, ('sub-quadratic monomial', d)
        s.maxdiss = max(s.maxdiss, p); s.n_reduced += 1
        exc = tot - 2; dd = dict(d); fac = arb(1)
        for a in s.order:
            if exc <= 0: break
            if a in dd and dd[a] > 0:
                cut = min(dd[a], exc); dd[a] -= cut; exc -= cut; fac = fac * s.pw(s.box[a], cut)
        assert exc == 0
        dd = {a: e for a, e in dd.items() if e > 0}
        norm = arb(1)
        for a, e in dd.items(): norm = norm * s.pw(s.cap[a], e / 2)
        for a, e in dd.items():
            w = c * fac * arb(e.numerator) / arb(2 * e.denominator) * s.cap[a] / norm
            if a == 'SL': s.lam(w)
            else: s.emit(w, a, a)
    def out(s): return [[up(v), a, b] for (a, b), v in sorted(s.terms.items())]


# ============================================================================================ base data
_JC = {}
def _jets(s, jdir):
    if (s, jdir) not in _JC: _JC[(s, jdir)] = T6R.jets(s, jdir)
    return _JC[(s, jdir)]


def t6_hulls(slabs, lo, hi, keys, jdir):
    """hull over the slabs and over all T6 sub-cells meeting [lo, hi] (certified T6 jets, widened 2 ulps by the reader)"""
    out = {k: None for k in keys}
    for s in slabs:
        J = _jets(s, jdir); X = J['x']
        m = (X[..., 1] >= lo) & (X[..., 0] <= hi)
        for k in keys:
            H = J[k][m]; a, b = float(H[:, 0].min()), float(H[:, 1].max())
            if out[k] is None: out[k] = [a, b]
            else: out[k] = [min(out[k][0], a), max(out[k][1], b)]
    return {k: arb((v[0] + v[1]) / 2, (v[1] - v[0]) / 2) + arb(0, 1e-300) for k, v in out.items()}


def frenet(k, p0, q0, n0, t0):
    """(T8 n3far.frenet) Frenet jets of p~, q~, n~, t1~ up to order 4"""
    from math import comb
    P = [p0]; Q = [q0]; N = [n0]; T1 = [t0]
    for m in range(4):
        P.append(-sum((comb(m, l) * k[l] * Q[m - l] for l in range(m + 1)), ZERO))
        Q.append((arb(1) if m == 0 else ZERO) + sum((comb(m, l) * k[l] * P[m - l] for l in range(m + 1)), ZERO))
        N.append(-sum((comb(m, l) * k[l] * T1[m - l] for l in range(m + 1)), ZERO))
        T1.append(sum((comb(m, l) * k[l] * N[m - l] for l in range(m + 1)), ZERO))
    return P, Q, N, T1


def g_bounds(Vb, eta1, th, xbar, cabs, babs):
    """majorants of g = U/alpha + 5 th k k', U = (x1 + th c) q + (x2 + th b) t1 (|x_i| <= xbar), with the box |e|, |e'| <= eta1
    and e'' = e''' = 0 (T8 q_strong with Z = (eta1, eta1, 0, 0)); the e''- and e'''-dependence is affine and returned
    separately: g'' += G2e e'', g''' += G3e e'' + Gtop e'''.   Returns G0..G3 (flat), G2e, G3e, Gtop."""
    from math import comb
    k = [Vb[f'k{m}'] for m in range(6)]
    P, Q, N, T1 = frenet(k, Vb['p'], Vb['q'], Vb['n'], Vb['t1'])
    Ve = {f'k{m}': k[m] for m in range(5)}
    Ve.update(e0=arb(0, eta1), e1=arb(0, eta1), e2=ZERO, e3=ZERO, e4=ZERO)
    inv = n3qforms.invalpha(Ve)
    dc = arb(xbar) + th * cabs; dt = arb(xbar) + th * babs
    Ek = lambda F, kk: sum((comb(kk, l) * F[l] * inv[kk - l] for l in range(kk + 1)), ZERO)
    Gs = []
    for kk in range(4):
        Ekap = sum((comb(kk, l) * k[l] * k[kk - l + 1] for l in range(kk + 1)), ZERO)
        Gs.append(ab(dc * Ek(Q, kk) + dt * Ek(T1, kk) + 5 * th * Ekap))
    al_min = 1 - ab(k[0]) * eta1
    Uabs = dc * ab(Q[0]) + dt * ab(T1[0]); U1abs = dc * ab(Q[1]) + dt * ab(T1[1])
    Gtop = Uabs * ab(k[0]) / al_min ** 2
    a1 = (ab(k[1]) + ab(k[0])) * eta1                               # |alpha'|
    G2e = Gtop                                                       # d/de'' (1/alpha)'' = k/alpha^2
    G3e = 3 * U1abs * ab(k[0]) / al_min ** 2 + Uabs * (3 * ab(k[1]) / al_min ** 2 + 6 * a1 * ab(k[0]) / al_min ** 3)
    return Gs, G2e, G3e, Gtop


def dn_bounds(H, eta1):
    """D(N) monomial coefficient bounds on the box |e|, |e'| <= eta1 (TD2 enclosures): order 0 -> sup|c|, order 1 ->
    sup|d0 c| eta1 + sup|d1 c| eta1, order 2 (constant term) -> c <= Hm eta1 w  (|c| <= Hm w^2/2, w = |e| + |e'| <= 2 eta1)"""
    V = {k: H[k] for k in ('k0', 'k1', 'k2', 'k3', 'k4', 'k5', 'r1', 'r2')}
    E = arb(0, eta1); V['e0'] = TD2.var(E, 0); V['e1'] = TD2.var(E, 1)
    out = G.dn(V); res = {}
    ORD0 = {(1, 1, 0, 0), (2, 0, 0, 0), (2, 1, 0, 0)}
    for m, c in zip(G.DN_MONOS, out):
        if not isinstance(c, TD2): c = TD2(c)
        if m == (0, 0, 0, 0):
            Hm = max(ab(c.h00), ab(c.h01), ab(c.h11), key=lambda z: float(z.upper())); res[m] = Hm * eta1   # coefficient on w
        elif m in ORD0: res[m] = ab(c.v)
        else: res[m] = (ab(c.g0) + ab(c.g1)) * eta1
    return res


def psi_bounds(H, eta1):
    """weak (M) form: per m (jet s_m of the test function): H0 (Hessian of P0 in (e, e')), L1 (gradient of P1), V2 (|P2|)"""
    V = {k: H[k] for k in ('k0', 'k1', 'k2')}
    E = arb(0, eta1); V['e0'] = TD2.var(E, 0); V['e1'] = TD2.var(E, 1)
    o = [c if isinstance(c, TD2) else TD2(c) for c in G.psi(V)]
    mx = lambda *v: max((ab(z) for z in v), key=lambda z: float(z.upper()))
    return [dict(H0=mx(o[3 * m].h00, o[3 * m].h01, o[3 * m].h11), L1=mx(o[3 * m + 1].g0, o[3 * m + 1].g1), V2=ab(o[3 * m + 2].v))
            for m in range(3)]


def phi_sups():
    """sup over [0, XE] of |s_m / rho| for psi = rho phi^_j (phi^_j = phi_j/|phi_j|_rho): s0/rho = f, s1/rho = f' - phi' f,
    s2/rho = f'' - 2 phi' f' + (phi'^2 - phi'') f; max over j.  Fixed functions of cert7 (coarse mesh, Taylor models)."""
    import base7, certdata as CD, cert7
    base7.set_base(os.path.join('..', 't1', 'polished_v2'))
    Sc = base7.Slab(100); Bc = Sc.base(arb(-1), CD.XE + 1); BL = base7.late_base(CD.XE + 1)
    nodes = CD.mesh(Bc, 0.1); FX = cert7.Fixed(cert7.coarse_nodes(Bc, nodes, 0.5), BL)
    nphi = [sum(((fx['e1'] * fx['phis'][j][0] * fx['phis'][j][0]).integral01() * fx['h'] for fx in FX.data), ZERO).sqrt() for j in range(3)]
    S = [ZERO] * 3
    for fx in FX.data:
        p1, p2 = fx['d'][0], fx['d'][1]
        for j in range(3):
            f = fx['phis'][j]
            vals = [f[0].sup(), (f[1] - p1 * f[0]).sup(), (f[2] - 2 * (p1 * f[1]) + (p1 * p1 - p2) * f[0]).sup()]
            for m in range(3): S[m] = max(S[m], vals[m] / nphi[j], key=lambda z: float(z.upper()))
    rint = sum(((fx['e1']).integral01() * fx['h'] for fx in FX.data), ZERO)          # int_0^XE rho
    return S, rint


# ============================================================================================ per interval
def interval(cert, C, N3, R6, box, LB, PHS):
    MU = U(C['meta']['closure']['mu_e']); c1W = U(C['meta']['weight_consts']['c1W']); PSI1 = U(C['meta']['weight_consts']['Psi1'])
    sl = [s for s in range(cert['slabs'][0], cert['slabs'][1] + 1) if str(s) in N3['per_slab']]
    if not sl: return None
    S7 = [C['slabs'][str(s)] for s in sl]
    Y = arb(box['Ybox']); eta1 = LB['ce1_sup'] * Y                   # |e|, |e'| <= eta1 on [0, XE] (Lemma 4b; c_e <= c_e1)
    XE, X1 = LB['XE'], LB['X1']
    # ---- box values (Lambda <= Ybox^2)
    s_ = (1 / (2 * EPS.sqrt()) + C1).sqrt()
    kap1 = (1 / (2 * MU * EPSI.sqrt()) + c1W / MU).sqrt(); kap2 = 1 / (MU * EPSI).sqrt()
    W34 = Wf(XE); W16 = Wf(XS)
    cSW = ((kap1 ** 2 / 18 + 2 * kap1 * kap2) / W34).sqrt()           # S_W <= cSW sqrt(Lambda)  (Lemma 4a iii)
    SWbox = cSW * Y
    Kreg = 1 / rho(XS); kR = Wf(X0R) / rho(XS)                       # kappa_R = sup_R W/rho
    sig = (rho(X1) / W34).sqrt()                                     # strip [X1, XE]: rho <= sig^2 W
    sig0 = LB['ce_sup'] * (3 * rho(X1)).sqrt()                       # (int_strip rho e^2)^1/2 <= sig0 sqrt(Lambda)
    Zs = [arb('2.1875') / 6, arb('7.5300') / 36, None, None]        # sup|zeta'|, |zeta''|: S' max 35/16, S'' max <= 7.53
    # sup|S'''|, |S''''| of S(u) = u^4 (35 - 84u + 70u^2 - 20u^3) on [0,1] by dense ball evaluation
    from flint import arb_poly
    Sp = arb_poly([0, 0, 0, 0, 35, -84, 70, -20])
    D3, D4p = Sp.derivative().derivative().derivative(), Sp.derivative().derivative().derivative().derivative()
    D2 = Sp.derivative().derivative()
    sup_on = lambda P: max((ab(P(arb((i + 0.5) / 200, 0.5 / 200))) for i in range(200)), key=lambda z: float(z.upper()))
    Zs = [sup_on(Sp.derivative()) / 6, sup_on(D2) / 36, sup_on(D3) / 216, sup_on(D4p) / 1296]
    Z1, Z2, Z3, Z4 = Zs
    # ---- modulation: |x^| <= xbar, |theta| <= thbar
    cphi = None
    for c in C['certificates']:
        v = c['c_phi']; cphi = v if cphi is None else [max(a, b) for a, b in zip(cphi, v)]
    CPHI = sum((U(v) ** 2 for v in cphi), ZERO).sqrt()
    R2 = max(U(R6[s]['norms']['R2']) for s in sl) * arb('1.001')
    mbar = max(U(e['C3']['mbar']) for e in S7); Cx = max(U(e['C3']['Cx_full']) for e in S7)
    lbar = max(U(e['C3']['lbar']) for e in S7)
    nV = min(U(e['C3']['nV_min']) for e in S7)
    xr = mbar * CPHI * R2
    rc = [R6[s]['scalars']['range_c'] for s in sl]; rb = [R6[s]['scalars']['range_b'] for s in sl]
    cabs = max(max(abs(v) for v in r) for r in rc); babs = max(max(abs(v) for v in r) for r in rb)
    cabs, babs = U(cabs), U(babs)
    xnl = box['xnl']                                                 # (M) nonlinear part at the box (fixed point, see main)
    xbar = Cx * Y + xr + xnl; thbar = xbar / nV; f1 = 1 + thbar
    dc = xbar + thbar * cabs
    fx = max(xbar / U(N3['per_slab'][str(s)]['mod']['xbar']) for s in sl)      # scale T8's x-dependent majorants
    fx = max(fx, arb(1), key=lambda z: float(z.upper()))
    # ---- tip constants (nl7, T5 box |e|, |e'| <= 1e-5 >= eta1) and base sups
    assert float(eta1.upper()) <= 1e-5
    tipk = {k: max(U(e['nonlinear']['tip'][k]) for e in S7) for k in ('n4', 'n3_0', 'm3', 'n2_0', 'n21', 'm2', 'm23')}
    n4, n30, n20, n21, m2, m23 = (tipk[k] for k in ('n4', 'n3_0', 'n2_0', 'n21', 'm2', 'm23'))
    jdir = JDIR
    keys = ['k0', 'k1', 'k2', 'k3', 'k4', 'k5', 'r1', 'r2', 'r3', 'p', 'q', 'n', 't1']
    def sups(lo, hi):
        Hh = t6_hulls(sl, lo, hi, keys, jdir); return Hh, {k: ab(v) for k, v in Hh.items()}
    Htip, St = sups(0.0, float(XE.upper()))
    K0, K1, K2, K3, K4 = (St[f'k{m}'] for m in range(5)); R1s, R2s = St['r1'], St['r2']
    alm = 1 - K0 * eta1; alM = 1 + K0 * eta1
    m3a = 10 * eta1 / alm ** 6                                       # sup|A31| = sup 10|e'|/J^6 (Lemma 4a i), J >= alpha
    C0 = 3 * K1 ** 2 + 4 * K0 * K2; C0p = 10 * K1 * K2 + 4 * K0 * K3; KP = K3 + R1s
    Gt, G2e_t, G3e_t, Gtop_t = g_bounds(Htip, eta1, thbar, xbar, cabs, babs)
    # ---- atoms
    A, B, D4, V0, V2, V4 = (Pl.a(n) for n in ('A', 'B', 'D4', 'V0', 'V2', 'V4'))
    W1, W3, V1, V3, SL, S, S3, S3z = (Pl.a(n) for n in ('W1', 'W3', 'V1', 'V3', 'SL', 'S', 'S3', 'S3z'))
    Nk = [A + SL * sig0, W1 + V0 * sig, B + V1 * sig, W3 + V2 * sig, D4 + V3 * sig]   # full rho-norms of e^(k) on [0, XE]
    parts = {}
    # =============================== (E0)  <chi(N + Q), u>_rho
    E0 = (A * Nk[4] * (2 * n4 * eta1) + A * Nk[3] * (2 * n30 * eta1) + Nk[2] * Nk[3] * (m3a * eta1)
          + (Nk[0] + Nk[1]) * (Nk[0] + Nk[1]) * (n20 * eta1) + (Nk[0] + Nk[1]) * Nk[2] * (n21 * eta1) + Nk[2] * Nk[2] * (m2 * eta1)
          + S * Nk[2] * Nk[2] * (m23 * eta1)) * f1
    E0 = E0 + A * (Nk[4] + Nk[2] * K0 ** 2 + Nk[0] * (C0 + KP * eta1 / alm)) * thbar + A * Nk[1] * Gt[0] + A * Nk[0] * dc
    parts['E0'] = E0
    # =============================== (E2)  eps <chi(N + Q), D_rho u>_rho <= eps |D_rho u| |chi(N + Q)|
    Delta = D4 + W3 * (2 * LB['PHI1']) + B * PHI2
    NQ = (Nk[4] * (2 * n4 * eta1) + Nk[3] * (2 * n30 * eta1) + S * Nk[3] * m3a + (Nk[0] + Nk[1]) * (2 * n20 * eta1)
          + Nk[2] * (2 * n21 * eta1) + S * Nk[2] * m2 + S * S * Nk[2] * m23) * f1
    NQ = NQ + (Nk[4] + Nk[2] * K0 ** 2 + Nk[0] * (C0 + KP * eta1 / alm)) * thbar + Nk[1] * Gt[0] + Nk[0] * dc
    parts['E2'] = Delta * NQ * EPS
    # =============================== region data for E_I
    def tilde(Sx, eta1, e2):
        """analytic majorants of At = 1 - J^-4 and its derivatives on a region with sups Sx, box |e|,|e'| <= eta1, |e''| <= e2"""
        k0, k1 = Sx['k0'], Sx['k1']; am = 1 - k0 * eta1; aM = 1 + k0 * eta1; J6 = am ** 6
        a4e = 4 * k0 * aM / J6; a4e1 = 4 * eta1 / J6; a4X = 4 * k1 * eta1 * aM / J6
        return dict(A=(a4e + a4e1) * eta1, a4e=a4e, a4e1=a4e1, a4X=a4X, A1f=a4X + a4e * eta1, am=am)
    HR, SR = sups(10.25, 16.25); HM, SM = sups(16.25, 30.25)
    tR = tilde(SR, eta1, None); tM = tilde(SM, eta1, SWbox)
    om1M = 4 * XS / (1 + XS ** 2); om2M = (20 * XS ** 2 - 4) / (1 + XS ** 2) ** 2
    intWM = 1 / (3 * XS ** 3)                                        # int_M W <= int_16.25^oo x^-4
    tau1R = PSI1 + 2 * Z1; tau2R = 2 * c1W + 4 * PSI1 * Z1 + 2 * Z1 ** 2 + 2 * Z2
    GR, G2e_R, G3e_R, Gtop_R = g_bounds(HR, eta1, thbar, xbar, cabs, babs)
    GM, G2e_M, G3e_M, Gtop_M = g_bounds(HM, eta1, thbar, xbar, cabs, babs)
    ce = LB['ce_sup']
    # =============================== (I0) on R (T = W zeta^2; jets via rho-norms, u = e on R)
    AR = tR['A'] * f1 + thbar
    I0R = (W3 * W3 * AR + B * W3 * (2 * tau1R * AR + (tR['A1f'] + 2 * n30 * eta1) * f1) + S * B * W3 * ((tR['a4e1'] + m3a) * f1)
           + W1 * W3 * (tau2R * AR + tau1R * (tR['A1f'] + 2 * n30 * eta1) * f1) + S * W1 * W3 * (tau1R * (tR['a4e1'] + m3a) * f1))
    CT = ((A + W1) * (A + W1) * n20 + (A + W1) * B * n21 + B * B * m2 + S * B * B * m23) * (tau1R * eta1 * f1)
    CTT = ((A + W1) * B * (2 * n20 * eta1) + B * B * (2 * n21 * eta1) + S * B * B * m2 + S * S * B * B * m23) * f1
    CTh = (W1 * tau1R + B) * (B * SR['k0'] ** 2 + A * (C0 + KP * eta1 / alm)) * thbar
    I0R = (I0R + CT + CTT + CTh + W1 * W1 * ((GR[1] + GR[0] * tau1R) / 2 + dc)) * kR
    # =============================== (I0) on M (T = W, zeta = 1, |e''| <= S_W box)
    AM = tM['A'] * f1 + thbar; A1M = (tM['A1f'] + tM['a4e1'] * SWbox) * f1; B0M = (2 * n30 * eta1 + m3a * SWbox) * f1
    gaE = (2 * n20 * eta1) * f1 + thbar * (C0 + KP * eta1 / alm); ga0 = (2 * n20 * eta1) * f1
    ga1 = (2 * n21 * eta1 + m2 * SWbox + m23 * SWbox ** 2) * f1 + thbar * SM['k0'] ** 2
    I0M = (V2 * V2 * AM + V1 * V2 * (2 * om1M * AM + A1M + B0M) + V0 * V2 * (om2M * AM + om1M * (A1M + B0M))
           + V0 * V0 * (om1M * ga0 + (GM[1] + GM[0] * om1M) / 2 + dc) + V0 * V1 * (om1M * ga1 + ga0) + V1 * V1 * ga1
           + SL * V0 * (om1M * gaE * ce * intWM.sqrt()) + SL * V1 * (gaE * ce * intWM.sqrt()))
    # =============================== (I0), (I2) on the far core x >= 30.249 (T8 n3.json per range, T8's certified boxes)
    def far_slab(s):
        r = N3['per_slab'][str(s)]; th8 = U(r['mod']['thbar'])
        acc = {}
        def mx(k, v): acc[k] = max(acc.get(k, ZERO), v, key=lambda z: float(z.upper()))
        eint = dict(I0V0=ZERO, I0V1=ZERO, I2V4=ZERO, I2V3=ZERO, I2V2=ZERO)
        for g in r['ranges']:
            if not g['core_box']: continue
            xlo, xhi = arb(g['x'][0]), arb(g['x'][1])
            o1 = min(PSI1, 8 / xlo, key=lambda z: float(z.upper())); o2 = min(2 * c1W, 2264 / xlo ** 2, key=lambda z: float(z.upper()))
            Atp_e = f1 * U(g['Atp']) + thbar * (1 + U(g['A'])); A_e = f1 * U(g['A']) + thbar
            A1e, B0e = f1 * U(g['A1']), f1 * U(g['B0'])
            kr = U(g['ke']) / U(g['Z'][0])
            c3r = A1e + B0e; c2r = f1 * (U(g['B1']) + U(g['dCe'])) + thbar * kr ** 2
            b0 = [U(v) * fx for v in g['beta']['0']]; b2 = [U(v) * fx for v in g['beta']['2']]
            Gg = [U(v) * fx for v in g['G']]
            mx('top0', Atp_e); mx('I0_12', 2 * o1 * A_e + A1e + B0e); mx('I0_02', o2 * A_e + o1 * (A1e + B0e))
            mx('I0_00', o1 * b0[1] + (Gg[1] + Gg[0] * o1) / 2 + dc); mx('I0_01', o1 * b0[2] + b0[1]); mx('I0_11', b0[2])
            be1 = b2[1] + dc
            mx('I2_44', Atp_e); mx('I2_43', c3r + 2 * o1 * A_e); mx('I2_42', c2r + o2 * A_e); mx('I2_41', b2[2]); mx('I2_40', be1)
            mx('I2_33', 2 * o1 * c3r); mx('I2_32', 2 * o1 * c2r + o2 * c3r); mx('I2_31', 2 * o1 * b2[2]); mx('I2_30', 2 * o1 * be1)
            mx('I2_22', o2 * c2r + (5 * Gg[1] + Gg[0] * o1) / 2); mx('I2_21', o2 * b2[2] + 3 * Gg[2]); mx('I2_20', o2 * be1 + Gg[3])
            K2x = LB['KS'] ** 2 + LB['H2'](xhi) / LB['MU']
            IW = (xhi - xlo) / (1 + xlo ** 2 / 4) ** 2 * K2x            # int_r W K^2 (W <= (1 + x^2/4)^-2 with flattening)
            eint['I0V0'] += o1 ** 2 * b0[0] ** 2 * IW; eint['I0V1'] += b0[0] ** 2 * IW
            eint['I2V4'] += b2[0] ** 2 * IW; eint['I2V3'] += 4 * o1 ** 2 * b2[0] ** 2 * IW; eint['I2V2'] += o2 ** 2 * b2[0] ** 2 * IW
        for k, v in eint.items(): acc['e_' + k] = sq(v)
        acc['sup_Atp'] = U(r['sup_Atp'])
        return acc
    FAR = {}
    for s in sl:
        for k, v in far_slab(s).items(): FAR[k] = max(FAR.get(k, ZERO), v, key=lambda z: float(z.upper()))
    I0F = (V2 * V2 * FAR['top0'] + V1 * V2 * FAR['I0_12'] + V0 * V2 * FAR['I0_02'] + V0 * V0 * FAR['I0_00'] + V0 * V1 * FAR['I0_01']
           + V1 * V1 * FAR['I0_11'] + SL * V0 * FAR['e_I0V0'] + SL * V1 * FAR['e_I0V1'])
    I2F = (V4 * V4 * FAR['I2_44'] + V4 * V3 * FAR['I2_43'] + V4 * V2 * FAR['I2_42'] + V4 * V1 * FAR['I2_41'] + V4 * V0 * FAR['I2_40']
           + V3 * V3 * FAR['I2_33'] + V3 * V2 * FAR['I2_32'] + V3 * V1 * FAR['I2_31'] + V3 * V0 * FAR['I2_30']
           + V2 * V2 * FAR['I2_22'] + V2 * V1 * FAR['I2_21'] + V2 * V0 * FAR['I2_20']
           + SL * V4 * FAR['e_I2V4'] + SL * V3 * FAR['e_I2V3'] + SL * V2 * FAR['e_I2V2'])
    # =============================== (I2) on M and R: D(N) monomials (TD2 enclosures) per 2-wide sub-range
    def dn_reg(lo, hi):
        res = None; x = lo
        while x < hi - 1e-9:
            Hh = t6_hulls(sl, x, min(x + 2.0, hi), keys, jdir); d = dn_bounds(Hh, eta1)
            res = d if res is None else {m: max(res[m], d[m], key=lambda z: float(z.upper())) for m in d}
            x += 2.0
        return res
    CM = dn_reg(16.25, 30.25); CR = dn_reg(10.25, 16.25)
    def lin(Cm, e2, p_e5, p_e4, p_e3, p_e2, p_e1, p_e, e3sq, k0s, k1s, k3s, k4s):
        """|zeta (1-th) dN + th (...)| as a Pl, given the jet norms (Pl) of zeta e^(k) and a pointwise e'' bound e2 (Pl or const)"""
        P = (p_e5 * Cm[(0, 0, 0, 1)] + p_e4 * (Cm[(0, 0, 1, 0)] + e2 * Cm[(1, 0, 1, 0)])
             + p_e3 * (Cm[(0, 1, 0, 0)] + e2 * Cm[(1, 1, 0, 0)] + e2 * e2 * Cm[(2, 1, 0, 0)]) + e3sq * Cm[(0, 2, 0, 0)]
             + p_e2 * (Cm[(1, 0, 0, 0)] + e2 * Cm[(2, 0, 0, 0)] + e2 * e2 * Cm[(3, 0, 0, 0)] + e2 * e2 * e2 * Cm[(4, 0, 0, 0)])
             + (p_e + p_e1) * Cm[(0, 0, 0, 0)]) * f1
        Pth = (p_e5 + p_e3 * k0s ** 2 + p_e2 * (2 * k0s * k1s + KP * eta1 / alm) + p_e1 * (C0 + KP * eta1 / alm)
               + p_e * (C0p + (k4s + R2s) * eta1 / alm + KP * (k1s + k0s) * eta1 ** 2 / alm ** 2)) * thbar + p_e1 * dc
        return P + Pth
    # M: zeta = 1, e^(k+1) = v^(k); |e''| <= S_W (box); e'''^2 <= S3 |v''|; |e| <= c_e sqrt(Lambda) (M inside [0, XE])
    LM = lin(CM, Pl.c(SWbox), V4, V3, V2, V1, V0, SL * (ce * intWM.sqrt()), S3 * V2, SM['k0'], SM['k1'], SM['k3'], SM['k4'])
    I2M_strong = (V4 + V3 * (2 * om1M) + V2 * om2M) * LM
    I2M_tr = (V2 * V2 * ((5 * GM[1] + GM[0] * om1M) / 2 + Gtop_M * eta1) + V2 * V1 * (3 * (GM[2] + G2e_M * SWbox))
              + V2 * V0 * (GM[3] + G3e_M * SWbox))
    # R: one zeta per monomial, assigned to the highest jet; e'' <= S; zeta e''' sup S3z; rho-norms with kappa_R
    sk = kR.sqrt()
    zE5 = V4 + (D4 * (4 * Z1) + W3 * (6 * Z2) + B * (4 * Z3) + W1 * Z4) * sk
    zE4 = V3 + (W3 * (3 * Z1) + B * (3 * Z2) + W1 * Z3) * sk
    zE3 = V2 + (B * (2 * Z1) + W1 * Z2) * sk
    zE2 = V1 + W1 * (Z1 * sk)
    LR = lin(CR, S, zE5, zE4, zE3, zE2, W1 * sk, A * sk, S3z * W3 * sk, SR['k0'], SR['k1'], SR['k3'], SR['k4'])
    Theta = V4 + V3 * (2 * PSI1) + V2 * (2 * c1W)
    I2R_strong = Theta * LR
    I2R_tr = (W3 * W3 * ((5 * GR[1] + GR[0] * tau1R) / 2 + Gtop_R * eta1) + W3 * B * (3 * GR[2]) + S * W3 * B * (3 * G2e_R)
              + W3 * W1 * GR[3] + S * W3 * W1 * G3e_R) * kR
    d3 = W1 * GR[3] + S * W1 * G3e_R + W3 * (Gtop_R * eta1) + B * (3 * GR[2]) + S * B * (3 * G2e_R) + W3 * (3 * GR[1]) + D4 * GR[0]
    com1 = (B * (2 * Z1) + W1 * Z2) * d3 * kR
    com2 = V2 * ((W1 * GR[2] + S * W1 * G2e_R + B * (2 * GR[1]) + W3 * GR[0]) * (2 * Z1) + (W1 * GR[1] + B * GR[0]) * Z2) * sk
    I2R = I2R_tr + com1 + com2
    parts['I0'] = (I0R + I0M + I0F) * MU
    parts['I2'] = (I2M_strong + I2M_tr + I2F + I2R_strong + I2R) * (MU * EPSI)
    # =============================== (M): nonlinear and E_ij parts of the modulation
    Hps = None
    xg = [0.0] + [float(v) for v in range(1, 34)] + [float(XE.upper())]
    for a_, b_ in zip(xg[:-1], xg[1:]):
        Hh = t6_hulls(sl, a_, b_, ['k0', 'k1', 'k2'], jdir); pb = psi_bounds(Hh, eta1)
        Hps = pb if Hps is None else [{k: max(o[k], n[k], key=lambda z: float(z.upper())) for k in o} for o, n in zip(Hps, pb)]
    Sph, rint = PHS
    w01 = Nk[0] + Nk[1]
    bNL = sum(((w01 * w01 * (Hps[m]['H0'] / 2) + w01 * Nk[2] * Hps[m]['L1'] + Nk[2] * Nk[2] * Hps[m]['V2']) * Sph[m] for m in range(3)), Pl())
    bNL = bNL + Nk[0] * Nk[1] * (Sph[0] * KP / alm)                  # P = (k''' - r') e e'/alpha
    ngm = [min(U(e['farg']['norm_chi_g_min'][i]) for e in S7) for i in range(2)]
    assert all(float(v.lower()) >= 1 for v in ngm), ngm               # |x_i| <= |x^_i| for i = 1, 2 (used in U)
    Qt = St['q']; Tt = St['t1']
    E1 = (Nk[1] * (Qt / alm) + Nk[0]) * (1 / ngm[0]); E2 = Nk[1] * (Tt / alm / ngm[1])
    cphim = max(U(v) for v in cphi)
    E3 = (Nk[2] * (cphim + K0 ** 2) + Nk[0] * (C0 + KP * eta1 / alm + cabs) + Nk[1] * (5 * K0 * K1 + (cabs * Qt + babs * Tt) / alm)
          + bNL) * (1 / nV)
    eR = cphim * R2 / nV                                              # theta part of (1 - theta) <S2, rho^-1(rho phi)''>
    def at_box(P):
        """value of a Pl at the box (atoms at their box values; S, S3, S3z by their box parts are not needed here)"""
        bv = dict(A=Y, B=Y / EPS.sqrt(), V0=Y / MU.sqrt(), V2=Y / (MU * EPSI).sqrt(), SL=Y, W1=s_ * Y, V1=kap1 * Y)
        tot = ZERO
        for k, c in P.d.items():
            v = c
            for a_, e in k: v = v * (bv[a_].log() * arb(e.numerator) / arb(e.denominator)).exp()
            tot += v
        return ub(tot)
    Ebox = arb(3).sqrt() * sq((at_box(E1) + eR) ** 2 + at_box(E2) ** 2 + (at_box(E3) + eR) ** 2)
    bbox = arb(3).sqrt() * at_box(bNL)
    if not float((mbar * Ebox).upper()) < 0.5:
        raise AssertionError(('mbar |E| >= 1/2', up(Ebox), up(mbar)))
    fM = mbar / (1 - mbar * Ebox)
    xNL_box = mbar * bbox + fM * Ebox * (Cx * Y + xr + mbar * bbox)
    # closure terms: |x_NL| <= mbar |b_NL| + fM |E| (Cx A + xr + mbar |b_NL|), |E| <= sqrt3 (sum_i |E_i|^2)^1/2 <= sqrt3 sum_i |E_i|
    Epl = (E1 + E2 + E3) * arb(3).sqrt()
    bpl = bNL * arb(3).sqrt()
    # |E| <= Epl + sqrt12 eR; the e-independent product sqrt12 eR xr fM is a residual increment of x_r (reported, not z-dependent)
    xNL = bpl * mbar + Epl * (A * Cx + xr + bpl * mbar) * fM + (A * Cx + bpl * mbar) * (eR * arb(12).sqrt() * fM)
    xr_inc = eR * arb(12).sqrt() * fM * xr
    g2E = max(U(e['farg']['sqrt_lmax_G2E']) for e in S7); g0 = max(U(e['farg']['sqrt_lmax_G0']) for e in S7)
    g2 = max(U(e['farg']['sqrt_lmax_G2']) for e in S7)
    parts['M'] = xNL * (A * arb(3).sqrt() + B * (EPS * g2E) + V0 * (MU * g0) + Theta * (MU * EPSI * g2))
    # =============================== reduce
    expand = dict(S=SL * cSW + Pl({(('B', Fr(3, 4)), ('D4', Fr(1, 4))): (2 * Kreg).sqrt()}) + B * ((2 * Kreg).sqrt() * C1 ** (arb(1) / 4)),
                  S3=(V2 * (1 / arb(18).sqrt() + arb(2).sqrt() * c1W ** (arb(1) / 4)) + Pl({(('V2', Fr(3, 4)), ('V4', Fr(1, 4))): arb(2).sqrt()})) * (1 / W34.sqrt()),
                  W1=Pl({(('A', Fr(1, 2)), ('B', Fr(1, 2))): arb(1)}) + A * C1.sqrt(),
                  W3=Pl({(('B', Fr(1, 2)), ('D4', Fr(1, 2))): arb(1)}) + B * C1.sqrt(),
                  V1=Pl({(('V0', Fr(1, 2)), ('V2', Fr(1, 2))): arb(1)}) + V0 * c1W.sqrt(),
                  V3=Pl({(('V2', Fr(1, 2)), ('V4', Fr(1, 2))): arb(1)}) + V2 * c1W.sqrt())
    def basic(P):
        for nm in ('W1', 'W3', 'V1', 'V3'): P = P.subs(nm, expand[nm])
        return P
    P1 = basic(V2 + (B * (2 * Z1) + W1 * Z2) * sk); P2 = basic(V3 + (W3 * (2 * Z1) + B * (3 * Z2) + W1 * Z3) * sk)
    expand['S3z'] = sqrtprod(P1, P2) * (2 / W16).sqrt()
    red = Reducer(MU, c1W, box['Ybox'], expand)
    if os.environ.get('NLDEBUG'):
        for k_, P in parts.items():
            big = sorted(P.d.items(), key=lambda kv: -float(kv[1].upper()))[:8]
            print('DEBUG', k_, [(dict(k), '%.2e' % float(v.upper())) for k, v in big])
    summ = {}
    for k_, P in parts.items():
        r1 = Reducer(MU, c1W, box['Ybox'], expand); r1.add(P)
        summ[k_] = dict(n_terms=len(r1.terms), max_coef=max((up(v) for v in r1.terms.values()), default=0.0),
                        diag={a: up(v) for (a, b), v in r1.terms.items() if a == b})
        red.add(P)
    top = dict(I0_V2sq=up(FAR['top0'] * MU), I2_V4sq=up(FAR['I2_44'] * MU * EPSI), sup_Atp=up(FAR['sup_Atp']))
    return dict(id=cert['id'], slabs=[sl[0], sl[-1]], terms=red.out(), max_dissipation_degree=float(red.maxdiss), n_reduced=red.n_reduced,
                parts=summ, top=top, M=dict(xNL_box=up(xNL_box), E_box=up(Ebox), b_box=up(bbox), mbar=up(mbar), Cx=up(Cx), xr=up(xr), xr_increment=up(xr_inc),
                                             slack_1e12=up(Cx * (arb('1e-12') - arb(box['Ycert'])))),
                box=dict(eta1=up(eta1), S_W=up(SWbox), xbar=up(xbar), thbar=up(thbar), f_x=up(fx), kappa_R=up(kR), sigma=up(sig), sigma0=up(sig0),
                         Z=[up(z) for z in Zs], tip=dict((k, up(v)) for k, v in tipk.items())))


JDIR = os.path.join(S4, 'p2', 't6c', 'jets')        # T6 jets of the (SE)-corrected base (review B1; equal to t6/jets for x <= 40)


def main():
    global JDIR
    ap = argparse.ArgumentParser()
    ap.add_argument('--Ybox', type=float, default=1.7e-12); ap.add_argument('--Ycert', type=float, default=6.3e-13)
    ap.add_argument('--n3', default=os.path.join(S4, 'p2', 't8se', 'constants', 'n3.json'))
    ap.add_argument('--work', default='.', help='directory with constants.json (relative to t7/); output WORK/nl_terms.json')
    ap.add_argument('--out', default=None); ap.add_argument('--ids', default='1-25')
    ap.add_argument('--res', default=os.path.join(S4, 'p2', 't6c', 'residuals.json'), help='T6 residuals of the corrected base (B1)')
    ap.add_argument('--jets', default=JDIR, help='T6 jets of the corrected base (B1)')
    A_ = ap.parse_args(); t0 = time.time()
    WORK = A_.work if os.path.isabs(A_.work) else os.path.normpath(os.path.join(HERE, A_.work))
    if A_.out is None: A_.out = os.path.join(WORK, 'nl_terms.json')
    C = json.load(open(os.path.join(WORK, 'constants.json'))); N3 = json.load(open(A_.n3))
    JDIR = A_.jets
    R6 = {}
    resfiles = [os.path.join(S4, 't6', 'residuals.json'), A_.res]      # max with the polished values (conservative)
    for f in resfiles:
        if os.path.exists(f):
            for r in json.load(open(f))['per_slab']:
                o = R6.get(r['slab'])
                if o is None: R6[r['slab']] = r
                else:
                    o['norms']['R2'] = max(o['norms']['R2'], r['norms']['R2'])
                    for k in ('range_c', 'range_b'):
                        o['scalars'][k] = [min(o['scalars'][k][0], r['scalars'][k][0]), max(o['scalars'][k][1], r['scalars'][k][1])]
    ell_min = min(r['scalars']['range_ell'][0] for r in R6.values())
    MU = U(C['meta']['closure']['mu_e']); c1W = U(C['meta']['weight_consts']['c1W'])
    tb, LB = lemma4b(C, MU, c1W, ell_min)
    print('Lemma 4b: sup|e| <= %.5e, sup|e\'| <= %.5e (x sqrt Lambda_e); eta1 at Ybox %.2e: %.4e (<= 1e-5: %s)' %
          (tb['c_e'], tb['c_e1'], A_.Ybox, tb['c_e1'] * A_.Ybox, tb['c_e1'] * A_.Ybox <= 1e-5), flush=True)
    PHS = phi_sups()
    print('phi^ sups (s0, s1, s2)/rho:', [up(v) for v in PHS[0]], 'int rho', up(PHS[1]), '(%.0fs)' % (time.time() - t0), flush=True)
    a_, b_ = [int(v) for v in A_.ids.split('-')]
    certs = [c for c in C['certificates'] if a_ <= c['id'] <= b_]
    box = dict(Ybox=A_.Ybox, Ycert=A_.Ycert, xnl=arb(0))
    import provenance as PV
    sl_all = sorted({s for c in certs for s in range(c['slabs'][0], c['slabs'][1] + 1) if str(s) in N3['per_slab']})
    here = lambda *p: os.path.join(HERE, *p)
    args = dict(argv=sys.argv[1:], Ybox=A_.Ybox, Ycert=A_.Ycert, ids=A_.ids, work=PV.rel(WORK),
                inputs=dict(constants=PV.entry(os.path.join(WORK, 'constants.json')), n3=PV.entry(A_.n3), n3_args=N3.get('args'),
                            residuals=[PV.entry(f) for f in resfiles if os.path.exists(f)], jets=PV.jets_entry(JDIR, sl_all)),
                code=PV.code([here('nlterms7.py'), here('nlsym', 'nlforms_gen.py'), here('nlsym', 'td2.py'), here('nlsym', 'nlgen.py'),
                              here('nlsym', 'rampsym.py'), here('nlsym', 'rampforms.json'), os.path.join(S4, 'p2', 't8se', 'n3qforms.py'),
                              os.path.join(S4, 't6', 'reader.py'), here('cert7.py'), here('certdata.py'), here('base7.py'), here('provenance.py')]))
    out = dict(Ybox=A_.Ybox, Ycert=A_.Ycert, tip_box=tb, intervals=[], n3=os.path.relpath(A_.n3, S4), args=args)
    for c in certs:
        # (M) fixed point: xbar contains x_NL, which is computed with xbar (two passes; x_NL is tiny)
        box['xnl'] = arb(0); r = interval(c, C, N3, R6, box, LB, PHS)
        if r is None: print(f"interval {c['id']}: no n3 slabs, skipped", flush=True); continue
        box['xnl'] = arb(r['M']['xNL_box']) * arb('1.01'); r2 = interval(c, C, N3, R6, box, LB, PHS)
        assert r2['M']['xNL_box'] <= float(box['xnl'].upper()), (r2['M']['xNL_box'], float(box['xnl'].upper()))
        out['intervals'].append(r2)
        print(f"interval {c['id']:2d} slabs {r2['slabs']}: {len(r2['terms'])} terms, max diss. degree {r2['max_dissipation_degree']:.2f}; "
              f"top I0 {r2['top']['I0_V2sq']:.3e} I2 {r2['top']['I2_V4sq']:.3e} (sup At+ {r2['top']['sup_Atp']:.3e}); x_NL(box) {r2['M']['xNL_box']:.2e} "
              f"(slack {r2['M']['slack_1e12']:.2e}) ({time.time() - t0:.0f}s)", flush=True)
    json.dump(out, open(A_.out, 'w'), indent=1)
    print('wrote', A_.out)


if __name__ == '__main__':
    main()
