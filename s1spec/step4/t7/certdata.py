"""Rigorous (arb) data of the projected form on a FEM mesh of [0, XE], XE = 32 - STIP (end of supp chi_32):
per FEM element: exact polynomials of the base coefficients (b2, b1, b0 and derivatives), Taylor models (tm.TM) of the
w-form coefficients C_ab, of e^{-phi/2} chi g_i (base directions at tau_c), e^{-phi/2} phi_j, e^{-phi/2} ell_j (constraints, strong form)
and of e^{-phi} chi g_i phi_j.  Constraint functions (fixed, tau-independent):
   phi_1 = chi p~_L, phi_2 = chi n~_L  (p~, n~ of the T1 base at its last node, k = 2.36e4: S1 surrogates, explicit)
   phi_3 = chi psi3~  (C^4 slow eigenfunction of step 3, c4fun.py: exact polynomial pieces on [0,28], explicit extension beyond)."""
import sys, numpy as np
from math import comb, factorial
from flint import arb, arb_poly, arb_series, ctx
ctx.prec = 128
import os
_HERE = os.path.dirname(os.path.abspath(__file__)); _STEP4 = os.path.dirname(_HERE)
for _p in (os.path.join(_STEP4, 't1'), os.path.join(os.path.dirname(_STEP4), 'step3nl', 'agent_int', 'v3'), _HERE):
    if _p not in sys.path: sys.path.insert(0, _p)
import tm
from tm import TM, from_series, poly_eval_series, HALF, T2S
from arbops import ArbBasis, ArbOp, A as Aexact
from rep import K0, DS1, STIP
import c4fun

STIPA = arb(STIP); X1 = arb(29) - STIPA; XE = arb(32) - STIPA
D_SER = 24
PSI = c4fun.PSI; PM = PSI.M

# ------------------------------------------------------------------ base in arb (per T1 element exact polynomials in u)
class ArbBase:
    def __init__(self, i, theta, xmax):
        from base2 import Base                         # prototype only (float T1 base)
        B = Base('P11'); R, D = B.rep(i); tau = B.tau0[i] + theta * B.dt[i]; Y = B.Y(tau, i)
        self.k = float(B.k(tau)); self.i, self.theta = i, theta
        AB = ArbBasis(R.B, ctx.prec); op = ArbOp(R, AB)
        Ya = [arb(float(v)) for v in Y]; Yd = [arb(0)] * len(Y)
        els, gl = op.elements(Ya, Yd)
        self.a = gl['a']; self.c = arb(1) / 4 - self.a; self.b = gl['k2t'] + self.c * arb(DS1)
        self.els = []; th = -arb.pi() / 2
        for el in els:
            if float(el['x'].lower()) > float(xmax.upper()): break
            kp = el['kp']; tpoly = kp.integral() * el['h'] + th
            self.els.append(dict(x=el['x'], h=el['h'], kp=kp, Q=el['Q'], th=tpoly))
            th = tpoly(arb(1))
        self.xk = [e['x'] for e in self.els] + [self.els[-1]['x'] + self.els[-1]['h']]
    def local(self, E, a, h):
        """exact polys in s of kappa^(m) (m<=4), Q, theta on [a, a+h] inside T1 element E"""
        el = self.els[E]; ih = 1 / el['h']
        U = arb_poly([(a - el['x']) * ih, h * ih])
        kp = el['kp']; K = []
        for m in range(5):
            K.append(kp(U) * ih ** m); kp = kp.derivative()
        return dict(K=K, Q=el['Q'](U), th=el['th'](U))

def coeff_polys(L, ab, a, h):
    """b2, b2', b2'', b1, b1', b0 as exact polys in s (T5 linearisation; b1 = -6 k k' + Q - (1/4 - a) x)"""
    k0, k1, k2 = L['K'][0], L['K'][1], L['K'][2]; c = ab.c; x = arb_poly([a, h])
    b2 = -(k0 * k0); b2p = -2 * k0 * k1; b2pp = -2 * k1 * k1 - 2 * k0 * k2
    b1 = -6 * k0 * k1 + L['Q'] - x * c; b1p = -5 * k1 * k1 - 6 * k0 * k2 - arb_poly([c])
    b0 = arb_poly([arb(1) / 4 - ab.a]) - 3 * k1 * k1 - 4 * k0 * k2
    return b2, b2p, b2pp, b1, b1p, b0

def trig_tms(Lth, a, h):
    """TMs of cos theta, sin theta (theta exact poly in s)"""
    Tt = Lth(arb_poly([HALF, 1]))                    # theta(1/2 + t)
    xm = a + h * HALF
    mk = lambda f: (lambda X: getattr(poly_eval_series(Tt, (X - xm) * (1 / h), D_SER + 1), f)())
    return from_series(mk('cos'), a, h, D_SER)[0], from_series(mk('sin'), a, h, D_SER)[0]

class Geo:
    """running geometry (Gamma accumulated from the tip) of one base along the FEM mesh"""
    def __init__(self, ab):
        self.ab = ab; self.Gr = -arb(DS1); self.Gi = arb(0)
    def element(self, E, a, h, nder=0):
        L = self.ab.local(E, a, h); cs, sn = trig_tms(L['th'], a, h)
        Gr = cs.antideriv(h) + self.Gr; Gi = sn.antideriv(h) + self.Gi
        self.Gr = self.Gr + h * cs.integral01(); self.Gi = self.Gi + h * sn.integral01()
        p = Gi * cs - Gr * sn; q = Gr * cs + Gi * sn; n = -1 * sn; t1 = cs
        out = dict(L=L, cos=cs, sin=sn, p=p, q=q, n=n, t1=t1)
        if nder:
            K = [TM(k) for k in L['K']]
            P = [p]; Qd = [q]; N = [n]; T = [t1]
            for m in range(1, nder + 1):          # p' = -k q, q' = 1 + k p, n' = -k t1, t1' = k n  (Leibniz)
                P.append(sum((-comb(m - 1, j)) * K[j] * Qd[m - 1 - j] for j in range(m)))
                Qd.append((1 if m == 1 else 0) + sum(comb(m - 1, j) * K[j] * P[m - 1 - j] for j in range(m)))
                N.append(sum((-comb(m - 1, j)) * K[j] * T[m - 1 - j] for j in range(m)))
                T.append(sum(comb(m - 1, j) * K[j] * N[m - 1 - j] for j in range(m)))
            out.update(P=P, N=N)
        return out

# ------------------------------------------------------------------ fixed functions
def chi_jets(a, h, nder=4):
    """chi_32 and x-derivatives as exact polys in s (1 for x <= X1)"""
    if float(a.upper()) <= float(X1.lower()) and float((a + h).upper()) <= float(X1.upper()) + 1e-30:
        return [arb_poly([1])] + [arb_poly([0])] * nder
    S = arb_poly([0, 0, 0, 0, 35, -84, 70, -20]); U = arb_poly([(XE - a) / 3, -h / 3])
    out = []; Sd = S
    for m in range(nder + 1):
        out.append(Sd(U) * (-arb(1) / 3) ** m); Sd = Sd.derivative()
    return out

def psi3_piece(e, a, h, n=5):
    """exact polys in s of psi3~ and derivatives on [a, a+h] inside c4fun piece e (x <= 28)"""
    P = PSI.pieces[e]; X = arb_poly([a, h])
    yy = (X - P.a) * P.sc + 1
    raw = [c4fun.F3.clen(P.cfs[d], yy) * P.sc ** d for d in range(n)]
    Lp = P.b - P.a; hv = (X - P.a) * (1 / Lp); out = []
    for k in range(n):
        v = arb_poly([0])
        for d in range(c4fun.NJET):
            for side, H in (('L', PSI.Hl[e][d]), ('R', PSI.Hr[e][d])):
                if H == 0: continue
                cf = c4fun.HB[(side, d)]; m = arb_poly([0])
                for ii in range(k, len(cf)):
                    f = 1
                    for t in range(k): f *= (ii - t)
                    if f: m = m + (arb(cf[ii].p * f) / arb(cf[ii].q)) * hv ** (ii - k)
                v = v + m * (H * Lp ** d)
        out.append(raw[k] + v * (1 / Lp ** k))
    return [TM(q) for q in out]

def psi3_ext(a, h, n=5):
    """TMs of the extension c1 g(s) + c2 E(s) (x >= 28) and derivatives"""
    E0 = 3 * ((PM.XC + PM.STIP) / 4) ** (arb(4) / 3)
    def F(X):
        S = X + PM.STIP; v = 0 * S
        for j in range(4): v = v + PM.AJ[j] * S ** (PM.A - 4 * j)
        return PM.C1 * v + PM.C2 * (-3 * (S * (1 / arb(4))) ** (arb(4) / 3) + E0).exp()
    return from_series(F, a, h, D_SER, nder=n - 1)

def phi_tms(a, h):
    """phi' .. phi'''' (TMs), e^{-phi/2}, e^{-phi}"""
    F = lambda X: (1 + X * X) ** (arb(2) / 3) / 4
    d = from_series(F, a, h, D_SER, nder=4)[1:]
    em = from_series(lambda X: (-(1 + X * X) ** (arb(2) / 3) / 8).exp(), a, h, D_SER)[0]
    e1 = from_series(lambda X: (-(1 + X * X) ** (arb(2) / 3) / 4).exp(), a, h, D_SER)[0]
    return d, em, e1

def leib(c, f, m):
    return sum(comb(m, k) * c[k] * f[m - k] for k in range(m + 1))

def ell_strong(f, b, d):
    """rho^-1 L^dag(rho f) from TMs f = (f..f''''), b = (b2,b2',b2'',b1,b1',b0) polys, d = (phi'..phi'''')"""
    p1, p2, p3, p4 = d
    r1 = -1 * p1; r2 = p1 * p1 - p2; r3 = -1 * (p1 * p1 * p1) + 3 * (p1 * p2) - p3
    r4 = p1 * p1 * p1 * p1 - 6 * (p1 * p1 * p2) + 3 * (p2 * p2) + 4 * (p1 * p3) - p4
    f0, f1, f2, f3, f4 = f; b2, b2p, b2pp, b1, b1p, b0 = [TM(v) for v in b]
    t4 = f4 + 4 * (r1 * f3) + 6 * (r2 * f2) + 4 * (r3 * f1) + r4 * f0
    h0 = b2 * f0; h1 = b2p * f0 + b2 * f1; h2 = b2pp * f0 + 2 * (b2p * f1) + b2 * f2
    t2 = h2 + 2 * (r1 * h1) + r2 * h0
    t1 = (b1p * f0 + b1 * f1) + r1 * (b1 * f0)
    return -1 * t4 + t2 - t1 + b0 * f0

# ------------------------------------------------------------------ mesh
def mesh(ab, htarget):
    bp = [arb(0)] + [x for x in ab.xk if 0 < float(x.mid()) < float(XE.mid())]
    bp += [p.a for p in PSI.pieces if 0 < float(p.a.mid())] + [PM.XC, X1, XE]
    bp = sorted(bp, key=lambda v: float(v.mid())); pts = [bp[0]]
    for v in bp[1:]:
        if float(v.mid()) - float(pts[-1].mid()) > 1e-9: pts.append(v)
    nodes = [pts[0]]
    for p, q in zip(pts[:-1], pts[1:]):
        L = float((q - p).mid()); m = max(1, int(np.ceil(L / htarget - 1e-9)))
        for k in range(1, m): nodes.append(arb(float((p + (q - p) * k / m).mid())))
        nodes.append(q)
    return nodes

def t1_index(ab, x):
    xf = float(x.mid()); xs = [float(v.mid()) for v in ab.xk]
    return max(0, min(int(np.searchsorted(xs, xf + 1e-12, side='right') - 1), len(ab.els) - 1))

def psi_index(x):
    xf = float(x.mid())
    for e, P in enumerate(PSI.pieces):
        if float(P.a.mid()) - 1e-12 <= xf < float(P.b.mid()) - 1e-12: return e
    return None
