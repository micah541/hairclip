"""Base geometry of slab SLAB at a t_c-ball T, from the T6 engine (bivariate (s, t_c) enclosures):
per element e: theta (exact BP), kappa, kappa_x, kappa_xx, kappa_xxx (exact BPs), position Gamma = Pg * G (pointwise CBPs).
Gamma(x=0) = (-2.500806276619802, 0), theta(0) = -pi/2 (T6/T1 convention)."""
import sys, time, math, pickle, os
S4 = '/home/claude/s1spec/step4'
sys.path.insert(0, f'{S4}/t6'); sys.path.insert(0, f'{S4}/t1')
from flint import arb, acb, arb_poly, ctx
from base import SlabBase
import tp as T
from engine import element_a, element_rot
ctx.prec = 128
ONE, ZERO = arb(1), arb(0)

def build(slab, xmax, poldir=f'{S4}/t1/polished_v2', verbose=True):
    t0 = time.time()
    B = SlabBase(slab, D=20, poldir=poldir); ctx.prec = 128
    Z = T.tconst(ZERO)
    st = dict(Qe=Z, rte=Z, the=T.tconst(-arb.pi() / 2), ze=(Z, Z), Ge=(Z, T.tconst(-arb(2.500806276619802))), Pe=(Z, T.tconst(-ONE)))
    els = []
    for e in range(B.N):
        xe0 = B.x[e](arb(0))
        if float(xe0.mid()) > xmax: break
        A_ = element_a(B, e, st['Qe'], st['rte'], st['the'])
        R_ = element_rot(B, e, A_, st['ze'], st['Ge'], st['Pe'], DS=30, DTE=8, DJ=8)
        els.append(dict(e=e, x=B.x[e], h=B.h[e], th=A_['th'], kap=A_['kap'], kx=A_['kx'], kxx=A_['kxx'], kxxx=A_['kxxx'], G=R_['G'], Pg=R_['Pg']))
        st['Qe'], st['rte'], st['the'] = A_['Q'].at(ONE), A_['rt'].at(ONE), A_['th'].at(ONE)
        st['ze'] = R_['z'].at(ONE); st['Ge'] = R_['G'].at(ONE); st['Pe'] = R_['Pg'].at(ONE)
        if verbose and e % 10 == 0: print(f'  element {e}: x_e {float(xe0.mid()):.3f} ({time.time() - t0:.0f}s)', flush=True)
    return B, els

class Geo:
    def __init__(self, slab, xmax, T_):
        self.B, self.els = build(slab, xmax)
        self.T = T_
        self.xe = [el['x'](T_) for el in self.els]; self.hh = [el['h'](T_) for el in self.els]
        self.xe_mid = [float(v.mid()) for v in self.xe]
    def elem(self, x):
        """element index for a float x (by midpoints)"""
        i = 0
        while i + 1 < len(self.xe_mid) and self.xe_mid[i + 1] <= x: i += 1
        return i
    def at(self, i, s):
        """quantities at local s (arb) in element i and t_c = T: theta, Gamma (acb), kappa, kappa_x, kappa_xx, kappa_xxx, x"""
        el = self.els[i]; T_ = self.T
        sh = arb_poly([arb(T_.mid()), arb(1)]); dl = arb(0, T_.rad())
        def tv(tp): return tp(sh)(dl)                     # centred (Taylor-shifted) evaluation on the t_c-ball
        ev = lambda bp: tv(bp.at(s))
        th = ev(el['th']); k0, k1, k2, k3 = ev(el['kap']), ev(el['kx']), ev(el['kxx']), ev(el['kxxx'])
        Gr, Gi = el['G'].at(s); Pr, Pi = el['Pg'].at(s)
        G = acb(tv(Gr), tv(Gi)); P = acb(tv(Pr), tv(Pi))
        x = self.xe[i] + self.hh[i] * (1 + s) / 2
        return dict(th=th, Gam=P * G, k0=k0, k1=k1, k2=k2, k3=k3, x=x, P=P)
