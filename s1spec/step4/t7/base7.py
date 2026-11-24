"""T7 production: the polished base and the fixed constraint functions.

Base on slab i = T6's definition (t6/base.py, copied unchanged as t6base.py/tp6.py, only the import name patched):
  tau = tau0 + dt (1 + tc)/2, tc in [-1, 1]; every stored unknown is the degree-8 polynomial in tc through the stored
  stage balls (Radau IIA nodes as exact dyadics); per element kappa = sum_q d_q(tc) N_q(s), outer knots as Taylor models in tc.
Sharp evaluation at a tc-ball: every quantity is first formed as a polynomial in tc (TP/BP of T6) and only then evaluated
at the ball, so ball radii reflect the true variation over the ball (evaluating the DOFs first loses a factor ~600 through
the cancellations of the Hermite basis).

Constraint functions (fixed, tau-independent; T9 must use exactly these):
  phi_1 = chi_32 p_L, phi_2 = chi_32 n_L, phi_3 = chi_32 psi3^ (c4fun.PSI of step 3).
  p_L = <Gamma_L, nu_L>, n_L = <e1, nu_L> for the curve Gamma_L with curvature kappa_L, theta_L(0) = -pi/2, Gamma_L(0) = -DS1 e1,
  arclength x from the tip; kappa_L = T1's piecewise polynomial on the slab-187 mesh (Hermite-C4 basis, nodal scaling
  (h_e/H_j)^m, H_j = sqrt(h_{j-1} h_j)) with the DOFs = the exact midpoints of the stored stage-8 row of slab 187
  (tc = +1, k = 2.36e4), stored losslessly in phi_store.json.  kappa_L is exactly C^4.
All paths are relative to this file."""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__))
STEP4 = os.path.dirname(HERE)
T1DIR = os.path.join(STEP4, 't1')
POLDIR = os.path.join(T1DIR, 'polished')
for p in (HERE, T1DIR):
    if p not in sys.path: sys.path.insert(0, p)
from flint import arb, arb_poly, ctx
import polstore
from rep import DS1, K0
import t6base as TB
import tp6 as T

PREC = 128
PHI_STORE = os.path.join(HERE, 'phi_store.json')
PHI_SLAB = 187
S_OF_U = arb_poly([-1, 2])          # s = 2u - 1


def bp_at(bp, t):
    """BP (list of s-polys by tc-power) evaluated at tc = t (arb, point or ball) -> arb_poly in s"""
    acc = bp.P[-1]
    for m in range(len(bp.P) - 2, -1, -1): acc = acc * t + bp.P[m]
    return acc


def set_base(path):
    """select the polished store (a directory of slab_XXX.pkl); relative paths are taken from t7/ (e.g. ../t1/polished_v2)"""
    global POLDIR
    p = path if os.path.isabs(path) else os.path.normpath(os.path.join(HERE, path))
    if not os.path.isdir(p): raise FileNotFoundError(p)
    POLDIR = p; TB.POLDIR = p
    return p


def base_arg(ap):
    ap.add_argument('--base', default=os.path.join('..', 't1', 'polished'), help='polished store (relative to t7/ or absolute)')


def store_dir(path=None):
    if path is None: return POLDIR
    return path if os.path.isabs(path) else os.path.normpath(os.path.join(HERE, path))


def manifest(path=None):
    """the store's manifest.json: per-slab records (slab, tau0, tau1, k0, k1, ...), sorted by slab"""
    return sorted(json.load(open(os.path.join(store_dir(path), 'manifest.json'))), key=lambda r: r['slab'])


def slab_range(path=None):
    M = manifest(path)
    return M[0]['slab'], M[-1]['slab']


# certificate intervals: group size by k at the group start (reproduces the v1 plan: slabs 5-18 / 19-58 / 59-128 / 129-187)
K_GROUPS = ((40.0, 2), (185.0, 4), (2500.0, 10), (float('inf'), 20))


def schedule(path=None, first=None, last=None):
    """certificate plan on the store's own slab / tau axis: consecutive slabs [lo, hi] covering slabs first..last (default: the
    whole store), size 2/4/10/20 by k0(lo) (thresholds 40, 185, 2500); a size-20 group absorbs a remainder < 10 at the end;
    tau_c = start of slab slab_c = lo + (hi - lo + 1) // 2.  (v1 store with first = 5 gives the v1 plan.)"""
    M = manifest(path); k0 = {r['slab']: r['k0'] for r in M}
    first = M[0]['slab'] if first is None else first; last = M[-1]['slab'] if last is None else last
    out = []; s = first
    while s <= last:
        size = next(n for kk, n in K_GROUPS if k0[s] < kk)
        hi = min(last, s + size - 1)
        if size == 20 and last - hi < size // 2: hi = last
        out.append((s, hi)); s = hi + 1
    return [dict(id=i, slabs=(lo, hi), slab_c=lo + (hi - lo + 1) // 2) for i, (lo, hi) in enumerate(out)]


class Slab:
    """polished slab i (T6 SlabBase): sharp states at tc points or balls"""
    def __init__(self, i, D=20):
        self.SB = SB = TB.SlabBase(i, D); ctx.prec = PREC
        self.i, self.n, self.N, self.D = i, SB.n, SB.N, D
        self.tau0 = float(SB.tau0); self.dtf = float(SB.d['dt']); self.k0 = float(SB.k)
        self._kbp = {}

    def tau(self, tc): return self.tau0 + self.dtf * (1 + float(tc)) / 2

    def kbp(self, e):
        if e not in self._kbp:
            d, dd = self.SB.local_coeffs(e)
            self._kbp[e] = (TB.expand(d, self.SB.Ns, self.D), TB.expand(dd, self.SB.Ns, self.D))
        return self._kbp[e]

    def base(self, tc, xmax=None, emax=None):
        return BaseData.from_slab(self, arb(tc) if not isinstance(tc, arb) else tc, xmax, emax)


class BaseData:
    """base at one tc (point or ball): per T1 element exact polys (in u in [0,1]) kp, Q, th, kt (kappa_tau at fixed scaling,
    tip elements only exact); scalars a, c = 1/4 - a, b, ell, lnlam, k"""
    @staticmethod
    def from_slab(S, t, xmax, emax):
        B = BaseData(); SB = S.SB
        B.a = SB.a(t); B.c = SB.c(t); B.b = SB.b(t); B.ell = SB.ell(t); B.lnlam = SB.lnlam(t)
        B.ellt = SB.elld(t); B.at = SB.a.derivative()(t) * (2 / SB.dt); B.ct = -B.at; B.bt = SB.b.derivative()(t) * (2 / SB.dt)
        B.k = float((arb(K0) * (-B.lnlam).exp()).mid())
        B.els = []; th = -arb.pi() / 2; Qe = arb(0)
        for e in range(S.N):
            x = SB.x[e](t); h = SB.h[e](t)
            if xmax is not None and float(x.lower()) > float(xmax.upper()): break
            if emax is not None and e >= emax: break
            kb, kdb = S.kbp(e)
            kp = bp_at(kb, t)(S_OF_U); ku = kp.derivative()
            Q = arb_poly([Qe]) + (ku * ku).integral() * (1 / h)
            tpoly = kp.integral() * h + th
            B.els.append(dict(x=x, h=h, kp=kp, Q=Q, th=tpoly, kt=bp_at(kdb, t)(S_OF_U), e=e))
            th = tpoly(arb(1)); Qe = Q(arb(1))
        B.xk = [el['x'] for el in B.els] + [B.els[-1]['x'] + B.els[-1]['h']]
        return B

    @staticmethod
    def from_state(op, Y, Yd, xmax, n):
        """exact (non-ball) state through T1's ArbOp (used for the stored constraint curve)"""
        B = BaseData(); els, gl = op.elements(Y, Yd)
        B.a = gl['a']; B.c = arb(1) / 4 - B.a; B.b = gl['k2t'] + B.c * arb(DS1); B.ell = gl['ell']; B.lnlam = Y[n + 1]
        B.k = float((arb(K0) * (-B.lnlam).exp()).mid())
        B.els = []; th = -arb.pi() / 2
        for el in els:
            if xmax is not None and float(el['x'].lower()) > float(xmax.upper()): break
            tpoly = el['kp'].integral() * el['h'] + th
            B.els.append(dict(x=el['x'], h=el['h'], kp=el['kp'], Q=el['Q'], th=tpoly, kt=el['kt'], e=el['e']))
            th = tpoly(arb(1))
        B.xk = [el['x'] for el in B.els] + [B.els[-1]['x'] + B.els[-1]['h']]
        return B

    def local(self, E, a, h):
        """exact polys in s in [0,1] of kappa^(m) (m <= 4), Q, theta on [a, a+h] inside T1 element E (fixed-knot elements)"""
        el = self.els[E]; ih = 1 / el['h']
        U = arb_poly([(a - el['x']) * ih, h * ih])
        kp = el['kp']; K = []; kt = el['kt']; Kt = []
        for m in range(5):
            K.append(kp(U) * ih ** m); kp = kp.derivative()
            Kt.append(kt(U) * ih ** m); kt = kt.derivative()
        return dict(K=K, Q=el['Q'](U), th=el['th'](U), Kt=Kt)


# ------------------------------------------------------------------------------------------ stored constraint base (phi_1, phi_2)
def make_phi_store():
    SP, d = polstore.load_slab(os.path.join(POLDIR, 'slab_%03d.pkl' % PHI_SLAB)); row = SP.Yn[len(SP.Yn) - 1]
    mids = [polstore.enc(arb(v.mid())) for v in row]
    json.dump(dict(slab=PHI_SLAB, row='Yn[8] (tc = +1)', note='exact midpoints (man, exp, 0, 0) of the stored balls; '
                   'kappa_L = T1 piecewise polynomial on the slab-187 mesh with these DOFs (t7/base7.py late_base)', Y=mids),
              open(PHI_STORE, 'w'))


_LATE = None
def late_base(xmax):
    """BaseData of the stored constraint curve kappa_L (phi_store.json holds the DOFs and the mesh: independent of --base)"""
    global _LATE
    if _LATE is None:
        d = json.load(open(PHI_STORE)); M_ = d['mesh']
        from rep import Rep, Mesh
        from arbops import ArbBasis, ArbOp
        m = Mesh(JT=M_['JT'], h0=M_['h0'], cT=M_['cT'], C=[float.fromhex(v) for v in M_['C']]); R = Rep(M_['P'], m)
        ctx.prec = PREC; op = ArbOp(R, ArbBasis(R.B, PREC))
        _LATE = (op, [polstore.dec(tuple(t)) for t in d['Y']], M_['n'])
    op, Y, n = _LATE; ctx.prec = PREC
    return BaseData.from_state(op, Y, [arb(0)] * len(Y), xmax, n)


if __name__ == '__main__':
    import time
    if not os.path.exists(PHI_STORE): make_phi_store(); print('wrote', PHI_STORE)
    t0 = time.time(); S = Slab(slab_range()[0]); print('slab load', time.time() - t0)
    t0 = time.time(); Bb = S.base(arb(0, 1)); print('ball base', time.time() - t0, 'a', Bb.a, 'b', Bb.b, 'ell', Bb.ell, 'nel', len(Bb.els))
    Bm = S.base(arb(-1)); Bp = S.base(arb(1))
    for e in (0, 10, 30, 45, 80, 110):
        kb = Bb.els[e]['kp']; rb = max(float(c.rad()) for c in kb.coeffs())
        tv = max(abs(float((a - b).mid())) for a, b in zip(Bm.els[e]['kp'].coeffs(), Bp.els[e]['kp'].coeffs()))
        print(e, 'x', Bb.els[e]['x'], 'ball rad', rb, 'true var', tv, 'ratio', rb / max(tv, 1e-300))
    L = late_base(arb(40)); print('late k', L.k, 'a', L.a, 'nel', len(L.els))
