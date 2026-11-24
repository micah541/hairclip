"""T6 layer 2: the stored polished base of one slab as exact objects in arb.

Definition of the base on slab i (used by every T6 bound; T7/T8 should use the same):
  * time: tau = tau0 + dtau (1 + tc)/2, tc in [-1, 1], dtau = the stored float dt (exact dyadic, as in T1's polish);
  * every unknown (DOF, ell, ln lam, c, t) is the degree-8 polynomial in tc interpolating the stored stage balls Yn[j]
    at the nodes nu_0 = 0, nu_1..nu_8 = the Radau IIA nodes of radau.py (60-digit mpmath binary values, taken as exact
    dyadic rationals; nu_8 = 1 exactly, so slabs join exactly);
  * per element e the curvature is sum_q d_q(tc) N_q(s), s = 2u - 1, with T1's normalised Hermite-C4 + bubble basis,
    d_q = y_q (h_e/H_j)^m (T1's nodal scaling), knots as in rep.Mesh / arbops.ArbOp (outer knots XT (1 + C ell^3)^(1/3)).
Outer knots depend on ell(tc) algebraically; they are carried as Taylor models in tc (tp.tpow)."""
import os, sys
from fractions import Fraction
from flint import arb, arb_poly, fmpq, fmpq_poly, ctx
HERE = os.path.dirname(os.path.abspath(__file__))
T1DIR = os.path.normpath(os.path.join(HERE, '..', 't1'))
if T1DIR not in sys.path: sys.path.insert(0, T1DIR)
from rep import Rep, Mesh, K0, DS1, basis
from arbops import hermite_exact
from radau import Radau
from polstore import load_slab
import tp as T

PREC = 128
POLDIR = os.path.join(T1DIR, 'polished')
_RADAU = {}


def _mpf_to_fmpq(v):
    s, man, ex, bc = v._mpf_
    q = fmpq(int(man)) * (fmpq(2) ** int(ex)) if ex >= 0 else fmpq(int(man), 2 ** (-int(ex)))
    return -q if s else q


def radau_nodes(s, RD=None):
    if s not in _RADAU:
        RD = Radau(s) if RD is None else RD
        _RADAU[s] = [fmpq(0)] + [_mpf_to_fmpq(c) for c in RD.cm]
    return _RADAU[s]


def lagrange_tc(nodes):
    """Lagrange basis on tc_i = 2 nu_i - 1 as exact fmpq_polys in tc"""
    tcs = [2 * v - 1 for v in nodes]
    out = []
    for j, tj in enumerate(tcs):
        p = fmpq_poly([1])
        for k, tk in enumerate(tcs):
            if k != j: p = p * fmpq_poly([-tk, 1]) * (1 / (tj - tk))
        out.append(p)
    return out


def q2a(q):
    return arb(q) if isinstance(q, fmpq) else arb(fmpq(q.numerator, q.denominator))


def dec(t):
    m, e, r, re = t
    if r == 0: return arb(mid=arb((m, e)))
    k = max(0, 32 - r.bit_length())
    return arb(mid=arb((m, e)), rad=arb(((r << k) - 1, re - k)))


_SBASIS = {}
def sbasis(P):
    """normalised local basis N_q(s), s = 2u - 1, exact rationals -> arb_polys; plus nrm (dyadic Fractions)"""
    if P in _SBASIS: return _SBASIS[P]
    FB = basis(P)
    nrm = [Fraction(float(v)) for v in FB.nrm]
    cols = hermite_exact(P)
    u_of_s = fmpq_poly([fmpq(1, 2), fmpq(1, 2)])
    Ns = []
    for j in range(P + 1):
        pu = fmpq_poly([fmpq(c.numerator, c.denominator) / fmpq(nrm[j].numerator, nrm[j].denominator) for c in cols[j]])
        ps = pu(u_of_s)
        Ns.append(arb_poly([arb(c) for c in ps.coeffs()]))
    _SBASIS[P] = (Ns, nrm)
    return _SBASIS[P]


class SlabBase:
    def __init__(self, i, D=20, poldir=None):
        SP, d = load_slab(os.path.join(poldir or POLDIR, f'slab_{i:03d}.pkl'))   # T1's loader (sets ctx.prec = 128)
        ctx.prec = PREC
        assert d['prec'] == PREC
        self.d = d; self.i = i; self.D = D
        cfg = d['cfg']
        self.R = R = SP.R
        self.mesh = m = R.mesh
        self.n, self.N, self.P = R.n, R.N, cfg['P']
        self.tau0, self.k = d['tau0'], d['k']
        self.dt = arb(float(d['dt']))
        Yn = SP.Yn
        self.Yn = Yn
        nodes = radau_nodes(cfg['s'], SP.RD)
        Ls = [arb_poly([arb(c) for c in p.coeffs()]) for p in lagrange_tc(nodes)]
        ntot = len(Yn[0])
        self.Y = []
        for k in range(ntot):
            acc = arb_poly([0])
            for j, L in enumerate(Ls): acc = acc + L * Yn[j][k]
            self.Y.append(acc)
        two_dt = 2 / self.dt
        self.Yd = [p.derivative() * two_dt for p in self.Y]      # d/dtau
        self.Ns, self.nrm = sbasis(self.P)
        self._mesh_tms()
        self._scalars()

    # -------------------------------------------------------------------------------------------- mesh TMs
    def _mesh_tms(self):
        D = self.D; m = self.mesh; n = self.n
        ell = self.Y[n]
        self.ell, self.elld = ell, self.Yd[n]
        XT = arb(m.XT)
        L2 = ell * ell; L3 = L2 * ell
        xs, dxs = [], []
        for v in m.xtip:
            xs.append(T.tconst(arb(float(v)))); dxs.append(T.tconst(arb(0)))
        third = arb(1) / 3
        for Cj in m.C:
            Ca = arb(float(Cj))
            g = T.tconst(arb(1)) + L3 * Ca
            xs.append(T.tpow(g, third, D) * XT)
            dxs.append(T.tmul(L2 * (XT * Ca), T.tpow(g, -2 * third, D), D))
        xs.append(ell); dxs.append(T.tconst(arb(1)))
        N = self.N
        assert len(xs) == N + 1
        h = [xs[e + 1] - xs[e] for e in range(N)]
        dh = [dxs[e + 1] - dxs[e] for e in range(N)]
        ih = [T.tpow(he, -1, D) for he in h]
        g = [T.tmul(dh[e], ih[e], D) for e in range(N)]
        # nodal scale ratios h_e/H_j (left node j=e, right node j=e+1)
        half = arb(1) / 2
        rl, rr = [], []
        for e in range(N):
            if e == 0: rl.append(T.tconst(arb(1)))
            else: rl.append(T.tpow(T.tmul(h[e], ih[e - 1], D), half, D))          # h_e / sqrt(h_{e-1} h_e)
            if e == N - 1: rr.append(T.tconst(arb(1)))
            else: rr.append(T.tpow(T.tmul(h[e], ih[e + 1], D), half, D))
        dHH = [g[0]] + [(g[j - 1] + g[j]) * half for j in range(1, N)] + [g[N - 1]]
        self.x, self.dx, self.h, self.dh, self.ih, self.gl = xs, dxs, h, dh, ih, g
        self.rl, self.rr, self.dHH = rl, rr, dHH

    def _scalars(self):
        n = self.n; R = self.R; nrm = self.nrm
        h0 = self.h[0]            # fixed (tip) element: constant TP
        ih0 = self.ih[0]
        i02, i04 = R.nidx[0, 2], R.nidx[0, 4]
        ih0c = ih0.coeffs()[0]
        self.k2t = self.Y[i02] * (ih0c ** 2 / q2a(nrm[2]))
        self.k4t = self.Y[i04] * (ih0c ** 4 / q2a(nrm[4]))
        K0a = arb(K0)
        self.a = T.tconst(arb(1) / 4) + (self.k4t + self.k2t * K0a ** 2) * (1 / K0a)
        self.c = T.tconst(arb(1) / 4) - self.a
        self.b = self.k2t + self.c * arb(DS1)
        self.lnlam = self.Y[n + 1]

    # -------------------------------------------------------------------------------------------- per element data
    def local_coeffs(self, e):
        """lists d_q(tc), dd_q(tc) (local coefficients and their tau-derivatives at fixed scaling), node index, m"""
        D = self.D; R = self.R
        nb = self.P + 1
        d, dd = [], []
        pw_l = [T.tconst(arb(1))]; pw_r = [T.tconst(arb(1))]
        for mm in range(1, 5):
            pw_l.append(T.tmul(pw_l[-1], self.rl[e], D)); pw_r.append(T.tmul(pw_r[-1], self.rr[e], D))
        for q in range(nb):
            g = R.L[e, q]
            if q < 5: sc = pw_l[q]
            elif q < 10: sc = pw_r[q - 5]
            else: sc = None
            if g >= 0:
                y, yd = self.Y[g], self.Yd[g]
            else:
                y, yd = T.tconst(arb(float(R.fixval[e, q]))), T.tconst(arb(0))
            if sc is not None:
                d.append(T.tmul(y, sc, D)); dd.append(T.tmul(yd, sc, D))
            else:
                d.append(y); dd.append(yd)
        return d, dd


def expand(tps, Ns, D):
    """sum_q tps[q](tc) N_q(s) -> BP"""
    deg = min(D, max(len(p.coeffs()) for p in tps) - 1)
    out = []
    cfs = [p.coeffs() for p in tps]
    for j in range(deg + 1):
        acc = arb_poly([0])
        for q, c in enumerate(cfs):
            if j < len(c): acc = acc + Ns[q] * c[j]
        out.append(acc)
    bp = T.BP(out, True)
    # fold degrees > D (only if some tp is longer)
    extra = [j for j in range(deg + 1, max(len(c) for c in cfs))]
    if extra:
        hi = []
        for j in extra:
            acc = arb_poly([0])
            for q, c in enumerate(cfs):
                if j < len(c): acc = acc + Ns[q] * c[j]
            hi.append(acc)
        bp.P[0] = T._fold(bp.P[0], hi)
    return bp
