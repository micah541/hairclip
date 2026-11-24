"""Seed of the T1 base from P1's exact 128-bit centre (t3/p1/ckpt/slab_NNNNN.json), entirely in arb (no float rounding).
P1: theta(alpha) = -pi/2 + alpha + sum_n u_n sin(n alpha) (sine class, exact dyadic u_n), X(alpha) = X0 + (L/2pi) int_0^alpha e^{i theta};
tip at alpha = pi, far symmetry point X0 at alpha = 0.  Base frame: Gamma~ = (X - c)/lam (translation + scaling, no reflection),
x = arclength from the tip, alpha(x) = pi - c_L lam x, kappa~^(m)(x) = (-1)^(m+1) (lam c_L)^(m+1) theta^(m+1)(alpha).
Projection onto the C^4 space: nodal kappa~^(0..4) exact at the knots; bubbles by Legendre moments of the remainder
(so int_e kappa_h = int_e kappa on every element: turning exactly -pi).  Far point put on the axis by a zero-mean C^inf bump
(physical s ~ 12) solved in arb.  usage: python3 seed_p1.py SLAB [NG] -> seed_p1_SLAB.pkl"""
import sys, json, pickle, time, numpy as np
from flint import arb, acb, arb_poly, acb_poly, arb_mat, ctx
from fractions import Fraction
from rep import Rep, Mesh, K0, DS1
from ops import Dyn
from arbops import ArbBasis, ArbOp, A, legendre_shift
from polstore import enc

PREC = 192

def ldarb(v):                       # a2j format [man, exp, rad_man, rad_exp] -> midpoint (exact)
    return arb(int(v[0])) * arb(2) ** int(v[1])

class P1Curve:
    """exact P1 centre; theta^(m)(alpha) evaluated in blocks: sum_n c_n e^{i n alpha} = sum_q e^{i qB alpha} sum_r c_{qB+r} e^{i r alpha}
    (every e^{i k alpha} from its own sin_cos, inner sums as arb_mat products): no long recurrences, so no wrapping blow-up
    (acb Horner on |z| = 1 inflates radii by (|cos|+|sin|)^n, i.e. 1e152 at n = 2800)."""
    B = 64
    def __init__(self, path):
        rec = json.load(open(path)); self.rec = rec
        self.u = [arb(int(m)) * arb(2) ** int(e) for m, e in rec['xbar']]; self.M = len(self.u)
        nz = max(i for i, (m, e) in enumerate(rec['xbar']) if int(m) != 0) + 1          # last nonzero mode
        self.L = ldarb(rec['L']); self.X0 = ldarb(rec['X0']); self.t = ldarb(rec['t'])
        self.cL = 2 * arb.pi() / self.L
        B = self.B; self.Q = Q = (nz + 1 + B - 1) // B
        cf = [arb(0)] + self.u[:nz] + [arb(0)] * (Q * B - nz - 1)                       # c_n, n = 0 .. QB-1
        self.Cm = [arb_mat(Q, B, [arb(n) ** m * cf[n] for n in range(Q * B)]) for m in range(7)]

    def theta_batch(self, als, ms, chunk=1500):
        """{m: [theta^(m)(alpha) for alpha in als]} (arb, rigorous)"""
        B, Q = self.B, self.Q; out = {m: [] for m in ms}; pi = arb.pi()
        for i0 in range(0, len(als), chunk):
            al = als[i0:i0 + chunk]; P = len(al)
            sr = [[None] * P for _ in range(B)]; cr = [[None] * P for _ in range(B)]
            sq = [[None] * P for _ in range(Q)]; cq = [[None] * P for _ in range(Q)]
            for p, a in enumerate(al):
                for r in range(B): sr[r][p], cr[r][p] = (r * a).sin_cos()
                for q in range(Q): sq[q][p], cq[q][p] = (q * B * a).sin_cos()
            Cr = arb_mat(B, P, [v for row in cr for v in row]); Sr = arb_mat(B, P, [v for row in sr for v in row])
            for m in ms:
                Ic = self.Cm[m] * Cr; Is = self.Cm[m] * Sr
                for p in range(P):
                    if m % 2:          # cosine sum  sum_q cq Ic - sq Is
                        v = sum((cq[q][p] * Ic[q, p] - sq[q][p] * Is[q, p] for q in range(Q)), arb(0))
                    else:              # sine sum    sum_q sq Ic + cq Is
                        v = sum((sq[q][p] * Ic[q, p] + cq[q][p] * Is[q, p] for q in range(Q)), arb(0))
                    v = v if m % 4 < 2 else -v
                    if m == 0: v = v + al[p] - pi / 2
                    if m == 1: v = v + 1
                    out[m].append(v)
        return out

    def theta_der(self, al, m):
        return self.theta_batch([al], [m])[m][0]

    def kap_batch(self, xs, ms, lam):
        """{m: kappa~^(m)(x)} = (-1)^(m+1) (lam c_L)^(m+1) theta^(m+1)(pi - c_L lam x)"""
        s = lam * self.cL; pi = arb.pi()
        th = self.theta_batch([pi - s * x for x in xs], [m + 1 for m in ms])
        return {m: [(-1) ** (m + 1) * s ** (m + 1) * v for v in th[m + 1]] for m in ms}

def main():
    slab = int(sys.argv[1]); NG = int(sys.argv[2]) if len(sys.argv) > 2 else 80
    ctx.prec = PREC; t0 = time.time()
    C1 = P1Curve(f'/home/claude/s1spec/step4/t3/p1/ckpt/slab_{slab:05d}.json')
    pi = arb.pi(); cL = C1.cL
    kt = -cL * C1.theta_der(pi, 1)                    # physical tip curvature (strand orientation)
    lam = A(K0) / kt
    ctx.prec = 128; ell = (+(C1.L / (2 * lam))).mid(); ctx.prec = PREC     # stored length: exact 128-bit dyadic (radius 0); alpha(ell) = O(2^-128)
    print(f'P1 slab {slab}: M = {C1.M}, Q = {C1.Q}, k = {kt.str(20)}, lam = {float(lam.mid()):.6e}, ell = {float(ell.mid()):.6f}  ({time.time()-t0:.1f}s)', flush=True)
    P = 11; mesh = Mesh(JT=40, h0=0.35, cT=0.08, NO=80, c=0.12, beta=0.02, ell0=float(ell.mid()))
    R = Rep(P, mesh); D = Dyn(R); AB = ArbBasis(R.B, PREC); op = ArbOp(R, AB)
    x, dx, h, dh, H, dHH = op.geometry(ell)
    N = R.N; nr = AB.nrm
    y = [arb(0)] * (R.n + 4)
    kn = C1.kap_batch(x, range(5), lam); Dn = [[kn[m][j] for m in range(5)] for j in range(N + 1)]
    for (j, m), v in R.fixed.items():
        err = float(abs(Dn[j][m] - A(v)).mid()); assert err < 1e-30, (j, m, err)      # sine class: exact symmetry, kappa(0) = K0
    for j in range(N + 1):
        for m in range(5):
            g = R.nidx[j, m]
            if g >= 0: y[g] = Dn[j][m] * H[j] ** m * A(nr[m])
    print(f'  N = {N}, n = {R.n}; nodal values done ({time.time()-t0:.1f}s)', flush=True)
    gr = [arb.legendre_p_root(NG, k, weight=True) for k in range(NG)]
    gu = [(g + 1) / 2 for g, w in gr]; gw = [w / 2 for g, w in gr]
    nb = AB.nbub; Pm = [arb_poly([A(c) for c in legendre_shift(m)]) for m in range(nb)]
    Pmu = [[Pm[m](u) for u in gu] for m in range(nb)]
    Binv = arb_mat(nb, nb, [(AB.Nq[10 + i] * Pm[m]).integral()(arb(1)) for m in range(nb) for i in range(nb)]).inv()
    kq_all = C1.kap_batch([x[e] + h[e] * u for e in range(N) for u in gu], [0], lam)[0]
    print(f'  {len(kq_all)} Gauss values done ({time.time()-t0:.1f}s)', flush=True)
    turn = arb(0)
    for e in range(N):
        d = op.local(y, h, H, e)
        for i in range(nb): d[10 + i] = arb(0)
        herm = op.poly(d); kq = kq_all[e * NG:(e + 1) * NG]
        mom = [sum((w * k * pu for w, k, pu in zip(gw, kq, Pmu[m])), arb(0)) - (herm * Pm[m]).integral()(arb(1)) for m in range(nb)]
        b = (Binv * arb_mat(nb, 1, mom)).entries()
        for i, g in enumerate(R.bidx[e]): y[g] = b[i]
        turn += h[e] * sum((w * k for w, k in zip(gw, kq)), arb(0))
    ctx.prec = 128; lnlam = (+lam.log()).mid(); ctx.prec = PREC
    y[R.n] = ell; y[R.n + 1] = lnlam; y[R.n + 2] = arb(0); y[R.n + 3] = C1.t
    seed = dict(slab=slab, P=P, C=mesh.C.copy(), cfg=dict(P=P, JT=40, h0=0.35, cT=0.08, NO=80, c=0.12, beta=0.02, rg=1.5, s=8, dt=0.15, XT=40.0,
                seed=f'P1 slab {slab}', prec=PREC, NG=NG), k=float(kt.mid()), lam=lam, L=C1.L, X0=C1.X0, t=C1.t)
    return C1, R, D, AB, op, y, seed, turn, NG, t0, lam

def closure(R, AB, op, y, NGc=48):
    """arb: int kappa_h (turning), Gamma~(ell) (far point) for state y"""
    ell = y[R.n]; x, dx, h, dh, H, dHH = op.geometry(ell)
    gr = [arb.legendre_p_root(NGc, k, weight=True) for k in range(NGc)]; gu = [(g + 1) / 2 for g, w in gr]; gw = [w / 2 for g, w in gr]
    th = -arb.pi() / 2; G = acb(-A(DS1)); thv = []
    for e in range(R.N):
        kp = op.poly(op.local(y, h, H, e)); tp = arb_poly([th]) + kp.integral() * h[e]
        for u, w in zip(gu, gw):
            tu = tp(u); G += h[e] * w * acb(tu.cos(), tu.sin())
        th = tp(arb(1))
    return th + arb.pi() / 2, G

if __name__ == '__main__':
    C1, R, D, AB, op, y, seed, turn, NG, t0, lam = main()
    tq, G = closure(R, AB, op, y)
    print(f'  projected: turning + pi = {(tq + arb.pi()).str(3, radius=True)} (exact-kappa quadrature: {(turn + arb.pi()).str(3, radius=True)}),'
          f' far point Gamma(ell) = {float(G.real.mid()):.9f} + {G.imag.str(3, radius=True)} i  ({time.time()-t0:.1f}s)', flush=True)
    # far point on the axis: kappa_h += eps * w'(x), w = exp(-((x - xc)/sg)^2), zero mean to e^-64, at physical s = 12 (width 1.5)
    xc = A(12.0) / lam; sg = A(1.5) / lam
    def wder(xv, m):
        t = (xv - xc) / sg; Hm = [arb(1), 2 * t]
        for k in range(1, m + 1): Hm.append(2 * t * Hm[-1] - 2 * k * Hm[-2])
        return (-1) ** m * Hm[m] * (-(t * t)).exp() / sg ** m
    x, dx, h, dh, H, dHH = op.geometry(y[R.n])
    gr = [arb.legendre_p_root(NG, k, weight=True) for k in range(NG)]; gu = [(g + 1) / 2 for g, w in gr]; gw = [w / 2 for g, w in gr]
    nb = AB.nbub; Pm = [arb_poly([A(c) for c in legendre_shift(m)]) for m in range(nb)]
    Binv = arb_mat(nb, nb, [(AB.Nq[10 + i] * Pm[m]).integral()(arb(1)) for m in range(nb) for i in range(nb)]).inv()
    def bump_vec(eps):
        """projection of eps * w' (same projection as the seed) as a state increment"""
        dv = [arb(0)] * len(y)
        for j in range(R.N + 1):
            for m in range(5):
                g = R.nidx[j, m]
                if g >= 0: dv[g] = eps * wder(x[j], m + 1) * H[j] ** m * A(AB.nrm[m])
        for e in range(R.N):
            if float(abs(((x[e] - xc) / sg).mid())) > 12 and float(abs(((x[e + 1] - xc) / sg).mid())) > 12: continue
            d = op.local(dv, h, H, e)
            for q in range(AB.nb):
                if op.L[e, q] < 0: d[q] = arb(0)          # increment: fixed DOF (kappa(0) = K0, ...) contribute 0, not their value
            for i in range(nb): d[10 + i] = arb(0)
            herm = op.poly(d); kq = [eps * wder(x[e] + h[e] * u, 1) for u in gu]
            mom = [sum((w * k * Pm[m](u) for w, k, u in zip(gw, kq, gu)), arb(0)) - (herm * Pm[m]).integral()(arb(1)) for m in range(nb)]
            b = (Binv * arb_mat(nb, 1, mom)).entries()
            for i, g in enumerate(R.bidx[e]): dv[g] = b[i]
        return dv
    one = bump_vec(arb(1))
    _, G1 = closure(R, AB, op, [a + arb('1e-10') * b for a, b in zip(y, one)])
    slope = (G1.imag - G.imag) / arb('1e-10'); eps = -G.imag / slope
    for it in range(2):
        y1 = [a + eps * b for a, b in zip(y, one)]
        tq1, G1 = closure(R, AB, op, y1)
        print(f'  closure it {it}: eps = {float(eps.mid()):.3e}  ImGamma(ell) = {G1.imag.str(3, radius=True)}  turning+pi = {(tq1 + arb.pi()).str(3, radius=True)}', flush=True)
        eps = eps - G1.imag / slope
    y = y1
    ctx.prec = 128
    yr = [(+v.mid()).mid() for v in y]                     # exact 128-bit dyadic midpoints (radius 0): the seed is a point, not a ball
    rad192 = max(float(v.rad()) for v in y)
    tq, G = closure(R, AB, op, yr)
    seed.update(y_enc=[enc(v) for v in yr], y_float=np.array([float(v.mid()) for v in yr]), bump_eps=float(eps.mid()), rad192=rad192,
                turn_defect=float((tq + arb.pi()).mid()), far_point=(float(G.real.mid()), float(G.imag.mid())))
    for kk in ('lam', 'L', 'X0', 't'): seed[kk] = enc(+seed[kk])
    seed['turn_defect_arb'] = enc(tq + arb.pi()); seed['far_imag_arb'] = enc(G.imag)
    print(f'  final (128-bit dyadic seed; max 192-bit radius {rad192:.1e}): turning+pi = {(tq + arb.pi()).str(3, radius=True)}, far point = {G.real.mid().str(20)} + {G.imag.str(3, radius=True)} i;'
          f' bump eps = {seed["bump_eps"]:.3e}  ({time.time()-t0:.1f}s)', flush=True)
    # projection error against the exact curve at random points: arb (128-bit state) and float state
    ctx.prec = PREC
    rng = np.random.default_rng(1); ellf = float(yr[R.n].mid())
    xs = sorted(A(float(v)) for v in np.r_[rng.uniform(0, 40, 150), rng.uniform(40, ellf, 150)])
    ke = C1.kap_batch(xs, [0], lam)[0]
    xk = [float(v.mid()) for v in x]
    ea, ef = [], []
    for xv, kv in zip(xs, ke):
        e = min(int(np.searchsorted(xk, float(xv.mid())) - 1), R.N - 1); e = max(e, 0)
        kp = op.poly(op.local(yr, h, H, e)); ea.append(abs(float((kp((xv - x[e]) / h[e]) - kv).mid())))
        ef.append(abs(R.eval_at(seed['y_float'][:R.n], ellf, np.array([float(xv.mid())]))[0] - float(kv.mid())))
    ea, ef = np.array(ea), np.array(ef); xf = np.array([float(v.mid()) for v in xs])
    seed['proj_err'] = dict(tip=float(ea[xf <= 40].max()), far=float(ea[xf > 40].max()), float_state=float(ef.max()))
    print(f'  projection |kappa_h - kappa_P1| (300 random x): x<=40 {ea[xf <= 40].max():.2e}, x>40 {ea[xf > 40].max():.2e}; float state {ef.max():.2e}')
    pickle.dump(seed, open(f'seed_p1_{seed["slab"]}.pkl', 'wb'), protocol=4)
    print(f'  wrote seed_p1_{seed["slab"]}.pkl  ({time.time()-t0:.1f}s)')
