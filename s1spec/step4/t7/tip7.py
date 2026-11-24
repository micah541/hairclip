"""T7 (C4), tip side, per slab (sup over the slab: nsub tc-balls x cells in x; sharp polynomial-first base).
E2 level (rho weight, x in [0, 34.25]):  (L u)'' = M u'' + e1 u' + e0 u,  M = -d^4 + b2 d^2 + (2 b2' + b1) d + (b2'' + 2 b1' + b0),
   G22 >= sup (-mu~2)^+ (Garding function of M in rho: Re<M v, v>_rho <= -1/2 ||v''||^2 - int rho mu~2 v^2),
   G12 >= sup |e1|, e1 = b1'' + 2 b0';  G02 >= sup |e0|, e0 = b0''.
Strip (chi_32, sigma in [29, 32], E0 level, rho weight): -<[L, chi] f, chi f>_rho = int A f^2 + B f'^2 with g_j = rho chi chi^(j),
   A = -2 g1''' + 3 g2'' - 2 g3' + g4 - b1 g1 + (b2 g1)' - b2 g2,  B = 6 g1' - 6 g2   (step 3 coupling_design / agent_coup, plus b2);
   sB = sup B+/W,  IA = int A+,  sA = 2 int A+ Phi_H (Phi_H(x) = int_16.25^x W^-1),  sA0 = 2 IA / rho(16.25);
   E2 strip: d_j = sup |coef_j| sqrt(rho/W), [L, chi] e = sum_j coef_j e^(j):  coef_3 = -4 chi', coef_2 = -6 chi'',
   coef_1 = -4 chi''' + 2 b2 chi', coef_0 = -chi'''' + b2 chi'' + b1 chi'   (reported as d3, d2, d1, d0).
Ramp (zeta on [10.25, 16.25], W = (1+x^2)^-2, operator L1 of v = e': c2 = -k^2, c1 = -8 k k' + Q - c x):
   Az = -2 g1''' + 3 g2'' - 2 g3' + g4 - c1 g1 + (c2 g1)' - c2 g2, Bz = 6 g1' - 6 g2 (g_j = W zeta zeta^(j));
   a_f = sup Az+/rho, b_f = sup Bz+/rho;  c_j = sup |coef_j| sqrt(W/rho) for [L1, zeta] v:
   c3 = 4|zeta'|, c2 = 6|zeta''|, c1 = |4 zeta''' - 2 c2 zeta'|, c0 = |zeta'''' - c2 zeta'' - c1 zeta'|.
Weight constants: Psi1 = sup_{x >= 10.25} |W'/W| = psi', c1W = sup (psi'^2 - psi'')/2  (exact functions of x, the same for all slabs).
Modulation in (E2), AB form: G_x <= C_x sqrt(lambda_max Gamma2) with Gamma2 the rho-Gram of (chi g_i)'' / ||chi g_i||_rho
   (normalised directions; C_x from c2r3.py, merged in merge.py); here gx = sqrt(trace Gamma2) over a tc-ball covering the slab.
usage: python3 tip7.py --slabs a-b --out DIR [--nsub 2]"""
import os, sys, json, time, argparse, math
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path: sys.path.insert(0, HERE)
from flint import arb, arb_poly, arb_series, ctx
import base7
from rep import STIP
ctx.prec = 128
XRHO = arb('34.25'); XR0, XR1 = arb('10.25'), arb('16.25')
SIG0, SIG1 = arb(29), arb(32)


def jets_fn(F, x, n=6):
    """derivatives 0..n-1 of F at the ball x via arb_series"""
    S = F(arb_series([x, 1], prec=n)); cf = S.coeffs(); cf = cf + [arb(0)] * (n - len(cf))
    return [cf[m] * math.factorial(m) for m in range(n)]


phiF = lambda X: (1 + X * X) ** (arb(2) / 3) / 4
WlogF = lambda X: (1 + X * X).log() * 2                     # psi = -ln W
rhoF = lambda X: (-(1 + X * X) ** (arb(2) / 3) / 4).exp()
WF = lambda X: (1 + X * X) ** (-2)


def S_ser(u):
    return u ** 4 * (35 - 84 * u + 70 * u * u - 20 * u ** 3)


def chi_jets(x, n=6):
    """chi_32 (1 for sigma <= 29, 0 for sigma >= 32) and derivatives at the ball x (cells inside the strip or outside)"""
    sig = x + arb(STIP)
    if float(sig.upper()) <= 29: return [arb(1)] + [arb(0)] * (n - 1)
    if float(sig.lower()) >= 32: return [arb(0)] * n
    return jets_fn(lambda X: S_ser((SIG1 - (X + arb(STIP))) / 3), x, n)


def zeta_jets(x, n=6):
    if float(x.upper()) <= 10.25: return [arb(0)] * n
    if float(x.lower()) >= 16.25: return [arb(1)] + [arb(0)] * (n - 1)
    return jets_fn(lambda X: S_ser((X - XR0) / 6), x, n)


def gmu(ps, c2, c2p, c2pp, c1, c1p, c0):
    p1, p2, p3, p4 = ps[1], ps[2], ps[3], ps[4]
    w2 = p1 * p1 - p2; w2p = 2 * p1 * p2 - p3; w2pp = 2 * p2 * p2 + 2 * p1 * p3 - p4
    rr2 = lambda f, fp, fpp: fpp - 2 * p1 * fp + w2 * f
    A = 2 * w2 - c2; Ap = 2 * w2p - c2p; App = 2 * w2pp - c2pp
    hlf = arb(1) / 2
    V = c0 - hlf * (c1p - p1 * c1) + hlf * rr2(c2, c2p, c2pp) - hlf * rr2(w2, w2p, w2pp)
    return -hlf * (A * A) - hlf * rr2(A, Ap, App) - V


def leib(f, g, m):
    return sum((math.comb(m, k) * f[k] * g[m - k] for k in range(m + 1)), arb(0))


def ibp_AB(wt, cut, b1, b2, b2p):
    """A, B of -<[L, cut] f, cut f>_wt = int A f^2 + B f'^2 (wt, cut: jet lists; b1, b2, b2' balls)"""
    g = [[leib(wt, cut, 0) * 0 for _ in range(4)] for _ in range(5)]
    wc = [leib(wt, cut, m) for m in range(5)]                          # (wt cut)^(m)
    gj = {}
    for j in range(1, 5):                                              # g_j = wt cut cut^(j): jets 0..4-j+... up to order 3
        cj = [cut[j + m] if j + m < len(cut) else arb(0) for m in range(4)]
        gj[j] = [leib(wc, cj, m) for m in range(4)]
    A = -2 * gj[1][3] + 3 * gj[2][2] - 2 * gj[3][1] + gj[4][0] - b1 * gj[1][0] + (b2p * gj[1][0] + b2 * gj[1][1]) - b2 * gj[2][0]
    B = 6 * gj[1][1] - 6 * gj[2][0]
    return A, B


class PhiH:
    """Phi_H(x) = int_16.25^x W^-1 = int (1+t^2)^2 dt (exact polynomial)"""
    @staticmethod
    def up(x):
        F = lambda t: t + 2 * t ** 3 / 3 + t ** 5 / 5
        return F(x) - F(XR1)


def base_cell(el, U, B):
    h = el['h']; ih = 1 / h; x = el['x'] + h * U
    kp = el['kp']; K = []
    for m in range(5):
        K.append(kp(U) * ih ** m); kp = kp.derivative()
    k0, k1, k2, k3, k4 = K; Q = el['Q'](U); c = B.c
    b2 = -k0 * k0; b2p = -2 * k0 * k1; b2pp = -2 * k1 * k1 - 2 * k0 * k2
    b1 = -6 * k0 * k1 + Q - c * x; b1p = -5 * k1 * k1 - 6 * k0 * k2 - c; b1pp = -16 * k1 * k2 - 6 * k0 * k3
    b0 = c - 3 * k1 * k1 - 4 * k0 * k2; b0p = -10 * k1 * k2 - 4 * k0 * k3; b0pp = -10 * k2 * k2 - 14 * k1 * k3 - 4 * k0 * k4
    return dict(x=x, K=K, Q=Q, b2=b2, b2p=b2p, b2pp=b2pp, b1=b1, b1p=b1p, b1pp=b1pp, b0=b0, b0p=b0p, b0pp=b0pp,
                c2v=b2, c1v=-8 * k0 * k1 + Q - c * x)


def tm_range(T_, nsc=16):
    """[lower, upper] of a TM on [0,1]: Horner on nsc sub-cells of the polynomial part, plus/minus the remainder"""
    lo, hi = math.inf, -math.inf; r = float(T_.r.upper())
    for i in range(nsc):
        v = T_.p(arb((i + 0.5) / nsc, 0.5 / nsc)); lo = min(lo, float(v.lower())); hi = max(hi, float(v.upper()))
    return lo - r, hi + r


def cut_jets(a, h, x0, L, up_=True, n=5):
    """exact polys in s of S((x - x0)/L) (up) or S((x0 - x)/L) (down) and x-derivatives; caller ensures the element is inside"""
    S = arb_poly([0, 0, 0, 0, 35, -84, 70, -20])
    U = arb_poly([(a - x0) / L, h / L]) if up_ else arb_poly([(x0 - a) / L, -h / L])
    out = []; Sd = S; f = arb(1) / L if up_ else -arb(1) / L
    for m in range(n):
        out.append(TM(Sd(U) * f ** m)); Sd = Sd.derivative()
    return out


def lb(f, g, m):
    return sum((math.comb(m, k) * f[k] * g[m - k] for k in range(m + 1)), TM(arb_poly([0])))


def ibp_tm(wt, cut, b1, b2, b2p):
    wc = [lb(wt, cut, m) for m in range(5)]
    gj = {}
    for j in range(1, 5):
        cj = [cut[j + m] if j + m < len(cut) else TM(arb_poly([0])) for m in range(4)]
        gj[j] = [lb(wc, cj, m) for m in range(4)]
    A = -2 * gj[1][3] + 3 * gj[2][2] - 2 * gj[3][1] + gj[4][0] - b1 * gj[1][0] + (b2p * gj[1][0] + b2 * gj[1][1]) - b2 * gj[2][0]
    B = 6 * gj[1][1] - 6 * gj[2][0]
    return A, B


def el_mesh(a, b, h):
    n = max(1, int(math.ceil((b - a) / h))); return [(arb(a) + (arb(b) - arb(a)) * k / n, (arb(b) - arb(a)) / n) for k in range(n)]


def run_slab(s, nsub=2, verbose=True):
    import r3core as CD
    from tm import TM as _TM, from_series, sup01
    global TM; TM = _TM
    t0 = time.time(); S = base7.Slab(s)
    acc = dict(G22=0.0, G12=0.0, G02=0.0, sB=0.0, IA=0.0, sA=0.0, d3=0.0, d2=0.0, d1=0.0, d0=0.0,
               a_f=0.0, b_f=0.0, c3=0.0, c2=0.0, c1=0.0, c0=0.0)
    # meshes (split at T1 knots so that the base is one polynomial per element)
    for j in range(nsub):
        t = arb((2 * j + 1) / nsub - 1, 1 / nsub); B = S.base(t, XRHO + 1)
        knots = [float(v.mid()) for v in B.xk]
        def elems(a, b, h):
            pts = sorted(set([a, b] + [k_ for k_ in knots if a < k_ < b]))
            out = []
            for u, v in zip(pts[:-1], pts[1:]): out += el_mesh(u, v, h)
            return out
        def bpolys(a, h):
            E = CD.t1_index(B, a); L = B.local(E, a, h); K = L['K']; x = arb_poly([a, h]); c = B.c
            k0, k1, k2, k3, k4 = K
            return dict(b2=-k0 * k0, b2p=-2 * k0 * k1, b2pp=-2 * k1 * k1 - 2 * k0 * k2, b1=-6 * k0 * k1 + L['Q'] - x * c,
                        b1p=-5 * k1 * k1 - 6 * k0 * k2 - arb_poly([c]), b1pp=-16 * k1 * k2 - 6 * k0 * k3,
                        b0=arb_poly([c]) - 3 * k1 * k1 - 4 * k0 * k2, b0p=-10 * k1 * k2 - 4 * k0 * k3,
                        b0pp=-10 * k2 * k2 - 14 * k1 * k3 - 4 * k0 * k4, c1v=-8 * k0 * k1 + L['Q'] - x * c)
        # E2 Garding on [0, 34.25]
        for a, h in elems(0.0, 34.25, 0.25):
            P = bpolys(a, h); d, em, e1 = CD.phi_tms(a, h); ps = [None] + list(d)
            T_ = {k: TM(v) for k, v in P.items()}
            mu2 = gmu(ps, T_['b2'], T_['b2p'], T_['b2pp'], 2 * T_['b2p'] + T_['b1'], 2 * T_['b2pp'] + T_['b1p'], T_['b2pp'] + 2 * T_['b1p'] + T_['b0'])
            acc['G22'] = max(acc['G22'], -tm_range(mu2)[0])
            acc['G12'] = max(acc['G12'], float(sup01(P['b1pp'] + 2 * P['b0p']).upper()))
            acc['G02'] = max(acc['G02'], float(sup01(P['b0pp']).upper()))
        # strip sigma in [29, 32]  (x = sigma - STIP)
        IA = 0.0; sA = 0.0
        for a, h in elems(29 - STIP, 32 - STIP, 0.25):
            P = bpolys(a, h); T_ = {k: TM(v) for k, v in P.items()}
            rho = from_series(rhoF, a, h, 24, nder=4); chi = cut_jets(a, h, arb(32) - arb(STIP), arb(3), up_=False)
            A, Bc = ibp_tm(rho, chi, T_['b1'], T_['b2'], T_['b2p'])
            Wi = from_series(lambda X: (1 + X * X) ** 2, a, h, 24)[0]
            acc['sB'] = max(acc['sB'], max(tm_range(Bc * Wi)[1], 0.0))
            Au = max(tm_range(A)[1], 0.0); hf = float(h.upper())
            IA += Au * hf; sA += Au * hf * float(PhiH.up(arb(float((a + h).upper()))).upper())
            srw = from_series(lambda X: ((-(1 + X * X) ** (arb(2) / 3) / 4).exp() * (1 + X * X) ** 2).sqrt(), a, h, 24)[0]
            coef = [-1 * chi[4] + T_['b2'] * chi[2] + T_['b1'] * chi[1], -4 * chi[3] + 2 * (T_['b2'] * chi[1]), -6 * chi[2], -4 * chi[1]]
            for jj in range(4): acc[f'd{jj}'] = max(acc[f'd{jj}'], float((coef[jj] * srw).sup().upper()))
        acc['IA'] = max(acc['IA'], IA); acc['sA'] = max(acc['sA'], 2 * sA)
        # zeta ramp x in [10.25, 16.25]
        for a, h in elems(10.25, 16.25, 0.25):
            P = bpolys(a, h); T_ = {k: TM(v) for k, v in P.items()}
            W = from_series(WF, a, h, 24, nder=4); z = cut_jets(a, h, XR0, arb(6), up_=True)
            Az, Bz = ibp_tm(W, z, T_['c1v'], T_['b2'], T_['b2p'])
            irho = from_series(lambda X: ((1 + X * X) ** (arb(2) / 3) / 4).exp(), a, h, 24)[0]
            acc['a_f'] = max(acc['a_f'], max(tm_range(Az * irho)[1], 0.0)); acc['b_f'] = max(acc['b_f'], max(tm_range(Bz * irho)[1], 0.0))
            swr = from_series(lambda X: ((1 + X * X) ** (arb(2) / 3) / 4).exp().sqrt() / (1 + X * X), a, h, 24)[0]
            c2v = T_['b2']
            coef = [z[4] - c2v * z[2] - T_['c1v'] * z[1], 4 * z[3] - 2 * (c2v * z[1]), 6 * z[2], 4 * z[1]]
            for jj in range(4): acc[f'c{jj}'] = max(acc[f'c{jj}'], float((coef[jj] * swr).sup().upper()))
    acc['sA0'] = 2 * acc['IA'] / float(rhoF(XR1).lower())
    out = dict(slab=s, tau=[S.tau0, S.tau0 + S.dtf], k=S.k0, nsub=nsub, t=time.time() - t0, **acc)
    if verbose:
        print(f"slab {s} (k {S.k0:.1f}): G22 {acc['G22']:.3f} G12 {acc['G12']:.3f} G02 {acc['G02']:.3f} | strip sB {acc['sB']:.3e} sA {acc['sA']:.3e} "
              f"sA0 {acc['sA0']:.2e} d3..0 {acc['d3']:.3e} {acc['d2']:.3e} {acc['d1']:.3e} {acc['d0']:.3e} | ramp a_f {acc['a_f']:.4f} b_f {acc['b_f']:.4f} "
              f"c3..0 {acc['c3']:.3f} {acc['c2']:.3f} {acc['c1']:.3f} {acc['c0']:.3f} ({out['t']:.0f}s)", flush=True)
    return out


def weight_consts():
    """Psi1 = sup_{x >= 10.25} psi', c1W = sup (psi'^2 - psi'')/2 for W = (1+x^2)^-2 (monotone in x beyond 10.25: values at 10.25,
    checked on a grid of balls up to 1e6)"""
    p1m = 0.0; c1m = 0.0; x = 10.25
    while x < 1e6:
        xb = arb(x * 1.0025, x * 0.0025); ps = jets_fn(WlogF, xb)
        p1m = max(p1m, float(abs(ps[1]).upper())); c1m = max(c1m, float(((ps[1] ** 2 - ps[2]) / 2).upper())); x *= 1.005
    return dict(Psi1=p1m, c1W=c1m)



def parse_slabs(spec):
    """'a-b' or comma-separated slabs / ranges, e.g. '5-20' or '5,13,21-23'"""
    out = []
    for part in spec.split(','):
        if '-' in part: a_, b_ = [int(v) for v in part.split('-')]; out += list(range(a_, b_ + 1))
        elif part.strip(): out.append(int(part))
    return out

if __name__ == '__main__':
    import base7
    ap = argparse.ArgumentParser(); ap.add_argument('--slabs', default='5-5'); ap.add_argument('--out', default=os.path.join(HERE, 'out_tip'))
    ap.add_argument('--nsub', type=int, default=2)
    base7.base_arg(ap)
    A = ap.parse_args(); os.makedirs(A.out, exist_ok=True); base7.set_base(A.base)
    wf = os.path.join(A.out, 'weight_consts.json')
    if not os.path.exists(wf): json.dump(weight_consts(), open(wf, 'w'), indent=1)
    
    for s in parse_slabs(A.slabs):
        fn = os.path.join(A.out, 'tip_%03d.json' % s)
        if os.path.exists(fn): continue
        json.dump(run_slab(s, A.nsub), open(fn, 'w'), indent=1)
