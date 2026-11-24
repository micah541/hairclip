"""T7 (C4), modulation Gram matrices per slab (normalised directions g^_i = g_i / ||chi g_i||_rho, g = (p, n, V_s)):
  Gamma0  = int_{10.25}^{ell} Wbar (zeta g^_i')(zeta g^_j'),     Gamma2 = int Wbar (zeta g^_i')''(zeta g^_j')''   ((I0), (I2)),
  Gamma2E = int_0^{34.25} rho (chi g^_i)''(chi g^_j)''                                                          ((E2), AB form),
with Wbar = W + (1 + ell^2/4)^-2 >= the flattened weight (X >= ell/2 on [ell/2, ell]).  Outputs sqrt(trace) >= sqrt(lambda_max);
merge.py forms C_gamma0 = C_x sqrt(lmax Gamma0), C_gamma2 = C_x sqrt(lmax Gamma2), G_x = C_x sqrt(lmax Gamma2E)
(|x| <= C_x ||u||_rho in normalised coordinates, C_x from c2r3.py).
Directions on the whole curve from T6's bivariate base (p, q, n, t1 sharp; Frenet derivatives):
  p' = -k q, p'' = -k' q - k (1 + k p), p''' = -k'' q - 2 k' (1 + k p) - k (k' p + k p');  n' = -k t1, n'' = -k' t1 - k^2 n,
  n''' = -k'' t1 - 3 k k' n - k^2 n';  V_s^(m) = -k^(m+2) + c p^(m) + b n^(m).
usage: python3 farg7.py --slabs a-b --out DIR [--nsub 1]      then  --refine 0.1 --nsub 4  (sharpens the peak slabs in place)"""
import os, sys, json, time, argparse, math
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path: sys.path.insert(0, HERE)
from flint import arb, arb_poly, ctx
import base7
from tm import TM, from_series, up
ctx.prec = 128
S_OF_U = arb_poly([-1, 2])
XR0, XR1, XRHO = arb('10.25'), arb('16.25'), arb('34.25')


def cutp(a, h, x0, L, up_, n=4):
    S = arb_poly([0, 0, 0, 0, 35, -84, 70, -20])
    U = arb_poly([(a - x0) / L, h / L]) if up_ else arb_poly([(x0 - a) / L, -h / L])
    out = []; Sd = S; f = arb(1) / L if up_ else -arb(1) / L
    for m in range(n):
        out.append(TM(Sd(U) * f ** m)); Sd = Sd.derivative()
    return out


def run_slab(s, nsub=1, verbose=True):
    import c2r3 as R3
    from rep import STIP
    t0 = time.time(); S = base7.Slab(s); SB = R3.SlabBP(S, arb(10 ** 7)); B = SB.B
    R3.BP = R3.BP if hasattr(R3, 'BP') else None
    G0 = [[arb(0)] * 3 for _ in range(3)]; G2 = [[arb(0)] * 3 for _ in range(3)]; G2E = [[arb(0)] * 3 for _ in range(3)]
    nrm = [arb(0)] * 3
    res = None
    for j in range(nsub):
        t = arb((2 * j + 1) / nsub - 1, 1 / nsub); sc = {k: v(t) for k, v in SB.sc.items()}; c, b = sc['c'], sc['b']
        ell = B.ell(t); Cl = (1 + ell * ell / 4) ** (-2)
        g0 = [[arb(0)] * 3 for _ in range(3)]; g2 = [[arb(0)] * 3 for _ in range(3)]; g2e = [[arb(0)] * 3 for _ in range(3)]; nr = [arb(0)] * 3
        for e, el in enumerate(SB.els):
            xe, he = B.x[e](t), B.h[e](t)
            if float(xe.lower()) > float(ell.upper()): break
            ih = 2 / he                                                         # ds/dx
            xa0, xb0 = float(xe.lower()), float((xe + he).upper())
            cuts = [0.0, 1.0]
            if e < B.mesh.JT:                                                   # fixed knots: split at the cutoff breakpoints
                for bp in (10.25, 16.25, 29 - STIP, 32 - STIP):
                    u = (bp - float(xe.mid())) / float(he.mid())
                    if 1e-12 < u < 1 - 1e-12: cuts.append(u)
            cuts = sorted(cuts)
            kb = el['k0']; ks = [kb]
            for m in range(5): ks.append(ks[-1].ds())
            ksb = [base7.bp_at(ks[m], t) for m in range(6)]; pqnt = [base7.bp_at(el[k_], t) for k_ in ('p', 'q', 'n', 't1')]
            # V_s and its x-derivatives polynomial-first in tc (they are small differences of O(1) terms at late times)
            D_ = B.D; T6 = R3.T; ihT = T6.tconst(ih)
            ihP = [T6.tconst(arb(1))]                                          # (ds/dx)^m = (2/h_e)^m as TPs in tc (moving outer knots)
            for m in range(1, 6): ihP.append(T6.tmul(ihP[-1], B.ih[e] * 2, D_))
            kx_ = [ks[0]] + [ks[m].tmul(ihP[m], D_) for m in range(1, 6)]
            pB, qB, nB, tB = el['p'], el['q'], el['n'], el['t1']
            oneB = R3.BP([arb_poly([1])], True)
            p1B = -(kx_[0].mul(qB, D_)); q1B = oneB + kx_[0].mul(pB, D_)
            p2B = -(kx_[1].mul(qB, D_)) - kx_[0].mul(q1B, D_)
            q2B = kx_[1].mul(pB, D_) + kx_[0].mul(p1B, D_)
            p3B = -(kx_[2].mul(qB, D_)) - kx_[1].mul(q1B, D_).scale(arb(2)) - kx_[0].mul(q2B, D_)
            n1B = -(kx_[0].mul(tB, D_)); t11B = kx_[0].mul(nB, D_)
            n2B = -(kx_[1].mul(tB, D_)) - kx_[0].mul(t11B, D_)
            t12B = kx_[1].mul(nB, D_) + kx_[0].mul(n1B, D_)
            n3B = -(kx_[2].mul(tB, D_)) - kx_[1].mul(t11B, D_).scale(arb(2)) - kx_[0].mul(t12B, D_)
            cT, bT = SB.sc['c'], SB.sc['b']
            VB = [(-kx_[m + 2]) + PB.tmul(cT, D_) + NB.tmul(bT, D_) for m, (PB, NB) in enumerate(((pB, nB), (p1B, n1B), (p2B, n2B), (p3B, n3B)))]
            Vsb = [base7.bp_at(v, t) for v in VB]
            for u0, u1 in zip(cuts[:-1], cuts[1:]):
                if e < B.mesh.JT and u0 > 0: u0a = (arb(10.25) if abs(float(xe.mid()) + u0 * float(he.mid()) - 10.25) < 1e-9 else None)
                U0 = arb(u0); U1 = arb(u1)
                S_ = arb_poly([2 * U0 - 1, 2 * (U1 - U0)]); a = xe + he * U0; h = he * (U1 - U0)
                K = [TM(ksb[m](S_) * ih ** m) for m in range(6)]
                p, q, n, t1 = [TM(v(S_)) for v in pqnt]
                k0, k1, k2, k3, k4, k5 = K
                one = TM(arb_poly([1]))
                p1 = -1 * (k0 * q); p2 = -1 * (k1 * q) - k0 * (one + k0 * p); p3 = -1 * (k2 * q) - 2 * (k1 * (one + k0 * p)) - k0 * (k1 * p + k0 * p1)
                n1 = -1 * (k0 * t1); n2 = -1 * (k1 * t1) - (k0 * k0) * n; n3 = -1 * (k2 * t1) - 3 * (k0 * k1 * n) - (k0 * k0) * n1
                P = [p, p1, p2, p3]; N = [n, n1, n2, n3]
                V = [TM(Vsb[m](S_)) for m in range(4)]                        # polynomial-first in tc
                G = [P, N, V]
                xa, xb = float(a.lower()), float((a + h).upper()); xm = float((a + h / 2).mid())
                if xm < float(XRHO.mid()):
                    rho = from_series(lambda X: (-(1 + X * X) ** (arb(2) / 3) / 4).exp(), a, h, 24)[0]
                    if xm + STIP < 29: chi = [one, TM(arb_poly([0])), TM(arb_poly([0]))]
                    else: chi = cutp(a, h, arb(32) - arb(STIP), arb(3), False, 3)
                    for i in range(3):
                        f = chi[0] * G[i][0]; nr[i] += (rho * f * f).integral01() * h
                    d2 = [chi[2] * G[i][0] + 2 * (chi[1] * G[i][1]) + chi[0] * G[i][2] for i in range(3)]
                    for i in range(3):
                        for k in range(3): g2e[i][k] += (rho * d2[i] * d2[k]).integral01() * h
                if xm > 10.25:
                    z = cutp(a, h, XR0, arb(6), True, 3) if xm < 16.25 else [one, TM(arb_poly([0])), TM(arb_poly([0]))]
                    W = from_series(lambda X: (1 + X * X) ** (-2), a, h, 24)[0] + Cl
                    f1 = [z[0] * G[i][1] for i in range(3)]
                    f2 = [z[2] * G[i][1] + 2 * (z[1] * G[i][2]) + z[0] * G[i][3] for i in range(3)]
                    for i in range(3):
                        for k in range(3):
                            g0[i][k] += (W * f1[i] * f1[k]).integral01() * h; g2[i][k] += (W * f2[i] * f2[k]).integral01() * h
        # normalise and take traces (sup over the balls)
        nl = [arb(float(v.lower())) for v in nr]
        tr = lambda M: sum((M[i][i] / nl[i] for i in range(3)), arb(0))
        r = dict(tr0=float(tr(g0).upper()), tr2=float(tr(g2).upper()), tr2E=float(tr(g2e).upper()), nrm=[float(v.sqrt().lower()) for v in nl],
                 d0=[float((g0[i][i] / nl[i]).upper()) for i in range(3)], d2=[float((g2[i][i] / nl[i]).upper()) for i in range(3)],
                 d2E=[float((g2e[i][i] / nl[i]).upper()) for i in range(3)])
        mx = lambda a_, b_: [max(x, y) for x, y in zip(a_, b_)] if isinstance(a_, list) else max(a_, b_)
        res = r if res is None else {k: (mx(res[k], r[k]) if k != 'nrm' else [min(x, y) for x, y in zip(res[k], r[k])]) for k in r}
    out = dict(slab=s, tau=[S.tau0, S.tau0 + S.dtf], k=S.k0, sqrt_lmax_G0=math.sqrt(res['tr0']), sqrt_lmax_G2=math.sqrt(res['tr2']),
               sqrt_lmax_G2E=math.sqrt(res['tr2E']), norm_chi_g_min=res['nrm'], nsub=nsub, t=time.time() - t0,
               sqrt_diag_G0=[math.sqrt(v) for v in res['d0']], sqrt_diag_G2=[math.sqrt(v) for v in res['d2']], sqrt_diag_G2E=[math.sqrt(v) for v in res['d2E']])
    if verbose:
        print(f"slab {s} (k {S.k0:.1f}): sqrt tr Gamma0 {out['sqrt_lmax_G0']:.4e}  Gamma2 {out['sqrt_lmax_G2']:.4e}  Gamma2E {out['sqrt_lmax_G2E']:.4e} "
              f"(||chi g|| >= {[round(v, 4) for v in res['nrm']]}) ({out['t']:.0f}s)", flush=True)
    return out



def parse_slabs(spec):
    """'a-b' or comma-separated slabs / ranges, e.g. '5-20' or '5,13,21-23'"""
    out = []
    for part in spec.split(','):
        if '-' in part: a_, b_ = [int(v) for v in part.split('-')]; out += list(range(a_, b_ + 1))
        elif part.strip(): out.append(int(part))
    return out

def combine(a, b):
    """two valid runs of the same slab -> the sharper bound per entry (upper bounds: min; the norm lower bounds: max)"""
    out = dict(b)
    for k in ('sqrt_lmax_G0', 'sqrt_lmax_G2', 'sqrt_lmax_G2E'): out[k] = min(a[k], b[k])
    for k in ('sqrt_diag_G0', 'sqrt_diag_G2', 'sqrt_diag_G2E'): out[k] = [min(x, y) for x, y in zip(a[k], b[k])]
    out['norm_chi_g_min'] = [max(x, y) for x, y in zip(a['norm_chi_g_min'], b['norm_chi_g_min'])]
    out['runs'] = [{k: r[k] for k in ('nsub', 'sqrt_lmax_G0', 'sqrt_lmax_G2', 'sqrt_lmax_G2E', 't')} for r in (a, b)]
    return out


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--slabs', default='5-5'); ap.add_argument('--out', default=os.path.join(HERE, 'out_farg'))
    ap.add_argument('--nsub', type=int, default=1)
    ap.add_argument('--refine', type=float, default=None, help='rerun with --nsub the slabs whose existing output (smaller nsub) has '
                    'sqrt_lmax_G0 > REFINE (the tc-hull of one slab-wide ball overestimates there, before outer-element insertions)')
    base7.base_arg(ap)
    A = ap.parse_args(); os.makedirs(A.out, exist_ok=True); base7.set_base(A.base)
    for s in parse_slabs(A.slabs):
        fn = os.path.join(A.out, 'farg_%03d.json' % s)
        if os.path.exists(fn):
            if A.refine is None: continue
            old = json.load(open(fn))
            if old.get('nsub', 1) >= A.nsub or old['sqrt_lmax_G0'] <= A.refine: continue
            new = combine(old, run_slab(s, A.nsub)); tmp = fn + '.tmp'
            json.dump(new, open(tmp, 'w'), indent=1); os.replace(tmp, fn)
            continue
        json.dump(run_slab(s, 1 if A.refine is not None else A.nsub), open(fn, 'w'), indent=1)
