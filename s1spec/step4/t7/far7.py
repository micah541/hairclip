"""T7 (C4), far curve, per slab (whole far curve x in [10.25, ell~]): Garding functions of the E_I levels.
  v = zeta e', (I0) operator L1 = -d^4 + c2 d^2 + c1 d + c0, c2 = -k^2, c1 = -8 k k' + Q - c x, c0 = -8 k'^2 - 10 k k'' (c = 1/4 - a~);
  (I2) operator L1^(2): c1 -> c1 + 2 c2', c0 -> c0 + 2 c1' + c2''.
  Weight: W = (1 + X^2)^-2, X = the far-end flattening of T8 (weights.Yl_jet, s0 = 1/2) in x: X = ell Phi(x/ell), X = x for x <= ell/2.
  mu = -A^2/2 - (W A)''/(2W) - V + psi_tau/2  (eicon.gaarding_mu, psi = -ln W; the last term is the tau-derivative of the
  time-dependent flattened weight, (1/2) d/dtau int W v^2 = int W v v_tau - (1/2) int W psi_tau v^2; psi_tau >= 0 since
  ell_tau > 0 and Phi - s Phi' >= 0; it restores the rate -c that the flat weight loses to the drift divergence).
  Hardy term (T5 (I0), x >= 16.25): hardy = sup |x b0'| 2/(2 P_eff + 1), P_eff = max(x psi'/2, 1/2), b0' = -10 k' k'' - 4 k k'''.
  mu1 = inf_{x >= 10.25} mu - hardy;  split at mu_f:  -int W mu v^2 <= -mu_f V0^2 + a_mu W1^2 with
  a_mu = sup over the cells where mu - hardy < mu_f of W (mu_f - mu + hardy) zeta^2 / rho (these must lie where chi_32 = 1).
Rigour: per T1 element, per tc-ball (sharp polynomial-first base), cells in u; all quantities by ball arithmetic on the cell;
the derivative jets of psi by arb_series at the cell ball (hull of both branches on cells straddling x = ell/2);
cells failing the target are bisected (depth <= 8).
With --se the cells on x >= 0.6 ell_lo use the base jets (k^(0..3), Q) widened to the (SE)-corrected base (se7.py).
usage: python3 far7.py --slabs a-b --out DIR [--muf 0.15] [--nsub 4] [--se]"""
import os, sys, json, time, argparse, math
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path: sys.path.insert(0, HERE)
from flint import arb, arb_poly, arb_series, ctx
import base7
import tp6 as T6
ctx.prec = 128
X0R, X1R = arb('10.25'), arb('16.25')
XCHI = arb(29) - arb(base7.STIP) if hasattr(base7, 'STIP') else None
from rep import STIP
XCHI = arb(29) - arb(STIP)          # chi_32 = 1 for x <= 29 - STIP
PSI = (arb(21) / 4, -arb(49) / 8, arb(15) / 8)
S0 = arb(1) / 2


def hull(a, b):
    lo = min(float(a.lower()), float(b.lower())); hi = max(float(a.upper()), float(b.upper()))
    return arb((lo + hi) / 2, (hi - lo) / 2 * (1 + 1e-15) + 1e-300)


def psi_jets(x, s, ell, ellt, n=5):
    """psi^(m)(x), m = 0..4, at x = ell s, psi = 2 ln(1 + X^2), X = ell Phi(x/ell); plus psi_tau = 4 X X_tau / (1 + X^2),
    X_tau = ell_tau (Phi - s Phi')  (>= 0).  s (sharp ball), ell, ell_tau balls; hull of both branches if s straddles 1/2."""
    def ser(branch):
        if branch == 'lin': X = arb_series([x, 1], prec=n)
        else: sS = arb_series([s, 1 / ell], prec=n)
        if branch == 'flat':
            u = (sS - S0) * 2
            a, b, c = PSI
            Pint = u ** 6 * (a / 6) + u ** 7 * (b / 7) + u ** 8 * (c / 8)
            X = (S0 + (u - Pint) * S0) * ell
        F = (1 + X * X).log() * 2
        cf = F.coeffs(); cf = cf + [arb(0)] * (n - len(cf))
        xc = X.coeffs(); X0, X1 = xc[0], xc[1]
        Xt = ellt * (X0 / ell - s * X1) if branch == 'flat' else arb(0)
        return [cf[m] * math.factorial(m) for m in range(n)] + [4 * X0 * Xt / (1 + X0 * X0)]
    s_lo = float(s.lower()); s_hi = float(s.upper())
    if s_hi <= 0.5: return ser('lin')
    if s_lo >= 0.5: return ser('flat')
    A, B = ser('lin'), ser('flat'); return [hull(p, q) for p, q in zip(A, B)]


def gmu(ps, c2, c2p, c2pp, c1, c1p, c0):
    p1, p2, p3, p4 = ps[1], ps[2], ps[3], ps[4]
    w2 = p1 * p1 - p2; w2p = 2 * p1 * p2 - p3; w2pp = 2 * p2 * p2 + 2 * p1 * p3 - p4
    rr2 = lambda f, fp, fpp: fpp - 2 * p1 * fp + w2 * f
    A = 2 * w2 - c2; Ap = 2 * w2p - c2p; App = 2 * w2pp - c2pp
    V = c0 - (c1p - p1 * c1) / 2 + rr2(c2, c2p, c2pp) / 2 - rr2(w2, w2p, w2pp) / 2
    return -A * A / 2 - rr2(A, Ap, App) / 2 - V


def smooth(u):
    """S(u) = u^4 (35 - 84u + 70u^2 - 20u^3) clipped to [0,1] (ball)"""
    lo, hi = float(u.lower()), float(u.upper())
    if hi <= 0: return arb(0)
    if lo >= 1: return arb(1)
    uc = arb((max(lo, 0) + min(hi, 1)) / 2, (min(hi, 1) - max(lo, 0)) / 2)
    return uc ** 4 * (35 - 84 * uc + 70 * uc * uc - 20 * uc ** 3)


def cell(el, U, B):
    """all needed quantities on the cell u in U (ball) of element el of base B (el['se'], el['sg']: sharp x_e/ell, h_e/ell)"""
    h = el['h']; ih = 1 / h
    if el['fixed']: x = el['x'] + h * U; s = x / B.ell            # tip elements: exact knots
    else: s = el['se'] + el['sg'] * U; x = B.ell * s               # outer knots scale with ell: x/ell sharp
    kp = el['kp']; K = []
    for m in range(4):
        K.append(kp(U) * ih ** m); kp = kp.derivative()
    Q = el['Q'](U); c = B.c
    if SEW is not None and float(x.upper()) >= SEW['xc']:             # (SE)-corrected base on x >= 0.6 ell_lo (se7.py)
        K = [K[m] + arb(0, SEW['k'][m]) for m in range(4)]; Q = Q + arb(0, SEW['Q'])
    k, k1, k2, k3 = K
    c2 = -k * k; c2p = -2 * k * k1; c2pp = -2 * k1 * k1 - 2 * k * k2
    c1 = -8 * k * k1 + Q - c * x; c1p = -7 * k1 * k1 - 8 * k * k2 - c; c0 = -8 * k1 * k1 - 10 * k * k2
    ps = psi_jets(x, s, B.ell, B.ellt); pst = ps[5]
    mu = gmu(ps, c2, c2p, c2pp, c1, c1p, c0) + pst / 2                    # + (1/2) psi_tau: d/dtau of the flattened weight
    mu2 = gmu(ps, c2, c2p, c2pp, c1 + 2 * c2p, c1p + 2 * c2pp, c0 + 2 * c1p + c2pp) + pst / 2
    b0p = -10 * k1 * k2 - 4 * k * k3
    Peff = x * ps[1] / 2; Pl = max(float(Peff.lower()), 0.5)
    hardy = abs(x * b0p) * 2 / (2 * arb(Pl) + 1)
    return dict(x=x, mu=mu, mu2=mu2, hardy=hardy)


SEW = None


def run_slab(s, muf=0.15, nsub=4, depth=8, verbose=True, se=False):
    global SEW
    t0 = time.time(); S = base7.Slab(s)
    if se:
        import se7; SEW = se7.widths(s)
    else: SEW = None
    res = dict(slab=s, tau=[S.tau0, S.tau0 + S.dtf], k=S.k0, muf=muf, nsub=nsub, se=bool(se))
    mu_min = math.inf; mu_min_at = None; mu2_min = math.inf; hardy = 0.0; ncell = 0
    bad = []                                # (x_lo, x_hi, mu_lo) of cells below the floor before the Hardy term
    good_min = math.inf
    for j in range(nsub):
        t = arb((2 * j + 1) / nsub - 1, 1 / nsub); B = S.base(t, None)
        SB = S.SB; il = T6.tpow(SB.ell, -1, SB.D)
        for el in B.els:
            el['se'] = T6.tmul(SB.x[el['e']], il, SB.D)(t); el['sg'] = T6.tmul(SB.h[el['e']], il, SB.D)(t)
            el['fixed'] = el['e'] < SB.mesh.JT
            xe_hi = float((el['x'] + el['h']).upper())
            if xe_hi < float(X0R.lower()): continue
            xr = max(float(el['x'].lower()), 1.0); hr = float(el['h'].upper())
            m = max(4, int(math.ceil(64 * hr / xr)))
            stack = [(k_ / m, (k_ + 1) / m, 0) for k_ in range(m)]
            while stack:
                a, b, d = stack.pop()
                U = arb((a + b) / 2, (b - a) / 2)
                C = cell(el, U, B); ncell += 1
                xl, xh = float(C['x'].lower()), float(C['x'].upper())
                if xh < float(X0R.lower()): continue
                if xl < float(X0R.lower()):                       # cell straddles 10.25: restrict u (x >= 10.25 only matters)
                    ua = float(((X0R - el['x']) / el['h']).lower())
                    if ua > a and d < depth: stack.append((max(a, ua - 1e-12), b, d + 1)); continue
                ml = float(C['mu'].lower()); m2l = float(C['mu2'].lower())
                if xl >= float(X1R.lower()) - 1e-9: hardy = max(hardy, float(C['hardy'].upper()))
                if ml < muf + 0.02 and d < depth and (b - a) > 1e-6:          # refine near the floor
                    mid = (a + b) / 2; stack += [(a, mid, d + 1), (mid, b, d + 1)]; continue
                if ml < mu_min: mu_min, mu_min_at = ml, (xl + xh) / 2
                mu2_min = min(mu2_min, m2l)
                if ml < muf + 0.02: bad.append((xl, xh, ml))
                else: good_min = min(good_min, ml)
    # split: cells with mu - hardy < muf need x <= 29 - STIP; a_mu = sup W (muf - mu + hardy) zeta^2 / rho
    a_mu = 0.0; bad_hi = 0.0; split_ok = True
    for xl, xh, ml in bad:
        if ml - hardy >= muf: continue
        bad_hi = max(bad_hi, xh)
        if xh > float(XCHI.lower()): split_ok = False; continue
        x = arb((xl + xh) / 2, (xh - xl) / 2)
        W = (1 + x * x) ** (-2); rho = (-(1 + x * x) ** (arb(2) / 3) / 4).exp(); z = smooth((x - X0R) / 6)
        a_mu = max(a_mu, float((W * (muf - ml + hardy) * z * z / rho).upper()))
    res.update(mu_raw_min=mu_min, mu_raw_min_at=mu_min_at, hardy=hardy, mu1_raw=mu_min - hardy, mu2=mu2_min,
               mu_min_outside_split=good_min - hardy, split_ok=split_ok, a_mu=a_mu, split_zone_hi=bad_hi,
               mu1=(muf if split_ok and good_min - hardy >= muf else min(muf, good_min - hardy)), ncell=ncell, t=time.time() - t0)
    if verbose:
        print(f"slab {s} (k {S.k0:.1f}): inf mu = {mu_min:.4f} at x = {mu_min_at:.2f}; hardy {hardy:.2e}; mu1 = {res['mu1']:.4f} "
              f"(split at {muf}: zone <= {bad_hi:.2f}, a_mu {a_mu:.2e}, ok {split_ok}); mu outside {good_min - hardy:.4f}; mu2 = {mu2_min:.4f} "
              f"({ncell} cells, {res['t']:.0f}s)", flush=True)
    return res



def parse_slabs(spec):
    """'a-b' or comma-separated slabs / ranges, e.g. '5-20' or '5,13,21-23'"""
    out = []
    for part in spec.split(','):
        if '-' in part: a_, b_ = [int(v) for v in part.split('-')]; out += list(range(a_, b_ + 1))
        elif part.strip(): out.append(int(part))
    return out

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--slabs', default='5-5'); ap.add_argument('--out', default=os.path.join(HERE, 'out_far'))
    ap.add_argument('--muf', type=float, default=0.15); ap.add_argument('--nsub', type=int, default=4)
    ap.add_argument('--se', action='store_true', help='far cells on the (SE)-corrected base (se7.py)')
    base7.base_arg(ap)
    A = ap.parse_args(); os.makedirs(A.out, exist_ok=True); base7.set_base(A.base)
    for s in parse_slabs(A.slabs):
        fn = os.path.join(A.out, 'far_%03d.json' % s)
        if os.path.exists(fn): continue
        json.dump(run_slab(s, A.muf, A.nsub, se=A.se), open(fn, 'w'), indent=1)
