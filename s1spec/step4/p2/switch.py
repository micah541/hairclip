"""P2 switch (task A), certified in arb: the true curve at P1's checkpoint against the v2 base at an interior time tc_a.

True curve = P1 centre (exact dyadic sine coefficients u_n, length L, tip-normalised scale lam_P = K0/kappa_tip)
           + P1 error set (|dtheta^(j)| <= |e|_j, j <= 4, recomputed here outward from G, w, box, tail T n^-6)
           + length uncertainty (scale factor 1 + d_L, |d_L| <= rL/(L - rL)).
Base = T6's exact definition (t6/base.py SlabBase on t1/polished_v2) at tc_a.
Bounds, per cell of the base mesh (x = arclength from the tip, cells of width <= HMAX):
  dk_k  >= sup |kappa_true^(k) - kappa_base^(k)| at equal arclength (k <= 3): P1-centre Taylor model (order K, rigorous
           remainder sum_n |u_n| n^(k+1) (n da)^K / K!) minus the base polynomial, + ball + length part;
  phi   = theta_true - theta_base (equal arclength): |e|_0 + |signed cumulative int (kP1 - kbase)| + length part;
  |D(x)| <= int_0^x |phi|  (D = Gamma_true - Gamma_base at equal arclength);
  normal graph e over the base: |e| <= sup|D|, foot shift |sigma| <= 1.01 sup|D|, |e'| <= (1+k|e|) tan(phi + k sigma),
  e'' .. e'''' from the exact curvature identity of a normal graph,
      e'' = [k1 S (S-1)/(sqrt S + 1) + Dk S^(3/2) - e'(k1' e + k1 e')]/(1 - k1 e),  S = (1-k1 e)^2 + e'^2,
      Dk = kappa_true(s(xi)) - k1(xi) = k1' sigma + k1'' sigma^2/2 + dk(xi) + dk'(xi) sigma + R,  sigma' = (S-1)/(sqrt S + 1),
  solved as truncated power series (arb_series) with interval coefficients valid on the whole cell.
Energies (T5 1.1, T8 F_j with crude weight majorants) -> Y_raw, F_j raw; modulation by Newton-Kantorovich (see kantorovich()).
usage: python3 switch.py CKPT_JSON BASE_SLAB TC_A [HMAX] [K]  -> out/switch_<slab>.json"""
import sys, os, json, math, time
S4 = '/home/claude/s1spec/step4'
sys.path.insert(0, f'{S4}/t6'); sys.path.insert(0, f'{S4}/t1')
import numpy as np
from flint import arb, arb_poly, arb_mat, arb_series, ctx
from base import SlabBase
import tp as T
from rep import K0

POLDIR = f'{S4}/t1/polished_v2'
STIP = arb('-2.249077651121112'); X0 = arb('2.2490776515')
EPS, EPSI = arb('1e-5'), arb('3e-5')
NS = 6                                      # series length: e, e', ..., e^(5) coefficients (we use j <= 4)


def ub(x): return math.nextafter(float(arb(x).abs_upper()), math.inf)   # outward
def ball(r): return arb(0, r) if not isinstance(r, arb) else arb(0, r.abs_upper())


def S_(u):  # smoothstep and derivatives (u may be an arb ball; clipping by cases)
    return u ** 4 * (35 - 84 * u + 70 * u ** 2 - 20 * u ** 3)
def S1_(u): return 140 * u ** 3 * (1 - u) ** 3
def S2_(u): return 420 * u ** 2 * (1 - u) ** 2 * (1 - 2 * u)


def step_hull(a, b, x0, w):
    """hulls of S((x - x0)/w) and its first two x-derivatives for x in [a, b] (S clipped to [0, 1])"""
    ua, ub_ = (a - x0) / w, (b - x0) / w
    if float(ub_.upper()) <= 0: return arb(0), arb(0), arb(0)
    if float(ua.lower()) >= 1: return arb(1), arb(0), arb(0)
    lo, hi = max(0.0, float(ua.lower())), min(1.0, float(ub_.upper()))
    u = arb(lo).union(arb(hi))
    return S_(u), S1_(u) / w, S2_(u) / w ** 2


# ------------------------------------------------------------------------------------------------ P1 side
class P1:
    B = 64
    def __init__(self, path, MMAX):
        rec = json.load(open(path)); self.rec = rec
        j2a = lambda v: arb(int(v[0])) * arb(2) ** int(v[1])
        def j2b(v):
            mid = j2a(v); rad = arb(int(v[2])) * arb(2) ** int(v[3]); return arb(mid, rad) if int(v[2]) else mid
        self.u = [j2a(v) for v in rec['xbar']]; self.M = len(self.u)
        self.L, self.Lmid = j2b(rec['L']), j2a(rec['L']); self.t = j2b(rec['t'])
        nz = max(i for i, (m, e) in enumerate(rec['xbar']) if int(m) != 0) + 1
        B = self.B; self.Q = Q = (nz + B) // B
        cf = [arb(0)] + self.u[:nz] + [arb(0)] * (Q * B - nz - 1)
        self.cf = cf
        self.Cm = [arb_mat(Q, B, [arb(n) ** m * cf[n] for n in range(Q * B)]) for m in range(MMAX + 1)]
        z = np.load(path[:-5] + '.npz'); self.G, self.w, self.box = z['G'], np.asarray(z['w'], float), np.asarray(z['box'], float)
        self.T = float(rec['T'])
        self.absu = [abs(v) for v in cf]

    def theta_batch(self, als, ms):
        """{m: [theta^(m)(alpha)]}, theta = -pi/2 + alpha + sum u_n sin(n alpha) (m >= 1 only)"""
        B, Q = self.B, self.Q; out = {m: [] for m in ms}; P = len(als)
        sr = [[None] * P for _ in range(B)]; cr = [[None] * P for _ in range(B)]
        sq = [[None] * P for _ in range(Q)]; cq = [[None] * P for _ in range(Q)]
        for p, a in enumerate(als):
            for r in range(B): sr[r][p], cr[r][p] = (r * a).sin_cos()
            for q in range(Q): sq[q][p], cq[q][p] = (q * B * a).sin_cos()
        Cr = arb_mat(B, P, [v for row in cr for v in row]); Sr = arb_mat(B, P, [v for row in sr for v in row])
        for m in ms:
            Ic = self.Cm[m] * Cr; Is = self.Cm[m] * Sr
            for p in range(P):
                if m % 2: v = sum((cq[q][p] * Ic[q, p] - sq[q][p] * Is[q, p] for q in range(Q)), arb(0))
                else: v = sum((sq[q][p] * Ic[q, p] + cq[q][p] * Is[q, p] for q in range(Q)), arb(0))
                v = v if m % 4 < 2 else -v
                if m == 1: v = v + 1
                out[m].append(v)
        return out

    def ball_norms(self):
        """|e|_j = sum_n n^j sup|du_n| (j = 0..4): float sums of nonnegative terms inflated by (1 + 1e-10), tail T n^-6"""
        col = np.abs(self.G) @ self.w; M = len(col); Nc = len(self.box)
        n = np.arange(1, Nc + 1, dtype=float); tot = self.box.copy(); tot[:M] += col
        out = []
        for j in range(5):
            s = float(np.sum(n ** j * tot)) * (1 + 1e-10)
            tail = self.T * Nc ** (j - 5) / (5 - j) * (1 + 1e-10)
            out.append(arb(repr(s)) + arb(repr(tail)))
        return [arb(float(v.upper())) for v in out]

    def taylor_rem(self, k, da, K):
        """sum_n |u_n| n^(k+1) (n da)^K / K!  (remainder of the K-term Taylor model of theta^(k+1))"""
        s = arb(0); fK = arb(math.factorial(K))
        for nn, a in enumerate(self.absu):
            if nn and a != 0: s += a * arb(nn) ** (k + 1) * (arb(nn) * da) ** K
        return s / fK


# ------------------------------------------------------------------------------------------------ main computation
def run(ck, slab, tca, HMAX=0.25, K=30, verbose=True, mod=None):
    t0 = time.time(); ctx.prec = 128
    SB = SlabBase(slab, poldir=POLDIR); ctx.prec = 128
    tc = arb(tca); n = SB.n
    ell = SB.ell(tc); tb = SB.Y[n + 3](tc); lnlam = SB.Y[n + 1](tc)
    P = P1(ck, 5 + K + 1)
    eb = P.ball_norms()
    pi = arb.pi(); cL = 2 * pi / P.Lmid
    th1pi = P.theta_batch([pi], [1])[1][0]
    ktip = -cL * th1pi; lamP = arb(K0) / ktip; s = lamP * cL           # alpha = pi - s x
    rL = P.L.rad(); dL = rL / (P.Lmid - rL)                             # |c_t/c_m - 1| <= dL
    D0 = arb(0); dL0 = dL; m1 = m2 = Sc = arb(0)
    if mod:                                                             # modulation balls: scale e^{mu1}, translation mu2, clock s
        dL = (1 + dL) * (1 + arb('1.01') * mod['m1']) - 1               # combined scale factor deviation (curvature level)
        m1, m2, Sc = mod['m1'], mod['m2'], mod['S']
    import reader as RD
    T6D = os.environ.get('T6DIR', f'{S4}/p2/t6')               # env T6DIR: one T6 data set for jets and residuals
    J = RD.jets(slab, d=f'{T6D}/jets'); cmax = arb(max(abs(v) for v in [r for r in json.load(open(f'{T6D}/residuals.json'))['per_slab'] if r['slab'] == slab][0]['scalars']['range_c']))
    def eh(name, e): a = J[name][e]; return arb(float(np.max(np.abs(a))))
    jel = []
    for e in range(SB.N):
        g = lambda nm: np.max(np.abs(J[nm][e]), axis=-1)            # per sub-cell |.| hull
        th_t = g('k3') + g('k0') ** 2 * g('k1') + g('k0') * g('Q') + float(cmax.upper()) * g('x') * g('k0') + g('r1') + g('k0') * g('A')
        jel.append(dict(P=eh('p', e), N=eh('n', e), V=eh('V', e), W=eh('varpi', e), Th=arb(float(np.max(th_t)) * (1 + 1e-12))))
    thmax = sum((arb(nn) * a for nn, a in enumerate(P.absu)), arb(0)) + 1   # sup |theta'|
    if verbose:
        print(f'P1 {os.path.basename(ck)}: k_tip {ktip.str(15)}, t {P.t.str(15)}; base slab {slab} tc {tca}: t~ {tb.str(15)}, '
              f'k~ {float((arb(K0) / lnlam.exp()).mid()):.10f}, ell {float(ell.mid()):.6f}; |e|_j = ' +
              ', '.join(f'{ub(v):.3e}' for v in eb), flush=True)
    # ---- cells on the base mesh at tc
    Ns, nrm = SB.Ns, SB.nrm
    cells = []
    for e in range(SB.N):
        xa, xb = SB.x[e](tc), SB.x[e + 1](tc); h = xb - xa
        d, _ = SB.local_coeffs(e)
        pe = arb_poly([0])
        for q, dq in enumerate(d): pe = pe + Ns[q] * dq(tc)                 # kappa_base(s), s in [-1, 1]
        nsub = max(1, math.ceil(ub(h) / HMAX))
        for j in range(nsub):
            ca = xa + h * j / nsub; cb = xa + h * (j + 1) / nsub; xc = (ca + cb) / 2; hc = (cb - ca) / 2
            sc = (xc - xa) * 2 / h - 1
            pxi = pe(arb_poly([sc, 2 / h]))                                  # kappa_base(xc + xi)
            cells.append(dict(e=e, a=ca, b=cb, xc=xc, hc=hc, p=pxi))
    if verbose: print(f'  {len(cells)} cells ({time.time() - t0:.1f}s)', flush=True)
    # ---- P1 centre Taylor coefficients at the cell centres
    MT = 5 + K
    ths = P.theta_batch([pi - s * c['xc'] for c in cells], list(range(1, MT + 1)))
    if verbose: print(f'  P1 jets done ({time.time() - t0:.1f}s)', flush=True)
    rem_cache = {}
    def rem(k, hc):
        hk = math.ceil(ub(hc) * 1e6) / 1e6; key = (k, hk)
        if key not in rem_cache: rem_cache[key] = s ** (k + 1) * P.taylor_rem(k, s * arb(repr(hk)), K)
        return rem_cache[key]
    fac = [arb(math.factorial(j)) for j in range(MT + 2)]
    phi_det = arb(0); Dcum = D0; out_cells = []
    for i, c in enumerate(cells):
        hc = c['hc']; pb = c['p']; xi = arb(0, hc.upper())
        dk = []; kb = []
        for k in range(5):
            a = arb_poly([(-1) ** (k + 1) * s ** (k + 1) * ths[k + 1 + j][i] * (-s) ** j / fac[j] for j in range(K)])
            diff = a - pb
            dk.append(sum((abs(cc) * hc.upper() ** j for j, cc in enumerate(diff.coeffs())), arb(0)) + rem(k, hc))
            if k == 0:
                I = diff.integral(); cint = I(hc) - I(-hc) + 2 * hc * ball(rem(0, hc))
                aint = sum((abs(cc) * hc.upper() ** (j + 1) * 2 / (j + 1) for j, cc in enumerate(diff.coeffs())), arb(0)) + 2 * hc * rem(0, hc)
            pb = pb.derivative()
        pp = c['p']
        for j in range(6):
            kb.append(pp(xi).abs_upper()); pp = pp.derivative()
        if mod:                                                         # clock: base at tau_a + s, |s| <= S
            Kt = mod['Ktau'][c['e']]
            dk = [v + mod['S'] * Kt[k] for k, v in enumerate(dk)]
            kb = [v + (mod['S'] * Kt[j]).abs_upper() for j, v in enumerate(kb)]
        c.update(dk=dk, kb=kb, phi_left=phi_det, aint=aint)
        phi_det = phi_det + cint
    # ---- per-cell bounds: angle, displacement (rest + explicit modulation pieces), normal component, series
    res = []
    DS1a = arb('2.500806276619802')
    for i, c in enumerate(cells):
        nb = [cells[j] for j in (i - 1, i, i + 1) if 0 <= j < len(cells)]
        xr = arb(float(c['b'].upper()))
        Kb = [max((cc['kb'][j] for cc in nb), key=ub) for j in range(6)]
        Dk = []
        for k in range(4):
            det = max((cc['dk'][k] for cc in nb), key=ub)
            detn = max((cc['dk'][k + 1] for cc in nb), key=ub)
            Lk = dL * ((k + 1) * (Kb[k] + det) + xr * (Kb[k + 1] + detn)) * (1 + dL) ** (k + 2)
            bk = s ** (k + 1) * eb[k + 1] * (1 + dL) ** (k + 1)
            Dk.append(det + Lk + bk)
        je = jel[c['e']]
        Kloc = Kb[0] + max((cc['dk'][0] for cc in nb), key=ub)
        Prest = eb[0] + abs(c['phi_left']).upper() + c['aint']
        Pang = Prest + arb('1.02') * (dL0 + m1) * xr * Kloc + Sc * je['Th']
        Dcum = Dcum + 2 * c['hc'] * Prest
        Dtot = Dcum + arb('1.01') * dL0 * xr + arb('1.02') * m1 * (2 * xr + DS1a) + arb('1.01') * m2 + Sc * (je['V'] + je['W'])
        em = (Dcum + arb('1.01') * dL0 * (je['P'] + DS1a) + arb('1.02') * m1 * (je['P'] + Dtot) + arb('1.01') * m2 * je['N'] + Sc * je['V']
              + Sc * Sc * je['W'] * je['Th'] + 2 * (dL0 ** 2 + m1 ** 2) * xr ** 2 * Kloc + arb('1.02') * (m1 + dL0) * xr * Pang
              + Dtot * Sc * je['Th'])
        c.update(Pang=Pang, Dtot=Dtot, em=em, Dk=Dk, Kb=Kb)
    for i, c in enumerate(cells):
        nb = [cells[j] for j in (i - 1, i, i + 1) if 0 <= j < len(cells)]
        Dt = max((cc['Dtot'] for cc in nb), key=ub); em = max((cc['em'] for cc in nb), key=ub)
        Pang = max((cc['Pang'] for cc in nb), key=ub)
        Kb = c['Kb']; Dk = [max((cc['Dk'][k] for cc in nb), key=ub) for k in range(4)]
        kbar = Kb[0]
        assert ub(3 * kbar * Dt) < 0.0099, 'foot point lemma needs 3 k |D| < 0.0099'
        sig = arb('1.01') * Dt
        E0 = em + 2 * kbar * Dt ** 2 + Dt * (Pang + kbar * Dt)     # normal graph vs normal component of D (second order)
        ang = Pang + kbar * sig
        E1 = (1 + kbar * E0) * ang.tan() if ub(ang) < 0.5 else None
        E = series_bounds(Kb, Dk, E0, E1, sig)
        c['E'] = E
        res.append(dict(x=[float(c['a'].lower()), float(c['b'].upper())], E=[ub(v) for v in E], Dk=[ub(v) for v in Dk]))
    if verbose: print(f'  bounds done ({time.time() - t0:.1f}s); max |e| {max(r["E"][0] for r in res):.3e}, '
                      f'max |e\'| {max(r["E"][1] for r in res):.3e}, max |dk| {max(r["Dk"][0] for r in res):.3e}', flush=True)
    en = energies(cells, ell)
    meta = dict(ckpt=ck, slab=slab, tc=tca, tau=SB.tau0 + SB.d['dt'] * (1 + tca) / 2, k_tip=ub(ktip), lamP=ub(lamP), t_P1=P.t.str(20),
                t_base=tb.str(20), ell=float(ell.mid()), ball=[ub(v) for v in eb], dL=ub(dL), K=K, HMAX=HMAX, ncells=len(cells),
                phi_det_end=ub(phi_det), cpu=time.time() - t0)
    return meta, en, res, cells, SB


def bp_sup(bp):
    """upper bound of |bp(s, tc)| for s, tc in [-1, 1]"""
    return sum((sum((abs(c) for c in p.coeffs()), arb(0)) for p in bp.P), arb(0))


def ktau_bounds(SB):
    """per element: sup over element x slab of |d^k/dx^k kappa_tau|_x| (k = 0..5), from T6 engine.element_a (BP exact in s)"""
    import engine as EN
    out = []
    z = T.tconst(arb(0))
    for e in range(SB.N):
        A_ = EN.element_a(SB, e, z, z, z)
        kt, f2 = A_['kt'], A_['f2']; f2s = T.tsup(f2)
        row = []
        for k in range(6):
            row.append(bp_sup(kt) * f2s ** k); kt = kt.ds()
        out.append(row)
    return out


def series_bounds(Kb, Dk, E0, E1, sig):
    """enclosures of |e^(j)| (j = 0..4) on a cell from the exact normal-graph curvature identity (truncated series)"""
    L = NS
    def ser(c): return arb_series(c + [arb(0)] * (L - len(c)), prec=L)
    k1 = ser([ball(Kb[j]) / math.factorial(j) for j in range(L)])
    dk = ser([ball(Dk[j]) / math.factorial(j) for j in range(min(L, 4))])
    k1p = ser([ball(Kb[j + 1]) / math.factorial(j) for j in range(L - 1)])
    k1pp = ser([ball(Kb[j + 2]) / math.factorial(j) for j in range(L - 2)])
    dkp = ser([ball(Dk[j + 1]) / math.factorial(j) for j in range(3)])
    # remainder of the Taylor expansions in sigma (|sigma| <= sig, sigma derivatives first order small): crude ball per coefficient
    SP = arb('1e-6')                                       # |sigma'|, |sigma''| <= 1e-6 (asserted below via S - 1)
    Rr = 4 * ((Kb[3] + Kb[4] + Kb[5]) * sig ** 2 * (sig + SP) + (Dk[2] + Dk[3]) * sig * (1 + SP))
    ec = [ball(E0), ball(E1)] + [arb(0)] * (L - 2)
    sc = [ball(sig)] + [arb(0)] * (L - 1)
    for order in range(L - 2):
        e = ser(ec)
        ep = ser([ec[j + 1] * (j + 1) for j in range(L - 1)])
        sgm = ser(sc)
        a = 1 - k1 * e
        S = a * a + ep * ep
        Sm1 = k1 * k1 * e * e - 2 * k1 * e + ep * ep
        rS = S.sqrt()
        Dser = k1p * sgm + k1pp * sgm * sgm / 2 + dk + dkp * sgm + ser([ball(Rr)] * L)
        rhs = (k1 * S * Sm1 / (rS + 1) + Dser * S * rS - ep * (k1p * e + k1 * ep)) / a
        sp = Sm1 / (rS + 1)                                   # sigma'
        rc = rhs.coeffs() + [arb(0)] * L; spc = sp.coeffs() + [arb(0)] * L
        # e'' series = rhs: coefficient of t^order gives ec[order + 2]
        ec[order + 2] = rc[order] / ((order + 1) * (order + 2))
        if order + 1 < L: sc[order + 1] = spc[order] / (order + 1)
    e = ser(ec); ep = ser([ec[j + 1] * (j + 1) for j in range(L - 1)])
    a = 1 - k1 * e; S = a * a + ep * ep; Sm1 = k1 * k1 * e * e - 2 * k1 * e + ep * ep
    spc = (Sm1 / (S.sqrt() + 1)).coeffs() + [arb(0)] * L
    assert ub(spc[0]) < 1e-6 and ub(2 * spc[1]) < 1e-6, 'sigma derivative bound used in Rr'
    return [ec[j].abs_upper() * math.factorial(j) for j in range(5)]


def energies(cells, ell):
    """A, B (rho, chi), V0, V2 (W_inf, zeta_inf on [10.25, ell]), F_j (j <= 3, weights majorised: T_j <= 28 y^(2j-2), T_0 <= 112/y^2)"""
    A2 = B2 = V02 = V22 = arb(0); F = [arb(0)] * 4
    X34 = arb('34.25')
    for c in cells:
        a, b = c['a'], c['b']; h = b - a; E = c['E']
        aL = arb(math.nextafter(float(a.lower()), -math.inf)); bU = arb(math.nextafter(float(b.upper()), math.inf))
        if float(aL.mid()) < 34.25:
            rho = (-(1 + aL * aL) ** (arb(2) / 3) / 4).exp()          # rho decreasing: sup at the left end
            chi, chi1, chi2 = step_hull(32 - (bU + STIP), 32 - (aL + STIP), arb(0), arb(3))
            chi, chi1, chi2 = chi.abs_upper(), chi1.abs_upper(), chi2.abs_upper()
            hh = min(b, X34, key=ub) - a
            A2 += rho * hh * (chi * E[0]) ** 2
            B2 += rho * hh * (chi2 * E[0] + 2 * chi1 * E[1] + chi * E[2]) ** 2
        if float(bU.mid()) > 10.25:
            W = (1 + aL * aL) ** -2
            z, z1, z2 = step_hull(aL, bU, arb('10.25'), arb(6))
            z, z1, z2 = z.abs_upper(), z1.abs_upper(), z2.abs_upper()
            V02 += W * h * (z * E[1]) ** 2
            V22 += W * h * (z2 * E[1] + 2 * z1 * E[2] + z * E[3]) ** 2
        y = aL - X0
        if float(bU.mid() - X0.mid()) > 14:
            yl = max(arb(14), y, key=lambda v: float(v.mid()))
            yu = bU - X0
            for j in range(4):
                wgt = 112 / (yl * yl) if j == 0 else 28 * yu ** (2 * j - 2)
                F[j] += wgt * h * E[j + 1] ** 2
    return dict(A=A2.sqrt(), B=B2.sqrt(), V0=V02.sqrt(), V2=V22.sqrt(), F=F)


def report(meta, en, mue):
    Y = (en['A'] ** 2 + EPS * en['B'] ** 2 + mue * (en['V0'] ** 2 + EPSI * en['V2'] ** 2)).sqrt()
    return dict(meta=meta, A=ub(en['A']), B=ub(en['B']), V0=ub(en['V0']), V2=ub(en['V2']), Y_raw=ub(Y), F_raw=[ub(v) for v in en['F']])


def _main():
    ck, slab, tca = sys.argv[1], int(sys.argv[2]), float(sys.argv[3])
    HMAX = float(sys.argv[4]) if len(sys.argv) > 4 else 0.25; K = int(sys.argv[5]) if len(sys.argv) > 5 else 30
    T7 = json.load(open(f'{S4}/t7/constants.json')); mue = arb(repr(T7['meta']['closure']['mu_e']))
    meta, en, res, cells, SB = run(ck, slab, tca, HMAX, K)
    out = report(meta, en, mue); out['cells'] = res
    json.dump(out, open(f'{S4}/p2/out/switch_raw_s{slab}.json', 'w'), indent=1)
    print({k: v for k, v in out.items() if k not in ('cells', 'meta')})


# ------------------------------------------------------------------------------------------------ modulation (IFT) at the switch
def chi_rho_norm(cells):
    s = arb(0)
    for c in cells:
        aL = arb(math.nextafter(float(c['a'].lower()), -math.inf))
        if float(aL.mid()) >= 34.25: break
        s += (-(1 + aL * aL) ** (arb(2) / 3) / 4).exp() * (min(c['b'], arb('34.25'), key=ub) - c['a'])
    return s.sqrt()


def ift(SB, slab, cells, A_raw, Ktau, apriori, verbose=True):
    """contraction G(mu^) = mu^ - A^{-1} F(mu^), A = -M(tau_a) (T7 normalisation, |M^{-1}| <= mbar):
       |G(0)| <= beta = mbar sqrt3 A_raw;  |DG| <= q = mbar (|M(tau_a + s) D(tau_a)^-1 - M(tau_a)| + |E^|) on the a-priori ball;
       unique zero with |mu^| <= R^ = beta/(1-q) if the ball B(0, R^) lies in the a-priori ball."""
    sys.path.insert(0, f'{S4}/t6'); import reader as RD
    T7 = json.load(open(f'{S4}/t7/constants.json')); s7 = T7['slabs'][str(slab)]
    fg = json.load(open(f'{S4}/t7/run_v2/out_farg/farg_{slab:03d}.json'))
    T6D = os.environ.get('T6DIR'); T6 = json.load(open(f'{T6D}/residuals.json' if T6D else f'{S4}/t6/residuals.json')); r6 = [r for r in T6['per_slab'] if r['slab'] == slab][0]
    mbar = arb(repr(s7['C3']['mbar']))
    gmin = [arb(repr(v)) for v in fg['norm_chi_g_min'][:2]] + [arb(repr(min(s7['C3']['nV_min'], fg['norm_chi_g_min'][2])))]
    X34 = 34.25
    def H(name): lo, hi = RD.hull(slab, name, 0.0, X34, d=(f'{T6D}/jets' if T6D else None)); return arb(max(abs(lo), abs(hi)))
    Pb, Qb, Nb, T1b, Vb, Wb, rb, k0b = (H(n) for n in ('p', 'q', 'n', 't1', 'V', 'varpi', 'r', 'k0'))
    # Phi = theta_tau (fixed x) on [0, 34.25]: |Phi| <= int_0^x sup|kappa_tau|
    Phib = arb(0)
    for e in range(SB.N):
        xa = float(SB.x[e](arb(0)).mid())
        if xa >= X34: break
        Phib += Ktau[e][0] * T.tsup(SB.h[e])
    Kt2 = max((Ktau[e][2] for e in range(SB.N) if float(SB.x[e](arb(0)).mid()) < X34), key=ub)
    two_dt = 2 / SB.dt; K0a = arb(K0); DS1a = arb('2.500806276619802')
    k2d, k4d = SB.k2t.derivative() * two_dt, SB.k4t.derivative() * two_dt
    ctau = -(k4d + k2d * K0a ** 2) * (1 / K0a); btau = k2d + ctau * DS1a
    ctb, btb = T.tsup(ctau), T.tsup(btau)
    cb = arb(max(abs(v) for v in r6['scalars']['range_c'])); bb = arb(max(abs(v) for v in r6['scalars']['range_b']))
    nchi = chi_rho_norm(cells)
    dg = [Vb + Phib * Qb, Phib * T1b, Kt2 + ctb * Pb + cb * (Vb + Phib * Qb) + btb * Nb + bb * Phib * T1b]
    Lcol = [nchi * dg[i] / gmin[i] for i in range(3)]
    dM = apriori['S'] * (3 * sum((L * L for L in Lcol), arb(0))).sqrt()
    # epsilon terms on the chi-region from the modulated cell bounds (valid on the whole a-priori ball)
    tip = [c for c in cells if float(c['a'].lower()) < X34]
    E0 = max((c['E'][0] for c in tip), key=ub); E1 = max((c['E'][1] for c in tip), key=ub)
    tpsi = E1 / (1 - k0b * E0)
    eps = [E0 + Qb * tpsi, T1b * tpsi + arb('1.01') * apriori['m1'] * Nb, rb + (Wb + E0 * Phib) * tpsi]
    Eh = (3 * sum(((nchi * eps[i] / gmin[i]) ** 2 for i in range(3)), arb(0))).sqrt()
    q = mbar * (dM + Eh)
    beta = mbar * arb(3).sqrt() * A_raw
    assert ub(q) < 0.5, f'contraction factor {ub(q)}'
    R = beta / (1 - q)
    m1, m2, S = R / gmin[0], R / gmin[1], R / gmin[2]
    ok = ub(m1) <= ub(apriori['m1']) and ub(m2) <= ub(apriori['m2']) and ub(S) <= ub(apriori['S'])
    out = dict(mbar=ub(mbar), gmin=[float(v.lower()) for v in gmin], chi_rho=ub(nchi), Phi_bar=ub(Phib), dg=[ub(v) for v in dg],
               dM=ub(dM), E_hat=ub(Eh), q=ub(q), beta=ub(beta), R_hat=ub(R), mu1=ub(m1), mu2=ub(m2), s=ub(S), inside_apriori=ok,
               apriori={k: ub(v) for k, v in apriori.items() if k != 'Ktau'})
    if verbose: print('IFT:', json.dumps(out), flush=True)
    return out


def certify(ck, slab, tca, HMAX=0.25, K=30, fac=1.02):
    """raw run -> a-priori modulation ball (fac x the raw IFT radius) -> modulated run -> IFT check -> Y(tau0), F_j(tau0)"""
    T7 = json.load(open(f'{S4}/t7/constants.json')); mue = arb(repr(T7['meta']['closure']['mu_e']))
    meta, en, res, cells, SB = run(ck, slab, tca, HMAX, K)
    raw = report(meta, en, mue)
    Ktau = ktau_bounds(SB)
    s7 = T7['slabs'][str(slab)]; fg = json.load(open(f'{S4}/t7/run_v2/out_farg/farg_{slab:03d}.json'))
    gm = [fg['norm_chi_g_min'][0], fg['norm_chi_g_min'][1], min(s7['C3']['nV_min'], fg['norm_chi_g_min'][2])]
    R0 = arb(repr(s7['C3']['mbar'])) * arb(3).sqrt() * en['A'] * fac
    apr = dict(m1=R0 / arb(repr(gm[0])), m2=R0 / arb(repr(gm[1])), S=R0 / arb(repr(gm[2])))
    meta2, en2, res2, cells2, _ = run(ck, slab, tca, HMAX, K, mod=dict(apr, Ktau=Ktau))
    mod = report(meta2, en2, mue)
    I = ift(SB, slab, cells2, en['A'], Ktau, apr)
    out = dict(raw={k: v for k, v in raw.items()}, modulated={k: v for k, v in mod.items() if k != 'meta'}, ift=I,
               Ktau_max=[max(ub(r[k]) for r in Ktau) for k in range(6)], cells_mod=res2)
    json.dump(out, open(f'{S4}/p2/out/switch_s{slab}.json', 'w'), indent=1, default=str)
    print('raw:', {k: v for k, v in raw.items() if k != 'meta'}); print('modulated:', out['modulated'])
    return out


if __name__ == '__main__':
    ck, slab, tca = sys.argv[1], int(sys.argv[2]), float(sys.argv[3])
    HMAX = float(sys.argv[4]) if len(sys.argv) > 4 else 0.25; K = int(sys.argv[5]) if len(sys.argv) > 5 else 30
    certify(ck, slab, tca, HMAX, K)
