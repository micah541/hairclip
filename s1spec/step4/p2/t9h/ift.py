"""T9 H1 part: alpha by the step-3 IFT for the base at a t_c-ball of slab 180 (constraints <chi w, f_j>_rho = 0 of step 3,
f_j from the step-3 IBP hull boxes, psi3^/h2^ from c4fun, S1 from s1chart).  Linear IFT about a float reference
(xc_f, lam_f, alpha_f): delta^T M = -G with G_j = <chi w0, f_j>, M_ij = <mode_i, g_j> (modes: translation, log-scale, alpha;
exact first-order expressions with interval slack), plus a second-order allowance.  Outputs alpha, xc, lam as balls."""
import sys, os, json, time, math
sys.path.insert(0, '/home/claude/s1spec/step4/p2/t9h'); sys.path.insert(0, '/home/claude/s1spec/step3nl/agent_mod/ibp')
from flint import arb, acb, arb_mat, ctx
ctx.prec = 128
import basegeo, s1chart
import ibp_boxes as IB
import c4fun as C4
U = lambda v: math.nextafter(float(abs(v).upper()), math.inf)
I = acb(0, 1)
TC = arb(1 - 1 / 1024, 1 / 1024)                    # t_c in [1 - 1/512, 1]: tau in [27.15 - 1.46e-4, 27.15]
if os.environ.get('T9TC'): TC = arb(*os.environ['T9TC'].split(','))   # a t_c ball 'mid,rad' (slab-wide hand-off, review A1)
REF = dict(xc=1.0079131539037775e-4, lam=0.9999950987042785, alpha=-1.1305273069129123e-4)   # float fit at k = 23700

def setup():
    t0 = time.time()
    g = basegeo.Geo(180, 40.0, TC); S = s1chart.S1()
    print(f'geometry + S1 chart ready ({time.time() - t0:.0f}s); Theta = {S.Theta.mid()} +- {S.Theta.rad()}', flush=True)
    return g, S

def foot(g, xb, Pfun, lam):
    """foot point on the base: element i and an s-ball containing the root of Re(Delta) (interval Newton, verified)"""
    for _ in range(40):
        i = g.elem(xb); h = float(g.hh[i].mid()); s = 2 * (xb - g.xe_mid[i]) / h - 1
        q = g.at(i, arb(s)); D = Pfun(q['Gam'], q['th'])
        c = math.cos(float(q['th'].mid()) - float(Pfun.th1.mid()))
        dxb = -float(D.real.mid()) * float(lam.mid()) / c; xb += dxb
        if abs(dxb) < 1e-14: break
    i = g.elem(xb); h = float(g.hh[i].mid()); s = 2 * (xb - g.xe_mid[i]) / h - 1
    S0 = arb(s); r0 = Pfun(g.at(i, S0)['Gam'], g.at(i, S0)['th']).real
    rad = 4 * U(r0) * float(lam.mid()) / h + 1e-15
    Sb = arb(s, rad)
    qb = g.at(i, Sb); dmin = (qb['th'] - Pfun.th1).cos() * g.hh[i] / (2 * lam)
    assert float(dmin.lower()) > 0
    Snew = S0 - r0 / dmin
    assert float(Snew.lower()) >= float(Sb.lower()) and float(Snew.upper()) <= float(Sb.upper()), 'foot point not enclosed'
    return i, Snew

class Pf:
    """Delta(Gamma~, theta~) = ((Gamma~ - xc)/lam - Gamma1) e^{-i theta1} (tip chart) or the ray-aligned form (far chart)"""
    def __init__(self, c, xc, lam, S, far=False):
        self.c, self.xc, self.lam, self.far, self.S = c, xc, lam, far, S
        self.th1 = c['th']
    def __call__(self, Gam, th):
        if not self.far:
            return ((Gam - self.xc) / self.lam - self.c['G']) * (-I * self.th1).exp()
        Dl = ((Gam - self.xc) / self.lam - self.S.G8) * (-I * self.S.Theta).exp() - self.c['tot']
        return Dl * (-I * self.c['D']).exp()

def ujets(g, S, x, xc, lam, far=False):
    """u, u', u'', u''' at the S1 point x (arb), and (i, S-ball) of the foot point"""
    c = S.far_chart(x) if far else S.tip_chart(x)
    P = Pf(c, xc, lam, S, far)
    i, Sb = foot(g, float((x * lam).mid()), P, lam)
    q = g.at(i, Sb)
    s0 = arb(Sb.mid()); q0 = g.at(i, s0); D0 = P(q0['Gam'], q0['th'])
    dIm = g.hh[i] / (2 * lam) * (q['th'] - P.th1).sin()            # d Im(Delta)/ds on the s-ball (mean-value form)
    u = D0.imag + arb(0, U(dIm) * U(Sb - s0))
    k1, k1p, k1pp = c['k'], c['k1'], c['k2']
    a = 1 - k1 * u; d = q['th'] - c['th']
    up = a * d.tan()
    J = (a * a + up * up).sqrt()
    kap = lam * q['k0']; kapx = lam * lam * J * q['k1']
    ap = -k1p * u - k1 * up
    F = kap * J ** 3 + up * ap - k1 * (a * a + up * up); upp = F / a
    Jp = (a * ap + up * upp) / J; app = -k1pp * u - 2 * k1p * up - k1 * upp
    Fp = kapx * J ** 3 + 3 * kap * J * J * Jp + upp * ap + up * app - k1p * (a * a + up * up) - k1 * (2 * a * ap + 2 * up * upp)
    uppp = (Fp * a - F * ap) / (a * a)
    return dict(u=u, u1=up, u2=upp, u3=uppp, i=i, S=Sb, q=q, c=c, d=d)

def run(g, S, xc, lam, al, verbose=True, limit=None):
    t0 = time.time(); xc, lam, al = arb(xc), arb(lam), arb(al)
    G = [arb(0)] * 3; M = [[arb(0)] * 3 for _ in range(3)]
    x8 = float(s1chart.TT.lower())
    nb = 0; worst = 0.0
    for b in IB.boxes():
        nb += 1
        if limit and nb > limit: break
        X = b['X']; xc_ = float(X.mid()); hX = arb(0, X.rad())
        gj = [b['g'][j][0] for j in (1, 2, 3)]
        far = float(X.upper()) > x8 - 1e-3
        if not far:
            J = ujets(g, S, arb(xc_), xc, lam)
            ps = C4.PSI.jets_any(arb(xc_), xc_, 5); hh = C4.H2.jets_any(arb(xc_), xc_, 5)
            w = [J['u'] - al * ps[0] - al * al * hh[0], J['u1'] - al * ps[1] - al * al * hh[1], J['u2'] - al * ps[2] - al * al * hh[2],
                 J['u3'] - al * ps[3] - al * al * hh[3]]
            # crude |u''''| over the box: |lam^3 kappa~''| J^2 + |k1''| + 10 (generous lower-order allowance)
            qb = g.at(J['i'], J['S']); u4 = U(lam ** 3 * qb['k2']) * 1.1 + 1.0 + 10.0
            w4 = arb(u4) + abs(al) * abs(b['ps'][4]) + al * al * abs(b['hh'][4])
            h2 = hX * hX; h3 = h2 * hX; h4 = h2 * h2                      # (arb ** on balls containing 0 gives nan)
            wbox = w[0] + w[1] * hX + w[2] * h2 / 2 + w[3] * h3 / 6 + arb(0, U(w4)) * h4 / 24
            dth = U(J['d']) + U(hX) * (U(J['u2']) + 1e-6)
        else:
            xs = max(xc_, x8 + 1e-9)
            J = ujets(g, S, arb(xs), xc, lam, far=True)
            u1b = U(J['u1']) + U(hX) * (U(J['q']['k0'] * lam) + U(J['c']['k']) + 1e-3)          # crude sup|u'| on the box
            wbox = J['u'] + arb(0, u1b * (U(hX) + abs(xs - xc_))) - al * b['ps'][0] - al * al * b['hh'][0]
            dth = U(J['d']) + U(hX) * 1e-3
        worst = max(worst, U(wbox) * U(b['w']) * U(gj[2]))
        for j in range(3): G[j] += b['w'] * wbox * gj[j]
        if not math.isfinite(float(G[2].rad())):
            print('terms', [v.str(3) for v in w] if not far else 'far', 'ps', [v.str(3) for v in ps] if not far else '', 'hh', [v.str(3) for v in hh] if not far else '', 'w4', w4.str(3) if not far else '', 'u4', u4 if not far else '')
            print('NaN at box', nb, 'J', {k: J[k].str(3) for k in ('u','u1','u2','u3')}, 'X', X.str(5), 'far', far, 'wbox', wbox.str(3), 'w', b['w'].str(3), 'g3', gj[2].str(3), 'ps0', b['ps'][0].str(3), 'u', J['u'].str(3)); raise SystemExit
        # modes (first order, exact up to the stated slack): translation, log-scale, alpha
        tdt = arb(0, U(arb(dth).tan()))
        m_tr = (-b['b2'][0] + arb(0, 1) * tdt) / lam
        m_sc = -b['b1'][0] - arb(0, U(J['u'])) + arb(0, U(b['q'])) * tdt
        m_al = -(b['ps'][0] + 2 * al * b['hh'][0])
        for i_, m_ in enumerate((m_tr, m_sc, m_al)):
            for j in range(3): M[i_][j] += b['w'] * m_ * gj[j]
        if verbose and nb % 500 == 0: print(f'  box {nb} ({time.time() - t0:.0f}s) G3 so far {G[2].str(5)}', flush=True)
    return G, M, nb

if __name__ == '__main__':
    g, S = setup()
    G, M, nb = run(g, S, REF['xc'], REF['lam'], REF['alpha'])
    Q = arb(0, 1e-12)                                                    # second-order allowance (|delta|^2 C, C <= 1e4)
    Mt = arb_mat([[M[i][j] for i in range(3)] for j in range(3)])        # rows j: sum_i delta_i M_ij = -G_j
    rhs = arb_mat([[-(G[j] + Q)] for j in range(3)])
    dl = Mt.solve(rhs)
    al = arb(REF['alpha']) + dl[2, 0]; xc = arb(REF['xc']) + dl[0, 0]; lam = arb(REF['lam']) * dl[1, 0].exp()
    out = dict(boxes=nb, G=[g_.str(6) for g_ in G], M=[[M[i][j].str(5) for j in range(3)] for i in range(3)],
               delta=[dl[k, 0].str(6) for k in range(3)], alpha=[float(al.lower()), float(al.upper())], xc=[float(xc.lower()), float(xc.upper())],
               lam=[float(lam.lower()), float(lam.upper())], alpha_rel_halfwidth=float(al.rad() / abs(al.mid())), Theta=[float(S.Theta.lower()), float(S.Theta.upper())],
               tc_ball=[1 - 1 / 512, 1.0], slab=180)
    json.dump(out, open('/home/claude/s1spec/step4/p2/out/t9_ift.json', 'w'), indent=1)
    print(json.dumps(out)[:1500])
