"""H1 uniformly on a t_c ball B of slab 180 (review A1: the good-time window is the whole slab): sup_B ||dw/dt_c||_Lam.
w = u - alpha psi3^ - alpha^2 h2^, p = (xc, log lam, alpha) from the step-3 IFT.  Along B (implicit function theorem):
  dw/dt = u_t + pdot_xc u_xc + pdot_l u_l - pdot_a (psi3^ + 2 alpha h2^),   pdot = -(M^T)^-1 G_t,
where u_t, u_xc, u_l are the exact first variations of the normal graph u and of u', u'', u''' (dual numbers through the h1 formulas)
under the base motion at fixed s (T6 bivariate polynomials, d/dt_c), the translation Z = -1/lam and the scaling Z = -Gamma~.
M_ij = <chi m_i, f_j>, G_t,j = <chi u_t, f_j> on the 3425 IBP boxes (third-order Taylor in x; the third derivative on the box hull).
The parameter box is P_B = p(t_k) + [-Rp, Rp]; the run asserts r_B |pdot| <= Rp - rad p(t_k), so p(t) stays in P_B on B.
Then for t in B:  sqrt Lam(t) <= sqrt Lam(t_k) + |t - t_k| sup_B ||w_t||_Lam   (sqrt Lam is a seminorm; mean value along B).
usage: python3 -I h1dt.py TC RB H1_JSON_AT_TC [RP_FACTOR]     (env COHERENT, S1PROF, S1TT as for h1.py)"""
import sys, os, json, time, math
TC, RB, HJ = sys.argv[1], float(sys.argv[2]), sys.argv[3]; RPF = float(sys.argv[4]) if len(sys.argv) > 4 else 3.0
sys.argv = ['h1.py', 'm']
sys.path.insert(0, '/home/claude/s1spec/step4/p2/t9h')
import h1
from flint import arb, acb, arb_mat, arb_poly
U, I = h1.U, h1.I

class D:
    """dual number with 3 directions (t_c, xc, log lam); value and derivatives arb or acb"""
    __slots__ = ('v', 'd')
    def __init__(s, v, d=(0, 0, 0)): s.v = v; s.d = list(d)
    def __add__(s, o):
        if isinstance(o, D): return D(s.v + o.v, [a + b for a, b in zip(s.d, o.d)])
        return D(s.v + o, s.d)
    __radd__ = __add__
    def __neg__(s): return D(-s.v, [-a for a in s.d])
    def __sub__(s, o): return s + (-o)
    def __rsub__(s, o): return (-s) + o
    def __mul__(s, o):
        if isinstance(o, D): return D(s.v * o.v, [a * o.v + s.v * b for a, b in zip(s.d, o.d)])
        return D(s.v * o, [a * o for a in s.d])
    __rmul__ = __mul__
    def inv(s):
        r = 1 / s.v; return D(r, [-a * r * r for a in s.d])
    def __truediv__(s, o): return s * o.inv() if isinstance(o, D) else D(s.v / o, [a / o for a in s.d])
    def __rtruediv__(s, o): return s.inv() * o
    def sqrt(s):
        r = s.v.sqrt(); return D(r, [a / (2 * r) for a in s.d])
    def tan(s):
        t = s.v.tan(); return D(t, [a * (1 + t * t) for a in s.d])

def at_dt(g, i, s):
    """t_c-derivatives at fixed local s (over the t_c ball): theta, kappa, kappa_x, Gamma"""
    el = g.els[i]; T_ = g.T; sh = arb_poly([arb(T_.mid()), arb(1)]); dl = arb(0, T_.rad())
    tv = lambda tp: tp(sh)(dl); tvd = lambda tp: tp.derivative()(sh)(dl)
    Gr, Gi = el['G'].at(s); Pr, Pi = el['Pg'].at(s)
    G = acb(tv(Gr), tv(Gi)); P = acb(tv(Pr), tv(Pi)); Gd = acb(tvd(Gr), tvd(Gi)); Pd = acb(tvd(Pr), tvd(Pi))
    return dict(th=tvd(el['th'].at(s)), k0=tvd(el['kap'].at(s)), k1=tvd(el['kx'].at(s)), Gam=Pd * G + P * Gd)

def ujets_d(g, S, x, xc, lam):
    """u, u', u'', u''' at S1 arclength x (exact point or ball) as duals: d[0] = d/dt_c, d[1] = d/dxc, d[2] = d/dlog lam"""
    far = float(x.upper()) > S.x13
    ch = S.far_chart(x) if far else (S.chart(x) if x.rad() == 0 else S.chart_ball(x))
    P = h1.Pf(ch, xc, lam, S, far)
    i, Sb = h1.foot(g, float((x * lam).mid()), P, lam)
    q = g.at(i, Sb); qt = at_dt(g, i, Sb); s0 = arb(Sb.mid()); q0 = g.at(i, s0); D0 = P(q0['Gam'], q0['th'])
    om = ch['th']; dIm = g.hh[i] / (2 * lam) * (q['th'] - om).sin()
    u0 = D0.imag + arb(0, U(dIm) * U(Sb - s0))
    d = q['th'] - om; td = d.tan(); hh2 = g.hh[i] / 2
    eom = (-I * om).exp(); Gt = (q['Gam'] - xc) / lam
    Z = [qt['Gam'] / lam * eom, -eom / lam, -Gt * eom]                 # rotated first variations of Gamma~ at fixed s
    sp_ = [-(z.real) / (hh2 / lam * d.cos()) for z in Z]               # foot-point motion s' (foot condition Re P = 0)
    U0 = D(u0, [z.imag - td * z.real for z in Z])
    TH = D(q['th'], [qt['th'] + hh2 * q['k0'] * sp_[0]] + [hh2 * q['k0'] * v for v in sp_[1:]])
    K0 = D(q['k0'], [qt['k0'] + hh2 * q['k1'] * sp_[0]] + [hh2 * q['k1'] * v for v in sp_[1:]])
    K1 = D(q['k1'], [qt['k1'] + hh2 * q['k2'] * sp_[0]] + [hh2 * q['k2'] * v for v in sp_[1:]])
    LAM = D(lam, [0, 0, lam])
    k1c, k1p, k1pp = ch['k'], ch['k1'], ch['k2']
    a = 1 - k1c * U0; u1 = a * (TH - om).tan()
    J = (a * a + u1 * u1).sqrt(); kap = LAM * K0; kapx = LAM * LAM * J * K1
    ap = -k1p * U0 - k1c * u1
    F = kap * J * J * J + u1 * ap - k1c * (a * a + u1 * u1); u2 = F / a
    Jp = (a * ap + u1 * u2) / J; app = -k1pp * U0 - 2 * k1p * u1 - k1c * u2
    Fp = kapx * J * J * J + 3 * kap * J * J * Jp + u2 * ap + u1 * app - k1p * (a * a + u1 * u1) - k1c * (2 * a * ap + 2 * u1 * u2)
    u3 = (Fp * a - F * ap) / (a * a)
    return [U0, u1, u2, u3]

if __name__ == '__main__':
    t0 = time.time()
    HK = json.load(open(HJ)); tcB = arb(TC, RB)
    pk = [arb(sum(HK[k]) / 2, (HK[k][1] - HK[k][0]) / 2) for k in ('xc', 'lam', 'alpha')]
    drift = [float(v) for v in os.environ.get('PDOT', '1e-5,2e-6,1e-5').split(',')]   # expected |pdot| (per t_c); checked below
    TREF = float(os.environ.get('TREF', TC)); dist = abs(float(TC) - TREF) + RB         # p is known (Krawczyk) at TREF
    rp = [float(pk[0].rad()) + RPF * dist * drift[0] + 1e-13, float(pk[1].rad()) / float(pk[1].mid()) + RPF * dist * drift[1] + 1e-13,
          float(pk[2].rad()) + RPF * dist * drift[2] + 1e-13]
    xcB = arb(pk[0].mid(), rp[0]); lamB = arb(pk[1].mid(), rp[1] * (1 + rp[1]) * float(pk[1].mid())); alB = arb(pk[2].mid(), rp[2])
    g = h1.basegeo.Geo(180, 64.0, tcB); S = h1.s1chart2.S1(); S.x13 = float((h1.s1chart2.T0 + 5).lower()) - 1e-6
    print(f'ball t_c {TC} +/- {RB}; P_B radii {rp} ({time.time() - t0:.0f}s)', flush=True)
    # ---- IBP boxes: the exact Jacobian M (modes = first variations) and G_t; defect correction for pdot
    M = [[arb(0)] * 3 for _ in range(3)]; Gt = [arb(0)] * 3; store = []
    for b in h1.IB.boxes():
        X = b['X']; xm = float(X.mid()); hX = arb(0, X.rad()); gj = [b['g'][j][0] for j in (1, 2, 3)]
        xs = arb(max(xm, 0.0)); Xb = arb(max(xm, 0.0), float(X.rad()))
        jd = ujets_d(g, S, xs, xcB, lamB); jh = ujets_d(g, S, Xb, xcB, lamB)
        ps = h1.C4.PSI.jets_any(xs, xm, 5); hh = h1.C4.H2.jets_any(xs, xm, 5); psb = h1.C4.PSI.jets_any(Xb, xm, 5); hhb = h1.C4.H2.jets_any(Xb, xm, 5)
        h2 = hX * hX; h3 = h2 * hX
        cols = [[jd[k].d[dr] for k in range(3)] + [jh[3].d[dr]] for dr in range(3)] + [[-(ps[k] + 2 * alB * hh[k]) for k in range(3)] + [-(psb[3] + 2 * alB * hhb[3])]]
        box = [c[0] + c[1] * hX + c[2] * h2 / 2 + c[3] * h3 / 6 for c in cols]       # u_t, u_xc, u_l, d w/d alpha on the box
        for j in range(3):
            Gt[j] += b['w'] * box[0] * gj[j]
            for i_ in range(3): M[i_][j] += b['w'] * box[i_ + 1] * gj[j]
        store.append((b['w'], gj, hX, h2, h3, cols))
    Mt = arb_mat([[M[i_][j] for i_ in range(3)] for j in range(3)])
    pd0 = Mt.solve(arb_mat([[-Gt[j]] for j in range(3)])); q = [arb(pd0[i_, 0].mid()) for i_ in range(3)]
    # defect: G~ = <chi (u_t + sum q_i m_i), f_j> has a small integrand, hence a small quadrature width; pdot = q + delta
    Gd = [arb(0)] * 3
    for (wb, gj, hX, h2, h3, cols) in store:
        c = [cols[0][k] + q[0] * cols[1][k] + q[1] * cols[2][k] + q[2] * cols[3][k] for k in range(4)]
        bx = c[0] + c[1] * hX + c[2] * h2 / 2 + c[3] * h3 / 6
        for j in range(3): Gd[j] += wb * bx * gj[j]
    for it_ in range(2):                                                # two defect corrections
        if it_: q = [arb(v.mid()) for v in pdv]; Gd = [arb(0)] * 3
        if it_:
            for (wb, gj, hX, h2, h3, cols) in store:
                c = [cols[0][k] + q[0] * cols[1][k] + q[1] * cols[2][k] + q[2] * cols[3][k] for k in range(4)]
                bx = c[0] + c[1] * hX + c[2] * h2 / 2 + c[3] * h3 / 6
                for j in range(3): Gd[j] += wb * bx * gj[j]
        dl_ = Mt.solve(arb_mat([[-Gd[j]] for j in range(3)]))
        pdv = [q[i_] + dl_[i_, 0] for i_ in range(3)]
    print(f'  G_t = {[v.str(3) for v in Gt]}; defect {[v.str(3) for v in Gd]}; pdot = {[v.str(4) for v in pdv]} ({time.time() - t0:.0f}s)', flush=True)
    dlam_rel = pdv[1]                                                   # d log lam / dt
    need = [U(abs(pdv[0])) * dist, U(abs(dlam_rel)) * dist, U(abs(pdv[2])) * dist]
    room = [rp[0] - float(pk[0].rad()), rp[1] - float(pk[1].rad()) / float(pk[1].mid()), rp[2] - float(pk[2].rad())]
    assert all(need[i_] <= room[i_] for i_ in range(3)), ('parameter box too small', need, room)
    # ---- Lambda norm of w_t on cells (coarse dyadic grid; third-order Taylor, third derivative on the cell hull)
    A2 = B2 = V02 = V22 = arb(0); sig_tip = 8 - h1.s1chart2.TT; ncell = 0
    grid = [(k / 64, 1 / 64) for k in range(656)] + [(10.25 + k / 32, 1 / 32) for k in range(768)] + [(34.25 + k / 16, 1 / 16) for k in range(412)]
    for x, h in grid:
        xm = x + h / 2; xs = arb(xm); Xb = arb(xm, h / 2)
        jd = ujets_d(g, S, xs, xcB, lamB); jh = ujets_d(g, S, Xb, xcB, lamB)
        ps = h1.C4.PSI.jets_any(xs, xm, 5); hh = h1.C4.H2.jets_any(xs, xm, 5); psb = h1.C4.PSI.jets_any(Xb, xm, 5); hhb = h1.C4.H2.jets_any(Xb, xm, 5)
        wt = [jd[j].d[0] + pdv[0] * jd[j].d[1] + pdv[1] * jd[j].d[2] - pdv[2] * (ps[j] + 2 * alB * hh[j]) for j in range(3)]
        w3 = jh[3].d[0] + pdv[0] * jh[3].d[1] + pdv[1] * jh[3].d[2] - pdv[2] * (psb[3] + 2 * alB * hhb[3])
        t = arb(0, h / 2); t2 = t * t; t3 = t2 * t
        W0 = wt[0] + wt[1] * t + wt[2] * t2 / 2 + w3 * t3 / 6; W1 = wt[1] + wt[2] * t + w3 * t2 / 2; W2 = wt[2] + w3 * t; W3 = w3
        Xs = arb(x).union(arb(x + h))
        if x < 34.25:
            sig = Xs + sig_tip; ch0, ch1, ch2 = h1.smooth((32 - sig) / 3); ch1 = -ch1 / 3; ch2 = ch2 / 9
            rho = (-(1 + arb(x) * arb(x)) ** (arb(2) / 3) / 4).exp()
            sa, sb = arb(U(ch0 * W0)), arb(U(ch2 * W0 + 2 * ch1 * W1 + ch0 * W2))
            A2 += arb(h) * rho * sa * sa; B2 += arb(h) * rho * sb * sb
        if x + h > 10.25:
            za = h1.smooth((Xs - arb('10.25')) / 6); zb = h1.smooth((60 - Xs) / 3)
            z0 = za[0] * zb[0]; z1 = za[1] / 6 * zb[0] - za[0] * zb[1] / 3; z2 = za[2] / 36 * zb[0] - 2 * za[1] / 6 * zb[1] / 3 + za[0] * zb[2] / 9
            q1 = 1 + arb(x) * arb(x); Wt_ = 1 / (q1 * q1)
            s0_, s2_ = arb(U(z0 * W1)), arb(U(z2 * W1 + 2 * z1 * W2 + z0 * W3))
            V02 += arb(h) * Wt_ * s0_ * s0_; V22 += arb(h) * Wt_ * s2_ * s2_
        assert math.isfinite(float((A2 + B2 + V02 + V22).rad())), ('non-finite cell', x)
        ncell += 1
    L2 = A2 + h1.EPS * B2 + h1.MU * (V02 + h1.EPSI * V22)
    out = dict(args=dict(TC=TC, RB=RB, H1_JSON=HJ, RP_FACTOR=RPF, drift=drift, env={k: os.environ[k] for k in ('COHERENT', 'S1PROF', 'S1TT', 'H1TAG', 'H1NSUB', 'H1REF', 'TREF', 'PDOT', 'T9TC', 'T9IFT', 'WINOUT', 'SUBW', 'T9JSON', 'TOPLOSS', 'N3JSON') if k in os.environ}), tc=[float(tcB.lower()), float(tcB.upper())], RB=RB, TREF=TREF, dist=dist, P_B=dict(xc=[float(xcB.lower()), float(xcB.upper())], lam=[float(lamB.lower()), float(lamB.upper())],
               alpha=[float(alB.lower()), float(alB.upper())]), pdot=[[float(v.lower()), float(v.upper())] for v in pdv], cells=ncell,
               A=U(A2.sqrt()), B=U(B2.sqrt()), V0=U(V02.sqrt()), V2=U(V22.sqrt()), wt_norm=U(L2.sqrt()), lip_term=math.nextafter(U(L2.sqrt()) * dist, math.inf),
               cpu=time.time() - t0)
    json.dump(out, open(f'/home/claude/s1spec/step4/p2/out/h1slab/dt_{TC}_{RB}.json', 'w'), indent=1); print(json.dumps(out))
