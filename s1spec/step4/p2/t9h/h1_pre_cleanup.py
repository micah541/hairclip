"""T9 H1 with the tight S1 chart (prof_tight + ext813): for t_c sub-ball K (K = 0..7 of [1 - 1/512, 1] in slab 180):
(1) the step-3 IFT (constraints <chi w, f_j>_rho = 0 on the 3425 step-3 IBP boxes), two linear passes;
(2) Lambda = A^2 + eps B^2 + mu (V0^2 + epsI V2^2) of w = u - alpha psi3^ - alpha^2 h2^ at the IFT solution (step-3 definitions:
    chi = chi_32, zeta = S((x - 10.25)/6) S((60 - x)/3), rho = exp(-(1 + x^2)^(2/3)/4), W <= (1 + x^2)^-2), by cells with
    third-order Taylor forms and an exact u'''' (sympy) on cell hulls.  usage: python3 -I h1.py K"""
import sys, json, time, math, os
sys.path.insert(0, '/home/claude/s1spec/step4/p2/t9h'); sys.path.insert(0, '/home/claude/s1spec/step3nl/agent_mod/ibp')
import sympy as sp
from flint import arb, acb, arb_mat, ctx
ctx.prec = 128
import basegeo, s1chart2
import ibp_boxes as IB
import c4fun as C4
from ift import foot, REF
U = lambda v: math.nextafter(float(abs(v).upper()), math.inf)
I = acb(0, 1)
KA = sys.argv[1] if len(sys.argv) > 1 else '0'
NSUB = int(os.environ.get('H1NSUB', '8'))                  # W' = [1 - 1/512, 1] covered by NSUB t_c balls
RADK = 1.0 / (1024 * NSUB)
if KA == 'm': K = 'm'; TCK = arb(1 - 1.0 / 1024)            # the midpoint of W' = [1 - 1/512, 1] (exact point)
elif KA.startswith('p'): K = KA; TCK = arb(KA[1:])          # an exact t_c point of slab 180, e.g. p-0.5 (t_c = -0.5)
else: K = int(KA); TCK = arb(1 - (2 * K + 1) * RADK, RADK)
EPS, EPSI, MU = arb('1e-5'), arb('3e-5'), arb('0.08')

# ---- exact u'''' from the normal-graph curvature identity (sympy)
u = sp.symbols('u0:5'); Kg = sp.symbols('K0:3'); c = sp.symbols('c0:4')
def ddx(e):
    out = sum(sp.diff(e, u[i]) * u[i + 1] for i in range(4)) + sum(sp.diff(e, Kg[i]) * Kg[i + 1] for i in range(2))
    return out + sum(sp.diff(e, c[i]) * c[i + 1] for i in range(3))
a0 = 1 - c[0] * u[0]; a1 = -c[1] * u[0] - c[0] * u[1]; Jx = sp.sqrt(a0 ** 2 + u[1] ** 2)
U2 = (Kg[0] * Jx ** 3 + u[1] * a1 - c[0] * (a0 ** 2 + u[1] ** 2)) / a0
U4 = ddx(ddx(U2))
from sympy.printing.pycode import PythonCodePrinter
class SafePrinter(PythonCodePrinter):
    """integer and half-integer powers as products (arb ** on balls containing 0 gives nan)"""
    def _print_Pow(self, expr, rational=False):
        b, e = expr.base, expr.exp
        if e.is_Rational and (e.q in (1, 2)):
            n = abs(int(e.p)); unit = '(' + self._print(b) + ')' if e.q == 1 else 'sqrt(' + self._print(b) + ')'
            prod = '*'.join([unit] * n) if n > 0 else '1'
            return '(' + prod + ')' if e > 0 else '(1/(' + prod + '))'
        return super()._print_Pow(expr)
mods = [{'sqrt': lambda a: a.sqrt()}, 'math']
fU4 = sp.lambdify(list(u[:4]) + list(Kg) + list(c), U4, modules=mods, printer=SafePrinter)

class Pf:
    def __init__(self, ch, xc, lam, S, far):
        self.ch, self.xc, self.lam, self.S, self.far = ch, xc, lam, S, far; self.th1 = ch['th']
    def __call__(self, Gam, th):
        if not self.far: return ((Gam - self.xc) / self.lam) * (-I * self.th1).exp() - self.ch['qp']     # theta1 enters once
        Dl = ((Gam - self.xc) / self.lam - self.S.G13) * (-I * self.S.Theta).exp() - self.ch['tot']
        return Dl * (-I * self.ch['D']).exp()

def ujets(g, S, x, xc, lam):
    far = float(x.upper()) > S.x13
    ch = S.far_chart(x) if far else S.chart(x)
    P = Pf(ch, xc, lam, S, far)
    i, Sb = foot(g, float((x * lam).mid()), P, lam)
    q = g.at(i, Sb); s0 = arb(Sb.mid()); q0 = g.at(i, s0); D0 = P(q0['Gam'], q0['th'])
    dIm = g.hh[i] / (2 * lam) * (q['th'] - P.th1).sin()
    u0 = D0.imag + arb(0, U(dIm) * U(Sb - s0))
    k1, k1p, k1pp = ch['k'], ch['k1'], ch['k2']
    a = 1 - k1 * u0; d = q['th'] - ch['th']; u1 = a * d.tan()
    J = (a * a + u1 * u1).sqrt(); kap = lam * q['k0']; kapx = lam * lam * J * q['k1']
    ap = -k1p * u0 - k1 * u1
    F = kap * J * J * J + u1 * ap - k1 * (a * a + u1 * u1); u2 = F / a
    Jp = (a * ap + u1 * u2) / J; app = -k1pp * u0 - 2 * k1p * u1 - k1 * u2
    Fp = kapx * J * J * J + 3 * kap * J * J * Jp + u2 * ap + u1 * app - k1p * (a * a + u1 * u1) - k1 * (2 * a * ap + 2 * u1 * u2)
    u3 = (Fp * a - F * ap) / (a * a)
    return dict(u=[u0, u1, u2, u3], i=i, S=Sb, q=q, ch=ch, d=d, far=far, J=J)

def u4_cell(g, S, Jc, h, lam):
    """hull of u'''' over the cell |x - x_c| <= h (exact expression on hulls)"""
    uu = [Jc['u'][j] + arb(0, h) * (U(Jc['u'][j + 1]) + 1.0) for j in range(3)] + [Jc['u'][3] + arb(0, h)]
    i = Jc['i']; sw = h * 2 / float(g.hh[i].mid()) * 1.01 + 1e-12
    qb = g.at(i, arb(Jc['S'].mid(), U(Jc['S'] - Jc['S'].mid()) + sw))
    ch = Jc['ch']
    if Jc['far']: k1h, k1ph = ch['k'], ch['k1']
    else:
        tc0 = s1chart2.T0 - arb(float(Jc['x'].mid())); H0 = S.hull(float((tc0 - h).lower()) - 1e-9, float((tc0 + h).upper()) + 1e-9)
        k1h, k1ph = H0[2], H0[3]
    # review B5: exact chain rule on the cell hull (no typed slack): sigma-derivatives of the base curvature at the foot,
    # K1 = lam^2 J k_x, K2 = lam^2 J' k_x + lam^3 J^2 k_xx, J = sqrt((1 - k1 u)^2 + u'^2), J' = (a a' + u' u'')/J, a' = -k1' u - k1 u'
    a_ = 1 - k1h * uu[0]; Jh = (a_ * a_ + uu[1] * uu[1]).sqrt(); aph = -k1ph * uu[0] - k1h * uu[1]; Jph = (a_ * aph + uu[1] * uu[2]) / Jh
    assert U(Jh) <= 1.01                                       # the foot s-range sw above assumes dx/dsigma <= 1.01
    K0 = lam * qb['k0']; K1 = lam * lam * Jh * qb['k1']
    K2 = lam * lam * Jph * qb['k1'] + lam * lam * lam * Jh * Jh * qb['k2']
    if Jc['far']:
        cc = [ch['k'], ch['k1'], ch['k2'], ch['k3']]
    else:
        tc_ = s1chart2.T0 - arb(float(Jc['x'].mid())); Hs = S.hull(float((tc_ - h).lower()) - 1e-9, float((tc_ + h).upper()) + 1e-9)
        p_, q_, k_, m_ = Hs[0], Hs[1], Hs[2], Hs[3]
        cc = [k_, m_, p_ / 4, -k_ * q_ / 4]
    u4 = fU4(*uu, K0, K1, K2, *cc)
    assert U(u4) < 1.0, ('u4 bootstrap', U(u4))              # the hulls uu assume |u''''| <= 1 on the cell (continuity in the cell radius)
    return u4

def smooth(uv):
    """S, S', S'' on an interval (clipped to [0, 1])"""
    lo, hi = float(uv.lower()), float(uv.upper())
    def poly(z):
        z2 = z * z; z3 = z2 * z; y = 1 - z
        return [z2 * z2 * (35 - 84 * z + 70 * z2 - 20 * z3), 140 * z3 * y * y * y, 420 * z2 * y * y * (1 - 2 * z)]
    if hi <= 0: return [arb(0)] * 3
    if lo >= 1: return [arb(1), arb(0), arb(0)]
    z = arb(max(lo, 0.0)).union(arb(min(hi, 1.0))); vals = poly(z)
    if lo < 0: vals = [vals[0].union(arb(0)), vals[1].union(arb(0)), vals[2].union(arb(0))]
    if hi > 1: vals = [vals[0].union(arb(1)), vals[1].union(arb(0)), vals[2].union(arb(0))]
    return vals

def ift_pass(g, S, xc, lam, al):
    G = [arb(0)] * 3; M = [[arb(0)] * 3 for _ in range(3)]
    for b in IB.boxes():
        X = b['X']; xc_ = float(X.mid()); hX = arb(0, X.rad()); gj = [b['g'][j][0] for j in (1, 2, 3)]
        xs = arb(max(xc_, 0.0))
        Jc = ujets(g, S, xs, xc, lam); Jc['x'] = xs
        ps = C4.PSI.jets_any(xs, xc_, 5); hh = C4.H2.jets_any(xs, xc_, 5)
        w = [Jc['u'][j] - al * ps[j] - al * al * hh[j] for j in range(4)]
        u4 = fU4_safe(g, S, Jc, U(hX), lam)
        w4 = abs(u4) + abs(al) * abs(b['ps'][4]) + al * al * abs(b['hh'][4])
        h2 = hX * hX; h3 = h2 * hX; h4 = h2 * h2
        wbox = w[0] + w[1] * hX + w[2] * h2 / 2 + w[3] * h3 / 6 + arb(0, U(w4)) * h4 / 24
        for j in range(3): G[j] += b['w'] * wbox * gj[j]
        uu = Jc['u']; u4b = arb(0, U(u4))
        U0b = uu[0] + uu[1] * hX + uu[2] * h2 / 2 + uu[3] * h3 / 6 + u4b * h4 / 24; U1b = uu[1] + uu[2] * hX + uu[3] * h2 / 2 + u4b * h3 / 6
        tdt = arb(0, U(abs(U1b) / (1 - abs(b['k']) * abs(U0b))))     # |tan d| on the box (u' = (1 - k1 u) tan d)
        modes = ((-b['b2'][0] + tdt) / lam, -b['b1'][0] + arb(0, U(U0b)) + arb(0, U(b['q'])) * tdt,
                 -(b['ps'][0] + 2 * al * b['hh'][0]))
        for i_, m_ in enumerate(modes):
            for j in range(3): M[i_][j] += b['w'] * m_ * gj[j]
    return G, M

def fU4_safe(g, S, Jc, h, lam):
    return u4_cell(g, S, Jc, h, lam)                       # review B3: no fallback, any failure is a hard error

def solve(G, M, Q):
    Mt = arb_mat([[M[i][j] for i in range(3)] for j in range(3)])
    rhs = arb_mat([[-(G[j] + arb(0, Q))] for j in range(3)])
    return Mt.solve(rhs)

def Lambda(g, S, xc, lam, al):
    st = S.st_tip if hasattr(S, 'st_tip') else None
    A2 = B2 = V02 = V22 = arb(0); sig_tip = 8 - s1chart2.TT
    out = dict(cells=0)
    grid = [(k / 128, 1 / 128) for k in range(1312)] + [(10.25 + k / 64, 1 / 64) for k in range(1536)] + [(34.25 + k / 32, 1 / 32) for k in range(824)]
    for x, h in grid:                                              # exact dyadic cells tiling [0, 60] (review B9)
        xm = x + h / 2; xs = arb(xm); hb = arb(0, h / 2)
        Jc = ujets(g, S, xs, xc, lam); Jc['x'] = xs
        ps = C4.PSI.jets_any(xs, xm, 5); hh = C4.H2.jets_any(xs, xm, 5)
        Xb = arb(xm, h / 2); psb = C4.PSI.jets_any(Xb, xm, 5); hhb = C4.H2.jets_any(Xb, xm, 5)
        w = [Jc['u'][j] - al * ps[j] - al * al * hh[j] for j in range(4)]
        w4 = abs(fU4_safe(g, S, Jc, h / 2, lam)) + abs(al) * abs(psb[4]) + al * al * abs(hhb[4]); w4 = arb(0, U(w4))
        t = hb; t2 = t * t; t3 = t2 * t; t4 = t2 * t2
        W0 = w[0] + w[1] * t + w[2] * t2 / 2 + w[3] * t3 / 6 + w4 * t4 / 24
        W1 = w[1] + w[2] * t + w[3] * t2 / 2 + w4 * t3 / 6
        W2 = w[2] + w[3] * t + w4 * t2 / 2
        W3 = w[3] + w4 * t
        Xs = arb(x).union(arb(x + h))
        if x < 34.25:
            sig = Xs + sig_tip; ch0, ch1, ch2 = smooth((32 - sig) / 3); ch1 = -ch1 / 3; ch2 = ch2 / 9
            rho = (-(1 + arb(x) * arb(x)) ** (arb(2) / 3) / 4).exp()       # rho decreasing in x
            a_ = ch0 * W0; b_ = ch2 * W0 + 2 * ch1 * W1 + ch0 * W2
            sa, sb = arb(U(a_)), arb(U(b_))                               # sup|.| (upper bounds; squares of balls via sups)
            A2 += arb(h) * rho * sa * sa
            B2 += arb(h) * rho * sb * sb
        if x + h > 10.25:
            za = smooth((Xs - arb('10.25')) / 6); zb = smooth((60 - Xs) / 3)
            z0 = za[0] * zb[0]; z1 = za[1] / 6 * zb[0] - za[0] * zb[1] / 3; z2 = za[2] / 36 * zb[0] - 2 * za[1] / 6 * zb[1] / 3 + za[0] * zb[2] / 9
            q1 = 1 + arb(x) * arb(x); Wt = 1 / (q1 * q1)                   # W <= (1 + x^2)^-2, decreasing
            v0 = z0 * W1; v2 = z2 * W1 + 2 * z1 * W2 + z0 * W3
            s0_, s2_ = arb(U(v0)), arb(U(v2))
            V02 += arb(h) * Wt * s0_ * s0_
            V22 += arb(h) * Wt * s2_ * s2_
        if not math.isfinite(float((A2 + B2 + V02 + V22).rad())):
            print('NaN cell', x, 'w', [v.str(3) for v in w], 'w4', w4.str(3), 'W0', W0.str(3), 'W2', W2.str(3)); raise SystemExit
        out['cells'] += 1
    Lam = A2 + EPS * B2 + MU * (V02 + EPSI * V22)
    out.update(A=U(A2.sqrt()), B=U(B2.sqrt()), V0=U(V02.sqrt()), V2=U(V22.sqrt()), sqrtLam=U(Lam.sqrt()))
    return out

if __name__ == '__main__':
    t0 = time.time()
    g = basegeo.Geo(180, 64.0, TCK); S = s1chart2.S1(); S.x13 = float((s1chart2.T0 + 5).lower()) - 1e-6
    print(f'K = {K}: t_c ball {TCK.str(10)}; setup {time.time() - t0:.0f}s', flush=True)
    ref = (arb(REF['xc']), arb(REF['lam']), arb(REF['alpha']))
    if os.environ.get('H1REF'): ref = tuple(arb(v) for v in os.environ['H1REF'].split(','))   # starting point (e.g. a neighbouring t_c)
    I3 = arb_mat([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    def pt(Mx): return arb_mat([[arb(Mx[i_, j_].mid()) for j_ in range(Mx.ncols())] for i_ in range(Mx.nrows())])
    for ps_ in range(5):                                           # Newton steps from the midpoint of G (floats); stop at the G-width level
        G, M = ift_pass(g, S, *ref)
        C = pt(pt(arb_mat([[M[i_][j_] for i_ in range(3)] for j_ in range(3)])).inv())
        st = [float((-(C * arb_mat([[arb(G[j_].mid())] for j_ in range(3)])))[i_, 0].mid()) for i_ in range(3)]
        noise = [sum(abs(float(C[i_, j_].mid())) * float(G[j_].rad()) for j_ in range(3)) for i_ in range(3)]
        print(f'  pass {ps_}: G = {[v.str(3) for v in G]}, step = {[f"{v:.3e}" for v in st]}, noise = {[f"{v:.2e}" for v in noise]} ({time.time() - t0:.0f}s)', flush=True)
        if all(abs(st[i_]) <= 4 * noise[i_] + 1e-15 for i_ in range(3)) or ps_ == 4: break
        ref = (arb(ref[0].mid() + st[0]), arb((ref[1] * arb(st[1]).exp()).mid()), arb(ref[2].mid() + st[2]))
    # Krawczyk (review B3) centred at p~ = p^ + step, delta = (dxc, dloglam, dalpha):  G(p~) in G(p^) + J(P) (p~ - p^)  (mean value, P ⊇ [p^, p~]),
    # K = -C G(p~) + (I - C J(P)) Delta, J(P) = M(P)^T, P = p~ + Delta.  K strictly inside Delta => exactly one zero in P, and it lies in p~ + K.
    Rk = [abs(st[i_]) + 4 * noise[i_] + 1e-15 for i_ in range(3)]; Dl = arb_mat([[arb(0, r)] for r in Rk])
    sv = arb_mat([[arb(v)] for v in st])
    pc = (ref[0] + sv[0, 0], ref[1] * sv[1, 0].exp(), ref[2] + sv[2, 0])
    Pb = (pc[0] + Dl[0, 0], pc[1] * Dl[1, 0].exp(), pc[2] + Dl[2, 0])
    _, MP = ift_pass(g, S, *Pb)
    JP = arb_mat([[MP[i_][j_] for i_ in range(3)] for j_ in range(3)])
    Gc = arb_mat([[G[j_]] for j_ in range(3)]) + JP * sv
    Kk = -(C * Gc) + (I3 - C * JP) * Dl
    kok = all(U(abs(Kk[i_, 0])) < Rk[i_] for i_ in range(3))
    print(f'  Krawczyk: {"ok" if kok else "FAIL"}; K = {[Kk[i_, 0].str(3) for i_ in range(3)]} in Delta = {[f"{r:.2e}" for r in Rk]} ({time.time() - t0:.0f}s)', flush=True)
    assert kok, 'Krawczyk test failed'
    xc = pc[0] + Kk[0, 0]; lam = pc[1] * Kk[1, 0].exp(); al = pc[2] + Kk[2, 0]
    L = Lambda(g, S, xc, lam, al)
    if S.C_rad > 0:                                                # rotation-mode bound: |dsqrtLam/dC| <= ||Pi q||_Lam
        gq = [arb(0)] * 3; N2 = {n: [arb(0)] * 2 for n in ('q', 'tr', 'sc', 'al')}
        for b in IB.boxes():
            k_, m_, p_, q_ = b['k'], b['m'], b['p'], b['q']; ch = b['chi']; c1 = arb(0, 1)
            for j in range(3): gq[j] += b['w'] * q_ * b['g'][j + 1][0]
            jets = dict(q=(q_, 1 + k_ * p_, m_ * p_ - k_ * k_ * q_), tr=(-b['b2'][0], k_ * c1, m_ * c1 + k_ * k_ * c1),
                        sc=(p_, -k_ * q_, -m_ * q_ - k_ * (1 + k_ * p_)), al=(b['ps'][0] + 2 * al * b['hh'][0], b['ps'][1] + 2 * al * b['hh'][1], b['ps'][2] + 2 * al * b['hh'][2]))
            for n, (g0, g1, g2) in jets.items():
                a_ = arb(U(ch[0] * g0)); bb = arb(U(ch[2] * g0 + 2 * ch[1] * g1 + ch[0] * g2))
                N2[n][0] += b['w'] * a_ * a_; N2[n][1] += b['w'] * bb * bb
        Mt = arb_mat([[M[i][j] for i in range(3)] for j in range(3)]); dp = Mt.solve(arb_mat([[gq[j]] for j in range(3)]))
        # E_I part (W <= (1+x^2)^-2 on [10.25, 60]): |g'| <= G1 and |(zeta g')''| <= G3 there (crude, generous): q' <= 1.01, tr', sc' <= 1e-3,
        # psi3^'+2a h2^' <= 1.6 x^1.03 (c4fun ext) -> sum W g'^2 <= g1max^2 * int_{10.25}^oo W; V2 with zeta'' <= 420/36, factor 50
        intW = arb(1) / (3 * arb('10.25') ** 3)
        gEI = dict(q=arb('1.01'), tr=arb('1e-3'), sc=arb('1e-3'), al=arb('1.6') * arb(60) ** arb('1.03'))
        nrm = {}
        for n in N2:
            e2 = N2[n][0] + EPS * N2[n][1] + MU * (gEI[n] * gEI[n] * intW) * (1 + EPSI * 2500)
            nrm[n] = e2.sqrt()
        dC = nrm['q'] + abs(dp[0, 0]) * nrm['tr'] / lam + abs(dp[1, 0]) * nrm['sc'] + abs(dp[2, 0]) * nrm['al']
        L['rot_mode_norm'] = U(dC); L['C_rad'] = S.C_rad; L['norms'] = {n: U(v) for n, v in nrm.items()}; L['dp'] = [U(dp[i, 0]) for i in range(3)]
        L['sqrtLam_total'] = L['sqrtLam'] + S.C_rad * U(dC) * 1.01
    L['krawczyk'] = dict(ok=kok, Delta=Rk, K=[U(abs(Kk[i, 0])) for i in range(3)])
    out = dict(args=dict(argv=sys.argv[1:], env={k: os.environ[k] for k in ('COHERENT', 'S1PROF', 'S1TT', 'H1TAG', 'H1NSUB', 'H1REF', 'TREF', 'PDOT', 'T9TC', 'T9IFT', 'WINOUT', 'SUBW', 'T9JSON', 'TOPLOSS', 'N3JSON') if k in os.environ}, grid='dyadic 1/128 on [0,10.25], 1/64 on [10.25,34.25], 1/32 on [34.25,60]'),
               K=K, tc=[float(TCK.lower()), float(TCK.upper())], alpha=[float(al.lower()), float(al.upper())], xc=[float(xc.lower()), float(xc.upper())],
               lam=[float(lam.lower()), float(lam.upper())], G=[v.str(4) for v in G], Theta=[float(S.Theta.lower()), float(S.Theta.upper())],
               Lambda=L, cpu=time.time() - t0)
    json.dump(out, open(f'/home/claude/s1spec/step4/p2/out/t9_h1_{K}' + ('_coh' if S.C_rad > 0 else '') + os.environ.get('H1TAG', '') + '.json', 'w'), indent=1)
    print(json.dumps(out))
