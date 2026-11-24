"""v3.3 task 6: copy of p2/t8se/h4loc.py (closed v3.2; read only) for the FIXED-L cone family of step 3 v3.3: unit L = 50, X = 3,
centres r_c in [300, 380] (rows of width 10), cylinder [r_c - 150, r_c + 150]; support [r_lo - 151, r_hi + 151] (h4loc.py: [r_lo/4 - 1,
7 r_hi/4 + 1], L = r_hi/4), M = 0.12, threshold M^2 = 0.0144.  Only cylinders() and the support line, the paths (HERE -> p2/t8se for
the imports; output default p2/out_v33/h4loc_v33.json) and the inputs' sha256 in args changed; the far family is not rerun (centres >= 490
are a subset of v3.2's certified family).  Env V33_VALIDATE=1: v3.2 cone rows, for a fidelity comparison in scratch.
Original docstring follows.
Review A1: Lemma O (H4), int (h0^(5))^2 dX <= M^2 per cylinder (cylinder units), from the good-time level-2 dissipation.

Localized level-2 dissipation on a cylinder support S = [x_lo, x_hi] (y_0 := x_lo - X0 - 1 >= 22, so T_2 = 28 Y_l^2 >= 28 y_0^2
on S and S lies in the core of every level):
    int_S w_4^2 <= D_S / (28 y_0^2),   D_S := int_S T_2 w_4^2 <= D_2,
and at the good time tau* of the window W (LEMMA_F C5, camp handoff):  D_2(tau*) <= [F_2(t_a) + int_W (2 g_2 sqrt F_2 + h_2)
+ inflow] / (theta_2 |W|).  A localized energy E_S = int T_2 psi^2 w_2^2 (psi = 1 on S) has the same leading term: its value at
t_a is bounded only by F_2(t_a), and its inflow (the ramp, or the transport through an inner cut-off) is that of F_2.  So the
window length, not the localization, decides; the localized forcing is <= the global one, which is used here.
In cylinder units (X = xi/L, h0 = h/L): int (h0^(5))^2 dX = L^7 int_{|xi| <= XL} (h^(5))^2 dxi, h the graph of the TRUE curve over the
cylinder's own centre tangent.  Chart (review O6b; derived and computed in arb as p2/t9chart.py and t9's KCH, no typed factor): with
gamma = Gamma~ + e N, x the base arclength, phi = theta~(x) - theta_c (theta_c the true centre tangent), xi' = M_t = (1 - k e) cos phi
- e' sin phi, and Q^t_4(x) := h^(5) at the true point of base parameter x, Q^b_4(x) the same for the base (e = 0) in the same frame:
    int (h^(5))^2 dxi = int (Q^t_4)^2 M_t dx,   Q^t_4 = Q^b_4 + [Q^t_4 - Q^b_4]   (at the same x: no second bracket is needed in L^2),
    |Q^t_4 - Q^b_4| <= A5 |e^(5)| + R,   A5 = |l_5| + r_5,   R = y0^-4 sum_{i <= 4} (|l_i| + r_i) B_i,
  by Taylor mode to order 4 in u = (x - x0)/y0 with the jets eta_i = y0^(i-1) e^(i), i = 0..5 (t9chart.E: exact linear coefficients l,
  nonlinear remainder r charged with the boxes; eta_4, eta_5 enter linearly and need no box).  So, by Minkowski,
    sqrt(int (h0^(5))^2 dX) <= L^3.5 sqrt(M_hi) [A5 sqrt(D_2(tau*) / (28 y_0^2)) + sqrt(|S|) R + sqrt(|S|) sup|h_base^(5)|].
  A5 sqrt(M_hi) is the chart factor (it replaces the typed 1.001); R carries the lower-order chain-rule terms.  Inputs: phi in [-P, P],
  P = sup|k| (X L + 2) 1.01 + beta, beta <= b0/(1 - sup|k| E0) the angle between the true and the base centre tangents; k^(m) balls
  [-sup|k^(m)|, sup|k^(m)|] (m <= 4, T6 hulls of the support cells); boxes |eta_0| <= E0/y0 (E0 = c_e Y + sqrt(F0)_Y (y^3/84)^(1/2) at
  the cells' x_hi, Y = 1.7e-12), |eta_1| <= b0, |eta_2| <= min(b1, y0 e2), |eta_3| <= min(b2, y0^2 e3), |eta_4| <= y0^3 e4 (n3 boxes
  b0, b1, b2, c_e and sqrt(F0)_Y as t9 recorded them in T9JSON 'b5'; e2, e3, e4 the C4 pointwise bounds at y0 with the amp's F_j and D2,
  decreasing in y; y0 <= x on S).  Base part as before (sympy majorant) with the same P.
  The arc over |xi| <= XL has |x - x_c| <= XL/M_lo; the excess over XL, XL (1/M_lo - 1), plus the hand-off frame's scale and shift
  ((1 - lam_f) x_hi + |x_c|, T9JSON 'frame') must fit in the supports' slack 1 (asserted).  lam_f <= 1 (asserted): the base-unit L^7 is
  conservative.
Cylinders as in p2/t9.py and the review: cone r_c in [200, 380] (rows of width 10, L = r_hi/4, support [r_lo/4 - 1, 7 r_hi/4 + 1],
M = 0.12); far centres 376 and 602 (L = 113, half-width 2L, M = 1; env FAR_M); far cells by radius or by arclength (x hull).
usage: python3 -I h4loc.py AMP_JSON [AMP_JSON ...]     (each a camp output; its handoff window, F_j and D2(tau*) are used;
       env T9JSON: t9_cone_kch.json, default p2/out/t9_cone_kch.json)"""
import os, sys, json, math
HERE = '/home/claude/s1spec/step4/p2/t8se'                      # v33: the original's directory (imports, relative defaults)
import hashlib
sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()
VAL = os.environ.get('V33_VALIDATE') == '1'
sys.path.insert(0, os.path.join(HERE, '..', '..', 't6'))
import sympy as sp
from flint import arb, ctx
import reader as RD
ctx.prec = 128
sys.path.insert(0, os.path.join(HERE, '..')); import t9chart as TC                       # E arithmetic (Taylor mode with jet boxes), as t9's KCH
T9J = os.environ.get('T9JSON') or os.path.normpath(os.path.join(HERE, '..', 'out', 't9_cone_kch.json'))   # env T9JSON: frame, n3 boxes
JD = os.path.join(os.path.normpath(os.path.abspath(os.environ['T6DIR'])) if os.environ.get('T6DIR')
                  else os.path.normpath(os.path.join(HERE, '..', 't6')), 'jets')          # env T6DIR as cconst
X0 = arb('2.2490776515')
U = lambda v: math.nextafter(float(v.upper()), math.inf)
def iv(lo, hi): return arb(lo).union(arb(hi))
def absb(a): return arb(U(abs(a)))
Lw = lambda v: math.nextafter(float(v.lower()), -math.inf)


def err_c4(F, D2, y):
    """C4 pointwise bounds |e''|, |e'''|, |e''''| at y (decreasing in y), from F_j and D2 at the good time (as p2/t9.py, p2/t9_far.py)"""
    r = arb('1.5') / arb('0.8828'); c2 = r * r; c3 = c2 * c2; y = arb(y)
    e2 = ((arb(3) / 28 * F[1] + 2 * r / 28 * (F[1] * F[2]).sqrt()) / y).sqrt()
    e3 = ((3 * c2 / 28 * F[2] + 2 * (c2 * c3).sqrt() / 28 * (F[2] * F[3]).sqrt()) / y ** 3).sqrt()
    e4 = (3 * c3 * F[3] / (28 * y ** 5) + 2 * (c2 * c3).sqrt() / 28 * (F[3] * D2).sqrt() / y ** 3).sqrt()
    return e2, e3, e4


def chart_h5(y0, P, ks, boxes):
    """Q^t_4 = y0^4 h^(5) of the true curve as a graph over a line at angle theta_c, at one point x0 of the support, as an E
    (t9chart): c (= Q^b_4) + sum_i l_i eta_i + NL, |NL| <= sum_i r_i |eta_i|, eta_i = y0^(i-1) e^(i) (i <= 5).  phi(x0) in [-P, P],
    k^(m)(x0) in [-ks[m], ks[m]] (m <= 4); boxes for eta_0..eta_3.  Slope F = N_t/M_t, N_t = a sin phi + e' cos phi, M_t = a cos phi
    - e' sin phi, a = 1 - k e (as t9chart.true_Q); D' = (1/M_t) d/du = y0 d/dxi, Q_j = D'^j F.  Returns (Q_4 at u = 0, M_t at u = 0)."""
    TC.BOX = list(boxes) + [None] * (6 - len(boxes))
    n = 4
    ph, kh, S, C, Yb = TC.base_series(y0, arb(0, P), [arb(0, k) for k in ks], arb(1), n=5)
    eh = [TC.E(0, {i: TC.ONE / TC.FACT[i]}) for i in range(n + 1)]              # e/y0 = sum eta_i u^i / i!
    ed = [TC.E(0, {i + 1: TC.ONE / TC.FACT[i]}) for i in range(n + 1)]          # e'   = sum eta_{i+1} u^i / i!
    Se, Ce, Ke = [TC.E(v) for v in S[:n + 1]], [TC.E(v) for v in C[:n + 1]], [TC.E(v) for v in kh[:n + 1]]
    a = TC.s_add([TC.E(1)] + [TC.E(0)] * n, TC.s_mul(Ke, eh, n), n, -1)
    Mt = TC.s_add(TC.s_mul(a, Ce, n), TC.s_mul(ed, Se, n), n, -1)
    Nt = TC.s_add(TC.s_mul(a, Se, n), TC.s_mul(ed, Ce, n), n)
    iM = TC.s_recip(Mt, n)
    Q = TC.s_mul(Nt, iM, n)
    for j in range(4): Q = TC.s_mul(iM, TC.s_der(Q), n - 1 - j)
    assert set(Q[0].l) | set(Q[0].r) <= set(range(6))
    return Q[0], Mt[0]


def h5_form():
    """h^(5) of a graph h(x) over a line, from k^(<=3) (arclength derivatives) and the angle psi (h' = tan psi); sympy gives
    h^(5) = [-60 k^4 S + 105 k^4 S/C^2 - 86 k^2 k1 C + 105 k^2 k1/C + 15 k k2 S + 10 k1^2 S + k3 C]/C^7  (S = sin psi, C = cos psi),
    checked below. Returned as a termwise majorant in K_m = sup|k^(m)| >= 0 and Pm = sup|psi| < pi/2."""
    s = sp.symbols('s'); k = sp.Function('k')(s); psi = sp.Function('psi')(s)
    D = lambda f: sp.diff(f, s).subs(sp.Derivative(psi, s), k) / sp.cos(psi)
    h = sp.tan(psi)
    for _ in range(4): h = D(h)
    K = sp.symbols('K0:4'); P = sp.symbols('P')
    rep = {sp.Derivative(k, (s, 3)): K[3], sp.Derivative(k, (s, 2)): K[2], sp.Derivative(k, s): K[1]}
    h = h.subs(rep).subs(k, K[0]).subs(psi, P)
    Sn, Cs = sp.sin(P), sp.cos(P)
    ref = (-60 * K[0] ** 4 * Sn + 105 * K[0] ** 4 * Sn / Cs ** 2 - 86 * K[0] ** 2 * K[1] * Cs + 105 * K[0] ** 2 * K[1] / Cs
           + 15 * K[0] * K[2] * Sn + 10 * K[1] ** 2 * Sn + K[3] * Cs) / Cs ** 7
    assert sp.simplify(h - ref) == 0

    def bound(K0, K1, K2, K3, Pm):
        sn, cs = Pm.sin(), Pm.cos()
        return (60 * K0 ** 4 * sn + 105 * K0 ** 4 * sn / cs ** 2 + 86 * K0 ** 2 * K1 + 105 * K0 ** 2 * K1 / cs + 15 * K0 * K2 * sn
                + 10 * K1 ** 2 * sn + K3) / cs ** 7
    return bound


def cells_of(slabs):
    out = []
    for s in slabs:
        J = RD.jets(s, d=JD)
        for e in range(J['x'].shape[0]):
            for c in range(J['x'].shape[1]):
                g = lambda nm: iv(*J[nm][e, c])
                out.append(dict(k=[g('k0'), g('k1'), g('k2'), g('k3'), g('k4')], p=g('p'), q=g('q'), x=(float(J['x'][e, c, 0]), float(J['x'][e, c, 1]))))
    for c in out: c['r'] = (c['p'] ** 2 + c['q'] ** 2).sqrt()
    return out


def cylinders():
    if VAL: return [dict(name=f'cone r_c {r}-{r + 10}', rlo=r, rhi=r + 10, L=arb(r + 10) / 4, X=3, M=arb('0.12')) for r in range(200, 380, 10)]
    return [dict(name=f'cone r_c {r}-{r + 10}', rlo=r, rhi=r + 10, L=arb(50), X=3, M=arb('0.12')) for r in range(300, 380, 10)]   # v33: fixed L = 50


def main():
    amps = [json.load(open(p)) for p in sys.argv[1:]]
    h5 = h5_form()
    T9 = json.load(open(T9J)); fr = T9['frame']; bx = T9['b5']['boxes']
    LAM = iv(*fr['lam']); XCF = arb(max(abs(v) for v in fr['xc']))
    assert LAM.upper() <= 1, 'frame scale lam_f > 1'
    YBOX = arb(bx['Y']); B0, B1, B2 = arb(repr(bx['b0'])), arb(repr(bx['b1'])), arb(repr(bx['b2'])); CE = arb(repr(bx['c_e'])); SF0 = arb(repr(bx['sqrtF0_at_Y']))
    rows = {}
    for path, A in zip(sys.argv[1:], amps):
        h = A['handoff']; D2 = arb(repr(h['at_Ybar']['D2_goodtime'])); win = h['window_slabs']; parts = h.get('D2_parts', {})
        F = [arb(repr(v)) ** 2 for v in h['at_Ybar']['sqrtF']]
        cells = cells_of(range(win[0], win[1] + 1))
        res = []
        for cy in cylinders():
            L, X = cy['L'], cy['X']
            if cy['name'].startswith('cone'):                             # union over the row of [r_c/4, 7 r_c/4] (as p2/t9.py)
                xlo, xhi = (0.25 * cy['rlo'] - 1, 1.75 * cy['rhi'] + 1) if VAL else (cy['rlo'] - 150 - 1, cy['rhi'] + 150 + 1)   # v33: [r_c - 3L, r_c + 3L] +- 1
            else:
                xlo = float((arb(cy['rlo']) - X * L).lower()) - 1; xhi = float((arb(cy['rhi']) + X * L).upper()) + 1
            y0 = arb(xlo) - X0 - 1
            assert y0 >= 28                                               # C4 pointwise bounds (far core)
            far = not cy['name'].startswith('cone')
            sel = [c for c in cells if (float(c['r'].upper()) >= xlo - 1 and float(c['r'].lower()) <= xhi + 1)
                   or (far and c['x'][1] >= xlo - 1 and c['x'][0] <= xhi + 1)]
            ks = [max((absb(c['k'][m]) for c in sel), key=U) for m in range(5)]
            E0 = max((CE * YBOX + SF0 * ((arb(c['x'][1]) - X0) ** 3 / 84).sqrt() for c in sel), key=U)   # |e| on the support (Y = 1.7e-12)
            e2, e3, e4 = err_c4(F, D2, y0)
            B = [E0 / y0, B0, min(B1, y0 * e2, key=U), min(B2, y0 ** 2 * e3, key=U), y0 ** 3 * e4]
            beta = B0 / (1 - ks[0] * E0)                                  # true vs base centre tangent: |beta| <= |tan beta| = |e'|/(1 - k e)
            P = ks[0] * (X * L + 2) * arb('1.01') + beta                  # |phi| on the arc: sup|k| (half-width + 2) * 1.01 + beta
            q, Mq = chart_h5(y0, arb(U(P)), [U(v) for v in ks], B[:4])
            A5 = arb(U(abs(q.l.get(5, arb(0))) + q.r.get(5, arb(0))))           # upper bounds as points (outward)
            R = arb(U(sum(((abs(q.l.get(i, arb(0))) + q.r.get(i, arb(0))) * B[i] for i in range(5)), arb(0)) / y0 ** 4))
            Mt = TC.enclose(Mq); mlo, mhi = arb(Lw(Mt)), arb(U(Mt)); assert mlo > 0
            ext = X * L * (1 / mlo - 1) + (1 - LAM.lower()) * arb(xhi) + XCF      # excess arclength + frame scale and shift
            assert ext <= 1 and X * L / mlo <= (X * L + 2) * arb('1.01'), (cy['name'], U(ext))
            b5 = h5(ks[0], ks[1], ks[2], ks[3], arb(U(P)))
            S_ = arb(xhi) - arb(xlo)
            base = L ** 7 * mhi * S_ * b5 * b5
            err = L ** 7 * mhi * A5 * A5 * D2 / (28 * y0 * y0)
            errR = L ** 7 * mhi * S_ * R * R
            tot = (err.sqrt() + errR.sqrt() + base.sqrt()) ** 2
            res.append(dict(name=cy['name'], L=U(L), y0=float(y0.lower()), error=U(err), error_lower=U(errR), base=U(base), total=U(tot),
                            M2=U(cy['M'] ** 2), ratio=U(tot / cy['M'] ** 2),
                            chart=dict(factor=U(A5 * mhi.sqrt()), A5=U(A5), l5=U(abs(q.l.get(5, arb(0)))), r5=U(q.r.get(5, arb(0))), M_t=[float(mlo.mid()), float(mhi.mid())],
                                       P=U(P), beta=U(beta), R=U(R), boxes=[U(v) for v in B], ext=U(ext))))
        worst = max(res, key=lambda r: r['ratio'])
        rows[path] = dict(window=win, window_length=parts.get('window_length'), D2=U(D2), D2_parts=parts, cylinders=res, worst=worst)
        print(f"{os.path.basename(path)}: window slabs {win}, |W| = {parts.get('window_length', float('nan')):.3e}, D2(tau*) <= {U(D2):.3e} "
              f"(Ybar-free {parts.get('zero', float('nan')):.2e}, linear {parts.get('linear', float('nan')):.2e}, "
              f"quadratic {parts.get('quadratic', float('nan')):.2e})")
        for r in res:
            if r['name'].startswith('far') or r['name'].endswith(('200-210', '290-300', '370-380')):
                print(f"   {r['name']:>18s}: L {r['L']:6.1f}  y0 {r['y0']:6.1f}  chart {r['chart']['factor']:.5f}  error {r['error']:.3e}  "
                      f"lower-order {r['error_lower']:.1e}  base {r['base']:.1e}  total/M^2 = {r['ratio']:.4g}")
        print(f"   worst: {worst['name']} total/M^2 = {worst['ratio']:.4g};  chart factor <= {max(r['chart']['factor'] for r in res):.5f} (was the typed 1.001)")
    env = {k: os.environ[k] for k in ('FAR_M', 'SUBW', 'T9JSON', 'T9WIN', 'T6DIR', 'H4OUT', 'TOPLOSS', 'NLFAC', 'BASE_T1', 'YBAR', 'T6RES', 'SLAB0', 'YOUT') if k in os.environ}
    out = dict(args=dict(argv=sys.argv[1:], amp=sys.argv[1:], jets=JD, chart_factor='computed per cylinder (review O6b; t9chart Taylor mode)',
                         t9json=T9J, frame=dict(lam=fr['lam'], xc=fr['xc']), boxes=bx, env=env), windows=rows)
    fo = os.environ.get('H4OUT') or '/home/claude/s1spec/step4/p2/out_v33/h4loc_v33.json'   # v33 default output
    worst_all = max((r for w in rows.values() for r in w['cylinders']), key=lambda r: r['ratio'])
    out['check'] = dict(item='Lemma O H4, fixed-L cone cylinders', worst=worst_all['name'], value_over_M2=worst_all['ratio'], total=worst_all['total'], threshold=0.0144,
                        margin=Lw(arb(1) / arb(repr(worst_all['ratio']))), ok=bool(arb(repr(worst_all['ratio'])) < 1))
    out['args']['inputs'] = {p: sha(p) for p in [os.path.abspath(__file__), T9J, f'{HERE}/../t9chart.py', os.path.join(JD, '..', 'residuals.json')] + sys.argv[1:]
                             + [os.path.join(JD, f'jets_{s:03d}.npz') for A in amps for s in range(A['handoff']['window_slabs'][0], A['handoff']['window_slabs'][1] + 1)]}
    if VAL:
        ref = json.load(open('/home/claude/s1spec/step4/p2/out/h4loc_A1_w180.json'))['windows']; rv = list(ref.values())[0]['cylinders']; mv = list(rows.values())[0]['cylinders']
        print('VALIDATE identical cone rows:', all(a == b for a, b in zip(mv, rv[:18])), len(mv)); sys.exit(0)
    out['args']['out'] = os.path.abspath(fo)
    json.dump(out, open(fo + '.tmp', 'w'), indent=1); os.replace(fo + '.tmp', fo)
    print(f"H4 fixed-L cone: worst {out['check']['worst']} total/M^2 = {out['check']['value_over_M2']:.4g}, margin x{out['check']['margin']:.3f}")
    assert out['check']['ok'], 'H4 fails'


if __name__ == '__main__':
    main()
