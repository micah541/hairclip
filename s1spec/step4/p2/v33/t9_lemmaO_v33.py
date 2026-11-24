"""v3.3 task 6: copy of p2/t9_lemmaO.py (closed v3.2; read only) for Lemma O H1 (graph) on the FIXED-L cone family of step 3 v3.3:
unit L = 50, X = 3 (half-width 150), centres r_c in [300, 380] in rows of width 10 (centre arcs [rlo + 2.0, rlo + 12.6], as t9_lemmaO.py),
C0 sqrt(t1) = 0.00227 (lemmaO_v32_consts.json, case 'S1 cone, X=3', deviation C0 at the needed ratio 1.6001e-7: 0.00226524, rounded up).
Changes: the cone rows (L = 50 for every row; t9_lemmaO.py: L = r_hi/4, centres [200, 380]); the far family is not rerun (centres >= 490
are a subset of v3.2's family, certified arc by arc by t9_lemmaO_far.py); the H1 graph test is made explicit: besides cos(theta - theta_c) > 0
on the whole arc [x_lo - 1.3 XL, x_hi + 1.3 XL] (xi monotone; t9_lemmaO.py: 1.1 XL, too short for the end test with the
width-10 centre rows and the slab-time motion at L = 50), it asserts xi < -XL at the cells containing the arc's start and xi > XL at
the cells containing its end (t9_lemmaO.py computed these cells but did not test them), so the arc covers [-XL, XL] once and is a graph
there; sup|h0| = sup|eta|/L over the arc cells with |xi| <= XL (+ |e|/L), H := sup|h0| + C0 sqrt(t1), so sup|h0| <= H - C0 sqrt(t1).
H8 (other arms) is computed as in t9_lemmaO.py and reported, informational only (not a requirement: step 3 applies Lemma O arc by arc).
Output p2/out_v33/t9_lemmaO_v33.json.  Env V33_VALIDATE=1: the v3.2 cone rows (C0 0.0023) for a fidelity comparison in scratch.
Original docstring follows.
Lemma O (lemmaO_v3 H1, H8) on the hand-off slab, per cylinder, from the T6 engine positions over the whole slab (uniform in tau*) plus the
position error |e| (normal displacement, T7 Lemma 4b tip box and the core integral of LEMMA_F C4).
For a cylinder (centre arc x_c, unit L, half-width X L): positions and angles from the T6 engine on s-balls of every element over the whole slab, the centre frame
(Gamma_c, theta_c) from the cells at x_c; xi + i eta = (Gamma - Gamma_c) e^{-i theta_c}.
  graph arc: the cells meeting [x_c - 1.1 X L, x_c + 1.1 X L] (on the doubled curve near ell~: mirror cells by x -> 2 ell~ - x);
  checks:   cos(theta - theta_c) > 0 on the graph arc (xi monotone) and |xi| > X L at its ends, so [-X, X] is covered once;
            H1: sup|h0| = sup |eta| / L over graph cells with |xi| <= X L (+ |e|/L), and H := sup|h0| + C0 sqrt(t1);
            H8: every other cell, and every mirror cell (the second arm of the clip), has distance >= 0.1 L from the box
                [-X L, X L] x [-H L, H L] (so from the graph): dist = |(|xi| - X L)^+, (|eta| - H L)^+|.
usage: python3 -I t9_lemmaO.py AMP_JSON [SLAB]   (env LO_NS: s-balls per T6 element, default 48)"""
import sys, os, json, math
sys.path.insert(0, '/home/claude/s1spec/step4/t6')
import numpy as np
from flint import arb, acb, ctx
import reader as RD
ctx.prec = 96
VAL = os.environ.get('V33_VALIDATE') == '1'
C0C = 0.0023 if VAL else 0.00227                                  # v33: C0 sqrt(t1) at L = 50 (lemmaO_v32_consts.json, cone X=3)
import hashlib
sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()
SLAB = int(sys.argv[2]) if len(sys.argv) > 2 else 180
amp = json.load(open(sys.argv[1])); sF0 = arb(repr(amp['handoff']['at_Ybar']['sqrtF'][0])); YB = arb('1e-12')
JD = os.path.join(os.environ.get('T6DIR', '/home/claude/s1spec/step4/p2/t6'), 'jets'); J = RD.jets(SLAB, d=JD); X0 = 2.2490776515   # env T6DIR
U = lambda v: math.nextafter(float(v.upper()), math.inf); Lo = lambda v: math.nextafter(float(v.lower()), -math.inf)
iv = lambda a: arb(float(a[0])).union(arb(float(a[1])))
# positions: T6 engine (basegeo) over the whole slab (t_c in [-1, 1]), each element split into NS s-balls (rigorous cover)
sys.path.insert(0, '/home/claude/s1spec/step4/p2/t9h'); import basegeo
NS = int(os.environ.get('LO_NS', '48'))
# correlated in time: Gamma(s, t) = Gamma(s, 0) + t Gamma_t(s, xi), |t| <= 1 (t_c), so relative positions in a cylinder frame move by
# t (Gamma_t(s) - Gamma_t(s_c)) only (the far field expands with the rescaling; slab-wide hulls of positions would not see the cylinder move along)
from flint import arb_poly
def at_dt(g, i, s):                                              # t_c-derivative of the position at fixed s (as t9h/h1dt.py)
    el = g.els[i]; T_ = g.T; sh = arb_poly([arb(T_.mid()), arb(1)]); dl = arb(0, T_.rad())
    tv = lambda tp: tp(sh)(dl); tvd = lambda tp: tp.derivative()(sh)(dl)
    Gr, Gi = el['G'].at(s); Pr, Pi = el['Pg'].at(s)
    return dict(Gam=acb(tvd(Pr), tvd(Pi)) * acb(tv(Gr), tv(Gi)) + acb(tv(Pr), tv(Pi)) * acb(tvd(Gr), tvd(Gi)))
g0 = basegeo.Geo(SLAB, 4.0e5, arb(0)); gB = basegeo.Geo(SLAB, 4.0e5, arb(0, 1)); G = []; GT = []; TH = []; XS = []
for i in range(len(g0.els)):
    for j in range(NS):
        sb = arb(-1 + (2 * j + 1) / NS, 1 / NS); q0 = g0.at(i, sb); qB = gB.at(i, sb); qt = at_dt(gB, i, sb)
        G.append(q0['Gam']); GT.append(qt['Gam']); TH.append(qB['th']); xb = qB['x']; XS.append((Lo(xb), U(xb)))
xs = np.array(XS); n = len(xs); ell = float(xs[-1, 1])
def ebox(x):                                                     # |e| at arclength x (T7 Lemma 4b on x <= 34.25; LEMMA_F C4 core)
    y = max(x - X0, 28.0)
    return U(arb('1.6352e4') * YB + sF0 * (arb(y) ** 3 / 84).sqrt())
EB = [ebox(float(xs[i, 1])) for i in range(n)]
def check(name, xlo, xhi, L, X, C0):
    XL = X * L; L_ = arb(L)
    cen = [i for i in range(n) if xs[i, 1] >= xlo and xs[i, 0] <= xhi]
    Gc = G[cen[0]]; tc = TH[cen[0]]
    for i in cen[1:]: Gc = acb(Gc.real.union(G[i].real), Gc.imag.union(G[i].imag)); tc = tc.union(TH[i])
    rot = acb(0, -tc).exp(); GTc = GT[cen[0]]
    for i in cen[1:]: GTc = acb(GTc.real.union(GT[i].real), GTc.imag.union(GT[i].imag))
    tb = arb(0, 1)
    def frame(i, mirror=False):
        g = G[i].conjugate() if mirror else G[i]; gt = GT[i].conjugate() if mirror else GT[i]
        z = ((g - Gc) + tb * (gt - GTc)) * rot; r = EB[i] + 1e-12
        return arb(z.real.mid(), z.real.rad() + r), arb(z.imag.mid(), z.imag.rad() + r)
    ARC = 1.1 if VAL else 1.3                                     # v33: arc [x_lo - ARC XL, x_hi + ARC XL] (1.1 too short for the end test at L = 50)
    a0, a1 = xlo - ARC * XL, xhi + ARC * XL
    graph = [(i, False) for i in range(n) if xs[i, 1] >= a0 and xs[i, 0] <= a1]
    if a1 > ell: graph += [(i, True) for i in range(n) if 2 * ell - xs[i, 1] <= a1]       # doubled curve through the far end
    gset = set(graph)
    for i, m in graph:                                           # xi monotone on the arc
        assert Lo((TH[i] - tc).cos()) > 0 if not m else True, (name, 'tangent')
    ends = [i for i in range(n) if xs[i, 0] <= a0 <= xs[i, 1]] + [i for i in range(n) if xs[i, 0] <= a1 <= xs[i, 1]]
    e0 = [i for i in range(n) if xs[i, 0] <= a0 <= xs[i, 1]]; e1 = [i for i in range(n) if xs[i, 0] <= a1 <= xs[i, 1]]
    assert e0 and e1 and a1 <= ell, (name, 'arc ends')                # v33: explicit end test (cone arcs stay below ell)
    xi_start = max(U(frame(i)[0]) for i in e0); xi_end = min(Lo(frame(i)[0]) for i in e1)
    ends_ok = xi_start < -XL and xi_end > XL
    hmax = 0.0
    for i, m in graph:
        xi, eta = frame(i, m)
        if Lo(abs(xi)) <= XL: hmax = max(hmax, U(abs(eta)))
    H = hmax / L + C0; HL = H * L
    dmin = math.inf; wh = None
    for i in range(n):
        for m in (False, True):
            if (i, m) in gset: continue
            xi, eta = frame(i, m)
            dx = max(Lo(abs(xi)) - XL, 0.0); dy = max(Lo(abs(eta)) - HL, 0.0); d = math.hypot(dx, dy)
            if d < dmin: dmin, wh = d, (int(i), m)
    return dict(name=name, L=L, X=X, sup_h0=hmax / L, H=H, C0=C0, H1_graph=bool(ends_ok), xi_start_upper=xi_start, xi_end_lower=xi_end, XL=XL,
                H8_dist=dmin / L, worst=wh, ok=(dmin / L >= 0.1), H8_informational=True)
out = []
for rlo in (range(200, 380, 10) if VAL else range(300, 380, 10)):   # cone: r_c in [rlo, rlo + 10], X = 3; x ~ r + 2.26 (as t9_lemmaO.py)
    out.append(check(f'cone r_c {rlo}-{rlo + 10}', rlo + 2.0, rlo + 10 + 2.6, (rlo + 10) / 4 if VAL else 50.0, 3, C0C))
if VAL:
    ref = json.load(open('/home/claude/s1spec/step4/p2/out/t9_lemmaO_v32.json'))['rows']
    print([(a['sup_h0'], a['H'], a['H8_dist']) for a in out[:3]]); print('VALIDATE identical sup_h0, H, H8:', all(a['sup_h0'] == b['sup_h0'] and a['H'] == b['H'] and a['H8_dist'] == b['H8_dist'] for a, b in zip(out, ref)), len(out)); sys.exit(0)
summ = dict(item='Lemma O H1 (graph, arc by arc), fixed-L cone cylinders', label='certified',
            H1_all=all(r['H1_graph'] for r in out), cylinders=len(out), sup_h0_max=max(r['sup_h0'] for r in out), H_max=max(r['H'] for r in out),
            C0_sqrt_t1=C0C, H8_min_dist_informational=min(r['H8_dist'] for r in out), fails=[r['name'] for r in out if not r['H1_graph']], rows=out,
            args=dict(jets=JD, argv=sys.argv[1:], slab=SLAB, C0=dict(cone=C0C), H8='informational', err='T7 4b + C4 core',
                      env={k: os.environ[k] for k in ('T6DIR', 'LO_NS') if k in os.environ},
                      inputs={p: sha(p) for p in [os.path.abspath(__file__), sys.argv[1], os.path.join(JD, f'jets_{SLAB:03d}.npz'),
                                                  '/home/claude/s1spec/step4/p2/t9h/basegeo.py', f'/home/claude/s1spec/step4/t1/polished_v2/slab_{SLAB:03d}.pkl']}))
fo = '/home/claude/s1spec/step4/p2/out_v33/t9_lemmaO_v33.json'
json.dump(summ, open(fo + '.tmp', 'w'), indent=1); os.replace(fo + '.tmp', fo)
print(json.dumps({k: v for k, v in summ.items() if k != 'rows'}))
assert summ['H1_all'], 'H1 graph fails'
