"""Lemma O (lemmaO_v3 H1, H8) on the hand-off slab, per cylinder, from the T6 engine positions over the whole slab (uniform in tau*) plus the
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
V32 = '--v32' in sys.argv; sys.argv = [a for a in sys.argv if a != '--v32']   # --v32: cone C0 sqrt(t) = 0.0023 (ratio 1.6e-7), output t9_lemmaO_v32.json
C0C = 0.0023 if V32 else 0.0021
SLAB = int(sys.argv[2]) if len(sys.argv) > 2 else 180
amp = json.load(open(sys.argv[1])); sF0 = arb(repr(amp['handoff']['at_Ybar']['sqrtF'][0])); YB = arb('1e-12')
J = RD.jets(SLAB, d='/home/claude/s1spec/step4/p2/t6/jets'); X0 = 2.2490776515
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
    a0, a1 = xlo - 1.1 * XL, xhi + 1.1 * XL
    graph = [(i, False) for i in range(n) if xs[i, 1] >= a0 and xs[i, 0] <= a1]
    if a1 > ell: graph += [(i, True) for i in range(n) if 2 * ell - xs[i, 1] <= a1]       # doubled curve through the far end
    gset = set(graph)
    for i, m in graph:                                           # xi monotone on the arc
        assert Lo((TH[i] - tc).cos()) > 0 if not m else True, (name, 'tangent')
    ends = [i for i in range(n) if xs[i, 0] <= a0 <= xs[i, 1]] + [i for i in range(n) if xs[i, 0] <= a1 <= xs[i, 1]]
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
    return dict(name=name, L=L, X=X, sup_h0=hmax / L, H=H, H8_dist=dmin / L, worst=wh, ok=(dmin / L >= 0.1))
out = []
for rlo in range(200, 380, 10):                                    # cone: r_c in [rlo, rlo + 10], L = r/4 (as t9.py), X = 3; x ~ r + 2.26
    out.append(check(f'cone r_c {rlo}-{rlo + 10}', rlo + 2.0, rlo + 10 + 2.6, (rlo + 10) / 4, 3, C0C))
c = 376.0
while c < ell:
    w = max(56.5, 0.05 * c); out.append(check(f'far centre {c:.0f}', c, min(c + w, ell), 113.0, 2, 0.05)); c += w
bad = [r for r in out if not r['ok']]
summ = dict(args=dict(argv=sys.argv[1:], slab=SLAB, C0=dict(cone=C0C, far=0.05), v32=V32, H8_min=0.1, err='T7 4b + C4 core'),
            cylinders=len(out), H8_min_dist=min(r['H8_dist'] for r in out), worst=min(out, key=lambda r: r['H8_dist']),
            sup_h0_max=dict(cone=max(r['sup_h0'] for r in out if r['name'].startswith('cone')), far=max(r['sup_h0'] for r in out if r['name'].startswith('far'))),
            fails=[r['name'] for r in bad], rows=out)
json.dump(summ, open('/home/claude/s1spec/step4/p2/out/t9_lemmaO_v32.json' if V32 else '/home/claude/s1spec/step4/p2/out/t9_lemmaO.json', 'w'), indent=1)
print(json.dumps({k: v for k, v in summ.items() if k != 'rows'}))
