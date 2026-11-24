"""Lemma O H1 and H8 for the far family (centres x_c >= 376, L = 113, X = 2) on slab 180, uniformly in tau*, by a pairwise
separation test that holds for every centre at once.
Positions: Gamma(s, t) = Gamma(s, 0) + t Gamma_t(s, xi), |t| <= 1 (t_c), from the T6 engine on NS s-balls per element (as t9_lemmaO.py).
The mirror arm is conj(Gamma), with doubled-curve parameter x' = -x (tip) and 2 ell~ - x (far end).
  Graph (H1): over any arc of base length <= 2 G, with G = 1.1 X L, the tangent turns by at most K0 (2 G + 2), where K0 is t9_far's true
     curvature bound. So the arc is a graph over its centre tangent with slope <= tan(that), and sup|h0| <= X^2/2 L h2 + |e|/L.
  Separation (H8 and the rest of H1): every pair of points at doubled-curve distance > G + (cell lengths) is at distance
     >= D = rho_B + 0.1 L + |e|_P + |e|_Q at every t.  Here rho_B = sqrt((X L)^2 + (H L)^2) and H = sup|h0| + C0 sqrt(t1), C0 sqrt(t1) <= 0.05.
     Then no point of the curve outside the graph arc is in, or within 0.1 L of, any far box.
  Pairs: d(t) >= |P0 - Q0| - |P_t - Q_t| (hulls, outward), with cells bucketed on a grid so that only nearby pairs are tested.
usage: python3 -I t9_lemmaO_far.py AMP_JSON [SLAB]   (env LO_NS, default 48)"""
import sys, os, json, math
sys.path.insert(0, '/home/claude/s1spec/step4/p2/t9h')
from flint import arb, acb, arb_poly, ctx
import basegeo
ctx.prec = 96
SLAB = int(sys.argv[2]) if len(sys.argv) > 2 else 180; NS = int(os.environ.get('LO_NS', '48'))
amp = json.load(open(sys.argv[1])); sF0 = arb(repr(amp['handoff']['at_Ybar']['sqrtF'][0])); YB = arb('1e-12')
far = json.load(open('/home/claude/s1spec/step4/p2/out/t9_far.json')); X0 = 2.2490776515
U = lambda v: math.nextafter(float(v.upper()), math.inf); Lo = lambda v: math.nextafter(float(v.lower()), -math.inf)
L, X = 113.0, 2.0; XL = X * L; G = 1.1 * XL
def at_dt(g, i, s):
    el = g.els[i]; T_ = g.T; sh = arb_poly([arb(T_.mid()), arb(1)]); dl = arb(0, T_.rad())
    tv = lambda tp: tp(sh)(dl); tvd = lambda tp: tp.derivative()(sh)(dl)
    Gr, Gi = el['G'].at(s); Pr, Pi = el['Pg'].at(s)
    return acb(tvd(Pr), tvd(Pi)) * acb(tv(Gr), tv(Gi)) + acb(tv(Pr), tv(Pi)) * acb(tvd(Gr), tvd(Gi))
g0 = basegeo.Geo(SLAB, 4.0e5, arb(0)); gB = basegeo.Geo(SLAB, 4.0e5, arb(0, 1))
cells = []                                                    # (cx, cy, r0, vx, vy, rv, xlo, xhi, e)
for i in range(len(g0.els)):
    for j in range(NS):
        sb = arb(-1 + (2 * j + 1) / NS, 1 / NS); P0 = g0.at(i, sb)['Gam']; Pt = at_dt(gB, i, sb); xb = gB.at(i, sb)['x']
        xl, xh = Lo(xb), U(xb); y = max(xh - X0, 28.0)
        e = U(arb('1.6352e4') * YB + sF0 * (arb(y) ** 3 / 84).sqrt())
        r0 = U(abs(P0 - acb(P0.real.mid(), P0.imag.mid()))) * (1 + 1e-9); rv = U(abs(Pt - acb(Pt.real.mid(), Pt.imag.mid()))) * (1 + 1e-9)
        cells.append((float(P0.real.mid()), float(P0.imag.mid()), r0, float(Pt.real.mid()), float(Pt.imag.mid()), rv, xl, xh, e))
ell = max(c[7] for c in cells)
K0 = U(arb(repr(far['h2'])))                                   # |h''| bound (>= true curvature) on the far region
turn = K0 * (2 * G + 2); assert turn < 0.5, turn
h0 = X * X / 2 * L * far['h2']; emax_rel = max(c[8] for c in cells if c[7] >= 150) / L
H = h0 + emax_rel + 0.05; rhoB = math.hypot(XL, H * L)
pts = [(c, False) for c in cells] + [(c, True) for c in cells]
def xpar(c, m):                                               # doubled-curve parameter interval
    if not m: return c[6], c[7]
    return (2 * ell - c[7], 2 * ell - c[6])
B = 1.2 * (rhoB + 0.1 * L + 2 * max(c[8] for c in cells) + 2 * max(c[2] + c[5] for c in cells))
grid = {}
for k, (c, m) in enumerate(pts):
    cy = -c[1] if m else c[1]; grid.setdefault((int(c[0] // B), int(cy // B)), []).append(k)
worst = (math.inf, None); npairs = 0
for k, (c, m) in enumerate(pts):
    if m or c[6] < 150 - G: continue                            # P: far-region points of the arm (by symmetry this covers the mirror)
    cy = c[1]; bx, by = int(c[0] // B), int(cy // B); a0, a1 = xpar(c, m)
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for k2 in grid.get((bx + dx, by + dy), ()):
                c2, m2 = pts[k2]; b0, b1 = xpar(c2, m2)
                sep = max(b0 - a1, a0 - b1, 0.0)
                if m2: sep = min(sep, max(-c2[7] - a1, a0 + c2[6], 0.0))        # mirror through the tip: x' = -x
                if sep <= G + (a1 - a0) + (b1 - b0): continue                    # same graph arc (or adjacent): handled by H1
                q = (c2[0], -c2[1] if m2 else c2[1]); qv = (c2[3], -c2[4] if m2 else c2[4])
                d0 = math.hypot(c[0] - q[0], c[1] - q[1]) - c[2] - c2[2]; v = math.hypot(c[3] - qv[0], c[4] - qv[1]) + c[5] + c2[5]
                d = (d0 - v) * (1 - 1e-9) - (rhoB + 0.1 * L + c[8] + c2[8]); npairs += 1
                if d < worst[0]: worst = (d, (round(c[6]), round(c2[6]), m2))
out = dict(args=dict(argv=sys.argv[1:], slab=SLAB, NS=NS, L=L, X=X, G=G, C0=0.05, err='T7 4b + C4 core', t='t_c in [-1, 1], correlated'),
           turn=turn, sup_h0=h0 + emax_rel, H=H, rho_B=rhoB, pairs=npairs, min_slack=worst[0], worst=worst[1], ok=worst[0] > 0)
json.dump(out, open('/home/claude/s1spec/step4/p2/out/t9_lemmaO_far.json', 'w'), indent=1); print(json.dumps(out))
