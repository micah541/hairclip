"""Lemma O H1 (graph) for the far family, arc by arc, step 3 v3.2: X = 2, unit L = 113 (lambda_1 units), (M, K) = (0.15, 0.3), every centre
on the arm (and, by the reflection symmetry, on the mirror arm) from the first far centre to the far end; certified pass/fail, arb, outward.
H1: in the frame of the centre with the curve's OWN tangent there horizontal, the arc of gamma through the centre meets [-X, X] x [-H, H]
in the graph of h0 over [-X, X], with sup|h0| <= H - C0 sqrt(t1).  Other arcs through the box are irrelevant: the curve is immersed and
step 3 applies Lemma O arc by arc (v3.1 l.97, l.105, l.111; amended (H5)); H8 (separation from other arms) is not required and not tested.
Data (FAR_JSON = t9_far_v32.json, slab-180 hulls at the good time; t9_far.py): on base arclength x >= XK = 148 (t9_far: k hulls of the
cells with x_hi >= 148, C4 error bounds at y = 147 - X0 and below them decreasing in y, |e| and |e'| at the far end and below increasing)
  |kappa_true| <= K0 <= h2 (t9_far's h'' = K0 sec^3 psi >= K0),  |k~| <= k0,  |e| <= e0,  |e'| <= e1.
True curve gamma = Gamma~ + e N over the base (x = base arclength): |gamma'| = J = sqrt((1 - k~ e)^2 + e'^2) in [Jmin, Jmax],
  Jmin = 1 - k0 e0,  Jmax = sqrt((1 + k0 e0)^2 + e1^2).
For a centre x_c and |x - x_c| <= G (doubled curve: through the far end the arc continues on the mirror arm, with the same bounds):
  turning: |theta_gamma(x) - theta_gamma(x_c)| <= K0 Jmax G =: phi < pi/2, so xi' = J cos(theta - theta_c) >= Jmin cos(phi) =: m > 0;
  xi is strictly increasing and |xi(x_c +- G)| >= G m > X L, so the arc meets {|xi| <= X L} in one sub-arc, a graph over [-X L, X L];
  h(0) = h'(0) = 0 (own frame) and |h''| = |kappa| sec^3 <= K0 sec^3(phi), so sup|h0| = sup|h| / L <= (X^2 / 2) L K0 sec^3(phi).
  The arc stays where the bounds hold: x_c - G >= XK for every far centre (asserted at the first centre, below).
H := sup|h0| + C0 sqrt(t1), rounded up to 4 decimals; C0 sqrt(t1) is the v3.2 far value (CONSTS_JSON far_v32.case,
  deviation_upper_at_needed_time.C0, at the needed ratio; (M, K) = (0.15, 0.3), X = 2, L = 113 asserted).  Then the arc stays in
  |eta| <= (H - C0 sqrt(t1)) L < H L on |xi| <= X L and leaves the box through its vertical sides: H1 holds with this H.
Frame (T9JSON = t9_cone_kch.json, block 'frame'): hand-off lengths are base lengths / lam_f, lam_f <= 1 (asserted), so the base-unit
  computation is conservative (kappa_h = lam_f kappa, X L in base units = lam_f X L <= X L); the frame shift x_c(frame) is a translation.
First far centre, both readings of 'centres >= 376' (hand-off units): arclength x_c >= 376 is base arclength >= 376 lam_f (asserted
  >= XK + G); radius |gamma - x_c(frame)| >= 376 lam_f: every slab-180 cell (T6 jets, p, q hulls) with x_lo < XK + G has
  r = |Gamma~| < 376 lam_f - |x_c(frame)| - |e| there (|e| <= c_e Ybar + sqrt(F0) (y^3/84)^(1/2) at its x_hi, as t9_lemmaO.py).
usage: python3 -I t9_lemmaO_far.py AMP_JSON FAR_JSON CONSTS_JSON   (env T6DIR: jets of slab 180; env T9JSON: t9_cone_kch.json)
output: p2/out/t9_lemmaO_far.json"""
import sys, os, json, math
from fractions import Fraction
from flint import arb, ctx
ctx.prec = 128
sys.path.insert(0, '/home/claude/s1spec/step4/t6')
import reader as RD
P2 = '/home/claude/s1spec/step4/p2'
U = lambda v: math.nextafter(float(v.upper()), math.inf); Lo = lambda v: math.nextafter(float(v.lower()), -math.inf)
amp, far, lo = (json.load(open(p)) for p in sys.argv[1:4])
T9J = os.environ.get('T9JSON', f'{P2}/out/t9_cone_kch.json'); fr = json.load(open(T9J))['frame']
JD = os.path.join(os.environ.get('T6DIR', f'{P2}/t6c'), 'jets'); SLAB = 180
# ---- step-3 v3.2 far family (constants file) and step-4 far data (t9_far --v32)
fc = lo['far_v32']['case']; L = arb(lo['far_v32']['L']); X = arb(Fraction(fc['X']).numerator) / Fraction(fc['X']).denominator
assert lo['far_v32']['L'] == 113 and Fraction(fc['X']) == 2 and Fraction(fc['M']) == Fraction(3, 20) and Fraction(fc['K2']) == Fraction(fc['K3']) == Fraction(3, 10)
C0S = arb(repr(fc['deviation_upper_at_needed_time']['C0']))                        # C0 sqrt(t1) (upper bound, outward at 1e-6)
assert far['args']['v32'] and far['args']['slab_hulls'] == SLAB and far['args']['y_err'] == 147
XK = 148.0                                                                          # t9_far.py: hulls of the cells with x_hi >= 148, C4 at y = 147 - X0
K0 = arb(repr(far['h2'])); k0 = arb(repr(far['base']['k0'])); e0 = arb(repr(far['err']['e0'])); e1 = arb(repr(far['err']['e1']))
LAM = arb(fr['lam'][0]).union(arb(fr['lam'][1])); XCF = max(abs(fr['xc'][0]), abs(fr['xc'][1]))
assert LAM.upper() <= 1, 'frame scale lam_f > 1: the base-unit computation would not be conservative'
XL = X * L
# ---- arc length G (base arclength on each side of the centre): smallest multiple of 1/8 with G Jmin cos(K0 Jmax G) > X L
Jmin = 1 - k0 * e0; Jmax = ((1 + k0 * e0) ** 2 + e1 * e1).sqrt(); assert Jmin > 0
g = float(XL.mid())
for _ in range(50): g = float(XL.mid()) / (float(Jmin.mid()) * math.cos(float((K0 * Jmax).mid()) * g))
G = arb(math.ceil(g * 8 + 1) / 8)
phi = K0 * Jmax * G; m = Jmin * phi.cos(); sec = 1 / phi.cos()
ok_turn = bool(phi < 1)                                                             # < pi/2 with room: xi strictly increasing on the arc
ok_cover = bool(G * m > XL)                                                         # the arc leaves |xi| <= X L at both ends
sup_h0 = X * X / 2 * L * K0 * sec ** 3                                              # cylinder units
H = arb(math.ceil(U(sup_h0 + C0S) * 1e4) / 1e4)
ok_height = bool(sup_h0 <= H - C0S)
# ---- first far centre: the arc must stay in x >= XK
xfirst = arb(XK) + G
ok_arcl = bool(xfirst <= 376 * LAM)                                                 # arclength reading: base x_c >= 376 lam_f >= XK + G
# radius reading.  Anchor: a slab-180 cell (x hull [a, b], r = |Gamma~| <= r_hi; slab-wide hulls, so at every time) near the tip.  At any
# time, with z a point of that cell, |Gamma~(x)| <= |Gamma~(z)| + |x - x(z)| <= r_hi + max(x - a, b) (arclength >= chord), and
# |gamma| <= |Gamma~| + |e|.  So every point with base x < XK + G has |gamma| <= r_hi + max(XK + G - a, b) + |e|(XK + G) =: R1, and the
# radius reading is covered if R1 + |x_c(frame)| < 376 lam_f.  |e| <= c_e Ybar + sqrt(F0) (y^3/84)^(1/2), increasing in y (as t9_lemmaO.py).
J = RD.jets(SLAB, d=JD); YB = arb('1e-12'); sF0 = arb(repr(amp['handoff']['at_Ybar']['sqrtF'][0])); X0 = arb('2.2490776515')
anc = None
for e_ in range(J['x'].shape[0]):
    for c_ in range(J['x'].shape[1]):
        a_, b_ = (float(v) for v in J['x'][e_, c_])
        if b_ > 40: continue
        p = arb(float(J['p'][e_, c_, 0])).union(arb(float(J['p'][e_, c_, 1]))); q = arb(float(J['q'][e_, c_, 0])).union(arb(float(J['q'][e_, c_, 1])))
        rh = U((p * p + q * q).sqrt())
        if anc is None or a_ - rh > anc[0] - anc[3]: anc = (a_, b_, (e_, c_), rh)
eloc = arb('1.6352e4') * YB + sF0 * ((xfirst - X0) ** 3 / 84).sqrt()
R1 = arb(anc[3]) + max(xfirst - anc[0], arb(anc[1]), key=lambda v: float(v.upper())) + eloc
ok_rad = bool(R1 + XCF < 376 * LAM)                                                 # radius reading: centres with r >= 376 have x >= XK + G
checks = dict(graph_monotone=ok_turn, graph_covers_X=ok_cover, height=ok_height, first_centre_arclength=ok_arcl, first_centre_radius=ok_rad)
out = dict(args=dict(argv=sys.argv[1:], slab=SLAB, jets=JD, t9json=T9J, reading='arc by arc (H8 not required, not tested)',
                     family=dict(X=fc['X'], L=113, M=fc['M'], K2=fc['K2'], K3=fc['K3'], needed=fc['needed'], source='far_v32.case')),
           data=dict(K0_le_h2=U(K0), k0=U(k0), e0=U(e0), e1=U(e1), XK=XK, lam_frame=[Lo(LAM), U(LAM)], xc_frame_abs=XCF),
           arc=dict(G=float(G.mid()), Jmin=Lo(Jmin), Jmax=U(Jmax), turn_max=U(phi), xi_dot_min=Lo(m), xi_end_min=Lo(G * m), XL=float(XL.mid())),
           first_centre=dict(x_needed=U(xfirst), x_arclength_reading=Lo(376 * LAM), anchor=dict(cell=list(anc[2]), x=[anc[0], anc[1]], r_hi=anc[3]),
                             e_loc=U(eloc), r_max_before=U(R1 + XCF), r_radius_reading=Lo(376 * LAM)),
           C0_sqrt_t1=U(C0S), sup_h0=U(sup_h0), H=float(H.mid()), H_minus_C0_sqrt_t1=Lo(H - C0S), slack_height=Lo(H - C0S - sup_h0),
           checks=checks, ok=all(checks.values()))
json.dump(out, open(f'{P2}/out/t9_lemmaO_far.json', 'w'), indent=1)
print(f"far H1 (v3.2, X = 2, L = 113, arc by arc): arc G = {float(G.mid())} each side, turning <= {U(phi):.5f} rad, xi(+-G) >= {Lo(G * m):.3f} > X L = 226; "
      f"sup|h0| <= {U(sup_h0):.5f}, C0 sqrt(t1) <= {U(C0S):.3e}, H = {float(H.mid())}; first centre x >= {U(xfirst):.3f} "
      f"(arclength 376 lam_f = {Lo(376 * LAM):.3f}; radius: points before have r <= {U(R1 + XCF):.3f}); checks {checks}: {'PASS' if out['ok'] else 'FAIL'}")
assert out['ok'], checks
