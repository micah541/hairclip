"""T9 (part): certified base cone jets (H3 K4, H5 cone cylinders) and the chart factor KCH, in arb, from the corrected T6 jet
hulls (p2/t6/jets; each hull is a sup over its sub-cell and over the whole slab), slabs 177-180 (the good-time window).
Cone cylinders (Lemma O, H5): centre radius r_c in [200, 380], L = r_c/4, X = 3 (half-width 0.75 r_c), graph over the
  centre tangent: h'' = k sec^3, h''' = sec^4 (k_s + 3 k^2 t), h'''' = sec^5 (k_ss + 10 k k_s t + 12 k^3 t^2 + 3 k^3 sec^2);
  |psi| <= sup|k| * (0.75 r_c + 2) * 1.01; K_j = L^(j-1) sup|h^(j)|, N3 = L^1.5 ||h'''||_L2(|x|<=XL) <= sqrt(6) K3.
  r = |Gamma~| = sqrt(p^2 + q^2) per cell; any frame with |x_c| <= 0.5, |lambda - 1| <= 1e-3 is covered by the slack 1.
  Error parts at Ybar (C4 at the cylinder's inner end y = 0.25 r_c - X0, constants as in boxes.py), with D2(tau*) from camp.
KCH and the e-terms (LEMMA_F C3; review B5, no typed allowance or factor): the exact decomposition (*) of p2/t9chart.py,
  D^j f_t(y) - D^j f_b(y) = [Q^t_j - Q^b_j](X_t(y)) + [Q^b_j(X_t(y)) - Q^b_j(X_b(y))],
  first bracket by Taylor mode in arb (exact linear part = D^j(sec^2 phi e' + sec^2 phi tan phi k e) minus the second bracket's,
  checked against sympy in p2/t9chart_check.py; the nonlinear remainder charged to the jets eta_i = y0^(i-1) e^(i) with the n3 boxes),
  second bracket <= sup|dQ^b_j/dx| |e| sup|sin phi|/inf cos phi.  E_j = int_[yw,1.5yw] (28/y^2)(D^j f)^2 dy, dy = M_t dx,
  T_m >= 28 y_T8^(2m-2) (m >= 1), T_0 >= 28/y_T8^2 off the flattening:
    sqrt(E_j[err]) <= sum_m K_jm sqrt(F_m) + sqrt(28/(3 yw)) sup_w(e-coefficient) sup_w|e|,  K_jm = sup (|l_{m+1}| + r_{m+1}) sqrt(M_t)/Yhat_t,
    |e(x)| <= c_e Ybar + sqrt(F0) ((x - X0)^3/84)^(1/2)   (T7 Lemma 4b on x <= 34.249, i.e. at y = 28; LEMMA_F C4 core on [28, y]).
  KCH = max_{(j,m) = (0,0), m >= 1} K_jm/S(j,m);  c0_j = K_j0 (j >= 1) + e-part per sqrt(F0);  c0tip_j = e-part per unit Ybar.
  Boxes (valid for Ybar <= n3's Y = 1.7e-12, n3boxes): |e'| <= b0, x|e''| <= b1, x^2|e'''| <= b2 (n3.json args), |e| as above with
  sqrt(F0) <= aF0 Y + bF0 at Y = 1.7e-12 (this amp, slabs SLABS).
  y_ray/y_T8 (T9RATIO=tip, default): y_ray(x) = y_ray(0) + x - int_0^x (1 - cos phi), the integral bounded on the arrangement of the
  cells' x-hulls (the true sub-cells tile the curve at every time), y_ray(0) = p(0) cos(Theta) by the base's symmetry about the real
  axis (q(0) = 0, theta(0) = -pi/2, asserted against the tip hulls).  The former q - x enclosure [q_lo - a, q_hi - b] (T9RATIO=qx)
  needs q - x monotone at fixed time AND static over the slab; the hulls are slab-wide, so it is not certified.
Inputs used as stated: Theta = THINF = 1.8617364190 +/- 1e-3 (CAPD gives int_tip^8 k_S1; the tail int_8^oo k_S1 = 1.13e-4 is the
  float value, covered by the allowance), X0 = 2.2490776515 (T8's y = x - X0).
usage: python3 -I t9.py AMP_JSON"""
import sys, os, json, math
import numpy as np
from flint import arb, ctx
ctx.prec = 128
sys.path.insert(0, '/home/claude/s1spec/step4/t6')
import reader as RD
sys.path.insert(0, '/home/claude/s1spec/step4/p2'); import curvjets as CJ
import t9chart as TC
N3P = os.environ.get('N3JSON', '/home/claude/s1spec/step4/p2/t8se/constants/n3.json')   # the A1 chain's n3far run (checked by n3boxes)
N3A = json.load(open(N3P))['args']
B0 = N3A['b0']                                                                          # |e'| box
T7C = '/home/claude/s1spec/step4/t7/constants.json'
JD = '/home/claude/s1spec/step4/p2/t6/jets'
SLABS = range(177, 181)
X0 = arb('2.2490776515'); TH = arb('1.8617364190179404', '1e-3')
amp = json.load(open(sys.argv[1])); hb = amp['handoff']['at_Ybar']
sF = [arb(repr(v)) for v in hb['sqrtF']]; F = [v * v for v in sF]; D2 = arb(repr(hb['D2_goodtime']))
U = lambda v: math.nextafter(float(v.upper()), math.inf)
def iv(lo, hi): return arb(lo).union(arb(hi))
def absb(a): return arb(U(abs(a)))

# ---- cells over the window slabs
cells = []
for s in SLABS:
    J = RD.jets(s, d=JD)
    for e in range(J['x'].shape[0]):
        for c in range(J['x'].shape[1]):
            g = lambda nm: iv(*J[nm][e, c])
            cells.append(dict(slab=s, a=float(J['x'][e, c, 0]), b=float(J['x'][e, c, 1]), k0=g('k0'), k1=g('k1'), k2=g('k2'), k3=g('k3'),
                              th=g('theta'), p=g('p'), q=g('q')))
for c in cells: c['r'] = (c['p'] ** 2 + c['q'] ** 2).sqrt()

# ---- error pointwise bounds at y (C4; T_0 >= 28/(1.5y)^2, T_m >= 28 (0.8828 y/1.5)^(2m-2))
rr = arb('1.5') / arb('0.8828'); c2 = rr * rr; c3 = c2 * c2
def err(y):
    y = arb(y)
    e2 = ((arb(3) / 28 * F[1] + 2 * rr / 28 * (F[1] * F[2]).sqrt()) / y).sqrt()
    e3 = ((3 * c2 / 28 * F[2] + 2 * (c2 * c3).sqrt() / 28 * (F[2] * F[3]).sqrt()) / y ** 3).sqrt()
    e4 = (3 * c3 * F[3] / (28 * y ** 5) + 2 * (c2 * c3).sqrt() / 28 * (F[3] * D2).sqrt() / y ** 3).sqrt()
    return e2, e3, e4

# ---- cone jets
out = dict(cone=[])
worst = dict(K2=0.0, K3=0.0, K4=0.0, N3=0.0, slope=0.0)
for rlo in range(200, 380, 10):
    rhi = rlo + 10; L = arb(rhi) / 4
    lo, hi = 0.25 * rlo - 1, 1.75 * rhi + 1
    sel = [c for c in cells if float(c['r'].upper()) >= lo and float(c['r'].lower()) <= hi]
    k0 = max((absb(c['k0']) for c in sel), key=U); k1 = max((absb(c['k1']) for c in sel), key=U); k2 = max((absb(c['k2']) for c in sel), key=U)
    e2, e3, e4 = err(max(0.25 * rlo - 1 - 2.25, 40.0))
    k3 = max((absb(c['k3']) for c in sel), key=U)
    yhi = arb(1.75 * rhi + 1) ; E0 = arb(repr(N3A['eta1'])) + F[0].sqrt() * (yhi ** 3 / 84).sqrt()   # |e| <= tip box eta1 (n3boxes vs T7 4b) + core integral (LEMMA_F C4)
    K0, K1, K2_ = (arb(v) for v in CJ.kmajor([arb(0, U(k0)), arb(0, U(k1)), arb(0, U(k2)), arb(0, U(k3))], [U(E0), B0, U(e2), U(e3), U(e4)], 3))
    psi = K0 * (arb('0.75') * rhi + 2) * arb('1.01'); t = psi.tan(); sec = 1 / psi.cos()
    assert U(sec) <= 1.01                                  # true half-length <= 1.01 (0.75 r + 2), so psi <= K0 * that
    h2 = K0 * sec ** 3; h3 = sec ** 4 * (K1 + 3 * K0 ** 2 * t)
    h4 = sec ** 5 * (K2_ + 10 * K0 * K1 * t + 12 * K0 ** 3 * t * t + 3 * K0 ** 3 * sec * sec)
    row = dict(r=[rlo, rhi], L=U(L), slope=U(t), K2=U(L * h2), K3=U(L ** 2 * h3), K4=U(L ** 3 * h4), N3=U(arb(6).sqrt() * L ** 2 * h3),
               base_k=[U(k0), U(k1), U(k2)], err=[U(e2), U(e3), U(e4)])
    out['cone'].append(row)
    for k in worst: worst[k] = max(worst[k], row[k])
out['cone_max'] = worst
h3k4 = out['cone'][0]                                       # r_c in [200, 210]: y in [50, 350] (H3's K4)
print(f"cone (r_c 200-380, base + error at Ybar): K2 <= {worst['K2']:.3e}, K3 <= {worst['K3']:.3e}, K4 <= {worst['K4']:.3e}, "
      f"N3 <= {worst['N3']:.3e}, slope <= {worst['slope']:.4f}  (H5: 0.12, 0.12, 0.12, 0.14);  H3 K4 (r_c = 200): {h3k4['K4']:.3e} (<= 0.06)")

# ---- KCH and the e-terms (review B5): exact decomposition (*) of p2/t9chart.py, no typed allowance or factor
S = {(1, 1): 1, (2, 1): 1, (2, 2): 1, (3, 1): 1, (3, 2): 3, (3, 3): 1, (0, 0): 1}
YBOX = arb(str(N3A['Y']))                                   # the n3 bootstrap Ybar: boxes hold for every Ybar <= YBOX (n3boxes)
TIPB = json.load(open(T7C))['meta']['tip_box']
assert TIPB['holds'] and abs(TIPB['Ybox'] - float(YBOX.mid())) <= 1e-6 * float(YBOX.mid()), 'T7 tip box at another Ybar'
CE = arb(repr(TIPB['c_e']))                                 # |e(x)| <= c_e Lambda_e^1/2 <= c_e Ybar on x <= 34.249 (T7 Lemma 4b)
_ps, _sl = amp['per_slab'], amp['slabs']
SF0BOX = max((arb(repr(_ps['0']['aF'][_sl.index(s)])) * YBOX + arb(repr(_ps['0']['bF'][_sl.index(s)])) for s in SLABS), key=U)
BX = [None, arb(repr(N3A['b0'])), arb(repr(N3A['b1'])), arb(repr(N3A['b2'])), None]   # |eta_1| <= b0, |eta_2| <= b1 y0/x, |eta_3| <= b2 (y0/x)^2
def ebox(yT8):                                              # |e| at x = X0 + yT8 >= 30.25, at Ybar <= YBOX
    return CE * YBOX + SF0BOX * (arb(yT8) ** 3 / 84).sqrt()
YW_LO, YW_HI = 50.0, 1.5 * 253                              # union of the windows [y, 1.5 y], y in [50, 253]
RATIO = os.environ.get('T9RATIO', 'tip')                    # y_ray/y_T8: 'tip' (rigorous, below; default) or 'qx' (q - x hull: static monotone argument, not implied by the time hulls)
def tip_anchor(slab):
    """y_ray(x) = y_ray(0) + x - int_0^x (1 - cos phi) at every time of the slab: the true sub-cells tile the curve, so for s in an elementary
    interval of the cell arrangement, 1 - cos(phi(s)) lies between the min lower / max upper bound over the cells covering it.  Anchor: the
    base is symmetric about the real axis with the tip on it (q(0) = 0, theta(0) = -pi/2, asserted against the tip hulls), so y_ray(0) = p(0) cos(Theta)."""
    cs = [c for c in cells if c['slab'] == slab]
    t0 = min(cs, key=lambda c: c['a']); assert t0['a'] <= 0 and t0['q'].contains(0) and t0['th'].contains(-arb.pi() / 2)
    yr0 = t0['p'] * TH.cos()
    xm = max(c['b'] for c in cs if c['a'] <= 2 * YW_HI); bp = sorted(set([c['a'] for c in cs] + [c['b'] for c in cs])); bp = [v for v in bp if v <= xm]
    v = [(c['a'], c['b'], 1 - (c['th'] - TH).cos()) for c in cs]
    Il, Ih = [arb(0)], [arb(0)]; lb, ub = [], []
    for k in range(len(bp) - 1):
        cov = [w for (a, b, w) in v if a <= bp[k] and b >= bp[k + 1]]; assert cov, (slab, bp[k])
        lb.append(min((w.lower() for w in cov), key=lambda z: float(z))); ub.append(max((w.upper() for w in cov), key=lambda z: float(z)))
        h = arb(bp[k + 1]) - arb(bp[k]); Il.append(Il[-1] + h * lb[-1]); Ih.append(Ih[-1] + h * ub[-1])
    return dict(yr0=yr0, bp=bp, Il=Il, Ih=Ih, lb=lb, ub=ub)
def I_at(T, x, hi):                                         # int_0^x (1 - cos phi) lower (hi=False) / upper (hi=True) bound, x >= 0
    import bisect
    assert 0 <= x <= T['bp'][-1]; k = min(max(0, bisect.bisect_right(T['bp'], x) - 1), len(T['lb']) - 1)
    return (T['Ih'][k] + (arb(x) - arb(T['bp'][k])) * T['ub'][k]) if hi else (T['Il'][k] + (arb(x) - arb(T['bp'][k])) * T['lb'][k])
TIPA = {s: tip_anchor(s) for s in SLABS} if RATIO == 'tip' else {}
PIECE = float(os.environ.get('T9PIECE', '0.5'))          # x is exact: evaluate on pieces of a cell's x-hull (the jet hulls stay those of the cell)
rows = []
for c in cells:
    if c['a'] > 2 * YW_HI or c['b'] < 40: continue          # windows: on the arm near the tip (Y_b' = cos phi > 0), and y_ray <= y_ray(0) + x <= 0.91 + x
    a_, b_ = arb(c['a']), arb(c['b'])
    phi = c['th'] - TH
    qmx = iv(float((c['q'].lower() - a_).lower()), float((c['q'].upper() - b_).upper()))   # q - x (monotone in x)
    npc = max(1, math.ceil((c['b'] - c['a']) / PIECE)); xs = [c['a'] + (c['b'] - c['a']) * i / npc for i in range(npc)] + [c['b']]
    for xlo, xhi in zip(xs[:-1], xs[1:]):
        yT8 = iv(xlo, xhi) - X0
        if RATIO == 'tip':                                  # y_ray = y_ray(0) + x - I(x): ratio = 1 + (y_ray(0) + X0 - I(x)) / y_T8
            T = TIPA[c['slab']]
            num = iv(float((T['yr0'] + X0 - I_at(T, xhi, True)).lower()), float((T['yr0'] + X0 - I_at(T, xlo, False)).upper()))
            ratio = 1 + num / yT8
        else:
            ratio = phi.cos() * (1 + (qmx + X0) / yT8) - c['p'] * phi.sin() / yT8        # y_ray / y_T8
        yr = yT8 * ratio
        assert yr.is_finite() and yT8 > 0
        if yr.upper() < YW_LO - 1 or yr.lower() > YW_HI + 1: continue                      # slack 1 >> |e| + frame shift
        ph_, kh, Sn, Cs, Yb = TC.base_series(yT8, phi, [c['k0'], c['k1'], c['k2'], c['k3']], ratio)
        Qb = TC.base_Q(Sn, Cs, Yb)
        TC.BOX[:] = [ebox(yT8.upper()) / yT8.lower(), BX[1], BX[2], BX[3], None]
        Qt, Mt, Yt = TC.true_Q(kh, Sn, Cs, Yb)
        wM = TC.enclose(Mt).sqrt() / TC.enclose(Yt)         # sqrt(M_t)/Yhat_t
        row = dict(slab=c['slab'], a=xlo, ca=c['a'], yr=(float(yr.lower()), U(yr)), yT8hi=U(yT8), K={}, Kl={}, ea={}, eb={},
                   tan=U(abs(phi.sin())) / float(phi.cos().lower()))
        for j in range(4):
            q = Qt[j]
            assert set(q.l) | set(q.r) <= set(range(j + 2)) and Qb[j][0].overlaps(q.c)
            for i in range(1, j + 2):
                li = TC.absu(q.l.get(i, arb(0)))
                row['K'][(j, i - 1)] = U((li + q.r.get(i, arb(0))) * wM); row['Kl'][(j, i - 1)] = U(li * wM)
            row['ea'][j] = U((TC.absu(q.l.get(0, arb(0))) + q.r.get(0, arb(0))) / yT8)     # first bracket, per |e|
            row['eb'][j] = U(TC.absu(Qb[j][1]) / yT8)                                        # second bracket, per |e| and per tan
        rows.append(row)
TANS = max(r['tan'] for r in rows)                          # sup|sin phi| / inf cos phi over the window cells
kch = kchl = 0.0; worstw = None; c0ch = {}
for r in rows:
    for (j, m), v in r['K'].items():
        if (j, m) in S:
            if v / S[(j, m)] > kch: kch, worstw = math.nextafter(v / S[(j, m)], math.inf), (r['slab'], r['a'], j, m)
            kchl = max(kchl, math.nextafter(r['Kl'][(j, m)] / S[(j, m)], math.inf))
        else:
            c0ch[j] = max(c0ch.get(j, 0.0), v)
# e-part per window bin: window starts in [yw, 1.05 yw] lie in [yw, 1.575 yw]; sqrt(int (28/y^2) dy) <= sqrt(28/(3 yw))
c0e, c0tip = {j: 0.0 for j in range(4)}, {j: 0.0 for j in range(4)}; yw = YW_LO; nbin = 0
while yw <= 253:
    sel = [r for r in rows if r['yr'][1] >= yw - 1 and r['yr'][0] <= min(1.575 * yw, YW_HI) + 1]
    assert sel
    wf = (arb(28) / (3 * arb(yw))).sqrt(); H = (arb(max(r['yT8hi'] for r in sel)) ** 3 / 84).sqrt()
    for j in range(4):
        co = wf * (arb(max(r['ea'][j] for r in sel)) + arb(max(r['eb'][j] for r in sel)) * arb(TANS))
        c0tip[j] = max(c0tip[j], U(co * CE)); c0e[j] = max(c0e[j], U(co * H))
    yw *= 1.05; nbin += 1
out['KCH'] = kch; out['KCH_lin'] = kchl; out['KCH_nl_factor'] = U(arb(kch) / arb(kchl)); out['KCH_where'] = worstw
out['c0_chain'] = dict(c0ch); out['c0_e'] = dict(c0e)
out['c0'] = {0: c0e[0], **{j: U(arb(c0ch[j]) + arb(c0e[j])) for j in range(1, 4)}}       # m = 0 chain term + e-term, per sqrt(F0)
out['c0tip'] = dict(c0tip)                                                                # e-term via the tip box, per unit Ybar
out['phi_max'] = max(U(abs(c['th'] - TH)) for c in cells if any(r['slab'] == c['slab'] and r['ca'] == c['a'] for r in rows))
out['b5'] = dict(ratio=RATIO, pieces=len(rows), piece=PIECE, bins=nbin, tan_star=TANS, boxes=dict(Y=str(N3A['Y']), b0=N3A['b0'], b1=N3A['b1'], b2=N3A['b2'],
                 c_e=TIPB['c_e'], sqrtF0_at_Y=U(SF0BOX)))
print(f"KCH <= {out['KCH']:.6f} (linear part {kchl:.6f}, computed nonlinear factor {out['KCH_nl_factor']:.6f}) at (slab, x, j, m) = {worstw}; "
      f"sup|phi| = {out['phi_max']:.4f};  c0_j = " + ', '.join(f'{j}: {v:.3e}' for j, v in sorted(out['c0'].items())) +
      ";  c0tip_j = " + ', '.join(f'{j}: {v:.3e}' for j, v in sorted(c0tip.items())))
out['args'] = dict(amp=sys.argv[1], window=amp['handoff']['window_slabs'], SLABS=list(SLABS), b0=B0, n3=N3P, t7const=T7C, env={k: os.environ[k] for k in ('COHERENT', 'S1PROF', 'S1TT', 'H1TAG', 'H1NSUB', 'H1REF', 'TREF', 'PDOT', 'T9TC', 'T9IFT', 'WINOUT', 'SUBW', 'T9JSON', 'TOPLOSS', 'N3JSON', 'T9OUT') if k in os.environ})
json.dump(out, open(os.environ.get('T9OUT', '/home/claude/s1spec/step4/p2/out/t9_cone_kch.json'), 'w'), indent=1)
