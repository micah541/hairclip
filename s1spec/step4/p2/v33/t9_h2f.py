"""v3.3 task 6, (H2f) stage 3: the error-transfer constants of p2/t9.py (KCH, c0, c0tip; review B5 block, copied verbatim except the
y-range and the e-term summation) on y_h in [60, 700] instead of the (H3) window union [50, 379.5].
Copied from p2/t9.py (closed v3.2 step 4; read only, not imported because it runs at module level): the T6 jet cells of slabs 177-180,
the n3 boxes, T7's tip box c_e, the hand-off frame (32 IFT balls + clam's true-curve shift), the tip anchor for y_ray/y_T8, the exact
decomposition (*) of p2/t9chart.py with Taylor mode, and per piece K_jm, ea_j, eb_j.  Changed:
  * YW_LO, YW_HI = 60, 700 (t9.py: 50, 379.5); the cell filter x <= 2 YW_HI and the tip anchor's range follow it.
  * KCH_f = max over the pieces with y_h in [59, 701] of K_jm / S(j, m), (j, m) in S; c0ch_f[j] = max K_j0 (j >= 1) on the same pieces
    (C3 of LEMMA_F holds on any y-interval: F_m are the global half-curve energies, T_m >= 28 y^(2m-2) up to y_l ~ 3.6e5 >> 700).
  * e-term: one interval [60, 700] instead of sliding windows, so the bins are summed in L^2 (Minkowski) rather than maximised:
      bins [y_k, y_k+1], y_k = 60 * 1.05^k (last bin cut at 700); per bin co_j = sup ea_j + sup eb_j * tan*, H_k = (yT8hi^3/84)^(1/2);
      sqrt(int_60^700 (28/y^2) (co |e|)^2 dy) <= c0tip_f[j] Ybar + c0e_f[j] sqrt(F0),
      c0tip_f[j] = c_e (sum_k 28 (1/y_k - 1/y_k+1) co_kj^2)^(1/2),  c0e_f[j] = (sum_k 28 (1/y_k - 1/y_k+1) co_kj^2 H_k^2)^(1/2)
    using |e| <= c_e Ybar + sqrt(F0) (yT8^3/84)^(1/2) (T7 Lemma 4b + LEMMA_F C4 core; y >= 28, valid up to y_l).
  * Outputs only these constants (the cone jets and shells of t9.py are not recomputed here).
Validity range of the inputs (checked for x <= 730, as CONTRACTS_r2 task 6 asks): the n3 boxes |e'| <= b0, x|e''| <= b1, x^2|e'''| <= b2
are asserted by n3boxes.py with y up to y_l = ell - X0 (ell ~ 3.6e5, constants_b1); the |e| bound is the core integral on [28, y], any y;
the (SE) correction changes the base only on x >= 0.6 ell ~ 2.1e5.  So every box holds on [60, 700].
usage: python3 -I t9_h2f.py AMP_JSON        (env T6DIR, N3JSON as t9.py; T9OUT: output, default p2/out_v33/t9_h2f.json)"""
import sys, os, json, math, hashlib, glob, bisect
import numpy as np
from flint import arb, arb_mat, ctx
ctx.prec = 128
sys.path.insert(0, '/home/claude/s1spec/step4/t6')
import reader as RD
sys.path.insert(0, '/home/claude/s1spec/step4/p2'); import t9chart as TC
sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()
N3P = os.environ.get('N3JSON', '/home/claude/s1spec/step4/p2/t8se/constants/n3_b1.json')
N3A = json.load(open(N3P))['args']
T7C = '/home/claude/s1spec/step4/t7/constants.json'
JD = os.path.join(os.environ.get('T6DIR', '/home/claude/s1spec/step4/p2/t6c'), 'jets')
SLABS = range(177, 181)
X0 = arb('2.2490776515'); TH = arb('1.8617364190179404', '1e-3')
amp = json.load(open(sys.argv[1]))
U = lambda v: math.nextafter(float(v.upper()), math.inf)
def iv(lo, hi): return arb(lo).union(arb(hi))

cells = []
for s in SLABS:
    J = RD.jets(s, d=JD)
    for e in range(J['x'].shape[0]):
        for c in range(J['x'].shape[1]):
            g = lambda nm: iv(*J[nm][e, c])
            cells.append(dict(slab=s, a=float(J['x'][e, c, 0]), b=float(J['x'][e, c, 1]), k0=g('k0'), k1=g('k1'), k2=g('k2'), k3=g('k3'),
                              th=g('theta'), p=g('p'), q=g('q')))

S = {(1, 1): 1, (2, 1): 1, (2, 2): 1, (3, 1): 1, (3, 2): 3, (3, 3): 1, (0, 0): 1}
YBOX = arb(str(N3A['Y']))
TIPB = json.load(open(T7C))['meta']['tip_box']
assert TIPB['holds'] and abs(TIPB['Ybox'] - float(YBOX.mid())) <= 1e-6 * float(YBOX.mid()), 'T7 tip box at another Ybar'
CE = arb(repr(TIPB['c_e']))
_ps, _sl = amp['per_slab'], amp['slabs']
SF0BOX = max((arb(repr(_ps['0']['aF'][_sl.index(s)])) * YBOX + arb(repr(_ps['0']['bF'][_sl.index(s)])) for s in SLABS), key=U)
BX = [None, arb(repr(N3A['b0'])), arb(repr(N3A['b1'])), arb(repr(N3A['b2'])), None]
def ebox(yT8):
    return CE * YBOX + SF0BOX * (arb(yT8) ** 3 / 84).sqrt()
YW_LO, YW_HI = 60.0, 700.0                                  # (H2f) range (t9.py: 50, 379.5)
VAL = os.environ.get('V33_VALIDATE') == '1'                 # scratch fidelity run: t9.py's range and levels 0-3 (compare KCH, c0_chain with t9_cone_kch.json)
if VAL: YW_LO, YW_HI = 50.0, 1.5 * 253
JN = 4 if VAL else 3
_ift = [json.load(open(f)) for f in sorted(glob.glob('/home/claude/s1spec/step4/p2/out/h1slab/ift_*.json'))]; assert len(_ift) == 32
_cl = json.load(open('/home/claude/s1spec/step4/p2/out/t9_h1_clam_s180.json'))
_Mi = arb_mat([[arb(s_) for s_ in row] for row in _cl['M']]).transpose().inv()
_fn = [arb(repr(v)) for v in _cl['f_norms']]; _tt = arb(repr(_cl['tail_t']))
CTR = [sum((abs(_Mi[i, j]) * _fn[j] for j in range(3)), arb(0)) * (1 + _tt) for i in range(3)]
DSH = 2 * max(CTR[0], CTR[1], key=U) * YBOX
XC = iv(min(d['xc'][0] for d in _ift), max(d['xc'][1] for d in _ift)) + DSH * arb(0, 1)
FL = iv(min(d['lam'][0] for d in _ift), max(d['lam'][1] for d in _ift)) + DSH * arb(0, 1)
assert all(TH.contains(iv(*d['Theta'])) for d in _ift)
FC = XC * TH.cos()
def tip_anchor(slab):
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
def I_at(T, x, hi):
    assert 0 <= x <= T['bp'][-1]; k = min(max(0, bisect.bisect_right(T['bp'], x) - 1), len(T['lb']) - 1)
    return (T['Ih'][k] + (arb(x) - arb(T['bp'][k])) * T['ub'][k]) if hi else (T['Il'][k] + (arb(x) - arb(T['bp'][k])) * T['lb'][k])
TIPA = {s: tip_anchor(s) for s in SLABS}
PIECE = float(os.environ.get('T9PIECE', '0.5'))
rows = []
for c in cells:
    if c['a'] > 2 * YW_HI or c['b'] < 40: continue
    phi = c['th'] - TH
    npc = max(1, math.ceil((c['b'] - c['a']) / PIECE)); xs = [c['a'] + (c['b'] - c['a']) * i / npc for i in range(npc)] + [c['b']]
    for xlo, xhi in zip(xs[:-1], xs[1:]):
        yT8 = iv(xlo, xhi) - X0
        T = TIPA[c['slab']]
        num = iv(float((T['yr0'] + X0 - I_at(T, xhi, True)).lower()), float((T['yr0'] + X0 - I_at(T, xlo, False)).upper()))
        ratio = 1 + num / yT8
        yr = (yT8 * ratio - FC) / FL
        assert yr.is_finite() and yT8 > 0
        if yr.upper() < YW_LO - 1 or yr.lower() > YW_HI + 1: continue
        ph_, kh, Sn, Cs, Yb = TC.base_series(yT8, phi, [c['k0'], c['k1'], c['k2'], c['k3']], ratio - FC / yT8)
        Qb = TC.base_Q(Sn, Cs, Yb)
        TC.BOX[:] = [ebox(yT8.upper()) / yT8.lower(), BX[1], BX[2], BX[3], None]
        Qt, Mt, Yt = TC.true_Q(kh, Sn, Cs, Yb)
        wM = (FL * TC.enclose(Mt)).sqrt() / TC.enclose(Yt)
        row = dict(slab=c['slab'], a=xlo, ca=c['a'], yr=(float(yr.lower()), U(yr)), yT8hi=U(yT8), K={}, Kl={}, ea={}, eb={},
                   tan=U(abs(phi.sin())) / float(phi.cos().lower()))
        for j in range(JN):                                  # (H2f): levels 0-2
            q = Qt[j]
            assert set(q.l) | set(q.r) <= set(range(j + 2)) and Qb[j][0].overlaps(q.c)
            for i in range(1, j + 2):
                li = TC.absu(q.l.get(i, arb(0)))
                row['K'][(j, i - 1)] = U((li + q.r.get(i, arb(0))) * wM); row['Kl'][(j, i - 1)] = U(li * wM)
            row['ea'][j] = U((TC.absu(q.l.get(0, arb(0))) + q.r.get(0, arb(0))) / yT8)
            row['eb'][j] = U(TC.absu(Qb[j][1]) / yT8)
        rows.append(row)
assert min(r['yr'][0] for r in rows) <= YW_LO and max(r['yr'][1] for r in rows) >= YW_HI, 'pieces do not cover [60, 700]'
TANS = max(r['tan'] for r in rows)
kch = kchl = 0.0; worstw = None; c0ch = {}
for r in rows:
    for (j, m), v in r['K'].items():
        if (j, m) in S:
            if v / S[(j, m)] > kch: kch, worstw = math.nextafter(v / S[(j, m)], math.inf), (r['slab'], r['a'], j, m)
            kchl = max(kchl, math.nextafter(r['Kl'][(j, m)] / S[(j, m)], math.inf))
        else:
            c0ch[j] = max(c0ch.get(j, 0.0), v)
# e-term on the single interval [60, 700]: bins summed in L^2
s_tip = {j: arb(0) for j in range(3)}; s_f0 = {j: arb(0) for j in range(3)}; yw = YW_LO; nbin = 0; bins = []
while yw < YW_HI:
    yn = min(yw * 1.05, YW_HI)
    sel = [r for r in rows if r['yr'][1] >= yw - 1 and r['yr'][0] <= yn + 1]
    assert sel
    wt = arb(28) * (1 / arb(yw) - 1 / arb(yn)); H = (arb(max(r['yT8hi'] for r in sel)) ** 3 / 84).sqrt()
    cos_ = []
    for j in range(3):
        co = arb(max(r['ea'][j] for r in sel)) + arb(max(r['eb'][j] for r in sel)) * arb(TANS)
        s_tip[j] += wt * co * co; s_f0[j] += wt * co * co * H * H; cos_.append(U(co))
    bins.append(dict(y=[yw, yn], co=cos_, H=U(H))); yw = yn; nbin += 1
if VAL:
    ref = json.load(open('/home/claude/s1spec/step4/p2/out/t9_cone_kch.json'))
    print('VALIDATE KCH', kch, ref['KCH'], kch == ref['KCH'], 'KCH_lin', kchl == ref['KCH_lin'], 'c0_chain', {j: c0ch[j] for j in c0ch}, ref['c0_chain']); sys.exit(0)
c0tip = {j: U(CE * s_tip[j].sqrt()) for j in range(3)}
c0e = {j: U(s_f0[j].sqrt()) for j in range(3)}
out = dict(label='certified (arb); constants of the (H2f) error transfer on y_h in [60, 700]',
           KCH=kch, KCH_lin=kchl, KCH_nl_factor=U(arb(kch) / arb(kchl)), KCH_where=worstw, c0_chain={j: c0ch.get(j, 0.0) for j in range(3)},
           c0_e=c0e, c0={0: c0e[0], **{j: U(arb(c0ch[j]) + arb(c0e[j])) for j in (1, 2)}}, c0tip=c0tip,
           tan_star=TANS, pieces=len(rows), bins=nbin, y_range=[YW_LO, YW_HI],
           frame=dict(xc=[float(XC.lower()), U(XC)], lam=[float(FL.lower()), U(FL)], shift=U(DSH)),
           boxes=dict(Y=str(N3A['Y']), b0=N3A['b0'], b1=N3A['b1'], b2=N3A['b2'], c_e=TIPB['c_e'], sqrtF0_at_Y=U(SF0BOX)), bin_table=bins,
           args=dict(argv=sys.argv[1:], SLABS=list(SLABS), jets=JD, n3=N3P, t7const=T7C, piece=PIECE,
                     env={k: os.environ[k] for k in ('T6DIR', 'N3JSON', 'T9OUT', 'T9PIECE') if k in os.environ},
                     inputs={p: sha(p) for p in [os.path.abspath(__file__), sys.argv[1], N3P, T7C, '/home/claude/s1spec/step4/p2/out/t9_h1_clam_s180.json',
                                                 '/home/claude/s1spec/step4/p2/t9chart.py'] + sorted(glob.glob('/home/claude/s1spec/step4/p2/out/h1slab/ift_*.json'))
                             + [f'{JD}/jets_{s:03d}.npz' for s in SLABS]}))
fo = os.environ.get('T9OUT', '/home/claude/s1spec/step4/p2/out_v33/t9_h2f.json')
json.dump(out, open(fo + '.tmp', 'w'), indent=1); os.replace(fo + '.tmp', fo)
print(f"(H2f) chart: KCH <= {kch:.6f} (linear {kchl:.6f}) at {worstw}; c0 = {out['c0']}; c0tip = {c0tip}; pieces {len(rows)}, bins {nbin}")
