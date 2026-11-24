"""v3.3 task 6, fixed-L cone family: certified cone jets (Lemma O H3, (H5) data) on the fixed-L cylinders of step 3 v3.3 (Lemma U,
task 5's Z5): unit L = 50, X = 3, centres r_c in [300, 380], cylinder [r_c - 150, r_c + 150] (half-width X L = 150), in L = 50 units on the
whole cylinder.  These replace v3.2's cone cylinders (unit r_c/4, centres [200, 380], cylinder [r_c/4, 7 r_c/4]).
Copied from p2/t9.py's cone block (closed v3.2; read only): the T6 jet cells of slabs 177-180, r = |Gamma~| per cell (any frame with
|x_c| <= 0.5, |lambda - 1| <= 1e-3 is covered by the slack 1), the C4 pointwise error bounds e2, e3, e4 at the cylinder's inner end, |e| from
the n3 tip box eta1 and the core integral, curvjets.kmajor (exact jets of base + e N), and the graph formulas over the centre tangent
  h'' = k sec^3, h''' = sec^4 (k_s + 3 k^2 t), h'''' = sec^5 (k_ss + 10 k k_s t + 12 k^3 t^2 + 3 k^3 sec^2),  |psi| <= sup|k| (XL + 2) 1.01,
  K_j = L^(j-1) sup|h^(j)|, N3 = L^1.5 ||h'''||_L2(|x| <= XL) <= sqrt(6) K3, slope = tan(psi).
Changed: rows r_c in [rlo, rlo + 10], rlo = 300, ..., 370; L = 50 (t9.py: r_hi/4); cells with r in [rlo - 150 - 1, rlo + 10 + 150 + 1]
(t9.py: [r_lo/4 - 1, 7 r_hi/4 + 1]); half-length X L + 2 = 152 (t9.py: 0.75 r_hi + 2); error bounds at y = rlo - 151 - 2.25 (t9.py:
r_lo/4 - 1 - 2.25); |e| at y = rlo + 161.  Checks (thresholds from Z5_REPORT / CONTRACTS_r2, values as v3.2): (H5) slope <= 0.08,
K2 = sup|h0''| <= 0.03, N3 = ||h0'''||_L2 <= 0.004; Lemma O H3: K2, K3, K4 <= M = 0.12.  margin = threshold / bound (rounded down).
usage: python3 -I t9_cone_v33.py AMP_JSON   (env T6DIR, N3JSON as t9.py; T9OUT, default p2/out_v33/t9_cone_v33.json)"""
import sys, os, json, math, hashlib
from flint import arb, ctx
ctx.prec = 128
sys.path.insert(0, '/home/claude/s1spec/step4/t6')
import reader as RD
sys.path.insert(0, '/home/claude/s1spec/step4/p2'); import curvjets as CJ
sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()
N3P = os.environ.get('N3JSON', '/home/claude/s1spec/step4/p2/t8se/constants/n3_b1.json')
N3A = json.load(open(N3P))['args']; B0 = N3A['b0']
JD = os.path.join(os.environ.get('T6DIR', '/home/claude/s1spec/step4/p2/t6c'), 'jets')
SLABS = range(177, 181)
amp = json.load(open(sys.argv[1])); hb = amp['handoff']['at_Ybar']
sF = [arb(repr(v)) for v in hb['sqrtF']]; F = [v * v for v in sF]; D2 = arb(repr(hb['D2_goodtime']))
U = lambda v: math.nextafter(float(v.upper()), math.inf); Lw = lambda v: math.nextafter(float(v.lower()), -math.inf)
def iv(lo, hi): return arb(lo).union(arb(hi))
def absb(a): return arb(U(abs(a)))
cells = []
for s in SLABS:
    J = RD.jets(s, d=JD)
    for e in range(J['x'].shape[0]):
        for c in range(J['x'].shape[1]):
            g = lambda nm: iv(*J[nm][e, c])
            cells.append(dict(slab=s, k0=g('k0'), k1=g('k1'), k2=g('k2'), k3=g('k3'), p=g('p'), q=g('q')))
for c in cells: c['r'] = (c['p'] ** 2 + c['q'] ** 2).sqrt()
rr = arb('1.5') / arb('0.8828'); c2 = rr * rr; c3 = c2 * c2
def err(y):                                                  # as p2/t9.py
    y = arb(y)
    e2 = ((arb(3) / 28 * F[1] + 2 * rr / 28 * (F[1] * F[2]).sqrt()) / y).sqrt()
    e3 = ((3 * c2 / 28 * F[2] + 2 * (c2 * c3).sqrt() / 28 * (F[2] * F[3]).sqrt()) / y ** 3).sqrt()
    e4 = (3 * c3 * F[3] / (28 * y ** 5) + 2 * (c2 * c3).sqrt() / 28 * (F[3] * D2).sqrt() / y ** 3).sqrt()
    return e2, e3, e4
LU, XU = 50, 3; XL = LU * XU                                  # fixed-L family (v3.3)
out = dict(cone=[], geometry=dict(L=LU, X=XU, centres=[300, 380], cylinder='[r_c - 150, r_c + 150]', rows='r_c in [rlo, rlo + 10], rlo = 300..370'))
worst = dict(K2=0.0, K3=0.0, K4=0.0, N3=0.0, slope=0.0)
VAL = os.environ.get('V33_VALIDATE') == '1'                  # scratch fidelity run: v3.2 geometry, compare rows with t9_cone_kch.json
for rlo in (range(200, 380, 10) if VAL else range(300, 380, 10)):
    rhi = rlo + 10; L = arb(rhi) / 4 if VAL else arb(LU)
    lo, hi = (0.25 * rlo - 1, 1.75 * rhi + 1) if VAL else (rlo - XL - 1, rhi + XL + 1)
    HL = arb('0.75') * rhi if VAL else arb(XL)
    sel = [c for c in cells if float(c['r'].upper()) >= lo and float(c['r'].lower()) <= hi]
    assert sel and (VAL or min(float(c['r'].lower()) for c in sel) <= lo and max(float(c['r'].upper()) for c in sel) >= hi)
    k0 = max((absb(c['k0']) for c in sel), key=U); k1 = max((absb(c['k1']) for c in sel), key=U); k2 = max((absb(c['k2']) for c in sel), key=U)
    e2, e3, e4 = err(max(lo - 2.25, 40.0))
    k3 = max((absb(c['k3']) for c in sel), key=U)
    yhi = arb(hi); E0 = arb(repr(N3A['eta1'])) + F[0].sqrt() * (yhi ** 3 / 84).sqrt()
    K0, K1, K2_ = (arb(v) for v in CJ.kmajor([arb(0, U(k0)), arb(0, U(k1)), arb(0, U(k2)), arb(0, U(k3))], [U(E0), B0, U(e2), U(e3), U(e4)], 3))
    psi = K0 * (HL + 2) * arb('1.01'); t = psi.tan(); sec = 1 / psi.cos()
    assert U(sec) <= 1.01
    h2 = K0 * sec ** 3; h3 = sec ** 4 * (K1 + 3 * K0 ** 2 * t)
    h4 = sec ** 5 * (K2_ + 10 * K0 * K1 * t + 12 * K0 ** 3 * t * t + 3 * K0 ** 3 * sec * sec)
    row = dict(r=[rlo, rhi], support_r=[lo, hi], L=U(L), slope=U(t), K2=U(L * h2), K3=U(L ** 2 * h3), K4=U(L ** 3 * h4), N3=U(arb(6).sqrt() * L ** 2 * h3),
               base_k=[U(k0), U(k1), U(k2), U(k3)], err=[U(e2), U(e3), U(e4)], E0=U(E0), err_at_y=max(lo - 2.25, 40.0))
    out['cone'].append(row)
    for k in worst: worst[k] = max(worst[k], row[k])
out['cone_max'] = worst
if VAL:
    ref = json.load(open('/home/claude/s1spec/step4/p2/out/t9_cone_kch.json'))
    print('VALIDATE rows identical:', all(all(a[k] == b[k] for k in ('slope', 'K2', 'K3', 'K4', 'N3', 'L')) for a, b in zip(out['cone'], ref['cone'])), len(out['cone']), len(ref['cone']), worst, ref['cone_max']); sys.exit(0)
TH = dict(slope=('0.08', 'H5 s0'), K2=('0.03', 'H5 m2 (also H3 M = 0.12)'), N3=('0.004', 'H5 n3'), K3=('0.12', 'H3 M'), K4=('0.12', 'H3 M'))
chk = {}
for k, (thr, name) in TH.items():
    v, tt = arb(repr(worst[k])), arb(thr)
    chk[k] = dict(item=name, value=worst[k], threshold=float(tt.mid()), margin=Lw(tt / v), ok=bool(v < tt))
chk['K2_vs_M'] = dict(item='H3 M (j = 2)', value=worst['K2'], threshold=0.12, margin=Lw(arb('0.12') / arb(repr(worst['K2']))), ok=bool(arb(repr(worst['K2'])) < arb('0.12')))
out['checks'] = chk; out['all_ok'] = all(v['ok'] for v in chk.values())
out['args'] = dict(argv=sys.argv[1:], SLABS=list(SLABS), b0=B0, n3=N3P, jets=JD, env={k: os.environ[k] for k in ('T6DIR', 'N3JSON', 'T9OUT') if k in os.environ},
                   inputs={p: sha(p) for p in [os.path.abspath(__file__), sys.argv[1], N3P, '/home/claude/s1spec/step4/p2/curvjets.py'] + [f'{JD}/jets_{s:03d}.npz' for s in SLABS]})
fo = os.environ.get('T9OUT', '/home/claude/s1spec/step4/p2/out_v33/t9_cone_v33.json')
json.dump(out, open(fo + '.tmp', 'w'), indent=1); os.replace(fo + '.tmp', fo)
print('fixed-L cone (L = 50, X = 3, r_c 300-380): ' + ', '.join(f"{k} {v['value']:.4e} (x{v['margin']:.2f})" for k, v in chk.items()))
assert out['all_ok'], 'fixed-L cone check FAILS'
