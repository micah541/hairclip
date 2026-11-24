"""v3.3 task 6, (H2f) stage 1: parameter balls for the (H2f) base sums, widened so that the base value covers the TRUE curve's hand-off
parameters.  For each of the 32 slab-180 t_c balls (p2/out/h1slab/ift_m.json, pinned v3.2, read only):
  xc, lam  widened by DSH = 2 max(C_tr, C_sc) YBOX   (exactly t9.py's true-curve frame shift: clam's transfer M^-T, f_norms, tail t;
                                                      factor 2 covers lam vs ln lam);
  alpha    widened by C_alpha YBOX                   (clam: |alpha - alpha~| <= C_alpha Ybar);
  Theta    unchanged (the S1 ray angle; not a curve parameter).
YBOX = the n3 bootstrap Ybar (n3_b1.json args Y = 1.7e-12) >= Ybar = 1e-12: the widened balls hold for every Ybar <= YBOX.
So the base sums of h2f_win.py, evaluated on these balls, enclose G_j of the base curve in every frame/alpha the true curve can have;
the error part (h2f_camp.py) then only transfers e at a fixed frame (as t9.py/camp.py do for the (H3) windows).
Writes p2/out_v33/h2f/ift_m.json (m = 0..31) and p2/out_v33/h2f/ift_widen.json (args block).
usage: python3 -I h2f_ift.py"""
import sys, os, json, math, hashlib, glob
from flint import arb, arb_mat, ctx
ctx.prec = 128
P2 = '/home/claude/s1spec/step4/p2'; O = f'{P2}/out'; OV = f'{P2}/out_v33/h2f'
os.makedirs(OV, exist_ok=True)
sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()
U = lambda v: math.nextafter(float(v.upper()), math.inf); Lo = lambda v: math.nextafter(float(v.lower()), -math.inf)
CL = f'{O}/t9_h1_clam_s180.json'; N3 = f'{P2}/t8se/constants/n3_b1.json'
_cl = json.load(open(CL)); YBOXS = str(json.load(open(N3))["args"]["Y"]); YBOX = arb(YBOXS)
_Mi = arb_mat([[arb(s_) for s_ in row] for row in _cl['M']]).transpose().inv()                 # as p2/t9.py (frame block)
_fn = [arb(repr(v)) for v in _cl['f_norms']]; _tt = arb(repr(_cl['tail_t']))
CTR = [sum((abs(_Mi[i, j]) * _fn[j] for j in range(3)), arb(0)) * (1 + _tt) for i in range(3)]
DSH = 2 * max(CTR[0], CTR[1], key=U) * YBOX
CA = arb(repr(_cl['C_alpha'])); DAL = CA * YBOX
ins = [CL, N3]
for m in range(32):
    src = f'{O}/h1slab/ift_{m}.json'; d = json.load(open(src)); ins.append(src)
    w = dict(alpha=[Lo(arb(repr(d['alpha'][0])) - DAL), U(arb(repr(d['alpha'][1])) + DAL)],
             xc=[Lo(arb(repr(d['xc'][0])) - DSH), U(arb(repr(d['xc'][1])) + DSH)],
             lam=[Lo(arb(repr(d['lam'][0])) - DSH), U(arb(repr(d['lam'][1])) + DSH)],
             Theta=d['Theta'], tc_ball=d['tc_ball'],
             args=dict(ball=m, source=src, source_sha256=sha(src), widened_by=dict(frame=U(DSH), alpha=U(DAL)), YBOX=YBOXS))
    json.dump(w, open(f'{OV}/ift_{m}.json', 'w'), indent=1)
out = dict(DSH=U(DSH), C_tr_sc=[U(CTR[0]), U(CTR[1])], C_alpha=_cl['C_alpha'], dalpha=U(DAL), YBOX=YBOXS,
           args=dict(argv=sys.argv[1:], env={}, inputs={p: sha(p) for p in ins}, outputs=[f'{OV}/ift_{m}.json' for m in range(32)]))
json.dump(out, open(f'{OV}/ift_widen.json', 'w'), indent=1)
print(f'widened 32 balls: frame +-{U(DSH):.3e}, alpha +-{U(DAL):.3e} (YBOX {YBOXS})')
