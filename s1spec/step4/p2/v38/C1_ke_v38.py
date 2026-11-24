#!/usr/bin/env python3
"""C1_ke.py -- [cert] the normal-graph factor a = 1 - k~ e of step 4's P2 representation stays >= 1/2 on the P2 boxes (fix round 14).
[v38 copy of <s1spec>/writeup/lemmaU/bulk/fix15_review/c/C1_ke.py for the step-4 v38 re-pin (p2/v38/final_tail_v38.py); only change: output path out_v38/iface/C1_ke.json]

Used in C1_identification.tex, Lemma "curvature bound on P2's interval": under the conclusions of Theorem P2 (boxes B u B' hold), the true
curve is the normal graph G~ + e nu~ with |e| <= eta_1 = 1e-5 for y < 28 (and on the tip x <= 34.25), and |e| <= 1.7 E(x),
E(x) = c_e Ybar + sqrt(F0) (y^3/84)^(1/2) for y >= 28 (step4_proof.tex, Objects, "Boxes"), y = x - X0, sqrt(F0) <= aF Ybar + bF per slab
(amp_A1_w180.json, the per-slab amplification of Theorem P2's far bound).  This script bounds sup |k~| |e| over every stored cell of the
corrected base (p2/t6c/jets, hulls over the whole slab, widened outward by 2 ulps as t6/reader.py does), slabs 2..180, with
Ybar = Ybar_box = 1.7e-12 (the box level of N3', above the Ybar = 1e-12 in use: conservative).
Read-only: reads step-4 files, writes only C1_ke.json next to this script.  Run: nice -n 19 python3 -B C1_ke.py"""
import sys, os, json, hashlib, platform, math
sys.dont_write_bytecode = True
import numpy as np
from flint import arb, ctx
ctx.prec = 128
HERE = os.path.dirname(os.path.abspath(__file__))
S4 = '/home/claude/s1spec/step4'
JETS = f'{S4}/p2/t6c/jets'
AMP = f'{S4}/p2/t8se/constants/amp_A1_w180.json'
U = lambda v: math.nextafter(float(v.upper()), math.inf)


def sha(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()


def widen(v):
    v = v.copy(); v[..., 0] = np.nextafter(np.nextafter(v[..., 0], -np.inf), -np.inf); v[..., 1] = np.nextafter(np.nextafter(v[..., 1], np.inf), np.inf)
    return v


amp = json.load(open(AMP)); ps = amp['per_slab']; slabs = amp['slabs']
Y = arb('1.7e-12'); CE = arb('1.6031e4'); X0 = arb('2.2490776515'); ETA1 = arb('1e-5'); YCORE = 28.0 + 2.2490776515
inputs = {AMP: sha(AMP)}
worst = dict(ke=0.0); worst_e = 0.0; per = []
for i, s in enumerate(slabs):
    path = f'{JETS}/jets_{s:03d}.npz'; inputs[path] = sha(path)
    z = np.load(path); xh = widen(z['x']).reshape(-1, 2); kh = widen(z['k0']).reshape(-1, 2)
    sF0 = arb(repr(ps['0']['aF'][i])) * Y + arb(repr(ps['0']['bF'][i]))
    kmax = np.maximum(np.abs(kh[:, 0]), np.abs(kh[:, 1]))
    best = 0.0
    for (xl, xu), km in zip(xh, kmax):
        eb = ETA1
        if xu >= YCORE:                                   # cell meets y >= 28: the far box (y at its upper end, E increasing in y)
            yu = arb(xu) - X0
            Ef = arb('1.7') * (CE * Y + sF0 * (yu * yu * yu / 84).sqrt())
            eb = Ef if xl >= YCORE else arb(max(U(Ef), U(ETA1)))
        ke = U(arb(km) * eb); best = max(best, ke); worst_e = max(worst_e, U(eb))
        if ke > worst['ke']: worst = dict(ke=ke, slab=s, x=[float(xl), float(xu)], k=float(km), e=U(eb))
    per.append([s, best])
ok = worst['ke'] < 0.5
print('sup |k~ e| over slabs 2..180:', worst, ' sup |e|_box:', worst_e, ' ok (< 1/2):', ok)
me = os.path.abspath(__file__)
out = dict(label='[cert]', ok=bool(ok), sup_ke=worst['ke'], worst=worst, sup_e_box=worst_e, per_slab_sup_ke=per,
           statement='sup over slabs 2..180, all stored cells x in [0, ell~], of |k~| * (box bound of |e|) < 1/2; hence a = 1 - k~ e >= 1/2 and J >= 1/2',
           args=dict(argv=sys.argv, env={}, python=platform.python_version(), Ybar=str(Y), c_e=str(CE), X0=str(X0), eta1=str(ETA1),
                     script=dict(path=me, sha256=sha(me)), inputs=inputs))
json.dump(out, open('/home/claude/s1spec/step4/p2/out_v38/iface/C1_ke.json', 'w'), indent=1)
