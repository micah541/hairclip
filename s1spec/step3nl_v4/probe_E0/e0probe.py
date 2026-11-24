"""[probe, est] effect of the E0 beta-loop constants on the tail anchor and first cell (side finding, not part of v4).
PI3: 0.01607 -> n3*BETA0_IBP = 7.52*(lead+dlead); f572: 5.72 -> 1 + 7.52*BETAF_IBP.  usage: python3 e0probe.py a_star L_star a1 L1lo L1hi"""
import sys, os, json, time
sys.path.insert(0, '/home/claude/s1spec/step3nl_v4/verifier')
import constants as K, rhs_iv as R, bb
from flint import arb
a_s, L_s, a1, L1lo, L1hi = map(float, sys.argv[1:6])
out = {}
for variant in ('v4', 'E0fix'):
    if variant == 'E0fix':
        K.PI3 = float(7.52 * (K.BETA0)); R._CC = None; c = R.consts(); c.cl['f572'] = R._c(1 + 7.52 * K.BETAF)
    t = time.time()
    r1 = bb.bound_cell(0.0, a_s, L_s, L_s * (1 + K.ETA_ANCHOR), tol=5e-4, rtol=0.005, floors=False, cem_iters=25)
    r2 = bb.bound_cell(0.0, a1, L1lo, L1hi, tol=5e-4, rtol=0.005, floors=True, cem_iters=25)
    out[variant] = dict(PI3=K.PI3, f572=float(R.consts().cl['f572'].hi), tail_anchor_U=r1['U'], first_cell_U=r2['U'], secs=round(time.time() - t, 1))
    print(variant, out[variant], flush=True)
json.dump(out, open('e0probe.json', 'w'), indent=1)
