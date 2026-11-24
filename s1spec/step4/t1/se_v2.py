"""(SE) far-end check of a stored base (polstore format), in arb at 192 bits on the stored midpoints:
tilt = theta~(ell~) + 3pi/2 = int_0^ell kappa~ + pi,  offset = Im Gamma~(ell~)  (far point vs the axis), at every stored row
(relative times 0, c_1..c_8) of every slab.  usage: python3 se_v2.py DIR [i0:i1] -> se_<DIR>_<i0>_<i1>.json"""
import sys, os, json, pickle
from flint import arb, ctx
from polstore import dec
from rep import Rep, Mesh
from arbops import ArbBasis, ArbOp
from seed_p1 import closure
DIR = sys.argv[1]; man = json.load(open(os.path.join(DIR, 'manifest.json')))
i0, i1 = (map(int, sys.argv[2].split(':')) if len(sys.argv) > 2 else (0, len(man)))
AB = None; out = []
for me in man[i0:i1]:
    d = pickle.load(open(os.path.join(DIR, me['file']), 'rb')); cfg = d['cfg']
    R = Rep(cfg['P'], Mesh(JT=cfg['JT'], h0=cfg['h0'], cT=cfg['cT'], C=d['C']))
    ctx.prec = 192
    if AB is None: AB = ArbBasis(R.B, 192)
    op = ArbOp(R, AB); rec = dict(slab=d['slab'], k0=me['k0'], k1=me['k1'])
    radmax = 0.0
    rows = []
    for row in range(len(d['Yn'])):
        Y = [dec(t) for t in d['Yn'][row]]; radmax = max(radmax, max(float(v.rad()) for v in Y))
        y = [v.mid() for v in Y]
        tq, G = closure(R, AB, op, y)
        lam = y[R.n + 1].exp()
        rows.append(dict(theta=d['nodes_rel'][row] if len(d['nodes_rel']) == len(d['Yn']) else ([0.0] + d['nodes_rel'])[row],
                         tilt=float((tq + arb.pi()).mid()), tilt_rad=float((tq + arb.pi()).rad()), offset=float(G.imag.mid()),
                         offset_rad=float(G.imag.rad()), offset_phys=float((G.imag * lam).mid()), ReG=float(G.real.mid()), ell=float(y[R.n].mid())))
    rec['start'], rec['end'] = rows[0], rows[-1]; rec['rows'] = rows
    rec['tilt_stage_max'] = max(abs(r['tilt']) for r in rows[1:-1]); rec['offset_max'] = max(abs(r['offset']) for r in rows)
    rec['Yn_rad_max'] = max(radmax, max(float(dec(t).rad()) for row in d['Yn'] for t in row))
    out.append(rec); print(json.dumps(rec), flush=True)
json.dump(out, open(f'se_{os.path.basename(os.path.normpath(DIR))}_{i0}_{i1}.json', 'w'), indent=1)
