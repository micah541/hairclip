"""does the per-slab tilt drift follow the polish tolerance?  slab 0 of P11v2 from the exact seed, Newton to em < TOL"""
import sys, pickle, time, numpy as np
from flint import arb, ctx
from base2 import Base
from arbops import SlabPolish, ArbBasis, ArbOp
from polstore import dec
from seed_p1 import closure
TOL = float(sys.argv[1])
B = Base('P11v2'); RD = B.RD; R, D = B.rep(0); q = B.slabs[0]
SP = SlabPolish(R, D, RD, q['nodes'], q['dt'])
ctx.prec = 128; sd = pickle.load(open('seed_p1_95.pkl', 'rb')); SP.Yn[0] = [dec(t) for t in sd['y_enc']]
J = D.jac(q['nodes'][len(q['nodes']) // 2]); RD.factor(J, q['dt']); t0 = time.time()
for it in range(200):
    em, dz = SP.step(J)
    if em < TOL or dz < 1e-40: break
ctx.prec = 192; AB = ArbBasis(R.B, 192); op = ArbOp(R, AB)
out = []
for i in range(len(SP.Yn)):
    tq, G = closure(R, AB, op, [v.mid() for v in SP.Yn[i]]); out.append((float((tq + arb.pi()).mid()), float(G.imag.mid())))
print(f'TOL {TOL:.0e}: it {it+1} em {em:.2e} dz {dz:.2e} ({time.time()-t0:.0f}s); tilt at nodes ' + ' '.join(f'{a:.2e}' for a, b in out) + '; offset end %.4e' % out[-1][1])
