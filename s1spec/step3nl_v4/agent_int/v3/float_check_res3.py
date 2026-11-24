"""float check of F1b: central finite difference (step 1e-4, 300-bit) of the full residuals R(x) (tight-tail profile at sigma = x + ST.mid) vs (i) symbolic R'(x)
and (ii) the rigorous box bound of the derivative-level residual (boxes of the final run)."""
import sys, json; sys.argv = ['x', '1']
from int_common import *
import resid_sym as RS, int_res3 as M
from c4fun import PSI, H2
from flint import ctx
def R_at(w, x, order):
    with ctx.workprec(300):
        X = arb(x); y = PSI.jets_any(X, x, 7); h = H2.jets_any(X, x, 7)
        sg = arb(x) + ST.mid(); p, q, k, m = [t.mid() for t in HR.tail_state(sg)]
        return RS.call(w, order, y, h, p, q, k, m, M.LAMS[w]).mid()
res = json.load(open(sys.argv[1] if False else 'v3_res_boxes_0.5.json'))
for x in (10.5, 14.0, 18.0, 22.5, 27.9, 33.0):
    row = []
    for w in ('psi', 'h2'):
        d = 1e-4
        fd = (R_at(w, x + d, 0) - R_at(w, x - d, 0)) / (2 * d); sym = R_at(w, x, 1)
        bd = [r['r'][w][1] for r in res if r['xa'] <= x <= r['xb']][0]
        row.append(f"{w}: FD {float(fd):+.4e} sym {float(sym):+.4e} box-bound {bd:.3e}")
    print(f"x={x}: " + " | ".join(row))
