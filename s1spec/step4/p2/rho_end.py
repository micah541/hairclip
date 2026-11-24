"""rigorous sup over each slab of |rho_theta(ell)| (angle residual at the far end; T6 engine, exact in s, TP in tc) and of
|theta(ell) + pi| (turning defect) for the v2 base: propagates Q, rho_theta, theta through all elements as T6's slab.py does.
usage: python3 rho_end.py a-b   -> out/rhoend/slab_NNN.json (restartable)"""
import sys, os, json, time
S4 = '/home/claude/s1spec/step4'
sys.path.insert(0, f'{S4}/t6'); sys.path.insert(0, f'{S4}/t1')
from flint import arb, ctx
import tp as T
from base import SlabBase
from engine import element_a
ONE = arb(1)
a, b = map(int, sys.argv[1].split('-'))
for i in range(a, b + 1):
    out = f'{S4}/p2/out/rhoend/slab_{i:03d}.json'
    if os.path.exists(out): continue
    t0 = time.process_time()
    B = SlabBase(i, poldir=f'{S4}/t1/polished_v2'); ctx.prec = 128
    Z = T.tconst(arb(0)); Qe, rte, the = Z, Z, T.tconst(-arb.pi() / 2)
    for e in range(B.N):
        A_ = element_a(B, e, Qe, rte, the)
        Qe, rte, the = A_['Q'].at(ONE), A_['rt'].at(ONE), A_['th'].at(ONE)
    sup_rt = T.tsup_sub(rte, K=32); turn = the + arb.pi() / 2          # theta(ell) - pi/2 (mod 2pi; -3pi/2 expected)
    json.dump(dict(slab=i, sup_rho_theta_end=float(sup_rt.upper()), rte_coeffs=[c.str(5, radius=True) for c in rte.coeffs()[:3]],
                   turn_end_range=T.trange(turn).str(8, radius=True), cpu=time.process_time() - t0), open(out, 'w'))
    print(i, f'{float(sup_rt.upper()):.3e}', f'{time.process_time() - t0:.1f}s', flush=True)
