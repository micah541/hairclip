import sys, time, numpy as np
from base2 import Base; from arbops import *; import resid
tag = sys.argv[1]; i = int(sys.argv[2])
B = Base(tag); R, D = B.rep(i); RD = B.RD; q = B.slabs[i]; nodes = q['nodes']; dt = q['dt']
print(f'slab {i}: k = {B.k(q["tau0"]):.2f}, n = {R.n}, N = {R.N}')
SP = SlabPolish(R, D, RD, nodes, dt)
def show(lbl, thetas):
    rows = []
    for th in thetas:
        o, _ = resid.norms(resid.samples(SP, th)); rows.append(o)
    keys = ['R2', 'R1', 'Rs', 'RIs', 'RF0', 'RF1', 'RF2', 'RF3', 'RF4', 'supS2_tip', 'supr_tip']
    print(lbl, ' '.join(f"{k}={max(r[k] for r in rows):.2e}" for k in keys), flush=True)
    return rows
t0 = time.time(); show('float base, nodes  :', [0.0] + list(RD.c)); print('  (%.1fs)' % (time.time() - t0))
J = D.jac(nodes[len(nodes) // 2]); RD.factor(J, dt)
for it in range(40):
    em, dz = SP.step(J)
    if em < 1e-33 or dz < 1e-36: break
print(f'polished: {it+1} iterations, last max|E| = {em:.1e}')
show('polished, nodes    :', [0.0] + list(RD.c))
th = np.linspace(0, 1, 21); rows = show('polished, 21 times :', th)
print('   R2 vs theta:', ' '.join(f'{r["R2"]:.1e}' for r in rows))
print('   Rs vs theta:', ' '.join(f'{r["Rs"]:.1e}' for r in rows))
