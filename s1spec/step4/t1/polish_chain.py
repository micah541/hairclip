"""polish a chain of consecutive slabs (same mesh) in arb; the first initial value is the float base, then exact continuity.
usage: python3 polish_chain.py TAG I0 NCHAIN [NTHETA]   -> prints residual norms per slab (max over theta samples and nodes)"""
import sys, time, json, numpy as np
from base2 import Base; from arbops import *; import resid

def polish_chain(B, i0, nchain, tol=1e-26, verbose=True):
    SPs = []; prev = None
    for i in range(i0, i0 + nchain):
        R, D = B.rep(i); RD = B.RD; q = B.slabs[i]
        SP = SlabPolish(R, D, RD, q['nodes'], q['dt'])
        if prev is not None:
            if len(B.slabs[i - 1]['C']) == len(q['C']): SP.Yn[0] = list(prev.Yn[-1])
            else:
                Y = list(prev.Yn[-1]); Rp = prev.R
                while Rp.N < R.N:                       # successive single insertions
                    from rep import Mesh, Rep
                    Cn = q['C']; Cp = Rp.mesh.C
                    newC = [c for c in Cn if not np.any(np.isclose(Cp, c, rtol=0, atol=0))]
                    m1 = Mesh(JT=Rp.mesh.JT, XT=Rp.mesh.XT, h0=Rp.mesh.h0, cT=Rp.mesh.cT, C=np.sort(np.r_[Cp, newC[0]]))
                    R1 = Rep(Rp.P, m1); Y = arb_refine(Y, Rp, R1, SP.AB); Rp = R1
                SP.Yn[0] = Y
        J = D.jac(q['nodes'][len(q['nodes']) // 2]); RD.factor(J, q['dt'])
        t0 = time.time()
        for it in range(60):
            em, dz = SP.step(J)
            if em < tol or dz < 1e-37: break
        if verbose: print(f'  slab {i}: polished in {it+1} it, max|E| {em:.1e} ({time.time()-t0:.0f}s)', flush=True)
        SPs.append(SP); prev = SP
    return SPs

KEYS = ['R2', 'R1', 'Rs', 'RIs', 'RF0', 'RF1', 'RF2', 'RF3', 'RF4']
def slab_norms(SP, nth=11):
    th = sorted(set(list(np.linspace(0, 1, nth)) + [0.0] + list(SP.RD.c)))
    rows = [resid.norms(resid.samples(SP, t))[0] for t in th]
    return {k: max(r[k] for r in rows) for k in KEYS}, rows, th

if __name__ == '__main__':
    tag, i0, nc = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]); nth = int(sys.argv[4]) if len(sys.argv) > 4 else 11
    B = Base(tag); SPs = polish_chain(B, i0, nc)
    for j, SP in enumerate(SPs):
        mx, rows, th = slab_norms(SP, nth)
        print(f'slab {i0+j} k={B.k(B.tau0[i0+j]):.1f}: ' + ' '.join(f'{k}={mx[k]:.1e}' for k in KEYS), flush=True)
    json.dump(dict(theta=th, rows=rows, k=B.k(B.tau0[i0 + nc - 1])), open(f'res_{tag}_{i0+nc-1}.json', 'w'))
