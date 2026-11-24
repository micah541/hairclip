import sys, numpy as np
sys.dont_write_bytecode=True
from run2 import *
from run_hairpin import axis_diag
def ncross_general(S, up=2, maxn=4096):
    """number of transversal self-intersections of the (resampled) closed polyline (all, not only on axis)"""
    al, th, Z, clo = S.curve(up)
    # resample down to <= maxn pts uniformly in arclength (alpha is arclength)
    step = max(1, len(Z)//maxn); Z = Z[::step]; n = len(Z)
    A = Z; B = np.roll(Z, -1)
    cnt = 0; pts=[]
    bb = 0
    for i in range(0, n, 256):
        a = A[i:i+256, None]; b = B[i:i+256, None]; c = A[None, :]; d = B[None, :]
        def cr(u, v): return (np.conj(u) * v).imag
        d1 = cr(b - a, c - a); d2 = cr(b - a, d - a); d3 = cr(d - c, a - c); d4 = cr(d - c, b - c)
        m = (d1 * d2 < 0) & (d3 * d4 < 0)
        ii, jj = np.nonzero(m)
        for p, q in zip(ii + i, jj):
            if q > p + 1 and not (p == 0 and q == n - 1):
                cnt += 1
                t = d3[p - i, q] / (d3[p - i, q] - d4[p - i, q]); pts.append(A[p] + t * (B[p] - A[p]))
    return cnt, pts
if __name__ == '__main__':
    for t in sys.argv[1:]:
        S = load_state(f'{OUT}/state_t{t}.npz'); d = axis_diag(S); c, pts = ncross_general(S)
        print(t, 'axis ncross', len(d['cross']), 'general', c, [complex(round(p.real,3), round(p.imag,3)) for p in pts][:12])
