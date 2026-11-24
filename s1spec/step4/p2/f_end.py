"""float: far-end point of the base, Re Gamma(ell)/ell, to compare the far-point angle with the turning defect"""
import sys, numpy as np
sys.path.insert(0, '/home/claude/s1spec/step4/t1')
import polstore
from rep import DS1
for s in map(int, sys.argv[1].split(',')):
    SP, d = polstore.load_slab(f'/home/claude/s1spec/step4/t1/polished/slab_{s:03d}.pkl'); R = SP.R; n0 = R.n
    y = np.array([float(v.mid()) for v in SP.Yn[-1]]); ell = y[n0]
    xk, _ = R.mesh.knots(ell); g, wq = np.polynomial.legendre.leggauss(40)
    xs = np.concatenate([0.5 * (a + b) + 0.5 * (b - a) * g for a, b in zip(xk[:-1], xk[1:])])
    ws = np.concatenate([0.5 * (b - a) * wq for a, b in zip(xk[:-1], xk[1:])])
    k = R.eval_at(y[:n0], ell, xs, 0)
    # theta at quadrature points: integrate kappa element-wise exactly enough (Gauss on sub-intervals)
    th = np.empty_like(xs); acc = 0.0; i = 0
    for a, b in zip(xk[:-1], xk[1:]):
        xx = 0.5 * (a + b) + 0.5 * (b - a) * g
        for j, xj in enumerate(xx):
            sub = 0.5 * (a + xj) + 0.5 * (xj - a) * g
            th[i + j] = -np.pi / 2 + acc + (0.5 * (xj - a) * wq * R.eval_at(y[:n0], ell, sub, 0)).sum()
        acc += (0.5 * (b - a) * wq * R.eval_at(y[:n0], ell, xx, 0)).sum(); i += len(g)
    Z = -DS1 + (ws * np.exp(1j * th)).sum()
    print(f'slab {s}: ell = {ell:.6f}, Gamma(ell) = {Z.real:.6f} + {Z.imag:.3e} i, Re/ell = {Z.real/ell:.5f}, turning-pi = {acc-np.pi:.3e}')
