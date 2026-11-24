import sys, time, numpy as np
from rep import *
from seed import HLSHalf

def project(R, F, ell):
    """nodal DOF from exact jets, bubbles by L2 projection of the remainder per element"""
    x, _ = R.mesh.knots(ell); h = np.diff(x); N = R.N; B = R.B
    D = np.array([F(x, m) for m in range(5)]).T
    y = R.to_y(D, np.zeros((N, B.nbub)), h)
    if B.nbub:
        d = R.local(y, h)                         # bubbles zero so far
        xq = x[:-1, None] + h[:, None] * B.uq[None, :]
        fq = F(xq.ravel(), 0).reshape(xq.shape)
        rem = fq - d @ B.V[0].T
        Vb = B.V[0][:, 10:]; G = Vb.T @ (B.wq[:, None] * Vb)
        bb = np.linalg.solve(G, (Vb.T * B.wq[None, :]) @ rem.T).T
        for e in range(N): y[R.bidx[e]] = bb[e]            # already in normalised units
    return y

if __name__ == '__main__':
    S = HLSHalf(sys.argv[1]); print('k', S.kmax, 'lam', S.lam, 'ell', S.ell, 'N', S.N)
    for P in (9, 11, 13):
        for (JT, h0, cT, NO, c, beta) in [(40, 0.35, 0.08, 80, 0.12, 0.02), (60, 0.25, 0.06, 100, 0.10, 0.015)]:
            R = Rep(P, Mesh(JT=JT, h0=h0, cT=cT, NO=NO, c=c, beta=beta))
            y = project(R, S.kap, S.ell)
            xs = np.r_[np.linspace(0, 40, 8001), np.linspace(40, S.ell, 20001)[1:]]
            e = np.abs(R.eval_at(y, S.ell, xs) - S.kap(xs))
            e2 = np.abs(R.eval_at(y, S.ell, xs, 2) - S.kap(xs, 2))
            x, _ = R.mesh.knots(S.ell)
            print(f"P={P} JT={JT} NO={NO} n={R.n}: max err k: tip[0,10] {e[xs<=10].max():.2e} [10,40] {e[(xs>10)&(xs<=40)].max():.2e} out {e[xs>40].max():.2e} | k'' tip {e2[xs<=10].max():.2e} out {e2[xs>40].max():.2e} | hmax_out {np.diff(x)[JT:].max():.2f} hmin_out {np.diff(x)[JT:].min():.2f}")
