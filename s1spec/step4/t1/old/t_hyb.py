import numpy as np
from rep import *; from ops import *
def hybrid_mass(R, ell):
    B = R.B; x, _ = R.mesh.knots(ell); h = np.diff(x); sc = R.lscale(h); N = R.N; n = R.n
    nm = B.nmom - 1
    M = np.zeros((n, n)); row = 0
    MV = B.Mom[:nm] @ B.V[0]                     # (nm, nb)
    G0 = (B.V[0].T * B.wq[None, :]) @ B.V[0]      # Gram on [0,1] (nb, nb)
    for e in range(N):
        for m in range(nm):
            for q in range(B.nb):
                g = R.L[e, q]
                if g >= 0: M[row, g] += MV[m, q] * sc[e, q]
            row += 1
    for j in range(1, N + 1):
        # test function = basis function of yhat_{j,4} (local index 4 in element j (left end) and 9 in element j-1 (right end)), scaled
        for e, qt in ((j - 1, 9), (j, 4)):
            if e >= N: continue
            for q in range(B.nb):
                g = R.L[e, q]
                if g >= 0: M[row, g] += h[e] * G0[qt, q] * sc[e, q] * sc[e, qt]
        row += 1
    assert row == n
    return M
for P in (9, 11):
    R = Rep(P, Mesh()); D = Dyn(R)
    for ell in (314.5, 4e4, 3.7e5):
        M = hybrid_mass(R, ell)
        # row equilibration (Galerkin rows have element-length factors)
        Ms = M / np.linalg.norm(M, axis=1)[:, None]
        s = np.linalg.svd(Ms, compute_uv=False); s0 = np.linalg.svd(D.mass(ell).toarray(), compute_uv=False)
        print(f"P={P} ell={ell:9.1f}: hybrid cond {s[0]/s[-1]:.2e} (smallest {s[-3:]})   pure moments cond {s0[0]/s0[-1]:.2e}")
