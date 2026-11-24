"""Float check: numerical abscissa of L on {chi psi1*, chi psi2*, chi psi3*}^perp (cut-off adjoint constraints,
as a certification would use) versus the exact spectral complement."""
import numpy as np, scipy.linalg as sl
import slowman as S
x = S.xg; sig = np.abs(x) + S.__dict__.get('STIP', -2.249)  # sigma = |x| + sigma_tip
from evans import STIP
sig = np.abs(x) + STIP
for lo, hi in [(29, 32), (19, 22), (12, 15)]:
    uu = np.clip((hi - sig) / (hi - lo), 0, 1); chi = uu**4 * (35 - 84 * uu + 70 * uu**2 - 20 * uu**3)
    chie = S.B.T @ chi                     # not used directly; cut in function space below
    F = []
    for i in range(3):
        fun = (S.B @ S.Lf[:, i]) / S.r     # adjoint eigenfunction on the full grid
        F.append(S.B.T @ (S.r * chi * fun))
    F = np.column_stack(F)
    Q, _ = np.linalg.qr(F, mode='complete'); Z = Q[:, 3:]
    w = sl.eigh(Z.T @ S.He @ Z, eigvals_only=True, subset_by_index=[Z.shape[1] - 1, Z.shape[1] - 1])[0]
    print(f"cutoff on sigma in [{lo},{hi}]: abscissa on the complement = {w:+.4f}")
print(f"exact spectral complement: {S.ws:+.4f}")
