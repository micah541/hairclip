import numpy as np, scipy.linalg as sl
from compress import build
from proto import cutoff
beta, X, h = 0.25, 30.0, 0.04
x, B, A, L, r, psi1, psi2 = build(beta, X, h)
Ae = B.T @ A @ B; He = (Ae + Ae.T)/2
ev, vl, vr = sl.eig(Ae, left=True, right=True)
i3 = int(np.argmin(abs(ev + 0.2563)))
lam3 = ev[i3].real
psi3s = np.real(vl[:, i3])                 # conjugated adjoint eigenvector (= r * psi3^*)
chi = cutoff(x, 30.0)
F = np.column_stack([B.T @ (r*chi*psi1), B.T @ (r*chi*psi2), psi3s])
Qm, _ = np.linalg.qr(F, mode='complete'); Z = Qm[:, 3:]
top = sl.eigh(Z.T @ He @ Z, eigvals_only=True, subset_by_index=[Z.shape[1]-2, Z.shape[1]-1])[::-1]
print("lam3 =", lam3, " compressed top on perp(p, -sin th, psi3*):", np.round(top, 4))
for omega in [0.03, 0.05, 0.07]:
    Aop = -He - omega*np.eye(He.shape[0])
    nA = int(np.sum(np.linalg.eigvalsh(Aop) < 0))
    G = F.T @ np.linalg.solve(Aop, F); G = (G + G.T)/2
    D = np.sqrt(abs(np.diag(G)))
    print(f"omega={omega}: n(A)={nA} eig(G)={np.round(np.linalg.eigvalsh(G),4)} normalized={np.round(np.linalg.eigvalsh(G/np.outer(D,D)),4)}")
# tail behaviour of the adjoint eigenfunction: rho * psi3^* (= r*(conj vector)) along x
n = (len(x)-1)//2
v = (B @ psi3s)[n:]; xs = x[n:]
for xx in [0, 2, 5, 10, 15, 20, 25]:
    i = int(round(xx/h)); print(f"x={xx:4.1f}  r*psi3* = {v[i]:+.3e}   psi3* = {v[i]/r[n+i]:+.3e}")
