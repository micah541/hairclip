"""Numerical (non-rigorous): the slow-mode structure of the modulated rescaled flow at S1, even class.
Spectral decomposition u = c psi3 + v, v in the spectral complement of {psi1, psi2, psi3}, spectral modulation
<u, psi_i*>_rho = 0 (i = 1, 2).  Then, exactly,
   c' = lam3 c + <N(u), psi3*>,   v' = L v + P N(u),   a = <N(u), psi1*>,  b = <N(u), psi2*>   (quadratic!)
Second order: c' = lam3 c + g2 c^2,  v ~ c^2 h2 with (L - 2 lam3) h2 = -P N2(psi3, psi3)."""
import sys, numpy as np, scipy.linalg as sl
sys.path.insert(0, '/home/claude/s1spec/step3')
from polyweight import build
from nlflow import Flow
X, h = 34.0, 0.05
xg, B, A, r, ps1, ps2 = build(0.25, X, h, kind='gauss')
n = int(round(X / h))
Ae = B.T @ A @ B                               # even part, r-scaled coordinates (Euclidean = L^2(rho))
ev, vl, vr = sl.eig(Ae, left=True, right=True)
lams = [1.0, 0.25, -0.2563]
idx = [int(np.argmin(abs(ev - l))) for l in lams]
R = np.real(vr[:, idx]); Lf = np.real(vl[:, idx])
Lf = Lf @ np.linalg.inv(R.T @ Lf).T            # biorthonormal: Lf^T R = I
print("eigenvalues:", np.round(ev[idx].real, 5), " biorth check:", np.round(Lf.T @ R, 8).tolist())
# half-line representation helpers: even vector e (len n+1, r-scaled) -> function on x>=0
rh = r[n:]
def to_fun(e): return (B @ e)[n:] / rh
def from_fun(f):                                 # f on x >= 0 (even) -> even r-scaled vector
    full = np.concatenate([f[::-1][:-1], f]) * r
    return B.T @ full
psi = [to_fun(R[:, i]) for i in range(3)]; s3 = psi[2][0]
psi = [p / (p[0] if i == 2 else 1.0) for i, p in enumerate(psi)]
R[:, 2] /= s3; Lf[:, 2] *= s3                   # normalise psi3(0) = 1, keep biorthonormality
# spectral-complement numerical abscissa
He = (Ae + Ae.T) / 2
Q, _ = np.linalg.qr(Lf, mode='complete'); Z = Q[:, 3:]
ws = sl.eigh(Z.T @ He @ Z, eigvals_only=True, subset_by_index=[Z.shape[1] - 1, Z.shape[1] - 1])[0]
Q2, _ = np.linalg.qr(Lf[:, :2], mode='complete'); Z2 = Q2[:, 2:]
w2 = sl.eigh(Z2.T @ He @ Z2, eigvals_only=True, subset_by_index=[Z2.shape[1] - 1, Z2.shape[1] - 1])[0]
print(f"numerical abscissa on the spectral complement of psi1,psi2,psi3: {ws:+.4f}  (of psi1,psi2 only: {w2:+.4f})")
# second variation of the rescaled CDF right-hand side along psi3
fl = Flow(X=X, h=h)
def R0(u): return fl.rhs_parts(u)[0]
eps = 1e-3
N2 = (R0(eps * psi[2]) + R0(-eps * psi[2])) / (2 * eps**2)
N3 = (R0(2 * eps * psi[2]) - 2 * R0(eps * psi[2]) + 2 * R0(-eps * psi[2]) - R0(-2 * eps * psi[2])) / (12 * eps**3)
# the far end of the finite domain is not meaningful for psi3 (growth |x|^2): cut the forcing beyond x = 25
cut = np.where(fl.x < 22, 1.0, np.where(fl.x < 25, np.cos(np.pi / 6 * (fl.x - 22))**2, 0.0))
N2e = from_fun(N2 * cut); N3e = from_fun(N3 * cut)
g2 = Lf[:, 2] @ N2e; a2 = Lf[:, 0] @ N2e; b2 = Lf[:, 1] @ N2e
g3 = Lf[:, 2] @ N3e
PN2 = N2e - R @ (Lf.T @ N2e)
h2e = np.linalg.solve(Ae - 2 * lams[2] * np.eye(len(Ae)), -PN2)
h2 = to_fun(h2e)
nr = lambda e: np.sqrt(e @ e)
print(f"c' = lam3 c + g2 c^2 + g3 c^3 + ...:  g2 = {g2:+.4f}, g3 = {g3:+.4f};   modulation a = {a2:+.4f} c^2, b = {b2:+.4f} c^2")
print(f"||psi3||_rho = {nr(R[:,2]):.3f},  ||P N2||_rho = {nr(PN2):.3f},  ||h2||_rho = {nr(h2e):.3f}")
for xx in [0, 1, 2, 3, 4, 6, 8, 10]:
    j = int(round(xx / h)); print(f"   x={xx:4.1f}: psi3={psi[2][j]:+9.3f}  h2={h2[j]:+10.3f}  N2={N2[j]:+10.3f}")
np.savez("slowman.npz", x=fl.x, psi3=psi[2], h2=h2, N2=N2, psi1=psi[0], psi2=psi[1],
         adj=np.array([to_fun(Lf[:, i]) * 1.0 for i in range(3)]), g2=g2, g3=g3, a2=a2, b2=b2, ws=ws)
