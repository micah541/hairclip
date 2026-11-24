import numpy as np, subprocess, sys
sys.path.insert(0, '/home/claude/s1spec'); sys.path.insert(0, '/home/claude/s1spec/lin')
from gen_evans import write
import evans as E
from block2 import G0
from block import mu
from scipy.integrate import solve_ivp
def frozen_plane(lam, S=8.0):
    """Gamma = 0 plane at sigma=S in the block basis: columns D_mu e+, D_mu e- (numerical)."""
    c = 0.25 - lam; m = mu(S); eps = c/m**4
    z = G0.copy()
    for t in np.linspace(0, 1, 401)[1:]:
        for _ in range(40): z = z - (z**4 + z - t*eps)/(4*z**3 + 1)
    z2, z3 = z; sp_ = z2 + z3; pp = z2*z3; s = -sp_; p = sp_**2 - pp
    ep = np.array([1, s/2, (s*s - 2*p)/2, (s**3 - 3*s*p)/2]); em = np.array([0, 1, s, s*s - p])
    D = np.diag([1, m, m*m, m**3])
    ya, yb = D @ ep, D @ em
    return np.array([ya[i]*yb[j] - ya[j]*yb[i] for (i, j) in [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]])
def numeric_tip(lam, W8):
    sol = solve_ivp(lambda s, w: E.rhs(s, w, lam) - E.psi_prime(s)*w, [8.0, E.STIP], W8.astype(complex),
                    method='DOP853', rtol=1e-12, atol=1e-16)
    return sol.y[:, -1]
lam = complex(sys.argv[1]) if len(sys.argv) > 1 else 0.3+0.2j
W8 = frozen_plane(lam)
Wn = numeric_tip(lam, W8)
# rescaling factor e^{psi}: psi' = mu(max(s,1)) integrated from 8 down to STIP (W~ = e^{psi}W, psi(8)=0)
from scipy.integrate import quad
psi = -quad(lambda s: E.mu(max(s, 1.0)), E.STIP, 8.0)[0]
Wn_true = Wn*np.exp(-psi)
Wre = []
for w in W8: Wre += [w.real, w.imag]
eps = 1e-15
write("test_in.txt", (lam.real, lam.real, lam.imag, lam.imag),
      [v - eps*abs(v) for v in Wre], [v + eps*abs(v) for v in Wre], [], [], [], [0.0]*12, [0.0]*12)
out = subprocess.run(["./evans_cap", "test_in.txt", "20"], capture_output=True, text=True).stdout
print(out)
print("numeric W13 =", Wn_true[4], "  W02 =", Wn_true[1])
