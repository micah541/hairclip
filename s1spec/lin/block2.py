import numpy as np
from block import mu, lognorm
G0 = np.array([np.exp(1j*np.pi/3), np.exp(-1j*np.pi/3)])
def scan(c, S=8.0, sig_max=4000.0, n=1500, w=1.0):
    sigs = np.exp(np.linspace(np.log(sig_max), np.log(S), n))   # decreasing: eps grows from ~0
    z = G0.copy(); rows = []
    for s_ in sigs:
        m = mu(s_); eps = c/m**4
        for _ in range(60):
            z = z - (z**4 + z - eps)/(4*z**3 + 1)
        z2, z3 = z
        sp_ = z2 + z3; pp = z2*z3; s = -sp_; p = sp_**2 - pp
        dz = 1/(4*z**3 + 1); dsp = dz[0] + dz[1]; dpp = dz[0]*z3 + z2*dz[1]
        ds = -dsp; dp = 2*sp_*dsp - dpp
        v = lambda zz: np.array([1, zz, zz*zz, zz**3]); dv = lambda zz, d: np.array([0, d, 2*zz*d, 3*zz*zz*d])
        ep = np.array([1, s/2, (s*s - 2*p)/2, (s**3 - 3*s*p)/2]); em = np.array([0, 1, s, s*s - p])
        dep = np.array([0, ds/2, (2*s*ds - 2*dp)/2, (3*s*s*ds - 3*ds*p - 3*s*dp)/2]); dem = np.array([0, 0, ds, 2*s*ds - dp])
        U = np.column_stack([ep, em, v(z2), v(z3)]); dU = np.column_stack([dep, dem, dv(z2, dz[0]), dv(z3, dz[1])])
        Ui = np.linalg.inv(U)
        Mt = Ui @ np.diag([0, 1, 2, 3]) @ U - 4*eps*Ui @ dU
        Bs = np.array([[s/2, 1], [(s*s - 4*p)/4, s/2]])
        Bm = -Mt/(3*s_); Bm[:2, :2] += m*Bs; Bm[2, 2] += m*z2; Bm[3, 3] += m*z3
        Wm = np.diag([1, w, 1, 1]); Bm = np.linalg.inv(Wm) @ Bm @ Wm
        lS = lognorm(Bm[:2, :2]); mG = -lognorm(-Bm[2:, 2:])
        rows.append((s_, lS, mG, mG - lS, np.linalg.norm(Bm[2:, :2], 2), np.linalg.norm(Bm[:2, 2:], 2), np.linalg.cond(U)))
    return np.array(rows)
def eta_fp(R):
    a, b, g = R[:, 4], R[:, 5], R[:, 3]
    if np.any(g <= 0): return np.inf, np.inf
    eta = 0.0
    for _ in range(200):
        eta_n = np.max((a + eta**2*b)/g)
        if eta_n > 10: return np.inf, np.inf
        eta = eta_n
    return eta, np.max(2*eta*b/g)
if __name__ == "__main__":
    for lam in [1.0, 0.25, -0.125, 2.0, 1.8+0.4j, -0.125+4j, 1+2j, 0.5+4.5j, -0.125+0.0001j]:
        c = 0.25 - lam
        out = []
        for w in [0.5, 1.0, 2.0, 4.0]:
            R = scan(c, w=w); e, L = eta_fp(R)
            i = np.argmin(R[:, 3])
            out.append(f"w={w}: eta={e:.4f} Lip={L:.3f} mingap={R[i,3]:.3f}@{R[i,0]:.1f}")
        R = scan(c, w=1.0)
        print(f"lam={lam}: " + " | ".join(out) + f" | cond(U) max {R[:,6].max():.1f}", flush=True)
