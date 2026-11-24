"""Frozen-root dichotomy along the tail for the linear eigen-ODE (model part):
   nu^4 + (sigma/4) nu - c = 0,  c = 1/4 - lam,   nu = mu*zeta,  zeta^4 + zeta = eps = c/mu^4.
Track roots continuously from eps=0 (0, -1, e^{+-i pi/3}) to eps(S); report min gap
   min_G Re zeta_G - max_S Re zeta_S   (relative to mu)."""
import numpy as np
def roots_path(eps_end, n=400):
    z = np.array([0, -1, np.exp(1j*np.pi/3), np.exp(-1j*np.pi/3)], dtype=complex)
    gaps = []
    for t in np.linspace(0, 1, n)[1:]:
        e = t*eps_end
        r = np.roots([1, 0, 0, 1, -e])
        # match continuously
        new = np.empty(4, dtype=complex); used = set()
        for i in range(4):
            j = min((j for j in range(4) if j not in used), key=lambda j: abs(r[j]-z[i]))
            new[i] = r[j]; used.add(j)
        z = new
        gaps.append(min(z[2].real, z[3].real) - max(z[0].real, z[1].real))
    return z, min(gaps)
mu = lambda s: (s/4)**(1/3)
for S in [8, 20, 30, 40, 60]:
    m4 = mu(S)**4
    out = []
    for lam in [2.7, -0.125, -0.125+5j, 2.7+5j, -0.125+9j, 2.7+9j, -0.125+15j, 2.7+15j, 1+12j]:
        c = 0.25 - lam
        z, g = roots_path(c/m4)
        out.append(f"{lam}: gap*mu={g*mu(S):.3f}")
    print(f"S={S} (mu={mu(S):.3f}, mu^4={m4:.2f}):", "; ".join(out))
