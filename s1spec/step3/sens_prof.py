import numpy as np, sys
from scipy.integrate import solve_ivp
import gram_float as GF
import gen_hs2
from gen_evans import profile_part
tau1, xlo, xhi, Ct, r0lo, r0hi = profile_part()
xm = np.array([(float(a) + float(b)) / 2 for a, b in zip(xlo, xhi)])
Cm = np.array([[float(Ct[i][j]) for j in range(4)] for i in range(4)])
r0 = np.array([float(b) for b in r0hi])
f = lambda s, y: [-y[2]*y[1], 1 + y[2]*y[0], y[3], y[0]/4]
orig_prof, orig_stip = GF.prof, GF.STIP
def make(y8):
    ev = lambda s, y: y[1]
    ev.direction = 1
    sol = solve_ivp(f, [8.0, -3.0], y8, method='DOP853', rtol=1e-13, atol=1e-16, dense_output=True, events=ev)
    stip = sol.t_events[0][0]
    def prof(s, freeze_tail=True):
        if s <= 8.0: return sol.sol(s)
        return orig_prof(s)
    return prof, stip
def G3for(prof, stip):
    GF.prof = prof; GF.STIP = stip
    sol, up = GF.run(sig_end=stip); G4, G3, tip = GF.gram(sol, up)
    GF.prof, GF.STIP = orig_prof, orig_stip
    Gs = (G3 + G3.T) / 2
    g2 = Gs[:2, :2]; g = Gs[:2, 2]
    return Gs, Gs[2, 2] - g @ np.linalg.solve(g2, g)
p0, s0 = make(xm); G0, S0 = G3for(p0, s0)
print("base stip", s0, "S", S0, "G33", G0[2, 2])
for k in range(4):
    for sg in [+1, -1]:
        y8 = xm + sg * Cm[:, k] * r0[k]
        pk, sk = make(y8); Gk, Sk = G3for(pk, sk)
        print(f"dir {k} {sg:+d}: stip shift {sk - s0:+.2e}  dS {Sk - S0:+.2e}  dG33 {Gk[2,2]-G0[2,2]:+.2e}  max|dG| {np.abs(Gk-G0).max():.2e}")
