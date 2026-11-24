"""DIAGNOSIS ONLY (floating point, not part of any certificate).
Linearised propagation of the certified sigma = 8 parallelepiped (../in_br2.txt: x + C r0) to the tip,
to see what drives the prof_tight widths there: the eps bracket (column 0), d = Z_nl,1(8) (columns 1, 2),
the q remainder (column 3), or the integrator (CAPD width / linear width).
Reversed time tau = 8 - sigma, field (p, q, k, m, th)' = (k q, -1 - k p, -m, -p/4, -k)  (as s1prof.cpp)."""
import numpy as np, json, os
from scipy.integrate import solve_ivp
HERE = os.path.dirname(os.path.abspath(__file__)); TIGHT = os.path.dirname(HERE)
v = [float(t) for t in open(os.path.join(TIGHT, "in_br2.txt")).read().split()]
x0 = np.array([(v[1 + 2 * i] + v[2 + 2 * i]) / 2 for i in range(5)])
C = np.array(v[11:36]).reshape(5, 5)
r0lo = np.array([v[36 + 2 * j] for j in range(5)]); r0hi = np.array([v[37 + 2 * j] for j in range(5)])
rad = (r0hi - r0lo) / 2; off = (r0hi + r0lo) / 2          # column-j parameter in off_j +- rad_j
x0 = x0 + C @ off
def f(t, y):
    p, q, k, m, th = y[:5]; V = y[5:].reshape(5, 5)
    J = np.array([[0, k, q, 0, 0], [-k, 0, -p, 0, 0], [0, 0, 0, -1, 0], [-0.25, 0, 0, 0, 0], [0, 0, -1, 0, 0]])
    return np.concatenate([[k * q, -1 - k * p, -m, -p / 4, -k], (J @ V).ravel()])
def qzero(t, y): return y[1]
qzero.direction = 1
y0 = np.concatenate([x0, C.ravel()])
TC = 10.248977652713181                                     # centre of the last prof_tight row
sol = solve_ivp(f, [0, 10.4], y0, method="DOP853", rtol=1e-13, atol=1e-22, dense_output=True,
                events=qzero)
names = ["p", "q", "k", "m", "th"]
yt = sol.sol(TC); V = yt[5:].reshape(5, 5)
lin = np.abs(V) * rad                                        # half-width contribution of each direction
L = open(os.path.join(TIGHT, "prof", "prof_tight.txt")).readlines()[-2].split()
capd = np.array([(float(L[4 + 2 * i]) - float(L[3 + 2 * i])) / 2 for i in range(5)])
out = {"tau_c": TC}
print(f"half-widths at tau = {TC} (last prof_tight row centre):")
print("comp   eps-dir    Re d       Im d       q-rem      linear sum  CAPD       CAPD/linear")
for i, n in enumerate(names):
    print(f"{n:4s} " + " ".join(f"{c:.3e}" for c in lin[i, :4]) + f"  {lin[i].sum():.3e}   {capd[i]:.3e}  {capd[i] / lin[i].sum():.3f}")
out["halfwidth_by_direction"] = {n: lin[i, :4].tolist() for i, n in enumerate(names)}
out["capd_halfwidth"] = capd.tolist()
# the tip: first q = 0 crossing (MinusPlus in tau), M = m there; Poincare-map sensitivities
tt = sol.t_events[0][0]; ye = sol.sol(tt); Ve = ye[5:].reshape(5, 5); fe = f(tt, ye)[:5]
dtau = -Ve[1, :] / fe[1]                                     # d tau_tip / d r_j
dx = Ve + np.outer(fe, dtau)                                 # d (state at crossing) / d r_j
M_r = dx[3, :]
print(f"\ntip: tau_tip = {tt:.13f}, M = m(tip) = {ye[3]:.3e} (float centre), dq/dtau = {fe[1]:.4f}, dm/dtau = {fe[3]:.4f}")
print("d tau_tip / d r_j:", " ".join(f"{a:.3e}" for a in dtau[:4]))
print("dM / d r_j       :", " ".join(f"{a:.3e}" for a in M_r[:4]))
Mw_fixed_eps = np.sum(np.abs(M_r[1:4]) * rad[1:4])
print(f"M half-width at fixed eps (d, q-rem): {Mw_fixed_eps:.3e}  (certified M enclosures: ~1.13e-10 half-width)")
print(f"eps half-width forced by it: {Mw_fixed_eps / abs(M_r[0]):.3e}   (bracket half-width used: {rad[0]:.3e})")
print(f"tau_tip half-width (bracket set, linear): {np.sum(np.abs(dtau[:4]) * rad[:4]):.3e}  (certified: 4.7e-10)")
# conditioned on M = 0: eliminate eps through M(eps, d, qrem) = 0
G = V[:, 1:4] - np.outer(V[:, 0], M_r[1:4] / M_r[0])          # d x(TC) / d (Re d, Im d, qrem) along M = 0
cond = np.sum(np.abs(G) * rad[1:4], axis=1)
print("\nhalf-widths at tau_c conditioned on M = 0 (eps eliminated):")
for i, n in enumerate(names):
    print(f"{n:4s} {cond[i]:.3e}   (vs linear box {lin[i].sum():.3e}, CAPD {capd[i]:.3e})")
out.update(tau_tip=tt, dtau=dtau.tolist(), dM=M_r.tolist(), M_halfwidth_fixed_eps=Mw_fixed_eps,
           cond_halfwidth=cond.tolist(), rad=rad.tolist())
json.dump(out, open(os.path.join(HERE, "diag_float.json"), "w"), indent=1)
