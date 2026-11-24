"""Independent float64 replica of the final certificate (scipy; z_a, z_b integrated directly, BCs solved directly)."""
import numpy as np, json, pickle
from scipy.integrate import solve_ivp
import gram_float as GF
from gram_float import coef, Mfwd, MLT, SIG0
import gen_hs3 as G3
STIP = GF.STIP; T0 = G3.T0; c_ = (8 - STIP) / T0; LB = G3.LB
GF.LAMT = G3.LAMT
T = pickle.load(open('tail_mid.pkl', 'rb'))['Tm']
za0, zb0 = T[:, 2].copy(), T[:, 3].copy(); w0 = np.array(G3.wstar_init())
def bumps(X):
    s = X / (c_ * LB)
    if s >= 1: return 0.0, 0.0
    z = (1 - s * s) ** 4
    return LB * s * z, LB**3 * s ** 3 * z / 6
NG = 6
def rhs(tau, Y):
    sg = SIG0 - tau; c = coef(sg, -GF.OMEGA); M = Mfwd(c)
    C, S = Y[0], Y[1]; ws = Y[2:6]; za, zb = Y[6:10], Y[10:14]
    Ia, Ib = Y[14:14 + NG], Y[14 + NG:14 + 2 * NG]
    yy = Y[14 + 2 * NG:14 + 6 * NG].reshape(NG, 4); II = Y[14 + 6 * NG:].reshape(NG, NG)
    u = min(max((SIG0 - sg) / 3.0, 0.0), 1.0); chi = u**4 * (35 - 84*u + 70*u*u - 20*u**3)
    e1, e3 = bumps(sg - STIP)
    fv = np.array([chi * c['p'], chi * C, chi * S, chi * ws[0], e1, e3]); gv = fv * c['R']
    out = np.empty_like(Y)
    out[0], out[1] = -S * c['k'], C * c['k']
    out[2:6] = -(MLT(c) @ ws) - c['P1'] * ws
    out[6:10] = -(M @ za); out[10:14] = -(M @ zb)
    out[14:14 + NG] = za[0] * gv; out[14 + NG:14 + 2 * NG] = zb[0] * gv
    dyy = -(yy @ M.T); dyy[:, 3] += fv
    out[14 + 2 * NG:14 + 6 * NG] = dyy.ravel()
    out[14 + 6 * NG:] = np.outer(yy[:, 0], gv).ravel()
    return out
Y0 = np.zeros(14 + 6 * NG + NG * NG); Y0[0] = 1; Y0[2:6] = w0; Y0[6:10] = za0; Y0[10:14] = zb0
# integrate in two pieces so the bump support start is a mesh point
tb = SIG0 - (STIP + c_ * LB)
s1 = solve_ivp(rhs, [0, tb], Y0, method='DOP853', rtol=1e-12, atol=1e-15)
s2 = solve_ivp(rhs, [tb, SIG0 - STIP], s1.y[:, -1], method='DOP853', rtol=1e-12, atol=1e-15)
Y = s2.y[:, -1]
C, S = Y[0], Y[1]; ws = Y[2:6]; za, zb = Y[6:10], Y[10:14]; Ia, Ib = Y[14:14 + NG], Y[14 + NG:14 + 2 * NG]
yy = Y[14 + 2 * NG:14 + 6 * NG].reshape(NG, 4); II = Y[14 + 6 * NG:].reshape(NG, NG)
G6 = np.zeros((NG, NG))
for j in range(NG):
    al, be = np.linalg.solve(np.array([[za[1], zb[1]], [za[3], zb[3]]]), -yy[j, [1, 3]])
    for i in range(NG): G6[i, j] = II[j, i] + al * Ia[i] + be * Ib[i]
a = c_ * ws[1]; b = c_**3 * (ws[3] + ws[1]) + 24 / LB**2 * c_ * ws[1]
Rm = np.array([[1, 0, 0, 0, 0, 0], [0, C, S, 0, 0, 0], [0, 0, 0, 1, -a, -b]])
G3f = Rm @ ((G6 + G6.T) / 2) @ Rm.T
g2, g = G3f[:2, :2], G3f[:2, 2]
Sch = G3f[2, 2] - g @ np.linalg.solve(g2, g)
cert = json.load(open('certificate_step3_linear.json'))
np.set_printoptions(precision=8, linewidth=160)
print("float: C_tip, S_tip =", C, S, "  a, b =", a, b, "  asym(G6) =", np.abs(G6 - G6.T).max())
print("float G3:\n", G3f)
print("float Schur:", Sch, "   rigorous enclosure:", cert["S"], " contained:", cert["S"][0] <= Sch <= cert["S"][1])
d = np.abs(G3f - np.array(cert["G3_mid"])); r = np.array(cert["G3_rad"])
print("max |float - rigorous mid| / rigorous rad:", (d / r).max())
