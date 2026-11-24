"""Float replica of the final Gram certificate for the theta-modified form
   int rho((1-theta) v''^2 - A v'^2 - (V + omega) v^2)  (coefficients A, V - lam divided by 1 - theta),
   with the same constraints f1, f2, f3. Prints the Schur complement and the eigenvalues of G3, plus the
   number of zeros of W01 on the plane run (Morse count)."""
import numpy as np, json, pickle, sys
from scipy.integrate import solve_ivp
import gram_float as GF
import gen_hs3 as G3
STIP = GF.STIP; T0 = G3.T0; c_ = (8 - STIP) / T0; LB = G3.LB; SIG0 = GF.SIG0
GF.LAMT = G3.LAMT
def bumps(X):
    s = X / (c_ * LB)
    if s >= 1: return 0.0, 0.0
    z = (1 - s * s) ** 4
    return LB * s * z, LB**3 * s ** 3 * z / 6
NG = 6
def tail_plane(th, om):
    c = GF.coef(SIG0, -om); c['A'] /= (1 - th); c['VML'] /= (1 - th)
    M = GF.Mfwd(c); ev, V = np.linalg.eig(M)
    i = np.argsort(ev.real)[:2]                 # two most negative: decaying modes in x (admissible)
    v = V[:, i[0]]; v = v / v[0]
    return np.column_stack([v.real, v.imag]), ev
def run(th, om, w0=None):
    Tp, ev = tail_plane(th, om); za0, zb0 = Tp[:, 0].copy(), Tp[:, 1].copy()
    if w0 is None: w0 = np.array(G3.wstar_init())
    def rhs(tau, Y):
        sg = SIG0 - tau; c = GF.coef(sg, -om); c['A'] /= (1 - th); c['VML'] /= (1 - th); M = GF.Mfwd(c)
        C, S = Y[0], Y[1]; ws = Y[2:6]; za, zb = Y[6:10], Y[10:14]
        Ia, Ib = Y[14:14 + NG], Y[14 + NG:14 + 2 * NG]
        yy = Y[14 + 2 * NG:14 + 6 * NG].reshape(NG, 4)
        u = min(max((SIG0 - sg) / 3.0, 0.0), 1.0); chi = u**4 * (35 - 84*u + 70*u*u - 20*u**3)
        e1, e3 = bumps(sg - STIP)
        fv = np.array([chi * c['p'], chi * C, chi * S, chi * ws[0], e1, e3]); gv = fv * c['R'] / (1 - th)
        out = np.empty_like(Y)
        out[0], out[1] = -S * c['k'], C * c['k']
        out[2:6] = -(GF.MLT(c) @ ws) - c['P1'] * ws
        out[6:10] = -(M @ za); out[10:14] = -(M @ zb)
        out[14:14 + NG] = za[0] * gv; out[14 + NG:14 + 2 * NG] = zb[0] * gv
        dyy = -(yy @ M.T); dyy[:, 3] += fv
        out[14 + 2 * NG:14 + 6 * NG] = dyy.ravel()
        out[14 + 6 * NG:] = np.outer(yy[:, 0], gv).ravel()
        return out
    Y0 = np.zeros(14 + 6 * NG + NG * NG); Y0[0] = 1; Y0[2:6] = w0; Y0[6:10] = za0; Y0[10:14] = zb0
    tb = SIG0 - (STIP + c_ * LB)
    # Morse count: zeros of W01 = za0 zb1 - za1 zb0 along the run
    nz = [0]; last = [None]
    def ev_fun(tau, Y): return Y[6] * Y[11] - Y[7] * Y[10]
    ev_fun.terminal = False
    s1 = solve_ivp(rhs, [0, tb], Y0, method='DOP853', rtol=1e-11, atol=1e-14, events=ev_fun)
    s2 = solve_ivp(rhs, [tb, SIG0 - STIP], s1.y[:, -1], method='DOP853', rtol=1e-11, atol=1e-14, events=ev_fun)
    nzero = len(s1.t_events[0]) + len(s2.t_events[0])
    Y = s2.y[:, -1]
    C, S = Y[0], Y[1]; ws = Y[2:6]; za, zb = Y[6:10], Y[10:14]; Ia, Ib = Y[14:14 + NG], Y[14 + NG:14 + 2 * NG]
    yy = Y[14 + 2 * NG:14 + 6 * NG].reshape(NG, 4); II = Y[14 + 6 * NG:].reshape(NG, NG)
    W01 = za[0] * zb[1] - za[1] * zb[0]; W13 = za[1] * zb[3] - za[3] * zb[1]
    G6 = np.zeros((NG, NG))
    for j in range(NG):
        al, be = np.linalg.solve(np.array([[za[1], zb[1]], [za[3], zb[3]]]), -yy[j, [1, 3]])
        for i in range(NG): G6[i, j] = II[j, i] + al * Ia[i] + be * Ib[i]
    a = c_ * ws[1]; b = c_**3 * (ws[3] + ws[1]) + 24 / LB**2 * c_ * ws[1]
    Rm = np.array([[1, 0, 0, 0, 0, 0], [0, C, S, 0, 0, 0], [0, 0, 0, 1, -a, -b]])
    G3f = Rm @ ((G6 + G6.T) / 2) @ Rm.T
    g2, g = G3f[:2, :2], G3f[:2, 2]
    Sch = G3f[2, 2] - g @ np.linalg.solve(g2, g)
    nA = nzero + (1 if W01 * W13 < 0 else 0)
    return Sch, np.linalg.eigvalsh(G3f), nA, np.linalg.eigvalsh(g2)
if __name__ == "__main__":
    for th, om in [(0, 3/64), (0.05, 3/64), (0.05, 1/32), (0.08, 1/32), (0.1, 1/32), (0.1, 1/64), (0.12, 1/64), (0.15, 1/64)]:
        Sch, eg, nA, eg2 = run(th, om)
        print(f"theta={th:5.3f} omega={om:.5f}: Schur={Sch:+.4e}  eig(G3)={np.array2string(eg, precision=4)}  n(A_even)={nA}  eig(g2)={np.array2string(eg2, precision=4)}", flush=True)
