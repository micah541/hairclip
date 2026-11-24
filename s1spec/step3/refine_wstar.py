"""Float Newton for (alpha, lamt): adjoint solution w~ from initial data alpha*wa - wb at sigma=32 with
w~'(tip) = ws1(tip) = 0 and ws3(tip) = 0 (then F'(0) = F'''(0) = 0)."""
import numpy as np, json
from scipy.integrate import solve_ivp
import gram_float as GF
from gram_float import coef, SIG0
STIP = GF.STIP
def lt_frame(lamt):
    roots = np.roots([1, 0, 0, -SIG0 / 4, -(0.5 - lamt)])
    z = [r for r in roots if r.real < -0.5][0]; v = np.array([1, z, z * z, z ** 3]); return v.real, v.imag
def tip(alpha, lamt):
    wa, wb = lt_frame(lamt); w0 = alpha * wa - wb
    def f(tau, Y):
        c = coef(SIG0 - tau, -GF.OMEGA); k, m, q, P1 = c['k'], c['m'], c['q'], c['P1']
        M = np.array([[0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1], [(0.5 - lamt), (k * m + q / 4), -k * k, 0]])
        return -(M @ Y) - P1 * Y
    s = solve_ivp(f, [0, SIG0 - STIP], w0, method='DOP853', rtol=1e-13, atol=1e-18)
    return s.y[:, -1]
al, lt = json.load(open('wstar_coef.json'))[0], -0.25626
for it in range(6):
    y = tip(al, lt); r = np.array([y[1], y[3]])
    h1, h2 = 1e-7, 1e-7
    J = np.column_stack([(np.array([tip(al + h1, lt)[1], tip(al + h1, lt)[3]]) - r) / h1,
                         (np.array([tip(al, lt + h2)[1], tip(al, lt + h2)[3]]) - r) / h2])
    d = np.linalg.solve(J, -r); al, lt = al + d[0], lt + d[1]
    print(f"it {it}: alpha = {al:.16f}  lamt = {lt:.16f}   |res| before = {np.abs(r).max():.2e}")
y = tip(al, lt); print("final tip ws:", y)
json.dump({"alpha": al, "lamt": lt}, open('wstar_refined.json', 'w'))
