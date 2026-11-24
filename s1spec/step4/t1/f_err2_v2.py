"""[t1/f_err2_v2.py: copy of p2/f_err2.py; only change: store directory from env STORE (default t1/polished_v2)]
float probe 3 (not certified): first-order normal error e of the P1 curve over the base on the WHOLE half curve.
Tip-normalised placement, clock tau* and (p,n)-coefficients chosen to minimise ||chi e||_rho; then A, B, V0 (T5 1.1).
e' = dtheta - kappa w, w' = kappa e (first order), e'' = dkappa - kappa' w - kappa^2 e."""
import sys, json, numpy as np
sys.path.insert(0, '/home/claude/s1spec/step4/t1')
import polstore
from rep import K0, DS1
from scipy.optimize import minimize_scalar
STIP = -2.249077651121112
ck, s = sys.argv[1], int(sys.argv[2])
rec = json.load(open(ck)); u = np.array([int(m) * 2.0 ** int(e) for m, e in rec['xbar']]); n = np.arange(1, len(u) + 1)
L = int(rec['L'][0]) * 2.0 ** int(rec['L'][1]); c = 2 * np.pi / L
def thd(al, m):   # d^m theta / d alpha^m, m >= 1
    ph = np.outer(al, n); f = [np.cos, lambda z: -np.sin(z), lambda z: -np.cos(z), np.sin][(m - 1) % 4]
    out = np.zeros(len(al))
    for i in range(0, len(al), 4000): out[i:i+4000] = (f(ph[i:i+4000]) * (n ** m * u)).sum(1)
    return out + (1.0 if m == 1 else 0.0)
lam = K0 / (-c * thd(np.array([np.pi]), 1)[0]); ellP = L / 2 / lam
import os
SP, d = polstore.load_slab(os.path.join(os.environ.get('STORE', '/home/claude/s1spec/step4/t1/polished_v2'), f'slab_{s:03d}.pkl'))
R = SP.R; nd = np.array(d['nodes_rel']); nd = nd if len(nd) == len(SP.Yn) else np.r_[0.0, nd]; n0 = R.n
Y = np.array([[float(v.mid()) for v in row] for row in SP.Yn])
def lag(tr):
    w = np.ones(len(nd))
    for i in range(len(nd)):
        for j in range(len(nd)):
            if i != j: w[i] *= (tr - nd[j]) / (nd[i] - nd[j])
    return w
ell_b = (lag(0.5) @ Y)[n0]; xe = min(ellP, ell_b) - 1e-9
x = np.r_[np.linspace(0, 40, 16001), np.linspace(40, xe, 40001)[1:]]; dx = np.diff(x)
kP = [(-1) ** (m + 1) * (lam * c) ** (m + 1) * thd(np.pi - c * lam * x, m + 1) for m in range(3)]
def cumtrap(f): return np.r_[0, np.cumsum(0.5 * (f[1:] + f[:-1]) * dx)]
rho = np.where(x <= 34.25, np.exp(-(1 + x ** 2) ** (2 / 3) / 4), 0.0)
uu = np.clip((32 - (x + STIP)) / 3, 0, 1); chi = uu ** 4 * (35 - 84 * uu + 70 * uu ** 2 - 20 * uu ** 3)
zz = np.clip((x - 10.25) / 6, 0, 1); zeta = zz ** 4 * (35 - 84 * zz + 70 * zz ** 2 - 20 * zz ** 3)
W = (1 + x ** 2) ** -2.0
def integ(f): return float(cumtrap(f)[-1])
def state(tr):
    y = lag(tr) @ Y; ell = y[n0]
    kb = [R.eval_at(y[:n0], ell, x, m) for m in range(3)]
    dk = [kP[m] - kb[m] for m in range(3)]
    thb = -np.pi / 2 + cumtrap(kb[0]); G = -DS1 + cumtrap(np.exp(1j * thb)); T = np.exp(1j * thb); nu = 1j * T
    dth = cumtrap(dk[0])
    # e' = dth - k w, w' = k e : trapezoid-implicit steps
    e = np.zeros(len(x)); w = np.zeros(len(x)); k = kb[0]
    for i in range(len(x) - 1):
        h = dx[i]; a = np.array([[1, h / 2 * k[i + 1]], [-h / 2 * k[i + 1], 1]])
        rhs = np.array([e[i] + h / 2 * (dth[i] - k[i] * w[i] + dth[i + 1]), w[i] + h / 2 * k[i] * e[i]])
        e[i + 1], w[i + 1] = np.linalg.solve(a, rhs)
    e1 = dth - k * w; e2 = dk[0] - kb[1] * w - k ** 2 * e
    # directions p = <G,nu>, n = <e1,nu> with derivatives (a' = -k b, b' = [Gamma] + k a)
    dirs = []
    for Z in (G, np.ones_like(G)):
        a0 = np.real(Z * np.conj(nu)); b0 = np.real(Z * np.conj(T)); gam = 1.0 if Z is G else 0.0
        a1 = -k * b0; b1 = gam + k * a0; a2 = -kb[1] * b0 - k * b1
        dirs.append((a0, a1, a2))
    return y, e, e1, e2, dirs, dk
def proj(tr, full=False):
    y, e, e1, e2, dirs, dk = state(tr)
    Bs = [chi * g[0] for g in dirs]; Gm = np.array([[integ(rho * a * b) for b in Bs] for a in Bs])
    cf = np.linalg.solve(Gm, [integ(rho * a * chi * e) for a in Bs])
    E = [e - cf @ np.array([g[i] for g in dirs]) for i, e in enumerate((e, e1, e2))] if False else \
        [e - sum(cf[j] * dirs[j][0] for j in range(2)), e1 - sum(cf[j] * dirs[j][1] for j in range(2)), e2 - sum(cf[j] * dirs[j][2] for j in range(2))]
    A = np.sqrt(integ(rho * (chi * E[0]) ** 2))
    if not full: return A
    chi1 = np.gradient(chi, x); chi2 = np.gradient(chi1, x)
    Bn = np.sqrt(integ(rho * (chi2 * E[0] + 2 * chi1 * E[1] + chi * E[2]) ** 2))
    V0 = np.sqrt(integ(W * (zeta * E[1]) ** 2))
    far = x > 34.25
    return dict(A=A, B=Bn, V0=V0, cf=cf, tau=d['tau0'] + tr * d['dt'], k=K0 / np.exp(y[n0 + 1]), t=y[n0 + 3],
                esup_tip=np.abs(E[0][~far]).max(), esup_far=np.abs(E[0][far]).max(), e1sup_far=np.abs(E[1][far]).max(),
                dksup=np.abs(dk[0]).max(), ell=(ellP, y[n0]))
g = np.linspace(-0.1, 1.1, 49); A = [proj(t) for t in g]; j = int(np.argmin(A))
r = minimize_scalar(proj, bounds=(g[max(j - 1, 0)], g[min(j + 1, len(g) - 1)]), method='bounded', options=dict(xatol=1e-12))
o = proj(r.x, True); mue, eps = 0.036, 1e-5
o['Y'] = np.sqrt(o['A'] ** 2 + eps * o['B'] ** 2 + mue * o['V0'] ** 2)
print(f"P1 {ck[-15:-5]} (k {K0/lam:.6f}) vs base slab {s}: " + ', '.join(f'{k}={v:.4g}' if np.isscalar(v) else f'{k}={v}' for k, v in o.items()))
