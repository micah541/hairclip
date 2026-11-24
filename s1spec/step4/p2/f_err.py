"""float probe 2: normal error e of the tip-normalised P1 curve over the base (one slab, exact Lagrange in tau),
projected off {chi p, chi n} in L2_rho; A(tau) = ||chi e||_rho minimised over the clock."""
import sys, json, numpy as np
sys.path.insert(0, '/home/claude/s1spec/step4/t1')
import polstore
from rep import K0, DS1
STIP = -2.249077651121112
ck = sys.argv[1]; s = int(sys.argv[2])
rec = json.load(open(ck))
u = np.array([int(m) * 2.0 ** int(e) for m, e in rec['xbar']]); n = np.arange(1, len(u) + 1)
L = int(rec['L'][0]) * 2.0 ** int(rec['L'][1]); c = 2 * np.pi / L
def tha(al): return 1 + (np.cos(np.outer(al, n)) * (n * u)).sum(1)
lam = K0 / (-c * tha(np.array([np.pi]))[0])
x = np.linspace(0, 34.25, 6851); dx = x[1] - x[0]
def cumtrap(f): return np.r_[0, np.cumsum(0.5 * (f[1:] + f[:-1]) * dx)]
def curve(k):
    th = -np.pi / 2 + cumtrap(k); return th, -DS1 + cumtrap(np.exp(1j * th))
kP = -lam * c * tha(np.pi - c * lam * x); thP, GP = curve(kP)
rho = np.exp(-(1 + x ** 2) ** (2 / 3) / 4)
uu = np.clip((32 - (x + STIP)) / 3, 0, 1); chi = uu ** 4 * (35 - 84 * uu + 70 * uu ** 2 - 20 * uu ** 3)
SP, d = polstore.load_slab(f'/home/claude/s1spec/step4/t1/polished/slab_{s:03d}.pkl')
R = SP.R; nd = np.array(d['nodes_rel']); nd = nd if len(nd) == len(SP.Yn) else np.r_[0.0, nd]; n0 = R.n
Y = np.array([[float(v.mid()) for v in row] for row in SP.Yn])
def lag(tr):  # Lagrange weights at relative time tr
    w = np.ones(len(nd))
    for i in range(len(nd)):
        for j in range(len(nd)):
            if i != j: w[i] *= (tr - nd[j]) / (nd[i] - nd[j])
    return w
def ip(f, g): return (rho * f * g).sum() * dx
def evalA(tr, full=False):
    y = lag(tr) @ Y; ell = y[n0]
    kb = R.eval_at(y[:n0], ell, x, 0); thb, Gb = curve(kb)
    nu = 1j * np.exp(1j * thb)
    e = np.real((GP - Gb) * np.conj(nu))
    p = np.real(Gb * np.conj(nu)); nn = np.real(np.conj(nu))
    B = np.array([chi * p, chi * nn]); G = np.array([[ip(a, b) for b in B] for a in B])
    ce = chi * e; coef = np.linalg.solve(G, [ip(a, ce) for a in B]); r = ce - coef @ B
    A = np.sqrt(ip(r, r)); A0 = np.sqrt(ip(ce, ce))
    return (A, tr, A0, np.abs(e).max(), np.abs(kP - kb).max(), ell, y[n0+1:n0+4]) if full else A
from scipy.optimize import minimize_scalar
g = np.linspace(-0.2, 1.2, 141); j = int(np.argmin([evalA(t) for t in g]))
res = minimize_scalar(evalA, bounds=(g[max(j-1,0)], g[min(j+1,len(g)-1)]), method="bounded", options=dict(xatol=1e-13))
A, tr, A0, esup, ksup, ell, sc = evalA(res.x, True)
print(f'slab {s}: best tau = {d["tau0"] + tr * d["dt"]:.6f} (rel {tr:.3f}); A (proj. off p,n) = {A:.3e}; A raw = {A0:.3e}; '
      f'sup|e|[0,34.25] = {esup:.3e}; sup|dkappa| = {ksup:.3e}; base k = {K0/np.exp(sc[0]):.8g}, t = {sc[2]:.12g}; P1 k = {K0/lam:.8g}')
