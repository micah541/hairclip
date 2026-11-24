"""Float 8x8 Gram with ODE-exact derivatives (no finite differences of p, sin theta), no mask; compare with rigorous mid values.
Uses profile grid prof_float.npz (from check_frame.py) and psi3md2_14.json (psi3^) as in mod_form2_float.py."""
import sys, json, numpy as np, sympy as sp
from math import comb
from numpy.polynomial import chebyshev as Ch
sys.path.insert(0, '/home/claude/s1spec'); import evans as E
d = np.load('prof_float.npz'); x, p, q, k, m, th = [d[n] for n in ('x', 'p', 'q', 'k', 'm', 'th')]
h = x[1] - x[0]; sig = x + E.STIP
xs = x <= 31.0 + 1e-9
P1, P2, P3, P4, K, M, Q, Pp, c, s = sp.symbols('P1 P2 P3 P4 K M Q Pp c s'); dd = sp.symbols('d0:5')
Em = [sp.Integer(1), -P1, P1**2 - P2, -P1**3 + 3*P1*P2 - P3, P1**4 - 6*P1**2*P2 + 3*P2**2 + 4*P1*P3 - P4]
der = lambda n: sum(comb(n, j) * Em[n - j] * dd[j] for j in range(n + 1))
Ls = sp.expand(-der(4) - K**2 * der(2) + (K*M + Q/4) * der(1) + dd[0] / 2)
aS = [sp.expand(Ls.coeff(dd[n])) for n in range(5)]
rules = {Pp: -K*Q, Q: 1 + K*Pp, K: M, M: Pp/4, c: -K*s, s: K*c}
def D(e): return sp.expand(sum(sp.diff(e, v) * r for v, r in rules.items()))
def ders(e):
    o = [e]
    for _ in range(4): o.append(D(o[-1]))
    return o
args = (P1, P2, P3, P4, K, M, Q, Pp, c, s)
lam = lambda e: sp.lambdify(args, e, 'numpy')
phi = (1 + x**2) ** (2 / 3) / 4
E_ = (1 + x*x) ** (-1/3); X2 = x*x; E2 = E_**2; E4 = E2**2; E7 = E4*E2*E_; E10 = E7*E2*E_
PP = [x*E_/3, E_/3 - 2*X2*E4/9, -2*x*E4/3 + 16*X2*x*E7/27, -2*E4/3 + 32*X2*E7/9 - 224*X2*X2*E10/81]
rho = np.exp(-phi); w8 = np.full(len(x), h); w8[0] = w8[-1] = h / 2; w8 = w8 * rho
ip = lambda a, b: np.sum(w8 * a * b)
cth, sth = np.cos(th), np.sin(th)
base = (*PP, k, m, q, p, cth, sth)
a = [lam(e)(*base) * np.ones_like(x) for e in aS]
dp = [lam(e)(*base) * np.ones_like(x) for e in ders(Pp)]
d2 = [lam(e)(*base) * np.ones_like(x) for e in ders(-s)]
u = np.clip((32 - sig) / 3, 0, 1); S = lambda t: t**4 * (35 - 84*t + 70*t**2 - 20*t**3)
t_ = sp.symbols('t'); Ssym = t_**4 * (35 - 84*t_ + 70*t_**2 - 20*t_**3)
chi = [sp.lambdify(t_, sp.diff(Ssym, t_, j), 'numpy')(u) * (-1/3)**j * ((sig > 29) & (sig < 32)) if j else S(u) for j in range(5)]
def fder(ps): return [sum(comb(n, j) * chi[j] * ps[n - j] for j in range(n + 1)) for n in range(5)]
g1, g2 = fder(dp), fder(d2)
out = {}
for key, g in (("1", g1), ("2", g2)):
    out["Lf" + key] = sum(a[n] * g[n] for n in range(5))
    out["Df" + key] = g[4] - 2 * PP[0] * g[3] + (PP[0]**2 - PP[1]) * g[2]
    out["f" + key] = g[0]
Jf = json.load(open('/home/claude/s1spec/step3nl/psi3md2_14.json'))
def ev(xv, dn):
    o = np.zeros_like(xv)
    for pc in Jf['pieces']:
        a_, b_ = pc['a'], pc['b']; msk = (xv >= a_) & (xv <= b_); cc = np.array(pc['coef'])
        for _ in range(dn): cc = Ch.chebder(cc)
        o[msk] = Ch.chebval(1 - 2 * (xv[msk] - a_) / (b_ - a_), cc) * (-2 / (b_ - a_)) ** dn
    return o
y0, y1 = ev(x, 0), ev(x, 1)
out["gA"] = chi[0] * (y0 - q * y1); out["gB"] = chi[0] * y1 * cth
NAMES = ["Lf1", "Lf2", "gA", "gB", "Df1", "Df2", "f1", "f2"]
Gf = np.array([[ip(out[a_] * xs, out[b_] * xs) for b_ in NAMES] for a_ in NAMES])
np.save('gram8_float.npy', Gf)
if __name__ == "__main__":
    np.set_printoptions(linewidth=200, precision=5, suppress=True)
    print(Gf)
    R = json.load(open('gram8_raw.json'))['G']; Gm = np.array([[(R[i][j][0] + R[i][j][1]) / 2 for j in range(8)] for i in range(8)])
    print("rigorous mid - float:"); print(Gm - Gf)
