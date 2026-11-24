"""Float: quadratic forcing Q2(psi3^) of the graph equation (Q(alpha psi3^) = alpha^2 Q2 + O(alpha^3)) and its
rho-norm, plus <chi Q2, f3> (needs F: use the float adjoint eigenfunction from refine data if available)."""
import numpy as np, sympy as sp, pickle, json, sys
sys.path.insert(0, '/home/claude/s1spec/step3nl')
G = pickle.load(open('graph_eq.pkl', 'rb'))
quad = G['quad']
s = sp.symbols('sigma'); U = sp.Function('U')(s)
p, q, k, m = [sp.Function(n)(s) for n in ('p', 'q', 'k', 'm')]
Ud = [U] + [sp.Derivative(U, (s, n)) for n in range(1, 5)]
J = sp.symbols('u0:5'); P_ = sp.symbols('pp qq kk mm')
expr = quad
for n in range(4, -1, -1): expr = expr.subs(Ud[n], J[n])
expr = expr.subs({p: P_[0], q: P_[1], k: P_[2], m: P_[3]})
f = sp.lambdify(list(J) + list(P_), expr, 'numpy')
# psi3^ jets on a grid from the float Chebyshev pieces + extension
from numpy.polynomial import chebyshev as Ch
Jd = json.load(open('psi3md2_14.json'))
def jets(x):
    for pc in Jd['pieces']:
        a, b = pc['a'], pc['b']
        if a <= x <= b:
            c = np.array(pc['coef']); sc = -2 / (b - a); y = 1 + sc * (x - a)
            return [Ch.chebval(y, Ch.chebder(c, n) if n else c) * sc**n for n in range(5)]
    # extension beyond X = 28: c1 g(s) + c2 E(s), s = x + stip
    A = 1 - 4 * Jd['lam']; aj = [1.0]
    for jj in range(1, 4): aj.append(aj[-1] * (A - 4*jj + 4) * (A - 4*jj + 3) * (A - 4*jj + 2) * (A - 4*jj + 1) / jj)
    sv = x + Jd['stip']; out = []
    for d in range(5):
        v = 0.0
        for j in range(4):
            e = A - 4 * j; fct = 1.0
            for i in range(d): fct *= (e - i)
            v += aj[j] * fct * sv ** (e - d)
        out.append(Jd['c_poly'] * v)
    return out
T0 = 10.249073; st = Jd['stip']
Dm = np.loadtxt('mc_prof.txt')[::5]
xs = (8 - (-2.249073)) * (1 - (Dm[:, 0] + Dm[:, 1]) / 2 / T0)
Pm = lambda i: (Dm[:, 2 + 2 * i] + Dm[:, 3 + 2 * i]) / 2
vals = []
for i, x in enumerate(xs):
    jt = jets(x); vals.append(f(*jt, Pm(0)[i], Pm(1)[i], Pm(2)[i], Pm(3)[i]))
vals = np.array(vals)
xx = np.linspace(10.25, 34.25, 4000)
vb = np.array([f(*jets(x), 0.0, x + st, 0.0, 0.0) for x in xx])
rho = lambda x: np.exp(-(1 + x * x) ** (2 / 3) / 4)
o = np.argsort(xs)
n2 = np.trapezoid(rho(xs[o]) * vals[o] ** 2, xs[o]) + np.trapezoid(rho(xx) * vb ** 2, xx)
print(f"||Q2(psi3^)||_rho (half line) = {np.sqrt(n2):.4f};  |Q2| at x=0,2,5,10,20: ", [float(f(*jets(t), *( [Pm(0)[np.argmin(abs(xs-t))], Pm(1)[np.argmin(abs(xs-t))], Pm(2)[np.argmin(abs(xs-t))], Pm(3)[np.argmin(abs(xs-t))]] if t<=10 else [0, t+st, 0, 0]))) for t in [0.01, 2, 5, 10, 20]])
