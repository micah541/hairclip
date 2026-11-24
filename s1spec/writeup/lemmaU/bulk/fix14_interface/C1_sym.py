#!/usr/bin/env python3
"""C1_sym.py -- exact sympy identities used by C1_identification.tex (fix round 14, Nov 20 2026).  New file; reads nothing.

 Y1  normal flow kinematics, general parametrisation X(t,p) = (x, y), X_t = V N, N = J X_p/|X_p| (J = +pi/2 rotation):
       theta_t|_p = V_p/|X_p| (= V_s)   and   d/dt |X_p| = -k V |X_p|,   k = (x_p y_pp - y_p x_pp)/|X_p|^3.
 Y2  arclength-proportional variables: theta = -pi/2 + a + th(a), d/ds = (2 pi/L) d/da  =>  k = (2pi/L)(1 + v), k_s = (2pi/L)^2 w,
       k_ss = (2pi/L)^3 w', k_sss = (2pi/L)^4 w'',  v = th', w = th''.
 Y3  the derivative identity d/da[(1 + v) w] = (1 + v) w' + w^2 (so int_0^a (1+v) w' = (1+v)w - (1+v(0))w(0) - int_0^a w^2, and over a
       period int (1+v) w' = -int w^2).
 Y4  tangential drift in the sigma clock (d sigma = (2pi/L)^4 dt):  d_sigma a|_p = int_0^a (1+v) w' + a <w^2>,
       assembled from d_t s|_p = -int_0^s k V ds',  L_t = -int k V ds,  V = -k_ss  (symbolic in the antiderivative symbols).
 Y5  the angle equation:  theta_sigma|_a = -th'''' - (1+v)^2 w + (1+v) int_0^a (w^2 - <w^2>)   (using w(0) = 0)
       and  (1+v)^2 w = th'' + d/da (v^2 + v^3/3),  so  th_sigma = -th'''' - th'' + N(th),  N = -(v^2 + v^3/3)' + (1+v) P(w^2 - <w^2>).
 Y6  for an odd trigonometric polynomial th = sum_{n<=3} u_n sin(n a) (symbolic u_n):  w^2 - <w^2> is even, its antiderivative from 0 is
       odd and has mean 0, i.e. equals the zero-mean antiderivative P of step 4;  <w^2> = (1/2) sum n^4 u_n^2.
 Y7  slaved quantities:  d ln L/d sigma = -(1/2 pi) int w^2 = -(1/2) sum n^4 u_n^2;  X0' := d_sigma X(a=0) = -(L/2pi) th'''(0) = (L/2pi) sum n^3 u_n,
       N(a=0) = 1 (theta(0) = -pi/2);  t' = (L/2pi)^4.
 Y8  normal-graph curvature (step 4 curvjets.py / nl7.py): X = G + e N over a unit-speed base (T' = kN, N' = -kT), a = 1 - k e,
       J^2 = a^2 + e'^2:  kappa = (X' x X'')/J^3 = (a e'' - e' a' + k (a^2 + e'^2))/J^3.
 Y9  tube lemma identities: c = S(X) + n nu(X) over a unit-speed S (T_S' = k nu, nu' = -k T_S):  <c', T_S(X)> = X'(1 - k n),
       <c', nu(X)> = n';  and for fixed c, d/dX <c - S(X), T_S(X)> = -1 + k <c - S(X), nu(X)>.
 Y10 symmetry: theta(-a) = -pi - theta(a) for odd th, hence e^{i theta(-a)} = -e^{-i theta(a)}.
Output: C1_sym.json (args block, all_pass).  Run: PYTHONHASHSEED=0 nice -n 19 python3 -B C1_sym.py
"""
import sys, os, json, hashlib, platform
sys.dont_write_bytecode = True
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
res = {}


def check(name, expr):
    z = sp.simplify(sp.expand(expr))
    ok = (z == 0)
    res[name] = bool(ok)
    print('%-72s %s' % (name, ok), flush=True)
    assert ok, (name, z)


# ---------------- Y1
t, p = sp.symbols('t p', real=True)
x = sp.Function('x')(t, p); y = sp.Function('y')(t, p); V = sp.Function('V')(t, p)
xp, yp = sp.diff(x, p), sp.diff(y, p); r = sp.sqrt(xp ** 2 + yp ** 2)
Nx, Ny = -yp / r, xp / r                                  # N = J T
xt, yt = V * Nx, V * Ny                                    # X_t = V N
xpt, ypt = sp.diff(xt, p), sp.diff(yt, p)                  # (X_t)_p = (X_p)_t
theta_t = (xp * ypt - yp * xpt) / r ** 2                   # d/dt atan2(y_p, x_p)
check('Y1a theta_t|_p - V_p/|X_p|', theta_t - sp.diff(V, p) / r)
kk = (xp * sp.diff(y, p, 2) - yp * sp.diff(x, p, 2)) / r ** 3
rt = (xp * xpt + yp * ypt) / r
check('Y1b d_t|X_p| + k V |X_p|', rt + kk * V * r)

# ---------------- Y2
a, L = sp.symbols('a L', positive=True)
th = sp.Function('th')(a)
theta = -sp.pi / 2 + a + th
c = 2 * sp.pi / L
Ds = lambda f: c * sp.diff(f, a)
v, w = sp.diff(th, a), sp.diff(th, a, 2)
k = Ds(theta)
check('Y2a k - (2pi/L)(1+v)', k - c * (1 + v))
check('Y2b k_s - (2pi/L)^2 w', Ds(k) - c ** 2 * w)
check('Y2c k_ss - (2pi/L)^3 w\'', Ds(Ds(k)) - c ** 3 * sp.diff(w, a))
check('Y2d k_sss - (2pi/L)^4 w\'\'', Ds(Ds(Ds(k))) - c ** 4 * sp.diff(w, a, 2))

# ---------------- Y3
check('Y3 d/da[(1+v)w] - (1+v)w\' - w^2', sp.diff((1 + v) * w, a) - (1 + v) * sp.diff(w, a) - w ** 2)

# ---------------- Y4: drift.  Symbols: F(a) = int_0^a (1+v) w' da', Wm = <w^2> = (1/2pi) int_0^{2pi} w^2.
F, Wm = sp.symbols('F Wm', real=True)
kV_ds_partial = -c ** 3 * F            # int_0^s k V ds' with V = -k_ss:  k V ds = -(2pi/L)(1+v)(2pi/L)^3 w' (L/2pi) da
Lt = c ** 3 * (-2 * sp.pi * Wm)        # L_t = -int k V ds = (2pi/L)^3 int_0^{2pi} (1+v) w' = -(2pi/L)^3 int w^2  (Y3 over a period)
ds_t = -kV_ds_partial                  # d_t s|_p
s_ = a * L / (2 * sp.pi)
dalpha_dt = 2 * sp.pi * ds_t / L - 2 * sp.pi * s_ * Lt / L ** 2
check('Y4 d_sigma a|_p - (F + a <w^2>)', dalpha_dt / c ** 4 - (F + a * Wm))
# check the sign of the kV integrand directly: k * V with V = -k_ss, ds = (L/2pi) da
check('Y4b k V ds/da + (2pi/L)^3 (1+v) w\'', k * (-Ds(Ds(k))) * L / (2 * sp.pi) + c ** 3 * (1 + v) * sp.diff(w, a))

# ---------------- Y5: angle equation (theta_t|_p = V_s = -k_sss; theta_t|_a = theta_t|_p - (d_t a|_p) theta_a)
I0 = sp.symbols('I0', real=True)       # int_0^a w^2
w0 = sp.symbols('w0', real=True)       # w(0) (= 0 for odd th)
Fexpr = (1 + v) * w - (1 + v.subs(a, 0)) * w0 - I0      # Y3 integrated
theta_sig = -Ds(Ds(Ds(k))) / c ** 4 - (Fexpr + a * Wm) * sp.diff(theta, a)
target = -sp.diff(th, a, 4) - (1 + v) ** 2 * w + (1 + v) * (I0 - a * Wm)
check('Y5a theta_sigma|_a (w(0)=0) - [-th\'\'\'\' - (1+v)^2 w + (1+v)(I0 - a<w^2>)]', (theta_sig - target).subs(w0, 0))
check('Y5b (1+v)^2 w - th\'\' - (v^2 + v^3/3)\'', (1 + v) ** 2 * w - w - sp.diff(v ** 2 + v ** 3 / 3, a))

# ---------------- Y6, Y7 on an odd trigonometric polynomial
u1, u2, u3 = sp.symbols('u1 u2 u3', real=True); al = sp.symbols('alpha', real=True)
U = [u1, u2, u3]
thp = sum(U[n - 1] * sp.sin(n * al) for n in (1, 2, 3))
vp, wp = sp.diff(thp, al), sp.diff(thp, al, 2)
mean_w2 = sp.integrate(sp.expand(wp ** 2), (al, 0, 2 * sp.pi)) / (2 * sp.pi)
check('Y6a <w^2> - (1/2) sum n^4 u_n^2', mean_w2 - sp.Rational(1, 2) * sum(n ** 4 * U[n - 1] ** 2 for n in (1, 2, 3)))
f = sp.expand(sp.expand_trig(wp ** 2) - mean_w2)
check('Y6b w^2 - <w^2> even', sp.expand_trig(f.subs(al, -al) - f))
z = sp.symbols('z', real=True)
Pf = sp.integrate(sp.expand(f.subs(al, z)), (z, 0, al))
check('Y6c int_0^a (w^2-<w^2>) odd', sp.expand_trig(Pf.subs(al, -al) + Pf))
check('Y6d int_0^a (w^2-<w^2>) has mean 0', sp.integrate(sp.expand(Pf), (al, 0, 2 * sp.pi)))
check('Y6e w(0) = 0 for odd th', wp.subs(al, 0))
Lsym = sp.symbols('L', positive=True)
dlnL = (-(2 * sp.pi / Lsym) ** 3 * 2 * sp.pi * mean_w2) / Lsym * (Lsym / (2 * sp.pi)) ** 4
check('Y7a d ln L/d sigma + (1/2) sum n^4 u_n^2', dlnL + sp.Rational(1, 2) * sum(n ** 4 * U[n - 1] ** 2 for n in (1, 2, 3)))
kss0 = (2 * sp.pi / Lsym) ** 3 * sp.diff(wp, al).subs(al, 0)
X0dot = -kss0 * (Lsym / (2 * sp.pi)) ** 4
check('Y7b X0\' - (L/2pi) sum n^3 u_n', X0dot - Lsym / (2 * sp.pi) * sum(n ** 3 * U[n - 1] for n in (1, 2, 3)))
check('Y7c N(0) = i e^{i theta(0)} = 1', sp.I * sp.exp(sp.I * (-sp.pi / 2)) - 1)

# ---------------- Y8 normal-graph curvature
xs = sp.symbols('xs', real=True)
phi = sp.Function('phi')(xs); e = sp.Function('e')(xs)
kb = sp.diff(phi, xs)
T = sp.Matrix([sp.cos(phi), sp.sin(phi)]); N = sp.Matrix([-sp.sin(phi), sp.cos(phi)])
Gp = T                                                     # unit-speed base: G' = T (G itself not needed)
Xp = Gp + sp.diff(e, xs) * N + e * sp.diff(N, xs)
Xpp = sp.diff(Xp, xs)
cross = Xp[0] * Xpp[1] - Xp[1] * Xpp[0]
A_ = 1 - kb * e; J2 = A_ ** 2 + sp.diff(e, xs) ** 2
check('Y8a |X\'|^2 - (a^2 + e\'^2)', Xp.dot(Xp) - J2)
check('Y8b X\'xX\'\' - (a e\'\' - e\' a\' + k(a^2+e\'^2))',
      cross - (A_ * sp.diff(e, xs, 2) - sp.diff(e, xs) * sp.diff(A_, xs) + kb * J2))

# ---------------- Y9 tube lemma identities
xi = sp.symbols('xi', real=True)
X = sp.Function('X')(xi); n = sp.Function('n')(xi)
psi = sp.Function('psi')                                   # tangent angle of S as a function of its arclength
TS = lambda q: sp.Matrix([sp.cos(psi(q)), sp.sin(psi(q))]); nuS = lambda q: sp.Matrix([-sp.sin(psi(q)), sp.cos(psi(q))])
kS = lambda q: sp.diff(psi(q), q)
Sx = sp.Function('Sx'); Sy = sp.Function('Sy')
# c = S(X) + n nu(X), with dS/dq = T_S(q): differentiate using the chain rule, replacing S'(X) by T_S(X)
cprime = TS(X) * sp.diff(X, xi) + sp.diff(n, xi) * nuS(X) + n * sp.diff(nuS(X), xi)
check('Y9a <c\', T(X)> - X\'(1 - k n)', sp.simplify(cprime.dot(TS(X)) - sp.diff(X, xi) * (1 - kS(X) * n)))
check('Y9b <c\', nu(X)> - n\'', sp.simplify(cprime.dot(nuS(X)) - sp.diff(n, xi)))
q = sp.symbols('q', real=True); cx, cy = sp.symbols('cx cy', real=True)
Sq = sp.Matrix([sp.Function('S1x')(q), sp.Function('S1y')(q)])
g_ = (sp.Matrix([cx, cy]) - Sq).dot(TS(q))
dg = sp.diff(g_, q).subs({sp.diff(Sq[0], q): sp.cos(psi(q)), sp.diff(Sq[1], q): sp.sin(psi(q))})
check('Y9c d/dX <c - S(X), T(X)> - (-1 + k <c - S, nu>)', dg - (-1 + kS(q) * (sp.Matrix([cx, cy]) - Sq).dot(nuS(q))))

# ---------------- Y10 symmetry of the angle
tho = sp.Function('tho')
tht = lambda q: -sp.pi / 2 + q + tho(q)
check('Y10a theta(-a) + pi + theta(a) for odd th', (tht(-a) + sp.pi + tht(a)).subs(tho(-a), -tho(a)))
check('Y10b e^{i(-pi - th)} + e^{-i th}', sp.exp(sp.I * (-sp.pi - xs)) + sp.exp(-sp.I * xs))

me = os.path.abspath(__file__)
out = dict(all_pass=all(res.values()), checks=res,
           args=dict(argv=sys.argv, env={k_: os.environ[k_] for k_ in ('PYTHONHASHSEED',) if k_ in os.environ},
                     python=platform.python_version(), sympy=sp.__version__,
                     script=dict(path=me, sha256=hashlib.sha256(open(me, 'rb').read()).hexdigest()), inputs={}))
json.dump(out, open(os.path.join(HERE, 'C1_sym.json'), 'w'), indent=1)
print('all_pass', out['all_pass'])
