#!/usr/bin/env python3
"""Exact (sympy) checks of the integration-by-parts identities used in EI_SPEC.md.
For each identity  integrand - reduced  must be a total x-derivative: checked by the Euler operator (variational derivative w.r.t. the
unknown function(s)) being identically 0, or the identity is pointwise (difference simplifies to 0).
No floating point.  Writes ei_identities.json (args block with sha256 of this script).  Run: nice -n 19 python3 ei_identities.py"""
import sympy as sp, json, sys, hashlib, os, time
t0 = time.time()
x = sp.symbols('x')
psi, c, Q, zeta, g, v, w = [sp.Function(n)(x) for n in ('psi', 'c', 'Q', 'zeta', 'g', 'v', 'w')]
W = sp.exp(-psi)
d = lambda f, n=1: sp.diff(f, x, n)

def euler(Lag, f, order=6):
    """Euler operator sum_k (-1)^k d^k/dx^k dLag/df^(k)"""
    E = 0
    for kk in range(order + 1):
        fk = f if kk == 0 else sp.Derivative(f, (x, kk))
        E += (-1) ** kk * sp.diff(sp.diff(Lag, fk), x, kk)
    return sp.simplify(sp.expand(E.doit()) / W)

res = {}
# (I1) level 0, top order:  W v (c v''')'  ~  W c v''^2 + W (c' - 2 psi' c) v' v'' - W ((psi' c)' - psi'^2 c) v v''
lhs = W * v * d(c * d(v, 3))
rhs = W * c * d(v, 2) ** 2 + W * (d(c) - 2 * d(psi) * c) * d(v) * d(v, 2) - W * (d(d(psi) * c) - d(psi) ** 2 * c) * v * d(v, 2)
res['I1_level0_top'] = str(euler(lhs - rhs, v))
# (I2) commutator (pointwise): zeta (c g''')' = (c (zeta g)''')' - c K4 - c' K3
K3 = 3 * d(zeta) * d(g, 2) + 3 * d(zeta, 2) * d(g) + d(zeta, 3) * g
K4 = 4 * d(zeta) * d(g, 3) + 6 * d(zeta, 2) * d(g, 2) + 4 * d(zeta, 3) * d(g) + d(zeta, 4) * g
res['I2_quasilinear_commutator'] = str(sp.simplify(sp.expand(zeta * d(c * d(g, 3)) - (d(c * d(zeta * g, 3)) - c * K4 - d(c) * K3))))
# (I3) D2v := (W v'')''/W = v'''' - 2 psi' v''' + (psi'^2 - psi'') v''   (pointwise)
res['I3_D2v'] = str(sp.simplify(sp.expand(d(W * d(v, 2), 2) / W - (d(v, 4) - 2 * d(psi) * d(v, 3) + (d(psi) ** 2 - d(psi, 2)) * d(v, 2)))))
# (I4) level-2 pairing: W v'' f'' ~ (W v'')'' f  for a generic f (both sides as functionals of (v,f)): Euler operators in v and in f vanish
f = sp.Function('f')(x)
diff_ = W * d(v, 2) * d(f, 2) - d(W * d(v, 2), 2) * f
res['I4_level2_pairing_Ev'] = str(euler(diff_, v)); res['I4_level2_pairing_Ef'] = str(euler(diff_, f))
# (I5) transport, level 0:  W v (Q v' + (Q'-1) v) ~ W (psi' Q/2 + Q'/2 - 1) v^2
res['I5_transport0'] = str(euler(W * v * (Q * d(v) + (d(Q) - 1) * v) - W * (d(psi) * Q / 2 + d(Q) / 2 - 1) * v ** 2, v))
# (I6) transport, level 2: (W v'')'' (Q v' + (Q'-1) v) ~ W[(psi' Q/2 + 5Q'/2 - 1) v''^2 + 2 Q'' v' v'' - Q'' v v''' + psi' Q'' v v'']
lhs = d(W * d(v, 2), 2) * (Q * d(v) + (d(Q) - 1) * v)
rhs = W * ((d(psi) * Q / 2 + sp.Rational(5, 2) * d(Q) - 1) * d(v, 2) ** 2 + 2 * d(Q, 2) * d(v) * d(v, 2) - d(Q, 2) * v * d(v, 3) + d(psi) * d(Q, 2) * v * d(v, 2))
res['I6_transport2'] = str(euler(lhs - rhs, v))
# (I7) translation-type term (coefficient C, no -v):  level 0: W v (C v' + C' v) ~ W (psi' C/2 + C'/2) v^2 ; level 2 as (I6) without the "-1"
C = sp.Function('C')(x)
res['I7_b0'] = str(euler(W * v * (C * d(v) + d(C) * v) - W * (d(psi) * C / 2 + d(C) / 2) * v ** 2, v))
lhs = d(W * d(v, 2), 2) * (C * d(v) + d(C) * v)
rhs = W * ((d(psi) * C / 2 + sp.Rational(5, 2) * d(C)) * d(v, 2) ** 2 + 2 * d(C, 2) * d(v) * d(v, 2) - d(C, 2) * v * d(v, 3) + d(psi) * d(C, 2) * v * d(v, 2))
res['I7_b2'] = str(euler(lhs - rhs, v))
# (I8) GN identity: W v''^4 = (W v' v''^3)' - 3 W v' v''^2 v''' + W psi' v' v''^3   (pointwise)
res['I8_GN'] = str(sp.simplify(sp.expand(W * d(v, 2) ** 4 - (d(W * d(v) * d(v, 2) ** 3) - 3 * W * d(v) * d(v, 2) ** 2 * d(v, 3) + W * d(psi) * d(v) * d(v, 2) ** 3))))
# (I9) (L w)' = L1 w' + c0' w  with L = -d^4 - k^2 d^2 - (5km+q/4) d + c0, c0 = 1/4 - 3m^2 - kp, L1 = -d^4 - k^2 d^2 - (7km+q/4) d - 8m^2 - (5/2)kp
k, m, p, q = [sp.Function(n)(x) for n in ('k', 'm', 'p', 'q')]
ode = {sp.Derivative(p, x): -k * q, sp.Derivative(q, x): 1 + k * p, sp.Derivative(k, x): m, sp.Derivative(m, x): p / 4}
c0 = sp.Rational(1, 4) - 3 * m ** 2 - k * p
Lw = -d(w, 4) - k ** 2 * d(w, 2) - (5 * k * m + q / 4) * d(w) + c0 * w
gg = d(w)
L1g = -d(gg, 4) - k ** 2 * d(gg, 2) - (7 * k * m + q / 4) * d(gg) + (-8 * m ** 2 - sp.Rational(5, 2) * k * p) * gg
e = sp.expand((d(Lw) - L1g - d(c0) * w).doit().subs(ode).doit().subs(ode))
res['I9_L1'] = str(sp.simplify(e))
# (I10) forcing / N pairing, level 0:  W v zeta N' ~ -(W zeta v)' N  (Euler in N and v of the difference vanish)
N = sp.Function('N')(x)
diff_ = W * v * zeta * d(N) + d(W * zeta * v) * N
res['I10_EN'] = str(euler(diff_, N)); res['I10_Ev'] = str(euler(diff_, v))
# (I11) A1 split: A1(u) = -u + u' q/(1-k u);  A1(u0+w) - A1(u0) = -w + w' Qu + E_A w,  Qu = q/(1-k(u0+w)), E_A = u0' q k/((1-k(u0+w))(1-k u0))
u0 = sp.Function('u0')(x)
A1 = lambda u: -u + d(u) * q / (1 - k * u)
Qu = q / (1 - k * (u0 + w)); EA = d(u0) * q * k / ((1 - k * (u0 + w)) * (1 - k * u0))
res['I11_A1_split'] = str(sp.simplify(A1(u0 + w) - A1(u0) - (-w + d(w) * Qu + EA * w)))
ok = all(vv == '0' for vv in res.values())
out = dict(args=dict(argv=sys.argv, script_sha256=hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest(), sympy=sp.__version__),
           residuals=res, all_zero=ok, runtime_s=round(time.time() - t0, 1))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ei_identities.json'), 'w'), indent=1)
print(json.dumps(res, indent=1)); print('ALL ZERO' if ok else 'FAIL')
