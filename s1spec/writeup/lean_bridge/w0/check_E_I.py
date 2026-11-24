"""W0 (SCOPE.md): sympy check of (E1)-(E4) and the energy identity (I) for curve diffusion flow.
Curve gamma(t,theta)=(x,y), normal N = i T, T = gamma_theta/v, v=|gamma_theta|, kappa = <g'', i g'>/v^3.
Flow gamma_t = V N with V arbitrary smooth (E1-E4), then V = -kappa_ss (I)."""
import sympy as sp
th = sp.symbols('theta')
K = 7
x = sp.Function('x')(th); y = sp.Function('y')(th); V = sp.Function('V')(th)
xs = [sp.diff(x, th, k) for k in range(K)]
ys = [sp.diff(y, th, k) for k in range(K)]
v = sp.sqrt(xs[1]**2 + ys[1]**2)
Tx, Ty = xs[1]/v, ys[1]/v
Nx, Ny = -Ty, Tx                      # N = i T
kap = (xs[1]*ys[2] - ys[1]*xs[2]) / v**3
ds = lambda f: sp.diff(f, th) / v
# time derivative: gamma_t = V N
xt, yt = V*Nx, V*Ny
def Dt(F):
    out = 0
    for k in range(1, K):
        out += sp.diff(F, xs[k]) * sp.diff(xt, th, k) + sp.diff(F, ys[k]) * sp.diff(yt, th, k)
    out += sp.diff(F, x) * xt + sp.diff(F, y) * yt
    return out
def z(e): return sp.simplify(sp.together(e)) == 0
res = {}
res['E1 v_t = -k V v'] = z(Dt(v) + kap*V*v)
f = sp.Function('f')
# commutator on a generic function of the jet: use f = kappa (any jet-expression works the same way); test on x and kappa
for name, F in [('x', x), ('kappa', kap)]:
    res[f'E2 [d_t,d_s]{name} = k V {name}_s  (time-indep. dep. through jet)'] = z(Dt(ds(F)) - ds(Dt(F)) - kap*V*ds(F))
res['E3a T_t = V_s N (x-comp)'] = z(Dt(Tx) - ds(V)*Nx)
res['E3b T_t = V_s N (y-comp)'] = z(Dt(Ty) - ds(V)*Ny)
res['E3c N_t = -V_s T (x-comp)'] = z(Dt(Nx) + ds(V)*Tx)
res['E4 k_t = V_ss + k^2 V'] = z(Dt(kap) - ds(ds(V)) - kap**2*V)
for k_, val in res.items(): print(('PASS' if val else 'FAIL'), k_)

# (I): in arclength parametrisation, k = k(s). Pointwise derivative of the integrand k_s^2 v:
# d/dt(k_s^2 v) = (2 k_s (k_s)_t + k_s^2 (-k V)) v, (k_s)_t = (k_t)_s + k V k_s, V = -k_ss, k_t = -k_ssss - k^2 k_ss
s = sp.symbols('s'); k = sp.Function('k')(s)
d = lambda e, n=1: sp.diff(e, s, n)
Vs = -d(k, 2)
kt = d(Vs, 2) + k**2*Vs
integrand = 2*d(k)*(d(kt) + k*Vs*d(k)) - k*Vs*d(k)**2
rhs = -2*d(k, 3)**2 + 2*k**2*d(k, 2)**2 + sp.Rational(1, 3)*d(k)**4
diff = sp.expand(integrand - rhs)
# Euler operator: diff is a total derivative iff its variational derivative vanishes
from sympy.calculus.euler import euler_equations
ee = euler_equations(diff, [k], s)
ok = all((e is sp.true) or (getattr(e,'lhs',None) is not None and sp.simplify(e.lhs - e.rhs) == 0) for e in ee) if ee else True
print(('PASS' if ok else 'FAIL'), '(I) integrand - RHS is a total s-derivative (Euler operator = 0)', ee)
# explicit antiderivative (the integration-by-parts steps used in the Lean proof)
Q = -2*d(k)*d(k,4) + 2*d(k,2)*d(k,3) - 2*k**2*d(k)*d(k,2) - k*d(k)**3/3
print(('PASS' if sp.expand(d(Q) - diff) == 0 else 'FAIL'), '(I) explicit: d/ds[-2k_s k_ssss + 2k_ss k_sss - 2k^2 k_s k_ss - k k_s^3/3] = integrand - RHS')
