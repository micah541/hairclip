"""Symbolic checks (sympy) of the T5 sign conventions.
Frame: x = arclength of the base from the tip, T = Gamma_x, nu = i T, kappa = theta_x, T_x = kappa nu, nu_x = -kappa T.
Base functions p = <G,nu>, q = <G,T>, n = <e1,nu>, t1 = <e1,T> (Frenet: p' = -k q, q' = 1 + k p, n' = -k t1, t1' = k n).
Base motion G_tau|_x = V nu + w T (w = transport varpi, w' = k V), theta_tau = V' + k w.
True curve Gamma = G + eps e nu, rescaled CDF: <Gamma_tau, nu_Gamma> = -kappa_ss + (1/4 - a)<Gamma, nu> + b <e1, nu>,
with a = a~, b = b~ and the base identity V = -k'' + (1/4 - a~) p + b~ n + r (r = eps r1)."""
import sympy as sp

x, eps, A, Bt = sp.symbols('x epsilon atil btil', real=True)
k = sp.Function('k')(x); e = sp.Function('e')(x); et = sp.Function('e_tau')(x)
p, q, n, t1, w, r1 = [sp.Function(s)(x) for s in ('p', 'q', 'n', 't1', 'w', 'r1')]
V = -sp.diff(k, x, 2) + (sp.Rational(1, 4) - A) * p + Bt * n + eps * r1
rules = {sp.Derivative(p, x): -k * q, sp.Derivative(q, x): 1 + k * p, sp.Derivative(n, x): -k * t1,
         sp.Derivative(t1, x): k * n, sp.Derivative(w, x): k * V}

def D(f, m=1):
    for _ in range(m):
        f = sp.diff(f, x)
        for _ in range(3): f = f.subs(rules)   # repeated, since substitutions create new derivatives
        f = sp.expand(f)
    return f

def Dv(v):                       # derivative of a frame vector [A along T, B along nu]
    a_, b_ = v; return [D(a_) - k * b_, D(b_) + k * a_]

def dot(u, v): return u[0] * v[0] + u[1] * v[1]
def ser(f, o=1): return sp.series(f, eps, 0, o + 1).removeO()

# true curve
Gam = [q, p + eps * e]
Gx = Dv(Gam)                                  # = [(1 - eps k e), eps e']
J = sp.sqrt(Gx[0]**2 + Gx[1]**2)
T = [Gx[0] / J, Gx[1] / J]; Nu = [-T[1], T[0]]
kap = dot(Dv(T), Nu) / J
Ds = lambda f: D(f) / J
kss = Ds(Ds(kap))
thtau = D(V) + k * w
Gt = [w - eps * e * thtau, V + eps * et]      # Gamma_tau at fixed x (nu~_tau = -theta_tau T~)
lhs = dot(Gt, Nu)
rhs = -kss + (sp.Rational(1, 4) - A) * dot(Gam, Nu) + Bt * dot([t1, n], Nu)
eq = ser(lhs - rhs, 1)
eq0 = sp.simplify(eq.coeff(eps, 0))
eq1 = sp.expand(eq.coeff(eps, 1))
sol = sp.solve(sp.Eq(eq1, 0), et)[0]
b2 = -k**2
b1 = -5 * k * sp.diff(k, x) - (sp.Rational(1, 4) - A) * q - Bt * t1 + w
b0 = sp.Rational(1, 4) - A - 3 * sp.diff(k, x)**2 - 4 * k * sp.diff(k, x, 2)
L = -sp.diff(e, x, 4) + b2 * sp.diff(e, x, 2) + b1 * sp.diff(e, x) + b0 * e
print('order eps^0 of the error equation (must be 0):', eq0)
print('e_tau - (L e - r1) =', sp.simplify(sp.expand(sol - (L - r1))))

# curvature form of the base flow used by T1: kappa_tau|_x = V'' + (k w)' with w = -k k' + Q + (1/4-a)(q - x) + b t1, Q' = k'^2
Q = sp.Function('Q')(x); rules[sp.Derivative(Q, x)] = sp.diff(k, x)**2
V0 = V.subs(eps, 0)
w_closed = -k * sp.diff(k, x) + Q + (sp.Rational(1, 4) - A) * (q - x) + Bt * t1
print("w' - k V (must be 0, so w = varpi when Q(0)=q(0)=t1(0)=0):", sp.simplify(D(w_closed) - sp.expand(k * V0)))
ktau = D(D(V0)) + D(k * w_closed)
ktau = ktau.subs(sp.Derivative(w, x), k * V0)
local = -sp.diff(k, x, 4) - k**2 * sp.diff(k, x, 2) - k * sp.diff(k, x)**2 + sp.diff(k, x) * Q - (sp.Rational(1, 4) - A) * (k + x * sp.diff(k, x))
print('kappa_tau - local form (must be 0):', sp.simplify(sp.expand(ktau - local)))
Phi = -sp.diff(k, x, 3) - k**2 * sp.diff(k, x) + k * Q - (sp.Rational(1, 4) - A) * x * k
print("Phi' - local form (flux form, must be 0):", sp.simplify(D(Phi) - local))
print('Phi - theta_tau (must be 0):', sp.simplify(sp.expand(Phi - (D(V0) + k * w_closed))))

# far end: Gamma~(ell(tau)) symmetry point: t1(ell) = 0, q(ell) = 0 for all tau, kappa'(ell) = 0.
# d/dtau t1(ell(tau),tau) = theta_tau n + ell_tau k n = 0 ; d/dtau q = (w + theta_tau p) + ell_tau (1 + k p) = 0
lt, th_t, pp, nn, kk, ww = sp.symbols('ell_tau theta_tau p_l n_l k_l w_l')
s = sp.solve([th_t * nn + lt * kk * nn, ww + th_t * pp + lt * (1 + kk * pp)], [lt, th_t], dict=True)[0]
print('far end: ell_tau =', s[lt], ' theta_tau =', s[th_t])
print('b1(ell) with q = t1 = kappa\' = 0 is w(ell) = varpi(ell) = -ell_tau:', sp.simplify(b1.subs({q: 0, t1: 0, sp.diff(k, x): 0}) - w))
# energy identity: d/dtau int_0^ell W e^2 contains the transport term int W 2 e b1 e' -> boundary W b1 e^2 at ell, plus moving-domain term ell_tau W e^2
print('boundary terms at ell: W e^2 (b1(ell) + ell_tau) =', sp.simplify(w + s[lt]))

# tip: Gamma(0) = -DS1 (on the axis), theta(0) = -pi/2: T(0) = -i, nu(0) = 1, p(0) = -DS1, n(0) = 1, q(0) = t1(0) = 0.
# tip fixed <=> V(0) = 0 <=> r(0) = 0:  0 = -k''(0) + (1/4-a)(-DS1) + b  =>  b = k''(0) + (1/4-a) DS1 ; kappa(0) = K0 fixed => a
DS1, k2_0, k4_0, K0 = sp.symbols('D_S1 k2_0 k4_0 K_0')
bsol = sp.solve(sp.Eq(-k2_0 + (sp.Rational(1, 4) - A) * (-DS1) + Bt, 0), Bt)[0]
print('tip: b =', bsol, ';  S1 (k2_0 = -DS1/4, a = 0):', sp.simplify(bsol.subs({k2_0: -DS1 / 4, A: 0})))
asol = sp.solve(sp.Eq(-k4_0 - K0**2 * k2_0 - (sp.Rational(1, 4) - A) * K0, 0), A)[0]
print('tip: a =', asol)
