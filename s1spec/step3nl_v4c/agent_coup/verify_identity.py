"""Verify (symbolically) -<[L,chi]w, chi w>_rho = int A w^2 + B w'^2 (difference is a total derivative: Euler operator = 0),
and that the compiled polynomials A/Wt, B/Wt of coup_common reproduce the stated A, B (checked numerically with exact sympy values)."""
import sympy as sp, random
from sympy.calculus.euler import euler_equations
x = sp.symbols('x')
w, rho, chi, kk, cc = [sp.Function(n)(x) for n in ('w', 'rho', 'chi', 'kk', 'cc')]
d = lambda e, n=1: sp.diff(e, x, n)
comm = -4*d(chi)*d(w,3) - 6*d(chi,2)*d(w,2) - (4*d(chi,3) + 2*kk**2*d(chi))*d(w) - (d(chi,4) + kk**2*d(chi,2) + cc*d(chi))*w
lhs = -rho*chi*w*comm
g1, g2, g3, g4 = rho*chi*d(chi), rho*chi*d(chi,2), rho*chi*d(chi,3), rho*chi*d(chi,4)
A = -2*d(g1,3) + 3*d(g2,2) - 2*d(g3) + g4 - d(rho*kk**2*chi*d(chi)) + rho*chi*(kk**2*d(chi,2) + cc*d(chi))
B = 6*d(g1) - 6*g2
diff_ = sp.expand(lhs - A*w**2 - B*d(w)**2)
eq = euler_equations(diff_, w, x)[0].lhs
print("Euler operator of (integrand - (A w^2 + B w'^2)) =", sp.simplify(sp.expand(eq.doit())))
# polynomial form vs generic A,B with rho = exp(-phi), k,m,p,q ODE functions
import coup_common as C
for a in (5, 7):
    Ap, Bp = C.AB_poly(a)
    phi = sp.Function('phi')(x); K, M, Pp, Q = [sp.Function(n)(x) for n in 'KMPQ']
    ode = {d(K): M, d(M): Pp/4, d(Pp): -K*Q, d(Q): 1 + K*Pp}
    def red(e):
        for _ in range(8):
            e = e.subs({sp.Derivative(K, (x, n)): d(sp.diff(K, x), n-1) for n in range(2, 5)}) if False else e
            e = e.subs({sp.Derivative(Pp, x): -K*Q, sp.Derivative(Q, x): 1 + K*Pp, sp.Derivative(K, x): M, sp.Derivative(M, x): Pp/4}).doit()
        return e
    rr = sp.exp(-phi); ch = chi
    # use cc = a K M + Q/4 (functions of x through K,M,Q), kk = K
    sub = {rho: rr, kk: K, cc: a*K*M + Q/4}
    Ag = A.subs(sub).doit(); Bg = B.subs(sub).doit()
    # replace chi, phi by explicit functions for numeric test
    chi_e = sp.exp(sp.sin(x)) ; phi_e = x**3/7 + sp.cos(x)
    def num(e):
        e = red(e).subs({chi: chi_e, phi: phi_e}).doit()
        return e
    Ag, Bg = num(Ag), num(Bg)
    syms = {C.cs[i]: sp.diff(chi_e, x, i) for i in range(5)}
    syms.update({C.Ps[i]: sp.diff(phi_e, x, i+1) for i in range(3)})
    syms.update({C.k_: K, C.m_: M, C.p_: Pp, C.q_: Q})
    Ap_e = Ap.subs(syms) * sp.exp(-phi_e); Bp_e = Bp.subs(syms) * sp.exp(-phi_e)
    for pt in (0.3, 1.1, -0.7):
        vals = {x: pt, K: 0.37, M: -0.21, Pp: 0.55, Q: 1.9}
        e1 = (Ag - Ap_e).subs(vals).evalf(30); e2 = (Bg - Bp_e).subs(vals).evalf(30)
        print(f"a={a} x={pt}: |A_generic - rho*Apoly| = {float(abs(e1)):.2e} (|A|={float(abs(Ag.subs(vals).evalf())):.2e}),  B diff = {float(abs(e2)):.2e}")
