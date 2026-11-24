import sympy as sp
s=sp.symbols('s')
p,q,k,m,th=[sp.Function(n)(s) for n in ['p','q','k','m','th']]
# shrinker ODE: p'=-kq, q'=1+kp, k'=m, m'=p/4, th'=k
rules={p.diff(s):-k*q, q.diff(s):1+k*p, k.diff(s):m, m.diff(s):p/4, th.diff(s):k}
def D(e):  # derivative along the ODE
    return sp.expand(e.diff(s).subs(rules))
def Dn(e,n):
    for _ in range(n): e=D(e)
    return e
def L(u):
    kp=m; kpp=p/4
    return sp.simplify(-Dn(u,4) - k**2*Dn(u,2) - (5*k*kp + q/4)*Dn(u,1) + (sp.Rational(1,4) - 3*kp**2 - 4*k*kpp)*u)
print("L q            =", sp.simplify(L(q)))
print("L p - p        =", sp.simplify(L(p)-p))
print("L(-sin th) + sin(th)/4 =", sp.simplify(L(-sp.sin(th)) + sp.sin(th)/4))
print("L(cos th) - cos(th)/4  =", sp.simplify(L(sp.cos(th)) - sp.cos(th)/4))
