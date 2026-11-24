import sympy as sp
X = sp.symbols('X', real=True)
phi = sp.Rational(1, 4)*(1 + X**2)**sp.Rational(2, 3)
E = (1 + X**2)**sp.Rational(-1, 3)
claims = {1: X*E/3, 2: E/3 - sp.Rational(2, 9)*X**2*E**4, 3: -sp.Rational(2, 3)*X*E**4 + sp.Rational(16, 27)*X**3*E**7,
          4: -sp.Rational(2, 3)*E**4 + sp.Rational(32, 9)*X**2*E**7 - sp.Rational(224, 81)*X**4*E**10}
for j, c in claims.items():
    print(f"P{j} check:", sp.simplify(sp.diff(phi, X, j) - c))
print("E' = -(2/3) X E^4:", sp.simplify(sp.diff(E, X) + sp.Rational(2, 3)*X*E**4))
# transpose eigen-ODE: L^T w = -w'''' - (k^2 w)'' + (b w)' + c w, b = 5km + q/4, c = 1/4 - 3m^2 - kp
x = sp.symbols('x'); w = sp.Function('w')(x)
k, m, p, q = [sp.Function(n)(x) for n in 'kmpq']
prof = {sp.Derivative(k, x): m, sp.Derivative(m, x): p/4, sp.Derivative(p, x): -k*q, sp.Derivative(q, x): 1 + k*p}
def D(e, n=1):
    for _ in range(n): e = sp.diff(e, x).subs(prof)
    return e
lam = sp.symbols('lam')
b = 5*k*m + q/4; c = sp.Rational(1, 4) - 3*m**2 - k*p
LT = -sp.diff(w, x, 4) - D(k**2*w, 2) + D(b*w) + c*w
claim = -sp.diff(w, x, 4) - k**2*sp.diff(w, x, 2) + (k*m + q/4)*sp.diff(w, x) + sp.Rational(1, 2)*w
print("L^T w = -w'''' - k^2 w'' + (km + q/4) w' + w/2 :", sp.simplify(sp.expand(LT - claim)))
# and L^T is the transpose of L: int w L u = int u L^T w for compactly supported
u = sp.Function('u')(x)
Lu = -sp.diff(u, x, 4) - k**2*sp.diff(u, x, 2) - b*sp.diff(u, x) + c*u
expr = w*Lu - u*LT
# exact derivative test via variational derivatives in u and w
def vd(E_, f):
    out = 0
    for j in range(0, 5):
        dE = sp.diff(E_, sp.diff(f, x, j)) if j > 0 else sp.diff(E_, f)
        out += (-1)**j*D(dE, j) if j > 0 else dE
    return sp.simplify(sp.expand(out.subs(prof)))
print("w Lu - u L^T w exact derivative:", vd(expr, u), vd(expr, w))
