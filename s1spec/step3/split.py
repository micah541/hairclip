"""Quadratic form Q(u) = Re<Lu,u>_rho written as a local form
   Q(u) = int rho [ -|u''|^2 + Acoef |u'|^2 + V |u|^2 ] dx   (no boundary terms on the whole line)
Check: (1) top of Q/|u|^2 on even functions agrees with the operator computation;
       (2) far region [R, X] with FREE end at R;  (3) near region [0,R] (even at 0, free at R)
           compressed against cut-off constraint functions (p, -sin th, psi3 restricted to [0,R])."""
import sys, numpy as np, scipy.linalg as sl
from compress import profile_on, build
from scipy.integrate import cumulative_trapezoid
beta = 0.25
def coefs(x):
    p, q, k, m = profile_on(x)
    s = 1 + x*x
    f1 = (4*beta/3)*x*s**(-1/3)
    f2 = (4*beta/3)*(s**(-1/3) - (2/3)*x*x*s**(-4/3))
    f3 = (4*beta/3)*(-(2/3)*x*s**(-4/3)*3 + (8/9)*x**3*s**(-7/3))
    # f4 by finite difference of f3 (smooth)
    hh = 1e-4
    def F3(z):
        sz = 1 + z*z
        return (4*beta/3)*(-2*z*sz**(-4/3) + (8/9)*z**3*sz**(-7/3))
    f4 = (F3(x+hh) - F3(x-hh))/(2*hh)
    r2 = f1**2 - f2                                  # rho''/rho
    r4 = f1**4 - 6*f1**2*f2 + 3*f2**2 + 4*f1*f3 - f4  # rho''''/rho
    k2 = k*k; k2p = 2*k*m; k2pp = 2*m*m + k*p/2
    b = 5*k*m + q/4; bp = 5*m*m + 5*k*p/4 + (1 + k*p)/4
    c = 0.25 - 3*m*m - k*p
    Acoef = 2*r2 + k2
    V = -0.5*r4 - 0.5*(k2pp - 2*k2p*f1 + k2*r2) + 0.5*(bp - b*f1) + c
    rho = np.exp(-beta*s**(2/3))
    return rho, Acoef, V, p, k
def form(a, bnd, h, even0):
    x = np.arange(a, bnd + h/2, h); n = len(x)
    rho, Ac, V, p, k = coefs(x)
    D1 = np.zeros((n, n)); D2 = np.zeros((n, n))
    for i in range(n):
        if 0 < i < n-1:
            D1[i, i-1], D1[i, i+1] = -1/(2*h), 1/(2*h)
            D2[i, [i-1, i, i+1]] = np.array([1, -2, 1])/h**2
    if even0:   # ghost u_{-1} = u_1
        D1[0, :] = 0; D2[0, [0, 1]] = np.array([-2, 2])/h**2
    else:
        D1[0, [0, 1, 2]] = np.array([-3, 4, -1])/(2*h); D2[0, [0, 1, 2]] = np.array([1, -2, 1])/h**2
    D1[-1, [-3, -2, -1]] = np.array([1, -4, 3])/(2*h); D2[-1, [-3, -2, -1]] = np.array([1, -2, 1])/h**2
    w = rho*h; w[0] *= (0.5 if not even0 else 0.5); w[-1] *= 0.5
    Qm = -(D2.T*w) @ D2 + (D1.T*(w*Ac)) @ D1 + np.diag(w*V)
    M = np.diag(w)
    return x, Qm, M, p, k
def top(Qm, M, C=None, nt=3):
    Mh = np.sqrt(np.diag(M)); H = Qm/np.outer(Mh, Mh)
    if C is not None:
        Cm = C*Mh[:, None]
        Qq, _ = np.linalg.qr(Cm, mode='complete'); Z = Qq[:, C.shape[1]:]
        H = Z.T @ H @ Z
    return sl.eigh(H, eigvals_only=True, subset_by_index=[H.shape[0]-nt, H.shape[0]-1])[::-1]
h = float(sys.argv[1]) if len(sys.argv) > 1 else 0.04
# (1) whole half-line, even
x, Qm, M, p, k = form(0.0, 30.0, h, True)
print("even half-line top of form:", np.round(top(Qm, M, nt=5), 4))
# eigenfunction psi3 from the operator discretization
xf, B, A, L, r, psi1, psi2 = build(beta, 30.0, h)
Ae = B.T @ A @ B
ev, vr = sl.eig(Ae)
i3 = int(np.argmin(abs(ev + 0.2563)))
v3 = np.real(vr[:, i3]); n0 = (len(xf)-1)//2
psi3_half = (B @ v3)[n0:]/r[n0:]          # back to un-conjugated function on x >= 0
psi3_half /= np.max(abs(psi3_half[:int(8/h)]))
theta = -np.pi/2 + cumulative_trapezoid(k, x, initial=0.0)
for R in [6.0, 8.0, 10.0, 12.0]:
    xa, Qa, Ma, pa, ka = form(R, 30.0, h, False)
    far = top(Qa, Ma, nt=2)
    xn, Qn, Mn, pn, kn = form(0.0, R, h, True)
    th = -np.pi/2 + cumulative_trapezoid(kn, xn, initial=0.0)
    C3 = np.column_stack([pn, -np.sin(th), psi3_half[:len(xn)]])
    near2 = top(Qn, Mn, C3[:, :2], nt=2); near3 = top(Qn, Mn, C3, nt=2)
    print(f"R={R:5.1f}  far(free at R) top={np.round(far,4)}   near: perp(p,sin)={np.round(near2,4)}  perp(p,sin,psi3)={np.round(near3,4)}", flush=True)
