"""Float cross-check of the IBP identities (task item 4).
Direct side uses F straight from graph_struct.pkl (numeric Ji = 1/sqrt(al^2+U1^2)), N(w) = F(a+w)-F(a)-L0 w.
IBP side: int w0*Lam  (linear part) and int QR (quadratic part) from ibp_sym.
Profile: ODE k'=m, m'=p/4, p'=-kq, q'=1+kp from the S1 tip data (any ODE solution would do).
u0, w even smooth functions, Psi even smooth decaying (so no boundary terms: w'(0)=w'''(0)=0)."""
import numpy as np, sympy as sp, sys, time, json
from scipy.integrate import solve_ivp
import ibp_sym as S
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
# profile
d = np.load('/home/claude/s1spec/step3nl/agent_mod/prof_float.npz'); k0, p0 = d['k'][0], d['p'][0]
sol = solve_ivp(lambda x, y: [y[1], y[3] / 4, -y[0] * y[2 + 0 + 1 - 1 + 0] if False else 0, 0], [0, 1], [0, 0, 0, 0]) if False else None
def rhs(x, y):
    kk, mm, pp, qq = y; return [mm, pp / 4, -kk * qq, 1 + kk * pp]
LMAX = 14.0
sol = solve_ivp(rhs, [0, LMAX], [k0, 0.0, p0, 0.0], method='DOP853', rtol=1e-13, atol=1e-14, dense_output=True)
X = sp.Symbol('X')
def jets_fn(expr, n):          # even smooth analytic function -> jets 0..n-1 as numpy callables
    return [sp.lambdify(X, sp.diff(expr, X, i), 'numpy') for i in range(n)]
ALPHA = 0.3
c = rng.normal(size=4) * 0.5
u0e = ALPHA * (X**2 / (1 + X**2 / 7)) * sp.exp(-X**2 / 30) * (1 + c[0] * sp.cos(X / 2) * X**2 / 5)
we = (c[1] + c[2] * X**2 + c[3] * X**4) * sp.exp(-X**2 / 8) * 0.05 + 0.02 * sp.cos(X) * sp.exp(-X**2 / 12) * (1 + X**2)
psie = sp.exp(-X**2 / 6) * (1 + 0.3 * X**2 + 0.05 * X**4)
au = jets_fn(u0e, 6); wu = jets_fn(we, 6); pu = jets_fn(psie, 6)
# quadrature
npan, ng = 300, 20
gx, gw = np.polynomial.legendre.leggauss(ng)
edges = np.linspace(0, LMAX, npan + 1)
xs = np.concatenate([(edges[i] + edges[i + 1]) / 2 + (edges[i + 1] - edges[i]) / 2 * gx for i in range(npan)])
ws_ = np.concatenate([(edges[i + 1] - edges[i]) / 2 * gw for i in range(npan)])
prof = sol.sol(xs); K, M, P, Q = prof
A = [f(xs) for f in au]; W = [f(xs) for f in wu]; Yv = [f(xs) for f in pu]
# make u0 and w even-extension smooth jets: ok analytic.  Lambdify QR, Lam, and the direct F
args = [*S.Y[:6], S.k, S.m, S.p, S.q, *S.a[:7], *S.w[:6], S.alf, S.Jf, S.alb, S.Jb]
fQR = sp.lambdify(args, S.QR, 'numpy', cse=True); fLam = sp.lambdify(args, S.Lam, 'numpy', cse=True)
fRd = sp.lambdify(args, S.Rdir, 'numpy', cse=True)
# direct F(a+w) - F(a) - L0 w
fFf = sp.lambdify(args, S.Ff, 'numpy', cse=True); fFb = sp.lambdify(args, S.Fb, 'numpy', cse=True)
alf = 1 - K * (A[0] + W[0]); Jf = 1 / np.sqrt(alf**2 + (A[1] + W[1])**2)
alb = 1 - K * A[0]; Jb = 1 / np.sqrt(alb**2 + A[1]**2)
z = np.zeros_like(xs)
val = [*Yv[:6], K, M, P, Q, *A[:6], z, *W[:6], alf, Jf, alb, Jb]
L0w = sum(c_.subs(S.aL and {}) if False else 0 for c_ in [])
aLn = [(-S.k**2), -(5 * S.k * S.m + S.q / 4), sp.Rational(1, 4) - 3 * S.m**2 - S.k * S.p]
L0 = (-W[4] - K**2 * W[2] - (5 * K * M + Q / 4) * W[1] + (0.25 - 3 * M**2 - K * P) * W[0])
Nw = fFf(*val) - fFb(*val) - L0
# linear / remainder split directly
lin_dir = sum(Yv[0] * 0 for _ in [0])
dCf = [sp.lambdify(args, S.dC[i], 'numpy', cse=True) for i in range(5)]
lin_integrand = Yv[0] * sum(dCf[i](*val) * W[i] for i in range(5))
Rint = fRd(*val)
# check N = Lin + R pointwise
print("pointwise |N - Lin - R| max:", np.max(np.abs(Nw - (lin_integrand / np.where(Yv[0] == 0, 1, Yv[0]) if False else sum(dCf[i](*val) * W[i] for i in range(5)) + Rint))))
I = lambda f: float(np.sum(ws_ * f))
Ilin_dir = I(lin_integrand); Ilin_ibp = I(W[0] * fLam(*val))
Iq_dir = I(Yv[0] * Rint); Iq_ibp = I(fQR(*val))
Ind = I(Yv[0] * Nw)
print(f"linear : direct {Ilin_dir:.15e}  ibp {Ilin_ibp:.15e}  diff {abs(Ilin_dir-Ilin_ibp):.2e}  rel {abs(Ilin_dir-Ilin_ibp)/abs(Ilin_dir):.2e}")
print(f"quad   : direct {Iq_dir:.15e}  ibp {Iq_ibp:.15e}  diff {abs(Iq_dir-Iq_ibp):.2e}  rel {abs(Iq_dir-Iq_ibp)/abs(Iq_dir):.2e}")
print(f"total  : direct <Psi,N> {Ind:.15e}  ibp {Ilin_ibp+Iq_ibp:.15e}  diff {abs(Ind-Ilin_ibp-Iq_ibp):.2e}")
print("w2 cubic coefficient of QR (should be identically 0): ", S.sp.expand(S.sp.diff(S.QR, S.w[2], 3)) == 0)
json.dump(dict(seed=int(sys.argv[1]) if len(sys.argv) > 1 else 1, lin_direct=Ilin_dir, lin_ibp=Ilin_ibp, quad_direct=Iq_dir, quad_ibp=Iq_ibp, total_direct=Ind,
               total_ibp=Ilin_ibp + Iq_ibp, abs_diff_lin=abs(Ilin_dir - Ilin_ibp), abs_diff_quad=abs(Iq_dir - Iq_ibp), abs_diff_total=abs(Ind - Ilin_ibp - Iq_ibp)),
          open(f'float_check_{sys.argv[1] if len(sys.argv)>1 else 1}.json', 'w'), indent=1)
