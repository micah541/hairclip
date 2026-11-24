import numpy as np
from scipy.integrate import solve_ivp
# 1) coupling matrix C = V^{-1} E V
w = np.exp(2j*np.pi*np.arange(3)/3)
V = np.array([[(-w[j])**i for j in range(3)] for i in range(3)])
C = np.linalg.solve(V, np.diag([0,1,2]) @ V)
print("C =\n", np.round(C, 6)); print("|C_jl| off-diag:", np.round(abs(C[~np.eye(3,dtype=bool)]),6))
print("predicted C01 =", -0.5-1j*np.sqrt(3)/6, " C10 =", -0.5+1j*np.sqrt(3)/6)
# 2) weighted norms of the true conical end, from an accurate backward integration
S = 8.0; a = 0.75*4**(-1/3); mu = lambda s: (s/4)**(1/3); E = lambda s: a*(s**(4/3)-S**(4/3))
f = lambda s,y: [-y[2]*y[1], 1+y[2]*y[0], y[3], y[0]/4]
def lin_state(s, A):   # pure decaying mode at s
    u = mu(s); lam = -u*w
    k = A; kp = (lam[0]*A).real; kpp = (lam[0]**2*A).real
    return [4*kpp, s, k, kp]
sfar = 24.0
A = 1.0e-15
sol = solve_ivp(f, [sfar, S], lin_state(sfar, A), method='DOP853', rtol=1e-13, atol=1e-300, dense_output=True)
def Z(s):
    p,q,k,m = sol.sol(s); u = mu(s); lam = -u*w
    Vs = np.vstack([np.ones(3), lam, lam**2]); return np.linalg.solve(Vs, np.array([k, m, p/4]))
eps = Z(S)[0].real
ss = np.linspace(S, 20, 2000)
wgt = lambda s: np.exp(E(s))*(s/S)**(1/3); g = lambda s: (s/S)**(4/3)
z0 = max(abs(Z(s)[0].real - eps*np.exp(-E(s))*(S/s)**(1/3))*wgt(s) for s in ss)
z1 = max(abs(Z(s)[1])*wgt(s)*g(s) for s in ss)
print(f"eps={eps:.4e}: sup w|Z0 - eps Phi0| = {z0:.3e} (bound r0 for eps_max=1.7e-4: 6.35e-7, scaled to this eps: {6.352e-7*eps/1.7e-4:.3e})")
print(f"            sup w g |Z1|        = {z1:.3e} (bound r1 scaled: {2.2e-6*eps/1.7e-4:.3e})")
print(f"            Z1(S)/eps = {Z(S)[1]/eps:.6f},  gamma1 = -0.010356+0.000333j,  |diff| = {abs(Z(S)[1]/eps-(-0.010356215835809902+0.000332793623115475j)):.2e} (bound delta1/eps = 2.12e-4)")
