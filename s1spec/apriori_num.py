"""Non-rigorous evaluation of the a priori (energy) bound potentials along S1.
Re lam ||u||^2 <= -(1-th)||u''||^2 + int U_th |u|^2 rho,
U_th = V + (1/2)(G rho)''/rho + G^2/(4 th),  G = 2 g + k^2,  g = phi'^2 - phi'',
V = -(1/2)(g rho)''/rho + (1/2)(k^2 phi' rho)'/rho - (1/2)(e' - e phi') + d,
e = -3 k k' - q/4,  d = 1/4 - 3 k'^2 - k p,  rho = exp(-phi),  phi = beta (1+s^2)^(2/3), s = sigma - sigma_tip."""
import numpy as np, sympy as sp
from evans import PROF, STIP
s, b = sp.symbols('s b', real=True)
phi = b*(1+s**2)**sp.Rational(2,3)
dphi = [sp.diff(phi, s, n) for n in range(6)]
fd = [sp.lambdify((s,b), d, 'numpy') for d in dphi]
def potentials(sig, beta, th):
    ss = sig - STIP
    p,q,k,m = PROF(sig) if sig <= 22 else (0.0, sig, 0.0, 0.0)
    kpp = p/4; kppp = -k*q/4                     # k'' = p/4, k''' = p'/4 = -kq/4
    P1,P2,P3,P4 = [fd[n](ss,beta) for n in (1,2,3,4)]
    g = P1**2 - P2; g1 = 2*P1*P2 - P3; g2 = 2*P2**2 + 2*P1*P3 - P4
    G = 2*g + k*k; G1 = 2*g1 + 2*k*m; G2 = 2*g2 + 2*(m*m + k*kpp)
    e = -3*k*m - q/4; e1 = -3*(m*m + k*kpp) - (1 + k*p)/4
    d = 0.25 - 3*m*m - k*p
    rho_term = lambda f, f1, f2: f2 - 2*f1*P1 + f*(P1**2 - P2)       # (f rho)''/rho
    kk = k*k; kk1 = 2*k*m
    V = -0.5*rho_term(g,g1,g2) + 0.5*(kk1*P1 + kk*P2 - kk*P1**2) - 0.5*(e1 - e*P1) + d
    U = V + 0.5*rho_term(G,G1,G2) + G*G/(4*th)
    return U, V, G
for beta in [0.05, 0.1, 0.15, 0.2, 0.25]:
    for th in [0.5, 1.0]:
        sig = np.linspace(STIP, 200, 20000)
        vals = np.array([potentials(x, beta, th)[0] for x in sig])
        i = np.argmax(vals)
        print(f"beta={beta:.2f} theta={th}: sup U = {vals.max():.4f} at s={sig[i]-STIP:.3f};  U(s=50)={potentials(STIP+50,beta,th)[0]:.3f}, U(s=190)={potentials(STIP+190,beta,th)[0]:.2f}")
