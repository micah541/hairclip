"""Near-optimal certificates for the numerical-range bounds of L in L^2(rho), rho = exp(-phi),
phi = beta (1+x^2)^(2/3).  For angle th we seek real functions w11, w22 (odd), wc, wd (even), w12 = wc + i wd,
and Lam such that pointwise on x >= 0 (parity covers x < 0):
    P(x) := H_th(x) + N_W(x) - Lam E33  <=  -tau I        (3x3 Hermitian, jet J = (u'', u', u))
where  Re(e^{i th} <Lu,u>_rho) = int J^* H_th J rho  and  int J^* N_W J rho = 0.
Then  Re(e^{i th} lam) <= Lam for every eigenvalue lam (and for the whole numerical range)."""
import numpy as np, cvxpy as cp, sys, json
sys.path.insert(0, '/home/claude/s1spec')
from evans import PROF, STIP

def profile(x):
    sig = x + STIP
    if sig <= 22.0:
        p, q, k, m = PROF(sig)
    else:
        p, q, k, m = 0.0, sig, 0.0, 0.0
    return p, q, k, m

def weight(x, beta):
    u = 1 + x*x
    phi = beta*u**(2/3)
    phi1 = beta*(4/3)*x*u**(-1/3)
    phi2 = beta*(4/3)*(u**(-1/3) - (2/3)*x*x*u**(-4/3))
    return phi, phi1, phi2

class Basis:
    """w(x) = sum theta_n b_n(x); returns value and derivative matrices on a grid."""
    def __init__(self, nodes, s):
        self.nodes, self.s = nodes, s
    def funcs(self, kind):
        s = self.s; F = []
        g = lambda x, c: np.exp(-(x - c)**2/(2*s*s))
        dg = lambda x, c: -(x - c)/(s*s)*np.exp(-(x - c)**2/(2*s*s))
        if kind == 'odd1':     # w11 ~ x^{1/3}
            F.append((lambda x: x*(1 + x*x)**(-1/3), lambda x: (1 + x*x)**(-1/3) - (2/3)*x*x*(1 + x*x)**(-4/3)))
        if kind == 'odd3':     # w22 ~ x, + x^{1/3}
            F.append((lambda x: x, lambda x: 1.0 + 0*x))
            F.append((lambda x: x*(1 + x*x)**(-1/3), lambda x: (1 + x*x)**(-1/3) - (2/3)*x*x*(1 + x*x)**(-4/3)))
        if kind == 'even2':    # w12 ~ x^{2/3}, + const
            F.append((lambda x: (1 + x*x)**(1/3), lambda x: (2/3)*x*(1 + x*x)**(-2/3)))
            F.append((lambda x: 1.0 + 0*x, lambda x: 0*x))
        for c in self.nodes:
            if kind.startswith('odd'):
                if c == 0: continue
                F.append((lambda x, c=c: g(x, c) - g(x, -c), lambda x, c=c: dg(x, c) - dg(x, -c)))
            else:
                F.append((lambda x, c=c: g(x, c) + g(x, -c), lambda x, c=c: dg(x, c) + dg(x, -c)))
        return F

def build_terms(xs, beta, th):
    """Return H (n,3,3) complex and the data needed for N."""
    n = len(xs); H = np.zeros((n, 3, 3), dtype=complex); PH1 = np.zeros(n)
    e = np.exp(1j*th)
    for i, x in enumerate(xs):
        p, q, k, m = profile(x)
        phi, phi1, phi2 = weight(x, beta)
        g = phi1**2 - phi2
        b = -(5*k*m + q/4); d = 0.25 - 3*m*m - k*p
        Hm = np.zeros((3, 3), dtype=complex)
        Hm[0, 0] = -np.cos(th)
        Hm[1, 0] = e*phi1;            Hm[0, 1] = np.conj(Hm[1, 0])
        Hm[2, 0] = -e*(g + k*k)/2;    Hm[0, 2] = np.conj(Hm[2, 0])
        Hm[2, 1] = e*b/2;             Hm[1, 2] = np.conj(Hm[2, 1])
        Hm[2, 2] = np.cos(th)*d
        H[i] = Hm; PH1[i] = phi1
    return H, PH1

def solve(beta, th, tau=1e-3, X1=26.0, h=0.05, nodes_step=0.5, s=0.6, verbose=False):
    xs = np.concatenate([np.arange(0, X1, h), np.exp(np.linspace(np.log(X1), np.log(1e6), 120))])
    H, PH1 = build_terms(xs, beta, th)
    B = Basis(np.arange(0, X1 - 2, nodes_step), s)
    F11, F22, Fc = B.funcs('odd1'), B.funcs('odd3'), B.funcs('even2')
    def mats(F):
        V = np.array([[f(x) for (f, df) in F] for x in xs]); D = np.array([[df(x) for (f, df) in F] for x in xs])
        return V, D
    V11, D11 = mats(F11); V22, D22 = mats(F22); Vc, Dc = mats(Fc)
    t11 = cp.Variable(len(F11)); t22 = cp.Variable(len(F22)); tc = cp.Variable(len(Fc)); td = cp.Variable(len(Fc))
    Lam = cp.Variable()
    w11 = V11 @ t11; dw11 = D11 @ t11; w22 = V22 @ t22; dw22 = D22 @ t22
    wc = Vc @ tc; dwc = Dc @ tc; wd = Vc @ td; dwd = Dc @ td
    cons = []
    for i in range(len(xs)):
        f1 = PH1[i]
        # real and imaginary parts of P = H + N - Lam E33 + tau I
        Ar = [[None]*3 for _ in range(3)]; Ai = [[None]*3 for _ in range(3)]
        Hr, Hi = H[i].real, H[i].imag
        Nr = [[0, w11[i], wc[i]],
              [w11[i], dw11[i] + 2*wc[i] - f1*w11[i], dwc[i] + w22[i] - f1*wc[i]],
              [wc[i], dwc[i] + w22[i] - f1*wc[i], dw22[i] - f1*w22[i]]]
        Ni = [[0, 0, wd[i]],
              [0, 0, dwd[i] - f1*wd[i]],
              [-wd[i], -(dwd[i] - f1*wd[i]), 0]]
        for a in range(3):
            for bb in range(3):
                Ar[a][bb] = Hr[a, bb] + Nr[a][bb] + (tau if a == bb else 0) - (Lam if a == bb == 2 else 0)
                Ai[a][bb] = Hi[a, bb] + Ni[a][bb]
        R = cp.bmat([[Ar[a][bb] for bb in range(3)] + [-Ai[a][bb] for bb in range(3)] for a in range(3)] +
                    [[Ai[a][bb] for bb in range(3)] + [Ar[a][bb] for bb in range(3)] for a in range(3)])
        cons.append((R + R.T)/2 << 0)
    prob = cp.Problem(cp.Minimize(Lam), cons)
    prob.solve(solver=cp.CLARABEL, verbose=verbose)
    sol = dict(beta=beta, th=th, tau=tau, X1=X1, nodes=list(B.nodes), s=s, Lam=float(Lam.value),
               t11=list(map(float, t11.value)), t22=list(map(float, t22.value)),
               tc=list(map(float, tc.value)), td=list(map(float, td.value)), status=prob.status)
    return sol

if __name__ == "__main__":
    beta = float(sys.argv[1]); th = float(sys.argv[2]); tau = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-3
    sol = solve(beta, th, tau=tau)
    print(f"beta={beta} th={th} tau={tau}: status={sol['status']}  Lam={sol['Lam']:.5f}")
    json.dump(sol, open(f"cert_b{beta}_th{th}_tau{tau}.json", "w"))
