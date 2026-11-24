"""Numerical (non-rigorous): numerical abscissa of the modulated linearization in a POLYNOMIAL weight
W = (1+x^2)^(-p).  Prediction: the essential spectrum is Re(lambda) <= (3-2p)/8, and for p < 2.525 the
lambda~ = -0.256 mode (growth |x|^2.02) is not an eigenvalue.  Question: is the form coercive on the
complement of the two symmetry modes alone?"""
import sys, numpy as np, scipy.linalg as sl
sys.path.insert(0, '/home/claude/s1spec')
from compress import profile_on, top_compressed
from scipy.integrate import cumulative_trapezoid

def build(pw, X, h, kind='poly'):
    n = int(round(X/h)); x = h*np.arange(-n, n+1); N = len(x)
    p, q, k, m = profile_on(x)
    a2 = -k*k; a1 = -(5*k*m + q/4); a0 = 0.25 - 3*m*m - k*p
    I = np.eye(N)
    D1 = (np.eye(N, k=1) - np.eye(N, k=-1))/(2*h)
    D2 = (np.eye(N, k=1) - 2*I + np.eye(N, k=-1))/h**2
    D4 = (np.eye(N, k=2) - 4*np.eye(N, k=1) + 6*I - 4*np.eye(N, k=-1) + np.eye(N, k=-2))/h**4
    L = -D4 + a2[:, None]*D2 + a1[:, None]*D1 + np.diag(a0)
    r = (1 + x*x)**(-pw/2) if kind == 'poly' else np.exp(-pw*(1+x*x)**(2/3)/2)
    A = (r[:, None]*L)/r[None, :]
    B = np.zeros((N, n+1)); B[n, 0] = 1.0
    for j in range(1, n+1):
        B[n+j, j] = B[n-j, j] = 1/np.sqrt(2)
    th = cumulative_trapezoid(k, x, initial=0.0); th = -np.pi/2 + th - th[n]
    return x, B, A, r, p, -np.sin(th)

if __name__ == "__main__":
    X = float(sys.argv[1]); h = float(sys.argv[2])
    pws = [float(a) for a in sys.argv[3:]]
    for pw in pws:
        x, B, A, r, psi1, psi2 = build(pw, X, h)
        Ae = B.T @ A @ B; He = (Ae + Ae.T)/2
        w0 = sl.eigh(He, eigvals_only=True, subset_by_index=[He.shape[0]-1, He.shape[0]-1])[0]
        C2 = B.T @ np.column_stack([r*psi1, r*psi2])
        wc = top_compressed(He, C2)
        ev, vl, vr = sl.eig(Ae, left=True, right=True)
        idx = [int(np.argmin(abs(ev - 1.0))), int(np.argmin(abs(ev - 0.25)))]
        ws = top_compressed(He, np.real(vl[:, idx]))
        top = ev[np.argsort(-ev.real)][:6]
        print(f"p={pw:.3f} X={X} h={h}: pred_ess={(3-2*pw)/8:+.4f}  top eig={np.round(top,4)}\n"
              f"      w0={w0:.4f}  orth(psi1,psi2)={np.round(wc,4)}  spec-compl={np.round(ws,4)}", flush=True)
