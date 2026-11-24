"""Numerical: coercivity of L in W=(1+x^2)^(-p) for functions supported in [R1, X] (the cone region).
Prediction: abscissa -> (3-2p)/8 as R1 grows (corrections O(p^2/R1^2))."""
import sys, numpy as np, scipy.linalg as sl
sys.path.insert(0, '/home/claude/s1spec')
from compress import profile_on
X, h = 200.0, 0.1
for pw in [1.75, 2.0, 2.5]:
    out = []
    for R1 in [3, 5, 8, 12, 20]:
        x = np.arange(R1, X + h/2, h); N = len(x)
        p, q, k, m = profile_on(x)
        a2 = -k*k; a1 = -(5*k*m + q/4); a0 = 0.25 - 3*m*m - k*p
        I = np.eye(N)
        D1 = (np.eye(N, k=1) - np.eye(N, k=-1))/(2*h)
        D2 = (np.eye(N, k=1) - 2*I + np.eye(N, k=-1))/h**2
        D4 = (np.eye(N, k=2) - 4*np.eye(N, k=1) + 6*I - 4*np.eye(N, k=-1) + np.eye(N, k=-2))/h**4
        L = -D4 + a2[:, None]*D2 + a1[:, None]*D1 + np.diag(a0)
        r = (1 + x*x)**(-pw/2)
        A = (r[:, None]*L)/r[None, :]
        H = (A + A.T)/2
        out.append(sl.eigh(H, eigvals_only=True, subset_by_index=[N-1, N-1])[0])
    print(f"p={pw}: pred {(3-2*pw)/8:+.4f}   R1=3,5,8,12,20 -> {np.round(out, 4)}", flush=True)
