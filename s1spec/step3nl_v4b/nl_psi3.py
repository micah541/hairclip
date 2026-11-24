"""basin test along the slow mode psi3 (lambda~ = -0.256), whose tail grows like |x|^2.02"""
import sys, numpy as np, scipy.linalg as sl
sys.path.insert(0, '/home/claude/s1spec/step3')
from polyweight import build
from nlflow import Flow, run
X, h = 34.0, 0.05
xg, B, A, r, psi1, psi2 = build(0.25, X, h, kind='gauss')
Ae = B.T @ A @ B
ev, vr = sl.eig(Ae)
i = int(np.argmin(abs(ev + 0.2563))); print("eig:", ev[i])
n = int(round(X / h))
vec = B @ np.real(vr[:, i])            # in r-scaled coords on the full grid
u = vec / r                             # actual eigenfunction
half = u[n:]                            # x >= 0
half = half / half[0]
fl = Flow(X=X, h=h); x = fl.x
for xx in [0, 2, 5, 10, 15, 20]:
    print(f"psi3({xx}) = {half[int(round(xx/h))]:+.4e}")
for Rc in [12.0, 20.0]:
    cut = np.where(x < Rc - 3, 1.0, np.where(x < Rc, np.cos(np.pi / 6 * (x - Rc + 3))**2, 0.0))
    for d in [float(a) for a in sys.argv[1].split(',')]:
        H, v = run(fl, d * half * cut, dt=2e-3, T=15.0, every=1.0, verbose=False)
        last = H[-1]
        ok = np.isfinite(last[1]) and last[0] >= 15 - 1e-9
        print(f"cut {Rc}: delta={d:+.2f} |u0|={H[0,1]:.3e} sup0(x<10)={H[0,2]:.3e} -> tau={last[0]:.1f} |u|={last[1]:.3e} {'ok' if ok else 'FAIL'}", flush=True)
