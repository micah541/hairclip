"""Is the hairpin's deviation from S1 the slow mode?  Fit u ~ c3 psi3 + a1 psi1 + a2 psi2 on [0, X]."""
import sys, numpy as np, scipy.linalg as sl
sys.path.insert(0, '/home/claude/s1spec/step3')
from polyweight import build
from measure import load, strand, graph, xs, K0, D_S1, p1, th1
Xb, h = 34.0, 0.05
xg, B, A, r, ps1, ps2 = build(0.25, Xb, h, kind='gauss')
ev, vr = sl.eig(B.T @ A @ B)
i = int(np.argmin(abs(ev + 0.2563))); n = int(round(Xb / h))
psi3 = ((B @ np.real(vr[:, i])) / r)[n:]; psi3 /= psi3[0]
psi3f = np.interp(xs, np.arange(n + 1) * h, psi3)
for kk in [100, 200, 400, 800]:
    S, z = load(f"snap_fine_{kk}.npz")
    al, th = S.theta(8); N = len(al)
    kt = abs((2 * np.pi / S.L) * np.gradient(th, al)[N // 2])
    Zs, als = strand(S); tip = Zs[0].real
    lam0 = K0 / kt; xc0 = tip + D_S1 * lam0
    u = graph(Zs, xc0, lam0, 34.0)
    ok = np.isfinite(u); m = len(u)
    for X in [6.0, 10.0]:
        sel = ok & (xs[:m] <= X)
        M = np.column_stack([psi3f[:m][sel], p1[:m][sel], -np.sin(th1[:m][sel])])
        c, *_ = np.linalg.lstsq(M, u[sel], rcond=None)
        res = u[sel] - M @ c
        print(f"k={kk:4d} fit on [0,{X:.0f}] (chart to {xs[:m][ok][-1]:.1f}): c3={c[0]:+.4f} a1={c[1]:+.4f} a2={c[2]:+.4f}  "
              f"max|u|={np.max(np.abs(u[sel])):.3f} max|resid|={np.max(np.abs(res)):.4f}  c3*k={c[0]*kk:+.1f}")
