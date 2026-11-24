"""scale-invariant checks against the notes' table: (c2 - tip) k_tip, |X| k_tip^2, and the tip-fit deviation."""
import sys, glob, numpy as np
from measure import load, strand, graph, xs, wq, K0, D_S1, HX
tag = sys.argv[1]
files = sorted(glob.glob(f"snap_{tag}_*.npz"), key=lambda f: int(f.split('_')[-1][:-4]))
print(" kmax    ktip    (c2-tip)k   |X|k^2   chart_end  sup|u|(x<=2)  sup|u|(x<=chart)  |u|_rho(chart)")
for f in files:
    S, z = load(f)
    al, th = S.theta(8); N = len(al)
    kt = abs((2 * np.pi / S.L) * np.gradient(th, al)[N // 2])
    Zs, als = strand(S)
    tip = Zs[0].real
    lam0 = K0 / kt; xc0 = tip + D_S1 * lam0
    u = graph(Zs, xc0, lam0, 34.0)
    ok = np.isfinite(u); n = len(u); xx = xs[:n]
    xend = xx[ok][-1] if ok.any() else 0
    mx = lambda a: np.max(np.abs(a)) if a.size else np.nan
    s2 = mx(u[(xx <= 2) & ok]); sall = mx(u[ok])
    nr = np.sqrt(np.sum(wq[:n][ok] * u[ok]**2))
    print(f"{float(z['kmax']):7.1f} {kt:8.2f}  {(float(z['c2'])-tip)*kt:8.4f}  {abs(float(z['AX']))*kt**2:8.4f}  {xend:6.2f}     {s2:.3e}     {sall:.3e}     {nr:.3e}")
print("S1:              2.5720    3.9760")
