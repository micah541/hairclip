"""Hand-off data for a run: scale-invariant checks, S1-chart extent, slow-mode coefficient c3 (psi3(0)=1)."""
import sys, re, glob, numpy as np, scipy.linalg as sl
sys.path.insert(0, '/home/claude/s1spec/step3')
from polyweight import build
from measure import load, strand, graph, xs, K0, D_S1, p1, th1, wq
Xb, h = 34.0, 0.05
xg, B, A, r, ps1, ps2 = build(0.25, Xb, h, kind='gauss')
ev, vr = sl.eig(B.T @ A @ B)
i = int(np.argmin(abs(ev + 0.2563))); n = int(round(Xb / h))
psi3 = ((B @ np.real(vr[:, i])) / r)[n:]; psi3 /= psi3[0]
psi3f = np.interp(xs, np.arange(n + 1) * h, psi3)
tag = sys.argv[1]
files = sorted([f for f in glob.glob(f"snap_{tag}_*.npz") if re.fullmatch(rf"snap_{re.escape(tag)}_\d+\.npz", f)],
               key=lambda f: int(f.split('_')[-1][:-4]))
print(f"{tag}:  kmax  (c2-tip)k  |X|k^2  chart_end  c1_resc   c3      c3*k   resid[0,6]/max|u|   |u|_rho(chart)")
rows = []
for f in files:
    S, z = load(f)
    al, th = S.theta(8); N = len(al)
    kt = abs((2 * np.pi / S.L) * np.gradient(th, al)[N // 2])
    Zs, als = strand(S); tip = Zs[0].real
    lam0 = K0 / kt; xc0 = tip + D_S1 * lam0
    u = graph(Zs, xc0, lam0, 34.0)
    ok = np.isfinite(u); m = len(u)
    if ok.sum() < 50: continue
    xend = xs[:m][ok][-1]
    sel = ok & (xs[:m] <= min(6.0, xend))
    M = np.column_stack([psi3f[:m][sel], p1[:m][sel], -np.sin(th1[:m][sel])])
    c, *_ = np.linalg.lstsq(M, u[sel], rcond=None)
    res = np.max(np.abs(u[sel] - M @ c)) / np.max(np.abs(u[sel]))
    nr = np.sqrt(np.sum(wq[:m][ok] * u[ok]**2))
    row = (float(z['kmax']), (float(z['c2']) - tip) * kt, abs(float(z['AX'])) * kt**2, xend, (float(z['c1']) - xc0) / lam0, c[0], c[0] * kt, res, nr)
    rows.append(row)
    print(f"      {row[0]:7.1f}  {row[1]:7.4f}  {row[2]:7.4f}  {row[3]:6.2f}   {row[4]:7.2f}  {row[5]:+.4f}  {row[6]:+6.1f}    {row[7]:.1e}           {row[8]:.3f}")
print("S1:            2.5720   3.9760")
np.save(f"handoff_{tag}.npy", np.array(rows))
