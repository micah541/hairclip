"""Summarize shape-search runs: outcome, T, lens margin, slow-mode amplitude c3*k at k ~ 200 and 400."""
import sys, re, glob, numpy as np, scipy.linalg as sl
sys.path.insert(0, '/home/claude/s1spec/step3')
from polyweight import build
from measure import load, strand, graph, xs, K0, D_S1, p1, th1
Xb, h = 34.0, 0.05
xg, B, A, r, ps1, ps2 = build(0.25, Xb, h, kind='gauss')
ev, vr = sl.eig(B.T @ A @ B)
i = int(np.argmin(abs(ev + 0.2563))); n = int(round(Xb / h))
psi3 = ((B @ np.real(vr[:, i])) / r)[n:]; psi3 /= psi3[0]
psi3f = np.interp(xs, np.arange(n + 1) * h, psi3)

def c3k(f):
    S, z = load(f)
    al, th = S.theta(8); N = len(al)
    kt = abs((2 * np.pi / S.L) * np.gradient(th, al)[N // 2])
    Zs, als = strand(S); tip = Zs[0].real
    lam0 = K0 / kt; xc0 = tip + D_S1 * lam0
    u = graph(Zs, xc0, lam0, 34.0)
    ok = np.isfinite(u); m = len(u)
    if ok.sum() < 50: return np.nan, 0.0
    xend = xs[:m][ok][-1]
    sel = ok & (xs[:m] <= min(6.0, xend))
    M = np.column_stack([psi3f[:m][sel], p1[:m][sel], -np.sin(th1[:m][sel])])
    c, *_ = np.linalg.lstsq(M, u[sel], rcond=None)
    return c[0] * kt, xend

tags = sys.argv[1:] if len(sys.argv) > 1 else sorted(set(re.sub(r"^hp_|\.log$", "", f) for f in glob.glob("hp_[srqk]*_*.log")))
print(f"{'run':22s} {'outcome':9s} {'T':>10s} {'x_T':>8s} {'|X|0':>6s} {'|Lam|(T)':>9s} {'c3k@200':>8s} {'c3k@400':>8s} {'chart@400':>9s}")
for tag in tags:
    rows = []; first = None; done = False
    for line in open(f"hp_{tag}.log"):
        if line.startswith("DONE"): done = True
        if 'kmax' not in line: continue
        t = float(re.search(r"t=([\d.e+-]+)", line).group(1)); k = float(re.search(r"kmax=\s*([\d.e+-]+)", line).group(1))
        lam = re.search(r"\|Lam\|=([\d.e+-]+)", line); X = re.search(r"\|X\|=([\d.e+-]+)", line)
        tip = float(re.search(r"tipx=([\d.e+-]+)", line).group(1))
        rows.append((t, k, float(lam.group(1)) if lam else np.nan, float(X.group(1)) if X else np.nan, tip))
        if X and first is None: first = float(X.group(1))
    rows = np.array(rows)
    if not done:
        print(f"{tag:22s} running (t={rows[-1,0]:.4f}, kmax={rows[-1,1]:.1f})"); continue
    late = np.isfinite(rows[:, 3]) & (rows[:, 3] < 0.3 * (first if first else 1)); sel = (rows[:, 1] > 60) & (rows[:, 1] < 450) & late
    blow = rows[-1, 1] > 400 and np.isfinite(rows[-1, 3]) and rows[-1, 3] < 1e-3
    if blow:
        Af = np.polyfit(rows[sel, 0], rows[sel, 1]**-4, 1); T = -Af[1] / Af[0]
        snaps = {int(f.split('_')[-1][:-4]): f for f in glob.glob(f"snap_{tag}_*.npz") if re.fullmatch(rf"snap_{re.escape(tag)}_\d+\.npz", f)}
        c2, _ = c3k(snaps[200]) if 200 in snaps else (np.nan, 0)
        c4, xe = c3k(snaps[400]) if 400 in snaps else (np.nan, 0)
        print(f"{tag:22s} {'blow-up':9s} {T:10.5f} {rows[-1,4]:8.3f} {first:6.3f} {rows[-1,2]:9.4f} {c2:8.2f} {c4:8.2f} {xe:9.1f}")
    else:
        last2 = rows[np.isfinite(rows[:, 3])][-1]
        print(f"{tag:22s} {'re-embeds':9s} {'t=%.4f' % last2[0]:>10s} {'':8s} {first if first else float('nan'):6.3f} {'X=%.4f' % last2[3]:>9s}")
