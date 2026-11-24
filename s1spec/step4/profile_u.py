import numpy as np
from measure import load, strand, graph, xs, K0, D_S1
for kk in [100, 400, 800]:
    f = f"snap_fine_{kk}.npz" if kk != 800 else "snap_fine_800.npz"
    S, z = load(f)
    al, th = S.theta(8); N = len(al)
    kt = abs((2 * np.pi / S.L) * np.gradient(th, al)[N // 2])
    Zs, als = strand(S); tip = Zs[0].real
    lam0 = K0 / kt; xc0 = tip + D_S1 * lam0
    u = graph(Zs, xc0, lam0, 34.0)
    vals = []
    for x in [0.5, 1, 1.5, 2, 3, 4, 5, 6, 8, 10, 12, 14, 16]:
        i = int(round(x / 0.01))
        vals.append(f"{u[i]:+.3f}" if i < len(u) and np.isfinite(u[i]) else "  --  ")
    print(f"k={kk:4d} lam={lam0:.3e}: u at x=0.5,1,1.5,2,3,4,5,6,8,10,12,14,16:", " ".join(vals))
