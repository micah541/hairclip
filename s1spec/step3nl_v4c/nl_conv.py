import numpy as np
from nlflow import Flow, run
cut = lambda x: np.where(x < 8, 1.0, np.where(x < 10, np.cos(np.pi / 4 * (x - 8))**2, 0.0))
for X, h, dt in [(34.0, 0.05, 2e-3), (34.0, 0.04, 1e-3), (44.0, 0.05, 2e-3)]:
    fl = Flow(X=X, h=h); x = fl.x
    u0 = 0.4 * (1 - (x / 10)**2)**2 * (x < 10) * cut(x)
    H, v = run(fl, u0, dt=dt, T=15.0, every=5.0, verbose=False)
    print(f"X={X} h={h} dt={dt}: " + "  ".join(f"tau={r[0]:.0f}: |u|={r[1]:.4e} loglam={r[3]:+.4f} x0={r[4]:+.4e}" for r in H), flush=True)
