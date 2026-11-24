import numpy as np, sys
from nlflow import Flow, run
fl = Flow(X=34.0, h=0.05); x = fl.x
cut = lambda x: np.where(x < 8, 1.0, np.where(x < 10, np.cos(np.pi / 4 * (x - 8))**2, 0.0))
d = float(sys.argv[1]); dt = float(sys.argv[2])
u0 = d * np.exp(-x**2 / 4) * cut(x)
H, v = run(fl, u0, dt=dt, T=float(sys.argv[3]), every=0.25, verbose=False)
for row in H: print(f"tau={row[0]:6.2f} |u|={row[1]:.3e} sup(x<10)={row[2]:.3e} loglam={row[3]:+.3f}")
