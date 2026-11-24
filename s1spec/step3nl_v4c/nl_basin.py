"""delta-sweep: does the modulated rescaled flow return to S1?"""
import numpy as np, time, sys
from nlflow import Flow, run
fl = Flow(X=34.0, h=0.05)
x = fl.x
cut = lambda x: np.where(x < 8, 1.0, np.where(x < 10, np.cos(np.pi / 4 * (x - 8))**2, 0.0))
shapes = {
    'tip':  np.exp(-x**2 / 4),
    'loop': np.exp(-(x - 4.0)**2 / 2) + np.exp(-(x + 4.0)**2 / 2),
    'wide': (1 - (x / 10)**2)**2 * (x < 10),
}
for name in sys.argv[1].split(','):
    for d in [float(a) for a in sys.argv[2].split(',')]:
        u0 = d * shapes[name] * cut(x)
        H, v = run(fl, u0, dt=2e-3, T=15.0, every=1.0, verbose=False)
        last = H[-1]
        status = 'ok' if np.isfinite(last[1]) and last[0] >= 15 - 1e-9 else f'fail at tau={last[0]}'
        print(f"{name:5s} delta={d:6.3f}  |u0|_rho={H[0,1]:.3e} sup0={H[0,2]:.3e}  ->  tau={last[0]:5.2f} |u|_rho={last[1]:.3e} "
              f"sup={last[2]:.3e} max sup={H[:,2].max():.3e}  {status}", flush=True)
