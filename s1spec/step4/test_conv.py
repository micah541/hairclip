import numpy as np
from hls import HLS, from_points
th, L, X0 = from_points(lambda p: 2*np.cos(p) + 1j*np.sin(p), 1, 512)
res = {}
for h in [8e-6, 4e-6, 2e-6, 1e-6]:          # tau' steps; (2pi/L)^4 ~ 1.7e-1 ... dt = h (L/2pi)^4
    s = HLS(th, L, X0, omega=1, filt=False, tail_tol=1.0)
    A0 = s.area()
    n = int(round(2e-3 / h))
    for i in range(n): s.step(h)
    res[h] = (s.theta()[1], s.L, s.area() / A0 - 1, s.t)
    print(f"h'={h}: t={s.t:.6f} L={s.L:.12f} area err {res[h][2]:+.3e}", flush=True)
ref = res[1e-6][0]
for h in [8e-6, 4e-6, 2e-6]:
    print(f"h'={h}: max|theta - theta_ref| = {np.max(np.abs(res[h][0] - ref)):.3e}  L-Lref={res[h][1]-res[1e-6][1]:+.3e}")
