import numpy as np
from hls import HLS, from_points
th, L, X0 = from_points(lambda p: 2*np.cos(p) + 1j*np.sin(p), 1, 256)
for cdt, filt, tol in [(2e-2, True, 1e-11), (5e-3, True, 1e-11), (2e-2, False, 1e-11), (2e-2, True, 1e-14)]:
    s = HLS(th, L, X0, omega=1, cdt=cdt, filt=filt, tail_tol=tol)
    A0 = s.area(); L0 = s.L
    s.run_to(tend=0.5)
    print(f"cdt={cdt} filt={filt} tol={tol}: N={s.N} rel area change {s.area()/A0-1:+.2e}  L {L0:.6f}->{s.L:.6f}", flush=True)
