import numpy as np
from hls import HLS, from_points
th, L, X0 = from_points(lambda p: 2*np.cos(p) + 1j*np.sin(p), 1, 512)
for cdt in [2e-2, 1e-2, 5e-3, 2.5e-3]:
    s = HLS(th, L, X0, omega=1, cdt=cdt, tail_tol=1.0)
    A0 = s.area(); out = []
    for te in [0.001, 0.01, 0.1, 0.5]:
        s.run_to(tend=te); out.append(s.area() / A0 - 1)
    print(f"cdt={cdt}: area err at t=1e-3,1e-2,0.1,0.5: " + " ".join(f"{e:+.2e}" for e in out), flush=True)
