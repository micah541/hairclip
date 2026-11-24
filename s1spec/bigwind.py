import numpy as np, sys
import evans as E
from scipy.integrate import solve_ivp
# extend profile beyond sigma=22 by the exact ray (k=p=m=0, q=sigma): error ~ e^{-E(22)} negligible
_PROF = E.PROF
def PROF(s):
    return _PROF(s) if s <= 22.0 else np.array([0.0, s, 0.0, 0.0])
E.PROF = PROF
def evans_big(lam):
    sfar = max(22.0, (40*abs(lam))**0.75)
    return E.evans(lam, sfar=sfar)
from winding import winding
which = sys.argv[1]; idx = 0 if which=='even' else 1
f = lambda lam: evans_big(lam)[idx]
for (x0,x1,y0,y1) in [(-0.2,3.0,4.0,25.0),(3.0,25.0,-25.0,25.0),(-0.2,3.0,-25.0,-4.0),(-0.2,3.0,-4.0,4.0)]:
    w,n = winding(f,[complex(x0,y0),complex(x1,y0),complex(x1,y1),complex(x0,y1)], n0=40)
    print(f"{which}: Re[{x0},{x1}] x Im[{y0},{y1}] winding = {w:.4f} ({n} evals)", flush=True)
