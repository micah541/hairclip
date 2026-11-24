import numpy as np, sys
from evans import evans

def winding(fun, corners, n0=40, maxref=12):
    """argument principle: winding number of fun along the polygon 'corners' (closed), with adaptive refinement
    so that the phase increment between consecutive samples stays below 0.5 rad."""
    total = 0.0
    pts = list(corners) + [corners[0]]
    nevals = 0
    for a, b in zip(pts[:-1], pts[1:]):
        ts = list(np.linspace(0, 1, n0 + 1))
        vals = {t: fun(a + (b - a) * t) for t in ts}; nevals += len(ts)
        i = 0
        while i < len(ts) - 1:
            t0, t1 = ts[i], ts[i + 1]
            d = np.angle(vals[t1] / vals[t0])
            if abs(d) > 0.5 and (t1 - t0) > 1e-6:
                tm = 0.5 * (t0 + t1); vals[tm] = fun(a + (b - a) * tm); nevals += 1
                ts.insert(i + 1, tm)
                continue
            total += d; i += 1
    return total / (2 * np.pi), nevals

if __name__ == "__main__":
    which = sys.argv[1]  # 'even' or 'odd'
    idx = 0 if which == 'even' else 1
    f = lambda lam: evans(lam)[idx]
    for (x0, x1, y0, y1) in [(-0.5, 3.0, -4.0, 4.0), (-1.5, 3.0, -4.0, 4.0), (0.6, 3.0, -4.0, 4.0),
                             (0.1, 0.6, -4.0, 4.0), (-0.1, 0.1, -1.0, 1.0)]:
        w, n = winding(f, [complex(x0, y0), complex(x1, y0), complex(x1, y1), complex(x0, y1)])
        print(f"{which}: rectangle Re[{x0},{x1}] x Im[{y0},{y1}]: winding = {w:.4f}  ({n} evaluations)", flush=True)
