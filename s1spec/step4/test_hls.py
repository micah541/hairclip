import numpy as np, time
from hls import HLS, from_points
# 1) ellipse relaxes to a circle, area conserved
th, L, X0 = from_points(lambda p: 2*np.cos(p) + 1j*np.sin(p), 1, 256)
s = HLS(th, L, X0, omega=1, cdt=2e-2)
A0 = s.area(); t0 = time.time()
s.run_to(tend=2.0)
al, thh, Z, clo = s.curve()
print(f"ellipse: t={s.t:.3f} N={s.N} area {A0:.12f} -> {s.area():.12f} (rel {abs(s.area()/A0-1):.1e}), "
      f"kmax={s.kmax():.6f} (circle {np.sqrt(np.pi/A0):.6f}), closure {clo:.1e}, {time.time()-t0:.1f}s")
# 2) limacon r = 1 + 1.7 sin(phi): inner loop collapses at T = 0.024322 (notes)
th, L, X0 = from_points(lambda p: (1 + 1.7*np.sin(p))*np.exp(1j*p), 2, 1024)
for cdt in [1e-2, 5e-3]:
    s = HLS(th, L, X0, omega=2, cdt=cdt)
    rec = []
    t0 = time.time()
    def cb(S): rec.append((S.t, S.kmax()))
    s.run_to(kstop=300.0, every=5, cb=cb)
    rec = np.array(rec)
    # type-I fit: kmax^-4 linear in t near the end
    sel = rec[:, 1] > 60
    A = np.polyfit(rec[sel, 0], rec[sel, 1]**-4, 1)
    T = -A[1] / A[0]
    print(f"limacon cdt={cdt}: N={s.N} t={s.t:.9f} kmax={s.kmax():.1f}  T_fit={T:.7f}  -d/dt k^-4={-A[0]:.4f}  ({time.time()-t0:.1f}s)", flush=True)
