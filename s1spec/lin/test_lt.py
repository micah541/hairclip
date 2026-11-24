import time, sys
from lintail import admissible_W, eta_bound
boxes = [(1.0, 1.0, 0.0, 0.0), (-0.125, -0.125, 0.0, 0.0), (-0.125, -0.125, 1.0, 1.0), (-0.125, -0.12, 3.9, 3.95),
         (1.75, 1.8, 0.3, 0.35), (0.3, 0.31, 0.2, 0.21), (-0.125, -0.125, 0.5, 0.5), (-0.125, -0.125, 2.0, 2.0)]
for box in boxes:
    t0 = time.time()
    try:
        R = admissible_W(box)
        print(f"box {box}: pieces={eta_bound.npieces} eta={float(R['eta'].upper()):.4f} Lip={float(R['Lip'].upper()):.3f} "
              f"gmin={R['gmin']:.3f} ({time.time()-t0:.1f}s)", flush=True)
    except Exception as e:
        print(f"box {box}: FAILED {e} ({time.time()-t0:.1f}s)", flush=True)
