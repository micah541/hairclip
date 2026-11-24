import numpy as np, sys
from numrange import build, Lam
delta = 0.125
beta = float(sys.argv[1])
for X in [80.0, 160.0, 260.0]:
    x, A = build(beta, X=X, h=0.05)
    row = []
    for th in [0.0, 0.25, 0.35, 0.42, 0.47, 0.50]:
        L = Lam(A, th)
        row.append(f"th={th}: Lam={L:.3f}" + ("" if th == 0 else f" -> Im<={(L+delta*np.cos(th))/np.sin(th):.2f}"))
    print(f"beta={beta} X={X}: " + "; ".join(row), flush=True)
