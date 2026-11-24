import sys; sys.path.insert(0, '/home/claude/s1spec/step3nl')
import os, sys, importlib, numpy as np
sys.argv = ['x', '1e-5', '1e-3', '1e-6', '1e-6']
res = []
for MU in (0.005, 0.02, 0.08):
    for EPSI in (3e-5, 1e-4, 3e-4, 1e-3):
        os.environ['MU'] = str(MU); os.environ['EPSI'] = str(EPSI)
        import inner_closure_float2 as M; importlib.reload(M)
        row = []
        for Ly in (1e-13, 1e-12, 3e-12, 1e-11, 3e-11, 1e-10):
            b, pb, i = M.worst(0.0, Ly); row.append(b / Ly)
        cap = None
        Ls = np.logspace(-14, -9, 26)
        for Ly in Ls:
            if M.worst(0.0, Ly)[0] > 0: cap = Ly; break
        print(f"MU={MU} EPSI={EPSI}: rates(a=0) at 1e-13..1e-10: " + " ".join(f"{r:+.4f}" for r in row) + f"   first positive at Lyap={cap}", flush=True)
