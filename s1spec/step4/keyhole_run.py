import os, sys, subprocess, numpy as np
R = sys.argv[1]; tag = sys.argv[2]
os.environ['CHAMBER'] = R
import hairpin
def defect(d):
    os.environ['CHDW'] = repr(float(d))
    s, k, L = hairpin.curvature_fn(1 << 20, 24, 0.025, 1.0, 0.05)
    ds = L / len(s); th = -np.pi / 2 + np.concatenate([[0], np.cumsum(0.5 * (k[1:] + k[:-1]) * ds)])
    return float((np.sum(np.exp(1j * th)) * ds).imag)
d0, d1 = 0.0, 0.003; f0, f1 = defect(d0), defect(d1)
for _ in range(6):
    d2 = d1 - f1 * (d1 - d0) / (f1 - f0); d0, f0, d1, f1 = d1, f1, d2, defect(d2)
env = dict(os.environ, CHAMBER=R, CHDW=repr(float(d1)), ARM="24", WID="1", SLOT="0.05", AMOLL="0.025")
subprocess.run(["python3", "run_hairpin.py", tag, "1e-2", "16384", "450"], env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=1500)
