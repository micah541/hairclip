import numpy as np
from handoff2 import *
def W(f):
    d = analyse(f); u = d['u']; n = len(u); return d, u - d['alpha']*PSI3[:n] - d['beta2']*H2[:n]
d0, w0 = W('snap_extA_1600.npz'); n0 = len(w0); nrm0 = np.sum(wh[:n0]*w0*w0)
out = open('numfloor_stable.txt', 'w'); out.write('coarse (cdt1e-2,tail1e-10) vs fine (cdt5e-3,tail1e-12): r2 with unstable phi removed\n kmax alpha c r2_removed r2_removed/|a|^3\n')
for f in ['snap_repA_400.npz', 'snap_fineB_400.npz', 'snap_extA_800.npz', 'snap_extB_800.npz']:
    d, w = W(f); n = min(len(w), n0); c = np.sum(wh[:n]*w[:n]*w0[:n])/nrm0
    rw = np.sqrt(np.sum(wh[:n]*(w[:n]-c*w0[:n])**2)); l = f"{f} {d['kmax']:8.1f} {d['alpha']:+.4e} {c:+.4e} {rw:.3e} {rw/abs(d['alpha'])**3:8.2f}"; print(l, flush=True); out.write(l+'\n')
