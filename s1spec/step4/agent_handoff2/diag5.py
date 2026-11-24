import numpy as np
from handoff2 import *
def W(f):
    d = analyse(f); u = d['u']; n = len(u); return d, u - d['alpha']*PSI3[:n] - d['beta2']*H2[:n]
d0, w0 = W('snap_extA_1600.npz'); n0 = len(w0); nrm0 = np.sum(wh[:n0]*w0*w0)
out = open('unstable_removed_k3_L24_ch0.3.txt', 'w')
out.write('ch0.3: r2 before / after removing the fitted unstable-mode component c*phi from w (phi = w at kmax 1608)\n kmax alpha r2 c r2_without_phi r2w/|a|^3\n')
for f in ['../snap_k3_L24_ch0.3_%d.npz' % k for k in (100, 200, 280, 400)] + ['snap_extA_%d.npz' % k for k in (800, 1600, 3200)]:
    d, w = W(f); n = min(len(w), n0); c = np.sum(wh[:n]*w[:n]*w0[:n])/nrm0
    rw = np.sqrt(np.sum(wh[:n]*(w[:n]-c*w0[:n])**2))
    l = f"{d['kmax']:8.1f} {d['alpha']:+.3e} {d['r2']:.3e} {c:+.3e} {rw:.3e} {rw/abs(d['alpha'])**3:9.2f}"; print(l, flush=True); out.write(l+'\n')
