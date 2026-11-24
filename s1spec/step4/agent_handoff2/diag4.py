import numpy as np
from handoff2 import *
def W(f):
    d = analyse(f); u = d['u']; n = len(u); return d, u - d['alpha']*PSI3[:n] - d['beta2']*H2[:n]
d0, w0 = W('snap_extA_1600.npz'); n0 = len(w0)
nrm0 = np.sum(wh[:n0]*w0*w0)
out = open('unstable_bracket.txt', 'w')
out.write('unstable-mode coefficient c = <w,phi>/<phi,phi>, phi = w of ch0.3 at kmax=1608 (so c(ch0.3,1608)=1); also c*(1608/kmax) removes the k-growth\n')
out.write('chamber kmax alpha*kmax r2  c  c*1608/kmax\n')
for ch, files in [('0.2', ['../snap_k3_L24_ch0.2_400.npz', 'snap_ext_ch0.2_1120.npz', 'snap_ext_ch0.2_1600.npz']),
                  ('0.3', ['../snap_k3_L24_ch0.3_400.npz', 'snap_extA_1120.npz', 'snap_extA_1600.npz']),
                  ('0.4', ['../snap_k3_L24_ch0.4_400.npz', 'snap_ext_ch0.4_1120.npz', 'snap_ext_ch0.4_1600.npz'])]:
    for f in files:
        d, w = W(f); n = min(len(w), n0); c = np.sum(wh[:n]*w[:n]*w0[:n])/nrm0
        line = f"{ch} {d['kmax']:8.1f} {d['alpha']*d['kmax']:+7.3f} {d['r2']:.3e} {c:+.4e} {c*1608/d['kmax']:+.4e}"
        print(line, flush=True); out.write(line + '\n')
