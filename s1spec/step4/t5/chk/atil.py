import numpy as np, sys
sys.path.insert(0,'/home/claude/s1spec/step4/feas')
K0=1.5398
names=['k7','k10','k14','k20','k28','k40','k56','k80','k113','k160']
d=[]
for n in names:
    z=np.load(f'/home/claude/s1spec/step4/feas2/snaps/p2_{n}.npz'); d.append((float(z['t']),float(z['kmax'])))
print(d[:3])
for (t0,k0),(t1,k1),n in zip(d[:-1],d[1:],names):
    # exact for k^-4 linear in t: d(k^-4)/dt = -4 k^-4 dlnk/dt ; a = 1/4 - lam^4 dlnk/dt = 1/4 + (K0^4/4) d(k^-4)/dt
    a=0.25+(K0**4/4)*(k1**-4-k0**-4)/(t1-t0)
    print(f"{n}->: k {k0:.2f}->{k1:.2f}  a_tilde(mean) = {a:+.4f}   (T-t)k^4 at start = {(4.5988639-t0)*k0**4:.4f}")
