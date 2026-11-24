import sys, numpy as np
sys.path.insert(0,'/home/claude/s1spec/step4')
from numpy.fft import irfft
def halfcurve(path):
    z=np.load(path); N=int(z['N']); L=float(z['L']); h=z['th_hat']
    k=np.arange(N//2+1)
    d1=1+irfft(1j*k*h,N); 
    kap=(2*np.pi/L)*d1
    al=2*np.pi*np.arange(N)/N
    idx=N//2-np.arange(N//2+1)  # from tip (alpha=pi) down to alpha=0
    s=(L/N)*np.arange(N//2+1)
    return s, -kap[idx], L, float(z['kmax'])
for nm in sys.argv[1:]:
    s,kap,L,km=halfcurve(nm)
    lam=1.5398/abs(kap[0]); print(nm, 'kmax',km,'ktip',kap[0],'L',L,'lam',lam,'ell',L/2/lam)
    # physical curvature profile: print at selected s
    for ss in [0,0.05,0.1,0.2,0.5,1,2,3,4,6,8,10,12,14,16,18,20,22,23,24,24.2]:
        i=np.argmin(abs(s-ss)); print(f"  s={s[i]:7.3f} x={s[i]/lam:9.2f} kphys={kap[i]: .5e} kresc={kap[i]*lam: .4e}")
