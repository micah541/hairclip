# float check: physical-frame (HLS Fourier, theta) linearisation at macro snapshots.
# counts non-contractive directions and growth rates -> size of a Lohner block.  Read-only on snapshots.
import numpy as np, sys
sys.path.insert(0,'/home/claude/s1spec/step4')
from hls import HLS
from numpy.fft import irfft
def load(p):
    d=np.load(p); N=int(d['N']); th=irfft(d['th_hat'],N)
    return HLS(th,float(d['L']),d['X0'],cdt=1e-2), d
for name in sys.argv[1:]:
    h,d=load('/home/claude/s1spec/step4/feas2/snaps/%s.npz'%name)
    a=np.abs(h.th_hat); nz=np.nonzero(a>1e-13*a[1:].max())[0]; Neff=nz.max()
    m=min(260,h.N//2-1)
    lin=-(h.kap**4-h.om**2*h.kap**2)
    def F(u): return lin*u+h.nonlin(u)[0]
    u0=h.th_hat.copy(); sc=h.N  # rfft scaling
    eps=1e-6*sc
    J=np.zeros((2*m,2*m))
    for j in range(1,m+1):
        for c,(dz) in enumerate((1.0,1j)):
            up=u0.copy(); um=u0.copy(); up[j]+=eps*dz; um[j]-=eps*dz
            dF=(F(up)-F(um))/(2*eps)
            J[:m,2*(j-1)+c]=dF[1:m+1].real; J[m:,2*(j-1)+c]=dF[1:m+1].imag
    # reorder columns to (re..., im...)
    idx=np.r_[np.arange(0,2*m,2),np.arange(1,2*m,2)]; J=J[:,idx]
    fac=(2*np.pi/h.L)**4
    ev=np.linalg.eigvals(J).real*fac
    sv=np.linalg.eigvalsh((J+J.T)/2)*fac
    # H^1-weighted (theta_alpha) inner product: D J D^-1
    D=np.r_[np.arange(1,m+1),np.arange(1,m+1)].astype(float)
    J1=D[:,None]*J/D[None,:]; s1=np.linalg.eigvalsh((J1+J1.T)/2)*fac
    f=lambda x:(int((x>0).sum()),float(x.max()),float(np.sort(x)[::-1][:60].clip(0).sum()))
    print('%-12s t=%.4f k=%.3g L=%.3f N=%d Neff(1e-13)=%d | eig(J): n+=%d max=%.3g sum+=%.3g | sym L2: n+=%d max=%.3g | sym H1: n+=%d max=%.3g  (rates per unit t)'%
          ((name,float(d['t']),float(d['kmax']),h.L,h.N,Neff)+f(ev)+f(sv)[:2]+f(s1)[:2]),flush=True)
