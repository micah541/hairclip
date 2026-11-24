import sys, numpy as np
sys.dont_write_bytecode=True
from run2 import *
from cross import ncross_general
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
ts=['0.5','1.5','3.0','4.5']
fig,ax=plt.subplots(4,1,figsize=(14,16))
for a,t in zip(ax,ts):
    S=load_state(f'{OUT}/state_t{t}.npz'); al,th,Z,clo=S.curve(4); N=len(al)
    ell=1.5398/S.kmax(); tip=Z[N//2]
    c,pts=ncross_general(S)
    # colour by sigma from tip (arclength/ell) along both strands
    sig=np.abs(al-np.pi)*S.L/(2*np.pi)/ell
    sc=a.scatter(Z.real,Z.imag,c=np.log10(np.exp(-0.25*(1+sig**2)**(2/3))+1e-300),s=1,cmap='viridis',vmin=-12,vmax=0)
    a.plot(tip.real,tip.imag,'ro'); 
    for p in pts: a.plot(p.real,p.imag,'kx')
    a.set_title(f't={t} kmax={S.kmax():.3f} ell={ell:.3f} L={S.L:.1f} crossings x={[round(p.real,2) for p in pts]}'); a.set_aspect('equal'); a.set_xlim(tip.real-3,tip.real+22)
    plt.colorbar(sc,ax=a,label='log10 rho(sigma)')
plt.tight_layout(); plt.savefig('geomzoom.png',dpi=60)
