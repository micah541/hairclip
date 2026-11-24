import sys, numpy as np
sys.dont_write_bytecode=True
from run2 import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
ts=['0','0.3','1.5','3.0','4.0','4.5','4.59']
fig,ax=plt.subplots(2,4,figsize=(20,9)); ax=ax.ravel()
for a,t in zip(ax,ts):
    S=load_state(f'{OUT}/state_t{t}.npz'); al,th,Z,clo=S.curve(2)
    a.plot(Z.real,Z.imag,lw=.6); a.set_title(f't={t} kmax={S.kmax():.2f} L={S.L:.2f} N={S.N}'); a.set_aspect('equal')
    N=len(al); a.plot(Z.real[N//2],Z.imag[N//2],'ro',ms=3); a.plot(Z.real[0],Z.imag[0],'go',ms=3)
plt.savefig('geom.png',dpi=70)
