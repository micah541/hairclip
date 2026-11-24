import sys, time; sys.path.insert(0,'/home/claude/s1spec/step4'); sys.path.insert(0,'/home/claude/s1spec'); sys.path.insert(0,'.')
import numpy as np, curves
from hls import HLS
def trial(Rd,w,Rc,phi,Lw,l0=0.3,tmax=30,wall=40):
    th,L,X0,b,e=curves.disk_theta(8192,Rd,w,Rc,phi,Lw,l0)
    S=HLS(th,L,X0,omega=1,cdt=1e-2,tail_tol=1e-10); t0=time.time(); kmin=1e9; n=0
    while S.t<tmax and time.time()-t0<wall:
        S.run_to(maxsteps=20); km=S.kmax(); kmin=min(kmin,km)
        if km>150 and S.t>1e-3: return ('SING',S.t,km,L)
    return ('smooth',S.t,S.kmax(),L)
for args in [(2.5,0.05,0.3,0,2.0),(2.5,0.05,0.3,0,3.5),(2.0,0.05,0.5,0,2.5),(3,0.05,0.3,0,5.0),(2.5,0.05,0.3,np.radians(25),1.8)]:
    print(args, trial(*args), flush=True)
