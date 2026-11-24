import sys, numpy as np, time
sys.path.insert(0,'/home/claude/s1spec/step4/feas'); sys.path.insert(0,'/home/claude/s1spec/step4/agent_handoff2'); sys.path.insert(0,'/home/claude/s1spec/step4')
t=time.time()
import chart as CH
H=CH.H
print("import",time.time()-t, "K0",H.K0,"D_S1",H.D_S1, len(H.xs), H.xs[-1])
from measure import STIP, load
print("STIP",STIP)
S,z=load('snap_main_2000.npz')
al,th=S.theta(1); N=len(al)
kt=abs((2*np.pi/S.L)*np.gradient(th,al)[N//2]); lam0=H.K0/kt
print("N",N,"L",S.L,"kt",kt,"lam0",lam0,"ds/lam",S.L/N/lam0, "kmax",S.kmax())
