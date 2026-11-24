import numpy as np
from scipy.integrate import solve_ivp
def f(s,y):
    p,q,k,m,th=y
    return [-k*q, 1+k*p, m, p/4, k]
D,k0=2.500806,1.539821
sol=solve_ivp(f,[0,14],[-D,0,k0,0,-np.pi/2],method='DOP853',rtol=1e-13,atol=1e-15,dense_output=True,max_step=0.01)
s=sol.t; p,q,k,m,th=sol.y
G=np.exp(1j*th)*(q+1j*p)
# report
for S in [1,2,3,4,5,6,7,8,9,10,11,12,13]:
    i=np.argmin(abs(s-S)); print(f"s={S:5.1f} k={k[i]: .3e} p={p[i]: .3e} q-s={q[i]-s[i]: .4f} |G|={abs(G[i]):.4f} arg={np.degrees(np.angle(G[i])):.3f} theta={np.degrees(th[i]):.3f}")
# axis crossings of the branch (Im G changes sign), excluding the tip
y=G.imag; idx=np.nonzero(np.sign(y[1:])!=np.sign(y[:-1]))[0]
print("axis crossings at s:",[round(s[i],4) for i in idx], "x:",[round(G[i].real,4) for i in idx])
# zeros of q
zq=np.nonzero(np.sign(q[1:])!=np.sign(q[:-1]))[0]; print("zeros of q at s:",[round(s[i],4) for i in zq])
