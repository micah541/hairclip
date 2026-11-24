import numpy as np
from scipy.integrate import solve_ivp
S=8.0; mu=(S/4)**(1/3)
w1=complex(-0.5,np.sqrt(3)/2); w2=w1.conjugate()
g=complex(-0.010356215835809902,0.000332793623115475)
def state(eps):
    Z1=eps*g
    k=eps+2*Z1.real; kp=-mu*(eps+2*(w1*Z1).real); kpp=mu**2*(eps+2*(w2*Z1).real)
    return [4*kpp, S+4*k*kp, k, kp, 0.0]
f=lambda s,y:[-y[2]*y[1],1+y[2]*y[0],y[3],y[0]/4,y[2]]
def M(eps,full=False):
    ev=lambda s,y:y[1]; ev.direction=0
    sol=solve_ivp(f,[S,-6],state(eps),method='DOP853',rtol=1e-13,atol=1e-16,events=ev,dense_output=True)
    z=sol.t_events[0]; y=sol.y_events[0]
    return (y[1][3], z, y, sol) if full else y[1][3]
for e in [1.40e-4,1.45e-4,1.50e-4,1.52e-4,1.53e-4,1.54e-4,1.55e-4,1.56e-4,1.60e-4,1.65e-4]:
    m,z,y,sol=M(e,True); print(f"eps={e:.3e}  M={m:+.5f}  zeros of q at sigma={np.round(z,4)}  q at sigma=-0.5: {sol.sol(-0.5)[1]:+.4f}")
from scipy.optimize import brentq
es=brentq(M,1.50e-4,1.56e-4,xtol=1e-20,rtol=1e-14); print("eps* =",es)
