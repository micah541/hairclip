import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
a=0.75*4**(-1/3)
mu=lambda s:(s/4)**(1/3)
w=np.exp(2j*np.pi*np.arange(3)/3)          # omega_j
def f(s,y):
    p,q,k,m,th=y
    return [-k*q, 1+k*p, m, p/4, k]
def state_from_Z(s,Z):
    """Z=(Z0,Z1,Z2) -> (p,q,k,m) ignoring the O(eps^2) q correction; Y=(k,k',k'')"""
    u=mu(s); lam=-u*w
    k=np.sum(Z).real; kp=np.sum(lam*Z).real; kpp=np.sum(lam**2*Z).real
    return np.array([4*kpp, s+4*k*kp, k, kp])
def M_of(A, sfar=22.0, retfull=False):
    # pure decaying mode at sfar, amplitude A
    y0=np.r_[state_from_Z(sfar,np.array([A,0,0])),0.0]
    ev=lambda s,y: y[1]; ev.direction=0
    sol=solve_ivp(f,[sfar,-6],y0,method='DOP853',rtol=1e-13,atol=1e-300,events=ev,dense_output=True)
    zs=sol.t_events[0]; ys=sol.y_events[0]
    if len(zs)<2: return np.nan if not retfull else (np.nan,sol)
    M=ys[1][3]  # m at second zero of q
    return (M,sol) if retfull else M
for A in [1e-14,3e-14,1e-13,3e-13,1e-12]:
    print(A, M_of(A))

Astar=brentq(M_of,3e-14,1e-13,xtol=1e-30,rtol=1e-14)
M,sol=M_of(Astar,retfull=True)
zs=sol.t_events[0]; ys=sol.y_events[0]
stip=zs[1]; ytip=ys[1]
print("A*=",Astar," M=",M)
print("tip sigma (tail normalization) =",stip, " D=-p=",-ytip[0], " k0=",ytip[2], " q'=",1+ytip[2]*ytip[0])
print("closest-approach zero of q at sigma=",zs[0])
# Z coordinates at S=8
def Z_from_state(s,y):
    p,q,k,m=y[:4]; u=mu(s); lam=-u*w
    V=np.vstack([np.ones(3),lam,lam**2]); Y=np.array([k,m,p/4])
    return np.linalg.solve(V,Y)
for S in [6.0,8.0,10.0]:
    Z=Z_from_state(S,sol.sol(S)); print(f"S={S}: Z0={Z[0].real:.6e}  gamma=Z1/Z0={Z[1]/Z[0]:.6f}  k(S)={sol.sol(S)[2]:.4e}")
