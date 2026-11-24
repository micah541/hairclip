import numpy as np
from scipy.integrate import quad
phi=lambda x:(1+x*x)**(2/3)/4
d1=lambda x:x/3*(1+x*x)**(-1/3)
d2=lambda x:(1+x*x)**(-1/3)/3-2/9*x*x*(1+x*x)**(-4/3)
x=np.linspace(0,34.25,3_000_001)
print('c1',0.5*max(0,(d1(x)**2-d2(x)).max()),'Phi1',d1(34.25))
print('K',[np.exp(phi(r)) for r in (1,2,4,8,12,20,34.25)])
S=lambda t:t**4*(35-84*t+70*t**2-20*t**3)
def sp(x):
    t=np.clip((x-30.25)/8,0,1);return S(t)
def Sd(x):
    t=(x-30.25)/8; return np.where((t>0)&(t<1),140*t**3*(1-t)**3/8,0)
A=lambda x:6*x/(1+x*x); B=d1
def p1(x): s=sp(x); return (1-s)*A(x)+s*B(x)
def p2(x):
    s=sp(x);Ap=6*(1-x*x)/(1+x*x)**2
    return (1-s)*Ap+s*d2(x)+Sd(x)*(B(x)-A(x))
for R in (40,60,80):
    y=np.linspace(6,R,3_000_001);print('c1W',R,0.5*max(0,(p1(y)**2-p2(y)).max()),'Psi1',p1(y).max())
psi=lambda r:3*np.log(1+30.25**2)+quad(p1,30.25,r)[0] if r>30.25 else 3*np.log(1+r*r)
print('KW',[np.exp(psi(r)) for r in (20,30,40,60,80)])
# numeric derivative check of p2
h=1e-6;print(max(abs((p1(y+h)-p1(y-h))/(2*h)-p2(y)).max() for y in [np.linspace(6,80,100001)]))
t=np.linspace(0,1,2_000_001)
import sympy as sy
T=sy.symbols('t');Ss=T**4*(35-84*T+70*T**2-20*T**3)
print([ (np.abs(sy.lambdify(T,sy.diff(Ss,T,j))(t))).max()/3**j for j in range(1,5)])
