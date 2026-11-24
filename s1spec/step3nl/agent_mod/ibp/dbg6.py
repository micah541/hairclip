import sys, itertools
from ibp_gen import remf as RemF
from flint import arb, ctx
ctx.prec=150
import ibp_consts as K, ibp_boxes as BX
al=5e-4
for b in BX.boxes():
    if b['tag']=='MC' and abs(float(b['X'].mid())-0.975)<0.006: break
f=lambda z: arb(float(z.mid()))
mag=lambda z: float(abs(z).upper())
def run(point_params, eta, ALB=0.0):
    ff=(lambda z: f(z)) if point_params else (lambda z: z)
    k_,m_,p_,q_=map(ff,(b['k'],b['m'],b['p'],b['q'])); ps=[ff(z) for z in b['ps']]; hh=[ff(z) for z in b['hh']]
    AL=arb(0,ALB) if ALB else arb(al)
    a=[AL*ps[i]+AL*AL*hh[i] for i in range(5)]
    alb=1-k_*a[0]; Jb=1/(alb*alb+a[1]*a[1]).sqrt()
    W=arb(0,eta); w0=W; w1=W
    alf=1-k_*(a[0]+w0); Jf=1/(alf*alf+(a[1]+w1)*(a[1]+w1)).sqrt()
    tt=RemF(k_,m_,p_,q_,*a,w0,w1,alf,Jf,alb,Jb)
    i=2
    Z=[sum((K.H[i][j]*ff(b['psih'][j+1][r]) for j in range(3)),arb(0)) for r in range(3)] if not point_params else [sum((f(K.H[i][j])*f(b['psih'][j+1][r]) for j in range(3)),arb(0)) for r in range(3)]
    return [max(mag(sum((Z[r]*tt[f't3_{s}_{r}'] for r in range(3)),arb(0))) for s in range(4)), max(mag(sum((Z[r]*tt[f't2_{s}_{r}'] for r in range(3)),arb(0))) for s in range(3))]
for pp in (False,True):
    for eta in (0.02,0.01,0.005,0.0025,0.00125,1e-5):
        print('point_params',pp,'eta',eta,run(pp,eta))
