import sys, itertools
from ibp_gen import qrf as QRF, coeff as CoefF
from flint import arb, ctx
ctx.prec=150
import ibp_consts as K, ibp_boxes as BX
al=5e-4
for b in BX.boxes():
    if b['tag']=='MC' and abs(float(b['X'].mid())-0.975)<0.006: break
f=lambda z: arb(float(z.mid()))
k_,m_,p_,q_=map(f,(b['k'],b['m'],b['p'],b['q'])); ps=[f(z) for z in b['ps']]; hh=[f(z) for z in b['hh']]
a=[al*ps[i]+al*al*hh[i] for i in range(5)]
print('x',b['X'].mid(),'a',[float(z) for z in a],'k',k_)
alb=1-k_*a[0]; Jb=1/(alb*alb+a[1]*a[1]).sqrt()
cc=CoefF(k_,m_,p_,q_,*a,alb,Jb)
i=2
Z=[sum((K.H[i][j]*f(b['psih'][j+1][r]) for j in range(3)),arb(0)) for r in range(3)]
print('Z',[float(z) for z in Z])
def Q0(w0,w1,w2=0):
    alf=1-k_*(a[0]+w0); Jf=1/(alf*alf+(a[1]+w1)*(a[1]+w1)).sqrt()
    return QRF(*Z,k_,m_,p_,q_,*a,arb(w0),arb(w1),arb(w2),alf,Jf,alb,Jb)['QR']
def quad(w0,w1,w2=0):
    g=lambda c: sum((Z[r]*cc[f'c{c}_{r}'] for r in range(3)),arb(0))
    return g('00')*w0*w0+g('01')*w0*w1+g('11')*w1*w1+g('02')*w0*w2+g('12')*w1*w2+g('22')*w2*w2
for (w0,w1,w2) in [(0.02,0,0),(0,0.02,0),(0.02,0.02,0),(-0.02,0.02,0),(0.005,0.005,0),(0.02,0,0.02),(0,0.02,0.02),(0.02,0.02,0.02),(-0.02,-0.02,0.02),(0.0,0.0,0.02)]:
    q=Q0(w0,w1,w2); qq=quad(w0,w1,w2)
    print((w0,w1,w2),'exact %.4e quad %.4e'%(float(q.mid()),float(qq.mid())))
