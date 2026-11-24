import sys, itertools
sys.argv=['x','5e-4','0.02']
from ibp_gen import linf as LinF, coeff as CoefF, remf as RemF
from flint import arb, ctx
ctx.prec=150
import ibp_consts as K, ibp_boxes as BX
AL=arb(0,5e-4); ETAB=arb(0,0.02)
mag=lambda z: arb(abs(z).upper())
best={}
for b in BX.boxes():
    if b['strip']: continue
    k_,m_,p_,q_=b['k'],b['m'],b['p'],b['q']; ps,hh=b['ps'],b['hh']
    a=[AL*ps[i]+AL*AL*hh[i] for i in range(5)]
    alb=1-k_*a[0]; Jb=1/(alb*alb+a[1]*a[1]).sqrt()
    w0=ETAB;w1=ETAB; alf=1-k_*(a[0]+w0); Jf=1/(alf*alf+(a[1]+w1)*(a[1]+w1)).sqrt()
    tt=RemF(k_,m_,p_,q_,*a,w0,w1,alf,Jf,alb,Jb)
    i=2
    Z=[sum((K.H[i][j]*b['psih'][j+1][r] for j in range(3)),arb(0)) for r in range(3)]
    for pre,n in (('t3',4),('t2',3),('t1',2)):
        for s in range(n):
            for r in range(3):
                v=float(mag(Z[r]*tt[f'{pre}_{s}_{r}']).upper())
                key=(pre,s,r)
                if v>best.get(key,(0,0))[0]: best[key]=(v,float(b['X'].mid()))
            v=float(mag(sum((Z[r]*tt[f'{pre}_{s}_{r}'] for r in range(3)),arb(0))).upper())
            key=(pre,s,'sum')
            if v>best.get(key,(0,0))[0]: best[key]=(v,float(b['X'].mid()))
for k_,v in sorted(best.items(), key=lambda t:-t[1][0])[:16]: print(k_, "%.4g at x=%.3f"%v)
