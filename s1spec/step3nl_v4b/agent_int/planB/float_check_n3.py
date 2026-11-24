import sys; sys.argv=['x','1e-3','1']
from pb_common import *
import pb_n3 as M
from c4fun import PSI,H2
for xv in (20.0,60.0,120.0,200.0,250.0):
    X=arb(xv); y=PSI.jets_any(X,xv,6); h=H2.jets_any(X,xv,6)
    sg=arb(xv)+ST.mid(); z=arb(0)
    from flint import arb_series
    u0s=[arb_series([z,y[j],h[j]],prec=8) for j in range(5)]
    c=(M.N0(u0s,arb(0),sg,arb(0),arb(0))).coeffs()+[arb(0)]*8
    print(xv,'y1,y2',float(y[1].mid()),float(y[2].mid()),'c2..c7',[float('%.3g'%float(abs(c[i]).mid())) for i in range(2,8)], ' alpha^j|c_j|:', [float('%.2g'%(1e-3**i*float(abs(c[i]).mid()))) for i in range(2,8)])
