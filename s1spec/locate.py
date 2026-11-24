import numpy as np, sys
from evans import evans
from winding import winding
def newton(f, z, tol=1e-12, it=40):
    for _ in range(it):
        h = 1e-6*(1+abs(z))
        fz = f(z); d = (f(z+h)-f(z-h))/(2*h)
        dz = fz/d; z = z - dz
        if abs(dz) < tol*(1+abs(z)): break
    return z
which = sys.argv[1]; idx = 0 if which=='even' else 1
f = lambda lam: evans(lam)[idx]
def boxes(x0,x1,y0,y1,depth=0):
    w,_ = winding(f,[complex(x0,y0),complex(x1,y0),complex(x1,y1),complex(x0,y1)], n0=16)
    w = int(round(w))
    if w == 0: return []
    if (x1-x0) < 0.3 and (y1-y0) < 0.3 or depth > 6:
        z = newton(f, complex((x0+x1)/2,(y0+y1)/2)); return [(z, w)]
    if (x1-x0) >= (y1-y0):
        xm=(x0+x1)/2; return boxes(x0,xm,y0,y1,depth+1)+boxes(xm,x1,y0,y1,depth+1)
    ym=(y0+y1)/2; return boxes(x0,x1,y0,ym,depth+1)+boxes(x0,x1,ym,y1,depth+1)
# shift box edges slightly off the symmetric/known points to avoid zeros on edges
res = boxes(-1.537, 3.013, -4.021, 4.017)
for z,w in sorted(res, key=lambda t:-t[0].real):
    print(f"{which}: zero at lambda = {z.real:+.10f} {z.imag:+.10f}i  (multiplicity in box {w})", flush=True)
