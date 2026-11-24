import sys; sys.argv=['x','1e-3','1e-2','1e-2','1']
src=open('pb_nlin.py').read().split("if __name__")[0]
exec(src)
bl = make_segments(BPS, [w * 1 for w in (0.02, 0.02, 0.025, 0.05, 0.1, 0.2, 0.2, 0.2)])
sa,sb=bl[2859]; b=make_box(sa,sb)
y = PSI.jets_any(b.X, b.xc, 5); h = H2.jets_any(b.X, b.xc, 5)
uj = [A * y[j] + A * A * h[j] for j in range(5)]
U0,U1=uj[0]+E0b,uj[1]+E1b
print('U0',U0,'U1',U1)
for kk in (arb(0),arb(0,1e-100),arb(0,3.9e-294),arb(0,1e-290)):
    for mm in (arb(0), arb(0,1e-290)):
        v=FU0(U0,U1,uj[2],uj[3],uj[4],arb(0,1e-290),b.q,kk,mm); print(kk,mm,v.is_nan(), str(v)[:40])
open('dU0_str.txt','w').write(str(dU0))
print('---')
for U0t in (arb(2300), arb(0,42), arb(2300,42)):
  for U1t in (arb(0.37), arb(0,0.386), arb(0.37,0.01)):
    for aa in ((uj[2],uj[3],uj[4]), (arb(0),arb(0),arb(0))):
        v=FU0(U0t,U1t,*aa,arb(0),b.q,arb(0),arb(0)); print(U0t,U1t,'zero' if aa[0]==0 else 'jets',v.is_nan())
print([str(t) for t in uj[2:]])
print('=== terms')
import sympy as sp
args = sp.Add.make_args(dU0)
bad=[]
for t in args:
    f = mk(t)
    v = f(arb(2300), arb(0,0.387), arb(0,1.8e-3), arb(0,1e-6), arb(0,1e-9), arb(0), b.q, arb(0), arb(0))
    if v.is_nan(): bad.append(t)
print(len(args), len(bad))
for t in bad[:3]: print(str(t)[:400])
print('=== pieces')
print(_Pr().doprint(bad[0]))
U1b=arb(0,0.387); one=arb(1)
den=(U1b*U1b+one*one); print(den, (den*den*den*den*den))
print(1/(den*den*den*den*den))
