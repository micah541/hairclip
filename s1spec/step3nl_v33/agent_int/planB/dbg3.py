import sys; sys.argv=['x','1e-3','1e-2','1e-2','1']
src=open('pb_nlin.py').read().split("if __name__")[0]
exec(src)
bl = make_segments(BPS, [w * 1 for w in (0.02, 0.02, 0.025, 0.05, 0.1, 0.2, 0.2, 0.2)])
found=None
for n,(sa,sb) in enumerate(bl):
    b = make_box(sa, sb); y = PSI.jets_any(b.X, b.xc, 5); h = H2.jets_any(b.X, b.xc, 5)
    d = box_vals(y, h, b.p, b.q, b.k, b.m)
    bad=[k_ for k_,v in d.items() if v.is_nan()]
    if bad: found=(n,b,y,h,bad); print(n,b.xa,bad); break
n,b,y,h,bad=found
uj = [A * y[j] + A * A * h[j] for j in range(5)]; U0,U1=uj[0]+E0b,uj[1]+E1b
print('U0',U0,'U1',U1,'k',b.k,'m',b.m,'p',b.p,'q',b.q)
for key in bad:
    e = {'A0':dU0,'A1':dU1}.get(key, coef.get(key) if not isinstance(key,str) else None)
    if e is None:
        e = coef[eval(key)]
    for t in sp.Add.make_args(sp.expand(e) if False else e):
        v = mk(t)(U0,U1,uj[2],uj[3],uj[4],b.p,b.q,b.k,b.m)
        if v.is_nan(): print('BAD term', key, str(t)[:300])
t=[t for t in sp.Add.make_args(dU0) if str(t).startswith('-36*U1**2*a2**3')][0]
print(_Pr().doprint(t)[:600])
k=b.k; a2=uj[2]
den=(U1*U1+ (-U0*k+1)*(-U0*k+1)); print('den',den, den.lower()>0, pw(den,5))
print('num', -arb(36)*(U1*U1)*pw(a2,3)*k*(-U0*k+1))
