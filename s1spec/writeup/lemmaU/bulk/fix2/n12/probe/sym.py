import sympy as sp, time
t0=time.time()
s=sp.Symbol('s', positive=True)
dl=sp.Function('dl')(s); u=sp.Function('u')(s); e=sp.Function('e')(s)
k=dl.diff(s); b=u.diff(s)
eta=2*sp.sin(dl/2)**2 + b*sp.sin(dl) + u*k*sp.cos(dl)
eps=(e/s+eta)/(1-eta)
E=sp.Function('E')(s)
N=(1-k*u+b**2)*sp.sin(dl)+k*u*b*sp.cos(dl); Dn=(1-k*u)*sp.cos(dl)-b*sp.sin(dl); Phi=N/Dn
T=lambda F:(1+E)*s*F.diff(s); Dop=lambda F:s*F.diff(s)
def subs_all(ex,K):
    # E^(i) -> eps derivatives; then e^(n) -> -(eta)^(n-1)
    for i in range(K,-1,-1):
        ex=ex.subs(E.diff(s,i) if i else E, eps.diff(s,i) if i else eps)
    for n in range(6,0,-1):
        ex=ex.subs(e.diff(s,n), -eta.diff(s,n-1))
    d=sp.symbols('d0:6'); us=sp.symbols('u0 b0 b1 b2 b3 b4'); e0=sp.Symbol('e0')
    for n in range(5,-1,-1):
        ex=ex.subs(u.diff(s,n) if n else u, us[n])
    for n in range(5,-1,-1):
        ex=ex.subs(dl.diff(s,n) if n else dl, d[n])
    ex=ex.subs(e,e0)
    return ex
out={}
Fb=b; Fp=Phi; Db=b
for kk in range(4):
    if kk>0:
        Fb=T(Fb); Fp=T(Fp); Db=Dop(Db)
    diffb=sp.expand(Fb-Db)
    chi=subs_all(diffb,kk)+subs_all(Fp,kk)
    out[kk]=chi
    fs=chi.free_symbols
    print(kk,'free',sorted(map(str,fs)),'ops',sp.count_ops(chi), round(time.time()-t0,1),flush=True)
b2,b3,b4=sp.symbols('b2 b3 b4')
for kk in range(4):
    c=out[kk]
    print(kk,'d2/db2^2',sp.simplify(sp.diff(c,b2,2))==0,'d/db2db3',sp.simplify(sp.diff(c,b2,b3))==0,'d2/db3^2',sp.simplify(sp.diff(c,b3,2))==0,'b4 in', b4 in c.free_symbols)
import pickle; pickle.dump(out,open('chi.pkl','wb'))
