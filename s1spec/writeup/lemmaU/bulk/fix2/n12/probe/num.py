import sys, pickle, sympy as sp, math
sys.path.insert(0,'/home/claude/s1spec/step3')
from flint import arb, ctx
ctx.prec=128
import hstail as HT, lintail2 as LT
out=pickle.load(open('chi.pkl','rb'))
d=sp.symbols('d0:6'); u0,b0,b1,b2,b3=sp.symbols('u0 b0 b1 b2 b3'); e0=sp.Symbol('e0'); s=sp.Symbol('s',positive=True)
for kk in range(4):
    c=out[kk]
    z={d[i]:0 for i in range(6)}; z[e0]=0
    print(kk,'zero at delta=0,e=0:', sp.simplify(c.subs(z))==0)
def ipow(a,n):
    n=int(n)
    if n<0: return 1/ipow(a,-n)
    r=arb(1)
    for _ in range(n): r=r*a
    return r
IP=sp.Function('ipow')
def fixpow(ex): return ex.replace(lambda q: q.is_Pow and q.exp.is_Integer, lambda q: IP(q.base, q.exp))
mods=[{'sin':lambda a:a.sin(),'cos':lambda a:a.cos(),'ipow':ipow}]
args=(s,e0,*d[:5],u0,b0,b1)
F={}
for kk in range(4):
    c=out[kk]; A=c.subs({b2:0,b3:0}); B2=sp.diff(c,b2); B3=sp.diff(c,b3)
    F[kk]=[sp.lambdify(args,fixpow(ex),modules=mods) for ex in (A,B2,B3)]
def box(x): return arb(0,x)
for s0 in (14,16,20,28):
    S=arb(s0)
    kb,mb,pb,qs=HT.profile_bounds(S); mu=LT.mu(S)
    dl=kb/mu
    U=0.03 if s0<20 else 0.09
    eb=qs+S*dl*dl/2+(pb+U)*dl
    vals=[box(float(dl.upper())),box(float(kb.upper())),box(float(mb.upper())),box(float((pb/4).upper())),box(float((kb*(S+qs)/4).upper()))]
    a=(S,box(float(eb.upper())),*vals,box(U),box(0.007),box(0.02))
    print(s0,'delta',float(dl.upper()),'e0',float(eb.upper()), ' '.join('k%d A %.3e B2 %.3e B3 %.3e'%(kk,*[float(abs(arb(0)+f(*a)).upper()) for f in F[kk]]) for kk in range(4)))
