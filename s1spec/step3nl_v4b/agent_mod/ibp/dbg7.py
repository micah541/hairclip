import sympy as sp, time
import ibp_sym as S
from ibp_codegen import make_func
w0,w1,w2=S.w[0],S.w[1],S.w[2]
V1=sp.Symbol('V1')
# QR r=0 , w2^0 part
Qr=sp.expand(sp.diff(S.QR,S.Y[0])); Q0=sp.expand(sp.Poly(Qr,w2).coeff_monomial(1))
print('terms Q0',len(Q0.args))
def reduce_expr(e):
    e=sp.expand(e.subs(w1,V1-S.a[1]))
    # replace V1^n (n>=2) : V1^2 -> 1/Jf^2 - alf^2
    P=sp.Poly(e,V1); out=0
    for (n,),c in P.terms():
        r,s=divmod(n,2)
        out+=c*(1/S.Jf**2-S.alf**2)**r*V1**s
    return sp.expand(out)
def D0(e): return S.dw(e,0)
# derivative operators in terms of V1: dw1 = d/dV1 + chain
def dV1(e):
    return sp.expand(sp.diff(e,V1)+sp.diff(e,S.Jf)*(-S.Jf**3*V1))
def dW0(e): return sp.expand(sp.diff(e,w0)+sp.diff(e,S.alf)*(-S.k)+sp.diff(e,S.Jf)*(S.Jf**3*S.alf*S.k))
R=reduce_expr(Q0); print('reduced terms',len(R.args))
def red(e): return reduce_expr_V(e)
def reduce_expr_V(e):
    P=sp.Poly(sp.expand(e),V1); out=0
    for (n,),c in P.terms():
        r,s=divmod(n,2); out+=c*(1/S.Jf**2-S.alf**2)**r*V1**s
    return sp.expand(out)
T_old=S.dw(S.dw(S.dw(Q0,1),1),1)
T_new=reduce_expr_V(dV1(reduce_expr_V(dV1(reduce_expr_V(dV1(R))))))
print('T3 w1w1w1: terms old',len(T_old.args),'new',len(T_new.args))
from flint import arb, ctx
import ibp_consts as K, ibp_boxes as BX
args=[S.k,S.m,S.p,S.q,*S.a[:5],w0,w1,S.alf,S.Jf,S.alb,S.Jb]
fo,_=make_func({'T':T_old},args,'fo')
args2=[S.k,S.m,S.p,S.q,*S.a[:5],w0,V1,S.alf,S.Jf,S.alb,S.Jb]
fn,_=make_func({'T':T_new},args2,'fn')
al=5e-4
for b in BX.boxes():
    if b['tag']=='MC' and abs(float(b['X'].mid())-0.975)<0.006: break
f=lambda z: arb(float(z.mid()))
k_,m_,p_,q_=map(f,(b['k'],b['m'],b['p'],b['q'])); ps=[f(z) for z in b['ps']]; hh=[f(z) for z in b['hh']]
AL=arb(al); a=[AL*ps[i]+AL*AL*hh[i] for i in range(5)]
alb=1-k_*a[0]; Jb=1/(alb*alb+a[1]*a[1]).sqrt()
for eta in (0.02,0.005,1e-6):
    W=arb(0,eta); alf=1-k_*(a[0]+W); Jf=1/(alf*alf+(a[1]+W)*(a[1]+W)).sqrt()
    vo=fo(k_,m_,p_,q_,*a,W,W,alf,Jf,alb,Jb)['T']; vn=fn(k_,m_,p_,q_,*a,W,a[1]+W,alf,Jf,alb,Jb)['T']
    print(eta,'old',float(abs(vo).upper()),'new',float(abs(vn).upper()))
