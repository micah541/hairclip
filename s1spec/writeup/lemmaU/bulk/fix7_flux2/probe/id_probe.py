# IDENTITY PROBE for FLUX2_DESIGN.md (design stage, Nov 20 2026): exact sympy checks of five hand identities. Not a proof-tree script; stage X1 redoes them.
import sympy as sp
from sympy.calculus.euler import euler_equations
x=sp.Symbol('x'); g=sp.Function('g')(x); W=sp.Function('W')(x); Z=sp.Function('Z')(x); j=sp.Function('j')(x)
def td(e):
    e=sp.expand(e)
    if e==0: return True
    eq=euler_equations(e,[g],x)
    return (not eq) or eq[0]==True or sp.simplify(sp.expand(eq[0].lhs))==0
G=[g.diff(x,i) if i else g for i in range(6)]
expr=(G[1]*G[2]).diff(x)*G[4]
print('residue: (g\'g\'\')\' g\'\'\'\' + (5/2) g\'\' g\'\'\'^2 is a total derivative:', td(expr+sp.Rational(5,2)*G[2]*G[3]**2), '; with 0:', td(expr), '; with 2:', td(expr+2*G[2]*G[3]**2))
kap=4*W*Z*Z.diff(x)*j
print('level-0 hand reduction:', td(kap*G[0]*G[3]-(-kap.diff(x)*G[0]*G[2]+sp.Rational(1,2)*kap.diff(x)*G[1]**2)))
psi=sp.Function('psi')(x); v=Z*g; Wv=sp.exp(-psi)
D2=sp.expand(sp.simplify((Wv*v.diff(x,2)).diff(x,2)/Wv))
p1=psi.diff(x); p2=p1**2-psi.diff(x,2); Zd=lambda i: Z.diff(x,i)
tgt=Z*G[4]+(4*Zd(1)-2*p1*Z)*G[3]+(6*Zd(2)-6*p1*Zd(1)+p2*Z)*G[2]+(4*Zd(3)-6*p1*Zd(2)+2*p2*Zd(1))*G[1]+(Zd(4)-2*p1*Zd(3)+p2*Zd(2))*G[0]
print('D2v coefficients:', sp.simplify(sp.expand(D2-tgt))==0)
K3=3*Zd(1)*G[2]+3*Zd(2)*G[1]+Zd(3)*G[0]; K4=4*Zd(1)*G[3]+6*Zd(2)*G[2]+4*Zd(3)*G[1]+Zd(4)*G[0]
print('j-commutator:', sp.simplify(sp.expand(-Z*(j*G[3]).diff(x)-(-(j*(Z*g).diff(x,3)).diff(x)+j*K4+j.diff(x)*K3)))==0)
s=sp.Symbol('s',positive=True); f=sp.Function('f')(s); Dop=lambda e: s*sp.diff(e,s)
D1=Dop(f);D2_=Dop(D1);D3=Dop(D2_);D4=Dop(D3)
print('D4 conv:', sp.simplify(s**4*f.diff(s,4)-(D4-6*D3+11*D2_-6*D1))==0, ' D3 conv:', sp.simplify(s**3*f.diff(s,3)-(D3-3*D2_+2*D1))==0)
