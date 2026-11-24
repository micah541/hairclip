"""Symbolic integration by parts of the modulation pairings <chi N(w), f_j>_rho  (task F3).

Notation (all symbols are *jets*; Dx is the total x-derivative, a derivation on the polynomial ring):
  k,m,p,q       S1 profile, ODE  k'=m, m'=p/4, p'=-kq, q'=1+kp
  a0..a6        jets of the background u0 = alpha psi3 + alpha^2 h2
  w0..w5        jets of w
  Y0..Y5        jets of the TEST function Psi = chi f_j rho   (Y_r = Psi^(r), no weight divided out)
  alf,Jf        1-k(a0+w0), 1/sqrt(alf^2+(a1+w1)^2)  (full state)     alb,Jb same at w=0 (background)
F(u0..u4) = sum_mon c_mon(Ji,al,u0,u1) u2^i u3^j u4^l + (p+u0-u1 q/al)/4  (graph_struct.pkl, same as nlin_E2).
  N(w)   = F(a+w) - F(a) - L0 w                      L0 = linearisation at u=0 (alpha=0)
         = Lin(w) + R(w),   Lin = DF(a)[w] - L0 w,   R = F(a+w)-F(a)-DF(a)[w]
Results (exact identities, verified with residual 0 by `verify()`):
  Y0*Lin(w) = w0*Lam + Dx(Blin)              Lam = sum_i (-1)^i Dx^i(dC_i Y0)
  Y0*R(w)   = QR(w0,w1,w2) + Dx(Bq)          QR free of w3,w4;  QR = O(|w|^2)
"""
import pickle, sympy as sp, time, sys
from math import comb
t00 = time.time()
GS = pickle.load(open('/home/claude/s1spec/step3nl/graph_struct.pkl', 'rb'))
k, m, p, q = sp.symbols('k m p q')
a = sp.symbols('a0:7'); w = sp.symbols('w0:6'); Y = sp.symbols('Y0:6')
alf, Jf, alb, Jb = sp.symbols('alf Jf alb Jb')
oal, oJi, ou0, ou1 = sp.symbols('al Ji u0 u1')
U = [a[i] + w[i] for i in range(5)]
alf_p = -m * U[0] - k * U[1]                      # d/dx alf
alb_p = -m * a[0] - k * a[1]
rules = {k: m, m: p / 4, p: -k * q, q: 1 + k * p, alf: alf_p, Jf: -Jf**3 * (alf * alf_p + U[1] * U[2]),
         alb: alb_p, Jb: -Jb**3 * (alb * alb_p + a[1] * a[2])}
for i in range(6): rules[a[i]] = a[i + 1]; rules[w[i]] = w[i + 1] if i < 5 else sp.Symbol('w6'); rules[Y[i]] = Y[i + 1] if i < 5 else sp.Symbol('Y6')
def Dx(e):
    e = sp.sympify(e)
    return sp.expand(sum(sp.diff(e, s) * r for s, r in rules.items() if e.has(s)))
# partial derivatives with respect to jets (chain rule through alf/Jf resp. alb/Jb)
def dw(e, i):     # d/dw_i , full state
    e = sp.sympify(e); r = sp.diff(e, w[i])
    if i == 0: r += sp.diff(e, alf) * (-k) + sp.diff(e, Jf) * (Jf**3 * alf * k)
    if i == 1: r += sp.diff(e, Jf) * (-Jf**3 * U[1])
    return sp.expand(r)
def da(e, i):     # d/da_i at fixed w ... for the background function Fb(a) (depends on alb,Jb only)
    e = sp.sympify(e); r = sp.diff(e, a[i])
    if i == 0: r += sp.diff(e, alb) * (-k) + sp.diff(e, Jb) * (Jb**3 * alb * k)
    if i == 1: r += sp.diff(e, Jb) * (-Jb**3 * a[1])
    return sp.expand(r)
def Fexpr(al_, Ji_, u):          # u = list of 5 jets
    tot = 0
    for (i, j, l), cf in GS.items():
        tot += cf.subs({oal: al_, oJi: Ji_, ou0: u[0], ou1: u[1]}, simultaneous=True) * u[2]**i * u[3]**j * u[4]**l
    return sp.expand(tot + (p + u[0] - u[1] * q / al_) / 4)
Ff = Fexpr(alf, Jf, U)
Fb = Fexpr(alb, Jb, list(a[:5]))
aL = {4: -1, 3: 0, 2: -k**2, 1: -(5 * k * m + q / 4), 0: sp.Rational(1, 4) - 3 * m**2 - k * p}
L0w = sum(aL[i] * w[i] for i in range(5))
# ---- decomposition of F:  c4 u4 + c3 u3 + G
def split(cfmap, al_, Ji_, u):
    c4 = sp.expand(sum(cf.subs({oal: al_, oJi: Ji_, ou0: u[0], ou1: u[1]}, simultaneous=True) for (i, j, l), cf in cfmap.items() if (i, j, l) == (0, 0, 1)))
    c3 = sp.expand(sum(cf.subs({oal: al_, oJi: Ji_, ou0: u[0], ou1: u[1]}, simultaneous=True) * u[2]**i for (i, j, l), cf in cfmap.items() if j == 1 and l == 0))
    G = sp.expand(sum(cf.subs({oal: al_, oJi: Ji_, ou0: u[0], ou1: u[1]}, simultaneous=True) * u[2]**i for (i, j, l), cf in cfmap.items() if j == 0 and l == 0) + (p + u[0] - u[1] * q / al_) / 4)
    assert all((j == 0 and l == 0) or (j == 1 and l == 0) or (i, j, l) == (0, 0, 1) for (i, j, l) in cfmap), "not quasilinear"
    return c4, c3, G
c4f, c3f, Gf = split(GS, alf, Jf, U)
c4b, c3b, Gb = split(GS, alb, Jb, list(a[:5]))
assert sp.expand(c4f * U[4] + c3f * U[3] + Gf - Ff) == 0
assert not c4f.has(w[2]) and not c3f.has(w[3])
# ---- linear part
dC = {}
for i in range(5):
    dC[i] = sp.expand(da(Fb, i) - aL[i])
def verify():
    out = {}
    # L0 = DF(0): d_i Fb at a=0
    z0 = {s: 0 for s in a}; z0.update({alb: 1, Jb: 1})
    for i in range(5):
        v = sp.expand(da(Fb, i).subs(z0, simultaneous=True) - aL[i]); assert v == 0, (i, v)
        assert sp.expand(dC[i].subs(z0, simultaneous=True)) == 0
    out['L0=DF(0)'] = 0
    # Lin identity
    lin = sum(Y[0] * dC[i] * w[i] for i in range(5))
    Lam = 0; Bl = 0
    for i in range(5):
        g = Y[0] * dC[i]
        gs = [g]
        for s in range(i): gs.append(Dx(gs[-1]))
        Lam += (-1)**i * gs[i]
        Bl += sum((-1)**s * gs[s] * w[i - 1 - s] for s in range(i))
    Lam = sp.expand(Lam)
    res = sp.expand(lin - w[0] * Lam - Dx(Bl)); out['lin_identity_residual'] = res; assert res == 0
    return Lam, Bl, out
Lam, Bl, vout = verify()
print(f"[{time.time()-t00:.0f}s] linear part ok; Lam terms {len(Lam.args)}", flush=True)
# ---- nonlinear remainder
D4 = sp.expand(c4f - c4b); D3 = sp.expand(c3f - c3b)
Rdir = sp.expand(Ff - Fb - sum(da(Fb, i) * w[i] for i in range(5)))              # direct definition (uses d_iFb)
Rdir = sp.expand(Ff - Fb - sum(da(Fb, i) * w[i] for i in range(5)))
Rrest = sp.expand(D4 * a[4] + D3 * a[3] + Gf - Gb - sum(da(Fb, i) * w[i] for i in range(3)))
# note: d_3 Fb = c3b, d_4 Fb = c4b  -> R = D4 w4 + D3 w3 + Rrest
assert sp.expand(da(Fb, 3) - c3b) == 0 and sp.expand(da(Fb, 4) - c4b) == 0
assert sp.expand(Rdir - (D4 * w[4] + D3 * w[3] + Rrest)) == 0
Th = sp.expand(Y[0] * D3 - Dx(Y[0] * D4))
assert not Th.has(w[3])
th0 = Th.subs(w[2], 0); th1 = sp.expand((Th - th0) / w[2]); assert not th1.has(w[2])
P = sp.expand(th0 * w[2] + th1 * w[2]**2 / 2)
assert sp.expand(sp.diff(P, w[2]) - Th) == 0
Bq = sp.expand(Y[0] * D4 * w[3] + P)
Pex = sp.expand(Dx(P) - Th * w[3])
assert not Pex.has(w[3]) and not Pex.has(w[4])
QR = sp.expand(Y[0] * Rrest - Pex)
assert not QR.has(w[3]) and not QR.has(w[4]) and not QR.has(w[5])
resq = sp.expand(Y[0] * Rdir - QR - Dx(Bq)); assert resq == 0
vout['quad_identity_residual'] = resq
print(f"[{time.time()-t00:.0f}s] nonlinear identity ok; QR terms {len(QR.args)}", flush=True)
# order-of-vanishing at w=0 (alf->alb, Jf->Jb)
def at0(e): return sp.expand(e.subs({alf: alb, Jf: Jb, w[0]: 0, w[1]: 0, w[2]: 0}, simultaneous=True))
assert at0(QR) == 0
for i in range(3): assert at0(dw(QR, i)) == 0, i
vout['QR_vanishes_to_second_order'] = True
# structure in w2: QR = sum_j w2^j Q_j(w0,w1)
QRp = sp.Poly(QR, w[2]); degw2 = QRp.degree(); print("deg_w2 =", degw2, flush=True)
Qj = {j: sp.expand(QRp.coeff_monomial(w[2]**j)) for j in range(degw2 + 1)}
# split into Y_r parts
def Ypart(e, r): return sp.expand(sp.diff(e, Y[r]))
for r in range(6):
    if r > 2: assert not QR.has(Y[r]), r
if __name__ == "__main__":
    print(vout)
    print("Lam has Y:", [r for r in range(6) if Lam.has(Y[r])], " a-jets:", [i for i in range(7) if Lam.has(a[i])])
    print("QR has Y:", [r for r in range(6) if QR.has(Y[r])], " a-jets:", [i for i in range(7) if QR.has(a[i])])
    import json; json.dump({k_: str(v) for k_, v in vout.items()}, open('verify_symbolic.json', 'w'), indent=1)
    print(f"[{time.time()-t00:.0f}s] done")
