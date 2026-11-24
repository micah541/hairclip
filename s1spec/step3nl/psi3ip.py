"""Rigorous inner products (half line, L^2(rho)) of psi3^ with the certified constraints of the linear certificate:
f1 = chi p, f2 = chi psi2 (psi2 = Ct C + St S), f3 = chi F - a eta1 - b eta3 (F = ws0), from the hulls of the
Lyapunov run outA5.txt (3425 boxes); plus ||psi3^|| on a fine grid. Each box contributes len * hull(integrand)."""
import sys, json
from flint import arb
sys.path.insert(0, '/home/claude/s1spec/step3'); sys.path.insert(0, '/home/claude/s1spec/step3nl')
import psi3fun as F
from gen_hs2 import chunks
from gram_post7 import DSet
import gen_hs3 as G3
U = arb(0, 1)
def A(lo, hi): return arb((lo + hi) / 2) + U * ((hi - lo) / 2 * (1 + 1e-15) + 1e-300)
def hull(a, b):
    lo = min(a.lower(), b.lower()); hi = max(a.upper(), b.upper()); return (lo + hi) / 2 + U * ((hi - lo) / 2)
def rho_iv(X): return (-(1 + X * X) ** (arb(2) / 3) / 4).exp()
def psi_iv(X):
    xc = float(X.mid()); out = None
    for P in F.PIECES:
        if X.upper() < P.a.lower() or X.lower() > P.b.upper(): continue
        jc = F.jet_poly(P, arb(xc), 3); jX = F.jet_poly(P, X, 3); hh = X - arb(xc)
        v = jc[0] + jc[1] * hh + jX[2] * hh * hh / 2 + F.Hcorr(P, X, 1)[0]
        out = v if out is None else out.union(v)
    if X.upper() > F.XC.lower():
        v = F.ext_jet(X, 1)[0]; out = v if out is None else out.union(v)
    return out
# smoothstep chi on sigma in [29, 32]
def S_(u): return u ** 4 * (35 - 84 * u + 70 * u * u - 20 * u ** 3)
T0 = arb(G3.T0); CH = chunks()
fin = DSet('/home/claude/s1spec/step3/outA5.txt', [0, 1, 2, 3, 4])
C8h, S8h, Ct2h, St2h = fin.hull[5], fin.hull[6], fin.hull[7], fin.hull[8]
CtB = C8h * Ct2h - S8h * St2h; StB = S8h * Ct2h + C8h * St2h
I = {"f1": arb(0), "f2": arb(0), "F": arb(0), "absbump": arb(0)}
fl = {"f1": 0.0, "f2": 0.0, "F": 0.0}
for line in open('/home/claude/s1spec/step3/outA5.txt'):
    tag = line[:2]
    if tag not in ("MB", "MC"): continue
    v = list(map(float, line.split()[1:])); ta, tb = v[0], v[1]
    h = [A(v[2 + 2 * i], v[3 + 2 * i]) for i in range((len(v) - 2) // 2)]
    if tag == "MB":
        st, C_, S_v, ws0 = h[0], h[1], h[2], h[5]
        Ct, St = CtB, StB
        c = min(int(ta + 1e-9), len(CH) - 1); pb = CH[c][2][2]
        pv = A(pb[0], pb[1])
        sig = A(32 - tb, 32 - ta); X = sig - st; ln = arb(tb - ta)
        chi = S_(A(ta / 3, tb / 3)) if tb <= 3.0 + 1e-12 else arb(1)
    else:
        pv, st, C_, S_v, ws0 = h[0], h[4], h[7], h[8], h[9]
        Ct, St = Ct2h, St2h
        cc = (8 - st) / T0
        X = hull(cc * (T0 - arb(tb)), cc * (T0 - arb(ta))); ln = cc * arb(tb - ta)
        chi = arb(1)
    Xc = arb(max(float(X.lower()), 0.0), 0) if False else X
    lo = max(float(X.lower()), 0.0); X = arb((lo + float(X.upper())) / 2, (float(X.upper()) - lo) / 2)
    u = psi_iv(X); w = ln * rho_iv(X) * u * chi
    psi2 = Ct * C_ + St * S_v
    for key, val in (("f1", pv), ("f2", psi2), ("F", ws0)):
        I[key] += w * val; fl[key] += float((w * val).mid())
    if tag == "MC" and ta >= float(T0) - G3.LB - 1e-12:
        I["absbump"] += ln * rho_iv(X) * abs(u)
a_, b_ = arb(0, 9.1e-8), arb(0, 2.6e-7)
bump = I["absbump"] * (abs(a_) * arb("0.833") + abs(b_) * arb("0.425"))   # |eta1| <= 4 max s zeta = 0.8324, |eta3| <= (32/3)max s^3 zeta = 0.4245
I3 = I["F"] + arb(0, float(bump.upper()))
# ||psi3^||^2 on [0, 40] with boxes of width 1e-3, plus tail (psi^2 <= (1.1 c1 s^A)^2)
nrm = arb(0); N = 40000
for i in range(N):
    X = arb((i + 0.5) * 1e-3, 0.5e-3); u = psi_iv(X); nrm += arb(1e-3) * rho_iv(X) * u * u
x1 = arb(40); s1 = x1 + F.STIP
tailn = (abs(F.C1) * 1.1) ** 2 * x1 ** (2 * F.A) * (-(x1 ** (arb(4) / 3)) / 4).exp() / (x1 ** (arb(1) / 3) / 3 - 2 * F.A / x1)
nrm += arb(0, float(tailn.upper()))
res = {k: [float(v.lower()), float(v.upper())] for k, v in (("psi3_f1", I["f1"]), ("psi3_f2", I["f2"]), ("psi3_chiF", I["F"]), ("psi3_f3", I3), ("psi3_norm2", nrm))}
res["float_mid"] = fl
for k, v in res.items(): print(f"  {k:12s} {v}")
json.dump(res, open('psi3ip.json', 'w'), indent=1)
