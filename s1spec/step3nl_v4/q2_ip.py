"""Rigorous <chi Q2(psi3^), f_j>_rho (half line), j = 1, 2, 3, on the 3425 hull boxes of outA5.txt, and the norm
||chi Q2||.  Q2 = quadratic part of the graph equation at u = psi3^ (exact sympy expression, graph_eq.pkl)."""
import sys, json, pickle, sympy as sp
from flint import arb
sys.path.insert(0, '/home/claude/s1spec/step3'); sys.path.insert(0, '/home/claude/s1spec/step3nl')
import psi3fun as F
from gen_hs2 import chunks
from gram_post7 import DSet
import gen_hs3 as G3
G = pickle.load(open('graph_eq.pkl', 'rb')); quad = G['quad']
s = sp.symbols('sigma'); U = sp.Function('U')(s)
pf, qf, kf, mf = [sp.Function(n)(s) for n in ('p', 'q', 'k', 'm')]
Ud = [U] + [sp.Derivative(U, (s, n)) for n in range(1, 5)]
J = sp.symbols('u0:5'); P_ = sp.symbols('pp qq kk mm')
expr = quad
for n in range(4, -1, -1): expr = expr.subs(Ud[n], J[n])
expr = sp.expand(expr.subs({pf: P_[0], qf: P_[1], kf: P_[2], mf: P_[3]}))
class _Pr(sp.printing.str.StrPrinter):
    def _print_Rational(self, e): return f"(arb({e.p})/arb({e.q}))"
    def _print_Integer(self, e): return f"arb({e.p})"
    def _print_Pow(self, e, rational=False):
        b, ex = e.as_base_exp(); n = int(ex); bb = f"({self._print(b)})"
        return "(" + "*".join([bb] * n) + ")" if n > 0 else "(1/(" + "*".join([bb] * (-n)) + "))"
Q2f = eval("lambda u0,u1,u2,u3,u4,pp,qq,kk,mm: " + _Pr().doprint(expr), {"arb": arb})
U_ = arb(0, 1)
def A(lo, hi): return arb((lo + hi) / 2) + U_ * ((hi - lo) / 2 * (1 + 1e-15) + 1e-300)
def hull(a, b):
    lo = min(a.lower(), b.lower()); hi = max(a.upper(), b.upper()); return (lo + hi) / 2 + U_ * ((hi - lo) / 2)
def rho_iv(X): return (-(1 + X * X) ** (arb(2) / 3) / 4).exp()
def jets_iv(X):
    xc = float(X.mid()); out = None
    for P in F.PIECES:
        if X.upper() < P.a.lower() or X.lower() > P.b.upper(): continue
        jc = F.jet_poly(P, arb(xc), 7); jX = F.jet_poly(P, X, 7); hh = X - arb(xc); Hc = F.Hcorr(P, X, 5)
        v = [jc[d] + jc[d + 1] * hh + jX[d + 2] * hh * hh / 2 + Hc[d] for d in range(5)]
        out = v if out is None else [o.union(t) for o, t in zip(out, v)]
    if X.upper() > F.XC.lower():
        v = F.ext_jet(X, 5); out = v if out is None else [o.union(t) for o, t in zip(out, v)]
    return out
def S_(u): return u ** 4 * (35 - 84 * u + 70 * u * u - 20 * u ** 3)
T0 = arb(G3.T0); CH = chunks()
fin = DSet('/home/claude/s1spec/step3/outA5.txt', [0, 1, 2, 3, 4])
C8h, S8h, Ct2h, St2h = fin.hull[5], fin.hull[6], fin.hull[7], fin.hull[8]
CtB = C8h * Ct2h - S8h * St2h; StB = S8h * Ct2h + C8h * St2h
I = {"f1": arb(0), "f2": arb(0), "F": arb(0), "nrm": arb(0)}
for line in open('/home/claude/s1spec/step3/outA5.txt'):
    tag = line[:2]
    if tag not in ("MB", "MC"): continue
    v = list(map(float, line.split()[1:])); ta, tb = v[0], v[1]
    h = [A(v[2 + 2 * i], v[3 + 2 * i]) for i in range(len(v) // 2 - 1)]
    if tag == "MB":
        st, C_, S_v, ws0 = h[0], h[1], h[2], h[5]
        c = min(int(ta + 1e-9), len(CH) - 1); bnd = CH[c][2]
        kv, mv, pv = [A(bnd[i][0], bnd[i][1]) for i in (0, 1, 2)]
        sig = A(32 - tb, 32 - ta); X = sig - st; qv = sig + A(bnd[3][0], bnd[3][1]); ln = arb(tb - ta)
        chi = S_(A(ta / 3, tb / 3)) if tb <= 3.0 + 1e-12 else arb(1); Ct, St = CtB, StB
    else:
        pv, qv, kv, mv, st, C_, S_v, ws0 = h[0], h[1], h[2], h[3], h[4], h[7], h[8], h[9]
        cc = (8 - st) / T0
        X = hull(cc * (T0 - arb(tb)), cc * (T0 - arb(ta))); ln = cc * arb(tb - ta); chi = arb(1); Ct, St = Ct2h, St2h
    lo = max(float(X.lower()), 0.0); X = arb((lo + float(X.upper())) / 2, (float(X.upper()) - lo) / 2)
    j = jets_iv(X)
    q2 = Q2f(*j, pv, qv, kv, mv)
    w = ln * rho_iv(X) * chi * q2
    I["f1"] += w * chi * pv; I["f2"] += w * chi * (Ct * C_ + St * S_v); I["F"] += w * chi * ws0
    I["nrm"] += ln * rho_iv(X) * chi * chi * q2 * q2
res = {k: [float(v.lower()), float(v.upper())] for k, v in I.items()}
for k, v in res.items(): print(f"  <chi Q2, {k}>: {v}")
json.dump(res, open('q2_ip.json', 'w'), indent=1)
