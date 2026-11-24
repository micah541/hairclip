"""Rigorous Gram matrix from hsrun3 doubleton outputs with SHARED parameter symbols.

Every run represents its final set as  x + C r0 + B r  (CAPD doubleton, no reorganization: C is a point matrix
transported exactly, r0 is the ORIGINAL parameter box).  The columns of C that belong to physical parameters
common to all runs are kept as symbols s_k:
    k = 0..3 : step-1 profile parameters at sigma = 8       (all runs)
    k = 4    : tip parameter st                               (all runs)
    k = 5..9 : cone coordinates of the admissible plane at x0 (Z runs only; same plane in every Z run)
The true solution of every run equals x + C s* + B r' with the SAME s* (true parameter values) and a run-specific
r' in [r].  Products of affine forms keep first-order cancellations between runs; second-order terms and all
run-specific parts go into the remainder.
G[i][j] = <A^{-1} f_j, f_i>_rho (half line) = l_i . (y_j, I_ji),   l_i = (0, m3_i/W13, 0, -m1_i/W13, 1)."""
import sys, json
import numpy as np
from flint import arb, ctx
sys.path.insert(0, '/home/claude/s1spec/step3')
from gen_hs3 import parse3, aux_for, QSHIFT
ctx.prec = 128
GEN = ["p", "C", "S", "w"]
SHARED_ALL = [0, 1, 2, 3, 4]
SHARED_Z = [5, 6, 7, 8, 9]
UNIT = arb(0, 1)

def A(lo, hi):
    assert lo <= hi, (lo, hi)
    return arb((lo + hi) / 2) + UNIT * ((hi - lo) / 2 * (1 + 1e-15) + 1e-300)
def mag(a): return abs(a).upper()          # arb upper bound of |a|

class Aff:
    """c + sum_k d[k] s_k + E,  s_k in S[k] (centred boxes), E in [-e, e];  c, d[k] arbs (balls)"""
    def __init__(self, c, d=None, e=arb(0)):
        self.c = c; self.d = dict(d or {}); self.e = e
    def __add__(self, o):
        if not isinstance(o, Aff): o = Aff(arb(o))
        d = dict(self.d)
        for k, v in o.d.items(): d[k] = d.get(k, arb(0)) + v
        return Aff(self.c + o.c, d, self.e + o.e)
    __radd__ = __add__
    def __neg__(self): return self * arb(-1)
    def __sub__(self, o): return self + (-o)
    def scale(self, a):
        return Aff(self.c * a, {k: v * a for k, v in self.d.items()}, self.e * mag(a))
    def lin_bound(self, S):
        return sum((mag(v) * mag(S[k]) for k, v in self.d.items()), arb(0))
    def __mul__(self, o):
        if not isinstance(o, Aff): return self.scale(arb(o))
        d = {}
        for k, v in self.d.items(): d[k] = d.get(k, arb(0)) + v * o.c
        for k, v in o.d.items(): d[k] = d.get(k, arb(0)) + self.c * v
        La, Lb = self.lin_bound(SYM), o.lin_bound(SYM)
        rem = La * Lb + self.e * (mag(o.c) + Lb + o.e) + o.e * (mag(self.c) + La)
        return Aff(self.c * o.c, d, rem)
    __rmul__ = __mul__
    def interval(self):
        tot = self.c
        for k, v in self.d.items(): tot += v * SYM[k]
        return tot + UNIT * self.e
    def tidy(self):
        """move the radius of the coefficient balls into the remainder"""
        d, e = {}, self.e
        for k, v in self.d.items():
            m = arb(v.mid()); d[k] = m; e += mag(v - m) * mag(SYM[k])
        cm = arb(self.c.mid()); e += mag(self.c - cm)
        return Aff(cm, d, e)

SYM = None
class DSet:
    def __init__(self, path, shared):
        global SYM
        D = parse3(open(path).read()); assert D["ok"], path
        self.x = [A(*v) for v in D["TX"]]
        self.C = [[A(*v) for v in row] for row in D["TC"]]
        self.B = [[A(*v) for v in row] for row in D["TB"]]
        self.r0 = [A(*v) for v in D["TR0"]]
        self.r = [A(*v) for v in D["TR"]]
        self.hull = [A(*v) for v in D["TIP"]]
        self.n = len(self.x); self.shared = shared
        syms = {k: self.r0[k] for k in shared}
        if SYM is None: SYM = {}
        for k, v in syms.items():                       # all runs must use the same parameter boxes
            if k in SYM: assert (SYM[k] - v).contains(0) and abs(float(SYM[k].rad()) - float(v.rad())) <= 1e-12 * float(v.rad()) + 1e-300, (path, k)
            else: SYM[k] = v
    def form(self, coef):
        """affine form of sum_i coef[i] z_i (coef: row -> arb)"""
        c = sum((a * self.x[i] for i, a in coef.items()), arb(0))
        d, e = {}, arb(0)
        for k in range(self.n):
            gC = sum((a * self.C[i][k] for i, a in coef.items()), arb(0))
            gB = sum((a * self.B[i][k] for i, a in coef.items()), arb(0))
            if k in self.shared: d[k] = gC
            else: e += mag(gC) * mag(self.r0[k])
            e += mag(gB) * mag(self.r[k])
        return Aff(c, d, e).tidy()
    def ratio(self, im, iw):
        """affine form of z_im / z_iw (mean-value form; gradient over the hull)"""
        Hm, Hw = self.hull[im], self.hull[iw]
        gm, gw = 1 / Hw, -Hm / (Hw * Hw)
        c = self.x[im] / self.x[iw]
        d, e = {}, arb(0)
        for k in range(self.n):
            gC = gm * self.C[im][k] + gw * self.C[iw][k]
            gB = gm * self.B[im][k] + gw * self.B[iw][k]
            if k in self.shared: d[k] = gC
            else: e += mag(gC) * mag(self.r0[k])
            e += mag(gB) * mag(self.r[k])
        return Aff(c, d, e).tidy()

def na(gs): return len(aux_for(gs))

def isect(a, b):
    lo = max(a.lower(), b.lower()); hi = min(a.upper(), b.upper())
    assert lo <= hi, "empty intersection: inconsistent enclosures"
    return (lo + hi) / 2 + UNIT * ((hi - lo) / 2)

def gram(prefix="out7", verbose=True):
    rho1, rho3 = {}, {}
    for gi in GEN:
        S = DSet(f"{prefix}Z{gi}.txt", SHARED_ALL + SHARED_Z); o = 5 + na([gi])
        rho1[gi] = S.ratio(o + 7, o + 4); rho3[gi] = S.ratio(o + 9, o + 4)
        if verbose: print(f"Z{gi}: m1/W13 = {rho1[gi].interval()}  m3/W13 = {rho3[gi].interval()}")
    G4f = [[None] * 4 for _ in range(4)]
    for a, gi in enumerate(GEN):
        for b, gj in enumerate(GEN):
            S = DSet(f"{prefix}Q{gj}{gi}.txt", SHARED_ALL); o = 5 + na([gj, gi])
            # l = (0, rho3, 0, -rho1, 1): centre part combined across rows first (keeps the cancellation of the
            # run-local errors along the admissible plane), small parts (l - l_centre) multiplied afterwards
            c1, c3 = rho1[gi].c, rho3[gi].c
            main = S.form({o + 1: c3, o + 3: -c1, o + 4: arb(1)})
            d1 = Aff(arb(0), rho1[gi].d, rho1[gi].e); d3 = Aff(arb(0), rho3[gi].d, rho3[gi].e)
            y1 = S.form({o + 1: arb(1)}); y3 = S.form({o + 3: arb(1)})
            G4f[a][b] = (main + d3 * y1 - d1 * y3).scale(arb(2) ** QSHIFT)
    Sc = DSet(f"{prefix}ZC.txt", SHARED_ALL + SHARED_Z)
    Ct, St = Sc.hull[5], Sc.hull[6]
    c1, c2 = arb(float(Ct.mid())), arb(float(St.mid()))
    Rm = [[arb(1), arb(0), arb(0), arb(0)], [arb(0), c1, c2, arb(0)], [arb(0), arb(0), arb(0), arb(1)]]
    G3f = [[None] * 3 for _ in range(3)]
    for a in range(3):
        for b in range(3):
            acc = Aff(arb(0))
            for i in range(4):
                for j in range(4):
                    w = Rm[a][i] * Rm[b][j]
                    if w != 0: acc = acc + G4f[i][j].scale(w)
            G3f[a][b] = acc
    G3 = [[isect(G3f[a][b].interval(), G3f[b][a].interval()) for b in range(3)] for a in range(3)]
    G4 = [[isect(G4f[a][b].interval(), G4f[b][a].interval()) for b in range(4)] for a in range(4)]
    if verbose:
        print("G4 (shared symbols, symmetric intersection):")
        for a in range(4): print("  " + "  ".join(f"{float(G4[a][b].mid()):+.7e}(+-{float(G4[a][b].rad()):.1e})" for b in range(4)))
        print("f2 = c1 cos K + c2 sin K,  c =", float(c1), float(c2), "  (C_tip, S_tip =", Ct, St, ")")
        for a in range(3): print("  G3:", [f"{float(G3[a][b].mid()):+.7e}(+-{float(G3[a][b].rad()):.1e})" for b in range(3)])
    g11, g12, g22, g13, g23, g33 = G3[0][0], G3[0][1], G3[1][1], G3[0][2], G3[1][2], G3[2][2]
    det2 = g11 * g22 - g12 * g12
    Sch = g33 - (g22 * g13 * g13 - 2 * g12 * g13 * g23 + g11 * g23 * g23) / det2
    ok = bool(g11 < 0 and det2 > 0 and Sch < 0)
    if verbose:
        print("G2: g11 =", g11, " det =", det2)
        print(f"Schur complement S in [{float(Sch.lower()):.4e}, {float(Sch.upper()):.4e}]")
        print("G3 NEGATIVE DEFINITE:", ok)
    return dict(G3=G3, G4=G4, S=Sch, ok=ok, c=(float(c1), float(c2)))

if __name__ == "__main__":
    gram(sys.argv[1] if len(sys.argv) > 1 else "out7")
