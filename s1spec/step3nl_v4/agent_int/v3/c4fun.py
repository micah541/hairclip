"""C^4 versions psi3~, h2~ of the certified C^3 functions psi3^ (psi3fun.py) and h2^ (h2fun.py).
The raw Chebyshev pieces P_e and the explicit far-field extension are those of psi3fun/h2fun; only the Hermite corrections change:
degree-9 corrections (exact rational basis, arb coefficients) remove the jet mismatches of orders 0..4 (the C^3 versions removed 0..3) at every
interface and at the extension point X = 28, and impose u'(0) = u'''(0) = 0.  Reason: the C^3 versions have u'''' jumps of 1e-7..4e-6 at every
interface, so (L psi3^)' etc. contain delta masses and the derivative-level (H^5) energy estimate is not defined; the C^4 versions have bounded
piecewise u^(5).  extension jets are available to any order (Bell recursion for (e^S)^(d)).  Jets: exact arb enclosures on balls."""
import sys, math
sys.path.insert(0, '/home/claude/s1spec/step3nl'); sys.path.insert(0, '/home/claude/s1spec/step3')
from flint import arb
import sympy as sp
import psi3fun as F3
import h2fun as H3
NJET = 5          # matched jets 0..4
def _hermite_basis():
    h = sp.symbols('h'); cs = sp.symbols('c0:%d' % (2 * NJET)); P = sum(cs[i] * h ** i for i in range(2 * NJET)); out = {}
    for side in ('L', 'R'):
        for d in range(NJET):
            eqs = []
            for k in range(NJET):
                vL = sp.diff(P, h, k).subs(h, 0); vR = sp.diff(P, h, k).subs(h, 1)
                eqs.append(vL - (1 if (side == 'L' and k == d) else 0)); eqs.append(vR - (1 if (side == 'R' and k == d) else 0))
            sol = sp.solve(eqs, cs); out[(side, d)] = [sp.Rational(sol[c]) for c in cs]
    return out
HB = _hermite_basis()
def mono_eval(cf, hval, k):
    tot = arb(0)
    for i in range(k, len(cf)):
        f = 1
        for t in range(k): f *= (i - t)
        if f: tot += arb(cf[i].p * f) / arb(cf[i].q) * hval ** (i - k)
    return tot
# ---- far-field extension jets of any order:  c1 g(s) + c2 E(s),  E = exp(-3 (s/4)^(4/3) + 3 (Xs/4)^(4/3)),  (e^S)^(d)/e^S = Ed_d
_NZ = 10
_zs = sp.symbols('z0:%d' % _NZ)
_P = [sp.Integer(1)]
for d in range(8):
    e = _zs[0] * _P[-1] + sum(sp.diff(_P[-1], _zs[i]) * _zs[i + 1] for i in range(_NZ - 1))
    _P.append(sp.expand(e))
class _Pr(sp.printing.str.StrPrinter):
    def _print_Rational(self, e): return f"(arb({e.p})/arb({e.q}))"
    def _print_Integer(self, e): return f"arb({e.p})"
    def _print_Pow(self, e, rational=False):
        b, ex = e.as_base_exp(); n = int(ex); return "(" + "*".join([f"({self._print(b)})"] * n) + ")"
EDF = [eval("lambda " + ",".join(map(str, _zs)) + ": " + _Pr().doprint(p_), {"arb": arb}) for p_ in _P]
def _zjets(s):
    out = []
    for j in range(_NZ):
        c = arb(1)
        for i in range(j): c = c * (arb(1) / 3 - i)
        out.append(-c / arb(4) ** j * (s / 4) ** (arb(1) / 3 - j))
    return out
def ext_jets(M, x, n):
    """derivatives 0..n-1 (n <= 9) of the extension of module M at x (arb)"""
    s = x + M.STIP; out = []
    for d in range(n):
        v = arb(0)
        for j in range(4):
            e = M.A - 4 * j; f = arb(1)
            for i in range(d): f = f * (e - i)
            v += M.AJ[j] * f * s ** (e - d)
        out.append(M.C1 * v)
    Xs = M.XC + M.STIP
    E = (-3 * (s / 4) ** (arb(4) / 3) + 3 * (Xs / 4) ** (arb(4) / 3)).exp()
    z = _zjets(s)
    return [o + M.C2 * E * EDF[d](*z) for d, o in enumerate(out)]
class C4:
    def __init__(self, M):
        self.M = M; self.pieces = M.PIECES; n = len(self.pieces)
        self.Hl = [[arb(0)] * NJET for _ in range(n)]; self.Hr = [[arb(0)] * NJET for _ in range(n)]
        j0 = self.raw(0, arb(0), NJET)
        self.Hl[0] = [arb(0), -j0[1], arb(0), -j0[3], arb(0)]
        self.jumps = []
        for e in range(n - 1):
            x = self.pieces[e + 1].a
            jl = self.raw(e, x, NJET); jr = self.raw(e + 1, x, NJET)
            self.Hl[e + 1] = [jl[k] - jr[k] for k in range(NJET)]
            self.jumps.append(max(float(abs(t).upper()) for t in self.Hl[e + 1]))
        jl = self.raw(n - 1, M.XC, NJET); je = ext_jets(M, M.XC, NJET)
        self.Hr[n - 1] = [je[k] - jl[k] for k in range(NJET)]
        self.jumps.append(max(float(abs(t).upper()) for t in self.Hr[n - 1]))
    def raw(self, e, x, n):
        P = self.pieces[e]; yy = P.y(x)
        return [F3.clen(P.cfs[d], yy) * P.sc ** d for d in range(n)]
    def corr(self, e, x, n):
        P = self.pieces[e]; L = P.b - P.a; hv = (x - P.a) / L; out = []
        for k in range(n):
            v = arb(0)
            for d in range(NJET):
                if self.Hl[e][d] != 0: v += self.Hl[e][d] * L ** d * mono_eval(HB[('L', d)], hv, k)
                if self.Hr[e][d] != 0: v += self.Hr[e][d] * L ** d * mono_eval(HB[('R', d)], hv, k)
            out.append(v / L ** k)
        return out
    def piece_jets(self, e, x, n):
        r = self.raw(e, x, n); c = self.corr(e, x, n); return [a + b for a, b in zip(r, c)]
    def u_box(self, X, xc, n):
        """n <= 7 jets on the ball X (Taylor form of the polynomial part about xc; corrections evaluated directly)"""
        out = None; xca = arb(xc); hh = X - xca
        for e, P in enumerate(self.pieces):
            if X.upper() < P.a.lower() or X.lower() > P.b.upper(): continue
            jc = self.raw(e, xca, n + 2); jX = self.raw(e, X, n + 2); Hc = self.corr(e, X, n)
            vals = [jc[d] + jc[d + 1] * hh + jX[d + 2] * hh * hh / 2 + Hc[d] for d in range(n)]
            out = vals if out is None else [o.union(v) for o, v in zip(out, vals)]
        return out
    def jets_any(self, X, xc, n):
        if X.upper() <= float(self.M.XC.lower()): return self.u_box(X, xc, n)
        if X.lower() >= float(self.M.XC.upper()): return ext_jets(self.M, X, n)
        return [a.union(b) for a, b in zip(self.u_box(X, xc, n), ext_jets(self.M, X, n))]
PSI = C4(F3); H2 = C4(H3)
if __name__ == "__main__":
    print("max jet mismatch removed per interface (psi3):", ["%.1e" % j for j in PSI.jumps])
    print("max jet mismatch removed per interface (h2):  ", ["%.1e" % j for j in H2.jumps])
    for name, Cf in (("psi3", PSI), ("h2", H2)):
        worst = 0.0
        for e in range(len(Cf.pieces) - 1):
            x = Cf.pieces[e + 1].a
            a = Cf.piece_jets(e, x, NJET); b = Cf.piece_jets(e + 1, x, NJET)
            worst = max(worst, max(float(abs(u - v).upper()) for u, v in zip(a, b)))
        a = Cf.piece_jets(len(Cf.pieces) - 1, Cf.M.XC, NJET); b = ext_jets(Cf.M, Cf.M.XC, NJET)
        worst = max(worst, max(float(abs(u - v).upper()) for u, v in zip(a, b)))
        j0 = Cf.piece_jets(0, arb(0), 4)
        print(name, "C^4 check: max jump of jets 0..4 at all interfaces + extension point:", worst, "; u'(0), u'''(0):", j0[1], j0[3])
        # compare with the C^3 function: jets 0..3 on an interior point and the 5th jet size
        x = arb(15.3); jn = Cf.piece_jets([i for i, P in enumerate(Cf.pieces) if P.af <= 15.3 <= P.bf][0], x, 6)
        old = F3.jet_poly if False else None
    # compare extension jets with module ext_jet
    X = arb(40.0)
    for M in (F3, H3):
        a = ext_jets(M, X, 5); b = M.ext_jet(X, 5)
        print("ext_jets vs ext_jet (orders 0..4) max diff:", max(float(abs(u - v).upper()) for u, v in zip(a, b)))
    # difference between C^4 and C^3 functions (order 0..4) at several points
    import h2res as HR
    for name, Cf, M in (("psi3", PSI, F3), ("h2", H2, H3)):
        d = []
        for xv in (3.3, 12.2, 21.7, 27.5):
            X = arb(xv); a = Cf.jets_any(X, xv, 5); b = HR.jets_any(M, X, xv)
            d.append(max(float(abs(u - v).upper()) / max(1.0, float(abs(v).upper())) for u, v in zip(a, b)))
        print(name, "max relative |C4 - C3| over jets 0..4 at x = 3.3, 12.2, 21.7, 27.5:", ["%.1e" % t for t in d])
