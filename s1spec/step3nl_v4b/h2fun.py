"""[h2fun: generated from psi3fun.py] The explicit second-order correction h2^ (same construction as psi3^ below,
with lamt -> 2 lamt, A = 1 - 8 lamt, data h2md_14.json; h2^(0) is NOT normalised).
Original docstring: The explicit approximate eigenfunction psi3^ (C^3, even, psi3^(0) = 1 up to a correction ~1e-12).
psi3^ = P_e + H_e on element e of psi3md2_14.json (float Chebyshev P_e in y = 1 - 2(x-a)/(b-a)), and
psi3^ = c1 g(x + stip) + c2 E(x + stip) + H_inf for x >= X, where
  g(s) = sum_{j<4} a_j s^(A-4j)   (A = 1 - 4 lamt, a_j by the recursion; R g = -4 a_4 s^(A-16) exactly),
  E(s) = exp(-3 (s/4)^(4/3) + 3 (Xs/4)^(4/3))   (Xs = X + stip; E' = -(s/4)^(1/3) E).
H_e are degree-7 Hermite corrections (exact real coefficients enclosed in arb) that remove the jet mismatches
(u, u', u'', u''') at every interface, impose u'(0) = u'''(0) = 0, and match the extension at X. H_inf = 0
(the far-end mismatch is absorbed in the last element). So psi3^ is C^3, piecewise smooth, in H^4_loc."""
import json, numpy as np
from flint import arb, arb_poly
from numpy.polynomial import chebyshev as Ch
import sys
sys.path.insert(0, '/home/claude/s1spec/step3')
from gen_hs3 import LAMT
J = json.load(open('/home/claude/s1spec/step3nl/h2md_14.json'))
STIP = arb(J["stip"]); XC = arb(J["X"]); LAM = 2 * arb(LAMT); A = 1 - 4 * LAM
C1, C2 = arb(J["c_poly"]), arb(J["c_dec"])
AJ = [arb(1)]
for j in range(1, 5): AJ.append(AJ[-1] * (A - 4 * j + 4) * (A - 4 * j + 3) * (A - 4 * j + 2) * (A - 4 * j + 1) / j)
NDER = 8
def clen(cf, y):
    b1, b2 = arb(0), arb(0)
    for ck in cf[:0:-1]: b1, b2 = 2 * y * b1 - b2 + ck, b1
    return y * b1 - b2 + cf[0]
class Piece:
    def __init__(s, a, b, coef):
        s.a, s.b = arb(a), arb(b); s.af, s.bf = a, b
        cf = np.array(coef); s.cfs = []
        for n in range(NDER + 1):
            s.cfs.append([arb(float(c)) for c in (Ch.chebder(cf, n) if n else cf)])   # chebder of floats: rounding!
        s.sc = -2 / (s.b - s.a); s.H = None
    def y(s, x): return 1 + s.sc * (x - s.a)
    def jet(s, x, n=4):
        """exact-ish derivatives at a point (arb x): uses the float chebder coefficients -> must be exact; see note"""
        yy = s.y(x); return [clen(s.cfs[d], yy) * s.sc ** d for d in range(n)]
# NOTE: numpy chebder on float coefficients is not exact. To make P_e a single well-defined polynomial we
# differentiate in arb instead (cheb derivative recursion), below.
def cheb_der_arb(c):
    n = len(c) - 1
    if n == 0: return [arb(0)]
    d = [arb(0)] * n
    d[n - 1] = 2 * n * c[n]
    if n >= 2: d[n - 2] = 2 * (n - 1) * c[n - 1]
    for k in range(n - 3, -1, -1): d[k] = d[k + 2] + 2 * (k + 1) * c[k + 1]
    d[0] = d[0] / 2
    return d
PIECES = []
for pc in J["pieces"]:
    P = Piece(pc["a"], pc["b"], pc["coef"])
    c0 = [arb(float(c)) for c in pc["coef"]]; P.cfs = [c0]
    for n in range(NDER): P.cfs.append(cheb_der_arb(P.cfs[-1]))
    PIECES.append(P)
# ---------------------------------------------------------------- extension
def ext_jet(x, n=5):
    """derivatives 0..n-1 of c1 g(s) + c2 E(s), s = x + stip (x: arb)"""
    s = x + STIP; out = []
    for d in range(n):
        v = arb(0)
        for j in range(4):
            e = A - 4 * j; f = arb(1)
            for i in range(d): f = f * (e - i)
            v += AJ[j] * f * s ** (e - d)
        out.append(C1 * v)
    Xs = XC + STIP
    E = (-3 * (s / 4) ** (arb(4) / 3) + 3 * (Xs / 4) ** (arb(4) / 3)).exp()
    z = -(s / 4) ** (arb(1) / 3); z1 = -(s / 4) ** (arb(-2) / 3) / 12; z2 = (s / 4) ** (arb(-5) / 3) / 72
    z3 = -5 * (s / 4) ** (arb(-8) / 3) / 864
    Ed = [arb(1), z, z1 + z * z, z2 + 3 * z * z1 + z ** 3,
          z3 + 4 * z * z2 + 3 * z1 * z1 + 6 * z * z * z1 + z ** 4]   # (e^S)^(d) / e^S
    return [o + C2 * E * Ed[d] for d, o in enumerate(out)]
# ---------------------------------------------------------------- Hermite corrections
from math import comb, factorial
def hermite_basis():
    """degree-7 polys phi_d (left jet e_d, right jet 0) and chi_d (left 0, right e_d) on h in [0,1], monomial coeffs (exact ints/rationals as arb)"""
    import sympy as sp
    h = sp.symbols('h'); cs = sp.symbols('c0:8'); P = sum(cs[i] * h ** i for i in range(8)); out = {}
    for side in ('L', 'R'):
        for d in range(4):
            eqs = []
            for k in range(4):
                vL = sp.diff(P, h, k).subs(h, 0); vR = sp.diff(P, h, k).subs(h, 1)
                eqs.append(vL - (1 if (side == 'L' and k == d) else 0)); eqs.append(vR - (1 if (side == 'R' and k == d) else 0))
            sol = sp.solve(eqs, cs)
            out[(side, d)] = [sp.Rational(sol[c]) for c in cs]
    return out
HB = hermite_basis()
def mono_eval(cf, hval, k):
    """k-th derivative of sum cf[i] h^i at hval (arb), cf exact rationals -> arb"""
    tot = arb(0)
    for i in range(k, len(cf)):
        f = 1
        for t in range(k): f *= (i - t)
        if f: tot += arb(cf[i].p * f) / arb(cf[i].q) * hval ** (i - k)
    return tot
def build_corrections():
    # interface mismatches
    for e, P in enumerate(PIECES): P.Hl = [arb(0)] * 4; P.Hr = [arb(0)] * 4
    j0 = PIECES[0].jet(arb(0))
    PIECES[0].Hl = [arb(0), -j0[1], arb(0), -j0[3]]          # make u'(0) = u'''(0) = 0 (u, u'' untouched)
    jumps = []
    for e in range(len(PIECES) - 1):
        L, R = PIECES[e], PIECES[e + 1]
        x = R.a
        jl = L.jet(x); jr = R.jet(x)
        dl = [jl[k] - jr[k] for k in range(4)]                # add to right piece's left jet
        R.Hl = dl; jumps.append(max(float(abs(t).upper()) for t in dl))
    Lp = PIECES[-1]; jl = Lp.jet(XC); je = ext_jet(XC, 4)
    Lp.Hr = [je[k] - jl[k] for k in range(4)]
    jumps.append(max(float(abs(t).upper()) for t in Lp.Hr))
    return jumps
JUMPS = build_corrections()
def Hcorr(P, x, n):
    """derivatives 0..n-1 (in x) of the correction on piece P at x (arb), via h = (x - a)/(b - a)"""
    L = P.b - P.a; hv = (x - P.a) / L; out = []
    for k in range(n):
        v = arb(0)
        for d in range(4):
            if P.Hl[d] != 0: v += P.Hl[d] * L ** d * mono_eval(HB[('L', d)], hv, k)
            if P.Hr[d] != 0: v += P.Hr[d] * L ** d * mono_eval(HB[('R', d)], hv, k)
        out.append(v / L ** k)
    return out
def piece_of(xf):
    for P in PIECES:
        if P.af <= xf <= P.bf: return P
    return None
def jet_poly(P, X, n=5):
    yy = P.y(X); return [clen(P.cfs[d], yy) * P.sc ** d for d in range(n)]
if __name__ == "__main__":
    print("max jet mismatch per interface (incl. tip handled separately, last = far end):", ["%.1e" % j for j in JUMPS])
    j0 = PIECES[0].jet(arb(0)); print("tip jet (u, u', u'', u'''):", [float(t.mid()) for t in j0])
    # check C^3 after correction
    worst = 0
    for e in range(len(PIECES) - 1):
        L, R = PIECES[e], PIECES[e + 1]; x = R.a
        jl = [a + b for a, b in zip(L.jet(x), Hcorr(L, x, 4))]; jr = [a + b for a, b in zip(R.jet(x), Hcorr(R, x, 4))]
        worst = max(worst, max(float(abs(a - b).upper()) for a, b in zip(jl, jr)))
    Lp = PIECES[-1]; jl = [a + b for a, b in zip(Lp.jet(XC), Hcorr(Lp, XC, 4))]; je = ext_jet(XC, 4)
    worst = max(worst, max(float(abs(a - b).upper()) for a, b in zip(jl, je)))
    j0c = [a + b for a, b in zip(PIECES[0].jet(arb(0)), Hcorr(PIECES[0], arb(0), 4))]
    print("after correction: max |jump| enclosure", worst, "; tip u', u''' =", j0c[1], j0c[3])
