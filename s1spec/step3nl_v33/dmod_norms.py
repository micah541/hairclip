"""Rigorous ||D_rho f_j||_rho (half line), D_rho g = g'''' - 2 phi' g''' + (phi'^2 - phi'') g'' = rho^-1 (rho g'')'',
for f1 = chi p and f2 = chi cos(theta - theta0) (any fixed theta0: |A cos + B sin| <= sqrt(A^2 + B^2)).
Derivatives from the profile ODE (theta' = k); chi = 1 for sigma <= 29; on sigma in [29, 32] the chi-derivatives are
included with |chi^(j)| bounds of the C^3 smoothstep (u = (sigma - 29)/3) and the step-1 chunk bounds of the profile.
These enter the E2-level modulation coupling  |a <D_rho f1, w~> + b <D_rho f2, w~>| <= Minv Lfn (||D f1||^2 + ||D f2||^2)^1/2 ||w~||^2."""
import sys, json, time
import sympy as sp
from flint import arb
sys.path.insert(0, '/home/claude/s1spec/step3nl'); sys.path.insert(0, '/home/claude/s1spec/step3')
import h2res as HR
k, m, p, q, c, s, X = sp.symbols('k m p q c s X')       # c = cos(theta - theta0), s = sin(theta - theta0)
rules = {p: -k * q, q: 1 + k * p, k: m, m: p / 4, c: -k * s, s: k * c}
def D(e): return sp.expand(sum(sp.diff(e, v) * r for v, r in rules.items()))
phi1 = X * (1 + X**2) ** sp.Rational(-1, 3) / 3
phi2 = sp.diff(phi1, X)
def Drho(g):
    d = [g]
    for _ in range(4): d.append(D(d[-1]))
    return sp.expand(d[4] - 2 * phi1 * d[3] + (phi1**2 - phi2) * d[2])
Dp = Drho(p)
Dc = Drho(c)                                           # = A c + B s
A_ = sp.expand(Dc).coeff(c); B_ = sp.expand(Dc).coeff(s)
assert sp.simplify(Dc - A_ * c - B_ * s) == 0
class _Pr(sp.printing.str.StrPrinter):
    def _print_Rational(self, e): return f"(arb({e.p})/arb({e.q}))"
    def _print_Integer(self, e): return f"arb({e.p})"
    def _print_Pow(self, e, rational=False):
        b, ex = e.as_base_exp(); bb = f"({self._print(b)})"
        if ex.is_Integer:
            n = int(ex); prod = "*".join([bb] * abs(n))
            return f"({prod})" if n > 0 else f"(1/({prod}))"
        return f"({bb}**(arb({ex.p})/arb({ex.q})))"
mk = lambda e: eval("lambda X, k, m, p, q: " + _Pr().doprint(e), {"arb": arb})
fDp, fA, fB = mk(Dp), mk(A_), mk(B_)
def rho_of(xlo): return HR.rho_of(xlo)
def core(GROUP=20):
    t0 = time.time(); n1 = arb(0); n2 = arb(0)
    lines = [l for l in open('/home/claude/s1_tight/prof/prof_tight.txt') if not l.startswith('END')]
    for g0 in range(0, len(lines), GROUP):
        grp = [list(map(float, l.split())) for l in lines[g0:g0 + GROUP]]
        ta, tb = min(v[0] for v in grp), max(v[1] for v in grp)
        pp, qq, kk, mm = [HR.I(min(v[13 + 2 * j] for v in grp), max(v[14 + 2 * j] for v in grp)) for j in range(4)]
        Xb = (HR.TT - arb(tb)).union(HR.TT - arb(ta))
        if Xb.upper() < 0: continue
        if Xb.lower() < 0: Xb = arb(0).union(arb(Xb.upper()))
        xlo, xhi = float(Xb.lower()), float(Xb.upper())
        v1 = abs(fDp(Xb, kk, mm, pp, qq)); a = fA(Xb, kk, mm, pp, qq); b = fB(Xb, kk, mm, pp, qq)
        v2sq = a * a + b * b
        w = arb(xhi - xlo) * rho_of(xlo)
        n1 += w * v1 * v1; n2 += w * v2sq
    # sigma in [8, 29] (chi = 1): chunk bounds
    N = 50
    for ci in range(8, 29):
        for i in range(N):
            sa = arb(ci) + arb(i) / N; sb = arb(ci) + arb(i + 1) / N
            kb, mb, pb, qs = HR.HS.profile_bounds(sa)
            Xb = (sa - HR.ST).union(sb - HR.ST)
            kk, mm, pp = arb(0, kb.upper()), arb(0, mb.upper()), arb(0, pb.upper()); qq = sa.union(sb) + arb(0, qs.upper())
            v1 = abs(fDp(Xb, kk, mm, pp, qq)); a = fA(Xb, kk, mm, pp, qq); b = fB(Xb, kk, mm, pp, qq)
            w = arb(float(Xb.upper()) - float(Xb.lower())) * rho_of(float(Xb.lower()))
            n1 += w * v1 * v1; n2 += w * (a * a + b * b)
    print(f"core: ||D f1||^2 <= {float(n1.upper()):.6g}, ||D f2||^2 <= {float(n2.upper()):.6g} [{time.time()-t0:.0f}s]")
    return n1, n2
def cutoff():
    """sigma in [29, 32]: (chi g)'''' etc.; |g^(j)| <= bounds (p: chunk bounds; cos/sin: <= 1 with derivatives
    bounded through k, m, p chunk bounds), |chi^(j)| <= S^(j)/3^j with |S'| <= 35/16*... (crude: max |S^(j)| on [0,1])."""
    Smax = [arb(1), arb(35) / 16 * 2, arb(10), arb(60), arb(400)]        # crude bounds for |S^(j)| on [0,1], j <= 4
    kb, mb, pb, qs = HR.HS.profile_bounds(arb(29))
    # |(chi g)^(j)| <= sum binom |chi^(i)| |g^(j-i)|; |g^(i)| <= G_i
    G1 = [pb, kb * (arb(32) + qs) + pb, arb(10) * (kb + mb + pb) * 33, arb(100) * (kb + mb + pb) * 33 * 33, arb(1000) * (kb + mb + pb) * 33 ** 3]
    G2 = [arb(1), kb, mb + kb * kb, arb(10) * (kb + mb + pb), arb(100) * (kb + mb + pb)]
    from math import comb
    def bnd(G):
        der = [sum((comb(j, i) * Smax[i] / arb(3) ** i * G[j - i] for i in range(j + 1)), arb(0)) for j in range(5)]
        x1 = arb(29) - HR.ST.upper(); ph1 = x1 * (1 + x1 * x1) ** (arb(-1) / 3) / 3 * 2   # generous |phi'| on the strip
        return der[4] + 2 * ph1 * der[3] + (ph1 * ph1 + 1) * der[2]
    x29 = float((arb(29) - HR.ST).lower()); x32 = float((arb(32) - HR.ST).upper())
    w = arb(x32 - x29) * rho_of(x29)
    c1, c2 = bnd(G1), bnd(G2)
    print(f"cutoff strip: |D(chi p)| <= {float(c1.upper()):.3g}, |D(chi cos)| <= {float(c2.upper()):.3g}, weight*len {float(w.upper()):.3g}")
    return w * c1 * c1, w * c2 * c2
if __name__ == "__main__":
    n1, n2 = core(); c1, c2 = cutoff()
    D1, D2 = (n1 + c1).sqrt(), (n2 + c2).sqrt()
    L2 = json.load(open('/home/claude/s1spec/step3/lyap_consts2.json'))
    coup = arb(L2["Minv"]) * arb(L2["Lf"]) * (D1 * D1 + D2 * D2).sqrt()
    out = {"Df1": float(D1.upper()), "Df2": float(D2.upper()), "E2_mod_coupling": float(coup.upper())}
    json.dump(out, open('dmod_norms.json', 'w'), indent=1)
    print(out)
