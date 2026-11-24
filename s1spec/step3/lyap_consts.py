"""Rigorous constants for the mixed triangular Lyapunov lemma (fourth supplement, Sec. 4).

Inputs: outA4.txt  (hsrun3 run of gen_hs3.build_A with MONALL=0.01: hulls of the state on subintervals,
        final doubleton).  Rows, stage B: st C S ws0..ws3 I*;  stage C: p q k m st C S ws0..ws3 I*.
On each subinterval everything is evaluated in ball arithmetic from the hull; integrals are bounded by
(length in X) * sup(integrand).  F^(n) = sum_j binom(n,j) ws_j B_{n-j}  (ws_j = w^(j)/rho, B_i = (e^phi)^(i)/e^phi).
Output: eps1, eps2, ||g||, ||phi3||, and the largest omega' with
     eps1 (||g|| + eps1 ||phi3||^2) <= (|lamt| - eps1 ||phi3|| - omega') (omega - eps2 - eps1 ||phi3|| - omega')."""
import sys, math, json
import sympy as sp
from math import comb
from flint import arb, ctx
sys.path.insert(0, '/home/claude/s1spec/step3')
import gen_hs3 as G3
from gen_hs2 import chunks
import hstail as H
from gram_post7 import Aff, DSet
ctx.prec = 100
U = arb(0, 1)
def A(lo, hi): return arb((lo + hi) / 2) + U * ((hi - lo) / 2 * (1 + 1e-15) + 1e-300)
def hull(a, b):                     # convex hull of two arbs
    lo = min(a.lower(), b.lower()); hi = max(a.upper(), b.upper()); return (lo + hi) / 2 + U * ((hi - lo) / 2)
def sq_sup(x): return (abs(x) * abs(x)).upper()          # sup |x|^2 (no ** : NaN for zero-midpoint balls)
def mag(x): return abs(x).upper()
R_ = lambda p, q: arb(p) / arb(q)

# ------------------------------------------------------------------ symbolic operators
X, P1, P2, P3, P4, k, m, q, p, lam = sp.symbols('X P1 P2 P3 P4 k m q p lam')
d = sp.symbols('d0:5'); ch = sp.symbols('c0:5')
Bp = [sp.Integer(1), P1, P1**2 + P2, P1**3 + 3*P1*P2 + P3, P1**4 + 6*P1**2*P2 + 3*P2**2 + 4*P1*P3 + P4]
Em = [sp.Integer(1), -P1, P1**2 - P2, -P1**3 + 3*P1*P2 - P3, P1**4 - 6*P1**2*P2 + 3*P2**2 + 4*P1*P3 - P4]
aL = {4: sp.Integer(-1), 3: sp.Integer(0), 2: -k**2, 1: -(5*k*m + q/4), 0: sp.Rational(1, 4) - 3*m**2 - k*p}
def der(n): return sum(comb(n, j) * Em[n - j] * d[j] for j in range(n + 1))   # e^phi (e^{-phi} psi)^(n)
Ls = sp.expand(-der(4) - k**2 * der(2) + (k*m + q/4) * der(1) + d[0] / 2)
aS = {n: sp.expand(Ls.coeff(d[n])) for n in range(5)}
assert sp.simplify(aS[4] + 1) == 0
class _Pr(sp.printing.str.StrPrinter):
    def _print_Rational(self, e): return f"R_({e.p},{e.q})"
    def _print_Integer(self, e): return f"arb({e.p})"
    def _print_Pow(self, e, rational=False):
        b, ex = e.as_base_exp()
        if ex.is_Integer:                       # integer powers as products (python-flint 0.9: ball**n is NaN if mid = 0)
            n = int(ex); bb = f"({self._print(b)})"
            prod = "*".join([bb] * abs(n)) if n != 0 else "arb(1)"
            return f"({prod})" if n > 0 else f"(1/({prod}))"
        return super()._print_Pow(e, rational)
def fn(expr, args):
    src = _Pr().doprint(sp.expand(expr))
    return eval("lambda " + ",".join(str(a) for a in args) + ": " + src, {"arb": arb, "R_": R_})
VARS = (P1, P2, P3, P4, k, m, q, p, lam)
cL = {n: fn(aL[n], VARS) for n in range(5)}
cS = {n: fn(aS[n], VARS) for n in range(5)}
# psi1 = p and the pieces C, S of psi2: sigma-derivatives via the profile ODE (X-derivatives = sigma-derivatives)
Cs, Ss = sp.symbols('Cs Ss')
rules = {p: -k*q, q: 1 + k*p, k: m, m: p/4, Cs: k*Ss, Ss: -k*Cs}
def D(e): return sp.expand(sum(sp.diff(e, v) * rv for v, rv in rules.items()))
def ders(e):
    out = [e]
    for _ in range(4): out.append(D(out[-1]))
    return out
PV = (k, m, q, p, Cs, Ss)
dp = [fn(e, PV) for e in ders(p)]; dC = [fn(e, PV) for e in ders(Cs)]; dS = [fn(e, PV) for e in ders(Ss)]

def weight(Xv):
    E = (1 + Xv * Xv) ** (arb(-1) / 3)
    X2 = Xv * Xv; E2 = E * E; E4 = E2 * E2; E7 = E4 * E2 * E; E10 = E7 * E2 * E
    P = [Xv * E / 3, E / 3 - 2 * X2 * E4 / 9, -2 * Xv * E4 / 3 + 16 * X2 * Xv * E7 / 27,
         -2 * E4 / 3 + 32 * X2 * E7 / 9 - 224 * X2 * X2 * E10 / 81]
    rho = (-(1 + Xv * Xv) ** (arb(2) / 3) / 4).exp()
    return P, rho
def bell(P, n):
    P1_, P2_, P3_, P4_ = P
    return [arb(1), P1_, P1_*P1_ + P2_, P1_*P1_*P1_ + 3*P1_*P2_ + P3_][n]
# smoothstep S(u) and derivatives
Su = sp.symbols('u'); Sm = Su**4 * (35 - 84*Su + 70*Su**2 - 20*Su**3)
dSm_pts = [fn(sp.diff(Sm, Su, j), (Su,)) for j in range(8)]          # exact at points; S has degree 7
def centered(ders_pts, I, j):
    """f^(j) on the ball I for a polynomial f with derivative list ders_pts (exact Taylor at the midpoint)."""
    m_ = arb(I.mid()); dl = I - m_; tot_ = arb(0); pw = arb(1); fac = 1
    for i in range(0, len(ders_pts) - j):
        if i > 0: pw = pw * dl; fac *= i
        tot_ += ders_pts[j + i](m_) * pw / fac
    return tot_
class _DS:
    def __getitem__(self, j): return lambda I: centered(dSm_pts, I, j)
dSm = _DS()
# bumps eta1 = LB s zeta(s), eta3 = LB^3 s^3 zeta/6, zeta = (1-s^2)^4, s = X/(c LB); derivatives in s
sv = sp.symbols('s'); zeta = (1 - sv**2) ** 4; LBs = sp.Rational(G3.LB).limit_denominator(10)
e1s = LBs * sv * zeta; e3s = LBs**3 * sv**3 * zeta / 6
de1_pts = [fn(sp.diff(e1s, sv, j), (sv,)) for j in range(10)]; de3_pts = [fn(sp.diff(e3s, sv, j), (sv,)) for j in range(12)]
de1 = [(lambda j: (lambda I: centered(de1_pts, I, j)))(j) for j in range(5)]
de3 = [(lambda j: (lambda I: centered(de3_pts, I, j)))(j) for j in range(5)]
def _pw(x, j):
    r = arb(1)
    for _ in range(j): r = r * x
    return r
def bump_ders(Xv, cc, which):
    s_ = Xv / (cc * arb(G3.LB)); f = de1 if which == 1 else de3
    return [f[j](s_) / _pw(cc * arb(G3.LB), j) for j in range(5)]

# ------------------------------------------------------------------ data
LAMT = arb(G3.LAMT); T0 = arb(G3.T0); OMEGA = arb(3) / 64
out = open('outA4.txt').read()
CH = chunks()
NA = len(G3.A_INT)
tot = {key: arb(0) for key in ["Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7", "Q8", "Q9", "P11", "P13", "P21", "P23",
                               "Lf1", "Lf2", "e1z", "e2z"]}
fin = DSet('outA4.txt', [0, 1, 2, 3, 4])
Ct, St = fin.hull[5], fin.hull[6]
nsub = 0
for line in out.splitlines():
    tag = line[:2]
    if tag not in ("MB", "MC"): continue
    v = list(map(float, line.split()[1:])); ta, tb = v[0], v[1]; h = [A(v[2 + 2 * i], v[3 + 2 * i]) for i in range((len(v) - 2) // 2)]
    nsub += 1
    if tag == "MB":
        st, C_, S_, ws = h[0], h[1], h[2], h[3:7]
        c = min(int(ta + 1e-9), len(CH) - 1); pv = CH[c][2]
        kk, mm, pp, dq = [A(lo, hi) for lo, hi in pv]
        sig = A(32 - tb, 32 - ta); Xv = sig - st
        kv, mv, pv_, qv = kk, mm, pp, sig + dq
        dX = arb(tb - ta)
        if tb <= 3.0 + 1e-12:
            uu = A(ta / 3, tb / 3)
            chi = [dSm[j](uu) * _pw(arb(-1) / 3, j) for j in range(5)]
        else:
            chi = [arb(1)] + [arb(0)] * 4
    else:
        pv_, qv, kv, mv, st, C_, S_, ws = h[0], h[1], h[2], h[3], h[4], h[5], h[6], h[7:11]
        cc = (8 - st) / T0
        Xv = hull(cc * (T0 - arb(tb)), cc * (T0 - arb(ta)))
        dX = cc.upper() * arb(tb - ta)
        chi = [arb(1)] + [arb(0)] * 4
    P, rho = weight(Xv)
    args = (P[0], P[1], P[2], P[3], kv, mv, qv, pv_, LAMT)
    aLv = {n: cL[n](*args) for n in range(5)}; aSv = {n: cS[n](*args) for n in range(5)}
    Fd = [sum((comb(n, j) * ws[j] * bell(P, n - j) for j in range(n + 1)), arb(0)) for n in range(4)]
    LmLs = sum(((aLv[n] - aSv[n]) * Fd[n] for n in range(4)), arb(0))           # (L - L*) F = (L - lamt) F
    w = dX * rho
    tot["Q1"] += w.upper() * sq_sup(chi[0] * LmLs) if True else 0
    if chi[1] != 0 or tag == "MB" and tb <= 3.0 + 1e-12:
        cs = sum((aSv[n] * sum((comb(n, j) * chi[j] * Fd[n - j] for j in range(1, n + 1)), arb(0)) for n in range(1, 5)), arb(0))
        cl = sum((aLv[n] * sum((comb(n, j) * chi[j] * Fd[n - j] for j in range(1, n + 1)), arb(0)) for n in range(1, 5)), arb(0))
        tot["Q2"] += (w * (abs(cs) * abs(cs))).upper(); tot["Q3"] += (w * (abs(cl) * abs(cl))).upper()
    # psi1 = p, psi2 = Ct C + St S and L*(chi psi_j)
    pvs = (kv, mv, qv, pv_, C_, S_)
    ps1 = [dp[n](*pvs) for n in range(5)]
    ps2 = [Ct * dC[n](*pvs) + St * dS[n](*pvs) for n in range(5)]
    for key, ps in (("Lf1", ps1), ("Lf2", ps2)):
        Lps = sum((aSv[n] * ps[n] for n in range(5)), arb(0))
        com = sum((aSv[n] * sum((comb(n, j) * chi[j] * ps[n - j] for j in range(1, n + 1)), arb(0)) for n in range(1, 5)), arb(0))
        tot[key] += (w * (abs(chi[0] * Lps + com) * abs(chi[0] * Lps + com))).upper()
    if tag == "MB" and tb <= 3.0 + 1e-12:                                          # (1 - chi) psi_i on sigma in [29, 32]
        tot["e1z"] += (w * (abs((1 - chi[0]) * ps1[0]) * abs((1 - chi[0]) * ps1[0]))).upper(); tot["e2z"] += (w * (abs((1 - chi[0]) * ps2[0]) * abs((1 - chi[0]) * ps2[0]))).upper()
    # tip bumps: support X <= LB c
    if tag == "MC" and ta >= float(T0) - G3.LB - 1e-12:
        for which, (kq, kl, kn) in ((1, ("Q4", "Q6", "Q8")), (3, ("Q5", "Q7", "Q9"))):
            ed = bump_ders(Xv, cc, which)
            es = sum((aSv[n] * ed[n] for n in range(5)), arb(0)) - LAMT * ed[0]
            el = sum((aLv[n] * ed[n] for n in range(5)), arb(0)) - LAMT * ed[0]
            tot[kq] += (w * (abs(es) * abs(es))).upper(); tot[kl] += (w * (abs(el) * abs(el))).upper(); tot[kn] += (w * (abs(ed[0]) * abs(ed[0]))).upper()
            tot[f"P1{which}"] += (w * abs(ps1[0]) * abs(ed[0])).upper(); tot[f"P2{which}"] += (w * abs(ps2[0]) * abs(ed[0])).upper()
print("subintervals:", nsub)
sq = lambda x: x.sqrt() if x > 0 else arb(0)
nr = {kk_: sq(vv) for kk_, vv in tot.items()}
# far tail sigma > 32 for (1 - chi) psi_i: |psi2| <= 1, |p| <= pb(32);  int_{x0}^oo rho <= exp(-x0^{4/3}/4) / (x0^{1/3}/3)
x0 = arb(32) - H.STIP.upper()
tail = (-(x0 ** (arb(4) / 3)) / 4).exp() / (x0 ** (arb(1) / 3) / 3)
kb32, mb32, pb32, qs32 = H.profile_bounds(32)
e1n = (tot["e1z"] + pb32 * pb32 * tail).sqrt(); e2n = (tot["e2z"] + tail).sqrt()
# ------------------------------------------------------------------ final values from the doubleton (shared affine forms)
nst = 5; iws1, iws3 = nst + 3, nst + 5
st_f = fin.form({4: arb(1)}); ws1 = fin.form({iws1: arb(1)}); ws3 = fin.form({iws3: arb(1)})
cc_f = (Aff(arb(8)) - st_f).scale(1 / T0)
a_ = (cc_f * ws1).interval(); b_ = (cc_f * cc_f * cc_f * (ws3 + ws1) + (cc_f * ws1).scale(arb(24) / arb(G3.LB) ** 2)).interval()
iI = {nm: 11 + j for j, nm in enumerate(G3.A_INT)}
Ctf, Stf = fin.form({5: arb(1)}), fin.form({6: arb(1)})
IFF = fin.form({iI["IFF"]: arb(1)}).interval()
psiF1 = fin.form({iI["IpF"]: arb(1)}).interval()
psiF2 = (Ctf * fin.form({iI["ICF"]: arb(1)}) + Stf * fin.form({iI["ISF"]: arb(1)})).interval()
M11 = fin.form({iI["Ipp"]: arb(1)}).interval()
M12 = (Ctf * fin.form({iI["IpC"]: arb(1)}) + Stf * fin.form({iI["IpS"]: arb(1)})).interval()
M22 = (Ctf * Ctf * fin.form({iI["ICC"]: arb(1)}) + (Ctf * Stf).scale(arb(2)) * fin.form({iI["ICS"]: arb(1)})
       + Stf * Stf * fin.form({iI["ISS"]: arb(1)})).interval()
tr, det = M11 + M22, M11 * M22 - M12 * M12
lmin = (tr - (tr * tr - 4 * det).sqrt()) / 2; lmax = (tr + (tr * tr - 4 * det).sqrt()) / 2
Minv = 1 / lmin.lower()
# ------------------------------------------------------------------ assemble
normF = IFF.sqrt()
nf3_hi = normF.upper() + mag(a_) * nr["Q8"] + mag(b_) * nr["Q9"]
nf3_lo = normF.lower() - mag(a_) * nr["Q8"] - mag(b_) * nr["Q9"]
r3 = nr["Q2"] + mag(a_) * nr["Q4"] + mag(b_) * nr["Q5"]
Lf3 = nr["Q1"] + nr["Q3"] + mag(a_) * nr["Q6"] + mag(b_) * nr["Q7"]
psif3 = [mag(psiF1) + mag(a_) * tot["P11"] + mag(b_) * tot["P13"], mag(psiF2) + mag(a_) * tot["P21"] + mag(b_) * tot["P23"]]
psif3n = (psif3[0] * psif3[0] + psif3[1] * psif3[1]).sqrt()
Lfn = (tot["Lf1"] + tot["Lf2"]).sqrt()
En = (e1n * e1n + e2n * e2n).sqrt()
eps1 = r3 + psif3n * Minv * Lfn
eps2 = En * Minv * Lfn
# ||Pi|| <= 1 + delta_Pi  (Pi = F*M^-1 F + E M^-1 F,  M = G_f + Delta,  ||Delta|| <= ||E|| sqrt(tr G_f))
trGf = tr + 2 * En * (M11.sqrt() + M22.sqrt()) + En * En
Dn = En * trGf.sqrt()
lmin_f = lmin.lower() - Dn
dPi = trGf * Minv * Dn / lmin_f + En * Minv * trGf.sqrt()
Pin = 1 + dPi
Nlo = nf3_lo * nf3_lo - Minv * psif3n * (psif3n + nf3_hi * En)
g_ = Pin * Lf3 / Nlo
phi3 = Pin * nf3_hi / Nlo
lamabs = abs(LAMT).lower()
def ok(wp):
    lhs = eps1 * (g_ + eps1 * phi3 * phi3)
    A0 = lamabs - eps1 * phi3 - wp; Bt = OMEGA.lower() - eps2 - eps1 * phi3 - wp
    return bool(A0 > 0 and Bt > 0 and lhs < A0 * Bt), lhs, A0 * Bt
lo, hi = 0.0, float(OMEGA)
for _ in range(40):
    mid = (lo + hi) / 2
    if ok(arb(mid))[0]: lo = mid
    else: hi = mid
res = dict(a=float(mag(a_)), b=float(mag(b_)), normF=float(normF.mid()), r3=float(r3), Lf3=float(Lf3),
           r3_far=float(nr["Q2"]), r3_bump=float(mag(a_) * nr["Q4"] + mag(b_) * nr["Q5"]),
           psif3=[float(x) for x in psif3], Minv=float(Minv), Lf=float(Lfn), E=float(En), eps1=float(eps1), eps2=float(eps2),
           Pi=float(Pin), N=float(Nlo), g=float(g_), phi3=float(phi3), rel_eps1=float(eps1 / normF.lower()),
           g_times_f3=float(g_ * normF.upper()), omega_prime=lo)
for kk_, vv in res.items(): print(f"  {kk_:12s} {vv}")
print("condition holds with omega' =", lo, ":", ok(arb(lo))[0])
json.dump(res, open("lyap_consts.json", "w"), indent=1)
