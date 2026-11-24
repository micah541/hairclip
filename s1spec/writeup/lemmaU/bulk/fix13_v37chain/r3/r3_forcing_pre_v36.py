#!/usr/bin/env python3
"""fix4_r3 (tag v36pre): the R3 correction of the forcing of Lemma S(A) (v3.3) and of Lemma U's c_dU.  Certified (python-flint arb, outward rounding;
every constant an exact rational).  Importable (functions below) and runnable (writes r3_forcing.json with an args block).

Defect (paper_v1/U2_TRACE.md, R3.2 and R3.2'): the scripts take lambda~ = -p/4, but p = 1.0246 exactly and lambda~ = -0.2561532224577463 (verifier:
LAMT = -0.2561532), so lambda~ + p/4 != 0; and the eps-part of the profile forcing was not bounded at all.

Derivation (exact).  (A)'s variable is u = f - U, U = c_U alpha(tau) y^p, c_U = 1.3056, s = ln y, eps = y^-4.  The far-field slope equation
(lemmaS.tex S1, lemmaS_symbolic.py: V_tau = -(1/4 - a) y V' + b cos(theta) V' - d_y^2 (A(V) V'' + B(V) V'^2), A(x) = (1+x^2)^-2, B(x) = -3x(1+x^2)^-3)
reads in s:  V_tau + c V_s = -eps (d^2 - 5d + 6) G(V),  G(V) = A(V)(V_ss - V_s) + B(V) V_s^2,  c = 1/4 - a - b cos(theta) e^-s
(d_y^2 (y^-2 G) = y^-4 (d^2 - 5d + 6) G).  With V = U + u, G(V) = G(U) + A_2 u_ss + A_1 u_s + A_0 u exactly, so the forcing is
  R = -U_tau - c U_s - eps (d^2 - 5d + 6) G(U).
(i) Non-eps part.  U_s = pU, U_tau = c_U alpha_tau y^p, alpha_tau = lambda~ alpha + beta, |beta| <= betabar (D1):
  -U_tau - c U_s = -c_U y^p [beta + (lambda~ + p/4) alpha - p a alpha] + c_U p b cos(theta) alpha y^(p-1).
  d_s^j y^p = p^j y^p and |d_s^j y^(p-1)| = |p-1|^j y^(p-1) <= p^j y^p / 14 on y >= 14.  Hence |d_s^j (non-eps part)| <= r0' p^j y^p with
  r0' = c_U (betabar + DL alpha_1) + c_U p alpha_1 (abar + bbar/14),   DL = max(|lambda~ + p/4|, |LAMT + p/4|)  (covers either reading of beta).
(ii) Eps part.  On U = c_U alpha y^p, d_s acts on any function of U as p*D with D = U d/dU (Euler operator), and d_s^j (eps h(U)) = eps (pD - 4)^j h(U).
  So d_s^j R_eps = -y^-4 Phi_j(U),  Phi_j = (pD - 4)^j (p^2 D^2 - 5pD + 6) g,  g(U) = G(U) = (p^2 - p) U (1+U^2)^-2 - 3 p^2 U^3 (1+U^2)^-3.
  Phi_j is odd; h_j = Phi_j / U is an exact rational function of w = U^2 (num and den with rational coefficients, computed by sympy, evaluated in arb).
  Power series: g = sum_k g_k U^(2k+1) gives Phi_j = sum_k g_k (p(2k+1) - 4)^j (p^2(2k+1)^2 - 5p(2k+1) + 6) U^(2k+1); the linear coefficient is
  c_eps,j = c_U |p-4|^j (p^2-5p+6)(p^2-p) (U2_TRACE R3.2), the cubic one is exported for the record.  The bound used is NOT the series: on each
  y-cell [ya, yb], |U| <= Um = c_U alpha_1 yb^p, so |d_s^j R_eps| <= c_U alpha_1 y^(p-4) k_cell, k_cell = sup_{0 <= w <= Um^2} |h_j(w)| (arb enclosure).
(iii) L^2 forcing.  F_j^2 = int omega |d_s^j R|^2 ds, omega = (y_ref/y)^gamma (gamma = 1), ds = dy/y.  Minkowski:
  F_j <= r0' p^j sqrt(I_int) + sqrt(Feps_j^2),  I_int = y_ref int y^(2p-2) dy,  Feps_j^2 = sum_cells (c_U alpha_1 k_cell)^2 y_ref int_cell y^(2p-10) dy.
(iv) Lemma U, delta-U lemma: d/dtau[alpha lambda^-p] = (beta + (lambda~ + p/4) alpha - p a alpha) lambda^-p, so
  c_dU = c_U (betabar + DL alpha_1 + p abar alpha_1) / (p (1/4 - abar))   (betabar' = betabar + DL alpha_1).
"""
import os, sys, json, hashlib, platform, math
sys.dont_write_bytecode = True
from fractions import Fraction
import sympy as sp
from flint import arb, ctx
ctx.prec = 200

def Aq(q):
    q = sp.Rational(q); return arb(int(q.p))/arb(int(q.q))
CU_Q = sp.Rational('1.3056'); P_Q = sp.Rational('1.0246')
LT_Q = sp.Rational('-0.2561532224577463')            # lambda~ (step3nl/h2md.py:13, psi3beta.py:5; LC:59)
LAMT_Q = sp.Rational('-0.2561532')                   # verifier constants.py:14
DL_Q = max(abs(LT_Q + P_Q/4), abs(LAMT_Q + P_Q/4))   # exact rational
CU, P, DL = Aq(CU_Q), Aq(P_Q), Aq(DL_Q)

# ---------------------------------------------------------------- exact symbolic Phi_j / U as num(w)/den(w), w = U^2
_U, _w = sp.symbols('U w', positive=True)
_g = (P_Q**2 - P_Q)*_U/(1 + _U**2)**2 - 3*P_Q**2*_U**3/(1 + _U**2)**3
def _Dop(e): return sp.together(_U*sp.diff(e, _U))
def _phi(j):
    e = sp.together(P_Q**2*_Dop(_Dop(_g)) - 5*P_Q*_Dop(_g) + 6*_g)
    for _ in range(j): e = sp.together(P_Q*_Dop(e) - 4*e)
    return sp.cancel(e)
_H = {}
def hpoly(j):
    """exact (num, den) coefficient lists in w = U^2 of h_j = Phi_j / U (lowest degree first)"""
    if j not in _H:
        n, d = sp.fraction(sp.cancel(_phi(j)/_U))
        n = sp.Poly(sp.expand(n), _U); d = sp.Poly(sp.expand(d), _U)
        assert all(m[0] % 2 == 0 for m in n.monoms()) and all(m[0] % 2 == 0 for m in d.monoms()), 'h_j must be even in U'
        cn = [sp.Rational(0)]*(n.degree()//2 + 1); cd = [sp.Rational(0)]*(d.degree()//2 + 1)
        for (k,), c in zip(n.monoms(), n.coeffs()): cn[k//2] = sp.Rational(c)
        for (k,), c in zip(d.monoms(), d.coeffs()): cd[k//2] = sp.Rational(c)
        assert cd[0] != 0
        _H[j] = (cn, cd)
    return _H[j]
def _horner(cs, x):
    r = arb(0)
    for c in reversed(cs): r = r*x + Aq(c)
    return r
def hsup(j, W):
    """certified upper bound of sup_{0 <= w <= W} |h_j(w)| (W an arb >= 0; ball [0, W.upper()])"""
    cn, cd = hpoly(j); wb = arb(0).union(arb(W.upper()))
    den = _horner(cd, wb); assert den.lower() > 0, 'denominator not bounded away from 0'
    return (abs(_horner(cn, wb))/den).upper()
def lin_coeff(j):
    """exact linear coefficient c_eps,j / c_U = (p-4)^j (p^2-5p+6)(p^2-p)  (= h_j(0)); checked against the symbolic h_j"""
    v = (P_Q - 4)**j*(P_Q**2 - 5*P_Q + 6)*(P_Q**2 - P_Q); cn, cd = hpoly(j); assert sp.Rational(cn[0], cd[0]) == v; return v
def cubic_coeff(j):
    """exact U^3 coefficient of Phi_j (record only): g_1 (3p-4)^j (9p^2-15p+6), g_1 = -2(p^2-p) - 3p^2"""
    g1 = -2*(P_Q**2 - P_Q) - 3*P_Q**2; v = g1*(3*P_Q - 4)**j*(9*P_Q**2 - 15*P_Q + 6)
    ser = sp.series(_phi(j), _U, 0, 4).removeO(); assert sp.simplify(ser.coeff(_U, 3) - v) == 0; return v

# ---------------------------------------------------------------- the corrected constants
def r0_honest(betabar, alpha1, abar, bbar):
    """r0' = c_U (betabar + DL alpha_1) + c_U p alpha_1 (abar + bbar/14)   (arb inputs = upper bounds)"""
    return CU*(betabar + DL*alpha1) + CU*P*alpha1*(abar + bbar/14)
def r0_old(betabar, alpha1, abar):
    return Aq('1.4')*(betabar + Aq('1.02')*alpha1*abar)
def feps2(j, ya, yb, yref, alpha1, n=400):
    """certified upper bound of int_{ya}^{yb} (yref/y) |d_s^j R_eps|^2 dy/y  (gamma = 1); n geometric y-cells"""
    ya, yb, yref = arb(ya), arb(yb), arb(yref); e = 2*P - 9; tot = arb(0); r = (yb/ya).log()
    for i in range(n):
        a_ = ya*(r*i/n).exp() if i else ya; b_ = ya*(r*(i + 1)/n).exp() if i + 1 < n else yb
        Um = CU*alpha1*(arb(b_.upper()).log()*P).exp()
        k = hsup(j, Um*Um)
        cell = (arb(a_.lower())**e - arb(b_.upper())**e)/(-e)          # int_cell y^(2p-10) dy, enlarged cell [a_.lower, b_.upper]
        tot += ((CU*alpha1*k)**2*yref*cell).upper()
    return arb(tot.upper())
def F_honest(j, r0n, Iint, fe2):
    """F_j <= r0' p^j sqrt(Iint) + sqrt(Feps_j^2)"""
    return r0n*P**j*Iint.sqrt() + fe2.sqrt()
def cdU_honest(betabar, abar, alpha1):
    return CU*(betabar + DL*alpha1 + P*abar*alpha1)/(P*(arb(1)/4 - abar))
def cdU_old(betabar, abar, alpha1):
    return CU*(betabar + P*abar*alpha1)/(P*(arb(1)/4 - abar))

# ---------------------------------------------------------------- standalone record
if __name__ == '__main__':
    import time
    HERE = os.path.dirname(os.path.abspath(__file__))
    BASE = '/home/claude/s1spec/writeup/lemmaU/bulk/fix4_r3/t2_A/base_ibp'
    FINP = BASE + '/lemmaS_F10p_M0.52.json'; VA = '/home/claude/s1spec/step3nl_v33/verifier/modsups_fin12f4/final_alpha_points_5e-8.json'
    a_ = json.load(open(FINP))['args']
    up = lambda x: Aq(repr(math.nextafter(float(x), math.inf)))       # same input rule as assembly_arb.UP
    AL, ABAR, BETA = up(a_['alpha1']), up(a_['abar']), up(a_['beta'])
    vs = json.load(open(VA))['1.5e-4']; ABV, BBV = Aq(str(vs['sup_ab'])), Aq(str(vs['sup_beta']))   # Lemma U reads these (t4lib.verifier_sups)
    G = arb(1); e = 2*P - G; Iint = arb(28)*((arb(700)**e - arb(14)**e)/e)
    r0o, r0n = r0_old(BETA, AL, ABAR), r0_honest(BETA, AL, ABAR, ABAR)
    res = dict(DL=str(DL_Q), DL_double=float(DL_Q), lambda_tilde=str(LT_Q), LAMT=str(LAMT_Q), p=str(P_Q), c_U=str(CU_Q),
               inputs_A=dict(alpha1=a_['alpha1'], abar=a_['abar'], beta=a_['beta'], rule='one ulp up (assembly_arb.UP)'),
               r0_old=r0o.str(20, radius=True), r0_new=r0n.str(20, radius=True), r0_ratio=(r0n/r0o).str(10, radius=True),
               r0_new_parts=dict(beta=(CU*BETA).str(12), lam=(CU*DL*AL).str(12), a=(CU*P*AL*ABAR).str(12), b=(CU*P*AL*ABAR/14).str(12)),
               Iint_14_700=Iint.str(20, radius=True), levels={})
    for j in range(4):
        fe = feps2(j, 14, 700, 28, AL); Fo = r0o*P**j*Iint.sqrt(); Fn = F_honest(j, r0n, Iint, fe)
        res['levels'][str(j)] = dict(c_eps_lin_over_cU=str(lin_coeff(j)), c_eps_lin=float(CU_Q*abs(lin_coeff(j))), cubic_coeff_over_cU=str(cubic_coeff(j)),
                                     hsup_at_700=float(hsup(j, (CU*AL*(arb(700).log()*P).exp())**2)), Feps2=fe.str(15, radius=True),
                                     F_old=Fo.str(15, radius=True), F_new=Fn.str(15, radius=True), F2_old=(Fo*Fo).str(15, radius=True), F2_new=(Fn*Fn).str(15, radius=True),
                                     ratio_F=float((Fn/Fo).upper()), label='cert')
    co, cn = cdU_old(BBV, ABV, Aq('1.5e-4')), cdU_honest(BBV, ABV, Aq('1.5e-4'))
    res['c_dU'] = dict(old=co.str(20, radius=True), new=cn.str(20, radius=True), ratio=float((cn/co).upper()), betabar_prime=(BBV + DL*Aq('1.5e-4')).str(15),
                       inputs=dict(sup_ab=str(vs['sup_ab']), sup_beta=str(vs['sup_beta'])), assert_lt_2e8=bool(cn.upper() < Aq('2e-8')), label='cert')
    assert res['c_dU']['assert_lt_2e8']
    ins = [os.path.abspath(__file__), FINP, VA]
    res['args'] = dict(argv=sys.argv, cwd=os.getcwd(), python=sys.version.split()[0], platform=platform.platform(), flint=__import__('flint').__version__,
                       sympy=sp.__version__, prec=ctx.prec, env={k: os.environ.get(k) for k in ('PYTHONHASHSEED', 'PYTHONDONTWRITEBYTECODE', 'OMP_NUM_THREADS')},
                       inputs={q: hashlib.sha256(open(q, 'rb').read()).hexdigest() for q in ins})
    out = os.path.join(HERE, 'r3_forcing.json')
    json.dump(res, open(out + '.tmp', 'w'), indent=1, sort_keys=True); os.replace(out + '.tmp', out)
    print(json.dumps({k: v for k, v in res.items() if k != 'args'}, indent=1))
