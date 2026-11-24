#!/usr/bin/env python3
"""Stage S1 (fin13v4): exact sympy checks of the identities that EI_SPEC.md marks [hand] (U-2), before any constant is computed.
 I12  chain rule for the coefficients:  D d = d^a + d^b g'  with d^b = d_U1 (exactly linear in g' = w''), for every d in {c, d010, d100, d200, d110, d300}
 I13  N decomposition: F(u0+w) - F(u0) - L w = c w'''' + d010 w''' + d100 w'' + d200 w''^2 + d110 w'' w''' + d300 w''^3 + d000, c = 1 - J^-4 (exact)
 I14  D(d000) = E1c g' + E0c(w, g),  E1c = d_U1 G(U) + (5km + q/4) (= the A1 expression),  E0c(0, 0) = 0
 I15  N_low' table (EI_SPEC 4.5) = D(N_low) exactly (generic coefficient functions)
 I16  v-equation: zeta * d/dx(w_tau) - (T1 + ... + T12) = 0 (EI_SPEC 2.3), profile ODE substituted
 FD   float cross-check (mpmath, 40 digits) of D(N_low) against a numerical x-derivative of N_low built from F itself, at 6 random points,
      with the profile ODE integrated numerically and random smooth u0, w.
Writes out/s1_identities.json.  Run: PYTHONHASHSEED=0 nice -n 19 python3 s1_identities.py"""
import sys, os, json, time, pickle, random
import sympy as sp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from v4common import args_block, dump, OUT
t0 = time.time()
def Z(e):
    """exact zero test: expand, then together/cancel (rational functions in sqrt(J^2)); returns '0' or the residual"""
    t_ = time.time(); r = sp.expand(e)
    print('  Z expand', round(time.time() - t_, 1), 'zero' if r == 0 else 'nonzero', flush=True)
    if r == 0: return '0'
    r = sp.cancel(sp.together(r))
    return '0' if r == 0 else str(r)[:200]
GSP = '/home/claude/s1spec/step3nl_v4b/graph_struct.pkl'
GS = pickle.load(open(GSP, 'rb'))
Ji, al, k, m, p, q, u0, u1 = sp.symbols('Ji al k m p q u0 u1')
U0, U1, a1, a2, a3, a4, a5, w2, w3, w4 = sp.symbols('U0 U1 a1 a2 a3 a4 a5 w2 w3 w4')
us, w, g, g1, ct, st = sp.symbols('us w g g1 ct st')
J = sp.sqrt((1 - k * U0) ** 2 + U1 ** 2)
sub = {Ji: 1 / J, al: 1 - k * U0, u0: U0, u1: U1}
Fpert = sum(cf.subs(sub) * (a2 + w2) ** i * (a3 + w3) ** j * (a4 + w4) ** l for (i, j, l), cf in GS.items())
Fpert += (p + U0 - U1 * q / (1 - k * U0)) / 4                     # as in agent_int/planB/pb_nlin.py
Lw = -w4 - k ** 2 * w2
Pw = sp.Poly(sp.expand(Fpert - Lw - Fpert.subs({w2: 0, w3: 0, w4: 0})), w2, w3, w4)
coef = {mon: c for mon, c in zip(Pw.monoms(), Pw.coeffs())}
res = {}
print('built', round(time.time() - t0, 1), flush=True)
res['monomials'] = str(sorted(coef))
assert set(coef) == {(0, 0, 1), (0, 1, 0), (1, 0, 0), (1, 1, 0), (2, 0, 0), (3, 0, 0)}, sorted(coef)
# I13: c = 1 - J^-4
res['I13_c_is_1_minus_Jm4'] = Z(coef[(0, 0, 1)] - (1 - 1 / J ** 4))
res['I13_d110'] = Z(coef[(1, 1, 0)] - 10 * U1 / J ** 6)
res['I13_d300'] = Z(coef[(3, 0, 0)] - 3 / J ** 6 * (1 - 6 * U1 ** 2 / J ** 2)); print("res['I13_d300']", res['I13_d300'][:60], round(time.time() - t0), flush=True)
# total derivative operator (rules: U0' = U1, U1' = a2 + g1, a_j' = a_{j+1}, profile ODE, theta' = k)
RULES = {U0: U1, U1: a2 + g1, a1: a2, a2: a3, a3: a4, a4: a5, k: m, m: p / 4, p: -k * q, q: 1 + k * p, us: a1, w: g, g: g1, ct: -k * st, st: k * ct}
def Dx(e): return sum(sp.diff(e, v) * r for v, r in RULES.items() if e.has(v))
NAMES = {(0, 0, 1): 'c', (0, 1, 0): 'd010', (1, 0, 0): 'd100', (2, 0, 0): 'd200', (1, 1, 0): 'd110', (3, 0, 0): 'd300'}
# I12
i12 = {}
for mon, nm in NAMES.items():
    d = coef[mon]; Dd = sp.expand(Dx(d))
    da = Dd.subs(g1, 0); db = sp.diff(d, U1)
    i12[nm] = Z(Dd - da - db * g1); print('I12', nm, i12[nm][:40], round(time.time() - t0), flush=True)
res['I12_chain_rule_linear_in_g1'] = i12
# I14 (structural, generic G = G(U0, U1, a2, a3, a4, k, m, p, q) as an undefined function, so the check is small):
#   D(d000) = E1c g' + E0c,  E1c = G_U1(U) + (5km + q/4),  E0c = H(U) - H(u0) - c0' w - c0 g + (5km+q/4)' g,  H := D G with U1' = a2 (no g'),
#   hence B0 = sup|H_U0(U) - c0'|, B1 = sup|H_U1(U) - c0 + (5km+q/4)'| (mean value), and E0c(0,0) = 0.
Gf = sp.Function('G')
c0 = sp.Rational(1, 4) - 3 * m ** 2 - k * p
xx_ = sp.Symbol('xx_')
def Dfun(e, RR):
    """total derivative: every symbol s in RR becomes s(xx_) with s' = RR[s]; works through undefined functions (chain rule by sympy)"""
    fs = {s_: sp.Function(str(s_) + '_f')(xx_) for s_ in RR}
    E = e.subs(fs, simultaneous=True); dE = sp.diff(E, xx_)
    dE = dE.subs({sp.Derivative(fs[s_], xx_): RR[s_].subs(fs, simultaneous=True) for s_ in RR}).doit()
    return dE.subs({fs[s_]: s_ for s_ in fs}, simultaneous=True)
RF = {us: a1, a1: a2, w: g, g: g1, a2: a3, a3: a4, a4: a5, k: m, m: p / 4, p: -k * q, q: 1 + k * p}
d000 = Gf(us + w, a1 + g, a2, a3, a4, k, m, p, q) - Gf(us, a1, a2, a3, a4, k, m, p, q) - c0 * w + (5 * k * m + q / 4) * g
Dd0 = sp.expand(Dfun(d000, RF))
E1c = sp.diff(Dd0, g1); E0c = Dd0.subs(g1, 0)
HX0, HX1 = sp.symbols('HX0 HX1')
GU1 = sp.diff(Gf(HX0, HX1, a2, a3, a4, k, m, p, q), HX1).subs({HX0: us + w, HX1: a1 + g}, simultaneous=True)
res['I14_E1c_eq_GU1_plus_drift'] = Z(E1c - (GU1 + 5 * k * m + q / 4))
res['I14_linear_in_g1'] = Z(sp.diff(Dd0, g1, 2))
res['I14_E0c_at_0'] = Z(E0c.subs({w: 0, g: 0}))
RH = {HX0: HX1, HX1: a2, a2: a3, a3: a4, a4: a5, k: m, m: p / 4, p: -k * q, q: 1 + k * p}
Hgen = Dfun(Gf(HX0, HX1, a2, a3, a4, k, m, p, q), RH)
def H_of(X0, X1): return Hgen.subs({HX0: X0, HX1: X1}, simultaneous=True)
c0p = -6 * m * (p / 4) - m * p - k * (-k * q)
dq4p = 5 * (m * m + k * p / 4) + (1 + k * p) / 4
res['I14_c0prime'] = Z(c0p - Dfun(c0, RH))
res['I14_drift_prime'] = Z(dq4p - Dfun(5 * k * m + q / 4, RH))
E0c_form = H_of(us + w, a1 + g) - H_of(us, a1) - c0p * w - c0 * g + dq4p * g
res['I14_E0c_form'] = Z(E0c - E0c_form)
print("I14", res['I14_E1c_eq_GU1_plus_drift'], res['I14_linear_in_g1'], res['I14_E0c_at_0'], res['I14_E0c_form'], res['I14_c0prime'], res['I14_drift_prime'], res['I14_c0prime'], res['I14_drift_prime'], round(time.time() - t0), flush=True)
# I15 N_low' table with generic coefficient functions of x
x = sp.symbols('x')
F_ = lambda n: sp.Function(n)(x)
D010, D100, D200, D110, D300, D000, G = [F_(n) for n in ('d010', 'd100', 'd200', 'd110', 'd300', 'd000', 'g')]
dx = lambda f, n=1: sp.diff(f, x, n)
Nlow = D010 * dx(G, 2) + D100 * dx(G) + D200 * dx(G) ** 2 + D110 * dx(G) * dx(G, 2) + D300 * dx(G) ** 3 + D000
table = (D010 * dx(G, 3) + (dx(D010) + D100) * dx(G, 2) + dx(D100) * dx(G) + 2 * D200 * dx(G) * dx(G, 2) + dx(D200) * dx(G) ** 2
         + D110 * (dx(G, 2) ** 2 + dx(G) * dx(G, 3)) + dx(D110) * dx(G) * dx(G, 2) + 3 * D300 * dx(G) ** 2 * dx(G, 2) + dx(D300) * dx(G) ** 3 + dx(D000))
res['I15_Nlow_prime_table'] = Z(dx(Nlow) - table); print("res['I15_Nlow_prime_table']", res['I15_Nlow_prime_table'][:60], round(time.time() - t0), flush=True)
# I16 v-equation
kf, mf, pf, qf, th = [F_(n) for n in ('k', 'm', 'p', 'q', 'theta')]
ode = {sp.Derivative(pf, x): -kf * qf, sp.Derivative(qf, x): 1 + kf * pf, sp.Derivative(kf, x): mf, sp.Derivative(mf, x): pf / 4, sp.Derivative(th, x): kf}
wf, u0f, zf, cf_, NL, PhiF, ps3, h2f = [F_(n) for n in ('w', 'u0', 'zeta', 'c', 'Nlow', 'PhiF', 'psi3', 'h2')]
aa, bb, be, alp = sp.symbols('a b beta alpha')
c0f = sp.Rational(1, 4) - 3 * mf ** 2 - kf * pf
L = lambda f: -dx(f, 4) - kf ** 2 * dx(f, 2) - (5 * kf * mf + qf / 4) * dx(f) + c0f * f
L1 = lambda f: -dx(f, 4) - kf ** 2 * dx(f, 2) - (7 * kf * mf + qf / 4) * dx(f) + (-8 * mf ** 2 - sp.Rational(5, 2) * kf * pf) * f
uf = u0f + wf
alg = lambda u: 1 - kf * u
A1 = lambda u: -u + dx(u) * qf / alg(u)
B1 = lambda u: dx(u) * sp.cos(th) / alg(u)
M = -be * ps3 - 2 * alp * be * h2f - aa * pf - bb * (-sp.sin(th)) + aa * A1(uf) + bb * B1(uf)
Nw = cf_ * dx(wf, 4) + NL
wtau = L(wf) + PhiF + Nw + M
gf = dx(wf); v = zf * gf
Qu = qf / alg(uf); Cu = sp.cos(th) / alg(uf)
dQ = qf * kf * uf / (1 - kf * uf); dC = sp.cos(th) * kf * uf / (1 - kf * uf)
EA = dx(u0f) * qf * kf / (alg(uf) * alg(u0f)); EB = dx(u0f) * kf * sp.cos(th) / (alg(uf) * alg(u0f))
K3 = 3 * dx(zf) * dx(gf, 2) + 3 * dx(zf, 2) * dx(gf) + dx(zf, 3) * gf
K4 = 4 * dx(zf) * dx(gf, 3) + 6 * dx(zf, 2) * dx(gf, 2) + 4 * dx(zf, 3) * dx(gf) + dx(zf, 4) * gf
T = [L1(v),                                                     # T1
     -(L1(v) - zf * L1(gf)),                                    # T2 = -[L1, zeta] g
     dx(cf_ * dx(v, 3)),                                        # T3
     -cf_ * K4 - dx(cf_) * K3,                                  # T4
     zf * dx(NL),                                               # T5
     zf * dx(c0f) * wf,                                         # T6
     zf * dx(PhiF),                                             # T7
     aa * (qf * dx(v) + (dx(qf) - 1) * v) + bb * (sp.cos(th) * dx(v) - kf * sp.sin(th) * v),     # T8
     -dx(zf) * (aa * Qu + bb * Cu) * gf,                        # T9
     aa * zf * dx(A1(u0f)) + bb * zf * dx(B1(u0f)),             # T10
     -be * zf * (dx(ps3) + 2 * alp * dx(h2f)) - aa * zf * dx(pf) - bb * zf * dx(-sp.sin(th)),   # T11
     aa * (dx(dQ * v) + zf * dx(EA * wf)) + bb * (dx(dC * v) + zf * dx(EB * wf))]               # T12
resid = zf * dx(wtau) - sum(T)
resid = resid.doit().subs(ode).doit().subs(ode)
res['I16_v_equation'] = Z(resid); print("res['I16_v_equation']", res['I16_v_equation'][:60], round(time.time() - t0), flush=True)
# ---------------------------------------------------------------- FD cross-check (mpmath)
import mpmath as mp
mp.mp.dps = 40
random.seed(12345)
symsF = (U0, U1, a2, a3, a4, w2, w3, w4, k, m, p, q)
Ffun = sp.lambdify(symsF, Fpert, 'mpmath')
cfun = {nm: sp.lambdify((U0, U1, a2, a3, a4, k, m, p, q), coef[mon], 'mpmath') for mon, nm in NAMES.items()}
dafun = {nm: sp.lambdify((U0, U1, a2, a3, a4, a5, k, m, p, q), sp.expand(Dx(coef[mon])).subs(g1, 0), 'mpmath') for mon, nm in NAMES.items()}
dbfun = {nm: sp.lambdify((U0, U1, a2, a3, a4, k, m, p, q), sp.diff(coef[mon], U1), 'mpmath') for mon, nm in NAMES.items()}
Flow = Fpert.subs({w2: 0, w3: 0, w4: 0})
RULES_H = {U0: U1, U1: a2, a2: a3, a3: a4, a4: a5, k: m, m: p / 4, p: -k * q, q: 1 + k * p}
Hc = sum(sp.diff(Flow, v_) * r_ for v_, r_ in RULES_H.items())                 # H(U0, U1) = D G with U1' = a2 (concrete G)
Hfn = sp.lambdify((U0, U1, a2, a3, a4, a5, k, m, p, q), Hc, 'mpmath')
GU1f = sp.lambdify((U0, U1, a2, a3, a4, k, m, p, q), sp.diff(Flow, U1), 'mpmath')
def E1f(us_, a1_, w_, g_, a2_, a3_, a4_, K_, M_, P_, Q_): return GU1f(us_ + w_, a1_ + g_, a2_, a3_, a4_, K_, M_, P_, Q_) + 5 * K_ * M_ + Q_ / 4
def E0f(us_, a1_, w_, g_, a2_, a3_, a4_, a5_, K_, M_, P_, Q_):
    c0_ = mp.mpf(1) / 4 - 3 * M_ ** 2 - K_ * P_; c0p_ = -mp.mpf(5) / 2 * M_ * P_ + K_ ** 2 * Q_; dq_ = 5 * (M_ ** 2 + K_ * P_ / 4) + (1 + K_ * P_) / 4
    return Hfn(us_ + w_, a1_ + g_, a2_, a3_, a4_, a5_, K_, M_, P_, Q_) - Hfn(us_, a1_, a2_, a3_, a4_, a5_, K_, M_, P_, Q_) - c0p_ * w_ - c0_ * g_ + dq_ * g_

fd = []
for trial in range(3):
    y0 = [mp.mpf(random.uniform(-0.3, 0.3)) for _ in range(3)] + [mp.mpf(random.uniform(5, 15))]     # k, m, p, q
    sol = mp.odefun(lambda xx, y: [y[1], y[2] / 4, -y[0] * y[3], 1 + y[0] * y[2]], 0, y0)
    cu = [mp.mpf(random.uniform(-0.2, 0.2)) for _ in range(6)]; cw = [mp.mpf(random.uniform(-0.1, 0.1)) for _ in range(3)]
    u0x = lambda xx, n=0: mp.diff(lambda t: sum(cu[i] * t ** i for i in range(6)) + cu[5] * mp.sin(t), xx, n)
    wx = lambda xx, n=0: mp.diff(lambda t: cw[0] * mp.sin(1.3 * t + cw[1]) + cw[2] * t ** 2, xx, n)
    def Nlow_at(xx):
        K_, M_, P_, Q_ = sol(xx)
        a = [u0x(xx, j) for j in range(5)]; ww = [wx(xx, j) for j in range(5)]
        Fu = Ffun(a[0] + ww[0], a[1] + ww[1], a[2], a[3], a[4], ww[2], ww[3], ww[4], K_, M_, P_, Q_)
        F0 = Ffun(a[0], a[1], a[2], a[3], a[4], 0, 0, 0, K_, M_, P_, Q_)
        c0_ = mp.mpf(1) / 4 - 3 * M_ ** 2 - K_ * P_
        Lw_ = -ww[4] - K_ ** 2 * ww[2] - (5 * K_ * M_ + Q_ / 4) * ww[1] + c0_ * ww[0]
        cc_ = cfun['c'](a[0] + ww[0], a[1] + ww[1], a[2], a[3], a[4], K_, M_, P_, Q_)
        return Fu - F0 - Lw_ - cc_ * ww[4]
    for xx in (mp.mpf(random.uniform(0.2, 0.8)), mp.mpf(random.uniform(0.9, 1.5))):
        num = mp.diff(Nlow_at, xx)
        K_, M_, P_, Q_ = sol(xx)
        a = [u0x(xx, j) for j in range(6)]; ww = [wx(xx, j) for j in range(5)]
        Uv = (a[0] + ww[0], a[1] + ww[1], a[2], a[3], a[4]); gg = ww[1]; gp, gpp, gppp = ww[2], ww[3], ww[4]
        D = {nm: dafun[nm](*Uv, a[5], K_, M_, P_, Q_) + dbfun[nm](*Uv, K_, M_, P_, Q_) * gp for nm in NAMES.values()}
        C = {nm: cfun[nm](*Uv, K_, M_, P_, Q_) for nm in NAMES.values()}
        Dd000 = E1f(a[0], a[1], ww[0], gg, a[2], a[3], a[4], K_, M_, P_, Q_) * gp + E0f(a[0], a[1], ww[0], gg, a[2], a[3], a[4], a[5], K_, M_, P_, Q_)
        tab = (C['d010'] * gppp + (D['d010'] + C['d100']) * gpp + D['d100'] * gp + 2 * C['d200'] * gp * gpp + D['d200'] * gp ** 2
               + C['d110'] * (gpp ** 2 + gp * gppp) + D['d110'] * gp * gpp + 3 * C['d300'] * gp ** 2 * gpp + D['d300'] * gp ** 3 + Dd000)
        # c' (for T3/T4): total derivative of c along the solution vs numerical
        cnum = mp.diff(lambda t: cfun['c'](u0x(t) + wx(t), u0x(t, 1) + wx(t, 1), u0x(t, 2), u0x(t, 3), u0x(t, 4), *sol(t)), xx)
        fd.append(dict(x=float(xx), Nlow_prime_num=mp.nstr(num, 15), table=mp.nstr(tab, 15), rel=float(abs(num - tab) / max(abs(num), mp.mpf('1e-30'))),
                       c_prime_rel=float(abs(cnum - D['c']) / max(abs(cnum), mp.mpf('1e-30')))))
res_fd_ok = all(r['rel'] < 1e-7 and r['c_prime_rel'] < 1e-7 for r in fd)      # consistency test: agreement limited by the mpmath ODE solution (odefun) used for the profile; a missing or wrong term gives O(1)
flat = [res['I13_c_is_1_minus_Jm4'], res['I13_d110'], res['I13_d300'], res['I14_E1c_eq_GU1_plus_drift'], res['I14_linear_in_g1'], res['I14_E0c_at_0'], res['I14_E0c_form'], res['I14_c0prime'], res['I14_drift_prime'],
        res['I15_Nlow_prime_table'], res['I16_v_equation']] + list(i12.values())
ok = all(s == '0' for s in flat) and res_fd_ok
out = dict(args=args_block([GSP]), residuals=res, fd_check=fd, fd_ok=res_fd_ok, all_zero=ok, runtime_s=round(time.time() - t0, 1))
h = dump(os.path.join(OUT, 's1_identities.json'), out)
print(json.dumps(res, indent=1)); print('FD max rel', max(r['rel'] for r in fd), max(r['c_prime_rel'] for r in fd)); print('ALL ZERO' if ok else 'FAIL', h)
