#!/usr/bin/env python3
"""Stage X4 (FLUX2_DESIGN.md sec. 3, constants #9, #10): data on the outer ramp from the (A) chain v36 and the chart/profile data of X3.
Inputs are taken from out/inputs.json (stage X0; every file re-hashed here and compared).

Notation: S = sigma-shell of OR, sigma_min = 54.75 <= inf S, sigma_max = 57.76 >= sup S, y_max = 57.76 >= sup Y(S) (X3: e <= 1e-40);
ds = dsigma/sigma; y_k(tau) := ||D^k g||_{L^2(S, ds)}; Y_k := ||D^k g||_{L^2(S x I, ds dtau)}; M_l^2 = int_OR (d^l g)^2 dx (dx = dsigma on S).
[hand] steps (each one line):
 (conv) Minkowski in L^2(ds) with sigma^l d^l = sum_m s(l,m) D^m (X1 (e)):  M_0 <= sigma_max^{1/2} y_0, M_1 <= sigma_min^{-1/2} y_1,
        M_2 <= sigma_min^{-3/2} (y_2 + y_1), M_3 <= sigma_min^{-5/2} (y_3 + 3 y_2 + 2 y_1), M_4 <= sigma_min^{-7/2} (y_4 + 6 y_3 + 11 y_2 + 6 y_1);
        the same in L^2(S x I) for window norms mu_l := (int_I M_l^2)^{1/2}.
 (chart) D^k g = z^A_k o Y - D^k d - chi_k (chart lemma (d), every k; X3), so y_k <= (1 + rho)^{1/2} ||z^A_k||_{L^2(Y(S), dy/y)} + ||D^k d|| + ||chi_k||
        (rho = X3 rho_meas: dsigma/sigma = (1 + eps_c) dy/y);  |chi_k| <= Abar_k + sum_i Bbar_ki (|d^i g| + |P^(i)|) + Bbar_422 (2 g''^2 + 2 P''^2) (k = 4);
        ||d^i g||_{L^2(ds)} <= sigma_min^{-1/2} M_i.
 (energies) pointwise: int_{Y(S)} (z^A_k)^2 dy/y <= (y_max/28) beta_A Ehat_k, beta_A = 2, k <= 2 (core weight omega = 28/y on Y(S));
        windows (B Phi): int_I int_{Y(S)} (z^A_{k+2})^2 dy/y dtau <= (y_max^5/28) 2 (C0_k + C1_k T), k = 1, 2 (export weight y^-4 (28/y) >= 28 y_max^-5).
 (Agmon) s_2^2 <= M_2^2/3 + 2 M_2 M_3 on [57, 60] (length 3); int_I s_2^2 <= sig2c + sig2s T, sig2c = Mbar_2 m3c/lam_s,
        sig2s = Mbar_2^2/3 + Mbar_2 (lam_s + m3s/lam_s), lam_s = 1.
 (affine) a window norm bounded by sum_i w_i sqrt(a_i + b_i T) is squared by Cauchy-Schwarz with fixed rational weights theta_i > 0,
        sum theta_i <= 1: (sum x_i)^2 <= sum x_i^2/theta_i.  theta_i are chosen deterministically from the float sizes at T = 1 (free parameters).
 (circular) the relative chart parts (Bbar_22, Bbar_33, Bbar_44, Bbar_43 ~ 1e-36) are absorbed by 1/(1 - beta) factors.
Outputs (all [cert], at beta_A = 2 and the (B Phi) level): Mbar_l, Nbar_l (l <= 2), m3c, m3s, m4c, m4s, sig2c, sig2s; regressions against the v36 pipeline.
Output out/flux2_data.json."""
import sys, os, json, math
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import f2common as FC
from flint import arb, ctx
ctx.prec = 200
T_ = FC.Timer()
CHK = {}


def ok(name, val):
    CHK[name] = bool(val); print('%-80s %s' % (name, bool(val)), flush=True); assert val, name


INP = json.load(open(os.path.join(FC.OUT, 'inputs.json')))
for p, h in INP['pinned'].items():
    assert FC.sha(p) == h, 'input changed since X0: %s' % p
ASM = INP['v36']['asm_design_N2ref4']
a = json.load(open(ASM))['res']['tight']
X2 = json.load(open(os.path.join(FC.OUT, 'flux2_coef.json')))['res']
X3 = json.load(open(os.path.join(FC.OUT, 'chart4_shell.json')))['res']
up = lambda r: FC.from_raw_up(r['up'])
BETA = arb(2)
SMIN = arb('54.75'); SMAX = arb('57.76'); YMAX = arb('57.76')
ok('X3 shell is [54.75, 57.76] and X2 sigma range inside it', X3['shell'] == ['54.75', '57.76'] and X2['profile']['sigma_min']['dn'] >= 54.75 and X2['profile']['sigma_max']['up'] <= 57.76)
ok('X3 chart shift e <= 1e-30 (so Y(S) within [54.75 - 1e-30, 57.76 + 1e-30]; y_max = 57.76 + e used below)', X3['e']['up'] <= 1e-30)
YMAX = YMAX + up(X3['e'])
RHO = up(X3['rho_meas']); SQR = (1 + RHO).sqrt()
E = [FC.load_arb(e).upper() for e in a['E']]
E = [arb(e) for e in E]
EX = {k: dict(C0=arb(FC.load_arb(a['export'][k]['C0']).upper()), C1=arb(FC.load_arb(a['export'][k]['C1']).upper())) for k in ('1', '2')}
ok('export statements have the weight y^-4 (28/y) on [ln 28, ln 350] for z_3 (k=1) and z_4 (k=2)',
   'y^-4 (28/y) z_3^2' in a['export']['1']['statement'] and 'y^-4 (28/y) z_4^2' in a['export']['2']['statement'] and 'ln28}^{ln350}' in a['export']['2']['statement'])
ok('level-2 export retained fraction r_2 = 1/1000 (A_interface r2_export)', abs(float(FC.load_arb(a['export']['2']['r']).mid()) - 1e-3) < 1e-15)
C0_2 = EX['2']['C0']
ok('STOP RULE: C0_2 = %.4e < 4e3 (kill 4e4; design: raise lambda_2 above 4e3)' % float(C0_2.upper()), C0_2.upper() < 4e3)
# chart / profile data on S (X3)
Abar = [up(r) for r in X3['Abar']]; Bb = {k: up(r) for k, r in X3['Bbar'].items()}; B422 = up(X3['B422'])
Pint = [up(r) for r in X3['P_int_ds']]; P2sq = up(X3['P2sq_int_ds']); Dint = [up(r) for r in X3['Dd_int_ds']]; Sds = up(X3['S_ds'])
sq = lambda v: arb(v).sqrt()


# ------------------------------------------------------------------ pointwise levels (beta_A = 2)
def zpt(k): return SQR * sq((YMAX / 28) * BETA * E[k])                # (1+rho)^{1/2} ||z_k||_{L2(Y(S))}
ypt = {}
for k in (0, 1):
    ypt[k] = zpt(k) + sq(Dint[k]) + Abar[k] * sq(Sds)              # chi_0 = A_0, chi_1 = A_1 (B_0i = B_1i = 0, X3 asserted)
Mb = {}
Mb[0] = SMAX.sqrt() * ypt[0]
Mb[1] = ypt[1] / SMIN.sqrt()
# level 2: y_2 <= r_2 + Bbar_22 sigma_min^{-1/2} M_2,  M_2 <= sigma_min^{-3/2} (y_2 + y_1)
r2 = zpt(2) + sq(Dint[2]) + Abar[2] * sq(Sds) + Bb['22'] * sq(Pint[2])
b22 = Bb['22'] / SMIN.sqrt()
fac2 = SMIN ** arb(-1.5)
ok('level-2 chart absorption: sigma_min^{-3/2} Bbar_22 sigma_min^{-1/2} < 1e-30', (fac2 * b22).upper() < 1e-30)
Mb[2] = fac2 * (r2 + ypt[1]) / (1 - fac2 * b22)
ypt[2] = r2 + b22 * Mb[2]
Mb = {k: FC.aup(v) for k, v in Mb.items()}
W57 = up(X2['W57'])
Nb = {k: FC.aup((W57 * Mb[k] ** 2).sqrt()) for k in Mb}


# ------------------------------------------------------------------ windows: term lists  [(w, a, b)] meaning  sum w sqrt(a + b T)
def zwin(k):        # k = 3, 4: (1+rho)^{1/2} ||z_k||_{L2(Y(S) x I)} <= (1+rho)^{1/2} sqrt((y_max^5/28) 2 (C0 + C1 T))
    ex = EX[str(k - 2)]; f = YMAX ** 5 / 28 * 2
    return [(SQR, f * ex['C0'], f * ex['C1'])]


def rate(c):        # sqrt(T) * c  ==  1 * sqrt(0 + c^2 T)
    return [(arb(1), arb(0), c * c)]


def scale(L, s_): return [(w * s_, a_, b_) for (w, a_, b_) in L]


def floats(L, Tref=1.0): return [float((w * (a_ + b_ * Tref).sqrt()).upper()) for (w, a_, b_) in L]


def square(L, Tref=1.0):
    """(sum w_i sqrt(a_i + b_i T))^2 <= sum w_i^2 (a_i + b_i T)/theta_i, theta_i = rational, sum theta_i <= 1 (exact)"""
    f = floats(L, Tref); tot = sum(f)
    th = [Fraction(max(1, int(1e9 * x / tot)), 10**9 + len(f)) for x in f]
    assert sum(th) <= 1
    c = arb(0); s_ = arb(0)
    for (w, a_, b_), t in zip(L, th):
        ta = FC.frac_arb(t); c += w * w * a_ / ta; s_ += w * w * b_ / ta
    return FC.aup(c), FC.aup(s_), [str(t) for t in th]


# Y_k for k <= 2 as window norms: sqrt(T) * ybar_k
Ywin = {k: rate(ypt[k]) for k in (0, 1, 2)}
# level 3: Y_3 <= zwin(3) + sqrt(T) (||D^3 d|| + Abar_3 |S|^1/2 + Bbar_32 (sigma_min^-1/2 Mbar_2 + ||P''||) + Bbar_33 ||P'''||) + Bbar_33 sigma_min^-1/2 mu_3
c3 = sq(Dint[3]) + Abar[3] * sq(Sds) + Bb['32'] * (Mb[2] / SMIN.sqrt() + sq(Pint[2])) + Bb['33'] * sq(Pint[3])
b33 = Bb['33'] / SMIN.sqrt()
f3 = SMIN ** arb(-2.5)
ok('level-3 chart absorption < 1e-30', (f3 * b33).upper() < 1e-30)
# mu_3 <= f3 (Y_3 + 3 Y_2 + 2 Y_1)  and Y_3 <= R_3 + b33 mu_3  =>  mu_3 <= f3 (R_3 + 3 Y_2 + 2 Y_1)/(1 - f3 b33)
den3 = 1 - f3 * b33
L3 = scale(zwin(3) + rate(c3) + scale(Ywin[2], arb(3)) + scale(Ywin[1], arb(2)), f3 / den3)
m3c, m3s, th3 = square(L3)
# s_2 window bound
lam_s = arb(1)
sig2c = FC.aup(Mb[2] * m3c / lam_s); sig2s = FC.aup(Mb[2] ** 2 / 3 + Mb[2] * (lam_s + m3s / lam_s))
# level 4: Y_4 <= zwin(4) + sqrt(T) c4 + Bbar_44 s^-1/2 mu_4 + Bbar_43 s^-1/2 mu_3 + Bbar_42 s^-1/2 sqrt(T) Mbar_2 + 2 Bbar_422 (||g''^2|| + ||P''^2||)
#   ||g''^2||_{L2(S x I, ds)} <= sigma_min^{-1/2} Mbar_2 (int_I s_2^2)^{1/2}
c4 = sq(Dint[4]) + Abar[4] * sq(Sds) + Bb['42'] * (Mb[2] / SMIN.sqrt() + sq(Pint[2])) + Bb['43'] * sq(Pint[3]) + Bb['44'] * sq(Pint[4]) + 2 * B422 * sq(P2sq)
b44 = Bb['44'] / SMIN.sqrt(); b43 = Bb['43'] / SMIN.sqrt()
f4 = SMIN ** arb(-3.5)
ok('level-4 chart absorption < 1e-30 (B_44 circular term, sec. 3.4)', (f4 * b44).upper() < 1e-30)
den4 = 1 - f4 * b44
L4 = zwin(4) + rate(c4) + scale(L3, b43) + [(2 * B422 * Mb[2] / SMIN.sqrt(), sig2c, sig2s)]
# mu_4 <= f4 (Y_4 + 6 Y_3 + 11 Y_2 + 6 Y_1), Y_3 <= R_3 + b33 mu_3 with R_3 = zwin(3) + rate(c3); mu_3 <= L3
Y3 = zwin(3) + rate(c3) + scale(L3, b33)
L4full = scale(L4 + scale(Y3, arb(6)) + scale(Ywin[2], arb(11)) + scale(Ywin[1], arb(6)), f4 / den4)
m4c, m4s, th4 = square(L4full)
# W-norm versions
n3c, n3s, n4c, n4s = [FC.aup(W57 * v) for v in (m3c, m3s, m4c, m4s)]

# ------------------------------------------------------------------ regressions against the v36 pipeline (inputs read correctly)
CF = json.load(open(INP['v36']['lemmaS_cflux']))
W57o = FC.from_raw_up(CF['W57_upper']); k0o = [FC.from_raw_up(v) for v in CF['kappa_I0']]
conv = arb(60) / 28
Eu = [conv * BETA * e for e in E]
N2o = [W57o * arb('57.75') * Eu[0], W57o * Eu[1] / arb('54.75'), W57o * 2 * (Eu[2] + Eu[1]) / arb('54.75') ** 3]     # assembly_arb.f4b, N_l^2
K0w_old = sum((k0o[l] * N2o[l] for l in range(3)), arb(0))
K0w_pipe = 2 * FC.load_arb(a['flux']['B0_d']['K0_w'])                  # f4b on beta E = beta * K0_w (linear)
ok('regression: f4b formula on (E, beta = 2, lemmaS_cflux) reproduces the v36 pipeline K0_w at beta = 2 (rel <= 1e-9)',
   abs(float((K0w_old / K0w_pipe).mid()) - 1) <= 1e-9)
ratio_N = [FC.fup(Nb[l] ** 2 / (N2o[l] * W57 / W57o)) for l in range(3)]
ok('X4 Nbar_l^2 (with d, chart) <= 1.05 x old f4b N_l^2 rescaled to the same W(57) (sanity; old has no d/chart)', all(r <= 1.05 for r in ratio_N))
# m3 against the v36 "honest I3" (different rules; reported only)
I3c_pipe = FC.load_arb(a['flux']['B0_d']['I3c_w']); I3s_pipe = FC.load_arb(a['flux']['B0_d']['I3s_w'])
k2o3 = FC.from_raw_up(CF['kappa_I2'][3])
R = dict(beta_A='2', lam_s='1', sigma_min='54.75', sigma_max='57.76', y_max=FC.rec(YMAX, 'cert'),
         E_hat=[FC.rec(e, 'cert', 'v36 asm_design_N2ref4 res.tight.E[%d]' % k) for k, e in enumerate(E)],
         export={k: dict(C0=FC.rec(v['C0'], 'cert'), C1=FC.rec(v['C1'], 'cert')) for k, v in EX.items()},
         Mbar={str(k): FC.rec(v, 'cert', 'sup_tau M_%d at beta_A = 2, with d and chart' % k) for k, v in Mb.items()},
         Nbar={str(k): FC.rec(v, 'cert', 'sup_tau N_%d = (int_OR W (g^(%d))^2)^1/2 <= W(57)^1/2 Mbar_%d' % (k, k, k)) for k, v in Nb.items()},
         m3c=FC.rec(m3c, 'cert'), m3s=FC.rec(m3s, 'cert'), m4c=FC.rec(m4c, 'cert'), m4s=FC.rec(m4s, 'cert'),
         n3c=FC.rec(n3c, 'cert', 'int_I N_3^2 <= n3c + n3s T'), n3s=FC.rec(n3s, 'cert'), n4c=FC.rec(n4c, 'cert'), n4s=FC.rec(n4s, 'cert'),
         sig2c=FC.rec(sig2c, 'cert', 'int_I s_2^2 <= sig2c + sig2s T'), sig2s=FC.rec(sig2s, 'cert'),
         theta=dict(level3=th3, level4=th4),
         regression=dict(K0w_old_formula=FC.fup(K0w_old), K0w_pipeline_beta2=FC.fup(K0w_pipe), ratio_Nbar2_new_over_old_sameW=ratio_N,
                         I3c_pipeline_w=FC.fup(I3c_pipe), I3s_pipeline_w=FC.fup(I3s_pipe), I3c_new_equiv=FC.fup(k2o3 * W57o * m3c), I3s_new_equiv=FC.fup(k2o3 * W57o * m3s),
                         note='I3 of the pipeline = kappa_I2[3] W57_old N3-pipeline (v3.1 rules: 2 a_min conversion, no chart); reported only'))
out = dict(stage='X4 flux2_data', res=R, checks=CHK, all_pass=all(CHK.values()),
           args=FC.args_block([os.path.abspath(__file__), os.path.join(FC.CODE, 'f2common.py'), os.path.join(FC.OUT, 'inputs.json'), os.path.join(FC.OUT, 'flux2_coef.json'),
                               os.path.join(FC.OUT, 'chart4_shell.json'), ASM, INP['v36']['lemmaS_cflux']]))
out['timing'] = dict(runtime_s=T_.s())
h = FC.dump(os.path.join(FC.OUT, 'flux2_data.json'), out)
print('Mbar', {k: FC.fup(v) for k, v in Mb.items()}, 'Nbar', {k: FC.fup(v) for k, v in Nb.items()})
print('m3 %.4e + %.4e T ; m4 %.4e + %.4e T ; s2 %.4e + %.4e T' % tuple(FC.fup(v) for v in (m3c, m3s, m4c, m4s, sig2c, sig2s)))
print('regression K0w old-formula/pipeline', float((K0w_old / K0w_pipe).mid()), 'Nbar^2 new/old', ratio_N)
print('X4 done sha', h, '(%.1fs)' % T_.s())
