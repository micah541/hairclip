#!/usr/bin/env python3
"""Stage X5 (FLUX2_DESIGN.md sec. 4, 5, 7; constants #11, #12): flux lemma v2 assembled, and the absorption conditions C0-C9 checked
on the stored fin13v4b rate tables (both alpha points) and the v36 hjw_rep.  arb, outward rounding; margins are certified LOWER bounds.

Level 0 (pointwise):  Kbar0 = sum_l kappa0_l Nbar_l^2 + K0nl + K0w,
    K0nl = (Nbar_1 + KAP2 Nbar_0) (cN2 Nbar_2 + cN1 Nbar_1 + cN0 Nbar_0 + cw wbar omega_W),  K0w = Nbar_0 wbar omega_W cK0w.
Level 2:  frak D_pt = a_2 Nbar_2 + a_1 Nbar_1 + a_0 Nbar_0;  frak X_pt = X.N2 Nbar_2 + X.N1 Nbar_1 + X.N0 Nbar_0 + X.w wbar omega_W;
    K2pt = D_pt X_pt + (lam1/2) D_pt^2 + (lam2/2) X_pt^2;   int_I phi2 <= C'' + S'' T,
    C'' = aX [2 x3^2 n3c + 2 (d110 Nbar_2)^2 sig2c] + aD [2 n4c + 2 a3^2 n3c],  aX = 1/(2 lam1) + 1/2, aD = 1/(2 lam2) + 1/2;  S'' likewise.
Verifier form: mu = 0.08, eps_I = 3e-5 (decimal, as rhs_iv._A), F0 = mu (FLUX0 + eps_I FLUX2);  s = F0 - mu Kbar0 - mu eps_I K2pt;
    C' = mu eps_I C'', S' = mu eps_I S''.
Conditions (sec. 7.3; [hand] readings of closure.py: march cell k certifies tot/Lambda <= U_k = rate_U on [Ly_lo_k, Ly_hi_k], barrier
    L' = 2 (U_k + 1e-9 max(1, |U_k|)) L, L monotone on the step; tail: rate_s <= -rho on [L*, (1 + ETA) L*], L_F = 2 F0/rho):
  C0  s > 0;   C1  S' <= s;
  C2  for every step k: max(L_k, L_{k+1}) + 2C' <= Ly_hi_k   (L_0 = Ly1 (1 + 1e-12), L_{k+1} <= sqrtL_next_k^2);
  C3  for every step k: U_k^+ 2C' < 1e-9 max(1, |U_k|) Ly_lo_k;
  C4  2C' <= ETA L_F,  L_F >= 2 F0/|rate_s_U|;
  C5  alpha condition whole tail: lhs >= rhs_whole_tail (1 + 2C'/L_F);
  C6  C' <= (kappa_J - kappa_J_certified) X_const   (v36 hjw_rep; X_const of the a = 1.5e-4 table);
  C7  consumers: march ones read Ly_hi (covered by C2); the mod_sups tail row (L*, A <= sqrt L*) with Lambda <= L* + 2C': reserve of the
      6-digit rounding; floor consumers of sqrt_LF need relative reserve C'/L_F (listed, from V36_REPORT.md);
  C8  regularity U-1: OPEN (not checkable by script);  C9 bootstrap order: the constants above are at the (B Phi) level (doubled exports).
Stop rules (task): combined pointwise margin < 2, any of C0-C6 failing -> FAIL (non-zero exit after writing the JSON).
Output out/flux2_cert.json."""
import sys, os, json, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import f2common as FC
from flint import arb, ctx
ctx.prec = 200
T_ = FC.Timer()
INP = json.load(open(os.path.join(FC.OUT, 'inputs.json')))
for p, h in INP['pinned'].items():
    assert FC.sha(p) == h, 'input changed since X0: %s' % p
X2 = json.load(open(os.path.join(FC.OUT, 'flux2_coef.json')))['res']
X4 = json.load(open(os.path.join(FC.OUT, 'flux2_data.json')))['res']
up = lambda r: FC.from_raw_up(r['up'])
dn = lambda r: FC.from_raw_dn(r['dn'])
COND = {}; CHK = {}


def cond(name, val, margin=None, note=None):
    COND[name] = dict(pass_=bool(val), margin=margin, note=note)
    print('%-6s %-5s %s %s' % (name, 'PASS' if val else 'FAIL', ('margin x%.4g' % margin) if isinstance(margin, float) else '', note or ''), flush=True)


def ok(name, val):
    CHK[name] = bool(val); print('%-80s %s' % (name, bool(val)), flush=True); assert val, name


# ------------------------------------------------------------------ inputs
N = [up(X4['Nbar'][str(l)]) for l in range(3)]
kap = [up(r) for r in X2['kappa0']]
KAP2 = up(X2['KAP2']); WBAR = arb('0.5'); OMW = up(X2['omega_W'])
cn = X2['K0nl_coef']; cK0w = up(X2['K0w_coef'])
A_ = {int(k): up(v) for k, v in X2['frak_a'].items()}
XC = {k: up(v) for k, v in X2['X2_pt_coef'].items()}
x3 = up(X2['x3']); d110 = up(X2['d110'])
n3c, n3s, n4c, n4s, sig2c, sig2s = [up(X4[k]) for k in ('n3c', 'n3s', 'n4c', 'n4s', 'sig2c', 'sig2s')]
LAM1 = arb(1); LAM2 = arb(1)
MU = arb('0.08'); EPSI = arb('3e-5'); FLUX0 = arb('1.64e-20'); FLUX2 = arb('8.0e-20')
F0 = MU * (FLUX0 + EPSI * FLUX2)

# ------------------------------------------------------------------ level 0
K0lin = sum((kap[l] * N[l] ** 2 for l in range(3)), arb(0))
K0nl = (N[1] + KAP2 * N[0]) * (up(cn['N2']) * N[2] + up(cn['N1']) * N[1] + up(cn['N0']) * N[0] + up(cn['w']) * WBAR * OMW)
K0w = N[0] * WBAR * OMW * cK0w
K0 = FC.aup(K0lin + K0nl + K0w)
# ------------------------------------------------------------------ level 2
Dpt = A_[2] * N[2] + A_[1] * N[1] + A_[0] * N[0]
Xpt = XC['N2'] * N[2] + XC['N1'] * N[1] + XC['N0'] * N[0] + XC['w'] * WBAR * OMW
K2pt = FC.aup(Dpt * Xpt + LAM1 / 2 * Dpt ** 2 + LAM2 / 2 * Xpt ** 2)
aX = 1 / (2 * LAM1) + arb(1) / 2; aD = 1 / (2 * LAM2) + arb(1) / 2
Cpp = FC.aup(aX * (2 * x3 ** 2 * n3c + 2 * (d110 * N[2]) ** 2 * sig2c) + aD * (2 * n4c + 2 * A_[3] ** 2 * n3c))
Spp = FC.aup(aX * (2 * x3 ** 2 * n3s + 2 * (d110 * N[2]) ** 2 * sig2s) + aD * (2 * n4s + 2 * A_[3] ** 2 * n3s))
Cp = FC.aup(MU * EPSI * Cpp); Sp = FC.aup(MU * EPSI * Spp)
s_ = FC.adn(F0 - MU * K0 - MU * EPSI * K2pt)
margin_comb = FC.adn((FLUX0 + EPSI * FLUX2) / (K0 + EPSI * K2pt))
margin_K0 = FC.adn(FLUX0 / K0)
print('Kbar0 %.6e (lin %.6e, nl %.3e, w %.3e); K2pt %.6e; C\'\' %.4e S\'\' %.4e; C\' %.4e S\' %.4e; s %.6e; combined margin x%.4f' % tuple(
    float(v.upper()) for v in (K0, K0lin, K0nl, K0w, K2pt, Cpp, Spp, Cp, Sp, s_, margin_comb)), flush=True)
STOP = []
if not margin_comb.lower() >= 2: STOP.append('combined pointwise margin < x2')

# ------------------------------------------------------------------ conditions
two = 2 * Cp
cond('C0', s_.lower() > 0, float((s_ / F0).lower()), 's/F0 (s = F0 - mu Kbar0 - mu eps_I K2pt > 0)')
cond('C1', (Sp - s_).upper() <= 0, float((s_ / Sp).lower()), "S' <= s")
RTs = INP['rate_tables']; per = {}
for a_, p in sorted(RTs.items()):
    d = json.load(open(p)); st = d['steps']
    Lprev = FC.from_raw_up(d['Ly1']) * (1 + arb('1e-12')) * (1 + arb('1e-15'))
    c2m = None; c3m = None; rows = []
    for k, sk in enumerate(st):
        Lnext = FC.from_raw_up(sk['sqrtL_next']) ** 2
        Lmax = Lprev if Lprev.upper() >= Lnext.upper() else Lnext
        hi = arb(sk['Ly_hi']); lo = FC.from_raw_dn(sk['Ly_lo'])   # Ly_hi is the exact upper end of the certified cell
        m2 = (hi - Lmax) / two                                      # C2 margin: room / 2C'
        U = arb(sk['rate_U'])
        pad = arb('1e-9') * (abs(U) if abs(U).lower() >= 1 else arb(1))
        if U.upper() > 0:
            m3 = pad * lo / (abs(U) * two)                          # C3 margin: pad L / (U+ 2C')
        else:
            m3 = None
        rows.append(dict(k=k, Lmax=FC.fup(Lmax), Ly_hi=sk['Ly_hi'], Ly_lo=sk['Ly_lo'], rate_U=sk['rate_U'], C2_margin=FC.fdn(m2), C3_margin=(FC.fdn(m3) if m3 is not None else 'U<=0: no condition')))
        c2m = m2 if c2m is None or m2.lower() < c2m.lower() else c2m
        if m3 is not None: c3m = m3 if c3m is None or m3.lower() < c3m.lower() else c3m
        Lprev = Lnext
    tl = d['tail']; env = d['args']['env']; ETA = arb(env['ETA'])
    assert tl['rate_s_U'] < 0
    rho_hi = FC.from_raw_up(abs(tl['rate_s_U']))                   # rho <= |rate_s_U|
    LF_lo = FC.adn(2 * F0 / rho_hi)
    LF_tab = FC.from_raw_up(tl['sqrt_LF']) ** 2
    ok('alpha %s: L_F lower bound 2F0/|rate_s_U| <= sqrt_LF^2 (consistency)' % a_, LF_lo.upper() <= LF_tab.upper())
    m4 = ETA * LF_lo / two
    ac = tl['alpha_condition_whole_tail']
    lhs = FC.from_raw_dn(ac['lhs']); rhs = FC.from_raw_up(ac['rhs_whole_tail'])
    m5 = lhs / (rhs * (1 + two / LF_lo))
    per[a_] = dict(steps=rows, C2_min_margin=FC.fdn(c2m), C3_min_margin=(FC.fdn(c3m) if c3m is not None else None), L_F_lower=FC.fdn(LF_lo),
                   C4_margin=FC.fdn(m4), C5_margin=FC.fdn(m5), ETA=env['ETA'])
    cond('C2[%s]' % a_, c2m.lower() > 1, FC.fdn(c2m), 'min over steps of (Ly_hi - max L)/2C\'')
    cond('C3[%s]' % a_, c3m is None or c3m.lower() > 1, FC.fdn(c3m) if c3m is not None else None, 'min over steps with U>0 of pad L/(U 2C\')')
    cond('C4[%s]' % a_, m4.lower() >= 1, FC.fdn(m4), '2C\' <= ETA L_F')
    cond('C5[%s]' % a_, m5.lower() > 1, FC.fdn(m5), 'lhs >= rhs_whole_tail (1 + 2C\'/L_F)')
# C6
hj = json.load(open(INP['v36']['hjw_rep']))
kJ = FC.from_raw_dn(hj['kappa_J']); kJc = FC.from_raw_up(hj['kappa_J_certified'])
from fractions import Fraction
pad = FC.frac_arb(Fraction(hj['kappa_J_pad'])) if 'kappa_J_pad' in hj else (kJ - kJc)
Xc = FC.from_raw_dn(json.load(open(RTs['0.00015']))['unc']['X_const'])
ok('hjw X_const equals the rate-table unc.X_const', hj['X_const'] == json.load(open(RTs['0.00015']))['unc']['X_const'])
reserve = kJ - kJc                       # certified lower bound of kappa_J^used - kappa_J^cert (doubles moved one ulp outward)
ok('hjw kappa_J reserve > 0.99 x pad (pad applied; double rounding only)', (reserve - pad * arb('0.99')).lower() > 0)
m6 = reserve * Xc / Cp
cond('C6', m6.lower() > 1, FC.fdn(m6), "C' <= (kappa_J - kappa_J_cert) X_const (pad %s; reserve %.4e)" % (hj.get('kappa_J_pad'), float(reserve.lower())))
# C7 (scripted part): mod_sups tail row reserve.  The tail row is (CAB0 + MAB a*) sqrt(L*) + QAB qbl L* + PHIAB Fn(a*) + STRIP_AB; every L-dependent
# term has degree <= 1 in L, so the row at L* + 2C' is <= (1 + 2C'/L*) x row(L*).  Needed: (1 + 2C'/L*) row(L*) <= sup_ab (6-digit up),
# i.e. reserve := sup_ab/row(L*) - 1 >= 2C'/L*.  row(L*) is recomputed here from the constants recorded in final_alpha_points args.
MS = '/home/claude/s1spec/step3nl_v4b/runs/P/modsups_fin13v4b/final_alpha_points_5e-8.json'
ms = json.load(open(MS)); kc = ms['args']['constants']
sys.path.insert(0, os.path.join(FC.V4B, 'verifier'))
C7 = {}
try:
    for k_ in ('PYTHONHASHSEED', 'TAG', 'NEW', 'THETA', 'ADD', 'FIX', 'S0V', 'S1V', 'DELTA', 'MOD', 'MODA', 'ETA', 'FLUX0', 'FLUX2', 'LY1SQ', 'V4B'):
        v = json.load(open(RTs['0.00015']))['args']['env'].get(k_)
        if v is not None: os.environ[k_] = v
    import unc as UNC
    qbl = arb(repr(UNC.qbl_sup()))
    A__ = lambda v: arb(str(v))
    for key, al in (('1.25e-4', '0.000125'), ('1.5e-4', '0.00015')):
        d = json.load(open(RTs[al])); t = d['tail']; Ls = arb(repr(t['L_star'])); a_s = arb(repr(t['a_star']))
        Fn = A__(kc['FN'][0]) * a_s + A__(kc['FN'][1]) * a_s ** 2 + A__(kc['FN'][2]) * a_s ** 3
        row_ab = (A__(kc['CAB0']) + A__(kc['MAB']) * a_s) * Ls.sqrt() + A__(kc['QAB']) * qbl * Ls + A__(kc['PHIAB']) * Fn + A__(kc['STRIP_AB'])
        row_b = (A__(kc['BETA0']) + A__(kc['BETA1']) * a_s) * Ls.sqrt() + A__(kc['QB']) * qbl * Ls + A__(kc['BETAF']) * Fn + A__(kc['STRIP_BETA'])
        need = two / Ls
        rab = arb(ms[key]['sup_ab']) / row_ab - 1; rb = arb(ms[key]['sup_beta']) / row_b - 1
        C7['mod_sups_tail_row[%s]' % key] = dict(row_ab=FC.fup(row_ab), sup_ab=ms[key]['sup_ab'], reserve_ab=FC.fdn(rab), row_beta=FC.fup(row_b), sup_beta=ms[key]['sup_beta'],
                                                  reserve_beta=FC.fdn(rb), needed=FC.fup(need), ok=bool(rab.lower() >= need.upper() and rb.lower() >= need.upper()))
        cond('C7.mod[%s]' % key, C7['mod_sups_tail_row[%s]' % key]['ok'], FC.fdn(arb(min(float(rab.lower()), float(rb.lower()))) / need), 'mod_sups tail row: 6-digit reserve >= 2C\'/L*')
except Exception as ex:
    C7['mod_sups_tail_row'] = dict(error=str(ex)); cond('C7.mod', False, None, 'could not evaluate: %s' % ex)
# tail states: Lambda <= ell + 2C' <= L* + 2C' must lie in the last march cell (far_max, S_ramp_max, sqrtL_max are maxima over march steps, closure.py:153, 308-310)
for a_, p in sorted(RTs.items()):
    d = json.load(open(p)); last = d['steps'][-1]
    Ls_ = FC.from_raw_up(d['tail']['L_star']); mt = (arb(last['Ly_hi']) - Ls_) / two
    C7['tail_inside_last_march_cell[%s]' % a_] = dict(L_star=d['tail']['L_star'], Ly_hi_last=last['Ly_hi'], margin=FC.fdn(mt))
    cond('C7.tail[%s]' % a_, mt.lower() > 1, FC.fdn(mt), "L* + 2C' <= Ly_hi of the last march step (far sups, sqrtL_max cover the tail)")
C7['consumers_v36'] = [
    dict(consumer='v36cfg J_m, slopes, J4, sJ4 -> base_ibp/assembly_arb', reads='unc const/slope', treatment='up 5 digits; kappa_J pad', covered_by='C6 (constants gain C\'; slopes unchanged)'),
    dict(consumer='v36cfg sqLam (sqrtL_max)', reads='max_k sqrt(Ly_hi_k)', covered_by='C2 + C7.tail'),
    dict(consumer='v36cfg rho (tail)', reads='tail rho, down 4 digits', covered_by='unchanged: the tail decay rate of the barrier is the same (sec. 7.3 C4)'),
    dict(consumer='mod_sups sup_ab, sup_beta (abar, betabar) -> assembly_arb, t4lib, zone3, r3', reads='march rows at Ly_hi; tail row at L*', covered_by='C2 (march rows) + C7.mod (tail row reserve)'),
    dict(consumer='hjw_rep kappa_J', reads='Ly_hi, D, X_const, F0, sqrt_LF, L*', covered_by='C6 (pad 1e-9); W2 (L_F <= L*) unaffected'),
    dict(consumer='assembly_rep.margins (mod_sups M1, sI1), n1_chart (far_max)', reads='max over march steps', covered_by='C2 + C7.tail'),
    dict(consumer='assembly_arb.f4b PHI_S (raw sqrt_LF), PHI_C (raw sqrt L*)', reads='floor / anchor doubles, one ulp', covered_by='NOT cushioned; consumed only by I3c/I3s and v3.1-only rows (V36_REPORT sec. 5), which flux lemma v2 does not use: unused'),
]
floor_need = FC.aup(Cp / FC.adn(arb(min(per[a]['L_F_lower'] for a in per))))
C7['floor_relative_reserve_needed'] = dict(value=FC.fup(floor_need), note="consumers of sqrt_LF / alpha_F need relative reserve >= C'/L_F (sqrt: (1 + 2C'/L_F)^1/2 - 1 <= C'/L_F)")
C7['march_consumers'] = 'farsup, sramp_sup, box_assert, sqrtL_max, mod_sups march rows, unconditional r1+ Ly_hi read Ly_hi of the cell: covered by C2 [read closure.py]'
COND['C8'] = dict(pass_=None, note='regularity lemma U-1 (Lambda absolutely continuous, g\'\'\'\' in L^2(OR) a.e., level-2 identity with a C^3 cutoff): OPEN, not checkable by script')
COND['C9'] = dict(pass_=True, note='all window constants are at the (B Phi) level: exports doubled in X4 (2 (C0_k + C1_k T)); pointwise energies at beta_A = 2')
bad = [k for k, v in COND.items() if v['pass_'] is False]
R = dict(K0bar=FC.rec(K0, 'cert', 'level-0 pointwise bound of K0^OR'), K0_parts=dict(lin=FC.rec(K0lin, 'cert'), nl=FC.rec(K0nl, 'cert'), w=FC.rec(K0w, 'cert')),
         K2pt=FC.rec(K2pt, 'cert', 'level-2 pointwise part'), D_pt=FC.rec(Dpt, 'cert'), X_pt=FC.rec(Xpt, 'cert'), lam1='1', lam2='1',
         Cpp=FC.rec(Cpp, 'cert', 'int_I phi2 <= C\'\' + S\'\' T'), Spp=FC.rec(Spp, 'cert'), Cp=FC.rec(Cp, 'cert', "C' = mu eps_I C''"), Sp=FC.rec(Sp, 'cert'),
         two_Cp=FC.rec(two, 'cert'), F0=FC.rec(F0, 'cert'), slack_s=FC.rec_dn(s_, 'cert', 's = F0 - mu Kbar0 - mu eps_I K2pt (lower bound)'),
         margin_combined=FC.rec_dn(margin_comb, 'cert', '(FLUX0 + eps_I FLUX2)/(Kbar0 + eps_I K2pt)'), margin_K0_vs_FLUX0=FC.rec_dn(margin_K0, 'cert'),
         K2pt_vs_FLUX2_record=FC.fup(K2pt / FLUX2), strictness_pad_C3=per, C6=dict(kappa_J=hj['kappa_J'], kappa_J_certified=hj['kappa_J_certified'], pad=hj.get('kappa_J_pad'),
         X_const=FC.fdn(Xc), reserve=FC.fdn(reserve), Cp_over_Xconst=FC.fup(Cp / Xc)), C7=C7, stop_rules_triggered=STOP)
out = dict(stage='X5 flux2_assemble', res=R, conditions=COND, checks=CHK, verdict=('PASS' if not bad and not STOP else 'FAIL'), failing=bad,
           args=FC.args_block([os.path.abspath(__file__), os.path.join(FC.CODE, 'f2common.py'), os.path.join(FC.OUT, 'inputs.json'), os.path.join(FC.OUT, 'flux2_coef.json'),
                               os.path.join(FC.OUT, 'flux2_data.json'), INP['v36']['hjw_rep'], MS] + list(RTs.values())))
out['timing'] = dict(runtime_s=T_.s())
h = FC.dump(os.path.join(FC.OUT, 'flux2_cert.json'), out)
print('VERDICT', out['verdict'], 'failing', bad, 'stop', STOP, 'sha', h, '(%.1fs)' % T_.s())
sys.exit(0 if out['verdict'] == 'PASS' else 3)
