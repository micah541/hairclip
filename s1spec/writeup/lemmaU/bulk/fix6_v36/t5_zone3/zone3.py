#!/usr/bin/env python3
"""Task 5 (r2): zone III, Lemma R outer ramp, fixed-L cone family.  python-flint arb, outward rounding.  Writes out/zone3_interface.json.
Inputs: out/ramp_U.json (ramp_U.py), lemmaO_v32_consts.json (Lemma O, via z3lib), lemmaR_consts.json (Lemma R at L = 50),
        [optional, comparison only] t1_statements/r2_probe/ramp_out_253_379.5_box0.15_0.52.json, and a U_cert.json for Bhat (--ucert).
r_chi = 253.  Covering A: cone members L = 50, X = 3, centres [300, 380]; far members L = 113, X = 2, centres >= 376.
Formulas: z3lib.py (Z1/Z2), ramp_U.py (ramp coefficients), lemmaR.tex Sec. 4 (ramp norms), lemmaO_v32_h6.py (H6 structure); Z5_REPORT.md.
usage: PYTHONHASHSEED=0 nice -n 19 python3 -B zone3.py [--ucert PATH]"""
import sys, os, json, math, hashlib, platform, argparse
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from flint import arb, ctx
import z3lib as Z
Z.FAR_CMIN = 490                                                 # far family restricted to centres >= 490 (bases at r >= 264; Z5_REPORT.md)
ctx.prec = 128
ap = argparse.ArgumentParser(); ap.add_argument('--ucert', default=None); ap.add_argument('--ucert_full', default=None); args = ap.parse_args()
A, UP, LO, amax = Z.A, Z.UP, Z.LO, Z.amax
RAMPJ = os.path.join(HERE, 'out', 'ramp_U.json'); LR = '/home/claude/s1spec/writeup/lemmaO_v32/lemmaR_consts.json'
MECH = '/home/claude/s1spec/writeup/lemmaU/bulk/t1_statements/r2_probe/ramp_out_253_379.5_box0.15_0.52.json'
def sha(p): return hashlib.sha256(open(p, 'rb').read()).hexdigest()
def num(x, label='cert', **kw):
    d = dict(arb=x.str(20, radius=True) if isinstance(x, arb) else str(x), double=UP(x) if isinstance(x, arb) else x, label=label); d.update(kw); return d
def lnum(x, label='cert'): return num(x, label, rounding='down', double_lower=LO(x))
RCHI = arb(253); R2 = arb('379.5'); p = Z.p; U1 = Z.U1; ell15 = Z.ell
sys.path.insert(0, '/home/claude/s1spec/writeup/lemmaU/bulk/fix6_v36'); import v36cfg as CFG
_VS = json.load(open(CFG.MODS))['1.5e-4']   # M14 (M3): abar read from the verifier; v36: per V36_VER
abar = arb(str(_VS['sup_ab'])); bbar = arb(str(_VS['sup_ab'])); assert abar.upper() <= arb('2.19103e-7').upper()   # v36: <= the v35 value (fin13v4b: 2.14041e-7); q, cminus, Db are monotone in abar
q = 1/(1 - 4*abar); cminus = arb(1)/4 - abar
CONE, FAR, eC, eF = Z.CONE, Z.FAR, Z.eC, Z.eF
RJ = json.load(open(RAMPJ)); LRC = [c for c in json.load(open(LR))['cases'] if c['L_m'] == 50.0][0]
assert abs(LRC['t_needed']/1.6001e-7 - 1) < 1e-9 and LRC['n3'] in ('1/250', '0.004')
G_HI = arb(600000); KB = arb('0.001')
DPHI = (1 + G_HI/RCHI**4)*(KB/RCHI).exp()                       # sup D phi on the ramp (y >= 253); W <= 1
OUTJ = {}

# ======================================================================= Z1, Z2
print('Z1/Z2', flush=True)
z1, n1 = Z.scan(RCHI, R2, Z.fixedL(300))
b0, b1, bz = arb(repr(z1['b0']['value'])), arb(repr(z1['b1']['value'])), arb(repr(z1['bz']['value']))
dec = 'primary' if (z1['b0']['value'] <= 0.15 and z1['b1']['value'] <= 0.52 and z1['bz']['value'] <= 0.16) else \
      ('fallback' if (z1['b0']['value'] <= 0.20 and z1['b1']['value'] <= 0.60 and z1['bz']['value'] <= 0.20) else 'fail')
z2, n2 = Z.scan(R2, arb(700), Z.far_members, step=arb(1), chained_from=Z.RWEND.lower())
b0_700 = max(z1['b0']['value'], z2['b0']['value']); b1_700 = max(z1['b1']['value'], z2['b1']['value'])
print(' Z1', z1, dec, '\n Z2 beyond', z2, flush=True)

# ======================================================================= covering assertions (Z4/Z5)
half_c = 50*(1 - eC); half_f = 113*(1 - eF)
cov_cone = dict(lo=UP(300 - half_c), hi=LO(380 + half_c)); cov_far_lo = UP(Z.FAR_CMIN - half_f)
assert cov_cone['lo'] <= 253 and cov_cone['hi'] >= 379.5 and cov_far_lo <= 379.5
COVER = dict(cone_members='L = 50, X = 3, centres [300, 380]: conclusion |r - r_c| <= 50(1 - e_cone)', e_cone=UP(eC), e_far=UP(eF),
             cone_covers=[cov_cone['lo'], cov_cone['hi']], cone_covers_253_379_5=True,
             far_members='L = 113, X = 2, centres >= 490 (subset of the v3.2 family, centres >= 376)', far_covers_from=cov_far_lo, far_covers_379_5_700=True)

# ======================================================================= ramp norms N_k (Lemma R recipe, unit-50 members), pieces asserted
C1c, C2c, G3c = CONE['d1'], CONE['d2'], A(json.load(open(Z.CONSTS))['cases'][0]['deviation_upper_at_needed_time']['sqrtC3t'])
C1f, C2f = FAR['d1'], FAR['d2']; G3f = A(json.load(open(Z.CONSTS))['far_v32']['case']['deviation_upper_at_needed_time']['sqrtC3t'])
s0, m2c, n3c, M = arb('0.08'), arb('0.03'), arb('0.004'), arb('0.12'); m2f, n3f = arb('0.04'), arb('0.004')
R4 = A(LRC['R4_upper']); kap4 = arb(2).sqrt()*M + R4.sqrt()
Z1R = b1 + p*U1(R2)                                              # sup |z_1| on the ramp (z_1 = y f' - p U_1)
def piece_norms(a, b, L, fam, Fb, kap3extra=True, k3=True, b1cap=None, z0cap=None):
    if fam == 'cone': kap2 = m2c + C2c; kap3 = n3c + G3c; pm = (s0.atan() + C1c).tan()
    else: kap2 = m2f + C2f; kap3 = n3f + G3f; pm = (s0.atan() + C1f).tan()
    sp_, sf = 1 + pm*pm, 1 + Fb*Fb; ln = ((b - a)/L).sqrt(); el = (b/a).log(); Ub = U1(b)
    assert ((b - a)/L).upper() <= 2
    Ks = kap3 + 3*pm*kap2**2*ln
    D1 = sf**arb('1.5')*(b/L)*kap2
    if b1cap is not None and b1cap.upper() < D1.upper(): D1 = b1cap
    N2f = (b/L)**arb('1.5')*sf**2*(sp_**arb('0.25')*Ks + 3*Fb*kap2**2*ln)
    z0 = (Fb + Ub) if (z0cap is None or z0cap.upper() > (Fb + Ub).upper()) else z0cap
    n = [z0*el.sqrt(), (D1 + p*Ub)*el.sqrt(), N2f + (D1 + p*p*Ub)*el.sqrt()]
    if k3:
        N3f = (b/L)**arb('2.5')*sf**arb('2.5')*sp_**arb('0.25')*(kap4 + 10*(pm + Fb)*kap2*Ks + (6 + 15*(pm*pm + Fb*Fb))*kap2**3*ln)
        n.append(N3f + 3*N2f + (D1 + p**3*Ub)*el.sqrt())
    return n, dict(piece=[UP(a), UP(b)], L=UP(L), kap2=UP(kap2), kap3=UP(kap3), D1=UP(D1), N2f=UP(N2f), n=[UP(x) for x in n])
def ramp_N(cut):
    pcs = [(RCHI, cut, arb(300)), (cut, R2, arb(380))]
    N = [arb(0)]*4; rows = []
    for a, b, rc in pcs:
        assert a.lower() >= (rc - half_c).upper() and b.upper() <= (rc + half_c).lower(), ('piece outside conclusion region', UP(a), UP(b), UP(rc))
        n, row = piece_norms(a, b, arb(50), 'cone', b0, b1cap=b1, z0cap=bz); row['member_centre'] = UP(rc); row['inside_conclusion_region'] = True
        rows.append(row)
        for k in range(4): N[k] = N[k] + n[k]**2
    return N, rows

# ======================================================================= Delta^R assembly
K = {int(j): {int(k): A(v) for k, v in d.items()} for j, d in RJ['K'].items()}
KZ = {int(j): {tuple(eval(k)): A(v) for k, v in d.items()} for j, d in RJ['Kz'].items()}
SUPL = {0: bz, 1: Z1R}
Kf = {j: dict(K[j]) for j in K}; GEN = {j: {} for j in K}; FOLD = []
for j in (0, 1, 2):
    for key, w in KZ[j].items():
        idx = list(key); ww = w
        while len(idx) > 2 and any(x < 2 for x in idx):
            x = min([x for x in idx if x < 2], key=lambda x: float(SUPL[x].mid())); idx.remove(x); ww = ww*SUPL[x]
        if len(idx) > 2: GEN[j][tuple(sorted(idx))] = GEN[j].get(tuple(sorted(idx)), arb(0)) + ww; continue
        a_, b_ = sorted(idx)
        for kk, wk in (((a_, ww),) if a_ == b_ else ((a_, ww/2), (b_, ww/2))): Kf[j][kk] = Kf[j].get(kk, arb(0)) + wk
        FOLD.append(dict(level=j, cls=list(key), to=[a_, b_], weight=UP(ww)))
for j in GEN: assert set(GEN[j]) <= {(2, 2, 2), (2, 2, 3)}, GEN[j]
CXI = A(RJ['CXI']); FOOTW = {int(j): A(d['footw_out']) for j, d in RJ['foot_out'].items()}; W33 = {int(j): A(d['W33']) for j, d in RJ['foot_out'].items()}
pref = q*RCHI**-4*DPHI
def delta(N):
    N2, N3 = N[2].sqrt(), N[3].sqrt()
    sob = (N[2]/ell15 + 2*N2*N3).sqrt()                          # time average of sup_ramp |z_2| (Sobolev on s-length ln1.5, Jensen)
    C222 = N[2]*sob; C223 = N2*sob*N3
    D = []; parts = []
    for j in (0, 1, 2):
        pair = sum((Kf[j][k]*N[k] for k in Kf[j]), arb(0))
        cub = GEN[j].get((2, 2, 2), arb(0))*C222 + GEN[j].get((2, 2, 3), arb(0))*C223
        s33 = W33[j]*CXI/3*C222; foot = FOOTW[j]*N[2]
        D.append(pref*(pair + cub + s33 + foot))
        parts.append(dict(pairs={str(k): UP(pref*Kf[j][k]*N[k]) for k in sorted(Kf[j])}, cubic=UP(pref*cub), S33_foot=UP(pref*s33), L4_foot=UP(pref*foot)))
    Db = [8*bbar*Z.A(RJ['ramp']['chi_prime_max'])/(RCHI*(1 - 4*abar))*DPHI*N[j] for j in range(3)]
    return D, Db, parts
best = None
for cs in ('331.5', '333', '335', '337.5', '340', '345', '348.5'):
    N, rows = ramp_N(arb(cs)); D, Db, parts = delta(N)
    if best is None or D[2].upper() < best[2][2].upper(): best = (cs, N, D, Db, rows, parts)
CUT, NR, DR, DRb, ROWS, PARTS = best
print(' ramp cut', CUT, 'N^2', [UP(x) for x in NR], 'DeltaR', [UP(x) for x in DR], 'DeltaR_b', [UP(x) for x in DRb], flush=True)
mech = None
if os.path.exists(MECH):                                         # contract-mechanical number (omega = 28/y weight, conversion 379.5/28): comparison only
    co = json.load(open(MECH))['c']
    mech = [UP(q*arb('379.5')/28*sum((A(co[str(j)][str(k)])*NR[k] for k in range(j + 2)), arb(0))) for j in range(3)]

# ======================================================================= Z3: n_III_k
far_pcs = [(R2, arb('600.5')), (arb('600.5'), arb(700))]
pieces3 = []
for a, b, rc in ((RCHI, arb(CUT), 300), (arb(CUT), R2, 380)):
    n, row = piece_norms(a, b, arb(50), 'cone', b0, k3=False, b1cap=b1, z0cap=bz); pieces3.append((a, b, n, row))
k0 = arb('0.04')/113
for a, b in far_pcs:
    c = (a + b)/2; assert c.lower() >= Z.FAR_CMIN and ((b - a)/2).upper() <= half_f.lower()
    psi = Z.WIN[-1]['F0'].atan() + k0*(1 + arb('0.25')).sqrt()*(b - arb('379.5')); Fb = (psi + C1f).tan(); assert UP(Fb) < 0.5
    n, row = piece_norms(a, b, arb(113), 'far', Fb, k3=False); row['member_centre'] = UP(c); pieces3.append((a, b, n, row))
_CDU = json.load(open(args.ucert_full))['launch']['c_dU']; cdU = arb(_CDU['arb']); assert _CDU['label'] == 'cert' and cdU.upper() < arb('2e-8')   # M14 (M3): certified c_dU of U_cert (u34; was typed 2e-8 [hyp])
S = [arb(0)]*3; mu_lo = arb(253)/700; NMU = 400
for i in range(NMU):
    m_a = mu_lo + (1 - mu_lo)*i/NMU; m_b = mu_lo + (1 - mu_lo)*(i + 1)/NMU
    lo, hi = 350*m_a, 700*m_b
    tot = [arb(0)]*3
    for a, b, n, row in pieces3:
        if b.upper() >= max(float(lo.lower()), 253.0) and a.lower() <= hi.upper():
            for k in range(3): tot[k] = tot[k] + n[k]**2
    for k in range(3): S[k] = S[k] if S[k].upper() >= tot[k].upper() else arb(UP(tot[k]))
nIII = [(S[k].sqrt() + p**k*cdU*arb(700)**p*arb(2).log().sqrt())**2 for k in range(3)]
tauIII = (arb(700)/253).log()/cminus
print(' n_III', [UP(x) for x in nIII], 'tau_III', UP(tauIII), flush=True)

# ======================================================================= Z5: H6 for the fixed-L family and the far family; U-fine constants
wmin = 1/(1 + 2*arb('-1.2').exp() * 0 + arb('-1.2').exp() + arb('-1.2').exp()*(arb(2)/3)**4)     # a_W = 0.3: Q convex, max at sigma = 0, ln1.5
WINF = [(arb(150), arb(150)*3/2, arb(150), arb(253)/1.5), (arb(253)/1.5, arb(253), arb(253)/1.5, arb(253))]   # (rho, 1.5 rho, sub-range lo, hi)
assert WINF[1][1].upper() <= 253 + 1e-12
def ufine(Bh):
    Eu = [arb(x)/wmin for x in Bh]
    S0 = (Eu[0]/ell15 + 2*(Eu[0]*Eu[1]).sqrt()).sqrt(); S1 = (Eu[1]/ell15 + 2*(Eu[1]*Eu[2]).sqrt()).sqrt()
    ws = []
    for rho, yb, y0, y1 in WINF: ws.append(dict(rho=rho, yb=yb, y0=y0, y1=y1, F0=U1(yb) + S0, D1=p*U1(yb) + S1, Eu=Eu, Ub=U1(yb)))
    return ws, S0, S1
sq = lambda s: 1 + s*s
tan_add = lambda a, d: (a.atan() + d).tan()
def H6(Bh):
    ws, S0, S1 = ufine(Bh); F0s = amax([w['F0'] for w in ws])
    Lc, Lf = arb(50), arb(113)
    s_cov = tan_add(s0, amax([C1c, C1f]))
    # U-fine points p in [150, 253): psi(p, t), psi(p, t_1) in atan(U_1(p) +- S_0) (Lemma U bound holds at t_1 too), so the angle change is <= 2 S_0;
    # data slope in frame c at t_1 is <= s0 (H5, whole cylinder):  s_str = tan(atan s0 + 2 S_0)
    s_str = tan_add(s0, 2*S0)
    sm_c = tan_add(s0, C1c); sm_f = tan_add(s0, C1f)
    # ---- cone (fixed-L) cylinders
    k_cov = sq(s0)**arb('1.5')*m2c + amax([C2c, (Lc/Lf)*C2f]); K2_cov = sq(s_cov)**arb('1.5')*k_cov
    k_str = amax([Lc*w['D1']/w['y0'] for w in ws]); K2_str = sq(s_str)**arb('1.5')*k_str
    def strip_part(Lx, s_x):
        acc = arb(0)
        for w in ws:
            pw = ((w['y0']**-3 - w['y1']**-3)/3).sqrt()
            nw = w['y0']**arb('-1.5')*(w['Eu'][2].sqrt() + w['Eu'][1].sqrt()) + ((p*p - p)*w['Ub'] + 3*(w['F0'] + s_x)*w['D1']**2)*pw
            acc += nw**2
        return sq(s_x)**2*sq(F0s)**arb('0.25')*Lx**arb('1.5')*acc.sqrt()
    def piece(Lx, s_x, kx, Lm, n3m, g3m, sm, length):
        return sq(s_x)**2*(sq(sm)**arb('0.25')*(Lx/Lm)**arb('1.5')*(n3m + g3m) + 3*(s_x + sm)*kx**2*(arb(length)/Lx).sqrt())
    for (a, b, rc) in ((253, 350, arb('301.5')), (350, LO(380 + half_c), arb(380))):
        assert a >= UP(rc - half_c) and b <= LO(rc + half_c)
    partsC = dict(ufine=strip_part(Lc, s_str), fixedL_253_350=piece(Lc, s_cov, k_cov, Lc, n3c, G3c, sm_c, 97),
                  fixedL_350_428=piece(Lc, s_cov, k_cov, Lc, n3c, G3c, sm_c, 79), far_428_530=piece(Lc, s_cov, k_cov, Lf, n3f, G3f, sm_f, 102))
    K3c = sum((v*v for v in partsC.values()), arb(0)).sqrt(); K2c = amax([K2_cov, K2_str])
    # ---- far cylinders (centre >= 376, base [r_c - 226, r_c + 226]; inner points [150, 253) from U-fine)
    sf_cov = tan_add(s0, amax([C1c, C1f]))
    sf_str = tan_add(s0, 2*S0)                                   # only used if a far base meets [150, 253): not the case for centres >= 490 (asserted)
    kf_cov = sq(s0)**arb('1.5')*m2f + amax([C2f, (Lf/Lc)*C2c]); K2f_cov = sq(sf_cov)**arb('1.5')*kf_cov
    K2f_str = arb(0)
    assert Z.FAR_CMIN - 226 >= 253 and Z.FAR_CMIN - 226 >= UP(300 - half_c)    # far bases start in the fixed-L coverage
    partsF = dict(fixedL_253_350=piece(Lf, sf_cov, kf_cov, Lc, n3c, G3c, sm_c, 97),
                  fixedL_350_428=piece(Lf, sf_cov, kf_cov, Lc, n3c, G3c, sm_c, 79), far1=piece(Lf, sf_cov, kf_cov, Lf, n3f, G3f, sm_f, 226),
                  far2=piece(Lf, sf_cov, kf_cov, Lf, n3f, G3f, sm_f, 226))
    K3f = sum((v*v for v in partsF.values()), arb(0)).sqrt(); K2f = amax([K2f_cov, K2f_str])
    KC, KF = arb('0.18'), arb('0.3')
    out = dict(cone=dict(slope=dict(covered=UP(s_cov), ufine=UP(s_str), margin=LO(arb('0.25')/amax([s_cov, s_str])), against=0.25),
                         K2=dict(covered=UP(K2_cov), ufine=UP(K2_str), bound=UP(K2c), margin=LO(KC/K2c), against=0.18),
                         K3=dict(parts={k: UP(v) for k, v in partsC.items()}, bound=UP(K3c), margin=LO(KC/K3c), against=0.18)),
               far=dict(slope=dict(covered=UP(sf_cov), ufine=None, margin=LO(arb('0.25')/sf_cov), against=0.25),
                        K2=dict(covered=UP(K2f_cov), ufine=None, bound=UP(K2f), margin=LO(KF/K2f), against=0.3),
                        K3=dict(parts={k: UP(v) for k, v in partsF.items()}, bound=UP(K3f), margin=LO(KF/K3f), against=0.3)),
               ufine_windows=[dict(rho=UP(w['rho']), yb=UP(w['yb']), F0=UP(w['F0']), D1=UP(w['D1'])) for w in ws], S0=UP(S0), S1=UP(S1),
               K2_ufine=UP(amax([Lc*w['D1']/w['y0'] for w in ws])), K3_ufine_strip=UP(strip_part(Lc, s_str)))
    out['min_margin'] = min(out['cone']['slope']['margin'], out['cone']['K2']['margin'], out['cone']['K3']['margin'], out['far']['slope']['margin'], out['far']['K2']['margin'], out['far']['K3']['margin'])
    return out
BH_PROV = ('5e-6', '1e-3', '0.53'); bh_src = 'provisional (CONTRACTS.md)'
if args.ucert:
    uc = json.load(open(args.ucert)); BH_PROV = tuple(repr(x) for x in uc['supersolution']['Bhat']); bh_src = args.ucert
h6_ref = H6(BH_PROV)
BH_FINE_EST = tuple(repr(float(UP(wmin*arb(x)))) for x in ('1.6e-5', '8.22e-4', '0.0446'))      # gonogo T5 CLEAN-GO Eu(150) [est], converted to Bhat
h6_est = H6(BH_FINE_EST)
# break-even: largest Bhat_2 (levels 0, 1 at the reference) with every H6 margin >= 1.5
lo_, hi_ = 0.0, 0.53
for _ in range(40):
    mid = (lo_ + hi_)/2
    if H6((BH_PROV[0], BH_PROV[1], repr(mid)))['min_margin'] >= 1.5: lo_ = mid
    else: hi_ = mid
print(' H6 at Bhat_ref', BH_PROV, h6_ref['min_margin'], ' at fine est', h6_est['min_margin'], ' break-even Bhat_2 (margin 1.5):', lo_, flush=True)

# ======================================================================= Z6: hand-off thresholds
HO = [
    dict(name='(H2f)', region='y in [60, 700] at t_1, G_j = int (28/y)(d_s^j (f - U_1))^2 ds', value=[5e-7, 5e-6, 1e-5], new=True,
         step4_margin_if_known=[1.5, 4.9, 4.0], margin_label='probe (base curve only, r2_probe/far_summary.json); treat as est', used_by_task5=False),
    dict(name='(H3) windows', region='[y_w, 1.5 y_w], y_w = 50 ... 253.125 (to 379.69)', value=[2e-7, 1e-5, 1e-5], new=False,
         step4_margin_if_known=[5.7, 78.7, 26.6], used_by_task5='Z1 data (|f|, |v|, y|f\'| at t_1), Z3, Z4'),
    dict(name='alpha_1', region='hand-off', value=1.5e-4, new=False, step4_margin_if_known=1.27, used_by_task5='U_1 in Z1/Z3/Z4'),
    dict(name='(H5) cone data, fixed-L cylinders', region='every cylinder [r_c - 150, r_c + 150], r_c in [300, 380], L = 50 units, X = 3',
         value=dict(slope_s0=0.08, h2_m2=0.03, L2_h3_n3=0.004), new=True, step4_margin_if_known=None,
         note='same values as v3.2; NEW GEOMETRY (v3.2: units r_m/4, centres [200, 380]); step 4 must rerun t9_cone_kch for it'),
    dict(name='Lemma O H1-H4, H7 for fixed-L cylinders', region='same cylinders', value=dict(H1='sup|h0| <= H - 0.00227 (C0 sqrt t, L=50)', H3='|h0^(j)| <= 0.12, j = 2,3,4',
         H4='int (h0^(5))^2 <= 0.0144', H7='ratio 1.6001e-7 (= 50^-4 q), t* margin 2.30'), new=True, step4_margin_if_known=None,
         note='new geometry; rerun h4loc (cone cylinders) and the Lemma O H1 graph check'),
    dict(name='(H5) far data', region='far cylinders, L = 113', value=dict(slope=0.08, h2=0.04, L2_h3=0.004), new=False, step4_margin_if_known=[1.87, 1.93, 5.6]),
]
RERUNS = [
    dict(stage='t9_cone_kch', input='fixed-L geometry: L = 50, X = 3, centres [300, 380] (cell cover), cylinders [r_c - 150, r_c + 150]', check='slope <= 0.08, |h0\'\'| <= 0.03, ||h0\'\'\'||_L2 <= 0.004, jets <= 0.12 (L = 50 units)'),
    dict(stage='h4loc (cone cylinders)', input='same geometry', check='int (h0^(5))^2 <= M^2 = 0.0144'),
    dict(stage='Lemma O H1 (graph, C0) for cone cylinders', input='same geometry', check='graph over [-3, 3] with sup|h0| <= H - 0.00227'),
    dict(stage='t9_h6_v32', input='superseded by this task\'s H6 (zone3_interface.json cone_family.H6); keep the far-family data rows', check='—'),
    dict(stage='(H2f) window stage, y in [60, 700]', input='geometry to x = 730; needs Micah\'s approval (R2_DECISIONS.md)', check='(5e-7, 5e-6, 1e-5)'),
]

# ======================================================================= JSON
ACC = dict(decision=dec, box_margin_min=min(0.15/z1['b0']['value'], 0.52/z1['b1']['value'], 0.16/z1['bz']['value']),
           HO1_ok=bool(b0_700 <= 0.5 and b1_700 <= 1.0), DeltaR2=UP(DR[2]), DeltaR2_ok=bool(UP(DR[2]) <= 0.03),
           tstar_margin=json.load(open(Z.CONSTS))['cases'][0]['margin_factor_lower'], H6_min_margin_at_Bhat_ref=h6_ref['min_margin'],
           H6_min_margin_at_fine_est=h6_est['min_margin'], covering_asserts=True)
out = dict(
    r_chi=253, tau_III=num(tauIII), provisional=bool(args.ucert is None),
    box_III_380=dict(region=[253, 379.5], b0=num(b0, at_r=z1['b0']['at_r']), b1=num(b1, at_r=z1['b1']['at_r']), bz=num(bz, at_r=z1['bz']['at_r']), cells=n1,
                     margins=dict(b0=0.15/z1['b0']['value'], b1=0.52/z1['b1']['value'], bz=0.16/z1['bz']['value'], against=[0.15, 0.52, 0.16])),
    decision=dec,
    box_III_700=dict(region=[253, 700], b0=num(arb(repr(b0_700))), b1=num(arb(repr(b1_700))), beyond_379_5=z2, HO1_ok=ACC['HO1_ok'], against=[0.5, 1.0]),
    n_III=[num(x) for x in nIII], n_III_note='sup_tau int over {350 <= y <= 700, r >= 253 lambda_1} of (z^A_k)^2 ds; nonempty only for tau - tau_1 <= tau_III; '
          'z^A = z + d_s^k deltaU with c_deltaU <= 2e-8 [hyp, task 4 certifies]; pieces below', n_III_pieces=[r[3] for r in pieces3],
    DeltaR=[num(x) for x in DR], DeltaR_b=[num(x) for x in DRb],
    DeltaR_detail=dict(pieces=ROWS, cut=float(CUT), N_sq=[num(x) for x in NR], N_sq_meaning='k <= 2: sup_t int_ramp z_k^2 ds; k = 3: time average over the horizon',
                       pref=num(pref), pref_meaning='q r_chi^-4 sup(D phi), g <= 6e5', parts=PARTS, folded_classes=FOLD,
                       genuine_cubic={str(j): {str(k): UP(v) for k, v in d.items()} for j, d in GEN.items()},
                       DeltaR_b_coefficient='8 bbar (1.875/ln1.5)/(r_chi (1 - 4 abar)) sup(D phi) N_j^2', mechanical_ramp_out_omega_weight=mech,
                       mechanical_note='contract-mechanical rerun (omega = 28/y weight x 379.5/28, Umax 0.05, jet-free): comparison only, NOT a bound for the Lemma-U weight'),
    cone_family=dict(L=50, X=3, centres=[300, 380], covering=COVER, tstar_margin=ACC['tstar_margin'], R4=num(R4), R4_label='cert (lemmaR_consts.json, L = 50)',
                     H6=dict(at_Bhat_ref=h6_ref, Bhat_ref=list(BH_PROV), Bhat_source=bh_src, at_fine_est=h6_est, Bhat_fine_est=list(BH_FINE_EST),
                             Bhat2_break_even_for_margin_1_5=lo_),
                     slope_check='closure slope margins in H6; the step-4 data check (slope <= 0.08 on fixed-L cylinders) is task 6'),
    ufine=dict(rule='windows [rho, 1.5 rho], rho = 150 and 253/1.5 (none enters the ramp); Eu_j = Bhat_j/w_min; S_0, S_1 window Sobolev (s-length ln1.5); '
                    'K2 = max_w (50/rho_w)(p U_1(1.5 rho_w) + S_1); K3 strip part as lemmaO_v32_h6.py strip formula with L_c = 50; |f| <= U_1(1.5 rho) + S_0',
               w_min=num(wmin), K2_at_Bhat_ref=h6_ref['K2_ufine'], K3_strip_at_Bhat_ref=h6_ref['K3_ufine_strip'], K2_at_fine_est=h6_est['K2_ufine'],
               K3_strip_at_fine_est=h6_est['K3_ufine_strip'], Bhat_ref=list(BH_PROV)),
    handoff_thresholds=HO, step4_reruns=RERUNS,
    inputs_used=dict(box_III_380='B4 (fixed-L cone members, far members), D1 (abar for the horizon), (H3) thresholds, alpha_1 <= 1.5e-4, (H5) s0 = 0.08, m2 = 0.03 (cone), 0.04 (far)',
                     box_III_700='as box_III_380, plus (H5) far m2 = 0.04 chained beyond 379.69',
                     n_III='B4, D1, (H3), (H5), Lemma O C_j, c_deltaU [hyp]; not (H2f)',
                     DeltaR='B4, D1 (abar, bbar), (H3), (H5) incl. n3 = 0.004, Lemma R (R4, L = 50), zone-III box; global Y_{j+1} of the widened design; g <= 6e5',
                     cone_family='B4, D1, (H5), Lemma O/R constants; H6 also U-fine Bhat (task 4)', never='B2prime, B3U'),
    acceptance=ACC,
    args=dict(argv=sys.argv, env={k: os.environ.get(k) for k in ('PYTHONHASHSEED',) if os.environ.get(k)}, python=platform.python_version(), prec=ctx.prec,
              inputs={pth: sha(pth) for pth in [RAMPJ, LR, Z.CONSTS, os.path.join(HERE, 'z3lib.py'), os.path.abspath(__file__)] + ([MECH] if os.path.exists(MECH) else []) + ([args.ucert] if args.ucert else []) + [args.ucert_full]}))
json.dump(out, open(os.path.join(HERE, 'out', 'zone3_interface.json'), 'w'), indent=1, sort_keys=True)
print(json.dumps(ACC, indent=1))
