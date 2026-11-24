"""V4_REVIEW R-2 (strip hypothesis sup|w'''| <= 0.02 in STRIP_BETA / STRIP_AB): one integration by parts in the strip commutator pairing.
The modulation lemma bounds <[L,chi] w, f_j>_rho = int rho f_j sum_i c_i w^(i) (strip sigma in [29,32]), c = (chi'''' + k^2 chi'' + (5km + q/4) chi',
4 chi''' + 2 k^2 chi', 6 chi'', 4 chi')  (agent_mod/lemma/lem_integrals.py) by sum_i kappa_{j,i} S_i, S_i = sup_strip |w^(i)|, kappa_{j,i} = int rho |f_j| |c_i|.
The w''' term is integrated by parts once (chi'(29) = chi'(32) = 0, no boundary term):
    int rho f 4 chi' w''' = - int (4 rho f chi')' w''  =>  new w'' coefficient  rho [6 f chi'' - 4 (f' chi' + f chi'' - phi' f chi')] = rho [2 f chi'' - 4 f' chi' + 4 phi' f chi'],
so kappa'_{j,2} = int rho |2 f_j chi'' - 4 f_j' chi' + 4 phi' f_j chi'|, kappa'_{j,3} = 0, kappa'_{j,0}, kappa'_{j,1} unchanged.
f_1 = chi p, f_2 = chi psi_2 (frame), f_3 = chi F with F = ws0, F' = ws1 + phi' ws0 (outA5 rows; lyap_consts2.py:140); f' = chi' g + chi g'.
Then the U-3 strip form (v4cert/u3_strip_form.py: same feedback and linearisation) is recomputed with kappa' -> STRIP_BETA_R2, STRIP_AB_R2,
which need only S0, S1, S2 <= 0.02 (asserted along the march) and |alpha| <= 5e-4: no hypothesis on w'''.
Regression: the unchanged kappa_{j,0}, kappa_{j,1}, kappa_{j,2}, kappa_{j,3} recomputed here equal lem_integrals.json to rounding.  Output out/r2_strip_ibp.json."""
import sys, os, json, time
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from v4bcommon import OUT, args_block, dump, fup, V4B
from flint import arb, ctx
t0 = time.time()
sys.argv_saved = sys.argv; sys.argv = sys.argv[:1]
sys.path.insert(0, '/home/claude/s1spec/step3nl/agent_mod'); sys.path.insert(0, '/home/claude/s1spec/step3nl')
import mod_gram as MG
sys.argv = sys.argv_saved
ctx.prec = 100
A_, centered, dSm_pts, _pw, dp, dC, dS, weight, chunks = MG.A, MG.centered, MG.dSm_pts, MG._pw, MG.dp, MG.dC, MG.dS, MG.weight, MG.chunks
ab = lambda z: arb(abs(z).upper())
OUTA5 = '/home/claude/s1spec/step3/outA5.txt'
LEMI = os.path.join(V4B, 'agent_mod', 'lemma', 'lem_integrals.json'); LI = json.load(open(LEMI))
CtB, StB, Ct2h, St2h = MG.frame_data(); CH = chunks()
Z = lambda: arb(0)
kap_old = {j: [Z() for _ in range(4)] for j in (1, 2, 3)}; kap2_new = {j: Z() for j in (1, 2, 3)}; nb = 0
for line in open(OUTA5):
    if line[:2] != "MB": continue
    v = list(map(float, line.split()[1:])); ta, tb = v[0], v[1]
    if tb > 3.0 + 1e-12: continue                                   # strip only (t = 32 - sigma in [0, 3])
    h = [A_(v[2 + 2 * i], v[3 + 2 * i]) for i in range((len(v) - 2) // 2)]
    st, C_, S_ = h[0], h[1], h[2]; ws = h[5:9]
    c = min(int(ta + 1e-9), len(CH) - 1)
    kv, mv, pv_, dq = [A_(lo, hi) for lo, hi in CH[c][2]]
    sig = A_(32 - tb, 32 - ta); Xv = sig - st; qv = sig + dq; dX = arb(tb) - arb(ta)
    uu = A_(ta / 3, tb / 3); chi = [centered(dSm_pts, uu, j) * _pw(arb(-1) / 3, j) for j in range(5)]
    lo = max(float(Xv.lower()), 0.0); hi = float(Xv.upper()); Xv = arb((lo + hi) / 2, (hi - lo) / 2 * (1 + 1e-15) + 1e-300)
    P, rho = weight(Xv); W = ab(dX * rho)
    pvs = (kv, mv, qv, pv_, C_, S_)
    g = {1: [dp[0](*pvs), dp[1](*pvs)], 2: [CtB * dC[0](*pvs) + StB * dS[0](*pvs), CtB * dC[1](*pvs) + StB * dS[1](*pvs)], 3: [ws[0], ws[1] + P[0] * ws[0]]}
    k_ = ab(kv); ch = [ab(chi[i]) for i in range(5)]; k2 = k_ * k_
    cc = [ch[4] + k2 * ch[2] + (arb(5) * k_ * ab(mv) + ab(qv) / 4) * ch[1], arb(4) * ch[3] + 2 * k2 * ch[1], arb(6) * ch[2], arb(4) * ch[1]]
    for j in (1, 2, 3):
        f = chi[0] * g[j][0]; fp = chi[1] * g[j][0] + chi[0] * g[j][1]
        for i in range(4): kap_old[j][i] += W * ab(f) * cc[i]
        kap2_new[j] += W * ab(2 * f * chi[2] - 4 * fp * chi[1] + 4 * P[0] * f * chi[1])
    nb += 1
up = lambda z: float(z.upper())
reg = {j: dict(mine=[up(z) for z in kap_old[j]], lem_integrals=LI[f'kappa_{j}'], rel=[abs(up(z) / LI[f'kappa_{j}'][i] - 1) for i, z in enumerate(kap_old[j])]) for j in (1, 2, 3)}
assert all(r < 1e-9 for j in reg for r in reg[j]['rel']), reg
KAPN = {j: [arb(LI[f'kappa_{j}'][0]), arb(LI[f'kappa_{j}'][1]), arb(kap2_new[j].upper()), arb(0)] for j in (1, 2, 3)}
# ---- U-3 strip form with kappa' (same code path as v4cert/u3_strip_form.py; the commutator coefficients of mod_ibp.json are sum_j |H_ij| kappa_{j,r}, ibp_final.py:80)
sys.path.insert(0, os.path.join(V4B, 'agent_mod', 'ibp'))
import ibp_consts as K
mag = lambda z: arb(abs(z).upper())
IBP = os.path.join(V4B, 'agent_mod', 'ibp'); P1 = os.path.join(IBP, 'mod_ibp.json'); P2 = os.path.join(IBP, 'final_5e-4.json')
Mj = json.load(open(P1))['alpha_max_5e-4']; Fj = json.load(open(P2))
A = lambda v: arb(repr(v))
SB = arb('0.02'); ETA = A(Mj['eta']); assert Mj['eta'] == 0.02 and Mj['alpha1'] == 0.0005
comps = ['a', 'b', 'beta']
d = [A(v) for v in Mj['neumann']['d']]
gam = [A(Fj['final']['core'][c]['gamma_ref']) for c in comps]
chk = {}
cst = {}; lal = {}
for i, c in enumerate(comps):
    s = Mj['strip'][c]; gq = {k: A(v) for k, v in s['quad_gamma_L1'].items()}
    kap_json = [A(v) for v in s['commutator_S0_S3']]
    kap_re = [arb(sum((mag(K.H[i][j]) * arb(K.KAP[j][ip]) for j in range(3)), arb(0)).upper()) for ip in range(4)]
    chk[c] = [abs(float(kap_re[ip].upper()) / float(kap_json[ip].upper()) - 1) for ip in range(4)]
    kap = [arb(sum((mag(K.H[i][j]) * KAPN[j + 1][ip] for j in range(3)), arb(0)).upper()) for ip in range(4)]
    M1, M2, M3 = A(s['cubic_M1_L1']), A(s['cubic_M2_L1']), A(s['cubic_M3_L1'])
    cub01 = ETA * (2 * SB * M3 / 3 + SB * M2)
    cst[c] = [kap[0] + SB * (gq['00'] + gq['01'] + gq['02']) + cub01, kap[1] + SB * (gq['11'] + gq['12']) + cub01, kap[2] + SB * gq['22'] + 2 * ETA * SB * M1, kap[3]]
    lal[c] = A(s['linear_alpha_S0'])
assert all(r < 1e-9 for c in chk for r in chk[c]), chk
def mx(lst):
    m = lst[0]
    for v in lst[1:]: m = v if v.upper() > m.upper() else m
    return arb(m.upper())
fb_c = {}; fb_l = {}
for i, c in enumerate(comps):
    fb_c[c] = [cst[c][r] + gam[i] * mx([cst[k][r] / d[j] for j, k in enumerate(comps)]) for r in range(4)]
    fb_l[c] = lal[c] + gam[i] * mx([lal[k] / d[j] for j, k in enumerate(comps)])
ab_c = [(fb_c['a'][r] ** 2 + fb_c['b'][r] ** 2).sqrt() for r in range(4)]; ab_l = (fb_l['a'] ** 2 + fb_l['b'] ** 2).sqrt()
AL = arb('5e-4'); S = arb('0.02')
def const(cv, lv): return sum((cv[r] * S for r in range(3)), arb(0)) + AL * lv * S          # r = 0, 1, 2 only: the S3 coefficient is 0
assert all(fb_c[c][3] == 0 for c in comps)
U3 = json.load(open(os.path.join(V4B, 'v4cert', 'out', 'u3_strip_form.json')))
out = dict(description=__doc__, n_strip_boxes=nb, kappa_regression=reg, commutator_regression_rel=chk,
           kappa_new={j: [up(z) for z in KAPN[j]] for j in (1, 2, 3)},
           beta=dict(const_coef_S0_S3=[fup(v) for v in fb_c['beta']], alpha_coef_S0=fup(fb_l['beta'])),
           ab=dict(const_coef_S0_S3=[fup(v) for v in ab_c], alpha_coef_S0=fup(ab_l)),
           STRIP_BETA_R2=fup(const(fb_c['beta'], fb_l['beta'])), STRIP_AB_R2=fup(const(ab_c, ab_l)),
           v4_with_S3_hypothesis=dict(STRIP_BETA=U3['STRIP_BETA'], STRIP_AB=U3['STRIP_AB']),
           hypotheses=dict(S0_S1_S2_le=0.02, alpha_le=5e-4, S3='none'), runtime_s=round(time.time() - t0, 1),
           args=args_block([OUTA5, LEMI, P1, P2, os.path.join(V4B, 'v4cert', 'out', 'u3_strip_form.json'), '/home/claude/s1spec/step3nl/agent_mod/mod_gram.py']))
print(json.dumps({k: v for k, v in out.items() if k not in ('description', 'args')}, indent=1))
print('sha256', dump(os.path.join(OUT, 'r2_strip_ibp.json'), out))
