#!/usr/bin/env python3
"""Task 4, U1: certified (arb) weight constants of Theta*phi and Ghat_jk at the box, with phi = exp(K_b e^{-s}), K_b = 1e-3.

  PYTHONHASHSEED=0 nice -n 19 python3 -B u1_ghat.py [--ramp_U PATH]      -> out/u1_ghat.json

(1) Gamma_m, kappa~, kappa2~ of Theta (P1: weight_consts_arb.certify(4, 4, a in [0, 3/10], theta_a <= 3/5, y_L = 56), rerun here) and of Theta*phi
    by the exact Leibniz/triangle bounds [paper, evaluated in arb]: with beta := K_b e^{-s} in (0, K_b/28] on s >= ln 28,
      phi'/phi = -beta, phi''/phi = beta^2 + beta, phi'''/phi = -(beta^3 + 3 beta^2 + beta)        (beta' = -beta)
      Gamma_m(Theta phi) <= sum_k C(m,k) Gamma_k(Theta) phi_{m-k},  phi_0 = 1, phi_1 = b, phi_2 = b + b^2, phi_3 = b + 3b^2 + b^3, b = K_b/28
      L_phi = L - beta  =>  kappa~_phi <= kappa~ + b ;   L_phi^2 + L_phi' = L^2 + L' - 2 L beta + beta^2 + beta  =>  kappa2~_phi <= kappa2~ + 2 kappa~ b + b^2 + b.
(2) Ghat_jk: verbatim port of phase2a/U_ramp/ghat_design_arb.run (general-weight pipeline of phase1/ghat, Ab incl. U_j, dfrac 2/25, inner-ramp foot
    adds) with Theta*phi constants, PLUS (a) task 5's outer-foot Young adds to the global Y (t5_zone3/out/ramp_U.json, foot_out: level 1 Y_2,
    level 2 Y_2, Y_3) and (b) a separate sup|z_1| = M1z for the factors z_1 of the identity (pairs with a low index 1, the L^4 identity, the S33 IBP,
    the foot adds).  The A-coefficient box keeps |f_s| <= M1.  Reason for (b): z_1 = d_s(f - U_1) = y f_y - p U_1, so the box |y f_y| <= M1 alone gives
    |z_1| <= M1 + p sup U_1.  Primary M1z = 0.534 >= 0.52 + p U_1(60 lambda_1) + p c_dU 60^p on the core [28, 60] (asserted below); on zones I/II the
    assembly asserts S_1(2 Bhat) <= M1z (window Sobolev at the hypothesis level), and on zone III b_1 + p U_1(379.5) <= M1z (task 5's box).
    Also run: the contract convention (M1z = M1, no phi, no t5 adds) which must reproduce ghat_design_arb.json (asserted), and the fallback box.
(3) Design data56 with phi_max = e^{K_b/28} in (R2):  g := ceil(Ghat_22 phi/(1 - 4 abar - Ghat_22 phi/56^4)) (integer), then in arb
      (R2)  (1 - 4abar) g >= Ghat_22 (1 + g/56^4) phi_max        rate_margin := (1 - 4abar) 56^4/(Ghat_22 phi_max)
      (R<2) (1 - 4abar) g > Ghat_kk D(28) phi_max, k = 0, 1      R_lt2_margin_k := (1 - 4abar) g/(Ghat_kk D(28) phi_max)
      theta(28) = g/(28^4 + g) <= 3/5 (P1's validity range);  D(56), D(28).
abar = sup|a| from the verifier output (step3nl/verifier/final_alpha_points_5e-8.json, alpha_1 = 1.5e-4 row)."""
import sys, os, json, time, math, pickle, argparse
sys.dont_write_bytecode = True
T0 = time.time()
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import t4lib as L
from t4lib import A, Q_, arb, fup, flo, num
sys.path.insert(0, L.P1)
from flint import ctx
import ghat_arb as GA, weight_consts_arb as WC
ctx.prec = 128
ap = argparse.ArgumentParser(); ap.add_argument('--ramp_U', default=L.T5 + '/out/ramp_U.json'); cli = ap.parse_args()
GDA = L.ROOT + '/phase2a/U_ramp/ghat_design_arb.json'
VS = L.verifier_sups(); ABAR = VS['abar']
KB = A(Q_(1, 1000)); BMAX = KB/28; PHIMAX = (KB/28).exp()
T = pickle.load(open(GA.TABLES, 'rb'))
LN2 = arb(2).log(); etap = A(Q_(15, 8))/LN2; DEL = A(Q_(1, 10))
PQ = A(Q_(10246, 10000)); ALPHA1 = A(Q_(15, 100000))

# ---------------------------------------------------------------- (1) weight constants
r = WC.certify(4, 4, Q_(0), Q_(3, 10), Q_(3, 5), 56); b = r['bound']
assert not r['unresolved_cells']
Gam0 = {0: arb(1), 1: b['G1'], 2: b['G2'], 3: b['G3']}; kap0 = b['kap']; kap20 = b['kap2']
ph = {0: arb(1), 1: BMAX, 2: BMAX + BMAX**2, 3: BMAX + 3*BMAX**2 + BMAX**3}
GamP = {m: sum((math.comb(m, k)*Gam0[k]*ph[m - k] for k in range(m + 1)), arb(0)) for m in range(4)}
GamP[0] = arb(1)
kapP = kap0 + BMAX; kap2P = kap20 + 2*kap0*BMAX + BMAX**2 + BMAX
print('Theta    : Gamma', [fup(Gam0[m]) for m in (1, 2, 3)], 'kappa~', fup(kap0), 'kappa2~', fup(kap20))
print('Theta*phi: Gamma', [fup(GamP[m]) for m in (1, 2, 3)], 'kappa~', fup(kapP), 'kappa2~', fup(kap2P), flush=True)

# ---------------------------------------------------------------- (2) Ghat
t5 = json.load(open(cli.ramp_U)); FO = t5['foot_out']
YADD_T5 = {1: {2: arb(FO['1']['Y_add_out']['2']), 3: arb(FO['1']['Y_add_out']['3'])},
           2: {2: arb(FO['2']['Y_add_out']['2']), 3: arb(FO['2']['Y_add_out']['3'])}, 0: {}}

def run(M0q, M1q, Umq, M0zq, M1zq, Gam, kap, kap2, yadd, dfrac=Q_(2, 25)):
    """ghat_design_arb.run (foot=True, honest Ab), with M1z for the z_1 factors and extra Y adds; returns per level dict"""
    Asym = GA.build_At_sym(Umq); At = GA.eval_At(Asym, M0q, M1q); Ab = {k: d.get((), arb(0)) for k, d in At.items()}; out = []
    M1z = A(M1zq)
    for j in (0, 1, 2):
        YCb = GA.YC_base(T[j], j, Ab, Gam)
        add, dis, spare, ndfr, items = GA.l4_extra_YC_x(T, j, M0q, M1zq, At, Gam, kap, dfrac=dfrac, s33='ibp', M0zq=M0zq)
        fadd = {k: arb(0) for k in add}; extra = Q_(0)
        w33 = sum((it[1] for it in items if it[0] == 'S33'), arb(0))
        if any(it[0] == 'S33' for it in items):
            c1 = w33*M1z*etap; fadd[3] += c1*A(kap)*DEL; extra = Q_(1, 10)
        for it in items:
            if it[0] != 'L4': continue
            w = it[2]
            if it[1] == 3: fadd[3] += w*M1z*etap*DEL
            elif it[1] == 2: fadd[2] += w*M1z*etap
        ok, sl = GA.budget_ok(dis, ndfr, M0q, extra)
        YC = {k: YCb[k] + add.get(k, 0) + fadd.get(k, 0) + yadd.get(j, {}).get(k, arb(0)) for k in YCb}
        G = GA.ghat_x(YC, j, M0q, kap2)
        out.append(dict(G=G, YC=YC, ok=bool(ok), slack=sl, w33=w33))
    return out

def design(Gjj, label):
    e28 = arb(28)**4; e56 = arb(56)**4; G22 = Gjj[2]
    den = 1 - 4*ABAR - G22*PHIMAX/e56
    assert den > 0
    g0 = G22*PHIMAX/den; gq = Q_(math.ceil(fup(g0))); g = A(gq)
    R2 = (1 - 4*ABAR)*g - G22*(1 + g/e56)*PHIMAX
    D28 = 1 + g/e28; D56 = 1 + g/e56; th28 = g/(e28 + g)
    rlt2 = [(1 - 4*ABAR)*g/(Gjj[k]*D28*PHIMAX) for k in (0, 1)]
    d = dict(g=str(gq), g_float=float(gq), R2_holds=bool(R2 > 0), R2_slack=num(R2, label),
             rate_margin=num((1 - 4*ABAR)*e56/(G22*PHIMAX), label, 'lo', against='>= 10 (R2 rate margin)'),
             R_lt2_margin=[num(x, label, 'lo', against='>= 5') for x in rlt2], R_lt2_holds=[bool(x > 1) for x in rlt2],
             D56=num(D56, label), D28=num(D28, label), theta28=num(th28, label), theta28_le_3_5=bool(th28 < A(Q_(3, 5))),
             c_rs=num(G22*g/e28**2, label))
    return d, g

M1Z_CORE = Q_(267, 500)            # 0.534
UC60 = A(Q_(13056, 10000))*ALPHA1*(arb(60).log()*PQ).exp()
cdU_bound = A(Q_(2, 10**8))        # c_dU <= 2e-8 (certified in u34; asserted there)
lhs = A(Q_(52, 100)) + PQ*UC60 + PQ*cdU_bound*(arb(60).log()*PQ).exp()
assert lhs < A(M1Z_CORE), fup(lhs)
CASES = {
    'contract_noPhi_M1z=M1': dict(box=(Q_(15, 100), Q_(52, 100)), U=Q_(9, 100), M0z=Q_(4, 25), M1z=Q_(52, 100), phi=False, t5=False),
    'phi_only': dict(box=(Q_(15, 100), Q_(52, 100)), U=Q_(9, 100), M0z=Q_(4, 25), M1z=Q_(52, 100), phi=True, t5=False),
    'primary': dict(box=(Q_(15, 100), Q_(52, 100)), U=Q_(9, 100), M0z=Q_(4, 25), M1z=M1Z_CORE, phi=True, t5=True),
    'primary_M1z=M1': dict(box=(Q_(15, 100), Q_(52, 100)), U=Q_(9, 100), M0z=Q_(4, 25), M1z=Q_(52, 100), phi=True, t5=True),
    'fallback': dict(box=(Q_(20, 100), Q_(60, 100)), U=Q_(9, 100), M0z=Q_(1, 5), M1z=Q_(69, 100), phi=True, t5=True),
    'fallback_M1z=M1': dict(box=(Q_(20, 100), Q_(60, 100)), U=Q_(9, 100), M0z=Q_(1, 5), M1z=Q_(60, 100), phi=True, t5=True),
}
res = dict(weight=dict(Theta=dict(Gamma=[num(Gam0[m]) for m in (1, 2, 3)], kappa=num(kap0), kappa2=num(kap20), source='weight_consts_arb.certify(4,4,[0,3/10],3/5,56) rerun'),
                       Theta_phi=dict(Gamma=[num(GamP[m]) for m in (1, 2, 3)], kappa=num(kapP), kappa2=num(kap2P), rule='Leibniz/triangle with beta <= K_b/28 [paper], arb'),
                       rel_change=dict(Gamma=[fup(GamP[m]/Gam0[m] - 1) for m in (1, 2, 3)], kappa=fup(kapP/kap0 - 1), kappa2=fup(kap2P/kap20 - 1)),
                       K_b='1/1000', phi_max=num(PHIMAX), beta_max=num(BMAX)),
           abar=dict(value=VS['strings']['sup_ab'], src=VS['src']), M1z_core_check=dict(lhs=num(lhs), M1z=str(M1Z_CORE),
           rule='0.52 + p U_1(60) + p c_dU 60^p <= M1z (core [28, 60], B2prime |y f_y| <= 0.52; c_dU <= 2e-8 certified in u34)'),
           t5_Yadd=dict(src=cli.ramp_U, level1={str(k): fup(v) for k, v in YADD_T5[1].items()}, level2={str(k): fup(v) for k, v in YADD_T5[2].items()}),
           cases={})
GD = json.load(open(GDA))['cases']
for name, c in CASES.items():
    Gm, kp, kp2 = (GamP, kapP, kap2P) if c['phi'] else (Gam0, kap0, kap20)
    o = run(c['box'][0], c['box'][1], c['U'], c['M0z'], c['M1z'], Gm, kp, kp2, YADD_T5 if c['t5'] else {})
    Gjj = [o[j]['G'][j] for j in range(3)]
    dsg, g = design(Gjj, 'cert')
    row = dict(box=dict(M0=str(c['box'][0]), M1=str(c['box'][1]), Umax=str(c['U']), M0z=str(c['M0z']), M1z=str(c['M1z']), design='data56'),
               phi=c['phi'], t5_Yadd=c['t5'],
               Ghat=[[num(o[j]['G'][k]) for k in range(j + 1)] for j in range(3)], Ghat_jj=[num(x) for x in Gjj],
               Y_next=[num(o[j]['YC'][j + 1]) for j in range(3)], Y_full={str(j): {str(k): fup(v) for k, v in o[j]['YC'].items()} for j in range(3)},
               W33_level2=num(o[2]['w33']), budget_ok=[o[j]['ok'] for j in range(3)], budget_slack=[o[j]['slack'] for j in range(3)], design=dsg)
    res['cases'][name] = row
    print(f"{name:26s} Ghat_jj {[round(fup(x), 2) for x in Gjj]}  Y_next {[round(fup(o[j]['YC'][j + 1]), 3) for j in range(3)]}  budget {row['budget_ok']}  g {dsg['g']}"
          f"  rate x{dsg['rate_margin']['double']:.3f}  R<2 x{[round(x['double'], 2) for x in dsg['R_lt2_margin']]}  D56 {dsg['D56']['double']:.6f} D28 {dsg['D28']['double']:.5f}"
          f"  theta28 {dsg['theta28']['double']:.4f}", flush=True)
    assert all(row['budget_ok'])
    if name == 'contract_noPhi_M1z=M1':
        ref = GD['widened (0.15,0.52), U .09, M0z .16']['Ghat_jj']
        for j in range(3): assert abs(fup(Gjj[j])/ref[j][1] - 1) < 1e-9, (j, fup(Gjj[j]), ref[j])
        print('   contract convention reproduces ghat_design_arb.json Ghat_jj (rel < 1e-9)')
p, c = res['cases']['primary'], res['cases']['contract_noPhi_M1z=M1']
pm = res['cases']['primary_M1z=M1']
po = res['cases']['phi_only']
res['phi_effect_rel'] = dict(rule='phi_only vs contract (same box, M1z = M1, no t5 adds)', Ghat_jj=[po['Ghat_jj'][j]['double']/c['Ghat_jj'][j]['double'] - 1 for j in range(3)])
res['t5_add_effect_rel'] = dict(Ghat_jj=[pm['Ghat_jj'][j]['double']/po['Ghat_jj'][j]['double'] - 1 for j in range(3)])
res['M1z_effect_rel'] = dict(Ghat_jj=[p['Ghat_jj'][j]['double']/pm['Ghat_jj'][j]['double'] - 1 for j in range(3)])
res['args'] = L.args_block([GA.TABLES, GDA, cli.ramp_U, L.VERIF_ALPHA, os.path.join(L.P1, 'ghat_arb.py'), os.path.join(L.P1, 'weight_consts_arb.py'),
                            os.path.abspath(__file__), os.path.join(HERE, 't4lib.py')])
print('phi effect (rel):', res['phi_effect_rel']['Ghat_jj'], ' M1z effect (rel):', res['M1z_effect_rel']['Ghat_jj'])
L.dump(res, os.path.join(L.OUTD, 'u1_ghat.json'))
print('wrote out/u1_ghat.json, wall %.1f s' % (time.time() - T0))
