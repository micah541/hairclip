#!/usr/bin/env python3
"""Task 4, U2: certified (arb) inner-ramp constants of Lemma U on y in [28, 56] for the weight Theta*phi (P2 rerun with phi).

  PYTHONHASHSEED=0 nice -n 19 python3 -B u2_K.py [--u1 out/u1_ghat.json] [--case primary] [--maxeval N] [--tol x] [--fresh]
      -> out/u2_K.json  (resumable: out/u2_K.partial.json)

Port of phase1/ramp/ramp_K_cert.py (library phase1/ramp/ulib.py, imported read-only), with these changes:
  * weight Theta*phi, phi = exp(K_b e^{-s}), K_b = 1e-3: with beta = K_b e^{-s} = (K_b/28) 2^{-t} (exact per cell) and RT_k = Theta^(k)/Theta (P2 forms),
      RTp_1 = RT_1 - beta,  RTp_2 = RT_2 - 2 beta RT_1 + (beta + beta^2),  RTp_3 = RT_3 - 3 beta RT_2 + 3 (beta + beta^2) RT_1 - (beta + 3 beta^2 + beta^3)
    (Leibniz with phi'/phi = -beta, phi''/phi = beta + beta^2, phi'''/phi = -(beta + 3 beta^2 + beta^3)); L_phi = -4 + RTp_1.
      rem_m = sum_{l=1..m} C(m,l) (eta^2)^(l) RTp_{m-l};  rho^ = (1/2)|(eta^2)'' + 2 L_phi (eta^2)'|;  ceta = L_phi eta eta' + eta'^2 + eta eta''.
  * a_W = 3/10 fixed (the weight of STATEMENTS_r2 Sec. 3); g in [G_LO, G_HI] = [3e5, 6.4e5] as a branch-and-bound dimension (covers the primary g and
    the contract fallback g = 6.29e5); centres s_rho >= ln 56 (all centres, at every time).
  * K_jj rho-part with the GLOBAL Y_{j+1} of the chosen case of u1 (box (0.15, 0.52), U .09, M0z .16, M1z .534, phi, task-5 outer-foot adds).
  * A-coefficient box on the ramp: B2prime's core box (0.11, 0.52), U_ramp = 1/80 >= 1.3056 alpha_1 56^p (as P2; U_1 frozen at tau_1, lambda <= lambda_1).
  * Folding sups on the ramp: |z_0| <= 0.11 (B2prime: |v_A| <= 0.0153, |dU| <= c_dU 56^p), |z_1| <= M1z = 0.534 (P2 used 0.52; z_1 = y f_y - p U_1).
  * Foot / S33 weights (port of the P2 foot block) with the global box of the case, Gamma_m, kappa~ of Theta*phi, M1z for the z_1 factors.
  * Genuinely cubic classes kept apart: (2,2,2) -> c33 (S33 data form, no IBP); (2,2,3) -> c223 := sup|c| per eps*Theta*phi for the Cauchy-Schwarz
    route int int |c| eps Theta z_2^2 |z_3| <= c223 (int int eps Theta z_2^4)^{1/2} (int int eps Theta z_3^2)^{1/2} (no integration by parts, hence no
    point term from the jump of (eta^2)''' at y = 56; the m-orders of the class are recorded).  Lemma C's IBP route is not used.
Every number is an arb upper endpoint over a branch-and-bound cover (ulib.bb_sup); lower bounds (bb block) show the sups are attained within tol."""
import sys, os, json, math, time, argparse, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import t4lib as L
sys.path.insert(0, L.P2); sys.path.insert(0, L.P1)
import sympy as sp
from flint import arb, ctx
import ulib as U
import ghat_arb as GA
from t4lib import A, Q_, fup, flo, num
from ulib import hull, LN2, LN28, LN56
ctx.prec = 128
ap = argparse.ArgumentParser()
ap.add_argument('--u1', default=os.path.join(L.OUTD, 'u1_ghat.json')); ap.add_argument('--case', default='primary')
ap.add_argument('--maxeval', type=int, default=150000); ap.add_argument('--tol', type=float, default=0.003); ap.add_argument('--fresh', action='store_true')
ap.add_argument('--glo', default='300000'); ap.add_argument('--ghi', default='640000'); ap.add_argument('--tag', default=''); ap.add_argument('--nophi', action='store_true')
args = ap.parse_args()
OUT = os.path.join(L.OUTD, f'u2_K{args.tag}.json'); PART = os.path.join(L.OUTD, f'u2_K{args.tag}.partial.json')
G_LO, G_HI = Q_(args.glo), Q_(args.ghi); A_W = arb(3)/10   # M14 (M3): exact 3/10 as an arb ball (was the double 0.3)
KB = arb(0) if args.nophi else A(Q_(1, 1000))
R = {} if args.fresh or not os.path.exists(PART) else json.load(open(PART))
def save(): L.dump(R, PART)
def cached(key, fn):
    if key in R: return R[key]
    t0 = time.time(); v = fn(); R[key] = v; save()
    print(f'  {key}: ub {v["ub"]:.7g} (lb {v["lb"]:.7g}, {v["nev"]} cells, {time.time() - t0:.1f} s)', flush=True); return v
def bb(f, boxes):
    ub, lb, nev, left = U.bb_sup(f, boxes, tol=args.tol, maxeval=args.maxeval)
    return dict(ub=ub, lb=lb, nev=nev, left=left)

T = U.load_tables()
U1 = json.load(open(args.u1)); CASE = U1['cases'][args.case]
YNEXT = [arb(CASE['Y_next'][j]['double']) for j in range(3)]
BOXG = CASE['box']; M0G, M1G, UG, M0ZG, M1ZG = (Q_(BOXG[k]) for k in ('M0', 'M1', 'Umax', 'M0z', 'M1z'))
WPH = U1['weight']['Theta_phi']; GamP = {0: arb(1)}
for m in (1, 2, 3): GamP[m] = arb(WPH['Gamma'][m - 1]['double'])
KAPP = arb(WPH['kappa']['double'])

# ================================================================== ramp geometry with phi; box (t, w, a, g)
WT = 16.0
def ramp_q(box, point=False):
    tb = box[0]
    t = arb((tb[0] + tb[1])/2) if point else hull(tb[0], tb[1])
    E2 = [U.e2(t, l) for l in range(5)]
    (wlo, whi), (alo, ahi), (gl, gh) = box[1], box[2], box[3]
    a = arb((alo + ahi)/2) if point else hull(alo, ahi)
    if whi == math.inf:
        p, q = U.pq_lower_tail((arb(tb[1]) - 1)*LN2 - wlo if not point else (t - 1)*LN2 - wlo, a)
    else:
        w = arb((wlo + whi)/2) if point else hull(wlo, whi)
        p, q = U.pq_from_sigma((t - 1)*LN2 - w, a)
    Rk = U.Rk(p, q)
    gb = arb((gl + gh)/2) if point else hull(gl, gh)
    s = LN28 + t*LN2; th = U.theta_s(s, gb)
    RT = {0: arb(1), 1: Rk[1] - 4*th, 2: Rk[2] + th*(16 - 8*Rk[1]), 3: Rk[3] + th*(-12*Rk[2] + 48*Rk[1] - 64)}
    be = (KB/28)*(-t*LN2).exp()
    f1, f2, f3 = -be, be + be*be, -(be + 3*be*be + be*be*be)
    RP = {0: arb(1), 1: RT[1] + f1, 2: RT[2] + 2*RT[1]*f1 + f2, 3: RT[3] + 3*RT[2]*f1 + 3*RT[1]*f2 + f3}
    rem = {m: sum((math.comb(m, l)*E2[l]*RP[m - l] for l in range(1, m + 1)), arb(0)) for m in (1, 2, 3)}
    Lp = -4 + RP[1]
    rho = (E2[2] + 2*Lp*E2[1])/2
    et = [U.eta(t, l) for l in range(3)]
    ceta = Lp*et[0]*et[1] + et[1]*et[1] + et[0]*et[2]
    return dict(rem=rem, rho=rho, ceta=ceta, L=Lp, RP=RP, W1=Rk[1], th=th)
def ramp_boxes():
    ts = [(i/16, (i + 1)/16) for i in range(16)]
    ws = [(0.0, 0.25), (0.25, 0.5), (0.5, 1.0), (1.0, 2.0), (2.0, 4.0), (4.0, 8.0), (8.0, WT), (WT, math.inf)]
    return [(tt, ww, (A_W, A_W), (float(G_LO), float(G_HI))) for tt in ts for ww in ws]
def mkramp(coef=None, Y=None, which='K'):
    def f(box, point=False):
        q = ramp_q(box, point)
        if which == 'K':
            tot = arb(0)
            for m, c in coef.items(): tot += c*abs(q['rem'][m])
            if Y is not None: tot += Y*abs(q['rho'])
            return tot
        if which in ('rho', 'ceta', 'L', 'W1'): return q[which]
        if which.startswith('RP'): return q['RP'][int(which[2])]
        raise ValueError(which)
    f.scales = [0.6931, 1.0, 1.0, 1e-5]
    return f

# ================================================================== coefficient arrays (ramp box)
M0R, M1R, UR = Q_(11, 100), Q_(52, 100), Q_(1, 80)
alpha1 = A(Q_(15, 100000)); p_ = A(U.PQ)
Uramp = A(Q_(13056, 10000))*alpha1*(LN56*p_).exp(); assert Uramp.upper() <= A(UR)
def coef_arrays(Ab):
    C = {}
    for j in (0, 1, 2):
        n = j + 2; C[j] = {k: {} for k in range(j + 2)}
        for (a, b), mons in T[j].items():
            if (a, b) == (n, n): continue
            for k in range(j + 2):
                w = 0
                if a == k: w += Q_(1) if a == b else Q_(1, 2)
                elif b == k: w += Q_(1, 2)
                if w == 0: continue
                for (m, i, r), co in mons.items():
                    if m == 0: continue
                    C[j][k][m] = C[j][k].get(m, arb(0)) + A(w*co)*Ab[(i, r)]
    return C
def coef_Z(At):
    Kz = {}; morders = {}
    for j in (0, 1, 2):
        n = j + 2; Kz[j] = {}
        for (a, b), mons in T[j].items():
            if (a, b) == (n, n): continue
            for (m, i, r), co in mons.items():
                if m == 0: continue
                for X, v in At[(i, r)].items():
                    if not X: continue
                    key = tuple(sorted([a, b] + list(X))); Kz[j].setdefault(key, {}); Kz[j][key][m] = Kz[j][key].get(m, arb(0)) + A(co)*v
                    morders.setdefault(f'{j}|{key}', set()).add(m)
    return Kz, {k: sorted(v) for k, v in morders.items()}
AsymR = U.build_At_sym(UR); AtR = U.eval_At(AsymR, M0R, M1R); AbR = U.ab_from_At(AtR)
C = coef_arrays(AbR); KZC, MORD = coef_Z(AtR)

print(f'U2: case {args.case}, Y_next {[fup(y) for y in YNEXT]}, g in [{float(G_LO):.4g}, {float(G_HI):.4g}], a_W = {A_W}', flush=True)
GEOM = {w: cached(f'geom|{w}', lambda w=w: bb(mkramp(which=w), ramp_boxes())) for w in ('rho', 'ceta', 'L', 'W1', 'RP1', 'RP2', 'RP3')}
K = {}
for j in (0, 1, 2):
    K[j] = {}
    for k in range(j + 2):
        Y = YNEXT[j] if k == j else None
        K[j][k] = cached(f'K|{j}{k}', lambda j=j, k=k, Y=Y: bb(mkramp(C[j][k], Y), ramp_boxes()))['ub']
Kpair = {j: {k: (cached(f'Kpair|{j}{k}', lambda j=j, k=k: bb(mkramp(C[j][k], None), ramp_boxes()))['ub'] if k == j else K[j][k]) for k in range(j + 2)} for j in (0, 1, 2)}
KZ = {j: {str(key): cached(f'Kz|{j}|{key}', lambda cm=cm: bb(mkramp(cm, None), ramp_boxes()))['ub'] for key, cm in KZC[j].items()} for j in (0, 1, 2)}

# ================================================================== folding of the cutoff x jet classes (box sups on the ramp)
SUPL = {0: A(M0R), 1: A(M1ZG)}
def augment(Kbase):
    Ka = {j: {k: arb(v) for k, v in d.items()} for j, d in Kbase.items()}; cubic = {}; folded = []
    for j in (0, 1, 2):
        for key, v in KZ[j].items():
            idx = list(eval(key)); w = arb(v)
            while len(idx) > 2 and any(x < 2 for x in idx):
                x = min([x for x in idx if x < 2], key=lambda x: float(SUPL[x].mid())); idx.remove(x); w = w*SUPL[x]
            if len(idx) > 2: cubic[f'{j}|{key}'] = fup(arb(v)); continue
            a_, b_ = sorted(idx); folded.append(dict(level=j, cls=key, to=[a_, b_], weight=fup(w)))
            for kk, ww in (((a_, w),) if a_ == b_ else ((a_, w/2), (b_, w/2))):
                Ka[j][kk] = Ka[j].get(kk, arb(0)) + ww
    return {str(j): {str(k): fup(v) for k, v in sorted(d.items())} for j, d in Ka.items()}, cubic, folded
Kaug, CUBIC, FOLDED = augment(K)
assert set(CUBIC) <= {'2|(2, 2, 2)', '2|(2, 2, 3)'}, CUBIC
assert all(int(k) <= 3 for d in Kaug.values() for k in d)

# ================================================================== foot / S33 weights (global box of the case)
AsymG = GA.build_At_sym(UG); AtG = GA.eval_At(AsymG, M0G, M1G)
etap = A(Q_(15, 8))/LN2; DEL = A(Q_(1, 10)); M1Z = A(M1ZG); aminG = (1 + A(M0G)**2)**-2; FOOT = {}
for j in (0, 1, 2):
    add, dis, spare, ndfr, items = GA.l4_extra_YC_x(T, j, M0G, M1ZG, AtG, GamP, KAPP, dfrac=Q_(2, 25), s33='ibp', M0zq=M0ZG)
    fw = arb(0); w33 = sum((it[1] for it in items if it[0] == 'S33'), arb(0)); has33 = any(it[0] == 'S33' for it in items); fadd = {2: arb(0), 3: arb(0)}
    if has33:
        c1 = w33*M1Z*etap; fadd[3] += c1*KAPP*DEL; fw += c1*KAPP/DEL; dd = A(Q_(1, 10))*aminG; fw += c1**2/(4*dd)
    for it in items:
        if it[0] != 'L4': continue
        w = it[2]
        if it[1] == 3: fadd[3] += w*M1Z*etap*DEL; fw += w*M1Z*etap/DEL
        elif it[1] == 2: fadd[2] += w*M1Z*etap; fw += w*M1Z*etap
        else: raise ValueError(('foot r', it))
    ok, sl = GA.budget_ok(dis, ndfr, M0G, Q_(1, 10) if has33 else Q_(0)); assert ok
    for k in (2, 3):     # consistency: the foot adds are inside the u1 Y_full of the case
        if fadd[k] != 0: assert fup(fadd[k]) <= CASE['Y_full'][str(j)][str(k)]
    FOOT[j] = dict(footw=num(fw), W33=num(w33), foot_Young_add_to_Y={str(k): fup(v) for k, v in fadd.items()}, budget_ok=bool(ok), budget_slack=sl)
CETA = arb(GEOM['ceta']['ub']); LMAX = arb(GEOM['L']['ub'])
Kf_S33 = arb(FOOT[2]['W33']['double'])*CETA/3
W222 = arb(CUBIC['2|(2, 2, 2)']); W223 = arb(CUBIC['2|(2, 2, 3)'])
c33 = Kf_S33 + W222
RRAT = (LN2*LMAX).exp()

out = dict(case=args.case, phi=not args.nophi, g_range=[str(G_LO), str(G_HI)], a_W=str(A_W), ramp=[28, 56], centres='s_rho >= ln 56',
           ramp_box=dict(M0=str(M0R), M1=str(M1R), U=str(UR), fold_sup_z0=str(M0R), fold_sup_z1=str(M1ZG)),
           global_box=dict(BOXG), Y_next=[fup(y) for y in YNEXT],
           K={str(j): {str(k): v for k, v in d.items()} for j, d in K.items()}, K_pair_only={str(j): {str(k): v for k, v in d.items()} for j, d in Kpair.items()},
           K_aug=Kaug, folded_classes=FOLDED, cubic_remaining=CUBIC, cubic_classes_raw={str(j): d for j, d in KZ.items()},
           class_m_orders=MORD, geometry={k: v['ub'] for k, v in GEOM.items()},
           foot={str(j): v for j, v in FOOT.items()}, Kf_unnorm_S33=num(Kf_S33), C_eta=num(CETA), L_max=num(LMAX), R_ratio=num(RRAT),
           c33=num(c33, rule='Kf_S33 (= W33 C_eta/3, exact foot) + weight of (2,2,2); S33 data form'),
           c223=num(W223, rule='sup |c| of the (2,2,3) class per eps*Theta*phi; Cauchy-Schwarz route (no IBP); m-orders ' + str(MORD.get('2|(2, 2, 3)'))),
           bb={k: dict(lb=v['lb'], nev=v['nev'], left=v['left']) for k, v in R.items() if isinstance(v, dict) and 'nev' in v})
out['args'] = L.args_block([U.TABLES, args.u1, os.path.join(L.P2, 'ulib.py'), os.path.join(L.P1, 'ghat_arb.py'), os.path.abspath(__file__), os.path.join(HERE, 't4lib.py')],
                           dict(maxeval=args.maxeval, tol=args.tol))
L.dump(out, OUT)
print('K_aug', Kaug); print('cubic', CUBIC, 'm-orders (2,2,3):', MORD.get('2|(2, 2, 3)'))
print('footw', [FOOT[j]['footw']['double'] for j in range(3)], 'Kf_S33', fup(Kf_S33), 'c33', fup(c33), 'c223', fup(W223), 'C_eta', fup(CETA), 'L_max', fup(LMAX), 'R', fup(RRAT))
print('wrote', OUT)
