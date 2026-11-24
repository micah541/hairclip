#!/usr/bin/env python3
"""Lemma U, Phase 1, P2: certified (arb) inner-ramp constants of the Lemma-U energy on y in [28, 56] (rescaled).

  python3 -B ramp_K_cert.py [--maxeval N] [--tol x] [--fresh]          (single process; resumable: ramp_K_cert.partial.json)

Output ramp_K_cert.json (keys documented in P2_REPORT.md).  Main consumer keys:
  ["4.0"]["K"]           same shape as feas_P.json["4.0"]["K"]: {j: {k: K_jk}} (floats = arb upper endpoints rounded up), Lemma-U weight, crude box,
                          rho part with the global Y_{j+1} INCLUDING the general-weight L^4 extras;  t1/t5 can read it in place of ../feas_P.json.
  ["variants"][name]["K"] the same for the other variants (no-L4 Y, small (A) core box, power weight, feas_P reproduction).
Definitions (s = ln y, eps = e^{-4s}, eta = S(t), t = (s - ln28)/ln2, Theta = W(sigma) D(s)):
  rem_m  = [ (eta^2 Theta)^(m) - eta^2 Theta^(m) ] / Theta = sum_{l=1..m} C(m,l) (eta^2)^(l) Theta^(m-l)/Theta      (cutoff part of the weight derivative)
  rho^   = (1/2) | (eta^2)'' + 2 L (eta^2)' |,   L = (ln eps Theta)' = -4 + W'/W - 4 theta                            (interpolation remainder, per unit Y)
  K_jk   = sup_{ramp x centres x a x g} [ sum_{pairs (a,b) != (n,n)} w_k(a,b) sum_{(m,i,r), m>=1} |coef| Ab_(i,r) |rem_m| + [k = j] Y_{j+1} rho^ ]
           (w_k = 1 for a = b = k, 1/2 for one index = k; Ab = Z-free part of the A-coefficient bounds on the ramp box);  dimensionless (divided by eps Theta).
Honest: |coef| tables, |rem_m| enclosed per cell (no sign cancellation between monomials), Y_{j+1} = sum |coef| Gamma_m Ab with certified global Gamma_m."""
import sys, os, json, math, time, hashlib, argparse, platform
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import sympy as sp
from flint import arb, ctx
import ulib as U
from ulib import A, Q_, hull, LN2, LN28, LN56

FEASP = '/home/claude/s1spec/writeup/step3_gap_v32/feas_P.json'
ASSEMBLY = '/home/claude/s1spec/writeup/step3_gap_v32/gonogo/lemmaS_A_rerun/rigor/assembly_arb.json'
OUT = os.path.join(HERE, 'ramp_K_cert.json'); PART = os.path.join(HERE, 'ramp_K_cert.partial.json')

ap = argparse.ArgumentParser()
ap.add_argument('--maxeval', type=int, default=40000)
ap.add_argument('--tol', type=float, default=0.002)
ap.add_argument('--fresh', action='store_true')
args = ap.parse_args()

def sha(p): return hashlib.sha256(open(p, 'rb').read()).hexdigest()
ARGS = dict(argv=sys.argv, cwd=os.getcwd(), python=sys.version, platform=platform.platform(), flint=__import__('flint').__version__,
            env={k: os.environ.get(k) for k in ('PYTHONHASHSEED', 'OMP_NUM_THREADS') if os.environ.get(k) is not None},
            inputs={p: sha(p) for p in (U.TABLES, FEASP, ASSEMBLY, os.path.join(HERE, 'ulib.py'), os.path.abspath(__file__))},
            params=dict(mu=U.MU, nu=U.NU, a=[str(U.A_LO), str(U.A_HI)], g=[str(U.G_LO), str(U.G_HI)], ramp=[28, 56], centres='s_rho >= ln 56',
                        maxeval=args.maxeval, tol=args.tol, prec=ctx.prec))

def fup(x):
    """float >= arb upper endpoint"""
    v = x.upper() if isinstance(x, arb) else arb(x); f = float(v)
    return f if arb(f) >= v else math.nextafter(f, math.inf)
def flo(x):
    v = x.lower(); f = float(v)
    return f if arb(f) <= v else math.nextafter(f, -math.inf)

R = {} if args.fresh or not os.path.exists(PART) else json.load(open(PART))
def save():
    json.dump(R, open(PART + '.tmp', 'w'), indent=1); os.replace(PART + '.tmp', PART)
def cached(key, fn):
    if key in R: return R[key]
    t0 = time.time(); v = fn(); v['sec'] = round(time.time() - t0, 1); R[key] = v; save()
    print(f'  {key}: {v.get("ub", "")}  (lb {v.get("lb", "")}, {v.get("nev", "")} cells, {v["sec"]} s)', flush=True); return v
def bb(f, boxes):
    ub, lb, nev, left = U.bb_sup(f, boxes, tol=args.tol, maxeval=args.maxeval)
    return dict(ub=ub, lb=lb, nev=nev, left=left)

T = U.load_tables()
G_BALL = hull(A(U.G_LO), A(U.G_HI))

# ================================================================== 0. profile bound on the ramp
alpha1 = A(Q_(15, 100000)); p_ = A(U.PQ)
Uramp = A(Q_(13056, 10000))*alpha1*(arb(56).log()*p_).exp()
UR_Q = Q_(1, 80); assert Uramp.upper() <= A(UR_Q)
R['U_ramp'] = dict(value=fup(Uramp), used=str(UR_Q), note='[cert arith; paper: |U_rho| <= 1.3056 alpha1 y^p on y <= 56, |alpha| <= alpha1 = 1.5e-4, lambda(tau) <= lambda(tau_rho)]')

# ================================================================== 1. global weight constants over supp: sigma in R, s >= max(ln28, sigma + ln56)
ST = 14.0
def theta_hi(slo):
    """upper bound of theta on {s >= max(ln28, sigma + ln56), sigma >= slo}; theta decreases in s"""
    if slo == -math.inf or not (arb(slo) + LN2 > 0): smin = LN28          # slo <= -ln2 (or undecided): s = ln28 is admissible
    else: smin = arb(slo) + LN56
    return (1/(1 + (4*smin).exp()/A(U.G_HI))).upper()
def glob_q(box, point=False):
    (slo, shi), (alo, ahi) = box
    a = arb((alo + ahi)/2) if point else hull(alo, ahi)
    if shi == math.inf: p, q = U.pq_upper_tail(arb(slo), a)
    elif slo == -math.inf: p, q = U.pq_lower_tail(arb(shi), a)
    else:
        sg = arb((slo + shi)/2) if point else hull(slo, shi); p, q = U.pq_from_sigma(sg, a)
    th = hull(0, theta_hi(slo))
    Rk = U.Rk(p, q)
    B = [None, arb(-4), 16 - 8*Rk[1], -12*Rk[2] + 48*Rk[1] - 64]
    L = -4 + Rk[1] - 4*th; Lp = Rk[2] - Rk[1]*Rk[1] + 16*th*(1 - th)
    return Rk, B, th, L, Lp
def glob_boxes():
    cuts = [-ST, -8, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 8, ST]
    bx = [((-math.inf, -ST), (0.0, 0.3))] + [((cuts[i], cuts[i + 1]), (0.0, 0.3)) for i in range(len(cuts) - 1)] + [((ST, math.inf), (0.0, 0.3))]
    return bx
def mkglob(which):
    def f(box, point=False):
        Rk, B, th, L, Lp = glob_q(box, point)
        if which == 'kap': return L
        if which == 'kap2': return L*L + Lp
        if which == 'W1': return Rk[1]
        m = int(which[1]); return Rk[m] + th*B[m]
    f.scales = [1.0, 1.0]; return f
print('1. global weight constants (supp, centres s_rho >= ln56)', flush=True)
GL = {w: cached(f'glob_{w}', lambda w=w: bb(mkglob(w), glob_boxes())) for w in ('G1', 'G2', 'G3', 'kap', 'kap2')}
Gm = {0: arb(1)}
for m in (1, 2, 3): Gm[m] = arb(GL[f'G{m}']['ub'])
KAP = arb(GL['kap']['ub']); KAP2 = arb(GL['kap2']['ub'])

# ================================================================== 2. global Y_{j+1} (crude box), with and without the general-weight L^4 extras
M0C, M1C = Q_(11, 100), Q_(52, 100); M0Z = Q_(11, 100); UGLOB = Q_(11, 100)
print('2. global YC / Ghat (crude box)', flush=True)
AsymG = U.build_At_sym(UGLOB); AtG = U.eval_At(AsymG, M0C, M1C); AbG = U.ab_from_At(AtG)
YG = {}
for dfrac in (Q_(1, 10), Q_(2, 25)):
    for j in (0, 1, 2):
        YCb = U.YC_gen(T[j], j, AbG, Gm); add, dis, ndp, items = U.l4_extra_gen(T, j, M0C, M0Z, M1C, AtG, Gm, KAP, dfrac)
        YC = {k: YCb[k] + add.get(k, arb(0)) for k in YCb}; ok, sl = U.budget(dis, ndp, dfrac, M0C)
        YG[(str(dfrac), j)] = dict(YCb=YCb, YC=YC, Ghat_noL4=U.ghat_gen(YCb, j, M0C, KAP2)[j], Ghat_L4=U.ghat_gen(YC, j, M0C, KAP2)[j], budget_ok=ok, slack=sl, add=add)
R['global'] = dict(Gamma=[fup(Gm[m]) for m in (1, 2, 3)], kappa=fup(KAP), kappa2=fup(KAP2), box=[str(M0C), str(M1C)], M0z=str(M0Z), Umax=str(UGLOB),
                   Y={f'dfrac={d}': {str(j): dict(Y_next_noL4=fup(YG[(d, j)]['YCb'][j + 1]), Y_next_L4=fup(YG[(d, j)]['YC'][j + 1]),
                                              Ghat_jj_noL4=fup(YG[(d, j)]['Ghat_noL4']), Ghat_jj_L4=fup(YG[(d, j)]['Ghat_L4']), budget_ok=YG[(d, j)]['budget_ok'], slack=YG[(d, j)]['slack'])
                                      for j in (0, 1, 2)} for d in ('1/10', '2/25')})
print('  Gamma_m', R['global']['Gamma'], 'kappa~', R['global']['kappa'], 'kappa2~', R['global']['kappa2'])
for d in ('1/10', '2/25'): print('  dfrac', d, {j: (round(v['Y_next_noL4'], 2), round(v['Y_next_L4'], 2), '%.4g' % v['Ghat_jj_L4']) for j, v in R['global']['Y'][f'dfrac={d}'].items()})
save()

# ================================================================== 5. foot / S33 weights (port of weights_exact.ramp, kappa -> kappa~, gam^m -> Gamma_m), global crude box
print('5. ramp L4 foot weights', flush=True)
FOOT = {}
etap = A(Q_(15, 8))/LN2; DEL = A(Q_(1, 10)); M1 = A(M1C); amin = (1 + A(M0C)**2)**-2; 
for j in (0, 1, 2):
    dfrac = Q_(2, 25)
    YCb = U.YC_gen(T[j], j, AbG, Gm); add, dis, ndp, items = U.l4_extra_gen(T, j, M0C, M0Z, M1C, AtG, Gm, KAP, dfrac)
    addc = dict(add); fadd = {k: arb(0) for k in add}
    fw = arb(0); extra = Q_(0); w33 = sum((it[1] for it in items if it[0] == 'S33'), arb(0)); has33 = any(it[0] == 'S33' for it in items)
    if has33:
        c1 = w33*M1*etap; fadd[3] += c1*KAP*DEL; fw += c1*KAP/DEL
        dd = A(Q_(1, 10))*amin; extra = Q_(1, 10); fw += c1**2/(4*dd)
    for it in items:
        if it[0] != 'L4': continue
        w = it[2]
        if it[1] == 3: fadd[3] += w*M1*etap*DEL; fw += w*M1*etap/DEL
        elif it[1] == 2: fadd[2] += w*M1*etap; fw += w*M1*etap
        else: raise ValueError(('foot r', it))
    ok, sl = U.budget(dis, ndp, dfrac, M0C, extra)
    YCtot = {k: YCb[k] + addc.get(k, arb(0)) + fadd.get(k, arb(0)) for k in YCb}
    FOOT[j] = dict(footw=fup(fw), W33=fup(w33), foot_Young_add_to_Y={str(k): fup(v) for k, v in fadd.items() if v != 0},
                   budget_ok=ok, budget_slack_over_amin=sl, dfrac='2/25', Ghat_jj_with_foot_adds=fup(U.ghat_gen(YCtot, j, M0C, KAP2)[j]))
    print('  level', j, {k: v for k, v in FOOT[j].items() if not k.startswith('_')}, flush=True)
    FOOT[j]['_fadd'] = fadd; FOOT[j]['_Yfull_next'] = YCtot[j + 1]; FOOT[j]['_w33'] = w33
    FOOT[j]['Y_next_full'] = fup(YCtot[j + 1])


# ================================================================== 3. ramp geometry: box (t, w, a), sigma = (t - 1) ln2 - w  (w in [0, inf): all centres s_rho >= ln56)
WT = 16.0
def ramp_q(box, mode, point=False):
    """returns dict of enclosures on the box: rem[1..3], rho, ceta, Lr, RT[1..3], epsTheta log-derivative"""
    tb = box[0]
    t = arb((tb[0] + tb[1])/2) if point else hull(tb[0], tb[1])
    E2 = [U.e2(t, l) for l in range(4)]
    if mode == 'power':                      # Theta = y^mu (feas_P / power reference): Theta^(k)/Theta = mu^k, L = mu - 4 = 0
        Rk = [arb(1), arb(4), arb(16), arb(64)]; th = arb(0)
    else:
        (wlo, whi), (alo, ahi) = box[1], box[2]
        a = arb((alo + ahi)/2) if point else hull(alo, ahi)
        if whi == math.inf:
            sigT = (arb(tb[1]) - 1)*LN2 - wlo
            p, q = U.pq_lower_tail(sigT if not point else (t - 1)*LN2 - wlo, a)
        else:
            w = arb((wlo + whi)/2) if point else hull(wlo, whi)
            p, q = U.pq_from_sigma((t - 1)*LN2 - w, a)
        Rk = U.Rk(p, q)
        gl, gh = box[3]
        gb = arb((gl + gh)/2) if point else hull(gl, gh)
        s = LN28 + t*LN2; th = U.theta_s(s, gb)
    # linear-in-theta forms (theta occurs once in each expression)
    rem = {1: E2[1],
           2: (2*E2[1]*Rk[1] + E2[2]) + th*(-8*E2[1]),
           3: (3*E2[1]*Rk[2] + 3*E2[2]*Rk[1] + E2[3]) + th*(3*E2[1]*(16 - 8*Rk[1]) - 12*E2[2])}
    L0 = -4 + Rk[1]
    rho = (E2[2] + 2*L0*E2[1] + th*(-8*E2[1]))/2
    et = [U.eta(t, l) for l in range(3)]
    ceta = L0*et[0]*et[1] + et[1]*et[1] + et[0]*et[2] + th*(-4*et[0]*et[1])
    RT = {1: Rk[1] - 4*th, 2: Rk[2] + th*(16 - 8*Rk[1]), 3: Rk[3] + th*(-12*Rk[2] + 48*Rk[1] - 64)}
    return dict(rem=rem, rho=rho, ceta=ceta, L=L0 - 4*th, RT=RT, W1=Rk[1], th=th)
def ramp_boxes(mode):
    ts = [(i/16, (i + 1)/16) for i in range(16)]
    if mode == 'power': return [(tt,) for tt in ts]
    ws = [(0.0, 0.25), (0.25, 0.5), (0.5, 1.0), (1.0, 2.0), (2.0, 4.0), (4.0, 8.0), (8.0, WT), (WT, math.inf)]
    return [(tt, ww, (0.0, 0.3), (float(U.G_LO), float(U.G_HI))) for tt in ts for ww in ws]
def mkramp(mode, coef=None, Y=None, which='K'):
    """coef: {m: arb} (sum |coef| Ab per m);  Y: arb or None (rho part);  which: 'K' | 'rho' | 'ceta' | 'L' | 'RT1'.. | 'W1'"""
    def f(box, point=False):
        q = ramp_q(box, mode, point)
        if which == 'K':
            tot = arb(0)
            for m, c in coef.items(): tot += c*abs(q['rem'][m])
            if Y is not None: tot += Y*abs(q['rho'])
            return tot
        if which == 'rho': return q['rho']
        if which == 'ceta': return q['ceta']
        if which == 'L': return q['L']
        if which == 'W1': return q['W1']
        if which.startswith('RT'): return q['RT'][int(which[2])]
        raise ValueError(which)
    f.scales = [0.6931, 1.0, 1.0, 1e-5] if mode != 'power' else [1.0]
    return f

# ================================================================== 4. ramp coefficient arrays
def coef_arrays(Ab, gamsign=None):
    """C[j][k][m] = sum_{pairs} w_k(a,b) sum_{(m,i,r) m>=1} |coef| Ab_(i,r).   Also returns YC_core with weight factor gamsign^m (feas_P: -4, signed) if given."""
    C = {}; YCs = {}
    for j in (0, 1, 2):
        n = j + 2; C[j] = {k: {} for k in range(j + 2)}; YC = {k: arb(0) for k in range(n + 1)}
        for (a, b), mons in T[j].items():
            if (a, b) == (n, n): continue
            if gamsign is not None:
                gC = sum((A(co)*arb(gamsign)**m*Ab[(i, r)] for (m, i, r), co in mons.items()), arb(0))
                if a == b: YC[a] += gC
                else: YC[a] += gC/2; YC[b] += gC/2
            for k in range(j + 2):
                w = 0
                if a == k: w += Q_(1) if a == b else Q_(1, 2)
                elif b == k: w += Q_(1, 2)
                if w == 0: continue
                for (m, i, r), co in mons.items():
                    if m == 0: continue
                    C[j][k][m] = C[j][k].get(m, arb(0)) + A(w*co)*Ab[(i, r)]
        YCs[j] = YC
    return C, YCs
def coef_Z(At):
    """cubic (and higher) cutoff-commutator classes: Kz[j][key][m] = sum |coef| v_X, key = sorted([a, b] + X), X a nonempty Z-multiset of A_(i,r)"""
    Kz = {}
    for j in (0, 1, 2):
        n = j + 2; Kz[j] = {}
        for (a, b), mons in T[j].items():
            if (a, b) == (n, n): continue
            for (m, i, r), co in mons.items():
                if m == 0: continue
                for X, v in At[(i, r)].items():
                    if not X: continue
                    key = tuple(sorted([a, b] + list(X))); Kz[j].setdefault(key, {}); Kz[j][key][m] = Kz[j][key].get(m, arb(0)) + A(co)*v
    return Kz

VARIANTS = {}
# (a) Lemma-U weight, crude box on the ramp, deviation formulation (U = U_ramp), Ab = Z-free part of At
AsymR = U.build_At_sym(UR_Q)
for vname, (M0q, M1q) in (('crude', (M0C, M1C)), ('smallA', (Q_(168, 10000), Q_(108, 1000)))):
    AtR = U.eval_At(AsymR, M0q, M1q); AbR = U.ab_from_At(AtR); C, _ = coef_arrays(AbR)
    VARIANTS[vname] = dict(box=[str(M0q), str(M1q)], U=str(UR_Q), C=C, Kz=coef_Z(AtR), At=AtR)
# (b) feas_P reproduction (Ab of A_bounds_s with U = 1e-12, signed gam = -4 in YC)  and power-weight honest (|coef| 4^m)
AbF = U.ab_exact(M0C, M1C, Q_(1, 10**12))
CF, YCF = coef_arrays(AbF, gamsign=-4)
YCP = {j: U.YC_gen(T[j], j, AbF, {m: arb(4)**m for m in range(4)}) for j in (0, 1, 2)}

print('3. ramp constants', flush=True)
out = {}
# ---- feas_P reproduction and power-weight honest (t only)
for vname, YCx in (('feasP_repro', YCF), ('power_honest', YCP)):
    K = {}
    for j in (0, 1, 2):
        K[j] = {}
        for k in range(j + 2):
            Y = YCx[j][j + 1] if k == j else None
            K[j][k] = cached(f'K|{vname}|{j}{k}', lambda j=j, k=k, Y=Y: bb(mkramp('power', CF[j][k], Y), ramp_boxes('power')))['ub']
    out[vname] = dict(K={str(j): {str(k): v for k, v in d.items()} for j, d in K.items()},
                      Y_next={str(j): float(YCx[j][j + 1].mid()) for j in (0, 1, 2)})
# ---- Lemma-U weight
GEOM = {}
for w in ('rho', 'ceta', 'L', 'W1', 'RT1', 'RT2', 'RT3'):
    GEOM[w] = cached(f'geom|{w}', lambda w=w: bb(mkramp('U', which=w), ramp_boxes('U')))
for vname in ('crude', 'smallA'):
    C = VARIANTS[vname]['C']
    for yname in ('none', 'full', 'L4', 'noL4'):
        if vname == 'smallA' and yname in ('noL4', 'L4'): continue
        K = {}
        for j in (0, 1, 2):
            K[j] = {}
            for k in range(j + 2):
                if yname != 'none' and k != j: K[j][k] = R[f'K|{vname}|none|{j}{k}']['ub']; continue
                Y = None if yname == 'none' else arb(fup(FOOT[j]['_Yfull_next'] if yname == 'full' else YG[('1/10', j)]['YC' if yname == 'L4' else 'YCb'][j + 1]))
                K[j][k] = cached(f'K|{vname}|{yname}|{j}{k}', lambda C=C, j=j, k=k, Y=Y: bb(mkramp('U', C[j][k], Y), ramp_boxes('U')))['ub']
        out[f'U_{vname}_Y{yname}'] = dict(K={str(j): {str(k): v for k, v in d.items()} for j, d in K.items()})
# ---- cubic cutoff-commutator classes (not in feas_P / t1): crude box
KZ = {}
for j in (0, 1, 2):
    KZ[j] = {}
    for key, cm in VARIANTS['crude']['Kz'][j].items():
        KZ[j][str(key)] = cached(f'Kz|crude|{j}|{key}', lambda cm=cm: bb(mkramp('U', cm, None), ramp_boxes('U')))['ub']

CETA = arb(GEOM['ceta']['ub'])
for j in (0, 1, 2): FOOT[j]['Kf_unnorm_S33'] = fup(FOOT[j]['_w33']*CETA/3)
R['foot'] = {str(j): {k: v for k, v in d.items() if not k.startswith('_')} for j, d in FOOT.items()}; save()
# ---- fold the cubic cutoff classes into K where box sups reduce them to pairs (|z_0| <= M0z, |z_1| <= M1); remainder reported as genuinely cubic
SUPL = {0: A(M0Z), 1: A(M1C)}
def augment(Kbase):
    Ka = {j: {k: arb(v) for k, v in d.items()} for j, d in Kbase.items()}; cubic = {}
    for j in (0, 1, 2):
        for key, v in KZ[j].items():
            idx = list(eval(key)); w = arb(v)
            while len(idx) > 2 and any(x < 2 for x in idx):
                lows = sorted([x for x in idx if x < 2], key=lambda x: float(SUPL[x].mid())); x = lows[0]; idx.remove(x); w = w*SUPL[x]
            if len(idx) > 2: cubic[f'{j}|{key}'] = fup(arb(v)); continue
            a_, b_ = sorted(idx)
            for kk, ww in (((a_, w),) if a_ == b_ else ((a_, w/2), (b_, w/2))):
                Ka[j][kk] = Ka[j].get(kk, arb(0)) + ww
    return {str(j): {str(k): fup(v) for k, v in sorted(d.items())} for j, d in Ka.items()}, cubic
Kfull = {int(j): {int(k): v for k, v in d.items()} for j, d in out['U_crude_Yfull']['K'].items()}
Kaug, CUBIC = augment(Kfull); out['U_crude_Yfull_aug'] = dict(K=Kaug)
print('augmented K (crude, full Y):', Kaug); print('genuinely cubic ramp classes (weights, raw):', CUBIC)

# ================================================================== 6. comparison and effective KS against t1's inputs
FP = json.load(open(FEASP))['4.0']['K']
AS = json.load(open(ASSEMBLY)); EA = AS['E']; gpA = AS['gp']
pp = 1.0246; amin_t1 = (16/17)**2; eps28 = 28.**-4
dUmax = 1.3056*56**pp*(3.33811e-9 + pp*2.19103e-7*1.5e-4)/(pp*(0.25 - 2.19103e-7)); dU = [(pp**j*dUmax)**2*math.log(2) for j in range(3)]
Ek = [(math.sqrt(EA[k]) + math.sqrt(dU[k]))**2 for k in range(3)]
X = {1: eps28*Ek[1], 2: gpA[0]*EA[0]/(2*amin_t1) + eps28*pp**4*dU[0], 3: gpA[1]*EA[1]/(2*amin_t1) + eps28*pp**6*dU[0]}
def t1_ramp(K):
    return [sum(float(v)*(eps28*Ek[int(k)] if int(k) <= j else X[int(k)]) for k, v in K[str(j)].items()) for j in range(3)]
used = t1_ramp(FP); cmpres = {}
for name, v in out.items():
    hr = t1_ramp(v['K'])
    cmpres[name] = dict(ratio_per_entry={j: {k: v['K'][j][k]/FP[j][k] for k in FP[j]} for j in FP}, extra_entries={j: {k: x for k, x in v['K'][j].items() if k not in FP[j]} for j in FP}, t1_ramp_bracket=hr,
                        effective_KS_per_level=[h/u for h, u in zip(hr, used)])
R['comparison'] = dict(feasP_K=FP, t1_ramp_bracket_used=used, t1_inputs=dict(EA=EA, gp=gpA, dU=dU, X=X), variants=cmpres)
print('6. effective KS (honest t1 ramp bracket / bracket with feas_P K):')
for name, c in cmpres.items(): print('  ', name, ['%.4g' % x for x in c['effective_KS_per_level']])

res = {'4.0': dict(K=out['U_crude_Yfull_aug']['K'], note='Lemma-U weight (4,4,a in [0,0.3]), g in [3e5,4.5e5], crude box (0.11,0.52); rho part with the global Y_{j+1} incl. general-weight L4 extras and ramp-foot Young adds (dfrac 2/25); cubic cutoff classes folded in by box sups; EXTRA entry K[1][3] (z_3 at level 1, data X_3) is not read by t1/t5 as written; arb upper endpoints'),
       'K_t1_shape_no_cubic': out['U_crude_Yfull']['K'], 'cubic_remaining': CUBIC,
       'variants': out, 'geometry_ramp': {k: v for k, v in R.items() if k.startswith('geom|')}, 'global': R['global'], 'foot': R['foot'],
       'C_eta_Xi': GEOM['ceta']['ub'], 'cubic_cutoff_classes_crude': {str(j): d for j, d in KZ.items()}, 'U_ramp': R['U_ramp'],
       'comparison': R['comparison'], 'bb': {k: dict(lb=v['lb'], nev=v['nev'], left=v['left']) for k, v in R.items() if isinstance(v, dict) and 'nev' in v},
       'args': ARGS}
json.dump(res, open(OUT + '.tmp', 'w'), indent=1, default=str); os.replace(OUT + '.tmp', OUT)
print('wrote', OUT)
