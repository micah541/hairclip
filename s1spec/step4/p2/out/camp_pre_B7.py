"""T8 certified: Lemma J and Corollary C (LEMMA_F.md sections 4-5) in arb, on the certified per-slab constants.
    python camp.py [--consts constants/constants.json] [--out constants/amp.json] [--F0 0,0,0,0]

For levels j = 0..3 and the slabs of constants.json (v2 base: 0..180; v1: 5..187) it computes, for every slab m, upper bounds
valid for ALL tau in slab m:
    sqrt F_j(tau) <= Ybar * aF[j][m] + bF[j][m]
(aF: amplification per unit Ybar; bF: residual part (T6 norms) + initial data F_j(tau0)).  Then the hand-off items
over the good-time window = the trailing slabs with k >= 2.0e4 at their start (v2: slabs 177-180, k in [2.04e4, 2.37e4];
v1: slabs 184-187, k in [2.03e4, 2.36e4]; --window a-b overrides):
    item 3  (E3 window energy, Stirling C3)      -> Ybar_3
    item 4  (K4 via D2 at a good time tau*, C5)    -> Ybar_4 (bisection, final value re-verified)
    item 2  (shell energies E0..E2)               -> Ybar_2
Every comparison is done on arb balls; the reported Ybar values are certified for the stated inputs (T7 values and
placeholders as recorded in constants.json; T1 base window energies and T9 chart factor as in amp.py)."""
import os, sys, json, math, argparse
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path: sys.path.insert(0, HERE)
from flint import arb, ctx
ctx.prec = 64
import cconst
from cconst import U, L, iv

T6DIR = cconst.T6DIR
A_S3, NU = arb('0.01'), arb('0.001')
EPSS = arb('1e-3')                        # (u + v)^2 <= (1 + eps) u^2 + (1 + 1/eps) v^2 for the Ybar / residual split
YSTAR = arb('1e-12')                      # hypothesis: the majorant Ybar is run with Ybar >= 1e-12 (it is ~2.3e-12)
DELTA_I0 = arb('1e-3')                    # Young weight for the residual terms of (I0)
BASE_T1_V1 = dict(E=(1.45e-8, 5.2e-8, 1.7e-7, 4.6e-7), K4=2.1e-4, shell=(8.0e-11, 3.4e-10, 1.4e-9))   # T1 report, v1 base, k = 2.3e4
# v2 base (polished_v2, t1/handoff.py P11v2, float, provisional T9 input): values at k = 20400 <= window start k = 20414.8
# (all six samples k = 20400..23700 decrease in k), rounded up
BASE_T1_V2 = dict(E=(2.28e-8, 8.01e-8, 2.51e-7, 6.27e-7), K4=2.33e-4, shell=(1.31e-10, 5.58e-10, 2.31e-9))
# certified (arb) base values at the hand-off sub-window W' = slab 180, t_c in [1 - 1/512, 1] (p2/t9h/win.py with the H1 IFT alpha;
# K4 from the certified T6 hulls, p2/t9.py row r_c in [200, 210])
BASE_T1_CERT = dict(E=(1.436e-8, 5.192e-8, 1.702e-7, 4.696e-7), K4=3.163e-4, shell=(2.832e-12, 3.659e-12, 8.924e-11))  # typed (review B4)
# review B4: read the certified window/shell values from T9's output and round outward; K4 is T9's row r_c in [200, 210] with
# the error parts set to 0 (same cells and formulas as p2/t9.py, slabs 177-180): 3.16322500e-4, rounded up to 3.1633e-4
_T9WIN = os.path.join(HERE, '..', 'out', 't9_win_tight.json')
if os.path.exists(_T9WIN):
    _w = json.load(open(_T9WIN)); _up = lambda v: math.nextafter(float(v), math.inf)
    BASE_T1_CERT = dict(E=tuple(_up(v) for v in _w['windows_max']), K4=3.1633e-4, shell=tuple(_up(v) for v in _w['shells']))
BASE_T1 = {'v1': BASE_T1_V1, 'v2': BASE_T1_V2}.get(os.environ.get('BASE_T1'), BASE_T1_CERT)
SUBW = os.environ.get('SUBW')                 # hand-off sub-window length in tau inside the last slab (C5 good time in it)
THR = dict(E=(2e-7, 1e-5, 1e-5, 1e-5), K4=0.06, shell=(2e-7, 1e-4, 0.2))
KCH = arb('1.02')                         # T9 chart factor (placeholder; overridden by T9JSON)
C0J = {1: 0.0, 2: 0.0, 3: 0.0}           # T9: m = 0 chain-rule terms, sqrt E_j[err] += c0_j sqrt F_0 (j >= 1)
if os.environ.get('T9JSON'):
    _t9 = json.load(open(os.environ['T9JSON'])); KCH = arb(repr(_t9['KCH'])); C0J = {int(k): v for k, v in _t9['c0'].items()}
STIR = {0: {0: 1}, 1: {1: 1}, 2: {1: 1, 2: 1}, 3: {1: 1, 2: 3, 3: 1}}
K_HAND = 2.0e4                            # the good-time window: the trailing slabs whose start has k >= 2.0e4
WINDOW = None                             # (first, last) slab; None: from K_HAND (v2 base: (177, 180); v1 base: (184, 187))
DT = 0.15                                 # slab length in tau (both bases)


def window_of(C):
    """the good-time window (first, last): last = the last slab, first = the first slab of the trailing run with k >= K_HAND"""
    per = C['per_slab']; m = len(per)
    while m > 0 and per[m - 1]['k'] >= K_HAND: m -= 1
    assert m < len(per), 'no slab with k >= K_HAND'
    return (per[m]['slab'], per[-1]['slab'])


def amin(*xs): return arb(min(U(x) for x in xs))
def amax(*xs): return arb(max(U(x) for x in xs))
def posa(x): return cconst.pos(x)
def sq(x): return cconst.sqrtp(x)


def cascade(lam, dt, P, R1=0.0, RIs=0.0):
    """Lemma J in arb.  lam: certified lower rate (arb), dt = t - tau0 >= 0.  Returns J (dict m -> arb upper bound,
    per unit Ybar^2), with the residual terms of (I0) [R1] and (I2) [R^I_s] included via Ybar >= YSTAR."""
    g = {k: (arb(v) if isinstance(v, (int, float)) else v) for k, v in P.items() if not k.startswith('_') and k not in ('cI', 'nVmin')}
    c0, c1, c2, c3 = [arb(v) for v in P['cI']]
    w0 = (-lam * dt).exp()
    B = (arb('0.5') + posa(g['rho'] - lam / 2) / lam) / g['kD']
    JA = 1 / lam
    JB = amin(4 * B / g['th0'], 1 / (lam * g['eps']))
    JD4 = 8 * B / g['eps']
    JW1 = sq(JA * JB) + g['c1'] * JA
    JW3 = sq(JB * JD4) + g['c1'] * JB
    # (I0) with its residual term R1 (V2 + 2 Psi1 V1 + Psi2 V0), Psi2 := 2 c1W, absorbed by Young with weight DELTA_I0:
    #   R1 (...) <= d (V0^2 + V2^2) + R1^2 Cres / (4 d),  Cres = (1 + Psi1)^2 + (Psi1 (1 + 2 sqrt c1W) + Psi2)^2
    Psi2 = 2 * g['c1W']
    Cres = (1 + g['Psi1']) ** 2 + (g['Psi1'] * (1 + 2 * g['c1W'].sqrt()) + Psi2) ** 2
    Kres = (arb(R1) / YSTAR) ** 2 * Cres / (4 * DELTA_I0) / lam
    K = w0 / (2 * g['muE']) + g['af'] * JW1 + g['bf'] * JB + Kres
    mg = g['Cg0'] * sq(JA)
    a = g['mu1'] - DELTA_I0 - lam / 2
    half = arb('0.5') - DELTA_I0
    JV0_L = 1 / (lam * g['muE'])
    if float(a.lower()) > 0:
        s = (mg + sq(mg * mg + 4 * a * K)) / (2 * a)
        JV0 = amin(JV0_L, s * s)
    else:
        JV0 = JV0_L
    JV2_I0 = (K + mg * sq(JV0) + posa(lam / 2 + DELTA_I0 - g['mu1']) * JV0) / half
    JV2 = amin(JV2_I0, 8 * B / g['muE'])
    sqX = c3 * sq(JD4) + c2 * sq(JW3) + c1 * sq(JB) + c0 * sq(JW1) + g['Cg2'] * sq(JA) + arb(RIs) / YSTAR / lam.sqrt()
    den = g['mu2'] - lam / 2
    assert float(den.lower()) > 0
    best = None
    for q in range(121):
        kap = arb(math.exp(-3 + 6 * q / 120))
        al = 1 + g['Psi1'] / kap
        bh = g['Psi1'] * kap + 2 * g['Psi1'] * g['c1W'].sqrt() + 2 * g['c1W']
        v = U(4 * (al * al + bh * bh / (4 * den)))
        best = v if best is None else min(best, v)
    JV4_I2 = 2 * w0 / (g['muE'] * g['epsI']) + arb(best) * sqX * sqX
    JV4 = amin(JV4_I2, 8 * B / (g['muE'] * g['epsI']))
    J = {0: JV0, 2: JV2, 4: JV4}
    J[1] = sq(J[0] * J[2]) + g['c1W'] * J[0]
    J[3] = sq(J[2] * J[4]) + g['c1W'] * J[2]
    return J, dict(B=U(B), JA=U(JA), JB=U(JB), JD4=U(JD4), JW1=U(JW1), JW3=U(JW3), K=U(K), Kres=U(Kres),
                   route_V0='I0' if U(JV0) < U(JV0_L) else 'Lambda', route_V2='I0' if U(JV2_I0) < U(8 * B / g['muE']) else 'Lambda',
                   route_V4='I2' if U(JV4_I2) < U(8 * B / (g['muE'] * g['epsI'])) else 'Lambda',
                   share_D4=U((c3 * sq(JD4)) ** 2 * arb(best) / JV4) if U(JV4_I2) < U(8 * B / (g['muE'] * g['epsI'])) else None)


T7CONST = os.path.join(HERE, '..', '..', 't7', 'constants.json')


def xr_per_slab(slabs, by6):
    """review B6.  T5 (M): |x| <= mbar/(1 - mbar e(A+W1)) (lbar A + c_phi R2 + n_phi eta1 (A+W1+B) + strip).  The linear part is
    T7's C_x Lambda^(1/2); the residual part x_r = mbar c_phi R2 (T7 C3 mbar per slab, the largest c_phi of the T7 certificates,
    T6 R2 per slab) is added here to the residual (B) parts of g_j and h_j.  The nonlinear and E_ij parts are O(eta_1 Ybar) relative;
    they fit in C_x (Ybar - Y_cert), Y_cert <= 6.3e-13 the certified (Y) level of p2/ycomp.py, for every Ybar >= YSTAR = 1e-12."""
    t7c = json.load(open(T7CONST))
    cphi = arb(repr(max(max(c['c_phi']) for c in t7c['certificates'])))
    return {s: arb(U(arb(repr(t7c['slabs'][str(s)]['C3']['mbar'])) * cphi * arb(repr(by6[s]['norms']['R2'])))) for s in slabs}


def integrate(C, F0=(0.0, 0.0, 0.0, 0.0), t7=None, jmax=3):
    t7 = t7 or C['t7']
    R6 = json.load(open(os.path.join(T6DIR, 'residuals.json')))
    by6 = {r['slab']: r for r in R6['per_slab']}
    XR = xr_per_slab([p['slab'] for p in C['per_slab']], by6)          # review B6: residual part of (M), per slab
    R1 = max(r['norms']['R1'] for r in R6['per_slab'] if C['slabs_range'][0] <= r['slab'] <= C['slabs_range'][1])
    RIs = max(r['norms']['RIs'] for r in R6['per_slab'] if C['slabs_range'][0] <= r['slab'] <= C['slabs_range'][1])
    per = C['per_slab']; n = len(per)
    slabs = [p['slab'] for p in per]
    assert slabs == list(range(slabs[0], slabs[0] + n)), 'slabs must be contiguous'
    assert all(abs(per[m + 1]['tau0'] - per[m]['tau0'] - DT) < 1e-9 for m in range(n - 1)), 'slab length must be DT'
    dt = arb('0.15')
    tau0 = arb(per[0]['tau0'])
    Cx = arb(t7['Cx'])
    lam = {j: [arb(p['lev'][str(j)]['lam']) for p in per] for j in range(jmax + 1)}
    # j = 0: the e-self-coupling of (S5), 2 sqrt(F0) * A2 sqrt(F0), has no lower level to go to: it lowers the rate by 2 A2
    lam[0] = [arb(L(lam[0][m] - 2 * arb(p['lev']['0']['A2']))) for m, p in enumerate(per)]
    for j in lam: assert all(float(v.lower()) > 0 for v in lam[j])
    aF = {}; bF = {}; diag = {}
    Jcache = {}
    for j in range(jmax + 1):
        DSj = arb(C['design']['dS'][j])
        Lend = []; acc = arb(0)
        for m in range(n): acc = acc + lam[j][m] * dt; Lend.append(acc)
        Lstart = [arb(0)] + Lend[:-1]
        lstar = []; cur = None
        for m in range(n):
            cur = lam[j][m] if cur is None else arb(min(L(cur), L(lam[j][m])))       # running min (lower bounds)
            lstar.append(cur)
        hA, hB, gA, gB = [], [], [], []
        for m in range(n):
            lv = per[m]['lev'][str(j)]
            c_ = {int(k): arb(v) for k, v in lv['c'].items()}; d_ = {int(k): arb(v) for k, v in lv['d'].items()}
            gA.append(sum((c_.get(i, arb(0)) / 2 * aF[i][m] for i in range(j)), arb(0)) + Cx * arb(lv['gam']) + arb(lv['A1']) / 2)
            xr = XR[per[m]['slab']]
            gB.append(sum((c_.get(i, arb(0)) / 2 * bF[i][m] for i in range(j)), arb(0)) + xr * arb(lv['gam']))
            modc = arb(lv['gr']) ** 2 / DSj + 8 / A_S3 * arb(lv['gha']) ** 2 + arb(lv['ghb']) ** 2 / NU
            mod = Cx * Cx * modc
            hA.append(sum(((1 + EPSS) * d_.get(i, arb(0)) * aF[i][m] ** 2 for i in range(j)), arb(0)) + arb(lv['eq']) + (1 + EPSS) * mod)
            # N3 couplings: (Y a + b)^2 <= (1 + s) Y^2 a^2 + (1 + 1/s) b^2 with s = b/(YSTAR a) clipped to [1e-3, 1e3] (valid for every Y)
            dn = {int(k_): arb(v) for k_, v in lv.get('dn3', {}).items()}
            sn = {i: arb(min(1e3, max(1e-3, U(bF[i][m]) / max(L(YSTAR * aF[i][m]), 1e-300)))) for i in dn}
            hA[-1] = hA[-1] + sum(((1 + sn[i]) * dn[i] * aF[i][m] ** 2 for i in dn), arb(0))
            n6 = by6[per[m]['slab']]['norms']
            res = arb(n6[f'RF{j}_opt4n']) ** 2 / DSj + 8 / A_S3 * arb(n6[f'RFa{j}_opt4n']) ** 2 + arb(n6[f'RFb{j}_opt4n']) ** 2 / NU
            res = (sq(res) + xr * sq(modc)) ** 2                    # (R + x_r g^r)-type terms, by Minkowski (review B6)
            hB.append(sum(((1 + 1 / EPSS) * d_.get(i, arb(0)) * bF[i][m] ** 2 for i in range(j)), arb(0)) + (1 + 1 / EPSS) * res
                      + arb(lv.get('hf_abs', 0.0)) + sum(((1 + 1 / sn[i]) * dn[i] * bF[i][m] ** 2 for i in dn), arb(0)))
        aF[j] = []; bF[j] = []; dg = dict(infl=[], HdA=[], GA=[], HB=[], GB=[], lstar=[], J=[])
        Kmax = None; rmax = {}; remax = 0.0; Klow = {}
        for m in range(n):
            lv = per[m]['lev'][str(j)]
            Kmax = lv['k'][str(j)] if Kmax is None else max(Kmax, lv['k'][str(j)])
            for i, v in lv['r'].items(): rmax[int(i)] = max(rmax.get(int(i), 0.0), v)
            remax = max(remax, lv['re'])
            key = (L(lstar[m]), m)
            if key not in Jcache:
                Jcache[key] = cascade(lstar[m], arb(per[m]['tau0']) - tau0, t7, R1, RIs)
            J, Jd = Jcache[key]
            for mm, v in lv['k'].items():
                if int(mm) != j: Klow[int(mm)] = max(Klow.get(int(mm), 0.0), v)
            infl = arb(Kmax) * J[j] + sum((arb(rmax.get(i, 0.0)) * sq(J[j] * J[i]) for i in range(j)), arb(0)) + arb(remax) * sq(J[j] / lstar[m]) \
                + sum((arb(v) * J[mm] for mm, v in Klow.items()), arb(0))
            HdA = arb(0); HdB = arb(0); GA = arb(0); GB = arb(0)
            for q in range(m):
                dec = (-(Lstart[m] - Lend[q])).exp(); dech = (-(Lstart[m] - Lend[q]) / 2).exp()
                HdA += dec * dt * hA[q]; HdB += dec * dt * hB[q]; GA += dech * dt * gA[q]; GB += dech * dt * gB[q]
            HdA += dt * hA[m]; HdB += dt * hB[m]; GA += dt * gA[m]; GB += dt * gB[m]
            HB = HdB + (-Lstart[m]).exp() * arb(F0[j])
            aF[j].append(sq(infl + HdA) + GA)
            bF[j].append(sq(HB) + GB)
            dg['infl'].append(U(infl)); dg['HdA'].append(U(HdA)); dg['GA'].append(U(GA)); dg['HB'].append(U(HB)); dg['GB'].append(U(GB))
            dg['lstar'].append(L(lstar[m])); dg['J'].append({k: U(v) for k, v in J.items()})
        dg['Kmax'] = Kmax; dg['rmax'] = rmax; dg['remax'] = remax; dg['hA'] = hA; dg['hB'] = hB; dg['gA'] = gA; dg['gB'] = gB
        dg['Lend'] = Lend; dg['Lstart'] = Lstart; dg['Jd'] = Jcache[(L(lstar[-1]), n - 1)][1]
        diag[j] = dg
    return dict(aF=aF, bF=bF, diag=diag, slabs=slabs, R1=R1, RIs=RIs)


def handoff(C, R, window=None):
    per = C['per_slab']; slabs = R['slabs']
    window = window or WINDOW or ((C['per_slab'][-1]['slab'],) * 2 if SUBW else window_of(C))
    idx = [slabs.index(s) for s in range(window[0], window[1] + 1)]
    aF, bF = R['aF'], R['bF']
    A = {j: amax(*[aF[j][m] for m in idx]) for j in aF}
    Bv = {j: amax(*[bF[j][m] for m in idx]) for j in bF}
    out = dict(window_slabs=list(window), window_k_start=[per[idx[0]]['k'], per[idx[-1]]['k']], A={j: U(v) for j, v in A.items()},
               B={j: U(v) for j, v in Bv.items()}, lim={})
    # item 3: KCH (sqrt F3 + 3 sqrt F2 + sqrt F1) <= sqrt(1e-5) - sqrt(E3_base)
    for jj in range(4):
        amp_ = KCH * sum((s * A[mm] for mm, s in STIR[jj].items()), arb(0)) + arb(repr(C0J.get(jj, 0.0))) * A[0]
        bres = KCH * sum((s * Bv[mm] for mm, s in STIR[jj].items()), arb(0)) + arb(repr(C0J.get(jj, 0.0))) * Bv[0]
        allow = arb(THR['E'][jj]).sqrt() - arb(BASE_T1['E'][jj]).sqrt()
        out['lim'][f'item3_E{jj}'] = L((allow - bres) / amp_)
    for jj in range(3):
        amp_ = KCH * sum((s * A[mm] for mm, s in STIR[jj].items()), arb(0)) + arb(repr(C0J.get(jj, 0.0))) * A[0]
        bres = KCH * sum((s * Bv[mm] for mm, s in STIR[jj].items()), arb(0)) + arb(repr(C0J.get(jj, 0.0))) * Bv[0]
        allow = arb(THR['shell'][jj]).sqrt() - arb(BASE_T1['shell'][jj]).sqrt()
        out['lim'][f'item2_E{jj}'] = L((allow - bres) / amp_)
    # item 4 at a good time tau* in the window (C5): D2(tau*) <= (1/|W|) int_W D2,
    #   theta2 int_W D2 <= F2(t_a) + int_W (2 g2 sqrt F2 + h2) + e^{lam* |W|} Ybar^2 inflow_2(t_b)
    dt = arb(SUBW) if SUBW else arb('0.15'); Wl = dt * len(idx)
    d2 = R['diag'][2]
    th2 = arb(min(per[m]['lev']['2']['Dcoef'] for m in idx))
    lst = arb(d2['lstar'][idx[-1]])
    inflW = (lst * Wl).exp() * arb(d2['infl'][idx[-1]])
    def D2bound(Y):
        F2a = (Y * aF[2][idx[0]] + bF[2][idx[0]]) ** 2
        s = arb(0)
        for m in idx:
            s += dt * (2 * (Y * d2['gA'][m] + d2['gB'][m]) * (Y * aF[2][m] + bF[2][m]) + Y * Y * d2['hA'][m] + d2['hB'][m])
        return (F2a + s + Y * Y * inflW) / (th2 * Wl)
    ya, Il, Lc = arb(50), arb(25), arb(50)
    allow4 = arb(THR['K4']) - arb(BASE_T1['K4'])
    def K4(Y):
        F3 = (Y * A[3] + Bv[3]) ** 2
        D2 = D2bound(Y)
        return KCH * Lc ** 3 * sq(F3 / (28 * ya ** 4 * Il) + 2 * sq(F3) * sq(D2) / (28 * ya ** 3))
    lo, hi = 0.0, 1e-8
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if U(K4(arb(mid))) <= float(allow4.lower()): lo = mid
        else: hi = mid
    assert U(K4(arb(lo))) <= float(allow4.lower())
    out['lim']['item4_K4_goodtime'] = lo
    out['D2_goodtime_per_Ybar2'] = U(D2bound(arb(1)) - D2bound(arb(0)))
    YB = arb(os.environ.get('YBAR', '1e-12'))                  # the certified majorant (p2/ycomp.py), for T9 item 5 / boxes
    out['at_Ybar'] = dict(Ybar=U(YB), D2_goodtime=U(D2bound(YB)), K4=U(K4(YB)), sqrtF=[U(YB * A[j] + Bv[j]) for j in range(4)])
    # review A1: the parts of the good-time bound (it is quadratic in Ybar): the Ybar-free (residual) part D2(0), the part
    # linear in Ybar (cross terms of (Ybar a + b)^2) and the Ybar^2 part, at YB; and the window length used
    d0, d1, dm = D2bound(arb(0)), D2bound(YB), D2bound(-YB)
    out['D2_parts'] = dict(window_length=U(Wl), zero=U(d0), linear=U((d1 - dm) / 2), quadratic=U((d1 + dm) / 2 - d0), theta2=L(th2))
    out['Ymax'] = min(out['lim'].values())
    out['binding'] = min(out['lim'], key=out['lim'].get)
    # decay of the initial data from the first slab to the window: sqrt F_j(tau) gets e^{-Lambda_j/2} sqrt F_j(tau0)
    out['initial_data_factor'] = {j: U((-R['diag'][j]['Lstart'][idx[0]] / 2).exp()) for j in R['diag']}
    # breakdown of A_3 at the window end
    m = idx[-1]; dj = R['diag'][3]
    out['A3_parts'] = dict(sqrt_inflow=U(sq(arb(dj['infl'][m]))), sqrt_duhamel=U(sq(arb(dj['HdA'][m]))), G=dj['GA'][m])
    return out


def apply_n3(C, path):
    """Lemma N3' (LEMMA_F.md 6.1; n3far.py): add the certified nonlinear and modulation increments, level by level and slab by slab."""
    N3 = json.load(open(path))['per_slab']
    for p in C['per_slab']:
        o = N3.get(str(p['slab']))
        if o is None: raise KeyError(f"n3 file has no slab {p['slab']}")
        levs = o.get('lev') or {'3': dict(dtheta=o['dtheta3'], dsigma=o['dsigma3'], dk=o['dk33'], d=o['d3'], k=o['k3'], c={}, hf=o['hf3'])}
        for js, inc in levs.items():
            if js not in p['lev']: continue
            lv = p['lev'][js]; j = int(js)
            lv['lam'] = L(arb(lv['lam']) - arb(inc['dsigma'])); assert lv['lam'] > 0, (p['slab'], j, lv['lam'])
            lv['Dcoef'] = L(arb(lv['Dcoef']) - arb(inc['dtheta'])); assert lv['Dcoef'] > 0, (p['slab'], j, lv['Dcoef'])
            lv['k'][js] = U(arb(lv['k'][js]) + arb(inc['dk']))
            for m in range(j):
                ms = str(m)
                lv['k'][ms] = U(arb(lv['k'].get(ms, 0.0)) + arb(inc['k'].get(ms, 0.0)))
                lv['c'][ms] = U(arb(lv['c'].get(ms, 0.0)) + arb(inc['c'].get(ms, 0.0)))
                lv.setdefault('dn3', {})[ms] = inc['d'].get(ms, 0.0)   # N3 couplings: adaptive Ybar/residual split in integrate()
            lv['hf_abs'] = U(arb(lv.get('hf_abs', 0.0)) + arb(inc['hf']))
    return C

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--consts', default=os.path.join(HERE, 'constants', 'constants.json'))
    ap.add_argument('--out', default=os.path.join(HERE, 'constants', 'amp.json'))
    ap.add_argument('--F0', default='0,0,0,0', help='F_j(tau0) at the start of the first slab (v2: slab 0; T4 input); default 0')
    ap.add_argument('--window', default=None, help='a-b: good-time window slabs (default: trailing slabs with k >= 2.0e4)')
    ap.add_argument('--Cg0', type=float, default=None, help='override T7 C_gamma0 (sensitivity)')
    ap.add_argument('--n3', default=None, help="Lemma N3' increments from n3far.py (level-3 nonlinear terms on the far boxes)")
    a = ap.parse_args()
    C = json.load(open(a.consts))
    if a.n3: C = apply_n3(C, a.n3)
    t7 = dict(C['t7'])
    if a.Cg0 is not None: t7['Cg0'] = a.Cg0
    F0 = [float(v) for v in a.F0.split(',')]
    R = integrate(C, F0, t7)
    win = tuple(int(v) for v in a.window.split('-')) if a.window else None
    h = handoff(C, R, win)
    Jd = R['diag'][3]['Jd']; Jend = R['diag'][3]['J'][-1]
    out = dict(inputs=dict(t7={k: v for k, v in t7.items() if not k.startswith('_')}, t7_sources=t7.get('_sources'), F0=F0,
                           R1=R['R1'], RIs=R['RIs'], YSTAR=U(YSTAR), KCH=U(KCH), THR=THR, BASE_T1=BASE_T1, window=h['window_slabs']),
               handoff=h, J_end=Jend, J_details=Jd,
               per_slab={j: dict(aF=[U(v) for v in R['aF'][j]], bF=[U(v) for v in R['bF'][j]], lstar=R['diag'][j]['lstar'],
                                 infl=R['diag'][j]['infl'], GA=R['diag'][j]['GA']) for j in R['aF']},
               slabs=R['slabs'])
    tmp = a.out + '.tmp'
    with open(tmp, 'w') as f: json.dump(out, f, indent=1)
    os.replace(tmp, a.out)
    print('A_j (window sup, per unit Ybar):', {j: f'{v:.4e}' for j, v in h['A'].items()}, ' residual parts B_j:', {j: f'{v:.2e}' for j, v in h['B'].items()})
    print('limits:', {k: f'{v:.3e}' for k, v in h['lim'].items()})
    print(f"Ybar requirement {h['Ymax']:.3e} (binding {h['binding']}); J at the end: " + ', '.join(f'{Jend[k]:.3e}' for k in sorted(Jend)) +
          f"; routes {Jd['route_V0']}/{Jd['route_V2']}/{Jd['route_V4']}; A3 parts {h['A3_parts']}")


if __name__ == '__main__':
    main()
