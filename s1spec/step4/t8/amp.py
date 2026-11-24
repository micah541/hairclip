"""Integrate the triangular comparison system of Lemma F along the T1 base (LEMMA_F.md, Corollary C) per unit Ybar and derive
the Ybar the hand-off items 2-4 need.   usage: python3 amp.py [cfg ...]   (reads consts_<cfg>.json)
  sqrt(F_j) <= sqrt(H_j) + G_j,   H_j = e^{-Lam_j} F_j(tau0) + Inflow_j + int e^{-(Lam_j(t)-Lam_j(s))} (sum_i d_ji F_i + e_j + 2 R_j^2)
                                  G_j = int e^{-(Lam_j(t)-Lam_j(s))/2} (sum_i c_ji/2 sqrt(F_i) + Cx gam_j + a_j)
  Inflow_j (W-ramp data, Lemma J):   sum_m k_jm J_m + sum_i r_ji sqrt(J_j J_i) + re_j sqrt(J_j/lam)
  Inflow_j (nested, ramp in the core of level j-2):  kD_j iD_{j-2} + sum_i sD_ji sqrt(iD_{j-2} F_i/lam) + ...
  iD_s := int omega D_s <= (H_s(lam) + 2 sup(g_s sqrt F_s)/lam) / th_s."""
import json, sys, numpy as np
from math import comb
import jcas

BASE_T1 = dict(E=(1.45e-8, 5.2e-8, 1.7e-7, 4.6e-7), K4=2.1e-4, shell=(8.0e-11, 3.4e-10, 1.4e-9))   # T1 report, k = 2.3e4
THR = dict(E=(2e-7, 1e-5, 1e-5, 1e-5), K4=0.06, shell=(2e-7, 1e-4, 0.2))
KCH = 1.02           # chart factor (T5 2.9; T9 to certify)
STIR = {0: {0: 1}, 1: {1: 1}, 2: {1: 1, 2: 1}, 3: {1: 1, 2: 3, 3: 1}}

def load(cfg):
    d = json.load(open(f'/home/claude/s1spec/step4/t8/consts_{cfg}.json'))
    rows = d['rows']; tau1 = d['tau1']
    import fb
    B = fb.base()
    slabs = sorted(set(r['slab'] for r in rows))
    per = []
    for i in slabs:
        rs = [r for r in rows if r['slab'] == i]
        t0 = float(B.tau0[i]); t1 = t0 + float(B.dt[i])
        lev = {}
        for j in range(5):
            L = [r['lev'][str(j)] for r in rs]
            def mx(key, sub=None):
                vals = []
                for l in L:
                    v = l.get(key, 0.0)
                    if sub is not None: v = v.get(str(sub), 0.0) if isinstance(v, dict) else 0.0
                    vals.append(v)
                return max(vals)
            lev[j] = dict(lam=min(l['lam'] for l in L), nested=L[0]['nested'],
                          k={m: mx('k', m) for m in range(j + 1)}, kD=mx('kD'), Dcoef=min(l['Dcoef'] for l in L),
                          c={i_: mx('cji', i_) for i_ in range(j)}, cf={i_: mx('cji_far', i_) for i_ in range(j)},
                          d={i_: mx('dji', i_) for i_ in range(j)}, r={i_: mx('rji', i_) for i_ in range(j)},
                          sD={i_: mx('sDji', i_) for i_ in range(j)}, A1=mx('A1'), re=mx('ramp_e'), eq=mx('e_quad'),
                          gam=mx('gam'), gamc=mx('gam_cone'), gr=mx('gam_rough'), kat=max(l.get('k_at_y', 0.0) for l in L))
        per.append(dict(i=i, t0=t0, t1=t1, k=rs[0]['k'], lev=lev))
    return per, tau1

def integrate(per, tau1, P=None, Cx=4.0, gam_key='gam', D4pt=None, Wscale=None, res=None, jmax=4, tau_start=None, far=None):
    """returns sqrtF[j] arrays on the slab-end grid, and diagnostics at tau1 (unit Ybar; residual R_j absolute if given)"""
    Pj = dict(jcas.T7); Pj.update(P or {}); Pj['Cx'] = Cx
    if tau_start is not None: per = [p for p in per if p['t1'] > tau_start]
    n = len(per); te = np.array([p['t1'] for p in per]); dt = np.array([p['t1'] - p['t0'] for p in per])
    tau0 = per[0]['t0']
    sq = {}; Hh = {}; Gg = {}; diag = {}
    iD = {}
    lamall = {j: np.array([p['lev'][j]['lam'] for p in per]) for j in range(jmax + 1)}
    for j in range(jmax + 1):
        lam = lamall[j].copy()
        for jj in range(j + 1, jmax + 1):          # nested consumer of this level: integrate with the slower rate
            if per[0]['lev'][jj]['nested'] == j: lam = np.minimum(lam, lamall[jj])
        Lam = np.cumsum(lam * dt)                                       # Lam at slab ends
        lam_min = np.minimum.accumulate(lam)
        H = np.zeros(n); G = np.zeros(n); infl = np.zeros(n)
        Jc = [jcas.cascade(lam_min[m], te[m] - tau0, Pj, D4_pointwise=D4pt)['J'] for m in range(n)]
        for m in range(n):
            lv = per[m]['lev'][j]
            # running sups of the ramp constants
            ks = {mm: max(per[q]['lev'][j]['k'][mm] for q in range(m + 1)) for mm in range(j + 1)}
            if Wscale is not None: ks = {mm: v * Wscale(j, mm, per[m]['lev'][j]['kat']) for mm, v in ks.items()}
            J = Jc[m]
            if lv['nested'] is None:
                rs = {i_: max(per[q]['lev'][j]['r'][i_] for q in range(m + 1)) for i_ in range(j)}
                re = max(per[q]['lev'][j]['re'] for q in range(m + 1))
                inf_ = sum(ks[mm] * J[mm] for mm in ks) + sum(rs[i_] * np.sqrt(J[j] * J[i_]) for i_ in range(j)) + re * np.sqrt(J[j] / lam_min[m])
            else:
                s_ = lv['nested']
                kD = max(per[q]['lev'][j]['kD'] for q in range(m + 1))
                sD = {i_: max(per[q]['lev'][j]['sD'][i_] for q in range(m + 1)) for i_ in range(j)}
                re = max(per[q]['lev'][j]['re'] for q in range(m + 1))
                Fi_sup = {i_: np.max(sq[i_][:m + 1])**2 for i_ in range(j)}
                inf_ = kD * iD[s_][m] + sum(sD[i_] * np.sqrt(iD[s_][m] * Fi_sup[i_] / lam_min[m]) for i_ in range(j)) \
                       + re * np.sqrt(iD[s_][m] / lam_min[m]) + sum(ks.get(mm, 0.0) * J[mm] for mm in ks if mm < j)
            infl[m] = inf_
        # Duhamel parts
        hd = np.zeros(n); gd = np.zeros(n)
        hsrc = np.zeros(n); gsrc = np.zeros(n)
        for m in range(n):
            lv = per[m]['lev'][j]
            Fi = {i_: max(sq[i_][m], sq[i_][m - 1] if m else sq[i_][0])**2 for i_ in range(j)}
            hsrc[m] = sum(lv['d'][i_] * Fi[i_] for i_ in range(j)) + lv['eq'] + (2 * (res[j] + Cx * lv['gr'])**2 if res is not None else 2 * (Cx * lv['gr'])**2)
            if far is not None:            # localized energies: far-region couplings through the global energies (factor Om_far)
                Rg, Omf = far
                hsrc[m] += sum(lv['cf'][i_] * Omf * Rg['sq'][j][m] * Rg['sq'][i_][m] for i_ in range(j))
            gsrc[m] = sum(lv['c'][i_] / 2 * np.sqrt(Fi[i_]) for i_ in range(j)) + Cx * lv[gam_key] + lv['A1'] / 2
        for m in range(n):
            w = np.exp(-(Lam[m] - Lam[:m + 1]))          # decay from the end of slab q to the end of slab m
            hd[m] = np.sum(w * dt[:m + 1] * hsrc[:m + 1])
            gd[m] = np.sum(np.sqrt(w) * dt[:m + 1] * gsrc[:m + 1])
        H = infl + hd; G = gd
        sq[j] = np.sqrt(H) + G
        Hh[j], Gg[j] = H, G
        # time-weighted dissipation integral of this level (for nested sources), weight rate = lam of the next levels (~ same)
        th = min(p['lev'][j]['Dcoef'] for p in per)
        iD[j] = (H + 2 * np.array([np.max(gsrc[:m + 1]) for m in range(n)]) * sq[j] / np.maximum(lam_min, 1e-3)) / max(th, 1e-3)
        diag[j] = dict(infl=infl, hd=hd, gd=gd, lam_min=lam_min, Dcoef=th, J=Jc, hsrc=hsrc, gsrc=gsrc)
    return dict(sq=sq, te=te, diag=diag, H=Hh, G=Gg, iD=iD, n=n)

def handoff(R, per, tau1, window_dt=0.6, Cx=4.0):
    """Ybar allowed by items 2-4 at the hand-off (tau1 = k 2.3e4), unit-Ybar amplifications A_j = sqrt(F_j(tau1))"""
    te = R['te']; m1 = int(np.argmin(np.abs(te - tau1)))
    A = {j: float(R['sq'][j][m1]) for j in R['sq']}
    out = dict(A=A)
    # item 3 windows (deviation energies, omega-normalised): sqrt(E_j^err) <= KCH * sum_m S(j,m) sqrt(F_m)
    lim = {}
    for j in range(4):
        amp = KCH * sum(Sjm * A[m] for m, Sjm in STIR[j].items())
        allow = np.sqrt(THR['E'][j]) - np.sqrt(BASE_T1['E'][j])
        lim[f'item3_E{j}'] = allow / amp
    for j in range(3):
        amp = KCH * sum(Sjm * A[m] for m, Sjm in STIR[j].items())
        lim[f'item2_E{j}'] = (np.sqrt(THR['shell'][j]) - np.sqrt(BASE_T1['shell'][j])) / amp
    # item 4: K4 = L^3 sup|w3| on y in [50,350], L = 50; window I = [y, 1.5 y], y = 50
    ya, Ilen, L = 50.0, 25.0, 50.0
    allow4 = THR['K4'] - BASE_T1['K4']
    if 4 in A:
        s2 = KCH**2 * L**6 * (A[3]**2 / (28 * ya**4 * Ilen) + 2 * A[3] * A[4] / (28 * ya**5))
        lim['item4_K4_viaF4'] = allow4 / np.sqrt(s2)
    # good time: D_2(tau*) <= (1/Dt) int_{tau1-Dt}^{tau1} D_2 <= e^{lam Dt} iD_2(tau1) / Dt
    lam2 = R['diag'][2]['lam_min'][m1]
    D2s = np.exp(lam2 * window_dt) * R['iD'][2][m1] / window_dt
    s2 = KCH**2 * L**6 * (A[3]**2 / (28 * ya**4 * Ilen) + 2 * A[3] * np.sqrt(D2s) / (28 * ya**3))
    lim['item4_K4_goodtime'] = allow4 / np.sqrt(s2)
    out['D2_goodtime'] = float(D2s)
    out['lim'] = lim
    out['Ymax'] = min(v for k, v in lim.items() if not (k == 'item4_K4_viaF4' and 'nolevel4' in out))
    out['Ymax_goodtime'] = min(v for k, v in lim.items() if k != 'item4_K4_viaF4')
    out['binding'] = min(lim, key=lim.get)
    return out

if __name__ == '__main__':
    cfgs = sys.argv[1:] or ['T5']
    for cfg in cfgs:
        per, tau1 = load(cfg)
        for gk in ('gam', 'gamc'):
            R = integrate(per, tau1, gam_key=gk)
            h = handoff(R, per, tau1)
            print(f"{cfg} [{gk}]: A_j = " + ', '.join(f"{h['A'][j]:.3e}" for j in range(5)),
                  f"| Ymax {h['Ymax']:.2e} (binding {h['binding']}), with good time {h['Ymax_goodtime']:.2e}")
            print('     limits: ' + ', '.join(f"{k}={v:.2e}" for k, v in h['lim'].items()))

def breakdown(R, per, tau1, j, Cx=4.0, gam_key='gam'):
    """shares of sqrt(F_j(tau1)) (unit Ybar): inflow pieces, Duhamel couplings, G pieces"""
    te = R['te']; m1 = int(np.argmin(np.abs(te - tau1))); n = m1 + 1
    lv = per[m1]['lev'][j]; J = R['diag'][j]['J'][m1]; lam = R['diag'][j]['lam_min'][m1]
    ks = {mm: max(per[q]['lev'][j]['k'][mm] for q in range(n)) for mm in range(j + 1)}
    out = {}
    if lv['nested'] is None:
        out['inflow k_jj J_j'] = ks[j] * J[j]
        out['inflow k_jm J_m (m<j)'] = sum(ks[mm] * J[mm] for mm in ks if mm < j)
        rs = {i_: max(per[q]['lev'][j]['r'][i_] for q in range(n)) for i_ in range(j)}
        out['ramp cross r_ji sqrt(JjJi)'] = sum(rs[i_] * np.sqrt(J[j] * J[i_]) for i_ in range(j))
        re = max(per[q]['lev'][j]['re'] for q in range(n)); out['ramp e-terms'] = re * np.sqrt(J[j] / lam)
    else:
        s_ = lv['nested']; kD = max(per[q]['lev'][j]['kD'] for q in range(n))
        out[f'nested kD iD_{s_}'] = kD * R['iD'][s_][m1]
    out['Duhamel (d_ji F_i, e, rough mod.)'] = R['diag'][j]['hd'][m1]
    H = R['diag'][j]['infl'][m1] + R['diag'][j]['hd'][m1]
    sqH = np.sqrt(H)
    shares = {k: v / H * sqH for k, v in out.items()}       # share of sqrt(H), proportional
    G = R['G'][j][m1]
    # G pieces: couplings / modulation / e-forcing (ratio of their time-averaged sources)
    gs = R['diag'][j]['gsrc'][:n]
    cpl = np.array([sum(per[q]['lev'][j]['c'][i_] / 2 * R['sq'][i_][q] for i_ in range(j)) for q in range(n)])
    mod = np.array([Cx * per[q]['lev'][j][gam_key] for q in range(n)]); ee = np.array([per[q]['lev'][j]['A1'] / 2 for q in range(n)])
    tot = cpl + mod + ee; tot[tot == 0] = 1
    Lam = np.cumsum(np.array([per[q]['lev'][j]['lam'] for q in range(n)]) * np.array([per[q]['t1'] - per[q]['t0'] for q in range(n)]))
    w = np.sqrt(np.exp(-(Lam[-1] - Lam))) * np.array([per[q]['t1'] - per[q]['t0'] for q in range(n)])
    shares['G: couplings c_ji'] = np.sum(w * cpl); shares['G: modulation |x| gam'] = np.sum(w * mod); shares['G: e-forcing'] = np.sum(w * ee)
    shares['total sqrt(F)'] = R['sq'][j][m1]
    return shares, J
