#!/usr/bin/env python3
"""Lemma U, Phase 2a (agent B): [est] split of the z_{j+1} ramp column K_{j,j+1} into the top pair (j+1,j+1) and the mixed pairs, and the
size of the IBP + Young route for the top pair (U_ramp.tex Lemma Y).  Float grid on arb-evaluated coefficient sums; NOT a certified bound.

  PYTHONHASHSEED=0 nice -n 19 python3 -B top_pair_est.py      (single process, < 150 MB, < 30 s)

Coefficients: C_pair[m] = sum_{(m,i,r)} |coef| Ab_(i,r) (P2's ramp convention: Ab = Z-free part of At, box (0.11, 0.52), U_ramp = 1/80),
rem_m majorant |rem_m| <= sum_l C(m,l) |(eta^2)^(l)| RT_{m-l},  RT = sup|Theta^(k)/Theta| on the ramp = (1, 3.8816, 15.5275, 62.1263) [cert, P2 geom].
Young route (level j, top pair): int c eps Theta z_{j+1}^2 = - int c eps Theta z_j z_{j+2} + (1/2) int (c eps Theta)'' z_j^2
  <= delta a_min int eps Theta eta^2 z_{j+2}^2 + int [ c^2/(4 delta a_min eta^2) + (1/2)|(c eps Theta)''|/(eps Theta) ] eps Theta z_j^2,
with |(c eps Theta)''|/(eps Theta) <= |c''| + 2 Lmax |c'| + L2max |c|, Lmax = 3.5213, L2max = Lmax^2 + 16 (estimate), c', c'' of the majorant
(the A^{(r+u)} parts are approximated by scaling the majorant; the jet parts of A^{(r+u)} -- new cubic classes -- are NOT included).
Output: top_pair_est.json (args block), stdout -> top_pair_est.out."""
import sys, os, json, math, hashlib, time, platform
sys.dont_write_bytecode = True
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); RAMP = '/home/claude/s1spec/writeup/lemmaU/phase1/ramp'
sys.path.insert(0, RAMP)
import sympy as sp
import ulib as U
Q_ = sp.Rational
T = U.load_tables()
At = U.eval_At(U.build_At_sym(Q_(1, 80)), Q_(11, 100), Q_(52, 100)); Ab = U.ab_from_At(At)
RT = [1.0, 3.8815925399248954, 15.52750011955708, 62.12633722823463]
LN2 = math.log(2); Lmax = 3.5212668029853376; L2max = Lmax**2 + 16
t = np.linspace(0, 1, 4001); s = math.log(28) + t*LN2; y = np.exp(s)
P = np.polynomial.Polynomial; Sp_ = P([0, 0, 0, 10, -15, 6]); E2 = Sp_**2
de2 = [E2.deriv(k)(t)/LN2**k for k in range(6)]; eta = Sp_(t)
def remmaj(m, extra=0):
    """majorant of |d^extra/ds^extra rem_m| (extra derivatives fall on the eta part only, plus the RT growth by Leibniz, crude)"""
    tot = np.zeros_like(t)
    for l_ in range(1, m + 1):
        for e in range(extra + 1):
            tot += math.comb(m, l_)*math.comb(extra, e)*np.abs(de2[l_ + e])*RT[min(m - l_ + extra - e, 3)]*(4.0**max(0, m - l_ + extra - e - 3))
    return tot
res = dict(levels={})
for j in (1, 2):
    n = j + 2; k = j + 1; parts = {}
    for (a, b), mons in T[j].items():
        if (a, b) == (n, n) or k not in (a, b): continue
        w = 1.0 if a == b else 0.5
        Cm = {}
        for (m, i, r), co in mons.items():
            if m == 0: continue
            Cm[m] = Cm.get(m, 0.0) + float(co)*float(Ab[(i, r)].mid())*w
        if not Cm: continue
        prof = sum(c*remmaj(m) for m, c in Cm.items())
        parts[(a, b)] = dict(C=Cm, sup=float(prof.max()), prof=prof)
    tot = sum(v['sup'] for v in parts.values())
    top = parts.get((k, k))
    out = dict(column_sup_sum=tot, pairs={str(kk): dict(C={str(m): c for m, c in v['C'].items()}, sup=v['sup']) for kk, v in parts.items()})
    if top is not None: out['f_top'] = top['sup']/tot
    # mixed pairs (a, k), a < k: exact IBP to energies (U_ramp.tex Lemma Y(a)); data coefficient on z_j..z_0 (all put on z_j: estimate)
    mib = 0.0
    for kk, v in parts.items():
        if kk == (k, k): continue
        cprof = v['prof']; cder = sum(cc*remmaj(m, 1) for m, cc in v['C'].items())
        mib += float((0.5*cder).max()) + (float(cprof.max()) if kk[0] < j else 0.0)
    out['mixed_ibp_coef'] = mib
    print(f'level {j}, column k = {k}: per-pair sups of the majorant (w included): ' + ', '.join(f'{kk}: {v["sup"]:.1f}' for kk, v in parts.items()) + f'  (sum {tot:.1f}; P2 cert K_{j}{k} = 532/537 pair part)')
    if top is not None:
        c = top['prof']; Cm = top['C']
        c1 = sum(cc*remmaj(m, 1) for m, cc in Cm.items()); c2 = sum(cc*remmaj(m, 2) for m, cc in Cm.items())
        amin = (1 + 0.15**2)**-2
        for delta in (0.1, 0.25, 0.5):
            with np.errstate(divide='ignore', invalid='ignore'):
                young = np.where(eta > 1e-12, c**2/(4*delta*amin*np.maximum(eta, 1e-300)**2), 0.0)
            ibp = 0.5*(c2 + 2*Lmax*c1 + L2max*c)
            Yt = young + ibp
            out[f'delta={delta}'] = dict(sup_young=float(young.max()), sup_ibp=float(ibp.max()), sup_total=float(Yt.max()))
            print(f'   top pair ({k},{k}): sup c = {c.max():.1f};  Young route, delta = {delta}: data coefficient on z_{j} (per eps Theta): sup {Yt.max():.4g} '
                  f'(c^2/(4 delta a_min eta^2) {young.max():.4g}, IBP remainder {ibp.max():.4g});  Ghat cost: B_c 1.5 -> {1.5 - delta} a_min')
            # late centre, design data56 (g = 3.503e5), pointwise D: swept data vs the X route for the top pair alone
            g = 3.503e5; D = 1 + g/y**4; l15 = math.log(1.5); xg = np.linspace(0, 16, 3201)
            Wf = lambda x: 1/(np.exp(-4*(x + 0.3)) + 1 + np.exp(4*(x - l15 - 0.3)))
            sig = np.log(y/56)[None, :] - xg[:, None]
            def sweep(f): sup = (f[None, :]*Wf(sig)).max(axis=1); return float(np.sum(0.5*(sup[1:] + sup[:-1]))*(xg[1] - xg[0]))/0.25*1.005
            SY = sweep(Yt*(28/y)**3*D); SX = float((D*y/28*Wf(np.log(y/56))).max())
            out[f'delta={delta}'].update(sweep_Y=SY, SX=SX)
            EA = (1.708e-07, 9.237e-05, 0.06006); D1 = 2*EA[1]/(1.25*(16/17)**2)
            if j == 2:
                dY = 28.**-4*SY*EA[2]; dX = float(c.max())*SX*D1
                out[f'delta={delta}'].update(data_Y_default=dY, data_X_top_default=dX)
                print(f'      late centre, default inputs: Young-route data {dY:.4g}  vs  X-route (top pair only) {dX:.4g}  [est]')
    res['levels'][str(j)] = out
res['args'] = dict(argv=sys.argv, cwd=os.getcwd(), python=sys.version, platform=platform.platform(),
                   env={k: os.environ.get(k) for k in ('PYTHONHASHSEED', 'PYTHONDONTWRITEBYTECODE')},
                   inputs={p: hashlib.sha256(open(p, 'rb').read()).hexdigest() for p in (U.TABLES, os.path.join(RAMP, 'ulib.py'), os.path.abspath(__file__))})
json.dump(res, open(os.path.join(HERE, 'top_pair_est.json'), 'w'), indent=1, default=float)
