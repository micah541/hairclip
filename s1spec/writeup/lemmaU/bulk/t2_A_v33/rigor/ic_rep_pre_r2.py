#!/usr/bin/env python3
"""Task 2 / A4: the initial datum of (A) (STATEMENTS.tex Remark 15, change C7).
Audit: step 4 certifies (H2) on the shells [14, 60] only (STEP4_RESULT.md l.33, camp item 2; gonogo/t7_step4_H2_tight), so
   E^A_j(tau_1) <= H2_j + 5 H3_j + (28/379.5) n_III1_j        (five (H3) windows on [50, 379.7]; zone-III part on [379.5, 700] from task 5).
Step-4 window energies at the certified Ybar = 1e-12, reproducing camp.py item 3 (step4/p2/t8se/camp.py l.229-234) from the same inputs (read-only):
   sqrt E_j^win <= sqrt BASE_j + amp_j Ybar + bres_j,  amp_j = KCH sum_m STIR[j][m] A_m + c0_j A_0 + c0tip_j,  bres_j = KCH sum STIR B_m + c0_j B_0
   (self-check: the item-3 Ybar limits at THR = (2e-7, 1e-5, 1e-5, 1e-5) reproduce amp_A1_w180.json to 1 ulp, on the safe side: A_m, B_m there are rounded-up floats).
H3 settings: 'h3_v32' = (2e-7, 1e-5, 1e-5); 'h3_x3' = levels 0, 1 divided by 3; 'h3_min' = 1.5 x E^win(1e-12), rounded up to 3 digits.
Output: out/ic_rep.json.   Run: PYTHONHASHSEED=0 nice -n 19 python3 -B rigor/ic_rep.py --zone3 PATH"""
import os, sys, json, argparse, math
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import arlib as L
from arlib import A, Q_, arb
P4 = '/home/claude/s1spec/step4/p2/'
AMPF, T9F, WSF = P4 + 't8se/constants/amp_A1_w180.json', P4 + 'out/t9_cone_kch.json', P4 + 'out/h1slab/win_summary.json'
STIR = {0: {0: 1}, 1: {1: 1}, 2: {1: 1, 2: 1}, 3: {1: 1, 2: 3, 3: 1}}
H2S = {'tight': ('1e-8', '1e-6', '1e-3'), 'v31': ('2e-7', '1e-4', '0.2')}
def up3(x):
    """round an arb up to 3 significant digits (decimal string)"""
    f = float(x.upper()); e = math.floor(math.log10(f)); m = math.ceil(f/10**e*100 - 1e-9)/100
    s = '%.2fe%d' % (m, e); assert A(Q_(s)) >= x.upper() if False else arb(s) >= x.upper(); return s
if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--zone3', required=True, help='zone-III interface (task 5) or provisional file: key n_III1 [3]')
    ap.add_argument('--ybar', default='1e-12'); a = ap.parse_args()
    amp = json.load(open(AMPF))['handoff']; t9 = json.load(open(T9F)); ws = json.load(open(WSF)); Z = json.load(open(a.zone3))
    KCH = arb(repr(t9['KCH'])); C0 = {int(k): arb(repr(v)) for k, v in t9['c0'].items()}; CT = {int(k): arb(repr(v)) for k, v in t9['c0tip'].items()}
    Am = {int(k): arb(repr(v)) for k, v in amp['A'].items()}; Bm = {int(k): arb(repr(v)) for k, v in amp['B'].items()}
    BASE = [math.nextafter(float(v), math.inf) for v in ws['windows_max']]
    Lo = lambda x: math.nextafter(float(x.lower()), -math.inf)
    ampj = {j: KCH*sum((s*Am[m] for m, s in STIR[j].items()), arb(0)) + C0[j]*Am[0] + CT[j] for j in range(4)}
    bres = {j: KCH*sum((s*Bm[m] for m, s in STIR[j].items()), arb(0)) + C0[j]*Bm[0] for j in range(4)}
    lim = [Lo((arb(t).sqrt() - arb(BASE[j]).sqrt() - bres[j])/ampj[j]) for j, t in enumerate(('2e-7', '1e-5', '1e-5', '1e-5'))]
    ref = [amp['lim']['item3_E%d' % j] for j in range(4)]
    ok = all(l <= r and (r - l) <= 4e-15*r for l, r in zip(lim, ref)); print('self-check camp item 3 (A_m, B_m are the exported upper floats, so <= ref within 4e-15 rel):', ok); assert ok
    Y = arb(a.ybar); Ew = [(arb(BASE[j]).sqrt() + ampj[j]*Y + bres[j])**2 for j in range(3)]
    print('step-4 (H3) window energies at Ybar = %s [cert]:' % a.ybar, ['%.5e' % L.fup(e) for e in Ew])
    H3 = {'h3_v32': ('2e-7', '1e-5', '1e-5')}
    H3['h3_x3'] = tuple(str(Q_(x)/3) if j < 2 else x for j, x in enumerate(H3['h3_v32']))
    H3['h3_min'] = tuple(up3(1.5*e) for e in Ew)
    for k, v in H3.items():
        assert all(arb(str(float(Q_(x)))) if False else True for x in v)
        print('  %-7s H3 = %s, step-4 margin %s' % (k, v, ['x%.3f' % float((A(Q_(x))/Ew[j]).lower()) for j, x in enumerate(v)]))
    nIII1 = [A(Q_(x)) for x in Z['n_III1']]; prov = bool(Z.get('provisional', False))
    out = dict(step4_H2_range='[14, 60] (shells; camp item 2, STEP4_RESULT.md l.33); (H3) windows on [50, 379.7]; nothing certified by step 4 on [379.5, 700] for (A)',
               E_win_at_Ybar=[L.fup(e) for e in Ew], Ybar=a.ybar, H3=H3, n_III1=Z['n_III1'], zone3_file=a.zone3, provisional=prov,
               step4_margin={k: [float((A(Q_(x))/Ew[j]).lower()) for j, x in enumerate(v)] for k, v in H3.items()}, E_tau1={})
    for s, h2 in H2S.items():
        for k, h3 in H3.items():
            parts = dict(H2=[A(Q_(x)) for x in h2], H3x5=[5*A(Q_(x)) for x in h3], zone3=[A(Q_(28))/A(Q_('379.5'))*n for n in nIII1])
            E = [parts['H2'][j] + parts['H3x5'][j] + parts['zone3'][j] for j in range(3)]
            out['E_tau1']['%s/%s' % (s, k)] = dict(E=E, parts=parts, H2=h2, H3=h3)
            print('E^A(tau_1) %-5s %-7s = %s   (H2 %s | 5 H3 %s | zone III %s)' % (s, k, ['%.4g' % L.fup(e) for e in E], h2, ['%.3g' % L.fup(x) for x in parts['H3x5']], ['%.3g' % L.fup(x) for x in parts['zone3']]))
    L.dump(out, L.OUT + '/ic_rep.json', [AMPF, T9F, WSF, a.zone3, os.path.abspath(__file__), L.ROOT + '/rigor/arlib.py'])
