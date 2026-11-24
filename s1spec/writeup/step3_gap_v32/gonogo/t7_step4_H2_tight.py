"""CLEAN-GO task 2 (CERTIFIED, arb): does step 4 meet the tightened hand-off shell thresholds (H2) = (1e-8, 1e-6, 1e-3)?
Reproduces camp.py's item-2 limit (step4/p2/t8se/camp.py l.236-242) from the SAME inputs, read-only:
    Ybar_2,j = L( (sqrt THR_j - sqrt BASE_j - bres_j) / amp_j ),
    shell route (t9 'shell' block present): ks = SHELL.KCH = 0, c0s = SHELL.c0 = 0  =>  amp_j = SHELL.c0tip_j (t9: sqrt E_j^sh[err] per unit Ybar), bres_j = 0,
    BASE_j = nextafter(win_summary.shells_j, +inf)   (camp: T9WIN = out/h1slab/win_summary.json).
Self-check: with the v3.1 thresholds (2e-7, 1e-4, 0.2) it must reproduce amp_A1_w180.json['handoff']['lim']['item2_E*'] bit for bit.
Also reports the energy form at the certified Ybar = 1e-12 (ycomp asserts Ybar = 1e-12 >= Y(tau_a) = 6.21e-13):
    E_j^sh <= (sqrt BASE_j + amp_j * 1e-12)^2   (arb upper)."""
import json, math
from flint import arb, ctx
ctx.prec = 64
P = '/home/claude/s1spec/step4/p2/'
AMP = json.load(open(P + 't8se/constants/amp_A1_w180.json'))['handoff']
T9 = json.load(open(P + 'out/t9_cone_kch.json'))['shell']
WS = json.load(open(P + 'out/h1slab/win_summary.json'))['shells']
assert T9['KCH'] == 0.0 and all(T9['c0'][str(j)] == 0.0 for j in range(3))
BASE = [math.nextafter(float(v), math.inf) for v in WS]
L = lambda x: math.nextafter(float(x.lower()), -math.inf); U = lambda x: math.nextafter(float(x.upper()), math.inf)   # = cconst.L / cconst.U (outward 1 ulp)
def lims(THR):
    out = []
    for j in range(3):
        amp_ = arb(repr(T9['c0tip'][str(j)])); allow = arb(THR[j]).sqrt() - arb(BASE[j]).sqrt()
        out.append(L(allow / amp_))
    return out
THR31, THRt, YBAR, YREQ = ('2e-7', '1e-4', '0.2'), ('1e-8', '1e-6', '1e-3'), arb('1e-12'), AMP['Ymax']
l31 = lims(THR31); ref = [AMP['lim'][f'item2_E{j}'] for j in range(3)]
print('self-check (v3.1 thresholds): recomputed', l31, '\n                              camp      ', ref, '\n  identical:', l31 == ref)
assert l31 == ref
lt = lims(THRt)
res = dict(inputs=dict(amp_json='step4/p2/t8se/constants/amp_A1_w180.json', t9='step4/p2/out/t9_cone_kch.json[shell]', base='step4/p2/out/h1slab/win_summary.json[shells]',
                       c0tip=[T9['c0tip'][str(j)] for j in range(3)], BASE=BASE, Ybar_certified=1e-12, Ybar_binding_step4=YREQ, binding_item=AMP['binding']),
           v31=dict(THR=THR31, Ybar_limit=l31, margin_vs_Ybar=[x/1e-12 for x in l31]), tight=dict(THR=THRt, Ybar_limit=lt, margin_vs_Ybar=[x/1e-12 for x in lt]))
E = []
for j in range(3):
    e = (arb(BASE[j]).sqrt() + arb(repr(T9['c0tip'][str(j)]))*YBAR)**2; E.append(U(e))
res['E_shell_at_Ybar_1e-12'] = E
res['tight']['margin_energy'] = [float(t)/e for t, e in zip(THRt, E)]
res['v31']['margin_energy'] = [float(t)/e for t, e in zip(THR31, E)]
res['tight']['still_not_binding'] = min(lt) > YREQ
print('tightened (H2) = (1e-8, 1e-6, 1e-3): Ybar limits (arb lower)', ['%.10e' % x for x in lt])
print('  margins vs certified Ybar = 1e-12:', ['x%.4g' % (x/1e-12) for x in lt], '| step-4 binding limit stays', AMP['binding'], '%.6e' % YREQ, ':', min(lt) > YREQ)
print('shell energies at Ybar = 1e-12 (arb upper):', ['%.6e' % x for x in E])
print('  energy margins tightened:', ['x%.4g' % m for m in res['tight']['margin_energy']], ' v3.1:', ['x%.4g' % m for m in res['v31']['margin_energy']])
json.dump(res, open('/home/claude/s1spec/writeup/step3_gap_v32/gonogo/t7_step4_H2_tight.json', 'w'), indent=1)
