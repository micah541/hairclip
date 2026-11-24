"""GO/NO-GO test T3 (FLOAT): outer cutoff of the shell energies of Lemma S(A) on a ramp fixed in rescaled coordinates, y in [350, 700]
(decreasing quintic cutoff, omega = (y/28)^-1, core box M0 = 0.11, M1 = 0.52; tables of lemmaS_F1, no L^4 extras).
Prints c_jk = sup of the coefficient of int_ramp z_k^2 ds (unweighted) in dE_j/dtau and K_jk = sup of the same relative to eps*omega.
Two ways to feed level 2:
  (CS)    z_3 from the level-1 dissipation of the uniformly local energies:  rate = sum_{k<=2} c_2k N_k^2 + K_23 * D3,
          D3 = int int eps omega z_3^2 per unit tau <= 2 (omega/Theta) * (level-1 source rate)/(2 a_min)
  (Young) z_3 absorbed in the dissipation of (A):  rate = sum_{k<=2} c_2k N_k^2 + K_23^2/(4 delta a_min) eps(350) omega(350) N_2^2,  delta = 0.1
Steady contribution to E_2^A = rate / gamma'_2 (0.143).  N_k^2: unweighted crude window energies of T1 (times 1.71 windows).  run from writeup/lemmaS_v32 with python3 -B (no bytecode there).
Nov 13 dictfix: reads t1's dict (key 'rows'); adds the CLEAN-GO rows (certified E^A, gamma' from lemmaS_A_rerun/rigor/assembly_arb.json). Backup: t3_A_outer_pre_dictfix.py."""
import sys, json, math
sys.dont_write_bytecode = True
sys.path.insert(0, '.')
from flint import arb
import lemmaS_F1 as F, run_F1 as R
y0, y1, gam = 350.0, 700.0, 1.0
box = R.box(0, 0, 0, M0=0.11, M1=0.52); Ab = F.A_bound_set('quasi', box); Rm = F.Ramp(y0, y1, gam, increasing=False); T = F.build_tables(); c = {}; K = {}
for j in (0, 1, 2):
    n = j + 2; tab = T[j]; YC = {k: 0.0 for k in range(n + 1)}; pairs = []
    for (a, b), mons in tab.items():
        if (a, b) == (n, n): continue
        gC = sum(co*gam**m*Ab[(i, r)] for (m, i, r), co in mons.items())
        if a == b: YC[a] += gC
        else: YC[a] += gC/2; YC[b] += gC/2
        pairs.append((a, b, mons))
    c[j] = {}; K[j] = {}
    for k in range(j + 2):
        def Yrem(t, k=k, rel=False):
            tot = arb(0)
            for (a, b, mons) in pairs:
                w = 0.0
                if a == k: w += 1.0 if a == b else 0.5
                elif b == k: w += 0.5
                if w == 0.0: continue
                for (m, i, r), co in mons.items():
                    if m == 0: continue
                    tot += arb(w*co*Ab[(i, r)])*Rm.rem(t, m)
            tot = tot*Rm.eps(t)
            if k == j: tot += arb(YC[j + 1])*Rm.rho(t)
            return tot/(Rm.eps(t)*Rm.om(t)) if rel else tot
        c[j][k] = F.sup_grid(Yrem, n=400); K[j][k] = F.sup_grid(lambda t, k=k: Yrem(t, k, True), n=400)
    print(j, 'c', {k: '%.3g' % v for k, v in c[j].items()}, 'K', {k: '%.3g' % v for k, v in K[j].items()}, flush=True)
T1 = json.load(open('../step3_gap_v32/gonogo/t1_crude_shell_model.json'))['rows']; amin = (16/17)**2; res = {}   # dictfix (Nov 13): t1 now writes a dict
AS = json.load(open('../step3_gap_v32/gonogo/lemmaS_A_rerun/rigor/assembly_arb.json'))                         # certified E^A, gamma' (arb upper ends)
CASES = (('v3.1 thresholds (H2) = (2e-7, 1e-4, 0.2)', 2, (1.81e-4, 0.247), (0.246, 0.246, 0.143)),
         ('tightened (H2) = (1e-8, 1e-6, 1e-3), estimated floor', 2, (8.2e-5, 0.048), (0.246, 0.246, 0.143)),
         ('CLEAN-GO', 2, (AS['E'][1], AS['E'][2]), tuple(AS['gp'])),
         ('CLEAN-GO', 4, (AS['E'][1], AS['E'][2]), tuple(AS['gp'])))
for case, KS, (EA1, E2A), gp in CASES:
    r = [x for x in T1 if x['case'].startswith(case) and (x['mu'], x['nu'], x['a'], x['KS']) == (4, 4, 0.3, KS)][0]
    N = [1.71*b/r['wmin'] for b in r['B']]
    src1 = 1.70*r['Gw']*EA1*1.3; D3 = 2*(28/350)/r['wmin']*src1/(2*amin)
    base = sum(c[2][k]*N[k] for k in range(3))
    cs = (base + K[2][3]*D3)/gp[2]; yg = (base + K[2][3]**2/(4*0.1*amin)*350.**-4*(28/350)*N[2])/gp[2]
    lvl1 = (c[1][0]*N[0] + c[1][1]*N[1] + c[1][2]*N[2])/gp[1]; lvl0 = (c[0][0]*N[0] + c[0][1]*N[1])/gp[0]
    key = f'{case} KS={KS}'
    res[key] = dict(N=N, D3=D3, level2_CS=cs, level2_Young=yg, level1=lvl1, level0=lvl0, E1A=EA1, E2A=E2A, gp=gp, Gw=r['Gw'], wmin=r['wmin'])
    print(f'{key}: added to E^A: level 0 {lvl0:.2e}, level 1 {lvl1:.2e} (E_1^A = {EA1:.2e}); level 2: CS {cs:.3g}, Young {yg:.3g}  against E_2^A = {E2A:.4g}  [est]')
import hashlib, os
def _sha(p): return hashlib.sha256(open(p, 'rb').read()).hexdigest()
_inp = ['../step3_gap_v32/gonogo/t1_crude_shell_model.json', '../step3_gap_v32/gonogo/lemmaS_A_rerun/rigor/assembly_arb.json', 'lemmaS_F1.py', 'run_F1.py', '../step3_gap_v32/gonogo/t3_A_outer.py']
args = dict(argv=sys.argv, cwd=os.getcwd(), python=sys.version, inputs={os.path.abspath(p): _sha(p) for p in _inp})
json.dump(dict(args=args, c=c, K=K, res=res), open('../step3_gap_v32/gonogo/t3_A_outer.json', 'w'), indent=1, default=float)
