"""Phase 1 checkpoint (P4): t1 CLEAN-GO shell model + launch_model accounting with Phase-1 certified inputs.  Float arithmetic on certified numbers.
Labels: [cert] = value read from a certified JSON; [est] = float model output; [est, unbounded in proof] = placeholder for a term with no proof yet.
Nothing is derived here: formulas are copied from t1_crude_shell_model.py and launch/launch_model.py; inputs are swapped via the dict cfg.
Self-check first (cfg ORIG = t1's inputs): must reproduce t1 CLEAN-GO M1 (4,4,0.3) KS=1 (x1.72) and launch_model LATE KS=1 M1 (x2.13).
Run: PYTHONHASHSEED=0 nice -n 19 python3 -B t1_checkpoint.py"""
import sys, os, json, math, hashlib, time
sys.dont_write_bytecode = True
import numpy as np
T0 = time.time()
R = '/home/claude/s1spec/writeup'; GO = f'{R}/step3_gap_v32/gonogo'; PH = f'{R}/lemmaU/phase1'
INP = dict(assembly_ibp=f'{GO}/lemmaS_A_ibp/rigor/assembly_arb.json', assembly_rerun=f'{GO}/lemmaS_A_rerun/rigor/assembly_arb.json',
           t8=f'{GO}/t8_L4_cost.json', feas_P=f'{R}/step3_gap_v32/feas_P.json', t1_json=f'{GO}/t1_crude_shell_model.json',
           launch_json=f'{PH}/launch/launch_model.json', ghat=f'{PH}/ghat/ghat_cert.json', ramp=f'{PH}/ramp/ramp_K_cert.json')
J = {k: json.load(open(v)) for k, v in INP.items()}
amin = (16/17)**2; eps28 = 28.**-4; l = math.log(1.5); p = 1.0246; A1 = 1.5e-4
H2 = (1e-8, 1e-6, 1e-3); H3 = (2e-7, 1e-5, 1e-5)                 # [paper] hand-off hypotheses (LAUNCH.md): typed, no JSON holds them
bbar, abar = 3.33811e-9, 2.19103e-7                              # [cert] verifier sups: typed (as t1 / launch_model)
DT6 = (9.86e-08, 1.04e-06, 0.0161)                               # [est] t5_t6_fine.out (Lemma R ramp): typed (as t1 / launch_model)
c_dU = 1.3056*(bbar + p*abar*A1)/(p*(0.25 - abar))
dU = tuple((p**j*c_dU*56**p)**2*math.log(2) for j in range(3))
xs = np.linspace(-20, 20, 400001); dx = xs[1] - xs[0]

def Wf(mu, nu, a, x): return 1/(np.exp(-mu*(x + a)) + 1 + np.exp(nu*(x - l - a)))
def supW(mu, nu, a, s0, s1, fac=None, n=4001):
    x = np.linspace(s0, s1, n); w = Wf(mu, nu, a, x); w = w*fac(x) if fac else w; return float(w.max())*(1 + 1e-3)
def IW_in(mu, nu, a, yc): m = xs <= math.log(56/yc); return 4*float(Wf(mu, nu, a, xs[m]).sum()*dx)
_wc = {}
def weight(mu, nu, a):
    k = (mu, nu, a)
    if k not in _wc:
        W = Wf(mu, nu, a, xs); wmin = float(W[(xs >= 0) & (xs <= l)].min())
        CA = supW(mu, nu, a, math.log(28/56), math.log(350/56), fac=lambda x: 2*np.exp(x))
        q = sum(supW(mu, nu, a, math.log(350/56) + i*l, math.log(350/56) + (i + 1)*l)/wmin for i in range(60))
        s = np.linspace(math.log(350), math.log(350) + 30, 300001)
        base = float(np.sum(Wf(mu, nu, a, s - math.log(56))*(c_dU*np.exp(p*s))**2)*(s[1] - s[0]))
        _wc[k] = dict(IW=4*float(W.sum()*dx), IWin56=IW_in(mu, nu, a, 56.), wmin=wmin, CA=CA, q=q, Dd=[p**(2*j)*base for j in range(3)])
    return _wc[k]

# ---------------- configurations ----------------
def Kfeas():
    K = {0: {0: 1.0e3, 1: 1.0e3}, 1: {0: 2e-12, 1: 1270., 2: 725.}, 2: {0: 2.3e-11, 1: 582., 2: 3100., 3: 704.}}
    for j, row in J['feas_P']['4.0']['K'].items():
        for k, v in row.items(): K[int(j)][int(k)] = float(v)
    return K
def cfg_orig():
    AS = J['assembly_rerun']; c = J['t8']['cases']['crude|mu=4,nu=4,a=0.3']; rr = max(c['L4_extras_ratio_general_over_power'], 1.)
    return dict(name='ORIG (t1 inputs)', EA=tuple(AS['E']), gp=tuple(AS['gp']), K=Kfeas(), footw={1: 1.5707250797145307, 2: 2862.4867820255}, rr=rr,
                Kf=1.17e4*rr, G22=c['Ghat_L4'][2])
def cfg_p1(G22, tag, cubic223):
    AS = J['assembly_ibp']; r = J['ramp']; K = {int(j): {int(k): float(v) for k, v in row.items()} for j, row in r['4.0']['K'].items()}
    ce = r['C_eta_Xi']; f2 = r['foot']['2']
    coef = f2['Kf_unnorm_S33'] + r['cubic_remaining']['2|(2, 2, 2)']           # 26.917 + 72.299 = 99.2 [cert arithmetic]
    coef23 = r['cubic_remaining']['2|(2, 2, 3)']*ce/3 if cubic223 else 0.        # [est, unbounded in proof] placeholder, like the (2,3,3) class: weight*C_eta/3
    return dict(name=f'PHASE1 G22={tag}' + (' +(2,2,3)' if cubic223 else ' no(2,2,3)'), EA=tuple(AS['E']), gp=tuple(AS['gp']), K=K,
                footw={1: r['foot']['1']['footw'], 2: f2['footw']}, rr=1., Kf=(coef + coef23)*28.**2, coef=coef, coef23=coef23, G22=G22, cubic223=cubic223)

# ---------------- model pieces (copied from t1 / launch_model) ----------------
def foot_terms(cfg, IW, Gw):
    EA = cfg['EA']; F2 = eps28*Gw*IW*2*EA[2]; JG = 2*(cfg['gp'][1]*EA[1]/(2*amin))*IW*Gw
    s33 = cfg['Kf']*math.sqrt(eps28*Gw*2*EA[2])*(F2/math.sqrt(math.log(2)) + math.sqrt(2)*F2**0.75*JG**0.25)
    rr = cfg['rr']; return (0., cfg['footw'][1]*rr*F2, cfg['footw'][2]*rr*F2 + s33)
def evol(cfg, IW, Gw, KS):
    EA, gp, K = cfg['EA'], cfg['gp'], cfg['K']
    Ek = [(math.sqrt(EA[k]) + math.sqrt(dU[k]))**2 for k in range(3)]
    X = {1: eps28*Ek[1], 2: gp[0]*EA[0]/(2*amin) + eps28*p**4*dU[0], 3: gp[1]*EA[1]/(2*amin) + eps28*p**6*dU[0]}
    ft = foot_terms(cfg, IW, Gw); out = []
    for j in range(3):
        br = sum(K[j][k]*(eps28*Ek[k] if k <= j else X[k]) for k in K[j])     # k = j+1 and (new) k = j+2 use the z_k dissipation data X_k
        out.append(dict(inflow=1.70*IW*Gw*Ek[j], ramp=Gw*IW*2*KS*br, foot=ft[j], outer=DT6[j]))
    return out
def box(B, wmin, Uwin):
    Eu = [b/wmin for b in B]
    S0 = math.sqrt(Eu[0]/l + 2*math.sqrt(Eu[0]*Eu[1])); S1 = math.sqrt(Eu[1]/l + 2*math.sqrt(Eu[1]*(Eu[1] + Eu[2])))
    return dict(f=Uwin + S0, yfp=1.05*Uwin + S1, M0=0.11/(Uwin + S0), M1=0.52/(1.05*Uwin + S1))
def Gwof(cfg): return 1/(1 - cfg['G22']/28.**4)
def sweep(cfg, mu, nu, a, KS, withH3=False):                  # t1 model: whole weight, hand-off data (60/28) H2 [+ (H3) windows]
    w = weight(mu, nu, a); Gw = Gwof(cfg); ev = evol(cfg, w['IW'], Gw, KS)
    dH3 = sum(supW(mu, nu, a, math.log(yw/56), math.log(1.5*yw/56))*(1.5*yw/28) for yw in (50*1.5**k for k in range(5))) if withH3 else 0.
    B = [Gw*((60/28)*H2[j] + dH3*H3[j]) + sum(ev[j].values()) for j in range(3)]
    return dict(B=B, **box(B, w['wmin'], 1.3056*A1*84**p))
def late(cfg, mu, nu, a, KS, conservative=False):
    w = weight(mu, nu, a); Gw = Gwof(cfg); EA = cfg['EA']; q, Dd, CA = w['q'], w['Dd'], w['CA']
    ev = evol(cfg, w['IWin56'], Gw, KS)
    if conservative:
        evr = evol(cfg, w['IW'], Gw, KS); ev = [dict(ev[j], ramp=evr[j]['ramp'], foot=evr[j]['foot']) for j in range(3)]
    B = [0.]*3
    for it in range(200):
        D = [Gw*(CA*EA[j] + (math.sqrt(q*B[j]) + math.sqrt(Dd[j]))**2) for j in range(3)]
        Bn = [D[j] + sum(ev[j].values()) for j in range(3)]
        if max(abs(Bn[j] - B[j])/max(Bn[j], 1e-300) for j in range(3)) < 1e-14: B = Bn; break
        B = Bn
    return dict(B=B, D=D, ev=ev, nfac=1/(1 - Gw*q), **box(B, w['wmin'], 1.3056*A1*84**p))
YC = (56., 60., 70., 84., 100., 126.6, 150., 200., 253.)
def present(cfg, mu, nu, a, KS, conservative=False):
    w = weight(mu, nu, a); Gw = Gwof(cfg); rows = []
    for yc in YC:
        ev = evol(cfg, w['IW'] if conservative else IW_in(mu, nu, a, yc), Gw, KS)
        dH2 = supW(mu, nu, a, math.log(28/yc), math.log(60/yc), fac=lambda x, yc=yc: yc*np.exp(x)/28)
        dH3 = sum(supW(mu, nu, a, math.log(yw/yc), math.log(1.5*yw/yc))*(1.5*yw/28) for yw in (50*1.5**k for k in range(5)))
        B = [Gw*(dH2*H2[j] + dH3*H3[j]) + sum(ev[j].values()) for j in range(3)]
        rows.append(dict(yc=yc, B=B, **box(B, w['wmin'], 1.3056*A1*(1.5*yc)**p)))
    return rows
def present_worst(rows):                                       # zones I/II used by the lemma: y_c <= 126.7 (zone III, y_c >= 190, not used)
    r = [x for x in rows if x['yc'] <= 126.7]; return dict(M0=min(x['M0'] for x in r), M1=min(x['M1'] for x in r),
                                                              yc_M0=min(r, key=lambda x: x['M0'])['yc'], yc_M1=min(r, key=lambda x: x['M1'])['yc'], rows=rows)
def bis(fun, lo, hi):
    for _ in range(80):
        mid = math.sqrt(lo*hi); (lo, hi) = (mid, hi) if fun(mid) > 1 else (lo, mid)
    return lo

res = dict(selfcheck={}, rows=[], flips={})
args = dict(argv=sys.argv, cwd=os.getcwd(), python=sys.version, env={k: os.environ.get(k) for k in ('PYTHONHASHSEED', 'PYTHONDONTWRITEBYTECODE')},
            inputs={k: dict(path=v, sha256=hashlib.sha256(open(v, 'rb').read()).hexdigest()) for k, v in INP.items()},
            script_sha256=hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest())
res['args'] = args
# ---------------- self-check ----------------
co = cfg_orig(); ok = True
print('SELF-CHECK (ORIG t1 inputs)')
t1row = [x for x in J['t1_json']['rows'] if x['case'].startswith('CLEAN-GO') and (x['mu'], x['nu'], x['a'], x['KS']) == (4, 4, 0.3, 1)][0]
s = sweep(co, 4, 4, 0.3, 1); d1 = abs(s['M1'] - t1row['margin_M1'])/t1row['margin_M1']; d0 = abs(s['M0'] - t1row['margin_M0'])/t1row['margin_M0']
print(f"  t1 sweep (4,4,0.3) KS=1: M1 x{s['M1']:.4f} (t1 json x{t1row['margin_M1']:.4f}, rel dev {d1:.1e}); M0 x{s['M0']:.3f} (t1 x{t1row['margin_M0']:.3f}, dev {d0:.1e})")
ok &= d1 < 1e-6 and d0 < 1e-6 and round(s['M1'], 2) == 1.72
lt = late(co, 4, 4, 0.3, 1); lj = J['launch_json']['late']['mu=4,nu=4,a=0.3']['rows']['1']
d2 = abs(lt['M1'] - lj['M1'])/lj['M1']; d3 = abs(lt['M0'] - lj['M0'])/lj['M0']
lc = late(co, 4, 4, 0.3, 1, True); d4 = abs(lc['M1'] - lj['rev_variant']['M1'])/lj['rev_variant']['M1']
print(f"  launch LATE (4,4,0.3) KS=1: M1 x{lt['M1']:.4f} (launch_model.json x{lj['M1']:.4f}, dev {d2:.1e}); M0 x{lt['M0']:.3f} (x{lj['M0']:.3f}); conservative M1 x{lc['M1']:.3f} (x{lj['rev_variant']['M1']:.3f}, dev {d4:.1e})")
pr = present(co, 4, 4, 0.3, 1); pj = J['launch_json']['present']['mu=4,nu=4,a=0.3']['1']
d5 = max(abs(a_['M1'] - b_['M1'])/b_['M1'] for a_, b_ in zip([x for x in pr if x['yc'] != 253.] , [y for y in pj if abs(y['yc'] - 126.7) > 1e-9 and y['yc'] != 253.0 and abs(y['yc']-190/1.5)>1e-9]))
print(f'  launch PRESENT (y_c list) max rel dev of M1 vs launch_model.json: {d5:.1e}')
ok &= d2 < 1e-6 and d3 < 1e-6 and d4 < 1e-6 and round(lt['M1'], 2) == 2.13 and d5 < 1e-6
s3 = sweep(co, 4, 4, 0.3, 1, True); print(f"  t1 sweep with (H3): M1 x{s3['M1']:.3f} (LAUNCH.md: x1.666)"); ok &= round(s3['M1'], 3) == 1.666
res['selfcheck'] = dict(passed=bool(ok), t1_M1=s['M1'], t1_json_M1=t1row['margin_M1'], launch_late_M1=lt['M1'], launch_json_M1=lj['M1'], launch_cons_M1=lc['M1'], t1_H3_M1=s3['M1'])
print('  SELF-CHECK', 'PASSED' if ok else 'FAILED'); assert ok

g = J['ghat']; r = J['ramp']
G22_W = max(g['sensitivity_Umax_0.09']['Ghat_jj'][2][1], r['foot']['2']['Ghat_jj_with_foot_adds'])
import copy
print('\nREVIEWER STRESS ROWS [rev, est]  weight (4,4,0.3), worst Ghat22 = %.1f' % G22_W)
def rows(tag, fac, G22=G22_W):
    for KS in (1, 2, 4, 5.5, 8):
        cfg = cfg_p1(G22, 'x', True); cfg['EA'] = tuple(e*f for e, f in zip(cfg['EA'], fac))
        ll = late(cfg, 4, 4, 0.3, KS); lc = late(cfg, 4, 4, 0.3, KS, True); pl = present_worst(present(cfg, 4, 4, 0.3, KS)); sw = sweep(cfg, 4, 4, 0.3, KS, True)
        print('  %-34s KS=%-3s Gw %.2f: LATE launch M0 x%.2f M1 x%.3f (B2 %.3f) | LATE cons M1 x%.3f | PRESENT launch M0 x%.2f M1 x%.2f | t1-sweep+H3 M1 x%.3f' % (tag, KS, Gwof(cfg), ll['M0'], ll['M1'], ll['B'][2], lc['M1'], pl['M0'], pl['M1'], sw['M1']))
rows('E^A as certified', (1, 1, 1))
rows('E^A x(1.41,1.13,1.25)', (1.41, 1.13, 1.25))
rows('E^A x(1.41,1.13,1.65)', (1.41, 1.13, 1.65))
rows('E^A x(1.5,1.5,2.5)', (1.5, 1.5, 2.5))
rows('E^A x(1.41,1.13,1.65), Ghat 3.3e5', (1.41, 1.13, 1.65), 3.3e5)
rows('E^A x(1.41,1.13,1.65), Gw->1.03', (1.41, 1.13, 1.65), 28.**4*(1 - 1/1.03))
