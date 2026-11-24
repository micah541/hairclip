#!/usr/bin/env python3
"""Kill check 2 reviewer scratch [rev, float]: U_model2.py (unmodified source, exec'd up to its 'runs' section, self-check included) with
 (a) window-uniform dissipation data: iint_{tau_rho}^{tau_rho+T} <= C0 + C1*T for EVERY tau_rho (A_REPAIR R4 / assembly_rep 'export'),
     consumed by Abel summation against the tau-decreasing ramp weight:  X <= S_X*C0 + C1*Sigma_X,  Sigma_X = int_tau sup_ramp(Psi/omega);
 (b) live sensitivity rows on A's certified export (H-O placeholders x10, E2A at the rule limit, fallback box, KS).
Writes only stdout.  Run: PYTHONHASHSEED=0 nice -n 19 python3 -B kc2_model.py > kc2_model.out"""
import sys, os, json, math
sys.dont_write_bytecode = True
SRC = '/home/claude/s1spec/writeup/lemmaU/phase2a/U_ramp/U_model2.py'
EXP = '/home/claude/s1spec/writeup/lemmaU/phase2a/A_repair/export_for_U.json'
src = open(SRC).read()
head = src.split('# ================= runs =================')[0]
a = "        return (math.sqrt(S_X*Dk) + math.sqrt(eps28*S_E*p**(2*k)*dU[0]))**2"
b = "        return (math.sqrt(S_X*Dk + cfg.get('C1', (0., 0., 0.))[k - 2]*SIG_X) + math.sqrt(eps28*S_E*p**(2*k)*dU[0]))**2"
assert head.count(a) == 1; head = head.replace(a, b)
a = "    D = cfg['D']\n"; b = "    D = cfg['D']; SIG_X = sweep(Dr*yR/28, yc, mu, nu, a)\n"
assert head.count(a) == 1; head = head.replace(a, b)
a = "    return out, dict(Gw=Gw, S_in=S_in, S_E=S_E, S_X=S_X, Xi=Xi"; b = "    return out, dict(Gw=Gw, S_in=S_in, S_E=S_E, S_X=S_X, SIG_X=SIG_X, Xi=Xi"
assert head.count(a) == 1; head = head.replace(a, b)
head = head.replace("ap.parse_args()", "ap.parse_args(['--inputs', EXP])").replace("os.path.abspath(__file__)", "SRC")
g = dict(__name__='kc2', __file__=SRC, EXP=EXP, SRC=SRC); exec(compile(head, SRC, 'exec'), g)
cfg_new, late, present, present_worst, box_flag = g['cfg_new'], g['late'], g['present'], g['present_worst'], g['cfg_box_is_wide']
X = json.load(open(EXP)); EA = tuple(X['EA']); gp = tuple(X['gp']); r_am = [X['kappa'][k]*X['amin'] for k in range(3)]
D_tau1 = tuple(X['D']); C0 = tuple(X['D_from_any_tau_a_with_E_tau_a_le_EA']); C1 = tuple(X['D_slope'])
assert max(int(k) for k in g['K_cert'][2]) <= 3 and max(int(k) for k in g['K_cert'][1]) <= 3          # D[2] (=71) is never consumed
def run(tag, inp, KS=1, C1v=None, boxname='widened (0.15,0.52), U .09, M0z .16', pres=True):
    cfg = cfg_new(dict(inp), 'data56', boxname, KS=KS)
    if C1v: cfg['C1'] = C1v
    box_flag[0] = boxname.startswith('widened (0.15')
    lt = late(cfg); e2 = lt['ev'][2]; pw = present_worst(present(cfg)) if pres else None
    print(f"{tag:66s} KS={KS} g {g['gD'](cfg):.3g} LATE M1 x{lt['M1']:.3f} M0(.11) x{lt['M0']:.2f} B2 {lt['B'][2]:.4f} [launch {lt['D'][2]:.4f} inflow {e2['inflow']:.4f} rampX {e2['rampX']:.4f} outer {e2['outer']:.4f}]"
          + (f" | PRESENT M1 x{pw['M1']:.2f} M0(.11) x{pw['M0']:.2f}" if pw else '') + f"  S_X {lt['info']['S_X']:.3f} Sig_X {lt['info']['SIG_X']:.3f}", flush=True)
    return lt
base = dict(EA=EA, gp=gp)
print('\n== R0 reproduce U_model2_A1 headline (D from tau_1, slope ignored)')
for KS in (1, 2): run('R0 D = D(tau_1) = (1.40e-7, 7.59e-5), slope ignored', dict(base, D=D_tau1), KS)
print('\n== R1 window-uniform data: X <= S_X*C0 + C1*Sigma_X (C0 = restart value, C1 = slope at H-O placeholders)')
for KS in (1, 1.33, 2, 4): run('R1 C0 = (2.92e-7, 1.51e-4), C1 = (5.1e-9, 5.9e-9)', dict(base, D=C0), KS, C1v=C1, pres=KS in (1, 2))
print('\n== R2 R1 + H-O placeholders x10 (outer rates (6.2e-7, 1.31e-6, 9.06e-4) per unit tau)')
r1 = (6.165e-9, 7.120e-9, 1.2614e-7); r10 = (6.21e-7, 1.31e-6, 9.06e-4)
EA10 = tuple(EA[k] + (r10[k] - r1[k])/gp[k] for k in range(3))
C0_10 = tuple(C0[k] + (EA10[k] - EA[k])/r_am[k] for k in range(3)); C1_10 = tuple(C1[k] + (r10[k] - r1[k])/r_am[k] for k in range(3))
print('   E^A(x10) =', ['%.4g' % x for x in EA10], ' C0 =', ['%.4g' % x for x in C0_10[:2]], ' C1 =', ['%.4g' % x for x in C1_10[:2]])
for KS in (1, 2): run('R2 H-O x10', dict(EA=EA10, gp=gp, D=C0_10), KS, C1v=C1_10)
print('\n== R3 R1 with E2A at the rule limit 0.1 (E0, E1 as certified)')
for KS in (1, 2): run('R3 E2A = 0.1', dict(EA=(EA[0], EA[1], 0.1), gp=gp, D=C0), KS, C1v=C1, pres=False)
print('\n== R4 R1 with the fallback boxes (K NOT certified for these boxes / g; M1 still measured against 0.52)')
for bn in ('widened (0.15,0.60), U .09, M0z .16', 'widened (0.20,0.60), U .09, M0z .20'):
    for KS in (1, 2): run('R4 ' + bn, dict(base, D=C0), KS, C1v=C1, boxname=bn, pres=False)
print('\n== R5 joint stress: H-O x10, KS = 2, C0 x1.5, E2A = 0.1')
run('R5 joint', dict(EA=(EA10[0], EA10[1], 0.1), gp=gp, D=tuple(1.5*x for x in C0_10)), 2, C1v=C1_10, pres=False)
run('R5b joint, fallback box (0.20,0.60)', dict(EA=(EA10[0], EA10[1], 0.1), gp=gp, D=tuple(1.5*x for x in C0_10)), 2, C1v=C1_10, boxname='widened (0.20,0.60), U .09, M0z .20', pres=False)
print('\n== R6 the slope if consumed naively as D + D_slope*T (T from tau_1), late centre launched at tau_rho = tau_1 + T: level-1 datum')
for T in (0, 100, 1e3, 1e4, 1e5): print(f'   T = {T:g}: D_1 = {D_tau1[1] + C1[1]*T:.4g}  (x{(D_tau1[1] + C1[1]*T)/D_tau1[1]:.2f});  D_0 = {D_tau1[0] + C1[0]*T:.4g} (x{(D_tau1[0] + C1[0]*T)/D_tau1[0]:.2f})')
