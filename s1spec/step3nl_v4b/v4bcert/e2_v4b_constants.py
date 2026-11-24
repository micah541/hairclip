"""fin13v4b E2 fix (E2_AUDIT.md section 5, F1): verifier constants of the E2 level, from the certified outputs of e2_consts.py and
e2_nlin_regions.py (both promoted from fix5/e2_audit, run in this chain), raw -> one ulp up.  No typed literals except the hypotheses
(AMAX_E2 = 2.5e-4, ETA_TIP_E2 = 1e-4 on [0,8]), which are asserted in closure.box_assert.
Nonlinear constants (rhs_iv.py ND): nd_k := max over R1 ([0,8], box 1e-4) and R2a ([8,12], box 1e-3); nf_k := max over R2b ([12,16.25], 1e-3) and
R3 ([16.25,34.25], 0.02); columns (0,0,1) w'''' -> nd1/nf1, (0,1,0) w''' -> 2, (1,0,0) w'' -> 3, A0 -> nd4 (A), A1 -> nd5 (W1), nf7 = max(A0, A1) (A and W1),
(2,0,0) -> nd6/nf4, (1,1,0) -> nd7/nf5, (3,0,0) -> nd8/nf6.  Also: the regions' jet boxes and alpha level are asserted to be the stated hypotheses.
Output out/e2v4b_constants.json."""
import sys, os, json, math
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from v4bcommon import OUT, args_block, dump
P1 = os.path.join(OUT, 'e2_consts.json'); P2 = os.path.join(OUT, 'nl_regions_a2.5e-4.json')
E = json.load(open(P1)); N = json.load(open(P2))
assert N['alpha_max'] == 2.5e-4 and N['pad_c4'] == 1
R = {r['name']: r for r in N['regions']}
assert (R['R1']['xlo'], R['R1']['xhi'], R['R1']['eta0'], R['R1']['eta1']) == (0.0, 8.0, 1e-4, 1e-4)
assert (R['R2a']['xlo'], R['R2a']['xhi'], R['R2a']['eta0'], R['R2a']['eta1']) == (8.0, 12.0, 1e-3, 1e-3)
assert (R['R2b']['xlo'], R['R2b']['xhi'], R['R2b']['eta0'], R['R2b']['eta1']) == (12.0, 16.25, 1e-3, 1e-3)
assert (R['R3']['xlo'], R['R3']['xhi'], R['R3']['eta0'], R['R3']['eta1']) == (16.25, 34.25, 0.02, 0.02)
S = N['sup']; C = {}
def put(name, raw, src):
    r = float(raw); v = math.nextafter(r, math.inf); C[name] = dict(raw=r, value=v, dir='up', source=src)
col = {1: '(0, 0, 1)', 2: '(0, 1, 0)', 3: '(1, 0, 0)'}
for i, key in col.items():
    put(f'ND{i}', max(S['R1'][key], S['R2a'][key]), f'nl_regions max(R1,R2a) {key}'); put(f'NF{i}', max(S['R2b'][key], S['R3'][key]), f'nl_regions max(R2b,R3) {key}')
put('ND4', max(S['R1']['A0'], S['R2a']['A0']), 'nl_regions max(R1,R2a) A0'); put('ND5', max(S['R1']['A1'], S['R2a']['A1']), 'nl_regions max(R1,R2a) A1')
for i, key in ((6, '(2, 0, 0)'), (7, '(1, 1, 0)'), (8, '(3, 0, 0)')): put(f'ND{i}', max(S['R1'][key], S['R2a'][key]), f'nl_regions max(R1,R2a) {key}')
for i, key in ((4, '(2, 0, 0)'), (5, '(1, 1, 0)'), (6, '(3, 0, 0)')): put(f'NF{i}', max(S['R2b'][key], S['R3'][key]), f'nl_regions max(R2b,R3) {key}')
put('NF7', max(S['R2b']['A0'], S['R3']['A0'], S['R2b']['A1'], S['R3']['A1']), 'nl_regions max(R2b,R3) of A0, A1')
g = E['c_garding']; put('R143', g['r143']['raw'], 'garding_cert.json via e2_consts c_garding.r143')
for k in ('r163', 'r107'): assert g[k]['ok']
dm = E['d_modulation']
put('DAB', dm['dab_lead_minus_Q2lead'], 'e2_consts d_modulation |lemma lead - Minv l|/A')
put('NF12', dm['NF12'], 'e2_consts lambda_max^(1/2) Gram(f1,f2)')
put('N3E2', dm['n3_c4'], 'e2_consts ||chi psi3|| (C^4 padded)')
put('NH2R', dm['nh2_c4'], 'e2_consts ||h2|| (C^4 padded)')
put('SRW', E['e_strip']['SRW_sup_sqrt_rho_over_W'], 'e2_consts sup_strip (rho/W)^(1/2)')
assert E['b_Drho']['ok'] and E['c_modE2']['ok']
C['AMAX_E2'] = dict(raw=2.5e-4, value=2.5e-4, dir='exact', source='hypothesis H-E2a, asserted in closure.box_assert')
C['ETA_TIP_E2'] = dict(raw=1e-4, value=1e-4, dir='exact', source='hypothesis H-E2b: sup(|w|,|w\'|) <= 1e-4 on [0,8], asserted in closure.box_assert')
typed = dict(nd1=6.9e-3, nd2=1.2e-2, nd3=4.2e-2, nd4=0.54, nd5=0.62, nd6=4.75, nd7=0.136, nd8=3.05, nf1=4.3e-3, nf2=6.9e-4, nf3=2.9e-5, nf4=0.014, nf5=0.46, nf6=3.06, nf7=8.6e-6)
cmp = {k: dict(typed=v, new=C[k.upper()]['value'], typed_covers=bool(v >= C[k.upper()]['value'])) for k, v in typed.items()}
out = dict(constants=C, typed_fin13v4_vs_new=cmp, args=args_block([P1, P2]))
print(json.dumps({k: (v['raw'], v['value']) for k, v in C.items()}, indent=0)); print(json.dumps(cmp))
print('sha256', dump(os.path.join(OUT, 'e2v4b_constants.json'), out))
