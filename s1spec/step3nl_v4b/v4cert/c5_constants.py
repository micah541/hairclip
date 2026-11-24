#!/usr/bin/env python3
"""Stage C5 (fin13v4): build out/constants_v4.json from the raw certified JSONs.
Every entry: dict(value=<double used by the verifier>, raw=<raw value read>, dir='up'|'down', source=<path>, key=<json key>).
Rule: value = nextafter(raw, +inf) for upper bounds ('up'), nextafter(raw, -inf) for lower bounds ('down': MU1, MU2), so value >= raw (resp. <=)
also when the producer stored float(upper) rounded to nearest.  Entries produced by the v4 stages C1-C4, U3 are already outward floats
(nextafter of an arb upper endpoint); they get one more nextafter (harmless).  A test at the end asserts value >= raw (up) / <= raw (down)
for every entry and that no entry is NaN/inf.  The verifier's constants.py reads this file (NEW=1) and contains no typed certified literal
for any entry listed here.  Legacy E/E2-block literals (outside EI_SPEC) are listed in 'legacy_typed' with their fin12f4 provenance.
Run: PYTHONHASHSEED=0 python3 c5_constants.py"""
import sys, os, json, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from v4common import args_block, dump, OUT, V4
J = lambda p: json.load(open(p))
P = {'coup': os.path.join(V4, 'agent_coup/coupling_constants.json'), 'sob': os.path.join(V4, 'agent_sob/sob_constants_P2.json'),
     'gard': os.path.join(V4, 'agent_int/v2/v2_T4_P2_0.25.json'), 'v3': os.path.join(V4, 'agent_int/v3/int_constants_v3.json'),
     'S2': os.path.join(V4, 'agent_int/c4switch/c4switch_S2_1.json'), 'S3': os.path.join(V4, 'agent_int/c4switch/c4switch_S3_1e-3.json'),
     'mod': os.path.join(V4, 'agent_mod/ibp/final_5e-4.json'), 'c1': os.path.join(OUT, 'c1_cut_sups.json'), 'c2': os.path.join(OUT, 'c2_t3d.json'),
     'c3': os.path.join(OUT, 'c3_modx.json'), 'c4': os.path.join(OUT, 'c4_hardy.json'), 'u3': os.path.join(OUT, 'u3_strip_form.json')}
D = {k: J(p) for k, p in P.items()}
C = {}
def put(name, raw, src, key, dir_='up'):
    raw = float(raw); assert math.isfinite(raw) and raw >= 0 or dir_ == 'down', (name, raw)
    v = math.nextafter(raw, math.inf if dir_ == 'up' else -math.inf)
    if raw == 0.0: v = 0.0
    C[name] = dict(value=v, raw=raw, dir=dir_, source=P[src], key=key)
def get(d, path):
    for p in path.split('/'): d = d[p]
    return d
# ---- coupling (agent_coup, P=2 rigorous)
cp = D['coup']['P']['2']['rigorous']
for nm, k in (('AF', 'a_f'), ('BF', 'b_f'), ('CF0', 'c3'), ('CF1', 'c2'), ('CF2', 'c1'), ('CF3', 'c0'), ('SB', 'sB'), ('SA', 'sA'), ('IA', 'IA'),
              ('DS0', 'd3'), ('DS1', 'd2'), ('DS2', 'd1'), ('DS3', 'd0'), ('PHI34', 'Phi_34.25'), ('K_rho_16.25', 'K_rho_16.25')):
    put(nm, cp[k], 'coup', 'P/2/rigorous/' + k)
# ---- Sobolev (agent_sob, P=2)
sb = D['sob']
put('c1', sb['rho']['c1'], 'sob', 'rho/c1'); put('PHX', sb['rho']['Phi1_X'], 'sob', 'rho/Phi1_X'); put('K2', sb['rho']['K_b2'], 'sob', 'rho/K_b2')
put('K12', sb['rho']['K_b12'], 'sob', 'rho/K_b12')
for r in ('[0,1]', '[1,2]', '[2,4]', '[4,8]', '[8,12]'): put('KREG' + r, sb['rho']['K_region'][r], 'sob', 'rho/K_region/' + r)
put('KW', sb['W']['K_W_region']['[30,34.25]'], 'sob', 'W/K_W_region/[30,34.25]')
put('PSI1W', sb['W']['Psi1']['R=60']['from_10.25'], 'sob', 'W/Psi1/R=60/from_10.25'); put('C1W', sb['W']['c1W']['R=60'], 'sob', 'W/c1W/R=60')
# ---- Garding (v2_T4_P2_0.25): MU1, MU2 lower bounds = min over regions of -(sup(-mu)); E1, E0 = max over regions
g = D['gard']['sup']
put('MU1', min(-v['-mu1'] for v in g.values()), 'gard', 'sup/*/-mu1 (min over regions of -value)', 'down')
put('MU2', min(-v['-mu2'] for v in g.values()), 'gard', 'sup/*/-mu2 (min over regions of -value)', 'down')
put('E1', max(v['|e1|'] for v in g.values()), 'gard', 'sup/*/|e1|'); put('E0', max(v['|e0|'] for v in g.values()), 'gard', 'sup/*/|e0|')
# ---- forcing (C^4 functions)
put('FN0', D['S2']['psi']['r_half'], 'S2', 'psi/r_half'); put('FN1', D['S2']['h2']['r_half'], 'S2', 'h2/r_half'); put('C3', D['S3']['C3'], 'S3', 'C3')
assert D['S3']['alpha_max'] == 0.001
v3 = D['v3']['R=60']
put('FNW0', v3['F1b_psi3_residual_level0'], 'v3', 'R=60/F1b_psi3_residual_level0'); put('FNW1', v3['F1b_h2_residual_level0'], 'v3', 'R=60/F1b_h2_residual_level0')
put('FNW2', v3['F1a_C3W_level0_for_C4_functions'], 'v3', 'R=60/F1a_C3W_level0_for_C4_functions')
put('FNWd0', v3['F1b_psi3_residual_level1_zeta'], 'v3', 'R=60/F1b_psi3_residual_level1_zeta'); put('FNWd1', v3['F1b_h2_residual_level1_zeta'], 'v3', 'R=60/F1b_h2_residual_level1_zeta')
put('FNWd2', v3['F1a_C3W_deriv_zeta'], 'v3', 'R=60/F1a_C3W_deriv_zeta')
# ---- modulation lemma (final_5e-4.json core)
fc = D['mod']['final']['core']
put('CAB0', fc['ab']['lead'], 'mod', 'final/core/ab/lead'); put('MAB', fc['ab']['m_total'], 'mod', 'final/core/ab/m_total'); put('QAB', fc['ab']['q_total'], 'mod', 'final/core/ab/q_total')
put('PHIAB', fc['ab']['cPhi_total'], 'mod', 'final/core/ab/cPhi_total')
put('BETA0', fc['beta']['lead'], 'mod', 'final/core/beta/lead'); put('BETA1', fc['beta']['m_total'], 'mod', 'final/core/beta/m_total'); put('QB', fc['beta']['q_total'], 'mod', 'final/core/beta/q_total')
put('BETAF', fc['beta']['cPhi_total'], 'mod', 'final/core/beta/cPhi_total')
# dlead (strip-sup feedback coupling into lead) is folded into lead in ibp_final? (lead = K.lead; dlead reported separately): add it conservatively
put('BETA0_dlead', fc['beta']['dlead'], 'mod', 'final/core/beta/dlead'); put('CAB0_dlead', fc['ab']['dlead'], 'mod', 'final/core/ab/dlead')
# ---- U-3 strip forms
u3 = D['u3']
put('STRIP_BETA', u3['STRIP_BETA'], 'u3', 'STRIP_BETA'); put('STRIP_AB', u3['STRIP_AB'], 'u3', 'STRIP_AB')
for c in ('beta', 'ab'):
    for r in range(4): put(f'SF_{c}_{r}', u3[c]['const_coef_S0_S3'][r], 'u3', f'{c}/const_coef_S0_S3/{r}')
    put(f'SF_{c}_al', u3[c]['alpha_coef_S0'], 'u3', f'{c}/alpha_coef_S0')
# ---- C1 cutoff / weight sups
c1 = D['c1']
for nm in ('PSI1', 'two_PSI1', 'PSI2', 'RWR', 'KWIR'): put(nm, c1['weight'][nm], 'c1', 'weight/' + nm)
for nm in ('KAP_IR', 'KAP_Z1', 'KAP_OR', 'KAP'): put(nm, c1['kappa'][nm], 'c1', 'kappa/' + nm)
put('Z1i', c1['zeta']['Z1i'], 'c1', 'zeta/Z1i'); put('Z2i', c1['zeta']['Z2i'], 'c1', 'zeta/Z2i'); put('Z1o', c1['zeta']['Z1o'], 'c1', 'zeta/Z1o')
for j in range(4): put(f'k4_{j+1}', c1['k4'][j], 'c1', f'k4/{j}')
for j in range(3): put(f'k3_{j+1}', c1['k3'][j], 'c1', f'k3/{j}')
for nm in ('T0', 'T2', 'TC0', 'TC2', 'Q2B'): put(nm, c1['transport'][nm], 'c1', 'transport/' + nm)
# ---- C3 modulation norms (C^4)
c3 = D['c3']
for nm in ('NP3', 'NH2', 'NAB', 'NP12', 'C0P'): put(nm, c3[nm], 'c3', nm)
# ---- C2 zone tables, C4 Hardy integrals, per alpha-level
TQ = {}
for lev in ('0.00025', '0.001'):
    T = D['c2']['tables'][lev]
    for zn in ('IR', 'ZA', 'ZB'):
        for kk, v in T[zn]['sup'].items():
            if kk.startswith('sup|'): continue
            put(f'TQ[{lev}][{zn}][{kk}]', v, 'c2', f'tables/{lev}/{zn}/sup/{kk}')
    for kk, v in D['c4'][lev].items(): put(f'HD[{lev}][{kk}]', v, 'c4', f'{lev}/{kk}')
    from flint import arb
    kUi = arb(repr(C[f'TQ[{lev}][IR][kU]']['value'])); assert kUi < arb('0.01')
    put(f'QUF[{lev}]', float((1 / (1 - kUi)).upper()), 'c2', f'1/(1 - tables/{lev}/IR/sup/kU)')     # Q_u, C_u = (q, cos)/(1 - k u): factor on IR
put('C0Hb', D['c4']['C0Hb'], 'c4', 'C0Hb'); put('PHI34_c4', D['c4']['Phi_34_25_upper'], 'c4', 'Phi_34_25_upper')
assert D['c2']['regression_all_equal'], 'C2 regression against r60/planB tables failed'
# ---- test
bad = [n for n, e in C.items() if not (e['value'] >= e['raw'] if e['dir'] == 'up' else e['value'] <= e['raw']) or not math.isfinite(e['value'])]
assert not bad, bad
legacy = {'note': 'E/E2-block and alpha-equation literals of fin12f4 constants.py / rhs_iv.py that EI_SPEC does not cover; carried over unchanged, provenance per ASSEMBLY.md (not re-derived in v4)',
          'names': ['LAMT', 'PI3', 'MOD_ALPHA (env MODA=25)', 'MOD_E2', 'ALIN', 'nd1..nd8', 'nf1..nf7', 't22', 't12', 'r143', 'r163', 'r107', 'f572', 'qb752 factor 7.52', 'QPSI (unused in v4 E_I)', 'cubic_E0fine json (read)', 'ADD_STRIP_CUBIC 2.6e-3']}
out = dict(constants=C, n=len(C), all_outward=True, legacy_typed=legacy, args=args_block(list(P.values())))
h = dump(os.path.join(OUT, 'constants_v4.json'), out)
print(len(C), 'constants; sha256', h)
for n in ('AF', 'CF0', 'CF2', 'CF3', 'SB', 'MU1', 'FN0', 'C3', 'STRIP_BETA', 'STRIP_AB', 'T2', 'NP3'): print(n, C[n]['raw'], C[n]['value'])
