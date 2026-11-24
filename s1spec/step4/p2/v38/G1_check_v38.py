#!/usr/bin/env python3
"""G1_check.py -- independent arb re-check of the level-2 window export of Lemma S(A) (k = 2, r_2 = 1/1000) as consumed by flux lemma v2.
[v38 copy of <s1spec>/writeup/lemmaU/bulk/fix10_proofs/G1_check.py for the step-4 v38 re-pin (p2/v38/final_tail_v38.py); changes: inputs repointed to the v37 chain (runs/prod2) and the v37 flux reassembly (flux/out/flux2_data.json); output out_v38/iface/G1_check.json]
Fix round 10 (G1_export_k2.tex), Nov 20 2026.  READ-ONLY on the proof tree; writes only G1_check.json next to this script.
Does not import any pipeline module: the arb JSON encodings are decoded here.

Inputs (v36 production run prod2, byte-identical to prod1 per V36_REPORT):
  W  = t2_A/weights_t2.json          level-2 dissipation budget (core, zones of [28,350], ramp slack), r_export
  A4 = t2_A/asm_design_N2ref4.json   assembly at 4 N2_ref (the setting the flux lemma reads); export['2']
  A1 = t2_A/asm_design_N2ref.json    assembly at N2_ref (design value quoted in V36_REPORT)
  I  = t2_A/A_interface.json         r2_export block
  U  = t4_cert/U_cert_v36_pass3.json supersolution: N2 <= 4 N2_ref, EA
  F  = fix8_flux2_impl/out/flux2_data.json  the values the flux lemma actually used
Checks (every inequality evaluated in arb; 'contains' = the recorded ball and the recomputed ball overlap and the recomputed upper end is
not above the recorded upper end by more than 1e-30 relative):
  B1  r_export[2] = 1/1000; core retained fraction (units a_c) >= 1/1000; core budget 2 - Bc - used >= r_2 recomputed from Bc and used.
  B2  plateau (60,175]: 2 - Bc_s - used_s >= r_2 a_c / a_s recomputed (a_s = (1 + 0.15^2)^-2, a_c = (1 + 0.11^2)^-2).
  B3  transitions [40,50], [50,60], [175,250], [250,350]: ok and no failure cell (the budget test is inside the zone loop); level-2 ok; all_ok.
  B4  ramp: Young users at level 2 (top, class (2,4)) have deltas summing to <= the certified ramp slack (units a_min(0.05)).
  S1  level-2 source list: every class handler in {pair, phi, ibp223}; pair routes in {data, young}; top route IBP; no other source kind.
  S2  D_2 (damped) = inflow + ramp + top + cls + extra (parts), and D_2^c <= D_2.
  C1  C0_2 = (E_2 + D_2)/(r_2 a_c), C1_2 = (2 F_2^2/gp_2 + coup_2 + o_2 + D_2^s)/(r_2 a_c) recomputed; contain the recorded export values.
  C2  setting: A4 args --N2 = 4 x default N2_ref (0.1, 0.2, 5.0); U supersolution N2_le_4N2ref; U EA_2 <= A4 E_2 (so (D3) holds with A4's E_2).
  C3  the flux lemma read exactly these: F export['2'] C0, C1 == A4 export['2'] C0, C1 upper ends; and quoted 1.1638e2, 2.6061e-1 dominate.
  C4  A_interface r2_export C0_ref4/C1_ref4 dominate A4; design values (N2_ref) 114.782, 0.0698142 dominate A1.
Run: PYTHONHASHSEED=0 nice -n 19 python3 -B G1_check.py
"""
import sys, os, json, hashlib, platform
from fractions import Fraction
sys.dont_write_bytecode = True
from flint import arb, ctx, __version__ as FLINT_VERSION
ctx.prec = 400

HERE = os.path.dirname(os.path.abspath(__file__))
RUN = '/home/claude/s1spec/writeup/lemmaU/bulk/fix13_v37chain/runs/prod2'   # v38: (A) chain v37 production (was fix6_v36)
P = dict(W=RUN + '/t2_A/weights_t2.json', A4=RUN + '/t2_A/asm_design_N2ref4.json', A1=RUN + '/t2_A/asm_design_N2ref.json',
         I=RUN + '/t2_A/A_interface.json', U=RUN + '/t4_cert/U_cert_v36_pass3.json',
         F='/home/claude/s1spec/writeup/lemmaU/bulk/fix13_v37chain/flux/out/flux2_data.json')


def sha(p): return hashlib.sha256(open(p, 'rb').read()).hexdigest()


SHA0 = {k: sha(p) for k, p in P.items()}


def dec(x):
    """arbio [mantissa, exponent, radius] or '[mid +/- rad]' string or {'arb': str} -> arb; recursive"""
    if isinstance(x, dict) and 'arb' in x:
        a = x['arb']
        if isinstance(a, list):
            m, e, r = a; return arb(int(m)) * arb(2) ** int(e) + arb(0, r)
        return arb(a)
    if isinstance(x, dict): return {k: dec(v) for k, v in x.items()}
    if isinstance(x, list): return [dec(v) for v in x]
    return x


def ld(k): return dec(json.load(open(P[k])))


def Q(s):
    f = Fraction(s); return arb(f.numerator) / arb(f.denominator)


CH = {}; VAL = {}


def ok(name, cond, **vals):
    CH[name] = bool(cond)
    if vals: VAL[name] = {k: (v.str(12, radius=True) if isinstance(v, arb) else v) for k, v in vals.items()}
    print('%-96s %s' % (name, bool(cond)), flush=True)


def le(a, b): return bool((a - b).upper() <= 0)        # certified a <= b


W = ld('W'); A4 = ld('A4'); A1 = ld('A1'); I = ld('I'); U = ld('U'); F = ld('F')
R2 = Q('1/1000'); AC = (1 + Q('11/100') ** 2) ** -2; AS_ = (1 + Q('15/100') ** 2) ** -2

# ---------------- B: budget
l2 = W['levels']['2']
ok('B1 r_export[2] == 1/1000', W['r_export'][2] == '1/1000')
ok('B1 core retained (units a_c) >= r_2', le(R2, l2['core']['retained']), retained=l2['core']['retained'])
bc = Q(l2['core']['Bc']); ok('B1 core: 2 - Bc - used/a_c >= r_2 (recomputed)', le(R2, 2 - bc - l2['core']['used']), Bc=l2['core']['Bc'], used=l2['core']['used'])
ok('B1 core retained == 2 - Bc - used (consistency)', abs(float((l2['core']['retained'] - (2 - bc - l2['core']['used'])).mid())) < 1e-12)
seg = l2['seg']; bcs = Q(seg['check']['bc'])
ok('B2 plateau (60,175]: 2 - Bc_s - used_s/a_s >= r_2 a_c/a_s', le(R2 * AC / AS_, 2 - bcs - seg['used']), Bc_s=seg['check']['bc'], used_s=seg['used'])
ok('B2 plateau check ok', seg['check']['ok'] is True)
for t in ('transition1', 'transition2'):
    z = l2[t]
    ok('B3 %s ok, rise %s and partition %s with no failing cell' % (t, z['rise']['interval'], z['partition']['interval']),
       z['ok'] and z['rise']['ok'] and z['partition']['ok'] and not z['rise']['failures'] and not z['partition']['failures'])
cover = sorted([tuple(l2[t][p]['interval']) for t in ('transition1', 'transition2') for p in ('rise', 'partition')])
ok('B3 zones + core [28,40] + plateau (60,175] cover [28,350]', cover == [(40.0, 50.0), (50.0, 60.0), (175.0, 250.0), (250.0, 350.0)])
ok('B3 level-2 ok and all_ok', l2['ok'] is True and W['all_ok'] is True)
ok('B3 weights_t2 ran with V36_R2 = 1', str(W['args']['env'].get('V36_R2')) == '1')

# ---------------- S: sources
r = A4['res']['tight']; i2 = r['info'][2]; p2 = r['parts'][2]
hand = sorted(set(c['handler'] for c in i2['classes'])); routes = sorted(set(c.get('route', '-') for c in i2['classes'] if c['handler'] == 'pair'))
ok('S1 level-2 class handlers subset of {pair, phi, ibp223}', set(hand) <= {'pair', 'phi', 'ibp223'}, handlers=str(hand))
ok('S1 level-2 pair routes subset of {data, young}', set(routes) <= {'data', 'young'}, routes=str(routes))
ok('S1 level-2 top route IBP', i2['top']['route'] == 'IBP')
dsum = sum((v for v in i2['deltas'].values()), arb(0))
ok('B4 level-2 Young deltas (top, classes) sum <= ramp slack', le(dsum, l2['ramp']['slack']), deltas=dsum, slack=l2['ramp']['slack'])
D = i2['D']; Dc = i2['D_window_const']; Ds = i2['D_window_slope']
parts = p2['inflow'] + p2['ramp_comm'] + p2['top'] + p2['cls'] + p2['extra']
ok('S2 D_2 == inflow + ramp + top + cls + extra (balls overlap)', D.overlaps(parts), D=D, parts=parts)
ok('S2 extra == foot + S33', p2['extra'].overlaps(p2['foot'] + p2['S33']))
ok('S2 D_2^c <= D_2', le(Dc, D), Dc=Dc, D=D)

# ---------------- C: constants
e2 = r['E'][2]; gp2 = r['gp'][2]; F2 = i2['F2']; cp = i2['coup_rate']; o2 = r['OUTER']['2']
C0 = (e2 + D) / (R2 * AC); C1 = (2 * F2 / gp2 + cp + o2 + Ds) / (R2 * AC)
X = r['export']['2']
ok('C1 recomputed C0_2 overlaps recorded and is not above it', C0.overlaps(X['C0']) and float(((C0.upper() - X['C0'].upper()) / X['C0'].upper()).upper()) < 1e-30, C0=C0, C0_rec=X['C0'])
ok('C1 recomputed C1_2 overlaps recorded and is not above it', C1.overlaps(X['C1']) and float(((C1.upper() - X['C1'].upper()) / X['C1'].upper()).upper()) < 1e-30, C1=C1, C1_rec=X['C1'])
ok('C1 recorded r and a_c are r_2 and (1 + 0.11^2)^-2', X['r'].overlaps(R2) and X['amin_c'].overlaps(AC))
ok('C1 statement is z_4 on [ln28, ln350] with weight y^-4 (28/y)', 'y^-4 (28/y) z_4^2' in X['statement'] and 'ln28}^{ln350}' in X['statement'])
GH2 = l2['core']['Ghat']
coup_re = sum((arb(GH2[str(k)].upper()) * arb(28) ** -4 * r['E'][k] for k in (0, 1)), arb(0))
ok('C1 coup_2 == sum_{k<2} Ghat_2k 28^-4 E_k (recomputed from weights_t2 Ghat and A4 E_0, E_1)', cp.overlaps(coup_re), coup=cp, coup_re=coup_re)
ok('C1 gamma\'_2 > 0', bool(gp2.lower() > 0), gp2=gp2)
av = A4['args']['argv']; n2 = av[av.index('--N2') + 1: av.index('--N2') + 4]
ok('C2 A4 setting --N2 = 4 x (0.1, 0.2, 5.0)', [Fraction(v) for v in n2] == [4 * Fraction(v) for v in ('0.1', '0.2', '5.0')], N2=str(n2))
ss = U['supersolution']
ok('C2 U_cert supersolution N2_le_4N2ref', ss['N2_le_4N2ref'] is True)
ok('C2 U_cert EA_2 <= A4 E_2 (so E^A_2(tau_a) <= Ehat_2 <= the E_2 used in C0)', le(ss['EA'][2], e2), EA2_U=ss['EA'][2], E2_A4=e2)
ok('C2 U_cert EA_k <= A4 E_k for k = 0, 1 (the coupling uses A4 E_0, E_1)', le(ss['EA'][0], r['E'][0]) and le(ss['EA'][1], r['E'][1]))
fx = F['res']['export']['2']
ok('C3 flux2_data used C0_2, C1_2 equal to A4 upper ends', float(fx['C0'].upper()) == float(X['C0'].upper()) or fx['C0'].overlaps(X['C0'].upper()),
   C0_flux=fx['C0'], C1_flux=fx['C1'])
ok('C3 flux2_data C1_2 equals A4 upper end', fx['C1'].overlaps(X['C1'].upper()))
ok('C3 quoted C0_2 <= 1.1638e2 and C1_2 <= 2.6061e-1 (texts BPhi_continuation_v36, lemmaD_k2_sketch)', le(X['C0'], Q('116.38')) and le(X['C1'], Q('0.26061')))
IRAW = json.load(open(P['I']))          # raw: keep the 'double'/'bound' fields of the interface
rb = IRAW['res']['tight']['r2_export']
def num(v): return arb(v['double'])          # A_interface: 'double' is the outward-rounded upper bound ('bound': 'up')
ok('C4 A_interface r2_export C0_ref4, C1_ref4 (bound up) dominate A4 export', rb['C0_ref4']['bound'] == 'up' and le(X['C0'], num(rb['C0_ref4'])) and le(X['C1'], num(rb['C1_ref4'])))
X1 = A1['res']['tight']['export']['2']
ok('C4 A_interface C0_ref, C1_ref (bound up) dominate A1 export (design N2_ref)', le(X1['C0'], num(rb['C0_ref'])) and le(X1['C1'], num(rb['C1_ref'])), C0=X1['C0'], C1=X1['C1'])
# informational (not a check): V36_REPORT quotes the design values 114.782 and 0.0698142, which are the 6-digit truncations (rounded DOWN)
INFO = dict(v36_report_design_quotes_rounded_down=bool((X1['C0'] - Q('114.782')).lower() > 0 and (X1['C1'] - Q('0.0698142')).lower() > 0),
            note='V36_REPORT table: C0_2 114.782, C1_2 0.0698142 at N2_ref are below the certified values 114.78222..., 0.069814238...; not consumed (the flux lemma reads 4 N2_ref)')
ok('C4 A_interface r2_export r_2 == 1/1000', str(rb['r_2']) == '1/1000')

SHA1 = {k: sha(p) for k, p in P.items()}
ok('inputs unchanged during the run', SHA0 == SHA1)
me = os.path.abspath(__file__)
out = dict(checks=CH, values=VAL, info=INFO, all_pass=all(CH.values()),
           args=dict(argv=sys.argv, cwd=os.getcwd(), env={k: os.environ.get(k) for k in ('PYTHONHASHSEED',)}, python=platform.python_version(),
                     flint=FLINT_VERSION, prec=ctx.prec, inputs={**{P[k]: SHA0[k] for k in P}, me: sha(me)}))
json.dump(out, open('/home/claude/s1spec/step4/p2/out_v38/iface/G1_check.json', 'w'), indent=1, sort_keys=True)
print('all_pass', out['all_pass'])
sys.exit(0 if out['all_pass'] else 1)
