#!/usr/bin/env python3
"""fix15_review part B (adapted from fix11_review/b_check.py): independent checks of the (A) chain v37 (read-only on fix13_v37chain).  arb, outward.  Writes b_check.json next to itself."""
import os, sys, json, hashlib, decimal, math, glob
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
from flint import arb, ctx
ctx.prec = 200
N = '/home/claude/s1spec/writeup/lemmaU/bulk/fix13_v37chain'
VP = '/home/claude/s1spec/step3nl_v4c/runs/P'
RT = VP + '/rate_table_fin13v4c_a0.00015.json'; RT125 = VP + '/rate_table_fin13v4c_a0.000125.json'
MS = VP + '/modsups_fin13v4c/final_alpha_points_5e-8.json'
sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()
def Aq(s):
    f = Fr(str(s)); return arb(f.numerator)/arb(f.denominator)
def Ad(x):  # exact double
    f = Fr(float(x)); return arb(f.numerator)/arb(f.denominator)
def val(x):
    """certified-number export -> arb (several formats)"""
    if isinstance(x, dict):
        if 'arb' in x and isinstance(x['arb'], list):
            m, e, r = x['arb'][:3]; return None
        for k in ('up', 'double', 'value'):
            if k in x: return Ad(x[k])
    return Ad(x)
res = {}; fails = []
def chk(name, ok, info=None):
    res[name] = dict(ok=bool(ok), info=info)
    if not ok: fails.append(name)
    print(('PASS ' if ok else 'FAIL ') + name, '' if info is None else str(info)[:300])
# ---------------- 1. wiring: argv of lemmaS_F10p vs raw verifier JSON, independent rounding
rt = json.load(open(RT)); rt125 = json.load(open(RT125)); ms = json.load(open(MS)); u = rt['unc']
a = json.load(open(N + '/t2_A/base_ibp/lemmaS_F10p_M0.52.json'))['args']
print('F10p args keys', sorted(a.keys()) if isinstance(a, dict) else a)
D = decimal.Decimal
def exact(x): return D(repr(float(x))) if False else D(float(x))      # exact binary value
def sig_ok(used, raw, sig, up=True):
    du = D(str(used)); dr = D(float(raw))
    ok = (du >= dr) if up else (du <= dr)
    rel = abs(du - dr)/dr
    return ok and rel < D(10)**(1 - sig), float(rel)
# lemmaS_F2 stores its parsed argv (floats); the decimal strings typed by v36cfg are the shortest reprs of these floats
J = [repr(x) for x in a['J']]; SL = [repr(x) for x in a['slopes']]; J4 = repr(a['J4']); SJ4 = repr(a['sJ4']); SQ = repr(a['sqLam']); RHO = repr(a['rho']); AB = repr(a['abar']); BE = repr(a['beta'])
print('consumed:', J, SL, J4, SJ4, SQ, RHO, AB, BE, 'alpha1', a['alpha1'], 'handoff', a['handoff'])
for m, nm in enumerate(('V0sq', 'V1sq', 'V2sq')):
    chk('J%d up' % m, *sig_ok(J[m], u[nm]['const'], 5)); chk('slope%d up' % m, *sig_ok(SL[m], u[nm]['slope_per_tau'], 5))
chk('J4 up', *sig_ok(J4, u['V4sq']['const'], 5)); chk('sJ4 up', *sig_ok(SJ4, u['V4sq']['slope_per_tau'], 5))
chk('sqLam up', *sig_ok(SQ, rt['sqrtL_max'], 5)); chk('rho down', *sig_ok(RHO, rt['tail']['rho'], 4, up=False))
chk('abar == mod_sups', D(AB) == D(str(ms['1.5e-4']['sup_ab'])), AB); chk('beta == mod_sups', D(BE) == D(str(ms['1.5e-4']['sup_beta'])), BE)
# worst of the two alpha points
w = []
for nm in ('V0sq', 'V1sq', 'V2sq', 'V4sq'):
    w.append(rt125['unc'][nm]['const'] <= u[nm]['const'] and rt125['unc'][nm]['slope_per_tau'] <= u[nm]['slope_per_tau'])
w += [rt125['sqrtL_max'] <= rt['sqrtL_max'], rt125['tail']['rho'] >= rt['tail']['rho'], ms['1.25e-4']['sup_ab'] <= ms['1.5e-4']['sup_ab'], ms['1.25e-4']['sup_beta'] <= ms['1.5e-4']['sup_beta'],
      ms['1.25e-4']['M1'] <= ms['1.5e-4']['M1'], ms['1.25e-4']['sI1'] <= ms['1.5e-4']['sI1'], ms['1.25e-4']['sI2'] <= ms['1.5e-4']['sI2'], ms['1.25e-4']['S_ramp'] <= ms['1.5e-4']['S_ramp']]
chk('1.5e-4 row dominates the 1.25e-4 row (J, slopes, peak, rho, sup_ab, sup_beta, M1, sI1, sI2, S_ramp)', all(w), w)
# mod_sups values vs rate table (is sup_ab/sup_beta derivable?  only consistency of J/peak/rho here)
chk('mod_sups J/peak/rho == rate table', ms['1.5e-4']['peak'] == rt['sqrtL_max'] and ms['1.5e-4']['rho'] == rt['tail']['rho'] and all(ms['1.5e-4']['J'][k] == [u[k]['const'], u[k]['slope_per_tau']] for k in ('V0sq', 'V1sq', 'V2sq', 'V4sq')))
chk('rate table env FLUX0/FLUX2/S0V/S1V/V4B', [rt['args']['env'].get(k) for k in ('FLUX0', 'FLUX2', 'S0V', 'S1V', 'V4B', 'TAG')] == ['1.64e-20', '8.0e-20', '0.0153', '0.6', '1', 'fin13v4c'], [rt['args']['env'].get(k) for k in ('FLUX0', 'FLUX2', 'S0V', 'S1V', 'V4B', 'TAG', 'MODA', 'ETA')])
# ---------------- 2. sha pins in every args block of the in-place outputs
bad = []; npin = 0; files = glob.glob(N + '/t2_A/out/*.json') + glob.glob(N + '/t4_cert/out/*.json') + glob.glob(N + '/t5_zone3/out/*.json') + glob.glob(N + '/t2_A/base_ibp/*.json') + glob.glob(N + '/t2_A/base_ibp/rigor/*.json') + [N + '/r3/r3_forcing.json', N + '/n1/out/n1_chart_fin13v4c.json']
verif_seen = set()
def pins(d):
    if isinstance(d, dict):
        for k, v in d.items():
            if isinstance(k, str) and k.startswith('/') and isinstance(v, str) and len(v) == 64: yield k, v
            else: yield from pins(v)
    elif isinstance(d, list):
        for x in d: yield from pins(x)
noargs = []
for f in sorted(files):
    j = json.load(open(f))
    if not (isinstance(j, dict) and 'args' in j): noargs.append(os.path.relpath(f, N)); continue
    for p, h in pins(j['args']):
        npin += 1
        if 'step3nl' in p: verif_seen.add(p)
        if not os.path.exists(p): bad.append((os.path.relpath(f, N), p, 'missing'))
        elif sha(p) != h: bad.append((os.path.relpath(f, N), p, 'mismatch'))
chk('all sha256 pins in args blocks match disk (%d pins, %d files)' % (npin, len(files)), not bad, bad[:12])
res['files_without_args'] = noargs; print('no args block:', noargs)
res['verifier_paths_pinned'] = sorted(verif_seen); print('verifier paths pinned:', *sorted(verif_seen), sep='\n   ')
chk('no pinned verifier path outside step3nl_v4b (or verifier-independent step3nl modules)', all('/step3nl_v4c/' in p or '/step3nl/' in p for p in verif_seen), sorted(p for p in verif_seen if '/step3nl_v4c/' not in p))
# ---------------- 3. kappa_J recomputed
h = json.load(open(N + '/t2_A/out/hjw_rep.json'))
Ly1 = Ad(rt['Ly1']); Xc = Ad(u['X_const']); F0 = Ad(rt['tail']['floor_terms']['F0'])
exc = max((Ad(s['Ly_hi']) - Ly1 for s in rt['steps']), key=lambda z: float(z.mid())); Ds = sum((Ad(s['D']) for s in rt['steps']), arb(0))
kJ = (exc/2 + F0*Ds)/Xc
chk('kappa_J certified reproduced', abs(float(kJ.upper()) - h['kappa_J_certified']) < 1e-15 and float(kJ.upper()) <= h['kappa_J_certified']*(1 + 1e-15), (float(kJ.upper()), h['kappa_J_certified']))
chk('kappa_J used - certified >= 1e-9 (1 ulp)', Ad(h['kappa_J']) - kJ >= Aq('1e-9') - Aq('1e-17'), float((Ad(h['kappa_J']) - kJ).lower()))
Cp = Ad(json.load(open(N + '/flux/out/flux2_cert.json'))['Cprime']) if 'Cprime' in json.load(open(N + '/flux/out/flux2_cert.json')) else Aq('2.65232e-26')   # [read] fix8 C' (integrated constant), FLUX2_IMPL_REPORT / fix9 recompute 2.65224e-26
chk("C'/X_const <= 1e-9 pad", (Cp/Xc).upper() < Aq('1e-9').lower(), 'margin x%.1f' % float((Aq('1e-9')/(Cp/Xc)).lower()))
# F0 = MU*FLUX0 + MU*EPSI*FLUX2 consistency is the verifier's; X_const parts
px = u['parts_Xc']; chk('X_const >= sum of its parts, excess < 1e-9 relative (verifier-side rounding; info)', 0 <= u['X_const'] - sum(px.values()) < 1e-9*u['X_const'], (u['X_const'] - sum(px.values()))/u['X_const'])
# ---------------- 4. R3 constants and c_dU
CU = Aq('1.3056'); P = Aq('1.0246'); LT = Aq('-0.2561532224577463'); LAMT = Aq('-0.2561532')
DLq = max(abs(Fr('-0.2561532224577463') + Fr('1.0246')/4), abs(Fr('-0.2561532') + Fr('1.0246')/4)); DL = arb(DLq.numerator)/arb(DLq.denominator)
AL = Aq('1.5e-4'); ABV = Aq(str(ms['1.5e-4']['sup_ab'])); BBV = Aq(str(ms['1.5e-4']['sup_beta']))
cdU = CU*(BBV + DL*AL + P*ABV*AL)/(P*(arb(1)/4 - ABV))
u34 = json.load(open(N + '/t4_cert/out/u34_primary.json'))
print('u34 launch keys', list(u34['launch'].keys())[:20]); print('c_dU stored', u34['launch']['c_dU'])
res['c_dU_mine'] = cdU.str(12); print('c_dU mine', cdU.str(15), ' 2e-8/c_dU', float((Aq('2e-8')/cdU).lower()))
chk('c_dU < 2e-8', cdU.upper() < Aq('2e-8').lower(), cdU.str(12))
const_txt = open('/home/claude/s1spec/step3nl_v4c/verifier/constants.py').read()
chk('verifier LAMT literal == -0.2561532', 'LAMT = -0.2561532;' in const_txt)
r3 = json.load(open(N + '/r3/r3_forcing.json'))
r0 = CU*(Aq(BE) + DL*AL) + CU*P*AL*(Aq(AB) + Aq(AB)/14)
res['r0_mine'] = r0.str(12); print('r0 mine (no ulp inflation)', r0.str(15), ' stored', r3['r0_new'], ' c_dU stored r3', r3['c_dU']['new'])
# ---------------- 5. r_2 export recomputation
AI = json.load(open(N + '/t2_A/out/A_interface.json')); r2 = AI['res']['tight']['r2_export']
print('r2_export', json.dumps({k: v for k, v in r2.items() if k not in ('note', 'statement')})[:1500])
res['r2_export'] = {k: v for k, v in r2.items() if k not in ('note',)}

# ---------------- 6. v37 extras: stale-output check, snapshot identity, verifier-path audit, flux
import filecmp
V36P = '/home/claude/s1spec/writeup/lemmaU/bulk/fix6_v36'
def walk(root):
    for d, _, fs in os.walk(root):
        for f in fs: yield os.path.relpath(os.path.join(d, f), root)
p2 = N + '/runs/prod2'; p1 = N + '/runs/prod1'
same12 = [r for r in walk(p2) if r.endswith('.json') and not filecmp.cmp(p2 + '/' + r, p1 + '/' + r, shallow=False)]
chk('prod1 == prod2 for every JSON (bytes)', not same12, same12)
# in-place outputs equal snapshot prod2
mp = {'t2_A': 't2_A/out', 't4_cert': 't4_cert/out', 't5_zone3': 't5_zone3/out', 'base_ibp': 't2_A/base_ibp', 'r3': 'r3', 'n1': 'n1/out'}
dif = []
for r in walk(p2):
    top = r.split('/')[0]
    if top in mp and r.endswith('.json'):
        q = N + '/' + mp[top] + r[len(top):]
        if not os.path.exists(q) or not filecmp.cmp(p2 + '/' + r, q, shallow=False): dif.append(r)
chk('in-place outputs == runs/prod2', not dif, dif)
# files byte-identical between v37 prod2 and v36 prod2 (candidates for stale / verifier-independent)
ident = [r for r in walk(p2) if r.endswith('.json') and os.path.exists(V36P + '/runs/prod2/' + r) and filecmp.cmp(p2 + '/' + r, V36P + '/runs/prod2/' + r, shallow=False)]
res['json_identical_to_v36_prod2'] = ident; print('identical to v36 prod2:', ident)
# in-place output JSON not in snapshot (possible stale leftovers from the rsync copy), with mtime
t0 = os.path.getmtime(N + '/runs/prod2/env.txt') - 500
stale = []
for sub in ('t2_A', 't4_cert', 't5_zone3', 'n1', 'r3'):
    for r in walk(N + '/' + sub):
        if r.endswith(('.json', '.pkl', '.npz')) and os.path.getmtime(N + '/' + sub + '/' + r) < t0: stale.append(sub + '/' + r)
res['data_files_older_than_prod2_start'] = stale; print('older than prod2:', len(stale)); [print('   ', x) for x in stale]
# any pinned path anywhere in prod2 pointing to v4b / fix6_v36
badp = []
for r in walk(p2):
    if r.endswith('.json'):
        txt = open(p2 + '/' + r).read()
        for needle in ('step3nl_v4b', 'bulk/fix6_v36', 'fin13v4b'):
            if needle in txt: badp.append((r, needle, txt.count(needle)))
res['prod2_mentions_of_v4b_or_fix6'] = badp; print('mentions v4b/fix6:', badp)
json.dump(dict(res=res, fails=fails, args=dict(argv=sys.argv, env={k: os.environ.get(k) for k in ('PYTHONHASHSEED', 'PYTHONDONTWRITEBYTECODE')}, prec=ctx.prec,
          inputs={p: sha(p) for p in (RT, RT125, MS, N + '/t2_A/out/hjw_rep.json', N + '/t2_A/out/A_interface.json', N + '/t4_cert/out/u34_primary.json', N + '/r3/r3_forcing.json', os.path.abspath(__file__))})),
          open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'b_check.json'), 'w'), indent=1, default=str)
print('FAILS:', fails)
