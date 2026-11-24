#!/usr/bin/env python3
"""fix15_review part B: flux reassembly v37 checks (read-only on fix13_v37chain/flux).  arb outward.  Writes b_flux.json."""
import os, sys, json, hashlib, filecmp, re
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
from flint import arb, ctx
ctx.prec = 200
B = '/home/claude/s1spec/writeup/lemmaU/bulk'; F = B + '/fix13_v37chain/flux'; F8 = B + '/fix8_flux2_impl'
sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()
Ad = lambda x: (lambda f: arb(f.numerator)/arb(f.denominator))(Fr(float(x)))
Aq = lambda s: (lambda f: arb(f.numerator)/arb(f.denominator))(Fr(str(s)))
res = {}; fails = []
def chk(n, ok, info=None):
    res[n] = dict(ok=bool(ok), info=info); print(('PASS ' if ok else 'FAIL ') + n, '' if info is None else str(info)[:400])
    if not ok: fails.append(n)
J = [f for f in sorted(os.listdir(F + '/out')) if f.endswith(('.json', '.pkl')) and 'timing' not in f]
for s in ('v37_pass1', 'v37_pass2'):
    d = [f for f in J if not filecmp.cmp(F + '/out/' + f, F + '/snap/' + s + '/out/' + f, shallow=False)]
    chk('flux/out == snap/%s (%d files)' % (s, len(J)), not d, d)
# pins
def pins(d):
    if isinstance(d, dict):
        for k, v in d.items():
            if isinstance(k, str) and k.startswith('/') and isinstance(v, str) and len(v) == 64: yield k, v
            else: yield from pins(v)
    elif isinstance(d, list):
        for x in d: yield from pins(x)
bad = []; old = set(); n = 0; allp = set()
for f in J:
    if not f.endswith('.json'): continue
    j = json.load(open(F + '/out/' + f))
    for p, h in pins(j):
        n += 1; allp.add(p)
        if not os.path.exists(p): bad.append((f, p, 'missing'))
        elif sha(p) != h: bad.append((f, p, 'mismatch'))
        if 'step3nl_v4b' in p or '/fix6_v36/' in p or '/fix8_flux2_impl/' in p: old.add(p)
chk('flux pins match disk (%d pins)' % n, not bad, bad[:10])
chk('no flux pin points to step3nl_v4b / fix6_v36 / fix8', not old, sorted(old)[:10])
res['pinned_dirs'] = sorted({re.sub(r'/[^/]+$', '', p) for p in allp}); print(*res['pinned_dirs'], sep='\n  ')
# recompute margins
c = json.load(open(F + '/out/flux2_cert.json')); r = c['res']
VP = '/home/claude/s1spec/step3nl_v4c/runs/P'
rt = json.load(open(VP + '/rate_table_fin13v4c_a0.00015.json'))
env = rt['args']['env']; chk('verifier run FLUX0/FLUX2 unchanged', env['FLUX0'] == '1.64e-20' and env['FLUX2'] == '8.0e-20', (env['FLUX0'], env['FLUX2']))
FL0 = Aq('1.64e-20'); FL2 = Aq('8.0e-20')
K0 = Ad(r['K0bar']['up']); K2 = Ad(r['K2pt']['up']); F0 = Ad(r['F0']['up'])
ft = rt['tail']['floor_terms']; print('floor_terms', ft)
chk('F0 cert == rate table floor F0 (1 ulp)', abs(r['F0']['up'] - ft['F0']) <= 1e-15*ft['F0'], (r['F0']['up'], ft['F0']))
# F0 = mu (FLUX0 + eps_I FLUX2): find mu, eps_I in data
txt = open('/home/claude/s1spec/step3nl_v4c/verifier/constants.py').read()
m = {k: re.search(r'^%s\s*=\s*([^;#\n]+)' % k, txt, re.M) for k in ('MU', 'EPSI', 'EPS_I')}
res['mu_eps_src'] = {k: (v.group(1).strip() if v else None) for k, v in m.items()}; print(res['mu_eps_src'])
# solve-free check: with (mu, eps) unknown, check combined margin via the identity margin = F0 / (mu K0 + mu eps K2) and slack s = F0 - that
s_lo = Ad(r['slack_s']['dn']); muK = F0 - s_lo              # = mu K0 + mu eps K2 (upper)
res['margin_combined_from_F0_and_slack'] = float((F0/muK).lower()); print('combined via F0/(F0-s):', (F0/muK).str(12))
chk('combined margin stored <= F0/(F0-s)+1e-9', r['margin_combined']['dn'] <= float((F0/muK).upper())*(1 + 1e-9), r['margin_combined']['dn'])
chk('K0 margin FLUX0/K0bar recomputed', abs(float((FL0/K0).lower()) - r['margin_K0_vs_FLUX0']['dn']) < 1e-12, float((FL0/K0).lower()))
# eps_I from combined: (FL0 + e FL2)/(K0 + e K2) = M  -> e = (FL0 - M K0)/(M K2 - FL2)
M = Ad(r['margin_combined']['dn']); e = (FL0 - M*K0)/(M*K2 - FL2); mu = F0/(FL0 + e*FL2)
res['eps_I_implied'] = e.str(10); res['mu_implied'] = mu.str(10); print('implied eps_I', e.str(10), 'mu', mu.str(10))
Cp = Ad(r['Cp']['up']); Xc = Ad(rt['unc']['X_const'])
chk("C6: C'/X_const <= 1e-9", (Cp/Xc).upper() < 1e-9, float((Aq('1e-9')/(Cp/Xc)).lower()))
chk("C1: S' <= s", Ad(r['Sp']['up']) < s_lo, float((s_lo/Ad(r['Sp']['up'])).lower()))
chk('all conditions pass except C8 (None)', all(v['pass_'] for k, v in c['conditions'].items() if k != 'C8') and c['conditions']['C8']['pass_'] is None and not c['failing'])
# compare with fix8
c8 = json.load(open(F8 + '/out/flux2_cert.json')); cmpv = {}
for k in ('K0bar', 'K2pt', 'Cp', 'Sp'): cmpv[k] = (c8['res'][k]['up'], r[k]['up'], r[k]['up']/c8['res'][k]['up'])
for k in ('margin_combined', 'margin_K0_vs_FLUX0', 'slack_s'): cmpv[k] = (c8['res'][k]['dn'], r[k]['dn'], r[k]['dn']/c8['res'][k]['dn'])
for k in c['conditions']:
    if c['conditions'][k].get('margin') is not None: cmpv[k] = (c8['conditions'][k]['margin'], c['conditions'][k]['margin'])
res['fix8_vs_v37'] = cmpv; [print(' ', k, v) for k, v in cmpv.items()]
chk('no margin below fix8 by more than 1e-3 relative', all(v[1] >= v[0]*(1 - 1e-3) for k, v in cmpv.items() if k.startswith(('margin', 'C', 'slack')) and k not in ('Cp',)), None)
# regression snapshot vs fix8 out (ignoring args/timing)
def strip(d):
    if isinstance(d, dict): return {k: strip(v) for k, v in d.items() if k not in ('args', 'runtime_s', 'timing', 'wall_s')}
    if isinstance(d, list): return [strip(x) for x in d]
    return d
dr = []
for f in J:
    if f.endswith('.json') and os.path.exists(F8 + '/out/' + f):
        if strip(json.load(open(F + '/snap/reg_v36/out/' + f))) != strip(json.load(open(F8 + '/out/' + f))): dr.append(f)
chk('flux regression snapshot == fix8 out apart from args/timing', not dr, dr)
json.dump(dict(res=res, fails=fails, args=dict(argv=sys.argv, env={k: os.environ.get(k) for k in ('PYTHONHASHSEED',)}, prec=ctx.prec,
    inputs={p: sha(p) for p in (F + '/out/flux2_cert.json', F8 + '/out/flux2_cert.json', VP + '/rate_table_fin13v4c_a0.00015.json', os.path.abspath(__file__))})), open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'b_flux.json'), 'w'), indent=1, default=str)
print('FAILS', fails)
