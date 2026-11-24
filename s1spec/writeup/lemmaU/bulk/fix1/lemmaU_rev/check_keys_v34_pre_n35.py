#!/usr/bin/env python3
"""Checks for lemmaU_v34.tex (read-only; writes nothing):
 1. sha256 of every source in constants_needed_v34.json;
 2. every key of every row resolves in its JSON file;
 3. every \\key{...} cited in lemmaU_v34.tex resolves in at least one listed row's file (coverage);
 4. every quoted number of the 'numbers' list is rounded outward against the certified value (exact rational comparison);
 5. the derived side conditions used by the text (domination by the 4 x N2_ref run, (A)'s core |z1|, D^c <= D^src, K6, B2', pins).
usage: nice -n 19 python3 -I check_keys_v34.py   (exit status 1 on any failure)"""
import json, hashlib, re, os, sys
from fractions import Fraction as Fr
from decimal import Decimal
HERE = os.path.dirname(os.path.abspath(__file__))
CN = json.load(open(os.path.join(HERE, 'constants_needed_v34.json')))
TEX = open(os.path.join(HERE, 'lemmaU_v34.tex')).read()
fails = []; notes = []
def fail(m): fails.append(m); print('FAIL', m)
D = {}
for k, f in CN['files'].items():
    h = hashlib.sha256(open(f['path'], 'rb').read()).hexdigest()
    if h != f['sha256']: fail('sha256 %s %s' % (k, f['path']))
    D[k] = json.load(open(f['path']))
    print('file %-3s %s  %s' % (k, h[:8], f['path']))
def resolve(d, path):
    """path components separated by '.', but dict keys may themselves contain '.' or ' ' -> greedy longest match."""
    parts = path.split('.') if path else []
    i = 0; cur = d
    while i < len(parts):
        if isinstance(cur, list):
            cur = cur[int(parts[i])]; i += 1; continue
        if not isinstance(cur, dict): raise KeyError(path)
        for j in range(len(parts), i, -1):
            k = '.'.join(parts[i:j])
            if k in cur: cur = cur[k]; i = j; break
        else: raise KeyError('%s (at %s)' % (path, '.'.join(parts[i:])))
    return cur
# 2. rows
nk = 0
for r in CN['constants']:
    if 'file' not in r: continue
    for k in r['key'].split(' / '):
        try: resolve(D[r['file']], k); nk += 1
        except Exception as e: fail('row key %s:%s (%s)' % (r['file'], k, e))
print('row keys resolved:', nk)
# 3. coverage of \key{} in the text
ALIAS = {'.C1': ['U:supersolution.C1'], '.margin_FLUX2': ['A4:res.tight.flux.beta.2.margin_FLUX2'], 'launch.*': ['U:launch'],
         'supersolution.K1_ok,K2_ok,K3_ok': ['U:supersolution.K1_ok', 'U:supersolution.K2_ok', 'U:supersolution.K3_ok'],
         'margins.present_rows[].Ubar': ['U:margins.present_rows.0.Ubar'], 'res.tight.': ['A:res.tight']}
rowkeys = {}
for r in CN['constants']:
    if 'file' in r:
        for k in r['key'].split(' / '): rowkeys.setdefault(r['file'], set()).add(k)
cited = sorted(set(re.findall(r'\\key\{([^}]*)\}', TEX)))
uncovered = []
for c in cited:
    cands = ALIAS.get(c)
    ok = False
    if cands:
        ok = all(any(resolve(D[x.split(':')[0]], x.split(':', 1)[1]) is not None for _ in [0]) for x in cands)
    else:
        for f, ks in rowkeys.items():
            for pre in ('', 'res.tight.'):
                if pre + c in ks: ok = True
        if not ok:   # sub-key of a listed row key (e.g. a component the text names)
            for f, ks in rowkeys.items():
                for pre in ('', 'res.tight.'):
                    try:
                        resolve(D[f], pre + c)
                        if any((pre + c).startswith(k) for k in ks): ok = True
                    except Exception: pass
    if not ok: uncovered.append(c)
for c in uncovered: fail('cited key not covered by a row: %s' % c)
print('cited keys:', len(cited), 'uncovered:', len(uncovered))
# 4. quoted numbers, outward
def interval(v):
    if isinstance(v, dict) and 'arb' in v:
        a = v['arb']
        if isinstance(a, list):
            m = Fr(int(a[0])) * (Fr(2) ** a[1]); r = Fr(a[2]); return m - r, m + r
        s = a.strip('[]')
        if '+/-' in s:
            m, r = s.split('+/-'); m = Fr(Decimal(m.strip())); r = Fr(Decimal(r.strip())); return m - r, m + r
        x = Fr(Decimal(s)); return x, x
    if isinstance(v, (int, float)) and not isinstance(v, bool): x = Fr(v); return x, x
    if isinstance(v, str): x = Fr(Decimal(v)); return x, x
    raise TypeError(repr(v)[:80])
nn = 0
for q, f, k, d in CN['numbers']:
    lo, hi = interval(resolve(D[f], k)); Q = Fr(Decimal(q))
    ok = Q >= hi if d == 'up' else Q <= lo
    nn += 1
    if not ok: fail('number %s %s:%s certified [%.10g, %.10g]' % (q, f, k, float(lo), float(hi)))
print('quoted numbers checked:', nn)
# 5. side conditions
T = 'res.tight.'
def up(f, k): return interval(resolve(D[f], k))[1]
def lo(f, k): return interval(resolve(D[f], k))[0]
for i in range(3):
    if not lo('A', T + 'EA_ref4.%d' % i) >= up('U', 'supersolution.EA.%d' % i): fail('EA_ref4 does not dominate Ehat at %d' % i)
    if not up('A4', T + 'E.%d' % i) == up('A', T + 'EA_ref4.%d' % i) and abs(float(up('A4', T + 'E.%d' % i)) / float(up('A', T + 'EA_ref4.%d' % i)) - 1) > 1e-12:
        fail('A4 E differs from A EA_ref4 at %d' % i)
    n2ref = Fr(Decimal(D['A']['N2_ref']['values'][i]))
    if not up('U', 'supersolution.N2.%d' % i) <= 4 * n2ref: fail('Nbar^2 > 4 N2_ref at %d' % i)
if resolve(D['U'], 'supersolution.N2_le_4N2ref') is not True: fail('N2_le_4N2ref not true')
S1 = up('A4', T + 'S1')
if not 2 * S1 * S1 <= Fr(1475, 10000) ** 2: fail('sqrt2 S1v > 0.1475')
print('(A) core |z1| <= sqrt2 * %.6f = %.6f' % (float(S1), 2 ** 0.5 * float(S1)))
for i in range(2):
    if not up('A4', T + 'info.%d.D_window_const' % i) <= lo('A', T + 'Dsrc.%d' % i): fail('D^c > D^src at %d' % i)
    if not up('A4', T + 'info.%d.D_window_slope' % i) <= up('A', T + 'B0.D_window_slope.%d' % i) * (1 + Fr(1, 10**12)): fail('D^s differs A/A4 at %d' % i)
k6 = lo('U', 'supersolution.EA.0') / up('A', T + 'IC.E_tau1.0')
if not k6 >= Fr(1190, 1000): fail('K6 margin %.5f < 1.190' % float(k6))
for i in range(3):
    if not up('A', T + 'IC.E_tau1.%d' % i) < lo('U', 'supersolution.EA.%d' % i): fail('K6 at %d' % i)
b2 = Fr(153, 10000) / (up('A4', T + 'S0') + up('A', T + 'B0.sup_core.0'))
if not b2 >= Fr(2116, 1000): fail("B2' S0v+sup|d| margin %.5f" % float(b2))
print("K6 margin %.5f, B2' margin (S0v + sup_core|d| vs 0.0153) %.5f" % (float(k6), float(b2)))
pins = D['U']['args']['inputs']
for k, f in (('A', 'A'), ('Z', 'Z'), ('U1', 'U1')):
    if pins.get(CN['files'][f]['path']) != CN['files'][f]['sha256']: fail('U_cert pin of %s' % f)
if D['R'].get('u1_ghat_sha') != CN['files']['U1']['sha256']:
    notes.append('ramp_U.json pins u1_ghat.json %s..., current %s... (stale pin, gap O3)' % (D['R']['u1_ghat_sha'][:8], CN['files']['U1']['sha256'][:8]))
if [float(x) for x in D['R']['Y_next_primary_used']] != [x['double'] for x in D['U1']['cases']['primary']['Y_next']]:
    notes.append('ramp_U Y_next_primary_used differs from u1 primary Y_next as doubles (check)')
s2 = D['S4']['results']
for key in ('H2f', 'H2_tightened_fU1', 'cone_fixedL', 'H4_fixedL', 'H1_fixedL', 'H7_fixedL'):
    if s2[key].get('ok') is not True: fail('S4 %s not ok' % key)
for n in notes: print('NOTE', n)
print('RESULT:', 'PASS' if not fails else 'FAIL (%d)' % len(fails))
sys.exit(1 if fails else 0)
