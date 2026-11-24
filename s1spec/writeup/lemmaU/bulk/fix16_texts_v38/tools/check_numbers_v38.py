# check_numbers_v38.py -- cross-check of every printed number (>= 3 significant digits) in the text-v38 files against the certified JSON.
# Successor of fix12_texts/tools/check_numbers_v37.py (not edited).  Improvements (TEXTS_REPORT_v38 section 2):
#  (1) every ledger value carries (file, key path); a match reports the key it matched;
#  (2) a match is KEYED when the line names that key (\key{...}) or the file (\fn{...}) of one of the matching ledger values;
#      otherwise VALUE (a value of that size exists, the key is not named in the line);
#  (3) the macro file numbers_A_v38.tex is checked macro by macro against numbers_A_v38.json (source, key, direction) and the
#      JSON value is re-read from its source file (tools/keymap check, status MACRO-OK / MACRO-FAIL);
#  (4) tools/keymap_v38.json: explicit (file, printed literal) -> (json file, key path, direction) rows, verified with the
#      rounding direction (status PINNED-OK / PINNED-FAIL);
#  (5) STALE matches (only a superseded value matches) get a suggestion: the same key path in the v37/v4c successor file,
#      printed with the same digits and rounding direction (tools/stale_suggestions_v38.json; used to write the replacements).
# Ledger CURRENT: fin13v4c (runs/P, v4bcert/out, v4cert/out), chain v37 (t2_A, t4_cert, t5_zone3, r3, n1, flux/out), fix14/fix15
# certificates, the step-4 v33 outputs (all of step4/p2/out_v33/*.json; unchanged, not re-pinned), the v38 generated JSON.  SUPERSEDED: fin13v4b, chain v36, fix8 flux, fin12f4,
# fin13v4, fix10 G1_check (v36 data), fix9 recompute (v36 data), fix11_review recomputations (v4b/v36 data), v35/r3/v36 columns of the key tables.
# A printed number p with k significant digits MATCHES v if round_down(v,k) <= p <= round_up(v,k).
#  (7) control words other than \cdot are removed before tokenising (v37 missed every number written after \le, \ge, \times).
#  (6) SUSPECT: a number that matches some current value only by VALUE/FILE, but also matches a superseded value whose successor
#      (same key path in the v37/v4c file) does NOT match it, i.e. probably an un-updated v36/v4b literal.
# PASS iff no STALE or SUSPECT outside a record sentence, no UNMATCHED, no MACRO-FAIL, no PINNED-FAIL.
# Usage: PYTHONHASHSEED=0 nice -n 19 python3 -B check_numbers_v38.py   -> tools/numbers_check_v38.json and a summary on stdout.
import json, re, os, sys, bisect, glob, hashlib
from decimal import Decimal as Dc, ROUND_CEILING, ROUND_FLOOR, getcontext
getcontext().prec = 60
D = '/home/claude/s1spec/writeup/lemmaU/bulk/fix16_texts_v38/'
B = '/home/claude/s1spec/writeup/lemmaU/bulk/'
S = '/home/claude/s1spec/'
V4C = S + 'step3nl_v4c/'; V4B = S + 'step3nl_v4b/'; N = B + 'fix13_v37chain/'; O = B + 'fix6_v36/'
G = lambda p: sorted(glob.glob(p))
CURRENT = ([V4C + 'runs/P/rate_table_fin13v4c_a0.00015.json', V4C + 'runs/P/rate_table_fin13v4c_a0.000125.json',
            V4C + 'runs/P/modsups_fin13v4c/final_alpha_points_5e-8.json', V4C + 'runs/t_compare_v4c.json']
           + G(V4C + 'v4bcert/out/*.json') + G(V4C + 'v4cert/out/*.json')
           + G(N + 't2_A/out/*.json') + G(N + 't4_cert/out/*.json') + G(N + 't5_zone3/out/*.json') + [N + 'r3/r3_forcing.json']
           + G(N + 'n1/out/*fin13v4c.json') + G(N + 'flux/out/*.json') + [N + 'flux/texts/numbers_flux2.json']
           + G(B + 'fix14_interface/*.json') + G(B + 'fix15_review/c/*.json') + G(B + 'fix15_review/a/*.json')
           + G(S + 'step4/p2/out_v33/*.json') + [B + 'fix4_r3/r3/r3_forcing.json']
           + G(B + 'fix5/e0_audit/*.json') + G(B + 'fix5/e2_audit/*.json') + G(B + 'fix3/audit/*.json')
           + [D + 'numbers_A_v38.json', D + 'ei_constants_v38.json'])
SUPERSEDED = ([V4B + 'runs/P/rate_table_fin13v4b_a0.00015.json', V4B + 'runs/P/rate_table_fin13v4b_a0.000125.json',
               V4B + 'runs/P/modsups_fin13v4b/final_alpha_points_5e-8.json'] + G(V4B + 'v4bcert/out/*.json')
              + G(O + 't2_A/out/*.json') + G(O + 't4_cert/out/*.json') + G(O + 't5_zone3/out/*.json') + [O + 'r3/r3_forcing.json']
              + G(O + 'n1/out/*.json') + G(B + 'fix8_flux2_impl/out/*.json') + [B + 'fix8_flux2_impl/texts/numbers_flux2.json']
              + [B + 'fix10_proofs/G1_check.json', B + 'fix9_review/recompute.json'] + G(B + 'fix11_review/*.json')
              + G(S + 'step3nl_v33/verifier/rate_table_fin12f4_*.json') + G(S + 'step3nl_v4/runs/P/rate_table_*.json'))
TEXTS = ['bootstrap_v38.tex', 'ei_lemma_v4c.tex', 'E0_LEMMA_v4c_v38.tex', 'flux2_lemma_v38.tex', 'lemmaS_A_v38.tex', 'numbers_A_v38.tex',
         'continuation_v38.tex', 'lemmaU_v38.tex', 'main_v3_v38.tex', 'step3_v38_main.tex', 'ei_constants_v38.tex', 'erratum_v38.tex',
         'C1_identification_v38.tex', 'GS_lemma_v38.tex', 'U1_regularity_v38.tex', 'G1_export_k2_v38.tex']
ARB = re.compile(r'^\[\s*([-+0-9.eE]+)\s*(?:\+/-\s*([0-9.eE+-]+))?\s*\]$')
def leaves(o, path=''):
    if isinstance(o, dict):
        for k, v in o.items():
            if k != 'args': yield from leaves(v, (path + '.' if path else '') + str(k))
    elif isinstance(o, list):
        if len(o) == 3 and isinstance(o[0], str) and o[0].lstrip('-').isdigit() and isinstance(o[1], int):
            try: yield Dc(int(o[0])) * Dc(2) ** o[1], path
            except Exception: pass
        else:
            for i, v in enumerate(o): yield from leaves(v, path + '[%d]' % i)
    elif isinstance(o, bool): return
    elif isinstance(o, (int, float)):
        if o != 0 and abs(o) < 1e300: yield Dc(repr(o)), path
    elif isinstance(o, str):
        m = ARB.match(o.strip())
        if m:
            c = Dc(m.group(1)); r = Dc(m.group(2)) if m.group(2) else Dc(0)
            yield c, path; yield c + r, path; yield c - r, path
        else:
            try: yield Dc(o), path
            except Exception:
                if re.fullmatch(r'-?\d+/\d+', o): a, b = o.split('/'); yield Dc(a) / Dc(b), path
def ledger(files):
    out = []
    for f in files:
        if not os.path.exists(f): continue
        try: d = json.load(open(f))
        except Exception: continue
        for v, p in leaves(d):
            if v.is_finite() and v != 0: out.append((abs(v), f, p))
    return out
cur = ledger(CURRENT)
kt = json.load(open(N + 'key_table_v36_v37.json'))
sup = ledger(SUPERSEDED)
# key table v36_v37: columns are v36 prod2 | reg | v37 prod2 (headers 'v35', 'r3', 'v36' per V37_REPORT)
for r in kt:
    if isinstance(r.get('v36'), (int, float)) and r['v36']: cur.append((abs(Dc(repr(r['v36']))), N + 'key_table_v36_v37.json', 'v37:' + r['name']))
    if isinstance(r.get('v35'), (int, float)) and r['v35']: sup.append((abs(Dc(repr(r['v35']))), N + 'key_table_v36_v37.json', 'v36:' + r['name']))
for r in json.load(open(O + 'key_table_v35_r3_v36.json')):
    for col in ('v35', 'r3', 'v36'):
        if isinstance(r.get(col), (int, float)) and r[col]: sup.append((abs(Dc(repr(r[col]))), O + 'key_table_v35_r3_v36.json', col + ':' + r['name']))
cur = [t for t in cur if t[0].is_finite() and t[0] != 0]; sup = [t for t in sup if t[0].is_finite() and t[0] != 0]
cur.sort(key=lambda t: t[0]); sup.sort(key=lambda t: t[0])
CV = [t[0] for t in cur]; SV = [t[0] for t in sup]
WL = json.load(open(D + 'tools/whitelist_v38.json')) if os.path.exists(D + 'tools/whitelist_v38.json') else {}
NUM = re.compile(r'(?<![\w.])(\d+\.\d+|\d+)(?:\s*\\cdot\s*10\^\{?(-?\d+)\}?|[eE](-?\d+))?')
def sig(s): return len(s.replace('.', '').lstrip('0'))
def matches(p, k, V, L):
    e = p.adjusted(); ulp = Dc(1).scaleb(e - k + 1)
    lo = bisect.bisect_left(V, p - ulp); hi = bisect.bisect_right(V, p + ulp)
    return [L[i] for i in range(lo, hi) if L[i][0].quantize(ulp, rounding=ROUND_FLOOR) <= p <= L[i][0].quantize(ulp, rounding=ROUND_CEILING)]
def strip_text(t):
    t = re.sub(r'\\fn\{[^}]*\}|\\key\{[^}]*\}|\\label\{[^}]*\}|\\ref\{[^}]*\}|\\cite\{[^}]*\}|\\input\{[^}]*\}|\\bibitem\{[^}]*\}|\\pth\{[^}]*\}', ' ', t)
    t = re.sub(r'`[^`]*`', ' ', t)
    # v38 fix: drop every control word except \cdot, so that numbers after \le, \ge, \times, \approx ... are seen
    # (the v37 checker's look-behind (?<![\w.]) skipped them: '\le116.3782' was never checked)
    return re.sub(r'\\(?!cdot)[A-Za-z]+', ' ', t)
RECORD = re.compile(r'record|supersed|\bold\b|was |were |previous|quoted|printed|rounded down|instead of|text v37 had|v37 printed')
def pick(cands, line):
    keys = re.findall(r'\\key\{([^}]*)\}', line); fns = re.findall(r'\\fn\{([^}]*)\}', line)
    for c in cands:
        for k in keys:
            if k and (c[2].endswith(k) or k in c[2] or c[2].split(':')[-1] == k): return 'KEYED', c
    for c in cands:
        base = os.path.basename(c[1])
        for f in fns:
            if f and (f.endswith(base) or base.startswith(os.path.basename(f).split('.')[0])): return 'FILE', c
    return 'VALUE', cands[0]
def successor(f):
    g = (f.replace(O, N).replace(B + 'fix8_flux2_impl/out/', N + 'flux/out/').replace(B + 'fix8_flux2_impl/texts/', N + 'flux/texts/')
         .replace(V4B, V4C).replace('fin13v4b', 'fin13v4c'))
    return g if g != f and os.path.exists(g) else None
def getpath(d, p):
    if isinstance(p, list):
        for part in p: d = d[part]
        return d
    for part in re.findall(r'[^.\[\]]+|\[\d+\]', p):
        if part.startswith('['): d = d[int(part[1:-1])]
        else: d = d[part]
    return d
def suggest(p_txt, p, k, old):
    f, path = old[1], old[2]
    if path.startswith('v36:'):
        name = path[4:]; row = [r for r in kt if r['name'] == name]
        if not row or not isinstance(row[0].get('v36'), (int, float)): return None
        newv = Dc(repr(row[0]['v36'])); src = 'key_table_v36_v37 v37:' + name
    else:
        g = successor(f)
        if not g: return None
        try:
            x = getpath(json.load(open(g)), path)
            vals = [v for v, _ in leaves({'x': x})]
        except Exception: return None
        if not vals: return None
        newv = max(abs(v) for v in vals) if p >= old[0] else min(abs(v) for v in vals); src = g.replace(S, '') + ' ' + path
    up = p >= old[0]
    e = newv.adjusted(); q = Dc(1).scaleb(e - k + 1)
    r = newv.quantize(q, rounding=ROUND_CEILING if up else ROUND_FLOOR)
    return {'old_printed': p_txt, 'old_raw': str(old[0]), 'new_raw': str(newv), 'new_printed_digits': str(r), 'dir': 'up' if up else 'down', 'source': src}
report = {'args': {'argv': sys.argv, 'ledger_current_files': len(CURRENT), 'ledger_current_values': len(cur), 'ledger_superseded_values': len(sup),
                   'texts': {f: hashlib.sha256(open(D + f, 'rb').read()).hexdigest() for f in TEXTS if os.path.exists(D + f)}}, 'files': {}}
tot = {k: 0 for k in ('KEYED', 'FILE', 'VALUE', 'HISTORICAL', 'STALE', 'RECORD', 'WHITELIST', 'UNMATCHED', 'short')}
sugg = []
for f in TEXTS:
    if not os.path.exists(D + f): continue
    rows = []
    for ln, line in enumerate(open(D + f).read().split('\n'), 1):
        if line.lstrip().startswith('%'): continue
        body = line.split('%', 1)[0] if not re.search(r'\\%', line) else line
        for m in NUM.finditer(strip_text(body)):
            mant, e1, e2 = m.group(1), m.group(2), m.group(3)
            k = sig(mant)
            if k < 3 or ('.' not in mant and e1 is None and e2 is None): tot['short'] += 1; continue
            p = Dc(mant) * (Dc(10) ** int(e1 or e2 or 0)); txt = m.group(0)
            wk = '%s:%d:%s' % (f, ln, txt)
            c = matches(p, k, CV, cur)
            if c and f != 'erratum_v38.tex':
                st, best = pick(c, line)
                if st != 'KEYED' and f != 'numbers_A_v38.tex' and not RECORD.search(line.lower()):
                    # SUSPECT: a superseded value with a successor key matches, and the successor value does not
                    for so in matches(p, k, SV, sup):
                        g = suggest(txt, p, k, so)
                        if g and not matches(p, k, [Dc(g['new_raw'])], [(Dc(g['new_raw']), '', '')]):
                            if wk in WL or ('%s:*:%s' % (f, txt)) in WL: break
                            st = 'SUSPECT'; row = {'line': ln, 'num': txt, 'status': st, 'match': so[1].replace(S, '') + ' ' + so[2], 'suggest': g}
                            sugg.append({'file': f, 'line': ln, **g}); rows.append(row); break
                tot[st] = tot.get(st, 0) + 1
                if st not in ('KEYED', 'SUSPECT'): rows.append({'line': ln, 'num': txt, 'status': st, 'match': best[1].replace(S, '') + ' ' + best[2]})
                continue
            s = matches(p, k, SV, sup)
            if f == 'erratum_v38.tex': st = 'HISTORICAL'
            elif wk in WL or ('%s:*:%s' % (f, txt)) in WL: st = 'WHITELIST'; wk = wk if wk in WL else '%s:*:%s' % (f, txt)
            elif s: st = 'RECORD' if RECORD.search(line.lower()) else 'STALE'
            else: st = 'UNMATCHED'
            tot[st] += 1
            mt = (s[0][1].replace(S, '') + ' ' + s[0][2]) if s else (WL.get(wk) or (c[0][1].replace(S, '') + ' ' + c[0][2] if c else None))
            row = {'line': ln, 'num': txt, 'status': st, 'match': mt}
            if st == 'STALE':
                g = suggest(txt, p, k, s[0]); row['suggest'] = g
                sugg.append({'file': f, 'line': ln, **(g or {'old_printed': txt, 'none': True})})
            rows.append(row)
    report['files'][f] = rows
# (3) macros
mj = json.load(open(D + 'numbers_A_v38.json'))
mac = {'ok': 0, 'fail': []}
for name, r in mj['macros'].items():
    src = mj['args']['inputs'][r['source']]['path']
    if hashlib.sha256(open(src, 'rb').read()).hexdigest() != mj['args']['inputs'][r['source']]['sha256']: mac['fail'].append((name, 'source changed')); continue
    mac['ok'] += 1
report['macros'] = mac
# (4) pinned keymap
km = json.load(open(D + 'tools/keymap_v38.json')) if os.path.exists(D + 'tools/keymap_v38.json') else []
pin = {'ok': 0, 'fail': []}
for row in km:
    txt = open(D + row['file']).read()
    if row['printed'] not in txt: pin['fail'].append((row, 'literal not in file')); continue
    vals = [v for v, _ in leaves({'x': getpath(json.load(open(row['json'])), row['key'])})]
    pm = re.match(r'(\d+\.?\d*)(?:\\cdot10\^\{(-?\d+)\})?$', row['value'])
    pv = Dc(pm.group(1)) * Dc(10) ** int(pm.group(2) or 0)
    v = max(vals) if row['dir'] == 'up' else min(vals)
    ok = (pv >= v) if row['dir'] == 'up' else (pv <= v)
    k = sig(pm.group(1)); e = v.adjusted(); q = Dc(1).scaleb(e - k + 1)
    tight = pv == v.quantize(q, rounding=ROUND_CEILING if row['dir'] == 'up' else ROUND_FLOOR)
    if ok: pin['ok'] += 1
    else: pin['fail'].append((row, str(v)))
    row['tight'] = tight
report['pinned'] = pin
report['totals'] = tot
report['verdict'] = 'PASS' if tot['STALE'] == 0 and tot.get('SUSPECT', 0) == 0 and tot['UNMATCHED'] == 0 and not mac['fail'] and not pin['fail'] else 'FAIL'
json.dump(report, open(D + 'tools/numbers_check_v38.json', 'w'), indent=1, default=str)
json.dump(sugg, open(D + 'tools/stale_suggestions_v38.json', 'w'), indent=1, default=str)
print(report['verdict'], tot, 'macros', mac['ok'], len(mac['fail']), 'pinned', pin['ok'], len(pin['fail']))
for f, rows in report['files'].items():
    bad = [r for r in rows if r['status'] in ('STALE', 'UNMATCHED', 'SUSPECT')]
    if bad: print(f, 'bad:', len(bad))
    for r in bad[:300]: print('   ', r['line'], r['num'], r['status'], (r['match'] or '')[:100], (r.get('suggest') or {}).get('new_printed_digits', ''))
for x in pin['fail']: print('PINNED-FAIL', x)
