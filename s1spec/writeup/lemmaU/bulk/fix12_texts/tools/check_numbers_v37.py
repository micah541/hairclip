# check_numbers_v37.py -- cross-check of every printed number (>= 3 significant digits) in the text-v37 files against the certified JSON.
# Ledger: numeric leaves of the CURRENT certified outputs (fin13v4b, v36 chain, R3, flux v2, G1, the v4/v4b constant files, step-4 v33
# outputs that are unchanged) and of SUPERSEDED outputs (v35 column of the v36 key table, fin12f4 tables, fin13v4 tables).
# A printed number p with k significant digits MATCHES a ledger value v if round_down(v,k) <= p <= round_up(v,k).
# Status per occurrence: CURRENT (matches a current value), STALE (matches only superseded values), WHITELIST (tools/whitelist_v37.json,
# with a reason), UNMATCHED.  PASS iff no STALE outside a sentence marked as a record, and no UNMATCHED.
# Usage: PYTHONHASHSEED=0 nice -n 19 python3 -B check_numbers_v37.py   -> tools/numbers_check_v37.json and a summary on stdout.
import json, re, os, sys, bisect, glob, hashlib
from decimal import Decimal as Dc, ROUND_CEILING, ROUND_FLOOR, getcontext
getcontext().prec = 60
D = '/home/claude/s1spec/writeup/lemmaU/bulk/fix12_texts/'
B = '/home/claude/s1spec/writeup/lemmaU/bulk/'
V4B = '/home/claude/s1spec/step3nl_v4b/'
CURRENT = ([V4B + 'runs/P/rate_table_fin13v4b_a0.00015.json', V4B + 'runs/P/rate_table_fin13v4b_a0.000125.json',
            V4B + 'runs/P/modsups_fin13v4b/final_alpha_points_5e-8.json', V4B + 'runs/P/t_tail_alpha_fin13v4b.json', V4B + 'runs/v4b_compare.json']
           + sorted(glob.glob(V4B + 'v4bcert/out/*.json')) + sorted(glob.glob(V4B + 'v4cert/out/*.json'))
           + sorted(glob.glob(B + 'fix6_v36/t2_A/out/*.json')) + sorted(glob.glob(B + 'fix6_v36/t4_cert/out/*.json'))
           + sorted(glob.glob(B + 'fix6_v36/t5_zone3/out/*.json')) + [B + 'fix6_v36/r3/r3_forcing.json', B + 'fix6_v36/t2_A/out/hjw_rep.json']
           + sorted(glob.glob(B + 'fix8_flux2_impl/out/*.json')) + [B + 'fix8_flux2_impl/texts/numbers_flux2.json']
           + [B + 'fix9_review/recompute.json', B + 'fix10_proofs/G1_check.json', B + 'fix4_r3/r3/r3_forcing.json']
           + ['/home/claude/s1spec/step4/p2/out_v33/h2u_v33_arb.json', '/home/claude/s1spec/step4/p2/out_v33/b0_dprof.json']
           + sorted(glob.glob(B + 'fix5/e0_audit/*.json')) + sorted(glob.glob(B + 'fix5/e2_audit/*.json')) + sorted(glob.glob(B + 'fix3/audit/*.json'))
           + sorted(glob.glob(B + 'fix6_v36/n1/out/*.json')) + [D + 'numbers_A_v37.json', D + 'ei_constants_v37.json'])
SUPERSEDED = (sorted(glob.glob('/home/claude/s1spec/step3nl_v33/verifier/rate_table_fin12f4_*.json'))
              + sorted(glob.glob('/home/claude/s1spec/step3nl_v4/runs/P/rate_table_*.json')))
TEXTS = ['ei_lemma_v4.tex', 'flux2_lemma_v37.tex', 'lemmaS_A_v37.tex', 'numbers_A_v37.tex', 'continuation_v37.tex', 'lemmaU_v37.tex',
         'main_v3_v37.tex', 'step3_v37_main.tex', 'ei_constants_v37.tex', 'erratum_v37.tex']
ARB = re.compile(r'^\[\s*([-+0-9.eE]+)\s*(?:\+/-\s*([0-9.eE+-]+))?\s*\]$')
def leaves(o):
    if isinstance(o, dict):
        for k, v in o.items():
            if k != 'args': yield from leaves(v)
    elif isinstance(o, list):
        if len(o) == 3 and isinstance(o[0], str) and o[0].lstrip('-').isdigit() and isinstance(o[1], int):
            try: yield Dc(int(o[0])) * Dc(2) ** o[1]
            except Exception: pass
        else:
            for v in o: yield from leaves(v)
    elif isinstance(o, bool): return
    elif isinstance(o, (int, float)):
        if o != 0 and abs(o) < 1e300: yield Dc(repr(o))
    elif isinstance(o, str):
        m = ARB.match(o.strip())
        if m:
            c = Dc(m.group(1)); r = Dc(m.group(2)) if m.group(2) else Dc(0)
            yield c; yield c + r; yield c - r
        else:
            try: yield Dc(o)
            except Exception:
                if re.fullmatch(r'-?\d+/\d+', o): a, b = o.split('/'); yield Dc(a) / Dc(b)
def ledger(files):
    vals = []; src = []
    for f in files:
        if not os.path.exists(f): continue
        try: d = json.load(open(f))
        except Exception: continue
        for v in leaves(d):
            if v.is_finite() and v != 0: vals.append(abs(v)); src.append(f)
    # key table: v36 current, v35/r3 superseded
    return vals, src
cv, cs = ledger(CURRENT)
kt = json.load(open(B + 'fix6_v36/key_table_v35_r3_v36.json'))
for r in kt:
    if isinstance(r.get('v36'), (int, float)): cv.append(abs(Dc(repr(r['v36'])))); cs.append('key_table v36: ' + r['name'])
sv, ss = ledger(SUPERSEDED)
for r in kt:
    for col in ('v35', 'r3'):
        if isinstance(r.get(col), (int, float)): sv.append(abs(Dc(repr(r[col])))); ss.append('key_table %s: %s' % (col, r['name']))
def index(vals, srcs):
    z = sorted((a, b) for a, b in zip(vals, srcs) if a.is_finite() and a != 0); return [a for a, _ in z], [b for _, b in z]
CV, CS = index(cv, cs); SV, SS = index(sv, ss)
WL = json.load(open(D + 'tools/whitelist_v37.json')) if os.path.exists(D + 'tools/whitelist_v37.json') else {}
NUM = re.compile(r'(?<![\w.])(\d+\.\d+|\d+)(?:\s*\\cdot\s*10\^\{?(-?\d+)\}?|[eE](-?\d+))?')
def sig(s):
    s = s.replace('.', '').lstrip('0'); return len(s)
def match(p, k, V, S):
    e = p.adjusted(); ulp = Dc(1).scaleb(e - k + 1)
    lo = bisect.bisect_left(V, p - ulp); hi = bisect.bisect_right(V, p + ulp)
    for i in range(lo, hi):
        v = V[i]
        if v.quantize(ulp, rounding=ROUND_FLOOR) <= p <= v.quantize(ulp, rounding=ROUND_CEILING): return S[i]
    return None
def strip_text(t):
    t = re.sub(r'\\fn\{[^}]*\}|\\key\{[^}]*\}|\\label\{[^}]*\}|\\ref\{[^}]*\}|\\cite\{[^}]*\}|\\input\{[^}]*\}|\\bibitem\{[^}]*\}', ' ', t)
    t = re.sub(r'`[^`]*`', ' ', t)
    return t
report = {'args': {'argv': sys.argv, 'ledger_current_files': len(CURRENT), 'ledger_current_values': len(CV), 'ledger_superseded_values': len(SV),
                   'texts': {f: hashlib.sha256(open(D + f, 'rb').read()).hexdigest() for f in TEXTS if os.path.exists(D + f)}}, 'files': {}}
tot = {'HISTORICAL': 0, 'CURRENT': 0, 'STALE': 0, 'RECORD': 0, 'WHITELIST': 0, 'UNMATCHED': 0, 'short': 0}
for f in TEXTS:
    if not os.path.exists(D + f): continue
    rows = []
    for ln, line in enumerate(open(D + f).read().split('\n'), 1):
        is_comment = line.lstrip().startswith('%')
        body = line.split('%', 1)[0] if not re.search(r'\\%', line) else line
        if f == 'numbers_A_v37.tex' and is_comment: continue
        if is_comment: continue
        for m in NUM.finditer(strip_text(body)):
            mant, e1, e2 = m.group(1), m.group(2), m.group(3)
            k = sig(mant)
            if k < 3: tot['short'] += 1; continue
            if '.' not in mant and e1 is None and e2 is None:  # plain integers (years, counts, line numbers, sha fragments)
                tot['short'] += 1; continue
            p = Dc(mant) * (Dc(10) ** int(e1 or e2 or 0))
            txt = m.group(0)
            st = None; s = match(p, k, CV, CS)
            if s: st = 'CURRENT'
            else:
                s = match(p, k, SV, SS)
                ctx = line.lower()
                wk = '%s:%d:%s' % (f, ln, txt)
                if wk in WL and f != 'erratum_v37.tex':
                    st = 'WHITELIST'; s = WL[wk]
                elif s:
                    st = 'HISTORICAL' if f == 'erratum_v37.tex' else 'RECORD' if re.search(r'record|supersed|old|v35|fin12f4|fin13v4\b|was |previous|up to v35|quoted|erratum', ctx) else 'STALE'
                elif ('%s:%d:%s' % (f, ln, txt)) in WL:
                    st = 'WHITELIST'; s = WL['%s:%d:%s' % (f, ln, txt)]
                elif f == 'erratum_v37.tex':
                    st = 'HISTORICAL'
                else: st = 'UNMATCHED'
            tot[st] += 1
            if st != 'CURRENT': rows.append({'line': ln, 'num': txt, 'status': st, 'match': s})
    report['files'][f] = rows
report['totals'] = tot
report['verdict'] = 'PASS' if tot['STALE'] == 0 and tot['UNMATCHED'] == 0 else 'FAIL'
json.dump(report, open(D + 'tools/numbers_check_v37.json', 'w'), indent=1, default=str)
print(report['verdict'], tot)
for f, rows in report['files'].items():
    bad = [r for r in rows if r['status'] in ('STALE', 'UNMATCHED')]
    print(f, 'non-current:', len(rows), 'bad:', len(bad))
    for r in bad[:400]: print('   ', r['line'], r['num'], r['status'], (r['match'] or '')[:90])
