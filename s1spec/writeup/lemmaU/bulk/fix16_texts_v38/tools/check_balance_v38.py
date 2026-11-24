# check_balance_v38.py -- brace, \begin/\end, $ and \ref/\label checks for the text-v38 files (LaTeX is not installed).
# The assembled main is checked on its \input-expanded text.  Usage: python3 -B check_balance_v38.py  (writes balance_v38.json)
import re, os, json, sys
D = '/home/claude/s1spec/writeup/lemmaU/bulk/fix16_texts_v38/'
FILES = ['bootstrap_v38.tex', 'ei_lemma_v4c.tex', 'E0_LEMMA_v4c_v38.tex', 'flux2_lemma_v38.tex', 'lemmaS_A_v38.tex', 'numbers_A_v38.tex',
         'continuation_v38.tex', 'lemmaU_v38.tex', 'main_v3_v38.tex', 'erratum_v38.tex', 'step3_v38_main.tex', 'b12_preamble_v38.tex',
         'chart_lemma_body_v38.tex', 'ei_constants_v38.tex', 'C1_identification_v38.tex', 'GS_lemma_v38.tex', 'U1_regularity_v38.tex',
         'G1_export_k2_v38.tex']
def strip(t):
    out = []
    for line in t.split('\n'):
        line = re.sub(r'\\\\', '  ', line); line = re.sub(r'\\[{}%$]', '  ', line)
        i = line.find('%'); out.append(line if i < 0 else line[:i])
    return '\n'.join(out)
def expand(name, seen=()):
    p = D + (name if name.endswith('.tex') else name + '.tex')
    if not os.path.exists(p): return None
    t = open(p).read()
    def sub(m):
        e = expand(m.group(1), seen + (name,))
        return e if e is not None else m.group(0)
    return re.sub(r'^(?!%)\\input\{([^}]*)\}', lambda m: sub(m), t, flags=re.M)
def check(t):
    s = strip(t); r = {}
    d = 0; mn = 0
    for ch in s:
        d += (ch == '{') - (ch == '}'); mn = min(mn, d)
    r['brace_final'] = d; r['brace_min'] = mn
    st = []; bad = []
    for m in re.finditer(r'\\(begin|end)\{([^}]*)\}', s):
        if m.group(1) == 'begin': st.append(m.group(2))
        elif not st or st.pop() != m.group(2): bad.append(m.group(2))
    r['env_unclosed'] = st; r['env_mismatch'] = bad
    r['dollar_even'] = (s.count('$') % 2 == 0)
    labels = set(re.findall(r'\\label\{([^}]*)\}', s)); refs = set(re.findall(r'\\ref\{([^}]*)\}', s))
    r['refs_missing'] = sorted(refs - labels)
    r['ok'] = d == 0 and mn >= 0 and not st and not bad and r['dollar_even']
    return r
res = {}
for f in FILES:
    if not os.path.exists(D + f): res[f] = {'missing': True}; continue
    res[f] = check(open(D + f).read())
t = expand('step3_v38_main.tex'); res['step3_v38_main.tex [expanded]'] = check(t)
res['step3_v38_main.tex [expanded]']['duplicate_labels'] = sorted({l for l in re.findall(r'\\label\{([^}]*)\}', strip(t)) if strip(t).count('\\label{%s}' % l) > 1})
json.dump(res, open(D + 'tools/balance_v38.json', 'w'), indent=1)
for k, v in res.items(): print(('PASS' if v.get('ok') else 'FAIL'), k, {a: b for a, b in v.items() if a != 'ok' and b not in ([], 0, True)})
