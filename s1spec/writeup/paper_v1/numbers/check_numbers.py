#!/usr/bin/env python3
"""Number-ledger checker for paper_v1 (unit U0).  Reads only; writes nothing.
  1. sha256 of every source file against ledger.json;
  2. every key resolves; the exact bound is recomputed (Fractions, outward) from the raw JSON values;
  3. the printed string equals the exact outward rounding of that bound at the printed digit
     (upper: up; lower and margin: down; exact: equal) and numbers_paper.tex equals the ledger;
  4. 'doc' entries: the quoted snippet is still present in its document;
  5. given .tex files as arguments: undefined macros, unused macros, and literal numbers that are in the ledger
     (they must be macros) or look like unledgered certified numbers (warning).
usage: nice -n 19 python3 -I check_numbers.py [file.tex ...]     (exit 1 on any FAIL)
Convention: a JSON decimal literal is taken as the value it spells (the pipelines store doubles already rounded outward)."""
import json, hashlib, os, re, sys, ast
from fractions import Fraction as Fr
from decimal import Decimal, getcontext
getcontext().prec = 80
HERE = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------- exact helpers
def dec(x):
    return Fr(Decimal(repr(x))) if isinstance(x, float) else Fr(x)

def resolve(d, path):
    """path components separated by '.', but dict keys may contain '.' or ' ' (greedy longest match)."""
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

def interval(v):
    if isinstance(v, dict) and 'arb' in v:
        a = v['arb']
        if isinstance(a, list):
            m = Fr(int(a[0])) * (Fr(2) ** a[1]); r = Fr(a[2]); return m - r, m + r
        return interval(a)
    if isinstance(v, bool): raise TypeError('bool')
    if isinstance(v, (int, float)): x = dec(v); return x, x
    if isinstance(v, str):
        s = v.strip()
        m = re.match(r'^\[?\s*([-+0-9.eE]+)\s*\+/-\s*([-+0-9.eE]+)\s*\]?$', s)
        if m: c = Fr(Decimal(m.group(1))); r = Fr(Decimal(m.group(2))); return c - r, c + r
        if '/' in s and re.match(r'^-?\d+/\d+$', s): x = Fr(s); return x, x
        x = Fr(Decimal(s.strip('[]'))); return x, x
    raise TypeError(repr(v)[:80])

def hull(ivs):
    return min(a for a, b in ivs), max(b for a, b in ivs)

def wild(d, pattern):
    """expand a path with '*' components (dict values or list items); returns the list of values."""
    parts = pattern.split('.'); cur = [d]
    for p in parts:
        nxt = []
        for c in cur:
            if p == '*':
                nxt += list(c.values()) if isinstance(c, dict) else list(c)
            else:
                try: nxt.append(resolve(c, p))
                except (KeyError, TypeError, IndexError): pass
        cur = nxt
    return cur

def get_input(D, spec):
    """spec = {'source','key'}; key may be 'len:PATH', 'max:PATTERN', 'min:PATTERN' or plain path.  Returns (lo, hi, raw)."""
    d = D[spec['source']]; k = spec['key']
    if k.startswith('len:'):
        v = resolve(d, k[4:]); n = Fr(len(v)); return n, n, len(v)
    if k.startswith('max:') or k.startswith('min:'):
        vals = wild(d, k[4:]); ivs = [interval(v) for v in vals]
        f = max if k.startswith('max:') else min
        lo, hi = f(a for a, b in ivs), f(b for a, b in ivs)
        return lo, hi, ('%s of %d values' % (k[:3], len(vals)))
    v = resolve(d, k); lo, hi = interval(v); return lo, hi, v

# interval arithmetic for expressions: + - * / exp ; operands are (lo, hi) Fractions
def _exp_end(x, up):
    e = Decimal(x.numerator) / Decimal(x.denominator); r = e.exp()
    f = Fr(r); eps = abs(f) * Fr(1, 10 ** 50)
    return f + eps if up else f - eps
def ev(node, env):
    if isinstance(node, ast.Expression): return ev(node.body, env)
    if isinstance(node, ast.Constant): x = dec(node.value); return x, x
    if isinstance(node, ast.Name): return env[node.id]
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        a, b = ev(node.operand, env); return -b, -a
    if isinstance(node, ast.BinOp):
        (a, b), (c, d) = ev(node.left, env), ev(node.right, env)
        if isinstance(node.op, ast.Add): return a + c, b + d
        if isinstance(node.op, ast.Sub): return a - d, b - c
        if isinstance(node.op, ast.Mult):
            p = [a * c, a * d, b * c, b * d]; return min(p), max(p)
        if isinstance(node.op, ast.Div):
            if c <= 0 <= d: raise ZeroDivisionError('interval divisor contains 0')
            p = [a / c, a / d, b / c, b / d]; return min(p), max(p)
        if isinstance(node.op, ast.Pow) and isinstance(node.right, ast.Constant) and isinstance(node.right.value, int):
            n = node.right.value; p = [a ** n, b ** n]; return min(p), max(p)
    if isinstance(node, ast.Call) and node.func.id == 'exp':
        a, b = ev(node.args[0], env); return _exp_end(a, False), _exp_end(b, True)
    raise ValueError('unsupported expression node %s' % ast.dump(node)[:60])

def entry_bound(D, e):
    """exact outward bound (Fraction) of a non-doc entry; also returns the raw list."""
    env = {}; raws = []
    for inp in e['inputs']:
        lo, hi, raw = get_input(D, inp); env[inp['var']] = (lo, hi); raws.append(raw)
    lo, hi = ev(ast.parse(e['expr'], mode='eval'), env)
    return (hi if e['direction'] == 'upper' else lo if e['direction'] in ('lower', 'margin') else None), (lo, hi), raws

# ---------------------------------------------------------------- printing
def round_fmt(x, fmt, mode):
    """x Fraction; fmt = {'mode':'sig'|'dec','n':int,'style':'sci'|'plain'}; mode 'up'|'down'|'exact'.  Returns the printed string."""
    if x == 0: return '0'
    neg = x < 0
    ax = abs(x)
    if fmt['mode'] == 'sig':
        e = 0
        while ax >= Fr(10) ** (e + 1): e += 1
        while ax < Fr(10) ** e: e -= 1
        q = Fr(10) ** (e - fmt['n'] + 1)
    else:
        q = Fr(10) ** (-fmt['n']); e = None
    t = ax / q
    if mode == 'exact':
        if t.denominator != 1: raise ValueError('exact value %s not representable at the printed digit' % x)
        n = int(t)
    else:
        up = (mode == 'up') != neg
        n = -(-t.numerator // t.denominator) if up else t.numerator // t.denominator
    if fmt['mode'] == 'sig':
        if n >= 10 ** fmt['n']: n //= 10; q *= 10; e += 1
        digs = str(n)
        if fmt['style'] == 'sci':
            m = digs[0] + ('.' + digs[1:] if len(digs) > 1 else '')
            s = '%se%d' % (m, e)
        else:
            nd = max(0, fmt['n'] - 1 - e)
            s = _plain(n, q, nd)
    else:
        nd = fmt['n']; s = _plain(n, q, nd)
    return ('-' if neg else '') + s

def _plain(n, q, nd):
    v = Fr(n) * q
    if nd == 0: return str(int(v))
    sc = v * 10 ** nd; assert sc.denominator == 1
    s = str(int(sc)).rjust(nd + 1, '0'); return s[:-nd] + '.' + s[-nd:]

def parse_printed(s):
    """printed string -> Fraction"""
    return Fr(Decimal(s))

def tex_of(s):
    m = re.match(r'^(-?[0-9.]+)e(-?\d+)$', s)
    if m: return r'\ensuremath{%s\cdot10^{%d}}' % (m.group(1), int(m.group(2)))
    return r'\ensuremath{%s}' % s

MANIFESTS = ['/home/claude/s1spec/step4/p2/out/FINAL_manifest.json', '/home/claude/s1spec/step4/p2/out_v33/FINAL_manifest_v33.json']
def manifest_pins():
    pins = {}
    for mf in MANIFESTS:
        m = json.load(open(mf))
        items = []
        for st in m.get('stages', []): items += st.get('inputs', []) + st.get('outputs', [])
        items += m.get('upstream', [])
        for it in items:
            if isinstance(it, dict) and 'path' in it and 'sha256' in it:
                p = os.path.normpath(os.path.join('/home/claude/s1spec/step4', it['path']))
                pins.setdefault(p, set()).add(it['sha256'])
    return pins
def pin_status(path, h, pins):
    if path not in pins: return 'unpinned'
    return 'pinned' if h in pins[path] else 'PIN-MISMATCH'

def load_sources(L):
    D = {}; bad = []
    pins = manifest_pins()
    for tag, f in L['sources'].items():
        if f['path'] is None: continue
        if not os.path.exists(f['path']): bad.append('missing source %s' % f['path']); continue
        h = hashlib.sha256(open(f['path'], 'rb').read()).hexdigest()
        if h != f['sha256']: bad.append('sha256 %s %s: ledger %s file %s' % (tag, f['path'], f['sha256'][:12], h[:12]))
        ps = pin_status(f['path'], h, pins)
        if ps != f.get('pin_status'): bad.append('pin status of %s changed: ledger %s now %s' % (f['path'], f.get('pin_status'), ps))
        if ps == 'PIN-MISMATCH': print('NOTE source %s differs from its manifest pin (documented exception, R7)' % f['path'])
        if f['path'].endswith('.json'): D[tag] = json.load(open(f['path']))
    return D, bad

def macro_text(L):
    out = ['% numbers_paper.tex -- GENERATED from ledger.json by build_ledger.py; do not edit.',
           '% Rounding: upper bounds up, lower bounds and margins down at the printed digit; margin = threshold / bound.',
           '% Use as $\\times\\NName$ for margins.  kind: cert certified, est estimate, hand hypothesis constant, doc only a document quotes it.']
    for e in L['entries']:
        src = ('doc:' if e['kind'] == 'doc' else '') + (os.path.basename(e['sources'][0]['path']) if e.get('sources') else 'hand')
        keys = '; '.join(i['key'] for i in e['inputs']) if e.get('inputs') else e.get('doc_snippet', '')[:60]
        out.append('%% %s [%s, %s] %s :: %s' % (e['name'], e['kind'], e['direction'], src, keys))
        out.append('\\newcommand{\\%s}{%s}' % (e['name'], tex_of(e['printed'])))
    return '\n'.join(out) + '\n'

# ---------------------------------------------------------------- the checks
def main(argv):
    L = json.load(open(os.path.join(HERE, 'ledger.json')))
    fails = []; warns = []
    def fail(m): fails.append(m); print('FAIL', m)
    def warn(m): warns.append(m); print('WARN', m)
    D, bad = load_sources(L)
    for b in bad: fail(b)
    names = set(); nchk = 0; ndoc = 0; nhand = 0
    for e in L['entries']:
        n = e['name']
        if n in names: fail('duplicate name ' + n)
        names.add(n)
        if not re.match(r'^N[A-Za-z]+$', n): fail('name not letters only: ' + n)
        if e['kind'] not in ('cert', 'est', 'hand', 'lit', 'doc'): fail('%s: bad kind %s' % (n, e['kind']))
        fmt = e['fmt']
        if e['kind'] == 'doc':
            ndoc += 1
            p = e['doc_file']; t = open(p, errors='replace').read()
            if hashlib.sha256(open(p, 'rb').read()).hexdigest() != e['doc_sha256']: warn('%s: document changed (sha256)' % n)
            if e['doc_snippet'] not in t: fail('%s: snippet not found in %s' % (n, p))
            if e['doc_value'].replace('-', '\u2212') not in e['doc_snippet'] and e['doc_value'] not in e['doc_snippet']: fail('%s: raw value not inside snippet' % n)
            want = e['doc_value'] if e['fmt'] is None or e['printed'] == e['doc_value'] else round_fmt(Fr(Decimal(e['doc_value'])), e['fmt'], 'up' if e['direction'] == 'upper' else 'down')
            if want != e['printed']: fail('%s: doc printed %s but rounding of %s gives %s' % (n, e['printed'], e['doc_value'], want))
            if not e.get('flag'): fail('%s: doc entry must be flagged' % n)
            continue
        try:
            bnd, (lo, hi), raws = entry_bound(D, e)
        except Exception as ex:
            fail('%s: cannot evaluate (%s: %s)' % (n, type(ex).__name__, ex)); continue
        nchk += 1
        if e['kind'] == 'hand': nhand += 1
        direction = e['direction']
        if direction == 'exact':
            if lo != hi: fail('%s: exact entry has an interval' % n); continue
            want = round_fmt(lo, fmt, 'exact')
        elif direction == 'upper': want = round_fmt(hi, fmt, 'up')
        else: want = round_fmt(lo, fmt, 'down')
        if want != e['printed']: fail('%s: printed %s but outward rounding of the bound gives %s' % (n, e['printed'], want))
        Q = parse_printed(e['printed'])
        if direction == 'upper' and not Q >= hi: fail('%s: printed %s below upper end %.12g' % (n, e['printed'], float(hi)))
        if direction in ('lower', 'margin') and not Q <= lo: fail('%s: printed %s above lower end %.12g' % (n, e['printed'], float(lo)))
        if e.get('raw_ledger') is not None and [json.dumps(r, sort_keys=True) for r in raws] != e['raw_ledger']:
            fail('%s: raw value changed' % n)
    gen = macro_text(L)
    cur = open(os.path.join(HERE, 'numbers_paper.tex')).read()
    if gen != cur: fail('numbers_paper.tex differs from the ledger (rerun build_ledger.py)')
    print('entries: %d (evaluated %d, of which hand %d; doc %d); sources: %d' % (len(L['entries']), nchk, nhand, ndoc, len(L['sources'])))
    # ---- tex files
    texs = [a for a in argv if a.endswith('.tex')]
    if texs:
        defined = set(re.findall(r'\\newcommand\{\\(N[A-Za-z]+)\}', cur))
        used = {}
        printed = {}
        for e in L['entries']:
            if e['kind'] == 'doc' and False: continue
            printed.setdefault(norm(e['printed']), []).append(e['name'])
        for f in texs:
            raw = open(f, errors='replace').read()
            body = re.sub(r'(?<!\\)%.*', '', raw)
            for m in re.finditer(r'\\(N[A-Za-z]+)(?![A-Za-z])', body):
                used.setdefault(m.group(1), []).append(os.path.basename(f))
            for m in re.finditer(r'\\(?:cite[a-z]*|ref|eqref|label|input|include|fn|url|href)\{[^}]*\}', body): pass
            clean = re.sub(r'\\(?:cite[a-z]*|ref|eqref|label|input|include|fn|url|href|cref|Cref)\{[^}]*\}', ' ', body)
            clean = re.sub(r'\\newcommand\{\\N[A-Za-z]+\}', ' ', clean)
            for ln, line in enumerate(clean.split('\n'), 1):
                for m in LIT.finditer(line):
                    s = lit_to_printed(m)
                    if s is None: continue
                    k = norm(s)
                    if re.match(r'^(19|20)\d\d$', s): continue
                    if k in printed:
                        fail('%s:%d literal number %s is a ledger value, use \\%s' % (os.path.basename(f), ln, s, printed[k][0]))
                    elif sigdigits(s) >= 4:
                        warn('%s:%d literal %s with >= 4 significant digits is not in the ledger' % (os.path.basename(f), ln, m.group(0).strip()))
        for u in sorted(used):
            if u not in defined: fail('undefined macro \\%s (used in %s)' % (u, ','.join(sorted(set(used[u])))))
        unused = sorted(defined - set(used))
        print('macros used: %d of %d defined; unused: %d' % (len(set(used) & defined), len(defined), len(unused)))
        if len(unused) <= 25:
            for u in unused: print('INFO unused macro \\%s' % u)
    print('SUMMARY: %d FAIL, %d WARN' % (len(fails), len(warns)))
    print('PASS' if not fails else 'FAIL')
    return 1 if fails else 0

LIT = re.compile(r'(?<![\w.\\{])(\d{1,3}(?:\{,\}\d{3})+|\d+\.\d+|\d{4,})(?:\s*\\(?:cdot|times)\s*10\^\{?(-?\d+)\}?|[eE](-?\d+))?(?![\w.])')
def lit_to_printed(m):
    mant = m.group(1).replace('{,}', '')
    ex = m.group(2) or m.group(3)
    try: return mant + ('e' + ex if ex else '')
    except Exception: return None
def norm(s):
    """canonical (digits, exponent) form of a printed number: leading/trailing zeros ignored"""
    try: d = Decimal(s)
    except Exception: return s
    if d == 0: return '0'
    t = d.normalize(); sign, digs, ex = t.as_tuple()
    return '%s%s@%d' % ('-' if sign else '', ''.join(map(str, digs)), ex)
def sigdigits(s):
    d = Decimal(s).normalize(); return len(d.as_tuple().digits)

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
