"""fix8_flux2_impl common helpers: hashing, args block, atomic JSON dump, outward rounding.  No proof-tree writes.
All scripts of this directory: run as  PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1 nice -n 19 python3 -B <script>."""
import os, sys, json, hashlib, platform, math, time
from fractions import Fraction
from flint import arb

ROOT = '/home/claude/s1spec/writeup/lemmaU/bulk/fix13_v37chain/flux'   # v37: copy of fix8_flux2_impl/code (flux lemma v2 reassembly on fin13v4c + v37)
OUT = os.environ.get('F2_OUT', os.path.join(ROOT, 'out'))     # F2_OUT: dry runs only (scratch); production uses ROOT/out
CODE = os.path.join(ROOT, 'code')
B = '/home/claude/s1spec/writeup/lemmaU/bulk'
VTAG = os.environ.get('F2_VTAG', 'fin13v4c')     # v37: verifier tag; fin13v4b + fix6_v36 = regression against fix8_flux2_impl
assert VTAG in ('fin13v4b', 'fin13v4c'), VTAG
V4B = '/home/claude/s1spec/step3nl_v4c' if VTAG == 'fin13v4c' else '/home/claude/s1spec/step3nl_v4b'   # verifier directory (name kept)
N36 = os.path.join(B, 'fix13_v37chain' if VTAG == 'fin13v4c' else 'fix6_v36')                           # (A) chain directory (name kept)
CHAIN_REPORT = os.path.join(N36, 'V37_LOG.md' if VTAG == 'fin13v4c' else 'V36_REPORT.md')               # report that must name the (A) chain run
F1 = os.path.join(B, 'fix3', 'flux')
DESIGN = os.path.join(B, 'fix7_flux2')


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for blk in iter(lambda: f.read(1 << 20), b''):
            h.update(blk)
    return h.hexdigest()


def args_block(inputs, extra=None):
    """argv, env, versions, every input path with sha256 (timing goes in a separate 'timing' key)."""
    import sympy
    import flint
    d = dict(argv=sys.argv, cwd=os.getcwd(), python=sys.version.split()[0], platform=platform.platform(),
             flint=flint.__version__, sympy=sympy.__version__,
             env={k: os.environ.get(k) for k in ('PYTHONHASHSEED', 'PYTHONDONTWRITEBYTECODE', 'F2_NCELL') + (('F2_VTAG',) if VTAG != 'fin13v4b' else ())},
             inputs={os.path.abspath(p): sha(p) for p in sorted(set(inputs))})
    if extra:
        d.update(extra)
    return d


def dump(path, obj):
    """atomic write, sorted keys; returns sha256 of the written file"""
    if isinstance(obj, dict) and 'timing' in obj:          # timing goes to a sidecar, so the JSON itself is reproducible byte for byte
        obj = dict(obj); tm = obj.pop('timing')
        with open(path[:-5] + '.timing.json' if path.endswith('.json') else path + '.timing', 'w') as f:
            json.dump(tm, f, sort_keys=True)
    tmp = path + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(obj, f, indent=1, sort_keys=True)
        f.write('\n')
    os.replace(tmp, path)
    return sha(path)


# ---------------------------------------------------------------- outward rounding to doubles
def fup(a):
    """upper bound of an arb as a double, rounded outward by one ulp (0 stays 0 only if the ball is exactly <= 0)"""
    a = arb(a)
    if not a.is_finite():
        raise ValueError('non-finite ball %s' % a)
    v = float(a.upper())
    return math.nextafter(v, math.inf)


def fdn(a):
    a = arb(a)
    if not a.is_finite():
        raise ValueError('non-finite ball %s' % a)
    v = float(a.lower())
    return math.nextafter(v, -math.inf)


def aup(a):
    """arb point = outward-rounded double upper bound"""
    return arb(fup(a))


def adn(a):
    return arb(fdn(a))


def from_raw_up(v):
    """a raw JSON double that is an UPPER bound, rounded up by one more ulp (convention: constants from raw JSON, rounded up)"""
    return arb(math.nextafter(float(v), math.inf))


def from_raw_dn(v):
    return arb(math.nextafter(float(v), -math.inf))


def load_arb(j):
    """the {'arb': [mant, exp, rad]} or {'arb': '[mid +/- rad]'} encodings used in the (A) chain outputs"""
    a = j['arb']
    if isinstance(a, str):
        return arb(a)
    m, e, r = a
    x = arb(int(m)) * arb(2) ** int(e)
    return x + arb(0, r) if r else x


def amax(*xs):
    b = xs[0]
    for x in xs[1:]:
        if x.upper() > b.upper():
            b = x
    return b


def astr(a):
    return a.str(20, radius=True)


def rec(a, label, note=None):
    """a JSON record of an upper bound: outward double + arb string + label"""
    d = dict(up=fup(a), arb=astr(arb(a)), label=label)
    if note:
        d['note'] = note
    return d


def rec_dn(a, label, note=None):
    d = dict(dn=fdn(a), arb=astr(arb(a)), label=label)
    if note:
        d['note'] = note
    return d


def frac_arb(fr):
    fr = Fraction(fr)
    return arb(fr.numerator) / arb(fr.denominator)


class Timer:
    def __init__(self):
        self.t0 = time.time()

    def s(self):
        return round(time.time() - self.t0, 1)


def no_pycache_check(dirs):
    bad = []
    for d in dirs:
        for r, ds, fs in os.walk(d):
            if '__pycache__' in ds:
                bad.append(os.path.join(r, '__pycache__'))
    return bad
