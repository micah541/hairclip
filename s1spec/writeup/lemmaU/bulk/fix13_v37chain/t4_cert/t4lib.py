"""Task 4 (Lemma U arb certification): shared helpers.  Import-safe, no side effects.
Interface numbers are exported twice (arb string + outward-rounded double) with a label (cert / est / hyp)."""
import os, sys, json, math, hashlib, platform
sys.dont_write_bytecode = True
from flint import arb, ctx
import sympy as sp
Q_ = sp.Rational
HERE = os.path.dirname(os.path.abspath(__file__))
OUTD = os.path.join(HERE, 'out')
ROOT = '/home/claude/s1spec/writeup/lemmaU'
P1 = ROOT + '/phase1/ghat'; P2 = ROOT + '/phase1/ramp'
T5 = ROOT + '/bulk/fix13_v37chain/t5_zone3'; T2 = ROOT + '/bulk/fix13_v37chain/t2_A'   # fix6_v36 (v36) copy
sys.path.insert(0, '/home/claude/s1spec/writeup/lemmaU/bulk/fix13_v37chain'); import v36cfg as CFG   # fix6_v36 (v36): verifier source + switches V36_VER/V36_R3/V36_R2
VERIF_ALPHA = CFG.MODS   # v36: per V36_VER (g13 was fin12f4)


def A(q):
    if isinstance(q, arb): return q
    q = Q_(q); return arb(int(q.p))/arb(int(q.q))


def D(s):
    """decimal string -> arb ball containing the decimal exactly (outward)"""
    q = Q_(str(s)); return arb(int(q.p))/arb(int(q.q))


def fup(x):
    v = x.upper() if isinstance(x, arb) else arb(x); f = float(v)
    if f != f: return math.inf
    return f if arb(f) >= v else math.nextafter(f, math.inf)


def flo(x):
    v = x.lower() if isinstance(x, arb) else arb(x); f = float(v)
    if f != f: return -math.inf
    return f if arb(f) <= v else math.nextafter(f, -math.inf)


def num(x, label='cert', d='up', **kw):
    """export an arb as {arb, double (outward), label}"""
    x = x if isinstance(x, arb) else arb(x)
    o = dict(arb=x.str(20, radius=True), double=(fup(x) if d == 'up' else flo(x)), label=label)
    o.update(kw); return o


def ub(x): return arb(fup(x))          # replace an arb by an exactly representable upper bound (stops radius growth)
def amax(xs):
    m = None
    for x in xs:
        u = x.upper() if isinstance(x, arb) else arb(x)
        m = u if m is None or u > m else m
    return m


def sha(p): return hashlib.sha256(open(p, 'rb').read()).hexdigest()


def args_block(inputs, extra=None):
    d = dict(argv=sys.argv, cwd=os.getcwd(), python=sys.version.split()[0], platform=platform.platform(),
             flint=__import__('flint').__version__, sympy=sp.__version__, prec=ctx.prec,
             env={k: os.environ.get(k) for k in ('PYTHONHASHSEED', 'PYTHONDONTWRITEBYTECODE', 'OMP_NUM_THREADS') + CFG.ENVKEYS},
             inputs={p: sha(p) for p in sorted(set(inputs)) if p and os.path.exists(p)})
    if extra: d.update(extra)
    return d


def dump(obj, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    json.dump(obj, open(path + '.tmp', 'w'), indent=1, sort_keys=True, default=str); os.replace(path + '.tmp', path)


def verifier_sups():
    """certified sup|a| (= sup|b|) and sup|beta| from the verifier output (alpha_1 = 1.5e-4 row), as decimal strings -> arb"""
    d = json.load(open(VERIF_ALPHA))['1.5e-4']
    return dict(abar=D(d['sup_ab']), bbar=D(d['sup_ab']), betabar=D(d['sup_beta']), src=VERIF_ALPHA,
                strings=dict(sup_ab=str(d['sup_ab']), sup_beta=str(d['sup_beta'])))
