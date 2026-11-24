"""Common helpers for the v4b certification stages (tag fin13v4b): args blocks, sha256 pins, outward float conversion.
Everything is written under /home/claude/s1spec/step3nl_v4c/v4bcert/out (or V4BOUT).  Read-only inputs elsewhere are pinned by sha256."""
import os, sys, json, hashlib, math, platform
sys.dont_write_bytecode = True
V4B = '/home/claude/s1spec/step3nl_v4c'
CERT = os.path.join(V4B, 'v4bcert')
OUT = os.environ.get('V4BOUT', os.path.join(CERT, 'out'))
os.makedirs(OUT, exist_ok=True)
ENVKEYS = ('PYTHONHASHSEED', 'PYTHONDONTWRITEBYTECODE', 'V4BOUT')
SCRIPT = os.path.abspath(sys.argv[0])
def sha(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()
def args_block(inputs, extra=None):
    import flint, sympy
    ins = list(inputs) + [SCRIPT, os.path.abspath(__file__)]
    d = dict(argv=sys.argv, env={k: os.environ.get(k) for k in ENVKEYS}, python=platform.python_version(), flint=flint.__version__, sympy=sympy.__version__,
             inputs={p: sha(p) for p in sorted(set(ins))})
    if extra: d.update(extra)
    return d
def fup(z):
    assert z.is_finite(), 'non-finite ball'
    return math.nextafter(float(z.upper()), math.inf)
def fdn(z):
    assert z.is_finite(), 'non-finite ball'
    return math.nextafter(float(z.lower()), -math.inf)
def dump(path, obj):
    with open(path, 'w') as f: json.dump(obj, f, indent=1, sort_keys=True)
    return sha(path)
# read-only modules / data used by the jet and profile code (same list as v4cert/v4common.EXT)
EXT = ['/home/claude/s1spec/step3nl/psi3fun.py', '/home/claude/s1spec/step3nl/h2fun.py', '/home/claude/s1spec/step3nl/h2res.py',
       '/home/claude/s1spec/step3nl/psi3md2_14.json', '/home/claude/s1spec/step3nl/h2md_14.json', '/home/claude/s1spec/step3nl/q2_derivs.pkl',
       '/home/claude/s1spec/step3/hstail.py', '/home/claude/s1spec/step3/hs_uforms.pkl', '/home/claude/s1spec/step3/gen_hs3.py',
       '/home/claude/s1_tight/bracket2.json', '/home/claude/s1_tight/prof/prof_tight.txt']
