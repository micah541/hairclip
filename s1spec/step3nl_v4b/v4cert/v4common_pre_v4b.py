"""Common helpers for the v4 certification stages (tag fin13v4).  Args blocks, sha256 pins, outward float conversion."""
import os, sys, json, hashlib, math, platform
V4 = '/home/claude/s1spec/step3nl_v4'
CERT = os.path.join(V4, 'v4cert')
OUT = os.environ.get('V4OUT', os.path.join(CERT, 'out'))     # second run of the chain: V4OUT=.../out_B
os.makedirs(OUT, exist_ok=True)
ENVKEYS = ('PYTHONHASHSEED', 'OMP_NUM_THREADS', 'TAG', 'NEW', 'THETA', 'ADD', 'FIX', 'S0V', 'S1V', 'DELTA', 'MOD', 'MODA', 'ETA', 'FLUX0', 'FLUX2', 'LY1SQ')
SCRIPT = os.path.abspath(sys.argv[0])   # taken at import (before any chdir)
def sha(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()
def args_block(inputs, extra=None):
    """argv, env (selected keys), python/flint/sympy versions, sha256 of every input path (and of the running script)"""
    import flint, sympy
    ins = list(inputs) + [SCRIPT]
    d = dict(argv=sys.argv, env={k: os.environ.get(k) for k in ENVKEYS}, python=platform.python_version(), flint=flint.__version__, sympy=sympy.__version__,
             inputs={p: sha(p) for p in sorted(set(ins))})
    if extra: d.update(extra)
    return d
def fup(z):
    """arb ball -> float >= every point of the ball (outward)"""
    assert z.is_finite(), 'non-finite ball'
    return math.nextafter(float(z.upper()), math.inf)
def fdn(z):
    assert z.is_finite(), 'non-finite ball'
    return math.nextafter(float(z.lower()), -math.inf)
def dump(path, obj):
    with open(path, 'w') as f: json.dump(obj, f, indent=1, sort_keys=True)
    return sha(path)
# external read-only modules/data imported by the jet / profile code (pinned by sha256 in every args block that uses them)
EXT = ['/home/claude/s1spec/step3nl/psi3fun.py', '/home/claude/s1spec/step3nl/h2fun.py', '/home/claude/s1spec/step3nl/h2res.py',
       '/home/claude/s1spec/step3nl/psi3md2_14.json', '/home/claude/s1spec/step3nl/h2md_14.json', '/home/claude/s1spec/step3nl/q2_derivs.pkl',
       '/home/claude/s1spec/step3/hstail.py', '/home/claude/s1spec/step3/hs_uforms.pkl', '/home/claude/s1spec/step3/gen_hs3.py',
       '/home/claude/s1_tight/bracket2.json']
