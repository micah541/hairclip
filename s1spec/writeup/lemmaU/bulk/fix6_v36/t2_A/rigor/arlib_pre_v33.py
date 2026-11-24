"""Shared helpers for the Phase-2a repair of Lemma S(A) (new files only; v3.1 untouched).
- imports the exact-arithmetic machinery of the certified (A) chain from the COPY  A_repair/base_ibp/rigor/weights_exact.py
  (byte-identical copy of step3_gap_v32/gonogo/lemmaS_A_ibp/rigor/weights_exact.py; that module chdirs on import, we restore cwd);
- args block (argv, env, input paths with sha256) for every output JSON;
- honest jet-free coefficient bounds Ab_h = Z-free part of the symbolic Aterms (includes every profile jet U_j = p^j U, repair 3).
All numbers produced with these helpers are arb (python-flint, outward rounding)."""
import os, sys, json, hashlib, math, platform
sys.dont_write_bytecode = True
ROOT = '/home/claude/s1spec/writeup/lemmaU/bulk/t2_A'
BASE = ROOT + '/base_ibp'                      # copy of gonogo/lemmaS_A_ibp (reference chain, tightened (H2))
BASE31 = ROOT + '/base_v31H2'                  # copy of gonogo/lemmaS_A_ibp_v31H2 (reference chain, v3.1 (H2))
OUT = ROOT + '/out'
os.makedirs(OUT, exist_ok=True)
_cwd = os.getcwd()
sys.path.insert(0, BASE + '/rigor'); sys.path.insert(0, BASE)
import weights_exact as WX                      # noqa: E402  (chdirs to BASE)
import arbio                                    # noqa: E402
os.chdir(_cwd)
import sympy as sp                              # noqa: E402
from flint import arb, ctx                      # noqa: E402
ctx.prec = 200
Q_ = sp.Rational
A = WX.A
P_Q = Q_(10246, 10000)

def sha(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()

def args_block(inputs, extra=None):
    d = dict(argv=list(sys.argv), cwd=os.getcwd(), python=sys.version, platform=platform.platform(),
             env={k: os.environ.get(k) for k in ('PYTHONHASHSEED', 'PYTHONDONTWRITEBYTECODE', 'OMP_NUM_THREADS')},
             inputs={p: sha(p) for p in inputs})
    import flint
    d['flint'] = getattr(flint, '__version__', '?'); d['sympy'] = sp.__version__
    if extra: d.update(extra)
    return d

def dump(obj, path, inputs, extra=None):
    """atomic write with args block"""
    obj = dict(obj); obj['args'] = args_block(inputs, extra)
    tmp = path + '.tmp'
    json.dump(arbio.dump(obj), open(tmp, 'w'), indent=1, default=str)
    os.replace(tmp, path)

def load(path):
    return arbio.load(json.load(open(path)))

def up(x): return x.upper()
def fup(x):
    f = float(x.upper()); return f if arb(f) >= x.upper() else math.nextafter(f, math.inf)
def flo(x):
    f = float(x.lower()); return f if arb(f) <= x.lower() else math.nextafter(f, -math.inf)

_ASYM = {}
def asym(Uq):
    """symbolic Aterms with profile bound U = Uq (exact rational); cached"""
    if Uq not in _ASYM:
        old = WX.UQ; WX.UQ = Q_(Uq)
        try: _ASYM[Uq] = WX.build_At_sym(rmax=6)
        finally: WX.UQ = old
    return _ASYM[Uq]

def at_eval(M0q, M1q, Uq):
    """At[(i,r)][X] arb upper bounds (X = tuple of jet indices >= 2; () = jet-free incl. U_j)"""
    old = WX.UQ; WX.UQ = Q_(Uq)
    try: return WX.eval_At(asym(Uq), Q_(M0q), Q_(M1q), rmax=6)
    finally: WX.UQ = old

def ab_honest(At):
    """repair 3: jet-free coefficient bounds INCLUDING the profile jets U_j (j >= 2), = Z-free part of Aterms"""
    return {k: d.get((), arb(0)) for k, d in At.items()}

def ab_v31(M0q, M1q, Uq, rmax=6):
    old = WX.UQ; WX.UQ = Q_(Uq)
    try: return WX.ab_exact(Q_(M0q), Q_(M1q), rmax)
    finally: WX.UQ = old

def U_of_y(y, al=Q_(15, 100000)):
    """|U(y)| <= 1.3056*alpha_1*y^p (arb)"""
    return A(Q_(13056, 10000))*A(al)*(A(y).log()*A(P_Q)).exp() if not isinstance(y, arb) else A(Q_(13056, 10000))*A(al)*(y.log()*A(P_Q)).exp()
