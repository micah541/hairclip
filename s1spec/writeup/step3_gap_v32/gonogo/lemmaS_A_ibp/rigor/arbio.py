"""exact (mid, rad) serialisation of arb balls for rigor/*.json"""
import math
from flint import arb, ctx
def dump(x):
    if isinstance(x, arb):
        m, e = x.mid().man_exp()
        return {'arb': [str(m), int(e), math.nextafter(float(x.rad()), math.inf)]}
    if isinstance(x, dict): return {str(k): dump(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)): return [dump(v) for v in x]
    return x
def load(x):
    ctx.prec = max(ctx.prec, 300)
    if isinstance(x, dict) and 'arb' in x:
        m, e, r = x['arb']
        return arb(int(m)) * arb(2)**int(e) + arb(0, r)
    if isinstance(x, dict): return {k: load(v) for k, v in x.items()}
    if isinstance(x, list): return [load(v) for v in x]
    return x
def up(x):   # certified upper endpoint as an exact arb number
    return x.upper()
def lo(x):
    return x.lower()
def fl_up(x):   # python float >= x
    f = float(x.upper()); return f if arb(f) >= x.upper() else math.nextafter(f, math.inf)
def fl_lo(x):
    f = float(x.lower()); return f if arb(f) <= x.lower() else math.nextafter(f, -math.inf)
