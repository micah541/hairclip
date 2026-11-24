import sys, os, time
from concurrent.futures import ThreadPoolExecutor
from flint import arb
import gen_hs3 as G3, gen_hs2 as G2
from gram_post6 import DSet, ratio, na
from gen_evans import profile_part as _pp
def build(variant, path):
    G3.ST_RAD_SAVE = G3.ST_RAD
    if "nost" in variant: G3.ST_RAD = 0.0
    if "noprof" in variant:
        orig = G3.profile_part
        tau1, xlo, xhi, Ct, r0lo, r0hi = orig()
        G3.profile_part = lambda: (tau1, xlo, xhi, Ct, [0.0]*4, [0.0]*4)
    if "nocone" in variant:
        res, T = G2.tail_data(); res = dict(res); res['eta'] = arb(0)
        G3.tail_data = lambda: (res, T)
    G3.build_Z(path, "w")
for v in sys.argv[1:]:
    pass
