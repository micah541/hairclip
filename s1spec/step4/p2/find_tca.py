"""Find the base clock point tc_a in a base slab where the base physical time t~ equals the P1 checkpoint time
(dyadic tc_a, 2^-50 grid).  usage: python3 -I find_tca.py CKPT_JSON SLAB"""
import sys, json
S4 = '/home/claude/s1spec/step4'
sys.path.insert(0, f'{S4}/t6'); sys.path.insert(0, f'{S4}/t1')
from flint import arb, ctx
from base import SlabBase
from rep import K0
ctx.prec = 128
ck, slab = sys.argv[1], int(sys.argv[2])
rec = json.load(open(ck)); v = rec['t']
t1 = arb(arb(int(v[0])) * arb(2) ** int(v[1]), arb(int(v[2])) * arb(2) ** int(v[3]))
SB = SlabBase(slab, poldir=f'{S4}/t1/polished_v2'); ctx.prec = 128; n = SB.n
f = lambda tc: float((SB.Y[n + 3](arb(tc)) - t1).mid())
lo, hi = -1.0, 1.0
flo, fhi = f(lo), f(hi)
print(f'P1 t = {t1.str(18)}; slab {slab}: t~(-1) - t = {flo:.4e}, t~(1) - t = {fhi:.4e}')
assert flo < 0 < fhi, 'P1 time not inside this base slab'
for _ in range(60):
    mid = 0.5 * (lo + hi)
    if f(mid) < 0: lo = mid
    else: hi = mid
tca = round(lo * 2 ** 50) / 2 ** 50
tb = SB.Y[n + 3](arb(tca)); lnlam = SB.Y[n + 1](arb(tca))
print(f'tc_a = {tca!r} (= {round(lo * 2 ** 50)}/2^50); t~(tc_a) - t = {(tb - t1).str(6)}; k~ = {(arb(K0) / lnlam.exp()).str(15)}; '
      f'P1 slab {rec["slab"]}')
