#!/usr/bin/env python3
"""Stage C4 (fin13v4): Hardy-type weighted integrals on Z1 = [16.25, 57] (EI_SPEC 4.5, 4.7, 4.9 T12).
On Z1, g = v, so |w(x)| <= |w(16.25)| + Phi(x)^{1/2} V0 with Phi(x) = int_16.25^x 1/W  (Cauchy-Schwarz).  Hence for a box-sup function f(x):
   ||1_Z1 f w||_W <= h_f |w(16.25)| + H_f V0,   h_f = (int_Z1 W f^2)^{1/2},  H_f = (int_Z1 W f^2 Phi)^{1/2}.
f in {A0, B0, EA1} per alpha-level (per-box sups from stage C2, out/c2_t3d_boxes.json, zones ZA and ZB, boxes with x_lo < 57),
and C0Hb = (int_34.25^57 W c0'^2 (Phi - Phi(34.25)))^{1/2} (c0' = -(5/2) m p + k^2 q, profile bounds per box).
W = exp(-psi) of the P=2 / R=60 design (agent_int/v2/int_common.W_up, psi increasing => W(x) <= W(x_lo) on a box);
Phi(x) <= Phi_up(x) := sum over the cover of [16.25, x] of len * 1/W_lo(x_hi) (1/W increasing).   Also Phi(34.25) vs the typed PHI34.
Output out/c4_hardy.json.  Run: PYTHONHASHSEED=0 nice -n 19 python3 c4_hardy.py"""
import sys, os, json, time, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from v4common import args_block, dump, OUT, EXT, V4, fup
V2 = os.path.join(V4, 'agent_int', 'v2')
sys.path.insert(0, V2)
import int_common as IC
from flint import arb
T0 = time.time()
BX = os.path.join(OUT, 'c2_t3d_boxes.json')
B = json.load(open(BX))['boxes']
Z1HI = 57.0
# Phi table on a fine cover (h = 1/64) of [16.25, 60]
H = 1.0 / 64; NPH = int(round((60.0 - 16.25) / H))
xs = [16.25 + i * H for i in range(NPH + 1)]
Phi = [arb(0)]
for i in range(NPH):
    Phi.append(Phi[-1] + arb(H) / IC.W_lo('2', xs[i + 1]))
def Phi_up(x):
    """upper bound of Phi(x): first table node >= x"""
    i = min(NPH, max(0, math.ceil((x - 16.25) / H - 1e-12)))
    return Phi[i]
PHI34 = Phi_up(34.25)
def hardy(lev, key, lo=16.25, hi=Z1HI, shift=None):
    h2 = arb(0); H2_ = arb(0)
    for zn in ('ZA', 'ZB'):
        for b in B[f"{zn}_{lev}"]:
            if b['xb'] <= lo or b['xa'] >= hi: continue
            xa = max(b['xa'], lo); xb = min(b['xb'], hi)
            if xb <= xa: continue
            L = arb(xb) - arb(xa); Wv = IC.W_up('2', max(xa, 1.0)); f = arb(b[key])
            ph = Phi_up(xb) if shift is None else Phi_up(xb) - shift
            h2 += L * Wv * f * f; H2_ += L * Wv * f * f * ph
    return fup(h2.sqrt()), fup(H2_.sqrt())
out = dict(description=__doc__, Phi_34_25_upper=fup(PHI34), PHI34_typed=9.8547e6, PHI34_typed_covers=bool(PHI34 <= arb('9.8547e6')))
for lev in ('0.00025', '0.001'):
    o = {}
    for key in ('A0', 'B0', 'EA1'):
        h_, H_ = hardy(lev, key); o['h' + key] = h_; o['H' + key] = H_
    out[lev] = o
# C0Hb: int_34.25^57 W c0'^2 (Phi - Phi(34.25)), profile bounds per sigma-box (hstail, monotone in sigma)
sys.path.insert(0, '/home/claude/s1spec/step3')
import hstail as HS
TIPLO = -2.2490776515195
acc = arb(0); n = 23 * 8; hh = (Z1HI - 34.25) / n
Phi34lo = arb(0)    # Phi - Phi(34.25) <= Phi_up(x) - 0 is too crude near 34.25 only by Phi(34.25); use the difference of table sums
for i in range(n):
    xa, xb = 34.25 + i * hh, 34.25 + (i + 1) * hh
    sig = xa + TIPLO - 1e-9
    kb, mb, pb, qs = HS.profile_bounds(arb(sig))
    qmax = arb(xb) - arb(TIPLO) + arb(qs.upper())
    c0p = arb(5) / 2 * arb(mb.upper()) * arb(pb.upper()) + arb(kb.upper()) ** 2 * qmax
    dphi = Phi_up(xb) - Phi[max(0, math.floor((34.25 - 16.25) / H + 1e-12))]
    acc += arb(xb - xa) * IC.W_up('2', xa) * c0p * c0p * dphi
out['C0Hb'] = fup(acc.sqrt())
out['args'] = args_block([BX, os.path.join(V2, 'int_common.py')] + EXT, dict(Phi_step=H))
out['runtime_s'] = round(time.time() - T0, 1)
h = dump(os.path.join(OUT, 'c4_hardy.json'), out)
print(json.dumps({k_: v for k_, v in out.items() if k_ not in ('args', 'description')}, indent=1)); print('sha256', h)
