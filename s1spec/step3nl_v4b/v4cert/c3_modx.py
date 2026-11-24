#!/usr/bin/env python3
"""Stage C3 (fin13v4): W-norms of the modulation source terms of the v-equation for the C^4 slow-manifold functions (EI_SPEC 4.9, U-6).
Copy of agent_int/v2/int_v56.py (V5 part; same boxes, same pointwise formulas, P=2, R=60), with the jets of psi3~, h2~ taken from the C^4
functions of agent_int/planB/c4fun.py instead of the C^3 ones (h2res.jets_any).  Norms ||zeta f||_W^2 = sum len * W_up(x_lo) * sup zeta^2 * sup f^2.
Quantities: NP12 = (||zeta psi1'||, ||zeta psi2'||), NP3 = ||zeta psi3'||, NH2 = ||zeta h2'||, NAB = (||zeta A1(u0)'||/alpha, ||zeta B1(u0)'||/alpha)
(alpha ball [-1e-3, 1e-3], valid at every alpha-level <= 1e-3), C0P = ||zeta c0'||.   Also reports the C^3 (v2) values for comparison.
Output out/c3_modx.json.  Run: PYTHONHASHSEED=0 nice -n 19 python3 c3_modx.py"""
import sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from v4common import args_block, dump, OUT, EXT, V4
V2 = os.path.join(V4, 'agent_int', 'v2'); PB = os.path.join(V4, 'agent_int', 'planB')
sys.path.insert(0, PB); sys.path.insert(0, V2)
import int_common as IC
from c4fun import PSI, H2
from flint import arb
T0 = time.time()
AMAX = 1e-3
A = arb(0, AMAX)
R = 60.0
def Spoly_f(t):
    t = min(max(t, 0.0), 1.0); a = arb(t)
    return a * a * a * a * (35 - 84 * a + 70 * a * a - 20 * a * a * a)
def zeta_sup(xa, xb):
    zi = IC.up(Spoly_f((xb - 10.25) / 6)) if xb < 16.25 else 1.0
    zo = IC.up(Spoly_f((R - xa) / 3)) if xa > R - 3 else 1.0
    return min(1.0, zi) * min(1.0, zo)
def quantities(b, y, h):
    k, m, p, q = b.k, b.m, b.p, b.q
    yh = [y[j] + A * h[j] for j in range(3)]
    ag = 1 - k * A * yh[0]; assert ag.lower() > 0
    agp = -(m * A * yh[0] + k * A * yh[1])
    qa = q / ag; qap = (1 + k * p) / ag - q * agp / (ag * ag)
    A1 = -(yh[1] - yh[2] * qa - yh[1] * qap)
    B1 = (abs(yh[2]) + abs(yh[1]) * abs(k)) / ag + abs(yh[1]) * abs(agp) / (ag * ag)
    c0p = -arb(5) / 2 * m * p + k * k * q
    d = {"psi1p": abs(k * q), "psi2p": abs(k), "psi3p": abs(y[1]), "h2p": abs(h[1]), "A1p_over_alpha": abs(A1), "B1p_over_alpha": B1, "c0p": abs(c0p)}
    return {kk: IC.up(v) for kk, v in d.items()}
bl = IC.make_segments([('x', 10.25), ('x', 16.25), ('x', 20.0), ('x', 30.25), ('x', 40.0), ('x', 60.0)], [0.01, 0.01, 0.0125, 0.0125, 0.025])
acc = {}; sups = {}
for n, (sa, sb) in enumerate(bl):
    b = IC.make_box(sa, sb)
    y = PSI.jets_any(b.X, b.xc, 5); h = H2.jets_any(b.X, b.xc, 5)
    for v in y + h: IC.chk(v)
    if b.xb > R + 1e-6: continue
    qd = quantities(b, y, h)
    L = arb(IC.up(arb(b.xb) - arb(b.xa))); Wv = IC.W_up('2', b.xa); z = arb(zeta_sup(b.xa, b.xb))
    for nm, v in qd.items():
        acc[nm] = acc.get(nm, arb(0)) + L * Wv * z * z * arb(v) * arb(v); sups[nm] = max(sups.get(nm, 0.0), v)
    if n % 1000 == 0: print(n, len(bl), f"[{time.time()-T0:.0f}s]", flush=True)
norms = {nm: IC.up(v.sqrt()) for nm, v in acc.items()}
NP12 = IC.up((arb(norms['psi1p']) ** 2 + arb(norms['psi2p']) ** 2).sqrt())
NAB = IC.up((arb(norms['A1p_over_alpha']) ** 2 + arb(norms['B1p_over_alpha']) ** 2).sqrt())
v2 = json.load(open(os.path.join(V2, 'v2_V56_1.json')))['V5']['P=2_R=60']
out = dict(description=__doc__, norms_zeta_C4=norms, NP3=norms['psi3p'], NH2=norms['h2p'], NAB_components=[norms['A1p_over_alpha'], norms['B1p_over_alpha']], NAB=NAB,
           NP12_components=[norms['psi1p'], norms['psi2p']], NP12=NP12, C0P=norms['c0p'], sup_on_support=sups,
           v2_C3_norms_zeta={k_: v2[k_]['norm_zeta'] for k_ in v2}, rel_C4_vs_C3={k_: norms[k_] / v2[k_]['norm_zeta'] - 1 for k_ in norms},
           args=args_block([os.path.join(V2, 'int_common.py'), os.path.join(V2, 'int_v56.py'), os.path.join(V2, 'v2_V56_1.json'), os.path.join(PB, 'c4fun.py')] + EXT, dict(nboxes=len(bl))),
           runtime_s=round(time.time() - T0, 1))
h = dump(os.path.join(OUT, 'c3_modx.json'), out)
print(json.dumps({k_: out[k_] for k_ in ('NP3', 'NH2', 'NAB', 'NP12', 'C0P', 'rel_C4_vs_C3')}, indent=1)); print('sha256', h)
