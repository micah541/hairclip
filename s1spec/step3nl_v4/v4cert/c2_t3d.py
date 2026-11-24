#!/usr/bin/env python3
"""Stage C2 "t3d" (fin13v4): zone x alpha-level tables of the coefficient sups of N = c w'''' + N_low (EI_SPEC 2.2, 4.3-4.6, 4.9 T12, 6).
Zones (x = sigma - sigma_tip):  IR = [10.25,16.25] (jet box |w|,|w'| <= ETA_IR = 1e-3, boxes as agent_int/planB/pb_ramp.py),
  ZA = [16.25,34.25], ZB = [34.25,60] (box B-c |w'| <= 0.02, B-d |w| <= 0.02 on x < 26 and <= 0.5 beyond; boxes as agent_int/planB/r60_t3.py,
  with an extra breakpoint at 34.25).   alpha-levels AMAX in {2.5e-4, 1e-3} (u0 = alpha psi3~ + alpha^2 h2~, C^4 jets, a5 = u0^(5) from c4fun).
Per box (arb, 200 bits), sup of |.| of:
  base (via pb_nlin.box_vals, the code path of r60_t3.py / pb_ramp.py, so the fin12f4 tables are reproduced):
     c (=d001), d010, d100, d200, d110, d300, A0, A1
  new (sympy, total x-derivative with U0' = U1, U1' = a2 + g', a_j' = a_{j+1}, profile ODE; checked exactly in stage S1):
     X^a := D X |_{g'=0},  X^b := d X/dU1    for X in {c, d010, d100, d200, d110, d300}      (so |D X| <= X^a + X^b |g'|; CA = c^a, CB = c^b)
     B0 := |d E0c/dw|, B1 := |d E0c/dg| over the (w, g) box  (D(d000) = E1c g' + E0c(w, g), E0c(0,0) = 0, |E1c| <= A1)
     kU := |k U0|;  DQ0 := |(dQ, dC)|, DQ1 := |(dQ', dC')|, EA0 := |(E_A, E_B)|, EA1 := |(E_A', E_B')|  (EI_SPEC 4.9 T12; |cos|,|sin| <= 1)
  Expressions that contain q are split q-affinely with the q-coefficient simplified symbolically (as pb_nlin does for A1), so that the
  sigma-width of the q-ball does not enter through a cancelling difference.
Output out/c2_t3d.json (zone tables, sub-region tables, regression against r60_T3_boxw1_0.02.json and planB_ramp.json) + out/c2_t3d_boxes.json
(per-box sups of A0, B0, EA1, kU for the Hardy stage C4).  Checkpoint per (zone, level) in out/c2_ckpt/.
Run: PYTHONHASHSEED=0 nice -n 19 python3 c2_t3d.py"""
import sys, os, json, time, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from v4common import args_block, dump, OUT, EXT, V4, sha
PB = os.path.join(V4, 'agent_int', 'planB')
sys.path.insert(0, PB)
_argv = sys.argv[:]; sys.argv = ['x', '1e-3', '1e-2', '1e-2', '1']
os.chdir(PB)
import pb_nlin as M
sys.argv = _argv
from pb_common import make_box, make_segments, up, chk
from c4fun import PSI, H2
from flint import arb, ctx
import sympy as sp
ctx.prec = 200
T0 = time.time()
CK = os.path.join(OUT, 'c2_ckpt'); os.makedirs(CK, exist_ok=True)
GS = M.GS
Ji, al, k, m, p, q, u0, u1 = sp.symbols('Ji al k m p q u0 u1')
U0, U1, a1, a2, a3, a4, a5, w2, w3, w4 = sp.symbols('U0 U1 a1 a2 a3 a4 a5 w2 w3 w4')
us, w, g, g1, ct, st = sp.symbols('us w g g1 ct st')
# symbolic objects: reuse pb_nlin's Fpert / coef (identical construction, same symbol names) -- one expansion only (RSS budget)
Fpert = M.Fpert; coef = M.coef
import resource; print("after pb_nlin import maxrss MB", resource.getrusage(resource.RUSAGE_SELF).ru_maxrss // 1024, flush=True)
RULES = {U0: U1, U1: a2 + g1, a1: a2, a2: a3, a3: a4, a4: a5, k: m, m: p / 4, p: -k * q, q: 1 + k * p, us: a1, w: g, g: g1, ct: -k * st, st: k * ct}
def Dx(e): return sum(sp.diff(e, v) * r for v, r in RULES.items() if e.has(v))
NAMES = {(0, 0, 1): 'c', (0, 1, 0): 'd010', (1, 0, 0): 'd100', (2, 0, 0): 'd200', (1, 1, 0): 'd110', (3, 0, 0): 'd300'}
DQ_G000 = sp.Poly(GS[(0, 0, 0)], q).degree()          # degree in q of the only q-dependent graph coefficient
DQMAX = max(DQ_G000, 1) + 1                               # Flow has q-degree <= max(deg G000, 1); D raises it by <= 1 (p' = -k q)
def qsplit(e, dmax=None):
    """exact q-Taylor split e = sum_{j<=dmax} q^j e_j, e_j = (d/dq)^j e |_{q=0} / j!  (exact because deg_q e <= dmax).
    The q-ball is the sigma-box (width ~ box width); taking q out avoids its width entering through a cancelling difference."""
    dmax = DQMAX if dmax is None else dmax
    if not e.has(q): return [(0, e)]
    out = []; d = e
    for j in range(dmax + 1):
        cj = d.subs(q, 0) / sp.factorial(j)
        if cj != 0: out.append((j, cj))
        d = sp.diff(d, q)
    return out
ARGS = (U0, U1, us, a1, w, g, a2, a3, a4, a5, k, m, p, q, ct, st)
_Pr = M._Pr
def mk(e): return eval("lambda " + ", ".join(map(str, ARGS)) + ": " + _Pr().doprint(e), {"arb": arb, "pw": M.pw})
import gc, ctypes
def _free():
    sp.core.cache.clear_cache(); gc.collect()
    try: ctypes.CDLL('libc.so.6').malloc_trim(0)
    except Exception: pass
CF = {}
def comp(key, e):
    """compile one expression and drop it (one big expression alive at a time)"""
    CF[key] = [(pw_, mk(c)) for pw_, c in qsplit(e)]
    import resource; print("  compiled", key, "maxrss MB", resource.getrusage(resource.RUSAGE_SELF).ru_maxrss // 1024, flush=True)
    _free()
for mon, nm in NAMES.items():
    Dd = Dx(coef[mon])          # not expanded (RSS); g1 enters only through U1' = a2 + g1, so subs(g1, 0) is exact
    comp(nm + '_a', Dd.subs(g1, 0)); comp(nm + '_b', sp.diff(coef[mon], U1)); del Dd
c0 = sp.Rational(1, 4) - 3 * m ** 2 - k * p
# D(d000) = E1c g' + E0c (stage S1, I14):  E0c = H(U) - H(u0) - c0' w - c0 g + (5km+q/4)' g,  H := D G with U1' = a2;  mean value in (w, g):
#   B0 := sup|H_U0(U) - c0'|,  B1 := sup|H_U1(U) - c0 + (5km+q/4)'|   over the (U0, U1) box
RULES_H = {U0: U1, U1: a2, a2: a3, a3: a4, a4: a5, k: m, m: p / 4, p: -k * q, q: 1 + k * p}
Flow = Fpert.subs({w2: 0, w3: 0, w4: 0})
c0p = -sp.Rational(5, 2) * m * p + k ** 2 * q
dq4p = 5 * (m * m + k * p / 4) + (1 + k * p) / 4
for key, var, extra in (('B0', U0, -c0p), ('B1', U1, -c0 + dq4p)):
    e = sum(sp.diff(sp.diff(Flow, v_), var) * r_ + (sp.diff(r_, var) * sp.diff(Flow, v_)) for v_, r_ in RULES_H.items()) + extra     # d/dvar of H = sum G_v r_v
    comp(key, e); del e
del Flow
# T12 objects (u = U0 = us + w, u' = U1 = a1 + g): dQ = q k u/(1-ku), dC = cos k u/(1-ku), E_A = u0' q k/((1-ku)(1-k u0)), E_B = u0' k cos/(...)
dQ = q * k * U0 / (1 - k * U0); dC = ct * k * U0 / (1 - k * U0)
EA = a1 * q * k / ((1 - k * U0) * (1 - k * us)); EB = a1 * k * ct / ((1 - k * U0) * (1 - k * us))
for key, e in (('dQp', Dx(dQ)), ('dCp', Dx(dC)), ('EAp', Dx(EA)), ('EBp', Dx(EB))):
    assert not e.has(g1), key        # no u'' in the T12 x-derivatives
    comp(key, e)
M.Fpert = M.Pw = M.coef = M.Flow = None; Fpert = coef = None
_free()
print('q-degree of G000:', DQ_G000, 'DQMAX', DQMAX); print('compiled', {k_: len(v) for k_, v in CF.items()}, f"[{time.time()-T0:.0f}s]", flush=True)
def ev(key, vals, qq):
    tot = arb(0)
    for pw_, f in CF[key]:
        t = f(*vals); tot = tot + (t if pw_ == 0 else t * qq ** pw_)
    return chk(tot)
def absup(z): return up(abs(chk(z)))
ONE = arb(0, 1)
def box_new(y, h, b, AMAX, e0, e1):
    A = arb(0, AMAX)
    uj = [A * y[j] + A * A * h[j] for j in range(6)]
    W0, G0 = arb(0, e0), arb(0, e1)
    Ub0, Ub1 = uj[0] + W0, uj[1] + G0
    vals = (Ub0, Ub1, uj[0], uj[1], W0, G0, uj[2], uj[3], uj[4], uj[5], b.k, b.m, b.p, b.q, ONE, ONE)
    o = {}
    for nm in NAMES.values():
        o[nm + '_a'] = absup(ev(nm + '_a', vals, b.q)); o[nm + '_b'] = absup(ev(nm + '_b', vals, b.q))
    o['B0'] = absup(ev('B0', vals, b.q)); o['B1'] = absup(ev('B1', vals, b.q))
    kU = abs(b.k * Ub0); assert kU.upper() < 0.5
    o['kU'] = up(kU)
    den = 1 - kU.upper()
    o['DQ0'] = up(kU.upper() / den * (b.q * b.q + 1).sqrt())
    o['DQ1'] = up((arb(absup(ev('dQp', vals, b.q))) ** 2 + arb(absup(ev('dCp', vals, b.q))) ** 2).sqrt())
    ea = abs(uj[1] * b.k) / ((1 - kU.upper()) * (1 - abs(b.k * uj[0]).upper()))
    o['EA0'] = up(ea * (b.q * b.q + 1).sqrt())
    o['EA1'] = up((arb(absup(ev('EAp', vals, b.q))) ** 2 + arb(absup(ev('EBp', vals, b.q))) ** 2).sqrt())
    return o
BASEKEYS = {"c": (0, 0, 1), "d010": (0, 1, 0), "d100": (1, 0, 0), "d200": (2, 0, 0), "d110": (1, 1, 0), "d300": (3, 0, 0), "A0": "A0", "A1": "A1"}
def run_zone(zname, bl, e0f, e1f, AMAX, subregs):
    ck = os.path.join(CK, f"{zname}_{AMAX:g}.json")
    if os.path.exists(ck): return json.load(open(ck))
    t = time.time(); zs = {}; subs = {}; boxes = []
    for n, (sa, sb) in enumerate(bl):
        b = make_box(sa, sb)
        y = PSI.jets_any(b.X, b.xc, 6); h = H2.jets_any(b.X, b.xc, 6)
        for v in y + h: chk(v)
        M.A = arb(0, AMAX); M.E0b = arb(0, e0f(b.xa)); M.E1b = arb(0, e1f(b.xa))
        d = M.box_vals(y, h, b.p, b.q, b.k, b.m)
        o = {nm: up(d[key]) for nm, key in BASEKEYS.items()}
        o["sup|u0'|"] = up(d["|u0'|"]); o["sup|u0''|"] = up(d["|u0''|"])
        o.update(box_new(y, h, b, AMAX, e0f(b.xa), e1f(b.xa)))
        for kk, v in o.items(): zs[kk] = max(zs.get(kk, 0.0), v)
        for lo, hi in subregs:
            if b.xb > lo + 1e-6 and b.xa < hi - 1e-6:
                s = subs.setdefault(f"[{lo:g},{hi:g}]", {})
                for kk, v in o.items(): s[kk] = max(s.get(kk, 0.0), v)
        boxes.append(dict(xa=b.xa, xb=b.xb, A0=o['A0'], B0=o['B0'], B1=o['B1'], A1=o['A1'], EA1=o['EA1'], EA0=o['EA0'], kU=o['kU']))
        if n % 200 == 0: print(zname, AMAX, n, len(bl), f"[{time.time()-t:.0f}s]", flush=True)
    r = dict(zone=zs, sub=subs, boxes=boxes, nboxes=len(bl), secs=round(time.time() - t, 1))
    json.dump(r, open(ck, 'w'))
    return r
ETA_IR = 1e-3
ZONES = {
 'IR': dict(bl=make_segments([('x', 10.25), ('x', 16.25)], [0.005]), e0=lambda x: ETA_IR, e1=lambda x: ETA_IR, subs=[(10.25, 16.25)]),
 'ZA': dict(bl=make_segments([('x', a) for a in (16.25, 20.0, 26.0, 30.25, 34.25)], [0.025] * 4), e0=lambda x: 0.02 if x < 26 else 0.5, e1=lambda x: 0.02,
            subs=[(16.25, 20.0), (20.0, 26.0), (26.0, 30.25), (30.25, 34.25)]),
 'ZB': dict(bl=make_segments([('x', a) for a in (34.25, 40.0, 60.0)], [0.025, 0.025]), e0=lambda x: 0.5, e1=lambda x: 0.02, subs=[(34.25, 40.0), (40.0, 60.0), (57.0, 60.0)]),
}
LEVELS = (2.5e-4, 1e-3)
if os.environ.get('C2TEST'):          # timing test only: first n boxes per zone, no checkpoint, no output
    nT = int(os.environ['C2TEST'])
    for zn, Z in ZONES.items():
        t_ = time.time(); r_ = None
        CK = os.path.join(OUT, 'c2_ckpt_test'); os.makedirs(CK, exist_ok=True)
        for f_ in os.listdir(CK): os.remove(os.path.join(CK, f_))
        r_ = run_zone(zn, Z['bl'][:nT], Z['e0'], Z['e1'], 2.5e-4, Z['subs'])
        print('TEST', zn, nT, 'boxes', round(time.time() - t_, 2), 's; total boxes', len(Z['bl']), {k_: float('%.3g' % v_) for k_, v_ in r_['zone'].items()}, flush=True)
    sys.exit(0)
res = {}; allboxes = {}
for lev in LEVELS:
    for zn, Z in ZONES.items():
        r = run_zone(zn, Z['bl'], Z['e0'], Z['e1'], lev, Z['subs'])
        res.setdefault(f"{lev:g}", {})[zn] = dict(sup=r['zone'], sub=r['sub'], nboxes=r['nboxes'], secs=r['secs'])
        allboxes[f"{zn}_{lev:g}"] = r['boxes']
# ---------------- regression: fin12f4 tables reproduced by the same code path
reg = {}
r60 = json.load(open(os.path.join(PB, 'r60_T3_boxw1_0.02.json')))['sub_regions_new']
for rg in ('[16.25,20]', '[20,26]', '[26,30.25]', '[40,60]'):
    mine = res['0.001']['ZA']['sub'].get(rg) or res['0.001']['ZB']['sub'].get(rg)
    reg[rg] = {kk: dict(r60=r60[rg][kk], v4=mine[kk if kk != 'd001' else 'c'], equal=(r60[rg][kk] == mine[kk if kk != 'd001' else 'c']), v4_le_old=(mine[kk if kk != 'd001' else 'c'] <= r60[rg][kk]))
               for kk in ('d001', 'd010', 'd100', 'd200', 'd110', 'd300', 'A0', 'A1')}
rt = json.load(open(os.path.join(PB, 'planB_ramp.json')))['table']['alpha=0.00025,eta=0.001']
mineIR = res['0.00025']['IR']['sup']
reg['IR_planB_ramp_2.5e-4_eta1e-3'] = {kk: dict(planB=rt[kk], v4=mineIR[m_], equal=(rt[kk] == mineIR[m_]), v4_le_old=(mineIR[m_] <= rt[kk])) for kk, m_ in
                                       (('d001_w4', 'c'), ('d010_w3', 'd010'), ('d100_w2', 'd100'), ('A1_w1', 'A1'), ('A0_w0', 'A0'))}
reg_ok = all(v['v4_le_old'] for rg in reg.values() for v in rg.values())     # same code path (pb_nlin.box_vals); v4 evaluates at 200 bits, the fin12f4 tables at 53 bits
reg_n_equal = sum(v['equal'] for rg in reg.values() for v in rg.values()); reg_n = sum(len(rg) for rg in reg.values())
ins = [os.path.join(PB, f) for f in ('pb_nlin.py', 'pb_common.py', 'c4fun.py', 'r60_T3_boxw1_0.02.json', 'planB_ramp.json')] + [os.path.join(V4, 'graph_struct.pkl'), '/home/claude/s1spec/step3nl/graph_struct.pkl'] + EXT
out = dict(description=__doc__, levels=list(LEVELS), eta_IR=ETA_IR, tables=res, regression=reg, regression_all_equal=reg_ok, regression_n_equal=[reg_n_equal, reg_n],
           args=args_block(ins, dict(prec=ctx.prec)), runtime_s=round(time.time() - T0, 1))
h = dump(os.path.join(OUT, 'c2_t3d.json'), out)
h2 = dump(os.path.join(OUT, 'c2_t3d_boxes.json'), dict(boxes=allboxes, args=out['args']))
for lev in res:
    for zn in res[lev]: print(lev, zn, {kk: float('%.4g' % v) for kk, v in res[lev][zn]['sup'].items()})
print('regression v4 <= fin12f4 tables:', reg_ok, 'bitwise equal', reg_n_equal, 'of', reg_n); print('sha256', h, h2)
