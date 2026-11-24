#!/usr/bin/env python3
"""GS_boot.py -- [cert] strict improvement of (GS) under the step-3 bootstrap (fix round 14; GS_lemma.tex, Lemma GS-2).
[v38 copy of <s1spec>/writeup/lemmaU/bulk/fix15_review/c/GS_boot_v4c.py for the step-4 v38 re-pin (p2/v38/final_tail_v38.py); changes: output path out_v38/iface/GS_boot_v4c.json; cpu_s recorded as null (time printed to the log)]

Under BS_U at tau (with (GS) at tau, so that u and w = u - alpha psi3^ - alpha^2 h2^ are defined on |x| <= 60), the verifier fin13v4b asserts
in closure.box_assert, at every march step and at the tail anchor (closure.py:143, 162), the box hypotheses below (their max ratios along
the production run alpha1 = 1.5e-4 are read from runs/P/rate_table_fin13v4c_a0.00015.json and asserted <= 1 here):
   |alpha| <= AMAX_E0 = 2.5e-4                       (E0_alpha; also Bf_mod_alpha <= 5e-4)
   sup(|w|, |w'|) <= ETA_TIP_E2 = 1e-4 on [0, 8]     (E2_tip_jet, H-E2b)
   |w| <= WB = 0.02 on [0, 34.25]                     (E0_w_all34)
   sup(|w|, |w'|) <= ETA_IR = 1e-3 on [10.25, 16.25]  (Bb_IR_jet)
   |w'| <= BOX_W1 = 0.02 on [16.25, 60]               (Bc_w1_far_EI, Bc_w1_far_S0V)
   |w| <= BOX_W0B = 0.5 on [34.25, 60]                (Bd_w_34_60)
 and on [8, 10.25] (no asserted w' box): |w'(x)| <= |w'(10.25)| + (2.25)^(1/2) ||w''||_{L^2(8,10.25)} <= ETA_IR + 1.5 B / rho(10.25)^(1/2),
   B = ||w~''||_rho <= (Lambda/eps)^(1/2) <= sqrtL_max/eps^(1/2) (B1: Lambda below the verifier barrier, whose max is sqrtL_max), chi = 1 there,
   rho = exp(-(1 + x^2)^(2/3)/4) decreasing.
With u = alpha psi3^ + alpha^2 h2^ + w, per cell of width 1/64 of [0, 60] (aligned with 8, 10.25, 16.25, 34.25, 60), in arb:
   sup|k_S1 u| <= K_cell (A |psi3^| + A^2 |h2^| + W0) < 1/2,   sup|u'| <= A |psi3^'| + A^2 |h2^'| + W1 < 1   (asserted),
psi3^, h2^ the C^4 functions of c4fun.py (cell-ball jets), K_cell from the S1 rows (s1chart2.chart_ball, tip-position ball included) for
x <= x13 and the step-1 tail envelope hstail.profile_bounds (decreasing) beyond.  Evenness gives |x| <= 60.
Read-only.  Output GS_boot.json.  Run: nice -n 19 python3 -B GS_boot.py"""
import sys, os, json, hashlib, platform, math, re, time
sys.dont_write_bytecode = True
sys.path.insert(0, '/home/claude/s1spec/step4/p2/t9h'); sys.path.insert(0, '/home/claude/s1spec/step3')
sys.path.insert(0, '/home/claude/s1spec/step3nl_v4c/agent_int/v3')
from flint import arb, ctx
ctx.prec = 128
import s1chart2
import hstail as HT
import c4fun as C4
HERE = os.path.dirname(os.path.abspath(__file__))
V4B = '/home/claude/s1spec/step3nl_v4c'
RUN = f'{V4B}/runs/P/rate_table_fin13v4c_a0.00015.json'
E0C = f'{V4B}/v4bcert/out/e0v4b_constants.json'; E2C = f'{V4B}/v4bcert/out/e2v4b_constants.json'
CONSTPY = f'{V4B}/verifier/constants.py'; CLOSPY = f'{V4B}/verifier/closure.py'
U = lambda v: math.nextafter(float(abs(v).upper()), math.inf)
t0 = time.time()


def sha(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()


run = json.load(open(RUN)); e0 = json.load(open(E0C))['constants']; e2 = json.load(open(E2C))['constants']
mr = run['box_hypotheses']['max_ratio']
need = ('E0_alpha', 'E2_tip_jet', 'E0_w_all34', 'Bb_IR_jet', 'Bc_w1_far_EI', 'Bc_w1_far_S0V', 'Bd_w_34_60', 'Bf_mod_alpha') + ('CF_w_16_31', 'CF_w1_16_31', 'CF_alpha')
for n in need: assert mr[n] <= 1.0, (n, mr[n])
src = open(CONSTPY).read()
lit = {}
for name in ('ETA_IR', 'BOX_W1', 'BOX_W0B'):
    m = re.search(r'\b%s\s*=\s*([0-9.eE+-]+)' % name, src); lit[name] = float(m.group(1))
assert lit == dict(ETA_IR=1e-3, BOX_W1=0.02, BOX_W0B=0.5), lit
assert 'box_assert(S_sup, FS, aq, log, f"step {k}")' in open(CLOSPY).read() and '"tail anchor"' in open(CLOSPY).read()
AM = arb(repr(e0['AMAX_E0']['value'])); WB = arb(repr(e0['WB']['value'])); ETIP = arb(repr(e2['ETA_TIP_E2']['value']))
assert float(AM.mid()) == 2.5e-4 and float(WB.mid()) == 0.02 and float(ETIP.mid()) == 1e-4
EIR, BW1, BW0B = arb('1e-3'), arb('0.02'), arb('0.5')
SQL = arb(repr(run['sqrtL_max'])); EPS = arb('1e-5')
Bmax = SQL / EPS.sqrt()
rho1025 = (-(1 + arb('10.25') ** 2) ** (arb(2) / 3) / 4).exp()
W1_8_10 = EIR + arb('1.5') * Bmax / rho1025.sqrt()

S = s1chart2.S1(); X13 = float((s1chart2.T0 + 5).lower()) - 1e-6


def kcell(xl, xu):
    K = 0.0
    if xl < X13:
        K = max(K, U(S.chart_ball(arb(xl).union(arb(min(xu, X13))))['k']))
    if xu > X13:
        sig = arb(max(xl, X13)) - s1chart2.T0 + 8 - arb(1e-6)       # sigma (lower end; envelopes decrease)
        K = max(K, U(HT.profile_bounds(arb(max(13.0, float(sig.lower()))))[0]))
    return K


def boxes(xl, xu):
    if xu <= 8: return ETIP, ETIP
    if xu <= 10.25: return WB, W1_8_10
    if xu <= 16.25: return EIR, EIR
    if xu <= 34.25: return WB, BW1
    return BW0B, BW1


h = 1.0 / 64; n = int(round(60 / h)); worst = dict(ku=0.0, up=0.0); regions = {}
for j in range(n):
    xl, xu = j * h, (j + 1) * h
    X = arb(xl).union(arb(xu)); xm = (xl + xu) / 2
    ps = C4.PSI.jets_any(X, xm, 2); hh = C4.H2.jets_any(X, xm, 2)
    W0, W1 = boxes(xl, xu)
    U0 = AM * abs(ps[0]) + AM * AM * abs(hh[0]) + W0
    U1 = AM * abs(ps[1]) + AM * AM * abs(hh[1]) + W1
    K = kcell(xl, xu)
    ku = math.nextafter(K * U(U0), math.inf) * (1 + 1e-12); up = U(U1)
    assert ku < 0.5 and up < 1.0, (xl, xu, ku, up)
    reg = '0-8' if xu <= 8 else '8-10.25' if xu <= 10.25 else '10.25-16.25' if xu <= 16.25 else '16.25-34.25' if xu <= 34.25 else '34.25-60'
    r_ = regions.setdefault(reg, dict(ku=0.0, up=0.0, sup_u=0.0)); r_['ku'] = max(r_['ku'], ku); r_['up'] = max(r_['up'], up); r_['sup_u'] = max(r_['sup_u'], U(U0))
    worst['ku'] = max(worst['ku'], ku); worst['up'] = max(worst['up'], up)
ok = worst['ku'] < 0.5 and worst['up'] < 1.0
print('GS strictness under the bootstrap on |x| <= 60:', worst, regions, 'W1 on [8,10.25] <=', U(W1_8_10), 'ok', ok, f'({time.time() - t0:.0f}s)')
me = os.path.abspath(__file__)
inputs = {p_: sha(p_) for p_ in (RUN, E0C, E2C, CONSTPY, CLOSPY, '/home/claude/s1_tight/prof/prof_tight.txt', '/home/claude/s1_tight/ext813/ext813_prof.txt')}
for m in ('s1chart2', 'hstail', 'c4fun', 'psi3fun', 'h2fun'):
    mod = sys.modules[m]; inputs[mod.__file__] = sha(mod.__file__)
out = dict(label='[cert]', ok=bool(ok), sup_ku=worst['ku'], sup_uprime=worst['up'], regions=regions, W1_on_8_10=U(W1_8_10),
           box_max_ratios={k_: mr[k_] for k_ in need}, sqrtL_max=run['sqrtL_max'], B_max=U(Bmax),
           statement='under the asserted box hypotheses of fin13v4b and B1, every (GS) graph function u on |x| <= 60 satisfies sup|k u| <= sup_ku < 1/2 and sup|u\'| <= sup_uprime < 1',
           args=dict(argv=sys.argv, env={}, python=platform.python_version(), cell=h, script=dict(path=me, sha256=sha(me)), inputs=inputs, cpu_s=None))   # v38: wall time printed to the stage log only
print(f'cpu_s {time.time() - t0:.1f}')
json.dump(out, open('/home/claude/s1spec/step4/p2/out_v38/iface/GS_boot_v4c.json', 'w'), indent=1)
