#!/usr/bin/env python3
"""U-3 (fin13v4): the strip part of the IBP modulation lemma (agent_mod/ibp, alpha1 = 5e-4, eta = 0.02) as a LINEAR FORM in the strip sups
S_r = sup_{sigma in [29,32]} |w^(r)|, r = 0..3, read from mod_ibp.json key alpha_max_5e-4.strip (code: ibp_run.py, ibp_final.py), per component i in {a, b, beta}:
   s_i(S) = |alpha| Lalpha_i S0 + sum_r kappa_ir S_r + sum_{ab} gamma_i,ab S_a S_b + eta [ M3_i/3 (S0+S1)^2 + M2_i (S0+S1) S2 + 2 M1_i S2^2 ]
(linear term: kappa_j = d_alpha Lam is per unit alpha, as m_j; commutator: [L,chi]w pairing; quadratic/cubic forms as in ibp_final.forms_of).
Linearised with the box S0, S1, S2 <= SB = 0.02 on the second factor (S0, S1 <= 0.02 is the lemma's strip hypothesis; S2 <= 0.02 is a new
asserted hypothesis, certified along the march: sup|w''| on [16.25,34.25] <= sI2):
   c_i0 = |alpha| Lalpha_i + kappa_i0 + SB (g00 + g01 + g02) + eta (2 SB M3_i/3 + SB M2_i)
   c_i1 = kappa_i1 + SB (g11 + g12) + eta (2 SB M3_i/3 + SB M2_i)
   c_i2 = kappa_i2 + SB g22 + 2 eta SB M1_i
   c_i3 = kappa_i3
NEUMANN FEEDBACK (correction found in v4): the strip terms are source terms of x = H R + H E x like every other term, so they get the same
feedback as cPhi in ibp_final.finalize: x_i <= s_i + gamma_i max_k s_k/d_k.  As a linear form: cf_ir = c_ir + gamma_i max_k c_kr/d_k (coefficientwise max
is >= the max of the forms).  The lemma's json (and fin12f4's STRIP_BETA, STRIP_AB) omit this feedback.
(a, b): Euclidean norm, Minkowski: |(s_a, s_b)| <= sum_r sqrt(cf_ar^2 + cf_br^2) S_r.
Outputs: coefficient vectors with alpha-part separated (c_i0 = c_i0^const + |alpha| Lalpha_i^fb), and the constants STRIP_BETA, STRIP_AB at S_r = 0.02 (r = 0..3) and
|alpha| = 5e-4 (the form used by the alpha-equation and mod_sups, where S_3 <= 0.02 remains a HYPOTHESIS: SF-2, open).
Output out/u3_strip_form.json.  Run: PYTHONHASHSEED=0 python3 u3_strip_form.py"""
import sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from v4common import args_block, dump, OUT, V4, fup
from flint import arb, ctx
ctx.prec = 200
IBP = os.path.join(V4, 'agent_mod', 'ibp')
P1 = os.path.join(IBP, 'mod_ibp.json'); P2 = os.path.join(IBP, 'final_5e-4.json')
Mj = json.load(open(P1))['alpha_max_5e-4']; F = json.load(open(P2))
A = lambda v: arb(repr(v))
SB = arb('0.02'); ETA = A(Mj['eta']); assert Mj['eta'] == 0.02 and Mj['alpha1'] == 0.0005
comps = ['a', 'b', 'beta']
d = [A(v) for v in Mj['neumann']['d']]
gam = [A(F['final']['core'][c]['gamma_ref']) for c in comps]
assert all(F['final']['core'][c]['gamma_ref'] == Mj['final']['core'][c]['gamma_ref'] for c in comps)
cst = {}; lal = {}
for c in comps:
    s = Mj['strip'][c]; g = {k: A(v) for k, v in s['quad_gamma_L1'].items()}; kap = [A(v) for v in s['commutator_S0_S3']]
    M1, M2, M3 = A(s['cubic_M1_L1']), A(s['cubic_M2_L1']), A(s['cubic_M3_L1'])
    cub01 = ETA * (2 * SB * M3 / 3 + SB * M2)
    cst[c] = [kap[0] + SB * (g['00'] + g['01'] + g['02']) + cub01, kap[1] + SB * (g['11'] + g['12']) + cub01, kap[2] + SB * g['22'] + 2 * ETA * SB * M1, kap[3]]
    lal[c] = A(s['linear_alpha_S0'])
def mx(lst):
    m = lst[0]
    for v in lst[1:]: m = v if v.upper() > m.upper() else m
    return arb(m.upper())
# feedback, coefficientwise (constant part and alpha part separately; both nonnegative)
fb_c = {}; fb_l = {}
for i, c in enumerate(comps):
    fb_c[c] = [cst[c][r] + gam[i] * mx([cst[k][r] / d[j] for j, k in enumerate(comps)]) for r in range(4)]
    fb_l[c] = lal[c] + gam[i] * mx([lal[k] / d[j] for j, k in enumerate(comps)])
ab_c = [(fb_c['a'][r] ** 2 + fb_c['b'][r] ** 2).sqrt() for r in range(4)]; ab_l = (fb_l['a'] ** 2 + fb_l['b'] ** 2).sqrt()
AL = arb('5e-4'); S = arb('0.02')
def const(cv, lv): return sum((cv[r] * S for r in range(4)), arb(0)) + AL * lv * S
nofb_beta = const(cst['beta'], lal['beta']); nofb_ab = sum((((cst['a'][r] ** 2 + cst['b'][r] ** 2).sqrt()) * S for r in range(4)), arb(0)) + AL * (lal['a'] ** 2 + lal['b'] ** 2).sqrt() * S
out = dict(description=__doc__,
           beta=dict(const_coef_S0_S3=[fup(v) for v in fb_c['beta']], alpha_coef_S0=fup(fb_l['beta'])),
           ab=dict(const_coef_S0_S3=[fup(v) for v in ab_c], alpha_coef_S0=fup(ab_l)),
           a=dict(const_coef_S0_S3=[fup(v) for v in fb_c['a']], alpha_coef_S0=fup(fb_l['a'])),
           b=dict(const_coef_S0_S3=[fup(v) for v in fb_c['b']], alpha_coef_S0=fup(fb_l['b'])),
           STRIP_BETA=fup(const(fb_c['beta'], fb_l['beta'])), STRIP_AB=fup(const(ab_c, ab_l)),
           STRIP_BETA_without_feedback=fup(nofb_beta), STRIP_AB_without_feedback=fup(nofb_ab),
           fin12f4_typed=dict(STRIP_BETA=0.02 * 5.3e-8, STRIP_AB=0.02 * 2.76e-7),
           hypotheses=dict(S0_S1_S2_le=0.02, S3_le_for_constants=0.02, alpha_le=5e-4, note='S3 <= 0.02 is needed only for the constants STRIP_BETA/STRIP_AB (alpha-equation, mod_sups); the E_I/E0 state form uses S3 <= supv2'),
           args=args_block([P1, P2]))
h = dump(os.path.join(OUT, 'u3_strip_form.json'), out)
print(json.dumps({k: v for k, v in out.items() if k not in ('description', 'args')}, indent=1)); print('sha256', h)
