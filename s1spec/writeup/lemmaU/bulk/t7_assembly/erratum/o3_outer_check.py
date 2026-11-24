#!/usr/bin/env python3
"""Erratum (task 7a), item A8: the repaired Lemma S(A) evaluated with the outer data that Hypothesis (O) as worded implies.

(O3) |z_k| <= 1 on O = [350, 700] (s = ln y) gives  Nbar_k^2 := sup_tau int_O z_k^2 ds <= ln 2  for k = 0, 1, 2  [paper].
t2_A / t2_A_v33 evaluate the affine export at the reference N2_ref = (0.1, 0.2, 5.0) [est] only. This script evaluates the same certified
affine rule (A_interface.json res.<setting>.EA_of_N2) at Nbar^2 = (ln 2, ln 2, ln 2), then the flux constants K0, K2 (lemmaS.tex:81,
assembly_arb.py:163-179: K = W conv (k_0 57.75 E0 + k_1 E1/54.75 + k_2 2(E2 + E1)/54.75^3), all coefficients >= 0) and the margins
FLUX/(beta K) at beta = 2 (the factor of B2 in v3.1/v3.2) and 5/4.
Checks: the rule at N2_ref reproduces EA_ref, and the K formula at EA_ref reproduces core.K0_at_E / K2_at_E (asserted, both as enclosures).
Inputs are read only. Doubles labelled 'up' are moved one ulp up and 'lo' one ulp down before use (outward).
Run: nice -n 19 python3 -I -B o3_outer_check.py  ->  o3_outer_check.json (with an args block)."""
import json, math, os, sys, hashlib, platform
from flint import arb, ctx
import flint
ctx.prec = 128
HERE = os.path.dirname(os.path.abspath(__file__))
BULK = '/home/claude/s1spec/writeup/lemmaU/bulk'
TREES = {'t2_A (fin11 verifier, H2f design 5e-7)': (BULK + '/t2_A', dict(FLUX0='4.1e-21', FLUX2='2.0e-20')),
         't2_A_v33 (fin12f4 verifier, H2f design 1e-6)': (BULK + '/t2_A_v33', dict(FLUX0='1.64e-20', FLUX2='8.0e-20'))}

def sha(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()

def num(x):
    """certified export entry -> arb, outward by one ulp of its stored double"""
    d = float(x['double']); b = x.get('bound', 'up')
    if b == 'up': return arb(repr(math.nextafter(d, math.inf)))
    if b == 'lo': return arb(repr(math.nextafter(d, -math.inf)))
    raise ValueError(b)

def UP(d):
    return arb(repr(math.nextafter(float(d), math.inf)))

def EA_of(r, G, N2):
    EA0 = [num(x) for x in r['EA0']]; gp = [num(x) for x in r['gp']]
    ko = [[num(x) for x in row] for row in r['kappa_o']]; q = arb(28) ** -4; E = []
    for j in range(3):
        o = sum((ko[j][k] * N2[k] for k in range(3)), arb(0))
        cp = sum((G[j][k] * q * (E[k] - EA0[k]) for k in range(j)), arb(0))
        E.append(EA0[j] + (o + 2 * cp) / gp[j])
    return E

def K_of(E, c):
    W = UP(c['W57_upper']); k0 = [UP(x) for x in c['kappa_I0']]; k2 = [UP(x) for x in c['kappa_I2']]
    conv = arb(60) / arb(28); Eu = [conv * x for x in E]
    N = [W * arb('57.75') * Eu[0], W * Eu[1] / arb('54.75'), W * 2 * (Eu[2] + Eu[1]) / arb('54.75') ** 3]
    return sum((k0[i] * N[i] for i in range(3)), arb(0)), sum((k2[i] * N[i] for i in range(3)), arb(0))

def up(x): return math.nextafter(float(x.upper()), math.inf)
def lo(x): return math.nextafter(float(x.lower()), -math.inf)

inputs = {}; out = {}
L2 = arb(2).log()
for name, (tree, FL) in TREES.items():
    pA = tree + '/out/A_interface.json'; pc = tree + '/base_ibp/lemmaS_cflux.json'
    inputs[pA] = sha(pA); inputs[pc] = sha(pc)
    d = json.load(open(pA)); c = json.load(open(pc))
    G = [[num(x) for x in row] for row in d['interface28']['Ghat_core']]
    Nref = [arb(x) for x in d['N2_ref']['values']]
    res = {}
    for st in ('tight', 'v31'):
        r = d['res'][st]; core = r['core']
        Eref = EA_of(r, G, Nref); EAref = [num(x) for x in r['EA_ref']]
        for a, b in zip(Eref, EAref):
            assert abs(float(a.mid()) - float(b.mid())) <= 1e-12 * abs(float(b.mid())), ('rule does not reproduce EA_ref', st, a, b)
        K0r, K2r = K_of(EAref, c)
        assert abs(float(K0r.mid()) / float(core['K0_at_E']['double']) - 1) < 1e-12 and abs(float(K2r.mid()) / float(core['K2_at_E']['double']) - 1) < 1e-12, ('K formula', st)
        N2o = [L2, L2, L2]
        E3 = EA_of(r, G, N2o); K0, K2 = K_of(E3, c)
        F0, F2 = arb(FL['FLUX0']), arb(FL['FLUX2'])
        row = dict(H2=r['IC']['H2'], H2f=r['IC']['H2f'],
                   EA_at_N2ref=[up(x) for x in EAref], EA_at_O3=[up(x) for x in E3],
                   rel_increase_vs_N2ref=[up((a - b) / b) for a, b in zip(E3, EAref)],
                   K0_at_O3=up(K0), K2_at_O3=up(K2), K0_at_N2ref=up(K0r), K2_at_N2ref=up(K2r),
                   FLUX0=FL['FLUX0'], FLUX2=FL['FLUX2'],
                   margin_beta2=dict(FLUX0_over_2K0=lo(F0 / (2 * K0)), FLUX2_over_2K2=lo(F2 / (2 * K2))),
                   margin_beta5_4=dict(FLUX0_over_K0=lo(F0 / (arb(5) / 4 * K0)), FLUX2_over_K2=lo(F2 / (arb(5) / 4 * K2))),
                   margin_beta2_at_N2ref=dict(FLUX0_over_2K0=lo(F0 / (2 * K0r)), FLUX2_over_2K2=lo(F2 / (2 * K2r))))
        res[st] = row
        print(f"{name} [{st}] E(O3) = {[f'{up(x):.6g}' for x in E3]}  (N2ref: {[f'{up(x):.6g}' for x in EAref]});  "
              f"K0 = {up(K0):.4e}, K2 = {up(K2):.4e};  beta=2: FLUX0/2K0 = {lo(F0 / (2 * K0)):.4f}, FLUX2/2K2 = {lo(F2 / (2 * K2)):.4f}", flush=True)
    out[name] = res
me = os.path.abspath(__file__); inputs[me] = sha(me)
args = dict(argv=sys.argv, cwd=os.getcwd(), python=sys.version.split()[0], platform=platform.platform(), flint=flint.__version__,
            env={k: os.environ.get(k) for k in ('PYTHONHASHSEED', 'OMP_NUM_THREADS') if k in os.environ}, inputs=inputs)
doc = dict(labels='cert = arb with outward rounding on certified exports; the implication (O3) => Nbar_k^2 <= ln 2 is [paper]',
           note='Nbar^2 = (ln 2, ln 2, ln 2); everything else as in the source export (kappa_J applied, A5 interface, point terms, R5prime).',
           res=out, args=args)
tmp = os.path.join(HERE, 'o3_outer_check.json.tmp'); json.dump(doc, open(tmp, 'w'), indent=1); os.replace(tmp, os.path.join(HERE, 'o3_outer_check.json'))
