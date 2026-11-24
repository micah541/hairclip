#!/usr/bin/env python3
"""Lemma U, Phase 2a (agent B): certified (arb) key constants of the Ghat design study.

  PYTHONHASHSEED=0 nice -n 19 python3 -B ghat_design_arb.py        (single process, < 200 MB, ~1-3 min)

What is certified [cert]:
  Ghat_jj (j = 0,1,2) of the general-weight pipeline of phase1/ghat (functions ghat_arb.YC_base, l4_extra_YC_x, budget_ok, ghat_x, used
  unchanged), with
    - Gamma_m, kappa~, kappa2~ from weight_consts_arb.certify(4, 4, a in [0, 3/10], theta_a <= 3/5, y_L = 56) (P1; monotone in theta_a, so valid for
      every damping constant g with g/(28^4 + g) <= 3/5);
    - Ab INCLUDING the profile jets U_j (j >= 2) (CHECKPOINT finding 3: Ab := jet-free part of At, as P2 / Fable's reconciliation recon.py);
    - dfrac = 2/25 and the ramp-foot Young adds (P2 conventions; the port is Fable's recon.py `run(..., foot=True)`, copied verbatim below);
    - sup|z_0| <= M0z = 4/25 (>= 0.152, CHECKPOINT conventions), profile bound Umax = 9/100 (support value) or 1/20 (zones I/II only).
  For each box (M0, M1) the design quantities, in arb with outward rounding:
    design 'global28' (Phase 1): g = G/(1 - G/28^4) (requires G < 28^4);  Gw = D(28) = 1/(1 - G/28^4);  rate margin 28^4/G.
    design 'data56' (this note, U_ramp.tex Sec. 4):  g = G/(1 - 4 abar - G/56^4) (requires G < (1 - 4 abar) 56^4);
        D(56) = 1 + g/56^4,  D(28) = 1 + g/28^4,  rate margin (1 - 4 abar) 56^4 / G,
        ramp-share data coefficient  c_rs = G g / 28^8  (multiplies the sweep integral S_rs of U_model2.py; U_ramp.tex eq. (rs)).
abar = 2.19103e-7 rounded up to 2.192e-7 (D1 of v3.2, certified sup |a|).
Output: ghat_design_arb.json (with args block: argv, env, input sha256), ghat_design_arb.out (stdout)."""
import sys, os, json, time, math, hashlib, pickle, platform
sys.dont_write_bytecode = True
T0 = time.time()
HERE = os.path.dirname(os.path.abspath(__file__))
GH = '/home/claude/s1spec/writeup/lemmaU/phase1/ghat'
sys.path.insert(0, GH)
import sympy as sp
from flint import arb, ctx
import ghat_arb as GA, weight_consts_arb as WC
ctx.prec = 128
Q_ = sp.Rational; A = GA.A

def sha(p): return hashlib.sha256(open(p, 'rb').read()).hexdigest()
def fup(x):
    x = x.upper() if isinstance(x, arb) else arb(x); v = float(x)
    return v if arb(v) >= x else math.nextafter(v, math.inf)
def flo(x):
    x = x.lower() if isinstance(x, arb) else arb(x); v = float(x)
    return v if arb(v) <= x else math.nextafter(v, -math.inf)

T = pickle.load(open(GA.TABLES, 'rb'))
r = WC.certify(4, 4, Q_(0), Q_(3, 10), Q_(3, 5), 56); b = r['bound']
Gam = {0: arb(1), 1: b['G1'], 2: b['G2'], 3: b['G3']}; kap = b['kap']; kap2 = b['kap2']
LN2 = arb(2).log(); etap = A(Q_(15, 8))/LN2; DEL = A(Q_(1, 10))

# ---- verbatim port of Fable's recon.py run() (phase1/checkpoint/fable_scratch/recon.py), returning arb Ghat_jj
def run(M0q, M1q, Umq, M0zq, dfrac, foot, honest=True, s33='ibp'):
    Asym = GA.build_At_sym(Umq); Ab = GA.ab_exact(M0q, M1q, 4, Umq); At = GA.eval_At(Asym, M0q, M1q); out = []
    if honest: Ab = {k: d.get((), arb(0)) for k, d in At.items()}
    for j in (0, 1, 2):
        YCb = GA.YC_base(T[j], j, Ab, Gam); add, dis, spare, ndfr, items = GA.l4_extra_YC_x(T, j, M0q, M1q, At, Gam, kap, dfrac=dfrac, s33=s33, M0zq=M0zq)
        fadd = {k: arb(0) for k in add}; extra = Q_(0); M1 = A(M1q)
        if foot:
            w33 = sum((it[1] for it in items if it[0] == 'S33'), arb(0))
            if any(it[0] == 'S33' for it in items):
                c1 = w33*M1*etap; fadd[3] += c1*A(kap)*DEL; extra = Q_(1, 10)
            for it in items:
                if it[0] != 'L4': continue
                w = it[2]
                if it[1] == 3: fadd[3] += w*M1*etap*DEL
                elif it[1] == 2: fadd[2] += w*M1*etap
        ok, sl = GA.budget_ok(dis, ndfr, M0q, extra)
        YC = {k: YCb[k] + add.get(k, 0) + fadd.get(k, 0) for k in YCb}
        G = GA.ghat_x(YC, j, M0q, kap2)
        out.append(dict(G=G, ok=bool(ok), slack=sl))
    return out

ABAR = A(Q_(2192, 10**10))
def design(G22):
    e28 = arb(28)**4; e56 = arb(56)**4; d = {}
    d['global28'] = dict(feasible=bool((G22 < e28)), g=fup(G22/(1 - G22/e28)), Gw_D28=fup(1/(1 - G22/e28)), rate_margin=flo(e28/G22))
    den = 1 - 4*ABAR - G22/e56
    g = G22/den
    d['data56'] = dict(feasible=bool(den > 0), g=fup(g), D56=fup(1 + g/e56), D28=fup(1 + g/e28), rate_margin=flo((1 - 4*ABAR)*e56/G22),
                       c_rs=fup(G22*g/e28**2), theta56=fup(G22/e56), theta28=fup(g/(e28 + g)),
                       Gamma_valid=bool(den > 0 and g/(e28 + g) <= A(Q_(3, 5))))   # P1's Gamma_m, kappa~ certified for theta(28) <= 3/5 only
    d['global28']['Gamma_valid'] = bool(G22 < e28 and G22/e28 <= A(Q_(3, 5)))
    return d

CASES = [  # name, M0, M1, Umax, M0z
    ('crude (0.11,0.52), U .09, M0z .16', Q_(11, 100), Q_(52, 100), Q_(9, 100), Q_(4, 25)),
    ('crude (0.11,0.52), U .05 (zones I/II), M0z .16', Q_(11, 100), Q_(52, 100), Q_(1, 20), Q_(4, 25)),
    ('widened (0.15,0.52), U .09, M0z .16', Q_(15, 100), Q_(52, 100), Q_(9, 100), Q_(4, 25)),
    ('widened (0.15,0.60), U .09, M0z .16', Q_(15, 100), Q_(60, 100), Q_(9, 100), Q_(4, 25)),
    ('widened (0.20,0.60), U .09, M0z .20', Q_(20, 100), Q_(60, 100), Q_(9, 100), Q_(1, 5)),
    ('zone-III-like (0.30,0.80), U .09, M0z .30', Q_(30, 100), Q_(80, 100), Q_(9, 100), Q_(3, 10)),
    ('Lemma-O-like (0.50,1.00), U .09, M0z .50', Q_(50, 100), Q_(1), Q_(9, 100), Q_(1, 2)),
]
res = dict(consts=dict(G1=fup(Gam[1]), G2=fup(Gam[2]), G3=fup(Gam[3]), kap=fup(kap), kap2=fup(kap2), theta_a_max='3/5', abar_used='2.192e-7'), cases={})
print('Gamma_m, kappa~, kappa2~ (P1 certify, theta_a <= 3/5):', [fup(Gam[m]) for m in (1, 2, 3)], fup(kap), fup(kap2), flush=True)
for name, M0, M1, Um, M0z in CASES:
    o = run(M0, M1, Um, M0z, Q_(2, 25), True)
    G = [o[j]['G'][j] for j in range(3)]
    row = dict(box=[str(M0), str(M1)], Umax=str(Um), M0z=str(M0z), dfrac='2/25', foot=True, Ab='jet-free part of At (incl. U_j)',
               Ghat_jj=[[flo(x), fup(x)] for x in G], Ghat_full={str(j): {str(k): fup(o[j]['G'][k]) for k in range(j + 1)} for j in range(3)},
               budget_ok=[o[j]['ok'] for j in range(3)], budget_slack=[o[j]['slack'] for j in range(3)], design=design(G[2]))
    res['cases'][name] = row
    dd = row['design']
    print(f"{name:48s} Ghat_jj = {fup(G[0]):.1f} {fup(G[1]):.1f} {fup(G[2]):.1f}  budget {row['budget_ok']} slack {row['budget_slack']}\n"
          f"    global28: g {dd['global28']['g']:.4g}  Gw=D(28) {dd['global28']['Gw_D28']:.5f}  rate margin x{dd['global28']['rate_margin']:.4f}\n"
          f"    data56  : g {dd['data56']['g']:.4g}  D(56) {dd['data56']['D56']:.5f}  D(28) {dd['data56']['D28']:.5f}  rate margin x{dd['data56']['rate_margin']:.3f}  c_rs = G g/28^8 = {dd['data56']['c_rs']:.5f}  Gamma valid (theta28 <= 3/5): global28 {dd['global28']['Gamma_valid']}, data56 {dd['data56']['Gamma_valid']}", flush=True)
base = res['cases']['crude (0.11,0.52), U .09, M0z .16']['Ghat_jj'][2][1]; wid = res['cases']['widened (0.15,0.52), U .09, M0z .16']['Ghat_jj'][2][1]
res['ratio_widened_over_crude_G22'] = wid/base
print(f'Ghat22 widened (0.15) / crude (0.11) = {wid/base:.4f}  (reviewer [rev]: x1.26)')
res['args'] = dict(argv=sys.argv, cwd=os.getcwd(), python=sys.version, platform=platform.platform(), flint=__import__('flint').__version__, prec=ctx.prec,
                   env={k: os.environ.get(k) for k in ('PYTHONHASHSEED', 'PYTHONDONTWRITEBYTECODE')},
                   inputs={p: sha(p) for p in (GA.TABLES, os.path.join(GH, 'ghat_arb.py'), os.path.join(GH, 'weight_consts_arb.py'), os.path.abspath(__file__))})
res['wall_s'] = round(time.time() - T0, 1)
json.dump(res, open(os.path.join(HERE, 'ghat_design_arb.json.tmp'), 'w'), indent=1)
os.replace(os.path.join(HERE, 'ghat_design_arb.json.tmp'), os.path.join(HERE, 'ghat_design_arb.json'))
print('wall', res['wall_s'], 's')
