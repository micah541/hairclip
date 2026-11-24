#!/usr/bin/env python3
"""v38 copy of p2/v33/f12_sliver.py (step-4 v38 stage 20 f12_sliver, p2/v38/final_tail_v38.py): the same F1/F2 check, repointed to the
verifier fin13v4c production run (step3nl_v4c/runs/P; its cert inputs incl. v4bcert/out/strip_cubic_v4c.json, every rate-table pin checked
against the disk) and the (A) chain v37 (writeup/lemmaU/bulk/fix13_v37chain: t2_A, t4_cert, t5_zone3, n1).  v38 changes, all listed:
  paths (verifier, chain, n1 script/output, sob, constants.py); the M1 value (A) reads is recomputed as (A) does (4 digits, up: 0.4329 for
  fin13v4c, was the literal 0.4328); U_cert_for_zone3 via the chain's make_ucert_for_zone3.py (repro validated against U_cert_for_zone3_v36);
  the U/zone3 harness argv are those of the v37 chain (U pass 3, zone3 from pass 2); the stage-19 output out_v38/n1_chart.json must equal
  the chain's n1_chart_fin13v4c.json (all keys but args); quoted thresholds (x1.1847, x3.365/x5.955, K6 x1.189, B2' x2.116) are replaced by
  the pinned v37 values rounded DOWN to 4 significant digits (the text convention for margins); the lemmaU_v34.tex quoted-number check is
  kept as information only (that text is superseded; lemmaU_v37 has no constants_needed file); the macro check uses the v37 generator
  fix12_texts/tools/gen_numbers_A_v37.py with every source repointed to chain v37 / fin13v4c (A, A1, A4, U, Z swapped pinned -> new).
  Closed-tree snapshot extended to the v37 chain, step3nl_v4c verifier/runs/P/v4bcert/out/agent_sob, fix12 tools.
--- original v3.3 docstring follows ---
v3.3 review items F1, F2 (t8_review/CHECK_N1_N5.md): standalone pinned check, step-4 v33 stage f12_sliver.  No pinned output is changed;
every recomputation is written under p2/out_v38/f12/ only.  python-flint arb, outward rounding.

F1  The zeta < 1 sliver sups G_i (i = 0..3) of (A)'s ramp conversion (base_ibp/rigor/assembly_arb.py: sliver_G) were taken on y in [14, 16.25 - X0];
    after N1 the w-part is integrated over sigma in Y^{-1}([14, 28]) inside [14 - 2^-19, 28 + 2^-19].  Here: G_i on the EXTENDED sliver
    y in [14 - 2^-19, 16.25 - X0lo] (X0lo = 2.2490776505 <= -sigma_tip, so it also covers the S1-chart tip ball), same formula and cells as sliver_G
    (G_i = sum_m S(i,m) ymax^m sup|d^m g/dx^m|, g = 1/zeta, zeta = S7((x - 10.25)/6)); the old G is recomputed by the same code and asserted equal
    (bitwise, arb mid and radius) to sliver_G() inside every assembly child.  Also reported: G on the needed domain 14 - e_ramp (e_ramp from n1_chart).
F2  sup|w''| on the sliver x in [16.249, 16.25] (contains every x = sigma - sigma_tip with sigma >= 14 - 2^-19).  On the sliver w' = g v, v = zeta w'
    the verifier's cut-off variable, so |w''| <= sup|g| |v'| + sup|g_x| |v|.  The verifier's far sups sI1, sI2 are the (I2_W) bound
    sqrt(K_W (Psi1 V^2 + 2 V V')) of the full weighted norms of v, v' (closure.py farsup); (I2_W) holds at EVERY x in [6, R] for v supported in
    (6, R) with Psi1(x) = sup_[x,R] psi' and the factor 1/W(x).  On the sliver 1/W(x) <= (1 + 16.25^2)^2 <= K_W = 2.323e6 (W = (1+x^2)^-2 there,
    certified below) and Psi1(x) <= 1.305 (agent_sob sob_constants_P2.json, Psi1 R=60 from_10.25), so sup_sliver |v'| <= sI2, sup_sliver |v| <= sI1
    with the SAME numbers.  Hence sup_sliver |w''| <= (1 + G0x) sI2 + G1x sI1 =: sI2_sl, sI2' := max(sI2, sI2_sl) (rounded up to a double).
Consequences (no rerun of the pinned chain; scratch recomputations from pinned inputs only):
  n1_chart with W2SUP = sI2' (exec of the pinned script with that one line replaced) -> f12/new/n1_chart.json;
  (A) assembly_rep.py, the 7 runs of run_t2.sh, with sliver_G -> extended G and --chart f12/new/n1_chart.json -> f12/new/asm_*.json;
  export_t2 (asm paths and output dir redirected) -> f12/new/A_interface.json;  assemble_U.py --A that, --out f12 -> f12/new/U_cert.json;
  zone3.py (exec, output path replaced) with the new U-fine Bhat -> f12/new/zone3_interface.json;  assemble_U pass 2 with that zone3 (fixed point).
  Every harness is first validated: the same harness with the pinned G / sI2 reproduces the pinned outputs (all keys except args) -> f12/repro/.
  Verifier M1 = 28 sI2' stays <= the 4-digit value 0.4328 that (A) reads (asserted), so (A)'s ramp boxes use a valid number unchanged.
Checks: every claim / margin / flag of the asm, A_interface, U, zone3 outputs against the pinned ones; the quoted numbers and side conditions of
  check_keys_v34.py (lemmaU_v34.tex) and the macros of numbers_A_v33.tex (n35_gen_numbers.py) evaluated on the new files; K6; M3/M6 (n35_asserts);
  a generic leaf diff (every changed number, relative change, and whether its outward rounding to 3..6 significant digits changes).
usage (step-4 v38 driver, stage f12_sliver):  nice -n 19 python3 -I -B f12_sliver_v38.py      (children: same interpreter, -I -B, one at a time)"""
import os, sys, json, hashlib, math, subprocess, re
from fractions import Fraction as Fr
from decimal import Decimal, getcontext, ROUND_CEILING, ROUND_FLOOR
sys.dont_write_bytecode = True
getcontext().prec = 80
P2 = '/home/claude/s1spec/step4/p2'; OV = f'{P2}/out_v38'; F12 = f'{OV}/f12'; NEW = f'{F12}/new'; REP = f'{F12}/repro'; LOGS = f'{F12}/logs'
B = '/home/claude/s1spec/writeup/lemmaU/bulk'; NC = f'{B}/fix13_v37chain'; T2 = f'{NC}/t2_A'; RIG = f'{T2}/rigor'; BASE = f'{T2}/base_ibp'
T4 = f'{NC}/t4_cert'; T5 = f'{NC}/t5_zone3'
N1PY = f'{NC}/n1/n1_chart.py'; N1PIN = f'{NC}/n1/out/n1_chart_fin13v4c.json'; N1S19 = f'{OV}/n1_chart.json'   # v38: chain copy of n1_chart; stage-19 output must equal it
VC = '/home/claude/s1spec/step3nl_v4c'
RT = f'{VC}/runs/P/rate_table_fin13v4c_a0.00015.json'; RT125 = f'{VC}/runs/P/rate_table_fin13v4c_a0.000125.json'
FAP = f'{VC}/runs/P/modsups_fin13v4c/final_alpha_points_5e-8.json'
SOBP2 = f'{VC}/agent_sob/sob_constants_P2.json'
CONSTPY = f'{VC}/verifier/constants.py'
SCUB = f'{VC}/v4bcert/out/strip_cubic_v4c.json'; A2B = f'{VC}/v4bcert/out/e0_A2_budget_v4c.json'      # v38: D-3 strip cubic (H8), A-1 budget
ASM = ['design_N2ref', 'design_N2ref4', 'design_N20', 'lo_N2ref', 'hi_N2ref', 'probe_N2ref', 'design_pm3ibp']
AIPIN = f'{T2}/out/A_interface.json'; UPIN = f'{T4}/out/U_cert_v36_pass3.json'; UZPIN = f'{T4}/out/U_cert_for_zone3_v36.json'; ZPIN = f'{T5}/out/zone3_interface.json'
UP2 = f'{T4}/out/U_cert_v36_pass2.json'; MKUZ = f'{T4}/make_ucert_for_zone3.py'
CN = f'{B}/fix1/lemmaU_rev/constants_needed_v34.json'; GENPY = f'{B}/fix12_texts/tools/gen_numbers_A_v37.py'
X0LO = '2.2490776505'                     # <= -sigma_tip (s1chart.SIG_TIP ball: [2.24907765057, 2.24907765152]) and <= X0 = 2.2490776515 (asserted in the child)
XF = '16.249'                             # F2 sliver x in [16.249, 16.25]
sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()
PY = [sys.executable, '-I', '-B']; ME = os.path.abspath(__file__)

# =============================================================================================== certified sliver functions (arb)
def g_ders():
    import sympy as sp
    v = sp.Symbol('v'); h = v**4*(35 - 84*v + 70*v**2 - 20*v**3); g = 1/(1 - h)
    return [sp.lambdify(v, sp.diff(g, v, m), modules=[{}]) for m in range(4)]
def sliver_G_dom(vmax, ymax, nsub=64, der=None):
    """verbatim arithmetic of assembly_arb.sliver_G with the domain v in [0, vmax], y <= ymax as parameters"""
    from flint import arb
    der = der or g_ders(); sup = {m: arb(0) for m in range(4)}
    for c in range(nsub):
        ball = (vmax*c/nsub).union(vmax*(c + 1)/nsub)
        for m in range(4):
            val = abs(der[m](ball)).upper() if m else abs(der[0](ball) - 1).upper()
            val = val/6**m if m else val
            if val > sup[m]: sup[m] = val
    S2 = {0: {0: 1}, 1: {1: 1}, 2: {1: 1, 2: 1}, 3: {1: 1, 2: 3, 3: 1}}
    G = {0: sup[0]}
    for i in (1, 2, 3): G[i] = sum((arb(S2[i][m])*ymax**m*sup[m] for m in S2[i]), arb(0))
    return G, vmax, ymax, sup
def Dq(s):
    from flint import arb
    f = Fr(str(s)); return arb(f.numerator)/arb(f.denominator)
def doms(X0=None):
    """(vmax, ymax) of: old (assembly_arb; X0 = its module constant when given), new (extended, stated 2^-19), needed (e_ramp of n1_chart), F2 (x >= 16.249)"""
    from flint import arb
    X0 = Dq('2.2490776515') if X0 is None else X0; XL, WS = Dq(X0LO), arb(2)**-19
    e_need = Dq(repr(json.load(open(N1PIN))['ramp']['e']))
    return dict(old=((Dq('16.25') - (X0 + 14))/6, Dq('16.25') - X0),
                new=((Dq('16.25') - (XL + 14 - WS))/6, Dq('16.25') - XL),
                needed=((Dq('16.25') - (XL + 14 - e_need))/6, Dq('16.25') - XL),
                f2=((Dq('16.25') - Dq(XF))/6, None))
def same(a, b): return a.mid() == b.mid() and a.rad() == b.rad()

# =============================================================================================== children
def child_asm(mode, argv):
    """(A) assembly_rep.py as __main__ in this process, with sliver_G patched (mode new) or checked only (mode repro)"""
    import runpy
    os.chdir(T2); sys.path.insert(0, RIG)
    import arlib                                         # noqa  (inserts base_ibp paths; chdir restored)
    import assembly_arb as AS
    os.chdir(T2)
    ORIG = AS.sliver_G; der = g_ders()
    def patched(nsub=64):
        G0, vm0, ym0, sp0 = ORIG(nsub); D = doms(AS.X0)
        Gc, vmc, ymc, spc = sliver_G_dom(*D['old'], nsub, der)
        assert all(same(G0[i], Gc[i]) for i in range(4)) and same(vm0, vmc) and same(ym0, ymc), 'G harness differs from sliver_G'
        if mode == 'repro': return G0, vm0, ym0, sp0
        return sliver_G_dom(*D['new'], nsub, der)
    AS.sliver_G = patched
    sys.argv = ['rigor/assembly_rep.py'] + argv
    runpy.run_path(f'{RIG}/assembly_rep.py', run_name='__main__')
def child_export(asmdir, outdir):
    import importlib.util
    os.chdir(T2); sys.path.insert(0, RIG)
    spec = importlib.util.spec_from_file_location('export_t2', f'{RIG}/export_t2.py'); mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    if asmdir != 'pinned':
        for k, nm in (('asm_ref', 'design_N2ref'), ('asm_ref4', 'design_N2ref4'), ('asm_0', 'design_N20'), ('asm_lo', 'lo_N2ref'), ('asm_hi', 'hi_N2ref'),
                      ('asm_probe', 'probe_N2ref'), ('asm_pm3', 'design_pm3ibp')):
            assert mod.F[k].endswith(f'/asm_{nm}.json'); mod.F[k] = f'{asmdir}/asm_{nm}.json'
    mod.O = outdir; sys.argv = ['rigor/export_t2.py']; mod.main()
def child_exec(path, out, repl, argv, extra):
    """exec a pinned script with exact, asserted single-occurrence substitutions (scratch copy in memory)"""
    src = open(path).read()
    for a, b in repl:
        assert src.count(a) == 1, ('substitution target not unique', path, a); src = src.replace(a, b)
    os.chdir(os.path.dirname(path)); sys.path.insert(0, os.path.dirname(path)); sys.argv = argv
    g = dict(__name__='__main__', __file__=path, **extra); exec(compile(src, path, 'exec'), g)

if len(sys.argv) > 1 and sys.argv[1].startswith('--child-'):
    k = sys.argv[1]
    if k == '--child-asm': child_asm(sys.argv[2], sys.argv[3:])
    elif k == '--child-export': child_export(sys.argv[2], sys.argv[3])
    elif k == '--child-n1':
        si2 = None if sys.argv[2] == 'repro' else float.fromhex(sys.argv[2])
        child_exec(N1PY, sys.argv[3], [] if si2 is None else [("W2SUP = arb(UPd(FM['sI2']))", "W2SUP = arb(_F12_SI2)")], ['n1_chart.py', '--out', sys.argv[3]], dict(_F12_SI2=si2))
    elif k == '--child-zone3':
        child_exec(f'{T5}/zone3.py', sys.argv[4], [("json.dump(out, open(os.path.join(HERE, 'out', 'zone3_interface.json'), 'w'), indent=1, sort_keys=True)",
                                                    "json.dump(out, open(_F12_OUT, 'w'), indent=1, sort_keys=True)")],
                   ['zone3.py', '--ucert', sys.argv[2], '--ucert_full', sys.argv[3]], dict(_F12_OUT=sys.argv[4]))
    sys.exit(0)

# =============================================================================================== main stage
from flint import arb, ctx
ctx.prec = 200
for d in (F12, NEW, REP, LOGS): os.makedirs(d, exist_ok=True)
LOG = []
def say(*a):
    s = ' '.join(str(x) for x in a); LOG.append(s); print(s, flush=True)
def run(tag, args, cwd):
    with open(f'{LOGS}/{tag}.log', 'w') as lf:                 # child logs: not pinned (some print wall times)
        r = subprocess.run(PY + args, cwd=cwd, env=dict(os.environ), stdout=lf, stderr=subprocess.STDOUT)
    assert r.returncode == 0, (tag, open(f'{LOGS}/{tag}.log').read()[-3000:])
fup = lambda x: (lambda f: f if arb(f) >= x.upper() else math.nextafter(f, math.inf))(float(x.upper()))
flo = lambda x: (lambda f: f if arb(f) <= x.lower() else math.nextafter(f, -math.inf))(float(x.lower()))
UPd = lambda x: math.nextafter(float(x), math.inf)
J = lambda p: json.load(open(p))
def strip_args(x): return {k: v for k, v in x.items() if k != 'args'} if isinstance(x, dict) else x

# ---- closed trees: snapshot (hash of every file) before
TREES = [T2, T4, T5, f'{NC}/n1', f'{NC}/r3', f'{NC}/v36cfg.py', f'{NC}/runs/prod2', f'{NC}/flux/out', f'{B}/fix1/lemmaU_rev', f'{B}/fix12_texts/tools',
         f'{VC}/verifier', f'{VC}/agent_sob', f'{VC}/runs/P', f'{VC}/v4bcert/out', f'{VC}/v4cert/out']
def snapshot():
    out = {}
    for t in TREES:
        if os.path.isfile(t): out[t] = sha(t); continue
        for dp, dn, fn in os.walk(t):
            for f in fn: p = os.path.join(dp, f); out[p] = sha(p)
    return out
SNAP0 = snapshot()
say(f'closed-tree snapshot: {len(SNAP0)} files')
CHK = []
def check(name, ok, detail=None):
    CHK.append(dict(check=name, ok=bool(ok), detail=detail)); say(('PASS ' if ok else 'FAIL ') + name, '' if detail is None else json.dumps(detail, default=str)[:400])
# ---- v38: the verifier objects pinned are runs/P only; every input pin of its two rate tables must match the disk
for rtp in (RT, RT125):
    pins = J(rtp)['args']['inputs']; badp = {p: h for p, h in pins.items() if not os.path.exists(p) or sha(p) != h}
    check(f'v38 rate table {os.path.basename(rtp)}: all {len(pins)} input pins match the disk', not badp, badp)
    check(f'v38 rate table {os.path.basename(rtp)} pins strip_cubic_v4c.json (D-3, H8) at its current sha256', pins.get(SCUB) == sha(SCUB), pins.get(SCUB))
    check(f'v38 rate table {os.path.basename(rtp)}: TAG fin13v4c, V4C = 1, FLUX0/FLUX2 = 1.64e-20/8.0e-20',
          J(rtp)['args']['env'].get('TAG') == 'fin13v4c' and J(rtp)['args']['env'].get('V4C') == '1' and J(rtp)['args']['env'].get('FLUX0') == '1.64e-20' and J(rtp)['args']['env'].get('FLUX2') == '8.0e-20')
for rtp in (RT, RT125):
    BH = J(rtp)['box_hypotheses']['max_ratio']
    check(f'v38 rate table {os.path.basename(rtp)}: all {len(BH)} box-hypothesis max ratios <= 1 (incl. D-2a CF_*)', all(v <= 1 for v in BH.values()) and all(k in BH for k in ('CF_w_16_31', 'CF_w1_16_31', 'CF_alpha')), BH)
check('v38 modsups: all input pins match the disk (incl. both rate tables)', all(os.path.exists(p) and sha(p) == h for p, h in J(FAP)['args']['inputs'].items()) and all(p in J(FAP)['args']['inputs'] for p in (RT, RT125)))
_pd = []
for d in ('t2_A', 't4_cert', 't5_zone3'):
    for f in sorted(os.listdir(f'{NC}/{d}/out')):
        if os.path.isfile(f'{NC}/{d}/out/{f}') and sha(f'{NC}/{d}/out/{f}') != sha(f'{NC}/runs/prod2/{d}/{f}'): _pd.append(f'{d}/{f}')
check('v38 chain v37 in-place outputs (t2_A, t4_cert, t5_zone3 out/) == runs/prod2 byte for byte', not _pd, _pd)
check('v38 chain n1 (n1_chart_fin13v4c.json) == runs/prod2/n1 copy', sha(N1PIN) == sha(f'{NC}/runs/prod2/n1/n1_chart_fin13v4c.json'))
check('v38 stage 19 out_v38/n1_chart.json == chain n1_chart_fin13v4c.json (all keys but args)', strip_args(J(N1S19)) == strip_args(J(N1PIN)))

# ---- F1: G_i old / extended / needed
D = doms(); der = g_ders()
GO = sliver_G_dom(*D['old'], 64, der); GN = sliver_G_dom(*D['new'], 64, der); GD = sliver_G_dom(*D['needed'], 64, der)
assert (Dq(X0LO) - Dq('2.2490776515')).upper() < 0
import importlib.util as _iu
sys.path.insert(0, f'{P2}/t9h'); import s1chart                              # noqa
assert (Dq(X0LO) - (-s1chart.SIG_TIP)).upper() < 0, 'X0LO must be <= -sigma_tip'
F1 = dict(domains=dict(old='y in [14, 16.25 - 2.2490776515] (assembly_arb.sliver_G)', new=f'y in [14 - 2^-19, 16.25 - {X0LO}] (stated domain, used)',
                       needed=f'y in [14 - e_ramp, 16.25 - {X0LO}], e_ramp = n1_chart ramp.e (reported only)'),
          vmax={k: fup(D[k][0]) for k in ('old', 'new', 'needed')},
          G_old=[fup(GO[0][i]) for i in range(4)], G_new=[fup(GN[0][i]) for i in range(4)], G_needed=[fup(GD[0][i]) for i in range(4)],
          ratio_new_over_old=[fup(GN[0][i]/GO[0][i]) for i in range(4)], ratio_needed_over_old=[fup(GD[0][i]/GO[0][i]) for i in range(4)],
          label='certified (arb upper bounds; 64 cells in v, as sliver_G)')
for i in range(4): assert (GN[0][i] - GO[0][i]).lower() >= 0
say('F1 G_old', F1['G_old']); say('F1 G_new', F1['G_new'], 'ratio', F1['ratio_new_over_old'])

# ---- F2: sliver sup|w''|
SOB = J(SOBP2); assert SOB['W']['W_normalisation'].startswith('W = (1+x^2)^(-2) for x <= 30.25')
PSI1 = SOB['W']['Psi1']['R=60']; assert PSI1['from_10.25'] == 1.305 and PSI1['from_16.25'] == 1.305
assert SOB['W']['psi1_min_lower_bound_on_[6,80]'] > 0                          # psi' > 0: W decreasing on [6, 80]
cs = open(CONSTPY).read(); assert 'KW = 2.323e6; PSI1W = 1.305; C1W = 0.8477' in cs
KW, PSI1W = Dq('2.323e6'), Dq('1.305')
invW = (1 + Dq('16.25')**2)**2; assert (invW - KW).upper() < 0                  # 1/W(x) <= 1/W(16.25) on the sliver (W decreasing)
xb = Dq(XF).union(Dq('16.25')); assert (4*xb/(1 + xb**2) - PSI1W).upper() < 0   # psi' = 4x/(1+x^2) on the sliver: Psi1(x) <= max(psi' on sliver, Psi1 from 16.25)
vF = D['f2'][0]; supF = {0: arb(0), 1: arb(0)}
for c in range(64):
    ball = (vF*c/64).union(vF*(c + 1)/64)
    supF[0] = max(supF[0], abs(der[0](ball) - 1).upper()); supF[1] = max(supF[1], (abs(der[1](ball))/6).upper())
G0x, G1x = supF[0], supF[1]
RTJ = J(RT); FAPJ = J(FAP); FM = RTJ['far_max']
assert RTJ['alpha1'] == 0.00015 and FAPJ['1.5e-4']['sI2'] == FM['sI2'] and FAPJ['1.5e-4']['sI1'] == FM['sI1'] and FAPJ['1.5e-4']['M1'] == FM['M1_ramp']
F2rows = {}
for row in ('1.5e-4', '1.25e-4'):
    r = FAPJ[row]; s2, s1 = arb(UPd(r['sI2'])), arb(UPd(r['sI1']))
    sl = (1 + G0x)*s2 + G1x*s1; new = max(s2, sl)
    F2rows[row] = dict(sI2=r['sI2'], sI1=r['sI1'], sI2_sliver=fup(sl), sI2_prime=fup(new), rel_increase=fup(new/arb(r['sI2']) - 1),
                       M1_prime=fup(28*arb(fup(new))), M1_verifier=r['M1'])
SI2P = F2rows['1.5e-4']['sI2_prime']
assert SI2P > FM['sI2']
def _up4(x):                                                                    # verbatim rounding of assembly_rep.margins._up4 (v37)
    d = Decimal(repr(float(x))); q = Decimal(1).scaleb(d.adjusted() - 3)
    r = d.quantize(q, rounding=ROUND_CEILING); assert r >= Decimal(repr(float(x))); return str(r)
VM1A = _up4(FAPJ['1.5e-4']['M1'])                                               # v38: the M1 string (A) reads (0.4329 for fin13v4c)
asmlog = open(f'{T2}/out/asm_design_N2ref.log').read(); assert f'M14 verifier ramp sups (rounded up): M1 {VM1A} ' in asmlog, 'pinned (A) log does not show the M1 read'
M1P = 28*arb(SI2P); assert (M1P - Dq(VM1A)).upper() < 0                           # (A) reads M1 rounded up to 4 digits (assembly_rep.margins)
F2 = dict(sliver=f'x in [{XF}, 16.25]', lemma='(I2_W) of agent_sob/SOB_NOTES.md at every sliver point for v, v\' (v = zeta w\', the verifier variable); '
          '|w\'\'| <= sup|g| |v\'| + sup|g_x| |v|, g = 1/zeta', inv_W_16_25=fup(invW), K_W=2.323e6, psi_prime_sliver_max=fup(4*xb/(1 + xb**2)), Psi1=1.305,
          sup_g_minus_1=fup(G0x), sup_g_x=fup(G1x), rows=F2rows, sI2_prime_used=SI2P, sI2_prime_hex=float(SI2P).hex(),
          M1_prime=fup(M1P), M1_read_by_A=VM1A + ' (4 digits, rounded up, as assembly_rep.margins)', M1_prime_le_M1_read_by_A=True, label='certified')
say('F2', json.dumps(F2rows['1.5e-4']))

# ---- harness validation (repro) and recomputation (new)
run('n1_repro', [ME, '--child-n1', 'repro', f'{REP}/n1_chart.json'], P2)
run('n1_new', [ME, '--child-n1', float(SI2P).hex(), f'{NEW}/n1_chart.json'], P2)
n1p, n1r, n1n = J(N1PIN), J(f'{REP}/n1_chart.json'), J(f'{NEW}/n1_chart.json')
assert strip_args(n1p) == strip_args(n1r), 'n1_chart harness does not reproduce the pinned n1_chart.json'
assert n1n['verifier']['w2_sup'] == SI2P
say('n1 harness reproduces pinned; new sup_add', n1n['ramp']['sup_add'], 'pinned', n1p['ramp']['sup_add'])
ASMPIN = {nm: f'{T2}/out/asm_{nm}.json' for nm in ASM}
for nm in ASM:
    av = J(ASMPIN[nm])['args']['argv'][1:]
    i = av.index('--chart'); assert av[i + 1] == N1PIN; o = av.index('--out'); assert av[o + 1] == f'out/asm_{nm}.json'
    rep = list(av); rep[o + 1] = f'{REP}/asm_{nm}.json'
    new = list(av); new[i + 1] = f'{NEW}/n1_chart.json'; new[o + 1] = f'{NEW}/asm_{nm}.json'
    run(f'asm_{nm}_repro', [ME, '--child-asm', 'repro'] + rep, T2)
    assert strip_args(J(ASMPIN[nm])) == strip_args(J(f'{REP}/asm_{nm}.json')), ('asm harness does not reproduce', nm)
    run(f'asm_{nm}_new', [ME, '--child-asm', 'new'] + new, T2)
    say(f'asm {nm}: harness reproduces pinned; new run done')
run('export_repro', [ME, '--child-export', 'pinned', REP], T2)
assert strip_args(J(AIPIN)) == strip_args(J(f'{REP}/A_interface.json')), 'export harness does not reproduce A_interface'
run('export_new', [ME, '--child-export', NEW, NEW], T2)
ZREL = ZPIN                                                                    # v38: U pass 3 of chain v37 reads the final zone3_interface (absolute path)
assert J(UPIN)['args']['argv'] == ['assemble_U.py', '--A', AIPIN, '--zone3', ZREL, '--out', 'U_cert_v36_pass3.json']
run('U_repro', [f'{T4}/assemble_U.py', '--A', AIPIN, '--zone3', ZREL, '--out', f'{REP}/U_cert.json'], T4)
assert strip_args(J(UPIN)) == strip_args(J(f'{REP}/U_cert.json')), 'assemble_U harness does not reproduce U_cert_v36_pass3'
check('v38 U fixed point of chain v37: U pass 3 == pass 2 (all keys but args)', strip_args(J(UPIN)) == strip_args(J(UP2)))
run('U_new', [f'{T4}/assemble_U.py', '--A', f'{NEW}/A_interface.json', '--zone3', ZREL, '--out', f'{NEW}/U_cert.json'], T4)
Un = J(f'{NEW}/U_cert.json')
run('Uz_repro', [MKUZ, UP2, f'{REP}/U_cert_for_zone3.json'], T4)                 # v38: the chain's make_ucert_for_zone3.py, validated on pass 2
assert strip_args(J(UZPIN)) == strip_args(J(f'{REP}/U_cert_for_zone3.json')), 'make_ucert_for_zone3 harness does not reproduce U_cert_for_zone3_v36'
run('Uz_new', [MKUZ, f'{NEW}/U_cert.json', f'{NEW}/U_cert_for_zone3.json'], T4)
assert J(ZPIN)['args']['argv'] == ['zone3.py', '--ucert', '../t4_cert/out/U_cert_for_zone3_v36.json', '--ucert_full', UP2]
run('zone3_repro', [ME, '--child-zone3', '../t4_cert/out/U_cert_for_zone3_v36.json', UP2, f'{REP}/zone3_interface.json'], T5)
assert strip_args(J(ZPIN)) == strip_args(J(f'{REP}/zone3_interface.json')), 'zone3 harness does not reproduce zone3_interface'
run('zone3_new', [ME, '--child-zone3', f'{NEW}/U_cert_for_zone3.json', f'{NEW}/U_cert.json', f'{NEW}/zone3_interface.json'], T5)
Zp, Zn = J(ZPIN), J(f'{NEW}/zone3_interface.json')
ZKEYS = ['DeltaR', 'DeltaR_b', 'n_III', 'box_III_380', 'box_III_700', 'tau_III', 'DeltaR_detail', 'n_III_pieces']
assert all(Zp[k] == Zn[k] for k in ZKEYS), 'zone3 quantities consumed by U changed'
run('U_new_pass2', [f'{T4}/assemble_U.py', '--A', f'{NEW}/A_interface.json', '--zone3', f'{NEW}/zone3_interface.json', '--out', f'{NEW}/U_cert_pass2.json'], T4)
assert strip_args(J(f'{NEW}/U_cert_pass2.json')) == strip_args(Un), 'U pass 2 differs from pass 1'
say('harnesses reproduce n1, 7 asm, A_interface, U_cert_v36_pass3, U_cert_for_zone3_v36, zone3; zone3 inputs of U unchanged; U pass 2 == pass 1 (fixed point)')

# =============================================================================================== comparisons
def interval(v):
    if isinstance(v, dict) and 'arb' in v:
        a = v['arb']
        if isinstance(a, list): m = Fr(int(a[0]))*(Fr(2)**a[1]); r = Fr(a[2]); return m - r, m + r
        s = a.strip('[]')
        if '+/-' in s: m, r = s.split('+/-'); m = Fr(Decimal(m.strip())); r = Fr(Decimal(r.strip())); return m - r, m + r
        x = Fr(Decimal(s)); return x, x
    if isinstance(v, bool): raise TypeError
    if isinstance(v, (int, float)): x = Fr(v); return x, x
    if isinstance(v, str): x = Fr(Decimal(v)); return x, x
    raise TypeError(repr(v)[:80])
def resolve(d, path):
    parts = path.split('.') if path else []; i = 0; cur = d
    while i < len(parts):
        if isinstance(cur, list): cur = cur[int(parts[i])]; i += 1; continue
        for j in range(len(parts), i, -1):
            k = '.'.join(parts[i:j])
            if k in cur: cur = cur[k]; i = j; break
        else: raise KeyError(path)
    return cur
def leaves(x, pre=''):
    if isinstance(x, dict):
        if 'arb' in x and isinstance(x['arb'], (list, str)):
            yield pre, x; return
        for k, v in x.items():
            if k == 'args': continue
            yield from leaves(v, f'{pre}.{k}' if pre else str(k))
    elif isinstance(x, list):
        for i, v in enumerate(x): yield from leaves(v, f'{pre}.{i}')
    elif isinstance(x, (int, float)) and not isinstance(x, bool): yield pre, x
    else: yield pre, x
def rnd(x, n, up):
    x = Decimal(x.numerator)/Decimal(x.denominator)
    if x == 0: return Decimal(0)
    q = Decimal(1).scaleb(x.adjusted() - n + 1); return x.quantize(q, rounding=ROUND_CEILING if up else ROUND_FLOOR)
def diff(old, new, name):
    """every leaf: changed numbers (relative change of the bound endpoints), non-numeric changes; and rounding flips at 3..6 digits"""
    lo_ = dict(leaves(old)); ln = dict(leaves(new)); ch = []; other = []; flips = []
    assert set(lo_) == set(ln), (name, sorted(set(lo_) ^ set(ln))[:10])
    for k in lo_:
        a, b = lo_[k], ln[k]
        if a == b: continue
        try: (al, ah), (bl, bh) = interval(a), interval(b)
        except (TypeError, ValueError, ArithmeticError): other.append(k); continue
        ref = max(abs(ah), abs(al), Fr(1, 10**300))
        ch.append(dict(key=k, old=[float(al), float(ah)], new=[float(bl), float(bh)], rel_hi=float((bh - ah)/ref), rel_lo=float((bl - al)/ref)))
        for n in (3, 4, 5, 6):
            if (bh > ah and rnd(bh, n, True) != rnd(ah, n, True)) or (bl < al and rnd(bl, n, False) != rnd(al, n, False)):
                flips.append(dict(key=k, digits=n, old_up=str(rnd(ah, n, True)), new_up=str(rnd(bh, n, True)), old_dn=str(rnd(al, n, False)), new_dn=str(rnd(bl, n, False))))
    worse = [c for c in ch if c['rel_hi'] > 0 or c['rel_lo'] < 0]
    return dict(file=name, leaves=len(lo_), changed=len(ch), max_rel_change=max([max(abs(c['rel_hi']), abs(c['rel_lo'])) for c in ch], default=0.0),
                nonnumeric_changed=other, rounding_flips=flips, changed_rows=sorted(ch, key=lambda c: -max(abs(c['rel_hi']), abs(c['rel_lo'])))[:60], worse_count=len(worse))
lo = lambda v: interval(v)[0]; hi = lambda v: interval(v)[1]
DIFFS = {}
# ---- n1_chart (chi bounds)
DIFFS['n1_chart'] = diff(n1p, n1n, 'n1_chart')
check('n1 chart ok condition (ramp e < 2^-19, sup chi^add_0 <= 2e-7)', n1n['ramp']['e'] < 2**-19 and n1n['ramp']['sup_add'][0] <= 2e-7,
      dict(e=n1n['ramp']['e'], sup_add0=n1n['ramp']['sup_add'][0]))
check('n1 chart: beta22, beta_c, beta_c3 < 1e-6 (assembly_rep assert)', max(n1n['ramp']['beta22'], n1n['core']['beta22'], n1n['core']['beta3']) < 1e-6)
CHI = {g: {k: dict(pinned=n1p[g][k], new=n1n[g][k]) for k in n1p[g] if n1p[g][k] != n1n[g][k]} for g in ('ramp', 'core')}
# ---- asm runs
ASMRES = {}
for nm in ASM:
    P, Nw = J(ASMPIN[nm]), J(f'{NEW}/asm_{nm}.json'); DIFFS[f'asm_{nm}'] = diff(P, Nw, f'asm_{nm}'); ASMRES[nm] = {}
    for st in ('tight', 'v31'):
        p, n = P['res'][st], Nw['res'][st]
        fp = sorted(c[0] for c in p['claims_failing']); fn = sorted(c[0] for c in n['claims_failing'])
        check(f'asm {nm} {st}: no new failing claim', set(fn) <= set(fp), dict(pinned_failing=fp, new_failing=fn))
        mrow = {}
        for k, v in p['margins'].items():
            a, b = lo(v), lo(n['margins'][k]); mrow[k] = dict(pinned=float(a), new=float(b), same_4digits=str(rnd(a, 4, False)) == str(rnd(b, 4, False)))
            if a > 1: check(f'asm {nm} {st}: margin {k} > 1', b > 1, mrow[k])
        fl = {}
        for bq, v in p['flux']['beta'].items():
            for mk in ('margin_FLUX0', 'margin_FLUX2'):
                a, b = lo(v[mk]), lo(n['flux']['beta'][bq][mk]); fl[f'beta{bq}.{mk}'] = dict(pinned=float(a), new=float(b), same_4digits=str(rnd(a, 4, False)) == str(rnd(b, 4, False)))
                if a > 1: check(f'asm {nm} {st}: flux {bq} {mk} > 1', b > 1, fl[f'beta{bq}.{mk}'])
        E = [dict(pinned=float(hi(p['E'][j])), new=float(hi(n['E'][j])), rel=float((hi(n['E'][j]) - hi(p['E'][j]))/hi(p['E'][j]))) for j in range(3)]
        ASMRES[nm][st] = dict(E=E, margins=mrow, flux=fl, S0=dict(pinned=float(hi(p['S0'])), new=float(hi(n['S0']))), S1=dict(pinned=float(hi(p['S1'])), new=float(hi(n['S1']))),
                              S0w=dict(pinned=float(hi(p['verifier_sups']['S0w'])), new=float(hi(n['verifier_sups']['S0w']))))
# ---- A_interface
Ap, An = J(AIPIN), J(f'{NEW}/A_interface.json'); DIFFS['A_interface'] = diff(Ap, An, 'A_interface')
for st in ('tight', 'v31'): check(f'A_interface {st}: rule_reproduces_EA_ref', An['res'][st]['rule_reproduces_EA_ref'] is True)
for st in ('tight', 'v31'):                                                   # M3 (n35_asserts)
    for k in ('ramp M1', 'ramp z1', 'ramp M0'):
        m = An['res'][st]['core']['margins'][k]; assert m['bound'] == 'lo'
        check(f'M3 {st} {k} > 1', lo(m) > 1, dict(new=float(lo(m)), pinned=float(lo(Ap['res'][st]['core']['margins'][k]))))
# ---- U
Up = J(UPIN); DIFFS['U_cert'] = diff(Up, Un, 'U_cert')
ss = Un['supersolution']
for k in ('K1_ok', 'K2_ok', 'K3_ok', 'N2_le_4N2ref'): check(f'U supersolution {k}', ss[k] is True)
check('M6 Nbar^2/(4 N2_ref) < 1', max(ss['N2_over_4N2ref']) < 1 - 1e-9, ss['N2_over_4N2ref'])
check('U M1z hypothesis checks', ss['M1z_hyp_check']['zones_I_II']['ok'] and ss['M1z_hyp_check']['zone_III']['ok'])
check('U-fine requirement of task 5', Un['ufine']['ok'] is True, Un['ufine']['Bhat'])
UM = {}
for k in ('late_M1', 'late_M0', 'late_M0z'):
    a, b = lo(Up['margins'][k]), lo(Un['margins'][k]); UM[k] = dict(pinned=float(a), new=float(b), same_4digits=str(rnd(a, 4, False)) == str(rnd(b, 4, False)))
    check(f'U margin {k} > 1', b > 1, UM[k])
for k in ('present_M1_min', 'present_M0_min', 'present_M0z_min'):
    a, b = Fr(Up['margins'][k]['value']), Fr(Un['margins'][k]['value']); UM[k] = dict(pinned=float(a), new=float(b), same_4digits=str(rnd(a, 4, False)) == str(rnd(b, 4, False)))
    check(f'U margin {k} > 1', b > 1, UM[k])
UM['loop_gain_total_max'] = dict(pinned=max(g['total'] for g in Up['loop_gain']), new=max(g['total'] for g in Un['loop_gain']))
check('U loop gain < 1', UM['loop_gain_total_max']['new'] < 1, UM['loop_gain_total_max'])
# ---- K6 (check_keys_v34.py side condition) and the other side conditions, on the new files
A4n = J(f'{NEW}/asm_design_N2ref4.json'); T = 'res.tight.'; A4p = J(ASMPIN['design_N2ref4'])
k6n = lo(resolve(Un, 'supersolution.EA.0'))/hi(resolve(An, T + 'IC.E_tau1.0')); k6p = lo(resolve(Up, 'supersolution.EA.0'))/hi(resolve(Ap, T + 'IC.E_tau1.0'))
f4 = lambda x: Fr(rnd(Fr(x), 4, False))                                        # v38: quoted = pinned v37 value rounded DOWN to 4 digits
check(f'K6 margin >= {float(f4(k6p))} (pinned v37 K6 rounded down)', k6n >= f4(k6p), dict(pinned=float(k6p), new=float(k6n)))
check('K6 at every level', all(hi(resolve(An, T + f'IC.E_tau1.{i}')) < lo(resolve(Un, f'supersolution.EA.{i}')) for i in range(3)))
check('EA_ref4 dominates Ehat', all(lo(resolve(An, T + f'EA_ref4.{i}')) >= hi(resolve(Un, f'supersolution.EA.{i}')) for i in range(3)))
check('A4 E == A EA_ref4', all(hi(resolve(A4n, T + f'E.{i}')) == hi(resolve(An, T + f'EA_ref4.{i}')) or abs(float(hi(resolve(A4n, T + f'E.{i}')))/float(hi(resolve(An, T + f'EA_ref4.{i}'))) - 1) <= 1e-12 for i in range(3)))
check('Nbar^2 <= 4 N2_ref', all(hi(resolve(Un, f'supersolution.N2.{i}')) <= 4*Fr(Decimal(An['N2_ref']['values'][i])) for i in range(3)))
S1n = hi(resolve(A4n, T + 'S1'))      # v38: the v3.3 quote 0.1475 (lemmaU_v34, S1v <= 0.10429) is stale since chain v36; lemmaU_v37/v38.tex quote S1v <= 0.10435, sqrt2 S1v <= 0.1476 < 0.52
check('(A) core |z1|: S1v <= 0.10435 and sqrt2 S1v <= 0.1476 < 0.52 (lemmaU_v37.tex l.354, 476)', S1n <= Fr(10435, 100000) and 2*S1n*S1n <= Fr(1476, 10000)**2,
      dict(S1v_new=float(S1n), S1v_pinned=float(hi(resolve(A4p, T + 'S1')))))
check('D^c <= D^src', all(hi(resolve(A4n, T + f'info.{i}.D_window_const')) <= lo(resolve(An, T + f'Dsrc.{i}')) for i in range(2)))
check('D^s A/A4', all(hi(resolve(A4n, T + f'info.{i}.D_window_slope')) <= hi(resolve(An, T + f'B0.D_window_slope.{i}'))*(1 + Fr(1, 10**12)) for i in range(2)))
b2p = Fr(153, 10000)/(hi(resolve(A4p, T + 'S0')) + hi(resolve(Ap, T + 'B0.sup_core.0')))
b2n = Fr(153, 10000)/(hi(resolve(A4n, T + 'S0')) + hi(resolve(An, T + 'B0.sup_core.0'))); check(f"B2' margin >= {float(f4(b2p))} (pinned v37 rounded down)", b2n >= f4(b2p), dict(pinned=float(b2p), new=float(b2n)))
# ---- quoted numbers of lemmaU_v34.tex (constants_needed_v34.json 'numbers'), evaluated on the new files
CNJ = J(CN); FILES = {'A': An, 'A4': A4n, 'U': Un, 'Z': Zn}
for k, f in CNJ['files'].items():
    if k not in FILES: FILES[k] = J(f['path'])
QN = []
for q, f, k, d in CNJ['numbers']:
    l_, h_ = interval(resolve(FILES[f], k)); Q = Fr(Decimal(q)); ok = Q >= h_ if d == 'up' else Q <= l_
    if not ok: QN.append(dict(quoted=q, file=f, key=k, dir=d, new=[float(l_), float(h_)]))
QN_INFO = dict(note='v38: information only. lemmaU_v34.tex is superseded (texts v37/v38); its quotes are v35 values and are not expected to bound chain v37',
               n_numbers=len(CNJ['numbers']), n_not_outward_on_new=len(QN), rows=QN[:40])
QNP = []
FILESP = {'A': Ap, 'A4': A4p, 'U': Up, 'Z': J(ZPIN)}
for k, f in CNJ['files'].items():
    if k not in FILESP: FILESP[k] = J(f['path'])
for q, f, k, d in CNJ['numbers']:
    l_, h_ = interval(resolve(FILESP[f], k)); Q = Fr(Decimal(q)); ok = Q >= h_ if d == 'up' else Q <= l_
    if not ok: QNP.append(dict(quoted=q, file=f, key=k, dir=d))
QN_INFO['n_not_outward_on_pinned_v37'] = len(QNP)
QN_INFO['same_set_pinned_vs_new'] = sorted((r['file'], r['key']) for r in QNP) == sorted((r['file'], r['key']) for r in QN)
check('lemmaU_v34.tex quotes [information only]: F1/F2 change no outcome (same failing set on pinned v37 and on new)', QN_INFO['same_set_pinned_vs_new'], QN_INFO)
QN = []
# ---- macros of numbers_A_v37.tex (fix12_texts/tools/gen_numbers_A_v37.py, every source repointed to chain v37 / fin13v4c) on pinned and new files
def macros(ovr):
    src = open(GENPY).read()
    for a, b in (("N = B + 'fix6_v36/'", "N = B + 'fix13_v37chain/'"), ("V4B = '/home/claude/s1spec/step3nl_v4b/runs/P/'", "V4B = '/home/claude/s1spec/step3nl_v4c/runs/P/'"),
                 ("'RT15': V4B + 'rate_table_fin13v4b_a0.00015.json'", "'RT15': V4B + 'rate_table_fin13v4c_a0.00015.json'"),
                 ("'RT125': V4B + 'rate_table_fin13v4b_a0.000125.json'", "'RT125': V4B + 'rate_table_fin13v4c_a0.000125.json'"),
                 ("'MS': V4B + 'modsups_fin13v4b/final_alpha_points_5e-8.json'", "'MS': V4B + 'modsups_fin13v4c/final_alpha_points_5e-8.json'"),
                 ("'FX': B + 'fix8_flux2_impl/out/flux2_cert.json'", "'FX': N + 'flux/out/flux2_cert.json'"),
                 ("'FXN': B + 'fix8_flux2_impl/texts/numbers_flux2.json'", "'FXN': N + 'flux/texts/numbers_flux2.json'"),
                 ("def J(k): return json.load(open(SRC[k]))", "SRC.update(_OVR)\ndef J(k): return json.load(open(SRC[k]))")):
        assert src.count(a) == 1, a; src = src.replace(a, b)
    g = dict(__name__='gen_v37_f12', _OVR=ovr); exec(compile(src, GENPY, 'exec'), g)
    assert all(os.path.exists(p) and ('fix6_v36' not in p and 'step3nl_v4b' not in p and 'fix8_flux2_impl' not in p) for k, p in g['SRC'].items() if k != 'G1'), g['SRC']
    return g['M'], dict(g['SRC'])
MP, SRCP = macros(dict(A=AIPIN, A1=ASMPIN['design_N2ref'], A4=ASMPIN['design_N2ref4'], U=UPIN, Z=ZPIN))
MN, SRCN = macros(dict(A=f'{NEW}/A_interface.json', A1=f'{NEW}/asm_design_N2ref.json', A4=f'{NEW}/asm_design_N2ref4.json', U=f'{NEW}/U_cert.json', Z=f'{NEW}/zone3_interface.json'))
MD = {k: dict(pinned=MP[k], new=MN[k]) for k in MP if MP[k] != MN[k]}
check(f'numbers_A_v37 macros ({len(MP)}, gen_numbers_A_v37 repointed to chain v37) unchanged on the new files', not MD, MD)
# ---- zone3 H6
H6p, H6n = Zp['cone_family']['H6']['at_Bhat_ref']['min_margin'], Zn['cone_family']['H6']['at_Bhat_ref']['min_margin']
DIFFS['zone3_interface'] = diff(Zp, Zn, 'zone3_interface')
check(f'zone3 H6 min margin >= {float(f4(H6p))} (pinned v37 rounded down)', H6n >= f4(H6p), dict(pinned=H6p, new=H6n))
check('zone3 acceptance unchanged apart from H6 value', {k: v for k, v in Zn['acceptance'].items() if k != 'H6_min_margin_at_Bhat_ref'} == {k: v for k, v in Zp['acceptance'].items() if k != 'H6_min_margin_at_Bhat_ref'})
# ---- named claims of the review (thinnest margins): ramp M1, K6, H6; E-hat; flux; S0v
Rp, Rn = Ap['res']['tight'], An['res']['tight']
NAMED = dict(
    ramp_M1=dict(pinned=float(lo(Rp['core']['margins']['ramp M1'])), new=float(lo(Rn['core']['margins']['ramp M1'])), quoted='x%s' % rnd(lo(Rp['core']['margins']['ramp M1']), 4, False)),
    ramp_z1=dict(pinned=float(lo(Rp['core']['margins']['ramp z1'])), new=float(lo(Rn['core']['margins']['ramp z1']))),
    ramp_M0=dict(pinned=float(lo(Rp['core']['margins']['ramp M0'])), new=float(lo(Rn['core']['margins']['ramp M0']))),
    S0v=dict(pinned=float(lo(Rp['core']['margins']['S0v'])), new=float(lo(Rn['core']['margins']['S0v']))),
    S0v_vs_0p0153=dict(pinned=float(lo(Rp['core']['margins']['S0v_vs_0.0153'])), new=float(lo(Rn['core']['margins']['S0v_vs_0.0153']))),
    flux_beta2=dict(pinned=[float(lo(Rp['core']['flux_margin']['2'][k])) for k in ('FLUX0_over_K0', 'FLUX2_over_K2')],
                    new=[float(lo(Rn['core']['flux_margin']['2'][k])) for k in ('FLUX0_over_K0', 'FLUX2_over_K2')]),
    flux_beta2_4N2ref=dict(pinned=[ASMRES['design_N2ref4']['tight']['flux'][f'beta2.{k}']['pinned'] for k in ('margin_FLUX0', 'margin_FLUX2')],
                           new=[ASMRES['design_N2ref4']['tight']['flux'][f'beta2.{k}']['new'] for k in ('margin_FLUX0', 'margin_FLUX2')], quoted='pinned v37 rounded down'),
    K6=dict(pinned=float(k6p), new=float(k6n), quoted='x%s' % rnd(k6p, 4, False)), H6=dict(pinned=H6p, new=H6n, quoted='x%s' % rnd(Fr(H6p), 4, False)),
    EA_ref=dict(pinned=[float(hi(x)) for x in Rp['EA_ref']], new=[float(hi(x)) for x in Rn['EA_ref']]),
    EA_ref4=dict(pinned=[float(hi(x)) for x in Rp['EA_ref4']], new=[float(hi(x)) for x in Rn['EA_ref4']]),
    Ehat_U=dict(pinned=[float(hi(x)) for x in Up['supersolution']['EA']], new=[float(hi(x)) for x in Un['supersolution']['EA']]),
    Bhat_U=dict(pinned=Up['supersolution']['Bhat'], new=Un['supersolution']['Bhat']), Bhat_ufine=dict(pinned=Up['ufine']['Bhat'], new=Un['ufine']['Bhat']))
for k in ('ramp_M1', 'ramp_z1', 'ramp_M0', 'S0v', 'S0v_vs_0p0153'): check(f'named margin {k} > 1', NAMED[k]['new'] > 1, NAMED[k])
check('ramp M1: pinned v37 value rounded down to 4 digits still a lower bound', Fr(NAMED['ramp_M1']['new']) >= f4(lo(Rp['core']['margins']['ramp M1'])), NAMED['ramp_M1'])
check('flux (4 N2_ref, beta 2): pinned v37 values rounded down to 4 digits still lower bounds',
      all(Fr(NAMED['flux_beta2_4N2ref']['new'][i]) >= f4(Fr(NAMED['flux_beta2_4N2ref']['pinned'][i])) for i in (0, 1)), NAMED['flux_beta2_4N2ref'])
# ---- pinned NUMBERS exceeded (not only thresholds): every leaf that moved in the unfavourable direction, summarised per file
EXCEEDED = {k: dict(changed=v['changed'], worse=v['worse_count'], max_rel_change=v['max_rel_change'], rounding_flips=v['rounding_flips'],
                    nonnumeric_changed=v['nonnumeric_changed'][:20]) for k, v in DIFFS.items()}
FLIPS = [dict(file=k, **f) for k, v in DIFFS.items() for f in v['rounding_flips']]

# ---- closed trees unchanged
SNAP1 = snapshot()
chg = sorted(p for p in set(SNAP0) | set(SNAP1) if SNAP0.get(p) != SNAP1.get(p))
check('closed trees unchanged by this stage (every file hash, no new files)', not chg, chg[:20])
pyc = [os.path.join(dp, f) for t in TREES if os.path.isdir(t) for dp, dn, fn in os.walk(t) for f in fn if f.endswith('.pyc') and os.path.join(dp, f) not in SNAP0]
check('no new .pyc in closed trees', not pyc, pyc[:10])

OUTS = [f'{REP}/n1_chart.json', f'{NEW}/n1_chart.json'] + [f'{d}/asm_{nm}.json' for nm in ASM for d in (REP, NEW)] + \
       [f'{REP}/A_interface.json', f'{NEW}/A_interface.json', f'{REP}/U_cert.json', f'{NEW}/U_cert.json', f'{NEW}/U_cert_pass2.json', f'{NEW}/U_cert_for_zone3.json',
        f'{REP}/zone3_interface.json', f'{NEW}/zone3_interface.json', f'{REP}/U_cert_for_zone3.json']
INPUTS = sorted(set([ME, N1PY, N1PIN, N1S19, RT, RT125, FAP, SOBP2, CONSTPY, SCUB, A2B, f'{P2}/t9h/s1chart.py', AIPIN, UPIN, UP2, UZPIN, ZPIN, CN, GENPY, MKUZ, f'{NC}/v36cfg.py',
                     f'{RIG}/assembly_rep.py', f'{RIG}/arlib.py', f'{RIG}/export_t2.py', f'{BASE}/rigor/assembly_arb.py', f'{T4}/assemble_U.py', f'{T4}/t4lib.py',
                     f'{T5}/zone3.py', f'{T5}/z3lib.py'] + list(ASMPIN.values()) +
                    [p for o in OUTS if 'args' in J(o) for p in (J(o)['args'].get('inputs') or {}) if isinstance(p, str) and p.startswith('/') and os.path.exists(p)]))
ALL_OK = all(c['ok'] for c in CHK)
out = dict(item='v3.3 review items F1, F2 (t8_review/CHECK_N1_N5.md): standalone pinned check (no pinned output changed); v38: on verifier fin13v4c and chain v37', label='certified (arb, outward) unless stated',
           F1=F1, F2=F2, chart_changes=CHI, asm=ASMRES, named=NAMED, U_margins=UM, quoted_numbers_failing=QN, quoted_numbers_lemmaU_v34_info=QN_INFO, macro_changes=MD, macro_sources=dict(pinned=SRCP, new=SRCN), M1_read_by_A=VM1A, rounding_flips=FLIPS,
           exceeded_summary=EXCEEDED, diffs=DIFFS, checks=CHK, all_ok=ALL_OK, outputs={o: sha(o) for o in OUTS},
           harness=dict(n1='exec of n1_chart.py with the single line W2SUP = arb(UPd(FM[\'sI2\'])) -> arb(sI2\'); repro without substitution == pinned',
                        asm='runpy of assembly_rep.py as __main__ with assembly_arb.sliver_G replaced by the extended-domain G (old G recomputed and asserted '
                            'bitwise equal to sliver_G); argv = pinned argv with --chart / --out redirected; repro == pinned (all keys but args)',
                        export='export_t2.main with F[asm_*] and O redirected; repro == pinned A_interface', U='assemble_U.py --A <new> --out <abs f12 path>; repro == U_cert_v35',
                        zone3='exec of zone3.py with the output path replaced; repro == pinned zone3_interface; U-consumed keys unchanged, so U pass 2 == pass 1'),
           args=dict(argv=sys.argv, cwd=os.getcwd(), env={k: os.environ.get(k) for k in ('PYTHONHASHSEED', 'PYTHONDONTWRITEBYTECODE', 'OMP_NUM_THREADS', 'S1TT', 'S1PROF', 'FORM_MULT', 'LC_ALL', 'PATH', 'V36_VER', 'V36_R3', 'V36_R2')},
                     python=sys.version, inputs={p: sha(p) for p in INPUTS}))
tmp = f'{F12}/f12_sliver.json.tmp'; json.dump(out, open(tmp, 'w'), indent=1, default=str); os.replace(tmp, f'{F12}/f12_sliver.json')
say('F12 all checks ok:', ALL_OK, f'({sum(c["ok"] for c in CHK)}/{len(CHK)})')
assert ALL_OK, [c for c in CHK if not c['ok']]
