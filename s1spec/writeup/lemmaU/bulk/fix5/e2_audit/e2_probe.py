"""E2 audit probe: certified cell bounds (bb.bound_cell of the fin13v4 verifier; rhs_iv replaced by vp/rhs_iv.py = fin13v4 + switchable E2 terms,
see mk_vp.py) for the coded E2 block and for the corrected one.  [cert given the constants]: U is a rigorous upper bound of the rate over the
cell for the stated rhs.  Also: worst-case sup(|w|,|w'|) on x <= 8 and x <= 12 over all shapes with Lyap <= peak (box of the E2 coefficients).
usage: python3 -B e2_probe.py [tol]      env as verifier/run_v4.sh"""
import sys, os, json, time, hashlib
D = os.path.dirname(os.path.abspath(__file__)); V = '/home/claude/s1spec/step3nl_v4/verifier'
sys.path.insert(0, V); sys.path.insert(0, D + '/vp')
import numpy as np
import constants as K, rhs_iv as R, bb
assert R.__file__.startswith(D + '/vp') and hasattr(R, 'E2X') and K.__file__.startswith(V) and bb.__file__.startswith(V)
from ival import P, cc, sqrt, mx
TOL = float(sys.argv[1]) if len(sys.argv) > 1 else 5e-4
RT = json.load(open('/home/claude/s1spec/step3nl_v4/runs/P/rate_table_fin13v4_a0.00015.json'))
tl = RT['tail']; s0 = RT['steps'][0]; a_s, L_s = tl['a_star'], tl['L_star']
_c = R._c
# constants: e2_consts.json [cert] rounded up by hand below (asserted against the json); v4b values [read, fin13v4b in progress] rounded up
EC = json.load(open(D + '/e2_consts.json')); M = EC['d_modulation']
CST = dict(dab=0.013114, n3=7.69699, nh2=13.1085, ngab=8.02860, nf12=2.06197, ndf=12.34096, srw=4.4292e-3,
           cw0=1.00114,      # 1 + CA0 + CB0   (v4b: 7.1463e-4, 4.2378e-4)
           cw1=33.0322,      # CA1 + CB1       (v4b: 32.000, 1.032186)
           r2=32.4022,       # RA2 + RB2       (v4b: 26.90245, 5.49971)
           xwc=5.71e-6)      # ||chi' w||_rho <= (35/48) (int_strip rho)^(1/2) sup_strip|w|, int_strip rho <= 3 rho(31.24) = 6.1283e-11 [cert, crude]
assert CST['dab'] >= M['dab_lead_minus_Q2lead'] + M['dlead_ab'] and CST['n3'] >= M['n3_c4'] and CST['nh2'] >= M['nh2_c4'] and CST['ngab'] >= M['NGAB']
assert CST['nf12'] >= M['NF12'] and CST['ndf'] >= M['NDF_proj12'] and CST['srw'] >= EC['e_strip']['SRW_sup_sqrt_rho_over_W']
assert CST['xwc'] >= 35 / 48 * EC['e_strip']['IR_crude_3rho'] ** 0.5 * (1 + 1e-12)
AL = 2.5e-4                                                    # alpha level of the forcing factors (cells have alpha <= 1.6e-4)
n3e = CST['n3'] + 2 * AL * CST['nh2']; u0p = AL * (CST['ngab'] + AL * CST['r2'])
FE2_CS = (1 + M['BETAF'] * n3e + M['PHIAB'] * (CST['nf12'] + u0p)) * (1 + 1e-9)
FE2_IBP = (1 + M['BETAF'] * n3e + M['PHIAB'] * u0p) * (1 + 1e-9); FE2A = M['PHIAB'] * CST['ndf'] * (1 + 1e-9)
def mk(keys, **kw):
    d = {k: _c(CST[k]) for k in keys}; d.update({k: _c(v) for k, v in kw.items()}); return d
base = ('dab', 'n3', 'nh2', 'ngab', 'srw', 'cw0', 'cw1', 'r2', 'xwc')
VARS = {   # name: (r143, E2X)
 'coded_fin13v4':        (14.34, None),
 'r143_up_only':         (14.3431, None),
 'E2fix_CS':             (14.3431, lambda: mk(base + ('nf12',), fe2=FE2_CS)),
 'E2fix_IBP':            (14.3431, lambda: mk(base + ('ndf',), fe2=FE2_IBP, fe2A=FE2A)),
 'E2fix_CS_x3':          (14.3431, lambda: mk(base + ('nf12',), fe2=FE2_CS)),          # all new E2 terms times 3 (sensitivity)
}
CELLS = [('first_cell', s0['alpha_up'], s0['Ly_lo'], s0['Ly_hi'], True), ('tail_anchor', a_s, L_s, L_s * (1 + K.ETA_ANCHOR), False)]
out = dict(tol=TOL, cells=[c[:4] for c in CELLS], constants=CST, FE2_CS=FE2_CS, FE2_IBP=FE2_IBP, FE2A=FE2A, variants={})
_orig_extra = R.e2_extra
for name, (r143, fx) in VARS.items():
    R._CC = None; c = R.consts(); c.cl['r143'] = _c(r143)
    R.E2X = fx() if fx is not None else None
    if name.endswith('_x3'):
        R.e2_extra = lambda *a, **k: cc(3.0) * _orig_extra(*a, **k); R.E2X['fe2'] = _c(1 + 3 * (FE2_CS - 1))
    else: R.e2_extra = _orig_extra
    res = {}
    for cn, a, Llo, Lhi, fl in CELLS:
        t = time.time(); r = bb.bound_cell(0.0, a, Llo, Lhi, tol=TOL, rtol=0.005, floors=fl, cem_iters=25)
        res[cn] = dict(U=r['U'], lb=r['lb'], secs=round(time.time() - t, 1))
    out['variants'][name] = dict(r143=r143, cells=res)
    print(name, json.dumps(res), flush=True)
R.E2X = None; R.e2_extra = _orig_extra
ref = out['variants']['coded_fin13v4']['cells']
out['coded_reproduces_production'] = bool(ref['first_cell']['U'] == s0['rate_U'] and ref['tail_anchor']['U'] == tl['rate_s_U'])
# ---- jet box of the E2 coefficients: sup(|w|,|w'|) on [0,8] and [0,12], worst case over all shapes with Lyap <= s^2 (as closure.sramp_sup)
def box_sup(s_hi, Kc):
    c = R.consts(); n = 4000
    t0 = np.arange(n) / n; t1 = (np.arange(n) + 1) / n
    ct, st = R.circ(t0, t1); s = cc(float(np.nextafter(s_hi, np.inf)))
    Aa = s * ct; Bb = s * st * c.inv_sqrt_eps; W1 = sqrt(Aa * (Bb + c.c1 * Aa))
    s0_ = sqrt(Aa * (Kc * (c.PHX * Aa + c.cl['two'] * W1))); s1_ = sqrt(W1 * (Kc * (c.PHX * W1 + c.cl['two'] * Bb)))
    return float(np.max(mx(s0_, s1_).hi))
c = R.consts(); speak = float(np.nextafter(RT['sqrtL_max'], np.inf)) if isinstance(RT['sqrtL_max'], float) else None
if speak is None: speak = max(x['Ly_hi'] for x in RT['steps']) ** 0.5 * (1 + 1e-12)
out['box'] = dict(sqrtL_peak=speak, sup_w_w1_x_le_8=box_sup(speak, c.K["4<x<=8"]), sup_w_w1_x_le_12=box_sup(speak, c.K12), sup_w_w1_x_le_16p25=box_sup(speak, c.KR16),
                  alpha_max_of_run=max(x['alpha_up'] for x in RT['steps']))
print(json.dumps(out['box'])); print('coded reproduces production:', out['coded_reproduces_production'])
out['args'] = dict(argv=sys.argv, env={k: os.environ.get(k) for k in ("PYTHONHASHSEED", "NEW", "THETA", "ADD", "FIX", "S0V", "S1V", "DELTA", "MOD", "MODA", "ETA", "FLUX0", "FLUX2")},
                   inputs={p: hashlib.sha256(open(p, 'rb').read()).hexdigest() for p in [V + '/' + f for f in ('constants.py', 'rhs_iv.py', 'bb.py', 'ival.py')] +
                           [D + '/vp/rhs_iv.py', D + '/e2_consts.json', os.path.abspath(__file__), '/home/claude/s1spec/step3nl_v4/v4cert/out/constants_v4.json',
                            '/home/claude/s1spec/step3nl_v4/runs/P/rate_table_fin13v4_a0.00015.json']})
json.dump(out, open(f'e2_probe_tol{TOL:g}.json', 'w'), indent=1)
