"""E0 audit probe: certified cell bounds (bb.bound_cell of the fin13v4 verifier, unmodified code, read only) with the E0 constants
PI3, f572, MOD_ALPHA, ALIN replaced in memory.  [cert given the constants]: U is a rigorous upper bound of the rate over the cell for the
patched rhs, lb a rigorous lower bound of that sup (U - lb <= tol).  Writes e0_probe.json (cwd).  Env: as verifier/run_v4.sh.
usage: python3 -B e0_probe.py [tol]"""
import sys, os, json, time, hashlib
V = '/home/claude/s1spec/step3nl_v4/verifier'
sys.path.insert(0, V)
import constants as K, rhs_iv as R, bb
TOL = float(sys.argv[1]) if len(sys.argv) > 1 else 5e-4
RT = json.load(open('/home/claude/s1spec/step3nl_v4/runs/P/rate_table_fin13v4_a0.00015.json'))
tl = RT['tail']; s0 = RT['steps'][0]
a_s, L_s = tl['a_star'], tl['L_star']
ALIN0 = dict(K.ALIN); ALINF = dict(ALIN0, sq=59.33, BB=6.22)       # raw alin_E0_spec_0.001.json: spec01 59.32715, |a22| 6.21990 (typed 59.3 / 6.2 are rounded DOWN)
VARS = {   # name: (PI3, f572, MOD_ALPHA, ALIN)
 'coded_fin13v4':            (0.01607,              5.72,              25.0, ALIN0),
 'v4_probe_E0fix':           (0.018347067272215117, 5.937087030370929, 25.0, ALIN0),
 'coded+ALINup':             (0.01607,              5.72,              25.0, ALINF),
 'T1_truenorm_crude':        (0.01877887,           6.0533,            43.3, ALINF),
 'T2_projected_nosym':       (0.01574642,           5.72,              38.4, ALINF),
 'T3_sharp':                 (0.0104,               5.72,              38.4, ALINF),
 'T3_sharp_MODA40':          (0.0104,               5.72,              40.0, ALINF),
}
CELLS = [('first_cell', s0['alpha_up'], s0['Ly_lo'], s0['Ly_hi'], True),
         ('tail_anchor', a_s, L_s, L_s * (1 + K.ETA_ANCHOR), False)] + \
        [(f'anchor_alpha={a:g}', a, L_s, L_s * (1 + K.ETA_ANCHOR), False) for a in (2e-5, 1e-5, 1e-6)]
out = dict(tol=TOL, cells=[c[:4] for c in CELLS], variants={})
for name, (pi3, f572, moda, alin) in VARS.items():
    K.PI3 = pi3; K.MOD_ALPHA = moda; K.ALIN = alin; R._CC = None
    c = R.consts(); c.cl['f572'] = R._c(f572)
    assert abs(float(c.cl['mod'].hi) - moda) < 1e-9 and abs(float(c.cl['alin_sq'].hi) - alin['sq']) < 1e-9
    res = {}
    for cn, a, Llo, Lhi, fl in CELLS:
        t = time.time(); r = bb.bound_cell(0.0, a, Llo, Lhi, tol=TOL, rtol=0.005, floors=fl, cem_iters=25)
        res[cn] = dict(U=r['U'], lb=r['lb'], secs=round(time.time() - t, 1))
    out['variants'][name] = dict(PI3=pi3, f572=f572, MOD_ALPHA=moda, ALIN_sq=alin['sq'], ALIN_BB=alin['BB'], cells=res)
    print(name, json.dumps(res), flush=True)
out['args'] = dict(argv=sys.argv, env={k: os.environ.get(k) for k in ("PYTHONHASHSEED", "NEW", "THETA", "ADD", "FIX", "S0V", "S1V", "DELTA", "MOD", "MODA", "ETA", "FLUX0", "FLUX2")},
                   inputs={p: hashlib.sha256(open(p, 'rb').read()).hexdigest() for p in [V + '/' + f for f in ('constants.py', 'rhs_iv.py', 'bb.py', 'ival.py')] +
                           ['/home/claude/s1spec/step3nl_v4/v4cert/out/constants_v4.json', '/home/claude/s1spec/step3nl_v4/runs/P/rate_table_fin13v4_a0.00015.json']})
json.dump(out, open(f'e0_probe_tol{TOL:g}.json', 'w'), indent=1)
