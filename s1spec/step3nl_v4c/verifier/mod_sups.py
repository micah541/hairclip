#!/usr/bin/env python3
"""Certified modulation sups sup|(a,b)|, sup|beta| over all march steps and the tail, and the summary final_alpha_points_5e-8.json.

Bounds (IBP modulation lemma, agent_mod/ibp/MOD_IBP.md, mod_ibp.json; constants of constants.py with MOD=ibp, valid for |alpha|<=5e-4, A,W1<=1e-3):
   |beta|  <= (BETA0 + BETA1 |alpha|) A + QB  Q + BETAF  Fn(alpha) + STRIP_BETA
   |(a,b)| <= (CAB0  + MAB   |alpha|) A + QAB Q + PHIAB  Fn(alpha) + STRIP_AB          (Euclidean norm of (a,b); hence also sup|a|, sup|b|)
with  Q = (A+W1+B)^2 <= qbl_sup * Ly,   ||chi Phi_F|| <= Fn(alpha),  strip terms = 0.02 * (strip sums per unit S).
Per march step k: A <= far A_max_k (worst case over shapes with Lyap <= Ly_hi_k), alpha <= alpha_up_k, Ly <= Ly_hi_k (monotone increasing in all three).
Tail: alpha <= a*, Ly <= L* (L* <= peak^2), A <= sqrt(L*) (A = s cw ct <= s): one extra row.  Every row is evaluated in arb (120 bits) and the max is
rounded UP to 6 significant decimal digits.
usage (verifier environment, see run_final.sh):  python3 mod_sups.py TAG [OUTDIR]   -> OUTDIR/final_alpha_points_5e-8.json (default OUTDIR = cwd)"""
import json, os, sys, math
from flint import arb, ctx
ctx.prec = 120
import constants as K, unc
A = lambda v: arb(str(v)) if not isinstance(v, arb) else v
def up6(x):                                   # outward decimal: smallest 6-digit decimal >= upper endpoint
    u = float(x.upper()); e = math.floor(math.log10(u)) - 5; m = math.ceil(u / 10 ** e * (1 - 1e-12)); 
    while arb(m) * arb(10) ** e < x.upper(): m += 1
    return float(f"{m}e{e}")
def Fn(al): return A(K.FN[0]) * al + A(K.FN[1]) * al * al + A(K.FN[2]) * al * al * al
def bounds(al, L, Amax):
    Q = A(repr(unc.qbl_sup())) * L
    b = (A(K.BETA0) + A(K.BETA1) * al) * Amax + A(K.QB) * Q + A(K.BETAF) * Fn(al) + A(K.STRIP_BETA)
    ab = (A(K.CAB0) + A(K.MAB) * al) * Amax + A(K.QAB) * Q + A(K.PHIAB) * Fn(al) + A(K.STRIP_AB)
    return ab, b
def sups(d):
    ab_m = b_m = None
    for s in d['steps']:
        ab, b = bounds(A(repr(s['alpha_up'])), A(repr(s['Ly_hi'])), A(repr(s['far']['A_max'])))
        ab_m = ab if ab_m is None or ab.upper() > ab_m.upper() else ab_m; b_m = b if b_m is None or b.upper() > b_m.upper() else b_m
    t = d['tail']; Ls = A(repr(t['L_star']))
    ab, b = bounds(A(repr(t['a_star'])), Ls, Ls.sqrt())
    ab_m = ab if ab.upper() > ab_m.upper() else ab_m; b_m = b if b.upper() > b_m.upper() else b_m
    return up6(ab_m), up6(b_m)
def summary(tdir, tag):
    out = {}
    for key, al in (('1.25e-4', '0.000125'), ('1.5e-4', '0.00015')):
        d = json.load(open(os.path.join(tdir, f'rate_table_{tag}_a{al}.json'), encoding='utf-8')); t = d['tail']; u = d['unc']; fm = d['far_max']
        sab, sb = sups(d)
        out[key] = {'steps': len(d['steps']), 'peak': d['sqrtL_max'], 'tau_tail': d['tau_tail_start'], 'rho': t['rho'], 'sqrtLF': t['sqrt_LF'],
            'alphaF': t['alpha_F'], 'a_star': t['a_star'], 'S_ramp': d['S_ramp_max_along_trajectory'], 'sI1': fm['sI1'], 'sI2': fm['sI2'],
            'M1': fm['M1_ramp'], 'A': fm['A_max'], 'W1': fm['W1_max'], 'sup_ab': sab, 'sup_beta': sb, 'delta_tail': u['delta_tail'],
            'J': {n: [u[n]['const'], u[n]['slope_per_tau']] for n in ('V0sq', 'V1sq', 'V2sq', 'V4sq')},
            'ramp_levels': sorted({(K.v4_level(s['alpha_up']) if K.NEW else K.ramp_level(s['alpha_up'])) for s in d['steps']})}
    return out
if __name__ == '__main__':
    tag = sys.argv[1]; outdir = sys.argv[2] if len(sys.argv) > 2 else '.'
    tdir = os.environ.get('TABLE_DIR', '.')
    s = summary(tdir, tag); p = os.path.join(outdir, 'final_alpha_points_5e-8.json')
    import hashlib
    _ins = [os.path.abspath(os.path.join(tdir, f'rate_table_{tag}_a{al}.json')) for al in ('0.000125', '0.00015')] + [os.path.join(os.path.dirname(os.path.abspath(__file__)), f) for f in ('mod_sups.py', 'constants.py', 'unc.py')] + ([os.path.join(K.SRC, 'v4bcert', 'out', f) for f in ('e0v4b_constants.json', 'r2_strip_ibp.json', 'e2v4b_constants.json')] if getattr(K, 'V4B', False) else [])
    if K.NEW: _ins.append(os.path.join(K.SRC, 'v4cert', 'out', 'constants_v4.json'))
    s['args'] = dict(argv=sys.argv, env={k_: os.environ.get(k_) for k_ in ('PYTHONHASHSEED', 'NEW', 'MOD', 'FIX', 'TABLE_DIR') + (('V4B',) if getattr(K, 'V4B', False) else ())}, inputs={q_: hashlib.sha256(open(q_, 'rb').read()).hexdigest() for q_ in _ins},
                     constants=dict(BETA0=K.BETA0, BETA1=K.BETA1, QB=K.QB, BETAF=K.BETAF, STRIP_BETA=K.STRIP_BETA, CAB0=K.CAB0, MAB=K.MAB, QAB=K.QAB, PHIAB=K.PHIAB, STRIP_AB=K.STRIP_AB, FN=list(K.FN)))
    json.dump(s, open(p, 'w', encoding='utf-8'), indent=1); print('wrote', p); print({k: (v['sup_ab'], v['sup_beta']) for k, v in s.items() if k != 'args'})
