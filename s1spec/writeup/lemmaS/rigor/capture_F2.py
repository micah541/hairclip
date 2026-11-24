#!/usr/bin/env python3
"""Rerun of the FINAL lemmaS_F2.py command (tag F10p) with F1.assemble wrapped to record its inputs G (incl. the ramp tables rampd that F2 does
not store) and the other assemble arguments.  Writes rigor/cap_F10p.json (+ rigor/lemmaS_rigcap_M0.52.json = the rerun's own output; compare with lemmaS_F10p_M0.52.json).
Rerun:  cd /home/claude/s1spec/writeup/lemmaS && PYTHONHASHSEED=0 python3 rigor/capture_F2.py     (about 2 min; CPU heavy; run alone)"""
import os, sys, json, shutil
HERE = os.path.dirname(os.path.abspath(__file__)); PAR = os.path.dirname(HERE)
os.chdir(PAR); sys.path.insert(0, PAR)
import lemmaS_F1 as F, lemmaS_F2 as F2
from flint import arb
def ser(x):
    if isinstance(x, arb):
        m, e = x.mid().man_exp(); return {'arb': [str(m), int(e.__int__()) if hasattr(e, '__int__') else e, float(x.rad())]}
    if isinstance(x, dict): return {str(k): ser(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)): return [ser(v) for v in x]
    if isinstance(x, (int, float, str)) or x is None: return x
    try: return float(x)
    except Exception: return str(x)
calls = []; orig = F.assemble
def wrap(G, J, sqLam, gam, alpha1, **kw):
    calls.append(dict(G={str(j): {k: ser(v) for k, v in G[j].items()} for j in G}, J=ser(J), sqLam=sqLam, gam=gam, alpha1=alpha1, kw=ser(kw)))
    r = orig(G, J, sqLam, gam, alpha1, **kw); calls[-1]['result'] = ser(r); return r
F.assemble = wrap
rl = []; orig_rl = F2.ramp_level
def wrap_rl(j, T, Abr, Ghr, YCr, Rm, gam, tmin=1e-3):
    out = orig_rl(j, T, Abr, Ghr, YCr, Rm, gam, tmin=tmin)
    rl.append(dict(j=j, Abr={str(k): ser(v) for k, v in Abr.items()}, Ghr=ser(Ghr), YCr=ser(YCr), out=ser(out[0]))); return out
F2.ramp_level = wrap_rl
sys.argv = ['lemmaS_F2.py', '--J', '6.0273e-12', '8.4673e-12', '1.8708e-12', '--handoff', '2e-7', '1e-4', '0.2', '--sqLam', '5.3075e-8', '--alpha1', '1.5e-4',
            '--core', '0.11', '0.52', '--ramp', '0.05', '0.9', '--mhr', '0', '--slopes', '8.6369e-19', '1.2134e-18', '2.6808e-19', '--rho', '2.693e-3', '--M1s', '0.52',
            '--J4', '6.2362e-8', '--sJ4', '8.9363e-15', '--rampl4', '--rs33', 'ibp2', '--l4', '--abar', '2.19103e-7', '--beta', '3.33811e-9', '--tag', 'rigcap']
F2.main()
json.dump(dict(assemble=calls, ramp_level=rl), open(os.path.join(HERE, 'cap_F10p.json'), 'w'))
shutil.move(os.path.join(PAR, 'lemmaS_rigcap_M0.52.json'), os.path.join(HERE, 'lemmaS_rigcap_M0.52.json'))
