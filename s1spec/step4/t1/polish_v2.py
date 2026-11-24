"""polish_all.py for the P1-seeded base (v2): slab 0 starts from the EXACT 128-bit seed (seed_p1.py, polstore-encoded),
not from its float copy; store directory selectable.  Otherwise identical to polish_all.py.
usage: python3 polish_v2.py TAG I_END eval_slabs --seed SEED.pkl [--store DIR]   -> DIR/slab_XXX.pkl + DIR/manifest.json"""
import sys, os, time, json, numpy as np
from polstore import save_slab, dec
import pickle
from base2 import Base
from arbops import SlabPolish, arb_refine
from rep import Mesh, Rep
from polish_chain import slab_norms, KEYS

def opt(name):
    if name not in sys.argv: return None
    i = sys.argv.index(name); v = sys.argv[i + 1]; del sys.argv[i:i + 2]; return v
SEEDF = opt('--seed'); SDIR = opt('--store'); STORE = SDIR is not None; RES = opt('--resume')   # --resume I: continue from stored slab I-1
assert SEEDF, '--seed SEED.pkl required'
tag, iend = sys.argv[1], int(sys.argv[2]); ev = set(int(v) for v in sys.argv[3].split(',')) if len(sys.argv) > 3 else set()
if STORE: os.makedirs(SDIR, exist_ok=True); man = []   # --store: also every 10th slab gets norms; slab_XXX.pkl + manifest.json
if STORE: ev |= set(range(0, iend + 1, 10))
B = Base(tag); RD = B.RD
out = dict(slabs=[], norms=[]); prev = None; t0 = time.time(); start = 0
if RES is not None:                                                # resume after an interruption: reuse the stored slabs 0..I-1
    from polstore import load_slab
    start = int(RES); man = json.load(open(f'{SDIR}/manifest.json'))[:start]
    assert len(man) == start and man[-1]['slab'] == start - 1
    prev, _ = load_slab(f'{SDIR}/slab_{start - 1:03d}.pkl')        # lossless: same arb state as in the interrupted run
    for l in open(f'polish_{tag}.log'):
        if l.startswith('NORMS '):
            d = json.loads(l[6:])
            if d['slab'] < start: out['norms'].append(d)
        elif l.startswith('{'):
            d = json.loads(l)
            if d['slab'] < start: out['slabs'].append(d)
    log = open(f'polish_{tag}.log', 'a'); log.write(f'RESUME at slab {start} (stored slabs 0..{start - 1} reused; fresh float Jacobian)\n'); log.flush()
else:
    log = open(f'polish_{tag}.log', 'w')
for i in range(start, iend + 1):
    R, D = B.rep(i); q = B.slabs[i]
    SP = SlabPolish(R, D, RD, q['nodes'], q['dt'])
    if prev is not None:
        Y = list(prev.Yn[-1]); Rp = prev.R
        while Rp.N < R.N:
            Cp = Rp.mesh.C; newC = [c for c in q['C'] if not np.any(Cp == c)]
            m1 = Mesh(JT=Rp.mesh.JT, XT=Rp.mesh.XT, h0=Rp.mesh.h0, cT=Rp.mesh.cT, C=np.sort(np.r_[Cp, newC[0]]))
            R1 = Rep(Rp.P, m1); Y = arb_refine(Y, Rp, R1, SP.AB); Rp = R1
        SP.Yn[0] = Y
        dY0 = max(abs(float(SP.Yn[0][k].mid()) - q['nodes'][0][k]) for k in range(R.n))      # polished vs float state
    else:                                                          # slab 0: the exact seed (128-bit dyadics, radius 0)
        from flint import ctx
        ctx.prec = 128; sd = pickle.load(open(SEEDF, 'rb')); Y = [dec(t) for t in sd['y_enc']]; assert len(Y) == R.n + 4
        dY0 = max(abs(float(Y[k].mid()) - q['nodes'][0][k]) for k in range(R.n + 4)); SP.Yn[0] = Y
    if prev is None or R.n != prev.R.n or i == start:              # float Jacobian refreshed after insertions, on resume (and adaptively below)
        J = D.jac(q['nodes'][len(q['nodes']) // 2]); RD.factor(J, q['dt'])
    ts = time.time(); fresh = False
    for it in range(120):
        em, dz = SP.step(J)
        if em < 1e-26 or dz < 1e-37: break
        if it == 25 and not fresh:                                 # slow: refresh the float Jacobian at the mid-slab state
            J = D.jac(np.array([float(v.mid()) for v in SP.Yn[len(SP.Yn) // 2]])); RD.factor(J, q['dt']); fresh = True
    rec = dict(slab=i, k=float(B.k(q['tau0'])), n=R.n, it=it + 1, maxE=em, t=time.time() - ts, dY0_float=dY0)
    out['slabs'].append(rec)
    log.write(json.dumps(rec) + '\n'); log.flush()
    if i in ev:
        mx, rows, th = slab_norms(SP, 6); mx.update(slab=i, k=rec['k'], n=R.n)
        out['norms'].append(mx); log.write('NORMS ' + json.dumps(mx) + '\n'); log.flush()
    if STORE:
        save_slab(f'{SDIR}/slab_{i:03d}.pkl', SP, i, q['tau0'], rec['k'], B.cfg, q['C'], rec)
        me = dict(slab=i, file=f'slab_{i:03d}.pkl', tau0=q['tau0'], tau1=q['tau0'] + q['dt'], k0=rec['k'], k1=float(B.k(q['tau0'] + q['dt'])),
                  dof=R.n + 4, n_el=R.N, maxE=em, norms=(out['norms'][-1] if i in ev else None))
        man.append(me); json.dump(man, open(f'{SDIR}/manifest.json', 'w'), indent=1)
    prev = SP
    if i % 10 == 0 or i == iend: json.dump(out, open(f'polish_{tag}.json', 'w'), indent=1)
log.write(f'DONE wall {time.time()-t0:.0f}s\n'); log.close()
json.dump(out, open(f'polish_{tag}.json', 'w'), indent=1)
