"""Interval evaluation of the IBP'd modulation pairings over all hull boxes.   usage: python3 ibp_run.py ALPHA1 ETA  -> run_<alpha>_<eta>.json"""
import sys, json, time, pickle
from flint import arb, ctx
ctx.prec = 150
from ibp_gen import linf as LinF, coeff as CoefF, remf as RemF
import ibp_consts as K
import ibp_boxes as BX
ALPHA1 = float(sys.argv[1]); ETA = float(sys.argv[2]); NC = int(sys.argv[3]) if len(sys.argv) > 3 else 4; TAG = f"{ALPHA1:g}_{ETA:g}_n{NC}"
AL = arb(0, ALPHA1); ETAB = arb(0, ETA)
up = lambda z: float(z.upper())
mag = lambda z: arb(abs(z).upper())
def mx(a, b): return a if a.upper() >= b.upper() else b
CN = ['00', '01', '11', '02', '12', '22']
comps = ['a', 'b', 'beta']
def newagg(): return {i: {'gam': {c: arb(0) for c in CN}, 'M3': arb(0), 'M2': arb(0), 'M1': arb(0)} for i in range(3)}
core = newagg(); near = newagg(); far = newagg(); strip = newagg(); farL1 = newagg()
XSPLIT = 16.25
sq_j = [arb(0) for _ in range(3)]; sq_i = [arb(0) for _ in range(3)]; off_ab = arb(0)
L1k = [arb(0) for _ in range(3)]; L1k_i = [arb(0) for _ in range(3)]
where = {i: {c: (0.0, None) for c in CN} for i in range(3)}
supk = [0.0] * 3
t0 = time.time(); nb = 0
for b in BX.boxes():
    nb += 1
    k_, m_, p_, q_ = b['k'], b['m'], b['p'], b['q']
    ps, hh = b['ps'], b['hh']
    a = [AL * ps[i] + AL * AL * hh[i] for i in range(5)]
    alb = 1 - k_ * a[0]; Jb = 1 / (alb * alb + a[1] * a[1]).sqrt()
    cc = CoefF(k_, m_, p_, q_, *a, alb, Jb)
    cells = []
    for ci in range(NC):
        for cj in range(NC):
            w0 = arb(-ETA + (2 * ci + 1) * ETA / NC, ETA / NC * (1 + 1e-12)); w1 = arb(-ETA + (2 * cj + 1) * ETA / NC, ETA / NC * (1 + 1e-12))
            alf = 1 - k_ * (a[0] + w0); Jf = 1 / (alf * alf + (a[1] + w1) * (a[1] + w1)).sqrt()
            cells.append(RemF(k_, m_, p_, q_, *a, w0, w1, alf, Jf, alb, Jb))
    # linear part:  kappa_j = d_alpha Lam (Y -> psih_j)
    kap = []
    for j in (1, 2, 3):
        Y = b['psih'][j]
        kap.append(LinF(*Y, *a, *ps, *hh, AL, k_, m_, p_, q_, alb, Jb)['LamA'])
    isstrip = b['strip']; W = b['w']
    kap_i = [sum((K.H[i][j] * kap[j] for j in range(3)), arb(0)) for i in range(3)]
    if not isstrip:
        for j in range(3): sq_j[j] += W * mag(kap[j]) * mag(kap[j])
        for i in range(3): sq_i[i] += W * mag(kap_i[i]) * mag(kap_i[i])
        off_ab += W * mag(kap_i[0]) * mag(kap_i[1])
    else:
        for j in range(3): L1k[j] += W * mag(kap[j])
        for i in range(3): L1k_i[i] += W * mag(kap_i[i])
    for i in range(3):
        Z = [sum((K.H[i][j] * b['psih'][j + 1][r] for j in range(3)), arb(0)) for r in range(3)]
        isfar = (not isstrip) and float(b['X'].mid()) >= XSPLIT
        def upd(tgt, nm, v, c=None):
            if c is None:
                if v.upper() > tgt[nm].upper(): tgt[nm] = v
            else:
                if v.upper() > tgt['gam'][c].upper(): tgt['gam'][c] = v
        for c in CN:
            g = sum((Z[r] * cc[f'c{c}_{r}'] for r in range(3)), arb(0)); gm = mag(g)
            if isstrip: strip[i]['gam'][c] += W * gm
            else:
                upd(core[i], None, gm, c); upd(far[i] if isfar else near[i], None, gm, c)
                if isfar: farL1[i]['gam'][c] += W * gm
                if up(gm) > where[i][c][0]: where[i][c] = (up(gm), float(b['X'].mid()))
        def Mq(prefix, n):
            best = arb(0)
            for tt in cells:
                for s in range(n):
                    v = mag(sum((Z[r] * tt[f'{prefix}_{s}_{r}'] for r in range(3)), arb(0)))
                    best = mx(best, v)
            return best
        for nm, pre, n in (('M3', 't3', 4), ('M2', 't2', 3), ('M1', 't1', 2)):
            v = Mq(pre, n)
            if isstrip: strip[i][nm] += W * v
            else:
                upd(core[i], nm, v); upd(far[i] if isfar else near[i], nm, v)
                if isfar: farL1[i][nm] += W * v
    for j in range(3): supk[j] = max(supk[j], up(mag(kap[j])))
    if nb % 500 == 0: print(nb, f"[{time.time()-t0:.0f}s]", flush=True)
FM = lambda A_: {comps[i]: {'gam': {c: up(A_[i]['gam'][c]) for c in CN}, 'M3': up(A_[i]['M3']), 'M2': up(A_[i]['M2']), 'M1': up(A_[i]['M1'])} for i in range(3)}
out = dict(alpha1=ALPHA1, eta=ETA, nboxes=nb,
           core=FM(core), near=FM(near), far=FM(far), far_L1=FM(farL1), strip_L1=FM(strip), xsplit=XSPLIT,
           where_gam_max={comps[i]: {c: where[i][c] for c in CN} for i in range(3)},
           m_j=[up(x.sqrt()) for x in sq_j], m_i=dict(zip(comps, [up(x.sqrt()) for x in sq_i])), offdiag_ab=up(off_ab), L1k_strip_j=[up(x) for x in L1k], L1k_strip_i=dict(zip(comps, [up(x) for x in L1k_i])), sup_kappa=supk)
json.dump(out, open(f'run_{TAG}.json', 'w'), indent=1)
print("m_j (per alpha, coefficient of |alpha| A):", out['m_j'], " m_i:", out['m_i'], f"[{time.time()-t0:.0f}s]")
