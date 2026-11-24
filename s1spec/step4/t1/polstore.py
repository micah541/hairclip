"""Lossless storage of a polished slab: every arb as exact (mid mantissa, mid exponent, rad mantissa, rad exponent)."""
import pickle, numpy as np
from flint import arb, ctx

def enc(v):
    m, e = v.mid().man_exp(); r, re = v.rad().man_exp()      # rad() is an exact arb with zero radius
    return (int(m), int(e), int(r), int(re))
def dec(t):
    m, e, r, re = t
    if r == 0: return arb(mid=arb((m, e)))
    k = max(0, 32 - r.bit_length())   # the arb(mid, rad) constructor rounds the radius up by one mag-ulp: pre-subtract it
    return arb(mid=arb((m, e)), rad=arb(((r << k) - 1, re - k)))

def save_slab(path, SP, i, tau0, k, cfg, C, rec):
    d = dict(slab=i, tau0=tau0, dt=SP.dt, k=k, cfg=cfg, C=np.array(C), n=SP.R.n, N=SP.R.N, P=SP.R.P, s=SP.s,
             prec=ctx.prec, nodes_rel=[float(v) for v in SP.nd],
             Yn=[[enc(v) for v in row] for row in SP.Yn],   # (s+1) x (n+4): rows = relative times 0, c_1..c_s; cols = DOF[0:n], ell, lnlam, c, t
             rec=rec)
    pickle.dump(d, open(path, 'wb'), protocol=4)

def load_slab(path):
    """-> (SlabPolish with Yn set, dict). Rebuilds mesh/Rep/Dyn/Radau from the stored cfg and knot parameters C."""
    from rep import Rep, Mesh
    from ops import Dyn
    from radau import Radau
    from arbops import SlabPolish
    d = pickle.load(open(path, 'rb')); cfg = d['cfg']
    m = Mesh(JT=cfg['JT'], h0=cfg['h0'], cT=cfg['cT'], C=d['C']); R = Rep(cfg['P'], m); D = Dyn(R); RD = Radau(cfg['s'])
    ctx.prec = d['prec']
    Y = [[dec(t) for t in row] for row in d['Yn']]
    SP = SlabPolish(R, D, RD, [[float(v.mid()) for v in row] for row in Y], d['dt'], d['prec'])
    SP.Yn = Y
    return SP, d
