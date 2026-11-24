"""T1 float driver v2: as run2.py, but seeded from an exact P1 seed (seed_p1.py output; float copy y_float only as the
solver's start, the arb polish uses the exact 128-bit seed).  usage: python3 run3.py TAG KEND SEED.pkl
Output: slabs_<TAG>.pkl (run2 format, readable by base2.Base), run_<TAG>.log."""
import os, sys, time, json, pickle, numpy as np
os.environ.setdefault('OMP_NUM_THREADS', '1')
from rep import Rep, Mesh, K0, DS1
from ops import Dyn, NS
from radau import Radau
from run import geom
from run2 import refine


def main():
    tag, KEND, seedf = sys.argv[1], float(sys.argv[2]), sys.argv[3]
    sd = pickle.load(open(seedf, 'rb')); cfg = dict(sd['cfg']); cfg.update(tag=tag, KEND=KEND, seed=os.path.abspath(seedf))
    P, h0, c, beta, rg, s, dt = cfg['P'], cfg['h0'], cfg['c'], cfg['beta'], cfg['rg'], cfg['s'], cfg['dt']
    R = Rep(P, Mesh(JT=cfg['JT'], h0=h0, cT=cfg['cT'], C=sd['C'])); D = Dyn(R); RD = Radau(s)
    Y = np.array(sd['y_float'], float); assert len(Y) == R.n + NS
    g0 = geom(R, D, Y)
    log = open(f'run_{tag}.log', 'w')
    log.write(json.dumps(cfg) + '\n'); log.write(f"seed {cfg['seed']} k={sd['k']} ell={Y[R.n]} n={R.n} N={R.N} turn+pi={g0['turn']+np.pi:.3e} ImG(ell)={g0['G_end'].imag:.3e}\n"); log.flush()
    slabs = []; tau = 0.0; Zg = None; t0 = time.time(); k = K0 / np.exp(Y[R.n + 1]); ks = 0
    while k < KEND:
        ts = time.time(); ins = 0
        while True:
            ell = Y[R.n]; x, _ = R.mesh.knots(ell); h = np.diff(x)
            xm = 0.5 * (x[:-1] + x[1:]); ht = 1 / (1 / (h0 + c * xm) + 1 / (beta * ell))
            bad = [e for e in range(R.mesh.JT, R.N - 1) if h[e] > rg * ht[e]]
            if not bad: break
            e = bad[0]; R, ynew = refine(R, Y[:R.n], ell, xm[e]); Y = np.r_[ynew, Y[-NS:]]; ins += 1
        if ins: D = Dyn(R); Zg = None
        sc = np.ones(len(Y)); sc[R.n:] = np.maximum(np.abs(Y[R.n:]), 1.0)
        Yj = Y if Zg is None else Y + Zg[len(Zg) // 2]
        J = D.jac(Yj); tj = time.time() - ts
        RD.factor(J, dt)
        Ys, hist = RD.solve(D.F, Y, dt, Z0=Zg, scale=sc, tol=1e-14, maxit=60)
        nodes = np.vstack([Y, Ys])
        slabs.append(dict(tau0=tau, dt=dt, C=R.mesh.C.copy(), nodes=nodes))
        W = RD.interp_matrix(1.0 + RD.c); Zg = W @ nodes - Ys[-1]
        Y = Ys[-1].copy(); tau += dt; ks += 1
        k = K0 / np.exp(Y[R.n + 1])
        g = geom(R, D, Y); p = D.parts(Y)
        log.write(f"slab {ks:4d} tau={tau:8.4f} k={k:10.3f} a={p['a']:+.5e} ell={Y[R.n]:.6e} N={R.N} n={R.n} ins={ins} it={len(hist):2d} last={hist[-1]:.1e} "
                  f"turn+pi={g['turn']+np.pi:+.2e} ImG={g['G_end'].imag:+.2e} wall={time.time()-ts:.1f}s (jac {tj:.1f}s)\n"); log.flush()
        if ks % 10 == 0 or k >= KEND:
            pickle.dump(dict(cfg=cfg, slabs=slabs, c=RD.c), open(f'slabs_{tag}.pkl', 'wb'))
    log.write(f"DONE {ks} slabs, wall {time.time()-t0:.0f}s\n"); log.close()


if __name__ == '__main__':
    main()
