"""T1 float driver (Lagrangian outer mesh with exact knot insertion).
usage: python3 run2.py TAG KEND P JT h0 cT NO c beta rg s dt
Seed: even part of an HLS snapshot (SEED env), closure corrected; Radau IIA slabs in tau until k >= KEND.
Output: slabs_<TAG>.pkl = list of dicts(tau0, dt, xi, nodes[(s+1), n+4]) + cfg; run_<TAG>.log."""
import os, sys, time, json, pickle, numpy as np
os.environ.setdefault('OMP_NUM_THREADS', '1')
from rep import Rep, Mesh, K0, DS1
from ops import Dyn, NS
from radau import Radau
from seed import HLSHalf
from run import seed_closed, geom


def refine(R, y, ell, xnew):
    """exact transfer of the state y (without scalars) to the mesh with an outer knot inserted at xnew"""
    B = R.B; m0 = R.mesh; m1 = m0.inserted(xnew, ell); R1 = Rep(R.P, m1)
    x0, _ = m0.knots(ell); h0 = np.diff(x0); x1, _ = m1.knots(ell); h1 = np.diff(x1)
    e = int(np.searchsorted(x0, xnew) - 1)                       # old element containing xnew
    D0 = R.nodal(y, h0)
    Dn = np.array([R.eval_at(y, ell, np.array([xnew]), mm)[0] for mm in range(5)])
    D1 = np.vstack([D0[:e + 1], Dn[None, :], D0[e + 1:]])
    bub0 = np.array([y[R.bidx[k]] for k in range(R.N)]) if B.nbub else np.zeros((R.N, 0))
    bub1 = np.zeros((R1.N, B.nbub))
    bub1[:e] = bub0[:e]; bub1[e + 2:] = bub0[e + 1:]
    y1 = R1.to_y(D1, np.zeros((R1.N, B.nbub)), h1)              # bubbles zero, nodal exact
    if B.nbub:
        d_old = R.local(y, h0)[e]
        for k in (e, e + 1):
            u = (x1[k] + h1[k] * B.uq - x0[e]) / h0[e]
            vals = B.evalpts(u, 0) @ d_old
            dl = R1.local(y1, h1)[k]
            rem = vals - B.V[0] @ dl
            Vb = B.V[0][:, 10:]
            bub1[k] = np.linalg.lstsq(Vb, rem, rcond=None)[0]
        for k in range(R1.N): y1[R1.bidx[k]] = bub1[k]
    return R1, y1


def main():
    a = sys.argv[1:]
    tag, KEND = a[0], float(a[1])
    P, JT, h0, cT, NO, c, beta, rg, s, dt = int(a[2]), int(a[3]), float(a[4]), float(a[5]), int(a[6]), float(a[7]), float(a[8]), float(a[9]), int(a[10]), float(a[11])
    seedf = os.environ.get('SEED', '/home/claude/s1spec/step4/feas2/snaps/p2_k20.npz')
    cfg = dict(tag=tag, KEND=KEND, P=P, JT=JT, h0=h0, cT=cT, NO=NO, c=c, beta=beta, rg=rg, s=s, dt=dt, seed=seedf, XT=40.0)
    S = HLSHalf(seedf)
    R = Rep(P, Mesh(JT=JT, h0=h0, cT=cT, NO=NO, c=c, beta=beta, ell0=S.ell)); D = Dyn(R); RD = Radau(s)
    y = seed_closed(R, D, S)
    Y = np.r_[y, S.ell, np.log(S.lam), 0.0, S.t]
    g0 = geom(R, D, Y)
    log = open(f'run_{tag}.log', 'w')
    log.write(json.dumps(cfg) + '\n'); log.write(f"seed k={S.kmax} lam={S.lam} ell={S.ell} n={R.n} N={R.N} turn+pi={g0['turn']+np.pi:.3e} ImG(ell)={g0['G_end'].imag:.3e}\n"); log.flush()
    slabs = []; tau = 0.0; Zg = None; t0 = time.time(); k = S.kmax; ks = 0
    while k < KEND:
        ts = time.time()
        # element splitting in the outer zone (exact refinement, at slab boundaries only)
        ins = 0
        while True:
            ell = Y[R.n]; x, _ = R.mesh.knots(ell); h = np.diff(x)
            xm = 0.5 * (x[:-1] + x[1:]); ht = 1 / (1 / (h0 + c * xm) + 1 / (beta * ell))
            bad = [e for e in range(R.mesh.JT, R.N - 1) if h[e] > rg * ht[e]]
            if not bad: break
            e = bad[0]; R, ynew = refine(R, Y[:R.n], ell, xm[e]); Y = np.r_[ynew, Y[-NS:]]; ins += 1
        if ins: D = Dyn(R); Zg = None
        sc = np.ones(len(Y)); sc[R.n:] = np.maximum(np.abs(Y[R.n:]), 1.0)
        Yj = Y if Zg is None else Y + Zg[len(Zg) // 2]                  # Jacobian at the (extrapolated) mid-slab state
        J = D.jac(Yj); tj = time.time() - ts
        RD.factor(J, dt)
        Ys, hist = RD.solve(D.F, Y, dt, Z0=Zg, scale=sc, tol=1e-14, maxit=60)
        nodes = np.vstack([Y, Ys])
        slabs.append(dict(tau0=tau, dt=dt, C=R.mesh.C.copy(), nodes=nodes))
        W = RD.interp_matrix(1.0 + RD.c); Zg = W @ nodes - Ys[-1]
        Y = Ys[-1].copy(); tau += dt; ks += 1
        lam = np.exp(Y[R.n + 1]); k = K0 / lam
        g = geom(R, D, Y); p = D.parts(Y)
        log.write(f"slab {ks:4d} tau={tau:8.4f} k={k:10.3f} a={p['a']:+.5e} ell={Y[R.n]:.6e} N={R.N} n={R.n} ins={ins} it={len(hist):2d} last={hist[-1]:.1e} "
                  f"turn+pi={g['turn']+np.pi:+.2e} ImG={g['G_end'].imag:+.2e} wall={time.time()-ts:.1f}s (jac {tj:.1f}s)\n"); log.flush()
        if ks % 10 == 0 or k >= KEND:
            pickle.dump(dict(cfg=cfg, slabs=slabs, c=RD.c), open(f'slabs_{tag}.pkl', 'wb'))
    log.write(f"DONE {ks} slabs, wall {time.time()-t0:.0f}s\n"); log.close()


if __name__ == '__main__':
    main()
