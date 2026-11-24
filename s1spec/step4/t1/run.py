"""T1 float driver: seed from an HLS snapshot (even part), Radau IIA slabs in tau until k >= KEND.
usage: python3 run.py TAG [KEND] [P] [JT] [h0] [cT] [NO] [c] [beta] [s] [dt]"""
import os, sys, time, json, numpy as np
os.environ.setdefault('OMP_NUM_THREADS', '1')
from rep import Rep, Mesh, K0, DS1
from ops import Dyn, NS
from radau import Radau
from seed import HLSHalf
from projtest import project

def geom(R, D, Y):
    """closure diagnostics: int kappa - pi, Im Gamma(ell), Gamma(ell), theta(ell) + pi/2 ... (float quadrature)"""
    p = D.parts(Y); B = R.B; h = p['h']
    th_e = np.r_[0.0, np.cumsum(h * (p['k0'] @ B.wq))] - np.pi / 2
    th = th_e[:-1, None] + h[:, None] * (p['k0'] @ B.Icum.T)
    E = np.exp(1j * th)
    G_e = np.r_[0.0, np.cumsum(h * (E @ B.wq))] - DS1
    return dict(turn=th_e[-1] + np.pi / 2, G_end=G_e[-1], theta_end=th_e[-1])

def seed_closed(R, D, S, log=None):
    """projection of the even part of the HLS state plus eps * w'(x) (zero mean, macro region) so that Im Gamma(ell) = 0"""
    xc = 12.0 / S.lam; sg = 2.0 / S.lam               # Gaussian bump at physical s = 12, width 2
    def wder(x, m):                                   # m-th derivative of w(x) = exp(-((x-xc)/sg)^2)
        t = (np.asarray(x) - xc) / sg
        Hm = np.polynomial.hermite.hermval(t, np.r_[np.zeros(m), 1.0])
        return (-1.0) ** m * Hm * np.exp(-t * t) / sg ** m
    eps = 0.0
    for it in range(4):
        F = lambda x, m: S.kap(x, m) + eps * wder(x, m + 1)
        y = project(R, F, S.ell)
        Y = np.r_[y, S.ell, 0.0, 0.0, 0.0]; g = geom(R, D, Y)
        # d ImGamma(ell)/d eps = int w(x) cos(theta(x)) dx (theta from the projected state)
        p = D.parts(Y); B = R.B; h = p['h']
        th_e = np.r_[0.0, np.cumsum(h * (p['k0'] @ B.wq))] - np.pi / 2
        th = th_e[:-1, None] + h[:, None] * (p['k0'] @ B.Icum.T)
        dIm = np.sum(h[:, None] * B.wq[None, :] * wder(p['xq'], 0) * np.cos(th))
        if log: log.write(f"  closure it {it}: eps={eps:.3e} ImG={g['G_end'].imag:.3e} turn+pi={g['turn']+np.pi:.3e}\n")
        eps -= g['G_end'].imag / dIm
    return y

def main():
    a = sys.argv[1:]
    tag = a[0]; KEND = float(a[1]) if len(a) > 1 else 2.3e4
    P = int(a[2]) if len(a) > 2 else 9; JT = int(a[3]) if len(a) > 3 else 40
    h0 = float(a[4]) if len(a) > 4 else 0.35; cT = float(a[5]) if len(a) > 5 else 0.08
    NO = int(a[6]) if len(a) > 6 else 80; c = float(a[7]) if len(a) > 7 else 0.12; beta = float(a[8]) if len(a) > 8 else 0.02
    s = int(a[9]) if len(a) > 9 else 8; dt = float(a[10]) if len(a) > 10 else 0.15
    seedf = os.environ.get('SEED', '/home/claude/s1spec/step4/feas2/snaps/p2_k20.npz')
    cfg = dict(tag=tag, KEND=KEND, P=P, JT=JT, h0=h0, cT=cT, NO=NO, c=c, beta=beta, s=s, dt=dt, seed=seedf)
    R = Rep(P, Mesh(JT=JT, h0=h0, cT=cT, NO=NO, c=c, beta=beta)); D = Dyn(R); RD = Radau(s)
    S = HLSHalf(seedf)
    y = seed_closed(R, D, S)
    Y = np.r_[y, S.ell, np.log(S.lam), (S.X0[0] if False else 0.0), S.t]
    # physical centre c: tip_phys = c - lam*DS1 (axis along the real line of the HLS frame, mirrored orientation not needed for c)
    g0 = geom(R, D, Y)
    log = open(f'run_{tag}.log', 'w')
    log.write(json.dumps(cfg) + '\n'); log.write(f"seed k={S.kmax} lam={S.lam} ell={S.ell} n={R.n} closure turn+pi={g0['turn']+np.pi:.3e} ImG(ell)={g0['G_end'].imag:.3e}\n"); log.flush()
    sc = np.ones(len(Y)); sc[R.n:] = np.maximum(np.abs(Y[R.n:]), 1.0)
    slabs = []; tau = 0.0; Zg = None; t0 = time.time(); k = S.kmax; ks = 0
    while k < KEND:
        ts = time.time()
        F0 = D.F(Y); J = D.jac(Y, F0); tj = time.time() - ts
        RD.factor(J, dt)
        Ys, hist = RD.solve(D.F, Y, dt, Z0=Zg, scale=sc, tol=2e-11)
        nodes = np.vstack([Y, Ys])                                   # (s+1, ntot) at relative times (0, c)
        slabs.append(dict(tau0=tau, dt=dt, nodes=nodes))
        # extrapolated guess for the next slab
        W = RD.interp_matrix(1.0 + RD.c)
        Zg = W @ nodes - Ys[-1]
        Y = Ys[-1].copy(); tau += dt; ks += 1
        lam = np.exp(Y[R.n + 1]); k = K0 / lam
        g = geom(R, D, Y); p = D.parts(Y)
        log.write(f"slab {ks:4d} tau={tau:8.4f} k={k:10.3f} a={p['a']:+.5e} ell={Y[R.n]:.6e} it={len(hist):2d} last={hist[-1]:.1e} "
                  f"turn+pi={g['turn']+np.pi:+.2e} ImG={g['G_end'].imag:+.2e} wall={time.time()-ts:.1f}s (jac {tj:.1f}s)\n"); log.flush()
        if ks % 10 == 0 or k >= KEND:
            np.savez(f'slabs_{tag}.npz', cfg=json.dumps(cfg), tau0=np.array([q['tau0'] for q in slabs]), dt=np.array([q['dt'] for q in slabs]),
                     nodes=np.array([q['nodes'] for q in slabs]), c=RD.c)
    log.write(f"DONE {ks} slabs, wall {time.time()-t0:.0f}s\n"); log.close()

if __name__ == '__main__':
    main()
