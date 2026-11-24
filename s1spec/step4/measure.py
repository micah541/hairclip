"""Rescaled deviation of the hairpin's collapsing pocket from S1.
For a snapshot: Gamma = (gamma - x_c)/lambda, normal graph u over S1's half-branch, modulation (x_c, lambda) fixed
by <u, f1>_rho = <u, f2>_rho = 0 (f1 = chi p, f2 = -chi sin theta1, chi the C^3 cutoff on sigma in [29,32],
additionally truncated where the graph chart ends).  Reports ||u||_rho, sup|u|, and the rescaled positions of c1, c2."""
import sys, glob, numpy as np
sys.path.insert(0, '/home/claude/s1spec')
from evans import PROF, STIP
from hls import HLS

# ---------------- S1 half-branch geometry, x = sigma - sigma_tip in [0, XS]
XS, HX = 40.0, 0.01
xs = np.arange(0, XS + HX / 2, HX)
sig = xs + STIP
P = np.zeros((4, len(xs))); ins = sig <= 22
P[:, ins] = PROF(sig[ins]); P[1, ~ins] = sig[~ins]
p1, q1, k1, m1 = P
th1 = -np.pi / 2 + np.concatenate([[0], np.cumsum((k1[1:] + k1[:-1]) / 2) * HX])
tau1 = np.exp(1j * th1); nu1 = 1j * tau1
G1 = q1 * tau1 + p1 * nu1                                   # Gamma1 = q tau1 + p nu1
rho = np.exp(-0.25 * (1 + xs**2) ** (2 / 3))
wq = np.full(len(xs), HX); wq[0] = wq[-1] = HX / 2; wq = 2 * wq * rho
uu = np.clip((32 - sig) / 3, 0, 1); chi = uu**4 * (35 - 84 * uu + 70 * uu**2 - 20 * uu**3)
F1, F2 = chi * p1, -chi * np.sin(th1)
D_S1 = -G1[0].real; K0 = k1[0]


def load(path):
    z = np.load(path)
    S = HLS(np.zeros(int(z['N'])), float(z['L']), z['X0'], omega=1)
    S.th_hat = z['th_hat']; S.t = float(z['t'])
    return S, z


def strand(S, up=8):
    """pocket-tip strand: alpha decreasing from pi (S1's half-branch direction), as an ordered polyline."""
    al, th, Z, clo = S.curve(up)
    N = len(al); i_tip = N // 2
    return Z[i_tip::-1], al[i_tip::-1]          # from the tip backwards to alpha = 0


def graph(Zs, xc, lam, xmax):
    """u(x) for x in [0, xmax] (S1 grid), by intersecting S1 normals with the strand. Returns u, valid mask."""
    n = int(round(xmax / HX)) + 1
    u = np.full(n, np.nan)
    Pp = xc + lam * G1[:n]
    j = 0; M = len(Zs)
    u[0] = ((Zs[0] - Pp[0]) * np.conj(nu1[0])).real / lam       # the tip itself (on the axis)
    for i in range(1, n):
        t, nv = tau1[i], nu1[i]
        # advance j until the tangential coordinate of the strand point relative to P changes sign
        f = lambda jj: ((Zs[jj] - Pp[i]) * np.conj(t)).real
        jj = max(j - 5, 0)
        fj = f(jj)
        found = False
        while jj < M - 1:
            fn = f(jj + 1)
            if abs(fj) < 1e-7 * lam:
                fn = fj; fj = 0.0
            if fj <= 0 < fn or fj >= 0 > fn or fj == 0.0:
                a = fj / (fj - fn) if fj != fn else 0.0
                zz = Zs[jj] + a * (Zs[jj + 1] - Zs[jj])
                s = ((zz - Pp[i]) * np.conj(nv)).real
                if abs(s) < 0.5 * lam * (1 + xs[i]):  # plausibility: within half a rescaled unit (scaled)
                    u[i] = s / lam; j = jj; found = True
                break
            jj += 1; fj = fn
            if jj - j > 4000: break
        if not found: break
    return u


def residual(Zs, xc, lam, xmax):
    u = graph(Zs, xc, lam, xmax)
    ok = np.isfinite(u); n = len(u)
    w = wq[:n] * ok; uz = np.where(ok, u, 0.0)
    return np.array([np.sum(w * uz * F1[:n]), np.sum(w * uz * F2[:n])]), u, ok


def modulate(Zs, xc0, lam0, xmax, it=12):
    v = np.array([xc0, lam0])
    for _ in range(it):
        r, u, ok = residual(Zs, v[0], v[1], xmax)
        Jm = np.zeros((2, 2))
        for k, d in enumerate([1e-6 * v[1], 1e-6 * v[1]]):
            e = np.zeros(2); e[k] = d
            rp, _, _ = residual(Zs, *(v + e), xmax)
            Jm[:, k] = (rp - r) / d
        dv = np.linalg.solve(Jm, -r)
        v = v + dv
        if np.max(np.abs(dv)) < 1e-12 * v[1]: break
    r, u, ok = residual(Zs, v[0], v[1], xmax)
    return v, u, ok, r


def analyse(path, xmax=34.0):
    S, z = load(path)
    Zs, als = strand(S)
    tip = Zs[0]
    km_tip = abs(S.kmax())
    # tip curvature from theta at alpha = pi
    al, th = S.theta(8); N = len(al)
    kt = (2 * np.pi / S.L) * np.gradient(th, al)[N // 2]
    lam0 = K0 / abs(kt); xc0 = tip.real + D_S1 * lam0
    # chart extent: stop where the strand passes c1 (rescaled) or at xmax
    c1r = (float(z['c1']) - xc0) / lam0
    v, u, ok, r = modulate(Zs, xc0, lam0, xmax)
    xc, lam = v
    n = len(u); xx = xs[:n]
    good = ok
    nrm = np.sqrt(np.sum(wq[:n][good] * u[good]**2))
    sup10 = np.nanmax(np.abs(u[(xx <= 10) & good])) if np.any((xx <= 10) & good) else np.nan
    xend = xx[good][-1] if np.any(good) else 0.0
    return dict(t=S.t, kmax=float(z['kmax']), ktip=abs(kt), lam=lam, lam0=lam0, xc=xc, xc0=xc0, nrm=nrm, sup10=sup10,
                xend=xend, c1r=(float(z['c1']) - xc) / lam, c2r=(float(z['c2']) - xc) / lam, AX=abs(float(z['AX'])),
                AL=abs(float(z['AL'])), u=u, res=r)


if __name__ == "__main__":
    tag = sys.argv[1]
    files = sorted(glob.glob(f"snap_{tag}_*.npz"), key=lambda f: int(f.split('_')[-1][:-4]))
    out = []
    for f in files:
        try:
            d = analyse(f)
        except Exception as e:
            print(f, "failed:", e); continue
        out.append(d)
        print(f"kmax={d['kmax']:8.1f} t={d['t']:.10f} lam={d['lam']:.4e} (tip-fit {d['lam0']:.4e}) xc={d['xc']:.6f} "
              f"|u|_rho={d['nrm']:.3e} sup|u|(x<=10)={d['sup10']:.3e} chart to x={d['xend']:.1f}  "
              f"c1 at {d['c1r']:.2f}, c2 at {d['c2r']:.3f} (S1: {-0.8302:.3f})", flush=True)
    np.save(f"measure_{tag}.npy", out, allow_pickle=True)
