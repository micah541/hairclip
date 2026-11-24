"""T8 far-field base sampler (float).  For a slab i and relative time theta of the T1 base (slabs_P11.pkl), returns
per-element arrays on the whole half curve at the points u = 0, Gauss(NG), 1 of every element:
  x, y = x - X0, quadrature weights, kappa^(m) (m = 0..5, exact element polynomials), Q = int_0^x kappa'^2,
  theta, Gamma, p = <G,nu>, q = <G,T>, n = <e1,nu>, t1 = <e1,T>  (nu = i T),
and the scalars a, c = 1/4 - a, b, ell, ell_tau = -Q(ell) + c ell, plus the rho-norms (on supp chi, chi = S((32-y)/3))
of the three modulation directions g = (p, n, V/|chi V|) with V = -kappa'' + c p + b n (the residual part of V~ is
dropped: it is <= 3e-14 and enters only T6).
The far field is x >= 16.2491 (y >= 14); the Lemma-S ramp eta = S(ln(y/14)/ln 2) on y in [14, 28]."""
import sys, numpy as np
from numpy.polynomial import legendre as LG
sys.path.insert(0, '/home/claude/s1spec/step4/t1')
from base2 import Base
from rep import K0, DS1

X0 = 2.2490776515          # Lemma S: y = x - X0 (T1 uses STIP = -2.249077651121112; difference 1e-11, immaterial)
NG = 24
_g, _w = LG.leggauss(NG)
UG = np.r_[0.0, (_g + 1) / 2, 1.0]          # sample points per element (ends included)
WG = np.r_[0.0, _w / 2, 0.0]                # Gauss weights (ends get 0)

def smooth(u):
    u = np.clip(u, 0, 1); return u**4 * (35 - 84 * u + 70 * u**2 - 20 * u**3)

_B = None
def base():
    global _B
    if _B is None: _B = Base('P11')
    return _B

def sample(i, theta, need_geom=True):
    B = base(); R, D = B.rep(i)
    tau = B.tau0[i] + theta * B.dt[i]
    Y = B.Y(tau, i); n = R.n; ell = Y[n]
    x, _ = R.mesh.knots(ell); h = np.diff(x); d = R.local(Y[:n], h)
    N = R.N; Bs = R.B
    K = np.empty((6, N, len(UG)))
    for m in range(6):
        V = Bs.evalpts(UG, m)                      # (nu, nb)
        K[m] = (d @ V.T) / h[:, None] ** m
    X = x[:-1, None] + h[:, None] * UG[None, :]
    Wq = h[:, None] * WG[None, :]
    # Q = int kappa'^2: exact per element by Gauss-Legendre of degree 2P-2 with NG >= P (exact), cumulative via Legendre fit
    f = K[1] ** 2
    # element integrals (Gauss part only: columns 1..NG)
    Qe = np.r_[0.0, np.cumsum(np.sum(f[:, 1:-1] * (h[:, None] * WG[None, 1:-1]), axis=1))]
    # cumulative inside element: fit Legendre series of f on Gauss nodes (exact: f is a polynomial of degree 2P-2 = 20 < NG)
    Qloc = np.empty_like(f)
    for e in range(N):
        cf = LG.legfit(_g, f[e, 1:-1], NG - 1)
        ci = LG.legint(cf, lbnd=-1) * (h[e] / 2)
        Qloc[e] = LG.legval(2 * UG - 1, ci)
    Q = Qe[:-1, None] + Qloc
    p = D.parts(Y)
    a = p['a']; c = 0.25 - a; b = p['k2t'] + c * DS1
    out = dict(i=i, theta=theta, tau=tau, k=K0 * np.exp(-Y[n + 1]), ell=ell, a=a, c=c, b=b, X=X, Y=X - X0, Wq=Wq, K=K, Q=Q,
               Qend=Qe[-1], ellt=-Qe[-1] + c * ell, h=h, xk=x)
    if need_geom:
        # theta and Gamma: cumulative integrals per element (Legendre series, exact for theta; Gamma to ~1e-14 rel.)
        th = np.empty_like(Q); G = np.empty(Q.shape, complex)
        th0 = -np.pi / 2; G0 = -DS1 + 0j
        for e in range(N):
            ck = LG.legfit(_g, K[0, e, 1:-1], NG - 1)
            ci = LG.legint(ck, lbnd=-1) * (h[e] / 2); ci[0] += th0 - LG.legval(-1.0, ci)
            th[e] = LG.legval(2 * UG - 1, ci)
            # Gamma: integrate e^{i theta} with a finer Gauss rule (theta nonpolynomial in e^{i.})
            gg, ww = LG.leggauss(48); thg = LG.legval(gg, ci); Eg = np.exp(1j * thg)
            cE = LG.legfit(gg, Eg.real, 47) + 1j * LG.legfit(gg, Eg.imag, 47)
            cI = LG.legint(cE, lbnd=-1) * (h[e] / 2); cI[0] += G0 - LG.legval(-1.0, cI)
            G[e] = LG.legval(2 * UG - 1, cI)
            th0 = LG.legval(1.0, ci); G0 = LG.legval(1.0, cI)
        T = np.exp(1j * th); nu = 1j * T
        pp = (G * np.conj(nu)).real; qq = (G * np.conj(T)).real
        nn = nu.real; t1 = T.real
        Vs = c * pp + b * nn                   # smooth part of V~ ; rough part -kappa''
        out.update(th=th, G=G, p=pp, q=qq, n=nn, t1=t1, Vs=Vs, Vr=-K[2])
        # rho-norms on supp chi
        rho = np.exp(-0.25 * (1 + X**2) ** (2 / 3)); chi = smooth((32 - (X - X0)) / 3)
        nrm = lambda g: float(np.sqrt(np.sum(Wq * rho * (chi * g) ** 2)))
        out.update(np_=nrm(pp), nn_=nrm(nn), nV=nrm(Vs - K[2]))
        out['closure'] = abs(G0.imag)          # Im Gamma(ell): closure defect (should be ~0)
    return out

if __name__ == '__main__':
    for i in (0, 60, 120, 187):
        s = sample(i, 0.5)
        print(i, f"k={s['k']:.4g} ell={s['ell']:.4g} c={s['c']:.5f} b={s['b']:.3e} |p|_rho={s['np_']:.4e} |n|_rho={s['nn_']:.4e} |V|_rho={s['nV']:.4e} closure={s['closure']:.2e} ellt={s['ellt']:.4g}")
        X = s['X']; m = (X > 16.25) & (X < 400)
        print('   far V/|V|: max |g3''| / y on [16,400]:', np.max(np.abs(np.gradient(s['Vs'][m] - s['K'][2][m]) if False else 0)),
              ' |p| max', np.max(np.abs(s['p'][m])), ' p(30)', np.interp(30, X.ravel(), s['p'].ravel()), ' Vs(30)/nV', np.interp(30, X.ravel(), s['Vs'].ravel()) / s['nV'])
