import numpy as np
from base2 import Base; from arbops import SlabPolish; import resid; from handoff import Curve
from rep import DS1
B = Base('C11'); i = 12; R, D = B.rep(i); q = B.slabs[i]
SP = SlabPolish(R, D, B.RD, q['nodes'], q['dt']); sm = resid.samples(SP, 0.0)
X, H = sm['X'], sm['H']
th = -np.pi / 2 + resid.cumint(sm['k0'], H); G = -DS1 + resid.cumint(np.exp(1j * th), H)
C = Curve(R, q['nodes'][0]); m = X <= 30
print('theta: resid-cumint vs Curve: max diff', np.abs(th[m] - C.theta(X[m])).max())
print('Gamma: resid-cumint vs Curve: max diff', np.abs(G[m] - C.Gam(X[m])).max())
print('kappa: arb sample vs float eval: max diff', np.abs(sm['k0'][m] - C.kap(X[m])).max())
print('kappa_xx: arb vs float:', np.abs(sm['k2'][m] - C.kap(X[m], 2)).max())
print('S1 check: kappa_xx - p/4 at x=0 (base not S1):', sm['k2'][0, 0], -DS1 / 4)
tau0 = q['tau0']; dlt = 1e-3
def thx(t): return Curve(R, B.Y(t, i)).theta(X[m])
tfd = (8 * (thx(tau0 + dlt) - thx(tau0 - dlt)) - (thx(tau0 + 2 * dlt) - thx(tau0 - 2 * dlt))) / (12 * dlt)
print('theta_tau: arb tt vs FD(4th order, d=1e-3): max diff', np.abs(sm['tt'][m] - tfd).max(), ' max|tt|', np.abs(tfd).max())
dlt = 2e-3; tfd2 = (8 * (thx(tau0 + dlt) - thx(tau0 - dlt)) - (thx(tau0 + 2 * dlt) - thx(tau0 - 2 * dlt))) / (12 * dlt)
print('   FD consistency d=1e-3 vs 2e-3:', np.abs(tfd - tfd2).max())
xx = X[m]; j = np.argsort(np.abs(sm['tt'][m] - tfd))[-5:]
print('   worst at x =', xx[j], (sm['tt'][m] - tfd)[j])
