import numpy as np
from rep import *; from ops import *; from seed import HLSHalf; from run import seed_closed, geom
S = HLSHalf('../feas2/snaps/p2_k20.npz')
R = Rep(9, Mesh()); D = Dyn(R)
ell = S.ell
x, dx = R.mesh.knots(ell); xp, _ = R.mesh.knots(ell * (1 + 1e-6)); xm, _ = R.mesh.knots(ell * (1 - 1e-6))
print('dx/dell check: max |FD - analytic| =', np.max(np.abs((xp - xm) / (2e-6 * ell) - dx)), ' max dx', dx.max())
y = seed_closed(R, D, S); Y = np.r_[y, ell, np.log(S.lam), 0.0, S.t]
Fv = D.F(Y)
# d/dtau int kappa via finite difference along F
for e_ in (1e-6, 1e-7):
    gp = geom(R, D, Y + e_ * Fv); gm = geom(R, D, Y - e_ * Fv)
    print('d/dtau turn (FD):', (gp['turn'] - gm['turn']) / (2 * e_), '  d/dtau Im G(ell):', (gp['G_end'].imag - gm['G_end'].imag) / (2 * e_))
# moments residual P0 sum check: sum_e int_e G dx + kappa(ell) * elld
p = D.parts(Y); B = R.B
intG = np.sum(p['h'] * (p['G'] @ B.wq)); kl = p['k0'][-1, -1]
print('int G =', intG, ' -kappa(ell)*elld =', -R.eval_at(y, ell, np.array([ell * (1 - 1e-12)]))[0] * p['elld'])
h = p['h']; I = lambda f: np.sum(h * (f @ B.wq))
k0, k1, k2, k4, Q, xq = p['k0'], p['k1'], p['k2'], p['k4'], p['Q'], p['xq']
D_ = R.nodal(y, h)
print('nodal at ends: k\'(0), k\'\'\'(0), k\'(l), k\'\'\'(l):', D_[0, 1], D_[0, 3], D_[-1, 1], D_[-1, 3])
print('int k4 =', I(k4), '  k3 ends', D_[-1, 3] - D_[0, 3])
print('int (k0^2 k2 + k0 k1^2) =', I(k0**2 * k2 + k0 * k1**2), ' vs [k^2 k\'] + int k k\'^2:', I(k0 * k1**2))
print('int k1 Q =', I(k1 * Q), ' vs kQ(l) - int k k1^2 =', D_[-1, 0] * p['Qe'][-1] - I(k0 * k1**2))
print('int (k0 + x k1) =', I(k0 + xq * k1), ' vs ell*k(ell) =', ell * D_[-1, 0])
# continuity check of kappa''' across nodes from element polynomials
k3 = p['d'] @ B.E1[3] / h**3; k3l = p['d'] @ B.E0[3] / h**3
print('max jump kappa\'\'\' at interior nodes:', np.max(np.abs(k3[:-1] - k3l[1:])), ' kappa\'\'\' at last node from element:', k3[-1])
for kk in range(5):
    vr = p['d'] @ B.E1[kk] / h**kk; vl = p['d'] @ B.E0[kk] / h**kk
    print(kk, 'max jump', np.max(np.abs(vr[:-1] - vl[1:])), ' max |val|', np.max(np.abs(vr)))
