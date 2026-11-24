import sys, numpy as np
from base2 import Base; import resid
from polish_chain import polish_chain
B = Base('P11'); i0, nc = int(sys.argv[1]), int(sys.argv[2])
SPs = polish_chain(B, i0, nc, verbose=False); SP = SPs[-1]
o, f = resid.norms(resid.samples(SP, 0.5)); X = f['X']; Wq = f['Wq']; eta = f['eta']
print('slab', i0 + nc - 1, 'k', B.k(B.tau0[i0 + nc - 1]), {k: f'{o[k]:.2e}' for k in ('RF2', 'RF3', 'RF4')}, ' ell', X.max())
for j, d in ((3, f['r2']), (4, f['r3'])):
    I = 28 * Wq * eta**2 * X**(2 * j - 2) * d**2
    Ie = I.sum(1); tot = Ie.sum(); cs = np.cumsum(Ie) / tot
    e50, e90 = np.searchsorted(cs, 0.5), np.searchsorted(cs, 0.9)
    em = np.argmax(Ie)
    print(f'RF{j}: 50% of the integral by x = {X[e50,0]:.0f}, 90% by x = {X[e90,0]:.0f}; max element at x = {X[em,0]:.0f} (h = {SP.R.mesh.knots(X.max())[0][em+1]-SP.R.mesh.knots(X.max())[0][em]:.1f})')
    print('   |d| max in [100,ell]:', np.abs(d[X > 100]).max(), ' |rho_kappa| max there', np.abs(f['rk'][X > 100]).max(), ' |rho_kappa_x|', np.abs(f['rk1'][X > 100]).max())
for xa, xb in ((0, 40), (40, 200), (200, 1000), (1000, 1e9)):
    m = (X >= xa) & (X < xb)
    if m.any(): print(f'  x in [{xa},{xb}): max|r| {np.abs(f["r"][m]).max():.2e} max|r_xx| {np.abs(f["r2"][m]).max():.2e} max|r_xxx| {np.abs(f["r3"][m]).max():.2e} max|rho_k| {np.abs(f["rk"][m]).max():.2e}')
