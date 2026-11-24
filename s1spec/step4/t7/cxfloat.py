"""T7 float: modulation constant C_x = sup |x| / ||u||_rho, x = M^-T l(u) (normalised directions g_i/||chi g_i||_rho), along the base
-> cx.json  (also ||M^-1|| and lbar = (sum ||l_j||^2)^(1/2) with phi_j normalised in rho)"""
import json, numpy as np, scipy.linalg as sla
import bgeo, fem
B = bgeo.base(); st = fem.Setup(34.25, 0.1, 'clamped'); F = st.fem
Mh = sla.sqrtm(st.M).real; Mhi = np.linalg.inv(Mh)
out = []
for kk, i, th in [(20.4, 0, 0.5), (40, None, None), (200, None, None), (2e3, None, None), (2e4, None, None)]:
    if i is None: tau, i = B.tau_of_k(kk); th = (tau - B.tau0[i]) / B.dt[i]
    g = bgeo.BaseGeom(i, th); op = fem.operator(st, g)
    ng = np.array([np.sqrt(F.ip(st.chi * gi, st.chi * gi)) for gi in op['g']])
    nphi = np.array([np.sqrt(F.ip(f[0], f[0])) for f in st.phi])
    M3 = op['M3'] / ng[:, None] / nphi[None, :]; Lv = op['Lv'] / nphi[None, :]
    T = np.linalg.inv(M3).T @ Lv.T @ Mhi                       # x = T (M^1/2 c)
    Cx = float(np.linalg.norm(T, 2)); lbar = float(np.linalg.norm(Lv.T @ Mhi, 'fro')); mbar = float(np.linalg.norm(np.linalg.inv(M3), 2))
    Tu = np.linalg.inv(op['M3']).T @ op['Lv'].T @ Mhi           # x for the unnormalised g = (p, n, Vs)
    out.append(dict(k=g.k, Cx=Cx, lbar=lbar, mbar=mbar, mbar_lbar=mbar * lbar, TTt=(Tu @ Tu.T).tolist()))
    print(f"k={g.k:9.1f}: C_x = {Cx:.3f}   ||M^-1|| = {mbar:.4f}  lbar = {lbar:.3f}  mbar*lbar = {mbar*lbar:.3f}", flush=True)
json.dump(out, open('cx.json', 'w'), indent=1)
