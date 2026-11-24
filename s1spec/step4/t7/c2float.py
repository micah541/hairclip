"""T7 float calibration of (C2): change of the projected form between tau_c and tau_c + dtau along the T1 base.
 sharp:  s_sharp = max |w^T dB w| / (||u||^2 + kappa ||u''||^2) on {<u,phi_j> = 0}   (FEM generalised eigenvalue, kappa = 0.25)
 crude:  T5 2.6-type sup bound  s_crude = sup|db0| + 1/2 sup|(rho db1)'/rho| + sup|db2|/(2 sqrt(kappa)) + ||dPi3|| (rank-3 change, FEM norm)
(C2) needs s <= m_c/2 ~ 0.019 on each interval.  -> c2float.json"""
import json, numpy as np, scipy.linalg as sla
import bgeo, fem
B = bgeo.base(); st = fem.Setup(34.25, 0.1, 'clamped'); F = st.fem; x = F.xq
Z = sla.null_space(st.Phiv.T); kap = 0.25
Mz = Z.T @ (st.M + kap * st.K2) @ Z
def sym(A): return 0.5 * (A + A.T)
out = []
for kk in (20.4, 40, 200, 2e3, 2e4):
    if kk == 20.4: i, th = 0, 0.5
    else: tau, i = B.tau_of_k(kk); th = (tau - B.tau0[i]) / B.dt[i]
    tau_c = B.tau0[i] + th * B.dt[i]
    g0 = bgeo.BaseGeom(i, th); op0 = fem.operator(st, g0); B0, _ = fem.projected(st, op0)
    Pi0 = op0['Gv'] @ np.linalg.inv(op0['M3']).T @ op0['Lv'].T
    for dt in (0.15, 0.3, 0.6, 1.2, 2.4):
        t1 = tau_c + dt
        if t1 > B.tau0[-1] + B.dt[-1]: continue
        j = B.slab_of(t1); g1 = bgeo.BaseGeom(j, (t1 - B.tau0[j]) / B.dt[j]); op1 = fem.operator(st, g1); B1, _ = fem.projected(st, op1)
        dB = Z.T @ sym(B1 - B0) @ Z; w = sla.eigh(dB, Mz, eigvals_only=True); s_sharp = float(np.max(np.abs(w)))
        Pi1 = op1['Gv'] @ np.linalg.inv(op1['M3']).T @ op1['Lv'].T
        wP = sla.eigh(Z.T @ sym(Pi1 - Pi0) @ Z, Z.T @ st.M @ Z, eigvals_only=True); dPi = float(np.max(np.abs(wP)))
        db = [b1 - b0 for b0, b1 in zip(op0['b'], op1['b'])]          # (db2, db1, db0) at quadrature points
        r1 = F.rho1 / F.rho
        db1p = np.gradient(db[1].ravel(), x.ravel()).reshape(x.shape)
        s_crude = float(np.max(np.abs(db[2])) + 0.5 * np.max(np.abs(db1p + r1 * db[1])) + np.max(np.abs(db[0])) / (2 * np.sqrt(kap)) + dPi)
        out.append(dict(k=g0.k, dtau=dt, s_sharp=s_sharp, s_crude=s_crude, dPi=dPi, db1_far=float(np.max(np.abs(db[1][x > 25]))), dc=float(g1.c - g0.c)))
        print(f"k={g0.k:8.1f} dtau={dt:4.2f}: s_sharp={s_sharp:.4f}  s_crude={s_crude:.4f}  (rank-3 change {dPi:.4f}, sup|db1| x>25 {out[-1]['db1_far']:.3f})", flush=True)
json.dump(out, open('c2float.json', 'w'), indent=1)
