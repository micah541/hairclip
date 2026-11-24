"""P2: relative energy-form norm of  S_b - S_a  (sym_rho forms of two snapshots, tip frame lam = 1.5398/kappa_tip, operator A: S1-frame terms + transport tau)
on the rho-orthogonal complement of ncon S1 eigenvectors:  theta_ab = max |eig( S_a - S_b , -S_a )|  (generalised, restricted).
Certificate at s gives -S_a >= w_a on the complement.  Then -S_b >= (1-theta) w_a  =>  loss = theta*w_a <= margin  <=>  theta <= m/(omega+m).
usage: python3 p2.py out.json files...  (files ordered along the trajectory)"""
import sys, json, glob, numpy as np, scipy.linalg as sla
sys.dont_write_bytecode = True
sys.path.insert(0, '/home/claude/s1spec/step4/feas2')
from f3geom import *; from f3fem import *
OM, MARG = 0.05, 0.02; TT = 4.598863926611
out = sys.argv[1]; files = sys.argv[2:]
mesh = Mesh(35.0, 0.1); M = mesh.mass(); sp1 = s1_splines(Smax=40, dx=0.01)
w1, U1 = eig_sym(mesh.form(coeffs(mesh, sp1, use_tau=False)), M)
def sym(K): return 0.5 * (K + K.T)
snaps = []
for f in files:
    z = np.load(f); t = float(z['t']); km = float(z['kmax'])
    dsig = 0.01 if km < 300 else 0.02
    g = numeric_geometry(f, Smax=40, dsig=dsig); sp = as_splines(g)
    K = mesh.form(coeffs(mesh, sp, use_tau=True)); S = sym(K)
    s = -np.log(TT - t) if km < 5.2 else 4 * np.log(km / K0S1)
    snaps.append(dict(f=f, t=t, km=km, s=float(s), S=S, lam=float(g['lam'])))
    print(f"loaded {f.split('/')[-1]} t={t:.10f} kmax={km:.3f} s={s:.3f}", flush=True)
def restr(ncon):
    C = U1[:, :ncon]; Z = sla.null_space(C.T @ M); return Z
Zs = {3: restr(3), 4: restr(4)}
def coercive(S, ncon):
    Z = Zs[ncon]; A = -(Z.T @ S @ Z); Mr = Z.T @ M @ Z
    return float(sla.eigh(A, Mr, eigvals_only=True)[0])
res = []
for a, b in zip(snaps[:-1], snaps[1:]):
    row = dict(f_a=a['f'].split('/')[-1], f_b=b['f'].split('/')[-1], t_a=a['t'], t_b=b['t'], k_a=a['km'], k_b=b['km'], s_a=a['s'], s_b=b['s'], ds=b['s'] - a['s'])
    for ncon in (3, 4):
        Z = Zs[ncon]; Mr = Z.T @ M @ Z; Ar = -(Z.T @ a['S'] @ Z); Dr = Z.T @ (a['S'] - b['S']) @ Z
        lam0 = float(sla.eigh(Ar, Mr, eigvals_only=True)[0]); row[f'coer_a_nc{ncon}'] = lam0
        # energy norm: A_r if positive definite, else A_r + shift*Mr (shift=1; flagged)
        shift = 0.0 if lam0 > 0.02 else (0.05 - lam0)
        th = float(np.max(np.abs(sla.eigh(Dr, Ar + shift * Mr, eigvals_only=True)))); row[f'theta_nc{ncon}'] = th; row[f'shift_nc{ncon}'] = shift
    res.append(row)
    print(f"{row['t_a']:.6f}->{row['t_b']:.6f} k {row['k_a']:.2f}->{row['k_b']:.2f} ds={row['ds']:.4f} coer3={row['coer_a_nc3']:.3f} theta3={row['theta_nc3']:.3e} | coer4={row['coer_a_nc4']:.3f} theta4={row['theta_nc4']:.3e}", flush=True)
json.dump(res, open(out, 'w'), indent=1)
