"""T5 float check 2: projected dynamics e_tau = P L e, P = I - G (Phi^T M G)^-1 Phi^T M  (constraints Phi fixed, directions G = base modes {p, n1, kappa''}).
top eigenvalue of sym of the projected form  B(e,e) = <Le,e> - <e,G>(Phi^T M G)^-1 <Phi, L e>  on Phi-perp, for Phi = S1-sym eigvecs (SYM), Phi = S1 {p,n1,psi3} (S1GEO), Phi = G (GEO)."""
import sys, numpy as np, scipy.linalg as sla
sys.dont_write_bytecode = True
sys.path.insert(0, '/home/claude/s1spec/step4/feas'); sys.path.insert(0, '/home/claude/s1spec/step4/t5/chk')
from f3geom import *; from f3fem import *
Sigma, h = 40.0, 0.1
mesh = Mesh(Sigma, h); M = mesh.mass(); xn = np.arange(mesh.ne + 1) * mesh.h
def cut(x, a=22.0, b=30.0):
    t = np.clip((x - a) / (b - a), 0, 1); return 1 - t**3 * (10 - 15 * t + 6 * t * t)
def fem_vec(x, f):
    sp = CubicSpline(x, f * cut(x)); d = np.zeros(mesh.ndof); d[0::2] = sp(xn); d[1::2] = sp(xn, 1); return d[mesh.free]
def gs(C):
    Q = []
    for c in C.T:
        for q in Q: c = c - (q @ M @ c) * q
        Q.append(c / np.sqrt(c @ M @ c))
    return np.array(Q).T
def proj_top(K, Phi, G):
    Z = sla.null_space(Phi.T @ M)
    B = K - M @ G @ np.linalg.solve(Phi.T @ M @ G, Phi.T @ K)
    Bz = Z.T @ B @ Z; w = sla.eigh(0.5 * (Bz + Bz.T), Z.T @ M @ Z, eigvals_only=True)
    cond = np.linalg.cond(Phi.T @ M @ G)
    return w[::-1][:2], cond
sp1 = s1_splines(Smax=Sigma + 5, dx=0.01); K1 = mesh.form(coeffs(mesh, sp1, use_tau=False)); ws1, U1 = eig_sym(K1, M)
x1 = np.arange(0, Sigma + 5, 0.01); g1 = s1_geometry(x1)
th = CubicSpline(x1, g1['kap']).antiderivative()(x1)
w, vr = eig_nonsym(K1, M); psi3 = vr[:, 2].real
PhiS1 = gs(np.c_[fem_vec(x1, g1['p']), fem_vec(x1, -np.cos(th)), psi3]); PhiSYM = U1[:, :3]
print('S1: GEO', np.round(proj_top(K1, PhiS1, PhiS1)[0], 3), ' SYM/dir S1GEO', np.round(proj_top(K1, PhiSYM, PhiS1)[0], 3))
for nm in sys.argv[1:]:
    g = numeric_geometry(f'/home/claude/s1spec/step4/feas2/snaps/p2_{nm}.npz', Smax=Sigma + 5, dsig=0.01)
    K = mesh.form(coeffs(mesh, as_splines(g), use_tau=True))
    G = gs(np.c_[fem_vec(g['x'], g['p']), fem_vec(g['x'], (np.conj(1j * g['T'])).real), fem_vec(g['x'], g['kap2'])])
    a, ca = proj_top(K, G, G); b, cb = proj_top(K, PhiSYM, G); c, cc = proj_top(K, PhiS1, G)
    # distance of G from PhiS1 (principal angles)
    sv = np.linalg.svd(PhiS1.T @ M @ G, compute_uv=False)
    print(f"{nm:6s} Phi=G {np.round(a,3)} | Phi=SYM {np.round(b,3)} cond {cb:.1f} | Phi=S1GEO {np.round(c,3)} cond {cc:.2f} | cos angles(G,S1GEO) {np.round(sv,3)}", flush=True)
