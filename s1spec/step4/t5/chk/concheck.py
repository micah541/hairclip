"""T5 float check: top of sym_rho(L_Gamma + transport) on the rho-complement of base-adapted constraint sets.
sets: SYM = first 3 S1-sym eigenvectors (FEAS reference);  GEO = {p, n1, third}, third = kappa'' (numeric) or psi3 (S1: direct eigenfunction, ev -0.256);
      ADJ = {p, n1, adjoint eigenfunction of the Galerkin operator (ev nearest -0.256)}.   single process, float only."""
import sys, numpy as np, scipy.linalg as sla
sys.dont_write_bytecode = True
sys.path.insert(0, '/home/claude/s1spec/step4/feas')
from f3geom import *; from f3fem import *
Sigma, h = 40.0, 0.1
mesh = Mesh(Sigma, h); M = mesh.mass(); xn = np.arange(mesh.ne + 1) * mesh.h
def cut(x, a=22.0, b=30.0):
    t = np.clip((x - a) / (b - a), 0, 1); return 1 - t**3 * (10 - 15 * t + 6 * t * t)
def fem_vec(x, f):
    sp = CubicSpline(x, f * cut(x)); d = np.zeros(mesh.ndof); d[0::2] = sp(xn); d[1::2] = sp(xn, 1); return d[mesh.free]
def top_on_complement(S, C):
    Z = sla.null_space(C.T @ M); w = sla.eigh(Z.T @ S @ Z, Z.T @ M @ Z, eigvals_only=True); return w[::-1][:3]
def gs(C):
    Q = []
    for c in C.T:
        for q in Q: c = c - (q @ M @ c) * q
        Q.append(c / np.sqrt(c @ M @ c))
    return np.array(Q).T
def run(label, g, use_tau):
    sp = as_splines(g); K = mesh.form(coeffs(mesh, sp, use_tau=use_tau)); S = 0.5 * (K + K.T)
    x = g['x']; p = g['p']
    n1 = g['n1']
    out = {}
    out['SYM'] = top_on_complement(S, U1[:, :3])
    # direct / adjoint eigenvectors of the (nonsymmetric) Galerkin operator
    w, vr = eig_nonsym(K, M); wl, vl = eig_nonsym(K.T, M)
    i3 = np.argmin(np.abs(w - w.real[2])) if True else 2
    ev = w.real[:4]
    psi3 = vr[:, 2].real; adj3 = np.linalg.solve(M, vl[:, 2].real) if False else vl[:, 2].real
    # adjoint w.r.t. <.,.>_M : K^T y = lam M y  -> y is the coefficient vector of the adjoint eigenfunction
    P, N1 = fem_vec(x, p), fem_vec(x, n1)
    third = psi3 if g.get('third') is None else fem_vec(x, g['third'])
    C = gs(np.c_[P, N1, third]); out['GEO'] = top_on_complement(S, C)
    ang = [float(abs(C[:, 2] @ M @ psi3) / np.sqrt(psi3 @ M @ psi3))]
    C2 = gs(np.c_[P, N1, adj3]); out['ADJ'] = top_on_complement(S, C2)
    C3 = gs(np.c_[P, N1, psi3]); out['PSI'] = top_on_complement(S, C3)
    # two-constraint only (p, n1): how bad is the third direction
    out['2only'] = top_on_complement(S, gs(np.c_[P, N1]))
    # oblique-projection data: rank of Gram between GEO and SYM
    print(f"{label:10s} ev(L)={np.round(ev,4)}  SYM {np.round(out['SYM'],3)}  GEO {np.round(out['GEO'],3)}  PSI {np.round(out['PSI'],3)}  ADJ {np.round(out['ADJ'],3)}  2only {np.round(out['2only'],3)}  cos(third_perp,psi3)={ang[0]:.3f}", flush=True)
sp1 = s1_splines(Smax=Sigma + 5, dx=0.01); K1 = mesh.form(coeffs(mesh, sp1, use_tau=False)); ws1, U1 = eig_sym(K1, M)
x1 = np.arange(0, Sigma + 5, 0.01); g1 = s1_geometry(x1)
# S1 normal x-component: n1 = dp/dx-type: use q' = 1 + k p, p' = -k q ; translation mode psi2 = <e1,nu>; with theta: T=(cos th, sin th), th' = k
th = CubicSpline(x1, g1['kap']).antiderivative()(x1)   # tangent angle from the tip (tip tangent vertical up to a constant): n1 = <e1,N>
# at the tip of S1 the tangent is perpendicular to the axis: T(0) = (0,1)*sign -> T = (-sin th, cos th), N = iT = (-cos th, -sin th) ; n1 = -cos th (sign irrelevant)
g1['n1'] = -np.cos(th); g1['third'] = None
run('S1', g1, False)
for nm in sys.argv[1:]:
    g = numeric_geometry(f'/home/claude/s1spec/step4/feas2/snaps/p2_{nm}.npz', Smax=Sigma + 5, dsig=0.01)
    g['n1'] = (np.conj(1j * g['T'])).real; g['third'] = g['kap2']
    run(nm, g, True)
