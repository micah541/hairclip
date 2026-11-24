"""fix11 review, part D (3): float check of the beta decomposition and of the E0 pairing identity on random data [num].
Finite-dimensional model with the SAME algebraic structure as the tip-zone equation (functions = vectors in R^n, <f,g> = sum rho_i f_i g_i,
L, D = arbitrary matrices, chi/q/k/cos = diagonal multipliers, F = arbitrary smooth nonlinear map with DF(0) = L):
    u = u0 + w, u0 = alpha psi3 + alpha^2 h2,  u_tau = F(u) - a (psi1 - A1(u)) - b (psi2 - B1(u)),  alpha_tau = lamt alpha + beta,
    A1(u) = -(u - q Du/ag(u)),  B1(u) = Du cos/ag(u),  ag(u) = 1 - k u,   f1 = chi psi1, f2 = chi psi2, f3 arbitrary,  constraints <chi w, f_j> = 0.
What is tested (all to rounding): (T1) x = H R - H E^T x with the E, R of the lemma; (T2) beta = beta_3 + beta_12 + beta_Lin + beta_Q + beta_F + beta_S + beta_E;
(T3) (1/2) d/dtau ||chi w||^2 = <L w~,w~> - <[L,chi]w,w~> + M + <chi Phi_F,w~> + <chi N(w),w~> with M as in E0_LEMMA_v4b.tex (no psi1, psi2 terms);
(T4) <g0,w~> = <Pi g0,w~>;  (T5) oblique identity with <g0,h> = 1;  (T6) rank-one bound of the beta_3 term;  (T7) E_13 = alpha e13 + O(alpha^2), E_23 = alpha e23 + O(alpha^2),
E_33 = 2 alpha <chi h2, f3> exactly, with e13 = <gA,f3>, e23 = -<gB,f3> (signs of BLOOP.md);  (T8) the Q2 lead identity of the E2 block: -a0<f1,Dr w~> - b0<f2,Dr w~> = -a0<Dr f1,w~> - b0<Dr f2,w~> for a rho-symmetric Dr.
This checks the hand identification of the pieces and their signs; it does NOT check the certified constants or the actual profile functions.
Run: PYTHONHASHSEED=0 nice -n 19 python3 -B d_beta_float.py"""
import sys, os, json, hashlib
sys.dont_write_bytecode = True
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
res = dict(trials=[]); worst = {}
def rec(name, val): worst[name] = max(worst.get(name, 0.0), float(val))
for seed in range(40):
    rng = np.random.default_rng(seed); n = 24
    rho = np.exp(-np.linspace(0, 3, n)) * (1 + 0.1 * rng.random(n)); ip = lambda f, g: float(np.sum(rho * f * g)); nrm = lambda f: np.sqrt(ip(f, f))
    L = rng.normal(size=(n, n)); Dm = rng.normal(size=(n, n)) * 0.7
    chi = np.clip(1.5 - np.linspace(0, 2, n), 0, 1); q = np.linspace(0, 3, n); kf = 0.3 * rng.normal(size=n); cs = np.cos(rng.normal(size=n))
    psi1, psi2, psi3, h2, f3 = (rng.normal(size=n) for _ in range(5)); f1, f2 = chi * psi1, chi * psi2; fs = [f1, f2, f3]
    lamt = -0.2561532; alpha = 10 ** rng.uniform(-3, -1)
    Q1m = rng.normal(size=(n, n)) * 0.5; Q2m = rng.normal(size=(n, n)) * 0.5
    def F(u): return L @ u + (Q1m @ u) * (Dm @ u) + (Q2m @ u) ** 2 * u            # DF(0) = L
    def DF(u, v): return L @ v + (Q1m @ v) * (Dm @ u) + (Q1m @ u) * (Dm @ v) + 2 * (Q2m @ u) * (Q2m @ v) * u + (Q2m @ u) ** 2 * v
    ag = lambda u: 1 - kf * u
    A1 = lambda u: -(u - q * (Dm @ u) / ag(u)); B1 = lambda u: (Dm @ u) * cs / ag(u)
    u0 = alpha * psi3 + alpha ** 2 * h2
    # w with <chi w, f_j> = 0
    w = 0.05 * rng.normal(size=n)
    Cm = np.array([rho * chi * f for f in fs]); w = w - np.linalg.pinv(Cm) @ (Cm @ w)
    wt = chi * w; assert max(abs(ip(wt, f)) for f in fs) < 1e-13
    u = u0 + w
    PhiF = F(u0) - lamt * alpha * (psi3 + 2 * alpha * h2); Nw = F(u) - F(u0) - L @ w; Lin = DF(u0, w) - L @ w; Nq = Nw - Lin
    dirs = [psi1 - A1(u), psi2 - B1(u), psi3 + 2 * alpha * h2]
    Mat = np.array([[ip(chi * dirs[k_], fs[j]) for j in range(3)] for k_ in range(3)])
    G = np.array([[ip(chi * [psi1, psi2, psi3][k_], fs[j]) for j in range(3)] for k_ in range(3)]); E = Mat - G
    R = np.array([ip(chi * (L @ w + PhiF + Nw), f) for f in fs])
    x = np.linalg.solve(Mat.T, R); a, b, beta = x                                   # x^T (G + E) = R
    H = np.linalg.inv(G).T
    rec('T1_x_eq_HR_minus_HETx', np.max(np.abs(x - (H @ R - H @ E.T @ x))) / np.max(np.abs(x)))
    # E as in the lemma
    Elem = np.array([[-ip(chi * A1(u), fs[j]) for j in range(3)], [-ip(chi * B1(u), fs[j]) for j in range(3)], [2 * alpha * ip(chi * h2, fs[j]) for j in range(3)]])
    rec('T1b_E_matrix_definition', np.max(np.abs(E - Elem)) / np.max(np.abs(Elem)))
    # adjoint in <,>_rho: L* = rho^-1 L^T rho
    Ls = (L.T * rho[None, :]) / rho[:, None]
    comm = L @ (chi * w) - chi * (L @ w)                                           # [L, chi] w
    Rj = [ip(wt, Ls @ f) - ip(comm, f) + ip(chi * PhiF, f) + ip(chi * Nw, f) for f in fs]
    rec('T1c_R_j_form', np.max(np.abs(np.array(Rj) - R)) / np.max(np.abs(R)))
    u3 = Ls @ f3 - lamt * f3
    b3 = H[2, 2] * ip(wt, u3); b12 = H[2, 0] * ip(wt, Ls @ f1) + H[2, 1] * ip(wt, Ls @ f2)
    bLin = sum(H[2, j] * ip(chi * Lin, fs[j]) for j in range(3)); bQ = sum(H[2, j] * ip(chi * Nq, fs[j]) for j in range(3))
    bF = sum(H[2, j] * ip(chi * PhiF, fs[j]) for j in range(3)); bS = -sum(H[2, j] * ip(comm, fs[j]) for j in range(3))
    bE = -sum(H[2, j] * sum(E[k_, j] * x[k_] for k_ in range(3)) for j in range(3))
    rec('T2_beta_decomposition', abs(beta - (b3 + b12 + bLin + bQ + bF + bS + bE)) / abs(beta))
    # E0 pairing
    wtau = F(u) - (beta + lamt * alpha) * (psi3 + 2 * alpha * h2) - a * (psi1 - A1(u)) - b * (psi2 - B1(u))
    lhs = ip(chi * wtau, wt)
    g0 = chi * psi3
    Mterm = -beta * ip(g0, wt) - 2 * alpha * beta * ip(chi * h2, wt) + a * ip(chi * A1(u), wt) + b * ip(chi * B1(u), wt)
    rhs = ip(L @ wt, wt) - ip(comm, wt) + Mterm + ip(chi * PhiF, wt) + ip(chi * Nw, wt)
    rec('T3_E0_pairing_identity', abs(lhs - rhs) / (abs(lhs) + 1e-300))
    rec('T3b_constraints_conserved', max(abs(ip(chi * wtau, f)) for f in fs) / max(abs(R)))
    # projection onto {f1,f2,f3}^perp
    Fm = np.array(fs).T; Gf = np.array([[ip(fi, fj) for fj in fs] for fi in fs])
    Pi = lambda v: v - Fm @ np.linalg.solve(Gf, np.array([ip(f, v) for f in fs]))
    rec('T4_projection', abs(ip(g0, wt) - ip(Pi(g0), wt)) / abs(ip(g0, wt)))
    h = sum(H[2, j] * fs[j] for j in range(3)); rec('T5a_g0_h_eq_1', abs(ip(g0, h) - 1))
    bF0 = ip(chi * PhiF, h)
    rec('T5b_oblique', abs((ip(chi * PhiF, wt) - bF0 * ip(g0, wt)) - ip(chi * PhiF, wt - h * ip(g0, wt))) / abs(ip(chi * PhiF, wt)))
    nu, ng = nrm(Pi(u3)), nrm(Pi(g0)); bound = 0.5 * abs(H[2, 2]) * (nu * ng + abs(ip(Pi(u3), Pi(g0)))) * ip(wt, wt)
    rec('T6_rank_one_ratio_le_1', abs(b3 * ip(g0, wt)) / bound)
    # first order of E
    gA = chi * (psi3 - q * (Dm @ psi3)); gB = chi * (Dm @ psi3) * cs
    def E_u0(al): 
        uu = al * psi3 + al ** 2 * h2; return -ip(chi * A1(uu), f3), -ip(chi * B1(uu), f3)
    e13, e23 = ip(gA, f3), -ip(gB, f3)
    r_ = [abs(E_u0(al)[0] - al * e13) / al ** 2 for al in (1e-3, 1e-4)] + [abs(E_u0(al)[1] - al * e23) / al ** 2 for al in (1e-3, 1e-4)]
    rec('T7_E13_E23_first_order_remainder_ratio_dev', max(abs(r_[0] / r_[1] - 1), abs(r_[2] / r_[3] - 1)))          # O(alpha^2): the scaled remainders agree
    rec('T7b_E33_exact', abs(E[2, 2] - 2 * alpha * ip(chi * h2, f3)) / abs(E[2, 2]))
    # rho-symmetric D_rho
    S = rng.normal(size=(n, n)); Dr = (S + (S.T * rho[None, :]) / rho[:, None]) / 2
    a0, b0 = rng.normal(), rng.normal()
    v1 = -a0 * ip(f1, Dr @ wt) - b0 * ip(f2, Dr @ wt); v2 = -a0 * ip(Dr @ f1, wt) - b0 * ip(Dr @ f2, wt)
    rec('T8_Q2_lead_selfadjoint', abs(v1 - v2) / abs(v1))
res['worst_over_40_seeds'] = worst
tol = dict(T6_rank_one_ratio_le_1=1.0 + 1e-12, T7_E13_E23_first_order_remainder_ratio_dev=0.05)
res['pass'] = {k_: bool(v <= tol.get(k_, 1e-8)) for k_, v in worst.items()}
res['n_fail'] = sum(1 for v in res['pass'].values() if not v)
res['args'] = dict(argv=sys.argv, env={k_: os.environ.get(k_) for k_ in ('PYTHONHASHSEED', 'PYTHONDONTWRITEBYTECODE')}, numpy=np.__version__,
                   script_sha256=hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest())
json.dump(res, open(os.path.join(HERE, 'd_beta_float.json'), 'w'), indent=1, sort_keys=True)
print(json.dumps(res['worst_over_40_seeds'], indent=1)); print(res['pass']); print('n_fail', res['n_fail'])
