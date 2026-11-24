"""From the rigorous 8x8 Gram (gram8*.json) to rigorous bounds of the two modulation forms on {f1,f2}^perp.
  y_a = <w, Pi v_a>,  ell = (<w, Pi Lf1>, <w, Pi Lf2>, <w, Pi u3>, <w, Pi u4>),  u = (gA, gB) [Q1] or (Df1, Df2) [Q2].
  ha = Minv11 Lf1 + Minv21 Lf2, hb = Minv12 Lf1 + Minv22 Lf2 ;  Q = y^T C y,  y = T ell,  Q = ell^T C' ell,  C' = T^T C T.
  Image of the unit ball of F^perp under w -> ell is the ellipsoid G^{1/2} B (G = Gram of the Pi v_a = Schur complement) hence
        sup Q/|w|^2 = lambda_max(C' G) = lambda_max(R C' R^T)  (G = R^T R).
  Rigour w.r.t. the interval enclosure [G] of the true Gram: for every true G in [G] (PSD) we have G <= G+ := mid(G) + diag(sum_j rad_ij)
  (diagonal dominance), hence E(G) is inside E(G+) and sup <= lambda_max(C' G+); symmetrically G- = mid - diag(rad) <= G gives
  inf Q <= lambda_min(C' G-) when G- is positive definite (used only for the upper end of the L2min enclosure).
  lambda_max of the 4x4 symmetric interval matrix B = R C' R^T is bounded by proving c I - mid(B) - diag(rad row sums) positive definite
  with an arb Cholesky.  (No regularisation of a near singular matrix is needed: the Loewner majorant G+ is positive definite by
  construction and differs from G by ~1e-6 relative.)"""
import sys, json
import numpy as np
from flint import arb, arb_mat, ctx
ctx.prec = 200
def ball(lo, hi):
    m = (lo + hi) / 2; r = (arb(hi) - arb(lo)) / 2
    return arb(m, float(r.upper()) * (1 + 1e-14) + 1e-300) if True else None
def up(x): return arb(x.upper())
def chol_pd(Mat, n):
    """arb Cholesky of a symmetric ball matrix (list of lists); returns R (upper, R^T R = Mat) or None if some pivot is not provably > 0"""
    R = [[arb(0)] * n for _ in range(n)]
    for i in range(n):
        s = Mat[i][i] - sum((R[k][i] * R[k][i] for k in range(i)), arb(0))
        if not (s.lower() > 0): return None
        R[i][i] = s.sqrt()
        for j in range(i + 1, n):
            t = Mat[i][j] - sum((R[k][i] * R[k][j] for k in range(i)), arb(0))
            R[i][j] = t / R[i][i]
    return R
def rad_ub(x): return arb(x.rad()).upper() if True else None
def majorant(Gm, n, sign=+1):
    """mid +- diag(sum_j rad_ij) as ball matrix with exact entries (radius of the result ~ 0)"""
    out = [[None] * n for _ in range(n)]
    for i in range(n):
        D = sum((arb(Gm[i][j].rad()).upper() for j in range(n)), arb(0))
        for j in range(n):
            mid = arb(Gm[i][j].mid())
            out[i][j] = mid + (sign * D if i == j else arb(0))
    return out
def eig_bound_max(B, n, c0):
    """smallest c >= c0 (found by doubling the gap) with c I - mid(B) - diag(row rad) provably PD; B symmetric interval matrix"""
    Dr = [sum((arb(B[i][j].rad()).upper() for j in range(n)), arb(0)) for i in range(n)]
    gap = max(abs(c0) * 1e-12, 1e-12); lo_ok = None
    for _ in range(200):
        c = arb(c0 + gap)
        Mt = [[(c if i == j else arb(0)) - arb(B[i][j].mid()) - (Dr[i] if i == j else arb(0)) for j in range(n)] for i in range(n)]
        if chol_pd(Mt, n) is not None: return c0 + gap
        gap *= 1.5
    raise RuntimeError("no bound")
def forms(G8, which):
    """G8: 8x8 arb matrix list (V = Lf1, Lf2, gA, gB, Df1, Df2; F = f1, f2).  which = 'Q1' or 'Q2'"""
    names = ["Lf1", "Lf2", "gA", "gB", "Df1", "Df2", "f1", "f2"]
    idx = [0, 1, 2, 3] if which == 'Q1' else [0, 1, 4, 5]
    F = [6, 7]
    M = [[G8[i][j] for j in F] for i in F]
    det = M[0][0] * M[1][1] - M[0][1] * M[1][0]; assert det.lower() > 0
    Mi = [[M[1][1] / det, -M[0][1] / det], [-M[1][0] / det, M[0][0] / det]]
    # Schur complement on the 4 chosen vectors
    G = [[G8[a][b] - sum((G8[a][F[i]] * Mi[i][j] * G8[F[j]][b] for i in range(2) for j in range(2)), arb(0)) for b in idx] for a in idx]
    for i in range(4):
        for j in range(i): G[i][j] = G[j][i]
    # C' = T^T C T
    T = [[Mi[0][0], Mi[1][0], arb(0), arb(0)], [Mi[0][1], Mi[1][1], arb(0), arb(0)], [arb(0), arb(0), arb(1), arb(0)], [arb(0), arb(0), arb(0), arb(1)]]
    C = [[arb(0)] * 4 for _ in range(4)]
    C[0][2] = C[2][0] = arb(-1) / 2; s24 = arb(1) / 2 if which == 'Q1' else arb(-1) / 2; C[1][3] = C[3][1] = s24
    TC = [[sum((T[k][i] * C[k][l] for k in range(4)), arb(0)) for l in range(4)] for i in range(4)]
    Cp = [[sum((TC[i][l] * T[l][j] for l in range(4)), arb(0)) for j in range(4)] for i in range(4)]
    for i in range(4):
        for j in range(i): Cp[i][j] = (Cp[i][j] + Cp[j][i]) / 2 if False else Cp[i][j].union(Cp[j][i])
    res = {"Minv": Mi, "G4": G, "Cp": Cp}
    # float reference
    Gf = np.array([[float(G[i][j].mid()) for j in range(4)] for i in range(4)]); Cf = np.array([[float(Cp[i][j].mid()) for j in range(4)] for i in range(4)])
    ev = np.sort(np.linalg.eigvals(Cf @ Gf).real); res["float_eigs"] = ev.tolist()
    nrm = np.sqrt(np.diag(Gf)); res["G4_normalised_eigs"] = np.linalg.eigvalsh(Gf / np.outer(nrm, nrm)).tolist()
    def sup_bound(Cmat):
        """rigorous upper bound of sup ell^T Cmat ell over E(G)"""
        Gp = majorant(G, 4, +1); R = chol_pd(Gp, 4); assert R is not None, "G+ not PD"
        # B = R Cmat R^T
        RC = [[sum((R[i][k] * Cmat[k][l] for k in range(4)), arb(0)) for l in range(4)] for i in range(4)]
        B = [[sum((RC[i][l] * R[j][l] for l in range(4)), arb(0)) for j in range(4)] for i in range(4)]
        for i in range(4):
            for j in range(i): B[i][j] = B[i][j].union(B[j][i]); B[j][i] = B[i][j]
        Bf = np.array([[float(B[i][j].mid()) for j in range(4)] for i in range(4)])
        lam = np.linalg.eigvalsh((Bf + Bf.T) / 2).max()
        return eig_bound_max(B, 4, lam)
    neg = [[-x for x in row] for row in Cp]
    res["lam_max_upper"] = sup_bound(Cp)
    res["lam_min_lower"] = -sup_bound(neg)
    # inner enclosures: Rayleigh quotients of explicit test vectors w = sum_a beta_a Pi v_a:
    #   Q(w)/|w|^2 = beta^T G C' G beta / beta^T G beta, evaluated in ball arithmetic with G in [G]  (rigorous for every true G in [G])
    from scipy.linalg import eigh
    Cf_ = np.array([[float(Cp[i][j].mid()) for j in range(4)] for i in range(4)])
    wv, V = eigh(Gf @ Cf_ @ Gf, Gf + 1e-14 * np.trace(Gf) * np.eye(4))
    def rq(beta):
        b = [arb(float(x)) for x in beta]
        Gb = [sum((G[i][j] * b[j] for j in range(4)), arb(0)) for i in range(4)]
        den = sum((b[i] * Gb[i] for i in range(4)), arb(0))
        num = sum((Gb[i] * Cp[i][j] * Gb[j] for i in range(4) for j in range(4)), arb(0))
        r = num / den
        return float(r.lower()), float(r.upper())
    res["rq_min"] = rq(V[:, 0]); res["rq_max"] = rq(V[:, -1])
    res["lam_min_upper"] = res["rq_min"][1]; res["lam_max_lower"] = res["rq_max"][0]
    return res
def load(path):
    R = json.load(open(path))["G"]
    return [[ball(R[i][j][0], R[i][j][1]) for j in range(8)] for i in range(8)]
if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else 'gram8_g3.json'
    G8 = load(path); out = {}
    for which in ("Q1", "Q2"):
        r = forms(G8, which)
        out[which] = {"lam_max_upper": float(r["lam_max_upper"]), "lam_min_lower": float(r["lam_min_lower"]), "lam_min_upper": r["lam_min_upper"], "lam_max_lower": r["lam_max_lower"], "rq_min": r["rq_min"], "rq_max": r["rq_max"],
                      "float_eigs": r["float_eigs"], "G4_normalised_eigs": r["G4_normalised_eigs"],
                      "Minv": [[float(x.mid()) for x in row] for row in r["Minv"]], "Minv_rad": float(max(x.rad() for row in r["Minv"] for x in row))}
        print(which, json.dumps(out[which], indent=1))
    print("R1_upper =", max(out["Q1"]["lam_max_upper"], -out["Q1"]["lam_min_lower"]), " (lower witnesses: sup>=", out["Q1"]["lam_max_lower"], " inf<=", out["Q1"]["lam_min_upper"], ")")
    print("L2max_upper =", out["Q2"]["lam_max_upper"], " L2min in [", out["Q2"]["lam_min_lower"], ",", out["Q2"]["lam_min_upper"], "]")
    json.dump(out, open('forms_' + path, 'w'), indent=1)
