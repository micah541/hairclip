"""Final certificate for the constrained linear estimate (step 3), constraint functions in final form:
    f1 = chi * p
    f2 = chi * (C_tip cos K + S_tip sin K) = chi * (-sin theta)            [exact: C_tip, S_tip are the true tip values]
    f3 = chi * w~  - a eta1 - b eta3,   a = c ws1(tip),  b = c^3 (ws3 + ws1)(tip) + (24/LB^2) c ws1(tip),   c = (8 - st)/T0
         so that f3'(0) = f3'''(0) = 0 and the even extension of f3 is C^3 (f3 in D(L*)).
    eta1 = LB s zeta(s), eta3 = LB^3 s^3 zeta(s)/6, zeta = (1 - s^2)^4, s = X/(c LB) in [0,1] (0 beyond), LB = 2.
All three are smooth even functions on the curve.  G[i][j] = <A^{-1} f_j, f_i>_rho on the half line, A = -H_sym - 3/64.
Uses the doubleton outputs out8*.txt (run8.py) and the shared-symbol affine arithmetic of gram_post7."""
import sys, json
from flint import arb
import gen_hs4  # noqa: must precede gram_post7 (patched aux_for)
sys.path.insert(0, '/home/claude/s1spec/step3')
import gram_post7 as P
from gram_post7 import Aff, DSet, isect, na, SHARED_ALL, SHARED_Z
from gen_hs3 import QSHIFT, T0, LB

BASE = ["wa", "wb", "w"]; EB = ["e1", "e3"]; ALL = BASE + EB

def certify(prefix="out8", verbose=True):
    # ---- l_i = (0, m3_i/W13, 0, -m1_i/W13, 1) from the plane runs
    rho1, rho3 = {}, {}
    for gi in ALL:
        S = DSet(f"{prefix}Z{gi}.txt", SHARED_ALL + SHARED_Z); o = 5 + na([gi])
        rho1[gi] = S.ratio(o + 7, o + 4); rho3[gi] = S.ratio(o + 9, o + 4)
    def entry(gj, gi):
        """<A^{-1} f_gj, f_gi> = l_gi . (y_gj, I_gj,gi)"""
        S = DSet(f"{prefix}Q{gj}{gi}.txt", SHARED_ALL); o = 5 + na([gj, gi])
        c1, c3 = rho1[gi].c, rho3[gi].c
        main = S.form({o + 1: c3, o + 3: -c1, o + 4: arb(1)})
        d1 = Aff(arb(0), rho1[gi].d, rho1[gi].e); d3 = Aff(arb(0), rho3[gi].d, rho3[gi].e)
        y1 = S.form({o + 1: arb(1)}); y3 = S.form({o + 3: arb(1)})
        return (main + d3 * y1 - d1 * y3).scale(arb(2) ** QSHIFT)
    E = {}                                       # E[(i, j)] = <A^{-1} f_j, f_i>  (computed versions)
    for gi in BASE:
        for gj in BASE: E[(gi, gj)] = entry(gj, gi)
    for ea in EB:
        for gi in ALL: E[(gi, ea)] = entry(ea, gi)
    # ---- exact coefficients as affine forms in the shared symbols
    R = []; AB = {}
    for g in BASE:
        Zg = DSet(f"{prefix}Z{g}.txt", SHARED_ALL + SHARED_Z)
        st = Zg.form({4: arb(1)}); g1 = Zg.form({6: arb(1)}); g3 = Zg.form({8: arb(1)})
        cc = (Aff(arb(8)) - st).scale(1 / arb(T0))
        a = cc * g1
        b = cc * cc * cc * (g3 + g1) + (cc * g1).scale(arb(24) / arb(LB) ** 2)
        R.append({g: Aff(arb(1)), "e1": -a, "e3": -b}); AB[g] = (a.interval(), b.interval())
    Ct = St = Aff(arb(0))
    def G3entry(A_, B_, flip):
        acc = Aff(arb(0))
        for i, ri in R[A_].items():
            for j, rj in R[B_].items():
                key = (j, i) if flip else (i, j)
                if key not in E: key = (key[1], key[0])
                acc = acc + ri * rj * E[key]
        return acc
    G3 = [[isect(G3entry(A_, B_, False).interval(), G3entry(A_, B_, True).interval()) for B_ in range(3)] for A_ in range(3)]
    G3 = [[isect(G3[x][y], G3[y][x]) for y in range(3)] for x in range(3)]
    # affine (shared-symbol) Schur test: bound only at the end
    best = None
    for flip in (False, True):
        g = [[G3entry(A_, B_, flip) for B_ in range(3)] for A_ in range(3)]
        sym = lambda x, y: (g[x][y] + g[y][x]).scale(arb(1) / 2)
        a11, a12, a22, a13, a23, a33 = sym(0, 0), sym(0, 1), sym(1, 1), sym(0, 2), sym(1, 2), sym(2, 2)
        dA = a11 * a22 - a12 * a12
        NA = a22 * a13 * a13 - (a12 * a13 * a23).scale(arb(2)) + a11 * a23 * a23
        QA = a33 * dA - NA                         # S * det(G2); S < 0  <=>  QA < 0 when det(G2) > 0
        cand = (dA.interval(), QA.interval())
        if best is None or float(cand[1].upper()) < float(best[1].upper()): best = cand
    det2, Q = best
    g11 = G3[0][0]
    Sch = Q / det2
    if verbose: print("affine: det(G2) in", det2, "  S*det(G2) in", Q)
    ok = bool(g11 < 0 and det2 > 0 and Sch < 0)
    if verbose:
        for g in BASE: print(f"BC correction {g}: a =", AB[g][0], " b =", AB[g][1])
        for e in EB:
            print(f"  <A^-1 {e}, f>:", {gi: f"{float(E[(gi, e)].interval().mid()):+.4e}" for gi in ALL})
        for x in range(3): print("  G3:", [f"{float(G3[x][y].mid()):+.8e}(+-{float(G3[x][y].rad()):.1e})" for y in range(3)])
        print("G2: g11 =", g11, "  det(G2) =", det2)
        print(f"Schur complement  S = g33 - g^T G2^-1 g  in  [{float(Sch.lower()):.5e}, {float(Sch.upper()):.5e}]")
        print("G3 NEGATIVE DEFINITE:", ok)
    return dict(G3=G3, S=Sch, ok=ok, AB=AB)

if __name__ == "__main__":
    r = certify(sys.argv[1] if len(sys.argv) > 1 else "out9")
    json.dump({"ok": r["ok"], "S": [float(r["S"].lower()), float(r["S"].upper())],
               "G3_mid": [[float(v.mid()) for v in row] for row in r["G3"]]}, open("certificate_spectral_omega_aff.json", "w"), indent=1)
