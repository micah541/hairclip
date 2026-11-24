"""Rigorous alpha-level forms on {f1,f2}^perp and {f1,f2,f3}^perp, including the beta-loop term.
 y = (<w,ha>, <w,hb>, <w,gA>, <w,gB>, <w,g0>), g0 = chi psi3^.
 Signs: x^T(G+E) = R, x = H R - H E^T x. alpha-linear beta = -alpha H33 (e13 a + e23 b) [+ H31,H32 terms, negligible, bounded below],
 e13 = <gA,f3> (E_13 = -<chi A1(u),f3>, A1(psi3) = -(psi - q psi')), e23 = -<gB,f3> (E_23 = -<chi B1,f3>, B1 = psi' cos).
 Energy term -beta <g0,w> =  + alpha H33 (e13 <w,ha> + e23 <w,hb>) <w,g0>  (a = <w,ha>, b = <w,hb> to leading order, H_1j = Minv_j1)."""
import sys, json
import numpy as np
from flint import arb, ctx
ctx.prec = 200
sys.path.insert(0, '/home/claude/s1spec/step3nl/agent_mod')
import mod_forms as MF
ball = MF.ball
R7 = json.load(open('gram7_g3.json'))['G']          # Lf1 Lf2 gA gB g0 f1 f2
P = json.load(open('pair_f3.json')); LY = json.load(open('/home/claude/s1spec/step3/lyap_consts2.json')); LX = json.load(open('/home/claude/s1spec/step3nl/agent_mod/lemma/lyap_extra.json'))
REM = json.load(open('rem.json'))
iv = lambda p: ball(p[0], p[1])
G = [[ball(*R7[i][j]) for j in range(7)] for i in range(7)]
pair = [iv(P[n]) for n in ("Lf1", "Lf2", "gA", "gB", "g0")]
f13 = arb(0, LY['psif3'][0]); f23 = arb(0, LY['psif3'][1])          # |<f1,f3>|, |<f2,f3>| (lyap_consts2 psif3)
nlo, nhi = LX['nf3_lo'], LX['nf3_hi']; n33 = ball(nlo * nlo, nhi * nhi)
# 8x8
G8 = [[None] * 8 for _ in range(8)]
for i in range(7):
    for j in range(7): G8[i][j] = G[i][j]
for i in range(5): G8[i][7] = G8[7][i] = pair[i]
G8[5][7] = G8[7][5] = f13; G8[6][7] = G8[7][6] = f23; G8[7][7] = n33
def schur(Gm, idx, cons):
    """project idx vectors off the constraint vectors cons (sequentially, 1 at a time with ball arithmetic)"""
    Gm = [row[:] for row in Gm]; cons = list(cons)
    for c in cons:
        for i in range(len(Gm)):
            for j in range(len(Gm)):
                if i != c and j != c: pass
        n = len(Gm)
        new = [[Gm[i][j] - Gm[i][c] * Gm[c][j] / Gm[c][c] for j in range(n)] for i in range(n)]
        Gm = new
    return [[Gm[i][j] for j in idx] for i in idx]
def schur2(Gm, idx, c1, c2):
    # exact 2x2 block elimination for (f1,f2) = (5,6)
    M = [[Gm[c1][c1], Gm[c1][c2]], [Gm[c2][c1], Gm[c2][c2]]]; det = M[0][0] * M[1][1] - M[0][1] * M[1][0]
    Mi = [[M[1][1] / det, -M[0][1] / det], [-M[1][0] / det, M[0][0] / det]]; cs = [c1, c2]; n = len(Gm)
    return [[Gm[i][j] - sum((Gm[i][cs[a]] * Mi[a][b] * Gm[cs[b]][j] for a in range(2) for b in range(2)), arb(0)) for j in range(n)] for i in range(n)], Mi
G2, Mi = schur2(G8, range(8), 5, 6)                      # projected off f1,f2 (rows/cols 5,6 become 0)
G3p = [[G2[i][j] - G2[i][7] * G2[7][j] / G2[7][7] for j in range(8)] for i in range(8)]
G5_2 = [[G2[i][j] for j in range(5)] for i in range(5)]; G5_3 = [[G3p[i][j] for j in range(5)] for i in range(5)]
for Gm in (G5_2, G5_3):
    for i in range(5):
        for j in range(i): Gm[i][j] = Gm[j][i]
# H33 etc. from the lemma block inverse: G33 = <chi psi3, f3>
g33 = pair[4]; H33 = 1 / (g33 - arb(0, abs(float(G8[4][5].upper())) * 0))              # leading; corrections from H31, H32 below
e13 = pair[2]; e23 = -pair[3]
T = [[Mi[0][0], Mi[1][0], 0, 0, 0], [Mi[0][1], Mi[1][1], 0, 0, 0], [0, 0, 1, 0, 0], [0, 0, 0, 1, 0], [0, 0, 0, 0, 1]]
T = [[arb(x) if not isinstance(x, arb) else x for x in r] for r in T]
def Cp_of(c13, c24, c15, c25):
    C = [[arb(0)] * 5 for _ in range(5)]
    C[0][2] = C[2][0] = arb(c13) / 2; C[1][3] = C[3][1] = arb(c24) / 2; C[0][4] = C[4][0] = c15 / 2; C[1][4] = C[4][1] = c25 / 2
    TC = [[sum((T[k][i] * C[k][l] for k in range(5)), arb(0)) for l in range(5)] for i in range(5)]
    Cp = [[sum((TC[i][l] * T[l][j] for l in range(5)), arb(0)) for j in range(5)] for i in range(5)]
    for i in range(5):
        for j in range(i): Cp[i][j] = Cp[i][j].union(Cp[j][i]); Cp[j][i] = Cp[i][j]
    return Cp
def bounds(Gm, Cp):
    n = 5
    Gp = MF.majorant(Gm, n, +1); R = MF.chol_pd(Gp, n); assert R is not None
    def sup(Cm):
        RC = [[sum((R[i][k] * Cm[k][l] for k in range(n)), arb(0)) for l in range(n)] for i in range(n)]
        B = [[sum((RC[i][l] * R[j][l] for l in range(n)), arb(0)) for j in range(n)] for i in range(n)]
        for i in range(n):
            for j in range(i): B[i][j] = B[i][j].union(B[j][i]); B[j][i] = B[i][j]
        Bf = np.array([[float(B[i][j].mid()) for j in range(n)] for i in range(n)])
        return MF.eig_bound_max(B, n, np.linalg.eigvalsh((Bf + Bf.T) / 2).max())
    neg = [[-x for x in r] for r in Cp]
    Gf = np.array([[float(Gm[i][j].mid()) for j in range(n)] for i in range(n)]); Cf = np.array([[float(Cp[i][j].mid()) for j in range(n)] for i in range(n)])
    return {"sup_upper": sup(Cp), "inf_lower": -sup(neg), "float_eigs": np.sort(np.linalg.eigvals(Cf @ Gf).real).tolist(),
            "G5_min_normalised_eig": float(np.linalg.eigvalsh(Gf / np.outer(np.sqrt(np.diag(Gf)), np.sqrt(np.diag(Gf)))).min())}
cH13 = H33 * e13; cH23 = H33 * e23
print("H33 =", H33, " c*e13 =", cH13, " c*e23 =", cH23)
out = {"H33": [float(H33.lower()), float(H33.upper())], "e13": [float(e13.lower()), float(e13.upper())], "e23": [float(e23.lower()), float(e23.upper())],
       "H33*e13": [float(cH13.lower()), float(cH13.upper())], "H33*e23": [float(cH23.lower()), float(cH23.upper())],
       "E33_over_alpha=2<chi h2,f3>": [2 * P['h2c'][0], 2 * P['h2c'][1]], "alpha2_remainders": REM, "pair_f3": P}
cfg = {"Q1": (-1, 1, arb(0), arb(0)), "extra": (0, 0, cH13, cH23), "Qfull": (-1, 1, cH13, cH23)}
def verts(lo_hi):
    return [arb(lo_hi[0]), arb(lo_hi[1])]
V13 = verts((float(cH13.lower()), float(cH13.upper()))); V23 = verts((float(cH23.lower()), float(cH23.upper())))
# the sup of lambda_max over a box of the (linear) parameters (c13, c23) is attained at a vertex (convexity); same for -inf.
for sp_name, Gm in (("f12perp", G5_2), ("f123perp", G5_3)):
    for nm in ("Q1", "extra", "Qfull"):
        sups, infs, fl = [], [], []
        vs = [(arb(0), arb(0))] if nm == "Q1" else [(a_, b_) for a_ in V13 for b_ in V23]
        for (c_, d_) in vs:
            a, b = (-1, 1) if nm != "extra" else (0, 0)
            if nm == "Q1": c_ = d_ = arb(0)
            r = bounds(Gm, Cp_of(a, b, c_, d_)); sups.append(float(r["sup_upper"])); infs.append(float(r["inf_lower"])); fl += r["float_eigs"]
        out[f"{nm}_{sp_name}"] = {"sup_upper": max(sups), "inf_lower": min(infs), "abs_upper": max(max(sups), -min(infs)), "float_eig_range_over_vertices": [min(fl), max(fl)]}
        print(nm, sp_name, "sup<=%.4f inf>=%.4f" % (max(sups), min(infs)))
# non-rigorous midpoint-F version for comparison
Pm = P['_mid_F_nonrigorous']; mid = lambda n: (Pm[n][0] + Pm[n][1]) / 2
cm13 = mid('g0') ** -1 * mid('gA'); cm23 = -mid('g0') ** -1 * mid('gB')
rr = bounds(G5_3, Cp_of(-1, 1, arb(cm13), arb(cm23))); out["Qfull_f123perp_midF_heuristic"] = {"c_e13": cm13, "c_e23": cm23, "float_eigs": rr["float_eigs"], "sup_upper": float(rr["sup_upper"]), "inf_lower": float(rr["inf_lower"])}
print("mid-F heuristic: ce13, ce23", cm13, cm23, "sup<=%.3f inf>=%.3f" % (rr["sup_upper"], rr["inf_lower"]))
json.dump(out, open('bloop.json', 'w'), indent=1)
